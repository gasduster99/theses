/******************************************************************************
File     : CCPSO.cpp
Author   : L. Shawn Matott
Copyright: 2009, L. Shawn Matott

Computation Constrained Particle Swarm Optimization (CCPSO) - PSO on a budget!

Version History
07-18-09    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "CCPSO.h"
#include "Model.h"
#include "Exception.h"
#include "WriteUtility.h"
#include "StatUtility.h"
#include "Utility.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"
#include "LevenbergAlgorithm.h"
#include "LatinHypercube.h"
#include "SuperMuseUtility.h"
#include "SuperMUSE.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
CCPSO::CCPSO(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pSwarm = NULL;
   m_pStats = NULL;
   m_pInit = NULL;
   m_SwarmSize = 0;
   m_MaxGens = 0;
   m_Budget = 0;
   m_BestIdx = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumInit = 0;

   //MPI-parallel communication arrays
   m_pMyBuf  = NULL;
   m_pTmpBuf = NULL;
   m_pBigBuf = NULL;
   m_pBuf    = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
DTOR

Free up memory used by PSO and it's member variables.
******************************************************************************/
CCPSO::~CCPSO(void)
{
   Destroy();
}/* end DTOR()*/

/******************************************************************************
Destroy()

Free up memory used by the PSO and it's member variables.
******************************************************************************/
void CCPSO::Destroy(void)
{
   int i;
   for(i = 0; i < m_SwarmSize; i++)
   {
      delete [] m_pSwarm[i].v;
      delete [] m_pSwarm[i].x;
      delete [] m_pSwarm[i].b;
   }
   delete [] m_pSwarm;

   for(i = 0; i < m_NumInit; i++)
   {
      delete [] m_pInit[i];
   }
   delete [] m_pInit;
   
   delete m_pStats;
   m_SwarmSize = 0;
   m_Budget = 0;
   m_BestIdx = 0;

   delete [] m_pMyBuf;
   delete [] m_pTmpBuf;
   delete [] m_pBigBuf;
   delete [] m_pBuf;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using PSO.
******************************************************************************/
void CCPSO::Calibrate(void)
{ 
   int id;
   char fileName[DEF_STR_SZ];
   FILE * pFile;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   
   Optimize();

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   if(id == 0)
   {
      sprintf(fileName, "OstOutput%d.txt", id);

      //write statistics of best parameter set to output file
      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      //write statistics of best parameter set to output file
      m_pStats->WriteStats(stdout);
   }/* end if() */
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using PSO.
******************************************************************************/
void CCPSO::Optimize(void)
{
   StatusStruct pStatus;
   int num, g, nleft;
   int id, lvl, idx;
   int i, j;
   double upr, lwr, sgn;
   double x, pl, pg, r1, r2, v, vmin;
   double rval, init; //initial inertia
   double c1, c2, fl, fg;
   LatinHypercube * pLHS = NULL;
   ParameterGroup * pGroup;
   
   InitFromFile(GetInFileName());

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   nleft = m_Budget;

   if(id == 0)
   {
      WriteSetup(m_pModel, "Computation Constrained Particle Swarm Optimization");
      //write banner
      WriteBanner(m_pModel, "gen   best value     ", "Pct. Complete");
   }/* end if() */

   //allocate swarm
   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();

   NEW_PRINT("ParticleStruct", m_SwarmSize);
   m_pSwarm = new ParticleStruct[m_SwarmSize];
   MEM_CHECK(m_pSwarm);

   for(i = 0; i < m_SwarmSize; i++)
   {
      NEW_PRINT("double", num);
      m_pSwarm[i].x = new double[num];

      NEW_PRINT("double", num);
      m_pSwarm[i].v = new double[num];

      NEW_PRINT("double", num);
      m_pSwarm[i].b = new double[num];

      m_pSwarm[i].n = num;      
   }/* end for() */
   MEM_CHECK(m_pSwarm[i-1].b);

   //initialize swarm
   NEW_PRINT("LatinHypercube", 1);
   pLHS = new LatinHypercube(num, m_SwarmSize);
   MEM_CHECK(pLHS);

   for(j = 0; j < num; j++)
   { 
      lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
      upr = pGroup->GetParamPtr(j)->GetUprBnd();
      pLHS->InitRow(j, lwr, upr);
   }/* end for() */
 
   lvl = idx = 0;
   for(i = 0; i < m_SwarmSize; i++) //for each particle
   {
      //initial velocity is 0.00
      for(j = 0; j < num; j++){ m_pSwarm[i].v[j] = 0.00;}

      //LHS_INIT
      for(j = 0; j < num; j++)
      { 
         rval = pLHS->SampleRow(j);
         m_pSwarm[i].x[j] = rval;
         m_pSwarm[i].b[j] = rval;
      }/* end for() */
   }/* end for() */

   //seed swarm with pre-specified values
   for(i = 0; (i < m_NumInit) && (i < m_SwarmSize); i++)
   {
      for(j = 0; j < num; j++)
      {         
         m_pSwarm[i].x[j] = m_pInit[i][j];
         m_pSwarm[i].b[j] = m_pInit[i][j];
      }/* end for() */
   }/* end for() */

   pLHS->Destroy();

   //evaluate swarm, possibly in parallel
   EvaluateSwarm();
   nleft -= m_SwarmSize;

   //perform intermediate bookkeeping
   m_pModel->Bookkeep(false);

   for(i = 0; i < m_SwarmSize; i++)
   {
	   m_pSwarm[i].fb = m_pSwarm[i].fx;
	   m_pSwarm[i].cb = m_pSwarm[i].cx;
   }

   /* --------------------------------------------
   enable special parameters now that local best 
   is initialized for each particle
   -------------------------------------------- */
   pGroup->EnableSpecialParams();

   //determine the best particle
   m_BestIdx = 0;
   m_Best = m_pSwarm[0].fb;
   for(i = 0; i < m_SwarmSize; i++)
   {
      if(m_Best > m_pSwarm[i].fx)
      {
         m_Best = m_pSwarm[i].fx;
         m_BestIdx = i;
      }/* end if() */
   }/* end if() */

   if(id == 0)
   {
      //write initial config.
      pStatus.curIter = 0;
      pStatus.maxIter = m_MaxGens;
      pStatus.pct = (float)100.00*((float)1.00-(float)nleft/(float)m_Budget);
      pStatus.numRuns = m_pModel->GetCounter();
      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);
      WriteRecord(m_pModel, 0, m_Best, pStatus.pct);

      WriteStatus(&pStatus);
   }/* end if() */

   //assign initial weights
   c1 = m_c1; c2 = m_c2;
   init = m_Inertia;

   //main optimization loop   
   for(g = 0; g < m_MaxGens; g++)
   {
      pStatus.curIter = m_CurGen = g+1;
      if(IsQuit() == true){ break;}
      if(nleft <= 0){ pStatus.pct = 100.00; break;}

      if(id == 0)
      {
         //update velocities, parameters and objective functions
         for(i = 0; i < m_SwarmSize; i++) //for each particle
         {            
            fl = (double)GetParticleRank(i);
            fg = (double)(m_SwarmSize - 1);
            //assign cognitive and social weights               
            m_c1 = 2.0 - (fl/fg);
            m_c2 = 4.00 - m_c1;

            //FILE * pOut = fopen("OstWeights.txt", "a");
            //fprintf(pOut, "Gen[%d] Particle[%d] C1=%0.2lf C2=%0.2lf  chi=%0.2lf\n", g, i, m_c1, m_c2, m_Constrict);
            //fclose(pOut);

            for(j = 0; j < num; j++) //for each parameter
            {
               //intermediary variables
               x = m_pSwarm[i].x[j];
               pl = m_pSwarm[i].b[j];
               pg = m_pSwarm[m_BestIdx].b[j];

               //random weights
               r1 = (double)MyRand() / (double)MY_RAND_MAX;
               r2 = (double)MyRand() / (double)MY_RAND_MAX;

               //revised velocity
               v = m_pSwarm[i].v[j];
               v = m_Constrict*((m_Inertia*v) + m_c1*r1*(pl-x) + m_c2*r2*(pg-x));

   			   //assign minimum perturbation to prevent stagnation
	   		   if(strcmp(pGroup->GetParamPtr(j)->GetType(), "real") == 0)
		   		   vmin = (0.01*fabs(x))/(g+1); 
			      else
				      vmin = 0.50;

               if(fabs(v) < vmin)
               {
                  //adjust randomized minimum velocity
                  sgn = (double)MyRand() / (double)MY_RAND_MAX; //random direction
                  if(sgn >= 0.50) v = +((1.00+r1)*vmin);
                  else            v = -((1.00+r2)*vmin);
               }
               m_pSwarm[i].v[j] = v;
            
               //revised position
               m_pSwarm[i].x[j] = x + v;
            }/* end for() */

            /*-----------------------------------------------------------
            Constrain revised position to stay within parameter limits
            but be sure to preserve angle (i.e. direction of movement)
            -----------------------------------------------------------*/
            double dx_min, dx_old, dx_new, dx_frac, frac;
            double L, G, low, hi;
            dx_min = 1.00;
            fl = (double)GetParticleRank(i);
            fg = (double)(m_SwarmSize)-1;
            frac=(2.0+(fl/fg))*(1.00-0.99*(double)g/(double)m_MaxGens); 

            for(j = 0; j < num; j++) //for each parameter
            {
               lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
               upr = pGroup->GetParamPtr(j)->GetUprBnd();

               /* ---------------------------------------------------
               Adaptively restrict parameter bounds.
               --------------------------------------------------- */
               x = m_pSwarm[i].x[j] - m_pSwarm[i].v[j];
               L = m_pSwarm[i].b[j];
               G = m_pSwarm[m_BestIdx].b[j];

               low = MyMin(MyMin(x,L),G);
               low -= frac*fabs(low);
               lwr = MyMax(low,lwr);

               hi = MyMax(MyMax(x,L),G);
               hi += frac*fabs(hi);
               upr = MyMin(hi, upr);

               v = m_pSwarm[i].v[j];
               x = m_pSwarm[i].x[j] - v; //original position
               //compute fractional move, store most restrictive fraction
               if(m_pSwarm[i].x[j] > upr)
               {
                  dx_old = v;
                  dx_new = 0.5*(upr-x);
                  dx_frac = fabs(dx_new/dx_old); //relative change
                  if(dx_frac < dx_min) dx_min = dx_frac;
                  m_NumUprViols++;
               }
               if(m_pSwarm[i].x[j] < lwr)
               {
                  dx_old = v;
                  dx_new = 0.5*(lwr-x);
                  dx_frac = fabs(dx_new/dx_old); //relative change
                  if(dx_frac < dx_min) dx_min = dx_frac;
                  m_NumLwrViols++;
               }
            }/* end for() */

            for(j = 0; j < num; j++) //for each parameter
            {
               v = m_pSwarm[i].v[j];
               x = m_pSwarm[i].x[j] - v; //original position      
               m_pSwarm[i].v[j] *= dx_min; //revised velocity
               m_pSwarm[i].x[j] = x + (v*dx_min); //revised position
            }
         }/* end for() */
      } /* end if() */

      //evaluate swarm, possibly in parallel
      EvaluateSwarm();
      nleft -= m_SwarmSize;

      //adjust weights
      m_RedRate = (double)g/(double)m_MaxGens;
      m_Inertia = init*(1.00 - m_RedRate);
   
      //revise avg, and local and global best
      for(i = 0; i < m_SwarmSize; i++)
      {
         //revise local best
         if(m_pSwarm[i].fx < m_pSwarm[i].fb)
         {
            for(j = 0; j < num; j++){ m_pSwarm[i].b[j] = m_pSwarm[i].x[j];}
            m_pSwarm[i].fb = m_pSwarm[i].fx;
			   m_pSwarm[i].cb = m_pSwarm[i].cx;
         }/* end if() */
         //revise global best         
         if(m_Best > m_pSwarm[i].fx)
         {
            m_Best = m_pSwarm[i].fx;
            m_BestIdx = i;
         }/* end if() */
      }/* end for() */

      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);

      pStatus.pct = ((float)100.00*(float)(g+1))/(float)m_MaxGens;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      if(id == 0){ WriteRecord(m_pModel, (g+1), m_Best, pStatus.pct);}

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */

   m_Constrict = 1.00;
   m_Inertia = init; //reset weights
   m_c1 = c1;
   m_c2 = c2;

   //place model at optimal prameter set
   pGroup->WriteParams(m_pSwarm[m_BestIdx].b);
   m_pModel->Execute();

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(id == 0)
   { 
      WriteOptimal(m_pModel, m_Best);
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      //write algorithm metrics
      WriteAlgMetrics(this);
   }
} /* end Optimize() */

/******************************************************************************
GetParticleRank()

Rank the particle based on local best.
******************************************************************************/
int CCPSO::GetParticleRank(int p)
{
   int count = 0;
   for(int i = 0; i < m_SwarmSize; i++)
   {
      if((i != p) && (m_pSwarm[i].fb < m_pSwarm[p].fb))
      {
         count++;
      }
   }/* end for() */
   return count;
}/* end GetParticleRank() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void CCPSO::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Computation Constrained Particle Swarm Optimization\n");
   fprintf(pFile, "Budget                  : %d\n", m_Budget);
   fprintf(pFile, "Max Generations         : %d\n", m_MaxGens);
   fprintf(pFile, "Actual Generations      : %d\n", m_CurGen);
   fprintf(pFile, "Swarm Size              : %d\n", m_SwarmSize);
   fprintf(pFile, "Constriction Factor     : %.2lf\n", m_Constrict);  
   fprintf(pFile, "Init. Cognitive Weight  : %.2lf\n", m_c1);
   fprintf(pFile, "Init. Social Weight     : %.2lf\n", m_c2);
   fprintf(pFile, "Inertia Weight          : %.2lf\n", m_Inertia);
   
   fprintf(pFile, "Inertia Reduction Rate  : ");
   fprintf(pFile, "Linear reduction to zero\n");

   fprintf(pFile, "Cognition Reduction Rate  : ");
   fprintf(pFile, "Linear reduction to 1.00\n");

   fprintf(pFile, "Social Awareness Rate   : ");
   fprintf(pFile, "Linear increase to 3.00\n");

   fprintf(pFile, "Initialization Method   : ");
   if(m_InitType == RANDOM_INIT){ fprintf(pFile, "Random\n");}
   else if(m_InitType == QUAD_TREE_INIT){ fprintf(pFile, "Quad-Tree\n");}
   else if(m_InitType == LHS_INIT){ fprintf(pFile, "Latin Hypercube Sampling\n");}
   else { fprintf(pFile, "Unknown\n");}

   fprintf(pFile, "Upper Violations        : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations        : %d\n", m_NumLwrViols);

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
EvaluateSwarm()

Evaluates the objective function of each particle in the swarm.
******************************************************************************/
void CCPSO::EvaluateSwarm(void)
{   
   int i, n;   
   ParameterGroup * pGroup;
   double val;

   MPI_Comm_size(MPI_COMM_WORLD, &n);
   
   if(n == 1) //serial execution
   {
      if(IsSuperMUSE() == false)
      {
         WriteInnerEval(WRITE_PSO, m_SwarmSize, '.');
         pGroup = m_pModel->GetParamGroupPtr();
         for(i = 0; i < m_SwarmSize; i++) 
         { 
            WriteInnerEval(i+1, m_SwarmSize, '.');
            pGroup->WriteParams(m_pSwarm[i].x);

			//let special parameters know about local best
			pGroup->ConfigureSpecialParams(m_pSwarm[i].fb, m_pSwarm[i].cb);

            val = m_pModel->Execute();
            m_pSwarm[i].fx = val;
			   m_pSwarm[i].cx = pGroup->GetSpecialConstraint();
         }
         WriteInnerEval(WRITE_ENDED, m_SwarmSize, '.');
      }
      else
      {
         EvalSwarmSuperMUSE();
      }
   }/* end if() */
   else /* parallel execution */
   {
      BcastSwarm();
      EvalSwarmParallel();      
   }/* end else() */
} /* end EvaluateSwarm() */

/******************************************************************************
BcastSwarm()

When in parallel, only the master computes the swarm movement. All the other 
processors just compute the objeective functions. The BcastSwarm() routine is 
called upon to broadcast  the current particle swarm from the master processor 
to all of the slave processors.
******************************************************************************/
void CCPSO::BcastSwarm(void)
{
   int num_vars, pop_size, buf_size;
   int i, j, num_procs, id, idx;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix
   pop_size = m_SwarmSize;
   num_vars = m_pSwarm[0].n;

   buf_size = pop_size*num_vars;

   if(m_pBuf == NULL)
   {
      NEW_PRINT("double", buf_size);
      m_pBuf = new double[buf_size];
      MEM_CHECK(m_pBuf);
   }

   for(i = 0; i < buf_size; i++){ m_pBuf[i] = 999.99;}

   //fill up the flattened matrix
   for(i = 0; i < num_vars; i++)
   {
      for(j = 0; j < pop_size; j++)
      {
         idx = (num_vars)*j + i;
         m_pBuf[idx] = m_pSwarm[j].x[i];
      }/* end for() */
   }/* end for() */

   //broadcast the flattened matrix
   MPI_Bcast(m_pBuf, buf_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   //use the flattened matrix to fill swarm
   for(i = 0; i < num_vars; i++)
   {
      for(j = 0; j < pop_size; j++)
      {
         
         idx = (num_vars)*j + i;
         m_pSwarm[j].x[i] = m_pBuf[idx];
      }/* end for() */
   }/* end for() */
}/* end BcastSwarm() */

/******************************************************************************
EvalSwarmParallel()

Compute objective function of entire particle swarm in parallel. Each processor 
evaluates a predetermined number of swarm particles, based on their processor id.
******************************************************************************/
void CCPSO::EvalSwarmParallel(void)
{    
   int i ,j, num_procs, id, bufsize, idx;
   ParameterGroup * pGroup;

   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
   bufsize = (m_SwarmSize/num_procs) + 1;

   //allocate space for intermediate buffers, if necessary
   if(m_pMyBuf == NULL)
   {
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];

      NEW_PRINT("double", bufsize);
      m_pTmpBuf = new double[bufsize];

      NEW_PRINT("double", m_SwarmSize);
      m_pBigBuf = new double[m_SwarmSize];
      MEM_CHECK(m_pBigBuf);
   }

   //perform parallel evaluations
   j = 0;
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_SwarmSize; i++) 
   { 
      if((i % num_procs) == id)
      {          
         pGroup->WriteParams(m_pSwarm[i].x);

		//let special parameters know about local best
		pGroup->ConfigureSpecialParams(m_pSwarm[i].fb, m_pSwarm[i].cb);

         m_pMyBuf[j] = m_pModel->Execute();
         m_pTmpBuf[j] = m_pMyBuf[j];
         j++;
      }/* end if() */
   }/* end for() */

   //gather results
   for(i = 0; i < num_procs; i++)
   {
      //receive someones buf, this will clobber myBuf
      MPI_Bcast(m_pMyBuf, bufsize, MPI_DOUBLE, i, MPI_COMM_WORLD);

      for(j = 0; j < bufsize; j++)
      {
         idx = (num_procs * j) + i; //idx maps myBuf into bigBuf

         if(idx < m_SwarmSize)
         {
            m_pBigBuf[idx] = m_pMyBuf[j]; //gather into bigbuf
            m_pMyBuf[j] = m_pTmpBuf[j]; //restore myBuf...clobbered by bcast
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //stuff results into swarm
   for(i = 0; i < m_SwarmSize; i++)
   {
      m_pSwarm[i].fx = m_pBigBuf[i];
   }/* end for() */
}/* end EvalSwarmParallel() */

/******************************************************************************
EvalSwarmSuperMUSE()

Compute objective functions of the swarm using SuperMUSE. This routine 
interfaces with the RepeatTasker SuperMUSE program, which assigns model 
evaluations to SuperMUSE clients on a first-come-first-served basis.
******************************************************************************/
void CCPSO::EvalSwarmSuperMUSE(void)
{  
   double val;
   bool bOk; 
   int i;
   ParameterGroup * pGroup;
   SuperMUSE * pSMUSE = GetSuperMusePtr();

   /* ----------------------------------------------------------------
   Generate task file that describes the desired parallel evaluations.
   This is somewhat analogous to the BcastGrid() operation used
   for MPI-parallel operations. Write the parameter values of each 
   swarm member as entries in the task file.

   Entries are first accumlated into a temp file to prevent the 
   SuperMUSE RepeatTasker program from prematurely processing the task
   file.
   ---------------------------------------------------------------- */
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_SwarmSize; i++)
   {
      //stuff the parameter group with values
      pGroup->WriteParams(m_pSwarm[i].x);
         
      //pass group to supermuse
      pSMUSE->WriteTask(pGroup);
   }/* end for() */

   //Finish task file (this will cause RepeatTasker to begin processing the job)
   pSMUSE->FinishTaskFile();

   //wait for SuperMUSE to report back (via the success or error files)
   bOk = pSMUSE->WaitForTasker();

   if(bOk == false) //SuperMUSE failed
   {
      LogError(ERR_SMUSE, "Reverting to serial execution.");
      DisableSuperMUSE();
      EvaluateSwarm();
   }
   else //SuperMUSE was successful
   {
      for(i = 0; i < m_SwarmSize; i++)
      {
         /* -----------------------------------------------
         Stuff the parameter group with ith swarm 
         member. This ensures that each objective function 
         gets associated with the correct parameter values.
         ------------------------------------------------ */
         pGroup->WriteParams(m_pSwarm[i].x);

         //stuff i-th result into chromosome pool
         val = pSMUSE->GatherResult(i);
         m_pSwarm[i].fx = val;
		 m_pSwarm[i].cx = pGroup->GetSpecialConstraint();
      }/* end for() */
   }/* end else() */
}/* end EvalSwarmSuperMUSE() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void CCPSO::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   int i, j, k, num;
   char * pTok;
   char * line;
   char tmp[DEF_STR_SZ];

   //swarm size recommended by Clerc
   m_SwarmSize = 10 + 2*(int)sqrt((double)(m_pModel->GetParamGroupPtr()->GetNumParams()));
   //num. gens recommended by Shi and Eberhart
   m_MaxGens = 350;
   m_Budget = m_SwarmSize*m_MaxGens;
   m_Constrict = 1.00;
   m_c1 = 2.00;  
   m_c2 = 2.00;
   m_Inertia = 1.2;
   m_RedRate = 0.10;
   m_LinRedFlag = true; //linear reduction to zero (always)
   m_InitType = LHS_INIT; //hypercube initialization always

   //read in PSO configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open PSO config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginCCPSO", pFileName) == true)
   {
      FindToken(pFile, "EndCCPSO", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginCCPSO", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndCCPSO") == NULL)
      {         
         if(strstr(line, "Budget") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_Budget); 
            if(m_Budget < 1)
            {
               LogError(ERR_FILE_IO, "Invalid CCPSO budget. Defaulting to 100.");
               m_Budget = 100;
            }
         }/*end else if() */         
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   /* initialize some or all swarm members to specified values */
   rewind(pFile);
   if(CheckToken(pFile, "BeginInitParams", pFileName) == true)
   {
      FindToken(pFile, "EndInitParams", pFileName);
      rewind(pFile);

      //allocate space for the parameter list
      num = m_pModel->GetParamGroupPtr()->GetNumParams();

      //count the number of entries
      FindToken(pFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      m_NumInit = 0;
      while(strstr(line, "EndInitParams") == NULL)
      {
         m_NumInit++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */

      //allocate space for entries
      if(m_NumInit > 0)
      {
         NEW_PRINT("double *", m_NumInit);
         m_pInit = new double * [m_NumInit];
         MEM_CHECK(m_pInit);
         for(i = 0; i < m_NumInit; i++)
         { 
            NEW_PRINT("double", num);
            m_pInit[i] = new double[num];
            MEM_CHECK(m_pInit[i]);
         }
      }/* end if() */

      //read in entries
      rewind(pFile);
      FindToken(pFile, "BeginInitParams", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      i = 0;
      while(strstr(line, "EndInitParams") == NULL)
      {
         pTok = line;
         //extract values, one-by-one, making any necessary conversions
         for(k = 0; k < num; k++)
         {
            j = ExtractString(pTok, tmp);
            j = ValidateExtraction(j, k, num, "PSO::InitFromFile()");
            pTok += j;            
            m_pInit[i][k] = m_pModel->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
         }/* end for() */                  
         i++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
   }/* end if() */

   fclose(pFile);

   /* -------------------------------------------------------------------
   Adjust swarm size and max. gens to reflect user-defined budget
   ------------------------------------------------------------------- */
   if(m_Budget > m_SwarmSize*m_MaxGens)
   {
      //inc. max. gens
      m_MaxGens = m_Budget/m_SwarmSize;
   }
   else if(m_Budget < (m_SwarmSize*m_MaxGens))
   {
      if(m_Budget < (m_SwarmSize*3)) //revise swarm size
      {
         m_SwarmSize = (int)MyMax(m_Budget/3, 3);
         m_MaxGens = m_Budget/m_SwarmSize;
      }
      else //revise max gens
      {
         m_MaxGens = m_Budget/m_SwarmSize;
      }
   }

   if(m_SwarmSize*m_MaxGens < m_Budget) m_MaxGens++;

} /* end InitFromFile() */

/******************************************************************************
CCPSO_Program()

Calibrate or optimize the model using CCPSO.
******************************************************************************/
void CCPSO_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("CCPSO", 1);
   CCPSO * PSO = new CCPSO(model);
   MEM_CHECK(PSO);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { PSO->Calibrate(); }
   else { PSO->Optimize(); }

   delete PSO;
   model->Destroy();
} /* end CCPSO_Program() */

/******************************************************************************
CCPSO_LEVMAR_Program()

Calibrate the model using CCPSO-Levenberg-Marquardt hybrid.
******************************************************************************/
void CCPSO_LEVMAR_Program(int argC, StringType argV[])
{
   int id;
   char file1[DEF_STR_SZ], file2[DEF_STR_SZ];
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("CCPSO", 1);
   CCPSO * PSO = new CCPSO(model);
   MEM_CHECK(PSO);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) 
   { 
      PSO->Calibrate();
      delete PSO;

      MPI_Comm_rank(MPI_COMM_WORLD, &id);
      sprintf(file1, "OstOutputPSO%d.txt", id);
      sprintf(file2, "OstOutput%d.txt", id);
      remove(file1);
      rename(file2, file1);

      sprintf(file1, "OstModelPSO%d.txt", id);
      sprintf(file2, "OstModel%d.txt", id);
      remove(file1);
      rename(file2, file1);

      sprintf(file1, "OstErrorsPSO%d.txt", id);
      sprintf(file2, "OstErrors%d.txt", id);
      remove(file1);
      rename(file2, file1);

      if(id == 0)
      {
         sprintf(file1, "OstStatusPSO%d.txt", id);
         sprintf(file2, "OstStatus%d.txt", id);
         remove(file1);
         rename(file2, file1);
      }

      NEW_PRINT("LevenbergAlgorithm", 1);
      LevenbergAlgorithm * LA = new LevenbergAlgorithm(model, false);
      MEM_CHECK(LA);
      LA->Calibrate(); 
      delete LA;      
   }
   else 
   {
      printf("Hybrid GML-CCPSO algorithm can only be used for calibration.\n"); 
   }
   
   model->Destroy();
}/* end CCPSO_LEVMAR_Program() */
