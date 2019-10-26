/******************************************************************************
File     : ParticleSwarmDESC.cpp
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

Particle Swarm Optimization (PSO) applies concepts of social behavior to
solve optimization problems. The PSO algorithm starts with a 'swarm' of 
particles (solutions) and "flies" this population through the design space in
search of the optimal solution. At each iteration, a given particle uses it's
own prior best solution (cognitive behavior) along with the current best 
solution of all particles (social behavior) to decide where to go next.

The PSODESC algorithm is a modified PSO, known as PSO using Diversity-Enhanced
Shuffled Complexes.PSODESC borrows concepts from the Shuffled Complex 
Evolutionary and Niched-Genetic algorithms to encourage efficient, effective, 
and reliable swarm exploration. After initialization, the swarm is subdivided 
into complexes using a diversity-enhancement scheme in which best-fit complex 
members are well-distributed using a niche concept.  The complexes are 
allowed to independently explore the design space, and are periodically merged, 
shuffled, and reformed.

Version History
09-14-07    lsm   created from ParticleSwarm.h
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "ParticleSwarmDESC.h"
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
ParticleSwarmDESC::ParticleSwarmDESC(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pSwarm = NULL;
   m_pSD = NULL;
   m_pSortedSwarm = NULL;
   m_pStats = NULL;
   m_pTrees = NULL;
   m_pInit = NULL;
   m_pComplex = NULL;
   m_TreeSize = 0;
   m_SwarmSize = 0;
   m_BestIdx = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumInit = 0;
   m_NumComplexes = 5;   //split swarm into fifths
   m_ShuffleRate = 0.1;    //shuffle 10% of the swarm   

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
ParticleSwarmDESC::~ParticleSwarmDESC(void)
{
   Destroy();
}/* end DTOR()*/

/******************************************************************************
Destroy()

Free up memory used by the PSO and it's member variables.
******************************************************************************/
void ParticleSwarmDESC::Destroy(void)
{
   int i;
   for(i = 0; i < m_SwarmSize; i++)
   {
      delete [] m_pSwarm[i].v;
      delete [] m_pSwarm[i].x;
      delete [] m_pSwarm[i].b;
   }
   delete [] m_pSwarm;
   delete [] m_pSortedSwarm;
   delete [] m_pSD;

   for(i = 0; i < m_NumComplexes; i++)
   {
      delete [] (m_pComplex[i].pSwarm);
   }
   delete [] m_pComplex;

   for(i = 0; i < m_NumInit; i++)
   {
      delete [] m_pInit[i];
   }
   delete [] m_pInit;
   
   delete [] m_Fmedian;
   delete m_pStats;
   delete [] m_pTrees;
   m_SwarmSize = 0;
   m_BestIdx = 0;
   m_NumComplexes = 0;

   delete [] m_pMyBuf;
   delete [] m_pTmpBuf;
   delete [] m_pBigBuf;
   delete [] m_pBuf;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
CalcPSOMedian()

Calculates and returns the median objective function of the swarm. This 
parameter is used in the termination criteria of the PSO algorithm.
******************************************************************************/  
double ParticleSwarmDESC::CalcPSOMedian(void)
{
   double med;
   int i;
    
   for(i = 0; i < m_SwarmSize; i++) { m_Fmedian[i] = m_pSwarm[i].fx; }
   med = CalcMedian(m_Fmedian, m_SwarmSize);

   return med;
}/* end CalcPSOMEdian() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using PSO.
******************************************************************************/
void ParticleSwarmDESC::Calibrate(void)
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
void ParticleSwarmDESC::Optimize(void)
{
   //char msg[DEF_STR_SZ];
   StatusStruct pStatus;
   int num, g;
   int id, lvl, idx;
   int i, j, k;
   double upr, lwr, r, range, sgn;
   double avg, x, pl, pg, r1, r2, v, vmin, median;
   double rval, init; //initial inertia
   double * pVals;
   LatinHypercube * pLHS = NULL;
   ParameterGroup * pGroup;
   ParticleStruct * pParticle, * pBest;
   
   InitFromFile(GetInFileName());

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if(id == 0)
   {
      WriteSetup(m_pModel, "PSODESC");
      //write banner
      WriteBanner(m_pModel, "gen   best value     ", "Convergence Value");
   }/* end if() */

   //allocate swarm
   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();

   NEW_PRINT("ParticleStruct", m_SwarmSize);
   m_pSwarm = new ParticleStruct[m_SwarmSize];
   MEM_CHECK(m_pSwarm);

   NEW_PRINT("ParticleStruct *", m_SwarmSize);
   m_pSortedSwarm = new ParticleStruct *[m_SwarmSize];
   MEM_CHECK(m_pSortedSwarm);

   NEW_PRINT("ComplexStruct", m_NumComplexes);
   m_pComplex = new ComplexStruct[m_NumComplexes];
   MEM_CHECK(m_pComplex);

   NEW_PRINT("double", m_SwarmSize);
   m_Fmedian = new double[m_SwarmSize];
   MEM_CHECK(m_Fmedian);

   //create complexes
   for(i = 0; i < m_NumComplexes; i++)
   {
      m_pComplex[i].size = (m_SwarmSize/m_NumComplexes);
      m_pComplex[i].fbprev = NEARLY_HUGE;
      m_pComplex[i].nfail = 0;
      NEW_PRINT("ParticleStruct *", m_pComplex[i].size);
      m_pComplex[i].pSwarm = new ParticleStruct *[m_pComplex[i].size];
      MEM_CHECK(m_pComplex[i].pSwarm);
   }
   NEW_PRINT("double", num);
   m_pSD = new double[num];
   MEM_CHECK(m_pSD);

   for(i = 0; i < num; i++)
   { 
      lwr = pGroup->GetParamPtr(i)->GetLwrBnd();
      upr = pGroup->GetParamPtr(i)->GetUprBnd();
      //estimate std. dev. by assuming parameter range is a 99.5% CI
      m_pSD[i] = (upr - lwr)/5.00;
   }/* end for() */

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
   if(m_InitType == LHS_INIT)
   {
      NEW_PRINT("LatinHypercube", 1);
      pLHS = new LatinHypercube(num, m_SwarmSize);
      MEM_CHECK(pLHS);

      for(j = 0; j < num; j++)
      { 
         lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
         upr = pGroup->GetParamPtr(j)->GetUprBnd();
         pLHS->InitRow(j, lwr, upr);
      }/* end for() */
   }/* end if() */
 
   lvl = idx = 0;
   for(i = 0; i < m_SwarmSize; i++) //for each particle
   {
      //initial velocity is 0.00
      for(j = 0; j < num; j++){ m_pSwarm[i].v[j] = 0.00;}

      if(m_InitType == RANDOM_INIT)
      {
         for(j = 0; j < num; j++) //for each parameter
         {
            //generate a random between lower and upper bound
            lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
            upr = pGroup->GetParamPtr(j)->GetUprBnd();
            range = upr - lwr;
            r = (double)MyRand() / (double)MY_RAND_MAX;
            rval = (r * range) + lwr;
            m_pSwarm[i].x[j] = rval;
            m_pSwarm[i].b[j] = rval;         
         }/* end for() */
      }/* end if() */
      else if(m_InitType == QUAD_TREE_INIT)
      {
         //initialize quad trees if needed
         if(m_pTrees == NULL)
         {
            m_TreeSize = num;
            NEW_PRINT("QuadTree", m_TreeSize);
            m_pTrees = new QuadTree[m_TreeSize];
            for(j = 0; j < m_TreeSize; j++)
            { 
               lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
               upr = pGroup->GetParamPtr(j)->GetUprBnd();
               m_pTrees[j].Init(lwr, upr);
            }/* end for() */
         }/* end if() */

         pVals = GetTreeCombo(lvl, idx, m_pTrees, m_TreeSize);
         //expand tree if needed.
         if(pVals == NULL)
         {
            for(j = 0; j < m_TreeSize; j++){ m_pTrees[j].Expand();}
            lvl++;
            idx = 0;            
            pVals = GetTreeCombo(lvl, idx, m_pTrees, m_TreeSize);
         }
         idx++;
         for(j = 0; j < num; j++)
         {
            m_pSwarm[i].x[j] = pVals[j];
            m_pSwarm[i].b[j] = pVals[j];
         }/* end for() */
         delete [] pVals;
      }/* end else if(QUAD_TREE_INIT) */
      else //LHS_INIT
      {
         for(j = 0; j < num; j++)
         { 
            rval = pLHS->SampleRow(j);
            m_pSwarm[i].x[j] = rval;
            m_pSwarm[i].b[j] = rval;
         }/* end for() */
      }/* end else() */
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

   if(m_InitType == LHS_INIT){ pLHS->Destroy();}

   //evaluate swarm, possibly in parallel
   EvaluateSwarm();

   //perform intermediate bookkeeping
   m_pModel->Bookkeep(false);

   for(i = 0; i < m_SwarmSize; i++){m_pSwarm[i].fb = m_pSwarm[i].fx;}

   //determine the best particle and average value
   m_BestIdx = 0;
   avg = 0.00;
   m_Best = m_pSwarm[0].fb;
   for(i = 0; i < m_SwarmSize; i++)
   {
      avg += m_pSwarm[i].fx;
      if(m_Best > m_pSwarm[i].fx)
      {
         m_Best = m_pSwarm[i].fx;
         m_BestIdx = i;
      }/* end if() */
   }/* end if() */
   avg /= m_SwarmSize;
   median = CalcPSOMedian();
   //current convergence value
   m_CurStop = fabs((median - m_Best)/median);

   if(id == 0)
   {
      //write initial config.
      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);
      WriteRecord(m_pModel, 0, m_Best, m_CurStop);
      pStatus.curIter = 0;
      pStatus.maxIter = m_MaxGens;
      pStatus.pct = 0.00;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
   }/* end if() */

   init = m_Inertia;
   //main optimization loop   
   for(g = 0; g < m_MaxGens; g++)
   {
      pStatus.curIter = m_CurGen = g+1;
      if(IsQuit() == true){ break;}
      if(m_CurStop < m_StopVal){ pStatus.pct = 100.00; break;}

      if(id == 0)
      {
         //form and shuffle complexes
         ShuffleComplexes();

         //update velocities, parameters and objective functions
         for(i = 0; i < m_NumComplexes; i++)
         {
            for(j = 0; j < m_pComplex[i].size; j++) //for each particle
            {
               pParticle = m_pComplex[i].pSwarm[j];
               pBest = m_pComplex[i].pSwarm[0]; //best of complex is always first

               for(k = 0; k < num; k++) //for each parameter
               {
                  //intermediary variables
                  x = pParticle->x[k];
                  pl = pParticle->b[k];
                  pg = pBest->b[k]; 

                  //random weights
                  r1 = (double)MyRand() / (double)MY_RAND_MAX;
                  r2 = (double)MyRand() / (double)MY_RAND_MAX;

                  //revised velocity
                  v = pParticle->v[k];
                  v = m_Constrict*((m_Inertia*v) + m_c1*r1*(pl-x) + m_c2*r2*(pg-x));
                  vmin = (0.01*fabs(x))/(g+1); //minimum perturbation, prevents stagnation
                  if(fabs(v) < vmin)
                  {                     
                     //preserve direction, but adjust randomized minimum velocity
                     sgn = (double)MyRand() / (double)MY_RAND_MAX; //random direction
                     if(sgn >= 0.50) v = +((1.00+r1)*vmin);
                     else            v = -((1.00+r2)*vmin);
                     //report stagnation
                     //sprintf(msg, "particle velocity stagnation, revised velocity = %E", v);
                     //LogError(ERR_STALL, msg);
                  }
                  pParticle->v[k] = v;
            
                  //revised position
                  pParticle->x[k] = x + v;

                  //constrain revised position to stay within parameter limits
                  lwr = pGroup->GetParamPtr(k)->GetLwrBnd();
                  upr = pGroup->GetParamPtr(k)->GetUprBnd();
                  //move half the distance to the goal!
                  if(pParticle->x[k] > upr){pParticle->x[k] = (upr+x)/2.00; m_NumUprViols++;}
                  if(pParticle->x[k] < lwr){pParticle->x[k] = (x+lwr)/2.00; m_NumLwrViols++;}
               }/* end for(params in a particle) */
            }/* end for(particles in a complex) */
         }/* end for(each complex) */
      } /* end if() */

      //evaluate swarm, possibly in parallel
      EvaluateSwarm();

      //reduce inertia
      if(m_LinRedFlag == true) //linearly reducing to zero
      {
         m_Inertia = init;
         m_RedRate = (double)g/(double)m_MaxGens;
      }
      m_Inertia *= (1.00 - m_RedRate);
   
      //revise avg, and local and global best
      avg = 0.00;
      for(i = 0; i < m_SwarmSize; i++)
      {
         avg += m_pSwarm[i].fx;
         //revise local best
         if(m_pSwarm[i].fx < m_pSwarm[i].fb)
         {
            for(j = 0; j < num; j++){ m_pSwarm[i].b[j] = m_pSwarm[i].x[j];}
            m_pSwarm[i].fb = m_pSwarm[i].fx;
         }/* end if() */
         //revise global best
         if(m_Best > m_pSwarm[i].fx)
         {
            m_Best = m_pSwarm[i].fx;
            m_BestIdx = i;
         }/* end if() */
      }/* end if() */
      avg /= m_SwarmSize;
      median = CalcPSOMedian();
      //current convergence value
      m_CurStop = fabs((median - m_Best)/median);
      pGroup->WriteParams(m_pSwarm[m_BestIdx].b);

      if(id == 0){ WriteRecord(m_pModel, (g+1), m_Best, m_CurStop);}
      pStatus.pct = ((float)100.00*(float)(g+1))/(float)m_MaxGens;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */

   m_Inertia = init; //reset inertia

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
ShuffleComplexes()

Subdivide the swarm into a set of diversity-enhanced complexes.
******************************************************************************/
void ParticleSwarmDESC::ShuffleComplexes(void)
{
   // char msg[DEF_STR_SZ];
   int c1, c2; //random complexes
   int e1, e2; ///random complex entries
   int z1, z2; //complex size
   int nswaps; //number of shuffles
   int i, j, k, p;
   double fmax, fcur;
   int fidx;
   //double d;
   ParticleStruct * tmp;

   //fill in sorted swarm array
   for(i = 0; i < m_SwarmSize; i++)
   {
      m_pSortedSwarm[i] = &(m_pSwarm[i]);
   }   
   
   //if only one complex, degenerate into standard PSO
   if(m_NumComplexes == 1)
   {
      for(i = 0; i < m_SwarmSize; i++)
      {
         m_pComplex[0].pSwarm[i] = m_pSortedSwarm[i];
      }
      return;
   }

   //NULL out the complexes
   for(i = 0; i < m_NumComplexes; i++)
   {
      for(j = 0; j < m_pComplex[i].size; j++)
      {
         m_pComplex[i].pSwarm[j] = NULL;
      }
   }/* end for() */

   //fill in unshuffled complexes
   j = 0;
   for(i = 0; i < m_NumComplexes; i++)
   {
      //sort the swarm based on local-best objective function value
      SortSwarmByObjFunc(&(m_pSortedSwarm[j]), (m_SwarmSize-j));

      //assign best-fit member to the complex
      m_pComplex[i].pSwarm[0] = m_pSortedSwarm[j];

      //place best-fit member at best-fit location
      m_pComplex[i].pSwarm[0]->fx = m_pComplex[i].pSwarm[0]->fb;
      for(p = 0; p < m_pComplex[i].pSwarm[0]->n; p++)
      {
         m_pComplex[i].pSwarm[0]->x[p] = m_pComplex[i].pSwarm[0]->b[p];   
      }

      //update failure metrics
      if(m_pComplex[i].pSwarm[0]->fb >= m_pComplex[i].fbprev)
      {
         m_pComplex[i].nfail++;
      }
      else
      {
         m_pComplex[i].fbprev = m_pComplex[i].pSwarm[0]->fb;
         m_pComplex[i].nfail = 0;
      }

      //sort the swarm based on distance of local-best to best-fit member
      SortSwarmByDistance(&(m_pSortedSwarm[j]), (m_SwarmSize-j));

      z1 = m_pComplex[i].size;

      //assign nearest (sz-1) members to the complex
      for(k = 1; k < z1; k++)
      {
         j++;
         m_pComplex[i].pSwarm[k] = m_pSortedSwarm[j];
      }
      j++;
   }/* end for() */

   //in 'polishing' mode?
   if((m_MaxGens - m_CurGen) <= m_PolishGens)
   {
      //swap in global best, swap out local worst
      for(i = 1; i < m_NumComplexes; i++)
      {
         //find local worst
         fmax = m_pComplex[i].pSwarm[0]->fb;
         fidx = 0;
         for(j = 1; j < m_pComplex[i].size; j++)
         {
            fcur = m_pComplex[i].pSwarm[j]->fb;
            if(fcur > fmax)
            {
               fmax = fcur;
               fidx = j;
            }
         }/* end for() */

         tmp = m_pComplex[i].pSwarm[0]; //preserve local best
         m_pComplex[i].pSwarm[0] = m_pComplex[0].pSwarm[0]; //swap in global best
         m_pComplex[i].pSwarm[fidx] = tmp; //swap out local worst, rplace with local best
      }/* end for() */
   }/* end if() */
   else //perform shuffling
   {  
      nswaps = (int)(m_ShuffleRate*m_SwarmSize);      
      for(i = 0; i < nswaps; i++)
      {
         //randomly select two complexes and swap two random entries.
         c1 = MyRand() % m_NumComplexes;
         do{
            c2 = MyRand() % m_NumComplexes;
         }while(c1 == c2);

         z1 = m_pComplex[c1].size;
         z2 = m_pComplex[c2].size;

         //best-fit member of each complex is excluded
         e1 = 1 + (MyRand() % (z1-1));
         e2 = 1 + (MyRand() % (z2-1));

         //perform swap
         tmp = m_pComplex[c1].pSwarm[e1];
         m_pComplex[c1].pSwarm[e1] = m_pComplex[c2].pSwarm[e2];
         m_pComplex[c2].pSwarm[e2]= tmp;
      }/* end for() */
   }/* end else() */

   //if an individual complex has stalled, set it into polishing mode
   for(i = 1; i < m_NumComplexes; i++)
   {
      if(m_pComplex[i].nfail >= 5)
      {
         //log the stall
         // sprintf(msg, "Complex #%d stalled in gen #%d, switching to polishing mode", i, m_CurGen);
         // LogError(ERR_STALL, msg);

         //find local worst
         fmax = m_pComplex[i].pSwarm[0]->fb;
         fidx = 0;
         for(j = 1; j < m_pComplex[i].size; j++)
         {
            fcur = m_pComplex[i].pSwarm[j]->fb;
            if(fcur > fmax)
            {
               fmax = fcur;
               fidx = j;
            }
         }/* end for() */

         tmp = m_pComplex[i].pSwarm[0]; //preserve local best
         m_pComplex[i].pSwarm[0] = m_pComplex[0].pSwarm[0]; //swap in global best
         m_pComplex[i].pSwarm[fidx] = tmp; //swap out local worst, replace with local best
      }/* end if() */
   }/* end for() */

   //check that complexes have been filled
   //printf("Shuffled Complexes\n");
   for(i = 0; i < m_NumComplexes; i++)
   {
      //printf("Complex # %d\n", (i+1));
      for(j = 0; j < m_pComplex[i].size; j++)
      {
         if(m_pComplex[i].pSwarm[j] == NULL)
         {
            LogError(ERR_NULL_PTR, "shuffled complexes not formed correctly");
            ExitProgram(1);
         }
         //d = GetParticleDistance(m_pComplex[i].pSwarm[0], m_pComplex[i].pSwarm[j]);
         //printf("(%d) %E   %E\n", (j+1), m_pComplex[i].pSwarm[j]->fb, d);
      }
   }/* end for() */
   //printf("Done\n");
}/* end ShuffleComplexes() */

/******************************************************************************
SortSwarmByObjFunc()

Sort the swarm according to locally optimal objective function values.
******************************************************************************/
void ParticleSwarmDESC::SortSwarmByObjFunc(ParticleStruct ** pSwarm, int size)
{
   int i,j;
   ParticleStruct * pTmp;
   
   for(i = 0; i < size; i++)
   {
      for(j = (i+1); j < size; j++)
      {
         if(pSwarm[j]->fb < pSwarm[i]->fb)
         {
            pTmp = pSwarm[j];
            pSwarm[j] = pSwarm[i];
            pSwarm[i] = pTmp;
         }
      }
   }

   //printf("Sorted Swarm (By Obj. Func.)\n");
   //for(i = 0; i < size; i++)
   //{
      //printf("(%d) %E\n", i+1, pSwarm[i]->fb);
   //}
}/* end SortSwarmByObjFunc() */

/******************************************************************************
SortSwarmByDistance()

Sort the swarm according to distance between locally optimal objective function 
values.
******************************************************************************/
void ParticleSwarmDESC::SortSwarmByDistance(ParticleStruct ** pSwarm, int size)
{
   int i,j;
   double dmin, d;
   ParticleStruct * pTmp;
   
   for(i = 1; i < size; i++)
   {
      dmin = GetParticleDistance(pSwarm[0], pSwarm[i]);

      for(j = (i+1); j < size; j++)
      {
         d = GetParticleDistance(pSwarm[0], pSwarm[j]);

         if(d < dmin)
         {
            pTmp = pSwarm[j];
            pSwarm[j] = pSwarm[i];
            pSwarm[i] = pTmp;
            dmin = d;
         }
      }
   }

   //printf("Sorted Swarm (By Distance)\n");
   //for(i = 0; i < size; i++)
   //{
        //d = GetParticleDistance(pSwarm[0], pSwarm[i]);
      //printf("(%d) %E   %E\n", i+1, pSwarm[i]->fb, d);
   //}
}/* end SortSwarmByDistance() */

/******************************************************************************
GetParticleDistance()

Compute the multi-dimensional distance between the local-best of the 
two particles. Distance is normalized to an estimated std. dev. of each 
respective parameter.
******************************************************************************/
double ParticleSwarmDESC::GetParticleDistance(ParticleStruct * p1, ParticleStruct * p2)
{
   int i;
   double v1, v2, v3;
   double d = 0.00;   

   for(i = 0; i < p1->n; i++)
   {
      v1 = p1->b[i];
      v2 = p2->b[i];
      v3 = (v1 - v2)/m_pSD[i];
      d += (v3*v3);
   }
   d = sqrt(d);
   return d;
}/* GetParticleDistance() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void ParticleSwarmDESC::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : PSODESC\n");
   fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
   fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
   fprintf(pFile, "Max Generations         : %d\n", m_MaxGens);
   fprintf(pFile, "Actual Generations      : %d\n", m_CurGen);
   fprintf(pFile, "Polishing Generations   : %d\n", m_PolishGens);
   fprintf(pFile, "Swarm Size              : %d\n", m_SwarmSize);
   fprintf(pFile, "Number of Complexes     : %d\n", m_NumComplexes);
   fprintf(pFile, "Shuffle Rate            : %.2lf\n", m_ShuffleRate);
   fprintf(pFile, "Constriction Factor     : %.2lf\n", m_Constrict);  
   fprintf(pFile, "Cognitive Weight        : %.2lf\n", m_c1);
   fprintf(pFile, "Social Weight           : %.2lf\n", m_c2);
   fprintf(pFile, "Inertia Weight          : %.2lf\n", m_Inertia);
   
   fprintf(pFile, "Inertia Reduction Rate  : ");
   if(m_LinRedFlag == true) fprintf(pFile, "Linear reduction to zero\n");
   else                     fprintf(pFile, "%.2lf\n", m_RedRate);

   fprintf(pFile, "Initialization Method   : ");
   if(m_InitType == RANDOM_INIT){ fprintf(pFile, "Random\n");}
   else if(m_InitType == QUAD_TREE_INIT){ fprintf(pFile, "Quad-Tree\n");}
   else if(m_InitType == LHS_INIT){ fprintf(pFile, "Latin Hypercube Sampling\n");}
   else { fprintf(pFile, "Unknown\n");}

   //fprintf(pFile, "Total Evals             : %d\n", m_pModel->GetCounter());      
   fprintf(pFile, "Upper Violations        : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations        : %d\n", m_NumLwrViols);

   m_pModel->WriteMetrics(pFile);
   if(m_CurStop <= m_StopVal)
   {
      fprintf(pFile, "Algorithm successfully converged on a solution\n");
   }
   else
   {
      fprintf(pFile, "Algorithm failed to converge on a solution, more generations may be needed\n");
   }   
}/* end WriteMetrics() */

/******************************************************************************
EvaluateSwarm()

Evaluates the objective function of each particle in the swarm.
******************************************************************************/
void ParticleSwarmDESC::EvaluateSwarm(void)
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
            val = m_pModel->Execute();
            m_pSwarm[i].fx = val;
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
void ParticleSwarmDESC::BcastSwarm(void)
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
void ParticleSwarmDESC::EvalSwarmParallel(void)
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
void ParticleSwarmDESC::EvalSwarmSuperMUSE(void)
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
      }/* end for() */
   }/* end else() */
}/* end EvalSwarmSuperMUSE() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void ParticleSwarmDESC::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   int i, j, k, num;
   char * pTok;
   char * line;
   char tmp[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];

   m_StopVal = 0.001;
   m_NumComplexes = 5;
   m_ShuffleRate = 0.1;
   m_SwarmSize = 20;
   m_MaxGens = 50;
   m_PolishGens = 0;
   m_Constrict = 1.00;
   m_c1 = 2.00;
   m_c2 = 2.00;
   m_Inertia = 1.2;   
   m_RedRate = 0.10;
   m_LinRedFlag = false;
   m_InitType = RANDOM_INIT;

   //read in PSO configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open PSO config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginPSODESC", pFileName) == true)
   {
      FindToken(pFile, "EndPSODESC", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginPSODESC", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndPSODESC") == NULL)
      {         
         if(strstr(line, "SwarmSize") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_SwarmSize); 
         }/*end else if() */         
         else if(strstr(line, "NumGenerations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxGens);
         }
         else if(strstr(line, "PolishingGenerations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_PolishGens);
         }
         else if(strstr(line, "ConstrictionFactor") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Constrict);
         }
         else if(strstr(line, "CognitiveParam") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_c1);
         }
         else if(strstr(line, "SocialParam") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_c2);
         }
         else if(strstr(line, "InertiaWeight") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Inertia);
         }
         else if(strstr(line, "InertiaReductionRate") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            MyStrLwr(tmp2);
            if(strcmp(tmp2, "linear") == 0){ m_LinRedFlag = true;}
            else{sscanf(line, "%s %lf", tmp, &m_RedRate);}
         }
         else if(strstr(line, "InitPopulationMethod") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2);
            MyStrLwr(tmp2);
            if(strcmp(tmp2, "random") == 0) {m_InitType = RANDOM_INIT;}
            else if(strcmp(tmp2, "quadtree") == 0) {m_InitType = QUAD_TREE_INIT;}
            else if(strcmp(tmp2, "lhs") == 0) {m_InitType = LHS_INIT;}
         }
         else if(strstr(line, "ConvergenceVal") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_StopVal);
         }
         else if(strstr(line, "NumComplexes") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumComplexes);
         }
         else if(strstr(line, "ShuffleRate") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_ShuffleRate);
         }
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   

   /* initialize some or all swarm members to specied values */
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
            j = ValidateExtraction(j, k, num, "PSODESC::InitFromFile()");
            pTok += j;            
            m_pInit[i][k] = m_pModel->GetParamGroupPtr()->GetParamPtr(k)->ConvertInVal(atof(tmp));
         }/* end for() */                  
         i++;
         line = GetNxtDataLine(pFile, pFileName);
      }/* end while() */
   }/* end if() */

   /* make sure DESC parameters are consistent with PSO configuration */
   if((m_NumComplexes < 1) || (m_NumComplexes > m_SwarmSize))
   {
      LogError(ERR_BAD_ARGS, "Invalid number of complexes, defaulting to 1");
      m_NumComplexes = 1;
   }
   if((m_ShuffleRate < 0.00) || (m_ShuffleRate > 1.00))
   {
      LogError(ERR_BAD_ARGS, "Invalid shuffle rate, defaulting to 0.10");
      m_ShuffleRate = 0.10;
   }
   if((m_SwarmSize % m_NumComplexes) != 0)
   {
      LogError(ERR_BAD_ARGS, "Swarm Size and Number of Complexes are not evenly divisible");
   }

   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
PSODESC_Program()

Calibrate or optimize the model using PSODESC.
******************************************************************************/
void PSODESC_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("ParticleSwarmDESC", 1);
   ParticleSwarmDESC * PSODESC = new ParticleSwarmDESC(model);
   MEM_CHECK(PSODESC);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { PSODESC->Calibrate(); }
   else { PSODESC->Optimize(); }

   delete PSODESC;
   model->Destroy();
} /* end PSODESC_Program() */

/******************************************************************************
PSODESC_LEVMAR_Program()

Calibrate the model using PSODESC-Levenberg-Marquardt hybrid.
******************************************************************************/
void PSODESC_LEVMAR_Program(int argC, StringType argV[])
{
   int id;
   char file1[DEF_STR_SZ], file2[DEF_STR_SZ];
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("ParticleSwarmDESC", 1);
   ParticleSwarmDESC * PSODESC = new ParticleSwarmDESC(model);
   MEM_CHECK(PSODESC);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) 
   { 
      PSODESC->Calibrate();
      delete PSODESC;

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
      printf("Hybrid GML-PSODESC algorithm can only be used for calibration.\n"); 
   }
   
   model->Destroy();
}/* end PSO_LEVMAR_Program() */
