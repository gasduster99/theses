/******************************************************************************
File      : c.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

An implementation of a computation-constrained simulated annealing algorithm linked
with the Vanderbilt and Louie formulation.

Version History
10-31-09    lsm   Created
******************************************************************************/
#include "mpi_stub.h"
#include "CCVSA.h"
#include "ModelBackup.h"
#include "Utility.h"
#include "StatUtility.h"
#include "WriteUtility.h"
#include "Exception.h"
#include "Model.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************
CTOR

Initializes parameters, reading user-specified input, if available.
******************************************************************************/
CCVSA::CCVSA(ModelABC * pModel)
{
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   int numParams;
   int i, j;
   ParameterGroup * pGroup;

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();
   m_pStats = NULL;
   m_pMelts = NULL;
   m_dx     = NULL;
   m_Q      = NULL;
   m_QT     = NULL;
   m_u      = NULL;
   m_cov    = NULL;
   m_Shape  = NULL;
   m_x      = NULL;
   m_A      = NULL;

   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   numParams = pGroup->GetNumParams();   

   //initialize to reasonable defaults using the findings of:
   // Zhang and Wang (2003) Optimal Parameters Selection for Simulated Annealing with Limited Computational Effort
   // Wang and Wu (1999) A Parameter Set Design Procedure for the Simulated Annealing Algorithm Under the Computational Time Constraint
  
   //inner evaluations is 4 x np (num. parameters)
   m_MaxInner = 4*numParams;
   m_MaxOuter = 246; //T_final = 0.01, T_init = 3000, T_reduce = 0.95
   m_Budget = m_MaxInner*m_MaxOuter;
   //init. everything to reasonable defaults
   m_InitProb = m_CurProb = 0.99;
   m_NumOuter = 0;
   m_NumMelts = 10; //no melts 
   m_InitTemp = m_CurTemp = 3000.00;
   m_TempFactor = 0.95;
   m_MeltCount = 0;
   m_TransCount = 0;
   m_EquilCount = 0;
   m_NumAborts = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumUphill = 0;
   m_NumDownhill = 0;   

   NEW_PRINT("double", numParams);
   m_pTransPoint = new double[numParams];
   MEM_CHECK(m_pTransPoint);

   NEW_PRINT("ModelBackup", 1);
   m_pTransBackup = new ModelBackup(m_pModel);
   MEM_CHECK(m_pTransBackup);
   
   NEW_PRINT("double", numParams);
   m_pBest = new double[numParams];
   MEM_CHECK(m_pBest);

   NEW_PRINT("double", numParams);
   m_dx = new double[numParams];
   MEM_CHECK(m_dx);

   NEW_PRINT("double", numParams);
   m_A = new double[numParams];
   MEM_CHECK(m_A);

   NEW_PRINT("double", numParams);
   m_u = new double[numParams];
   MEM_CHECK(m_u);

   NEW_PRINT("double *", numParams);
   m_Q = new double * [numParams];
   MEM_CHECK(m_Q);

   NEW_PRINT("double *", numParams);
   m_QT = new double * [numParams];
   MEM_CHECK(m_QT);

   NEW_PRINT("double *", numParams);
   m_cov = new double * [numParams];
   MEM_CHECK(m_cov);

   NEW_PRINT("double *", numParams);
   m_Shape = new double * [numParams];
   MEM_CHECK(m_Shape);

   for(i = 0; i  < numParams; i++)
   {
      NEW_PRINT("double", numParams);
      m_Q[i] = new double[numParams];
      MEM_CHECK(m_Q[i]);

      NEW_PRINT("double", numParams);
      m_QT[i] = new double[numParams];
      MEM_CHECK(m_QT[i]);

      NEW_PRINT("double", numParams);
      m_cov[i] = new double[numParams];
      MEM_CHECK(m_cov[i]);

      NEW_PRINT("double", numParams);
      m_Shape[i] = new double[numParams];
      MEM_CHECK(m_Shape[i]);
   }/* end for() */

   for(i = 0; i  < numParams; i++)
   {
      for(j = 0; j  < numParams; j++)
      {
         if(i == j){ m_Q[i][j] = 1.00;}
         else      { m_Q[i][j] = 0.00;}
      }
   }
   
   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("CCVSA::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginSimulatedAlg", pFileName) == true)
   {
      FindToken(inFile, "EndSimulatedAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginSimulatedAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndSimulatedAlg") == NULL)
      {
         if(strstr(line, "Budget") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_Budget);
            if(m_Budget < 100)
            {
               LogError(ERR_FILE_IO, "Invalid CCSA budget. Defaulting to 100.");
               m_Budget = 100;
            }
         }
         line = GetNxtDataLine(inFile, pFileName);
      }/* end while() */
   } /* end if() */
   else
   {
      LogError(ERR_FILE_IO, "Using default algorithm setup.");
   }/* end else() */

   fclose(inFile);

   IncCtorCount();
}/* end default CTOR */

/******************************************************************************
DTOR
******************************************************************************/
CCVSA::~CCVSA(void)
{ 
   Destroy();
}/* end DTOR */

/******************************************************************************
Destroy()

Frees up memory used by the algorithm.
******************************************************************************/
void CCVSA::Destroy(void)
{
   int i;
   delete [] m_pBest;
   delete [] m_pTransPoint;
   delete [] m_pMelts;
   delete m_pStats;
   delete m_pTransBackup;

   delete [] m_dx;
   delete [] m_u;
   delete [] m_A;

   for(i = 0; i < m_pModel->GetParamGroupPtr()->GetNumParams(); i++)
   {
      delete [] m_Q[i];
      delete [] m_QT[i];
      delete [] m_cov[i];
      delete [] m_Shape[i];
	  if(m_x != NULL){ delete [] m_x[i];}
   }

   delete [] m_Q;
   delete [] m_QT;
   delete [] m_cov;
   delete [] m_Shape;
   delete [] m_x;
   delete [] m_Finner;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void CCVSA::WriteMetrics(FILE * pFile)
{   
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Simulated Annealing (Vanderbilt-Louie Implementation)\n");   
   fprintf(pFile, "Computational Budget    : %d\n", m_Budget);
   fprintf(pFile, "Outer Iterations        : %d\n", m_MaxOuter);
   fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
   fprintf(pFile, "Bookeeping Evals.       : %d\n", m_MaxInner);
   fprintf(pFile, "Temperature Reduction   : %.2lf%%\n", m_TempFactor*100.0);
   fprintf(pFile, "Initial Temperature     : %.2lf\n", m_InitTemp);
   fprintf(pFile, "Final Temperature       : %.2lf\n", m_CurTemp);
   fprintf(pFile, "Initial Pr[Acc]         : %.2lf%%\n", m_InitProb*100.0);
   fprintf(pFile, "Actual Final Pr[Acc]    : %.2lf%%\n", m_CurProb*100.0);
   fprintf(pFile, "Expected Final Pr[Acc]  : 50.00%%\n");
   fprintf(pFile, "Melting Evals           : %d\n", m_MeltCount);
   fprintf(pFile, "Transition Evals        : %d\n", m_TransCount);
   fprintf(pFile, "Equilibration Evals     : %d\n", m_EquilCount);
   //fprintf(pFile, "Total Evals             : %d\n", m_pModel->GetCounter());
   fprintf(pFile, "Rejected Transitions    : %d\n", m_NumAborts);
   fprintf(pFile, "Uphill Transitions      : %d\n", m_NumUphill);
   fprintf(pFile, "Downhill Transitions    : %d\n", m_NumDownhill);
   fprintf(pFile, "Upper Violations        : %d\n", m_NumUprViols);
   fprintf(pFile, "Lower Violations        : %d\n", m_NumLwrViols);

   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
StoreBest()

Saves the currently active parameter set into the m_pBest array.
******************************************************************************/
void CCVSA::StoreBest(void)
{
   m_pModel->GetParamGroupPtr()->ReadParams(m_pBest);
} /* end StoreBest() */

/******************************************************************************
RestoreBest()

Copies the parameter set stored in the m_pBest array into the model Parameter
Group, then the model is rerun, so that all constraints, response vars and
observations are consistent.
******************************************************************************/
void CCVSA::RestoreBest(void)
{
   m_pModel->GetParamGroupPtr()->WriteParams(m_pBest);
   m_pModel->Execute();
}/* end RestoreBest() */

/******************************************************************************
Melt()

'Melts' the design space to determine the initial temperature. The initial 
temperature is computed as the standard deviation of a user-specified number
of random moves.

Returns objective function value at starting location.
******************************************************************************/
double CCVSA::Melt(double initVal)
{
   int i; 
   double bestVal, avgVal, curVal, median;
   char c;

   //allocate space for melts
   if(m_pMelts == NULL)
   {
      NEW_PRINT("double", m_NumMelts); 
      m_pMelts = new double[m_NumMelts];
      MEM_CHECK(m_pMelts);
   }

   bestVal = curVal = initVal;

   WriteMelt(0, m_NumMelts, '.');
   avgVal = 0.00;
   for(i = 0; i < m_NumMelts; i++)
   {
      //make a random move
      GenerateRandomMove();
	   curVal = m_pModel->Execute();
      avgVal += curVal;
      m_MeltCount++;

      //accumulate values
      m_pMelts[i] = curVal;

      c = '+';
      if(curVal < bestVal)
      {         
         c = '-';
         StoreBest();
         bestVal = curVal;
      } /* end if() */

      WriteMelt(i+1, m_MaxInner, c);
   } /* end while() */
   WriteMelt(-1, -1, '.');

   //compute current convergence value (Eqn 18 of paper)
   avgVal /= (double)m_NumMelts;
   median = CalcMedian(m_pMelts, m_NumMelts);

   RestoreBest();   
   m_MeltCount++;

   m_CurTemp = m_InitTemp = CalcStdDev(m_pMelts, m_NumMelts);

   /* -------------------------------------------------------------------
   Adjust cooling schedule to reflect user-defined computational budget
   ------------------------------------------------------------------- */   
   while((m_Budget - m_NumMelts) < (m_MaxInner+1)*10) //need at least 10 temperature reductions
   {
      m_MaxInner--;
   }
   m_MaxOuter = (int)((m_Budget - m_NumMelts)/(m_MaxInner+1));

   m_TempFactor= exp((log(0.01) - log(m_InitTemp))/m_MaxOuter);

   if(m_MaxOuter*m_MaxInner < (m_Budget - m_NumMelts)) m_MaxOuter++;



   return (m_pModel->GetObjFuncVal());
}/* end Melt() */

/******************************************************************************
Optimize()

Optimize the objective function using the SA algorithm.
******************************************************************************/
void CCVSA::Optimize(void)
{
   int numEvals = 0;   
   StatusStruct pStatus;
   double curVal;
   int i;

   //write setup
   WriteSetup(m_pModel, "Computation-Constrained Simulated Annealing (Vanderbilt-Louie Implementation)");

   curVal = m_pModel->Execute();
   StoreBest();
   m_MeltCount++;

   curVal = Melt(curVal);
   numEvals = m_MeltCount;

   //write banner and initial result
   WriteBanner(m_pModel, "iter  obj. function  ", "Pct. Complete");
   WriteRecord(m_pModel, 0, curVal, (double)numEvals*100.00/(double)m_Budget);
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxOuter;
   pStatus.pct = (float)numEvals*(float)100.00/(float)m_Budget;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   for(i = 0; i < m_MaxOuter; i++)
   {
      if(IsQuit() == true){ break;}
      if(numEvals >= m_Budget){ pStatus.pct = 100.00; break;}

      curVal = Equilibrate(curVal);
      m_CurTemp *= m_TempFactor;
      numEvals = m_pModel->GetCounter();

      //write iteration result
      WriteRecord(m_pModel, (i+1), curVal, (double)numEvals*100.00/(double)m_Budget);
      pStatus.curIter = m_NumOuter = i+1;
      pStatus.pct = ((float)100.00*(float)(i+1))/(float)m_MaxOuter;
      pStatus.numRuns = numEvals;
      WriteStatus(&pStatus);

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   } /* end while() */

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   //write optimal results 
   WriteOptimal(m_pModel, curVal);
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);
   //write algorithm metrics
   WriteAlgMetrics(this);
}/* end Optimize() */

/******************************************************************************
Calibrate()

Calibrate the model using the SA algorithm.
******************************************************************************/
void CCVSA::Calibrate(void)
{ 
   FILE * pFile;
   char fileName[DEF_STR_SZ];
   int id;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
     
   Optimize();

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if(id == 0)
   {
      //compute statistics (variance and covariance)   
      m_pStats->CalcStats();

      sprintf(fileName, "OstOutput%d.txt", id);

      pFile = fopen(fileName, "a");   
      m_pStats->WriteStats(pFile);
      fclose(pFile);

      m_pStats->WriteStats(stdout);
   }
} /* end Calibrate() */

/******************************************************************************
Equilibrate()

Allows the system to come to equilibrium at the current temperature. This is 
accomplished by making a number of moves equal to m_MaxInner and choosing 
the best one. Returns the objective function value of the best move.
******************************************************************************/
double CCVSA::Equilibrate(double initVal)
{
   int i, j, m, n; 
   double bestVal, curVal, lastVal, avgVal, median;
   char c;
   ParameterGroup * pGroup;
   ParameterABC * pParam;
 
   //allocate m_x array, if needed
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   if(m_x == NULL)
   {
      NEW_PRINT("double *", n);
      m_x = new double *[n];
      MEM_CHECK(m_x);

      for(i = 0; i < n; i++)
      {
         NEW_PRINT("double", m_MaxInner);
         m_x[i] = new double[m_MaxInner];
         MEM_CHECK(m_x[i]);
      }/* end for() */

      //storage for inner loop f(x) calculations
      NEW_PRINT("double", m_MaxInner);
      m_Finner = new double [m_MaxInner];
      MEM_CHECK(m_Finner);
   }/* end if() */

   bestVal = curVal = initVal;

   WriteInnerEval(WRITE_SA, m_MaxInner, '.');

   //initialize probability metrics for this set of tranisitions
   m_NumProbTests = 0;
   m_TotProb      = 0.00;
   avgVal = 0.00;

   for(m = 0; m < m_MaxInner; m++)
   {
      lastVal = curVal;
      //curVal = Transition(curVal);
      curVal = GaussTransition(curVal);
      avgVal += curVal;
      m_Finner[m] = curVal;

      //store param values used in transition
      for(i = 0; i < n; i++)
      {
         pParam    = pGroup->GetParamPtr(i);
         m_x[i][m] = pParam->GetEstVal();
      }

      //update current best
      if(curVal < bestVal)
      {         
         StoreBest();
         bestVal = curVal;
      } /* end if() */

      if(curVal < lastVal)      { c = '-';}
      else if(curVal == lastVal){ c = '.';}
      else                      { c = '+';}
      WriteInnerEval(m+1, m_MaxInner, c);
   } /* end for() */

   //compute current convergence value (Eqn 18 of paper)
   avgVal /= (double)m_MaxInner;
   median = CalcMedian(m_Finner, m_MaxInner);

   //avg param transition (eqn 10 of the paper)
   for(i = 0; i < n; i++)
   {
      m_A[i] = 0.00;
      for(m = 0; m < m_MaxInner; m++)
      {
         m_A[i] += m_x[i][m];
      }
      m_A[i] /= (double)m_MaxInner;
   }

   //shape estimate (equation 11 of paper)
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         m_Shape[i][j] = 0.00; 
         for(m = 0; m < m_MaxInner; m++)
         {
            m_Shape[i][j] += (m_x[i][m] - m_A[i])*(m_x[j][m] - m_A[j]);
         }
         m_Shape[i][j] /= (double)m_MaxInner;
      }/* end for() */
   }/* end for() */

   //estimate cov of next iter (eqn 13 of paper)
   for(i = 0; i < n; i++)
   {
      for(j = 0; j < n; j++)
      {
         m_cov[i][j] = (3.00 * m_Shape[i][j])/(0.11*(double)m_MaxInner);
      }/* end for() */
   }/* end for() */

   //compute Q of next iter, using Cholesky decomposition (eqn 7 of paper)
   CholeskyDecomp(m_cov, m_Q, m_QT, n);

   //compute avg. probability of acceptance for the equilibration
   if(m_NumProbTests > 0)
   { 
      m_CurProb = m_TotProb / (double)m_NumProbTests;
      if(m_InitProb < 0.00){ m_InitProb = m_CurProb;}
   }

   WriteInnerEval(WRITE_ENDED, m_MaxInner, '.');

   RestoreBest();
   m_EquilCount++;

   return bestVal;
} /* end Equilibrate() */

/******************************************************************************
Transition()

Attempts to make a move from the current paramter set (location). For each 
move attempted, the resulting objective function value is tested against the
acceptance criteria (either the move reduces ob. func. or a randomly generated
number is less than the acceptance probability. Returns the value of the obj.
function at the revised location.
******************************************************************************/
double CCVSA::Transition(double initVal)
{
   int i, n;
   double increase, prob, r, curVal;
   double upr, lwr, val, adj;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

	//store initial Parameter values
   m_pTransBackup->Store();
  
   /*------------------------------------------
   Randomly initialize the 'u' vector. Uniformly
   dsitributed from -sqrt(3) to +sqrt(3).
   -------------------------------------------*/
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   for(i = 0; i < n; i++)
   {
      m_u[i] = (((double)MyRand()/(double)MY_RAND_MAX)*2.00*sqrt(3.00))-sqrt(3.00);
   }/* end for() */

   /*-----------------------------------------------
   Compute delta x, based on equation (5) of paper.
   ------------------------------------------------*/
   VectMult(m_Q, m_u, m_dx, n, n);
   
   /*---------------------------------------
   Make a move.
   ---------------------------------------*/	
   for(i = 0; i < n; i++)
   {
      pParam = pGroup->GetParamPtr(i);
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      val = pParam->GetEstVal();
      adj = val + m_dx[i];
      //if out of bounds, move half the distance
      if(adj > upr){adj = (upr+val)/2.00; m_NumUprViols++;}
      if(adj < lwr){adj = (val+lwr)/2.00; m_NumLwrViols++;}
      pParam->SetEstVal(adj);      
   }
	curVal = m_pModel->Execute();
   m_TransCount++;

   /*---------------------------------------
   Accept or reject move.
   ---------------------------------------*/	
	if(curVal <= initVal)
	{
      m_NumDownhill++;
   } /* end if() */
   else /* test against acceptance probability */
   {			
      increase = curVal - initVal;

      // generate a random between 0 and 1
      r = (double)MyRand() / (double)MY_RAND_MAX;

      //check if the increase is acceptable
      prob = exp(-increase / m_CurTemp);

      //probability metrics
      m_TotProb += prob;
      m_NumProbTests++;

      if(prob >= r)
      {
         // accept the move     
         m_NumUphill++;
      } /* end if() */
      else
      {
         m_pTransBackup->SemiRestore();         
         m_NumAborts++;
         curVal = initVal;
      } /* end else() */
   }/* end else() */

   return curVal;
} /* end Transition() */

/******************************************************************************
GaussTransition()

Attempts to make a move from the current paramter set (location). For each 
move attempted, the resulting objective function value is tested against the
acceptance criteria (either the move reduces ob. func. or a randomly generated
number is less than the acceptance probability. Returns the value of the obj.
function at the revised location.

The amount of parameter change induced by a transition is defined by the size
of the 'neighborhood' of adjacent solutions, which for continuous
SA is limited to +/-10% of the parameter range - centered on current parameter
location.
******************************************************************************/
double CCVSA::GaussTransition(double initVal)
{
   int i, n;
   double increase, prob, r;
   double upr, lwr, val, curVal;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

	//store initial Parameter values
   m_pTransBackup->Store();
  
   /*---------------------------------------
   Make a move.
   ---------------------------------------*/	
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   double sd = sqrt(MyMax(NEARLY_ZERO,fabs(initVal))/(double)n);

   for(i = 0; i < n; i++)
   {
      //perterb the parameter to a neighboring state (+/- 10% of range)
      pParam = pGroup->GetParamPtr(i);
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      val = curVal = pParam->GetEstVal();

      r = (double)MyRand() / (double)MY_RAND_MAX; //0 to 1
      // epsilon perturbation using normal distribution
      // centered on childVal with standard deviation estimated 
      // by the fitness value
      val = MyGaussRand(curVal, sd);

      //enforce parameter limits
      if(val > upr) val = curVal + (upr-curVal)*r;
      if(val < lwr) val = curVal - (curVal-lwr)*r;

      pParam->SetEstVal(val);
   }
	curVal = m_pModel->Execute();
   m_TransCount++;

   /*---------------------------------------
   Accept or reject move.
   ---------------------------------------*/	
	if(curVal <= initVal)
	{
      m_NumDownhill++;
   } /* end if() */
   else /* test against acceptance probability */
   {			
      increase = curVal - initVal;

      // generate a random between 0 and 1
      r = (double)MyRand() / (double)MY_RAND_MAX;

      //check if the increase is acceptable
      prob = exp(-increase / m_CurTemp);

      //probability metrics
      m_TotProb += prob;
      m_NumProbTests++;

      if(prob >= r)
      {
         // accept the move     
         m_NumUphill++;
      } /* end if() */
      else
      {
         m_pTransBackup->SemiRestore();         
         m_NumAborts++;
         curVal = initVal;
      } /* end else() */
   }/* end else() */

   return curVal;
} /* end GaussTransition() */

/******************************************************************************
GenerateRandomMove()

Adjusts a single parameter by a random displacement.
******************************************************************************/
void CCVSA::GenerateRandomMove(ParameterABC * pParam)
{
   double r; //random number
   double width; //width of move set
   double curVal; //current parameter value
   double upr; //upper limit of current move
   double lwr; //lower limit of current move
   double adjVal; //adjusted value 

   upr = pParam->GetUprBnd();
   lwr = pParam->GetLwrBnd();
   curVal = pParam->GetEstVal();
   width  = (upr - lwr);

   /*----------------------------------------------
   Generate a random displacement whose range will 
   meet side constraints
   ----------------------------------------------*/
   r = (((double)MyRand() / (double)MY_RAND_MAX) * width);   
   adjVal = r + lwr;

   //if out of bounds, move half the distance
   if(adjVal > upr){adjVal = (upr+curVal)/2.00; m_NumUprViols++;}
   if(adjVal < lwr){adjVal = (curVal+lwr)/2.00; m_NumLwrViols++;}
   
   pParam->SetEstVal(adjVal);
} /* end GenerateRandomMove() */

/******************************************************************************
GenerateRandomMove()

Adjusts parameters by a random displacement.
******************************************************************************/
void CCVSA::GenerateRandomMove(void)
{
   ParameterGroup * pGroup;
   ParameterABC * pParam;
   int numParams;
   int i;

   pGroup = m_pModel->GetParamGroupPtr();   
   numParams = pGroup->GetNumParams();
   for(i = 0; i < numParams; i++)
   {
      pParam = pGroup->GetParamPtr(i);
      GenerateRandomMove(pParam);
   }
} /* end GenerateRandomMove() */

/******************************************************************************
CCVSA_Program()

Calibrate the model using the SA algorithm.
******************************************************************************/
void CCVSA_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("VandSA", 1);
   CCVSA * VSA = new CCVSA(model);
   MEM_CHECK(VSA);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ VSA->Calibrate(); }
   else { VSA->Optimize(); }

   delete VSA;
   model->Destroy();
} /* end CCVSA_Program() */


