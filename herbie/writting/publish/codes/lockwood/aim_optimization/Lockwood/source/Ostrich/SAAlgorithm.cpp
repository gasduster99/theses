/******************************************************************************
File      : SAAlgorithm.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

An implementation of the simulated annealing algorithm for continuously 
varying parameters.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   adjusted melt operation
                  outputs obj. func. str.
07-08-04    lsm   switched to ParameterABC
                  WriteSetup() is now the first action of algorithm
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
09-01-04    lsm   Algorithm now based on: Vanderbilt and Louie. 1984. "A Monte 
                  Carlo Simulated Annealing Approach to Optimization over 
                  Continuous Variables". Journal of Computational Physics. 
                  vol. 56, pg. 259-271.
11-18-04    lsm   Added convergence criteria, based on median F() of inner loop.
11-30-04    lsm   Replaced calls to time() with MyTime()
10-21-05    lsm   Switched to homegrown implementation that is easier to study
                  than the Vanderbilt and Louie implementation
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC), added
                  model init. and bookkepping calls. Some statistics can now 
                  be calculated in parallel. 
******************************************************************************/
#include "mpi_stub.h"
#include "SAAlgorithm.h"
#include "Model.h"
#include "ModelBackup.h"
#include "Utility.h"
#include "StatUtility.h"
#include "WriteUtility.h"
#include "Exception.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************
CTOR

Initializes parameters, reading user-specified input, if available.
******************************************************************************/
SAAlgorithm::SAAlgorithm(ModelABC * pModel)
{
   FILE * inFile;
   char * line;
   char tmp[DEF_STR_SZ];
   int numParams;
   ParameterGroup * pGroup;

   RegisterAlgPtr(this);
   IroncladString pFileName = GetInFileName();
   m_pStats = NULL;
   m_pMelts = NULL;
   m_Finner = NULL;

   //init. everything to reasonable defaults
   m_InitProb = m_CurProb = -1.00;
   m_StopVal = 0.001;
   m_CurStop = 1.00;
   m_NumOuter = 0;
   m_MaxOuter = 20;
   m_MaxInner = 10;
   m_NumMelts = 100;
   m_InitTemp = m_CurTemp = 10.00;
   m_TempFactor = 0.9;
   m_MeltCount = 0;
   m_TransCount = 0;
   m_EquilCount = 0;
   m_NumAborts = 0;
   m_NumUprViols = 0;
   m_NumLwrViols = 0;
   m_NumUphill = 0;
   m_NumDownhill = 0;

   m_pModel = pModel;
   pGroup = m_pModel->GetParamGroupPtr();
   numParams = pGroup->GetNumParams();   

   NEW_PRINT("ModelBackup", 1);
   m_pTransBackup = new ModelBackup(m_pModel);
   MEM_CHECK(m_pTransBackup);
   
   NEW_PRINT("double", numParams);
   m_pBest = new double[numParams];
   MEM_CHECK(m_pBest);

   inFile = fopen(pFileName, "r");  
   if(inFile == NULL) 
   {
      FileOpenFailure("SAAlgorithm::CTOR", pFileName);
   }/* end if() */ 

   if(CheckToken(inFile, "BeginSimulatedAlg", pFileName) == true)
   {
      FindToken(inFile, "EndSimulatedAlg", pFileName);
      rewind(inFile);

      FindToken(inFile, "BeginSimulatedAlg", pFileName);
      line = GetNxtDataLine(inFile, pFileName);
      
      while(strstr(line, "EndSimulatedAlg") == NULL)
      {
         if(strstr(line, "NumInitialTrials") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumMelts);
         }
         else if(strstr(line, "TemperatureScaleFactor") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_TempFactor);
         }
         else if(strstr(line, "OuterIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxOuter);
         }
         else if(strstr(line, "InnerIterations") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxInner);
         }
         else if(strstr(line, "ConvergenceVal") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_StopVal);
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
SAAlgorithm::~SAAlgorithm(void)
{ 
   Destroy();
}/* end DTOR */

/******************************************************************************
Destroy()

Frees up memory used by the algorithm.
******************************************************************************/
void SAAlgorithm::Destroy(void)
{
   delete [] m_pBest;
   delete [] m_pMelts;
   delete m_pStats;
   delete m_pTransBackup;

   delete [] m_Finner;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void SAAlgorithm::WriteMetrics(FILE * pFile)
{
   double sd, p;
   sd = ((m_InitTemp/100.00)-m_dEavg)/3.00;
   p = exp(-m_dEavg/m_CurTemp);
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Simulated Annealing for Continuous Parameters\n");
   fprintf(pFile, "Desired Convergence Val : %E\n", m_StopVal);
   fprintf(pFile, "Actual Convergence Val  : %E\n", m_CurStop);
   fprintf(pFile, "Max Outer Iterations    : %d\n", m_MaxOuter);
   fprintf(pFile, "Actual Outer Iterations : %d\n", m_NumOuter);
   fprintf(pFile, "Inner Iterations        : %d\n", m_MaxInner);
   fprintf(pFile, "Temperature Reduction   : %.2lf%%\n", m_TempFactor*100.0);
   fprintf(pFile, "Initial Temperature     : %.2lf\n", m_InitTemp);
   fprintf(pFile, "Avg. Energy Change      : %.2lf\n", m_dEavg);
   fprintf(pFile, "Std. Dev. Energy Change : %.2lf\n", sd);
   fprintf(pFile, "Final Temperature       : %.2lf\n", m_CurTemp);
   fprintf(pFile, "Initial Pr[Acc]         : %.2lf%%\n", m_InitProb*100.0);
   fprintf(pFile, "Actual Final Pr[Acc]    : %.2lf%%\n", m_CurProb*100.0);
   fprintf(pFile, "Expected Final Pr[Acc]  : %.2lf%%\n", p*100.0);
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
   if(m_CurStop <= m_StopVal)
   {
      fprintf(pFile, "Algorithm successfully converged on a solution\n");
   }
   else
   {
      fprintf(pFile, "Algorithm failed to converge on a solution, more outer iterations may be needed\n");
   }      
}/* end WriteMetrics() */

/******************************************************************************
StoreBest()

Saves the currently active parameter set into the m_pBest array.
******************************************************************************/
void SAAlgorithm::StoreBest(void)
{
   m_pModel->GetParamGroupPtr()->ReadParams(m_pBest);
} /* end StoreBest() */

/******************************************************************************
RestoreBest()

Copies the parameter set stored in the m_pBest array into the model Parameter
Group, then the model is rerun, so that all constraints, response vars and
observations are consistent.
******************************************************************************/
void SAAlgorithm::RestoreBest(void)
{
   m_pModel->GetParamGroupPtr()->WriteParams(m_pBest);
   m_pModel->Execute();
}/* end RestoreBest() */

/******************************************************************************
Melt()

'Melts' the design space to determine the initial temperature. The initial 
temperature is computed so that a statistically large energy increase (dE) from a 
sample of random moves (melting trials) will be accepted with ~100% probability. 
That is:
   T0 = -dEmax/ln(0.99) ~ 100*dEmax
Where dEmax is computed as 3 std. deviations from the mean dE estimated from the
melting trials.

Returns objective function value at starting location.
******************************************************************************/
double SAAlgorithm::Melt(double initVal)
{
   int i; 
   double Ebest, Eprev, Ecur; //objective function values
   double Emed;               //melting trial metrics
   double dE, dEmax, dEavg;   //energy change (delta E)
   double * pdE;              //array of positive energy changes
   char c;

   //allocate space for melts (needed for convergence value)
   if(m_pMelts == NULL)
   {
      NEW_PRINT("double", m_NumMelts); 
      m_pMelts = new double[m_NumMelts];
      pdE      = new double[m_NumMelts];
      MEM_CHECK(m_pMelts);
   }

   Ebest = Ecur = initVal;
   dE = dEavg = 0.00;
   
   WriteMelt(0, m_NumMelts, '.');
   for(i = 0; i < m_NumMelts; i++)
   {
      //make a random move
      GenerateRandomMove();
      Eprev = Ecur;
	   Ecur  = m_pModel->Execute();
      dE = Ecur - Eprev;      
      m_MeltCount++;

      //accumulate values
      pdE[i] = fabs(dE);
      dEavg += pdE[i];
      m_pMelts[i] = Ecur;

      c = '+';      
      if(dE < 0.00) //decrease in energy
      {         
         c = '-';
         if(Ecur < Ebest) //new best solution?
         {
            StoreBest();
            Ebest = Ecur;
         }
      } /* end if() */      

      WriteMelt(i+1, m_MaxInner, c);
   } /* end while() */
   WriteMelt(-1, -1, '.');

   dEavg /= (double)m_NumMelts;
   //compute current convergence value
   Emed = CalcMedian(m_pMelts, m_NumMelts);   
   m_CurStop = fabs((Emed - Ebest)/Emed);

   RestoreBest();   
   m_MeltCount++;

   //assign initial temperatue
   dEmax = dEavg + 3.00*CalcStdDev(pdE, m_NumMelts);
   m_CurTemp = m_InitTemp = 100.00*dEmax;
   m_dEavg = dEavg; //store average energy change as a metric   

   //if user-supplied reduction rate is not valid, compute internally
   if((m_TempFactor >= 1.00) || (m_TempFactor <= 0.00))
   {
      LogError(ERR_BAD_ARGS, "Invalid temperature reduction rate; using internally calculated value");
      m_TempFactor = pow((1.00/m_InitTemp), 1.00/(double)m_MaxOuter);
   }

   delete [] pdE;

   return (m_pModel->GetObjFuncVal());
}/* end Melt() */

/******************************************************************************
Optimize()

Optimize the objective function using the SA algorithm.
******************************************************************************/
void SAAlgorithm::Optimize(void)
{   
   StatusStruct pStatus;
   double curVal;
   int i;

   //write setup
   WriteSetup(m_pModel, "Simulated Annealing for Continuous Parameters");

   curVal = m_pModel->Execute();
   StoreBest();
   m_MeltCount++;

   curVal = Melt(curVal);

   //write banner and initial result
   WriteBanner(m_pModel, "iter  obj. function  ", "Convergence Value");
   WriteRecord(m_pModel, 0, curVal, m_CurStop);
   pStatus.curIter = 0;
   pStatus.maxIter = m_MaxOuter;
   pStatus.pct = 0.00;
   pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&pStatus);

   for(i = 0; i < m_MaxOuter; i++)
   {
      if(IsQuit() == true){ break;}

      curVal = Equilibrate(curVal);
      m_CurTemp *= m_TempFactor;

      //write iteration result
      WriteRecord(m_pModel, (i+1), curVal, m_CurStop);
      pStatus.curIter = m_NumOuter = i+1;
      pStatus.pct = ((float)100.00*(float)(i+1))/(float)m_MaxOuter;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);

      //check for convergence
      if(m_CurStop <= m_StopVal){pStatus.pct = 100.00; break;}

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
void SAAlgorithm::Calibrate(void)
{ 
   FILE * pFile;
   char fileName[DEF_STR_SZ];
   int id;

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
double SAAlgorithm::Equilibrate(double initVal)
{
   int m,n; 
   double bestVal, curVal, lastVal, avgVal, median;
   char c;
   ParameterGroup * pGroup;
 
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();

   if(m_Finner == NULL)
   {
      //storage for inner loop f(x) calculations (used in convergence test)
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
   m_CurStop = fabs((median - bestVal)/median);

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

The amount of parameter change induced by a transition is defined by the size
of the 'neighborhood' of adjacent solutions, which for continuous
SA is limited to +/-10% of the parameter range.
******************************************************************************/
double SAAlgorithm::Transition(double initVal)
{
   int i, n;
   double increase, prob, r, curVal;
   double upr, lwr, range, val, adj;
   ParameterGroup * pGroup;
   ParameterABC * pParam;

	//store initial Parameter values
   m_pTransBackup->Store();
  
   /*---------------------------------------
   Make a move.
   ---------------------------------------*/	
   pGroup = m_pModel->GetParamGroupPtr();
   n = pGroup->GetNumParams();
   for(i = 0; i < n; i++)
   {
      //perterb the parameter to a neighboring state (+/- 10% of range)
      pParam = pGroup->GetParamPtr(i);
      upr = pParam->GetUprBnd();
      lwr = pParam->GetLwrBnd();
      range = upr - lwr;
      if((range >= 1.00) && (range < 10.00) && 
         (strcmp(pParam->GetType(), "real") != 0)) range = 10.00;
     
      val = pParam->GetEstVal()-lwr;
      
      //determine random adjustment
      adj = (((2.00*((double)MyRand() / (double)MY_RAND_MAX)) - 1.00) * range)/10.00;
      
      //perterb parameter value
      val += adj;
      range = upr - lwr; //reset range
      if((val < 0) || (val > range)) val -= (2.00*adj); 
      
      //reversing direction should prevent bounds violations
      if(val > range){val = range; m_NumUprViols++;}
      if(val < 0)    {val = 0;     m_NumLwrViols++;}
      pParam->SetEstVal(val+lwr);
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
GenerateRandomMove()

Adjusts a single parameter by a random displacement.
******************************************************************************/
void SAAlgorithm::GenerateRandomMove(ParameterABC * pParam)
{   
   double range; //width of move set
   double val; //current parameter value
   double upr; //upper limit of current move
   double lwr; //lower limit of current move
   double adj; //random adjustment 
   
   upr = pParam->GetUprBnd();
   lwr = pParam->GetLwrBnd();
   val = pParam->GetEstVal()-lwr;
   range  = (upr - lwr);
   if((range >= 1.00) && (range < 10.00) && 
      (strcmp(pParam->GetType(), "real") != 0)) range = 10.00;
   
   /*----------------------------------------------
   Generate a random displacement
   ----------------------------------------------*/
   adj = (((2.00*((double)MyRand() / (double)MY_RAND_MAX)) - 1.00) * range)/10.00;

   //perterb parameter value
   range = upr - lwr; //reset range
   val += adj;
   if((val < 0) || (val > range)) val -= (2.00*adj); 
   
   //reversing direction should prevent bounds violations
   if(val > range){val = range; m_NumUprViols++;}
   if(val < 0)    {val = 0;     m_NumLwrViols++;}
   pParam->SetEstVal(val+lwr);
} /* end GenerateRandomMove() */

/******************************************************************************
GenerateRandomMove()

Adjusts parameters by a random displacement.
******************************************************************************/
void SAAlgorithm::GenerateRandomMove(void)
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
double SAAlgorithm::GaussTransition(double initVal)
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
SA_Program()

Calibrate the model using the SA algorithm.
******************************************************************************/
void SA_Program(int argc, StringType argv[])
{  
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model;

   NEW_PRINT("SAAlgorithm", 1);
   SAAlgorithm * SA = new SAAlgorithm(model);
   MEM_CHECK(SA);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE){ SA->Calibrate(); }
   else { SA->Optimize(); }

   delete SA;
   model->Destroy();
} /* end SA_Program() */


