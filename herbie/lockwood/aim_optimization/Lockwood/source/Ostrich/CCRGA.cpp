/******************************************************************************
File     : CCRGA.cpp
Author   : L. Shawn Matott
Copyright: 2009, L. Shawn Matott

A computation-constrained real-coded GA.

Version History
07-26-09    lsm   created
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "CCRGA.h"
#include "Model.h"
#include "Exception.h"
#include "Utility.h"
#include "WriteUtility.h"
#include "ParameterGroup.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
CCRGA::CCRGA(ModelABC * pModel)
{
   RegisterAlgPtr(this);

   m_pModel = pModel;

   NEW_PRINT("ChromosomePool", 1);
   m_pPopulation = new ChromosomePool;

   MEM_CHECK(m_pPopulation);
   
   m_pStats = NULL;
   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
DTOR

Free up memory used by the GA and it's member variables.
******************************************************************************/
CCRGA::~CCRGA(void)
{
   Destroy();
}/* end DTOR()*/

/******************************************************************************
Destroy()

Free up memory used by the GA and it's member variables.
******************************************************************************/
void CCRGA::Destroy(void)
{
   delete m_pPopulation;
   delete m_pStats;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using the GA.
******************************************************************************/
void CCRGA::Calibrate(void)
{ 
   int id;
   char fileName[DEF_STR_SZ];
   FILE * pFile;
      
   Optimize();

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);

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

Minimize the objective function using the GA.
******************************************************************************/
void CCRGA::Optimize(void)
{
   int numEvals = 0;
   StatusStruct pStatus;
   int maxGens;
   int i, id;
   Chromosome * pBest;
        
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //initialize sampling with replacement
   int np = m_pModel->GetParamGroupPtr()->GetNumParams();
   SampleWithReplacement(-1, np);

   m_pPopulation->CreateComm(m_pModel);
   m_pPopulation->Initialize(&m_Budget); 
   maxGens = m_pPopulation->GetNumGens();

   if(id == 0)
   {
      //write setup
      WriteSetup(m_pModel, "Computation-Constrained Real-coded Genetic Algorithm (CCRGA)");   
      //write banner
      WriteBanner(m_pModel, "gen    best fitness   ", " Pct. Complete");
   }/* end if() */

   pStatus.maxIter = m_MaxGens = maxGens;
   for(i = 0; i <= maxGens; i++)
   {
      pStatus.curIter = m_CurGen = i;      
      if(IsQuit() == true){ break;}
      if(numEvals >= m_Budget){ break;}

      //evaluate current generation
      m_pPopulation->EvalFitness();
      numEvals += m_pPopulation->GetPoolSize();

      //compute best
      pBest = m_pPopulation->GetBestFit();

      if(id == 0)
      {         
         //write results
         m_pPopulation->ConvertChromosome(pBest);
         pStatus.pct  = ((float)100.00*(float)numEvals)/(float)m_Budget;
         pStatus.numRuns = m_pModel->GetCounter();
         WriteRecord(m_pModel, i, pBest->GetFitness(), pStatus.pct);
         WriteStatus(&pStatus);

         //create next generation
         if(i < maxGens)
         {
            m_pPopulation->CreateNxtGen((double)i/(double)maxGens);
         }
      }/* end if() */

      if(numEvals >= m_Budget)
      { 
         pBest = m_pPopulation->GetBestFit();
         m_pPopulation->ConvertChromosome(pBest);
         pStatus.pct = 100.00;
         break;
      }

      //perform intermediate bookkeeping
      m_pModel->Bookkeep(false);
   }/* end for() */

   m_pModel->Execute();

   //perform final bookkeeping
   m_pModel->Bookkeep(true);

   if(id == 0)
   {   
      WriteOptimal(m_pModel, pBest->GetFitness());
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      //write algorithm metrics
      WriteAlgMetrics(this);
   }

   //free up memory used by SampleWithReplacement()
   SampleWithReplacement(-3, np);
} /* end Optimize() */

/******************************************************************************
WriteMetrics()

Write out setup and metrics for the algorithm.
******************************************************************************/
void CCRGA::WriteMetrics(FILE * pFile)
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Computation-Constrained Real-coded Genetic Algorithm (CCRGA)\n");
   fprintf(pFile, "Max Generations         : %d\n", m_MaxGens);
   fprintf(pFile, "Actual Generations      : %d\n", m_CurGen);
   fprintf(pFile, "Computational Budget    : %d\n", m_Budget);
   m_pPopulation->WriteMetrics(pFile);
   //fprintf(pFile, "Total Evals             : %d\n", m_pModel->GetCounter());
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
CCRGA_Program()

Calibrate the model using the GA.
******************************************************************************/
void CCRGA_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("CCRGA", 1);
   CCRGA * GA = new CCRGA(model);

   MEM_CHECK(GA);
   
   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { GA->Calibrate(); }
   else { GA->Optimize(); }

   delete GA;
   model->Destroy();
} /* end CCRGA_Program() */


