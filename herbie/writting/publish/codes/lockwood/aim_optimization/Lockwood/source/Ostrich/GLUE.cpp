/******************************************************************************
File     : GLUE.cpp
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Generalized Likelihood Uncertainty Engine - GLUE

Version History
06-23-10    lsm   added copyright information and initial comments.
******************************************************************************/
#include "mpi_stub.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "GLUE.h"
#include "Model.h"
#include "Exception.h"
#include "WriteUtility.h"
#include "Utility.h"
#include "ParameterGroup.h"
#include "ParameterABC.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
GLUE::GLUE(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pBehavioral = NULL;
   m_pSamples = NULL;
   m_MaxSamples = 0;
   m_NumDesired = 0;
   m_NumFound = 0;
   m_SamplesPerIter = 0;
   m_Threshold = -1.00;

   //MPI-parallel communication arrays
   m_pMyBuf  = NULL;
   m_pTmpBuf = NULL;
   m_pBigBuf = NULL;
   m_pBuf    = NULL;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
DTOR

Free up memory used by GLUE and it's member variables.
******************************************************************************/
GLUE::~GLUE(void)
{
   Destroy();
}/* end DTOR()*/

/******************************************************************************
Destroy()

Free up memory used by the GLUE and it's member variables.
******************************************************************************/
void GLUE::Destroy(void)
{
   int i;
   for(i = 0; i < m_NumDesired; i++)
   {
      delete [] m_pBehavioral[i].x;
   }
   delete [] m_pBehavioral;

   for(i = 0; i < m_SamplesPerIter; i++)
   {
      delete [] m_pSamples[i].x;
   }
   delete [] m_pSamples;
   
   m_MaxSamples = 0;
   m_NumDesired = 0;
   m_NumFound = 0;
   m_SamplesPerIter = 0;
   m_Threshold = -1.00;

   delete [] m_pMyBuf;
   delete [] m_pTmpBuf;
   delete [] m_pBigBuf;
   delete [] m_pBuf;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using GLUE.
******************************************************************************/
void GLUE::Calibrate(void)
{    
   Optimize();
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using GLUE.
******************************************************************************/
void GLUE::Optimize(void)
{
   int num,id,i,j,k,m,g, numSamples;
   double range, r, rval, upr, lwr;
   int maxGens;
   StatusStruct pStatus;
   ParameterGroup * pGroup;
   
   InitFromFile(GetInFileName());

   maxGens = 1+(m_MaxSamples/m_SamplesPerIter);

   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if(id == 0)
   {
      WriteSetup(m_pModel, "Generalized Likelihood Uncertainty Engine");
      //write banner
      WriteBanner(m_pModel, "gen   best value     ", "Num Found");
   }/* end if() */

   //allocate list of behavioral samples
   pGroup = m_pModel->GetParamGroupPtr();
   num = pGroup->GetNumParams();

   NEW_PRINT("SampleStruct", m_NumDesired);
   m_pBehavioral = new SampleStruct[m_NumDesired];
   MEM_CHECK(m_pBehavioral);

   for(i = 0; i < m_NumDesired; i++)
   {
      NEW_PRINT("double", num);
      m_pBehavioral[i].x = new double[num];
      m_pBehavioral[i].fx = HUGE_VAL;
      m_pBehavioral[i].n = num;      
   }/* end for() */
   MEM_CHECK(&(m_pBehavioral[i-1]));

   //allocate list of random samples
   NEW_PRINT("SampleStruct", m_SamplesPerIter);
   m_pSamples = new SampleStruct[m_SamplesPerIter];
   MEM_CHECK(m_pSamples);

   for(i = 0; i < m_SamplesPerIter; i++)
   {
      NEW_PRINT("double", num);
      m_pSamples[i].x = new double[num];
      m_pSamples[i].fx = HUGE_VAL;
      m_pSamples[i].n = num;      
   }/* end for() */
   MEM_CHECK(&(m_pSamples[i-1]));
 
   //main optimization loop   
   for(g = 0; g < maxGens; g++)
   {
      pStatus.curIter = g+1;
      if(IsQuit() == true){ break;}
      if(m_NumFound == m_NumDesired){ pStatus.pct = 100.00; break;}
      numSamples = g*m_SamplesPerIter;
      if(numSamples >= m_MaxSamples){ pStatus.pct = 100.00; break;} 

      if(id == 0)
      {
         /*-----------------------------------------------------------
         Generate random samples
         -----------------------------------------------------------*/
         for(i = 0; i < m_SamplesPerIter;  i++)
         {
            for(j = 0; j < num; j++) //for each parameter
            {
               lwr = pGroup->GetParamPtr(j)->GetLwrBnd();
               upr = pGroup->GetParamPtr(j)->GetUprBnd();

               range = upr - lwr;
               r = (double)MyRand() / (double)MY_RAND_MAX;
               rval = (r * range) + lwr;
               m_pSamples[i].x[j] = rval;               
            }/* end for() */
         }/* end for() */
      } /* end if() */

      //evaluate samples, possibly in parallel
      EvaluateSamples();
   
      //revise list of best samples
      for(i = 0; i < m_SamplesPerIter; i++)
      {
         for(j = 0; j < m_NumDesired; j++)
         {
            if(m_pSamples[i].fx < m_pBehavioral[j].fx)
            {
               //shift previous list to make room
               for(k = m_NumDesired-1; k > j; k--)
               {
                  for(m = 0; m < num; m++)
                  {
                     m_pBehavioral[k].x[m] = m_pBehavioral[k-1].x[m];
                     m_pBehavioral[k].fx = m_pBehavioral[k-1].fx;
                  }
               }
               //insert new entry
               for(m = 0; m < num; m++)
               {
                  m_pBehavioral[j].x[m] = m_pSamples[i].x[m];
                  m_pBehavioral[j].fx = m_pSamples[i].fx;
               }
               break;
            }/* end if() */
         }/* end for() */
      }/* end for() */

      //count number of behavioral
      m_NumFound = 0;
      for(j = 0; j < m_NumDesired; j++)
      {
         if(m_pBehavioral[j].fx < m_Threshold)
         {                      
            m_NumFound++;
         }/* end if() */
      }/* end for() */
      
      //first entry is always the best
      pGroup->WriteParams(m_pBehavioral[0].x);

      pStatus.pct = ((float)100.00*(float)(g+1))/(float)maxGens;
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      if(id == 0){ WriteRecord(m_pModel, (g+1), m_pBehavioral[0].fx, m_NumFound);}
   }/* end for() */

   //place model at optimal prameter set
   pGroup->WriteParams(m_pBehavioral[0].x);
   m_pModel->Execute();

   if(id == 0)
   { 
      WriteOptimal(m_pModel, m_pBehavioral[0].fx);
      pStatus.numRuns = m_pModel->GetCounter();
      WriteStatus(&pStatus);
      //write algorithm metrics
      WriteAlgMetrics(this);
   }
} /* end Optimize() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void GLUE::WriteMetrics(FILE * pFile) 
{
   ParameterGroup * pGroup;
   pGroup = m_pModel->GetParamGroupPtr();

   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm               : Generalized Likelihood Uncertainty Estimation\n");
   fprintf(pFile, "Behavioral Threshold    : %E\n", m_Threshold);
   fprintf(pFile, "Max Samples             : %d\n", m_MaxSamples);
   fprintf(pFile, "Actual Num. Behavorial  : %d\n", m_NumFound);
   fprintf(pFile, "Desired Num. Behavorial : %d\n", m_NumDesired);

   fprintf(pFile,"Sample  obj.function  ");
   pGroup->Write(pFile, WRITE_BNR);
   fprintf(pFile,"\n");

   for(int i = 0; i < m_NumDesired; i++)
   {
	   fprintf(pFile, "%-4d  ", i);
      WritePreciseNumber(pFile, m_pBehavioral[i].fx);
      fprintf(pFile, "  ");
      pGroup->WriteParams(m_pBehavioral[i].x);
      pGroup->Write(pFile, WRITE_SCI);
      fprintf(pFile, "\n");
   }
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
EvaluateSamples()

Evaluates the objective function of each sample.
******************************************************************************/
void GLUE::EvaluateSamples(void)
{   
   int i, n;   
   ParameterGroup * pGroup;
   double val;

   MPI_Comm_size(MPI_COMM_WORLD, &n);
   
   if(n == 1) //serial execution
   {
      WriteInnerEval(WRITE_GLUE, m_SamplesPerIter, '.');
      pGroup = m_pModel->GetParamGroupPtr();
      for(i = 0; i < m_SamplesPerIter; i++) 
      { 
         WriteInnerEval(i+1, m_SamplesPerIter, '.');
         pGroup->WriteParams(m_pSamples[i].x);

         val = m_pModel->Execute();
         m_pSamples[i].fx = val;
      }
      WriteInnerEval(WRITE_ENDED, m_SamplesPerIter, '.');
   }/* end if() */
   else /* parallel execution */
   {
      BcastSamples();
      EvalSamplesParallel();      
   }/* end else() */
} /* end EvaluateSamples() */

/******************************************************************************
BcastSamples()

When in parallel, only the master computes the samples. All the other 
processors just compute the objeective functions. The BcastSamples() routine is 
called upon to broadcast  the current set of samples from the master processor 
to all of the slave processors.
******************************************************************************/
void GLUE::BcastSamples(void)
{
   int num_vars, pop_size, buf_size;
   int i, j, num_procs, id, idx;

   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   //size the flattened variable matrix
   pop_size = m_SamplesPerIter;
   num_vars = m_pSamples[0].n;

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
         m_pBuf[idx] = m_pSamples[j].x[i];
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
         m_pSamples[j].x[i] = m_pBuf[idx];
      }/* end for() */
   }/* end for() */
}/* end BcastSamples() */

/******************************************************************************
EvalSamplesParallel()

Compute objective function of entire set of samples in parallel. Each processor 
evaluates a predetermined number of samples, based on their processor id.
******************************************************************************/
void GLUE::EvalSamplesParallel(void)
{    
   int i ,j, num_procs, id, bufsize, idx;
   ParameterGroup * pGroup;

   //setup processor id and number of processors
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
   bufsize = (m_SamplesPerIter/num_procs) + 1;

   //allocate space for intermediate buffers, if necessary
   if(m_pMyBuf == NULL)
   {
      NEW_PRINT("double", bufsize);
      m_pMyBuf = new double[bufsize];

      NEW_PRINT("double", bufsize);
      m_pTmpBuf = new double[bufsize];

      NEW_PRINT("double", m_SamplesPerIter);
      m_pBigBuf = new double[m_SamplesPerIter];
      MEM_CHECK(m_pBigBuf);
   }

   //perform parallel evaluations
   j = 0;
   pGroup = m_pModel->GetParamGroupPtr();
   for(i = 0; i < m_SamplesPerIter; i++) 
   { 
      if((i % num_procs) == id)
      {          
         pGroup->WriteParams(m_pSamples[i].x);

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

         if(idx < m_SamplesPerIter)
         {
            m_pBigBuf[idx] = m_pMyBuf[j]; //gather into bigbuf
            m_pMyBuf[j] = m_pTmpBuf[j]; //restore myBuf...clobbered by bcast
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //stuff results into swarm
   for(i = 0; i < m_SamplesPerIter; i++)
   {
      m_pSamples[i].fx = m_pBigBuf[i];
   }/* end for() */
}/* end EvalSamplesParallel() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename.
******************************************************************************/
void GLUE::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   char * line;
   char tmp[DEF_STR_SZ];

   m_MaxSamples = 100;
   m_NumDesired = 10;
   m_NumFound = 0;
   m_SamplesPerIter = 10;
   m_Threshold = 1000.00;

   //read in GLUE configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open GLUE config. file. Using Defaults");      
      return;
   }/* end if() */   

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginGLUE", pFileName) == true)
   {
      FindToken(pFile, "EndGLUE", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginGLUE", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndGLUE") == NULL)
      {         
         if(strstr(line, "SamplesPerIter") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_SamplesPerIter); 
            if(m_SamplesPerIter < 1)
            {
               LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 10.");
               m_SamplesPerIter = 10;
            }
         }/*end if() */         
         else if(strstr(line, "NumBehavioral") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumDesired); 
            if(m_NumDesired < 1)
            {
               LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 10.");
               m_NumDesired = 10;
            }
         }/*end else if() */         
         else if(strstr(line, "MaxSamples") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MaxSamples); 
            if(m_NumDesired < 1)
            {
               LogError(ERR_FILE_IO, "Invalid GLUE setting. Defaulting to 100.");
               m_MaxSamples = 100;
            }
         }/*end else if() */         
         else if(strstr(line, "Threshold") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Threshold); 
         }/*end else if() */         
         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   
   fclose(pFile);
} /* end InitFromFile() */

/******************************************************************************
GLUE_Program()

Calibrate or optimize the model using GLUE.
******************************************************************************/
void GLUE_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("GLUE", 1);
   GLUE * pGLUE = new GLUE(model);
   MEM_CHECK(pGLUE);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { pGLUE->Calibrate(); }
   else { pGLUE->Optimize(); }

   delete pGLUE;
   model->Destroy();
} /* end GLUE_Program() */

