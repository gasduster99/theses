/******************************************************************************
File     : GLUE.h
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Generalized Likelihood Uncertainty Engine - GLUE

Version History
06-23-10    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef GLUE_H
#define GLUE_H

#include <stdio.h>
#include "AlgorithmABC.h"
#include "ModelABC.h"

/* define a structure to encapsulate each behavioral sample. */
typedef struct GLUE_SAMPLE_STRUCT
{
	double * x; //current parameters
	double fx;  //objective function value
   int n;      //number of parameters
}SampleStruct;

/******************************************************************************
class GLUE

******************************************************************************/
class GLUE : public AlgorithmABC
{
   public:
      GLUE(ModelABC * pModel);
      ~GLUE(void);
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);      

   private:
      void EvaluateSamples(void);
      void BcastSamples(void);
      void EvalSamplesParallel(void);

      ModelABC * m_pModel;
      SampleStruct * m_pSamples;
      SampleStruct * m_pBehavioral;
      int m_MaxSamples;
      int m_NumDesired;
      int m_NumFound;
      int m_SamplesPerIter;
      double m_Threshold;

      //buffers used in MPI-parallel communication
      double * m_pBuf;
      double * m_pMyBuf;
      double * m_pTmpBuf;
      double * m_pBigBuf;
}; /* end class GLUE */

extern "C" {
void GLUE_Program(int argC, StringType argV[]);
}

#endif /* GLUE_H */


