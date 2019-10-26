/******************************************************************************
File     : GridAlgorithm.h
Author   : L. Shawn Matott
Copyright: 2005, L. Shawn Matott

This algorithm will "grid" the objective function surface by evaluating an
exhaustive (or nearly exhaustive) set of parameter configurations.

Version History
01-25-05    lsm   created
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
017-18-07   lsm   Added SuperMUSE support
******************************************************************************/
#ifndef GRID_ALG_H
#define GRID_ALG_H

#include <stdio.h>
#include "AlgorithmABC.h"
#include "StatsClass.h"
#include "ModelABC.h"
#include "MyTypes.h"

/******************************************************************************
class GridAlgorithm

******************************************************************************/
class GridAlgorithm : public AlgorithmABC
{
   public:
      GridAlgorithm(ModelABC * pModel);
      ~GridAlgorithm(void);
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);

   private:
      void EvaluateGrid(void);
      void BcastGrid(void);
      void EvalGridParallel(void);
      void EvalGridSuperMUSE(void);
      double GetGridVal(int i, int j);

      ModelABC * m_pModel;
      StatsClass * m_pStats;
      GridStruct * m_pMini;

      int * m_pDims; //grid dimensions
      double * m_Lwr; //parameter lower bounds
      int * m_Rval;  //rollover counts
      int m_GridSize;
      int m_NumIters;
      int m_NumLeft;
      int m_CurIter;
      int m_MiniSize;      
      double * m_pBest;

      //buffers used in MPI-parallel communication
      double * m_pBuf;
      double * m_pMyBuf;
      double * m_pTmpBuf;
      double * m_pBigBuf;
}; /* end class GridAlgorithm */

extern "C" {
void GRID_Program(int argC, StringType argV[]);
}

#endif /* GRID_ALG_H */


