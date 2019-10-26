/******************************************************************************
File     : CCRGA.h
Author   : L. Shawn Matott
Copyright: 2009, L. Shawn Matott

A computation-constrained real-coded GA.

Version History
07-26-09    lsm   created
******************************************************************************/
#ifndef CCRGA_H
#define CCRGA_H

#include <stdio.h>
#include "AlgorithmABC.h"
#include "ChromosomePool.h"
#include "StatsClass.h"
#include "ModelABC.h"

/******************************************************************************
class CCRGA

******************************************************************************/
class CCRGA : public AlgorithmABC
{
   public:
      CCRGA(ModelABC * pModel);
      ~CCRGA(void);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);

   private:
      ModelABC * m_pModel;
      ChromosomePool * m_pPopulation;
      StatsClass * m_pStats;
      int m_Budget;
      int m_MaxGens;
      int m_CurGen;
}; /* end class CCRGA */

extern "C" {
void CCRGA_Program(int argC, StringType argV[]);
}

#endif /* CCRGA_H */


