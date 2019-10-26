/******************************************************************************
File     : BinaryGA.h
Author   : L. Shawn Matott
Copyright: 2005, L. Shawn Matott

A Genetic Algorithm applies concepts (namely survival of the fittest and 
natural selection) from evolutionary theory to optimization problems. The 
Genetic Algorithm starts with a population of coded solutions (ChromosomePool) 
and evolves this population using the processes of Selection, Crossover and 
Mutation such that each successive genration of solutions is an improvement 
(on average) over previous generations.

Version History
10-18-05    lsm   added copyright information and initial comments.
01-01-07    lsm   Algorithm uses abstract model base class (ModelABC)
******************************************************************************/
#ifndef BINARY_GA_H
#define BINARY_GA_H

#include <stdio.h>
#include "AlgorithmABC.h"
#include "ChromosomePool.h"
#include "StatsClass.h"
#include "ModelABC.h"

/******************************************************************************
class BinaryGA

******************************************************************************/
class BinaryGA : public AlgorithmABC
{
   public:
      BinaryGA(ModelABC * pModel);
      ~BinaryGA(void);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);

   private:
      ModelABC * m_pModel;
      ChromosomePool * m_pPopulation;
      StatsClass * m_pStats;
      double m_StopVal; //convergence criteria
      double m_CurStop; //current convergence val (compared against m_StopVal)
      int m_MaxGens;
      int m_CurGen;
}; /* end class BinaryGA */

extern "C" {
void BGA_Program(int argC, StringType argV[]);
}

#endif /* BINARY_GA_H */


