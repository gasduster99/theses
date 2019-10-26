/******************************************************************************
File      : PowellAlgorithm.h
Author    : L. Shawn Matott
Copyright : 2003, L. Shawn Matott

An implementation of Powell's optimization algorithm.

Version History
08-26-03    lsm   created 
02-17-04    lsm   switched over to Numerical Recipes implementation
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC

******************************************************************************/
#ifndef POWELL_ALGORITHM_H
#define POWELL_ALGORITHM_H

#include "AlgorithmABC.h"
#include "ModelABC.h"
#include <stdio.h>
#include "ParameterGroup.h"
#include "StatsClass.h"
#include "OptSearchClass.h"

/******************************************************************************
class PowellAlgorithm

Powell's method is a zero-order optimization algorithm, which utilizes the
concept of conjugate directions to determine the optimial search direction.
This implementation of Powell's Method is derived from the description given
by Press et.al. in Numerical Recipes in C, pages 415-418.
******************************************************************************/
class PowellAlgorithm : public AlgorithmABC
{
   public:
      PowellAlgorithm(ModelABC * pModel);
      ~PowellAlgorithm(void);
      void Destroy(void);
      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);

   private:
      //max. # of iterations, where one 1D search is counted as an iteration
      int m_MaxIter; 

      /* if difference between obj. function from the previous iteration 
      is less than the convergence value, the algorithm exits.*/
      double m_ConvVal; 

      int m_NumDirs; //numer of search directions

      double ** m_pSearchDirs; //array of search directions

      ModelABC * m_pModel;
      StatsClass * m_pStats; //calibration statistics
      OptSearchClass * m_pSearch; //1-dimensional search

      //metrics
      int m_AlgCount;
      int m_NumRestarts;
      int m_NumUprViols;
      int m_NumLwrViols;
      int m_CurIter;

}; /* end class PowellAlgorithm */

extern "C" {
void PWL_Program(int argc, StringType argv[]);
}

#endif /* POWELL_ALGORITHM_H */

