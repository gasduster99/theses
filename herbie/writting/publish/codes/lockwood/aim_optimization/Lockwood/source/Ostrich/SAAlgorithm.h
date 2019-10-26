/******************************************************************************
File      : SAAlgorithm.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

An implementation of the simulated annealing algorithm.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   adjusted melt operation
                  outputs obj. func. str.
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added support for user-requested program termination
09-01-04    lsm   Algorithm now based on: Vanderbilt and Louie. 1984. "A Monte 
                  Carlo Simulated Annealing Approach to Optimization over 
                  Continuous Variables". Journal of Computational Physics. 
                  vol. 56, pg. 259-271.
11-18-04    lsm   Added convergence criteria, based on median F() of inner loop.
10-21-05    lsm   Switched to homegrown implementation that is easier to study
                  than the Vanderbilt and Louie implementation
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
******************************************************************************/
#ifndef SA_ALGORITHM_H
#define SA_ALGORITHM_H

#include "AlgorithmABC.h"
#include "ModelABC.h"
#include "ModelBackup.h"
#include <stdio.h>
#include "ParameterGroup.h"
#include "StatsClass.h"

/******************************************************************************
class SAAlgorithm

******************************************************************************/
class SAAlgorithm : public AlgorithmABC
{
   public:
      SAAlgorithm(ModelABC * pModel);
      ~SAAlgorithm(void);
      void Destroy(void);

      void Optimize(void);
      void Calibrate(void);
      void WriteMetrics(FILE * pFile);

   private :
      void GenerateRandomMove(ParameterABC * pParam);
      void GenerateRandomMove(void);
      double Transition(double initVal);
      double GaussTransition(double initVal);
      double Equilibrate(double initVal);
      double Melt(double initVal);
      void StoreBest(void);
      void RestoreBest(void);

      int m_NumOuter; //current iteration number
      int m_MaxOuter; //maximum outer iterations
      int m_MaxInner; //maximum inner iterations
      double m_Hood;  //neighborhood size
      double m_InitTemp; //Initial temperature
      double m_CurTemp;  //Current temperature
      double m_TempFactor; //Temperature reduction factor
      double m_dEavg;
      double m_StopVal;  //convergence criteria
      double m_CurStop; //current convergence val (compared against m_StopVal)
      ModelABC * m_pModel;
      ModelBackup * m_pTransBackup;
      double * m_pMelts;
      double * m_Finner;
     
      int m_NumMelts; //num. of obj. func. evals. used to determine initial temperature

      double * m_pBest;/* array containg current best parameter set */
      StatsClass * m_pStats;

      //metrics
      int m_MeltCount;
      int m_TransCount;
      int m_NumAborts;
      int m_EquilCount;
      int m_NumUprViols;
      int m_NumLwrViols;
      int m_NumUphill;
      int m_NumDownhill;
      double m_CurProb;
      double m_InitProb;
      double m_TotProb;
      int m_NumProbTests;
}; /* end class SAAlgorithm */

extern "C" {
void SA_Program(int argc, StringType argv[]);
}

#endif /* SA_ALGORITHM_H */

