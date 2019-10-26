/******************************************************************************
File      : VandSA.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

An implementation of the simulated annealing algorithm.

Version History
09-01-04    lsm   Algorithm based on: Vanderbilt and Louie. 1984. "A Monte 
                  Carlo Simulated Annealing Approach to Optimization over 
                  Continuous Variables". Journal of Computational Physics. 
                  vol. 56, pg. 259-271.
01-01-07    lsm   Algorithm now uses abstract model base class (ModelABC).
******************************************************************************/
#ifndef VAND_SA_H
#define VAND_SA_H

#include "AlgorithmABC.h"
#include "ModelABC.h"
#include "ModelBackup.h"
#include <stdio.h>
#include "ParameterGroup.h"
#include "StatsClass.h"

/******************************************************************************
class VandSA

******************************************************************************/
class VandSA : public AlgorithmABC
{
   public:
      VandSA(ModelABC * pModel);
      ~VandSA(void);
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
      double m_InitTemp; //Initial temperature
      double m_CurTemp;  //Current temperature
      double m_TempFactor; //Temperature reduction factor
      double m_StopVal;  //convergence criteria
      double m_CurStop; //current convergence val (compared against m_StopVal)
      ModelABC * m_pModel;
      ModelBackup * m_pTransBackup;
      double * m_pMelts;
      double * m_Finner;

      //matrices and vectors needed for implementation of Vanderbilt and Louie (1984) 
      double *  m_dx;  //transition vector (delta x)
      double ** m_Q;  //step distribution matrix
      double ** m_QT; //Q transpose
      double *  m_u;   //random steps
      double ** m_cov; //covariance matrix ('s')
      double ** m_Shape; //random walk shape matrix ('S')
      double ** m_x; //parameter values used in transition
      double *  m_A; //avg parameter value during a series of transitions

      int m_NumMelts; //num. of obj. func. evals. used to determine initial temperature

      double * m_pTransPoint;/* transition parameter set */
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
}; /* end class VandSA */

extern "C" {
void VSA_Program(int argc, StringType argv[]);
}

#endif /* VAND_SA_H */

