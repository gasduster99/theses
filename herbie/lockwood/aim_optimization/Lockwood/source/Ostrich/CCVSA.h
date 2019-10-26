/******************************************************************************
File      : CCVSA.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

An implementation of a computation-constrained simulated annealing algorithm linked
with the Vanderbilt and Louie formulation.

Version History
10-31-09    lsm   Created
******************************************************************************/
#ifndef CC_VAND_SA_H
#define CC_VAND_SA_H

#include "AlgorithmABC.h"
#include "ModelABC.h"
#include "ModelBackup.h"
#include <stdio.h>
#include "ParameterGroup.h"
#include "StatsClass.h"

/******************************************************************************
class CCVSA

******************************************************************************/
class CCVSA : public AlgorithmABC
{
   public:
      CCVSA(ModelABC * pModel);
      ~CCVSA(void);
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
      int m_Budget;
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
}; /* end class CCVSA */

extern "C" {
void CCVSA_Program(int argc, StringType argv[]);
}

#endif /* CC_VAND_SA_H */

