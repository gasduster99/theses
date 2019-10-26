/******************************************************************************
File      : CCSA.h
Author    : L. Shawn Matott
Copyright : 2009, L. Shawn Matott

An implementation of a computation-constrained simulated annealing algorithm.

Version History
10-30-09    lsm   created from "standard" SA
******************************************************************************/
#ifndef CCSA_ALGORITHM_H
#define CCSA_ALGORITHM_H

#include "AlgorithmABC.h"
#include "ModelABC.h"
#include "ModelBackup.h"
#include <stdio.h>
#include "ParameterGroup.h"
#include "StatsClass.h"

/******************************************************************************
class CCSA

******************************************************************************/
class CCSA : public AlgorithmABC
{
   public:
      CCSA(ModelABC * pModel);
      ~CCSA(void);
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
      int m_Budget;
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
}; /* end class CCSA */

extern "C" {
void CCSA_Program(int argc, StringType argv[]);
}

#endif /* CCSA_ALGORITHM_H */

