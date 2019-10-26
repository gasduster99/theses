/******************************************************************************
File     : CCPSO.h
Author   : L. Shawn Matott
Copyright: 2009, L. Shawn Matott

Computation Constrained Particle Swarm Optimization (CCPSO) - PSO on a budget!

Version History
07-18-09    lsm   added copyright information and initial comments.
******************************************************************************/

#ifndef CCPSO_H
#define CCPSO_H

#include <stdio.h>
#include "AlgorithmABC.h"
#include "StatsClass.h"
#include "QuadTree.h"
#include "ModelABC.h"
#include "ParticleSwarm.h"

/******************************************************************************
class CCPSO

******************************************************************************/
class CCPSO : public AlgorithmABC
{
   public:
      CCPSO(ModelABC * pModel);
      ~CCPSO(void);
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);      

   private:
      void EvaluateSwarm(void);
      void BcastSwarm(void);
      void EvalSwarmParallel(void);
      void EvalSwarmSuperMUSE(void);
      int GetParticleRank(int p);

      ModelABC * m_pModel;
      ParticleStruct * m_pSwarm;
      StatsClass * m_pStats;
      int m_Budget;
      int m_SwarmSize;
      int m_MaxGens;
      int m_BestIdx;
      double m_Best;
      double m_Constrict; //constriction factor
      double m_c1; //cognitive weight
      double m_c2; //social weight
      double m_Inertia; //Inertia weight
      double m_RedRate;   
      PopInitType m_InitType;      
      bool m_LinRedFlag; //true: linearly reduce interia to zero
      int m_CurGen;
      double m_StopVal;  //convergence criteria
      double m_CurStop; //current convergence val (compared against m_StopVal)

      //buffers used in MPI-parallel communication
      double * m_pBuf;
      double * m_pMyBuf;
      double * m_pTmpBuf;
      double * m_pBigBuf;

      //buffer for initial parameter values
      int       m_NumInit;
      double ** m_pInit;

      //metrics
      int m_NumUprViols;
      int m_NumLwrViols;
}; /* end class CCPSO */

extern "C" {
void CCPSO_Program(int argC, StringType argV[]);
void CCPSO_LEVMAR_Program(int argC, StringType argV[]);
}

#endif /* CCPSO_H */


