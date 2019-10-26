/******************************************************************************
File     : ParticleSwarmDESC.h
Author   : L. Shawn Matott
Copyright: 2007, L. Shawn Matott

Particle Swarm Optimization (PSO) applies concepts of social behavior to
solve optimization problems. The PSO algorithm starts with a 'swarm' of 
particles (solutions) and "flies" this population through the design space in
search of the optimal solution. At each iteration, a given particle uses it's
own prior best solution (cognitive behavior) along with the current best 
solution of all particles (social behavior) to decide where to go next.

The PSODESC algorithm is a modified PSO, known as PSO using Diversity-Enhanced
Shuffled Complexes.PSODESC borrows concepts from the Shuffled Complex 
Evolutionary and Niched-Genetic algorithms to encourage efficient, effective, 
and reliable swarm exploration. After initialization, the swarm is subdivided 
into complexes using a diversity-enhancement scheme in which best-fit complex 
members are well-distributed using a niche-radius concept.  The complexes are 
allowed to independently explore the design space, and are periodically merged, 
shuffled, and reformed.

Version History
09-14-07    lsm   created from ParticleSwarm.h
******************************************************************************/

#ifndef PARTICLE_SWARM_DESC_H
#define PARTICLE_SWARM_DESC_H

#include <stdio.h>
#include "AlgorithmABC.h"
#include "StatsClass.h"
#include "QuadTree.h"
#include "ModelABC.h"
#include "ParticleSwarm.h" //for definition of ParticleStruct

/* define a structure to encapsulate each complex. */
typedef struct PSO_COMPLEX
{
   int size;
   int nfail;  //number of consective iterations for which fb did not improve
   double fbprev; //previous best objective function value   
   ParticleStruct ** pSwarm; //swarm members
}ComplexStruct;

/******************************************************************************
class ParticleSwarmDESC

******************************************************************************/
class ParticleSwarmDESC : public AlgorithmABC
{
   public:
      ParticleSwarmDESC(ModelABC * pModel);
      ~ParticleSwarmDESC(void);
	   void InitFromFile(IroncladString pFileName);
      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);      

   private:
      void ShuffleComplexes(void);
      void EvaluateSwarm(void);
      void BcastSwarm(void);
      void EvalSwarmParallel(void);
      void EvalSwarmSuperMUSE(void);
      double CalcPSOMedian(void);
      double GetParticleDistance(ParticleStruct * p1, ParticleStruct * p2);
      void SortSwarmByObjFunc(ParticleStruct ** pSwarm, int size);
      void SortSwarmByDistance(ParticleStruct ** pSwarm, int size);

      ModelABC * m_pModel;
      ParticleStruct * m_pSwarm;
      ParticleStruct ** m_pSortedSwarm;
      ComplexStruct * m_pComplex;
      StatsClass * m_pStats;
      QuadTree * m_pTrees;
      double * m_pSD; //parameter std. dev., used in distance calculations
      int m_NumComplexes; //number of complexes      
      int m_TreeSize;
      int m_SwarmSize;
      int m_MaxGens;
      int m_PolishGens; //number of 'polishing' generations, where DESC is disabled
      int m_BestIdx;
      double m_Best;
      double m_Constrict; //constriction factor
      double m_c1; //cognitive weight
      double m_c2; //social weight
      double m_Inertia; //Inertia weight
      double m_RedRate;   
      double m_ShuffleRate;  //fraction of swarm members that are shuffled
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
      double * m_Fmedian;
}; /* end class ParticleSwarm */

extern "C" {
void PSODESC_Program(int argC, StringType argV[]);
void PSODESC_LEVMAR_Program(int argC, StringType argV[]);
}

#endif /* PARTICLE_SWARM_DESC_H */


