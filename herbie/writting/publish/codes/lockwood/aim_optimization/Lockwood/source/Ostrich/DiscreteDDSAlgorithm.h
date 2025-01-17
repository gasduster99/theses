/******************************************************************************
File     : DiscreteDDS.h
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

The DiscreteDDS algorithm modifies to DDS algorithm so that parameter perturbations
are ensured to be at least +/- 1

Version History
03-25-10    lsm   created file
******************************************************************************/
#ifndef DDDS_ALGORITHM_H
#define DDDS_ALGORITHM_H

#include <stdio.h>
#include <math.h>
#include "AlgorithmABC.h"
#include "StatsClass.h"
#include "ModelABC.h"

/******************************************************************************
class DiscreteDDSAlgorithm
******************************************************************************/
class DiscreteDDSAlgorithm : public AlgorithmABC
{
   public:
      DiscreteDDSAlgorithm(ModelABC * pModel);
      ~DiscreteDDSAlgorithm(void);

      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);

   private:
      ModelABC   *m_pModel;  //Pointer to model being optimized
	   StatsClass *m_pStats;  //Pointer to statistics

		double  m_r_val;     //perturbation number 0<r<1
		int	  m_MaxIter;   //maximum number of iterations                                                                                            
		int	  m_UserSeed;	//random number generator seed
      //if true, then algorithm starts with users best guess (param->EstVal)
      //if false, random parameter set chosen		
		bool    m_UserSuppliedInit;		
      int     m_nCorr;  //number of times perturbations were corrected
	   double PerturbParam(const double &best_value, ParameterABC * pParam);  

}; /* end class DiscreteDDSAlgorithm */

extern "C" {
void DiscreteDDS_Program(int argC, StringType argV[]);
}

#endif /* DDDS_ALGORITHM_H */
