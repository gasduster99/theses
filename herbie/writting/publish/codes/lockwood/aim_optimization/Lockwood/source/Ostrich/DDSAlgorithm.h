/******************************************************************************
File     : DDS.h
Author   : James R. Craig and Bryan Tolson
Copyright: 2006, James R. Craig and Bryan Tolson

The DDS algorithm uses an intelligent greedy search to search the parameter 
space and find the global minimum of an optimization problem. The algorithm is 
discussed in Tolson and Shoemaker, Water Resources Research, 2007

Dynamically dimensioned Search (DDS) version 1.1 algorithm by Bryan Tolson
c++ version (original was coded in Matlab, translated to Fortran by Bryan Tolson)
Translated to c++ by James Craig (July-Aug 2006)
DDS is an n-dimensional continuous global optimization algorithm

Ostrich updates required in following files:
	MyTypes.h
	Ostrich.cpp
	Utilities.h
	Utilities.cpp

Version History
03-01-06    jrc   created file
******************************************************************************/
#ifndef DDS_ALGORITHM_H
#define DDS_ALGORITHM_H

#include <stdio.h>
#include <math.h>
#include "AlgorithmABC.h"
#include "StatsClass.h"
#include "ModelABC.h"

/******************************************************************************
class DDSAlgorithm
******************************************************************************/
class DDSAlgorithm : public AlgorithmABC
{
   public:
      DDSAlgorithm(ModelABC * pModel);
      ~DDSAlgorithm(void);

      void Optimize(void);
      void Calibrate(void);
      void Destroy(void);
      void WriteMetrics(FILE * pFile);

   private:
      ModelABC   *m_pModel;             //Pointer to model being optimized
			StatsClass *m_pStats;             //Pointer to statistics

			double			m_r_val;							//perturbation number 0<r<1
			int					m_MaxIter;						//maximum number of iterations                                                                                            
			int					m_UserSeed;					  //random number generator seed
			bool				m_UserSuppliedInit;		//if true, then algorithm starts with users best guess (param->EstVal)
																				//if false, random parameter set chosen		

			double PerturbParam(const double &best_value, ParameterABC * pParam);  

}; /* end class DDSAlgorithm */

extern "C" {
void DDS_Program(int argC, StringType argV[]);
}

#endif /* DDS_ALGORITHM_H */
