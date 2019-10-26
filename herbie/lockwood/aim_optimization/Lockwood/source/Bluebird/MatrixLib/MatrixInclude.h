#ifndef MATRIXINCLUDE_H
#define MATRIXINCLUDE_H

#include "GeneralInclude.h"
#include <math.h>
#include <iostream>

//Matrix Stuff
//************************************************************************
enum normtype  {
	VECTOR_MAGNITUDE_NORM,
	LARGEST_COMPONENT_NORM
};
enum cgm_type  {
	BASIC_CGM, 
	PRECONDITIONED_CGM,
	BICONJUGATE
};
const double gauss_min  =1e-40;         //really,really small number

//library exit strategy functions (functions located in MatrixSolvers.cpp)------

void   MatExitGracefully(char *statement, badcode code);

void   MatExitGracefullyIf(bool condition, char *statement, badcode code);


/********************************************************************************
  Matrix Algebra Function declarations (functions located in MatrixSolvers.cpp)
********************************************************************************/
void   MatVectMult    (Ironclad2DArray  A, 
											 Ironclad1DArray  x,
											 Writeable1DArray b, 
											 const int size, 
											 bool transpose);
void   MatVectMult    (Ironclad2DArray  A, 
											 Ironclad1DArray  x,
											 Writeable1DArray b, 
											 const int size1, 
											 const int size2);
void   MatMult        (Ironclad2DArray  A, 
											 Ironclad2DArray  B,
											 Writeable2DArray C, 
											 const int				N, 
											 const int				M, 
											 const int				P);
bool   Gauss          (Writeable2DArray A, 
											 Writeable1DArray b,
											 Writeable1DArray x, 
											 const int        size);
void   PCG            (double         **A,
											 double          *b,
											 double          *x, 
											 const int        size, 
											 const double     tolerance, 
											 int             &numiter, 
											 const cgm_type   type);
void   ThomasAlgorithm(Ironclad1DArray  e, 
											 Ironclad1DArray  f,
											 Ironclad1DArray  g, 
											 Ironclad1DArray  b, 
											 Writeable1DArray sol,
											 const int        size);
bool   SVD            (Writeable2DArray A, 
											 Writeable1DArray b,
											 Writeable1DArray x, 
											 const int        size);
double CalculateNorm  (Ironclad1DArray  b, 
											 int              size, 
											 const normtype   type);


#endif
