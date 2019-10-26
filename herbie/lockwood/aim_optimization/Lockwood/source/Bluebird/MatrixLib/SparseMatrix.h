//SparseMatrix.h
#ifndef SPARSEMATRIX
#define SPARSEMATRIX

#include "MatrixInclude.h"

const int    BCG_max_iter=60;   //Maximum iterations for sparse BCG solver

enum BCGtestparam{RELATIVE_RESIDUAL,RELATIVE_TRANS_RESIDUAL,NORM_ERROR,MAX_ERROR};

/****************************************************
 *  Class CSparseMatrix
 *  Sparse Matrix Data Abstraction
 ***************************************************/
class CSparseMatrix{ 
	//Adapted from Press et al c++, pg.81

 private: /*------------------------------------------------*/
	
	double *sA;
	int    *ija;
	int     size;     //actual size of original square matrix
	int     nNonZero; //number of non-zero elements
  int     buffersize;
	bool   *dirichlet;  //dirichlet entries (used for dirichlet conditions)

	double *p,*pp,*r,*rr,*z,*zz;

	double CalculateAdjustedNorm(Ironclad1DArray b, int size, const normtype type) const;
  void   SolveEst   (Ironclad1DArray b, Writeable1DArray x, const int size, const bool transpose) const;

 public: /*------------------------------------------------*/

  CSparseMatrix();
	CSparseMatrix(const CSparseMatrix &A);
	CSparseMatrix(Ironclad2DArray A, const int size, const double thresh);	
 ~CSparseMatrix();

	int    GetNumEntries() const;

	bool   K_to_IJ    (const int k, int &i, int &j) const;	
	int    IJ_to_K    (const int i, const int j) const;	
	double GetAij     (const int i, const int j) const;

	int    GetNumOffDiag    (const int row) const;
	int    GetNumRows       () const;
	double GetAvgDiagCoeff  () const;
	double GetMaxCoeff      () const;

	bool   AddToAij         (const double &a, const int i, const int j);
	bool   SetAij           (const double &a, const int i, const int j);

	void   MakeDirichlet    (const int *rows, const int ndirichlet);
	void   MakeEssential1   (const int j, const int k, const double ac, const double bc  );
	void   MakeEssential2   (const CSparseMatrix *A, const int j1, const int k1, const double a, const double b);

	void   DynamicInitialize(const int N);
  void   DynamicInitialize(Ironclad1DArray Adiag,  const int N);
  void   DynamicAdd       (const double &a, const int i, const int j);

	void   ScalarMult       (const double     w);
	void   MatVectMult      (Ironclad1DArray  x, 
		                       Writeable1DArray b, 
													 const int        size, 
													 const bool       transpose,
													 const bool       ignore) const;


	void   Print            (const bool full, int sample) const;

	void   BCG              (Ironclad1DArray     b, 
						               Writeable1DArray    x,
													 const int           size,
													 const BCGtestparam  BCGtype,
													 double             &err,
													 const double        tolerance) const;
};	

#endif