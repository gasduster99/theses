//SparseMatrix.cpp
#include "SparseMatrix.h"

//-----------------------------------------------------------------------
CSparseMatrix::CSparseMatrix(){
	sA      =NULL;
	ija     =NULL; 
	nNonZero=0;
	size    =0;	
	buffersize=0;
	dirichlet=NULL;
	p =NULL;
	pp=NULL;
	r =NULL;
	rr=NULL;
	z =NULL;
	zz=NULL;
}
//-----------------------------------------------------------------------
CSparseMatrix::CSparseMatrix(const CSparseMatrix &A){
	
	int j,k;
	this->size    =A.size;
	this->nNonZero=A.nNonZero;
	this->buffersize=this->nNonZero+1;
	//Allocate memory
	this->sA =new double [this->nNonZero+1]; 
	this->ija=new int    [this->nNonZero+1];

	for (k=0; k<=nNonZero;k++){
		this->sA [k]=A.sA [k];
		this->ija[k]=A.ija[k];
	}
	p        =new double [size];
	pp       =new double [size];
	r        =new double [size];
	rr       =new double [size];
	z        =new double [size];
	zz       =new double [size];
	dirichlet=new bool   [size];
	for (j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
//-----------------------------------------------------------------------
CSparseMatrix::~CSparseMatrix(){
	if (globaldebug){cout<<"DELETING SPARSE MATRIX"<<endl;}
	delete [] sA;
	delete [] ija; 
	nNonZero=0;
	size    =0;
	buffersize=0;
	delete [] p;
	delete [] pp;
	delete [] r;
	delete [] rr;
	delete [] z;
	delete [] zz;
	delete [] dirichlet;
}
/************************************************************************
	Translates A[N][N] to sparse matrix S
------------------------------------------------------------------------*/
CSparseMatrix::CSparseMatrix(Ironclad2DArray A, const int N, const double thresh){
	int i,j,k;
	size    =N;
	nNonZero=N;
	//count non-zero elements
	for (i=0; i<size; i++){
		for (j=0;j<size; j++){
			if ((fabs(A[i][j])>=thresh) && (i!=j)){nNonZero++;}
		}
	}
	//Allocate memory
	sA =new double [nNonZero+1]; 
	ija=new int    [nNonZero+1];
	if (ija==NULL){MatExitGracefully("CSparseMatrix::Constructor::Out of memory",OUT_OF_MEMORY);}
	buffersize=nNonZero+1;

	ija[0]=size+1;//index to row 0 of off-diagonal element
	for (j=0;j<size; j++){sA[j]=A[j][j];}
	k=size;
	for (i=0; i<size; i++){
		for (j=0;j<size; j++){
			if ((fabs(A[i][j])>thresh) && (i!=j)){
				if (k+1>nNonZero){MatExitGracefully("sprsTranslate: bad index calc",RUNTIME_ERR);}
				k++;
				//++k;
				sA [k]=A[i][j];
				ija[k]=j;	
			}
		}
		if (i>nNonZero){MatExitGracefully("sprsTranslate: bad index calc",RUNTIME_ERR);}
		ija[i+1]=k+1;
	}
	p       =new double [size];
	pp      =new double [size];
	r       =new double [size];
	rr      =new double [size];
	z       =new double [size];
	zz      =new double [size];
	dirichlet=new bool   [size];
	for (j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
//-----------------------------------------------------------------------
int CSparseMatrix::GetNumEntries() const{return nNonZero;}
/************************************************************************
 DynamicInitialize:
	Initializes dynamic sparse array with zero diagonal values 
------------------------------------------------------------------------*/
void CSparseMatrix::DynamicInitialize(const int N){
	if (size!=0) {
		MatExitGracefully("CSparseMatrix::DynamicInitialize::Already initialized",RUNTIME_ERR);}
	size    =N;
	nNonZero=N;
	//Allocate memory
	sA =new double [size+1]; 
	ija=new int    [size+1];
	if (ija==NULL){
		MatExitGracefully("CSparseMatrix::DynamicInitialize::Out of memory",OUT_OF_MEMORY);}
	buffersize=size+1;
	ija[0]=size+1;
	for (int i=0;i<size; i++){
		sA[i]   =0.0;
		ija[i+1]=size+1;
	}
	sA[size]=0.0; //Arbitrary value (doesn't matter)
	p       =new double [size];
	pp      =new double [size];
	r       =new double [size];
	rr      =new double [size];
	z       =new double [size];
	zz      =new double [size];
  dirichlet=new bool   [size];
	for (int j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
/************************************************************************
 DynamicInitialize:
	Initializes dynamic sparse array with diagonal values only
------------------------------------------------------------------------*/
void CSparseMatrix::DynamicInitialize(Ironclad1DArray Adiag,
																      const int       N){
	if (size!=0) {MatExitGracefully("CSparseMatrix::DynamicInitialize::Already initialized",RUNTIME_ERR);}
	size    =N;
	nNonZero=N;
	//Allocate memory
	sA =new double [size+1]; 
	ija=new int    [size+1];
	if (ija==NULL){MatExitGracefully("CSparseMatrix::DynamicInitialize::Out of memory",OUT_OF_MEMORY);}
	buffersize=size+1;

	ija[0]=size+1;
	for (int i=0;i<size; i++){
		sA [i]  =Adiag[i];
		ija[i+1]=size+1;
	}
	sA[size]=0.00001; //Arbitrary value (doesn't matter)
	p =new double [size];
	pp=new double [size];
	r =new double [size];
	rr=new double [size];
	z =new double [size];
	zz=new double [size];
  dirichlet=new bool [size];
	for (int j=0; j<size;j++){
		p[j]=pp[j]=r[j]=rr[j]=z[j]=zz[j]=0.0;
		dirichlet[j]=false;
	}
}
/************************************************************************
 DynamicAdd:
	Dynamically inserts new coefficient row=i, col=j, into sparse matrix
	could get costly with large matrices;
	but is helped with buffering of the matrices-rebuffering is called 4 times for each bandwidth increase
------------------------------------------------------------------------*/
void CSparseMatrix::DynamicAdd(const double &a, const int i, const int j){
	
	//messes up if j>i ( I dont think this is true any longer)
	int buffer=(size/4);

	int    *tmpija;
	double *tmpsA;
	int     knew, i1,k;
	
	if (size==0)                  {MatExitGracefully("CSparseMatrix::DynamicAdd::not initialized correctly",RUNTIME_ERR);}
	if ((sA==NULL) || (ija==NULL)){MatExitGracefully("CSparseMatrix::DynamicAdd: not initialized correctly",RUNTIME_ERR);}
	if ((i>=size)  || (j>=size))  {MatExitGracefully("CSparseMatrix::DynamicAdd: bad indices specified",RUNTIME_ERR);}
	if ((i<0)      || (j<0))      {MatExitGracefully("CSparseMatrix::DynamicAdd: bad indices specified(2)",RUNTIME_ERR);}
	
	if (i==j){sA[i]=a; return;}
  
	//find k index of inserted coefficient--------------------------------
	knew=ija[i+1];                           //default- guess last coefficient in row i
	for (k=ija[i]; k<ija[i+1]; k++){           //searches through row i
		if      (ija[k]==j){sA[k]=a;return;}     //if already exists, merely replace coefficient (same as SetAij)
		else if (ija[k]> j){knew=k; break;}      //otherwise, adopt knew
	}
	if (a==0.0){return;}                     //do not dynamically add if a zero entry (should this be here?}

	//copy all terms in array and insert---------------------------------- 
	nNonZero++;                              //item added to total array

	///cout <<"added ("<<i<<","<<j<<")"<<endl;
	//if (i>20){MatExitGracefully("CSparseMatrix::DynamicAdd: TMP DEBUG",RUNTIME_ERR);}
	if (knew>nNonZero){cout << "bad knew"<<endl;}

	if (buffersize<=nNonZero){
		//cout<<"*"<<endl;
		tmpija=new int    [nNonZero+1+buffer];				
		tmpsA =new double [nNonZero+1+buffer];
		if (tmpsA==NULL){
			MatExitGracefully("CSparseMatrix::DynamicAdd::Out of memory",OUT_OF_MEMORY);}
		buffersize=nNonZero+1+buffer;

		for (k=0; k<knew; k++){               //before inserted coefficient 
			tmpija[k]=ija[k];
			tmpsA [k]= sA[k];
		}	             
		for (k=nNonZero; k>=knew+1; k--){     //after inserted coefficient
			tmpija[k]=ija[k-1];
			tmpsA [k]= sA[k-1];
		}
		tmpija[knew]=j;                       //inserted coefficient
		tmpsA [knew]=a;

		delete [] ija;                        //delete old values
		delete [] sA;

		ija=tmpija;                           //redirect array pointers
		sA =tmpsA;
	}
	else{  //no need to update buffer

		for (k=nNonZero; k>=knew+1; k--){     //after inserted coefficient //still pretty slow!
			ija[k]=ija[k-1];
			sA [k]= sA[k-1];
		}
		ija[knew]=j;                          //inserted coefficient
		sA [knew]=a;
	}

	//cout << "inserted coeff k="<<knew<< ". Non zero entries="<<nNonZero<<endl; 
	//TMP DEBUG
	//for (k=0; k<=nNonZero; k++){
	//	cout << ija[k]<<"     "<<sA[k]<<endl;
	//	if (k==size){cout<<"---"<<endl;}
	//}
	//cout <<endl;

	for (i1=i; i1<size; i1++){ //update indexing information
		ija[i1+1]++;             //item added to row i (ija[i+1] is index of lowest term in row i+1)
	}


	return;
}
/************************************************************************
 MakeDirichlet:
	marks all rows and colums specified by rows[] as dirichlet 
	affects BCG, MatVectMultiply, SolveEst, Print and CalculateAdjustedNorm
------------------------------------------------------------------------*/
void CSparseMatrix::MakeDirichlet(const int *rows, const int ndirichlet){
	for (int k=0; k<size; k++){
		dirichlet[k]=false;
	}
	for (int j=0; j<ndirichlet; j++){
		dirichlet[rows[j]]=true;
	}
}
/************************************************************************
 MakeEssential1:
	revises matrix to meet essential conditions X_{j1}+a*X_{k1}=b
	affects matrix only- if condition is removed, matrix must be remade!!!!
	current implementation is pretty slow- CAN OPTIMIZE!!
	revised from Akin, pg. 403
------------------------------------------------------------------------*/
void CSparseMatrix::MakeEssential1(const int j1, const int k1, const double a, const double b){
	int i,j;
	if ((j1<0)      || (k1<0))      {
		MatExitGracefully("CSparseMatrix::MakeEssential1: bad indices specified",RUNTIME_ERR);}
	if ((j1==k1) && (a!=0.0)) {return;} //cannot force constraint on itself
  
	double Sjj=this->GetAij(j1,j1);
	double *Skrow=new double [this->size];
	double *Sjrow=new double [this->size];
	double *Sij=new double [this->size];
	double *Sik=new double [this->size];
	for (i=0; i<size; i++){
		Skrow[i]=this->GetAij(k1,i);
		Sjrow[i]=this->GetAij(j1,i);
    Sij[i]=this->GetAij(i,j1);
		Sik[i]=this->GetAij(i,k1);
	}

	for (i=0; i<size; i++){ //for each row in matrix
		
		if ((i==k1) && (a!=0.0)){ //special row (constraint for X_k1)
			for (j=0; j<size; j++){ //for each column
				double Sk1=Skrow[j];//this->GetAij(k1,j);
				double Sj1=Sjrow[j];//this->GetAij(j1,j);

				if      (j==j1){this->DynamicAdd(a,k1,j1);}
				else if (j==k1){this->DynamicAdd(Skrow[k1]+2.0*a*(Sjrow[k1]-Sjj)+a*a,k1,k1);}//
				else           {this->DynamicAdd(Skrow[j]-a*Sjrow[j],k1,j);}
			}
		}

		else if (i==j1){ //special row (constraint for X_j1)
			for (j=0; j<size; j++){ //for each column
				this->SetAij(0.0,i,j); //does nothing for most sparse matrix rows
				if (j==j1){this->SetAij(1.0,j1,j1);}
				if (a!=0.0){if (j==k1){this->DynamicAdd(a,j1,k1);}}
			}
		}

		else 	{		//All other rows
			if (a!=0.0){this->DynamicAdd(Sik[i]-a*Sij[i],i,k1);} 
			this->DynamicAdd(0.0,i,j1);
		}
	}
	delete [] Skrow;
	delete [] Sjrow;
	delete [] Sik;
	delete [] Sij;
}
/************************************************************************
 MakeEssential1:
	revises matrix to meet essential conditions X_{j1}+a*X_{k1}=b
	affects matrix only- if condition is removed, matrix must be remade!!!!
	current implementation is pretty slow- CAN OPTIMIZE!!
	revised from Akin, pg. 403
------------------------------------------------------------------------*/
void CSparseMatrix::MakeEssential2(const CSparseMatrix *A, const int j1, const int k1, const double a, const double b){
	int i,j;
	if ((j1<0)      || (k1<0))      {MatExitGracefully("CSparseMatrix::MakeEssential2: bad indices specified",RUNTIME_ERR);}
	if ((j1==k1) && (a!=0.0)) {return;} //cannot force constraint on itself
  

	//cout <<"CSparseMatrix::MakeEssential2: "<<a<<endl;
	for (i=0; i<size; i++){ //for each row in matrix
		
		if ((i==k1) && (a!=0.0)){ //special row (constraint for X_k1)
			for (j=0; j<size; j++){ //for each column

				if      (j==j1){this->DynamicAdd(a,k1,j1);}
				else if (j==k1){this->DynamicAdd(A->GetAij(k1,k1)-2.0*a*(A->GetAij(j1,k1)-A->GetAij(j1,j1))+a*a,k1,k1);}//
				else           {this->DynamicAdd(A->GetAij(k1,j)-a*A->GetAij(j1,j),k1,j);}
			}
		}

		else if (i==j1){ //special row (constraint for X_j1)
			for (j=0; j<size; j++){ //for each column
				this->SetAij(0.0,i,j); //does nothing for most sparse matrix rows
				if (j==j1){this->SetAij(1.0,j1,j1);}
				if (a!=0.0){if (j==k1){this->DynamicAdd(a,j1,k1);}}
			}
		}

		else 	{		//All other rows
			if (a!=0.0){this->DynamicAdd(A->GetAij(i,k1)-a*A->GetAij(i,j1),i,k1);} 
			this->DynamicAdd(0.0,i,j1);
		} 
	}

}
/************************************************************************
 IJ_TO_K:
	Translates i (row) and j (column) index of equivalent matrix A[size][size] 
	to single index [k] of sparse matrix sA[nNonZero]
	The k returned skips over k=size
------------------------------------------------------------------------*/
int CSparseMatrix::IJ_to_K   (const int i, const int j) const{
	if (i==j){return i;}
	else{
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			if (ija[k]==j){return k-1;}
		}
		return -1; //must always check
	}
}
/************************************************************************
 K_to_IJ:
	Translates k index of sparse matrix sA[nNonZero] to i,j of equivalent matrix A[size][size] 
------------------------------------------------------------------------*/
bool   CSparseMatrix::K_to_IJ    (const int k, int &i, int &j) const{
	i=-2;
	int ktmp;
	if ((k<0) || (k>nNonZero)){i=j=0; return false;}
	if (k<size){i=j=k;}
	else{
		ktmp=k+1;
    j=ija[ktmp]; 
		for (int itmp=0; itmp<=size; itmp++){ //searches for correct row i 
			if (ija[itmp]> ktmp){i=itmp-1;return true;} 
      if (ija[itmp]==ktmp){i=itmp;  return true;}
			i=itmp;//takes care of k==nNonZero
		}	
	}
  
	return true;
}

int CSparseMatrix::GetNumOffDiag (const int row) const{return ija[row+1]-ija[row];}
//-----------------------------------------------------------------------
int CSparseMatrix::GetNumRows    () const{return size;}
//-----------------------------------------------------------------------
double CSparseMatrix::GetAvgDiagCoeff  () const{
	double avg(0.0);
  for (int i=0; i<size;i++){
		avg+=sA[i];
	}
	return avg/size;
}
//-----------------------------------------------------------------------
double CSparseMatrix::GetMaxCoeff      () const{
	double max=0.0;
	for (int k=0; k<nNonZero+1; k++){
		if (k!=size){upperswap(max,fabs(sA[k]));}
	}
	return max;
}
/************************************************************************
 GetAij:
	Gets entry of equivalent matrix A, A[i][j]
	if A[i][j] not represented, returns 0.0;
	i=row, j=column
------------------------------------------------------------------------*/
double  CSparseMatrix::GetAij(const int i, const int j) const{
	if (i==j){return sA[i];}
	else{
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			//if (k==size){ExitGracefully("CSparseMatrix::GetAij: badly initialized matrix",RUNTIME_ERR);}
			if (ija[k]==j){return sA[k];}
		}
		return 0.0;
	}
}
/************************************************************************
 SetAij:
	Overwrites entry of equivalent matrix A, A[i][j] with a
	if A[i][j] not represented, returns false;
------------------------------------------------------------------------*/
bool CSparseMatrix::SetAij(const double &a, const int i, const int j){
	if (i==j){sA[i]=a;return true;}
	else{
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			if (ija[k]==j){sA[k]=a;return true;}
		}
	  if (a==0.0){return true;}
		else       {return false;} 
	}
}
/************************************************************************
 AddToAij:
	Updates entry of equivalent matrix A, A[i][j] by adding a
	if A[i][j] not represented, returns false;
------------------------------------------------------------------------*/
bool CSparseMatrix::AddToAij(const double &a, const int i, const int j){
	if (i==j){sA[i]+=a;return true;}
	else{
		for (int k=ija[i]; k<ija[i+1]; k++){ //searches through row i
			if (ija[k]==j){sA[k]+=a;return true;}
		}
	  if (a==0.0){return true;}
		else       {return false;} 
	}
}
/************************************************************************
 MatVectMult:
	Sparse Matrix vector multiply
------------------------------------------------------------------------*/
void CSparseMatrix::MatVectMult(Ironclad1DArray  x,
															  Writeable1DArray b,
																const int        N, 
																const bool       transpose,
																const bool       ignore) const {
	if (size!=N){
		MatExitGracefully("CSparseMatrix::MatVectMult: Cannot multiply matrix and vector of different sizes",BAD_DATA);}
	
	int i,k;	
	if (!transpose){
		for (i=0;i<size;i++){
			//if (!(ignore && dirichlet[i])){
				b[i]=sA[i]*x[i];
				for (k=ija[i];k<ija[i+1];k++){
					//if (!(ignore && dirichlet[ija[k]])){
						b[i]+=sA[k]*x[ija[k]];
					//}
				}
			//}
			//else{
			//	b[i]=x[i];
			//}
		}
	}
	else{
		int j;
		for (i=0;i<size;i++){
			if (!dirichlet[i]){
				b[i]=sA[i]*x[i];
			}
			else{
				b[i]=x[i];
			}
		}
		for (i=0;i<size;i++){
			//if (!dirichlet[i]){ //not sure if activity handled correctly here
				for (k=ija[i];k<ija[i+1];k++){
					j=ija[k];
					//if (!dirichlet[j]){
						b[j]+=sA[k]*x[i];
					//}
				}
			//}
		}
	}

}
/************************************************************************
	ScalarMult:
		multiplies sparse matrix by a scalar value w
------------------------------------------------------------------------*/
void CSparseMatrix::ScalarMult (const double w){
	for (int k=0; k<nNonZero+1; k++){
		if (k!=size){
			sA[k]*=w;
		}
	}
}
/************************************************************************
	Print
------------------------------------------------------------------------*/
void CSparseMatrix::Print(const bool full, int sample) const{
	double tmpA;
	if (sample>size){sample=size;} 
	for (int i=0;i<sample;i++){//for each row
		//if (ProgramAborted()){return;} 
		if (!dirichlet[i]){
			for (int j=0;j<sample;j++){//for each column
				tmpA=GetAij(i,j);
				if ((full) || (tmpA!=0.0)){cout <<tmpA<<"  ";}
			}
		}
		else{
			for (int j=0;j<sample;j++){
				if (i==j){tmpA=1.0;}else{tmpA=0.0;}
				if ((full) || (tmpA!=0.0)){cout <<tmpA<<"  ";}
			}
		}
		cout <<endl;
	}
}
/************************************************************************
	SolveEst:
		result,{x}, is {b} normalized by the diagonal of A
		S is sparse matrix
		same as asolve from Press et al
		can be changed to speed up convergence (using off diagonal terms?)
------------------------------------------------------------------------*/
void CSparseMatrix::SolveEst(Ironclad1DArray b, Writeable1DArray x, const int size, const bool transpose) const{

		for (int i=0; i<size;i++){
			if (!dirichlet[i]){x[i]=(sA[i]!=0.0 ? b[i]/sA[i] : b[i]);}
			else              {x[i]=b[i];                            }
		}

}

/************************************************************************
 CalculateAdjustedNorm:
	Calculates norm of a vector (dirichlet terms removed)
------------------------------------------------------------------------*/
double CSparseMatrix::CalculateAdjustedNorm(Ironclad1DArray b, int size, const normtype type) const{
	int i,isamax;
	static double norm;
	if      (type==VECTOR_MAGNITUDE_NORM){

		norm=0.0;for (i=0;i<size;i++){norm+=b[i]*b[i];}//if (!dirichlet[i]){norm+=b[i]*b[i];}}                           
		return sqrt(norm);
	}
	else if (type==LARGEST_COMPONENT_NORM){

		isamax=0;for (i=1;i<size;i++){if (fabs(b[i])>fabs(b[isamax])){isamax=i;}}//if (!dirichlet[i]){if (fabs(b[i])>fabs(b[isamax])){isamax=i;}}} 
		return fabs(b[isamax]);
	}
	else{                                                                  
		return 0.0;
	}
}
/************************************************************************
 BCG:
	Sparse Preconditioned Biconjugate Gradient Solver 
	S is sparse matrix representing A[size][size] matrix
	x[size] is output vector
	b[size] is RHS 
	BCGtype is testing criteria 
		if RELATIVE_RESIDUAL,       |Ax-b|/|b|         <BCG_tol
		if RELATIVE_TRANS_RESIDUAL, |At(Ax-b)|/|(At*b)|<BCG_tol
		if NORM_ERROR,              |xerr|/|x|         <BCG_tol
		if MAX_ERROR,               xerr_max/|x|_max   <BCG_tol

	MUST OPTIMIZE!!!!!
------------------------------------------------------------------------*/
void CSparseMatrix::BCG(Ironclad1DArray     b, 
												Writeable1DArray    x, 
												const int           size, 
												const BCGtestparam  BCGtype, 
												double             &err,
//												double              normalize,
												const double        BCG_tol) const{

	static double ak,akden;
	static double bk,bkden(1.0),bknum;
  static double bnorm,dxnorm,xnorm,zm1norm,znorm;
	int    j;	
	bkden=1.0;
	const double eps=1e-14;

	if (size!=this->size){MatExitGracefully("BCG: Matrix and vector not same size",BAD_DATA);}

	//double *bcopy=new double [size];
	/*for (j=0; j<size;j++){
		bcopy[j]=b[j];*/
		/*for (int i=0;i<size;i++){
			if (dirichlet[i]){
				bcopy[j]-=b[i]*GetAij(i,j); //do not want to do this!!
			}
		}*/
	/*}*/

	int iter=0;
	MatVectMult(x,r,size,false,true);
	
	for (j=0; j<size;j++){
		if (!dirichlet[j]){
			r [j]=b[j]-r[j];
			rr[j]=r[j];
		}
	}
	if (CalculateAdjustedNorm(r,size,VECTOR_MAGNITUDE_NORM)==0.0){
		cout <<"converged"<<endl;
		return;
		//ExitGracefully("BCG: residual norm is zero (already converged)",BAD_DATA);
	}

	//MatVectMult(r,rr,size,false,true); //minimum residual variant if uncommented (blows up dirichlet conditions)

	if      (BCGtype==RELATIVE_RESIDUAL){

		bnorm=CalculateAdjustedNorm(b,size,VECTOR_MAGNITUDE_NORM);
		if (bnorm==0.0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA);}

		SolveEst(r,z,size,false);
	}
	else if (BCGtype==RELATIVE_TRANS_RESIDUAL){

		SolveEst(b,z,size,false);
		bnorm=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA);}

		SolveEst(r,z,size,false);
	}
	else if (BCGtype==NORM_ERROR){

		SolveEst(b,z,size,false);
		bnorm=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA);}

		SolveEst(r,z,size,false);

		znorm=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM);
		if (znorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA);}
	}
	else if (BCGtype==MAX_ERROR){

		SolveEst(b,z,size,false);
		bnorm=CalculateAdjustedNorm(z,size,LARGEST_COMPONENT_NORM);
		if (bnorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA);}

		SolveEst(r,z,size,false);
		znorm=CalculateAdjustedNorm(z,size,LARGEST_COMPONENT_NORM);
		if (znorm==0){MatExitGracefully("BCG: bnorm is zero (all RHS elements zero)",BAD_DATA);}
	}
	else {MatExitGracefully("BCG: Illegal value for tolerance",RUNTIME_ERR);}

	//cout <<fixed <<setprecision(6);

	while (iter<=BCG_max_iter){
		iter++;
		SolveEst(rr,zz,size,true); //true indicates transpose matrix At
		
		bknum=0.0;
		for (j=0;j<size;j++){bknum+=z[j]*rr[j];}
		if (bknum==0.0){cout <<z[0]<<" "<<rr[0]<<endl;MatExitGracefully("BCG: beta numerator is zero",BAD_DATA);}
		
		//Calculate coeff bk and direction vectors p and pp--------------------------
		if (iter==1){
			for (j=0;j<size;j++){
				if (!dirichlet[j]){
					p [j]=z [j];
					pp[j]=zz[j];
				}
			}
		}
		else{
			bk=bknum/bkden;
			for (j=0;j<size;j++){
				if (!dirichlet[j]){
					p [j]=bk*p [j]+z [j];
					pp[j]=bk*pp[j]+zz[j];
				}
				else{
					p[j]=0.0;
					pp[j]=0.0;
				}
			}
		}

		//Calculate coeff ak, new iterate x and new residuals r and rr---------------
		bkden=bknum;
		
		MatVectMult(p,z,size,false,true);

		akden=0.0; 
		for (j=0;j<size;j++){if (!dirichlet[j]){akden+=z[j]*pp[j];}} //{S*p}*{S/rr}
		if (akden==0){MatExitGracefully("BCG: alpha denominator is zero",BAD_DATA);}

		ak=bknum/akden;

		MatVectMult(pp,zz,size,true,true);
		for (j=0;j<size;j++){
			if (!dirichlet[j]){
				x [j]+=ak*p [j];
				r	[j]-=ak*z [j];
				rr[j]-=ak*zz[j];
			}
			else{
				x[j]=b[j];
				r[j]=0.0;
				rr[j]=0.0;
			}
		}
		
		//Solve A*z=r and check stopping criteria------------------------------------
		SolveEst(r,z,size,false);

		if       (BCGtype==RELATIVE_RESIDUAL){
			
			err=CalculateAdjustedNorm(r,size,VECTOR_MAGNITUDE_NORM)/bnorm;
		}
		else if  (BCGtype==RELATIVE_TRANS_RESIDUAL){

			err=CalculateAdjustedNorm(z,size,VECTOR_MAGNITUDE_NORM)/bnorm;
		}
		else if ((BCGtype==NORM_ERROR) || 
			       (BCGtype==MAX_ERROR)){

			normtype ntype;
			if (BCGtype==NORM_ERROR){ntype=VECTOR_MAGNITUDE_NORM;}
			else                    {ntype=LARGEST_COMPONENT_NORM;}

			zm1norm=znorm;
			znorm  =CalculateAdjustedNorm(z,size,ntype);
			if (fabs(zm1norm-znorm)>eps*znorm){
				dxnorm=fabs(ak)*CalculateAdjustedNorm(p,size,ntype);
				err=znorm/fabs(zm1norm-znorm)*dxnorm;
			}
			else{
				err=znorm/bnorm;
				continue;
			}
			xnorm=CalculateAdjustedNorm(x,size,ntype);
			if (err<=0.5*xnorm){err/=xnorm;}
			else{
				err=znorm/bnorm;
				continue;
			}
		}
		//cout <<"BCG iter: " <<iter<<" "<<err<<endl;
		if (err<=BCG_tol){break;}

	}//end while
	//cout <<"BCG_iter:"<<iter<<endl;

	//cout <<"BCG Condition number"<<endl;
	if (iter>=BCG_max_iter){MatExitGracefully("BCG: Unable to solve system of equations",RUNTIME_ERR);}
}
