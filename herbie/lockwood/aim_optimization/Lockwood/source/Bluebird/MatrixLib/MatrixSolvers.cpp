//Matrix Solvers 

#include "MatrixInclude.h"

#include <math.h>
#include <iostream>

/************************************************************************
 MatMult:
	Multiplies square matrix and vector A*x. Returns b
-----------------------------------------------------------------------*/
void   MatVectMult    (Ironclad2DArray A, Ironclad1DArray x,Writeable1DArray b, const int size, bool transpose){
	if (!transpose){
		for(int i=0;i<size;i++){
			b[i]=0.0;  for(int j=0; j<size;j++){b[i]+=A[i][j]*x[j];}}
	}
	else{
		for(int i=0;i<size;i++){
			b[i]=0.0;  for(int j=0; j<size;j++){b[i]+=A[j][i]*x[j];}}	
	}
}
/************************************************************************
 MatMult:
	Multiplies matrix and vector A*x. Returns b
-----------------------------------------------------------------------*/
void   MatVectMult    (Ironclad2DArray A, Ironclad1DArray x,Writeable1DArray b, const int size1, const int size2){

	for(int i=0;i<size1;i++){
		b[i]=0.0;  for(int j=0; j<size2;j++){b[i]+=A[i][j]*x[j];}}

}
/************************************************************************
 MatMult:
	Multiplies NxM matrix A times MxP matrix B. Returns C (NxP)
-----------------------------------------------------------------------*/
void   MatMult(Ironclad2DArray A, Ironclad2DArray B,Writeable2DArray C, const int N, const int M, const int P){
	for (int n=0; n<N; n++){
		for (int p=0; p<P; p++){
			C[n][p]=0.0;
			for (int m=0; m<M; m++){
				C[n][p]+=A[n][m]*B[m][p];
			}
		}
	}
}
/************************************************************************
 Gauss:
	Performs Gauss Elimination on the matrix a and RHS b. Returns solution, x and rank
	this operation destroys the A matrix and RHS
	OPTIMIZED
-----------------------------------------------------------------------*/
bool Gauss(Writeable2DArray A, 
					 Writeable1DArray b,
					 Writeable1DArray x, 
					 const int        size){
  static int nrow,j,i; 
	static double Y,t,v,amax;

	/*cout <<endl;
	for (j=0;j<size;j++){  
		for (i=0;i<size;i++){  
			cout << A[i][j] <<" ";
		}
		cout << "| "<<b[j]<<" | "<<endl;
	}
	cout <<endl;*/


  for (j=0;j<size;j++){  
    //sort rows by size of A[j][j]
    amax=sqrt(pow(A[j][j],2)); nrow=j;
    for (int jrow=j+1; jrow<size;jrow++){ 
      if(sqrt(pow(A[jrow][j],2))> amax){
				amax=sqrt(pow(A[jrow][j],2));
				nrow=jrow;
			}
		}
    for(int irow=0; irow < size; irow++){
      t=A[j][irow]; 
			A[j][irow]=A[nrow][irow]; 
			A[nrow][irow]=t;
		}    
		t=b[j]; 
		b[j]=b[nrow]; 
		b[nrow]=t;
    //perform elimination on lower rows 
    for (i=(j+1);i<size;i++){ 
      v=A[i][j]/A[j][j];
      for (int k=j; k<size;k++){
				A[i][k]=A[i][k]-A[j][k]* v;
			}
      b[i]=b[i]-b[j]*v;
		}
	}

  for(i=0; i < size; i++){
    x[i]=0.0; 
    if(sqrt(pow(A[i][i],2))< gauss_min){  return false; }
	}
  //backsubstitution
  x[size-1]=b[size-1]/A[size-1][size-1];
  for (i=size-1; i>=0; i--){
    Y=0; 
		for (int l=size-1; l>=i; l--){
			Y+=x[l+1]*A[i][l+1];
		}
    x[i]=(b[i]-Y)/A[i][i];
	}
	/*for (j=0;j<size;j++){  
		cout << "| "<<x[j]<<" | "<<endl;
	}*/
	return true;
}
/************************************************************************
 PCG:
	Preconditioned Conjugate Gradient Solver (still some bugs!)
------------------------------------------------------------------------*/
void  PCG(double **A,double *b,double *x, const int size, const double tolerance, int &numiter, const cgm_type type){
	
	double *res;						//residual vector
	double *p;							//improving direction
	double *w,*z;           //temporary vectors
	double alpha,beta;      //step length
	int    iter;
	int    i;
  bool   done(false);
	double sumr(0),maxr(0);	
	double sumb(0),maxb(0);
	double sum1,sum2;     
  
	//BICONJUGATE parameters
	double sumrc,maxrc;
	double *pc;
	double *rc;
	double *wc;
	double  rho;
	pc=new double [size];
	rc=new double [size];
	wc=new double [size];


	res=new double [size];
	p  =new double [size];
	w  =new double [size];
	z  =new double [size];

	//****initialize***********************************************

	//prepare termination criteria (sumb, maxb)
	for (i=0;i<size; i++){
		sumb+=b[i]*b[i];
		if (fabs(b[i])>maxb){maxb=fabs(b[i]);}
	}

	//initial guess
	for (i=0; i<size; i++){
		x[i]=0.0;//b[i]/A[i][i];
	}

	//initial residual {r}={b}-[A]{x}, direction vector {p}={r}
	MatVectMult(A,x,w,size,false);
	for (i=0;i<size; i++){
		res[i]=b[i]-w[i];
		p  [i]=res[i];
	}
	if (type==BICONJUGATE){
		rho=0.0;
		for (i=0;i<size; i++){
			pc[i]=res[i];
			rc[i]=res[i];
			rho+=rc[i]*res[i];
		}
	}

	//****iterate**************************************************
	iter=0;
	do {
		cout << "PCJ iter"<<iter << endl;
		
		iter++;
		//get {w}----------------------------------------------------  
		MatVectMult(A,p,w,size,false);	  //{w}=[A]{p}
		if (type==BICONJUGATE){
		 MatVectMult(A,pc,wc,size,true);	//{wc}=[A]T{pc}
		}

		//get alpha--------------------------------------------------
		sum1=0.0; sum2=0,0;
	  if (type==BICONJUGATE){
			for (i=0; i<size; i++){    //alpha=({rc}*{r})/({pc}*{w})
				sum1+=rc[i]*res[i];
				sum2+=pc[i]*w[i];
			}
		}
		else{                        //alpha=({p}*{r})/({p}*{w})
			for (i=0; i<size; i++){
				sum1+=p[i]*res[i];
				sum2+=p[i]*w[i];
			}
		}
		alpha=sum1/sum2;
		cout <<"alpha:"<<alpha<<" sum1: "<< sum1 <<" sum2: "<<sum2<<endl;

		//improve solution guess ------------------------------------
		for (i=0; i<size; i++){
			x[i]+=alpha*p[i];          //{x}={x}+alpha*{p}
		}
		
		//calculate/estimate new residual ---------------------------			
		if ((iter%8)==0){
		  MatVectMult(A,x,w,size,false);
			for (i=0; i<size; i++){
				res[i]=b[i]-w[i];        //actual residual {r}={b}-[A]{x}
			}		
			if (type==BICONJUGATE){    //not sure about this
        MatVectMult(A,x,wc,size,true);
        for (i=0; i<size; i++){
					rc[i]=b[i]-wc[i];
				}
				/*for (i=0; i<size; i++){
					rc[i]-=alpha*wc[i];
				}*/
			}
		}
		else{
			for (i=0; i<size; i++){
				res[i]-=alpha*w[i];      //guess at residual {r}={r}-a*([A]{x})old 
			}	
			if (type==BICONJUGATE){
				for (i=0; i<size; i++){
					rc[i]-=alpha*wc[i];
				}
			}
		}

		//calculate termination criteria-------------------------------
		sumr=maxr=0;
		for (i=0; i<size; i++){ 
			sumr+=res[i]*res[i];                         // |r|=sum(r^2)
			if (fabs(res[i])>maxr){maxr=fabs(res[i]);}   // rmax=max(r)
		}	
		if (type==BICONJUGATE){
			sumrc=maxrc=0;
			for (i=0; i<size; i++){
				sumrc+=rc[i]*rc[i];
				if (fabs(rc[i])>maxrc){maxrc=fabs(rc[i]);}
			}
		}
		cout << "residual: " << sumr  << " " << maxr  <<endl;
		cout << "cresidual:" << sumrc << " " << maxrc <<endl;

		if      ((type!=BICONJUGATE) &&
			       (sumr <tolerance*sumb) && (maxr <tolerance*maxb))   {done=true;}
		else if ((type==BICONJUGATE) && 
			       (sumrc<tolerance*sumb) && (maxrc<tolerance*maxb) && 
						 (sumr <tolerance*sumb) && (maxr <tolerance*maxb))   {done=true;}
		else if (iter>numiter)                                       {done=true;}  //bad done
		else{
			//choose new improving direction------------------------------
			sum1=sum2=0.0;
			switch (type){

			case(BASIC_CGM):{                //-------------------------------
				sum1=sum2=0,0;
				for (i=0; i<size; i++){
					sum1+= res[i]*w[i];
					sum2+= p  [i]*w[i];
				}
				beta=-sum1/sum2;          //b=-({r}*{w})/({p}*{w})
				for (i=0; i<size; i++){
					p[i]=res[i]+beta*p[i];  // {p}={r}+b*{p}
				}						
			}
			case(PRECONDITIONED_CGM):{       //-------------------------------
					//fill in later
					/*btemp[0]=A[0][0]; 
					ctemp[0]=A[0][1];
					for (i=1; i<size-1; i++){
            atemp[i]=A[i][i-1];          
						btemp[i]=A[i][i];          
						ctemp[i]=A[i][i+1];           
					}
					atemp[size-1]=A[size  ][size-1]; 
					btemp[size-1]=A[size-1][size-1]; 
					TriDiagonal(atemp, btemp, ctemp, r, z, N, code);
					*/
					for (i=0; i<size; i++){
						z[i]=res[i]/A[i][i];
					}
					sum1=sum2=0,0;
					for (i=0; i<size; i++){
	          sum1+=z[i]*w[i];
	          sum2+=p[i]*w[i];
					}
	        beta=-sum1/sum2;
					for (i=0; i<size; i++){
						p[i]=z[i]+beta*p[i];
					}
					
				}
				case(BICONJUGATE):         //-------------------------------
				{
					double rhoold(rho);
					rho=0;
					for (i=0; i<size; i++){
						rho+=rc[i]*res[i];
					}			
					beta=rho/rhoold;
					for (i=0; i<size; i++){
						p[i]=res[i]+beta*p[i];
						pc[i]=rc[i]+beta*pc[i];
					}
				}
				default:{
				}
			}
		}
	} while (!done);

	numiter=iter;

	delete [] res;
	delete [] p;
	delete [] w;
	delete [] z;

	delete [] pc;
	delete [] rc;
	delete [] wc;
}
/************************************************************************
 CalculateNorm:
	Calculates norm of a vector
------------------------------------------------------------------------*/
double CalculateNorm(Ironclad1DArray b, int size, const normtype type){
	int i,isamax;
	static double norm;
	if      (type==VECTOR_MAGNITUDE_NORM){

		norm=0.0;for (i=0;i<size;i++){norm+=b[i]*b[i];}                           
		return sqrt(norm);
	}
	else if (type==LARGEST_COMPONENT_NORM){

		isamax=0;for (i=0;i<size;i++){if (fabs(b[i])>fabs(b[isamax])){isamax=i;}} 
		return fabs(b[isamax]);
	}
	else{                                                                  
		return 0.0;
	}
}


/************************************************************************
 ThomasAlgorithm:
	Solves Tridiagonal matrix with diagonals e,f, and g. RHS=b
	Returns solution vector x
------------------------------------------------------------------------*/
void ThomasAlgorithm(Ironclad1DArray  e,  Ironclad1DArray f,
										 Ironclad1DArray  g,  Ironclad1DArray b, 
										 Writeable1DArray x,const int     size){
	int i;
	double *f1=new double [size];
	double *b1=new double [size];

	//calculate f prime vector
	f1[0]=f[0];
	for (i=1; i<size; i++){
		f1[i]=f[i]-(e[i]*g[i-1]/f1[i-1]);
	}
	//calculate b prime vector
	b1[0]=b[0];
	for (i=1; i<size; i++){
		b1[i]=b[i]-(e[i]*b1[i-1]/f1[i-1]);
	}
	//solve for x using back substitution
	x[size-1]=b1[size-1]/f1[size-1];
	for (i=size-2; i>=0; i--){
		x[i]=(b1[i]-g[i]*x[i+1])/f1[i];
	}
	delete [] f1;
	delete [] b1;
}
inline double SIGN(const double &a, const double &b){
	return b>=0.0 ? (a>=0.0 ? a : -a): (a>=0.0 ? -a : a);
}
inline bool BADIND(const int i, const int k, const int n){
	if ((i<0) || (i>=n)){return true;}
	if ((k<0) || (k>=n)){return true;}
	return false;
}
/************************************************************************
 Singular Value Decomposition:
Returns solution, x and rank
	this operation destroys the A matrix
-----------------------------------------------------------------------*/
bool SVD(Writeable2DArray A, 
				 Writeable1DArray b,
				 Writeable1DArray x, 
				 const int size){
	bool flag;
	int i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,xx,y,z;

	int m=size;
	int n=size;
	double *tmp=new double [n];
	double *w  =new double [n];
	double **v =new double *[n];
	for (i=0;i<n;i++){
		v[i]=new double [n]; 
	}
	if (v==NULL){MatExitGracefully("SVD: out of memory",OUT_OF_MEMORY);}

	/*	cout <<endl;
		for (j=0;j<size;j++){  
			for (i=0;i<size;i++){  
				cout << A[i][j] <<" ";
			}
			cout << "| "<<b[j]<<" | "<<endl;
		}
		cout <<endl;*/

	g=scale=anorm=0.0;
	for (i=0; i<n; i++){
		l=i+2;
		tmp[i]=scale*g;
		g=s=scale=0.0;
		if (i<m){
			for (k=i;k<m;k++){scale+=fabs(A[k][i]);}
			if (scale!=0.0){
				for (k=i; k<m; k++){
					A[k][i]/=scale;
					s+=A[k][i]*A[k][i];
				}
				f=A[i][i];
				g=-SIGN(sqrt(s),f);
				h=f*g-s;
				A[i][i]=f-g;
				for (j=l-1;j<n;j++){
					s=0.0;
					for (k=i;k<m;k++){s+=A[k][i]*A[k][j];}
					f=s/h;
					for (k=i;k<m;k++){A[k][j]+=f*A[k][i];}
				}
				for (k=i;k<m;k++){A[k][i]*=scale;}
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if ((i+1<=m) && (i!=n)){
			for (k=l-1;k<n;k++){scale+=fabs(A[i][k]);}
			if (scale!=0.0){
				for (k=l-1;k<n;k++){
					A[i][k]/=scale;
					s+=A[i][k]*A[i][k];
				}
				f=A[i][l-1];
				g=-SIGN(sqrt(s),f);
				h=f*g-s;
				A[i][l-1]=f-g;
				for (k=l-1;k<n;k++){tmp[k]=A[i][k]/h;}
				for (j=l-1;j<m;j++){
					s=0.0;
					for (k=l-1;k<n;k++){s      +=A[j][k]*A[i][k];}
					for (k=l-1;k<n;k++){A[j][k]+=s*tmp[k];  }			
				}
				for (k=l-1;k<n;k++){A[i][k]*=scale;}
			}
		}
		upperswap(anorm,(fabs(w[i])+fabs(tmp[i])));
	}
	for (i=n-1; i>=0; i--){
		if (i<n-1){
			if (g!=0.0){
				for (j=l;j<n;j++){v[j][i]=(A[i][j]/A[i][l])/g;}
				for (j=l;j<n;j++){
					s=0.0;
					for (k=l;k<n;k++){s 		 +=A[i][k]*v[k][j];}
					for (k=l;k<n;k++){v[k][j]+=s*v[k][i];}
				}
			}
			for (j=l;j<n;j++){v[i][j]=v[j][i]=0.0; }
		}
		v[i][i]=1.0;
		g=tmp[i];
		l=i;
	}
	for (i=min(n,m)-1;i>=0;i--){
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++){A[i][j]=0.0;}
		if(g!=0.0){
			g=1.0/g;
			for(j=l;j<n;j++){
				s=0.0;
				for (k=l;k<m;k++){s+=A[k][i]*A[k][j];}
				f=(s/A[i][i])*g;
				for (k=i;k<m;k++){A[k][j]+=f*A[k][i];}
			}
			for(j=i;j<m;j++){A[j][i]*=g;}
		}
		else{
			for(j=i;j<m;j++){A[j][i]=0.0;}
		}
		++A[i][i];
		//A[i][i]++;
	}
	for (k=n-1;k>=0;k--){
		for(its=0; its<30;its++){
			//cout <<"SVD iter "<< its <<endl;
			flag=true;
			for (l=k;l>=0;l--){
				nm=l-1;
				if (fabs(tmp[l])+anorm==anorm){
					flag=false;
					break;
				}
				if (fabs(w[nm])+anorm==anorm){break;}
			}//end for l
			if (flag){
				c=0.0;
				s=1.0;
				for (i=l-1;i<k+1;i++){
					f=s*tmp[i];
					tmp[i]*=c;
					if (fabs(f)+anorm==anorm){break;}
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=-f*h;
					for (j=0;j<m;j++){
						y=A[j][nm];
						z=A[j][i];
						A[j][nm]=y*c+z*s;
						A[j][i]=z*c-y*s;
					}
				}
			}//end if flag
			z=w[k];
			if (l==k){
				if (z<0.0){
					w[k]=-z;
					for (j=0;j<n;j++){v[j][k]=-v[j][k];}
				}
				break;
			}
			if (its==29){MatExitGracefully("SVD unconverged",SINGMAT);}
			xx=w[l];
			nm=k-1;
			y=w[nm];
			g=tmp[nm];
			h=tmp[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((xx-z)*(xx+z)+h*((y/(f+SIGN(g,f)))-h))/xx;
			c=s=1.0;
			for (j=l; j<=nm; j++){
				i=j+1;
				g=tmp[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				tmp[j]=z;
				c=f/z;
				s=h/z;
				f=xx*c+g*s;
				g=g*c-xx*s;
				h=y*s;
				y*=c;
				for (jj=0;jj<n;jj++){
					xx=v[jj][j];
					z=v[jj][i];
					v[jj][j]=xx*c+z*s;
					v[jj][i]=z*c-xx*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if(z!=0.0){
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				xx=c*y-s*g;
				for (jj=0;jj<m;jj++){
					y=A[jj][j];
					z=A[jj][i];
					A[jj][j]=y*c+z*s;
					A[jj][i]=z*c-y*s;
				}
			}//end for j=l to nm
			tmp[l]=0.0;
			tmp[k]=f;
			w[k]=xx;
		}//end for its
	}
	//end sub svdcmp (Press et al)

	double maxw(-ALMOST_INF);
	double minw( ALMOST_INF);

	for (j=0; j<n; j++){upperswap(maxw,w[j]);}              //get condition number info
	for (j=0; j<n; j++){ 
		if (fabs(maxw/w[j])>1e8){w[j]=0.0;}                    //zero out bad terms
	}
	for (j=0; j<n; j++){if (w[j]!=0) lowerswap(minw,w[j]);} //reevaluate min

	//back substitute
	for (j=0; j<n; j++){
		s=0.0;
		if (w[j]!=0.0){
			for (i=0; i<m; i++){s+=A[i][j]*b[i];}
			s/=w[j];
		}
		tmp[j]=s;
	}

	for (j=0; j<n; j++){ 
		s=0.0;
		for (jj=0;jj<n;jj++){s+=v[j][jj]*tmp[jj];}
		x[j]=s;
	}

	/*	for (j=0;j<size;j++){  
			cout << "| "<<x[j]<<" | "<< w[j]<<endl;
		}*/

	//calculate error

	//cout << "n: "<<n<< " m: "<<m<<endl;
	delete [] tmp;
	delete [] w;
	for (i=0; i<n; i++){delete [] v[i];} delete [] v;

	//cout <<"SVD CONDITION #: " <<maxw/minw<<endl;
	return true; //should be dependent upon condition number (max wj/minwj) should be low

}
void   MatExitGracefully(char *statement, badcode code)
{
//TMP DEBUG
}
void   MatExitGracefullyIf(bool condition, char *statement, badcode code)
{
//TMP DEBUG
}