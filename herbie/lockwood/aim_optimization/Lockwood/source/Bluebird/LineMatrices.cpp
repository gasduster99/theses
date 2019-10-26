#include "StringElem.h"

/***********************************************************************
				STATIC MEMBER INITIALIZATION
***********************************************************************/
double**				CLineElem::CD=  new double *[MAX_FAR_FIELD];
double**        CLineElem::f=   new double *[MAX_ORDER+1];
double**        CLineElem::g=   new double *[MAX_ORDER+1];
double**        CLineElem::A=   new double *[MAX_ORDER+5];
double**        CLineElem::unit=new double *[MAXLINECONTROL]; 
cmplex*         CLineElem::cheby=new cmplex [MAX_ORDER+1];
double*         CLineElem::rhs= new double  [MAXLINECONTROL];
double*         CLineElem::b=   new double  [MAX_ORDER+5];
double*         CLineElem::sol= new double  [MAX_ORDER+5];
//double*         CLineElem::c[4];
double          CLineElem::c1=0;
double          CLineElem::c2=0;
double          CLineElem::c3=0;
double          CLineElem::c4=0;
double          CLineElem::ConstrntVal[3];
bool            CLineElem::ConstrntOn [3];

/***********************************************************************
				Prepare matrices for string solution
***********************************************************************/
void CLineElem::Prepare(){
	int i,j,k,l,m,n;
  //create storage for permanent matrices
  for(i=0;i< MAX_FAR_FIELD;  i++){  CD[i]=new double[MAX_ORDER+1];   }
  for(i=0;i< MAX_ORDER+1;    i++){   f[i]=new double[MAX_ORDER+1];   }
	for(i=0;i< MAX_ORDER+1;    i++){   g[i]=new double[MAX_ORDER+1];   }
  for(i=0;i< MAX_ORDER+5;    i++){   A[i]=new double[MAX_ORDER+5];   }
	for(i=0;i< MAXLINECONTROL; i++){unit[i]=new double[MAX_ORDER+1];   }

  double C[MAX_FAR_FIELD][MAX_ORDER+1];
	double D[MAX_ORDER+1][MAX_ORDER+1];
  double sum=0;

  //create D matrix 
  //---------------------------------------------------
  //initialize
  for (i=0; i<=MAX_ORDER; i++){
    for (k=0; k<=MAX_ORDER; k++){D[i][k]=0;}}
  //fill first row
  for (i=0; i<=MAX_ORDER; i+=2){
    if ((i%4)==0){D[0][i]=1;} else {D[0][i]=-1;}}
  D[1][1]=1;
  for (i=1; i<=MAX_ORDER; i++){
    for (k=i; k<=MAX_ORDER; k+=2){
      if (k!=1){D[i][k]=2*D[i-1][k-1]-D[i][k-2];}}}

  //create C matrix 
  //---------------------------------------------------
  for (i=0; i<MAX_FAR_FIELD; i++){
    for (k=0; k<=MAX_ORDER; k++){
      if (((i+k)%2)==0){C[i][k]=(2.0/(double(i+k+1)));}else {C[i][k]=0;}}}

  //create CD matrix 
  //---------------------------------------------------
  for (i=0; i<MAX_FAR_FIELD; i++){
    for (j=0; j<=MAX_ORDER; j++){		
      sum=0.0;
      for (k=0; k<=MAX_ORDER; k++){sum+=C[i][k]*D[k][j];}
	    CD[i][j]=sum*IOVER2PI;}}
	for (i=MAX_FAR_FIELD-1; i>0; i--){
    for (j=0; j<=MAX_ORDER; j++){	CD[i][j]=CD[i-1][j];}}

  //create f  matrix 
  //---------------------------------------------------
  for (m=0; m<=MAX_ORDER; m++){
    for (n=0; n<=MAX_ORDER; n++){
      if ((m<n) && (((m+n)%2)==1) && (m>0))    {f[n][m]= 4.0/double(n-m);}
      else if ((m==0) && ((n%2)==1) && (n!=0)) {f[n][m]= 2.0/double(n);}
      else                                     {f[n][m]= 0.0;}}}

  //create g  matrix
  //---------------------------------------------------
  for (n=0; n<=MAX_ORDER; n++){
    for (m=0; m<=MAX_ORDER; m++){
      sum=0.0;
      if ((m<n) && ((m+n)%2==0) && (m>0)){
				for (l=1;l<=((n-m)/2); l++){sum+=double(n-2*l+1)/double(2*l-1);}
				g[n][m]=8.0*sum;}
			else if((m==0) && ((n%2)==0) && (n!=0)){
				for (l=1;l<=(n/2)    ; l++){sum+=double(n-2*l+1)/double(2*l-1);}
				g[n][m]=4.0*sum;}
      else g[n][m]=0.0;}}
	
	//initialize A,b, unit matrices
	//---------------------------------------------------
  for (n=0; n<MAX_ORDER+5; n++){
		for (m=0; m<MAX_ORDER+5; m++){
			A[n][m]=0.0;
		}
		b[n]=0.0;
	}
	for (n=0; n<=MAX_ORDER; n++){
    for (m=0; m<MAXLINECONTROL; m++){
			unit[m][n]=0.0;
		}
	}
}
//-----------------------------------------------------------------------
void CLineElem::Destroy(){
	if (globaldebug){cout <<"DESTROYING LINE ELEMENT STATIC DATA"<<endl;}
	int i;
	delete [] rhs;
	delete [] b;
	delete [] sol;
  for(i=0;i< MAX_FAR_FIELD;  i++){delete []   CD[i];}
	for(i=0;i< MAX_ORDER+1;    i++){delete []    f[i];}  
	for(i=0;i< MAX_ORDER+1;    i++){delete []    g[i];}
  for(i=0;i< MAX_ORDER+5;    i++){delete []    A[i];}
	for(i=0;i< MAXLINECONTROL; i++){delete [] unit[i];} 
  delete [] CD;
	delete [] f;
	delete [] g;
	delete [] A;
	delete [] unit;
}




























/***********************************************************************
				MATRIX ACCESSOR FUNCTIONS
***********************************************************************
***********************************************************************/
/*cmplex CLineElem::h(const cmplex *pCoeff, const int n) const{
	if      (n>=2){ return (pCoeff[n-1]-pCoeff[n+1])/(double)(n+n);}
  else if (n==1){ return pCoeff[0]-0.5*pCoeff[2];}
	else          { return 0.0;}
}
//-----------------------------------------------------------------------
cmplex CLineElem::hh(const cmplex *pCoeff, const int n) const{
	if      (n>=3){ return (pCoeff[n-2]+2.0*pCoeff[n]-pCoeff[n+2])/(double)(2.0*n*n);}
  else if (n==2){ return 0.25*pCoeff[0]-pCoeff[2]/12.0-pCoeff[4]/24.0;}
	else          { return 0.0;}
}*/





/*double          CLineElem::Xsaved       =1000;
double          CLineElem::Xlogtermsaved=1;
cmplex          CLineElem::Zsaved       =cmplex(1000,0);
cmplex          CLineElem::Zlogtermsaved=cmplex(0,0);*/
//***********************************************************************
/*cmplex CLineElem::f(const int r,const int c,const cmplex Z){
  if (r!=c)																	{return pf[r][c];}
	if(Z!=Zsaved)															{Zsaved=Z; Zlogtermsaved=log(Zm1oZp1(Z));}
	//if(Z!=Zsaved)															{Zsaved=Z; Zlogtermsaved=log((Z-1.0)/(Z+1.0));}
  return Zlogtermsaved;
}
//***********************************************************************
double CLineElem::f(const int r,const int c,const double X){	
  if (r!=c)																	{return pf[r][c];}
  if (X!=Xsaved)														{Xsaved=X; Xlogtermsaved=log(fabs((X-1.0)/(X+1.0)));}
  return Xlogtermsaved;
}*/
//***********************************************************************
/*cmplex CLineElem::g(const int r,const int c,const cmplex Z){
  //if (Z!=Zsaved) {Zsaved=Z; Zlogtermsaved=log((Z-1.0)/(Z+1.0));}
	if(Z!=Zsaved)															{Zsaved=Z; Zlogtermsaved=log(Zm1oZp1(Z));}
  if ((c<r)&& (((r+c)%2)==1) && (c>0))      {return 2.0*double(r)*Zlogtermsaved;}
  else if ((c==0) && ((r%2)==1) && (r!=0))  {return     double(r)*Zlogtermsaved;}
  else if (r==c)                            {return 2.0/((Z-1.0)*(Z+1.0));}
  else                                      {return pg[r][c];}
}
//***********************************************************************
double CLineElem::g(const int r,const int c,const double X){
  if (X!=Xsaved)  {Xsaved=X; Xlogtermsaved=log(fabs((X-1.0)/(X+1.0)));}
  if      ((c<r) &&(((r+c)%2)==1) && (c>0))	{return 2.0*r*Xlogtermsaved;}
  else if ((c==0)&&((r%2)==1)     && (r!=0)){return r*Xlogtermsaved;}
  else if (r==c)														{return 2.0/((X-1.0)*(X+1.0));}
  else																			{return pg[r][c];}
}*/
//***********************************************************************
//double CLineElem::CD(const int r,const int c)					{return pCD[r][c];}
//***********************************************************************
//double CLineElem::D(const int r,const int c)					{return pD[r][c];}

/*double CLineElem::B(const int r,const int c) const{
  if ((r==0) && ((c%2)==1))									{return (2.0/(double)(c));}
  else if ((r<c) && (((r+c)%2)==1))					{return (4.0/(double)(c-r));}
  else																			{return 0;}
}	
//-----------------------------------------------------------------------
double CLineElem::d(const int n,const int j,double **a,const int Nj) const{
  if (n>=Nj)															  {return 0;}
  else if (n==0)														{return (d(2,j,a,Nj)/2.0)+a[j][1];}
  else																			{return d(n+2,j,a,Nj)+2.0*(n+1)*a[j][n+1];}
}*/
//-----------------------------------------------------------------------