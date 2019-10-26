
#include "MasterInclude.h"
#include "MatrixInclude.h"
/***********************************************************************
				Laurent
************************************************************************
Calculates Laurent Expansion about Z=0.0 using horner's rule (knuth's rule for the real version)
	should translate into assembly code
	MUST OPTIMIZE
-----------------------------------------------------------------------*/
cmplex Laurent(const cmplex &Z, const int N, const cmplex * const &cf ){
	static cmplex w,sum; 
	w=1.0/Z;
  sum=cf[N-1];
  for (int i=N-2; i>=0; i--){ sum=(sum*w)+cf[i]; }
  return sum; 
}
//-----------------------------------------------------------------------
cmplex LaurentRe(const cmplex &Z, const int N, const double * const &cf){
	static cmplex w;
	static double r,s,c,d,e;
	w=1.0/Z;
	r=2.0*w.real();
  s=(w*conj(w)).real();
	c=cf[N-1]; 
	d=cf[N-2];

  for (int i=N-3; i>=0; i--){
    e=c;
	  c=d+r*e;
	  d=cf[i]-s*e;
  }
  return -(c*w+d);
}
/***********************************************************************
				OutsideW 
***********************************************************************
calculates derivative of Laurent Expansion
MUST OPTIMIZE
-----------------------------------------------------------------------*/
cmplex LaurentDer(const cmplex &Z,const int N, const cmplex * const &cf, const double r){
	//Returns negative derivative (W=-dPhi/dz)
  static cmplex w,sum;
  if (abs(Z)>=1.0){ 
		w=1.0/Z;
    sum=cf[N-1]*(double)(N-1);
    for (int i=N-2; i>=0; i--){ sum=(sum*w)+(cf[i]*(double)(i)); }
    return sum*w/r;
  }
  else{            //Taylor Series
    return 0.0;
  }
}	
//-----------------------------------------------------------------------
cmplex LaurentDerRe(const cmplex &Z,const int N, const double * const &cf){
	//knuths rule- derivative
  static int i;
	static cmplex w;
	static double r,s,c,d,e;
	w=1.0/Z;
  r=2.0*w.real();
  s=(w*conj(w)).real();
  c=cf[N-1]*(N-1);
	d=cf[N-2]*(N-2);

  for (i=N-3; i>=0; i--){
    e=c;
	  c=d+r*e;
	  d=cf[i]*(double)(i)-s*e;
  }
  return -(c*w+d)*w;

}	
/***********************************************************************
				OutsideWder 
***********************************************************************
calculates 2nd derivative of Laurent Expansion in x direction
MUST OPTIMIZE
-----------------------------------------------------------------------*/
cmplex LaurentDxx(const cmplex &Z,const int N, const cmplex * const &cf, const double r){
  static cmplex w,sum;

  if (abs(Z)>=1.0){  //Laurent Series
		w=1.0/Z;
    sum=cf[N-1]*double((N-1)*(N-2));
    for (int i=N-2; i>=0; i--){ sum=(sum*w)+(cf[i]*(double)((i+1)*i)); }
    return sum*w*w/r/r;
  }
  else{           //Taylor series
		return 0.0; 
  }
}	
//-----------------------------------------------------------------------
cmplex LaurentDxxRe(const cmplex &Z,const int N, const double * const &cf){
	//knuths rule- 2nd derivative - correct!
  static int i;
	static cmplex w;
	static double r,s,c,d,e;
	w=1.0/Z;
  r=2.0*w.real();
  s=(w*conj(w)).real();
  c=cf[N-1]*(N-1)*(N-2);
	d=cf[N-2]*(N-2)*(N-3);

  for (i=N-3; i>=0; i--){
    e=c;
	  c=d+r*e;
	  d=cf[i]*(double)(i)*(double)(i+1)-s*e;
  }
  return -(c*w+d)*w*w;

}	
/***********************************************************************
				Taylor expansion
***********************************************************************/
cmplex Taylor(const cmplex &Z, const int N,const cmplex * const &cf){
	//horners rule
	static cmplex sum;
  if (N==0){return cf[0];}
  sum=cf[N-1];
  for (int n=N-2; n>=0; n--){ sum=(sum*Z)+cf[n]; }
  return sum;
}
//-----------------------------------------------------------------------
cmplex TaylorDer(const cmplex &Z,const int N,const cmplex * const &cf,const double &r ){
	//returns NEGATIVE derivative
	static cmplex sum;
  if (N==0){return 0.0;}
  sum=(double)(N-1)*cf[N-1];
  for (int n=N-2; n>=1; n--){ sum=(sum*Z)+(cf[n]*(double)(n)); }
  return -sum/r;
}
//-----------------------------------------------------------------------
cmplex TaylorDxx(const cmplex &Z,const int N,const cmplex * const &cf,const double &r ){
	static cmplex sum;
  if (N==0){return 0.0;}
  sum=(double)((N-1)*(N-2))*cf[N-1];
  for (int n=N-2; n>=2; n--){ sum=(sum*Z)+(cf[n]*(double)((n)*(n-1))); }
  return sum/r/r; 
}
//************************************************************************
bool InSquare(const cmplex &z, const cmplex &squarecen, const double w){
	if ((z.real()> squarecen.real() + (w/2.0)) ||
      (z.real()< squarecen.real() - (w/2.0)) ||
      (z.imag()> squarecen.imag() + (w/2.0)) ||
      (z.imag()< squarecen.imag() - (w/2.0)))     {return false;}
  return true;
}
//------------------------------------------------------------------------
bool InEllipse(const cmplex &z, const cmplex &z1, const cmplex &z2, const double dist){
	if ((abs(z-z1)+abs(z-z2))<dist){return true;}
	else                           {return false;}
}
//------------------------------------------------------------------------
bool InTriangle(const cmplex &z, const cmplex &z1, const cmplex &z2, const cmplex &z3){

	double A1(TriArea(z, z1, z2));
	double A2(TriArea(z, z2, z3));
	double A3(TriArea(z, z3, z1));

	return (((A1 >= 0.0) && (A2 >= 0.0) && (A3 >= 0.0)) || ((A1 <= 0.0) && (A2 <= 0.0) && (A3 <= 0.0)));
}
//------------------------------------------------------------------------
intercepttype Intersect(const cmplex &z1a,const cmplex &z1b,const cmplex &z2a,const cmplex &z2b,cmplex &zint){
	//Adapted from work by Paul D. Bourke
	double ua,ub,den;

	den=((z2b - z2a).imag() * (z1b - z1a).real() - (z2b - z2a).real() * (z1b - z1a).imag());

	if (den != 0.0) {
		ua = ((z2b - z2a).real() * (z1a - z2a).imag() - (z2b - z2a).imag() * (z1a - z2a).real()) / den;
		ub = ((z1b - z1a).real() * (z1a - z2a).imag() - (z1b - z1a).imag() * (z1a - z2a).real()) / den;
	}
	else{
		ua=0.0;
		ub=0.0;
	}
	zint = z1a + ua * (z1b - z1a);
	if ((ua>0.0) && (ua<1.0) && (ub>0.0) && (ub < 1.0)) {
		return INTERSECTED;
	}
	else if ((ub==0.0) || (ub==1.0) || (ua==0.0) || (ua==1.0)){
		return INTNODE;
	}
	else {
		return NO_INTERSECT;
	}
}
//------------------------------------------------------------------------
intercepttype LineIntersectCircle(const cmplex &z1, 
																	const cmplex &z2, 
																	const cmplex &zcen, 
																	const double  R, 
																	      cmplex &zint1, 
																	      cmplex &zint2){

  //Adapted from work by Paul D. Bourke
	double A,B,C,u1,u2,den;
	int    nint(0);

	//A=(x2 - x1) ^ 2 + (y2 - y1) ^ 2
	A = abs(z2-z1)*abs(z2-z1);

	//B = 2 * (((x2 - x1) * (x1 - Xc)) + ((y2 - y1) * (y1 - Yc)))
	B=2.0*(((z2.real()-z1.real())*(z1.real()-zcen.real()))+
		     ((z2.imag()-z1.imag())*(z1.imag()-zcen.imag())));

	//C = (Xc ^ 2) + (Yc ^ 2) + (x1 ^ 2) + (y1 ^ 2) - 2 * ((Xc * x1) + (Yc * y1)) - R ^ 2
	C = (zcen.real()*zcen.real()) + 
		  (zcen.imag()*zcen.imag()) + 
			(  z1.real()*  z1.real()) + 
			(  z1.imag()*  z1.imag()) - 2.0 * ((zcen.real() * z1.real()) + (zcen.imag() * z1.imag())) - R*R;
 
	den=B*B-4*A*C;

	if (den>=0.0 && (A!=0.0)){
		u1 = (-B + sqrt(den)) / (2.0*A);
    u2 = (-B - sqrt(den)) / (2.0*A);

		if ((u1>0.0) && (u1<1.0)){
      if (nint==0){zint1 = z1 + u1 * (z2 - z1);}
			nint++;
		}
    if ((u2>0.0) && (u2<1.0)){
      if (nint==0){zint1 = z1 + u2 * (z2 - z1);}
		  else        {zint2 = z1 + u2 * (z2 - z1);}
			nint++;
		}
	}
	if      (nint==0){zint1=zint2=0.0; return NO_INTERSECT;}
	else if (nint==1){      zint2=0.0; return INTERSECTED;}
	else if (nint==2){                 return TWO_INTERSECTIONS;}   
	else             {                 return NO_INTERSECT;}
	
}
/********************************************************************
IntegrateLine
---------------------------------------------------------------------
calculates positive and negative integral of line defined by (x1,y1) and (x2,y2)
positive is, of course >=0
negative is, of course <=0
returns total integral positive+negative
********************************************************************/
double IntegrateLine(const double &x1, const double &x2,
										 const double &y1, const double &y2, 
										 double &positive, double &negative){

	if (x1>x2){ExitGracefully("IntegrateLine: bad interval",RUNTIME_ERR);}
	if      ((y1>=0.0) && (y2>=0.0)){positive=(y1+y2)*(x2-x1)/2.0;negative=0.0;}//line wholly above x-axis
	else if ((y1<=0.0) && (y2<=0.0)){negative=(y1+y2)*(x2-x1)/2.0;positive=0.0;}//line wholly below x-axis
	else{
		double xint=x1-(x2-x1)*y1/(y2-y1); //find intercept-use area of triangles
		
		if (y1>y2){
			positive=y1*(xint-x1  )/2.0;
			negative=y2*(x2  -xint)/2.0;
		}
		else{
			positive=y2*(x2  -xint)/2.0;
			negative=y1*(xint-x1  )/2.0;
		}
	}

	return positive+negative;
}		
/********************************************************************
Polyfit
---------------------------------------------------------------------
calculates coefficients a0..aN-1 of polynomial that fits through points 
(xi,yi) where i=0..N-1
********************************************************************/
bool            PolyFit(Ironclad1DArray   x,
												Ironclad1DArray   y,
												Writeable1DArray  a,
												const int         N){
	int n,n2;
  double tmp;
	bool worked;
	ExitGracefullyIf(N<=0,"PolyFit::Must fit polynomial where N>=1",BAD_DATA);
	
	if (N==1){a[0]=y[0];return true;}
  else {
		//create, initialize matrices
		double **A  =new double *[N];
		double  *b  =new double  [N];
		for (n=0; n<N; n++){
			A[n]=new double [N];
			b[n]=y[n];
			for(n2=0; n2<N; n2++){A[n][n2]=0.0;}
		}

		for (n=0; n<N; n++){
			A[n][0]=1;
			tmp=1;
			for (n2=1; n2<N; n2++){
				tmp*=x[n];
				A[n][n2]=tmp;
			}
		}

		if (Gauss(A,b,a,N)){worked=true;}
		else               {worked=false;}

		for (n=0; n<N; n++){delete A[n];} delete A;
	}
  return worked;
}
/********************************************************************
PolyEval
---------------------------------------------------------------------
Evaluates polynomial with coefficients a0..aN-1 at point x
uses Knuth's rule
********************************************************************/
double         PolyEval(const double      x,
												Ironclad1DArray   a,
												const int         N){
	double sum,tmp;

	sum=a[0];
	tmp=1;
	for (int n=1; n<N; n++){
		tmp*=x;
		sum+=a[n]*tmp;
	}
	return sum;
}
/**************************************************************************
	TriArea:
		Computes Trangular area between 3 points
		positive if counterclockwise, negative if clockwise
-------------------------------------------------------------------------*/
double TriArea(const cmplex &z1, const cmplex &z2,const cmplex &z3){
	return 0.5 * ((z2.real()-z1.real())*(z3.imag()-z1.imag()) - 
		            (z2.imag()-z1.imag())*(z3.real()-z1.real()));
}
//**************************************************************************
bool DynArrayAppend(void **& pArr,void *xptr,int &size){  
	void **tmp=NULL;
  if (xptr==NULL){return false;}
	if ((pArr==NULL) && (size>0)) {return false;}
  size=size+1;																								//increment size
  tmp=new void *[size];					                              //allocate memory 
	if (tmp==NULL){ExitGracefully("DynArrayAppend::Out of memory",OUT_OF_MEMORY);}
  for (int i=0; i<(size-1); i++){                             //copy array
		tmp[i]=pArr[i];
		if (pArr[i]==NULL){ExitGracefully("DynArrayAppend::Bad existing array",BAD_DATA);}
	}	                             
  tmp[size-1]=xptr;																		        //add new pointer
  if (size>1){delete [] pArr; }										            //delete old array of pointers														
	pArr=tmp;																				            //redirect pointer
	return true;
}
//-------------------------------------------------------------------------
void   StraightSort(Writeable1DArray a, const int size){
	double tmp;
	int i;
	for (int j=1; j<size; j++){
		tmp=a[j];
		i=j;
		while (i>0 && (a[i-1]>tmp)){
			a[i]=a[i-1];
			i--;
		}
		a[i]=tmp;
	}
}
//*********************************************************************
void line_dxf(ofstream &DXF,
							double x1, 
							double y1, 
							double z1, 
							double x2, 
							double y2, 
							double z2, 
							char  *layer,
							int    color){

	DXF<<"0"   <<endl<<"LINE"<<endl;
  DXF<<"8"   <<endl<< layer<<endl;
	DXF<<"62"  <<endl<< color<<endl;
	DXF<<"10"  <<endl<< x1   <<endl;
	DXF<<"20"  <<endl<< y1   <<endl;
	DXF<<"30"  <<endl<< z1   <<endl;
	DXF<<"11"  <<endl<< x2   <<endl;
	DXF<<"21"  <<endl<< y2   <<endl;
	DXF<<"31"  <<endl<< z2   <<endl;
}
//******************************************************************************
void randinit(long init){
  idum=new long;
	(*idum)=init;
}
//******************************************************************************
double ran1(){
	//random number between 0 and 1 (from press et al)
	int j;
	long k;
	static long iy=false;
	static long iv[NTAB];
	float temp;

	if (((*idum)<=0) || (!iy)){         //initialization procedure
		if (-(*idum)<1){(*idum)=1;}
		else           {(*idum)=-(*idum);}
		for (j=NTAB+7; j>=0; j--){
			k=(*idum)/IQ;
			(*idum)=IA*((*idum)-k*IQ)-IR*k;
			if ((*idum)<0){(*idum)+=IMI;}
			if (j<NTAB){iv[j]=(*idum);}
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	(*idum)=IA*((*idum)-k*IQ)-IR*k;
	if ((*idum)<0){(*idum)+=IMI;}
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=(*idum);
	temp=AM*iy;

	if (temp>RNMAX){return RNMAX;}
	else           {return temp;}

}
//******************************************************************************
double gaussrand(){
  //returns random variable with mean of 0 and standard deviation of 1
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2,tmp;

	if ((*idum)<0){iset=0;}
	if (iset==0){
		do{
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			rsq=v1*v1+v2*v2;
		} while ((rsq>=1.0) || (rsq==0.0));
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		tmp= v2*fac;
	}
	else{
		iset=0;
		tmp= gset;
	}
	return tmp;
}
//******************************************************************************
double quickrandom (const double low,const double high){
	return (double)(rand())/(double)(RAND_MAX)*(high-low)+low;
}
//******************************************************************************

//******************************************************************************

void ExitGracefullyIf(bool condition, char *statement, badcode code){
	if (condition==true){ExitGracefully(statement,code);}
}
//******************************************************************************
void ExitGracefully(char *statement, badcode code){
	int thiscode(0);
	switch (code){	
		case(SINGMAT):      {thiscode=0;break;}
		case(COMPLETED):    {thiscode=1;break;}
		case(RUNTIME_ERR):  {thiscode=-1;break;}
		case(BAD_DATA):     {thiscode=-1;break;}
		case(TOO_MANY):	    {thiscode=-2;break;}
		case(BAD_SOL_READ): {thiscode=-6;break;}
		case(OTHER):        {thiscode=-11;break;}
    case(NO_INPUT_FILE):{thiscode=-15;break;}
		default:            {thiscode=-11;break;}
	}
	cout <<endl<<endl;
	cout <<"****************************************************"<<endl;
	cout <<"Exiting Gracefully: "<< statement <<endl;
	cout <<"****************************************************"<<endl;
	ofstream GLOBALDEBUG;
	GLOBALDEBUG.open("debug.out");
	GLOBALDEBUG<<thiscode<<endl<<statement<<endl;
	GLOBALDEBUG.close();
	bool pause(true);
	if (code==COMPLETED){pause=false;}
	GlobalDelete(pause);
	exit(0);
}
//******************************************************************************
bool ProgramAborted(){
	ifstream STOP;
	STOP.open("stop",ios::in);
	if (STOP.is_open())  {
		cout << "STOP file generated. Please Wait. "<<endl;
		STOP.close();
		return true;
	}
	else {return false;}
}

//******************************************************************************
int GetNumTriGaussPts(const trigausspoints gausspts){
	switch(gausspts){
		case(TG_1POINT):    return 1;  break;
		case(TG_3POINT):    return 3;  break;
		case(TG_3POINTB):   return 3;  break;
		case(TG_3ENDPOINT): return 3;  break;
		case(TG_4POINT):    return 4;  break;
		case(TG_6POINT):    return 6;  break;
		case(TG_7POINT):    return 7;  break;
		case(TG_9POINT):    return 9;  break;
		case(TG_13POINT):   return 13; break;
		case(TG_16POINT):   return 16; break;
	}
	ExitGracefully("GetNumTriGaussPts",RUNTIME_ERR);return 0;
}
//------------------------------------------------------------------------
int GetNumLinIntegrationPts(const lin_integrationtype gausspts){

	switch(gausspts){
		case(LG_1POINT):      return 1;  break;
		case(LG_2POINT):      return 3;  break;
		case(LG_3POINT):      return 3;  break;
		case(LG_4POINT):      return 4;  break;
		case(LG_5POINT):      return 5;  break;
		case(LG_TRAPEZOID1):  return 2;  break;
		case(LG_TRAPEZOID2):  return 3;  break;
		case(LG_TRAPEZOID3):  return 4;  break;
		case(LG_TRAPEZOID10): return 11; break;
		case(LG_TRAPEZOID100):return 101;break;
	}
	ExitGracefully("GetNumLinIntegrationPts: integration method is invalid",BAD_DATA);
	return 0;
}
//------------------------------------------------------------------------
void GetLinIntegrationInfo(const double &x1,
													 const double &x2,
													 const lin_integrationtype gausspts, 
													 Writeable1DArray x,  //integration points
													 Writeable1DArray w){
	//Gauss Legendre Weights & Points
	switch(gausspts){
		case(LG_1POINT):   //linear
		  w[0]=1.0;						x[0]=x1+0.5*(x2-x1);  
			break;
		//Gauss Legendre---------------------------------------
		case(LG_2POINT):   
			w[0]=0.5;						x[0]=x1+0.21132486541*(x2-x1);
			w[1]=0.5;           x[1]=x1+0.78867513459*(x2-x1);
			break;
		case(LG_3POINT):  
      w[0]=0.27777777778; x[0]=x1+0.11270166538*(x2-x1);
			w[1]=0.44444444444; x[1]=x1+0.50000000000*(x2-x1);
			w[2]=0.27777777778; x[2]=x1+0.88729833462*(x2-x1);
			break;
		case(LG_4POINT):  //
      w[0]=0.17392742257;	x[0]=x1+0.06943184421*(x2-x1);
			w[1]=0.32607257743;	x[1]=x1+0.33000947821*(x2-x1);
			w[2]=0.32607257743;	x[2]=x1+0.66999052179*(x2-x1);
			w[3]=0.17392742257;	x[3]=x1+0.93056815580*(x2-x1);
			break;
		case(LG_5POINT): //
      w[0]=0.11846344253;	x[0]=x1+0.04691007703*(x2-x1);
			w[1]=0.23931433525;	x[1]=x1+0.23076534495*(x2-x1);
			w[2]=0.28444444444;	x[2]=x1+0.50000000000*(x2-x1);
			w[3]=0.23931433525;	x[3]=x1+0.76923465505*(x2-x1);
			w[4]=0.11846344253;	x[4]=x1+0.95308992297*(x2-x1);
			break;
		//Trapezoid Rule----------------------------------------
		case(LG_TRAPEZOID1): //
			w[0]=0.5;					  x[0]=x1;
			w[1]=0.5;						x[1]=x2;
			break;
		case(LG_TRAPEZOID2): //
			w[0]=0.25;				  x[0]=x1;
			w[1]=0.50;					x[1]=x1+0.5*(x2-x1);
      w[2]=0.25;					x[2]=x2;
			break;
		case(LG_TRAPEZOID3): //
			w[0]=0.1666666667;	x[0]=x1;
			w[1]=0.3333333333;	x[1]=x1+0.33333333333*(x2-x1);
			w[2]=0.3333333333;	x[2]=x1+0.66666666667*(x2-x1);
      w[3]=0.1666666667;	x[3]=x2;
			break;
		case(LG_TRAPEZOID10): //
			w[0]=0.0500000000;	x[0]=x1;
			w[1]=0.1000000000;	x[1]=x1+0.1*(x2-x1);
			w[2]=0.1000000000;	x[2]=x1+0.2*(x2-x1);
			w[3]=0.1000000000;	x[3]=x1+0.3*(x2-x1);
			w[4]=0.1000000000;	x[4]=x1+0.4*(x2-x1);
			w[5]=0.1000000000;	x[5]=x1+0.5*(x2-x1);
			w[6]=0.1000000000;	x[6]=x1+0.6*(x2-x1);
			w[7]=0.1000000000;	x[7]=x1+0.7*(x2-x1);
			w[8]=0.1000000000;	x[8]=x1+0.8*(x2-x1);
			w[9]=0.1000000000;	x[9]=x1+0.9*(x2-x1);
      w[10]=0.050000000;	x[10]=x2;
			break;
		case(LG_TRAPEZOID100): //
			w[0]  =0.0050000000;	x[0]  =x1;
			for (int i=1; i<100;i++){
			  w[i]=0.0100000000;	x[1]  =x1+double(i)*0.01*(x2-x1);
			}
      w[100]=0.0050000000;	x[100]=x2;
		break;
	}
}
//------------------------------------------------------------------------
void    GetTriGaussInfo(const cmplex &z1,
												const cmplex &z2,
												const cmplex &z3,
												const trigausspoints gausspts, 
												cmplex       *z, 
												double       *w){

	static double a1,a2,a3,a4,b1,b2,b3,b4,c3,c4,t1,t2,t3,t4,t5;

/*   H R Schwarz,   Methode der Finiten Elemente, Teubner Studienbuecher, 1980.
     Strang and Fix, An Analysis of the Finite Element Method, Prentice Hall, 1973, page 184.
     Arthur H Stroud, Approximate Calculation of Multiple Integrals, Prentice Hall, 1971.
     O C Zienkiewicz, The Finite Element Method, McGraw Hill, Third Edition, 1977, page 201.*/

	switch(gausspts){
		case(TG_1POINT):   //linear
		  w[0]=1.0;        z[0]=TriLocalToGlobal(z1,z2,z3,1.0/3.0,1.0/3.0,1.0/3.0);  
			break;
		case(TG_3POINT):   //quadratic
			w[0]=1.0/3.0;    z[0]=TriLocalToGlobal(z1,z2,z3, 0.5, 0.5, 0.0);
			w[1]=1.0/3.0;    z[1]=TriLocalToGlobal(z1,z2,z3, 0.5, 0.0, 0.5);
			w[2]=1.0/3.0;    z[2]=TriLocalToGlobal(z1,z2,z3, 0.0, 0.5, 0.5);
			break;
		case(TG_3POINTB):  //quadratic
      w[0]=1.0/3.0;    z[0]=TriLocalToGlobal(z1,z2,z3,2.0/3.0, 1.0/6.0, 1.0/6.0);
			w[1]=1.0/3.0;    z[1]=TriLocalToGlobal(z1,z2,z3,1.0/6.0, 2.0/3.0, 1.0/6.0);
			w[2]=1.0/3.0;    z[2]=TriLocalToGlobal(z1,z2,z3,1.0/6.0, 1.0/6.0, 2.0/3.0);
			break;
		case(TG_3ENDPOINT)://linear-poor approximation!
      w[0]=1.0/3.0;    z[0]=TriLocalToGlobal(z1,z2,z3, 1.0, 0.0, 0.0);
			w[1]=1.0/3.0;    z[1]=TriLocalToGlobal(z1,z2,z3, 0.0, 1.0, 0.0);
			w[2]=1.0/3.0;    z[2]=TriLocalToGlobal(z1,z2,z3, 0.0, 0.0, 1.0);
			break;
		case(TG_4POINT):  //Strang and Fix formula #3.
		  w[0]=-9.0/16.0;  z[0]=TriLocalToGlobal(z1,z2,z3, 1.0/3.0, 1.0/3.0, 1.0/3.0);
		  w[1]=25.0/48.0;  z[1]=TriLocalToGlobal(z1,z2,z3, 3.0/5.0, 1.0/5.0, 1.0/5.0);
		  w[2]=25.0/48.0;  z[2]=TriLocalToGlobal(z1,z2,z3, 1.0/5.0, 3.0/5.0, 1.0/5.0);
		  w[3]=25.0/48.0;  z[3]=TriLocalToGlobal(z1,z2,z3, 1.0/5.0, 1.0/5.0, 3.0/5.0);
			break;
		case(TG_6POINT):  //Strang and Fix, formula #5.
			a1 = 0.816847572980459;    b1 = 0.091576213509771;
			a2 = 0.108103018168070;    b2 = 0.445948490915965;
			t1 = 0.109951743655322;    t2 = 0.223381589678011;		

			w[0]=t1;         z[0]=TriLocalToGlobal(z1,z2,z3, a1, b1, b1); 
			w[1]=t1;         z[1]=TriLocalToGlobal(z1,z2,z3, b1, a1, b1);
			w[2]=t1;         z[2]=TriLocalToGlobal(z1,z2,z3, b1, b1, a1);
			w[3]=t2;         z[3]=TriLocalToGlobal(z1,z2,z3, a2, b2, b2);
			w[4]=t2;         z[4]=TriLocalToGlobal(z1,z2,z3, b2, a2, b2);
			w[5]=t2;         z[5]=TriLocalToGlobal(z1,z2,z3, b2, b2, a2);
      break;
		case(TG_7POINT): //"quintic"
			a1=0.05961587;  b1=0.5*(1.0-a1);//0.47014206;
			
			a2=0.79742699;  b2=0.5*(1.0-a2);//0.10128651;
			t1=0.225;       t2=0.13239415;   t3=(1.0-t1-3.0*t2)/3.0;//0.12593918;
      //cout <<"QUINTIC"<<endl;
      w[0]=t1;         z[0]=TriLocalToGlobal(z1,z2,z3,1.0/3.0,1.0/3.0,1.0/3.0);
			w[1]=t2;         z[1]=TriLocalToGlobal(z1,z2,z3,a1,b1,b1);
			w[2]=t2;         z[2]=TriLocalToGlobal(z1,z2,z3,b1,a1,b1);
			w[3]=t2;         z[3]=TriLocalToGlobal(z1,z2,z3,b1,b1,a1);
			w[4]=t3;         z[4]=TriLocalToGlobal(z1,z2,z3,a2,b2,b2);
			w[5]=t3;         z[5]=TriLocalToGlobal(z1,z2,z3,b2,a2,b2);
			w[6]=t3;         z[6]=TriLocalToGlobal(z1,z2,z3,b2,b2,a2);
			break;
		case(TG_9POINT): //Strang and Fix formula #8.

			a1 = 0.124949503233232;  b1 = 0.437525248383384;
			a3 = 0.797112651860071;  b3 = 0.165409927389841;  c3 = 0.037477420750088;
      t1 = 0.205950504760887;  t3 = 0.063691414286223;

			w[0]=t1;         z[0]=TriLocalToGlobal(z1,z2,z3,a1,b1,b1);
			w[1]=t1;         z[1]=TriLocalToGlobal(z1,z2,z3,b1,a1,b1);
			w[2]=t1;         z[2]=TriLocalToGlobal(z1,z2,z3,b1,b1,a1);
			w[3]=t3;         z[3]=TriLocalToGlobal(z1,z2,z3,a3,b3,c3);
			w[4]=t3;         z[4]=TriLocalToGlobal(z1,z2,z3,b3,c3,a3);
			w[5]=t3;         z[5]=TriLocalToGlobal(z1,z2,z3,c3,a3,b3);
			w[6]=t3;         z[6]=TriLocalToGlobal(z1,z2,z3,a3,c3,b3);
			w[7]=t3;         z[7]=TriLocalToGlobal(z1,z2,z3,b3,a3,c3);
			w[8]=t3;         z[8]=TriLocalToGlobal(z1,z2,z3,c3,b3,a3);
      break;
		case(TG_13POINT)://Strang and Fix, formula #10.
			
			a1 = 0.479308067841923;  b1 = 0.260345966079038;
			a2 = 0.869739794195568;  b2 = 0.065130102902216;
			a3 = 0.638444188569809;  b3 = 0.312865496004875;  c3 = 0.048690315425316;
			t1 = 0.175615257433204;  t2 = 0.053347235608839;
			t3 = 0.077113760890257;  t4 =-0.149570044467670;

			w[0] =t4;        z[0] =TriLocalToGlobal(z1,z2,z3,1.0/3.0,1.0/3.0,1.0/3.0);
			w[1] =t1;        z[1] =TriLocalToGlobal(z1,z2,z3,a1,b1,b1);
			w[2] =t1;        z[2] =TriLocalToGlobal(z1,z2,z3,b1,a1,b1);
			w[3] =t1;        z[3] =TriLocalToGlobal(z1,z2,z3,b1,b1,a1);

			w[4] =t2;        z[4] =TriLocalToGlobal(z1,z2,z3,a2,b2,b2);
			w[5] =t2;        z[5] =TriLocalToGlobal(z1,z2,z3,b2,a2,b2);
			w[6] =t2;        z[6] =TriLocalToGlobal(z1,z2,z3,b2,b2,a2);
			
			w[7] =t3;        z[7] =TriLocalToGlobal(z1,z2,z3,a3,b3,c3);
			w[8] =t3;        z[8] =TriLocalToGlobal(z1,z2,z3,b3,c3,a3);
			w[9] =t3;        z[9] =TriLocalToGlobal(z1,z2,z3,c3,a3,b3);
			w[10]=t3;        z[10]=TriLocalToGlobal(z1,z2,z3,a3,c3,b3);
			w[11]=t3;        z[11]=TriLocalToGlobal(z1,z2,z3,b3,a3,c3);
			w[12]=t3;        z[12]=TriLocalToGlobal(z1,z2,z3,c3,b3,a3);
			break;
    case(TG_16POINT)://	http://www.math.niu.edu/~rusin/known-math/99/int_triangle    

			a1=0.0814148234145540; b1=0.4592925882927230;
			a2=0.6588613844964800; b2=0.1705693077517600;
			a3=0.8989055433659380; b3=0.0505472283170310;
			a4=0.0083947774099580; b4=0.2631128296346380; c4=0.7284923929554040;
			
			t1=0.0950916342672850;
			t2=0.1032173705347180;
			t3=0.0324584976231980;
			t4=0.0272303141744350;
			t5=0.1443156076777870;

			w[0] =t5;        z[0] =TriLocalToGlobal(z1,z2,z3,1.0/3.0,1.0/3.0,1.0/3.0);

	    w[1] =t1;        z[1] =TriLocalToGlobal(z1,z2,z3,a1,b1,b1);
			w[2] =t1;        z[2] =TriLocalToGlobal(z1,z2,z3,b1,a1,b1);
			w[3] =t1;        z[3] =TriLocalToGlobal(z1,z2,z3,b1,b1,a1);
      
			w[4] =t2;        z[4] =TriLocalToGlobal(z1,z2,z3,a2,b2,b2);
			w[5] =t2;        z[5] =TriLocalToGlobal(z1,z2,z3,b2,a2,b2);
			w[6] =t2;        z[6] =TriLocalToGlobal(z1,z2,z3,b2,b2,a2);
         
			w[7] =t3;        z[7] =TriLocalToGlobal(z1,z2,z3,a3,b3,b3);
			w[8] =t3;        z[8] =TriLocalToGlobal(z1,z2,z3,b3,a3,b3);
			w[9] =t3;        z[9] =TriLocalToGlobal(z1,z2,z3,b3,b3,a3);		  
       
			w[10]=t4;        z[10]=TriLocalToGlobal(z1,z2,z3,a4,b4,c4);
			w[11]=t4;        z[11]=TriLocalToGlobal(z1,z2,z3,b4,c4,a4);
			w[12]=t4;        z[12]=TriLocalToGlobal(z1,z2,z3,c4,a4,b4);
			w[13]=t4;        z[13]=TriLocalToGlobal(z1,z2,z3,a4,c4,b4);
			w[14]=t4;        z[14]=TriLocalToGlobal(z1,z2,z3,b4,a4,c4);
			w[15]=t4;        z[15]=TriLocalToGlobal(z1,z2,z3,c4,b4,a4);
			break;
	}
}
/***************************************************************************
                 TRIANGLE COORDINATE TRANSFORMATION
****************************************************************************
	Converts to barycentric coordinates
--------------------------------------------------------------------------*/
void   TriGlobalToLocal(const cmplex &z1 , const cmplex &z2,  const cmplex &z3,
												const cmplex &z, 
															double &xii,       double &xij,       double &xik){

	static cmplex zk,zloc;
	static double aloc;

	//correct element-local coordinate system:
	//zi=cmplex(0.0,0.0), zj=cmplex(1.0,0.0);
	//(13 adds/subtracts, 14 mults, 3 divisions)
	zk  =(z3-z1)/(z2-z1);
	zloc=(z -z1)/(z2-z1);
	aloc=0.5*zk.imag();

	//3 divides, 7 mults, 4 adds/subs
	xii=0.5/aloc*(zk.imag()-zk.imag()*zloc.real()+(zk.real()-1.0)*zloc.imag());
	xij=0.5/aloc*(zk.imag()*zloc.real()-zk.real()*zloc.imag());
	xik=0.5/aloc*(zloc.imag());

	//a bit faster? (one divide, 5 mults, 4 adds/subs), :
	//static double yk,yl;
	//yk=1.0/zk.imag();
	//yl=zloc.imag();
	//xii=(1.0-zloc.real())+yl*(zk.real()-1.0)*yk;
	//xij=zloc.real()-yl*zk.real()*yk;
	//xik=yl*yk;	  

	//even faster? bit faster (one divide, 3 mults, 3 adds/subs), :
	//static double yk,yl;
	//yk=1.0/zk.imag();
	//yl=zloc.imag();
	//xij=zloc.real()-yl*zk.real()*yk;
	//xik=yl*yk;	  
	//xii=1-xij-xik;
		
	 //Alternate slow formulation (36 adds/subtracts 24 mults 3 divisions)
	/*xii=((z-z3)*conj(z3-z2)).imag()/((z1-z3)*conj(z3-z2)).imag();
	xij=((z-z1)*conj(z1-z3)).imag()/((z2-z1)*conj(z1-z3)).imag();
	xik=((z-z2)*conj(z2-z1)).imag()/((z3-z2)*conj(z2-z1)).imag();*/
}
//******************************************************************************
cmplex arccosh(const cmplex &z){
	if (z.imag()>0){
		return log(z+pow(z*z-1.0,0.5));
	}
	else{
		return log(z-pow(z*z-1.0,0.5));
	}
}
//------------------------------------------------------------------------------
cmplex arccosz(const cmplex &z){
	double x=z.real();
	double y=z.imag();
	double a=0.5*sqrt((x+1.0)*(x+1.0)+y*y);
	double b=0.5*sqrt((x-1.0)*(x-1.0)+y*y);
	double X=a+b;
	
	return cmplex (acos(a-b),-y/fabs(y)*log(X+sqrt(X*X-1)));

}
//cmplex ChebyshevZ(const cmplex &z){
	//should only be used to evaluate single terms
//}
//------------------------------------------------------------------------------
double my_erfc(const double &x){
	double tmp(fabs(x)); //take abs so that we are always in positive quadrant.
	double fun;
	double f1;
	double tmp2;
	double tmp3;

	if(tmp > 3.0){
		f1  = (1.0 - 1.0/(2.0 * tmp * tmp) 
			       + 3.0/(4.0 * pow(tmp,4)) 
			       - 5.0/(6.0 * pow(tmp,6)));
		fun = f1 * exp(-tmp * tmp) / (tmp * sqrt(PI));
	} 
	else{
		tmp2 = 1.0 / (1.0 + (0.3275911 * tmp));
		tmp3 =   0.254829592  * tmp2       //5th order polynomial interpolation
			   - (0.284496736 * tmp2 * tmp2) 
			   + (1.421413741 * pow(tmp2,3))
			   - (1.453152027 * pow(tmp2,4))
			   + (1.061405429 * pow(tmp2,5));
		fun = tmp3 * exp(-tmp * tmp);
	} 
	if (tmp == x) {return fun;}
	else{return (2-fun);}
} 
//------------------------------------------------------------------------------
double my_erf(const double &x){
	return 1-my_erfc(x);
}
//------------------------------------------------------------------------------
double hantush(const double &u, const double b){
  //From Charbeneau, 2000 pg. 539
	double eps=0.0000001;
	double oldst=-1e30;
	double oldsum =-1e30;
	double sum,st,del,x;
	int i,j,it;
	double tnm;

	double v=u;

	if (u>(b/2.0)){v=0.25*b*b/u;}
	
	j=1;
	while (j<=40) {
		if (j==1){
			st=0.5*exp(-(v+0.25*b*b/v));
		}
		else{
			it=ipow(2,j-2);
			tnm=it;
			del=v/tnm;
			x=del/2.0;
			sum=0.0;
			i=1;
			while (i<=it){
				sum+=exp(-(x+0.25*b*b/x))/x;
				x  +=del;
				i++;
			}
			st=(st+v*sum/tnm)/2.0;
		}
		sum=(4.0*st-oldst)/3.0;
		if (fabs(sum-oldsum)<eps){break;}
		oldsum=sum;
		oldst=st;
		j++;
	}
	if (u==v){return 2.0*bessel_K0(b)-sum;}
	else     {return sum;}
}
//------------------------------------------------------------------------------
double bessel_I0(const double &x){
	//modified bessel function of the first kind, Ko
	//From Charbeneau, 2000 pg. 53
	double ax,y;

	ax=fabs(x);
	if (ax<3.75){
		y=x/3.75;
		y*=y;
		return 1.0+y*(3.5156229+
			         y*(3.0899424+
							 y*(1.2067492+
							 y*(0.2659732+
							 y*(0.0360768+
							 y*(0.0045813))))));}
	else{
		y=3.75/ax;
		return exp(ax)/sqrt(ax)*(0.39894228+y*( 0.01328592+
			                                  y*( 0.00225319+
																				y*(-0.00157565+
																				y*( 0.00916281+
																				y*(-0.02057706+
																				y*( 0.02635537+
																				y*(-0.01647633+
																				y*( 0.00392377)))))))));}
}
//----------------------------------------------------------------------------
double bessel_K0(const double &x){
	//modified bessel function of the second kind, Ko
	double y;

	if (x<=2.0){
		y=0.25*x*x;
		return (-log(0.5*x)*bessel_I0(x))+(-EULER+y*(0.42278420+
																						  y*(0.23069756+
																						  y*(0.03488590+
																						  y*(0.00262698+
																						  y*(0.00010750+
																						  y* 0.0000074  ))))));}
	else {
		y=2.0/x;
		return exp(-x)/sqrt(x)*(1.25331414+y*(-0.07832358+
																			 y*( 0.02189568+
																			 y*(-0.01062446+
																			 y*( 0.00587872+
																			 y*(-0.00251540+
																			 y*( 0.00053208)))))));}

}
//------------------------------------------------------------------------------
double expint(const double &x){
	//computes exponential integral for x>0 (revised from Press et. al. Numerical recipes)

	int k,maxiter(100);
	double fact,prev,sum(0.0),term;

	if (x<=0.0){
		ExitGracefully("Negative Argument in exponential integral",RUNTIME_ERR);}
	if      (x<REALSMALL){
		return log(x)+EULER;
	} 
	else if (x<=-log(REALSMALL)){
		fact=1.0;
		for (k=1;k<=maxiter;k++){
			fact*=x/k;
			term=fact/k;
			sum+=term;
			if (term<REALSMALL*sum){break;}
		}
		return sum+log(x)+EULER;
	}
	else {
		term=1.0;
		for (k=1;k<=maxiter;k++){
			prev=term;
			term*=k/x;
			if (term<REALSMALL){break;}
			if (term<prev)     {sum+=term;}
			else               {sum-=prev; break;}
		}
		return exp(x)*(1.0+sum)/x;
	}
	return 0.0;
}
//------------------------------------------------------------------------------
double MQInterpolate(Ironclad1DArray   valctrl,
									   Ironclad1DArray_z zctrl,
									   int               numctrl,
									   Writeable1DArray  MQcoeff,
									   Ironclad1DArray_z zbasis,
									   int               order, 
										 double            R,
									   double           &ave_val){
//returns max error
	int m,n,n2,s,rank;
	double rs,rn;

	if (order<=1){
    ave_val=0.0;
    for(m=0; m<numctrl; m++){ave_val+=valctrl[m]/double(numctrl);}
		if (order==0){return 0.0;}
    if (order==1){MQcoeff[0]=0;}
  }
  else{  
		//create, initialize matrices
		double **A  =new double *[order+1];
		double  *b  =new double  [order+1];
		double  *sol=new double  [order+1];
		for (n=0; n<=order; n++){
			A[n]=new double [order+1];
			for(n2=0; n2<=order; n2++){A[n][n2]=0.0;}
      b[n]=0.0;
			sol[n]=0.0;
		}
    for (m=0; m<numctrl; m++){
      for (s=0; s<order; s++){
	     	rs=sqrt(pow(abs(zbasis[s]-zctrl[m]),2)+R*R);
				for (n=0; n<order; n++){
					//rn=abs(zbasis[n]-zctrl[m]);
					rn=sqrt(pow(abs(zbasis[n]-zctrl[m]),2)+R*R);
					A[s][n]+=(rs*rn);
				}
				A[s][order]+=rs;
				b[s]       +=rs*valctrl[m]; 
      }
    }
    for (n=0; n<order; n++){
      A[order][n]=1.0;
    }
    A[order][order]=0;
    b[order]       =0;

		rank=0;
    if (!SVD(A,b,sol,(order+1))){ 
			delete [] sol; 
			delete [] b;
			for (n=0; n<=order; n++){delete [] A[n];} delete [] A; 
			ExitGracefully("MQInterpolate: SVD routine failed",SINGMAT);}
    
    //pick up solution
    for (n=0; n<order; n++){
			MQcoeff[n]=sol[n]; 
		}
		ave_val=sol[order]; 

		delete [] sol;
		delete [] b;
		for (n=0; n<=order; n++){delete [] A[n];} delete [] A; 
	}

	double maxerr(0.0),err,val,maxval(-ALMOST_INF),minval(ALMOST_INF);
  for (m=0; m<numctrl; m++){
		upperswap(maxval,valctrl[m]);
		lowerswap(minval,valctrl[m]);
		val=ave_val;
		for (n=0; n<order; n++){
			//rn=abs(zbasis[n]-zctrl[m]);
			rn=sqrt(pow(abs(zbasis[n]-zctrl[m]),2)+R*R);
			val+=MQcoeff[n]*rn;
		}
		err=fabs(valctrl[m]-val);
		upperswap(maxerr,err);
	}
	return maxerr/(maxval-minval);
}
