#include "StringElem.h"

/************************************************************************
                           CLineElem
************************************************************************	
	  		                   CONSTRUCTOR
************************************************************************/
CLineElem::CLineElem(){}
//------------------------------------------------------------------------
CLineElem::CLineElem(char						 *Name,
										 const CSingleLayerABC *pLay, 
										 const cmplex     end1, 
										 const cmplex     end2,
										 const linetype   LType,
		                 const int        prec)
					  :CAnalyticElem(Name,pLay){

	int n;
	double fold;

  //initialize variable values     
	ltype       =LType;
  z1=end1;
	z2=end2;

	SetPrecision(prec,order,fold);
	nlinecontrol=int(order*fold);

  c[0]=NULL;	   c[1]=NULL;      c[2]=NULL;       c[3]=NULL; 
	
  if (order>MAX_ORDER){ExitGracefully("StringElem::Order too large",TOO_MANY);}
  if (nlinecontrol>MAXLINECONTROL){ExitGracefully("StringElem::Too many control points on element",TOO_MANY);}

	//Create Dynamic Arrays 
	FForder=     MAX_FAR_FIELD;   
	disabled=    false;

	X=           new double   [nlinecontrol];
  zctrl=       new cmplex   [nlinecontrol]; 
  JumpCoeff=   new double   [order+1];
  FarFldCoeff= new double   [MAX_FAR_FIELD];

  //create matrices of control point location 
  Interpolate(X,nlinecontrol,-1,1);
  Interpolate(zctrl,nlinecontrol,z1,z2);

  //initialize all coeff to 0
  for (n=0; n<=order;  n++)  {   JumpCoeff[n]=0.0;}
  for (n=0; n<FForder; n++)  { FarFldCoeff[n]=0.0;}
}

//-----------------------------------------------------------------------
CLineElem::~CLineElem(){
  if (globaldebug){cout <<"   DESTROYING LINEELEM "<<name<<endl;}
	delete [] X;
  delete [] zctrl;
  delete [] JumpCoeff;
  delete [] FarFldCoeff;
} 
/***********************************************************************
				ACCESSOR FUNCTIONS
***********************************************************************/
bool    CLineElem::IsDisabled()            const {return disabled;}
//-----------------------------------------------------------------------
cmplex	CLineElem::GetZ(const int pointno) const {if (pointno==0){return z1;} else {return z2;}}
//-----------------------------------------------------------------------
bool    CLineElem::HasFlux()               const {return ltype==LINESINK;}   
/***********************************************************************
				MANIPULATOR FUNCTIONS
***********************************************************************/
void CLineElem::SetAsDisabled(){disabled=true;}

/***********************************************************************
				Get Discharge Potential
----------------------------------------------------------------------*/
cmplex CLineElem::GetDischargePotential(const cmplex &z,const double &t) const{

	if (disabled){return 0.0;}

  cmplex  omega(0.0),jump(0.0),logend1,logend2;
  cmplex  Z ((z-0.5*(z1+z2))/(0.5*(z2-z1)));
	double  modulus(abs(Z));
  int     n;
	cmplex  mult(1.0,0);

	if (ltype==DOUBLET){mult=-IM;}

	//far-field 
  if (modulus>ISFAR){ 
    omega=mult*LaurentRe(Z,(int)(FForder*pow(modulus/ISFAR,FF_POWER)),FarFldCoeff);
		omega+=GetLog(z,t);
  }  
  //near-field
  else if (((abs(Z-E1)>MOVE_DIST*abs(z2-z1)) && (abs(Z-E2)>MOVE_DIST*abs(z2-z1)))){
		cmplex logterm=log(Zm1oZp1(Z));
    for (n=0; n<=order; n++){
			jump+=JumpCoeff[n]*cClenshaw2(n,Z,f[n],logterm);
		}
    omega=mult*IOVER2PI*jump; 
		omega+=GetLog(z,t);
  }
  //singularities
  else {
		if (abs(Z-E1) == abs(Z-E2)){
			return 0.0;
			ExitGracefully("CLineElem::GetPotential: Bad (infinite) Input location",RUNTIME_ERR);}
		else {omega=GetDischargePotential(z+MOVE_DIST*abs(z1-z2),t);} //trouble for inhomogeneities-direction matters
	}

  return omega;
}

/***********************************************************************
				Get Discharge Function
----------------------------------------------------------------------*/
cmplex CLineElem::GetW(const cmplex &z,const double &t) const{
	if (disabled){return 0.0;}

  cmplex  W(0.0),jump(0.0);   
  cmplex  Z ((z-0.5*(z1+z2))/(0.5*(z2-z1)));
  cmplex  divisor(PI*(z2-z1));
	double  end1(0.0),end2(0.0),modulus(abs(Z)),sum;
  int     n,j;
	cmplex  mult(1.0,0);

  if (ltype==LINESINK){
		sum=-1.0;
    for (n=0; n<=order; n++){
			sum*=-1.0;
      end1+=JumpCoeff[n]; 
	    end2+=JumpCoeff[n]*sum; 
		}
  }
	else if (ltype==DOUBLET){mult=-IM;}

  //far-field
  if (modulus>ISFAR){ 
		W+=mult*2.0*PI*LaurentDerRe(Z,(int)(FForder*pow(modulus/ISFAR,FF_POWER)),FarFldCoeff)/divisor;
    if (ltype==LINESINK){
			W-=(end2/(Z+1.0)-end1/(Z-1.0))/divisor;
		}
  }  
  //near-field
  else if ((abs(Z-E1)>MOVE_DIST*abs(z1-z2)) && (abs(Z-E2)>MOVE_DIST*abs(z1-z2))){
		cmplex logterm(log(Zm1oZp1(Z)));
    for (n=0; n<=order; n++){
			for (j=0; j<n; j++)    {cheby[j]=g[n][j];}
			for (j=n-1; j>=1; j-=2){cheby[j]=2.0*(double)(n)*logterm;}
      if  ((n%2)==1)         {cheby[0]=    (double)(n)*logterm;}
			cheby[n]=2.0/((Z-1.0)*(Z+1.0));
      jump+=JumpCoeff[n]*cClenshaw(n,Z,cheby);
		}
		if (ltype==LINESINK){W-=(jump+end2/(Z+1.0)-end1/(Z-1.0))/divisor;}		
    else                {W-=mult*jump/divisor;                       }
  }
	
  //singularities
  else {
		if (abs(Z-E1) == abs(Z-E2)){
			return 0.0;
			ExitGracefully("CLineElem::GetW: Bad (infinite) Input location",RUNTIME_ERR);}
		else {W=GetW(z+cmplex(MOVE_DIST*abs(z1-z2),0),t);}
	}

  return W;
}
/***********************************************************************
				Get Discharge Derivative Function
----------------------------------------------------------------------*/
cmplex CLineElem::GetGx(const cmplex &z,const double &t) const{
	cmplex Gx(0.0);
	return Gx; //TMP DEBUG
}
/***********************************************************************
				Get Flux Through Face
----------------------------------------------------------------------*/
cmplex CLineElem::GetFluxThruFace(const cmplex &z1a,const cmplex &z2a,const double &t) const{ 
  static cmplex zint;
  //SHOULD OPTIMIZE -intersect may be too slow
	if (Intersect(z1a,z2a,z1,z2,zint)!=INTERSECTED){
		return GetDischargePotential(z1a,t).imag()-GetDischargePotential(z2a,t).imag();
    //return IM*conj(GetDischargePotential(z1a,t)-GetDischargePotential(z2a,t)); //complex flux thru face
	}
	else{
    return GetDischargePotential(z1a,t).imag()-GetDischargePotential(zint-(z2a-z1a)*REALSMALL,t).imag()-
           GetDischargePotential(z2a,t).imag()+GetDischargePotential(zint+(z2a-z1a)*REALSMALL,t).imag();
    /*return IM*conj(GetDischargePotential(z1a,t)-GetDischargePotential(zint-(z2a-z1a)*REALSMALL,t)-
                   GetDischargePotential(z2a,t)+GetDischargePotential(zint+(z2a-z1a)*REALSMALL,t));*/
	}
}
/***********************************************************************
				Get Leakage Function
----------------------------------------------------------------------*/
double CLineElem::GetLeakage (const cmplex &z, const double &t,const leak_type ltype) const{
	return 0.0;
}
/***********************************************************************
				Get Integrated Budget Function
----------------------------------------------------------------------*/
void CLineElem::GetIntegratedBudget(const cmplex &z1t, const cmplex &z2t, const cmplex &z3t, 
																		const double &t, double &inQ, double &outQ) const{
	inQ=outQ=0.0;
	if (ltype==LINESINK){
		cmplex zint1,zint2,zint;
		int  nint=0;		
		bool in1 =InTriangle(z1,z1t,z2t,z3t);
		bool in2 =InTriangle(z2,z1t,z2t,z3t);		

    //find intersections
		if (Intersect(z1t,z2t,z2,z1,zint)==INTERSECTED){zint1=zint;nint++;}
		if (Intersect(z2t,z3t,z2,z1,zint)==INTERSECTED){if (nint==0){zint1=zint;}else{zint2=zint;}nint++;}
		if (Intersect(z3t,z1t,z2,z1,zint)==INTERSECTED){if (nint==0){zint1=zint;}else{zint2=zint;}nint++;}
		
		if ((!in1) && (!in2) && (nint==0)){return;} //no element/triangle overlap

		double X1(-1.0),X2(1.0); //default- element totally contained
		
		if (nint>0){             //element partially contained
			if      ((nint==1) && (in1)){X1=abs(zint1-z1)/abs(z2-z1);}
			else if ((nint==1) && (in2)){X2=abs(zint1-z1)/abs(z2-z1);}
			else if  (nint==2)          {
				X1=abs(zint1-z1)/abs(z2-z1);
				X2=abs(zint2-z1)/abs(z2-z1);        
				if (X2<X1){double Xtmp=X1;X1=X2;X2=Xtmp;}
			}
		}
		//integrate cumulative extraction from X1 to X2
		//GetCumExtraction(i,X1,X2,inQ,outQ); //TMP DEBUG- not yet a function
	}
}
/**********************************************************************
				Get Segment Net Discharge Function
----------------------------------------------------------------------*/
double CLineElem::GetNetDischarge(const double &t) const{
	if (disabled){return 0.0;}
  //SHOULD BE SIMPLIFIED- test with superblocks
  double Discharge(0.0);
	double sum(0),end1(0),end2(0);
  if (ltype==LINESINK){
		sum=-1.0;
    for (int n=0; n<=order; n++){
			sum*=-1.0;
			end1+=JumpCoeff[n];  
			end2+=JumpCoeff[n]*sum;//sum=(-1)^n
		}
		Discharge=end2-end1;
		/*for (int n=0; n<=order; n++){
			sum*=-1.0; 
			Discharge+=JumpCoeff[n]*(sum-1.0);//sum=(-1)^n
		}*/
	}
	return Discharge;
}
/***********************************************************************
 Get Segment Logarithm Potential Function
	Returns the potential Omega=Phi + i Psi from the endpoint singularities of a linesink
----------------------------------------------------------------------*/
cmplex  CLineElem::GetLog(const cmplex &z, const double &t) const{
	if (disabled){return 0.0;}
	cmplex omega(0.0);

	if (ltype==LINESINK){
		cmplex  logend1(0.0),logend2(0.0);
	  cmplex  Z((z-0.5*(z1+z2))/(0.5*(z2-z1)));
    double  end1(0.0),end2(0.0);
		int     n;
	  double sum;

    cmplex  ZBH((pLayer->GetBlackHole()-0.5*(z1+z2))/(0.5*(z2-z1)));
    logend1=log((Z-1.0)/abs(ZBH-1.0));
    logend2=log((Z+1.0)/abs(ZBH+1.0));  

		sum=-1.0;
    for (n=0; n<=order; n++){
			sum*=-1.0;
			end1+=JumpCoeff[n];  
			end2+=JumpCoeff[n]*sum;//sum=(-1)^n
		}

		omega+=IOVER2PI*(-(end1*logend1)+(end2*logend2));

		if ((CAnalyticElem::BranchcutLocus) && (ltype==LINESINK)){
			double argz=arg((z               -z1));
			double arg1=arg(-(zBranchcutLocus-z1));//if negative, away from BCL, else towards
			if      ((arg1>=0) && (argz>0.0) && (argz>arg1)) {omega-=IM*(end2-end1);}
			else if ((arg1< 0) && (argz<0.0) && (argz<arg1)) {omega+=IM*(end2-end1);}
		}
  }

	return omega;
}

/***********************************************************************
				Jump Functions
************************************************************************
 Returns the jump in potential DOmega=DPhi + i DPsi across the element
----------------------------------------------------------------------*/
cmplex CLineElem::GetOmJump(double X,double t) const{
	if (disabled){return 0.0;}
  cmplex jump(0.0);
	cmplex mult(IM);
	if ((ltype==LINESINK) || (ltype==DIPOLE)){mult=cmplex(1.0,0);}
  for(int n=0; n <=order; n++){jump+=mult*JumpCoeff[n]*cos(n*acos(X));} //TMP DEBUG- need to change
  return -IM*jump;
}
/* Returns the jump in complex discharge DW=DQX - iDQY across the element
-----------------------------------------------------------------------*/
cmplex CLineElem::GetWJump(double X,double t) const{
	if (disabled){return 0.0;}
  double sum, length(abs(z2-z1));
	cmplex jump(0.0);
  cmplex mult(IM);
	if ((ltype==LINESINK) || (ltype==DIPOLE)){mult=cmplex(1.0,0);}
  for(int n=0; n <=order; n++){  
    sum=0.0;
    for(int k=n-1; k>=0; k-=2){
			if (k!=0){sum+=2.0*n*cos(k*acos(X));}
			else     {sum+=    n*cos(k*acos(X));}
		} 
		jump+=mult*JumpCoeff[n]*2.0/length*sum;
  }
  return -IM*jump;  
}
/* Returns the Cumulative Extraction of the element
-----------------------------------------------------------------------*/
double CLineElem::GetCumExtract(const double X, const leftright direct, double t) const{
	if (disabled){return 0.0;}
  double cumext(0),mult(1.0),sum(-1.0);
	int    n;
	if (direct==DIR_RIGHT){mult=-1.0;}
  if (ltype==LINESINK){
		for (n=0; n<=order; n++){
			sum*=-1.0;
			cumext+=JumpCoeff[n]*(2.0*sum+cos(n*acos(mult*X))-1.0);
		}
		return cumext;
	}
	return 0.0;		
}
/***********************************************************************
				PURELY VIRTUAL MEMBER FUNCTIONS 
************************************************************************/
void CLineElem::GetMatrixBuildInfo (MatrixInfo &info){
	ExitGracefully("CLineElem::GetSegMatrixBuildInfo: Virtual function. Should never get called.",VIRTUAL_ERROR);}

/***********************************************************************
				MEMBER FUNCTIONS
***********************************************************************/
//-----------------------------------------------------------------------
void CLineElem::WriteItself(ofstream &SOL, const double &t) const{
	// * [] [] [] [] [] [] [] [] ... ([] = coeff from a_0 to a_N)
  int n;
	SOL << " * ";
	for (n=0; n<=order; n++){
		SOL<< JumpCoeff[n] <<" ";
	}
  SOL<<endl;
}
//-----------------------------------------------------------------------
bool CLineElem::ReadItself(ifstream &SOL){
	// * [] [] [] [] [] [] [] [] ... ([] = coeff from a_0 to a_N)
	int Len(0);
	char *s[MAXINPUTITEMS];
	bool warm(true);
	int n;

	if (!TokenizeLine(SOL,s,Len)) {
		if ((strcmp(s[0],"*"))){return false;}

		for (n=0; n<=min(order,Len-2); n++){
			JumpCoeff[n]=s_to_d(s[n+1]);          
		}
		for (n=min(order,Len-2)+1; n<=order; n++){
			JumpCoeff[n]=0.0;          
		}		

    SetFarFieldCoeff();
	}
	else {warm=false;}

	//cout << "warm= " <<warm<<endl;
  return warm;

}
/************************************************************************
*************************************************************************
				EXPLICIT SOLVER FUNCTIONS 
*************************************************************************
************************************************************************/
int  CLineElem::GetDegreesOfFreedom() const {return order+1;}
//-----------------------------------------------------------------------
void CLineElem::GetUnitInfluences  (const int n, const cmplex *pts,
																		const int NumCtrl, double *uPhi, cmplex *uQ, const double t){
  int m,n2;

	//set pJumpCoeff[n] to 1, others to zero
	for (n2=0; n2<=order; n2++){
		JumpCoeff[n2]=0.0;
	}
	JumpCoeff[n]=1.0;
	SetFarFieldCoeff();

	//obtain unit influences	
	for (m=0;m<NumCtrl;m++){
		uPhi[m]=GetDischargePotential(pts[m],t).real();
		uQ  [m]=GetW                 (pts[m],t);
	}
}
//-----------------------------------------------------------------------
void CLineElem::SetCoeff           (double *coeff){
	//copies coefficients
	for (int n=0; n<=order; n++){JumpCoeff[n]=coeff[n];}
	SetFarFieldCoeff();
}
/************************************************************************
*************************************************************************
				Geometry Functions
*************************************************************************
************************************************************************/

/************************************************************************
				Interpolate Control Points
************************************************************************/
void CLineElem::Interpolate(double *Array, int size, double v1, double v2){
  for (int m=0; m<size; m++){
    Array[m]=((v2-v1)*0.5*CONTROLEND*cos(((double)(m+1)-0.5)/(double)(size)*PI)+(v1+v2)/2.0);
	}
}
//------------------------------------------------------------------------
void CLineElem::Interpolate(cmplex *Array, int size, cmplex v1, cmplex v2){
  for (int m=0; m<size; m++){
    Array[m]=((v2-v1)*0.5*CONTROLEND*(double)(cos(((double)(m+1)-.5)/(double)(size)*PI))+(v1+v2)/2.0);
	}
}
//*********************************************************************** 
void CLineElem::WriteGeometry(ofstream &BASEMAP) const              {
	BASEMAP << "\" StringElem \",  "<<-2<<endl;	 
	BASEMAP <<z1.real()<< " , " <<z1.imag()<<endl;
	BASEMAP <<z2.real()<< " , " <<z2.imag()<<endl;
}
//------------------------------------------------------------------------
cmplex CLineElem::Centroid() const{
  return 0.5*(z1+z2); 
}
//------------------------------------------------------------------------
bool CLineElem::IsInSquare(const cmplex &zc,const double w) const{
	if (z1.real() > zc.real() + (w/2.0)) {return false;}
  if (z1.real() < zc.real() - (w/2.0)) {return false;}
  if (z1.imag() > zc.imag() + (w/2.0)) {return false;}
  if (z1.imag() < zc.imag() - (w/2.0)) {return false;}
	if (z2.real() > zc.real() + (w/2.0)) {return false;}
  if (z2.real() < zc.real() - (w/2.0)) {return false;}
  if (z2.imag() > zc.imag() + (w/2.0)) {return false;}
  if (z2.imag() < zc.imag() - (w/2.0)) {return false;}
  return true;
}
//------------------------------------------------------------------------
bool CLineElem::IsInCircle(const cmplex &zc,const double r) const{
	if (abs(z1-zc)>r){return false;}
	if (abs(z2-zc)>r){return false;}
	return true;
}
//------------------------------------------------------------------------
bool CLineElem::SharesNode(const cmplex &zn) const{
	if (zn==z1) {return true;}
	if (zn==z2) {return true;}
	return false;
}
//------------------------------------------------------------------------
bool CLineElem::PartInCircle   (const cmplex &zc, const double r) const{
	double a,b,c;
  //solve quadratic
	a = pow(z2.real() - z1.real(),2)+ pow(z2.imag() - z1.imag(),2) ;

  b = 2.0*((z2.real() - z1.real())*(z1.real() - zc.real()) + 
					 (z2.imag() - z1.imag())*(z1.imag() - zc.imag())); 
      
	c = pow(     zc.real(),2) + pow(     zc.imag(),2) + 
			pow(z1.real(),2) + pow(z1.imag(),2) - 
			2.0*(zc.real()* z1.real() + 
					 zc.imag()* z1.imag()) - pow(r,2);

	if ((abs(z1-zc)<=r) || 
			(abs(z2-zc)<=r) ||
			(b*b-4.0*a*c>=0.0)){
    return true;
	}
	else {
		return false;
	}
} 

	/*
	Old GetSegment Leakage Function
		if ((LJumpCoeff!=NULL) && (LFarFldCoeff!=NULL)){
  
		double  logend1,logend2,leak(0.0),end1(0.0),end2(0.0);
	  cmplex  Z((z-0.5*(z1+z2))/(0.5*(z2-z1)));
		cmplex  jump(0.0);
		double  modulus(abs(Z));

		if (Leakltype==LINESINK){
			cmplex  ZBH((pLayer->GetBlackHole()-0.5*(z1+z2))/(0.5*(z2-z1)));
			logend1=log((Z-1.0)/abs(ZBH-1.0)).real();
			logend2=log((Z+1.0)/abs(ZBH+1.0)).real();
			for (int n=0; n<=order[i]; n++){
				end1+=LJumpCoeff[i][n].real();
				end2+=LJumpCoeff[i][n].real()*pow(-1.0,n);
			}
		}
		//far-field (without truncation)
		if (modulus>ISFAR){ 
			leak+=Outside(Z,(int)(FForder[i]*pow(modulus/ISFAR,FF_POWER)),LFarFldCoeff[i]).real();
			if (Leakltype==LINESINK){
				leak+=IOVER2PI*(-(end1*logend1)+(end2*logend2));
			}
		}  
		//near-field
		else if ((abs(Z-E1)>MOVE_DIST*abs(z1-z2)) && (abs(Z-E2)>MOVE_DIST*abs(z1-z2))){
			cmplex logterm(log(Zm1oZp1(Z)));
			for (int n=0; n<=order; n++){     
				jump+=conj(LJumpCoeff[i][n])*cClenshaw2(n,Z,f[n],logterm);
			}
      leak+=IOVER2PI*jump.real();
		  if (Leakltype==LINESINK) {
				leak+=IOVER2PI *(-(end1*logend1)+(end2*logend2));
			}
  
		}
		//singularities
		else {leak=GetSegmentLeakage(i,z+cmplex(MOVE_DIST*abs(z1-z2),0),c,t);}

		return leak/c;
	}
	else {return 0.0;}*/

/*
//***********************************************************************
//				Get Segment Complex Potential due to Leakage  
//***********************************************************************
cmplex CLineElem::GetSegmentLeakOm(const int i,const cmplex &z,const double &t) const{
	double  logend1,logend2,omega(0.0),modulus;
	double  end1(0.0),end2(0.0),end3(0.0),end4(0.0),end5(0.0),end6(0.0);
	cmplex  Z,ZBH,jump(0.0),jump2(0.0),jump3(0.0);
  cmplex *fZj;          
  int     n,j;
  double L(abs(z1-z2));
	if ((LJumpCoeff!=NULL) && (LFarFldCoeff!=NULL)){
    
		fZj=new cmplex [order+1];

		Z=  (           z-0.5*(z1+z2))/(0.5*(z2-z1));
		ZBH=(pLayer->GetBlackHole()-0.5*(z1+z2))/(0.5*(z2-z1));
		modulus=abs(Z);
		logend1=log((Z-1.0)/abs(ZBH-1.0)).real();
		logend2=log((Z+1.0)/abs(ZBH+1.0)).real();


		for (n=0; n<=order; n++){
			end1+=h(LJumpCoeff[i],n).real();
			end2+=h(LJumpCoeff[i],n).real()*pow(-1.0,n);
			end3+=h(LJumpCoeff[i],n).imag();
      end4+=h(LJumpCoeff[i],n).imag()*pow(-1.0,n);
			end5+=hh(LJumpCoeff[i],n).imag();
      end6+=hh(LJumpCoeff[i],n).imag()*pow(-1.0,n);  
		}

		//far-field (without truncation)
		if (modulus>ISFAR){ 
			//int nterms=(int)(min(FForder[i],
			//                 MIN_FAR_FIELD+(FForder[i]-MIN_FAR_FIELD)*(ISFAR/modulus)));
			//leak+=Outside(Z,FForder[i],LFarFldCoeff[i]).real();
			//if (Leakltype==LINESINK){leak+=IOVER2PI*(-(end1*logend1)+(end2*logend2));}
		}  
		//near-field
		else if ((abs(Z-E1)>MOVE_DIST*abs(z1-z2)) && (abs(Z-E2)>MOVE_DIST*abs(z1-z2))){
			for (n=0; n<=order; n++){
				for (j=0; j<n; j++){fZj[j]=f[n][j];}
        fZj[n]=log(Zm1oZp1(Z));      
        jump +=conj( h(LJumpCoeff[i],n)).imag()*cClenshaw(n,Z,fZj);   //doublet term (without log)
        jump2+=conj( h(LJumpCoeff[i],n)).real()*cClenshaw(n,Z,fZj);   //dipole term 1 (without logs)
				jump3+=conj(hh(LJumpCoeff[i],n)).real()*cClenshaw(n,Z,fZj);   //dipole term 2 (without logs)
			}
			//doublet
			omega=(L*(Z-conj(Z))*IOVER2PI/16.0*(jump-(end1*logend1)+(end2*logend2))).real(); 
      //dipole 
			omega =((Z-conj(Z))*L*IOVER2PI/16.0*(jump2-(end3*logend1)+(end4*logend2))).real();
			omega+=(            L*IOVER2PI/8.0 *(jump2-(end5*logend1)+(end6*logend2))).real();
      omega+=(           -L*IOVER2PI/8.0 *(end3*((Z-1.0)*logend1-(Z-1.0)))).real();
			omega+=(            L*IOVER2PI/8.0 *(end4*((Z+1.0)*logend2-(Z+1.0)))).real();

		  //add multipoles for linesinks 
			//add taylor series
		}
		//singularities
		//else {leak=GetSegmentLeakage(i,z+cmplex(MOVE_DIST*abs(z1-z2),0),c,t);}

		delete [] fZj;


		return omega;// /c;
	}
	else {return 0.0;}
}*/
//***********************************************************************
//				Get Segment Complex Discharge due to Leakage
//***********************************************************************
/*cmplex CLineElem::GetSegmentLeakW(const int i,const cmplex &z,const double &t) const{
	return 0.0;
}
*/
//-----------------------------------------------------------------------
/*double CLineElem::GetHeadJump(int i,double X,double t) const{
  double jump(0.0);
  for(int n=0; n <=order; n++){jump+=LJumpCoeff[i][n].imag()*cos(n*acos(X));}
  return jump;
}
//-----------------------------------------------------------------------
double CLineElem::GetHeadGradJump(int i,double X,double t) const{
  double sum, length(abs(z2-z1));
	cmplex jump(0);

  for(int n=0; n <=order; n++){  
    sum=0.0;
    for(int k=n-1; k>=0; k-=2){
			if (k!=0){sum+=2.0*n*cos(k*acos(X));}
			else     {sum+=    n*cos(k*acos(X));}
		} 
    jump+=conj(LJumpCoeff[i][n])*2.0/length*sum;
	}
  return jump.real();
}*/