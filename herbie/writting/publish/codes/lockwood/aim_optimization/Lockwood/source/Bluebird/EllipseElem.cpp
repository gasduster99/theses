//EllipseElem.cpp
#include "EllipseElem.h"
//************************************************************************
//************************************************************************
//                           CEllipseElem
//************************************************************************
//************************************************************************
//							 CONSTRUCTOR
//***********************************************************************
CEllipseElem::CEllipseElem(){
	nellcontrol=0;order=0;OutCoeff=NULL;InCoeff=NULL;zctrl=NULL;
}
//-----------------------------------------------------------------------
CEllipseElem::CEllipseElem(char			       *Name,      
													 const CSingleLayerABC *pLay,       
													 const cmplex		  zcen,         
													 const double     MajorAxis, 
													 const double		  MinorAxis,   
													 const double     theta,		
													 const int			  prec):
							CAnalyticElem(Name,pLay){
	int    n,m;
	double fold;

	zc=zcen;                      
	a=MajorAxis;
	b=MinorAxis;
	c=sqrt(a*a-b*b);
	angle=theta;

	if (b>a){ExitGracefully("CEllipseElem::Constructor: bad axes specified",BAD_DATA);}

	SetPrecision(prec,order,fold);   
	nellcontrol=(int)(fold*order);

  //set up and initialize dynamic arrays
  OutCoeff= new cmplex [order];
	InCoeff = new cmplex [order];
  zctrl=    new cmplex [nellcontrol];  

	eta0=GlobalToLocal(zc+a*exp(IM*angle)).real();

	for (n=0; n<order; n++){OutCoeff[n]=cmplex(0.0,0.0);InCoeff[n]=cmplex(0.0,0.0);}
	double psi;
  for (m=0; m<nellcontrol; m++){
    psi=-PI+2.0*PI*(double)(m)/(double)(nellcontrol);
    zctrl[m]=LocalToGlobal(cmplex(eta0,psi)); 
	}


}
//-----------------------------------------------------------------------
CEllipseElem::~CEllipseElem(){
  if (globaldebug){cout <<"    DESTROYING ELLIPSEELEM "<<name<<endl;}
  delete [] OutCoeff;
	delete [] InCoeff;
  delete [] zctrl;;  
}
/***********************************************************************
							 STATIC MEMBERS
***********************************************************************/
void CEllipseElem::SetPrecision(const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("Ellipse:SetPrecision::Improper precision level specified",BAD_DATA);}
	switch(Precision){
		case(0): {order=1;  fold=1.0; break;}
		case(1): {order=10; fold=1.5; break;}
		case(2): {order=20; fold=1.8; break;}
		case(3): {order=30; fold=2.0; break;}
	  case(4): {order=40; fold=2.0; break;}
	  case(5): {order=75; fold=2.0; break;}
	  case(9): {order=100;fold=3.0; break;}
	}
}
/***********************************************************************
							 ACCESSORS
***********************************************************************/
cmplex CEllipseElem::GetCenter()const {return zc;}
//-----------------------------------------------------------------------
void   CEllipseElem::GetAxes(double &Major, double &Minor) const {Major=a; Minor=b;}
//-----------------------------------------------------------------------
bool   CEllipseElem::HasFlux() const {return false;}
/***********************************************************************
							 Read & Write
***********************************************************************/
void   CEllipseElem::WriteItself(ofstream &SOL, const double &t) const{
  int n;
	SOL << name <<endl;
	SOL << Q    <<endl;
	SOL << " * ";
	for (n=0; n<order; n++){
		SOL<< OutCoeff[n].real() <<" "<<OutCoeff[n].imag() <<" ";
	}
	SOL <<endl<< " * ";
	for (n=0; n<order; n++){
		SOL<< InCoeff[n].real() <<" "<<InCoeff[n].imag() <<" ";
	}

  SOL<<endl;
}
//-----------------------------------------------------------------------
bool CEllipseElem::ReadItself(ifstream &SOL){
	int Len(0);
	char *s[MAXINPUTITEMS];
	bool warm(true);
	int n;

	if  (!TokenizeLine(SOL,s,Len)) {} //skip name of element
	if ((!TokenizeLine(SOL,s,Len)) && (Len==1)) {
		Q=s_to_d(s[0]);
	}
	else {warm=false;}

	if ((warm) && (!TokenizeLine(SOL,s,Len))) {
		if ((strcmp(s[0],"*"))){return false;}

		for (n=0; n<min(order,(Len-1)/2); n++){
			OutCoeff[n]=s_to_c(s[2*n+1],s[2*n+2]);          
		}
		for (n=min(order,(Len-1)/2)+1; n<order; n++){
			OutCoeff[n]=0.0;          
		}
	}
	else {return false;}//eof
	
	if ((warm) && (!TokenizeLine(SOL,s,Len))) {
		if ((strcmp(s[0],"*"))){return false;}

		for (n=0; n<min(order,(Len-1)/2); n++){
			InCoeff[n]=s_to_c(s[2*n+1],s[2*n+2]);          
		}
		for (n=min(order,(Len-1)/2)+1; n<order; n++){
			InCoeff[n]=0.0;          
		}
	}
	else {return false;}//eof
	return warm;
}
//-----------------------------------------------------------------------
void CEllipseElem::UpdateBlock(const double &t) const{
	if ((pBlock!=NULL) && (pBlock->IsOn())){
			pBlock->Update(myBlockID,NOT_A_STRING,t);
	}
}
/***********************************************************************
							 GlobalToLocal
***********************************************************************/
cmplex CEllipseElem::GlobalToLocal(const cmplex &z) const{
	cmplex Z=(z-zc)/exp(IM*angle);
	cmplex zeta;
	
	if (Z.real()>0.0){zeta=log(Z/c+sqrt(Z*Z/(c*c)-1.0));}
	else             {zeta=log(Z/c-sqrt(Z*Z/(c*c)-1.0));}

	return zeta;
}
cmplex CEllipseElem::LocalToGlobal(const cmplex &zeta) const{
	//return zc+Z;
	
	return zc+cmplex(c*cosh(zeta.real())*cos(zeta.imag()),c*sinh(zeta.real())*sin(zeta.imag()))*exp(IM*angle);
}

/***********************************************************************
							 GEOMETRY MEMBER FUNCTIONS
***********************************************************************/
void CEllipseElem::WriteGeometry(ofstream &BASEMAP) const              {
	//This will not draw ellipses with orientation!
	BASEMAP << "\" CircleElem \",  2" <<endl;
	BASEMAP << zc.real()<< " , "    << zc.imag()<<endl;
	BASEMAP << a                 << " , "    << b            <<endl;
}
//-----------------------------------------------------------------------
cmplex CEllipseElem::Centroid() const  {return zc;}
//-----------------------------------------------------------------------
bool   CEllipseElem::IsInside(const cmplex &z) const{
  if (GlobalToLocal(z).real()<=eta0)   {return true; } 
	else                                 {return false;}
} 
//-----------------------------------------------------------------------
bool   CEllipseElem::IsInSquare(const cmplex &zc,const double w) const{
	for (int m=0; m<nellcontrol; m++){
		if ((zctrl[m].real() < zc.real()-w) ||
				(zctrl[m].real() > zc.real()+w) ||
			  (zctrl[m].imag() < zc.imag()-w) ||
			  (zctrl[m].imag() > zc.imag()+w)) {return false;}
	}
	return true;
}
//-----------------------------------------------------------------------
bool   CEllipseElem::IsInCircle(const cmplex &zc,const double r) const{
	for (int m=0; m<nellcontrol; m++){
	  if (abs(zctrl[m]-zc)>r) {return false;} 
	}
	return true; 
}
//-----------------------------------------------------------------------
bool   CEllipseElem::PartInCircle(const cmplex &zc,const double r) const{
	for (int m=0; m<nellcontrol; m++){
	  if (abs(zctrl[m]-zc)<r) {return true;} 
	}
	return false;
}
//-----------------------------------------------------------------------
double CEllipseElem::GetArea() const     {return PI*a*b;}
//------------------------------------------------------------------------
bool CEllipseElem::SharesNode(const cmplex &zn) const{return false;}
//***********************************************************************
//							 GetDischargePotential
//***********************************************************************
cmplex CEllipseElem::GetDischargePotential (const cmplex &z,const double &t) const{
  cmplex Omega(0.0);
  cmplex zeta(GlobalToLocal(z));
	int n;
  if (zeta.real()<eta0){ //Inside ellipse
		Omega=2.0*InCoeff[0];

		for(n=1; n<order; n++){
	    Omega+=conj(InCoeff[n])*(exp((double)(n)*zeta)+exp(-(double)(n)*zeta));
		}
	}
	else{                 //Outside Ellipse
		for(n=1;n<order;n++){      
			if ((double)(n)*zeta.real()<MAXEXP) {
				Omega+=OutCoeff[n]*exp(-(double)(n)*zeta);
			}
		}
		Omega += Q*zeta - Q*GlobalToLocal(pLayer->GetBlackHole());
	}

  return Omega;

}
//***********************************************************************
//							 GetW
//***********************************************************************
cmplex CEllipseElem::GetW(const cmplex &z,const double &t) const{
  cmplex W(0.0);
  cmplex zeta(GlobalToLocal(z));
	cmplex sinhzeta;
	int n;
	sinhzeta=cmplex(sinh(zeta.real())*cos(zeta.imag()),cosh(zeta.real())*sin(zeta.imag()));

  if (zeta.real()<eta0){ //Inside ellipse
		for(n=1; n<order; n++){
			W -=conj(InCoeff[n])*(exp((double)(n)*zeta)-exp(-(double)(n)*zeta))*(double)(n)/(sinhzeta*c);
		}
	}
	else{                  //Outside Ellipse
		for(n=1;n<order;n++){      
			if ((double)(n)*zeta.real()<MAXEXP) {
				W+=OutCoeff[n]*exp(-(double)(n)*zeta)*(double)(n)/(sinhzeta*c);
			}
		}
	}
	
	W -=Q/(sinhzeta*c);
  W*=exp(-IM*(angle));

  return W;
}
//***********************************************************************
//							 GetGx
//***********************************************************************
cmplex CEllipseElem::GetGx(const cmplex &z,const double &t) const{
	//TMP DEBUG-must have Raghu Handle this
  cmplex Gx(0.0);
  cmplex zeta(GlobalToLocal(z));
  if (zeta.real()<eta0){ //Inside ellipse
	
	}
	else{                  //Outside Ellipse

	}

  return Gx;
}
/***********************************************************************
							 Get Flux Through Face
----------------------------------------------------------------------*/
cmplex CEllipseElem::GetFluxThruFace(const cmplex &z1, const cmplex &z2, const double &t) const{
	
	//if (Q!=0.0) && IntersectEllipse(z1,z2)... //TMP DEBUG- should hangle intersections
  
  return IM*conj(GetDischargePotential(z1,t)-GetDischargePotential(z2,t)); //complex flux thru face
}
//***********************************************************************
//							 GetLeakage
//***********************************************************************
double  CEllipseElem::GetLeakage (const cmplex &z, const double &t,const leak_type ltype) const{
	return 0.0; 
}
//***********************************************************************
//							 GetNetDischarge
//***********************************************************************
double  CEllipseElem::GetNetDischarge(const double &t) const{
	return Q; 
}
