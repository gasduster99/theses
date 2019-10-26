#include "Circle.h"
/***********************************************************************
************************************************************************
                           CCircularElem
***********************************************************************/
/***********************************************************************
							 CONSTRUCTOR
***********************************************************************/
CCircularElem::CCircularElem(){
	ncircontrol=0;order=0;OutCoeff=NULL;InCoeff=NULL;zctrl=NULL;
}
//-----------------------------------------------------------------------
CCircularElem::CCircularElem(char						 *Name, 
														 const CSingleLayerABC *pLay, 
														 const cmplex			zc, 
														 const double			rad, 
														 const int				prec):
							 CAnalyticElem(Name,pLay){
	int    n,m;
  double angle, fold;

	zcen       =zc;        
	R          =rad;         
	SetPrecision(prec,order,fold);   
	ncircontrol=(int)(fold*order);

  //set up and initialize dynamic arrays
  OutCoeff= new cmplex [order];
	InCoeff = new cmplex [order];
  zctrl=     new cmplex [ncircontrol];  

	for (n=0; n<order; n++){OutCoeff[n]=cmplex(0.0,0.0);InCoeff[n]=cmplex(0.0,0.0);}
  for (m=0; m<ncircontrol; m++){
    angle=2.0*PI*(double)(m)/(double)(ncircontrol);
    zctrl[m]= R*cmplex(cos(angle), sin(angle))+zcen; 
	}

}
//-----------------------------------------------------------------------
CCircularElem::~CCircularElem(){
  if (globaldebug){cout <<"    DESTROYING CIRELEM "<<name<<endl;}
  delete [] OutCoeff;
	delete [] InCoeff;
  delete [] zctrl;;  
}
/***********************************************************************
							 STATIC MEMBERS
***********************************************************************/
void CCircularElem::SetPrecision(const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("SetPrecision::Improper precision level specified",BAD_DATA);}
	switch(Precision){
		case(0): {order=3;  fold=3.0; break;}
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
cmplex CCircularElem::GetZcen()   const {return zcen;}
//-----------------------------------------------------------------------
double CCircularElem::GetRadius() const {return R;}
//-----------------------------------------------------------------------
bool   CCircularElem::HasFlux()   const {return false;}

/***********************************************************************
							 Read & Write
***********************************************************************/
void   CCircularElem::WriteItself(ofstream &SOL, const double &t) const{
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
bool CCircularElem::ReadItself(ifstream &SOL){

	// if old solution has different number of coefficients, then truncation is used

	int Len(0);
	char *s[MAXINPUTITEMS];
	bool warm(true);
	int n;

	if  (!TokenizeLine(SOL,s,Len)) {} //skip name of element
	if ((!TokenizeLine(SOL,s,Len)) && (Len==1)) {
		Q=s_to_d(s[0]);
	}
	else {return false;}

	if (!TokenizeLine(SOL,s,Len)) {
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
void CCircularElem::UpdateBlock(const double &t) const{
	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,NOT_A_STRING,t);
	}
}

/***********************************************************************
							 GEOMETRY MEMBER FUNCTIONS
***********************************************************************/
void CCircularElem::WriteGeometry(ofstream &BASEMAP) const              {
	BASEMAP << "\" CircleElem \",  2" <<endl;
	BASEMAP << zcen.real()<< " , "    <<zcen.imag()<<endl;
	BASEMAP << R          << " , 0.0"              <<endl;
}
//-----------------------------------------------------------------------
cmplex CCircularElem::Centroid() const {return zcen;}
//-----------------------------------------------------------------------
bool   CCircularElem::IsInside(const cmplex &z) const{
  if (abs(z-zcen)<=R)                  {return true; } 
	else                                 {return false;}
} 
//-----------------------------------------------------------------------
bool   CCircularElem::IsInSquare(const cmplex &zc,const double w) const{
  if ((zcen.real()-R > zc.real()-w) &&
      (zcen.real()+R < zc.real()+w) &&
      (zcen.imag()-R > zc.imag()-w) &&
      (zcen.imag()+R < zc.imag()+w)) {return true; }
	else                               {return false;}
}
//-----------------------------------------------------------------------
bool   CCircularElem::IsInCircle(const cmplex &zc,const double r) const{
	if   (abs(zcen-zc)<(r-R))            {return false;} 
	else                                 {return true; }
}
//-----------------------------------------------------------------------
bool   CCircularElem::PartInCircle(const cmplex &zc,const double r) const{
	if (abs(zcen-zc)<(r+R)) {return false;}
	else                    {return true;}
}
//-----------------------------------------------------------------------
double CCircularElem::GetArea() const     {return PI*R*R;}
//------------------------------------------------------------------------
bool CCircularElem::SharesNode(const cmplex &zn) const{
	return false;
}
/***********************************************************************
							 GetDischargePotential
----------------------------------------------------------------------*/
cmplex CCircularElem::GetDischargePotential (const cmplex &z,const double &t) const{
  cmplex Omega(0.0);
  cmplex Z((z-zcen)/R);
  if (abs(Z)>=1.0) {Omega=Laurent(Z,order,OutCoeff)+IOVER2PI*Q*log((z-zcen)/(pLayer->GetBlackHole()-zcen));}
  else             {Omega=Taylor (Z,order,InCoeff);}
  return Omega;
}

/***********************************************************************
							 Get W
----------------------------------------------------------------------*/
cmplex CCircularElem::GetW(const cmplex &z,const double &t) const{
  cmplex W(0.0);
  cmplex Z((z-zcen)/R);
  if (abs(Z)>1.0){ W=LaurentDer(Z,order,OutCoeff,R)-IOVER2PI*Q/(z-zcen);}
  else           { W=TaylorDer (Z,order,InCoeff, R);}
  return W;
}
/***********************************************************************
							 Get Gx
----------------------------------------------------------------------*/
cmplex CCircularElem::GetGx(const cmplex &z,const double &t) const{
  cmplex Gx(0.0);
  cmplex Z((z-zcen)/R);
  if (abs(Z)>1.0){ Gx=-(LaurentDxx(Z,order,OutCoeff,R))+(IOVER2PI*Q/(z-zcen)/(z-zcen));} 
  else           { Gx=-(TaylorDxx (Z,order,InCoeff ,R));}
  return Gx;
}
/***********************************************************************
							 Get Flux Through Face
----------------------------------------------------------------------*/
cmplex CCircularElem::GetFluxThruFace(const cmplex &z1, const cmplex &z2, const double &t) const{
	
	intercepttype tmp;
	cmplex zint1,zint2;
	tmp=LineIntersectCircle(z1,z2,zcen,R,zint1,zint2);
	if      (tmp==NO_INTERSECT){
		return GetDischargePotential(z1,t).imag()-GetDischargePotential(z2,t).imag();
		//return Im*conj(GetDischargePotential(z1,t)-GetDischargePotential(z2,t)); //complex flux thru face
	}
	else if (tmp==INTERSECTED){ //one intersection
		return GetDischargePotential(z1,t).imag()-GetDischargePotential(zint1-(z2-z1)*REALSMALL,t).imag()-
           GetDischargePotential(z2,t).imag()+GetDischargePotential(zint1+(z2-z1)*REALSMALL,t).imag();
		/*return IM*conj(GetDischargePotential(z1,t)-GetDischargePotential(zint1-(z2-z1)*REALSMALL,t)-
                   GetDischargePotential(z2,t)+GetDischargePotential(zint1+(z2-z1)*REALSMALL,t));*/
	}
	else if (tmp==TWO_INTERSECTIONS){
		if (abs(zint1-z1)<abs(zint2-z1)){
		  return GetDischargePotential(z1                     ,t).imag()-GetDischargePotential(zint1-(z2-z1)*REALSMALL,t).imag()+
             GetDischargePotential(zint1+(z2-z1)*REALSMALL,t).imag()-GetDischargePotential(zint2-(z2-z1)*REALSMALL,t).imag()+				
		         GetDischargePotential(zint2+(z2-z1)*REALSMALL,t).imag()-GetDischargePotential(z2                     ,t).imag();
			/*return IM*conj(GetDischargePotential(z1                     ,t)-GetDischargePotential(zint1-(z2-z1)*REALSMALL,t)+
					           GetDischargePotential(zint1+(z2-z1)*REALSMALL,t)-GetDischargePotential(zint2-(z2-z1)*REALSMALL,t)+				
			               GetDischargePotential(zint2+(z2-z1)*REALSMALL,t)-GetDischargePotential(z2                     ,t));*/ //complex flux thru face
		}
		else {
		  return GetDischargePotential(z1                     ,t).imag()-GetDischargePotential(zint2-(z2-z1)*REALSMALL,t).imag()+
             GetDischargePotential(zint2+(z2-z1)*REALSMALL,t).imag()-GetDischargePotential(zint1-(z2-z1)*REALSMALL,t).imag()+				
			       GetDischargePotential(zint1+(z2-z1)*REALSMALL,t).imag()-GetDischargePotential(z2                     ,t).imag();			
			/*return IM*conj(GetDischargePotential(z1                     ,t)-GetDischargePotential(zint2-(z2-z1)*REALSMALL,t)+
                     GetDischargePotential(zint2+(z2-z1)*REALSMALL,t)-GetDischargePotential(zint1-(z2-z1)*REALSMALL,t)+				
			               GetDischargePotential(zint1+(z2-z1)*REALSMALL,t)-GetDischargePotential(z2                     ,t));*/
		}
	}
	else {return 0.0;} //never happens
}
/***********************************************************************
							 Get Leakage
----------------------------------------------------------------------*/
double  CCircularElem::GetLeakage (const cmplex &z, const double &t,const leak_type ltype) const{
	return 0.0; 
}
/***********************************************************************
							 Get Net Discharge
----------------------------------------------------------------------*/
double  CCircularElem::GetNetDischarge(const double &t) const{
	return Q; 
}


