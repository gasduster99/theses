#include "FarField.h"
/**********************************************************************
                             FarField
***********************************************************************
														CONSTRUCTORS
***********************************************************************/
CFarField::CFarField():
           CAnalyticElem("FarField",NULL){

  zref=0;         
	RefHead=0;       
	Constant=0; 
	RefPoint=true; 
	alpha=0;   
	Qo=0;  
	for (int j=0;j<10;j++){OldC[j]=OldQ[j]=0;}
	given=true;
}
//----------------------------------------------------------------------
CFarField::CFarField(CSingleLayerABC *pLay, 
										 cmplex					  UnifFlow, 
										 double						NetExt):
           CAnalyticElem("FarField",pLay){
       
  Qo=UnifFlow;    
	RefHead=0;      
	Constant=0;      
	RefPoint=false;

	zref=cmplex(BHDIST,0); 
  if (Qo.real()!=0) {alpha=atan((Qo.imag())/(Qo.real()));}
  else              {alpha=PI/2.0;                       }
	for (int j=0;j<10;j++){OldC[j]=OldQ[j]=0;}
		given=true;
}
//----------------------------------------------------------------------
CFarField::CFarField(CSingleLayerABC *pLay, 
										 cmplex						UnifFlow, 
										 cmplex						Zref, 
										 double						refhead):
           CAnalyticElem("FarField",pLay){

  Qo=UnifFlow;    
	zref=Zref;      
	RefHead=refhead;   
	RefPoint=true; 
	Constant=0;

  if (Qo.real()!=0) {alpha=atan((Qo.imag())/(Qo.real()));}
  else              {alpha=PI/2.0;                       } 
	for (int j=0;j<10;j++){OldC[j]=OldQ[j]=0;}
	given=true;
}
//----------------------------------------------------------------------
CFarField::~CFarField(){
	if (globaldebug){cout <<"   DESTROYING FARFIELD"<<endl;}
}
/***********************************************************************
					ACCESSORS
***********************************************************************/
cmplex CFarField::GetZref() const             {return zref;}
//-----------------------------------------------------------------------
cmplex CFarField::GetUnifFlow() const         {return Qo;}
//-----------------------------------------------------------------------
bool   CFarField::IsReferencePt() const       {return RefPoint;}

/***********************************************************************
					MANIPULATORS
***********************************************************************/
void   CFarField::SetZref(cmplex z)           {zref=z;}

/***********************************************************************
					Get Discharge Potential
----------------------------------------------------------------------*/
cmplex CFarField::GetDischargePotential(const cmplex &z, const double &t) const{  
  return -conj(Qo)*(z-zref)+Constant;
}
/***********************************************************************
					Get W
----------------------------------------------------------------------*/
cmplex CFarField::GetW(const cmplex &z,const double &t) const{
	return conj(Qo);
}
/***********************************************************************
					Get Gx
----------------------------------------------------------------------*/
cmplex CFarField::GetGx(const cmplex &z,const double &t) const{
	return 0.0;
}
/***********************************************************************
					GetFlux Through Face
----------------------------------------------------------------------*/
cmplex CFarField::GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const double &t) const{
	
	//return GetDischargePotential(z1,t).imag()-GetDischargePotential(z2,t).imag();
  return IM*conj(GetDischargePotential(z1,t)-GetDischargePotential(z2,t)); //complex flux thru face
}
/***********************************************************************
					Get Leakage
----------------------------------------------------------------------*/
double CFarField::GetLeakage(const cmplex &z, const double &t,const leak_type ltype) const{
	return 0.0;
}

/***********************************************************************
					Get Net Discharge
----------------------------------------------------------------------*/
double          CFarField::GetNetDischarge(const double &t) const{return 0;}

/***********************************************************************
					SolveItself
***********************************************************************/
void CFarField::SolveItself(double &change,double &objective,const double t){
  double DesiredPot, RefPot,oldConst,T,K,B;
	double relax(1.0);

  if (RefPoint) {
    K=pLayer->GetCond (zref);
		B=pLayer->GetBase (zref);
    T=pLayer->GetThick(zref);
    RefPot=pLayer->GetDischargePotential(zref,t).real()-Constant;
		DesiredPot=ConvertToPotential(RefHead-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());

    oldConst=Constant;
		//if (CSingleLayer::fresh ){ relax=1.0; }else{relax=1.0;}
    Constant=oldConst+(DesiredPot-RefPot-oldConst)*relax;
		change=fabs(Constant-oldConst)/pLayer->GetDeltaPot();
  	objective=0.0;
		//cout <<RefHead<< " "<< B << " " << K << " " <<T<<endl;
		//cout <<"Global C ="<< Constant<< " Ref Pot " << RefPot << " Des. Pot " << DesiredPot <<endl;
  }
  else {
	  double NetQ=pLayer->GetNetDischarge(t);
		for (int j=8;j>=0;j--){
		  OldC[j+1]=OldC[j];
      OldQ[j+1]=OldQ[j];
		}
		OldC[0]=Constant;
		OldQ[0]=NetQ;
	  change=objective=0.0;
	} //net extraction
}
//----------------------------------------------------------------------
double CFarField::GetMaxError(const double &t) const{
	if (RefPoint){
		double B=pLayer->GetBase (zref);
		return (pLayer->GetHead(zref,t)-RefHead);//absolute head error
	//	return (pLayer->GetHead(zref,t)-RefHead)/max(RefHead-B,1.0);/relative head error
	}
	else {return (pLayer->GetNetDischarge(t)-NetExt)/max(NetExt,1.0);} 
}
/***********************************************************************
					EXPLICIT SOLVER:
***********************************************************************/
int  CFarField::GetDegreesOfFreedom   () const {return 1;}
//----------------------------------------------------------------------
void CFarField::GetMatrixBuildInfo    (MatrixInfo &info){
	
	if (!RefPoint){ExitGracefully("CFarField::GetMatrixBuildInfo-Cannot have explicit solver with net extraction condition",BAD_DATA);}
	cout << "CFarField: Sending Explicit Matrix info"<<endl;
	double K=pLayer->GetCond (zref);
  double B=pLayer->GetBase (zref);
  double T=pLayer->GetThick(zref);	  
	info.nctrl=1;
	info.zctrl[0]=zref;
	info.unit[0][0]=1.0;
	info.elemrhs[0]=ConvertToPotential(RefHead-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());
	info.phiCoeff=-1;
	info.QxCoeff=0.0;
	info.QyCoeff=0.0;
//Part must be sent as given
}
//----------------------------------------------------------------------
void CFarField::GetUnitInfluences     (const int n,const cmplex *pts,const int NumCtrl,double *uPhi, cmplex *uQ, const double t){
	
	double oldConst;  
	int    m;

	if (n>0)      {ExitGracefully("CFarField::GetUnitInfluences-n should be zero",BAD_DATA);}
	if (!RefPoint){ExitGracefully("CFarField::GetUnitInfluences-Cannot have explicit solver with net extraction condition",BAD_DATA);}

	oldConst=Constant;
	Constant=1.0;
	for (m=0;m<NumCtrl;m++){
		if (t>=0.0){//non-given (C)
			uPhi[m]=1.0;
			uQ  [m]=0.0;
		}
		else       {//Given (Unif. flow)
			uPhi[m]=GetDischargePotential(pts[m],t).real();
			uQ  [m]=GetW                 (pts[m],t);		
		}
	}
	Constant=oldConst;
	Constant=0.0;
}
//----------------------------------------------------------------------
void CFarField::SetCoeff              (double *coeff){Constant=coeff[0];}
/***********************************************************************
					OTHERS:
***********************************************************************/
void CFarField::WriteItself(ofstream &SOL, const double &t) const{
	SOL << "FarField"<<endl;
	SOL << Constant << endl;
}
//------------------------------------------------------------------------
bool CFarField::ReadItself(ifstream &SOL){
	int Len(0);
	char *s[MAXINPUTITEMS];

	TokenizeLine(SOL,s,Len); //skip name of element

	if ((!TokenizeLine(SOL,s,Len)) && (Len==1)) {
		Constant=s_to_d(s[0]);
	} 
	else {return false;}

	return true;
}
//------------------------------------------------------------------------
void   CFarField::WriteOutput (const double &t)                 const{}
cmplex CFarField::Centroid    () const                               {return cmplex(0,0);}
bool   CFarField::IsInSquare  (const cmplex &zc,const double w) const{return false;}
bool   CFarField::IsInCircle  (const cmplex &zc,const double r) const{return false;}
bool   CFarField::PartInCircle(const cmplex &zc,const double r) const{return false;}
double CFarField::GetArea     () const                               {return ALMOST_INF;}
bool   CFarField::SharesNode  (const cmplex &zn) const               {return false;}