#include "MLFarField.h"
/**********************************************************************
                             FarField
***********************************************************************
														CONSTRUCTORS
***********************************************************************/
CmlFarField::CmlFarField(CMultiLayerABC *pLay, 
												 cmplex						UnifFlow, 
												 cmplex						Zref, 
												 double						refhead):
						 CmlAnalyticElem("FarField",pLay){

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
CmlFarField::~CmlFarField(){
	if (globaldebug){cout <<"   DESTROYING FARFIELD"<<endl;}
}
/***********************************************************************
					ACCESSORS
***********************************************************************/
cmplex CmlFarField::GetZref() const             {return zref;}
//-----------------------------------------------------------------------
cmplex CmlFarField::GetUnifFlow() const         {return Qo;}
//-----------------------------------------------------------------------
bool   CmlFarField::IsReferencePt() const       {return RefPoint;}

/***********************************************************************
					MANIPULATORS
***********************************************************************/
void   CmlFarField::SetZref(cmplex z)           {zref=z;}

/***********************************************************************
					Get Discharge Potential
----------------------------------------------------------------------*/
double CmlFarField::GetCompPotential(const cmplex &z, const double &t) const{  
	static double pot;
	//TMP DEBUG : STUB
  return pot;
}
//----------------------------------------------------------------------
Ironclad1DArray CmlFarField::GetHelmPotentialVect(const cmplex &z, const double &t) const{  
	static double hpot[MAX_MLAYERS];
	//TMP DEBUG : STUB
  return hpot;
}
//----------------------------------------------------------------------
Ironclad1DArray CmlFarField::GetPotentialVect(const cmplex &z, const double &t) const{  
	static double pot[MAX_MLAYERS];
	//TMP DEBUG : STUB
  return pot;
}
/***********************************************************************
					Get W
----------------------------------------------------------------------*/
Ironclad1DArray_z CmlFarField::GetQxQyVect(const cmplex &z,const double &t) const{
	static cmplex W[MAX_MLAYERS];
	//TMP DEBUG : STUB
  return W;
}
/***********************************************************************
					Get Gx
----------------------------------------------------------------------*/
Ironclad1DArray_z CmlFarField::GetGxVect(const cmplex &z,const double &t) const{
	static cmplex Gx[MAX_MLAYERS];
	for (int l=0; l<pLayer->GetNumLayers(); l++){Gx[l]=0.0;}
  return Gx;
}
/***********************************************************************
					GetFlux Through Face
----------------------------------------------------------------------*/
cmplex CmlFarField::GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const int lev, const double &t) const{
	
	//return GetDischargePotential(z1,t).imag()-GetDischargePotential(z2,t).imag();
  return 0.0;//IM*conj(GetDischargePotential(z1,t)-GetDischargePotential(z2,t)); //complex flux thru face
}
/***********************************************************************
					Get Leakage
----------------------------------------------------------------------*/
Ironclad1DArray CmlFarField::GetLeakageVect(const cmplex &z, const double &t,const leak_type ltype) const{
	static double leak[MAX_MLAYERS];
	for (int l=0; l<pLayer->GetNumLayers(); l++){leak[l]=0.0;}
  return leak;
}

/***********************************************************************
					Get Net Discharge
----------------------------------------------------------------------*/
double CmlFarField::GetNetDischarge(const int    lev,const double &t) const{return 0.0;}
double CmlFarField::GetNetDischarge(const double &t) const{return 0.0;}

/***********************************************************************
					SolveItself
***********************************************************************/
void CmlFarField::SolveItself(double &change,double &objective,const double t){
  //double DesiredPot, RefPot,oldConst,T,K,B;
	//double relax(1.0);

  /*K=pLayer->GetCond (zref);
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
*/
}
//----------------------------------------------------------------------
double CmlFarField::GetMaxError(const double &t) const{

	return (pLayer->GetHeadVect(zref,t)[0]-RefHead);//absolute head error
	//double B=pLayer->GetBase (zref);
	//	return (pLayer->GetHead(zref,t)-RefHead)/max(RefHead-B,1.0);/relative head error
}
/***********************************************************************
					OTHERS:
***********************************************************************/
void CmlFarField::WriteItself(ofstream &SOL, const double &t) const{
	SOL << "FarField"<<endl;
	SOL << Constant << endl;
}
//------------------------------------------------------------------------
bool CmlFarField::ReadItself(ifstream &SOL){
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
void   CmlFarField::WriteOutput (const double &t)                 const{}
cmplex CmlFarField::Centroid    () const                               {return cmplex(0,0);}
bool   CmlFarField::IsInSquare  (const cmplex &zc,const double w) const{return false;}
bool   CmlFarField::IsInCircle  (const cmplex &zc,const double r) const{return false;}
bool   CmlFarField::PartInCircle(const cmplex &zc,const double r) const{return false;}
double CmlFarField::GetArea     () const                               {return ALMOST_INF;}
bool   CmlFarField::SharesNode  (const cmplex &zn) const               {return false;}