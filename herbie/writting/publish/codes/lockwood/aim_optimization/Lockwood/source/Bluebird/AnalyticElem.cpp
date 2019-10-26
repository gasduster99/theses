#include "AnalyticElem.h"

/************************************************************************
                           CAnalyticElem
													 CONSTRUCTORS
************************************************************************/
CAnalyticElem::CAnalyticElem(){nSubElems=0;}
//------------------------------------------------------------------------
CAnalyticElem::CAnalyticElem(char *Name, const CSingleLayerABC *pLay){
	if (Name!=NULL){
    name = new char[strlen(Name)+1];
    name=strcpy(name,Name);
	}
	nSubElems=0; 
	pSubElemArray=NULL; 
	pBlock=NULL; 	
	pLayer=pLay;	
	given=false;
	elemID=TotalElems;       
	
	//TotalElems++;  
	//pAllElements[TotalElems-1]=this;			//add to global array		

	if (!DynArrayAppend((void**&)(pAllElements),(void*)(this),TotalElems)){
		 ExitGracefully("CAnalyticElem::Constructor: creating NULL element",BAD_DATA);};
}
//------------------------------------------------------------------------
CAnalyticElem::CAnalyticElem(CAnalyticElem **p, const CSingleLayerABC *p2, int s){ 
		
	pLayer=p2; 
	nSubElems=s; 
	pSubElemArray=NULL;
	pBlock=NULL; 
	given=false;
  for (int i=0;i<nSubElems; i++){AddToContainer(p[i]);}

  elemID=TotalElems; 
	
	//TotalElems++; 
	//pAllElements[TotalElems-1]=this;			//add to global array		

	if (!DynArrayAppend((void**&)(pAllElements),(void*)(this),TotalElems)){
		 ExitGracefully("CAnalyticElem::Constructor: creating NULL element",BAD_DATA);};

}
//------------------------------------------------------------------------
CAnalyticElem::~CAnalyticElem(){ 
	if (globaldebug){cout <<"     DESTROYING ANALYTICELEM "<<name<<endl;}
	delete [] name; 
	if (nSubElems>0){delete [] pSubElemArray;}
}
/************************************************************************
                           STATIC INITIALIZATION 
************************************************************************/
int             CAnalyticElem::TotalElems      =0;
CAnalyticElem **CAnalyticElem::pAllElements    =NULL; 
int             CAnalyticElem::DefaultPrecision=3; 
bool            CAnalyticElem::Cosmetic        =false;
bool            CAnalyticElem::BranchcutLocus  =false;
cmplex          CAnalyticElem::zBranchcutLocus =0.0;
//------------------------------------------------------------------------
void						CAnalyticElem::DestroyAllElements   (){
	if (globaldebug){cout <<"DESTROYING ALL ELEMENTS"<<endl;}
	for (int i=0; i<TotalElems; i++){
		delete pAllElements[i];
	}
	delete [] pAllElements;
}
//------------------------------------------------------------------------
void            CAnalyticElem::SetDefaultPrecision(const int prec){
	CAnalyticElem::DefaultPrecision=prec;
}
//------------------------------------------------------------------------
void            CAnalyticElem::SetCosmetic         (const bool on){Cosmetic=on;}
//------------------------------------------------------------------------
void            CAnalyticElem::SetBranchCutLocus    (const bool on, const cmplex &z){
	BranchcutLocus=on; 
	zBranchcutLocus=z;
} 
/************************************************************************
                           ACCESSOR FUNCTIONS
************************************************************************/
int             CAnalyticElem::GetID()         const {return elemID;}
//------------------------------------------------------------------------
int             CAnalyticElem::GetNumSubElems()const {return nSubElems;}
//------------------------------------------------------------------------
char*           CAnalyticElem::GetName()       const {return name;}
//------------------------------------------------------------------------
bool						CAnalyticElem::IsGiven()       const {return given;}
//------------------------------------------------------------------------
bool            CAnalyticElem::HasFlux()       const {return false;}
//------------------------------------------------------------------------
CAnalyticElem*  CAnalyticElem::GetAllSubElems()const {return *pSubElemArray;} //must edit for recursive search of containers
//------------------------------------------------------------------------
int             CAnalyticElem::GetTotalElements()    {return TotalElems;}
/************************************************************************
                           ASSIGNMENT FUNCTIONS
************************************************************************/
void CAnalyticElem::SetBlockOwner(COwnerABC *BlockPtr, int seg, int IDinBlock){
	pBlock=BlockPtr;
	myBlockID=IDinBlock;
}
//------------------------------------------------------------------------
void CAnalyticElem::AddToContainer(CAnalyticElem *Elemptr){
  //elemType must be changed if type !=(aggregate)
	if (!DynArrayAppend((void**&)(pSubElemArray),(void*)(Elemptr),nSubElems)){
	  ExitGracefully("CAnalyticElem::AddToContainer: adding NULL element",BAD_DATA);};
} 
void CAnalyticElem::UpdateBlock   (const double &t) const{
	/*if ((pBlock!=NULL) && (pBlock->IsOn)){
			pBlock->Update(myBlockID,NOT_A_STRING,t);
	}*/
}
/************************************************************************
                           GetDischargePotential
************************************************************************/
cmplex CAnalyticElem::GetDischargePotential(const cmplex &z,const double &t) const{
  cmplex Omega(0,0);
  for (int i=0; i<nSubElems; i++){Omega+=pSubElemArray[i]->GetDischargePotential(z,t);}
  return Omega;
}
/************************************************************************
                           GetW
************************************************************************/
cmplex CAnalyticElem::GetW(const cmplex &z,const double &t) const{
  cmplex W(0,0);
  for (int i=0; i<nSubElems; i++){W+=pSubElemArray[i]->GetW(z,t);}
  return W;
}
/************************************************************************
                           GetGx
************************************************************************/
cmplex CAnalyticElem::GetGx(const cmplex &z, const double &t) const{
  cmplex Gx(0,0);
  for (int i=0; i<nSubElems; i++){Gx+=pSubElemArray[i]->GetGx(z,t);}
  return Gx;
}
/************************************************************************
                          GetLeakage
************************************************************************/
double CAnalyticElem::GetLeakage(const cmplex &z, const double &t,const leak_type ltype) const{
  double Leak(0.0);
  for (int i=0; i<nSubElems; i++){Leak+=pSubElemArray[i]->GetLeakage(z,t,ltype);}
  return Leak;
}
/************************************************************************
                          GetLeakage
************************************************************************/
double CAnalyticElem::GetCurl(const cmplex &z, const double &t) const{
  double Curl(0.0);
  for (int i=0; i<nSubElems; i++){Curl+=pSubElemArray[i]->GetCurl(z,t);}
  return Curl;
}
/************************************************************************
                           GetNetDischarge
************************************************************************/
double CAnalyticElem::GetNetDischarge(const double &t) const{
  double TD(0.0);
  for (int i=0; i<nSubElems; i++){TD+=pSubElemArray[i]->GetNetDischarge(t);}
  return TD;
}
/************************************************************************
                           GetFluxThruFace
Gets integrated normal and tangential flux between points z1 and z2 
returns Qflux=Qn + i Qt, where Qn is integrated normal flux through face and Qt is integrated tangential flux
Requires BranchCutLocus To be enabled
************************************************************************/
cmplex CAnalyticElem::GetFluxThruFace   (const cmplex &z1, const cmplex &z2, const double &t) const{
  cmplex Qflux(0.0);
  for (int i=0; i<nSubElems; i++){Qflux+=pSubElemArray[i]->GetFluxThruFace(z1,z2,t);}
  return Qflux;
}
/************************************************************************
                           GetIntegratedLeakage
Gets integrated leakage [L^3/T] into triangle defined by z1,z2,z3
************************************************************************/
double CAnalyticElem::GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, const leak_type ltype) const{
  double Q(0.0);
  for (int i=0; i<nSubElems; i++){Q+=pSubElemArray[i]->GetIntegratedLeakage(z1,z2,z3,t,ltype);}
  return Q;
}
/************************************************************************
                           GetIntegratedBudget
Gets influx and outflux [L^3/T] into triangle defined by z1,z2,z3
************************************************************************/
void CAnalyticElem::GetIntegratedBudget(const cmplex &z1,const cmplex &z2,const cmplex &z3, 
																				const double &t,double &inQ, double &outQ) const{
  double tmpin(0.0),tmpout(0.0);
	inQ=outQ=0.0;
  for (int i=0; i<nSubElems; i++){
		pSubElemArray[i]->GetIntegratedBudget(z1,z2,z3,t,tmpin,tmpout);
		inQ +=tmpin;
		outQ+=tmpout;
	}
}
/************************************************************************
                           Get Flux Distribution
	Gets integrated extraction between points z1 and z2 along line
	Linearly Distributes these fluxes to endpoints Q1 and Q2
	Used for finite element source distribution
************************************************************************/
void CAnalyticElem::GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const{
  double tmpQ1(0.0),tmpQ2(0.0);
	Q1=Q2=0.0;
  for (int i=0; i<nSubElems; i++){
		pSubElemArray[i]->GetFluxDistribution(z1,z2,t,tmpQ1,tmpQ2);
		Q1+=tmpQ1;
		Q2+=tmpQ2;
	}
}
/************************************************************************
                           Get Cumulative Flux
	Returns cumulative extraction ("base flow") along surface water network
	if this is an extraction element (e.g., linesink) this will return a value 
	if the point z is within NEAR_FEATURE of the element, otherwise, returns zero
************************************************************************/
double CAnalyticElem::GetCumulativeFlux    (const cmplex &z, const double &t) const{
  double tmpFlux(0.0);
  for (int i=0; i<nSubElems; i++){
		tmpFlux=max(pSubElemArray[i]->GetCumulativeFlux(z,t),tmpFlux);
	}
	return tmpFlux;
}
/************************************************************************
                           Get Cumulative Flux
	Returns complex strength of singularity at point z (e.g., Q/2pi for a well)
	if the point z is within NEAR_FEATURE of the element, otherwise, returns zero
************************************************************************/
cmplex CAnalyticElem::GetSingularStrength  (const cmplex &z, const double &t) const{
  cmplex S(0.0);
  for (int i=0; i<nSubElems; i++){
		S+=pSubElemArray[i]->GetSingularStrength(z,t);
	}
	return S;
}
/************************************************************************
                           SolveItself
************************************************************************/
void CAnalyticElem::SolveItself(double &change, double &objective,const double t) {

  double tempchange,maxchange(0.0),tempobj,maxobj(0.0);

  for (int i=0; i<nSubElems; i++){
    tempchange=tempobj=0.0;
    pSubElemArray[i]->SolveItself(tempchange,tempobj,t);
		upperswap(maxchange,tempchange);
		upperswap(maxobj   ,tempobj   );
  } 
  change   =maxchange;  
	objective=maxobj   ;
}
/************************************************************************
                           WriteItself
************************************************************************/
void	 CAnalyticElem::WriteOutput(const double &t)								const					{}
bool   CAnalyticElem::ReadItself(ifstream &SOL)                                 {return true;}
double CAnalyticElem::GetMaxError(const double &t)								const					{return 0.0;}
/************************************************************************
              EXPLICIT SOLVER FUNCTIONS (purely virtual) 
************************************************************************/
int  CAnalyticElem::GetDegreesOfFreedom   () const  {return 0;}
//------------------------------------------------------------------------
void CAnalyticElem::GetMatrixBuildInfo    (MatrixInfo &info){
	ExitGracefully("CAnalyticElem::GetMatrixBuildInfo:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void CAnalyticElem::GetUnitInfluences     (const int n,const cmplex *pts,const int NumCtrl,
																					 double *uPhi, cmplex *uQ, const double t){
	ExitGracefully("CAnalyticElem::GetUnitInfluences:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void CAnalyticElem::SetCoeff              (double *coeff){
	ExitGracefully("CAnalyticElem::SetCoeff:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void  CAnalyticElem::GetFluxMatrixBuildInfo (MatrixInfo &info){
	ExitGracefully("CAnalyticElem::GetFluxMatrixBuildInfo:this should never be called!",VIRTUAL_ERROR); 
}
//------------------------------------------------------------------------
void CAnalyticElem::GetFluxUnitInfluences  (const cmplex *pts,const int NumCtrl,double *uPhi, const double t){
	ExitGracefully("CAnalyticElem::GetFluxUnitInfluences:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void CAnalyticElem::SetFluxCoeff           (double coeff){
	ExitGracefully("CAnalyticElem::SetFluxCoeff:this should never be called!",VIRTUAL_ERROR);
}
/************************************************************************
                           GEOMETRY FUNCTIONS (all blank)
************************************************************************/
cmplex         CAnalyticElem::Centroid() const                                    {return 0;}
bool           CAnalyticElem::IsInside(const cmplex &z) const                     {return false;}
bool           CAnalyticElem::IsInSquare(const cmplex &zc,const double w) const   {return false;}
bool           CAnalyticElem::IsInCircle(const cmplex &zc,const double r) const   {return false;}
bool           CAnalyticElem::PartInCircle(const cmplex &zc,const double r) const {return false;}
bool           CAnalyticElem::SharesNode(const cmplex &zn) const                  {return false;} 
void           CAnalyticElem::WriteGeometry(ofstream &BASEMAP) const              {}
