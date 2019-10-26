#include "MLAnalyticElem.h"

/************************************************************************
                           CmlAnalyticElem
													 CONSTRUCTORS
************************************************************************/
CmlAnalyticElem::CmlAnalyticElem(char *Name, const CMultiLayerABC *pLay){
	if (Name!=NULL){
    name = new char[strlen(Name)+1];
    name=strcpy(name,Name);
	}
	nSubElems=0; 
	pSubElemArray=NULL; 
	pLayer=pLay;	
	given=false;
	elemID=TotalElems;       
	
	if (!DynArrayAppend((void**&)(pAllElements),(void*)(this),TotalElems)){
		 ExitGracefully("CmlAnalyticElem::Constructor: creating NULL element",BAD_DATA);};
}
//------------------------------------------------------------------------
CmlAnalyticElem::CmlAnalyticElem(CmlAnalyticElem **p, const CMultiLayerABC *p2, int s){ 
		
	pLayer=p2; 
	nSubElems=s; 
	pSubElemArray=NULL;
	given=false;
  for (int i=0;i<nSubElems; i++){AddToContainer(p[i]);}

  elemID=TotalElems; 
	
	if (!DynArrayAppend((void**&)(pAllElements),(void*)(this),TotalElems)){
		 ExitGracefully("CmlAnalyticElem::Constructor: creating NULL element",BAD_DATA);};

}
//------------------------------------------------------------------------
CmlAnalyticElem::~CmlAnalyticElem(){ 
	if (globaldebug){cout <<"     DESTROYING MULTILAYER ANALYTICELEM "<<name<<endl;}
	delete [] name; 
	if (nSubElems>0){delete [] pSubElemArray;}
}
/************************************************************************
                           STATIC INITIALIZATION 
************************************************************************/
int             CmlAnalyticElem::TotalElems      =0;
CmlAnalyticElem **CmlAnalyticElem::pAllElements    =NULL; 
int             CmlAnalyticElem::DefaultPrecision=3; 
//------------------------------------------------------------------------
void						CmlAnalyticElem::DestroyAllMLElements   (){
	if (globaldebug){cout <<"DESTROYING ALL ELEMENTS"<<endl;}
	for (int i=0; i<TotalElems; i++){
		delete pAllElements[i];
	}
	delete [] pAllElements;
}
//------------------------------------------------------------------------
void            CmlAnalyticElem::SetDefaultPrecision(const int prec){
	CmlAnalyticElem::DefaultPrecision=prec;
}
/************************************************************************
                           ACCESSOR FUNCTIONS
************************************************************************/
int             CmlAnalyticElem::GetID()         const {return elemID;}
//------------------------------------------------------------------------
int             CmlAnalyticElem::GetNumSubElems()const {return nSubElems;}
//------------------------------------------------------------------------
char*           CmlAnalyticElem::GetName()       const {return name;}
//------------------------------------------------------------------------
bool						CmlAnalyticElem::IsGiven()       const {return given;}
//------------------------------------------------------------------------
bool            CmlAnalyticElem::HasFlux()       const {return false;}
//------------------------------------------------------------------------
CmlAnalyticElem*  CmlAnalyticElem::GetAllSubElems()const {return *pSubElemArray;} //must edit for recursive search of containers
//------------------------------------------------------------------------
int             CmlAnalyticElem::GetTotalElements()    {return TotalElems;}
/************************************************************************
                           ASSIGNMENT FUNCTIONS
************************************************************************/
/*void CmlAnalyticElem::SetBlockOwner(COwnerABC *BlockPtr, int seg, int IDinBlock){
	pBlock=BlockPtr;
	myBlockID=IDinBlock;
}*/
//------------------------------------------------------------------------
void CmlAnalyticElem::AddToContainer(CmlAnalyticElem *Elemptr){
  //elemType must be changed if type !=(aggregate)
	if (!DynArrayAppend((void**&)(pSubElemArray),(void*)(Elemptr),nSubElems)){
	  ExitGracefully("CmlAnalyticElem::AddToContainer: adding NULL element",BAD_DATA);};
} 
//void CmlAnalyticElem::UpdateBlock   (const double &t) const{
	/*if ((pBlock!=NULL) && (pBlock->IsOn)){
			pBlock->Update(myBlockID,NOT_A_STRING,t);
	}*/
//}
/************************************************************************
                           Get Comprehensive Potential
************************************************************************/
double CmlAnalyticElem::GetCompPotential(const cmplex &z,const double &t) const{
  double pot;int i;
  pot=0.0;
  for (i=0; i<nSubElems; i++){
		pot+=pSubElemArray[i]->GetCompPotential(z,t);
	}
  return pot;
}
/************************************************************************
                           Get Helmholtz potential vector
************************************************************************/
Ironclad1DArray	CmlAnalyticElem::GetHelmPotentialVect (const cmplex &z, const double &t) const{
  static double hpot[MAX_MLAYERS];
	Unchangeable1DArray tmphpot;
	int i,l;

  for (l=0; l<pLayer->GetNumLayers(); l++){hpot[l]=0.0;}

  for (i=0; i<nSubElems; i++){
		tmphpot=pSubElemArray[i]->GetHelmPotentialVect(z,t);
		for (l=0; l<pLayer->GetNumLayers(); l++){
			hpot[l]+=tmphpot[l];
		}
	}
  return hpot;
}
/************************************************************************
                           Get Discharge Potential Vector
************************************************************************/
Ironclad1DArray CmlAnalyticElem::GetPotentialVect(const cmplex &z,const double &t) const{
  static double pot[MAX_MLAYERS];
	Unchangeable1DArray tmppot;
	int i,l;

  for (l=0; l<pLayer->GetNumLayers(); l++){pot[l]=0.0;}

  for (i=0; i<nSubElems; i++){
		tmppot=pSubElemArray[i]->GetPotentialVect(z,t);
		for (l=0; l<pLayer->GetNumLayers(); l++){
			pot[l]+=tmppot[l];
		}
	}
  return pot;
}
/************************************************************************
                           GetW
************************************************************************/
Ironclad1DArray_z CmlAnalyticElem::GetQxQyVect(const cmplex &z,const double &t) const{
  static cmplex W[MAX_MLAYERS];
	Unchangeable1DArray_z tmpW;
	int i,l;

  for (l=0; l<pLayer->GetNumLayers(); l++){W[l]=0.0;}

  for (i=0; i<nSubElems; i++){
		tmpW=pSubElemArray[i]->GetQxQyVect(z,t);
		for (l=0; l<pLayer->GetNumLayers(); l++){
			W[l]+=tmpW[l];
		}
	}
  return W;
}
/************************************************************************
                           GetGx
************************************************************************/
Ironclad1DArray_z CmlAnalyticElem::GetGxVect(const cmplex &z, const double &t) const{
  static cmplex Gx[MAX_MLAYERS];
	Unchangeable1DArray_z tmpGx;
	int i,l;

  for (l=0; l<pLayer->GetNumLayers(); l++){Gx[l]=0.0;}

  for (i=0; i<nSubElems; i++){
		tmpGx=pSubElemArray[i]->GetGxVect(z,t);
		for (l=0; l<pLayer->GetNumLayers(); l++){
			Gx[l]+=tmpGx[l];
		}
	}
  return Gx;
}
/************************************************************************
                          GetLeakage
************************************************************************/
Ironclad1DArray CmlAnalyticElem::GetLeakageVect(const cmplex &z, const double &t,const leak_type ltype) const{
  static double Leak[MAX_MLAYERS];
	Unchangeable1DArray tmpleak;
	int i,l;

  for (l=0; l<pLayer->GetNumLayers(); l++){Leak[l]=0.0;}

  for (i=0; i<nSubElems; i++){
		tmpleak=pSubElemArray[i]->GetLeakageVect(z,t,ltype);
		for (l=0; l<pLayer->GetNumLayers(); l++){
			Leak[l]+=tmpleak[l];
		}
	}
  return Leak;
}
/************************************************************************
                          Get Curl Vector
************************************************************************/
Ironclad1DArray CmlAnalyticElem::GetCurlVect(const cmplex &z, const double &t) const{
  static double Curl[MAX_MLAYERS];
	Unchangeable1DArray tmpcurl;
	int i,l;

  for (l=0; l<pLayer->GetNumLayers(); l++){Curl[l]=0.0;}

  for (i=0; i<nSubElems; i++){
		tmpcurl=pSubElemArray[i]->GetCurlVect(z,t);
		for (l=0; l<pLayer->GetNumLayers(); l++){
			Curl[l]+=tmpcurl[l];
		}
	}
  return Curl;
}
/************************************************************************
                           GetNetDischarge
************************************************************************/
double CmlAnalyticElem::GetNetDischarge(const int lev, const double &t) const{
  double TD(0.0);
  for (int i=0; i<nSubElems; i++){TD+=pSubElemArray[i]->GetNetDischarge(lev,t);}
  return TD;
}
//-----------------------------------------------------------------------
double CmlAnalyticElem::GetNetDischarge(const double &t) const{
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
cmplex CmlAnalyticElem::GetFluxThruFace   (const cmplex &z1, const cmplex &z2, const int lev, const double &t) const{
  cmplex Qflux(0.0);
  for (int i=0; i<nSubElems; i++){Qflux+=pSubElemArray[i]->GetFluxThruFace(z1,z2,lev,t);}
  return Qflux;
}
/************************************************************************
                           GetIntegratedLeakage
Gets integrated leakage [L^3/T] into triangle defined by z1,z2,z3
************************************************************************/
double CmlAnalyticElem::GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, const cmplex &z3, const int lev, const double &t, const leak_type ltype) const{
  double Q(0.0);
  for (int i=0; i<nSubElems; i++){Q+=pSubElemArray[i]->GetIntegratedLeakage(z1,z2,z3,lev,t,ltype);}
  return Q;
}
/************************************************************************
                           GetIntegratedBudget
Gets influx and outflux [L^3/T] into triangle defined by z1,z2,z3
************************************************************************/
void CmlAnalyticElem::GetIntegratedBudget(const cmplex &z1,const cmplex &z2,const cmplex &z3,  const int lev,
																				const double &t,double &inQ, double &outQ) const{
  double tmpin(0.0),tmpout(0.0);
	inQ=outQ=0.0;
  for (int i=0; i<nSubElems; i++){
		pSubElemArray[i]->GetIntegratedBudget(z1,z2,z3,lev,t,tmpin,tmpout);
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
void CmlAnalyticElem::GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const int lev, const double &t, double &Q1, double &Q2) const{
  double tmpQ1(0.0),tmpQ2(0.0);
	Q1=Q2=0.0;
  for (int i=0; i<nSubElems; i++){
		pSubElemArray[i]->GetFluxDistribution(z1,z2,lev,t,tmpQ1,tmpQ2);
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
double CmlAnalyticElem::GetCumulativeFlux    (const cmplex &z, const double &t) const{
  double tmpFlux(0.0);
  for (int i=0; i<nSubElems; i++){
		tmpFlux=max(pSubElemArray[i]->GetCumulativeFlux(z,t),tmpFlux);
	}
	return tmpFlux;
}
/************************************************************************
                           SolveItself
************************************************************************/
void CmlAnalyticElem::SolveItself(double &change, double &objective,const double t) {

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
void	 CmlAnalyticElem::WriteOutput(const double &t)								const					{}
bool   CmlAnalyticElem::ReadItself(ifstream &SOL)                                 {return true;}
double CmlAnalyticElem::GetMaxError(const double &t)								const					{return 0.0;}
/************************************************************************
              EXPLICIT SOLVER FUNCTIONS (purely virtual) 
************************************************************************/
/*int  CmlAnalyticElem::GetDegreesOfFreedom   () const  {return 0;}
//------------------------------------------------------------------------
void CmlAnalyticElem::GetMatrixBuildInfo    (MatrixInfo &info){
	ExitGracefully("CmlAnalyticElem::GetMatrixBuildInfo:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void CmlAnalyticElem::GetUnitInfluences     (const int n,const cmplex *pts,const int NumCtrl,
																					 double *uPhi, cmplex *uQ, const double t){
	ExitGracefully("CmlAnalyticElem::GetUnitInfluences:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void CmlAnalyticElem::SetCoeff              (double *coeff){
	ExitGracefully("CmlAnalyticElem::SetCoeff:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void  CmlAnalyticElem::GetFluxMatrixBuildInfo (MatrixInfo &info){
	ExitGracefully("CmlAnalyticElem::GetFluxMatrixBuildInfo:this should never be called!",VIRTUAL_ERROR); 
}
//------------------------------------------------------------------------
void CmlAnalyticElem::GetFluxUnitInfluences  (const cmplex *pts,const int NumCtrl,double *uPhi, const double t){
	ExitGracefully("CmlAnalyticElem::GetFluxUnitInfluences:this should never be called!",VIRTUAL_ERROR);
}
//------------------------------------------------------------------------
void CmlAnalyticElem::SetFluxCoeff           (double coeff){
	ExitGracefully("CmlAnalyticElem::SetFluxCoeff:this should never be called!",VIRTUAL_ERROR);
}*/
/************************************************************************
                           GEOMETRY FUNCTIONS (all blank)
************************************************************************/
cmplex         CmlAnalyticElem::Centroid() const                                    {return 0;}
bool           CmlAnalyticElem::IsInside(const cmplex &z) const                     {return false;}
bool           CmlAnalyticElem::IsInSquare(const cmplex &zc,const double w) const   {return false;}
bool           CmlAnalyticElem::IsInCircle(const cmplex &zc,const double r) const   {return false;}
bool           CmlAnalyticElem::PartInCircle(const cmplex &zc,const double r) const {return false;}
bool           CmlAnalyticElem::SharesNode(const cmplex &zn) const                  {return false;} 
void           CmlAnalyticElem::WriteGeometry(ofstream &BASEMAP) const              {}
