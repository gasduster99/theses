
#include "Layer.h"

void BuildMatrix(CAnalyticElem **pElems, int *segindices, int numelems, CAnalyticElem *OtherElems, const double t);
/*************************************************************************
                             CSingleLayer
************************************************************************
                           CONSTRUCTORS
************************************************************************
************************************************************************/
CSingleLayer::CSingleLayer(){ 
  black_hole=0;  level=0; 
	pLayerAbove=pLayerBeneath=NULL; 
	pFarField      =NULL;
  pElemArray     =NULL; nElems=0;
	pHSElems       =NULL; nHSElems=0;
	pGivenElems    =NULL; nGivenElems=0;
	pCondZoneArray =NULL; nCondZones=0;
	pBaseZoneArray =NULL; nBaseZones=0;
	pThickZoneArray=NULL; nThickZones=0;
	pPoroZoneArray =NULL; nPoroZones=0;
	pTaylorCirArray=NULL; nTaylorCirs=0;
	pInterpGrid    =NULL; InterpolationOn=false;
  QxInterp       =NULL;
	QyInterp       =NULL;

	DeltaPot=1.0;
	pMasterBlock   =NULL;

  Extents.w= ALMOST_INF; Extents.s= ALMOST_INF;
  Extents.e=-ALMOST_INF; Extents.n=-ALMOST_INF;
}
//------------------------------------------------------------------------
CSingleLayer::~CSingleLayer(){ 
	if (globaldebug){cout <<" DESTROYING LAYER "<<level<<endl;}
	delete [] pElemArray;     //just pointers: elements & propZones deleted globally
	delete [] pHSElems;
	delete [] pGivenElems;
	delete [] pCondZoneArray;
	delete [] pBaseZoneArray;
	delete [] pThickZoneArray;
	delete [] pPoroZoneArray;
	delete pInterpGrid;
	delete QxInterp;
	delete QyInterp;
	for (int i=0; i<nTaylorCirs; i++){delete pTaylorCirArray[i];} delete [] pTaylorCirArray;
	if ((pMasterBlock!=NULL) && (pMasterBlock->IsOn())){delete pMasterBlock;}
}
/************************************************************************
                 STATIC INITIALIZERS (default values)
*************************************************************************
************************************************************************/
bool   CSingleLayer::fresh             =true;
bool   CSingleLayer::solved            =true;
int    CSingleLayer::iter              =0;
bool   CSingleLayer::directflux        =false;
bool   CSingleLayer::GridWhileSolve    =false;
bool   CSingleLayer::FFeverytime       =false;
bool   CSingleLayer::FullAnalyze       =false;
window CSingleLayer::GWSW;
int    CSingleLayer::GWSresolution     =50; 
int    CSingleLayer::MinLayerIterations=1; 
int    CSingleLayer::MaxLayerIterations=300;
double CSingleLayer::LayerTolerance    =0.00001;
int    CSingleLayer::WriteInterval     =1000000;
double CSingleLayer::saltwaterSG       =1.0;
double CSingleLayer::sealevel          =-ALMOST_INF;
cmplex CSingleLayer::black_hole        =500.0;
/************************************************************************
************************************************************************
                           ACCESSOR FUNCTIONS
************************************************************************
************************************************************************/

/************************************************************************
 GetBlackHole
	returns location where logarithm terms go to zero (black hole)
-----------------------------------------------------------------------*/
cmplex          CSingleLayer::GetBlackHole()   const    {return CSingleLayer::black_hole;}
//------------------------------------------------------------------------
CSuperblock    *CSingleLayer::GetMasterBlock() const    {return pMasterBlock;}
//------------------------------------------------------------------------
CAnalyticElem  *CSingleLayer::GetElem(int i)   const    {if ((i<nElems) && (i>=0)){return pElemArray[i];} else{return NULL;}}   
//------------------------------------------------------------------------
int             CSingleLayer::GetNumElems()    const    {return nElems;}   
//------------------------------------------------------------------------
window          CSingleLayer::GetExtents()     const    {return Extents;}
//------------------------------------------------------------------------
double          CSingleLayer::GetDeltaPot()    const    {return DeltaPot;}
//------------------------------------------------------------------------
int             CSingleLayer::GetCurrentIter() const    {return iter;}
//------------------------------------------------------------------------
double          CSingleLayer::GetSeaLevel()    const    {return sealevel;}
//------------------------------------------------------------------------
double          CSingleLayer::GetSaltwaterSG() const		{return saltwaterSG;}
/************************************************************************
*************************************************************************
                           ASSIGNMENT FUNCTIONS
*************************************************************************
************************************************************************/
void CSingleLayer::SetLevelAbove (CAquicludeABC *above)  {pLayerAbove=above;}
//------------------------------------------------------------------------
void CSingleLayer::SetLevelBelow (CAquicludeABC *below)  {pLayerBeneath=below;}
//------------------------------------------------------------------------
void CSingleLayer::SetMasterBlock(CSuperblock *Master)  {pMasterBlock=Master;}
//------------------------------------------------------------------------
void CSingleLayer::SetInterpolationGrid (CMesh *pGrid)  {pInterpGrid=pGrid;}
//------------------------------------------------------------------------
void CSingleLayer::SetCoastalInfo(double sea_elev, double brineSG){
	sealevel=sea_elev; 
	if ((brineSG<=1.0) || (brineSG>5.0)){
		ExitGracefully("CSingleLayer::SetCoastalInfo: Bad value for saltwater density",BAD_DATA);}
	saltwaterSG=brineSG;
}

/************************************************************************
                           Set Black Hole
*************************************************************************
	places black hole where unif. flow potential is zero and far from centroid of domain
	this is where log terms and 1/z should all vanish
------------------------------------------------------------------------*/
void CSingleLayer::SetBlackHole(){
  cmplex zc;
  double alpha;

	if (pFarField==NULL){ExitGracefully("CSingleLayer::SetBlackHole: FarField doesn't exist", BAD_DATA);}

	if (pFarField->IsReferencePt()){black_hole=pFarField->GetZref();}
	else{
		cmplex unifflow=pFarField->GetUnifFlow();
		if (unifflow.real()!=0.0){alpha=(atan((unifflow.imag())/(unifflow.real())))+(PI*0.5);}
		else                     {alpha=0.0;}
		zc=cmplex((Extents.e+Extents.w)/2.0,(Extents.n+Extents.s)/2.0); //centroid of domain
		black_hole=cmplex(BHDIST*cos(alpha)+zc.real(), BHDIST*sin(alpha)+zc.imag());
		pFarField->SetZref(black_hole);
	}
}
/************************************************************************
                           SET BACKGROUND VALUES
*************************************************************************
	sets the background values for base, thickness, conductivity,porosity and level
-----------------------------------------------------------------------*/
void CSingleLayer::SetBackgroundValues(double B, double T, double K, double n, int L){
	if (K<0)            {ExitGracefully("CSingleLayer::SetBackgroundValues: negative conductivity",BAD_DATA);}
	if ((n<0) || (n>=1)){ExitGracefully("CSingleLayer::SetBackgroundValues: bad porosity value   ",BAD_DATA);}
  base_elev=B; thickness=T; conductivity=K; porosity=n; level=L;
}

/************************************************************************
                           CHANGE EXTENTS
*************************************************************************
	increases the extents of the domain window ("Extents") to include the point specified
------------------------------------------------------------------------*/
void CSingleLayer::UpdateExtents(const cmplex z){
	upperswap(Extents.e,z.real()); lowerswap(Extents.w,z.real());
	upperswap(Extents.n,z.imag()); lowerswap(Extents.s,z.imag());
	if (Extents.e<z.real()){Extents.e=z.real();}
}
/*------------------------------------------------------------------------
	increases the extents of the domain window ("Extents") to include the circle specified
------------------------------------------------------------------------*/
void CSingleLayer::UpdateExtents(const cmplex z,const double r){
	UpdateExtents(z+r);	
	UpdateExtents(z-r);	
	UpdateExtents(z+IM*r);	
	UpdateExtents(z-IM*r);
}
/************************************************************************
                           SET SOLVE DATA 
*************************************************************************
	sets the static solution data. 
  information does not have to be entered all at once
  if only one piece of information is required, send 0 as the other values
  these values are not currently tested for realism!!! TMP DEBUG 
------------------------------------------------------------------------*/
void CSingleLayer::SetSolveData(const int miniter, const int maxiter, const double tolerance){
	if (miniter>0)   {CSingleLayer::MinLayerIterations=miniter;  }
	if (maxiter>0)   {CSingleLayer::MaxLayerIterations=maxiter;  }
	if (tolerance>0) {CSingleLayer::LayerTolerance=    tolerance;}
}
//------------------------------------------------------------------------
void CSingleLayer::SetFFEveryTime() {CSingleLayer::FFeverytime =true;}
//------------------------------------------------------------------------
void CSingleLayer::SetFullAnalysis(){CSingleLayer::FullAnalyze =true;}
//------------------------------------------------------------------------
void CSingleLayer::SetWriteInterval(const int interval){
	if (interval>0){CSingleLayer::WriteInterval =interval;}
	else           {CSingleLayer::WriteInterval =1000;    }
}
/************************************************************************
                           Set Grid While Solve
*************************************************************************
	specifies gridding window and resolution for gridding every iteration
  this procedure takes place after any element is solved
------------------------------------------------------------------------*/
void CSingleLayer::SetGridWhileSolve(const double n, const double e, const double w, const double s, const int res){
	cout <<n<<" "<<e<<" "<<w<<" "<<s<<" "<<res<<endl;
	if ((n>s) && (e>w) && (res>0)){ 
		CSingleLayer::GWSW.n=n;
		CSingleLayer::GWSW.e=e;
		CSingleLayer::GWSW.w=w;
		CSingleLayer::GWSW.s=s;
		CSingleLayer::GWSresolution=res;
		CSingleLayer::GridWhileSolve =true;
	}
}
/*************************************************************************
                           ADD TO LAYER
**************************************************************************
	The following functions add elements, far field elements, 
  taylor series and property zones to the layer
------------------------------------------------------------------------*/

//	Adds Element to Layer
//------------------------------------------------------------------------
void CSingleLayer::AddToLayer(CAnalyticElem *Elemptr, elementtype type){

	if (!DynArrayAppend((void**&)(pElemArray),(void*)(Elemptr),nElems)){
		 ExitGracefully("CSingleLayer::AddToLayer(Elem): adding NULL element",BAD_DATA);};
	
	if (type==GIVEN_ELEM){
		if (!DynArrayAppend((void**&)(pGivenElems),(void*)(Elemptr),nGivenElems)){
		  ExitGracefully("CSingleLayer::AddToLayer(Elem): adding NULL given element",BAD_DATA);};
	}
	else if (type==HEAD_SPECIFIED_ELEM){
		if (!DynArrayAppend((void**&)(pHSElems),(void*)(Elemptr),nHSElems)){
		  ExitGracefully("CSingleLayer::AddToLayer(Elem): adding NULL head specified element",BAD_DATA);};
	}
} 

//	Adds FarField Element to Layer
//------------------------------------------------------------------------
void CSingleLayer::AddToLayer(CFarField  *FFptr){
 if (FFptr==NULL){ExitGracefully("CSingleLayer::AddToLayer- adding NULL far field",BAD_DATA);}
 if (pFarField!=NULL){ExitGracefully("CSingleLayer::AddToLayer- Two far field elements added!",BAD_DATA);}
 pFarField=FFptr;
}

//	Adds Taylor Series Circle to layer
//------------------------------------------------------------------------
void CSingleLayer::AddToLayer(CTaylorCir *Cirptr){
	if (!DynArrayAppend((void**&)(pTaylorCirArray),(void*)(Cirptr),nTaylorCirs)){
	  ExitGracefully("CSingleLayer::AddToLayer(Elem): adding NULL element",BAD_DATA);};
}

//	Adds Property Zone to layer
//------------------------------------------------------------------------
void CSingleLayer::AddToLayer(CPropZone *Propptr){
  if (Propptr==NULL){
		ExitGracefully("CSingleLayer::AddToLayer- adding NULL property zone",BAD_DATA);
	}
  if (Propptr->GetType()==hydraulic_conductivity){
		if (!DynArrayAppend((void**&)pCondZoneArray,(void*)(Propptr),nCondZones)){
			ExitGracefully("CSingleLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else if (Propptr->GetType()==base_elevation){
		if (!DynArrayAppend((void**&)pBaseZoneArray,(void*)(Propptr),nBaseZones)){
			ExitGracefully("CSingleLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else if (Propptr->GetType()==poro){
		if (!DynArrayAppend((void**&)pPoroZoneArray,(void*)(Propptr),nPoroZones)){
			ExitGracefully("CSingleLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else if (Propptr->GetType()==layer_thickness){
		if (!DynArrayAppend((void**&)pThickZoneArray,(void*)(Propptr),nThickZones)){
			ExitGracefully("CSingleLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else {ExitGracefully("unused prop type in layer",OTHER);}
}
/************************************************************************
                           GetCond
*************************************************************************
 Returns the hydraulic conductivity, K, at a given location
-----------------------------------------------------------------------*/
double CSingleLayer::GetCond(const cmplex &z) const {	
	static double Kloc;
  if (CSuperblock::SuperblocksOn()) {
		Kloc=pMasterBlock->GetCond(z);
		if (Kloc!=NO_VALUE){return Kloc;} 
		else               {return conductivity;}
	}
  else{
		return CPropZone::NestedSift(z,pCondZoneArray,nCondZones,conductivity);
	}
}
/************************************************************************
                           GetAnisotropy
*************************************************************************
 Returns the hydraulic conductivity anisotropy at a given location
-----------------------------------------------------------------------*/
anisotropy CSingleLayer::GetAnisotropy(const cmplex &z) const {	
	static anisotropy A;
	A.dir=0;
	A.ratio=0.0;
	return A;
}
/************************************************************************
                           GetBase
*************************************************************************
 Returns the layer base elevation, B, for any point and time
-----------------------------------------------------------------------*/
double CSingleLayer::GetBase(const cmplex &z) const {
	static double Bloc;
  if (CSuperblock::SuperblocksOn()) {
		Bloc=pMasterBlock->GetBase(z);
		if (Bloc!=NO_VALUE){return Bloc;     } 
		else               {return base_elev;}
	}
  else{
		return CPropZone::NestedSift(z,pBaseZoneArray,nBaseZones,base_elev);
	}
}
/************************************************************************
                           GetThick
*************************************************************************
 Returns the layer thickness, T for any point and time
-----------------------------------------------------------------------*/
double CSingleLayer::GetThick(const cmplex &z) const {
	static double Tloc;	
  if (CSuperblock::SuperblocksOn()) {
		Tloc=pMasterBlock->GetThick(z);
		if (Tloc!=NO_VALUE){return Tloc;     } 
		else               {return thickness;}
	}
  else{
		return CPropZone::NestedSift(z,pThickZoneArray,nThickZones,thickness);
	}
}
/************************************************************************
                           GetPoro
*************************************************************************
 Returns the porosity, n (or theta) for any point and time
-----------------------------------------------------------------------*/
double CSingleLayer::GetPoro(const cmplex &z) const{
	static double nloc;
	nloc=NO_VALUE;
  if (CSuperblock::SuperblocksOn()) {
		nloc=pMasterBlock->GetPoro(z);
		if (nloc!=NO_VALUE){return nloc;     } 
		else               {return porosity;}
	}
  else{
		return CPropZone::NestedSift(z,pPoroZoneArray,nPoroZones,porosity);
	}
}
/************************************************************************
                           GetDischargePotential
*************************************************************************
 Returns the discharge potential Omega= Phi + iPsi for any point and time
-----------------------------------------------------------------------*/
cmplex CSingleLayer::GetDischargePotential(const cmplex &z,const double &t) const {
  //static cmplex Omega;
	//"static" removed for false potential of area vortex (which has to know the potential to contribute a potential)
	cmplex Omega;
	Omega=0.0;
	if (CTaylorCir::Active){
		for (int i=0; i<nTaylorCirs; i++){
			if (pTaylorCirArray[i]->IsInside(z)){
			  return pTaylorCirArray[i]->GetDischargePotential(z,t); 
			}
		}
	}
  if (CSuperblock::SuperblocksOn()) {
	  Omega+=pMasterBlock->GetDischargePotential(z,t);
	}
  else {
		for (int i=0; i<nElems; i++){
	    Omega+=pElemArray[i]->GetDischargePotential(z,t);
		}
	}
	Omega+= pFarField->GetDischargePotential(z,t); 

	if (pLayerAbove   !=NULL) {Omega+=pLayerAbove  ->GetDischargePotential(BELOW,z,t);}
	if (pLayerBeneath !=NULL) {Omega+=pLayerBeneath->GetDischargePotential(ABOVE,z,t);}
	
  return Omega;
}

/************************************************************************
                           GetW
*************************************************************************
 Returns the complex discharge W = Qx - iQy for any point and time
-----------------------------------------------------------------------*/
cmplex CSingleLayer::GetW(const cmplex &z,const double &t) const {
  static cmplex W;
	W=0.0;
	
	if (!InterpolationOn){
		if (CTaylorCir::Active){
			for (int i=0; i<nTaylorCirs; i++){
				if (pTaylorCirArray[i]->IsInside(z)){
					return pTaylorCirArray[i]->GetW(z,t);
				}
			}
		}
		if (CSuperblock::SuperblocksOn()) {
			W+=pMasterBlock->GetW(z,t);
		}
		else {
			for (int i=0; i<nElems; i++)    {
				W+=pElemArray[i]->GetW(z,t);
			}
		}
		W+=pFarField->GetW(z,t);

		if (pLayerAbove   !=NULL) {W+=pLayerAbove->  GetW(BELOW,z,t);}
		if (pLayerBeneath !=NULL) {W+=pLayerBeneath->GetW(ABOVE,z,t);}
	}
	else if (InterpolationOn){
		if (pInterpGrid==NULL){ExitGracefully("CSingleLayer:GetW: no interpolation grid",BAD_DATA);}
		W=cmplex(pInterpGrid->InterpolateValue(c2Dto3D(z),QxInterp),
						 pInterpGrid->InterpolateValue(c2Dto3D(z),QyInterp));
	}
  return W;
}
/************************************************************************
                           GetGx
*************************************************************************
 Returns the discharge derivative Gx = dQx/dx - idQy/dy for any point and time
-----------------------------------------------------------------------*/
cmplex CSingleLayer::GetGx(const cmplex &z,const double &t) const {
	cmplex Gx(0.0);
	for (int i=0; i<nElems; i++)    {
		Gx+=pElemArray[i]->GetGx(z,t);
	}
	//if (pLayerAbove   !=NULL) {Gx+=pLayerAbove->  GetGx(BELOW,z,t);} //TMP DEBUG
	//if (pLayerBeneath !=NULL) {Gx+=pLayerBeneath->GetGx(ABOVE,z,t);}

	return Gx;
}
/************************************************************************
                           GetLeakage
*************************************************************************
 Returns the leakage N any point and time
  if ltype=FROMTOP ,         returns only leakage from above the layer (i.e. recharge, lakes, rivers, confining layer)
  if ltype=FROMBOTTOM,       returns only leakage from below the layer (i.e. through Aquiclude)
  if ltype=FROMTOPANDBOTTOM, returns net leakage
  leakage is POSITIVE if water is being EXTRACTED from the layer
------------------------------------------------------------------------*/
double CSingleLayer::GetLeakage(const cmplex &z, const double &t,const leak_type ltype) const {
  static double Leak;
	Leak=0.0;
  if (CSuperblock::SuperblocksOn()) {
	  Leak+=pMasterBlock ->GetLeakage(z,t,ltype);
	}
  else {
		for (int i=0; i<nElems; i++)  {
			Leak+=pElemArray[i]->GetLeakage(z,t,ltype);
		} 
	}	

	if ((ltype!=FROMBOTTOM) && (pLayerAbove  !=NULL)){Leak-=pLayerAbove->   GetLeakage(z,t);}
	if ((ltype!=FROMTOP   ) && (pLayerBeneath!=NULL)){Leak+=pLayerBeneath ->GetLeakage(z,t);}

  return Leak;
}
/************************************************************************
                           GetCurl
*************************************************************************
 Returns the Curl (beta) any point and time
------------------------------------------------------------------------*/
double CSingleLayer::GetCurl(const cmplex &z, const double &t) const {
  static double Curl;
	Curl=0.0;
  if (CSuperblock::SuperblocksOn()) {
	  //Curl+=pMasterBlock ->GetLeakage(z,t); //TMP DEBUG
	}
  else {
		for (int i=0; i<nElems; i++)  {
			Curl+=pElemArray[i]->GetCurl(z,t);
		} 
	}	

  return Curl;
}
/************************************************************************
                           GetNetDischarge
*************************************************************************
 Returns the net discharge Q from the layer at some time,t 
-----------------------------------------------------------------------*/
double CSingleLayer::GetNetDischarge(const double &t) const{
	static double Q;
	Q=0.0;
  if (CSuperblock::SuperblocksOn()) {
	  Q+=pMasterBlock->GetNetDischarge(t);
	}
  else {
		for (int i=0; i<nElems; i++)  {
			Q+=pElemArray[i]->GetNetDischarge(t);
		}
	}
	return Q;
}
/************************************************************************
                           GetFluxThruFace
*************************************************************************
 Returns Qn +iQt, the integrated normal & tangential flux between the points z1 and z2 at time t
 Qn=difference in stream function after branch cuts are reoriented
 Qt=difference in potential function after branch cuts are redirected
-----------------------------------------------------------------------*/
cmplex CSingleLayer::GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const double &t) const{
	static cmplex Qflx;

	CAnalyticElem::SetBranchCutLocus(true,0.5*(z1+z2));
	
  Qflx=0.0;
	for (int i=0; i<nElems; i++){
    Qflx+=pElemArray[i]->GetFluxThruFace(z1,z2,t);
	}

  Qflx+=pFarField->GetFluxThruFace(z1,z2,t);

	CAnalyticElem::SetBranchCutLocus(false,0.0);

	return Qflx;
}
/************************************************************************
                           Integrated Leakage
*************************************************************************
 Returns the integrated leakage [L^3/T] over a triangle defined by z1,z2,z3
 Returns positive value for leakage, negative for recharge
------------------------------------------------------------------------*/
double CSingleLayer::GetIntegratedLeakage(const cmplex    &z1, 
																		const cmplex    &z2, 
																		const cmplex    &z3, 
																		const double    &t, 
																		const leak_type  ltype) const {
	static double Q;
	Q=0.0;
  //if (CSuperblock::SuperblocksOn()) {
	  //Q+=pMasterBlock->GetIntegratedLeakage(z1,z2,z3,t,ltype);
	//}
  //else {
		for (int i=0; i<nElems; i++)  {
			Q+=pElemArray[i]->GetIntegratedLeakage(z1,z2,z3,t,ltype);
		}
	//}
	return Q;
}
/************************************************************************
                           Integrated Leakage
*************************************************************************
 Returns the integrated budget (net influx and outflux) over a triangle defined by z1,z2,z3
------------------------------------------------------------------------*/
void CSingleLayer::GetIntegratedBudget(const cmplex &z1, 
																 const cmplex &z2, 
																 const cmplex &z3, 
																 const double &t, 
																       double &inQ, 
																			 double &outQ) const{
  static double tmpin,tmpout;
	inQ=outQ=0.0;
  /*if (CSuperblock::SuperblocksOn()) {
	  Q+=pMasterBlock->GetIntegratedBudget(z1,z2,z3,t,tmpin,tmpout);
		inQ+=tmpin;
		outQ+=tmpout;
	}
  else {*/
		for (int i=0; i<nElems; i++)  {
			pElemArray[i]->GetIntegratedBudget(z1,z2,z3,t,tmpin,tmpout);
			inQ+=tmpin;
			outQ+=tmpout;
		}
	//}
}
/************************************************************************
                           Get Flux Distribution 
*************************************************************************
	Gets integrated extraction between points z1 and z2 along line
	Linearly Distributes these fluxes to endpoints Q1 and Q2
	Used for finite element source distribution
------------------------------------------------------------------------*/
void CSingleLayer::GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const{

  static double tmpQ1,tmpQ2;
	Q1=Q2=0.0;
  /*if (CSuperblock::SuperblocksOn()) {
	  pMasterBlock->GetFluxDistribution(z1,z2,t,tmpQ1,tmpQ2);
		Q1+=tmpQ1;
		Q2+=tmpQ2;
	}
  else {*/
		for (int i=0; i<nElems; i++)  {
			pElemArray[i]->GetFluxDistribution(z1,z2,t,tmpQ1,tmpQ2);
			Q1+=tmpQ1;
			Q2+=tmpQ2;
		}
	//}
}
/************************************************************************
                           Get Singularity Strength  
*************************************************************************
	Returns strength of singularity located at point z (e.g., Q/2pi for a well)
------------------------------------------------------------------------*/
cmplex CSingleLayer::GetSingularStrength  (const cmplex &z, const double &t) const{
  static cmplex S;
	S=0.0;
  /*if (CSuperblock::SuperblocksOn()) {
	  S+=pMasterBlock->GetSingularStrength(z,t);
	}
  else {*/
		for (int i=0; i<nElems; i++)  {
			S+=pElemArray[i]->GetSingularStrength(z,t);
		}
	//}
	return S;
}
/************************************************************************
                           Horizontal Seepage Velocity
*************************************************************************
 Returns the horizontal seepage velocity (vx, vy) at any complex point and time
-----------------------------------------------------------------------*/
cmplex CSingleLayer::GetVelocity2D(const cmplex &z, const double &t) const{
  static cmplex W;
  static double K,T,B,n,head,pot;

  K   =GetCond              (z);
  T   =GetThick             (z);
  n   =GetPoro              (z);
  B   =GetBase              (z);
  W   =GetW                 (z,t);                       

  pot =GetDischargePotential(z,t).real();  

	head=ConvertToHead(pot,K,T,sealevel-B,saltwaterSG);

	if ((abs(W)>=ALMOST_INF) || (T<=0.0) || (n<=0.0) || (head<=0.0)){return 0.0;}

	if (IsConfined(pot,K,T)){return conj(W)/(n*T   );}
	else                    {return conj(W)/(n*head);}
	

}

/************************************************************************
                           Seepage Velocity (3D)
*************************************************************************
 Returns the seepage velocity (vx, vy, vz) at any 3D point and time
 calculates vertical and horizontal velocity component
 MUST OPTIMIZE
------------------------------------------------------------------------*/
vector CSingleLayer::GetVelocity3D(const pt3D &pt, const double &t) const{
	
	static cmplex W,v,z;           //parameters declared static to reduce computational cost of call
	static double pot,head;
	static double leaktop,leakbot;
	static double B,K,T,n,v_z;

	z      =c3Dto2D(pt);

  K      =GetCond              (z); 
	T      =GetThick             (z); 
	n      =GetPoro              (z);
	B			 =GetBase              (z);

	W      =GetW                 (z,t);
	pot    =GetDischargePotential(z,t).real();

	leaktop=GetLeakage           (z,t,FROMTOP   ); //leakup & leakbot are negative if sources of water
	leakbot=GetLeakage           (z,t,FROMBOTTOM);

  head   =ConvertToHead(pot,K,T,sealevel-B,saltwaterSG);	//local head (relative to base)

	if ((abs(W)>=ALMOST_INF) || (T<=0.0) || (n<=0.0) || (head<=0.0)){return 0.0;}

	v_z=ConvertToVertVelocity(pt.z,K,T,B,n,head,pot,W,leaktop,leakbot);
	
	//if ((pt.z<B) || (pt.z>(T+B)) || (pt.z>(head+B))) {return 0.0;}

	if (IsConfined(pot,K,T)){v=conj(W/(n*T   ));}
	else                    {v=conj(W/(n*head));}
  
	return vector(v.real(),v.imag(),v_z);

}
/*************************************************************************
                           Vertical Seepage Velocity
**************************************************************************
 Returns the vertical seepage velocity (vz) at any complex point and time
------------------------------------------------------------------------*/
/*double CSingleLayer::vz(const pt3D &pt, const double &t) const{
	static cmplex W,z; 
	static double K,T,B,n,pot,leaktop,leakbot,head;
	z=c3Dto2D(pt);
  K=GetCond (z); 
	T=GetThick(z); 
	B=GetBase (z);
	n=GetPoro (z);
	pot    =GetDischargePotential(z,t).real();
	W			 =GetW                 (z,t);
	leaktop=GetLeakage           (z,t,FROMTOP   ); //leakup & leakbot are negative if sources of water
	leakbot=GetLeakage           (z,t,FROMBOTTOM); 
	head   =ConvertToHead        (pot,K,T,sealevel-B,saltwaterSG);

	if ((fabs(W.real())>=ALMOST_INF) ||(fabs(W.imag())>=ALMOST_INF)){W=0.0;}

	return ConvertToVertVelocity(pt.z,K,T,B,n,head,pot,W,leaktop,leakbot);
}*/
//----------------------------------------------------------------------------
double CSingleLayer::ConvertToVertVelocity(const double &z,const double &K,const double &T, 
																		 const double &B, const double &n,const double &head,
																		 const double &pot,const cmplex &W,const double &leaktop,const double &leakbot) const{
	if ((T==0.0) ||(n==0.0)){return 0.0;}
	if (head>0){
		if (IsConfined(pot,K,T)){ //confined
			return ((leakbot+leaktop)/T*(z-B)-leakbot)/n;
		}
		else{                     //unconfined
			return ((-1.0/head)*(0.5*abs(W)*abs(W)/pot-(leakbot+leaktop))*(z-B)-leakbot)/n; //revised
		//return (-(abs(W)*(1.0/hz-1.0/hz2)/step_size)-(leak/hz))*h+leakbot; //from Strack, 1984
		}
	}	
	else{											  //dry aquifer
		return 0.0; 
	}
}
/************************************************************************
                           Effective Velocity
*************************************************************************
 Returns the effective seepage velocity (vex, vey, vz) at any 3D point and time
 calculates horizontal velocity components only (for now)
 MUST OPTIMIZE
------------------------------------------------------------------------*/
vector CSingleLayer::GetEffVelocity3D(const pt3D &pt,  const double &t, const disp_info &disp) const{
	static cmplex W,v,z,Gx,Gy, nprime;           //parameters declared static to reduce computational cost of call
	static double pot,head,vex,vey;
	static double leaktop,leakbot;
	static double B,K,T,n,Dxx,Dyy,Dxy,Tx,Ty;

	z=c3Dto2D(pt);

  K      =GetCond              (z); 
	T      =GetThick             (z); 
	n      =GetPoro              (z);
	B			 =GetBase              (z);
	nprime=0.0;//GetPoroSlope(z);

	pot    =GetDischargePotential(z,t).real();
	W      =GetW                 (z,t);
	Gx		 =GetGx                (z,t);

	head   =ConvertToHead(pot,K,T,sealevel-B,saltwaterSG); //local head

	leaktop=GetLeakage           (z,t,FROMTOP   ); //leakup & leakbot are negative if sources of water
	leakbot=GetLeakage           (z,t,FROMBOTTOM); 
	Gy     =-Gx+(leaktop+leakbot);

	if (abs(W)>=ALMOST_INF){W=0.0;}
	//if ((fabs(W.real())>=ALMOST_INF) ||(fabs(W.imag())>=ALMOST_INF)){W=0.0;}
	
	
	CVector effv= c2Dto3D(ConvertToEffective(pot,head,K,T,n,nprime,Gx,Gy,W,disp));
	effv.z=ConvertToVertVelocity(pt.z,K,T,B,n,head,pot,W,leaktop,leakbot);
	//effv.z=ConvertToEffectiveVertical(...)

	return effv;
}
//---------------------------------------------------------------------------
cmplex CSingleLayer::GetEffVelocity2D(const cmplex &z, const double &t, const disp_info &disp) const{	
	static cmplex W,v,Gx,Gy, nprime;           //parameters declared static to reduce computational cost of call
	static double pot,head;
	static double leakage;
	static double B,K,T,n;

  K      =GetCond              (z); 
	T      =GetThick             (z); 
	n      =GetPoro              (z);
	B			 =GetBase              (z);
	nprime=0.0;//GetPoroSlope          (z);

	pot    =GetDischargePotential(z,t).real();
	W      =GetW                 (z,t);
	Gx		 =GetGx                (z,t);

	head   =ConvertToHead(pot,K,T,sealevel-B,saltwaterSG); //local head

	leakage=GetLeakage           (z,t,FROMTOPANDBOTTOM); //negative if sources of water
	Gy     =IM*Gx + IM*leakage; //+curl;

	if (abs(W )>=ALMOST_INF){W    =0.0;}
	if (abs(Gx)>=ALMOST_INF){Gx=Gy=0.0;}

	return ConvertToEffective(pot,head,K,T,n,nprime,Gx,Gy,W,disp);
	
}
//---------------------------------------------------------------------------
cmplex CSingleLayer::ConvertToEffective(const double &pot, const double &head, 
													        const double &K,   const double &T,  const double &n, const cmplex    &nprime, 
													        const cmplex &Gx,  const cmplex &Gy, const cmplex &W, const disp_info &disp) const{
	bool isconfined;
  static cmplex v;
	static double Dxx,Dyy,Dxy,vex,vey,Tx,Ty,vnrm,Wnrm,vx,vy,dhdx,dhdy;
	double aL=disp.al;
	double aT=disp.at;

	if (head>0){
		isconfined=IsConfined(pot,K,T);

		Tx=-nprime.real()/n;
		Ty=-nprime.imag()/n;

		if (isconfined){
			v=conj(W/(n*T)); 
			dhdx=dhdy=0.0;  
		}
		else           {
			v=conj(W/(n*head));
			Tx  += W.real()/head/head/K; 
			Ty  +=-W.imag()/head/head/K;
			dhdx =-W.real()/head     /K;
			dhdy = W.imag()/head     /K;
		}

		vx=v.real();
		vy=v.imag();
		vnrm=abs(v);
		Wnrm=abs(W);
		if (vnrm==0.0){return 0.0;}

		Dxx=(   aL)*vx*vx/vnrm+(aT)*vy*vy/vnrm+disp.D;
		Dyy=(   aT)*vx*vx/vnrm+(aL)*vy*vy/vnrm+disp.D;
		Dxy=(aL-aT)*vx*vy/vnrm					  		+disp.D;

		//Tx=Ty=0.0;
		//TMP DEBUG-plots only d|v|/dx, d|v|/dy-------------------------
		/*double vmagdx(0.0),vmagdy(0.0);
    if (abs(v)!=0.0){vmagdx=abs(v)*( (conj(W)*Gx).real()/Wnrm/Wnrm+Tx);} //WORKING!!!
    if (abs(v)!=0.0){vmagdy=abs(v)*( (conj(IM*W)*Gy).real()/Wnrm/Wnrm+Ty);}//WORKING!!!-CHANGE TO DERIVATION
		return cmplex(vmagdx,vmagdy);*/

		//TMP DEBUG-plots only dvx/dx, dvy/dx-------------------------
		/*double vxdx(0.0),vydx(0.0);
    if (vx!=0.0){vxdx=vx*(Gx.real()/W.real()+Tx);} //WORKING!!!
    if (vy!=0.0){vydx=vy*(Gx.imag()/W.imag()+Tx);} //WORKING!!! 
		return cmplex(vxdx,vydx);*/

		//TMP DEBUG-plots only dvx/dy, dvy/dy-------------------------
		/*double vxdy(0.0),vydy(0.0);
    if (vx!=0.0){vxdy=vx*(-Gy.imag()/W.real()+Ty);} //WORKING!!!-CHANGE TO DERIVATION
    if (vy!=0.0){vydy=vy*( Gy.real()/W.imag()+Ty);} //WORKING!!!-CHANGE TO DERIVATION
		return cmplex(vxdy,vydy);*/

		double DxxDx(0.0),DyyDy(0.0),DxyDx(0.0),DxyDy(0.0);

		if (vx   !=0.0){DxxDx+=aL*vx*vx/vnrm*(2.0*Gx.real()/W.real()-(conj(W)*Gx).real()/Wnrm/Wnrm+Tx);}//WORKING!!!
		if (vy   !=0.0){DxxDx+=aT*vy*vy/vnrm*(2.0*Gx.imag()/W.imag()-(conj(W)*Gx).real()/Wnrm/Wnrm+Tx);}//WORKING!!!

		if (vx   !=0.0){DyyDy+=aT*vx*vx/vnrm*( 2.0*Gy.imag()/W.real()-(conj(IM*W)*Gy).real()/Wnrm/Wnrm+Ty);}//WORKING!!!-CHANGE TO DERIVATION
		if (vy   !=0.0){DyyDy+=aL*vy*vy/vnrm*(-2.0*Gy.real()/W.imag()-(conj(IM*W)*Gy).real()/Wnrm/Wnrm+Ty);}//WORKING!!!-CHANGE TO DERIVATION
		
    if (vx*vy!=0.0){DxyDx=(aL-aT)*vx*vy/vnrm*(Gx.real()/W.real()+Gx.imag()/W.imag()-(conj(W)*Gx).real()/Wnrm/Wnrm+Tx);}//WORKING!!!

    if (vx*vy!=0.0){DxyDy=(aL-aT)*vx*vy/vnrm*(Gy.imag()/W.real()+(-Gy.real())/W.imag()-(conj(IM*W)*Gy).real()/Wnrm/Wnrm+Ty);}//WORKING!!!-CHANGE TO DERIVATION
	
		//TMP DEBUG-plots only DxxDx, DyyDy-----------------------------
		//return cmplex(DxxDx,DyyDy);
		//TMP DEBUG-plots only DxyDx, DxyDy-----------------------------		
		//return cmplex(DxyDx,DxyDy);

		//plus equals is actually v* =v (+) d/di(nhD dC/di), where (+) should be (-) for "true" effective velocity
		vex=vx;
		vex+=(DxxDx+DxyDy+Dxx/n*nprime.real()+Dxy/n*nprime.imag()+Dxx/head*dhdx+Dxy/head*dhdy);

		vey=vy;
		vey+=(DyyDy+DxyDx+Dyy/n*nprime.imag()+Dxy/n*nprime.real()+Dyy/head*dhdy+Dxy/head*dhdx);

		return cmplex(vex,vey);

	}	
	else{											  //dry aquifer
		return 0.0; 
	}
}

/************************************************************************
                           GetHead
*************************************************************************
 Returns the GLOBAL head value (phi_local+B) at any complex point and time
-----------------------------------------------------------------------*/
double CSingleLayer::GetHead(const cmplex &z, const double &t) const {
	static double K,B,T,head;
  K=GetCond(z);	
	B=GetBase(z);
	T=GetThick(z);
	head=ConvertToHead(GetDischargePotential(z,t).real(),K,T,sealevel-B,saltwaterSG);
  if (head<=0.0){return B;}
  else          {return head+B;}
}
/*************************************************************************
                           GetSaturatedThickness
**************************************************************************
 Returns the saturated thickness of the aquifer layer
------------------------------------------------------------------------*/
double CSingleLayer::GetSaturatedThickness(const cmplex &z, const double &t) const{
	return GetSaturatedThick(GetDischargePotential(z,t).real(),GetCond(z),GetThick(z));
}
/************************************************************************
                           GetHeadAndPotential
*************************************************************************
 Returns the GLOBAL head value (phi_local+B) and discharge potential, Omega, 
 at any complex point and time
-----------------------------------------------------------------------*/
void CSingleLayer::GetHeadAndPotential(const cmplex &z,cmplex &omega, double &head,const double &t) const {
	//returns GLOBAL head values, local potential
	static double K,B,T;
  K=GetCond(z);	
	B=GetBase(z);
	T=GetThick(z);
  omega=GetDischargePotential(z,t);
	head =ConvertToHead(omega.real(),K,T,sealevel-B,saltwaterSG);
  if (head<=0.0){head=B;}
  else          {head=head+B;}
}
/************************************************************************
                           GetBaseflow
*************************************************************************
   Returns the cumulative baseflow in any linesink or other extraction element
	 located within NEAR_FEATURE of the point
-----------------------------------------------------------------------*/
double CSingleLayer::GetBaseflow(const cmplex &z, const double t) const{
	double tmpBF=0;

	for (int i=0; i<nElems; i++){
	  tmpBF=max(pElemArray[i]->GetCumulativeFlux(z,t),tmpBF);
	}

  return tmpBF;
}

/************************************************************************
                           IterativeSolve
*************************************************************************
 Iteratively solves for the coefficients of each element in the layer
	 maxchange is the maximum change of any element during the current iteration
	 maxobj    is the maximum error of any element during the current iteration
	 PROGRESS  is the output progress file
------------------------------------------------------------------------*/
void CSingleLayer::IterativeSolve(double &maxchange, double &maxobj, const double &t, ofstream &PROGRESS){

	double         tempchange,tempobj;
  int            i;
	CAnalyticElem *WorstChange;
	CAnalyticElem *WorstObj;
	double         AvgObj;
	bool           stopped(false);
	ofstream       ANALYSIS;

  cout << "Solving Layer #" <<level+1 << " ("<<nElems << " element(s))..."<<endl;
	cout << "------------------------------------------------------------"<<endl;

  PROGRESS << "solve"<<endl;

	CSingleLayer::solved=false;
	CTaylorCir::Active=false;	
	
	iter=1; 

	if (CSingleLayer::FullAnalyze){
		ANALYSIS.open("iterative_analysis.csv");
    StartAnalysis(ANALYSIS);
	}

	//Solve given elements first
	for (i=0; i<nGivenElems; i++){
		pGivenElems[i]->SolveItself(tempchange,tempobj,0.0);
	}

	while ((!stopped) &&
		     (((iter<=CSingleLayer::MaxLayerIterations) && (maxchange>CSingleLayer::LayerTolerance)) ||
				   (iter<=CSingleLayer::MinLayerIterations)) ) { 
		
		tempchange =tempobj=0.0;	
		maxchange  =maxobj =0.0;
		WorstChange=NULL;
		WorstObj   =NULL;

		if (!stopped){

			AvgObj=0.0;
      IdentifyDeltaPotential(t);

			pFarField->SolveItself(tempchange,tempobj,t);
			upperswap(maxobj   ,tempobj   );if (tempobj   ==maxobj   ) {WorstObj   =pFarField;}
      upperswap(maxchange,tempchange);if (tempchange==maxchange) {WorstChange=pFarField;}  
			AvgObj+=tempobj;
			if (CSingleLayer::FullAnalyze){PrintAnalysis(ANALYSIS,iter,-1,tempchange,t);}

			//TMP DEBUG-A1REVISE
			//Solve for fluxes of explicit flux elements
			if (directflux){
				//BuildDirectFluxMatrix();
				//SolveDirectFluxMatrix();
			}

			for (i=0; i<nElems; i++){
				tempchange=tempobj=0.0;			

				if (ProgramAborted())  { stopped=true;i=nElems;}

				if (!stopped) {

					WriteEllipsisToScreen(i,nElems,40);

					/*if (CSingleLayer::FFeverytime){ //SHOULD PROBABLY REMOVE THIS - should be handled by explicit flux method
						pFarField->SolveItself(tempchange,tempobj,t); 
						upperswap(maxobj   ,tempobj   );if (tempobj   ==maxobj   ){WorstObj   =pFarField;}
            upperswap(maxchange,tempchange);if (tempchange==maxchange){WorstChange=pFarField;}  
						if (CSingleLayer::FullAnalyze){PrintAnalysis(ANALYSIS,iter,-1,tempchange,t);}
					}*/

					if (!pElemArray[i]->IsGiven()){ 
						pElemArray[i]->SolveItself(tempchange,tempobj,t);
						upperswap(maxobj   ,tempobj   );if (tempobj   ==maxobj   ) {WorstObj   =pElemArray[i];}
						upperswap(maxchange,tempchange);if (tempchange==maxchange) {WorstChange=pElemArray[i];}  
						AvgObj+=tempobj;
						if (CSingleLayer::FullAnalyze){PrintAnalysis(ANALYSIS,iter,i,tempchange,t);}
					}

 					if (CSingleLayer::GridWhileSolve) {
						CSingleLayer::solved=true; //so stream function is cleaned up
						Grid(CSingleLayer::GWSW,CSingleLayer::GWSresolution,iter,i,t);
						CSingleLayer::solved=false; 
					} 
				}
			} 
      if (nElems>40) {cout<<endl;}

			cout << "iteration "      << iter
				   << ",   maxchange: " << maxchange
				   << ",    max. obj: " << maxobj;
			if (WorstObj   !=NULL){cout << "  (worst obj: "   <<WorstObj   ->GetName()<<") ";}
			if (WorstChange!=NULL){cout << "  (worst change: "<<WorstChange->GetName()<<") ";}
			cout << ",   avg. obj: " << AvgObj/nElems<< endl;

			PROGRESS << iter << " "<<maxchange<<" "<<maxobj<<" "<<GetNetDischarge(t)<<endl;
			
			if (iter%WriteInterval==0){
				cout << "Writing Back-up Solution..."<<endl;
				ofstream SOL;
				SOL.open("solution.bbs");
				SOL.precision(10);
				WriteItself(SOL,t);
				SOL.close();
			}

			iter++;
			CSingleLayer::fresh=false;
		}
  }
	CSingleLayer    ::solved=true;
	CTaylorCir::Active=true;

	if (iter<CSingleLayer::MaxLayerIterations)  {cout << "...success!"                         <<endl;}
	else if (stopped)                     {cout << "...STOP file generated. Please wait."<<endl;}
  else                                  {cout << "...exceeded maximum layer iterations"<<endl;}
	PROGRESS<<"done"<<endl;
	if (CSingleLayer::FullAnalyze){ANALYSIS.close();}
}
/************************************************************************
                           ExplicitSolve
*************************************************************************
 Explicitly solves for the coefficients of each element in the layer
	 PROGRESS  is the output progress file
-----------------------------------------------------------------------*/
void CSingleLayer::ExplicitSolve        (ofstream &PROGRESS, const double &t){
	
	int i,l,count;
	CAnalyticElem **AllElems;
	int            *AllSegs;

	//sift through elements, cout segments
	count=0;
	for (i=0; i<nElems; i++){
		if (!pElemArray[i]->IsString()){count++;}
		else                           {count+=pElemArray[i]->GetNumSegs();}
	}
	count++; //for farfield

	AllElems=new CAnalyticElem *[count];
	AllSegs= new int            [count];

	count=0;
	for (i=0; i<nElems; i++){
		if (!pElemArray[i]->IsString()){
			AllElems[count]=pElemArray[i];
			AllSegs [count]=NOT_A_STRING;
			cout <<"ExplicitSolve: "<< AllSegs [count]<<endl;
		  count++;
		}
		else                           {
			for (l=0; l<pElemArray[i]->GetNumSegs(); l++){
				if (!pElemArray[i]->IsDisabled(i)){
					AllElems[count]=pElemArray[i];
					AllSegs [count]=l;
					count++;
				}
			}
		}

 	}
  AllElems[count]=pFarField;
	AllSegs[count]=NOT_A_STRING;

	cout << "Solving Explicitly..."<<endl;
	BuildMatrix(AllElems,AllSegs,count+1,NULL,t); //Builds and solves problem
	cout << "...Done solving explicitly"<<endl;

	delete [] AllElems; //deletes pointers only
	delete [] AllSegs;
}
/************************************************************************
                           IdentifyDeltaPotential
*************************************************************************
  Takes 16 points across domain to estimate potential range (max-min)
-----------------------------------------------------------------------*/
void CSingleLayer::IdentifyDeltaPotential(const double t){
	
	double minpot(ALMOST_INF), maxpot(-ALMOST_INF),pot(0.0);
	if ((Extents.e!=Extents.w) && (Extents.n!=Extents.s)){
		for   (double x=Extents.w; x<=Extents.e+(Extents.e-Extents.w)*0.1; x+=(Extents.e-Extents.w)/4.0){
			for (double y=Extents.s; y<=Extents.s+(Extents.n-Extents.s)*0.1; y+=(Extents.n-Extents.s)/4.0){
				pot=GetDischargePotential(cmplex(x,y),t).real();
				if ((pot<minpot) && (pot >-ALMOST_INF)){minpot=pot;}
				if ((pot>maxpot) && (pot <ALMOST_INF) ){maxpot=pot;}
			}
		}
	}
	DeltaPot=maxpot-minpot;
	if (DeltaPot<1.0){DeltaPot=1.0;}

}
/************************************************************************
                           InterpolateQxQy
*************************************************************************
  Obtain Qx, Qy values and store in memory for interpolation
	Enables interpolation of W (and thus velocity)
-----------------------------------------------------------------------*/
void CSingleLayer::InterpolateQxQy      (const double &t){
	int k;
	int nNodes;
	cmplex W;
	cout <<"Creating Interpolation Grid";
	if (pInterpGrid==NULL){ExitGracefully("CSingleLayer::InterpolateQxQy: NULL grid",BAD_DATA);}
	nNodes=pInterpGrid->GetNumNodes();
	QxInterp=new double [nNodes];
	QyInterp=new double [nNodes];
	InterpolationOn=false;
	for (k=0; k<nNodes; k++){
		WriteEllipsisToScreen(k,nNodes,20);
		W=GetW(c3Dto2D(pInterpGrid->GetNodeLocation(k)),t); 
		QxInterp[k]=W.real();
		QyInterp[k]=-W.imag();
	}
	cout <<endl;
	InterpolationOn=true;
}
/************************************************************************
                           Write Itself/Write Output
*************************************************************************
  CSingleLayer::WriteItself: Writes solution file
  CSingleLayer::WriteOutput: Writes generic output
-----------------------------------------------------------------------*/
void CSingleLayer::WriteItself(ofstream &SOL, const double &t) const{
	pFarField->WriteItself(SOL,t);
	for (int i=0; i<nElems; i++){				
		WriteEllipsisToScreen(i,nElems,30);
    pElemArray[i]->WriteItself(SOL,t);
	}
}
//-----------------------------------------------------------------------
void CSingleLayer::WriteOutput(const double &t) const{
	bool stopped(false);
	pFarField->WriteOutput(t);
	for (int i=0; i<nElems; i++){
		if (ProgramAborted())  {stopped=true;i=nElems;}
		WriteEllipsisToScreen(i,nElems,30);
		if (!stopped){pElemArray[i]->WriteOutput(t);}
	}
  if (CSuperblock::SuperblocksOn()) {cout<<endl;}
}
/************************************************************************
                           ANALYSIS TOOLS
*************************************************************************
 writes extensive output during iterative solve
-----------------------------------------------------------------------*/
void CSingleLayer::StartAnalysis(ofstream &ANALYSIS){
	ANALYSIS<< "iter,elem,change,FarField,";
	for (int i=0; i<nElems; i++){
		ANALYSIS << pElemArray[i]->GetName()<<",";
	}
	ANALYSIS << endl;
}
//-----------------------------------------------------------------------
void CSingleLayer::PrintAnalysis(ofstream &ANALYSIS, const int iter, const int i,const double change, const double &t){
	ANALYSIS<<iter<< ","<< i<<","<<change<<",";
	ANALYSIS<<pFarField->GetMaxError(t)<<",";
	for (int j=0; j<nElems; j++){
		ANALYSIS << pElemArray[j]->GetMaxError(t)<<",";
	}
	ANALYSIS << endl;
}
/************************************************************************
                           GRIDDING FUNCTIONS
************************************************************************/
void CSingleLayer::Grid(window W, int resolution, int iter, int i, double t) const{

  ostrstream hfile;
	ostrstream sfile;
	char filename1[1000];
	char filename2[1000];
	filename1[0]='\0'; 
	filename2[0]='\0';

	hfile << "head"<<iter<<"-"<<i<<".grd"<<'\0';
	sfile << "stream"<<iter<<"-"<<i<<".grd"<<'\0';

  strcpy(filename1,hfile.str());
  strcpy(filename2,sfile.str());
	SetPlotLayer(this);
	cGrid(filename1,filename2,W,resolution,grdHeadStream,t); 
 
}
/***********************************************************************
				SolveFluxMatrix
************************************************************************
Builds & solves explicit matrix for solution of head-specified Analytic element coefficients
-----------------------------------------------------------------------*/
void CSingleLayer::SolveFluxMatrix      (double &maxchange, double &maxobj, const double &t){
	double	 **AA;
	double		*BB;
	double		*SS;
	int				 i,j,m,col,row;
	MatrixInfo MI;

	//reserve memory for matrix
	AA        = new double *[nHSElems];
	BB        = new double  [nHSElems];
	SS        = new double  [nHSElems];	
	if (SS==NULL){ExitGracefully("BuildMatrix: Out of memory(2)",OUT_OF_MEMORY);}
	for (i=0;i<nHSElems;i++){
		AA[i]=new double  [nHSElems];
		for (j=0;j<nHSElems;j++){
			AA[i][j]=0.0;
		}
		BB[i]=0.0;
		SS[i]=0.0;
	}

	//go through column by column (a1 by a1)
	for (row=0; row<nHSElems;row++){
		pHSElems[row]->GetFluxMatrixBuildInfo(MI); 

		cout << "building matrix...row "<< row <<endl;
			
		//go thru column by column
		for (col=0; col<nHSElems; col++){	

			for (m=0; m<MI.nctrl; m++){
				AA[row][col]+=MI.unit[0][m]*MI.unit[0][m]; 
			}	
		}
		 
	}

	for (i=0;i<nHSElems;i++){delete [] AA[i];	}
	delete [] AA;
	delete [] BB;
	delete [] SS;

}

	//Q=GetDischargePotential(z1,t).imag()-GetDischargePotential(z2,t).imag();

  //identify branch cut discontinuity
	
	/*bool   discontinuous(true);
	double Psi1,Psi2,PsiLeft,PsiRight,tmpPsi1,tmpPsi2;
	double slope,slope1,slope2,slope3,maxslope,lastslope(0.0);
	cmplex zleft,zright;
	double eta=0.2;

	Psi1=GetDischargePotential(z1,t).imag();
	Psi2=GetDischargePotential(z2,t).imag();
	zleft =z1; PsiLeft =Psi1;
	zright=z2; PsiRight=Psi2;

	slope=(Psi2-Psi1)/abs(z2-z1); //dPhi/dX

	while ((abs(zright-zleft)>0.01*abs(z1-z2)) && (discontinuous) ){

		//calculate 2 intermediate Psi Values at X=-1/3, X=1/3
		tmpPsi1=GetDischargePotential((1.0*(zright-zleft)/3.0+zleft),t).imag();
		tmpPsi2=GetDischargePotential((2.0*(zright-zleft)/3.0+zleft),t).imag(); 
  
		//Calculate slopes dPhi/dX
		slope1=(tmpPsi1 -PsiLeft )/(abs(zright-zleft)/3.0);
		slope2=(tmpPsi2 -tmpPsi1 )/(abs(zright-zleft)/3.0);
		slope3=(PsiRight-tmpPsi2 )/(abs(zright-zleft)/3.0);

		maxslope=max(max(fabs(slope1),fabs(slope2)),fabs(slope3));

		if ((maxslope)<(2.0*(lastslope))){ 
			//if relative slope is not growing rapidly (should grow by factor of ~3), then it is not a discontinuity
			PsiLeft=PsiRight;
			discontinuous=false;
		}
		else if ((maxslope==0.0) ||                 //may not be neccessary
			((fabs(slope-slope1)/(maxslope)<eta) && 
			 (fabs(slope-slope2)/(maxslope)<eta) && 
			 (fabs(slope-slope3)/(maxslope)<eta))){
			PsiLeft=PsiRight;
			discontinuous=false;
		}
		else if (maxslope==fabs(slope1)){ //discontinuity in first third
			zright  =1.0*(zright-zleft)/3.0+zleft;
			PsiRight=tmpPsi1;
		}
		else if (maxslope==fabs(slope2)){ //discontinuity in middle third 
			zleft   =1.0*(zright-zleft)/3.0+zleft;
			zright  =2.0*(zright-zleft)/3.0+zleft;
			PsiRight=tmpPsi2;
			PsiLeft =tmpPsi1;
		}
		else if (maxslope==fabs(slope3)){   //discontinuity in last third
			zleft   =2.0*(zright-zleft)/3.0+zleft;
			PsiLeft =tmpPsi2;
		}
		else{
			//??
			cout <<"BADDD"<<endl;
			discontinuous=false;
		}

		lastslope=fabs(PsiRight-PsiLeft)/abs(zright-zleft);
	} //End While loop

	Q=(Psi1-PsiLeft)+(PsiRight-Psi2);

	if (discontinuous){
		//cout <<"discontinuous"<<endl;
		ofstream DISCONT;
		DISCONT.open("discontinuity.bna",ios::app);
		DISCONT << "\" Discontinuity \",  1" <<endl;
		DISCONT << zleft.real()<< " , "    <<zleft.imag()<<endl;
		DISCONT.close();
	}*/