
#include "MultiLayer.h"

/*************************************************************************
                             CMultiLayer
************************************************************************
                           CONSTRUCTORS
************************************************************************
************************************************************************/
CMultiLayer::CMultiLayer(){ 
  black_hole=0;  
	pLayerAbove=pLayerBeneath=NULL; 
	pFarField        =NULL;
  pElemArray       =NULL; nElems=0;
	pHSElems         =NULL; nHSElems=0;
	pGivenElems      =NULL; nGivenElems=0;
	pCondZoneArray   =NULL; nCondZones=0;
	pThickZoneArray  =NULL; nThickZones=0;
	pPoroZoneArray   =NULL; nPoroZones=0;
	pAqThickZoneArray=NULL; nAqThickZones=0;
	pAqPoroZoneArray =NULL; nAqPoroZones=0;
	conductivity		 =NULL;
  thickness				 =NULL;       
  porosity				 =NULL;
	aqtardporo			 =NULL;
	aqtardcond			 =NULL;
	aqtardthick			 =NULL;
	conductance			 =NULL;
	top_elev				 =0.0;

	DeltaPot=1.0;
	//pMasterBlock   =NULL;
	
	nLayers=1.0;

  Extents.w= ALMOST_INF; Extents.s= ALMOST_INF;
  Extents.e=-ALMOST_INF; Extents.n=-ALMOST_INF;
}
//------------------------------------------------------------------------
CMultiLayer::~CMultiLayer(){ 
	if (globaldebug){cout <<" DESTROYING MULTILAYER "<<endl;}
	delete [] pElemArray;     //just pointers: elements & propZones deleted globally
	delete [] pHSElems;
	delete [] pGivenElems;
	delete [] pCondZoneArray;
	delete [] pThickZoneArray;
	delete [] pPoroZoneArray;
	delete [] pAqCondZoneArray;
	delete [] pAqThickZoneArray;
	delete [] pAqPoroZoneArray;

	delete [] conductivity;
  delete [] thickness;       
  delete [] porosity;
	delete [] aqtardporo;
	delete [] aqtardcond;
	delete [] aqtardthick;
	delete [] conductance;
}
/************************************************************************
                 STATIC INITIALIZERS (default values)
*************************************************************************
************************************************************************/
bool   CMultiLayer::solved							=true;
int    CMultiLayer::iter								=0;
int    CMultiLayer::MinMLayerIterations	=1; 
int    CMultiLayer::MaxMLayerIterations	=500;
double CMultiLayer::MLayerTolerance			=0.00001;
int    CMultiLayer::WriteInterval				=1000000;
cmplex CMultiLayer::black_hole					=500.0;
/************************************************************************
************************************************************************
                           ACCESSOR FUNCTIONS
************************************************************************
************************************************************************/

/************************************************************************
 GetBlackHole
	returns location where logarithm terms go to zero (black hole)
-----------------------------------------------------------------------*/
cmplex						CMultiLayer::GetBlackHole()   const    {return CMultiLayer::black_hole;}
//------------------------------------------------------------------------
CmlAnalyticElem  *CMultiLayer::GetElem(int i)   const    {if ((i<nElems) && (i>=0)){return pElemArray[i];} else{return NULL;}}   
//------------------------------------------------------------------------
int								CMultiLayer::GetNumElems()    const    {return nElems;}   
//------------------------------------------------------------------------
window						CMultiLayer::GetExtents()     const    {return Extents;}
//------------------------------------------------------------------------
double						CMultiLayer::GetDeltaPot()    const    {return DeltaPot;}
//------------------------------------------------------------------------
int								CMultiLayer::GetCurrentIter() const    {return iter;}
//------------------------------------------------------------------------
int								CMultiLayer::GetNumLayers()   const    {return nLayers;}
/************************************************************************
*************************************************************************
                           ASSIGNMENT FUNCTIONS
*************************************************************************
************************************************************************/
void CMultiLayer::SetLevelAbove (CAquicludeABC *above)  {pLayerAbove=above;}
//------------------------------------------------------------------------
void CMultiLayer::SetLevelBelow (CAquicludeABC *below)  {pLayerBeneath=below;}


/************************************************************************
                           Set Black Hole
*************************************************************************
	places black hole where unif. flow potential is zero and far from centroid of domain
	this is where log terms and 1/z should all vanish
------------------------------------------------------------------------*/
void CMultiLayer::SetBlackHole(){
  cmplex zc;

	if (pFarField==NULL){
		ExitGracefully("CMultiLayer::SetBlackHole: FarField doesn't exist", BAD_DATA);}

	if (pFarField->IsReferencePt()){
		black_hole=pFarField->GetZref();
	}
	else{
		ExitGracefully("CMultiLayer::SetBlackHole():Non-reference point far-field",BAD_DATA);
	}
}
/************************************************************************
                           SET NUMBER OF LAYERS
*************************************************************************
	sets the background values for base, thickness, conductivity,porosity and level
	must be done before any other steps
-----------------------------------------------------------------------*/
void CMultiLayer::SetNumLayers(int N){

	if (N<=0){
		ExitGracefully("CMultiLayer::SetBackgroundValues: bad number of layers specified",BAD_DATA);}

	nLayers=N;

	conductivity=new double [nLayers];        //background conductivity
  thickness		=new double [nLayers];        //background aquifer thickness
  porosity		=new double [nLayers];        //background porosity
	aqtardporo	=new double [nLayers-1];        //background aquitard porosity
	aqtardcond	=new double [nLayers-1];        //background aquitard K
	aqtardthick	=new double [nLayers-1];        //background aquitard thickness
	conductance	=new double [nLayers-1];        //background aquitard conductance (K*t)
	for (int l=0; l<nLayers-1;l++){ //default values
		conductivity[l]=1.0;
		thickness		[l]=5.0;
		porosity		[l]=0.3;
		aqtardporo	[l]=0.3;
		aqtardcond	[l]=0.001;
		aqtardthick	[l]=0.1;
		conductance	[l]=0.01;
	}
	conductivity[nLayers-1]=1.0;
	thickness		[nLayers-1]=5.0;
	porosity		[nLayers-1]=0.3;
}
/************************************************************************
                           SET BACKGROUND VALUES
*************************************************************************
	sets the background values for base, thickness, conductivity,porosity and level
-----------------------------------------------------------------------*/
void CMultiLayer::SetBackgroundValues(double top,
																			Ironclad1DArray T, //size N
																			Ironclad1DArray K, //size N
																			Ironclad1DArray N, //size N
																			Ironclad1DArray t, //size N-1 
																			Ironclad1DArray k, //size N-1
																			Ironclad1DArray n){//size N-1
	int l;

	if (nLayers<=0){
		ExitGracefully("CMultiLayer::SetBackgroundValues: number of layers unknown",BAD_DATA);}
	total_transmissivity=0.0;
	top_elev=top;
	for (l=0; l<nLayers; l++){
		if ((K[l]<0) || (k[l]<0)){
			ExitGracefully("CMultiLayer::SetBackgroundValues: negative conductivity",BAD_DATA);}
		if ((N[l]<0) || (N[l]>=1) || (n[l]<0) || (n[l]>=1) ){
			ExitGracefully("CMultiLayer::SetBackgroundValues: bad porosity value   ",BAD_DATA);}
		if ((T[l]<=0) || (t[l]<0)){
			ExitGracefully("CMultiLayer::SetBackgroundValues: bad thickness value   ",BAD_DATA);}

		conductivity[l]=K[l];
		thickness		[l]=T[l];
		porosity		[l]=N[l];
		if (l<nLayers-1){
			aqtardporo	[l]=n[l];
			aqtardcond	[l]=k[l];
			aqtardthick	[l]=t[l];
			base_elev-=(T[l]);
		}
		else{
			base_elev-=(T[l]+t[l]);
		}
		if (t[l]!=0.0){
			conductance	[l]=k[l]/t[l];
		}
		else{
			conductance	[l]=0.0;//calculated below...
		}
		total_transmissivity+=K[l]*T[l]; //TMP DEBUG - is this correct
	}
	for (l=0; l<nLayers-1; l++){
		if (t[l]==0.0){
			//there are other options for calculating conductance
			conductance	[l]=k[l]/0.5*(T[l]+T[l+1]);
		}
	}
}

/************************************************************************
                           CHANGE EXTENTS
*************************************************************************
	increases the extents of the domain window ("Extents") to include the point specified
------------------------------------------------------------------------*/
void CMultiLayer::UpdateExtents(const cmplex z){
	upperswap(Extents.e,z.real()); lowerswap(Extents.w,z.real());
	upperswap(Extents.n,z.imag()); lowerswap(Extents.s,z.imag());
	if (Extents.e<z.real()){Extents.e=z.real();}
}
/*------------------------------------------------------------------------
	increases the extents of the domain window ("Extents") to include the circle specified
------------------------------------------------------------------------*/
void CMultiLayer::UpdateExtents(const cmplex z,const double r){
	UpdateExtents(z+r);	   UpdateExtents(z-r);	
	UpdateExtents(z+IM*r); UpdateExtents(z-IM*r);
}
/************************************************************************
                           SET SOLVE DATA 
*************************************************************************
	sets the static solution data. 
  information does not have to be entered all at once
  if only one piece of information is required, send 0 as the other values
  these values are not currently tested for realism!!! TMP DEBUG 
------------------------------------------------------------------------*/
void CMultiLayer::SetSolveData(const int miniter, const int maxiter, const double tolerance){
	if (miniter>0)   {CMultiLayer::MinMLayerIterations=miniter;  }
	if (maxiter>0)   {CMultiLayer::MaxMLayerIterations=maxiter;  }
	if (tolerance>0) {CMultiLayer::MLayerTolerance=    tolerance;}
}
/*************************************************************************
                           ADD TO LAYER
**************************************************************************
	The following functions add elements, far field elements, 
  taylor series and property zones to the layer
------------------------------------------------------------------------*/
//------------------------------------------------------------------------
//	Adds Element to Layer
//------------------------------------------------------------------------
void CMultiLayer::AddToLayer(CmlAnalyticElem *Elemptr, mlelementtype type){

	if (!DynArrayAppend((void**&)(pElemArray),(void*)(Elemptr),nElems)){
		 ExitGracefully("CMultiLayer::AddToLayer(Elem): adding NULL element",BAD_DATA);};
	
	if (type==GIVEN_MLELEM){
		if (!DynArrayAppend((void**&)(pGivenElems),(void*)(Elemptr),nGivenElems)){
		  ExitGracefully("CMultiLayer::AddToLayer(Elem): adding NULL given element",BAD_DATA);};
	}
	else if (type==HEAD_SPECIFIED_MLELEM){
		if (!DynArrayAppend((void**&)(pHSElems),(void*)(Elemptr),nHSElems)){
		  ExitGracefully("CMultiLayer::AddToLayer(Elem): adding NULL head specified element",BAD_DATA);};
	}
} 
//------------------------------------------------------------------------
//	Adds FarField Element to Layer
//------------------------------------------------------------------------
void CMultiLayer::AddToLayer(CmlFarField  *FFptr){
 if (FFptr==NULL){ExitGracefully("CMultiLayer::AddToLayer- adding NULL far field",BAD_DATA);}
 if (pFarField!=NULL){ExitGracefully("CMultiLayer::AddToLayer- Two far field elements added!",BAD_DATA);}
 pFarField=FFptr;
}
//------------------------------------------------------------------------
//	Adds Property Zone to layer
//------------------------------------------------------------------------
void CMultiLayer::AddToLayer(CmlPropZone *Propptr){
  if (Propptr==NULL){
		ExitGracefully("CMultiLayer::AddToLayer- adding NULL property zone",BAD_DATA);
	}
  if (Propptr->GetType()==hydraulic_conductivity){
		if (!DynArrayAppend((void**&)pCondZoneArray,(void*)(Propptr),nCondZones)){
			ExitGracefully("CMultiLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else if (Propptr->GetType()==poro){
		if (!DynArrayAppend((void**&)pPoroZoneArray,(void*)(Propptr),nPoroZones)){
			ExitGracefully("CMultiLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else if (Propptr->GetType()==layer_thickness){
		if (!DynArrayAppend((void**&)pThickZoneArray,(void*)(Propptr),nThickZones)){
			ExitGracefully("CMultiLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
  if (Propptr->GetType()==aq_hydraulic_conductivity){
		if (!DynArrayAppend((void**&)pAqCondZoneArray,(void*)(Propptr),nAqCondZones)){
			ExitGracefully("CMultiLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else if (Propptr->GetType()==aquitard_poro){
		if (!DynArrayAppend((void**&)pAqPoroZoneArray,(void*)(Propptr),nAqPoroZones)){
			ExitGracefully("CMultiLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else if (Propptr->GetType()==aquitard_thickness){
		if (!DynArrayAppend((void**&)pAqThickZoneArray,(void*)(Propptr),nAqThickZones)){
			ExitGracefully("CMultiLayer::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else {ExitGracefully("unused prop type in layer",OTHER);}
}
/************************************************************************
                           GetLevel
*************************************************************************
 Returns the index of the level this point is in
-----------------------------------------------------------------------*/
int	CMultiLayer::GetLevel(const cmplex &z, const double &elev, bool &inaquitard) const{

	double              tmpelev=top_elev;
	Unchangeable1DArray thick;
	Unchangeable1DArray aqthick;
	int                 l;

	thick  =GetThicknessVect(z);
	aqthick=GetAqThicknessVect(z);
	for (l=0;l<nLayers;l++){
		tmpelev-=thick[l];
		if (elev>tmpelev){inaquitard=false; return l;}

		tmpelev-=aqthick[l];
		if (elev>tmpelev){inaquitard=true;  return l;}
	}
	inaquitard=false;
	return 0; //TMP DEBUG
}
/************************************************************************
                           GetAnisotropy
*************************************************************************
 Returns the hydraulic conductivity anisotropy at a given location
-----------------------------------------------------------------------*/
anisotropy CMultiLayer::GetAnisotropy(const cmplex &z, const int lev) const {	
	static anisotropy A;
	A.dir=0;
	A.ratio=0.0;
	return A;
}
/************************************************************************
                           GetCondVect
*************************************************************************
 Returns the hydraulic conductivity vector, K, at a given location
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetCondVect(const cmplex &z) const {	
	return CmlPropZone::NestedSift(z,pCondZoneArray,nCondZones,conductivity);
}
/************************************************************************
                           GetThickVect
*************************************************************************
 Returns the layer thicknesses, T, at a given location
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetThicknessVect(const cmplex &z) const {	
	return CmlPropZone::NestedSift(z,pThickZoneArray,nThickZones,thickness);
}
/************************************************************************
                           GetPoroVect
*************************************************************************
 Returns the porosity, n (or theta) for any point and time
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetPoroVect(const cmplex &z) const{
	return CmlPropZone::NestedSift(z,pPoroZoneArray,nPoroZones,porosity);
}
/************************************************************************
                           GetAqCondVect
*************************************************************************
 Returns the hydraulic conductivity vector, K, at a given location
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetAqCondVect(const cmplex &z) const {	
	return CmlPropZone::NestedSift(z,pAqCondZoneArray,nAqCondZones,conductivity);
}
/************************************************************************
                           GetAqThickVect
*************************************************************************
 Returns the layer thicknesses, T, at a given location
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetAqThicknessVect(const cmplex &z) const {	
	return CmlPropZone::NestedSift(z,pAqThickZoneArray,nAqThickZones,thickness);
}
/************************************************************************
                           GetAqPoroVect
*************************************************************************
 Returns the porosity, n (or theta) for any x,y location
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetAqPoroVect(const cmplex &z) const{
	return CmlPropZone::NestedSift(z,pAqPoroZoneArray,nAqPoroZones,porosity);
}
/************************************************************************
                           GetBase
*************************************************************************
 Returns the layer base elevation, B, for any x,y location
-----------------------------------------------------------------------*/
double CMultiLayer::GetBase(const cmplex &z) const {
	return base_elev; //TMP DEBUG
	//return CmlPropZone::NestedSift(z,pBaseZoneArray,nBaseZones,base_elev);
}
/************************************************************************
                           GetTotalThickness
*************************************************************************
 Returns the total multilayer thickness, Ttot, for any x,y location
-----------------------------------------------------------------------*/
double CMultiLayer::GetTotalThickness 	 (const cmplex &z) const{
	double tot_thick;
	Unchangeable1DArray thick;
	Unchangeable1DArray aqthick;
	thick  =GetThicknessVect(z);
	aqthick=GetAqThicknessVect(z);
  tot_thick=0;
	for (int l=0; l<nLayers-1; l++){
		tot_thick+=thick[l]+aqthick[l];
	}
	tot_thick+=thick[nLayers-1];

	return tot_thick;
}
/************************************************************************
                           GetBottomVect
*************************************************************************
 Returns the layer base elevations, B, for any point and time
 elevations are the top elevation of the aquitard beneath them
-----------------------------------------------------------------------*/
Ironclad1DArray	CMultiLayer::GetBottomVect				 (const cmplex &z) const{
	static double       bottom[MAX_MLAYERS];
	Unchangeable1DArray thick;
	Unchangeable1DArray aqthick;
	int                 l;

	thick  =GetThicknessVect(z);
	aqthick=GetAqThicknessVect(z);

	bottom[0]=top_elev-thick[0];
	for (l=1;l<nLayers;l++){
		bottom[l]=bottom[0]-aqthick[l-1]-thick[l];
	}
	return bottom;
}
/************************************************************************
                           GetConductanceVector
*************************************************************************
 Returns the aquitard conductance vector (size:nLayers-1), for any x,y and time
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetConductVect       (const cmplex &z) const{
	static double       conductance[MAX_MLAYERS-1];
	Unchangeable1DArray thick;
	Unchangeable1DArray cond;
	int l;

  thick=CmlPropZone::NestedSift(z,pAqThickZoneArray,nAqThickZones,thickness);
  cond =CmlPropZone::NestedSift(z,pAqCondZoneArray ,nAqCondZones ,conductivity);
	for (l=1;l<nLayers-1;l++){
		if (thick[l]!=0){
			conductance[l]=cond[l]/thick[l];
		}
		else{
			conductance[l]=0; //TMP DEBUG
		}
	}
	return conductance;
}
/************************************************************************
                           Get Normalized Transmissivity Vector
*************************************************************************
 Returns the Normalized Transmissivity vector (size:nLayers), for any x,y and time
-----------------------------------------------------------------------*/
Ironclad1DArray	CMultiLayer::GetNormTransmissVect (const cmplex &z) const{
	static double trans[MAX_MLAYERS];
	double totaltrans=0;
	Unchangeable1DArray thick;
	Unchangeable1DArray cond;
	int l;

  thick=CmlPropZone::NestedSift(z,pAqThickZoneArray,nAqThickZones,thickness);
  cond =CmlPropZone::NestedSift(z,pAqCondZoneArray ,nAqCondZones ,conductivity);
	for (l=0; l<nLayers; l++){
		trans[l]=cond[l]*thick[l];//TMP DEBUG - correct?
		totaltrans+=trans[l];//TMP DEBUG - correct?
	}
	for (l=0; l<nLayers; l++){
		trans[l]/=totaltrans;
	}
	return trans; 
}
/************************************************************************
                           GetSystemEigenvector
*************************************************************************
 Returns the system eigenvector (size:nLayers), for any x,y and time
-----------------------------------------------------------------------*/
Ironclad1DArray	CMultiLayer::GetSystemEigenvector (const cmplex &z) const{
	static double eigen[MAX_MLAYERS];
	int l;
	for (l=0; l<nLayers; l++){
		eigen[l]=0.0; //TMP DEBUG
	}
	return eigen;
}
/************************************************************************
                           Get Comprehensive Potential 
*************************************************************************
 Returns the discharge potential, Phi, for any point and time
-----------------------------------------------------------------------*/
double CMultiLayer::GetCompPotential(const cmplex &z,const double &t) const {
  static double		Phi;
	int             i;

	Phi=0.0;
	for (i=0; i<nElems; i++){
		Phi+=pElemArray[i]->GetCompPotential(z,t);
	}
	Phi+=pFarField->GetCompPotential(z,t);

	//if (pLayerAbove   !=NULL) {Phi+=pLayerAbove  ->GetCompPotential(BELOW,z,t);}
	//if (pLayerBeneath !=NULL) {Phi+=pLayerBeneath->GetCompPotential(ABOVE,z,t);}
	
  return Phi;
}
/************************************************************************
                           Get Discharge Potential Vector
*************************************************************************
 Returns the discharge potential, Phi, for any point and time
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetPotentialVect(const cmplex &z,const double &t) const {
  static double				Phi[MAX_MLAYERS];
	Unchangeable1DArray ptemp;
	int             i,l;

	for (l=0;l<nLayers;l++){Phi[l]=0.0;}

	for (i=0; i<nElems; i++)
	{
		ptemp=pElemArray[i]->GetPotentialVect(z,t);
		for (l=0;l<nLayers;l++){
			Phi[l]+=ptemp[l];
		}
	}
	ptemp=pFarField->GetPotentialVect(z,t);

	for (l=0;l<nLayers;l++)
	{
		Phi[l]+=ptemp[l];
	} 

	//if (pLayerAbove   !=NULL) {Omega+=pLayerAbove  ->GetDischargePotential(BELOW,z,t);}
	//if (pLayerBeneath !=NULL) {Omega+=pLayerBeneath->GetDischargePotential(ABOVE,z,t);}
	
  return Phi;
}

/************************************************************************
                           GetW
*************************************************************************
 Returns the complex discharge W = Qx - iQy for any point and time
-----------------------------------------------------------------------*/
Ironclad1DArray_z CMultiLayer::GetQxQyVect(const cmplex &z,const double &t) const {
  static cmplex					QxQy[MAX_MLAYERS];
	Unchangeable1DArray_z ptemp;
	int             i,l;

	for (l=0;l<nLayers;l++){QxQy[l]=0.0;}

	for (i=0; i<nElems; i++)
	{
		ptemp=pElemArray[i]->GetQxQyVect(z,t);
		for (l=0;l<nLayers;l++){
			QxQy[l]+=ptemp[l];
		}
	}
	ptemp=pFarField->GetQxQyVect(z,t);

	for (l=0;l<nLayers;l++)
	{
		QxQy[l]+=ptemp[l];
	} 

	//if (pLayerAbove   !=NULL) {Omega+=pLayerAbove  ->GetDischargePotential(BELOW,z,t);}
	//if (pLayerBeneath !=NULL) {Omega+=pLayerBeneath->GetDischargePotential(ABOVE,z,t);}

  return QxQy;
}
/************************************************************************
                           GetGxVector
*************************************************************************
 Returns the discharge derivative Gx = dQx/dx - idQy/dy for any point and time
-----------------------------------------------------------------------*/
Ironclad1DArray_z CMultiLayer::GetGxVect(const cmplex &z,const double &t) const {
  static cmplex					Gx[MAX_MLAYERS];
	Unchangeable1DArray_z ptemp;
	int             i,l;

	for (l=0;l<nLayers;l++){Gx[l]=0.0;}

	for (i=0; i<nElems; i++)
	{
		ptemp=pElemArray[i]->GetGxVect(z,t);
		for (l=0;l<nLayers;l++){
			Gx[l]+=ptemp[l];
		}
	}
	ptemp=pFarField->GetGxVect(z,t);

	for (l=0;l<nLayers;l++)
	{
		Gx[l]+=ptemp[l];
	} 

	//if (pLayerAbove   !=NULL) {Omega+=pLayerAbove  ->GetDischargePotential(BELOW,z,t);}
	//if (pLayerBeneath !=NULL) {Omega+=pLayerBeneath->GetDischargePotential(ABOVE,z,t);}

  return Gx;
}
/************************************************************************
                           GetLeakage
*************************************************************************
 Returns the leakage N any point and time
  if ltype=FROMTOP ,         returns only leakage from above the layer (i.e. recharge, lakes, rivers, confining layer)
  if ltype=FROMBOTTOM,       returns only leakage from below the layer (i.e. through aquitard)
  if ltype=FROMTOPANDBOTTOM, returns net leakage
  leakage is POSITIVE if water is being EXTRACTED from the layer
------------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetLeakageVect(const cmplex &z, const double &t,const leak_type ltype) const {
  static double Leak[MAX_MLAYERS];
	Unchangeable1DArray ltemp;
	int             i,l;

	for (l=0;l<nLayers;l++){Leak[l]=0.0;}

	for (i=0; i<nElems; i++)
	{
		ltemp=pElemArray[i]->GetLeakageVect(z,t,ltype);
		for (l=0;l<nLayers;l++){
			Leak[l]+=ltemp[l];
		}
	}

	//if ((ltype!=FROMBOTTOM) && (pLayerAbove  !=NULL)){Leak-=pLayerAbove->   GetLeakage(z,t);}
	//if ((ltype!=FROMTOP   ) && (pLayerBeneath!=NULL)){Leak+=pLayerBeneath ->GetLeakage(z,t);}

  return Leak;
}
/************************************************************************
                           GetCurl
*************************************************************************
 Returns the Curl (beta) any point and time
------------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetCurlVect(const cmplex &z, const double &t) const {
  static double Curl[MAX_MLAYERS];
	int             l;
	for (l=0;l<nLayers;l++){
		Curl[l]=0.0;
	}
  return Curl;//TMP DEBUG
}
/************************************************************************
                           GetNetDischarge
*************************************************************************
 Returns the net discharge Q from the layer at some time,t 
-----------------------------------------------------------------------*/
double CMultiLayer::GetNetDischarge    (const int lev, const double &t) const{
	static double Q;
	Q=0.0;
	for (int i=0; i<nElems; i++)  {
		Q+=pElemArray[i]->GetNetDischarge(lev,t);
	}
	return Q;
}
//-----------------------------------------------------------------------
double CMultiLayer::GetNetDischarge(const double &t) const{
	static double Q;
	Q=0.0;
	for (int i=0; i<nElems; i++)  {
		Q+=pElemArray[i]->GetNetDischarge(t);
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
cmplex CMultiLayer::GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const int lev, const double &t) const{
	static cmplex Qflx;
	
  Qflx=0.0;
	for (int i=0; i<nElems; i++){
    Qflx+=pElemArray[i]->GetFluxThruFace(z1,z2,lev,t);
	}

  Qflx+=pFarField->GetFluxThruFace(z1,z2,lev,t);

	return Qflx;
}
/************************************************************************
                           Integrated Leakage
*************************************************************************
 Returns the integrated leakage [L^3/T] over a triangle defined by z1,z2,z3
 Returns positive value for leakage, negative for recharge
------------------------------------------------------------------------*/
double CMultiLayer::GetIntegratedLeakage(const cmplex    &z1, 
																				 const cmplex    &z2, 
																				 const cmplex    &z3, 
																				 const int        lev,
																				 const double    &t,
																				 const leak_type  ltype) const {
	static double Q;
	Q=0.0;
	for (int i=0; i<nElems; i++)  {
		Q+=pElemArray[i]->GetIntegratedLeakage(z1,z2,z3,lev,t,ltype);
	}
	return Q;
}
/************************************************************************
                           Integrated Leakage
*************************************************************************
 Returns the integrated budget (net influx and outflux) over a triangle defined by z1,z2,z3
------------------------------------------------------------------------*/
void CMultiLayer::GetIntegratedBudget(const cmplex &z1, 
																 const cmplex &z2, 
																 const cmplex &z3, 
																 const int     lev,
																 const double &t, 
																       double &inQ, 
																			 double &outQ) const{
  double tmpin,tmpout;
	inQ=outQ=0.0;
	for (int i=0; i<nElems; i++)  {
		pElemArray[i]->GetIntegratedBudget(z1,z2,z3,lev,t,tmpin,tmpout);
		inQ+=tmpin;
		outQ+=tmpout;
	}
}
/************************************************************************
                           Get Flux Distribution 
*************************************************************************
	Gets integrated extraction between points z1 and z2 along line
	Linearly Distributes these fluxes to endpoints Q1 and Q2
	Used for finite element source distribution
------------------------------------------------------------------------*/
void CMultiLayer::GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const int lev, const double &t, double &Q1, double &Q2) const{

  double tmpQ1,tmpQ2;
	Q1=Q2=0.0;
	for (int i=0; i<nElems; i++)  {
		pElemArray[i]->GetFluxDistribution(z1,z2,lev,t,tmpQ1,tmpQ2);
		Q1+=tmpQ1;
		Q2+=tmpQ2;
	}
}
/************************************************************************
                           Horizontal Seepage Velocity
*************************************************************************
 Returns the horizontal seepage velocity (vx, vy) at any complex point and time
-----------------------------------------------------------------------*/
cmplex CMultiLayer::GetVelocity2D(const cmplex &z, const int lev, const double &t) const{
  static cmplex W;
  static double K,T,B,n,head,pot;

  T   =GetThicknessVect(z)[lev];
  n   =GetPoroVect	   (z)[lev];
  W   =GetQxQyVect   (z,t)[lev];                       

	if ((abs(W)>=ALMOST_INF) || (T<=0.0) || (n<=0.0)){return 0.0;}

	return conj(W)/(n*T);

}

/************************************************************************
                           Seepage Velocity (3D)
*************************************************************************
 Returns the seepage velocity (vx, vy, vz) at any 3D point and time
 calculates vertical and horizontal velocity component
 MUST OPTIMIZE
------------------------------------------------------------------------*/
vector CMultiLayer::GetVelocity3D(const pt3D &pt, const double &t) const{
	
	static cmplex W,v,z;           //parameters declared static to reduce computational cost of call
	static double pot,head;
	static double leaktop,leakbot;
	static double B,K,T,n,v_z;

	int lev;
	bool   inaquitard;

	z      =c3Dto2D(pt);

	lev    =GetLevel(z,pt.z,inaquitard);

	if (!inaquitard){
		T   =GetThicknessVect(z)[lev];
		n   =GetPoroVect		 (z)[lev];
		W   =GetQxQyVect		 (z,t)[lev];

		leaktop=GetLeakageVect       (z,t,FROMTOP   )[lev]; //leakup & leakbot are negative if sources of water
		leakbot=GetLeakageVect       (z,t,FROMBOTTOM)[lev];

		head   =GetHeadVect(z,t)[lev];

		if ((abs(W)>=ALMOST_INF) || (T<=0.0) || (n<=0.0) || (head<=0.0)){return 0.0;}

		v_z=ConvertToVertVelocity(pt.z,K,T,B,n,head,pot,W,leaktop,leakbot);
	}
	else{
		n  =GetAqPoroVect(z)[lev];
		v_z=leakbot/n;
	}
	
	//if ((pt.z<B) || (pt.z>(T+B)) || (pt.z>(head+B))) {return 0.0;}

	//if (IsConfined(pot,K,T)){v=conj(W/(n*T   ));}
	//else                    {v=conj(W/(n*head));}
  v=0.0; //TMP DEBUG
	return vector(v.real(),v.imag(),v_z);

}
/*************************************************************************
                           Vertical Seepage Velocity
**************************************************************************
 Returns the vertical seepage velocity (vz) at any complex point and time
------------------------------------------------------------------------*/
double CMultiLayer::ConvertToVertVelocity(const double &z,const double &K,const double &T, 
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
/************************************************************************
                           GetHead
*************************************************************************
 Returns the GLOBAL head value at any complex point and time
-----------------------------------------------------------------------*/
Ironclad1DArray CMultiLayer::GetHeadVect(const cmplex &z, const double &t) const {
	static double head[MAX_MLAYERS];
	static Unchangeable1DArray K,T;
  int l;

  K=GetCondVect(z);	
	T=GetThicknessVect(z);

	for (l=0;l<nLayers;l++){
		head[l]=0.0;
	}
	//head=ConvertToHead(GetDischargePotential(z,t).real(),K,T,sealevel-B,saltwaterSG);
  //if (head<=0.0){return B;}
  //else          {return head+B;}
	return head; //TMP DEBUG
}
double CMultiLayer::GetHead(const pt3D &pt, const double &t) const{
	bool inaq;
	cmplex z=c3Dto2D(pt);
	int l=GetLevel(z,pt.z,inaq);
	if (inaq){
		return (GetHeadVect(z,t)[l]+GetHeadVect(z,t)[l+1])*0.5; //TMP DEBUG
	}
	else{
		return GetHeadVect(z,t)[l];
	}
}
/************************************************************************
                           IterativeSolve
*************************************************************************
 Iteratively solves for the coefficients of each element in the layer
	 maxchange is the maximum change of any element during the current iteration
	 maxobj    is the maximum error of any element during the current iteration
	 PROGRESS  is the output progress file
------------------------------------------------------------------------*/
void CMultiLayer::IterativeSolve(double &maxchange, double &maxobj, const double &t, ofstream &PROGRESS){

	double					 tempchange,tempobj;
  int							 i;
	CmlAnalyticElem *WorstChange;
	CmlAnalyticElem *WorstObj;
	double					 AvgObj;
	bool						 stopped(false);
	ofstream				 ANALYSIS;

  cout << "Solving Multilayer system ("<<nElems << " element(s))..."<<endl;
	cout << "------------------------------------------------------------"<<endl;

  PROGRESS << "solve"<<endl;

	CMultiLayer::solved=false;
	
	iter=1; 

	//Solve given elements first
	for (i=0; i<nGivenElems; i++){
		pGivenElems[i]->SolveItself(tempchange,tempobj,0.0);
	}

	while ((!stopped) &&
		     (((iter<=CMultiLayer::MaxMLayerIterations) && (maxchange>CMultiLayer::MLayerTolerance)) ||
				   (iter<=CMultiLayer::MinMLayerIterations)) ) { 
		
		tempchange =tempobj=0.0;	
		maxchange  =maxobj =0.0;
		WorstChange=NULL;
		WorstObj   =NULL;

		if (!stopped){

			AvgObj=0.0;
      CalculateDeltaPotential(t);

			pFarField->SolveItself(tempchange,tempobj,t);
			upperswap(maxobj   ,tempobj   );if (tempobj   ==maxobj   ) {WorstObj   =pFarField;}
      upperswap(maxchange,tempchange);if (tempchange==maxchange) {WorstChange=pFarField;}  
			AvgObj+=tempobj;


			for (i=0; i<nElems; i++){
				tempchange=tempobj=0.0;			

				if (ProgramAborted())  { stopped=true;i=nElems;}

				if (!stopped) {

					WriteEllipsisToScreen(i,nElems,40);

					if (!pElemArray[i]->IsGiven()){ 
						pElemArray[i]->SolveItself(tempchange,tempobj,t);
						upperswap(maxobj   ,tempobj   );if (tempobj   ==maxobj   ) {WorstObj   =pElemArray[i];}
						upperswap(maxchange,tempchange);if (tempchange==maxchange) {WorstChange=pElemArray[i];}  
						AvgObj+=tempobj;
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
			
			/*if (iter%WriteInterval==0){
				cout << "Writing Back-up Solution..."<<endl;
				ofstream SOL;
				SOL.open("solution.bbs");
				SOL.precision(10);
				WriteItself(SOL,t);
				SOL.close();
			}*/

			iter++;
			//CMultiLayer::fresh=false;
		}
  }
	CMultiLayer    ::solved=true;

	if (iter<CMultiLayer::MaxMLayerIterations)  {cout << "...success!"                         <<endl;}
	else if (stopped)                     {cout << "...STOP file generated. Please wait."<<endl;}
  else                                  {cout << "...exceeded maximum layer iterations"<<endl;}
	PROGRESS<<"done"<<endl;
}
/************************************************************************
                           IdentifyDeltaPotential
*************************************************************************
  Takes 16 points across domain to estimate potential range (max-min)
-----------------------------------------------------------------------*/
void CMultiLayer::CalculateDeltaPotential(const double t){
	
	double minpot(ALMOST_INF), maxpot(-ALMOST_INF),pot(0.0);
	if ((Extents.e!=Extents.w) && (Extents.n!=Extents.s)){
		for   (double x=Extents.w; x<=Extents.e+(Extents.e-Extents.w)*0.1; x+=(Extents.e-Extents.w)/4.0){
			for (double y=Extents.s; y<=Extents.s+(Extents.n-Extents.s)*0.1; y+=(Extents.n-Extents.s)/4.0){
				pot=GetCompPotential(cmplex(x,y),t);
				if ((pot<minpot) && (pot >-ALMOST_INF)){minpot=pot;}
				if ((pot>maxpot) && (pot <ALMOST_INF) ){maxpot=pot;}
			}
		}
	}
	DeltaPot=maxpot-minpot;
	if (DeltaPot<1.0){DeltaPot=1.0;}

}
/************************************************************************
                           Write Itself/Write Output
*************************************************************************
  CMultiLayer::WriteItself: Writes solution file
  CMultiLayer::WriteOutput: Writes generic output
-----------------------------------------------------------------------*/
void CMultiLayer::WriteItself(ofstream &SOL, const double &t) const{
	pFarField->WriteItself(SOL,t);
	for (int i=0; i<nElems; i++){				
    pElemArray[i]->WriteItself(SOL,t);
	}
}
//-----------------------------------------------------------------------
void CMultiLayer::WriteOutput(const double &t) const{
	bool stopped(false);
	pFarField->WriteOutput(t);
	for (int i=0; i<nElems; i++){
		if (ProgramAborted())  {stopped=true;i=nElems;}
		if (!stopped){pElemArray[i]->WriteOutput(t);}
	}
}

