//CAquiclude.cpp
#include "Aquitard.h"

//*************************************************************************
//                             CAquiclude
//************************************************************************
//                           CONSTRUCTORS
//************************************************************************
CAquiclude::CAquiclude(){ 
}
//------------------------------------------------------------------------
CAquiclude::~CAquiclude(){
	if (globaldebug){cout <<" DESTROYING AQUICLUDE "<<endl;}
  int i;
	for (i=0; i<nCondZones; i++)	{delete pCondZoneArray [i];} delete [] pCondZoneArray;
	delete [] pElemArray; //just pointer, elements deleted globally
}
//------------------------------------------------------------------------
CAquiclude::CAquiclude(CLayerSetABC	  *pLaybelow,
										   CLayerSetABC		*pLayabove, 
										   const double		 elev, 
										   const int			 lev){

	pLayerBelow=pLaybelow;  
	pLayerAbove=pLayabove;  
  nCondZones=0;
  nElems=0;
  level=lev;

	pCondZoneArray=NULL;
	pElemArray=NULL;

	conductance=0.0; //default value -should always be background
	elevation=elev;
	DeltaLeakage=1.0;

	Extents.n=-ALMOST_INF;
	Extents.e=-ALMOST_INF;
	Extents.w= ALMOST_INF;
	Extents.s= ALMOST_INF;
}
//************************************************************************
//                  STATIC INITIALIZERS (default values)
//************************************************************************
int    CAquiclude::MinAquicludeIterations=1; 
int    CAquiclude::MaxAquicludeIterations=300;
double CAquiclude::AquicludeTolerance    =0.00001;
//************************************************************************
//                           ACCESSOR FUNCTIONS
//************************************************************************
double CAquiclude::GetConductance(const cmplex &z) const{
	//this implementation may be silly- nested conductance inhoms will not be effective
	return CPropZone::NestedSift(z,pCondZoneArray,nCondZones,conductance);
}
//------------------------------------------------------------------------
double				CAquiclude::GetElevation (const cmplex &z) const{return elevation;  }
//------------------------------------------------------------------------
CLayerSetABC *CAquiclude::GetLayerAbove() const{return pLayerAbove;}
//------------------------------------------------------------------------
CLayerSetABC *CAquiclude::GetLayerBelow() const{return pLayerBelow;}
//------------------------------------------------------------------------
double				CAquiclude::GetDeltaLeakage() const{return DeltaLeakage;}
//************************************************************************
//                           ASSIGNMENT FUNCTIONS
//************************************************************************
void CAquiclude::UpdateExtents (const cmplex &z){
	upperswap(Extents.e,z.real()); lowerswap(Extents.w,z.real());
	upperswap(Extents.n,z.imag()); lowerswap(Extents.s,z.imag());
}
//------------------------------------------------------------------------
void CAquiclude::UpdateExtents (const cmplex &z,const double &r){
	upperswap(Extents.e,z.real()+r); lowerswap(Extents.w,z.real()-r);
	upperswap(Extents.n,z.imag()+r); lowerswap(Extents.s,z.imag()-r);
}
//************************************************************************
//                           AddToLayer
//************************************************************************
void CAquiclude::AddToLayer(CAnalyticElem *Elemptr){
	if (!DynArrayAppend((void**&)(pElemArray),(void*)(Elemptr),nElems)){
	  ExitGracefully("CAquiclude::AddToLayer(Elem): adding NULL element",BAD_DATA);};
} 
//-----------------------------------------------------------------------
void CAquiclude::AddToLayer(CPropZone *Propptr){
	if (Propptr==NULL)          {ExitGracefully("CAquiclude::AddToLayer- adding NULL prop zone",        BAD_DATA);}
  if (Propptr->GetType()==hydraulic_conductance){
		if (!DynArrayAppend((void**&)pCondZoneArray,(void*)(Propptr),nCondZones)){
			                         ExitGracefully("CAquiclude::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
	}
	else {
		ExitGracefully("Unused property type in Aquiclude",OTHER);}
}
//************************************************************************
//                           MEMBER FUNCTIONS
//************************************************************************
//************************************************************************
//                           GetDischargePotential
//************************************************************************
cmplex CAquiclude::GetDischargePotential(dir direct, const cmplex &z,const double &t) const{
  cmplex Omega(0.0,0.0);
	for (int i=0; i<nElems; i++){
	  Omega+=pElemArray[i]->GetDischargePotential(z,t); 
		//areasink-negative pot for positive leakage (above)
		//areasink-positive pot for positive leakage (below)
		//areasink-negative pot for positive leakage (if in layer)
	}
	if (direct==BELOW){Omega=-Omega;}

  return Omega;
}

//************************************************************************
//                           GetW
//************************************************************************
cmplex CAquiclude::GetW(dir direct, const cmplex &z,const double &t) const {
  cmplex W(0.0,0.0);
	for (int i=0; i<nElems; i++)    {
		W+=pElemArray[i]->GetW(z,t);
	}
	if (direct==BELOW){W=-W;}
  return W;
}
//************************************************************************
double CAquiclude::GetLeakage           (const cmplex &z,const double &t) const{
	//leakage is positive if water moving from upper to lower
  double Leak(0.0);
	for (int i=0; i<nElems; i++)    {
		Leak+=pElemArray[i]->GetLeakage(z,t,FROMTOPANDBOTTOM);
	}
  return Leak;
}
//************************************************************************
double CAquiclude::GetDesiredLeakage    (const cmplex &z,const double &t) const{
	double headup,headdown,cond;
	if (pLayerAbove==NULL){headup=  0.0;                      }
	else                  {headup=  0.0;}//pLayerAbove->GetHead(z,t);}//TMP DEBUG- ML MIGRATION
	if (pLayerBelow==NULL){headdown=0.0;                      }
	else                  {headdown=0.0;}//pLayerBelow->GetHead(z,t);}//TMP DEBUG- ML MIGRATION

	cond=GetConductance(z);
	//leakage from global head difference (actual leakage)
	//leakage is positive if water moving from upper to lower
	if   (headdown>=elevation){																									//bottom confined
		if (headup  >=elevation){	return cond*(headup   -headdown );} //top wet
		else                    {	return cond*(elevation-headdown );}}//top dry 
	else {																																			//bottom unconfined
		if (headup  >=elevation){	return cond*(headup   -elevation);} //top wet
		else										{	return 0.0;                       } //top dry	
	}
}
//************************************************************************
double CAquiclude::GetNetFlux (const double &t) const{
return 0.0;
}
//************************************************************************
void CAquiclude::SolveItself          (double &change, double &objective,const double t){
	double tempchange(0.0),maxchange(0.0);
	double tempobj(0.0),maxobj(0.0);
	int    iter(1);
	bool   stopped(false);
	cout << "Solving Aquiclude #"<<level+1<<" ("<<nElems<<" elements)..."<<endl;
	cout << "------------------------------------------------------------"<<endl;

	if (nElems>0){
	maxchange=0;
	while ((iter<MaxAquicludeIterations) && 
				 (maxchange>AquicludeTolerance) && 
				 (!stopped)){ 
		maxchange=0.0;
		maxobj=0.0;
		for(int i=0;i<nElems;i++){
			if (ProgramAborted()){stopped=true;i=nElems;}

			if (!stopped){
				pElemArray[i]->SolveItself(tempchange,tempobj,t);
				upperswap(maxobj,   tempobj   );
				upperswap(maxchange,tempchange);
			}
		}
		cout << "iteration "      << iter
				 << ",   maxchange: " << maxchange
				 << ",    max. obj: " << maxobj<<endl;
		iter++;
	}
	}

	change=maxchange;
	objective=maxobj;
}
//************************************************************************
void CAquiclude::WriteItself          (ofstream &SOL, const double &t) const{
	for (int i=0; i<nElems; i++){				
    pElemArray[i]->WriteItself(SOL,t);
	}
}
//************************************************************************
void CAquiclude::WriteOutput          (const double &t) const{
	for (int i=0; i<nElems; i++){				
    pElemArray[i]->WriteOutput(t);
	}

}
