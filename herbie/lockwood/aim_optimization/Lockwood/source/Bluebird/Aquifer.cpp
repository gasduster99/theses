#include "Aquifer.h"
/************************************************************************
                             CAquifer
*************************************************************************
                           CONSTRUCTORS
************************************************************************/
CAquifer::CAquifer(char * aqname){

 name=new char[strlen(aqname)+1]; 
 name=strcpy(name,aqname);
 
 nLevels=0;
 
	for (int l=0; l<MAXAQLEVELS; l++){
		pSLayerArray[l]			=NULL; //initially empty, layers added with AddSingleLayerToAquifer or AddMultiLayerToAquifer
		pMLayerArray[l]			=NULL; //initially empty, layers added with AddSingleLayerToAquifer or AddMultiLayerToAquifer
		LevelIsMultiLayer[l]=false;
		pAquicludeArray[l]	=NULL;
	}

}

//------------------------------------------------------------------------
CAquifer::CAquifer(char * aqname, CSingleLayer **pLayers, CAquiclude **pAquicludes, int numlevels){
 //TMP DEBUG: ML MIGRATION: Should remove entirely 

	name = new char[strlen(aqname)+1]; 
	name=strcpy(name,aqname);
 
	ExitGracefullyIf(pLayers==NULL, "CAquifer::Constructor: NULL Layer ",BAD_DATA);
	ExitGracefullyIf(pAquicludes==NULL, "CAquifer::Constructor: NULL Aquiclude",BAD_DATA);
	nLevels=numlevels;
  int L;
	for (L=0; L<nLevels; L++){
		pSLayerArray     [L]=pLayers  [L]; 
		pMLayerArray     [L]=NULL;
		LevelIsMultiLayer[L]=false;
		pAquicludeArray  [L]=pAquicludes[L];
	}
	for (L=nLevels; L<MAXAQLEVELS;L++){
		pSLayerArray		 [L]=NULL;
		pMLayerArray		 [L]=NULL;
		LevelIsMultiLayer[L]=false;
		pAquicludeArray	 [L]=NULL;
	}
}
//------------------------------------------------------------------------
CAquifer::~CAquifer(){
	if (globaldebug){cout <<"DESTROYING AQUIFER"<<endl;}
	//destroys levels, aquicludes
	for (int L=0;L<nLevels;L++){
		delete pSLayerArray   [L];
		delete pMLayerArray   [L];
		delete pAquicludeArray[L];
	}
	delete [] name;
}
/************************************************************************
                           STATIC INITIALIZERS
************************************************************************/
int    CAquifer::MinAqIterations=1; 
int    CAquifer::MaxAqIterations=50;
double CAquifer::AqTolerance    =0.00001;
bool   CAquifer::Multilevel     =false;
/************************************************************************
                           ACCESSOR FUNCTIONS
************************************************************************/
//CSingleLayerABC  *CAquifer::GetLayer    (int L) const {return (CSingleLayerABC*)(pSLayerArray[L]);}//TMP DEBUG-ML MIGRATION
//CLayerSetABC     *CAquifer::GetLayer    (int L) const {return (CLayerSetABC)(pLayerArray[L]);}
//------------------------------------------------------------------------
CAquiclude		*CAquifer::GetAquiclude(int L) const {return pAquicludeArray[L];} 
//------------------------------------------------------------------------
int						 CAquifer::GetNumLevels()      const {return nLevels;}
//------------------------------------------------------------------------
int						 CAquifer::GetNumLayers()      const {
  int sum=0;
	for (int l=0; l<nLevels; l++){
		if (LevelIsMultiLayer[l]){
			sum+=pMLayerArray[l]->GetNumLayers();
		}
		else {
			sum++;
		}
	}
	return sum;
}
//------------------------------------------------------------------------
CSingleLayer *CAquifer::GetLayer(const int L)      const {return pSLayerArray[0];}//TMP DEBUG-ML MIGRATION
//CLayerSet    **CAquifer::GetAllLayerSets()   const {return pSLayerArray;}
//------------------------------------------------------------------------
int        CAquifer::GetLevel    (const pt3D pt) const{
	if (nLevels==1){return 0;} //quick out for single level models
	cmplex z(pt.x,pt.y);
	int L=0; //top level is default
	while ((L<nLevels-1) && 
		     (pt.z<(pSLayerArray[L]->GetThick(z)+pSLayerArray[L]->GetBase(z)))){
		L++;
	}
	//if (pt.z>(pLayerArray[L]->GetThick(z)+pLayerArray[L]->GetBase(z))){return -1;}
	return L;
}
/************************************************************************
                           ASSIGNMENT FUNCTIONS
************************************************************************/
void CAquifer::SetSolveData(const int miniter, const int maxiter, const double tolerance){
	if (miniter!=0)  {CAquifer::MinAqIterations=miniter;  }
	if (maxiter!=0)  {CAquifer::MaxAqIterations=maxiter;  }
	if (tolerance>0) {CAquifer::AqTolerance=    tolerance;}
}
//------------------------------------------------------------------------
/*void CAquifer::AddLayerSetToAquifer(CLayerSet *pLayer){
	if (!DynArrayAppend((void**&)(pLayerArray),(void*)(pLayer),nLevels)){
	  ExitGracefully("CAquifer::AddLayerSetToAquifer(LayerSet): adding NULL Layer Set",BAD_DATA);};
}*/
//------------------------------------------------------------------------
void CAquifer::AddSingleLayer(CSingleLayer *pLayer, CAquiclude *pAquicludeBeneath){
	ExitGracefullyIf(pLayer==NULL,"CAquifer::AddSingleLayer: NULL Single Layer Added", BAD_DATA);
	ExitGracefullyIf(nLevels>=MAXAQLEVELS-1,"CAquifer::AddSingleLayer: Too many levels added to aquifer", BAD_DATA);
	
	pSLayerArray     [nLevels]=pLayer;
  pMLayerArray     [nLevels]=NULL;
	pAquicludeArray  [nLevels]=pAquicludeBeneath;
	LevelIsMultiLayer[nLevels]=false;
	nLevels++;
}
//------------------------------------------------------------------------
void CAquifer::AddMultiLayer(CMultiLayer  *pLayer, CAquiclude *pAquicludeBeneath){
	ExitGracefullyIf(pLayer==NULL,"CAquifer::AddSingleLayer: NULL MultiLayer Added", BAD_DATA);
	ExitGracefullyIf(nLevels>=MAXAQLEVELS-1,"CAquifer::AddSingleLayer: Too many levels added to aquifer", BAD_DATA);

	pSLayerArray     [nLevels]=NULL;
  pMLayerArray     [nLevels]=pLayer;
	pAquicludeArray  [nLevels]=pAquicludeBeneath;
	LevelIsMultiLayer[nLevels]=true;
	nLevels++;
}

/************************************************************************
                           IterativeSolve
************************************************************************/
void CAquifer::IterativeSolve(double &t,ofstream &PROGRESS){
	
	double maxchange(0.0),maxobjective(0.0),tempchange,tempobj;
  int    iter(1),L;
	bool   stopped(false);
  cout << "Solving Aquifer"<<endl;
	cout << "------------------------------------------------------------"<<endl;

  while (
		     (((iter<=CAquifer::MaxAqIterations) && (maxchange>CAquifer::AqTolerance)) ||
				  (iter<=CAquifer::MinAqIterations)) && 
				 (!stopped)
		    ){ 

		if (ProgramAborted()){stopped=true;}

		if (!stopped){
			//Solve Layers (may have to modify layer criteria on last iteration)
			for (L=0; L<nLevels; L++){
				tempchange=tempobj=0.0;
				if (LevelIsMultiLayer[L]){
					pMLayerArray[L]->IterativeSolve(tempchange,tempobj,t,PROGRESS);
				}
				else{
					pSLayerArray[L]->IterativeSolve(tempchange,tempobj,t,PROGRESS);
				}
				upperswap(maxchange   ,tempchange);
				upperswap(maxobjective,tempobj   );
			} 
			//Solve Leaky Layers
			if (Multilevel){
				cout <<"Solving for leakage between aquifer levels"<<endl;
				cout << "------------------------------------------------------------"<<endl;
				for (L=0; L<nLevels; L++){
					tempchange=tempobj=0.0;
					pAquicludeArray[L]->SolveItself(tempchange,tempobj,t);
					upperswap(maxchange   ,tempchange);
					upperswap(maxobjective,tempobj   );
				} 
			}	
			else {
				CAquifer::MaxAqIterations=1;
				CAquifer::MinAqIterations=0;
			}
		}
    cout <<endl;
		cout << "==============================================================="<<endl;
    cout << "Aquifer iteration "      << iter
				 << ",   maxchange: " << maxchange
				 << ",    max. obj: " << maxobjective << endl;
    cout <<"==============================================================="<<endl;
		
    iter++;
  }
  if ((iter<CAquifer::MaxAqIterations)){cout << "...success!" <<endl;}
  else if (!stopped)                   {cout << "...exceeded maximum iterations"<<endl;}
	else                                 {cout << "...stop file generated"<<endl;}

}
//----------------------------------------------------------------------
void CAquifer::IterExplicitSolve(double &t, ofstream &PROGRESS){
	ExitGracefully("CAquifer::IterExplicitSolve: Stub Function",BAD_DATA);
	return;
}
/************************************************************************
                           GetHead
************************************************************************/
double CAquifer::GetHead              (const pt3D &pt, const double &t) const{
	cmplex z(pt.x,pt.y);  
	int L(GetLevel(pt));
	if ((L>=0) && (L<nLevels)){
    if (LevelIsMultiLayer[L]){return pMLayerArray[L]->GetHead(pt,t);}
		else										 {return pSLayerArray[L]->GetHead(z,t); }
	}
	else                      {return 0.0;}
}

/************************************************************************
                           GetNetDischarge
************************************************************************/
double CAquifer::GetNetDischarge(const double &t) const    {
	double Q(0);
	for (int L=0; L<nLevels; L++){
		if (LevelIsMultiLayer[L]){Q+=pMLayerArray[L]->GetNetDischarge(t);}
		else										 {Q+=pSLayerArray[L]->GetNetDischarge(t);}
	}
	return Q;
}
/************************************************************************
                           Get Seepage Velocity (2D)
************************************************************************/
vector CAquifer::GetVelocity2D(const pt3D &pt,const double &t) const{
	cmplex z=c3Dto2D(pt);
	int L=GetLevel(pt);
	//TMP DEBUG - should handle aquiclude velocity as well
  if ((L>=0) && (L<nLevels)){
		if (LevelIsMultiLayer[L]){ExitGracefully("CAquifer::GetVelocity2D:STUB",BAD_DATA);return 0;}
		else										 {return c2Dto3D(pSLayerArray[L]->GetVelocity2D(z,t));}
	}
  else  {return 0.0; }
} 
/************************************************************************
                           Get Seepage Velocity (3D)
************************************************************************/
vector CAquifer::GetVelocity3D(const pt3D &pt,const double &t) const{
	int L=GetLevel(pt);
	//TMP DEBUG - should handle aquiclude velocity as well
  if ((L>=0) && (L<nLevels)){return pSLayerArray[L]->GetVelocity3D(pt,t);}
  else                      {return 0.0;}
}
/************************************************************************
                           Get Effective Velocity (2D)
************************************************************************/
vector CAquifer::GetEffVelocity2D(const pt3D &pt,  const double &t, const disp_info &disp) const{
	int L=GetLevel(pt);
	//TMP DEBUG - should handle aquiclude velocity as well
  if ((L>=0) && (L<nLevels)){return c2Dto3D(pSLayerArray[L]->GetEffVelocity2D(c3Dto2D(pt),t,disp));}
  else                      {return 0.0;}
}
/************************************************************************
                           Get Effective Velocity (3D)
************************************************************************/
vector CAquifer::GetEffVelocity3D(const pt3D &pt,  const double &t, const disp_info &disp) const{
	int L=GetLevel(pt);
	//TMP DEBUG - should handle aquiclude velocity as well
  if ((L>=0) && (L<nLevels)){return pSLayerArray[L]->GetEffVelocity3D(pt,t,disp);}
  else                      {return 0.0;}
}
/************************************************************************
                           GetConductivity
************************************************************************/
double     CAquifer::GetCond(const cmplex &z, int L) const {
	if ((L>=0) && (L<nLevels)){return pSLayerArray[L]->GetCond(z);}
	else											{return ALMOST_INF;}
}
//------------------------------------------------------------------------
double CAquifer::GetCond              (const pt3D &pt) const{
	cmplex z(pt.x,pt.y);  
	return GetCond(z,GetLevel(pt));
}
/************************************************************************
                           GetPorosity
************************************************************************/
double CAquifer::GetPoro(const cmplex &z,const int L) const{
  if ((L>=0) && (L<nLevels)){return pSLayerArray[L]->GetPoro(z);}
	else                      {return 1.0;}
}
//------------------------------------------------------------------------
double CAquifer::GetPoro              (const pt3D &pt) const{
	cmplex z(pt.x,pt.y);  
	return GetPoro(z,GetLevel(pt));
}
/************************************************************************
                           GetSaturatedThickness
************************************************************************/
double CAquifer::GetSaturatedThickness(const cmplex &z,const int L, const double &t) const{
  if ((L>=0) && (L<nLevels)){return pSLayerArray[L]->GetSaturatedThickness(z,t);}
	else                      {return 0.0;}
}
//------------------------------------------------------------------------
double CAquifer::GetSaturatedThickness(const pt3D &pt,const double &t) const{
	cmplex z(pt.x,pt.y);  
	return GetSaturatedThickness(z,GetLevel(pt),t);
}
/************************************************************************
                           Get Leakage
************************************************************************/
double CAquifer::GetLeakage (const pt3D &pt,const double &t,const leak_type ltype) const{
	int L=GetLevel(pt);
	cmplex z(pt.x,pt.y); 
  if ((L>=0) && (L<nLevels)){return pSLayerArray[L]->GetLeakage(z,t,ltype);}
	else                      {return 0.0;}
}
/*double CAquifer::GetLeakage(const cmplex &z,const int L,const double &t) const{
  if ((L>=0) && (L<nLevels)){return pAquicludeArray[L]->GetLeakage(z,t);}
	else                      {return 0.0;}
}*/
/************************************************************************
                           Get Baseflow
************************************************************************/
double CAquifer::GetBaseflow(const cmplex &z, int L, const double t) const{
	if ((L>=0) && (L<nLevels)){return pSLayerArray[L]->GetBaseflow(z,t);}
	else											{return 0.0;}
}
//------------------------------------------------------------------------
double CAquifer::GetBaseflow(const pt3D &pt, const double t) const{
	cmplex z(pt.x,pt.y);  
	return GetBaseflow(z,GetLevel(pt),t);
}
/************************************************************************
                           InterpolateQxQy
************************************************************************/
void   CAquifer::InterpolateQxQy      (const double &t){
int L;
	for (L=0; L<nLevels; L++){
		 pSLayerArray[L]->InterpolateQxQy(t);
	}
}
/************************************************************************
                           WriteItself
*************************************************************************
Calls all solution creation functions of contained layers & elements
-----------------------------------------------------------------------*/
void CAquifer::WriteItself(ofstream &SOL, const double &t) const{
	int L;
	for (L=0; L<nLevels; L++){
		if (LevelIsMultiLayer[L]){
			pMLayerArray[L]->WriteItself(SOL,t);
		}
		else{
			pSLayerArray[L]->WriteItself(SOL,t);
		}
		if (pAquicludeArray[L]!=NULL){
		   pAquicludeArray[L]->WriteItself(SOL,t);  
		}
	}
}
//************************************************************************
void CAquifer::WriteOutput(const double &t) const{
	int L;
	for (L=0; L<nLevels; L++){
		if (LevelIsMultiLayer[L]){
			pMLayerArray[L]->WriteOutput(t);
		}
		else{
			pSLayerArray[L]->WriteOutput(t);
		}	
		if (pAquicludeArray[L]!=NULL){
		  pAquicludeArray[L]->WriteOutput(t);  
		}
	}
}






