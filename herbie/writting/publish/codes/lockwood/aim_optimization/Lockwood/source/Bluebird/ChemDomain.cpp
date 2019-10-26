//CChemDomain
#include "ChemDomain.h"

/*********************************************************************
                       GENERIC TRANSPORT DOMAIN
**********************************************************************
                           CONSTRUCTORS
**********************************************************************/
CChemDomain::CChemDomain():
             C2DDomainABC(){
  												
	//Domain always totally empty at creation       
  pAq=NULL;                             //Initialized in SetAquifer

	poro         =NULL;                   //Initialized in InitializeGrid

	pSpeciesArray=NULL;	nspecies=0;       //initialized/dynamically modified in AddToDomain Routines      
	pSoilArray   =NULL; nSoils=0;

	StartTime    =0.0; //default
	OutputTimes  =NULL;	nOutputTimes=0;WriteInterval=1;
	ObsPoints    =NULL; nObsPoints=0;

	pConcGrid    =NULL;                   //initialized in SetGrid	

	pTransScheme =NULL;                   //initialized in SetTransportScheme
	pDispScheme  =NULL;										//initialized in SetDispersionScheme
  pReactSchemes=NULL;                   //initialized in SetReactionScheme
	nRxNs=0;
	
	C            =NULL;										//initialized in InitializeGrid
	Ctmp         =NULL;										//initialized in InitializeGrid 
	Cnew         =NULL;										//initialized in InitializeGrid	
	Csorb        =NULL;                   //initialized in InitializeGrid	 

	IsothermTable=NULL;                   //initialized in InitializeIsothermTable
	for (int s=0;s<MAX_SPECIES;s++){
		sorbing     [s]=false;
		UsesIsotherm[s]=false;
    Cback       [s]=0.0;
	}

	NodalSourceFluxes=NULL;               //initialized in InitializeSources
	pNodalSources    =NULL;          
	nNodalSources    =NULL;               
	
	srand(2);

	TimeStep				=0.1;									//default values
	AdvTimeStep			=0.01;
	AdvSpaceStep    =0.001;
	AdvectionType		=ADAPTIVE_TIME_STEP; 
	DispersionType  =EULERIAN; 
	alpha_L     		=0.0;
	alpha_TH     		=0.0;
	alpha_TV     		=0.0;

	UseInitialConditions=false;
}
//------------------------------------------------------------------------
CChemDomain::~CChemDomain(){

	if (globaldebug){cout <<"DESTROYING TRANSPORT DOMAIN... "<<endl;}
	int i,k,s;
	delete [] poro;

	for (s=0; s<nspecies;          s++){delete pSpeciesArray    [s];} delete [] pSpeciesArray;
	for (i=0; i<nSoils;            i++){delete pSoilArray       [i];} delete [] pSoilArray;
	for (i=0; i<nRxNs;             i++){delete pReactSchemes    [i];} delete [] pReactSchemes;
	delete pTransScheme;
	delete pDispScheme;


	if (NodalSourceFluxes!=NULL){for (k=0; k<pConcGrid->GetNumNodes(); k++){delete [] NodalSourceFluxes[k];} delete [] NodalSourceFluxes;}
	if (pNodalSources    !=NULL){for (k=0; k<pConcGrid->GetNumNodes(); k++){delete [] pNodalSources    [k];} delete [] pNodalSources    ;}         
	delete [] nNodalSources;           

	if (pConcGrid!=NULL){
		if (C    !=NULL){for (k=0; k<pConcGrid->GetNumNodes(); k++){delete [] C[k];    }  delete [] C;    }
		if (Ctmp !=NULL){for (k=0; k<pConcGrid->GetNumNodes(); k++){delete [] Ctmp[k]; }  delete [] Ctmp; }
		if (Cnew !=NULL){for (k=0; k<pConcGrid->GetNumNodes(); k++){delete [] Cnew[k]; }  delete [] Cnew; }
		if (Csorb!=NULL){for (k=0; k<pConcGrid->GetNumNodes(); k++){delete [] Csorb[k];}  delete [] Csorb;}                                      
		delete pConcGrid;
	}
	delete [] OutputTimes;
//	delete [] noneq_sorpt;
  
	if (IsothermTable!=NULL){
		for (i=0;i<nSoils; i++){
			for (s=0; s<nspecies; s++){delete IsothermTable[i][s];}
			delete [] IsothermTable[i];
		}
		delete [] IsothermTable;
	}
}
/************************************************************************
                           ACCESSOR/INQUIRY FUNCTIONS
************************************************************************/
//SAME FOR 2D and 3D 
int    CChemDomain::GetNumSpecies    ()            const{return nspecies;} 
//-------------------------------------------------------------------------
double CChemDomain::GetDiffusionCoeff(const int s) const{
	ExitGracefullyIf(((s<0) || (s>=nspecies)),"CChemDomain::GetDiffusionCoeff: Bad Species index",RUNTIME_ERR);
	return pSpeciesArray[s]->GetDiffusionCoeff();
}
//-------------------------------------------------------------------------
double CChemDomain::GetDispersivity(const pt3D &pt,disp_dir dir) const{
	switch (dir){ 
		case(LONGITUDINAL):	{return alpha_L; break;} 
		case(TRANSVERSE):		{return alpha_TH;break;}
		case(VERTICAL):			{return alpha_TV;break;}
	}
	return 0.0;
}
//-------------------------------------------------------------------------
void CChemDomain::GetDispersivities(const pt3D &pt, const int s, disp_info &disp) const{
	disp.al=alpha_L;
	disp.at=alpha_TH;
	disp.D=pSpeciesArray[s]->GetDiffusionCoeff();
}
/************************************************************************
                           GET DISPERSION COEFFICIENT
*************************************************************************
TMP DEBUG -re-orient vector v
          -now only works if alignment.z=0  (transformation in x,y plane only)
returns NON-DIFFUSIVE portion of dispersion coeff in cardinal directions specified
alignment vector is a UNIT vector that depicts the orientation of the x-axis
SAME FOR 2D AND 3D 
-----------------------------------------------------------------------*/
double CChemDomain::CalcDispersionCoeff(vector &v, orientation dir, const CVector &alignment) const{

  if (abs(v)==0.0){return 0.0;}

	
	//TMP DEBUG : TRANSFORM v -currently only two dimensional
	ExitGracefullyIf((alignment.z!=0.0),"CChemDomain::GetDispersionCoeff: cannot rotate into this plane",RUNTIME_ERR);
  
	double cosa,sina,d;
	d=abs(alignment);
	cosa=alignment.x/d;
	sina=alignment.y/d;
	
	v=CVector(v.x*cosa-v.y*sina,v.x*sina+v.y*cosa,v.z);

	double aL =alpha_L; 
	double aTH=alpha_TH;
	double aTV=alpha_TV;

	switch (dir){
		case(OR_XX):{return (aL*v.x*v.x+aTH*v.y*v.y+aTV*v.z*v.z)/(abs(v));break;}
		case(OR_YY):{return (aL*v.y*v.y+aTH*v.x*v.x+aTV*v.z*v.z)/(abs(v));break;}
		case(OR_ZZ):{return (aL*v.z*v.z+aTH*v.x*v.x+aTV*v.y*v.y)/(abs(v));break;}
		case(OR_XY):{return (aL-aTH)*v.x*v.y/(abs(v));break;}
		case(OR_YX):{return (aL-aTH)*v.x*v.y/(abs(v));break;}
		case(OR_XZ):{return (aL-aTV)*v.x*v.z/(abs(v));break;}
		case(OR_ZX):{return (aL-aTV)*v.x*v.z/(abs(v));break;}
		case(OR_YZ):{return (aL-aTV)*v.z*v.y/(abs(v));break;}
		case(OR_ZY):{return (aL-aTV)*v.z*v.y/(abs(v));break;}
	}
	return 0.0;
}
/************************************************************************
                           GET RETARDATION
*************************************************************************
returns retardation factor for equilibrium species, R=1.0 otherwise
-----------------------------------------------------------------------*/
double CChemDomain::GetRetardation (const pt3D &pt, const double &t, const int s) const{

	static double tmpconc[MAX_SPECIES];
	static int    k;
	static double pb;
	static double n;

	if (IsothermTable==NULL){return 1.0;}//shouldn't happen?
	else{
		if (UsesIsotherm[s]==false){return 1.0;}//quick out
	
		k =GetSoilType(pt);
		//if (IsothermTable[k][s].type==no_isotherm){return 1.0;}//quick out
		if (IsothermTable[k][s]==NULL){return 1.0;}//quick out

		pb=pSoilArray[k]->GetBulkDryDensity();
		n =pAq->GetPoro(pt);

		return IsothermTable[k][s]->GetRetardation(pb,n,GetConcentration(pt,t,s));
	}
	return 1.0;

}
//-------------------------------------------------------------------------
//SAME FOR 2D AND 3D 
double CChemDomain::TranslateToSorbed (const pt3D &pt, const double C, const int s) const{

  if (IsothermTable==NULL){return 0.0;}
	else{
		ExitGracefullyIf(!UsesIsotherm[s],"CChemDomain::TranslateToSorbed: Cannot translate species that does not use isotherm",RUNTIME_ERR);

		int    k =GetSoilType(pt);
		

		if (IsothermTable[k][s]==NULL){return 0.0;}

		return IsothermTable[k][s]->TranslateToSorbed(C);
	}
	return 0.0;
}
/************************************************************************
                           GET SPECIES DECAY
*************************************************************************
-----------------------------------------------------------------------*/
double CChemDomain::GetSpeciesDecay (const int s) const{
	return pSpeciesArray[s]->GetDecayRate();
}

/************************************************************************
                           GET CONCENTRATIONS
*************************************************************************
-----------------------------------------------------------------------*/
double CChemDomain::GetConcentration(const pt3D &pt, const double &t, const int s) const{

	static double concs[MAX_SPECIES];
	pConcGrid->InterpolateValues(pt,C,concs,nspecies);
	return concs[s];

}
//-------------------------------------------------------------------------
void CChemDomain::GetConcentrations(const pt3D &pt, const double &t, Writeable1DArray concs) const{

	pConcGrid->InterpolateValues(pt,C,concs,nspecies);

}
//-------------------------------------------------------------------------
double CChemDomain::GetConcentration(const int k,const double &t, const int s  ) const{

	if ((s>=0) && (s<nspecies) && (k>=0) && (k<pConcGrid->GetNumNodes())){
		return C[k][s];
	}
	else {ExitGracefully("CChemDomain::GetConcentration",RUNTIME_ERR);return 0.0;}

}
//-------------------------------------------------------------------------
void CChemDomain::GetConcentrations (const int k,const double &t, Writeable1DArray concs) const{
	
	ExitGracefullyIf((concs==NULL),"CChemDomain::GetConcentrations",RUNTIME_ERR);

	if ( (k>=0) && (k<pConcGrid->GetNumNodes())){
		for (int s=0;s<nspecies; s++){concs[s]=C[k][s];}
	}
	else{ExitGracefully("CChemDomain::GetConcentrations",RUNTIME_ERR);}

}
//-------------------------------------------------------------------------
double  CChemDomain::GetAmbientConc       (const int s) const{
	if ((s>=0) && (s<nspecies) ){return Cback[s];}
	else{ExitGracefully("CChemDomain::GetAmbientConc: bad species index",RUNTIME_ERR);return 0.0;}
}
/************************************************************************
                           GET SORBED CONCENTRATIONS
*************************************************************************
-----------------------------------------------------------------------*/
double  CChemDomain::GetSorbedConc     (const pt3D &pt, const double &t, const int s) const{

	static double concs[MAX_SPECIES];

	if (UsesIsotherm[s]){
		pConcGrid->InterpolateValues(pt,C,concs,nspecies);
		return TranslateToSorbed(pt,concs[s],s);
	}
	else{
		pConcGrid->InterpolateValues(pt,Csorb,concs,nspecies);
		return concs[s]; 
	}
}
//-------------------------------------------------------------------------
void CChemDomain::GetSorbedConcs(const pt3D &pt, const double &t, Writeable1DArray sorb) const {

	static double concs[MAX_SPECIES];
	
	pConcGrid->InterpolateValues(pt,Csorb,sorb,nspecies);

	pConcGrid->InterpolateValues(pt,C,concs,nspecies);
	for (int s=0;s<nspecies;s++){
		if (UsesIsotherm[s]){
			sorb[s]=TranslateToSorbed(pt,concs[s],s);
		}
	}

}
//-------------------------------------------------------------------------
/*void CChemDomain::GetSorbedConcs (const int k,const double &t, Writeable1DArray concs) const{
	
	if (concs==NULL){ExitGracefully("CChemDomain::GetSorbedConcs",RUNTIME_ERR);}
	if ( (k>=0) && (k<pConcGrid->GetNumNodes())){
		for (int s=0;s<nspecies; s++){concs[s]=C[k][s];}
	}
	else{ExitGracefully("CChemDomain::GetSorbedConcs",RUNTIME_ERR);}

}*/
//-------------------------------------------------------------------------
double CChemDomain::GetSorbedConc  (const int k, const double &t, const int s  ) const{

	if ((s>=0) && (s<nspecies) && (k>=0) && (k<pConcGrid->GetNumNodes())){
		pt3D pt=pConcGrid->GetNodeLocation(k);

		if (UsesIsotherm[s]){return TranslateToSorbed(pt,C[k][s],s);}
		else                {return Csorb[k][s];                    }
	}
	else {cout<<k<<" "<<s<<endl;ExitGracefully("CChemDomain::GetConcentration",RUNTIME_ERR);return 0.0;}

}

//-------------------------------------------------------------------------
double CChemDomain::GetNewSorbedConc  (const int k, const double &t, const int s  ) const{

	if ((s>=0) && (s<nspecies) && (k>=0) && (k<pConcGrid->GetNumNodes())){
		pt3D pt=pConcGrid->GetNodeLocation(k);
		if (UsesIsotherm[s]){return TranslateToSorbed(pt,Cnew[k][s],s);}
		else                {return Csorb[k][s];             }//TMP DEBUG -requires Csorbnew?
	}
	else {cout<<k<<" "<<s<<endl;ExitGracefully("CChemDomain::GetConcentration",RUNTIME_ERR);return 0.0;}

}
/************************************************************************
         WRITE CONCENTRATION TO FILE/READ CONCENTRATION FROM FILE
*************************************************************************
-----------------------------------------------------------------------*/
void CChemDomain::WriteConcentrationToFile(const int outputstep) const{

	//format:
	//double t, int numnodes, int nspecies
	char filename[FILENAME_SIZE];
	sprintf(filename,"Concentrations-%d.csv",outputstep);	

	ofstream CONCFILE;
	CONCFILE.open(filename);
	CONCFILE << OutputTimes[outputstep-1]<<","<<pConcGrid->GetNumNodes() << ","<<nspecies<<endl;
	for (int k=0; k<pConcGrid->GetNumNodes(); k++){
	  for (int s=0; s<nspecies; s++){
			CONCFILE << Cnew[k][s] << ",";
		}
		CONCFILE <<endl;
	}
	CONCFILE.close();
}
//-----------------------------------------------------------------------
void CChemDomain::SetInitialConditions (){
	cout << "CChemDomain::SetInitialConditions: Initial Conditions are being used!"<<endl;
	UseInitialConditions=true;
}
//-----------------------------------------------------------------------
void CChemDomain::ReadInitialConditionsFromFile(){
	int        Len,thisnumnodes,thisspecies;
	char     *s[MAXINPUTITEMS];

	ExitGracefullyIf((C==NULL),"CChemDomain::ReadInitialConditionsFromFile: Concentration Array not yet initialized",RUNTIME_ERR);
	ifstream CONCFILE;
	CONCFILE.open("InitialConcs.csv");
	//if (CONCFILE.bad()){ ExitGracefully("CChemDomain::ReadInitialConditionsFromFile: bad initial conditions file",BAD_DATA);}
	if (!TokenizeLine(CONCFILE,s,Len)){ //insure comma delimiter
		if (Len==3){
			StartTime   =s_to_d(s[0]);
			cout <<"START TIME : "<<StartTime<<endl;
			thisnumnodes=s_to_i(s[1]);
			thisspecies =s_to_i(s[2]);
			if (thisnumnodes!=pConcGrid->GetNumNodes()){
				cout <<thisnumnodes<<" "<<pConcGrid->GetNumNodes()<<endl;
				CONCFILE.close();
				ExitGracefully("Bad Initial conditions file: different grid used",BAD_DATA);}
		}
		else{CONCFILE.close();ExitGracefully("Bad Initial conditions file: incorrect length",BAD_DATA);}
		
		for (int k=0; k<pConcGrid->GetNumNodes();k++){
			if (!TokenizeLine(CONCFILE,s,Len)){ 
				if (Len==thisspecies){
					for (int m=0; m<thisspecies; m++){
						Cnew[k][m]=C[k][m]=s_to_d(s[m]);
					}
				}
				else{CONCFILE.close();ExitGracefully("Bad Initial conditions file: incorrect number of species",BAD_DATA);}
			}
			else{CONCFILE.close();ExitGracefully("Bad Initial conditions file: incorrect number of lines",BAD_DATA);}
		}
	}

	CONCFILE.close();
}
/************************************************************************
                           GET SOURCE CONCENTRATION
*************************************************************************
Averages Cin from multiple sources returns Cin=(Q1C1+Q2C2+C3C3...) / (Q1+Q2+Q3...)
-----------------------------------------------------------------------*/
double  CChemDomain::GetSourceConcentration(const int k, const double &t, const int s) const{
	double tmpflux(0.0);

	for (int m=0; m<nNodalSources[k]; m++){
    tmpflux+=NodalSourceFluxes[k][m]*pNodalSources[k][m]->GetConcentration(t,s);
	}
	ExitGracefullyIf((tmpflux<0.0),"CChemDomain::GetSourceConcentration: poor initialization of flux information",RUNTIME_ERR);

	return tmpflux;
}
//-------------------------------------------------------------------------
void CChemDomain::AddSourceFlux(const int k, CSourceSink *pSource, const double &Q){
	
	ExitGracefullyIf((pSource==NULL),"CChemDomain::AddSourceFlux: poor initialization of flux information",RUNTIME_ERR);
  ExitGracefullyIf((Q<0)          ,"CChemDomain::AddSourceFlux: poor identification of flux information",RUNTIME_ERR);
	ExitGracefullyIf(((k>=pConcGrid->GetNumNodes()) || (k<0)),
		                               "CChemDomain::AddSourceFlux: improper cell index"                    ,RUNTIME_ERR);
	
	NodalSourceFluxes[k][nNodalSources[k]]=Q;
  pNodalSources    [k][nNodalSources[k]]=pSource;
	nNodalSources    [k]++;
}

/************************************************************************
                           ASSIGNMENT/MANIPULATOR FUNCTIONS
************************************************************************/
//------------------------------------------------------------------------
void CChemDomain::SetGrid (CMesh *pGrid)         {

	ExitGracefullyIf((pGrid==NULL),"CChemDomain::SetGrid: NULL grid added",RUNTIME_ERR);
  //cout<<"GRID CONSTRUCTED AND ADDED"<<endl;double t;cin>>t;
	pConcGrid=pGrid;

}
//------------------------------------------------------------------------
void CChemDomain::SetOutputTimes(double *times,int numtimes, int interval){	

	ExitGracefullyIf((times==NULL),"CChemDomain:SetOutputTimes: NULL input",BAD_DATA);
	ExitGracefullyIf((numtimes<=0),"CChemDomain:SetOutputTimes: bad number of transport output times (<=0) specified",BAD_DATA);
	ExitGracefullyIf((interval<=0),"CChemDomain:SetOutputTimes: negative interval",BAD_DATA);

	nOutputTimes=numtimes;
	OutputTimes =new double[nOutputTimes];
  for (int i=0; i<nOutputTimes; i++){
		OutputTimes[i]=times[i];
	}
	StraightSort(OutputTimes,nOutputTimes);
	WriteInterval=interval; 

}
//------------------------------------------------------------------------
void CChemDomain::SetAmbientConcs(Ironclad1DArray BC,Ironclad1DArray BS){	

	ExitGracefullyIf((BC==NULL),"CChemDomain:SetAmbientConcs: NULL input",BAD_DATA);
	ExitGracefullyIf((BS==NULL),"CChemDomain:SetAmbientConcs: NULL input(2)",BAD_DATA);
  
	int s;
  for (s=0; s<nspecies; s++){
		Cback[s]=BC[s];
		Sback[s]=BS[s];
	}
	for (s=nspecies;s<MAX_SPECIES;s++){
		Cback[s]=0.0;
		Sback[s]=0.0;
	}

}
//------------------------------------------------------------------------
void CChemDomain::SetObservationPts	(pt3D *obspts, int nobs){

	ExitGracefullyIf((obspts==NULL),"CChemDomain:SetObservationPts: NULL input array",BAD_DATA);
	ExitGracefullyIf((nobs<=0),     "CChemDomain:SetObservationPts: bad number of transport observation points specified",BAD_DATA);
	
	nObsPoints=nobs;
	ObsPoints=new pt3D [nObsPoints];
	for (int i=0; i<nObsPoints; i++){
		ObsPoints[i]=obspts[i];
	}

}
//------------------------------------------------------------------------
//2D ONLY?
void CChemDomain::SetTransportScheme (C2DTransportScheme *pTScheme){

	ExitGracefullyIf((pTScheme==NULL),"CChemDomain::SetTransportScheme: NULL Transport Scheme",BAD_DATA);

	if (pTransScheme!=NULL){delete pTransScheme;}//if overwriting scheme
	pTransScheme=pTScheme;
}

//------------------------------------------------------------------------
//2D ONLY?
void CChemDomain::SetDispersionScheme	(C2DTransportScheme *pDScheme){

	ExitGracefullyIf((pDScheme==NULL),"CChemDomain::SetDispersionScheme: NULL Dispersion Scheme",BAD_DATA);

	if (pDispScheme!=NULL){delete pDispScheme;}//if overwriting existing scheme
	pDispScheme=pDScheme;

}
//------------------------------------------------------------------------
void CChemDomain::AddReactionScheme (CReactionScheme *pRScheme){

	if (!DynArrayAppend((void**&)pReactSchemes,(void*)(pRScheme),nRxNs)){
		ExitGracefully("CChemDomain::SetReactionScheme: adding NULL Reaction Scheme",BAD_DATA);}

}
//------------------------------------------------------------------------
void CChemDomain::SetTransportTimeStep (const double tstep){
	ExitGracefullyIf((tstep<0),"CChemDomain::SetTransportTimeStep: negative time step not allowed",BAD_DATA);
  TimeStep=tstep;
}
//------------------------------------------------------------------------
void CChemDomain::SetVolumeRatio          (const double V){
	ExitGracefullyIf((V<=0.0),"CChemDomain::SetVolumeRatio: negative or zero volume ratio not allowed",BAD_DATA);
	VolumeRatio=V;
}
//------------------------------------------------------------------------
void CChemDomain::SetAdvectionParams(advtype   Atype, double tstep, double sstep){

	ExitGracefullyIf((tstep<0),"CChemDomain::SetAdvectionParams: negative time step",BAD_DATA);
	ExitGracefullyIf((sstep<0),"CChemDomain::SetAdvectionParams: negative space step",BAD_DATA);
  if ((tstep==0) && (Atype==CONSTANT_TIME_STEP)){
		ExitGracefully("CChemDomain::SetAdvectionParams: zero time step not allowed",BAD_DATA);}
  if ((sstep==0) && (Atype==CONSTANT_SPACE_STEP)){
		ExitGracefully("CChemDomain::SetAdvectionParams: zero space step not allowed",BAD_DATA);}

	AdvectionType=Atype;
	AdvTimeStep  =tstep;
	AdvSpaceStep =sstep;

}
//------------------------------------------------------------------------
void CChemDomain::SetDispersionParams(const disptype  Dtype,const double aL, const double aTH, const double aTV){

	ExitGracefullyIf(((aL<0)||(aTH<0)||(aTV<0)),"CChemDomain::SetDispersionParams: negative dispersivity",BAD_DATA);

  DispersionType  =Dtype;	
	alpha_L     		=aL;
	alpha_TH     		=aTH;
	alpha_TV     		=aTV;

}
//------------------------------------------------------------------------
void CChemDomain::EnableEffectiveParams(){EffectiveParams=true;}

/************************************************************************
                           ADD TO DOMAIN FUNCTIONS
-------------------------------------------------------------------------
************************************************************************/
void CChemDomain::AddToDomain(CSpecies *spec){

	if (!DynArrayAppend((void**&)(pSpeciesArray),(void*)(spec),nspecies)){
	  ExitGracefully("CChemDomain::AddToDomain(Species): adding NULL species",BAD_DATA);};
}
//------------------------------------------------------------------------
void CChemDomain::AddToDomain	(CSoilType    *soil){

	if (!DynArrayAppend((void**&)(pSoilArray),(void*)(soil),nSoils)){
	  ExitGracefully("CChemDomain::AddToDomain(Soils): adding NULL soil",BAD_DATA);};
}
//------------------------------------------------------------------------
void CChemDomain::AddIsotherm(const int i, const int s, CIsotherm *iso){

	if (IsothermTable==NULL){
		ExitGracefully("CChemDomain::AddIsotherm: Isotherm Table uninitialized",RUNTIME_ERR);}
	if ((i<0) || (i>=nSoils) || (s<0) || (s>=nspecies)){
		ExitGracefully("CChemDomain::AddIsotherm: Bad soil or species index",BAD_DATA);}
	
	sorbing      [s]=true; 
	UsesIsotherm [s]=true;
	IsothermTable[i][s]=iso;
}


/************************************************************************
                           INITIALIZE ISOTHERM TABLE 
-------------------------------------------------------------------------
done once either in parse or in Transport inititialization
must be done before GridInitialize
VALID IN 2D AND 3D 
************************************************************************/
void CChemDomain::InitializeIsothermTable(){
	
	int i,s;

	if (nSoils  <=0)   {
		nSoils=1;
		pSoilArray      =new CSoilType*[1];
		pSoilArray[0]   =new CSoilType("Default Background",DEFAULT_DRY_DENSITY);
	}
	if (nspecies <=0)   {
		nspecies=1;
		pSpeciesArray   =new CSpecies*[1];
		pSpeciesArray[0]=new CSpecies("Default Species",0.0,1.0,0.0);
	}
	
	for (s=0; s<nspecies; s++){sorbing[s]=false;}

	/*for (s=0; s<nspecies; s++){
    if (false){//pSpeciesArray[s]->IsNonEquilibrium() //TMP DEBUG
			UsesIsotherm[s]=false;
		}
		else {
      UsesIsotherm[s]=true; //could be dangerous...
		}
	}*/

	//modify for surface reactions  (e.g., Cation exchange) //TMP DEBUG
	for (i=0;i<nRxNs;i++){
		if (pReactSchemes[i]->IsSorptionReaction()){
			for (s=0; s<nspecies; s++){
				sorbing     [s]=true;
				UsesIsotherm[s]=false;
			}
		}
	}
	if (IsothermTable==NULL){
		IsothermTable =new CIsotherm **[nSoils];
		for (i=0;i<nSoils; i++){
			IsothermTable[i]=new CIsotherm *[nspecies];
			for (s=0; s<nspecies; s++){
				IsothermTable[i][s]=NULL;	
			}		
		}
	}
	else{
		ExitGracefully("CChemDomain::InitializeIsothermTable: trying to re-initialize sorption isotherm table",BAD_DATA);
	}

}
/************************************************************************
                           INITIALIZE GRID (done once)
-------------------------------------------------------------------------
Creates arrays to store concentration information
initializes concentration values to background/initial conditions
2D ONLY : FLOW FIELD INFORMATION, BUDGETPOLY STATISTICS
2D/3D : CONCENTRATION ARRAYS 
************************************************************************/
/*void CChemDomain::InitializeGrid  (){

	int           k,s;
	double        junkC[MAX_SPECIES];
  double t(0.0);

	//Check quality of data---------------------------------------------------
	if (pConcGrid==NULL)  {ExitGracefully("CChemDomain::InitializeGrid: Concentration Grid not created.",BAD_DATA);}
	if (nspecies<=0)      {ExitGracefully("CChemDomain::InitializeGrid: No species to transport",BAD_DATA);}
	if (noneq_sorpt==NULL){ExitGracefully("CChemDomain::InitializeGrid: not properly initialized",RUNTIME_ERR);}

	int nNodes(pConcGrid->GetNumNodes());
  if (nNodes<=0)        {ExitGracefully("CChemDomain::InitializeGrid:no nodes to create",BAD_DATA);}

	if (ProgramAborted()){return;}



	cout <<"Initializing Concentration Grid/Mesh ("<<nNodes<<" nodes)..."<<endl;
	//--------------------------------------------------------------------------
  C    =new double *[nNodes];
	Ctmp =new double *[nNodes];
	Cnew =new double *[nNodes];	
  Csorb=new double *[nNodes];
	if (Csorb==NULL)     {ExitGracefully("CChemDomain::InitializeGrid:out of memory(1b).",OUT_OF_MEMORY);}

	for (k=0; k<nNodes; k++){ 
    C    [k]=new double [nspecies];
    Ctmp [k]=new double [nspecies];
		Cnew [k]=new double [nspecies];
    Csorb[k]=new double [nspecies];
	  if (Cnew[k]==NULL){ExitGracefully("CChemDomain::InitializeGrid:out of memory(2b).",OUT_OF_MEMORY);}
		
		//Initialize concentrations
		for (s=0; s<nspecies; s++){
			C    [k][s]=pSpeciesArray[s]->GetBckgrndConc();
			Ctmp [k][s]=0.0;
		  Cnew [k][s]=0.0;		
		  Csorb[k][s]=0.0;		
		}
	}
	
	cout <<"Initializing Flow field Information ("<<nNodes<<" nodes)..."<<endl;
	//-------------------------------------------------------------------------
	if (pAq==NULL)        {ExitGracefully("CChemDomain::InitializeGrid: Aquifer is NULL", RUNTIME_ERR);}
	poro    =new double [nNodes];
	satthick=new double [nNodes];
	if (satthick==NULL)     {ExitGracefully("CChemDomain::InitializeGrid:out of memory(flow).",OUT_OF_MEMORY);}
	
	for (k=0;k<nNodes;k++){
		poro    [k]=pAq->GetPoro              (pConcGrid->GetNodeLocation(k));
		satthick[k]=pAq->GetSaturatedThickness(pConcGrid->GetNodeLocation(k),0.0);
	}

	//Add Initial Conditions-------------------------------------------------
	AdjustForSourceTerms(t,0.0,junkC);
	for (k=0; k<nNodes; k++){ 
		for (s=0; s<nspecies; s++){
			C[k][s]=Cnew[k][s];	
		}
	}
	
	//Print Grid Geometry----------------------------------------------------
	pConcGrid  ->WriteGeometry();

	cout <<"...done initializing Concentration Grid."<<endl;
}  */


/************************************************************************
													TRANSPORT
*************************************************************************
	Transports contaminants using Operator splitting
	2D/3D- mostly generic with exception of 
		budgetpoly calculations, satthick in initialization of transport scheme 
-----------------------------------------------------------------------*/
/*
void    CChemDomain::Transport(ofstream &PROGRESS){

	int     i,s;
	double  t(0.0),simulation_end_time;
	cmplex  tmpz;
	pt3D    tmppt;
	bool    SIAdone(true);
	int     SIAiter;
	int     latesttimeplotted(0);
  bool    stopped(false);
  int     nNodes;

	double  Csorbed      [MAX_SPECIES];

	double  sys_mass     [MAX_SPECIES];  //All used for mass budget
	double  massadded    [MAX_SPECIES];
	double  adv_error    [MAX_SPECIES];
	double  disp_error   [MAX_SPECIES];
	double  boundary_loss[MAX_SPECIES];
	double  feature_loss [MAX_SPECIES]; 
	double  leakage_loss [MAX_SPECIES];
	double  rxn_change   [MAX_SPECIES];

	cout<<endl<<endl<<"_______________*_Transporting Contaminant_*______________________"<<endl<<endl;

  PROGRESS<<"transport"<<endl<<0<<endl;

	//check for bad data----------------------------------------------------
	if (pConcGrid==NULL)   {ExitGracefully("CChemDomain::Transport: Mesh is NULL"                           , BAD_DATA);}
	if (pTransScheme==NULL){ExitGracefully("CChemDomain::Transport: TransScheme is NULL"                    , BAD_DATA);}
	if (nconctimes<=0)     {ExitGracefully("CChemDomain::Transport: No concentration output times specified", BAD_DATA);}
	if (TimeStep<=0)       {ExitGracefully("CChemDomain::Transport: Invalid transport timestep specified"   , BAD_DATA);}
  
	//if ((DispersionType==EULERIAN) && (pDispScheme==NULL)){
	//	                     ExitGracefully("CChemDomain::Transport: Dispersion scheme is NULL", RUNTIME_ERR);}
  if ((DispersionType!=EULERIAN) && (pDispScheme!=NULL)){delete pDispScheme;pDispScheme=NULL;}



	//prepare output file for monitoring points-----------------------------
	ofstream OBSCONC;
	OBSCONC.open("observed_conc.csv");
	OBSCONC<<"t, ";
	for (i=0; i<nObsPoints; i++){
		for (s=0; s<nspecies; s++){
			OBSCONC<<"pt "<<i<<" "<<pSpeciesArray[s]->GetName();
			if ((i!=nObsPoints-1) || (s!=nspecies-1)){OBSCONC<<",";}
		}
	}
	OBSCONC<<endl;

	ofstream MB;
	MB.open("MassBalance.csv");
	MB<<"time,";
	for (s=0; s<nspecies; s++){
		MB<<"mass "<<s+1<<", source "<<s+1<<", adv "<<s+1<<", disp "<<s+1<<", boundary "<<s+1<<",feature "<<s+1<<", leakage "<<s+1<<", reaction"<<s+1;
		if (s!=nspecies-1){MB<<",";}
	}
	MB<<endl;

	//Initialize Isotherm Table (if not done previously)--------------------
	if (IsothermTable==NULL) {InitializeIsothermTable();} 

	//Initialize Concentration Grid-----------------------------------------
  InitializeGrid();
	if (ProgramAborted()){stopped=true;}

	//Initialize source/sink information------------------------------------
	InitializeSources();
	if (ProgramAborted()){stopped=true;}

	//Initialize Transport/Dispersion Schemes-------------------------------
	pTransScheme->SetAdvectionParams(AdvectionType, AdvTimeStep, AdvSpaceStep);

	if (pDispScheme!=NULL){ 
		pTransScheme->Initialize(ADVECTION_ONLY ,poro,satthick);//TMP DEBUG
		pDispScheme ->Initialize(DISPERSION_ONLY,poro,satthick);
	}
	else{
		pTransScheme->Initialize(ADV_AND_DISP   ,poro,satthick);
	}

	//Initialize Local Variables--------------------------------------------
  nNodes=pConcGrid->GetNumNodes();
	simulation_end_time=OutputTimes[nOutputTimes-1];

 	//Plot Initial Concentrations, write initial output (if required) ------
	if ((latesttimeplotted<nOutputTimes) && (OutputTimes[latesttimeplotted]==0)){
		latesttimeplotted++;
		for (s=0; s<nspecies; s++){
			PlotConcentration(latesttimeplotted,s);
			PlotSorbedConc   (latesttimeplotted,s);
		}
		pTransScheme->WriteOutput(latesttimeplotted);
	}

	SIAdone=false;
	SIAiter=0;

	//Transport Contaminant-------------------------------------------------
	do {
		if (ProgramAborted()){stopped=true;}
		if (!stopped){
	
			cout << endl << "Transport time: " <<t;

			OBSCONC<<t<<",";
			for (i=0; i<nObsPoints; i++){
				for (s=0; s<nspecies; s++){
				  OBSCONC<<GetConcentration(ObsPoints[i],t,s);
					if ((i!=nObsPoints-1) || (s!=nspecies-1)){OBSCONC<<",";}
				}
			}
			OBSCONC<<endl;
			
			do { //sequential iterative approach loop-not presently used

				for (s=0; s<nspecies; s++){
					massadded[s]=adv_error[s]=disp_error[s]=boundary_loss[s]=feature_loss[s]=leakage_loss[s]=0.0;
          sys_mass [s]=GetTotalSystemMass(s,t);
				}

				//Advective/Dispersive Transport of Contaminants-------------------------------------------- 
				CopyConcentrationArray(Cnew,Ctmp);
				pTransScheme->Transport(t,TimeStep,C,Ctmp,Cnew);  
				
				for (s=0; s<nspecies; s++){
					for (m=0; m<nlinesourcesinks; m++){
						feature_loss[s]+=0.0;//pLineSinks[m]->GetNetMassFlux(t,grdConc)*TimeStep; //TMP DEBUG
					}
					leakage_loss [s]=0.0;//TMP DEBUG
					adv_error    [s]+=GetNewSystemMass(s,t)-sys_mass[s]-boundary_loss[s];
					sys_mass     [s] =GetNewSystemMass(s,t);
				}
							

				//Eulerian Dispersive Transport of Contaminants---------------------------------------------
				if ((DispersionType==EULERIAN) && (pDispScheme!=NULL)){
          CopyConcentrationArray(Cnew,Ctmp);
					pDispScheme->Transport(t,TimeStep,C,Ctmp,Cnew);  
				
					for (s=0; s<nspecies; s++){
						disp_error[s]+=GetNewSystemMass(s,t)-sys_mass[s];
						sys_mass  [s] =GetNewSystemMass(s,t);
					}
				}

				//Reaction-----------------------------------------------------------------------
				//sift through cells, perform batch reactions at each
				CopyConcentrationArray(Cnew,Ctmp);
				if (pReactScheme!=NULL){
					for (k=0; k<nNodes; k++){			
						tmppt=pConcGrid->GetNodeLocation(k);
						//tmpz =c3Dto2D(tmppt);
						for (s=0; s<nspecies; s++){
							if (noneq_sorpt[s]){Csorbed[s]=Csorb[k][s];                          }
							else               {Csorbed[s]=TranslateToSorbed(tmppt,Ctmp[k][s],s);}
						}
						pReactScheme->React(C[k],    Ctmp[k], Cnew[k],
																Csorbed, Csorbed, Csorbed, //tmp debug
																GetSoilType(tmppt),
																pAq->GetPoro(tmppt),
																TimeStep);  
					}
					for (s=0; s<nspecies; s++){
						rxn_change[s]+=GetNewSystemMass(s,t)-sys_mass[s];
						sys_mass  [s] =GetNewSystemMass(s,t);
					}
				}
		
				//Adjust for source zones/points------------------------------------------------*/
				/*AdjustForSourceTerms(t,TimeStep,massadded);
				for (s=0; s<nspecies; s++){
					massadded [s]+=GetNewSystemMass(s,t)-sys_mass[s];
					sys_mass  [s] =GetNewSystemMass(s,t);
				}*/
				//test for SIA convergence-not presently used
				/*for (k=0; k<nNodes; k++){			
					for (s=0; s<nspecies; s++){
					  upperswap(SIAmaxchange,Cnew[k][s]-Csia[k][s]);
						max_conc=upperswap(max_conc,Cnew[k][s])
						SIAmaxchange/=max_conc;
					}
				}
				CopyConcentrationArray(Cnew,Csia);*/ //may wish to use relaxation 

				/*SIAdone=true;//TMP DEBUG
				SIAiter++;

			}while ((SIAdone==false) && (SIAiter<MAX_SIA_ITER));

      //total update of concentrations------------------------------------
		/*	CopyConcentrationArray(Cnew,C);

 			//Plot Concentrations, write output (if required) ------------------
			if ((latesttimeplotted<nOutputTimes) && (t+TimeStep>=OutputTimes[latesttimeplotted])){
				latesttimeplotted++;
				for (s=0; s<nspecies; s++){
					PlotConcentration       (latesttimeplotted,s);
					PlotSorbedConc          (latesttimeplotted,s);
					WriteConcentrationToFile(latesttimeplotted);
				}
				pTransScheme->WriteOutput(latesttimeplotted);	
			}
			
			//Plot Mass Balance information-------------------------------------
			MB<<t<<",";

			for (s=0;s<nspecies;s++){
				MB<<sys_mass     [s]<<","
					<<massadded    [s]<<","
					<<adv_error    [s]<<","
					<<disp_error   [s]<<","
					<<boundary_loss[s];
				if (s!=nspecies-1){MB<<",";}
			}
			MB<<endl;

			//increment time----------------------------------------------------
			t+=TimeStep;	
		
			PROGRESS<<t          <<" "
				      <<sys_mass[0]<<" "
							<<massadded[0]-boundary_loss[0]<<" "
							<<adv_error[0]<<" "<<disp_error[0]<<endl;

		}//if (!stopped)...

	}	while ((t<simulation_end_time) && (!stopped));

	PROGRESS<<"done"<<endl;
	cout << "...Transport calculations done"<<endl;

	OBSCONC.close();
	MB     .close();

}*/
//---------------------------------------------------------------------
//2D AND 3D 
void CChemDomain::CopyConcentrationArray(Ironclad2DArray  c1,
																				 Writeable2DArray cout){
	for (int k=0; k<pConcGrid->GetNumNodes(); k++){
		for (int s=0; s<nspecies; s++){
			cout[k][s]=c1[k][s];;
		}	
	}
}
