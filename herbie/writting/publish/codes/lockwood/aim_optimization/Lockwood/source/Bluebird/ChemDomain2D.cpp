//CChemDomain2D
#include "ChemDomain.h"

CChemDomain2D *pCurrDomain;
int 					 gridspecies;

bool CDomainABC 	::EffectiveParams=false;

/*********************************************************************
													 TRANSPORT DOMAIN
**********************************************************************
													 CONSTRUCTORS
**********************************************************************/
CChemDomain2D::CChemDomain2D():
							 CChemDomain(){
													
	//Domain always totally empty at creation 			

	pAq=NULL; 														//Initialized in SetAquifer

	satthick				 =NULL;

	pAreaSources		 =NULL; nsources=0; 	//initialized/dynamically modified in AddToDomain Routines	 
	pLineSourceSinks =NULL; nlinesourcesinks =0;	
	pPointSourceSinks=NULL; npointsourcesinks=0;
	pSoilZones			 =NULL; nSoilZones=0;  
	pKdZone 				 =NULL;

	pCurrDomain 		=this;								//for gridding routine
	gridspecies 		=0;
	MBEnabled 			=true;
	initialize_only =false;
	PlotSorbed      =false;
}
//------------------------------------------------------------------------
CChemDomain2D::~CChemDomain2D(){

	if (globaldebug){cout <<"DESTROYING 2D TRANSPORT DOMAIN... "<<endl;}
	int m;

	delete [] satthick;

	for (m=0; m<nsources; 				 m++){delete pAreaSources 		[m];} delete [] pAreaSources;
	for (m=0; m<npointsourcesinks; m++){delete pPointSourceSinks[m];} delete [] pPointSourceSinks;
	for (m=0; m<nlinesourcesinks;  m++){delete pLineSourceSinks [m];} delete [] pLineSourceSinks;

	delete [] pSoilZones; //individuals destroyed by PropertyZones:Destroy all

} 
void		CChemDomain2D::DisableMassBalance 		 (){MBEnabled=false;}
/************************************************************************
													 GET DISPERSION COEFFICIENT
*************************************************************************
Gets dispersion coefficient at point - should only use if velocity has not yet been obtained
-----------------------------------------------------------------------*/
double CChemDomain2D::GetDispersionCoeff(const pt3D &pt, const double &t, orientation dir, const CVector &alignment) const{

	vector v=pAq->GetVelocity3D(pt,t); 

	return CChemDomain::CalcDispersionCoeff(v,dir,alignment);

}
//-------------------------------------------------------------------------
int CChemDomain2D::GetSoilType(const pt3D &pt) const{

	//background aquifer soil type is always soilindex==0
	int i=(int)(CPropZone::NestedSift(c3Dto2D(pt),pSoilZones,nSoilZones,0));
	if ((i<0) || (i>=nSoils)){ExitGracefully("CChemDomain2D::GetSoilType: bad soil zone returned",RUNTIME_ERR);}
	return i;

}
//-------------------------------------------------------------------------
double	CChemDomain2D::GetBulkDryDensity			 (const pt3D &pt) const{
	//returns [kg/l]
	return pSoilArray[GetSoilType(pt)]->GetBulkDryDensity();
}
/************************************************************************
													 PLOT CONCENTRATION
*************************************************************************
ONLY GOOD FOR 2D SYSTEMS
-----------------------------------------------------------------------*/
void grdConcs(const cmplex &z,const double &t,double *concs, const int nspec) {
	
	// Acts as shell function for plotting routine
	if (pCurrDomain->GetGrid()->IsInside(c2Dto3D(z))){	
		pCurrDomain->GetConcentrations(c2Dto3D(z),t,concs);
		for (int s=0;s<nspec;s++){
			if (concs[s]<0){concs[s]=0.0;}//SURFERMAGIC;}//TMP DEBUG
		}
	}
	else{
		for (int s=0;s<nspec;s++){concs[s]=SURFERMAGIC;}
	}
}
//--------------------------------------------------------------------------
void grdSorbedConcs(const cmplex &z,const double &t,double *concs, const int nspec) {
	
	// Acts as shell function for plotting routine
	pCurrDomain->GetSorbedConcs(c2Dto3D(z),t,concs);
	for (int s=0;s<nspec;s++){
		if (concs[s]<0){concs[s]=SURFERMAGIC;}
	}
}
//--------------------------------------------------------------------------
double grdKd(const cmplex &z, const double &t){
	return pCurrDomain->GetKd(z);
}
//-------------------------------------------------------------------------
void CChemDomain2D::PlotConcentrations(const int outputstep) const{
	
	char **filenames;
	filenames=new char *[nspecies];
	
	window GridWindow;
	int GridRes,s;
	GridWindow=pConcGrid->GetBoundingBox();
	GridRes 	=max((int)(1.5*sqrt((double)(pConcGrid->GetNumNodes()))),100); 
	//TMP DEBUG 	- mt3d compare
	/*GridWindow.w=0;
	GridWindow.e=1000.0;
	GridWindow.s=0;
	GridWindow.n=1000.00;
	GridRes=100;*/
	GridRes=150; //TMP DEBUG
	
	for ( s=0;s<nspecies;s++){
		filenames[s]=new char[255];sprintf(filenames[s],"conc-%d-%d.grd",outputstep,s);
	}
	vGrid(filenames,GridWindow,GridRes,nspecies,grdConcs,outputstep);
	for ( s=0;s<nspecies;s++){delete [] filenames[s];}delete [] filenames;

}
//-------------------------------------------------------------------------
void CChemDomain2D::PlotSorbedConcs(const int outputstep) const{

	if (PlotSorbed)
	{
		char **filenames;
		filenames=new char *[nspecies];
		
		window GridWindow;
		int GridRes,s;
		GridWindow=pConcGrid->GetBoundingBox();
		GridRes 	=max((int)(1.5*sqrt((double)(pConcGrid->GetNumNodes()))),100); 
		//TMP DEBUG 	
		/*GridWindow.w=0;
		GridWindow.e=100.0;
		GridWindow.s=0;
		GridWindow.n=100.00;
		GridRes=100;*/
		GridRes=150; //TMP DEBUG

		for ( s=0;s<nspecies;s++){
			filenames[s]=new char[255];sprintf(filenames[s],"sorbed-%d-%d.grd",outputstep,s);
		}
		vGrid(filenames,GridWindow,GridRes,nspecies,grdSorbedConcs,outputstep);

		for ( s=0;s<nspecies;s++){delete [] filenames[s];}delete [] filenames;
	}
}

/************************************************************************
									 GET SOURCE CONCENTRATIONS
*************************************************************************
-----------------------------------------------------------------------*/
double	CChemDomain2D::GetDirichletConc(const pt3D &pt, const double &t, const int s) const{


	static double tmpconc[MAX_SPECIES];
	double thisconc(NO_INFLUENCE);
	cmplex z=c3Dto2D(pt);

	for (int m=0;m<nsources; m++){
		if (pAreaSources[m]->IsInside(z)){
			pAreaSources[m]->GetSpecifiedConcentrations(z,t,tmpconc); 
			thisconc=tmpconc[s];
		}
	}
	return thisconc;

}
//-------------------------------------------------------------------------
double	CChemDomain2D::GetRechargeConcentration(const pt3D &pt, const double &t, const int s) const{

	static double tmpconc[MAX_SPECIES];
	double thisconc(0.0);
	cmplex z=c3Dto2D(pt);

	for (int m=0;m<nsources; m++){
		if (pAreaSources[m]->IsInside(z)){
			pAreaSources[m]->GetRechargeConcentrations(z,t,tmpconc); 
			thisconc=tmpconc[s];
		}
	}
	//if (thisconc>0){
	//cout <<"Getting Recharge Concs!"<<thisconc<<endl;
	//}
	return thisconc;

}

//-------------------------------------------------------------------------
double	CChemDomain2D::GetSpecifiedMassFlux(const pt3D &pt, const double &t, const int s) const{
 
	static double tmpflux[MAX_SPECIES];
	double thisflux(0.0);
	cmplex z=c3Dto2D(pt);
	
	for (int m=0;m<nsources; m++){
		//if (pAreaSources[m]->GetType()==SPECIFIED_FLUX){
		if (pAreaSources[m]->IsInside(z)){
			pAreaSources[m]->GetSpecifiedMassFlux(z,t,tmpflux); //Mass flux per unit area [mg/L^2]
			thisflux=tmpflux[s];
		}
		//}
	}
	
	return thisflux;
}
//-------------------------------------------------------------------------
double CChemDomain2D::GetRetardation (const pt3D &pt, const double &t, const int s) const{
	if (pKdZone==NULL){return CChemDomain::GetRetardation(pt,t,s);}
	else {
		double pb=1.0; //TMP DEBUG
		double n=pLayer->GetPoro(c3Dto2D(pt));
		return 1.0+pb/n*pKdZone->GetValue(c3Dto2D(pt));
	}
}
//-------------------------------------------------------------------------
double CChemDomain2D::TranslateToSorbed (const pt3D &pt, const double C, const int s) const{
	//right now, Kd for single species only
	if (pKdZone==NULL){return CChemDomain::TranslateToSorbed(pt,C,s);}
	else {
		return pKdZone->GetValue(c3Dto2D(pt))*C;
	}
}

//-------------------------------------------------------------------------
void CChemDomain2D::PlotKd() const{
	if (pKdZone!=NULL){
	
		char filename[FILENAME_SIZE];
		sprintf(filename,"Kd.grd");

		window GridWindow;
		int GridRes;
		GridWindow=pConcGrid->GetBoundingBox();
		GridRes 	=max((int)(1.5*sqrt((double)(pConcGrid->GetNumNodes()))),150); 
		rGrid(filename,GridWindow,GridRes,grdKd,0);
	}
}
double CChemDomain2D::GetKd(const cmplex &z) const{
	return pKdZone->GetValue(z);
}


/************************************************************************
													 ASSIGNMENT/MANIPULATOR FUNCTIONS
************************************************************************/
//SHOULD BE "SETLAYER" IN 2D 
void CChemDomain2D::SetAquifer(const CAquifer  *pAquifer){

	if (pAquifer==NULL){
		ExitGracefully("CChemDomain2D::SetAquifer: NULL Aquifer added",RUNTIME_ERR);}

	pAq=pAquifer; 

}
//------------------------------------------------------------------------
void CChemDomain2D::SetLayer(const CSingleLayerABC	*pLay){

	if (pLay==NULL){
		ExitGracefully("CChemDomain2D::SetLayer: NULL Layer added",RUNTIME_ERR);}

	pLayer=pLay; 

}
//------------------------------------------------------------------------
/*void CChemDomain2D::SetGrid (CMesh *pGrid)				 {

	if (pGrid==NULL){
		ExitGracefully("CChemDomain2D::SetGrid: NULL grid added",RUNTIME_ERR);}

	pConcGrid=pGrid;
}*/

/************************************************************************
													 ADD TO DOMAIN FUNCTIONS
-------------------------------------------------------------------------
************************************************************************/
void CChemDomain2D::AddToDomain(CAreaSource *zone){

	if (!DynArrayAppend((void**&)(pAreaSources),(void*)(zone),nsources)){
		ExitGracefully("CChemDomain2D::AddToDomain(AreaSource): adding NULL source",BAD_DATA);};
}
//------------------------------------------------------------------------
void CChemDomain2D::AddToDomain(CPointSourceSink *ptsource){

	if (!DynArrayAppend((void**&)(pPointSourceSinks),(void*)(ptsource),npointsourcesinks)){
		ExitGracefully("CChemDomain2D::AddToDomain(PointSourceSink): adding NULL source/sink",BAD_DATA);};
}
//------------------------------------------------------------------------
void CChemDomain2D::AddToDomain(CLinearSourceSink *sink){

	if (!DynArrayAppend((void**&)(pLineSourceSinks),(void*)(sink),nlinesourcesinks)){
		ExitGracefully("CChemDomain2D::AddToDomain(LinearSourceSink): adding NULL source/sink",BAD_DATA);};
}
//------------------------------------------------------------------------
void CChemDomain2D::AddToDomain (CPropZone		*zone){
	if (zone->GetType()==soil_type){
		if (!DynArrayAppend((void**&)(pSoilZones),(void*)(zone),nSoilZones)){
			ExitGracefully("CChemDomain2D::AddToDomain(SoilZones): adding NULL soil zone",BAD_DATA);};
	} 
}
//------------------------------------------------------------------------
void CChemDomain2D::AddKdZone(CPropZone 				*kdzone){
	if (kdzone==NULL){ExitGracefully("CChemDomain2D::AddKdZone: adding NULL zone",BAD_DATA);}
	pKdZone=kdzone;
}
/************************************************************************
													 INITIALIZE GRID (done once)
-------------------------------------------------------------------------
Creates arrays to store concentration information
initializes concentration values to background/initial conditions
2D ONLY : FLOW FIELD INFORMATION, BUDGETPOLY STATISTICS
2D/3D : CONCENTRATION ARRAYS 
************************************************************************/
void CChemDomain2D::InitializeGrid	(){

	int 					k,s;
	double				junkC[MAX_SPECIES];
	double				t(0.0); //Should be starttime

	//Check quality of data---------------------------------------------------
	if (pConcGrid==NULL)	{ExitGracefully("CChemDomain2D::InitializeGrid: Concentration Grid not created.",BAD_DATA);}
	if (nspecies<=0)			{ExitGracefully("CChemDomain2D::InitializeGrid: No species to transport",BAD_DATA);}

	int nNodes(pConcGrid->GetNumNodes());
	if (nNodes<=0)				{ExitGracefully("CChemDomain2D::InitializeGrid:no nodes to create",BAD_DATA);}

	if (ProgramAborted()){return;}



	cout <<"Initializing Concentration Grid/Mesh ("<<nNodes<<" nodes)..."<<endl;
	//--------------------------------------------------------------------------
	C 	 =new double *[nNodes];
	Ctmp =new double *[nNodes];
	Cnew =new double *[nNodes]; 
	Csorb=new double *[nNodes];
	if (Cnew==NULL) 		{ExitGracefully("CChemDomain2D::InitializeGrid:out of memory(1b).",OUT_OF_MEMORY);}
	
	for (k=0; k<nNodes; k++){ 
		C 	 [k]=new double [nspecies];
		Ctmp [k]=new double [nspecies];
		Cnew [k]=new double [nspecies];
		Csorb[k]=new double [nspecies];
		if (Csorb[k]==NULL){ExitGracefully("CChemDomain2D::InitializeGrid:out of memory(2b).",OUT_OF_MEMORY);}
		
		//Initialize concentrations to background values
		for (s=0; s<nspecies; s++){
			C 	 [k][s]=Cback[s];
			Ctmp [k][s]=Cback[s];
			Cnew [k][s]=Cback[s]; 
			Csorb[k][s]=Sback[s]; 	
		}
	}

	cout <<"Initializing Flow field Information ("<<nNodes<<" nodes)..."<<endl;
	//-------------------------------------------------------------------------
	if (pAq==NULL){
		ExitGracefully("CChemDomain2D::InitializeGrid: Aquifer is NULL", RUNTIME_ERR);}
	poro		=new double [nNodes];
	satthick=new double [nNodes];

	if (satthick==NULL){
		ExitGracefully("CChemDomain2D::InitializeGrid:out of memory(flow).",OUT_OF_MEMORY);}
	
	for (k=0;k<nNodes;k++){
		poro		[k]=pLayer->GetPoro							 (c3Dto2D(pConcGrid->GetNodeLocation(k)));
		satthick[k]=pLayer->GetSaturatedThickness(c3Dto2D(pConcGrid->GetNodeLocation(k)),t);
		//cout <<"satthick[k]:"<<satthick[k]<<" poro[k] " <<poro[k]<< endl;
	}

	//Add Initial Conditions-------------------------------------------------
	AdjustForSourceTerms(t,0.0,junkC); //zero time step

	CopyConcentrationArray(Cnew,C);

	//Print Grid Geometry----------------------------------------------------
	//pConcGrid  ->WriteGeometry(); //Not neccesary- done when read from file or created

	cout <<"...done initializing Concentration Grid."<<endl;
}  
/************************************************************************
													InitializeSources
*************************************************************************
Creates Important Nodal source & Sink matrices
identifies flux for each node or cell, associates nodes with source 
2D ONLY 
-----------------------------------------------------------------------*/
void CChemDomain2D::InitializeSources (const double &t){
	
	int k,m;

	cout <<"Initializing Sinks and Sources..."<<endl;

	int nNodes(pConcGrid->GetNumNodes());
	if (nNodes<=0){ExitGracefully("CChemDomain2D::InitializeGrid:no nodes to create",BAD_DATA);}

	cout <<"  Number of point  source/sinks "<<npointsourcesinks<<endl;
	cout <<"  Number of linear source/sinks "<<nlinesourcesinks <<endl;

	for (m=0;m<npointsourcesinks; m++){
		pPointSourceSinks[m]->CalculateFluxes(t);
		//cout <<"pt "<<m<<"calculating fluxes"<<endl;
	}
	for (m=0;m<nlinesourcesinks; m++){
		cmplex z1a,z2a;
		pLineSourceSinks[m]->CalculateFluxes(t);
		//cout <<"ln "<<m<<"calculating fluxes"<<pLineSourceSinks[m]->GetInflux(t)<<" ---"<<pLineSourceSinks[m]->GetOutflux(t)<<endl;
		//pLineSourceSinks[m]->GetEndpoints(z1a,z2a);
	}

	//Initializing Source Arrays-------------------------------------------
	//---------------------------------------------------------------------
	NodalSourceFluxes=new double			 *[nNodes];
	pNodalSources 	 =new CSourceSink **[nNodes];
	nNodalSources 	 =new int 					[nNodes];
	if (nNodalSources==NULL)		 {
		ExitGracefully("CChemDomain2D::InitializeSources:out of memory(1).",OUT_OF_MEMORY);}
	for (k=0; k<nNodes;k++){
		NodalSourceFluxes[k]=new double 			[MAX_SOURCES_PER_CELL];
		pNodalSources 	 [k]=new CSourceSink *[MAX_SOURCES_PER_CELL];
		if (pNodalSources[k]==NULL) 		{
			ExitGracefully("CChemDomain2D::InitializeSources:out of memory(2).",OUT_OF_MEMORY);}
		nNodalSources 	 [k]=0;
		for (m=0; m<MAX_SOURCES_PER_CELL; m++){
			NodalSourceFluxes[k][m]=0.0;
			pNodalSources 	 [k][m]=NULL; 
		}
	}

	//===========================================================================================
	if (pConcGrid->GetType()==CELL_BASED){

		//linear sources/sinks------------------------------------------
		cmplex	z1a,z2a;
		int 		cellscrossed; 
		pt3D	 *ptint=new pt3D[302];
		int 	 *ids  =new int [300];

		for (m=0;m<nlinesourcesinks; m++){
			cellscrossed=300;
			pLineSourceSinks[m]->GetEndpoints(z1a,z2a);
			pConcGrid->GetMeshBoundaryIntercepts(c2Dto3D(z1a),c2Dto3D(z2a),ptint,ids,cellscrossed);
			for(int j=0; j<cellscrossed;j++){ 
				AddSourceFlux(ids[j],pLineSourceSinks[m],pLineSourceSinks[m]->GetOutfluxThruSegment(c3Dto2D(ptint[j]),c3Dto2D(ptint[j+1]),t));
			}

		}
		delete [] ptint;
		delete [] ids;

		//point sources/sinks-------------------------------------------------- 
		for (m=0;m<npointsourcesinks; m++){
			if (pConcGrid->GetCellIndex(c2Dto3D(pPointSourceSinks[m]->GetZ()),k)){
				AddSourceFlux(k,pPointSourceSinks[m],pPointSourceSinks[m]->GetOutflux(t));
			}
			else{
				//Not on grid
			}
		}
	}
	//===========================================================================================
	else{ //Node based 
		//point sources/sinks-------------------------------------------------- 
		for (m=0;m<npointsourcesinks; m++){
			//if (pConcGrid->GGetCellIndex(c2Dto3D(pPointSourceSinks[m]->GetZ()),k)){
			//	AddSourceFlux(k,pPointSourceSinks[m],pPointSourceSinks[m]->GetOutflux(t));
			//}
		} 
	}

	//ExitGracefully("dead on arrival",BAD_DATA);

	cout <<"...done initializing Sinks and Sources."<<endl;
	
}
/************************************************************************
													Adjust For Source Terms
*************************************************************************
Adjusts matrix Cnew to account for source influence
mass added is the total amount of species mass added by source terms
should only be used for LAGRANGIAN DISPERSION SCHEMES!!!
2D ONLY 
-----------------------------------------------------------------------*/
void CChemDomain2D::AdjustForSourceTerms(const double t, const double tstep, Writeable1DArray massadded){
	int 	 k,m,s;
	double tmpconc[MAX_SPECIES];
	double tmpflux[MAX_SPECIES];
	double vol,recharge;
//TMP DEBUG - FEREVISE

	for (s=0; s<nspecies; s++){massadded[s]=0.0;}

	int 	 nNodes=pConcGrid->GetNumNodes();
	if (nNodes==0){ExitGracefully("CChemDomain2D::AdjustForSourceTerms: zero nodes",RUNTIME_ERR);}

	//sift through areal source zones
	for (k=0; k<nNodes; k++){
		pt3D pt=pConcGrid->GetNodeLocation(k);
		cmplex z=c3Dto2D(pt);

		if (pConcGrid->GetType()==CELL_BASED){vol=pConcGrid->GetCellArea(k)*satthick[k]; }
		recharge=fabs(pAq->GetLeakage(pt,t,FROMTOP));
			
		//area sources---------------------------------------------------------
		for (m=0;m<nsources; m++){
			if (pAreaSources[m]->IsInside(z)){
				//cout <<"INSIDE SOURCE"<<endl;
				//pAreaSources[k]->GetType();
				pAreaSources[m]->GetSpecifiedMassFlux(z,t+(tstep/2.0),tmpflux);
				for (s=0; s<nspecies; s++){if (tmpconc[s]!=NO_INFLUENCE){Cnew[k][s]+=tmpflux[s]*tstep/(satthick[k]*poro[k]);}}
				
				pAreaSources[m]->GetSpecifiedConcentrations(z,t+tstep,tmpconc);
				for (s=0; s<nspecies; s++){if (tmpconc[s]!=NO_INFLUENCE){Cnew[k][s]=tmpconc[s];}}

				pAreaSources[m]->GetInitialSorbedConcentrations(z,t,tmpconc);
				for (s=0; s<nspecies; s++){if (tmpconc[s]!=NO_INFLUENCE){Csorb[k][s]=tmpconc[s];}}

				pAreaSources[m]->GetRechargeConcentrations(z,t+(tstep/2.0),tmpconc);
				for (s=0; s<nspecies; s++){if (tmpconc[s]!=NO_INFLUENCE){Cnew[k][s]+=recharge*(tmpconc[s]-Cnew[k][s])*tstep/(satthick[k]*poro[k]);}}
				
			}
		}

		//point source/sinks---------------------------------------------------------
		/*for (m=0;m<npointsourcesinks; m++){
			if (pConcGrid->IsInCell(c2Dto3D(pPointSources[m]->GetZ()),k)){

				pPointSources[m]->GetSpecifiedMassFlux(t+(tstep/2.0),tmpflux);
				for (s=0;s<nspecies; s++){Cnew[k][s]+=tmpflux[s]*tstep/(vol*poro[k]);}
				
				pPointSources[m]->GetSpecifiedConcentrations(t+tstep,tmpconc);
				for (s=0;s<nspecies; s++){if (tmpconc[s]!=no_influence){Cnew[k][s]=tmpconc[s];}}
			}
		}*/
	}

		//line sources/sinks-----------------------------------------------------------
		/*pt3D pt1,pt2,*ptint;
		int *cIds;
		int nint;
		for (m=0; m<nlinesources;m++){
			if (pConcGrid->GetMeshBoundaryIntercepts(pt1,pt2,ptint,cIds,nint)){
				for (int j=0; j<nint; j++){
				
				}
			}
		}*/
	//} //end if cell based
	//else if (pConcGrid->GetType()==NODE_BASED){ //finite element method
		
	//}
}
/************************************************************************
													TRANSPORT
*************************************************************************
	Transports contaminants using Operator splitting
	2D/3D- mostly generic with exception of 
		budgetpoly calculations, satthick in initialization of transport scheme 
-----------------------------------------------------------------------*/
void		CChemDomain2D::Transport(ofstream &PROGRESS){

	int 		i,k,s;
	double	t(0.0),simulation_end_time;
	cmplex	tmpz;
	pt3D		tmppt;
	bool		SIAdone(true);
	int 		SIAiter;
	int 		nextplottime(0);
	bool		stopped(false);
	int 		nNodes;
	double	M;
	double  LocTimeStep;

	double	Csnew 			 [MAX_SPECIES];

	double	sys_mass		 [MAX_SPECIES];  //All used for mass budget
	double	sorb_mass 	 [MAX_SPECIES];
	double	MB_error		 [MAX_SPECIES];
	double	boundary_loss[MAX_SPECIES];
	double	sink_loss 	 [MAX_SPECIES]; 
	double	rxn_change	 [MAX_SPECIES];
	double	source_gain  [MAX_SPECIES];
	double	init_mass 	 [MAX_SPECIES];
	ofstream MB;

	cout<<endl<<endl<<"_______________*_Transporting Contaminant_*______________________"<<endl<<endl;

	PROGRESS<<"transport"<<endl<<0<<endl;

	//check for bad data----------------------------------------------------
	if (pConcGrid==NULL)	 {ExitGracefully("CChemDomain2D::Transport: Grid or Mesh is NULL" 														, BAD_DATA);}
	if (pTransScheme==NULL){ExitGracefully("CChemDomain2D::Transport: TransScheme is NULL"											, BAD_DATA);}
	if (nOutputTimes<=0)	 {ExitGracefully("CChemDomain2D::Transport: No transport model output times specified", BAD_DATA);}
	if (TimeStep<=0)			 {ExitGracefully("CChemDomain2D::Transport: Invalid transport timestep specified" 		, BAD_DATA);}
	if ((DispersionType!=EULERIAN) && (pDispScheme!=NULL)){delete pDispScheme;pDispScheme=NULL;}	


	//PlotKd(); //TMP DEBUG 


	//Initialize Local Variables--------------------------------------------
	nNodes						 =pConcGrid->GetNumNodes();
	simulation_end_time=OutputTimes[nOutputTimes-1];

	//Initialize Isotherm Table (if not done previously)--------------------
	if (IsothermTable==NULL) {InitializeIsothermTable();} 

	//Initialize Concentration Grid/Read Initial Conditions-----------------
	InitializeGrid();
	if (ProgramAborted()){return;}

	if ((!stopped) && (UseInitialConditions)){
		ReadInitialConditionsFromFile();
		t=StartTime;
		for (i=0;i<nOutputTimes; i++){if (StartTime>OutputTimes[i]){nextplottime++;}}
	}
	//Update sorbed concentrations (for isotherm relationships)
	for (k=0; k<nNodes; k++){ 
		for (s=0; s<nspecies; s++){
			if (UsesIsotherm[s]){Csorb[k][s]=TranslateToSorbed(pConcGrid->GetNodeLocation(k),Ctmp[k][s],s);}
		}
	}
	if (ProgramAborted()){return;}

	//Initialize source/sink information------------------------------------
	InitializeSources(StartTime);
	if (ProgramAborted()){return;}

	//Initialize Transport/Dispersion Schemes-------------------------------
	pTransScheme->SetAdvectionParams(AdvectionType, AdvTimeStep, AdvSpaceStep);

	//pTransScheme->Initialize(ADV_AND_DISP ,poro,satthick); //For overhaul

	if (pDispScheme!=NULL){ 
		pTransScheme->Initialize(StartTime,ADVECTION_ONLY ,poro,satthick);//TMP DEBUG
		pTransScheme->SetSubScheme(pDispScheme);
		pDispScheme ->Initialize(StartTime,DISPERSION_ONLY,poro,satthick);
	}
	else{
		pTransScheme->Initialize(StartTime,ADV_AND_DISP 	,poro,satthick);
	}

	//cout <<"Initializing Reaction Schemes..."<<endl;
	//Initialize Reaction Schemes-------------------------------------------
	for (i=0; i<nRxNs; i++){
		pReactSchemes[i]->SetSoilArray(pSoilArray,nSoils);
		pReactSchemes[i]->SetSpeciesArray(pSpeciesArray,nspecies);
		pReactSchemes[i]->Initialize();
	}
	//cout <<"Writing Initial Output..."<<endl;
	
	//Plot Initial Concentrations, write initial output (if required) ------
	if (!initialize_only){
		if ((nextplottime<nOutputTimes) && (OutputTimes[nextplottime]==StartTime)){
			nextplottime++;
			PlotConcentrations			 (nextplottime);
			PlotSorbedConcs 				 (nextplottime);
			WriteConcentrationToFile (nextplottime);
			pTransScheme->WriteOutput(nextplottime);
		} 
	}
	
	if (ProgramAborted()){return;}
	for (s=0; s<nspecies; s++){
		MB_error [s]=boundary_loss[s]=sink_loss[s]=source_gain[s]=0.0;
		sys_mass [s]=pTransScheme->CalculateSystemMass(s,C,Csorb,sorb_mass[s]);
		init_mass[s]=sys_mass[s];
	}
	if (MBEnabled){ 
		MB.open("MassBalance.csv"); MB<<"time,";
		MB.precision(9);
		for (s=0; s<nspecies; s++){
			MB<<"mass "<<s+1<<", source "<<s+1<<", sink "<<s+1<<", bndry "<<s+1<<", MBerr "<<s+1<<",sorb "<<s+1;
			if (s!=nspecies-1){MB<<",";}
		}
		MB<<endl<<t<<",";
		for (s=0;s<nspecies;s++){

			//tot_mass_error=sys_mass_start+source_gain-sink_loss-boundary_loss+RxN_gain-RxN_loss-sys_mass_end
			MB<<sys_mass			[s]<<","<<source_gain 	[s]<<","<<-sink_loss		[s]<<","
				<<-boundary_loss[s]<<","<<MB_error			[s]<<","<<sorb_mass 		[s];
			if (s!=nspecies-1){MB<<",";}
		}
		MB<<endl; 	
	}
  if (initialize_only){PROGRESS<<"done"<<endl; return;}

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

	

	SIAdone=false;
	SIAiter=0;	
	int steps(0);
	LocTimeStep=TimeStep;
	//Transport Contaminant-------------------------------------------------
	do {
		if (ProgramAborted()){return;}

		LocTimeStep=min(TimeStep,OutputTimes[nextplottime]-t);

		if (!stopped){
	
			if ((steps%WriteInterval)==0){
				cout << endl << "Transport time: " <<t;
			}

			SIAiter=0;
			do //sequential iterative approach (SIA) loop-not presently used
			{ 

				//InitializeMassBalance(t);
				for (s=0; s<nspecies; s++){
					MB_error [s]=0.0;
					sys_mass [s]=pTransScheme->CalculateSystemMass(s,C,Csorb,sorb_mass[s]);
				}
				
				//Advective/Dispersive Transport of Contaminants-------------------------------------------- 
				CopyConcentrationArray(Cnew,Ctmp);
				pTransScheme->Transport(t,LocTimeStep,C,Ctmp,Cnew);	
				
				//Eulerian Dispersive Transport of Contaminants---------------------------------------------
				if ((DispersionType==EULERIAN) && (pDispScheme!=NULL))
				{
					CopyConcentrationArray(Cnew,Ctmp);
					pDispScheme->Transport(t,LocTimeStep,C,Ctmp,Cnew);  

				}
				//Adjust for source zones/points (for fully Lagrangian schemes)-----------------------------
				
				//Update sorbed concentrations (for isotherm relationships)
				for (k=0; k<nNodes; k++){ 
					for (s=0; s<nspecies; s++){
						if (UsesIsotherm[s]){Csorb[k][s]=TranslateToSorbed(pConcGrid->GetNodeLocation(k),Ctmp[k][s],s);}
					}
				}

				//UpdateMassBalance(adv/advdisp,t);
				for (s=0; s<nspecies; s++){
					M 								 =pTransScheme->CalculateSystemMass(s,Cnew,Csorb,sorb_mass[s]);
					if (MBEnabled){
						boundary_loss[s]+=pTransScheme->CalculateBorderFlux(s,C,Cnew,t,LocTimeStep);
						sink_loss 	 [s]+=pTransScheme->CalculateSinkLoss  (s,C,Cnew,t,LocTimeStep);
						source_gain  [s]+=pTransScheme->CalculateSourceGain(s,C,Cnew,t,LocTimeStep);
						MB_error		 [s]+=(M-init_mass[s])-source_gain[s]+(sink_loss[s]+boundary_loss[s]);//Cumulative error
						sys_mass		 [s] =M;
					}
				} 			

				//Reaction-----------------------------------------------------------------------
				//sift through cells, perform batch reactions at each
				for (i=0;i<nRxNs;i++){
					CopyConcentrationArray(Cnew,Ctmp);
					for (k=0; k<nNodes; k++){ 	
						tmppt=pConcGrid->GetNodeLocation(k);
						pReactSchemes[i]->React(C[k], 		Ctmp[k],	Cnew[k],
																		Csorb[k], Csorb[k], Csnew, //tmp debug
																		GetSoilType(tmppt),poro[k], LocTimeStep);
						
						for (s=0; s<nspecies; s++){
							if (!UsesIsotherm[s]){Csorb[k][s]=Csnew[s];}
						}
					}
				}
				//UpdateMassBalance(rxn, t);
				for (s=0; s<nspecies; s++){
					//RxN_mass_error=sys_mass_start+RxN_gain(pred)-RxN_loss(pred)-last_sys_mass
					M 						=pTransScheme->CalculateSystemMass(s,Cnew,Csorb,sorb_mass[s]);
					rxn_change[s]+=M-sys_mass[s];
					sys_mass	[s] =M;
				}
				
				//test for SIA convergence-not presently used
				/*for (k=0; k<nNodes; k++){ 		
					for (s=0; s<nspecies; s++){
						upperswap(SIAmaxchange,Cnew[k][s]-Csia[k][s]);
						max_conc=upperswap(max_conc,Cnew[k][s])
						SIAmaxchange/=max_conc;
					}
				}
				CopyConcentrationArray(Cnew,Csia);*/ //may wish to use relaxation 

				SIAdone=true;//TMP DEBUG
				SIAiter++;

			}while ((SIAdone==false) && (SIAiter<MAX_SIA_ITER));

			//total update of concentrations------------------------------------
			CopyConcentrationArray(Cnew,C);

			//Plot Concentrations, write Major output (if required) ------------
			if ((nextplottime<nOutputTimes) && (t+LocTimeStep>=OutputTimes[nextplottime]))
			{
				nextplottime++;
				PlotConcentrations				(nextplottime);
				PlotSorbedConcs 					(nextplottime);
				WriteConcentrationToFile	(nextplottime);
				pTransScheme->WriteOutput (nextplottime);	
			}
			

			//increment time----------------------------------------------------
			t+=LocTimeStep;	
			steps++;

			//Write Minor Output------------------------------------------------
			if ((steps%WriteInterval)==0)
			{
				if (MBEnabled)
				{
					MB<<t<<",";
					for (s=0;s<nspecies;s++){
						//tot_mass_error=sys_mass_start+source_gain-sink_loss-boundary_loss+RxN_gain-RxN_loss-sys_mass_end
						MB<<sys_mass			[s]<<","<<source_gain 	[s]<<","<<-sink_loss		[s]<<","
							<<-boundary_loss[s]<<","<<MB_error			[s]<<","<<sorb_mass 		[s];
						if (s!=nspecies-1){MB<<",";}
					}
					MB<<endl;
				}
			
				OBSCONC<<t<<",";
				for (i=0; i<nObsPoints; i++){for (s=0; s<nspecies; s++){
						OBSCONC<<GetConcentration(ObsPoints[i],t,s);
						if ((i!=nObsPoints-1) || (s!=nspecies-1)){OBSCONC<<",";}}}OBSCONC<<endl;
			
			
				PROGRESS.precision(9);
				PROGRESS<<t 							 <<" "
								<<sys_mass			[0]<<" "
								<<source_gain 	[0]<<" "
								<<-sink_loss		[0]<<" "
								<<-boundary_loss[0]<<" "
								<<MB_error			[0]<<" "
								<<sorb_mass 		[0]<<endl;
			} /* end if ((steps%WriteInterval)==0)...*/
		} /* end if (!stopped)...*/
	} while ((t<simulation_end_time) && (!stopped));

	PROGRESS<<"done"<<endl;
	cout << "...Transport calculations done"<<endl;

	OBSCONC.close();
	MB		 .close();

}
//-------------------------------------------------------------------
void CChemDomain2D::CalculateMoments(const int s, const double &t) const{ 
	double sum(0.0);
	pt3D	 pt,center(0.0);
	int 	 k;

	for (k=0; k<pConcGrid->GetNumNodes(); k++){
		pt=pConcGrid->GetNodeLocation(k);
		if (pConcGrid->GetType()==CELL_BASED){
			sum 		+=		 Cnew[k][s]*pConcGrid->GetCellArea(k)*satthick[k]*poro[k]/VolumeRatio;
			center.x =pt.x*Cnew[k][s]*pConcGrid->GetCellArea(k)*satthick[k]*poro[k]/VolumeRatio;
		}

	}
	center.x/=sum;
	center.y/=sum;
	center.z/=sum;
}