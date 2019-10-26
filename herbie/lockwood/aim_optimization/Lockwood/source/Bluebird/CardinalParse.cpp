//CardinalParse.cpp
#include "Bluebird.h"

/******************************************************************************
		CARDINAL FILE PARSING ROUTINE
-------------------------------------------------------------------------------
IMPORTANT: order of commands:
					-species and soil declarations before Isotherms
					-species declarations before sourcezone declaration
					-soiltype declarations before soilzone declaration (???)
					-concentration grid or mesh must be included
					-"Transport" must be set to on before grid will be created
					-TransportParams before DispersionParams before (ReadGrid or ReadMesh) (for Eulerian Advection)
					-TemporalWeighting before (ReadGrid or ReadMesh)
					-UpstreamWeighting before (ReadGrid or ReadMesh)
					-MMOCParams after ReadGrid
					-PrintGridDiagnostics after ReadGrid
******************************************************************************/
bool CardinalParse(char 						*filename,
									 CChemDomain2D	 *&TransDomain,
									 CAquifer 			 *&pAq,
									 CTriMesh        *pMesh){

	//parsing variables----------------------------------------------------
	char	*s[MAXINPUTITEMS];
	char	*thisname=new char	[255];

	int 	 Len,code,l(0),i;
	bool	 eof (false);
	bool	 done(false);
	bool	 noisy(parserdebug); //FOR PARSER DEBUG, turn parserdebug to true (in "MasterInclude.h")
	bool   initconc=false;

	//creation variables-------------------------------------------------
	CSpecies					 *pSpecies=NULL;
	CSoilType 				 *pSoil 	=NULL;
	CPropZone 				 *cast2 	=NULL;

	C2DTransportScheme *pTS 					=NULL;
	CStreamlineTrans	 *pStreamTrans	=NULL;
	C2DPatchSource		 *pPatchSource	=NULL;
	CRandomWalk 			 *pRandWalk 		=NULL;
	C2DAnalytic        *pAnalyticScheme=NULL;
	CMOC							 *pMOC					=NULL;
	CMMOC 						 *pMMOC 				=NULL;
	C2DPointSource		 *p2DPointSource=NULL;
	C2DFEEulerian 		 *pFEM					=NULL;
	C2DFDEulerian 		 *pFDM					=NULL;

	CRectGrid 				 *pGrid 				=NULL;

	CTimeSeries 			 *pTimeSeries[MAX_SPECIES];

	int 								NumSeries=0;
	int                 nspecies=0;

	double							aqdispersivity(0); 	 
	bool								elemdisabled;
	bool								eulerian_advect(false);
	bool								SUPG(false);
	bool								ConsistentFormulation(false);
	double							temporalweight(0.5); //Crank Nicholson by default
	double							upstreamweight(0.0); //central difference by default
	bool								DXFMeshOn(false),BNAMeshOn(false);

	char								blank[1]; blank[0]='\0';
	char								space[1]; space[0]=' ';

	TransDomain 				=new CChemDomain2D();		
	ExitGracefullyIf(TransDomain==NULL,"Parse::Transport Domain Creation- out of memory",OUT_OF_MEMORY);

	TransDomain->SetAquifer(pAq);
	TransDomain->SetLayer(pAq->GetLayer(0));//top layer by default

	for (int sp=0; sp<MAX_SPECIES; sp++){pTimeSeries[sp]=NULL;}

	ofstream BASEMAP;
	BASEMAP.open("TransportBasemap.bna");

	ifstream CDL(filename);  
	if (CDL.fail()){cout << "Cannot find file "<< filename <<endl; return false;}

	//-----------------------------------------------------------------
	cout << "Parsing File " << filename <<"..."<<endl;
	cout << "------------------------------------------------------------"<<endl;

	//--sift through file-----------------------------------------------
	while ((!TokenizeLine(CDL,s,Len)) && (!eof)){

		l++;	
		if (noisy){ cout << "reading line " << l << ": ";}

		if			(Len==0)																														{code=-1;}
		/*------------BASIC OPTIONS (100-200)------------------------------------------------*/
		else if ((!strcmp(s[0],"SpeciesList"					 )) ) 												{code=100;}
		else if ((!strcmp(s[0],"SoilList" 						 )) ) 												{code=101;}
		else if ((!strcmp(s[0],"OutputTimes"					 )) ) 												{code=102;}
		else if ((!strcmp(s[0],"BackgroundConcentrations")) ) 											{code=103;}
		else if ((!strcmp(s[0],"AmbientConcentrations" )) )													{code=103;}
		else if ((!strcmp(s[0],"LiterConversionFactor" )) ) 												{code=104;}
		else if ((!strcmp(s[0],"PlotSorbedConcentrations" )) ) 											{code=105;}

		/*------------SOURCE ZONES (200-300)-------------------------------------------------*/
		else if ((!strcmp(s[0],"ConcDiscTimeSeries" 	 )) ) 												{code=201;}
		else if ((!strcmp(s[0],"SourceZone" 					 )) ) 												{code=202;}
		else if ((!strcmp(s[0],"InitConcZone" 				 )) ) 												{code=202;}
		else if ((!strcmp(s[0],"EllSourceZone"				 )) ) 												{code=203;}
		else if ((!strcmp(s[0],"ContaminatedRecharge"  )) ) 												{code=204;}
		else if ((!strcmp(s[0],"StreamlineStartFace"	 )) ) 												{code=205;}
		else if ((!strcmp(s[0],"PointSourceSink"			 )) ) 												{code=206;}
		else if ((!strcmp(s[0],"LinearSourceSink" 		 )) ) 												{code=207;}
		else if ((!strcmp(s[0],"ImportInitialConditions")) )												{code=208;}
		else if ((!strcmp(s[0],"2DPlaneSource"         )) ) 												{code=209;}
		else if ((!strcmp(s[0],"2DInitialPointSource"  )) ) 												{code=210;}
		else if ((!strcmp(s[0],"2DConstantPointSource" )) ) 												{code=211;}

		/*------------SOIL ZONES/INFO (300-400)----------------------------------------------*/
		else if ((!strcmp(s[0],"SoilZone" 						 )) ) 												{code=301;}
		else if ((!strcmp(s[0],"EllSoilZone"					 )) ) 												{code=302;}
		else if ((!strcmp(s[0],"HeterogeneousKd"			 )) ) 												{code=303;}
		else if ((!strcmp(s[0],"SorptionIsotherms"		 )) ) 												{code=304;}

		/*------------TRANSPORT ALGORITHM (400-500)------------------------------------------*/
		else if ((!strcmp(s[0],"TransportParams"			 )) ) 												{code=401;}
		else if ((!strcmp(s[0],"AdvectionParams"			 )) ) 												{code=402;}
		else if ((!strcmp(s[0],"DispersionParams" 		 )) ) 												{code=403;}
		else if ((!strcmp(s[0],"EffectiveParams"			 )) ) 												{code=404;}
		else if ((!strcmp(s[0],"TemporalWeight" 			 )) ) 												{code=405;}
		else if ((!strcmp(s[0],"UpstreamWeight" 			 )) ) 												{code=406;}
		else if ((!strcmp(s[0],"StreamlineUpwind" 		 )) ) 												{code=407;}		
		else if ((!strcmp(s[0],"DisableMassBalance" 	 )) ) 												{code=408;}
		else if ((!strcmp(s[0],"ConsistentFormulation" )) ) 												{code=-1;}
		else if ((!strcmp(s[0],"DiagnosticsOnly"       )) ) 												{code=410;}

		else if ((!strcmp(s[0],"ConcObsPoints"				 )) ) 												{code=450;}

		/*------------SPECIFIC TRANSPORT ALGORITHM (500-600)---------------------------------*/
		else if ((!strcmp(s[0],"StreamlineParams" 		 )) ) 												{code=501;}
		else if ((!strcmp(s[0],"MOCParams"						 )) ) 												{code=502;}
		else if ((!strcmp(s[0],"RandomWalkParams" 		 )) ) 												{code=503;}
		else if ((!strcmp(s[0],"FiniteElementParams"	 )) ) 												{code=504;}
		else if ((!strcmp(s[0],"MMOCParams" 						)) )												{code=505;}

		else if ((!strcmp(s[0],"PatchSourceParams"		 )) ) 												{code=550;}
		else if ((!strcmp(s[0],"PointSourceParams"		 )) ) 												{code=551;}

		/*------------GRID / MESH (600-700)--------------------------------------------------*/
		else if ((!strcmp(s[0],"ConcGrid" 						 )) || 
						 (!strcmp(s[0],"ConcentGrid"					 )) || 
						 (!strcmp(s[0],"ConcentrationGrid"		 )) ) 												{code=601;}
		else if ((!strcmp(s[0],"ReadFDGridFromFile" 	 )) ) 												{code=602;}
		else if ((!strcmp(s[0],"ConcMesh" 						 )) || 
						 (!strcmp(s[0],"ConcentMesh"					 )) ) 												{code=603;}
		else if ((!strcmp(s[0],"GridGenerateOnly" 		 )) ) 												{code=604;}
		else if ((!strcmp(s[0],"ReadFEMeshFromFileOld" )) ) 												{code=605;}		
		else if ((!strcmp(s[0],"ReadFEMeshFromFile" 	 )) ) 												{code=606;}
		else if ((!strcmp(s[0],"GridSpacingFunction" 	 )) ) 												{code=607;}

		else if ((!strcmp(s[0],"MeshDXFOutputOn"       )))                          {code=650;}
		else if ((!strcmp(s[0],"MeshBNAOutputOn"       )))                          {code=651;}
		else if ((!strcmp(s[0],"PrintGridDiagnostics"  )))                          {code=652;}

		else if ((!strcmp(s[0],"CationExchangeRxN"		 )) ) 												{code=900;}

		//disabled or commented out
		else if ((!strcmp(s[0],"#"										 )) || 
						 (!strcmp(s[0],"rem"									 )) || 
						 (!strcmp(s[0],"*"										 )) || 
						 (!strcmp(s[0],"&"										 )) ) 												{code=-1;}
		//end file
		else if ((!strcmp(s[0],"end"									 )) || 
						 (!strcmp(s[0],"End"									 )) || 
						 (!strcmp(s[0],"EndInput" 						 )) ) 												{code=-1;eof=true;}
		//unrecognized
		else																																				{code=-3;}

		//disable 
		if ((code>2) && (Len>1) && 
				((!strcmp(s[1],"disable")) || (!strcmp(s[1],"Disable")))) 							{elemdisabled=true; code=-2;}
		else																																				{elemdisabled=false;}

		switch(code){
		case(100)://-----------------------------------------------------------------------------------
			{//SpeciesList
			 //{string species name, double diffusion coeff, double decay, double mol.wt}x numspecies
			 //&
			 if (noisy){cout<<"Species List"<<endl;}
			 done=false;																					 eof=TokenizeLine(CDL,s,Len); l++;
				do {
					if	((Len>3) && (strcmp(s[0],"&"))){
						thisname=strcpy(thisname,blank);
						for(i=0; i<Len-4; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
						pSpecies= new CSpecies(thisname, s_to_d(s[Len-3]),s_to_d(s[Len-1]),s_to_d(s[Len-2]) );
						TransDomain->CChemDomain::AddToDomain(pSpecies); eof=TokenizeLine(CDL,s,Len); l++;
						nspecies++;
					}
					else if ((Len==1) && (!strcmp(s[0],"&"))) { 
						done=true;
					}
					else if (Len==0) {eof=TokenizeLine(CDL,s,Len); l++;}
					else {ImproperFormat(s,l); break;}
				} while ((!done) && (!eof));
			 break;
			}
		case(101)://-----------------------------------------------------------------------------------
			{//SoilsList
			 //{string soil name, double dry density}x numsoils
			 //&
			 if (noisy){cout<<"Soiltype List"<<endl;}
			 done=false;																					 eof=TokenizeLine(CDL,s,Len); l++;
				do {
					if	((Len>=1) && (strcmp(s[0],"&"))){
						thisname=strcpy(thisname,blank);
						for(i=0; i<Len-1; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
						pSoil= new CSoilType(thisname,s_to_d(s[Len-1]));
						TransDomain->CChemDomain::AddToDomain(pSoil);
																														 eof=TokenizeLine(CDL,s,Len); l++;
					}
					else if ((Len==1) && (!strcmp(s[0],"&"))) { 
						done=true;
					}
					else if (Len==0) {eof=TokenizeLine(CDL,s,Len); l++;}
					else {ImproperFormat(s,l); break;}
				} while ((!done) && (!eof));
			 break;
			}
		case(102)://-----------------------------------------------------------------------------------
			{//Transport Model Output Times
			 //string "OutputTimes" [optional int WriteInterval]
			 //{double t}x(numtimes)
			 //&
				int ntimes(0);
				int interval(1);
				double *pTimes;
				pTimes=new double [MAX_CONC_TIMES]; 
				if (noisy){cout << "Transport Model Output Times"<<endl;}
				if (Len==2){interval=s_to_i(s[1]);}
				done=false; 																				 eof=TokenizeLine(CDL,s,Len); l++;
				do {
					if	((Len==1) && (strcmp(s[0],"&"))){
						if (ntimes>=MAX_CONC_TIMES){ExitGracefully("Parse: Too many transport model output times",BAD_DATA);}
						pTimes[ntimes]=s_to_d(s[0]); ntimes++;			 eof=TokenizeLine(CDL,s,Len); l++;
					}
					else if ((Len==1) && (!strcmp(s[0],"&"))) { 
						TransDomain->SetOutputTimes(pTimes, ntimes,interval);
						done=true;
					}
					else if (Len==0) {eof=TokenizeLine(CDL,s,Len); l++;}
					else {ImproperFormat(s,l); break;}
				} while ((!done) && (!eof));	
				delete [] pTimes;
				break;
			}
		case(103)://-----------------------------------------------------------------------------------
			{//Background Concentrations
			 //string "BackgroundConcentrations" 
			 //{double C double S}x(nspecies)
			 //&
				int nspec(0);
				double pC[MAX_SPECIES];
				double pS[MAX_SPECIES]; 
				if (noisy){cout << "Transport Model Ambient Concentrations"<<endl;}
				done=false; 																				 eof=TokenizeLine(CDL,s,Len); l++;
				do {
					if (Len==2){
						if (nspec>=MAX_SPECIES){ExitGracefully("Parse: Too many ambient species concentrations",BAD_DATA);}
						pC[nspec]=s_to_d(s[0]);pS[nspec]=s_to_d(s[1]);
						nspec++;																				 eof=TokenizeLine(CDL,s,Len); l++;
					}
					else if ((Len==1) && (!strcmp(s[0],"&"))) { 
						if (nspec!=TransDomain->GetNumSpecies()){ExitGracefully("Parse: Not enough background species concentrations specified",BAD_DATA);}
						TransDomain->SetAmbientConcs(pC, pS);
						done=true;
					}
					else if (Len==0) {eof=TokenizeLine(CDL,s,Len); l++;}
					else {ImproperFormat(s,l); break;}
				} while ((!done) && (!eof));	
				break;
			}
		case(104)://-----------------------------------------------------------------------------------
			{//Liter Conversion Ratio
			 //string "LiterConversionFactor", double V [liters/L^3]
				if (noisy){cout <<"Liter Conversion Factor"<<endl;}
				if (Len==2){TransDomain->SetVolumeRatio(s_to_d(s[1])); }
				else			 {ImproperFormat(s,l);}
			 break;
			}
		case(105)://-----------------------------------------------------------------------------------
			{//PlotSorbedConcentrations
			 //string "PlotSorbedConcentrations"
				if (noisy){cout <<"Plot Sorbed Concentrations"<<endl;}
				TransDomain->PlotSorbedMass(); 
			 break;
			}
		case(201)://-----------------------------------------------------------------------------------
			{//Concentration Time Series (discrete)
			 //string "ConcDiscTimeSeries" 
			 //{double t double C} x{num time series steps}
			 //&
				if (noisy){cout << "Discrete Concentration Time Series"<<endl;}
				ExitGracefullyIf((NumSeries>=nspecies),"Parse: Inappropriate number of time series",BAD_DATA);
				pTimeSeries[NumSeries]=CDiscreteTimeSeries::Parse(CDL,l);
				NumSeries++;
				break;
			}
		case(202)://-----------------------------------------------------------------------------------
			{//Source Zone
			 //string "SourceZone"
			 //double starttime
			 //double endtime
			 //{double x double y}x(numlines+1)
			 //&
			 //{double initconc}x(numspecies)
			 //&
			 if (noisy) {cout <<"Constant Concentration Source Zone"<<endl;}
			 CAreaSource *tmp;
			 tmp=NULL;
			 tmp=CPolyAreaSource::Parse(CDL,SPECIFIED_CONC,l,nspecies);
			 if (tmp!=NULL){
				 TransDomain->AddToDomain(tmp);
				 if (NumSeries!=0){tmp->SetTimeSeries(pTimeSeries,NumSeries);NumSeries=0;}
			 }
			 break;
			}
		case(203)://-----------------------------------------------------------------------------------
			{//Elliptical Source Zone
			 //string "EllSourceZone"
			 //double starttime
			 //double endtime
			 //double x1 double y1 int soilindex	double A double B double angle (degrees)
			 //&
			 //{double initconc}x(numspecies)
			 //&
			 if (noisy) {cout <<"Elliptical Contaminant Source Zone"<<endl;}
			 CAreaSource *tmp;
			 tmp=NULL;
			 tmp=CEllAreaSource::Parse(CDL,SPECIFIED_CONC,l,nspecies);
			 if (tmp!=NULL){TransDomain->AddToDomain(tmp);}
			 break;
			}
		case(204)://-----------------------------------------------------------------------------------
			{//Contaminated Recharge Zone
			 //string "ContaminatedRecharge"
			 //double starttime
			 //double endtime
			 //{double x double y}x(numlines+1)
			 //&
			 //{double initconc}x(numspecies)
			 //&
			 if (noisy) {cout <<"Contaminated Recharge zone"<<endl;}
			 CAreaSource *tmp;
			 tmp=NULL;
			 tmp=CPolyAreaSource::Parse(CDL,SPECIFIED_CONC_RECHARGE,l,nspecies);
			 if (tmp!=NULL){TransDomain->AddToDomain(tmp);}
			 break;
			}
		case(205)://-----------------------------------------------------------------------------------
			{//Streamline Starting Face (Injection Front)
			 //string "StreamlineStartFace",double x1,y1,x2,y2, int divisions
			 if (noisy){cout <<"Streamline Starting Face"<<endl;}
			 if (Len==6){
				 CTransect *tmp;
				 tmp=NULL;
				 tmp=new CTransect(-1,TransDomain->GetLayer(),s_to_c(s[3],s[4]),s_to_c(s[1],s[2]),s_to_i(s[5]));
				 if (tmp!=NULL){CStreamline::AddStartingFace(tmp);}
			 }
			 else 		 {ImproperFormat(s,l);}
			 break;
			}
		case(206)://------------------------------------------------------------------------------------
			{//Point Source/Sink
			 //string "PointSourceSink" 
			 //double x double y double r 
			 //&
			 //[conc[s]] x nspecies
			 //&
				CPointSourceSink *tmp=NULL;
				tmp=CPointSourceSink::Parse(CDL,l,TransDomain->GetLayer(),TransDomain);
				if (tmp!=NULL){TransDomain->AddToDomain(tmp);}
				break;	 
			}
		case(207)://------------------------------------------------------------------------------------
			{//Linear Source/Sink
			 //string "PointSourceSink" 
			 //double x1 double y1 double x2 y2 
			 //&
			 //[conc[s]] x nspecies
			 //&
				CLinearSourceSink *tmp=NULL;
				tmp=CLinearSourceSink::Parse(CDL,l,TransDomain->GetLayer(),TransDomain);
				if (tmp!=NULL){TransDomain->AddToDomain(tmp);}
				break;	 
			}
		case(208)://-----------------------------------------------------------------------------------
			{//Warm start on intial conditions
			 //string "ImportInitialConditions"
				if (noisy) {cout <<"Import Initial Conditions"<<endl;}
				if (Len==1){initconc=true;}
				else			 {ImproperFormat(s,l);}
			 break;
			}
		case(209)://-----------------------------------------------------------------------------------
			{//2D Plane Source 
			 //string "2DPlaneSource" 
			 //double x1 double y1 double x2 double y2 double deteriorationrate
			 //&
			 //{double conc}x(numspecies) ->if C<0, then NO_INFLUENCE
			 //&
				if (noisy) {cout <<"2D Analytic Plane Source"<< endl;}
				if (pAnalyticScheme!=NULL)
				{
					C2DPlaneSource *pAnalPlaneSrc=NULL;
					pAnalPlaneSrc=C2DPlaneSource::Parse(CDL,SPECIFIED_CONC,TransDomain,l);
					pAnalyticScheme->AddSolution(pAnalPlaneSrc);
				}
				else{ExitGracefully("Cardinal Parse: Analytic transport style must be declared before Plane Source is created",BAD_DATA);}
			 break;
			}
		case(210)://-----------------------------------------------------------------------------------
			{//2D Initial Point Source 
			 //string "2DInitialPointSource" 
			 //double x double y
			 //&
			 //{double conc}x(numspecies) ->if C<0, then NO_INFLUENCE
			 //&
				if (noisy) {cout <<"2D Analytic Point Source (initial mass)"<< endl;}
				if (pAnalyticScheme!=NULL)
				{
					C2DPtSource *pAnalPtSrc=NULL;
					pAnalPtSrc=C2DPtSource::Parse(CDL,INITIAL_CONCENTRATION,TransDomain,l);
					pAnalyticScheme->AddSolution(pAnalPtSrc);
				}
				else{ExitGracefully("Cardinal Parse: Analytic transport style must be declared before Point Source is created",BAD_DATA);}
			 break;
			}
		case(211)://-----------------------------------------------------------------------------------
			{//2D Constant Point Source 
			 //string "2DConstantPointSource" 
			 //double x double y
			 //&
			 //{double conc}x(numspecies) ->if C<0, then NO_INFLUENCE
			 //&
				if (noisy) {cout <<"2D Analytic Point Source (specified concentration)"<< endl;}
				if (pAnalyticScheme!=NULL)
				{
					C2DPtSource *pAnalPtSrc=NULL;
					pAnalPtSrc=C2DPtSource::Parse(CDL,SPECIFIED_CONC,TransDomain,l);
					pAnalyticScheme->AddSolution(pAnalPtSrc);
				}
				else{ExitGracefully("Cardinal Parse: Analytic transport style must be declared before Point Source is created",BAD_DATA);}
			 break;
			}
		case(301)://-----------------------------------------------------------------------------------
			{//Soil Zone
			 //string "SoilZone", string name (ignored) 
			 //int soilindex 
			 //{double x double y}x(numlines+1)
			 //&
			 if (noisy) {cout <<"Soil Zone"<<endl;}
			 cast2=CPolyPropZone::Parse(CDL,l,soil_type);
			 if (cast2!=NULL){TransDomain->AddToDomain(cast2);}			
			 break;
			}
		case(302)://-----------------------------------------------------------------------------------
			{//Elliptical Soil Zone
			 //string "EllSoilZone", string name (ignored) 
			 //double x1 double y1 int soilindex	double A double B double angle (degrees)
			 //& 
			 if (noisy) {cout <<"Elliptical Soil Zone"<<endl;}
			 cast2=CEllPropZone::Parse(CDL,l,soil_type);
			 if (cast2!=NULL){
				 TransDomain->AddToDomain(cast2);
			 }			
			 break;
			}
		case(303)://-----------------------------------------------------------------------------------
			{//HeterogeneousKd Field
			 //string "HeterogeneousKd"
			 //{double x double y}x(numlines+1)
			 //&
			 //{double x double y}x(numbasis+1)
			 //&
				if (noisy){cout <<"Heterogeneous Kd Field"<<endl;}
				if (Len==1){
					CMQPolyPropZone *pZone;
					pZone=CMQPolyPropZone::Parse(CDL,l,soil_type);
					TransDomain->AddKdZone(pZone); 
				}
				else			 {ImproperFormat(s,l);}
			 break;
			}
		case(304)://-----------------------------------------------------------------------------------
			{//Sorption Isotherms
			 //string "SorptionIsotherms"
			 //{int soilindex int species}x{numisotherms} (not all combinations neccesary) 
			 //&
			 if (noisy) {cout <<"Sorption Isotherms"<<endl;}
			 if (Len==1){ 		 
				CIsotherm *iso;
				int sp(0),k(0);
				TransDomain->InitializeIsothermTable();
				done=false; 																				 eof=TokenizeLine(CDL,s,Len); l++;
				do {																 
					if (Len>=3){
						k =s_to_i(s[0]);
						sp=s_to_i(s[1]);
						if ((Len==4) && (!strcmp(s[2],"Linear"))){
							iso=new CLinearIsotherm(s_to_d(s[3]));
							TransDomain->AddIsotherm(k,sp,iso); 					 eof=TokenizeLine(CDL,s,Len); l++;
						}
						else if ((Len==5) && (!strcmp(s[2],"Freundlich"))){
							iso=new CFreundlichIsotherm(s_to_d(s[3]),s_to_d(s[4]));
							cout <<"Added Freundlich!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
							TransDomain->AddIsotherm(k,sp,iso); 					 eof=TokenizeLine(CDL,s,Len); l++;						
						}
						else if ((Len==6) && (!strcmp(s[2],"Langmuir"))){
							iso=new CLangmuirIsotherm(s_to_d(s[3]),s_to_d(s[4]),s_to_d(s[5]));
							TransDomain->AddIsotherm(k,sp,iso); 					 eof=TokenizeLine(CDL,s,Len); l++;						
						}
						else if ((Len==5) && (!strcmp(s[2],"PDM"))){
							iso=new CIsotherm();//TMP DEBUG
							TransDomain->AddIsotherm(k,sp,iso); 					 eof=TokenizeLine(CDL,s,Len); l++;						
						}
						else if ((Len==3) && (!strcmp(s[2],"NonEquilibrium"))){
							//iso=new CNonEqIsotherm();
							//TransDomain->AddIsotherm(k,sp,iso); 					 eof=TokenizeLine(CDL,s,Len); l++;						
						}
						else{
							ExitGracefully("CardinalParse: Incorrect format of isotherm table",BAD_DATA);
						}
					}
					else if ((Len==1) && (!strcmp(s[0],"&"))) { 
						done=true;
					}
					else if (Len==0) {																 eof=TokenizeLine(CDL,s,Len); l++;}
					else {ImproperFormat(s,l); break;}
				} while ((!done) && (!eof));	
				break;
			 }
			 else {ImproperFormat(s,l);}
			 break;
			}
		case(401)://-----------------------------------------------------------------------------------
			{//Transport Parameters
			 //string "TransportParameters", string transtype, double timestep
				if (noisy){cout <<"Transport Params"<<endl;}

				if (Len==3){
					if			(!strcmp(s[1],"RandomWalk"))		{pTS=pRandWalk		  =new CRandomWalk 		 (TransDomain);}
					else if (!strcmp(s[1],"Streamline"))		{pTS=pStreamTrans   =new CStreamlineTrans(TransDomain);}
					else if (!strcmp(s[1],"MMOC"))					{pTS=pMMOC				  =new CMMOC 					 (TransDomain);}
					else if (!strcmp(s[1],"MOC")) 					{pTS=pMOC 				  =new CMOC						 (TransDomain);}
					else if (!strcmp(s[1],"2DPatchSource")) {pTS=pPatchSource   =new C2DPatchSource	 (TransDomain);}
					else if (!strcmp(s[1],"2DPointSource")) {pTS=p2DPointSource =new C2DPointSource	 (TransDomain);}
					else if (!strcmp(s[1],"Analytic"))      {pTS=pAnalyticScheme=new C2DAnalytic     (TransDomain);}
					else if (!strcmp(s[1],"Eulerian"))			{eulerian_advect=true;}
					else																		{pTS= 						NULL;}
					
					if (!eulerian_advect){TransDomain->SetTransportScheme(pTS);}
					TransDomain->SetTransportTimeStep(s_to_d(s[2]));
				}
				else {ImproperFormat(s,l);}
			 break;
			}		
		case(402)://-----------------------------------------------------------------------------------
			{//Advection Parameters
			 //string "AdvectionParameters", string advtype, double timestep, double spacestep 
				if (noisy){cout <<"Advection Params"<<endl;}
				advtype Atype(ADAPTIVE_TIME_STEP);//default value
				if (Len==2){
					if			(!strcmp(s[1],"variable")){Atype=ADAPTIVE_TIME_STEP;}
					else if (!strcmp(s[1],"Variable")){Atype=ADAPTIVE_TIME_STEP;}
					TransDomain->SetAdvectionParams(Atype,G_MIN_TIME_STEP,G_MIN_STEP); //highest possible precision
				}
				else if (Len==4){
					if			(!strcmp(s[1],"constant"))	{Atype=CONSTANT_TIME_STEP;}
					else if (!strcmp(s[1],"Constant"))	{Atype=CONSTANT_TIME_STEP;}
					else if (!strcmp(s[1],"variable"))	{Atype=ADAPTIVE_TIME_STEP;}
					else if (!strcmp(s[1],"Variable"))	{Atype=ADAPTIVE_TIME_STEP;}
					else if (!strcmp(s[1],"adaptive"))	{Atype=ADAPTIVE_TIME_STEP;}
					else if (!strcmp(s[1],"Adaptive"))	{Atype=ADAPTIVE_TIME_STEP;}
					else if (!strcmp(s[1],"ConstSpace")){Atype=CONSTANT_SPACE_STEP;}
					else if (!strcmp(s[1],"constspace")){Atype=CONSTANT_SPACE_STEP;}

					TransDomain->SetAdvectionParams(Atype,s_to_d(s[2]),s_to_d(s[3]));
				}
				else {ImproperFormat(s,l);}
			 break;
			}
		case(403)://-----------------------------------------------------------------------------------
			{//Dispersion Parameters/Dispersivities
			 //string "DispersionParams", Dispersion Type, double alpha_l, [double alpha TH], [double alphaTV]
				double aL(0.0),aTH(0.0),aTV(0.0);
				if (noisy) {cout<<"Dispersion Params"<<endl;}
				disptype Dtype;
				Dtype=EULERIAN; //Default
				if			(Len>=2){
					if			(!strcmp(s[1],"Eulerian"	)){Dtype=EULERIAN;}
					else if (!strcmp(s[1],"eulerian"	)){Dtype=EULERIAN;}
					else if (!strcmp(s[1],"Lagrangian")){Dtype=LAGRANGIAN;}
					else if (!strcmp(s[1],"lagrangian")){Dtype=LAGRANGIAN;}
					else if (!strcmp(s[1],"Analytic")  ){Dtype=ANALYTIC;}
					else if (!strcmp(s[1],"analytic")  ){Dtype=ANALYTIC;}
				}
				else			 {ImproperFormat(s,l);break;}
				if			(Len<=3){aL=s_to_d(s[2]);aTH=aTV=aL;}
				else if (Len==4){aL=s_to_d(s[2]);aTH=s_to_d(s[3]);aTV=aTH;}
				else if (Len==5){aL=s_to_d(s[2]);aTH=s_to_d(s[3]);aTV=s_to_d(s[4]);}
				else			 {ImproperFormat(s,l);break;}
				TransDomain->SetDispersionParams(Dtype,aL,aTH,aTV); 

			 break;
			}
		case(404)://-----------------------------------------------------------------------------------
			{//EffectiveParams
				TransDomain->EnableEffectiveParams();
				break;
			}
		case(405)://-----------------------------------------------------------------------------------
			{//Temporal Weighting
			 //string "TemporalWeighting", double w
				if (noisy){cout <<"Temporal Weighting"<<endl;}
				if (Len==2){temporalweight=s_to_d(s[1]);}
				else			 {ImproperFormat(s,l);}
			 break;
			}
		case(406)://-----------------------------------------------------------------------------------
			{//Upstream Weighting
			 //string "UpstreamWeighting", double w
				if (noisy){cout <<"Upstream Weighting"<<endl;}
				if (Len==2){upstreamweight=s_to_d(s[1]);}
				else			 {ImproperFormat(s,l);}
			 break;
			}
		case(407)://-----------------------------------------------------------------------------------
			{//StreamlineUpwind
			 //string "StreamlineUpwind"
				if (noisy) {cout <<"Use Streamline Upwind Petrov Galerkin"<< endl;}
				if (Len==1){SUPG=true;}
				else {ImproperFormat(s,l);}
			 break;
			}
		case(408)://-----------------------------------------------------------------------------------
			{//Disable Mass Balance
			 //string "DisableMassBalance"
				if (noisy){cout <<"Disabling Mass Balance"<<endl;}
				if (Len==1){TransDomain->DisableMassBalance();}
				else			 {ImproperFormat(s,l);}
			 break;
			} 
		case(410)://-----------------------------------------------------------------------------------
			{//Dont Transport- only do grid/mesh diagnostics / initialization
			 //string "DisableMassBalance"
				if (noisy){cout <<"DiagnosticsOnly"<<endl;}
				if (Len==1){TransDomain->InitializeOnly();}
				else			 {ImproperFormat(s,l);}
			 break;
			} 
		case(450)://-----------------------------------------------------------------------------------
			{//Concentration Observation Points
			 //string "ConcObsPoints"
			 //{x y} x (NumPoints)
			 //&
			 if (noisy) {cout <<"concentration observation points"<<endl;}
			 if (Len==1){ 		 
				int npts(0);
				pt3D *pObsPoints;
				pObsPoints=new pt3D [MAX_OBS_POINTS]; 
				done=false; 																				 eof=TokenizeLine(CDL,s,Len); l++;
				do {
					if	((Len==2) && (strcmp(s[0],"&"))){
						pObsPoints[npts]=c2Dto3D(s_to_c(s[0],s[1]));		 eof=TokenizeLine(CDL,s,Len); l++;
						npts++;
						ExitGracefullyIf(npts>MAX_OBS_POINTS,"Cardinal Parse:: Maximum number of observation points exceeded",BAD_DATA);
					}
					else if ((Len==1) && (!strcmp(s[0],"&"))) { 
						TransDomain->SetObservationPts(pObsPoints, npts);
						done=true;
					}
					else if (Len==0) {																 eof=TokenizeLine(CDL,s,Len); l++;}
					else {ImproperFormat(s,l); break;}
				} while ((!done) && (!eof));	
				delete [] pObsPoints;
				break;
			 }
			 else {ImproperFormat(s,l);}
			 break;
			}
		case(501)://-----------------------------------------------------------------------------------
			{//Streamline Parameters
			 //string "Streamline Params", int numstreamlines, string UpdateType, [string 1Dsolvetyp], [int precision]
				if (noisy) {cout <<"Streamline Parameters"<< endl;}
				str_update_type type=UPDATE_DIRECT;//default
				if (pStreamTrans!=NULL){
					if (Len==3){
						if (!strcmp(s[2],"GridUpdate" 			 )){type=UPDATE_FROM_UNEVEN; }
						if (!strcmp(s[2],"InterpolatedUpdate")){type=UPDATE_DIRECT;  }
						if (!strcmp(s[2],"NoUpdate" 				 )){type=NO_UPDATE; 		 }

						pStreamTrans->SetParameters(s_to_i(s[1]),type);
					}
					else {ImproperFormat(s,l);}
				}
				else{ExitGracefully("Cardinal Parse: Transport style must be declared before Streamline parameters are assigned",BAD_DATA);}
			 break;
			}
		case(502)://-----------------------------------------------------------------------------------
			{//MOC Parameters
			 //string "MOCParams", string distribution pattern, int npartlow, int nparthigh
				if (noisy) {cout <<"MOC Parameters"<< endl;}
				if (pMOC!=NULL){
					if (Len==4){
						if (!strcmp(s[1],"Random" )){pMOC->SetParameters(RANDOM_PATTERN,s_to_i(s[2]),s_to_i(s[3]));}
						if (!strcmp(s[1],"Pattern")){pMOC->SetParameters(FIXED_PATTERN, s_to_i(s[2]),s_to_i(s[3]));}
					}
					else {ImproperFormat(s,l);}
				}
				else{ExitGracefully("Cardinal Parse: Transport style must be declared before MOC parameters are assigned",BAD_DATA);}
			 break;
			}
		case(503)://-----------------------------------------------------------------------------------
			{//Random Walk Parameters
			 //string "RandomWalkParams", int numparticles
				if (noisy) {cout <<"Random Walk Parameters"<< endl;}
				bool EPVA(true);
				if (pRandWalk!=NULL){
					if (Len==3){
						if			(!strcmp(s[2],"EPVA"		)){EPVA=true;  }
						else if (!strcmp(s[2],"NoEPVA"	)){EPVA=false;	}
						pRandWalk->SetParameters(s_to_i(s[1]),EPVA);
					}
					else {ImproperFormat(s,l);}
				}
				else{ExitGracefully("Cardinal Parse: Transport style must be declared before Random Walk parameters are assigned",BAD_DATA);}
			 break;
			}
		case(504)://-----------------------------------------------------------------------------------
			{//Finite Element Parameters
			 //string "FiniteElementParams", string gausstype, string velocity calculation type, string mass matrix type
				if (noisy) {cout <<"Finite Element Parameters"<< endl;}
				if (pFEM!=NULL){
					if (Len>=3){
						trigausspoints gt=TG_7POINT;
						bool trad=false;bool lumped=true;

						if			(!strcmp(s[1],"1a"))		 {gt=TG_1POINT; 	}
						else if (!strcmp(s[1],"Linear")) {gt=TG_1POINT; 	}
						else if (!strcmp(s[1],"Quartic")){gt=TG_3POINT; 	}
						else if (!strcmp(s[1],"3a"))		 {gt=TG_3POINT; 	}
						else if (!strcmp(s[1],"3b"))		 {gt=TG_3POINTB;	}
						else if (!strcmp(s[1],"3c"))		 {gt=TG_3ENDPOINT;}
						else if (!strcmp(s[1],"4a"))		 {gt=TG_4POINT; 	}
						else if (!strcmp(s[1],"6a"))		 {gt=TG_6POINT; 	}
						else if (!strcmp(s[1],"7a"))		 {gt=TG_7POINT; 	}
						else if (!strcmp(s[1],"Quintic")){gt=TG_7POINT; 	}
						else if (!strcmp(s[1],"9a"))		 {gt=TG_9POINT; 	}
						else if (!strcmp(s[1],"13a")) 	 {gt=TG_13POINT;	}
						else if (!strcmp(s[1],"16a")) 	 {gt=TG_16POINT;	}

						if			(!strcmp(s[2],"Traditional")){trad=true;	}
						else if (!strcmp(s[2],"Continuous")) {trad=false; }

						if (Len>=4){
							if			(!strcmp(s[3],"Lumped"		)){lumped=true;  }
							else if (!strcmp(s[3],"lumped"		)){lumped=true;  }	
							else if (!strcmp(s[3],"Consistent")){lumped=false; }
							else if (!strcmp(s[3],"NotLumped" )){lumped=false; }							
						}

						pFEM->SetParameters(gt,trad,lumped);
					}
					else {ImproperFormat(s,l);}
				}
				else{ExitGracefully("Cardinal Parse: Transport style must be declared before Finite element parameters are assigned",BAD_DATA);}
			 break;
			}
		case(505)://-----------------------------------------------------------------------------------
			{//MMOC Parameters
			 //string "MMOCParams", string PollockOn/PollockOff, int npartlow, int nparthigh 
				if (noisy) {cout <<"MMOC Parameters"<< endl;}
				if (pMMOC!=NULL){
					if (Len==4){
						if (!strcmp(s[1],"PollockOn" )){pMMOC->SetParameters(s_to_i(s[2]), pGrid);}
						if (!strcmp(s[1],"PollockOff")){pMMOC->SetParameters(s_to_i(s[2]), NULL);}
					}
					else {ImproperFormat(s,l);}
				}
				else{ExitGracefully("Cardinal Parse: Transport style must be declared before MOC parameters are assigned",BAD_DATA);}
			 break;
			}
		case(550)://-----------------------------------------------------------------------------------
			{//Patch Source Parameters
			 //string "PatchSourceParams", double xcenter, double ycenter, double width, double Co, double deterioration 
				if (noisy) {cout <<"Patch Source Parameters"<< endl;}
				if (pPatchSource!=NULL){
					if (Len==7){
						pPatchSource->SetParameters(s_to_c(s[1],s[2]),s_to_c(s[3],s[4]),s_to_d(s[5]),s_to_d(s[6]));
					}
					else {ImproperFormat(s,l);}
				}
				else{ExitGracefully("Cardinal Parse: Transport style must be declared before Patch Source parameters are assigned",BAD_DATA);}
			 break;
			}
		case(551)://-----------------------------------------------------------------------------------
			{//Point Source Parameters
			 //string "PointSourceParams", double xcenter, double ycenter, double Mass, string type
				if (noisy) {cout <<"Point Source Parameters"<< endl;}
				if (p2DPointSource!=NULL){
					sourcetype stype=INITIAL_CONCENTRATION;
					if (Len==5){
						if			(!strcmp(s[4],"Constant")){stype=SPECIFIED_CONC;}
						p2DPointSource->SetParameters(s_to_c(s[1],s[2]),s_to_d(s[3]),stype);
					}
					else {ImproperFormat(s,l);}
				}
				else{ExitGracefully("Cardinal Parse: Transport style must be declared before Point Source parameters are assigned",BAD_DATA);}
			 break;
			}
		case(601)://-----------------------------------------------------------------------------------
			{//Concentration Grid
			 //string "ConcGrid", double north double south double east double west int resolutionX {int resolutionY}
/*			 if (noisy) {cout <<"concentration grid"<<endl;}
			 if ((Len==6) || (Len==7)) {
				 if (Len==7){
					 pGrid=new CRectGrid(s_to_d(s[1]),s_to_d(s[2]),s_to_d(s[3]),s_to_d(s[4]),s_to_i(s[5]),s_to_i(s[6]));}
				 else{
					 pGrid=new CRectGrid(s_to_d(s[1]),s_to_d(s[2]),s_to_d(s[3]),s_to_d(s[4]),s_to_i(s[5]),s_to_i(s[5]));
				 }
				 if (pGrid!=NULL){
					 if (BNAMeshOn){pGrid->WriteGeometry();}
					 TransDomain->SetGrid(pGrid);
					 pFDM=new C2DFDEulerian(TransDomain,pGrid,temporalweight,upstreamweight);
					 if (pFDM!=NULL){
						 TransDomain->SetDispersionScheme(pFDM); //TBR
						 pTS->SetSubScheme(pFDM); 
					 }
					 else{
						 ExitGracefully("CardinalParse: Unable to successfully create Eulerian finite difference scheme",BAD_DATA);
					 }
				 }
				 TransDomain->SetGrid(pGrid);
			 }
			 else {ImproperFormat(s,l);}*/ //NOT USED!!!!
			 break;
			}
		case(602)://------------------------------------------------------------------------------------
			{//Concentration file Attached
			 //string "ReadFDGridFromFile"
				if (noisy) {cout <<"Concentration Grid Read..."<< endl;}
				if (Len==1){
					ifstream CONCGRID;
					CONCGRID.open("concgrid.bbg");
					int cl(0);
					pGrid=NULL;
					pGrid=CRectGrid::ReadBBGFile(CONCGRID,cl);
					CONCGRID.close();
					if (pGrid!=NULL){
						TransDomain->SetGrid(pGrid);
						if (BNAMeshOn){pGrid->WriteGeometry();}
						pFDM=new C2DFDEulerian(TransDomain,pGrid,temporalweight,upstreamweight);
						ExitGracefullyIf(pFDM==NULL,"CardinalParse: Unable to successfully create Eulerian finite difference scheme",BAD_DATA);
						if (eulerian_advect){TransDomain->SetTransportScheme (pFDM);}
						else                {
							pTS->SetSubScheme(pFEM);
							TransDomain->SetDispersionScheme(pFDM);
						}
					}
					else{
						ExitGracefully("CardinalParse: Unable to successfully create rectangular grid",BAD_DATA);
					}
				}
				else {ImproperFormat(s,l);}
				break;	 
			}
		case(603)://-----------------------------------------------------------------------------------
			{//Concentration Mesh
			 //string "ConcMesh"
			 //MaxNodes (or <=zero if default)
			 //{x y F} x (NumPoints)
			 //&
			 //{p1 p2} x (NumBoundarySegs) 
			 //&
			 //{p1 p2} x (NumInternalSegs) 
			 //&
			 /*if (noisy) {cout <<"Concentration Mesh"<<endl;}
			 if (Len>=1){ 		
				 pMesh=CDelunayMesh::DynamicParse(CDL,l,"ConcMesh");
				 if (pMesh!=NULL){
					 pMesh->WriteGeometry(DXFMeshOn,BNAMeshOn);
					 TransDomain->SetGrid(pMesh);
					 if (eulerian_advect){
						 pFEM=new C2DFEEulerian(TransDomain,(CFEMesh*)(pMesh),temporalweight,upstreamweight);
						 if (pFEM!=NULL){
							 TransDomain->SetTransportScheme(pFEM);
							 //pTS->SetSubScheme(pFEM); //not necc.
						 }
					 }
					 else{
						 pFEM=new C2DFEEulerian(TransDomain,(CFEMesh*)(pMesh),temporalweight,upstreamweight);
						 if (pFEM!=NULL){
							 TransDomain->SetDispersionScheme(pFEM);//TBR
							 pTS->SetSubScheme(pFEM);
						 }
					 }
				 }
			 }
			 else {ImproperFormat(s,l);}*/ //NO LONGER USED
			 break;
			}
		case(604)://-----------------------------------------------------------------------------------
			{//GridGenerateOnly
				break;
			}		
		case(605)://-----------------------------------------------------------------------------------
			{//ReadFEMeshFromFileOld
				if (noisy) {cout <<"Concentration Mesh Read (Old Format)"<<endl;}
				
				ifstream NODE,ELEM,SIDE,DNOD;
				NODE.open("ConcMesh.n"); //TMP DEBUG- name should be provided
				ELEM.open("ConcMesh.e");
				SIDE.open("ConcMesh.s");
				DNOD.open("ConcMesh.d");

				CFEMesh *pMesh;
				pMesh=new CFEMesh("FiniteElementMesh");
				pMesh->ReadItself(NODE,ELEM,SIDE,DNOD);

				if (pMesh!=NULL){
					pMesh->WriteGeometry(DXFMeshOn,BNAMeshOn);
					TransDomain->SetGrid(pMesh);
				}
				pFEM=new C2DFEEulerian(TransDomain,pMesh,temporalweight,upstreamweight);
				ExitGracefullyIf(pFEM==NULL,"CardinalParse: Unable to successfully create Eulerian finite element scheme",BAD_DATA);
				
				if (eulerian_advect){TransDomain->SetTransportScheme(pFEM);}
				else{
					pTS->SetSubScheme(pFEM);
					TransDomain->SetDispersionScheme(pFEM);
				}
				NODE.close();
				ELEM.close();
				SIDE.close();
				DNOD.close();
				break;
			}
		case(606)://-----------------------------------------------------------------------------------
			{//ReadFEMeshFromFile 
				if (noisy) {cout <<"Concentration Mesh Read (BTM Format)"<<endl;}
				
				ifstream BTM;
				BTM.open("ConcMesh.btm"); //TMP DEBUG- name should be provided

				CFEMesh *pMesh;
				pMesh=new CFEMesh("FiniteElementMesh");
				pMesh->ReadItself2(BTM);

				if (pMesh!=NULL){
					pMesh->WriteGeometry(DXFMeshOn,BNAMeshOn);
					TransDomain->SetGrid(pMesh);
				}

				if (eulerian_advect){
					pFEM=new C2DFEEulerian(TransDomain,pMesh,temporalweight,upstreamweight);
					if (pFEM!=NULL){TransDomain->SetTransportScheme(pFEM);}
				}
				else{
					pFEM=new C2DFEEulerian(TransDomain,pMesh,temporalweight,upstreamweight);
					pTS->SetSubScheme(pFEM);
					if (pFEM!=NULL){TransDomain->SetDispersionScheme(pFEM);}
				}
				BTM.close();
				break;
			}	
		case(607)://-----------------------------------------------------------------------------------
			{//GridSpacingFunction 
				if (noisy) {cout <<"Spacing-function specified grid generation"<<endl;}
				//string "GridSpacingFunction", string label (unused) 
				//double x_bottom_left y_bottom_left
				//double pivot_angle
				//double x-length, double y-length
	      //double minspacing, double maxspacing
				//double nX (number of columns) double nY (number of rows)
				//{double X double Fx}xNX
				//& 
				//{double Y double Fy}xNY
				//& 
				pGrid=NULL;
				pGrid=CRectGrid::ParseSpacing(CDL,l);
				if (pGrid!=NULL){
					TransDomain->SetGrid(pGrid);
					if (BNAMeshOn){pGrid->WriteGeometry();}
					if (eulerian_advect){
						pFDM=new C2DFDEulerian(TransDomain,pGrid,temporalweight,upstreamweight);
						if (pFDM!=NULL){TransDomain->SetTransportScheme(pFDM);}
					}
					else{
						pFDM=new C2DFDEulerian(TransDomain,pGrid,temporalweight,upstreamweight);
						if (pFDM!=NULL){
							TransDomain->SetDispersionScheme(pFDM);
						}
					}
				}
				break;
			}	
		case(650)://-----------------------------------------------------------------------------------
			{//Mesh Generation
			 //string "MeshDXFOutputOn"
			 DXFMeshOn=true; 
		   break;
		  }
		case(651)://-----------------------------------------------------------------------------------
			{//Mesh Generation
			 //string "MeshBNAOutputOn"
			 BNAMeshOn=true; 
		   break;
		  }
		case(652)://-----------------------------------------------------------------------------------
			{ //PrintGridDiagnostics
				//string "PrintGridDiagnostics"
				if (pFDM!=NULL){pFDM->PrintGridDiagnostics();} 
				if (pFEM!=NULL){pFEM->PrintMeshDiagnostics();}
				break;
		  }
		case(900)://-----------------------------------------------------------------------------------
			{//Cation Exchange Reaction
			 //string "CationExchangeRxN"
			 //double capacity double anionactivity
			 //{int index int valence double selectivity}x(NumCations)
			 //& 
				if (noisy) {cout <<"Cation Exchange Reaction"<< endl;}
				CReactionScheme *pRxN;
				pRxN=CCationExchange::Parse(CDL,l);
				if (pRxN!=NULL){TransDomain->AddReactionScheme(pRxN);}

				break;
			}
		case(-1)://------------------------------------------------------------------------------------
			{//commented out
				if (noisy) {cout <<"##"<< endl;}
				break;	 
			}
		case(-2)://------------------------------------------------------------------------------------
			{//disabled
				if (noisy) {cout <<"disabled"<< endl;} 
				break;	
			}
		case(-3)://------------------------------------------------------------------------------------
			{//unrecognize
				cout <<"transport input file: command unrecognized at line "<< l << endl;
				break;	 
			}
		default: {
				if (noisy) {cout << "##"<<endl;} 
				break;
			}//{cout<< "BAD BBL FILE INPUT, line "<< l+1 <<endl;return false;}

		}/*end switch*/
	}/* end while Tokenize & !eof*/



	//*******************************************************************************************
	//											 DONE PARSING- CREATE & BUILD SYSTEM
	//*******************************************************************************************

	if ((SUPG) && (pFEM!=NULL)){
		pFEM->UseSUPG();
	}
	if (initconc){TransDomain->SetInitialConditions();}
	BASEMAP.close();
	CDL.close();

	cout<<"...done parsing transport input file"<<endl<<endl;
	//delete [] thisname;

	return true;
}
