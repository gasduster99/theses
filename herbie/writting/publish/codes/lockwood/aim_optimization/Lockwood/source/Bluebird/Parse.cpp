#include "Bluebird.h"
//-----------------------------------------------------------------------
bool DynAddElement(CAnalyticElem *xptr, elementtype etype, CAnalyticElem **&pElems,elementtype *&types, int &size);
//-----------------------------------------------------------------------
bool DynAddPZone  (CPropZone *xptr, CPropZone **&pZones, int &size);
//-----------------------------------------------------------------------

/***********************************************************************
							BLUEBIRD PARSER
***********************************************************************
 reads input file filename, creates complete domain

 IMPORTANT:
 order of commands:
		-"WarmStart" must be declared before elements are created
		-Uniform flow must be input before reference point for each layer
		-reference point must be input before any elements on each layer (if warm start is used)
	  -"ResistanceLayer" and "NextLayer" must be declared before layer elements are added
		-Layers must be added from the bottom up
-----------------------------------------------------------------------*/
bool Parse(char            *filename,                                  
				   CAquifer        *&aq, 
					 CPathline      **particles,  int &numpart,     //particle vars
					 EngineCommands   &eng,
					 CTriMesh        *&pMesh){                                                  

  //parsing variables----------------------------------------------------
  char  *s[MAXINPUTITEMS];
	char  *thisname=new char  [255];

  int    Len,nlines,code,l(0),i;
	bool   eof(false);
	bool   noisy(parserdebug); //FOR PARSER DEBUG, turn parserdebug to true (in "MasterInclude.h")
  bool   solutionexists(true);

  //creation variables---------------------------------------------------
  CAnalyticElem   *cast=NULL;
	CPropZone       *cast2=NULL;

	CSingleLayer    *layer         [MAXAQLEVELS];//Single layer arrays
  CAnalyticElem  **AllElems      [MAXAQLEVELS];
	elementtype     *ElemTypes     [MAXAQLEVELS];
  CAnalyticElem  **TransientElems[MAXAQLEVELS];
	CPropZone      **PropZones     [MAXAQLEVELS];
  CFarField       *ffp           [MAXAQLEVELS];
	int              elemsparsed   [MAXAQLEVELS];
	int              transelems    [MAXAQLEVELS]; 
  int              numpzones     [MAXAQLEVELS];

	CAquiclude      *aquicludes    [MAXAQLEVELS];//Aquiclude arrays
	CAnalyticElem  **AllLeakElems  [MAXAQLEVELS];
  CPropZone			 **LeakPropZones [MAXAQLEVELS];
	int              leakelems     [MAXAQLEVELS];
	int              numlpzones    [MAXAQLEVELS];

	CMultiLayer     *mlayer        [MAXAQLEVELS];//Multilayer arrays
  CmlAnalyticElem**AllMLElems    [MAXAQLEVELS];
	mlelementtype   *MLElemTypes   [MAXAQLEVELS];
	CmlPropZone    **MLPropZones   [MAXAQLEVELS];
  CmlFarField     *MLffp         [MAXAQLEVELS];
	int              mlnumlevels   [MAXAQLEVELS];
  int							 mlelemsparsed [MAXAQLEVELS];
	int              nummlpzones   [MAXAQLEVELS];

	int              numlevels;

  int              defaultprec, blockprec;
	bool						 elemdisabled;

  double					 aqbase        [MAXAQLEVELS];
	double					 aqthickness   [MAXAQLEVELS];
	double					 aqporosity    [MAXAQLEVELS];  
	double					 aqconductivity[MAXAQLEVELS];

  cmplex					 FFQxy				 [MAXAQLEVELS];
	cmplex           FFgradient;
	cmplex					 FFzref;
  double					 FFrefhead;

  int              cl(0);
	bool             singlelayer(true);

	int              SuperNest;    
	bool             blocknested(true);
  bool             blocksolve,blockimplicit,done;
 	char             blank[1]; blank[0]='\0';
	char             space[1]; space[0]=' ';

	spc_type         spacetype(USER_SPECIFIED);
	double           spacemin(0),spacemax(ALMOST_INF),spacerad(1),spaceparam(0);

  int L;
  for (L=0;L<MAXAQLEVELS; L++) {
		layer					[L]=NULL;
    AllElems			[L]=NULL;
    TransientElems[L]=NULL;
		PropZones     [L]=NULL;
		ffp           [L]=NULL;
		elemsparsed		[L]=0;
		numpzones			[L]=0;

		aquicludes    [L]=NULL;
		AllLeakElems  [L]=NULL;
		LeakPropZones [L]=NULL;
		leakelems     [L]=0;
		numlpzones    [L]=0;

		mlayer        [L]=NULL;
    AllMLElems    [L]=NULL;
	  MLElemTypes   [L]=NULL;
	  MLPropZones   [L]=NULL;
		MLffp         [L]=NULL;
		mlnumlevels   [L]=1;
		mlelemsparsed	[L]=0;
		nummlpzones		[L]=0;

		transelems    [L]=0;
	}

	//Create Layer 0 -others are added as layers are added ("NextLayer" Command, code=205)
	//-----------------------------------------------------------------
	layer[0]=new CSingleLayer();
  if (layer[0]==NULL){ExitGracefully("Parse::Layer Creation- out of memory",OUT_OF_MEMORY);}

	mlayer[0]=new CMultiLayer();
  if (mlayer[0]==NULL){ExitGracefully("Parse::Multilayer Creation- out of memory",OUT_OF_MEMORY);}

  //default values---------------------------------------------------
  for (L=0;L<MAXAQLEVELS; L++) {
		//for (l=0;l<MAX_MLAYERS;l++){
			aqbase        [L]=0+(double)(L)*100; 
			aqthickness   [L]=100.0;	
			aqconductivity[L]=0.001;
			aqporosity    [L]=0.3;	
			FFQxy         [L]=0; //TMP DEBUG- to be removed
		//}
	}
	FFgradient				 =0.0;
	FFzref						 =0.0; 
	FFrefhead					 =50.0;	

	numlevels					 =1;

	defaultprec				 =3; 
  CAnalyticElem  ::SetDefaultPrecision(defaultprec);
	CmlAnalyticElem::SetDefaultPrecision(defaultprec);

	SuperNest         =1; 
  blocksolve        =false; 
	blockimplicit     =true;
	blockprec         =defaultprec;

  eng.explicitsolve =false;
	eng.writesol      =true;
	eng.writeout      =true;
	eng.writesolvetime=false;
  eng.solve         =true;
	eng.grid          =true;
	eng.track         =true; 
	eng.transport     =false;
	eng.transient     =false;
	eng.warm          =false;
	eng.warmtime      =0.0;
	eng.debug         =false;
	eng.socket        =false;
	eng.interpolate   =false;
	eng.meshgenerate  =false;
	eng.DXFMeshOn     =false;
	eng.BNAMeshOn     =false;
	eng.obsfileexists =false;

	eng.Grid.head     =true;
	eng.Grid.stream   =true;
	eng.Grid.leak     =true;
	eng.Grid.QxQy     =false;
	eng.Grid.GxGy     =false;
	eng.Grid.GxGynum  =false;
	eng.Grid.pot      =false;
	eng.Grid.cond     =false;
	eng.Grid.base     =false;
	eng.Grid.satthick =false;
	eng.Grid.thick    =false;
	eng.Grid.top      =false;
	eng.Grid.conduct  =false;
	eng.Grid.elem     =false;
	eng.Grid.saltH20  =false;
	eng.Grid.vxvy     =false;
	eng.Grid.effvel   =false;
	eng.Grid.curl     =false;
	eng.Grid.DxxDx    =false;
	eng.Grid.DxyDx    =false;
	eng.Grid.VmagDer  =false;
	eng.Grid.vder     =false;
	eng.Grid.effvelnum=false;
	eng.Grid.res=50;
	eng.Grid.W.n=100; eng.Grid.W.s=0;
	eng.Grid.W.e=100; eng.Grid.W.w=0;
	eng.nTimes=1;
	eng.outtimes[0]=0;
	eng.timestep=ALMOST_INF;
	numpart=0; 


 	noisy=true;
	CSingleLayer::fresh=true;

  //-----------------------------------------------------------------
  cout << "Parsing Input File " << filename <<"..."<<endl;
	cout << "------------------------------------------------------------"<<endl;

  //Input/Output files-----------------------------------------------
  ofstream BASEMAP;
  BASEMAP.open("Elements.bna");

  ifstream BBD(filename);  
  if (BBD.fail()){cout << "Cannot find file "<<filename <<endl; return false;}

  ifstream SOLUTION("solution.bbs");  
  if (SOLUTION.fail()){cout << "Cannot find solution file solution.bbs" <<endl;solutionexists=false;} 

	//--sift through file-----------------------------------------------
	while ((!TokenizeLine(BBD,s,Len)) && (!eof)){

		l++;  
		if (noisy){ cout << "reading line " << l << ": ";}

		/*assign code for switch statement
		--------------------------------------------------------------------------------------
		<100         : ignored/special
		0   thru 100 : Elements/Element Info
		100 thru 200 : Solve Type Info
		200 thru 300 : Aquifer Info
		300 thru 400 : FarField Info
		400 thru 500 : Solve Info
		500 thru 600 : Grid Info
		600 thru 700 : Track info, particles
		700 thru 800 : Surface Water Info
		800 thru 900 : Contaminant Transport Info
		--------------------------------------------------------------------------------------
		*/
		code=-3; //not recognized is default

		if (code==-3){
    if      (Len==0)                                                            {code=-1;}
		//---------------------ELEMENTS        ------------------------------------------------
		if      ((!strcmp(s[0],"Well"									 )) || 
						 (!strcmp(s[0],"well"									 )) )													{code=4;  }
		else if ((!strcmp(s[0],"HSWell"								 )) || 
			       (!strcmp(s[0],"hswell"								 )) )													{code=5;  }
    else if ((!strcmp(s[0],"head"						 		   )) || 
			       (!strcmp(s[0],"Head"									 )) )													{code=6;  }
		else if ((!strcmp(s[0],"Inhom"								 )) || 
			       (!strcmp(s[0],"Inhomogeneity"				 )) ||
			       (!strcmp(s[0],"inhom"								 )) || 
						 (!strcmp(s[0],"inhomogeneity"				 )) )													{code=7;  }
		else if ((!strcmp(s[0],"leaky"								 )) || 
						 (!strcmp(s[0],"LeakyWall"						 )) )													{code=8;  }
		else if ((!strcmp(s[0],"drain"								 )) ||  
						 (!strcmp(s[0],"Drain"								 )) )													{code=9;  }
		else if ((!strcmp(s[0],"stage"								 )) || 
						 (!strcmp(s[0],"River"                 )) ||
						 (!strcmp(s[0],"river"                 )) ||
						 (!strcmp(s[0],"Stage"								 )) )													{code=10; }
    else if ((!strcmp(s[0],"areasink"							 )) || 
						 (!strcmp(s[0],"AreaSink"							 )) )												  {code=11; }
    else if ((!strcmp(s[0],"CircularInh"					 )) ||
						 (!strcmp(s[0],"CircularInhomogeneity" )) ||
						 (!strcmp(s[0],"CirInhom"							 )) )													{code=12; }
		else if ((!strcmp(s[0],"vortex"								 )) || 
						 (!strcmp(s[0],"Vortex"								 )) )													{code=13; }
		else if ((!strcmp(s[0],"dipole"								 )) || 
						 (!strcmp(s[0],"Dipole"								 )) )													{code=14; }
    else if ((!strcmp(s[0],"Qnorm"								 )) || 
						 (!strcmp(s[0],"qnorm"								 )) ||
						 (!strcmp(s[0],"NormalDischarge"			 )) )													{code=15; }
    else if ((!strcmp(s[0],"CirLake"							 )) || 
						 (!strcmp(s[0],"cirlake"							 )) )		   										{code=16; }
	  else if ((!strcmp(s[0],"BaseInhom"						 )) || 
						 (!strcmp(s[0],"BInhom"								 )) )		   										{code=17; }
		else if ((!strcmp(s[0],"Reservoir"						 )) )		   		                {code=18; }
		else if ((!strcmp(s[0],"HWell"								 )) )		   		                {code=19; }
		else if ((!strcmp(s[0],"VarStage"							 )) )		   	                  {code=20; }
		else if ((!strcmp(s[0],"Extraction"						 )) || 
						 (!strcmp(s[0],"extraction"						 )) )													{code=21; }
		else if ((!strcmp(s[0],"DryWell"							 )) || 
						 (!strcmp(s[0],"drywell"							 )) )									        {code=22; }
		else if ( !strcmp(s[0],"SmoothInhom"					 )) 									        {code=23; }
		else if ((!strcmp(s[0],"EllipticalInhom"			 )) || 
						 (!strcmp(s[0],"EllipticalInhomogeneity")) )									      {code=24; }
		else if ((!strcmp(s[0],"EllLake"			         )) || 
						 (!strcmp(s[0],"EllipticalLake"        )) )									        {code=25; }
		else if ((!strcmp(s[0],"SurfaceDrain"			     )) || 
						 (!strcmp(s[0],"SeepageFace"           )) )									        {code=26; }

		else if ((!strcmp(s[0],"ConductInhom"          )) || 
						 (!strcmp(s[0],"conductinhom"					 )) )													{code=70; }
		else if ( !strcmp(s[0],"ConductHole"           ))                           {code=71; } 

		//-------------------MULTILAYER ELEMENTS-----------------------------------------------
		else if ( !strcmp(s[0],"ML_Well"               ))                           {code=75; } 
		else if ( !strcmp(s[0],"ML_CirInhom"           ))                           {code=76; } 
		else if ( !strcmp(s[0],"ML_EllipticalInhom"    ))                           {code=77; }
		else if ( !strcmp(s[0],"ML_Inhomogeneity"			 ))														{code=78; }
		else if ( !strcmp(s[0],"ML_Head"							 ))														{code=79; }
		else if ( !strcmp(s[0],"ML_LeakyWall"					 ))														{code=80; }
		else if ( !strcmp(s[0],"ML_Drain"							 ))														{code=81; }
		else if ( !strcmp(s[0],"ML_RechargeZone"			 ))														{code=82; }
		else if ( !strcmp(s[0],"ML_LeakageZone"				 ))														{code=83; }

		else if ((!strcmp(s[0],"BoxOfInhomCells"			 )) )                         {code=90; }
		}

		//---------------------SOLVE OPTIONS  -------------------------------------------------
		if (code==-3){
		if      ((!strcmp(s[0],"SolveOnly"             )) || 
			       (!strcmp(s[0],"solveonly")))			                                  {code=101;} 
    else if ((!strcmp(s[0],"GridOnly"              )) || 
			       (!strcmp(s[0],"gridonly"              )) )				                  {code=102;}
    else if ((!strcmp(s[0],"TrackOnly"             )) || 
			       (!strcmp(s[0],"trackonly"             )) )												  {code=103;}
    else if ( !strcmp(s[0],"SolveAndGrid"          ))														{code=104;}
    else if ( !strcmp(s[0],"SolveAndTrack"				 ))														{code=105;}
    else if ( !strcmp(s[0],"SolveGridTrack"        ))														{code=106;}
    else if ((!strcmp(s[0],"WarmStart"             )) ||
						 (!strcmp(s[0],"warm"                  )) )													{code=107;}
	  else if ((!strcmp(s[0],"Superblock"            )) ||
			       (!strcmp(s[0],"Superblocks"           )) ||
						 (!strcmp(s[0],"superblock"            )) ||
						 (!strcmp(s[0],"superblocks"           )) )													{code=108;}
		else if ( !strcmp(s[0],"GridWhileSolve"        ))														{code=109;}
		else if ( !strcmp(s[0],"TransportOnly"         ))														{code=110;}
		else if ( !strcmp(s[0],"SolveGridTransport"    ))														{code=111;}
		else if ( !strcmp(s[0],"SolveAndTransport"     ))														{code=112;}
		else if ( !strcmp(s[0],"TaylorCir"             ))														{code=113;}
		else if ( !strcmp(s[0],"Transient"             ))														{code=114;}
		else if ( !strcmp(s[0],"TransientTimes"        ))														{code=115;}
		else if ( !strcmp(s[0],"Relax"                 ))				          					{code=116;}
		else if ( !strcmp(s[0],"FullyExplicit"         ))				          					{code=117;}
		else if ( !strcmp(s[0],"SolveFarFieldEveryTime"))				          				  {code=118;}
		else if ( !strcmp(s[0],"WriteInterval"         ))				          					{code=119;}
		else if ( !strcmp(s[0],"SocketOn"              ))				          					{code=120;}
		else if ( !strcmp(s[0],"FullAnalyze"           ))				          					{code=121;}
		else if ( !strcmp(s[0],"ClassicSuperblocks"    ))				          					{code=122;}
		else if ( !strcmp(s[0],"Formatted"             ))				          					{code= -1;}
		else if ( !strcmp(s[0],"AreaSinkRelax"         ))				          					{code=124;}
		else if ( !strcmp(s[0],"WriteSolveTime"        ))				          					{code=125;}
		}
		//---------------------LAYER OPTIONS---------------------------------------------------
		if (code==-3){
		if      ( !strcmp(s[0],"BaseAndThickness"      ))														{code=201;}
    else if ( !strcmp(s[0],"Conductivity"          ))														{code=202;}
    else if ( !strcmp(s[0],"Porosity"              ))														{code=203;}
		else if ( !strcmp(s[0],"MultiLayer"            ))														{code=204;}
		else if ( !strcmp(s[0],"NextLayer"             ))										  			{code=205;}
	  else if ( !strcmp(s[0],"Aquiclude"             ))														{code=207;}
		else if ( !strcmp(s[0],"PorosityZone"          ))                           {code=208;}
		else if ( !strcmp(s[0],"SingleLayer"           ))														{code=209;}

		else if ( !strcmp(s[0],"Conductivities"        ))														{code=210;}
		else if ( !strcmp(s[0],"Thicknesses"					 ))														{code=211;}
		else if ( !strcmp(s[0],"Porosities"						 ))														{code=212;}
		else if ( !strcmp(s[0],"MultiLayerProperties"	 ))														{code=213;}
		else if ( !strcmp(s[0],"KAndThicknessAndPorosity" ))												{code=213;}
		}

		//---------------------FAR FIELD    ---------------------------------------------------
		if (code==-3){
    if      ( !strcmp(s[0],"Uniform"               ))														{code=301;}
    else if ( !strcmp(s[0],"Gradient"							 ))														{code=302;}
    else if ( !strcmp(s[0],"ReferencePoint"				 ))														{code=303;}
    else if ( !strcmp(s[0],"NetExtraction"				 ))														{code=304;}
		}

		//---------------------ITERATIVE SOLVE OPTIONS-----------------------------------------
		if (code==-3){
	  if      ( !strcmp(s[0],"Precision"						 ))														{code=401;}
    else if ( !strcmp(s[0],"Silent"								 ))														{code=402;}
    else if ( !strcmp(s[0],"Noisy"                 ))														{code=-1; }
    else if ( !strcmp(s[0],"ColdStart"             ))														{code=-1; }
    else if ( !strcmp(s[0],"MaximumSolveIterations"))													  {code=406;}
    else if ( !strcmp(s[0],"MinimumSolveIterations"))													  {code=407;}
    else if ( !strcmp(s[0],"RelaxationCoefficient" ))														{code=-1; }
    else if ( !strcmp(s[0],"LocalSolveIterations"  ))														{code=-1; }
    else if ( !strcmp(s[0],"PotentialTolerance"		 ))														{code=408;}
    else if ( !strcmp(s[0],"CompareHeadOn"				 ))														{code=409;}
    else if ( !strcmp(s[0],"CompareHeadOff"				 ))													  {code=-1; }
    else if ( !strcmp(s[0],"Inverse"							 ))														{code=-1; }
    else if ( !strcmp(s[0],"WriteSolutionOn"			 ))														{code=413;}
    else if ( !strcmp(s[0],"WriteSolutionOff"			 ))														{code=414;}
    else if ( !strcmp(s[0],"WriteOutputOff"				 ))														{code=415;}
		else if ( !strcmp(s[0],"BluebirdSolution"      ))				          					{code= -1;}
		}

		//---------------------GRIDDING/OUTPUT OPTIONS-----------------------------------------
    if (code==-3){
		if      ( !strcmp(s[0],"Window"                ))                           {code=501;}
    else if ( !strcmp(s[0],"Window2"							 ))                           {code=502;}
    else if ( !strcmp(s[0],"PlotQxQy"							 ))                           {code=503;}
    else if ( !strcmp(s[0],"PlotBaseThickness"		 ))                           {code=504;}
    else if ( !strcmp(s[0],"PlotCond"		           ))                           {code=505;}
    else if ( !strcmp(s[0],"PlotElement"		       ))                           {code=506;}
    else if ( !strcmp(s[0],"AnalysisFace"		       ))                           {code=507;}
    else if ( !strcmp(s[0],"ZoneBudget"		         ))                           {code=508;}
		else if ( !strcmp(s[0],"PlotGxGy"		           ))                           {code=509;}
		else if ( !strcmp(s[0],"PlotVxVy"		           ))                           {code=510;}
		else if ( !strcmp(s[0],"PlotEffVelocity"		   ))                           {code=511;}
		else if ( !strcmp(s[0],"PlotCurl"		           ))                           {code=512;}
		else if ( !strcmp(s[0],"PlotDxxDxDyyDy"		     ))                           {code=513;}
		else if ( !strcmp(s[0],"PlotGxGyNumerical"		 ))                           {code=514;}
		else if ( !strcmp(s[0],"PlotVelocityMagnitudeDerivatives"))                 {code=515;}
		else if ( !strcmp(s[0],"PlotVelocityDerivatives"))                          {code=516;}
		else if ( !strcmp(s[0],"PlotDxyDxDxyDy"        ))                           {code=517;}
		else if ( !strcmp(s[0],"PlotEffVelNum"         ))                           {code=518;}
		else if ( !strcmp(s[0],"PlotSatThick"          ))                           {code=519;}

		else if ( !strcmp(s[0],"ObservationFileExists" ))				          					{code=520;}
		}
		//---------------------TRACKING OPTIONS------------------------------------------------
		if (code==-3){
    if			((!strcmp(s[0],"ForwardTracing"        )) )                         {code=-1; }//{code=601;}
    else if ((!strcmp(s[0],"ReverseTracing"        )) )                         {code=-1; }//{code=602;}
    else if ((!strcmp(s[0],"Time"                  )) )                         {code=603;}
    else if ((!strcmp(s[0],"Particle"              )) ||
			       (!strcmp(s[0],"particle"              )) )                         {code=604;}
    else if ((!strcmp(s[0],"CaptureZone"           )) || 
						 (!strcmp(s[0],"CircleOfParticles"     )) )                         {code=605;}
    else if ((!strcmp(s[0],"CaptureAll"            )) )													{code=-1; }//{code=606;}
    else if ((!strcmp(s[0],"BoxOfParticles"        )) )                         {code=607;}
    else if ((!strcmp(s[0],"TrackDuration"				 )) )                         {code=608;}
		else if ((!strcmp(s[0],"TrackResolution"			 )) )                         {code=609;}
		else if ((!strcmp(s[0],"TrackPrecision"				 )) )                         {code=610;}
		else if ((!strcmp(s[0],"InterpolationGridFile" )) )                         {code=611;}
		else if ((!strcmp(s[0],"InterpolationGrid"		 )) )                         {code=612;}
    else if ((!strcmp(s[0],"Particle3D"            )) )                         {code=613;}
		else if ((!strcmp(s[0],"PollocksMethodOn"      )) )                         {code=615;}
		}

		//---------------------SURFACE WATER FEATURES------------------------------------------
		if (code==-3){
    if      ((!strcmp(s[0],"Headwater"             )) || 
						 (!strcmp(s[0],"HDWater"               )))                          {code=701;}
    else if ((!strcmp(s[0],"Coastal"               )))                          {code=702;}
    else if ((!strcmp(s[0],"Pipeline"              )))                          {code=703;}

		//---------------------MESH GENERATION-------------------------------------------------
    else if ((!strcmp(s[0],"MeshGenerateOnly"      )))                          {code=801;}
    else if ((!strcmp(s[0],"MeshBoundaries"        )))                          {code=802;}
		else if ((!strcmp(s[0],"MeshDXFOutputOn"       )))                          {code=803;}
		else if ((!strcmp(s[0],"MeshBNAOutputOn"       )))                          {code=804;}
		else if ((!strcmp(s[0],"CourantMeshGenerate"   )))                          {code=805;}

		//---------------------TEMPORARY-------------------------------------------------------
    else if ((!strcmp(s[0],"Voronoi"               )))                          {code=901;}

		//disabled, commented out, or recognized but not important
    else if ((!strcmp(s[0],"#"										 )) || 
			       (!strcmp(s[0],"rem"									 )) || 
						 (!strcmp(s[0],"*"										 )) || 
						 (!strcmp(s[0],"&"										 )) )													{code=-1;}
		//end file
		else if ((!strcmp(s[0],"end"									 )) || 
						 (!strcmp(s[0],"End"									 )) || 
			       (!strcmp(s[0],"EndInput"							 )) )							            {code=-1;eof=true;}
		else if ((!strcmp(s[0],"Debug"								 )) || 
						 (!strcmp(s[0],"debug"								 )) )												  {code=-10;}
		//unrecognized
    else																																			  {code=-3;}
		}

		//-------------------------------------------------------------------------------------
		//disable (not currently implemented properly-should read coeff, but not do anything
		if ((code>2) && (code<100) && (Len>1) && 
			  ((!strcmp(s[1],"disable")) || (!strcmp(s[1],"Disable"))))					      {elemdisabled=true; if (noisy) {cout <<"disabled"<< endl;} }
		else if (code!=-3)                                                          {elemdisabled=false;}
		//-------------------------------------------------------------------------------------		

    switch(code){
	  case(4):  //------------------------------------------------------------------------------------
		  {//discharge well
       //string "Well", string name 
		   //double x double y double Q {double r (ignored)}
		   //&
				ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer well in multilayer aquifer",BAD_DATA);
        thisname=strcpy(thisname,blank);
			  for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
				cast=CDischargeWell::Parse(BBD,l,layer[cl],thisname);
				if (cast!=NULL){
					if (!elemdisabled){
            if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						DynAddElement(cast,GIVEN_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					}
				}
		   break;
		  }
    case(5): //------------------------------------------------------------------------------------
		  {//head specified well
       //string "HSWell", string name
		   //double x double y double Head double r
       //&  
				ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer well in multilayer aquifer",BAD_DATA);
        thisname=strcpy(thisname,blank);
			  for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
				cast=CHeadSpecifiedWell::Parse(BBD,l,layer[cl],thisname);
				if (cast!=NULL)
				{
					if (!elemdisabled){
            if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						DynAddElement(cast,HEAD_SPECIFIED_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);					
					}
				}
		   break;
		  }
	  case(6): //------------------------------------------------------------------------------------
		  {//River
		   //string "Head", string name 
		   //{double x double y double head}x(numlines+1)
		   //&[int precision]
				ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer river in multilayer aquifer",BAD_DATA);
        thisname=strcpy(thisname,blank);
			  for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
				cast=CRiver::Parse(BBD,l,layer[cl],thisname);
				if (cast!=NULL)
				{
					if (!elemdisabled){
            if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						DynAddElement(cast,HEAD_SPECIFIED_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);	
					}
				}
				else {cout << "Bad River Entry"<<endl;}
			  
				/*CRiverSegment **tmp;
			  int NumSegs(0);
			  thisname=strcpy(thisname,blank);
			  for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			  tmp=CRiverSegment::Parse(BBD,l,layer[cl],thisname,NumSegs);
			  for (i=0;i<NumSegs; i++){
				  cast=tmp[i];
				  if (cast!=NULL){
					  if (!elemdisabled){
              if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					    DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					  }	
				  }
			  }
			  delete [] tmp;
			  //*/
		   break;
		  }
	  case(7): //------------------------------------------------------------------------------------
		  {//Inhomogeneity
		  //string "Inhom", string name 
			//double K 
		  //{double x double y}x(numlines+1)
		  //&[int precision]
			ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer inhomogeneity in multilayer aquifer",BAD_DATA);
			CInhom * tmp;
			thisname=strcpy(thisname,blank);
			for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			cast=tmp=CInhom::Parse(BBD,l,layer[cl],thisname);
			if (cast!=NULL)
			{
			  if (!elemdisabled){
					if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					DynAddPZone(tmp->GetCondZone(),PropZones[cl],numpzones[cl]);
				}	
			}

			/*CInhomSegment **tmp;
			CPropZone *zone=NULL;
			int NumSegs(0);
			thisname=strcpy(thisname,blank);
			for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			tmp=CInhomSegment::Parse(BBD,l,layer[cl],thisname,NumSegs,zone);
			for (i=0;i<NumSegs; i++){
				cast=tmp[i];
				if (cast!=NULL){
					if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					if (!elemdisabled){
					  DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					}	
				}
			}
			delete [] tmp;
			if ((!elemdisabled) && (zone!=NULL)){
				DynAddPZone(zone,PropZones[cl],numpzones[cl]);
			}*/
		  break;
			}
	  case(8): //------------------------------------------------------------------------------------
		  {//Leaky wall
			 //string "LeakyWall", string name 
		   //double thick
			 //double K
		   //{double x double y}x(numlines+1)
       //& [int precision] 
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer leaky wall in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CLeakyWall::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
					if (!elemdisabled){
						if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);	
					}
				}
		   break;
		  }
	  case(9): //------------------------------------------------------------------------------------
		  {//Drainy Crack
		   //string "Drain", string name 
		   //double thick
			 //double K
		   //{double x double y}x(numlines+1)
       //& [int precision] 
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer drain in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CDrain::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);					
				 }
			 }
		   break;
		  }
	  case(10): //------------------------------------------------------------------------------------
		  {/*Stage-Specified Linesink
		     string "Stage", string name  
		     double cond
			   double thickness
			   double width
		     {double x double y double head double depth}x(numlines+1) 
		     & [int precsion]
				*/
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer river in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CStage::Parse(BBD,l,layer[cl],thisname);
				if (cast!=NULL){
					if (!elemdisabled){
						if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						DynAddElement(cast,HEAD_SPECIFIED_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					}
				}
		   break;
		  }
    case(11): //------------------------------------------------------------------------------------
		  {/*AreaSink
		     string AreaSink, string name 
		     {double x double y}x(numlines+1) 
		     & [int order] [double fold]
		     {double x double y double leakage}x(numcontrolpts+1)
		     & [int precsion]
				*/
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer area sink in multilayer aquifer",BAD_DATA);
			 thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CAreaSink::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL){
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,GIVEN_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }
		   break;
		  }
	  case(12): //------------------------------------------------------------------------------------
		  {/*Circular Inhomogeneity
		     string "CircularInhomogeneity", string name 
		     double x double y double cond double radius
		     & [int precision] 
				*/
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer circular inhomogeneity in multilayer aquifer",BAD_DATA);
			 CCirInhom * tmp;
			 thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=tmp=CCirInhom::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					 DynAddPZone(tmp->GetCondZone(),PropZones[cl],numpzones[cl]);
				 }	
			 }
		   break;
		  }
	  case(13):  //------------------------------------------------------------------------------------
		  {/*vortex
         string "Vortex" , string name 
		     double x double y double strength
		     &*/
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer vortex in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CVortex::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 { 
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }
		   break;
		  }
	  case(14):  //------------------------------------------------------------------------------------
		  {//Point Dipole
       //string "Dipole" , string name 
		   //double x double y double strength double orientation
		   //&
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single point dipole in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CPtDipole::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
				 	 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }
		   break;
		  }
	  case(15): //------------------------------------------------------------------------------------
		  {//Qnorm
       //string "Qnorm", string name 
		   //{double x double y double Q}x(numlines+1)
       //&[int precision]
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer normal-discharge element in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CQnorm::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }
		   break;
		  }
		case(16): //------------------------------------------------------------------------------------
		  {//Circular Lake
		   //string "CirLake", string name 
		   //double x double y double elev double radius
		   //& [int precision] 
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer circular lake in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CCirLake::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
				 	 DynAddElement(cast,HEAD_SPECIFIED_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }
		   break;
		  }
	  case(17): //------------------------------------------------------------------------------------
		  {//Inhomogeneity in Base & thickness
		   //string "BaseInhom", string name 
			 //double Base double thickness
		   //{double x double y}x(numlines+1)
		   //&[int precision]
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer base inhomogeneity in multilayer aquifer",BAD_DATA);
			 CBaseInhom *tmp;
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=tmp=CBaseInhom::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
				 	 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					 DynAddPZone(tmp->GetBaseZone(),PropZones[cl],numpzones[cl]);
					 DynAddPZone(tmp->GetThickZone(),PropZones[cl],numpzones[cl]);
				 }
			 }
		   break;
		  }
	  case(18): //------------------------------------------------------------------------------------
		  {//Reservoir
		   //string "Reservoir", string name 
			 //double K 
		   //{double x double y}x(numlines+1)
		   //&[int precision]
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer reservoir in multilayer aquifer",BAD_DATA);
			 int thisprec;
			 double thisQ;
			 cmplex	stringp [MAXLINES];
       cmplex	leakzp [1];
       double leakp  [1];
			 if (noisy) {cout << "Reservoir element"<<endl;}
			 thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
		   eof=TokenizeLine(BBD,s,Len); l++; done=false;
		   do{ 
          if      (Len==0) {eof=TokenizeLine(BBD,s,Len); l++;}
			    else if (Len==1) {thisQ=s_to_d(s[0]); done=true;eof=TokenizeLine(BBD,s,Len); l++;}
			    else             {ImproperFormat(s,l); break;}
		   } while ((!done) && (!eof));
		   done=false;  nlines=0;
       do {
		     if (Len==0) {eof=TokenizeLine(BBD,s,Len); l++;}
			   else if (Len==2){
           stringp[nlines]=s_to_c(s[0],s[1]);
			     layer[cl]->UpdateExtents(stringp[nlines]);
			     nlines++; eof=TokenizeLine(BBD,s,Len); l++;}
         else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			     stringp[nlines]=stringp[0]; nlines++;
			     BASEMAP << "\" Reservoir \",  "<<-(nlines)<<endl;	 
			     for (i=0; i<nlines; i++){BASEMAP <<stringp[i].real()<< " , " <<stringp[i].imag()<<endl;}
			     if  (Len==2){thisprec= s_to_i(s[1]);}
			     else        {thisprec= defaultprec;}
			     cast = new CDrain(thisname,
														 layer[cl],
														 stringp,
														 nlines,
														 1e7,
														 1.0,
														 thisprec);
					 if (!elemdisabled){
						 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					 }
           leakzp[0]=0;
					 int    j;
           double Area(0.0);
           for (i=0;i<nlines;i++) {
             j = (i+1)%nlines;
             Area+= (stringp[i].real()*stringp[j].imag()-stringp[i].imag()*stringp[j].real())/2.0;
					 }
           leakp[0]=thisQ/Area;
			     cast = new CAreaSink(thisname,
																layer[cl],
																stringp,
																nlines-1,
																leakp,
																leakzp,
																thisprec,
																1);
					 if (!elemdisabled){
						 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						 DynAddElement(cast,GIVEN_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					 }
			     done=true;}
         else {ImproperFormat(s,l); break;}
		   } while ((!done) && (!eof));
		   break;
		  }
	  case(19): //------------------------------------------------------------------------------------
		  {//Horizontal Well
		   //string "HWell", string name 
		   //double Pumping Rate
		   //{double x double y}x(numlines+1)
       //&[int precision] 
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer well in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CHorizontalWell::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }		 
		   break;
		  }
	  case(20): //------------------------------------------------------------------------------------
		  {//Variable Stage-Specified Linesink
			 //string "VarStage", string name  
	     //double cond
	     //double thickness
	     //double width
	     //double roughness
	     //{double x double y double base double depth double runoff}x(numlines+1) 
	     //& [int precsion]
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer river in multilayer aquifer",BAD_DATA);
			 thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CVariableStage::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
				 	 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }		 
		   break;
		  }
	  case(21): //------------------------------------------------------------------------------------
		  {//Extraction Specified
       //string "Extraction", string name 
		   //{double x double y double Q}x(numlines+1)
       //&
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer extraction river in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CQSpec::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,GIVEN_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }
		   break;
		  }
	  case(22): //------------------------------------------------------------------------------------
		  {//Dry Well
			 //string "DryWell", string name
			 //double x double y double Qpump double r
			 //&  
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer well in multilayer aquifer",BAD_DATA);
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=CDryWell::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }
			 }
		   break;
		  }
	  case(23): //------------------------------------------------------------------------------------
		  {//Smooth Inhomogeneity
		   //string "Inhom", string name 
			 //{double x double y double Kxy)x(numpoints)
			 //&
		   //{double x double y}x(numlines+1)
		   //&[int precision]
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer smooth inhomogeneity  in multilayer aquifer",BAD_DATA);
			 CSmoothInhom * tmp;
			 thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=tmp=CSmoothInhom::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					 DynAddPZone(tmp->GetCondZone(),PropZones[cl],numpzones[cl]);
				 }	
			 }
		   break;
		  }
	  case(24): //------------------------------------------------------------------------------------
		  {//Elliptical Inhomogeneity
		   //string "EllipticalInhomogeneity", string name 
		   //double x double y double cond double a double b double angle(degrees)
		   //& [int precision] 
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer elliptical inhomogeneity  in multilayer aquifer",BAD_DATA);
			 CEllipseInhom * tmp;
			 thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=tmp=CEllipseInhom::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					 DynAddPZone(tmp->GetCondZone(),PropZones[cl],numpzones[cl]);				 
				 }	
			 }
		   break;
		  }
	  case(25): //------------------------------------------------------------------------------------
		  {//Elliptical Lake
		   //string "EllipticalLake", string name 
		   //double x double y double elev double a double b double angle(degrees)
		   //& [int precision] 
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer elliptical lake in multilayer aquifer",BAD_DATA);
			 CEllLake * tmp;
			 thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=tmp=CEllLake::Parse(BBD,l,layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
					 if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,HEAD_SPECIFIED_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
				 }	
			 }
		   break;
		  }
	  case(26): //------------------------------------------------------------------------------------
		  {//Surface Drain
		   //string "SurfaceDrain", string name 
		   //{double x double y double head}x(numlines+1)
		   //&[int precision]
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer surface drain in multilayer aquifer",BAD_DATA);
        thisname=strcpy(thisname,blank);
			  for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
				cast=CSurfaceDrain::Parse(BBD,l,layer[cl],thisname);
				if (cast!=NULL)
				{	
					if (!elemdisabled){
						if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
						DynAddElement(cast,HEAD_SPECIFIED_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);	
					}
				}
				else {cout << "Bad Surface Drain/Seepage Face Entry"<<endl;}
			  
		   break;
		  }
	  case(70): //------------------------------------------------------------------------------------
		  {//Inhomogeneity in Conductance
			 //string "ConductInhom", string name
			 //double conductance 
		   //{double x double y}x(numlines+1)
		   //&[int precision]
			 CConductInhom *tmp;
       thisname=strcpy(thisname,blank);
			 for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
			 cast=tmp=CConductInhom::Parse(BBD,l,aquicludes[cl],layer[cl],thisname);
			 if (cast!=NULL)
			 {
				 if (!elemdisabled){
				   if (eng.warm){eng.warm=cast->ReadItself(SOLUTION);}
					 DynAddElement(cast,REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					 DynAddPZone(tmp->GetConductZone(),LeakPropZones[cl],numlpzones[cl]);	
				 }
			 }
		   break;
		  }
		case(90): //------------------------------------------------------------------------------------
			{//Box of inhomgeneity cells
			 //string "BoxOfInhomCells", double xmin, double xmax, double ymin double ymax int nx double mincond double maxcond
				 //ostream ost;
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot create single layer box of inhomogeneity cells in multilayer aquifer",BAD_DATA);
				if (noisy) {cout << "Box Of Inhom Cells"<<endl;}
				if (Len==8)
				{
					CInhom **Inhoms;
					Inhoms=CInhom::CreateBoxOfInhoms(layer[cl],
						                               s_to_d(s[1]),s_to_d(s[2]),s_to_d(s[3]),
																					 s_to_d(s[4]),s_to_i(s[5]),s_to_d(s[6]),
																					 s_to_d(s[7]));
          for (i=0; i<(s_to_i(s[5])*s_to_i(s[5])); i++){
						//cannot read solution as of now.
						eng.warm =false;
					  DynAddElement(Inhoms[i],REGULAR_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					  DynAddPZone(Inhoms[i]->GetCondZone(),PropZones[cl],numpzones[cl]);	
					}
					if (eng.warm){ExitGracefully("The Bluebird engine cannot currently read the solution for a box of inhomogeneity cells",BAD_DATA);}
				 }
				 else {ImproperFormat(s,l); break;}
				 break;
			}
	  case(101): //------------------------------------------------------------------------------------
			{//Solve Only 
			 if (noisy) {cout <<"SolveOnly"<<endl;}      
			 eng.solve=true; eng.grid=false; eng.track=false; eng.transport=false;
			 break;
			}            
	  case(102): //------------------------------------------------------------------------------------                      
			{//Grid Only
			 if (noisy) {cout <<"GridOnly"<<endl;}  
			 if (solutionexists==false){ExitGracefully("Cannot grid without solution",BAD_DATA);}
			 eng.solve=false; eng.grid=true;  eng.track=false; eng.warm=true; eng.transport=false;
			 break;
			}
    case(103): //------------------------------------------------------------------------------------
			{//Track only
			 if (noisy) {cout <<"TrackOnly"<<endl;}
			 if (solutionexists==false){ExitGracefully("Cannot track particles without solution",BAD_DATA);}
       eng.solve=false; eng.grid=false; eng.track=true;  eng.warm=true; eng.transport=false;  eng.transient=false;
			 break;
			} 
	  case(104): //------------------------------------------------------------------------------------
			{//Solve And Grid
			 if (noisy) {cout <<"Solve And Grid"<<endl;}
			 eng.solve=true; eng.grid=true;  eng.track=false; eng.transport=false; 
			 break;
			}						
	  case(105): //------------------------------------------------------------------------------------
			{//Solve And Track
			 if (noisy) {cout <<"Solve And Track"<<endl;}
			 eng.solve=true; eng.grid=false; eng.track=true; eng.transport=false; eng.transient=false;
			 break;
			}						
	  case(106): //------------------------------------------------------------------------------------
			{//Solve Grid Track
			 if (noisy) {cout <<"SolveGridTrack"<<endl;}
			 eng.solve=true; eng.grid=true;  eng.track=true; eng.transport=false;  eng.transient=false; 
			 break;
			}						
	  case(107): //------------------------------------------------------------------------------------
		  {//warm start 
				if (noisy) {cout <<"Warm Start"<<endl;}
			  if (solutionexists){
					if (Len==2) {eng.warmtime=s_to_d(s[1]);}
					eng.warm=true; CSingleLayer::fresh=false;
				}
				else {
					ExitGracefully("Parse::WarmStart used without solution file",BAD_DATA);
				}
				break;
			}                                     
	  case(108): //------------------------------------------------------------------------------------
			{//Superblocks
			 //string "Superblocks", [int nesting levels], [int precision], [str "explicit" or "implicit"]
				if (noisy) {cout << "SuperBlocks On" <<endl;}
				if      (Len==1){SuperNest=1;           blocksolve=true;blockprec=defaultprec; blockimplicit=true;}
			  else if (Len==2){SuperNest=s_to_i(s[1]);blocksolve=true;blockprec=defaultprec; blockimplicit=true;}
        else if (Len==3){SuperNest=s_to_i(s[1]);blocksolve=true;blockprec=s_to_i(s[2]);blockimplicit=true;}
        else if (Len==4){
					if      (!strcmp(s[3],"explicit")) {
						SuperNest=s_to_i(s[1]);blocksolve=true;blockprec=s_to_i(s[2]);blockimplicit=false;
					}
					else if (!strcmp(s[3],"implicit")) {
						SuperNest=s_to_i(s[1]);blocksolve=true;blockprec=s_to_i(s[2]);blockimplicit=true;
					}
					else					{ImproperFormat(s,l);}	
				}
        else            {ImproperFormat(s,l);}		 
			 break;
			}    
		case(109): //------------------------------------------------------------------------------------
			{//GridWhileSolve
			 //string "GridWhileSolve", double topbound, double bottombound, double leftbound, double rightbound, int resolution
				if (noisy) {cout << "GridWhileSolve On" <<endl;}
			  if (Len==6) {
					CSingleLayer::SetGridWhileSolve(s_to_d(s[1]),s_to_d(s[4]),s_to_d(s[3]),s_to_d(s[2]),s_to_i(s[5]));
				}
		    else {ImproperFormat(s,l);}
		   break;
			}
		case (110)://------------------------------------------------------------------------------------
			{//Transport Only
			 if (noisy) {cout<< "Transport Only"<<endl;}     
			 if (solutionexists==false){ExitGracefully("Cannot transport contaminant without solution",BAD_DATA);}
			 eng.solve=false; eng.grid=false;  eng.track=false; eng.warm=true; eng.transport=true; eng.transient=false;
			 break;
			}
		case (111)://------------------------------------------------------------------------------------
			{//Solve Grid Transport
			 if (noisy) {cout<< "Solve Grid Transport"<<endl;}
       eng.solve=true; eng.grid=true;  eng.track=false; eng.transport=true; eng.transient=false;
			 break;
			}
		case (112)://------------------------------------------------------------------------------------
			{//Solve and Transport
			 if (noisy) {cout<< "Solve and Transport"<<endl;}
			 eng.solve=true; eng.grid=false;  eng.track=false; eng.transport=true; eng.transient=false;
			 break;
			}
		case (113)://------------------------------------------------------------------------------------
			{//Taylor Series Circle
				CTaylorCir *tmp;
				tmp=CTaylorCir::Parse(BBD,l,layer[cl],defaultprec);
				if (tmp!=NULL){layer[cl]->AddToLayer(tmp);}
			 break;
			}
		case (114)://------------------------------------------------------------------------------------
			{//Transient Flow Solution
			 //string Transient, double timestep 
				if (Len==2){
					eng.timestep=s_to_d(s[1]);
					eng.transient=true;
				}
			 break;
			}
		case (115)://------------------------------------------------------------------------------------
			{//Transient Flow Output Times
				eng.nTimes=0;
				if (noisy){cout << "Transient Output Times"<<endl;}
				done=false;                                              eof=TokenizeLine(BBD,s,Len); l++;
			  do {
				  if  ((Len==1) && (strcmp(s[0],"&"))){
				 		eng.outtimes[eng.nTimes]=s_to_d(s[0]); eng.nTimes++; eof=TokenizeLine(BBD,s,Len); l++;
					}
					else if ((Len==1) && (!strcmp(s[0],"&"))) { 
						done=true;
					}
					else if (Len==0) {eof=TokenizeLine(BBD,s,Len); l++;}
					else {ImproperFormat(s,l); break;}
				} while ((!done) && (!eof) && (eng.nTimes<MAXTIMES));				
			 break;
			}
		case (116)://------------------------------------------------------------------------------------
			{//Relaxation
			 //string "Relax", string "Coeff", double relaxationcoeff 
				if (noisy) {cout << "River Relaxation"<<endl;}
				if (Len==3){CRiver::SetRelaxation(s_to_d(s[2]));}
			  else       {ImproperFormat(s,l);}
			 break;
			}
		case (117)://------------------------------------------------------------------------------------
			{//Fully Explicit Solve
			 //string FullyExplicit
				if (noisy) {cout << "Fully Explicit Solve"<<endl;}
				if (Len==1){eng.explicitsolve=true;}
			  else       {ImproperFormat(s,l);}
			 break;
			}
		case (118)://------------------------------------------------------------------------------------
			{//Solve Far Field Every Time
			 //string "SolveFarFieldEveryTime"
				if (noisy) {cout << "Solve FarField Every Time On"<<endl;}
				if (Len==1){CSingleLayer::SetFFEveryTime();}
				else       {ImproperFormat(s,l);}
			 break;
			}
		case (119)://------------------------------------------------------------------------------------
			{//Write interval
			 //string "WriteInterval", int number of iterations 
				if (noisy) {cout << "WriteInterval"<<endl;}
				if (Len==2){CSingleLayer::SetWriteInterval(s_to_i(s[1]));}
				else       {ImproperFormat(s,l);}
			 break;
			}
		case (120)://------------------------------------------------------------------------------------
			{//Socket On
			 //string "SocketOn"
				if (noisy) {cout << "Sockets On"<<endl;}
				if (Len==1){eng.socket=true;}
				else       {ImproperFormat(s,l);}
			 break;
			}
		case (121)://------------------------------------------------------------------------------------
			{//Full Analysis of iterative convergence
			 //string "FullAnalyze"
				if (noisy) {cout << "Full Analysis"<<endl;}
				if (Len==1){CSingleLayer::SetFullAnalysis();}
				else       {ImproperFormat(s,l);}
			 break;
			}

		case (122)://------------------------------------------------------------------------------------
			{//Non-nested (Classic) Superblocks
			 //string "ClassicSuperblocks"
				if (noisy) {cout << "Classic Superblocks On"<<endl;}
				if (Len==1){blocknested=false;}
				else       {ImproperFormat(s,l);}
				break;
			}
		case (123)://------------------------------------------------------------------------------------
			{ 
			break;
			}
		case (124)://------------------------------------------------------------------------------------
			{//Area Sink Relaxation  on
			 //string "AreaSinkRelax", int startiter, int enditer
				if (noisy) {cout << "Area Sink Relaxation  on"<<endl;}
				if (Len==3){CAreaSink::SetRelaxation(s_to_i(s[1]),s_to_i(s[2]));}
				else       {ImproperFormat(s,l);}
				break;
			}
		case (125)://------------------------------------------------------------------------------------
			{//Write Solve Time
			 //string "AreaSinkRelax", int startiter, int enditer
				if (noisy) {cout << "Write Solve Time  on"<<endl;}
				eng.writesolvetime =true;
				break;
			}
	  case(201): //------------------------------------------------------------------------------------
		  {//Single Layer Base And Thickness
				if (noisy) {cout <<"layer base and thickness"<<endl;}
			  ExitGracefullyIf(!singlelayer,"Parse::Cannot specify single layer base and thickness in multilayer aquifer",BAD_DATA);
				if (Len==3){aqbase[cl]=s_to_d(s[1]); aqthickness[cl]=s_to_d(s[2]);}
				else       {ImproperFormat(s,l);}
				break;
			}
    case(202): //------------------------------------------------------------------------------------
		  {//Single Layer Conductivity
			 if (noisy) {cout <<"layer conductivity"<<endl;}
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot specify single layer conductivity in multilayer aquifer",BAD_DATA);
		   if (Len==2){aqconductivity[cl]=s_to_d(s[1]);}
		   else       {ImproperFormat(s,l);}
       break; 
		  }
	  case(203): //------------------------------------------------------------------------------------
		  {//Single Layer Porosity 
			 if (noisy) {cout <<"layer porosity"<<endl;}
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot specify single layer porosity in multilayer aquifer",BAD_DATA);
		   if (Len==2){aqporosity[cl]=s_to_d(s[1]);}
		   else       {ImproperFormat(s,l);}
       break; 
		  }
	  case(204): //------------------------------------------------------------------------------------
			{//MultiLayer 
			//string "MultiLayer", int numlayers
			singlelayer=false;
			mlnumlevels[cl]=s_to_i(s[1]);
			mlayer[cl]->SetNumLayers(s_to_i(s[1]));
			break;
			}
	  case(205): //------------------------------------------------------------------------------------
		  {//New Layer Start
			 if (noisy) {cout <<"Next Layer Start-----------------------------------------"<<endl;}
			 if (Len==1){
				  cl++;
					numlevels++;
					ExitGracefullyIf(cl>=MAXAQLEVELS,"Parse::Exceeded maximum levels in aquifer",BAD_DATA);

					if (singlelayer)
					{
						layer[cl]=new CSingleLayer();
						if ((layer[cl]==NULL)){
							ExitGracefully("Parse::New SingleLayer Creation- out of memory",OUT_OF_MEMORY);}
					}
					else//multilayer
					{
						mlayer[cl]=new CMultiLayer();
						if ((mlayer[cl]==NULL)){
							ExitGracefully("Parse::New MultiLayer Creation- out of memory",OUT_OF_MEMORY);}
					}
					//TMP DEBUG- new farfield should be created here so that it may read its solution

			 }
		   else       {ImproperFormat(s,l);}
       break; 
		  }	  
		case(207): //------------------------------------------------------------------------------------
		  {//Aquiclude (for below layer)
			 //must come after aquifer base/thickness, layer information
			 //string "Aquiclude", double thick, double k {k ignored in BB}
			 if (noisy) {cout <<"Aquiclude"<<endl;}
			 if (Len==3){
					double thisthick=s_to_d(s[1]);
					double thiscond =s_to_d(s[2]); //ignored for now
					//TMP DEBUG- need function to check all bases and thicknesses for multilayer

					//TMP DEBUG- transition to multilayer
					//if (cl==0){aquicludes[cl]=new CAquiclude(NULL       ,layer[cl],aqbase[cl],cl);}
					//else      {aquicludes[cl]=new CAquiclude(layer[cl-1],layer[cl],aqbase[cl],cl);}
					
					if (aquicludes[cl]==NULL){
						ExitGracefully("Parse::New Aquiclude Creation- out of memory",OUT_OF_MEMORY);}
			 }
		   else       {ImproperFormat(s,l);}
       break; 
		  }
		case(208)://-----------------------------------------------------------------------------------
			{//Porosity Zone
		   //string "PorosityZone", string name (ignored) 
			 //double n 
		   //{double x double y}x(numlines+1)
		   //&
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot specify single layer porosity zone in multilayer aquifer",BAD_DATA);
			 cast2=CPolyPropZone::Parse(BBD,l,poro);
			 if (cast2!=NULL){
				 DynAddPZone(cast2,PropZones[cl],numpzones[cl]);
			 }			
		   break;
		  }
	  case(209): //------------------------------------------------------------------------------------
			{//SingleLayer 
			//string "SingleLayer"
			singlelayer=true;
			break;
			}
	  case(210): //------------------------------------------------------------------------------------
			{//Conductivities 
			//string "Conductivities"
			//{double K} x (numlevels)
			//&
			//{double k} x (numlevels-1)
			//&
			ExitGracefullyIf(singlelayer,"Parse::Cannot specify multiple conductivities for a single layer aquifer",BAD_DATA);
			
			break;
			}
	  case(211): //------------------------------------------------------------------------------------
			{//Thicknesses 
			//string "Thicknesses"
			//{double T} x (numlevels)
			//&
			//{double t} x (numlevels-1)
			//&
			ExitGracefullyIf(singlelayer,"Parse::Cannot specify multiple thicknesses for a single layer aquifer",BAD_DATA);
			
			break;
			}
	  case(212): //------------------------------------------------------------------------------------
			{//Porosities 
			//string "Porosities"
			//{double n} x (numlevels)
			//&
			//{double naq} x (numlevels-1)
			//&
			ExitGracefullyIf(singlelayer,"Parse::Cannot specify multiple thicknesses for a single layer aquifer",BAD_DATA);
			
			break;
			}
	  case(213): //------------------------------------------------------------------------------------
			{//MLProperties  
			//string "MultiLayerProperties"
			//{double K double T double n} x (numlevels)
			//&
			//{double k double t double naq} x (numlevels-1)
			//&
			ExitGracefullyIf(singlelayer,"Parse::Cannot specify multilayer properties for a single layer aquifer",BAD_DATA);
			
			break;
			}
	  case(301): //------------------------------------------------------------------------------------
		  {//Uniform Flow
			 //string "UniformFlow", Qx, Qy
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot specify uniform flow discharge in multilayer aquifer. Must specify a gradient",BAD_DATA);
			 if (noisy) {cout <<"uniform flow"<<endl;}
		   if (Len==3) {FFQxy[cl]=s_to_c(s[1],s[2]);}
		   else        {ImproperFormat(s,l);break;}
		   break;
		  }
	  case(302): //------------------------------------------------------------------------------------
		  {//Gradient
			 //string "Gradient", Ix, Iy
			 if (noisy) {cout <<"uniform flow (gradient)"<<endl;}
		   if (Len==3) {FFgradient=s_to_c(s[1],s[2]);}
		   else        {ImproperFormat(s,l);break;}
		   break;
		  }
	  case(303): //------------------------------------------------------------------------------------
		  {//Reference point
			//string "ReferencePoint" x, y, head
			if (noisy) {cout <<"Farfield: Reference Point"<<endl;}
		  if (Len==4) {
				FFzref   =s_to_c(s[1],s[2]);
				FFrefhead=s_to_d(s[3]);
				if (singlelayer){
					ffp[cl]=new CFarField(layer[cl],FFQxy[cl],FFzref,FFrefhead);
					if (ffp[cl]!=NULL){ 
					  if (eng.warm){eng.warm=ffp[cl]->ReadItself(SOLUTION);}
					}
					else{ExitGracefully("Parse::Unable to create Single Layer Farfield element",BAD_DATA);}
				}
				else if (!singlelayer){
					MLffp[cl]=new CmlFarField(mlayer[cl],FFgradient,FFzref,FFrefhead);
					if (MLffp[cl]!=NULL){ 
					  if (eng.warm){eng.warm=MLffp[cl]->ReadItself(SOLUTION);}
					}
					else{ExitGracefully("Parse::Unable to create Multi-Layer Farfield element",BAD_DATA);}					
				}
		  }
		  else{ImproperFormat(s,l);break;}
		  break;}		   
	  case(304): //------------------------------------------------------------------------------------
		  {//Net Extraction
			 if (noisy) {cout <<"farfield-net extraction"<<endl;}
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot specify net extraction condition in multilayer aquifer",BAD_DATA);
		   if (Len==2) {
         ffp[cl]=new CFarField(layer[cl],FFQxy[cl],s_to_d(s[1]));
				 if (ffp[cl]!=NULL){ 
					 if (eng.warm){eng.warm=ffp[cl]->ReadItself(SOLUTION);}
				 }
				 else{ExitGracefully("Parse::Unable to create Farfield element",BAD_DATA);}
		   }
		   else{ImproperFormat(s,l);break;}
		   break;
		  }	
	  case(401): //------------------------------------------------------------------------------------
		  {//Precision
			 if (noisy) {cout <<"Precision"<<endl;}
			 if (Len==2){defaultprec=s_to_i(s[1]);
									 CAnalyticElem::SetDefaultPrecision(defaultprec);}
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
		case(402): //------------------------------------------------------------------------------------
		  {//Silent
			 if (noisy) {cout <<"Silent"<<endl;}
		   if (Len==1){noisy=false;}
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
		case(403): //------------------------------------------------------------------------------------
		  {//Noisy
			 if (noisy) {cout <<"Noisy"<<endl;}
		   if (Len==1){noisy=true;}
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
	  case(406): //------------------------------------------------------------------------------------
		  {//MaximumSolveIterations
			 if (noisy) {cout <<"MaxIterations"<<endl;}
			 if (Len==2){
				 CSingleLayer::SetSolveData(0,s_to_i(s[1]), 0);
				 CMultiLayer ::SetSolveData(0,s_to_i(s[1]), 0);
			 }
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
		case(407): //------------------------------------------------------------------------------------
		  {//MinimumSolveIterations
			 if (noisy) {cout <<"MinIterations"<<endl;}
		   if (Len==2){
				 CSingleLayer::SetSolveData(s_to_i(s[1]), 0, 0);
				 CMultiLayer ::SetSolveData(s_to_i(s[1]), 0, 0);
			 }
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
		case(408): //------------------------------------------------------------------------------------
		  {//Potential Tolerance
			 if (noisy) {cout <<"Tolerance"<<endl;}
		   if (Len==2){
				 CSingleLayer::SetSolveData(0, 0,s_to_d(s[1]));
				 CMultiLayer ::SetSolveData(0, 0,s_to_d(s[1]));
			 }
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
		case(409): //------------------------------------------------------------------------------------
		  {//Compare Heads On
			  if (noisy) {cout <<"Compare Heads On"<<endl;}
			  eng.obsfileexists=true; //old format still accounted for
		   break;
		  }
		case(413): //------------------------------------------------------------------------------------
		  {//WriteSolutionOn
			 if (noisy) {cout <<"Write Solution On"<<endl;}
		   if (Len==1){eng.writesol=true;}
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
		case(414)://------------------------------------------------------------------------------------
		  {//WriteSolutionOff
			 if (noisy) {cout <<"Write Solution Off"<<endl;}
		   if (Len==1){eng.writesol=false;}
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
		case(415)://------------------------------------------------------------------------------------
		  {//WriteOutputOff
			 if (noisy) {cout <<"Write Output Off"<<endl;}
		   if (Len==1){eng.writeout=false;}
		   else       {ImproperFormat(s,l);break;}
		   break;
		  }
	  case(501): //------------------------------------------------------------------------------------
		  {//Gridding Window 
		   //window center x,y, cell size, resolution in x resolution in y
			 if (noisy) {cout <<"Gridding Window: type 1"<<endl;}
		   if (Len==6) {
			   double GridCenX=s_to_d(s[1]);
			   double GridCenY=s_to_d(s[2]);
				 eng.Grid.W.e = GridCenX + (s_to_d(s[3])*s_to_d(s[4])/ 2.0);
         eng.Grid.W.w = GridCenX - (s_to_d(s[3])*s_to_d(s[4])/ 2.0);
			   eng.Grid.W.n = GridCenY + (s_to_d(s[3])*s_to_d(s[5])/ 2.0);
         eng.Grid.W.s = GridCenY - (s_to_d(s[3])*s_to_d(s[5])/ 2.0);
			   eng.Grid.res=s_to_i(s[4]);
		   }
		   else {ImproperFormat(s,l);}
		   break;
		  }
	  case(502): //------------------------------------------------------------------------------------
		  {//Gridding Window 2  
		   //window2 double north double south double east double west int resolution 
			 if (noisy) {cout <<"Gridding Window: type 2"<<endl;}
		   if (Len==6) {
				 eng.Grid.W.n=s_to_d(s[1]);
				 eng.Grid.W.s=s_to_d(s[2]);
				 eng.Grid.W.e=s_to_d(s[3]); 
				 eng.Grid.W.w=s_to_d(s[4]);
			   eng.Grid.res=s_to_i(s[5]);
		   }
		   else {ImproperFormat(s,l);}
		   break;
		  }
		case(503): //------------------------------------------------------------------------------------
		  {//GridQxQy
		   //string "GridQxQy"
			 if (noisy)  {cout <<"Grid Qx and Qy on"<<endl;}
		   eng.Grid.QxQy=true;
		   break;
		  }
		case(504): //------------------------------------------------------------------------------------
		  {//GridBaseThickness
		   //string "GridBaseThickness"
			 if (noisy)  {cout <<"Grid Base and Thickness on"<<endl;}
		   eng.Grid.base =true;
			 eng.Grid.top  =true;
			 eng.Grid.thick=true;
		   break;
		  }
		case(505): //------------------------------------------------------------------------------------
		  {//GridCond
		   //string "GridCond"
			 if (noisy)  {cout <<"Grid Conductivity and Conductance"<<endl;}
		   if (Len==1) {eng.Grid.cond   =true;
										eng.Grid.conduct=true;}
		   else {ImproperFormat(s,l);}
		   break;
		  }
		case(506): //------------------------------------------------------------------------------------
		  {//GridElement
		   //string "GridElement"
			 //Grids Previously parsed element
			 if (noisy)  {cout <<"Grid Previous Element"<<endl;}
		   if (Len==1) {
				 if (cast!=NULL){
					 eng.Grid.elem        =true;
					 eng.Grid.pGrdElement =cast;
				 }
			 }
		   else {ImproperFormat(s,l);}
		   break;
		  }
		case(507): //------------------------------------------------------------------------------------
		  {//AnalysisFace
		   //string "AnalysisFace" int ID double x1 double y1 double x2 double y2 int numpoints
			 //creates output along line
			  if (noisy)  {cout <<"Analysis Face (Transect)"<<endl;}
			  if (Len==7) {
					ExitGracefullyIf(!singlelayer,"Parse::Cannot currently use transect in multilayer aquifer",BAD_DATA);
				  CTransect *tmp=new CTransect(s_to_i(s[1]),
																			 layer[cl],//should sent MLsublayer
																			 s_to_c(s[2],s[3]),
																			 s_to_c(s[4],s[5]),
																			 s_to_i(s[6]));
				}
			  else {ImproperFormat(s,l);}
		   break;
		  }
		case(508): //------------------------------------------------------------------------------------
		  {//Zone Budget
		   //string "ZoneBudget"
		   //{double x double y}x(numlines+1)
		   //&
			 ExitGracefullyIf(!singlelayer,"Parse::Cannot currently use zone budget in multilayer aquifer",BAD_DATA);
			 CZoneBudget * tmp;
			 tmp=CZoneBudget::Parse(BBD,l,layer[cl]);
		   break;
		  }
		case(509): //------------------------------------------------------------------------------------
		  {//PlotGxGy
		   //string "PlotGxGy"
			 if (noisy)  {cout <<"Plot Gx and Gy on"<<endl;}
		   eng.Grid.GxGy=true;
		   break;
		  }
		case(510): //------------------------------------------------------------------------------------
		  {//PlotVxVy
		   //string "PlotVxVy"
			 if (noisy)  {cout <<"Plot Velocity on"<<endl;}
		   eng.Grid.vxvy=true;
		   break;
		  }
		case(511): //------------------------------------------------------------------------------------
		  {//PlotEffVelocity
		   //string "PlotEffVelocity"
			 if (noisy)  {cout <<"Plot Effective Velocity on"<<endl;}
		   eng.Grid.effvel=true;
		   break;
		  }
		case(512): //------------------------------------------------------------------------------------
		  {//PlotCurl
		   //string "PlotCurl"
			 if (noisy)  {cout <<"Plot Curl on"<<endl;}
		   eng.Grid.curl=true;
		   break;
		  }
		case(513): //------------------------------------------------------------------------------------
		  {//PlotDxxDxDyyDy
		   //string "PlotDxxDxDyyDy"
			 if (noisy)  {cout <<"Plot Dxx/Dx-Dyy/Dy on"<<endl;}
		   eng.Grid.DxxDx=true;
		   break;
		  }
		case(514): //------------------------------------------------------------------------------------
		  {//PlotGxGyNumerical
		   //string "PlotGxGyNumerical"
			 if (noisy)  {cout <<"Plot Gx and Gy (numerical) on"<<endl;}
		   eng.Grid.GxGynum=true;
		   break;
		  }
		case(515): //------------------------------------------------------------------------------------
		  {//PlotVelocityMagnitudeDerivatives
		   //string "PlotVelocityMagnitudeDerivatives"
			 if (noisy)  {cout <<"Plot Velocity Magnitude (numerical) on"<<endl;}
		   eng.Grid.VmagDer=true;
		   break;
		  }
		case(516): //------------------------------------------------------------------------------------
		  {//PlotVelocityDerivatives
		   //string "PlotVelocityDerivatives"
			 if (noisy)  {cout <<"Plot Velocity Derivatives (numerical) on"<<endl;}
		   eng.Grid.vder=true;
		   break;
		  }
		case(517): //------------------------------------------------------------------------------------
		  {//PlotDxyDxDxyDy
		   //string "PlotDxyDxDxyDy"
			 if (noisy)  {cout <<"Plot Dxy/Dx-Dxy/Dy on"<<endl;}
		   eng.Grid.DxyDx=true;
		   break;
		  }
		case(518): //------------------------------------------------------------------------------------
		  {//PlotEffVelNum
		   //string "PlotEffVelNum"
			 if (noisy)  {cout <<"Plot Numerical Effective Velocity on"<<endl;}
		   eng.Grid.effvelnum=true;
		   break;
		  }
		case(519): //------------------------------------------------------------------------------------
		  {//PlotSatThick
		   //string "PlotSatThick"
			 if (noisy)  {cout <<"Plot Saturated Thickness on"<<endl;}
		   eng.Grid.satthick=true;
		   break;
		  }
		case (520)://------------------------------------------------------------------------------------
			{ //"ObservationFileExists"
			eng.obsfileexists=true;
			break;
			}
		case(603): //------------------------------------------------------------------------------------
		  {//Time
		   //"time", starttime(ignored), endtime 
			 if (noisy) {cout <<"Tracking Duration"<<endl;}
		   if (Len==3) {CPathline::SetDuration(s_to_d(s[1]));}
		   else        {ImproperFormat(s,l);                 }
		   break;
		  }
		case(604)://------------------------------------------------------------------------------------
		  {//Particle
		   //string "Particle", double x, double y
			 if (noisy) {cout <<"Particle"<<endl;}
			 if (Len==3) {
				 if (numpart<MAX_PARTICLES){
           particles[numpart]= new CPathline(numpart,aq,c2Dto3D(s_to_c(s[1],s[2])),TRACK_FORWARD);
				   numpart++;
				 }
				 else{
					 ExitGracefully("Parse: too many particles",BAD_DATA);
				 }
			 }
			 else {ImproperFormat(s,l);}
		   break;
		  }
		case(605)://------------------------------------------------------------------------------------
		  {//Capture Zone
		   //string "CaptureZone", double x, double y, int numparticles, [ double radius ]
			 if (noisy) {cout <<"Capture Zone"<<endl;}
			 double thisradius=5; //TEMP DEBUG -should find nearest well
			 cmplex loc;
			 double angle;
			 if (Len==4) {
				 for (i=0; i<s_to_i(s[3]); i++){
					 angle=2.0*PI*(double)(i)/(double)(s_to_i(s[3]));
					 loc=thisradius*cmplex(cos(angle), sin(angle))+s_to_c(s[1],s[2]);
					 if (numpart<MAX_PARTICLES){
				     particles[numpart]= new CPathline(numpart, 
																							 aq, 
																							 c2Dto3D(loc), 
																							 TRACK_BACKWARD);
						 numpart++;
					 }
				   else{
					   ExitGracefully("Parse: too many particles",BAD_DATA);
					 }
				 }
			 }
			 else if (Len==5) { //more of a circle of particles
				 thisradius=s_to_d(s[4]);
				 for (i=0; i<s_to_i(s[3]); i++){
					 angle=2.0*PI*(double)(i)/(double)(s_to_i(s[3]));
					 loc=thisradius*cmplex(cos(angle), sin(angle))+s_to_c(s[1],s[2]);
					 if (numpart<MAX_PARTICLES){
				     particles[numpart]= new CPathline(numpart, 
																							 aq, 
																							 c2Dto3D(loc), 
																							 TRACK_BACKWARD);
						 numpart++;
					 }
					 else{
					   ExitGracefully("Parse: too many particles",BAD_DATA);
					 }
				 }
			 }
			 else {ImproperFormat(s,l);}
		   break;
		  }
		case(607)://------------------------------------------------------------------------------------
		  {//Box Of Particles 
		   //string "BoxOfParticles", double xmin, double xmax, double ymin double ymax int nx int ny 
			 if (noisy) {cout <<"Box Of Particles"<<endl;}
			 if (Len==7) {
				 double xincrem=(s_to_d(s[2])-s_to_d(s[1]))/(s_to_i(s[5])-1); 
         double yincrem=(s_to_d(s[4])-s_to_d(s[3]))/(s_to_i(s[6])-1); 

				 for   (double xgrid=s_to_d(s[1]); xgrid<=s_to_d(s[2]); xgrid+=xincrem){
           for (double ygrid=s_to_d(s[3]); ygrid<=s_to_d(s[4]); ygrid+=yincrem){
						 if (numpart<MAX_PARTICLES){
						 	 particles[numpart]= new CPathline(numpart, 
																								 aq, 
																								 c2Dto3D(cmplex (xgrid,ygrid)), 
																								 TRACK_FORWARD);
							 numpart++;
						 }
						 else{
							 ExitGracefully("Parse: too many particles",BAD_DATA);
						 }
					 }
				 }
			 }
			 else {ImproperFormat(s,l);}
		   break;
		  }
		case(608)://-----------------------------------------------------------------------------------
			{//Tracking Duration
			 //string "TrackDuration", double time
       if (noisy)	 {cout <<"Tracking Duration"<<endl;}
			 if (Len==2) {CPathline::SetDuration(s_to_d(s[1]));}
			 else        {ImproperFormat(s,l);}
			 break;
			}
		case(609)://-----------------------------------------------------------------------------------
			{//Tracking Resolution
			 //string "TrackResolution", integer resolution
       if (noisy) {cout <<"Tracking Resolution"<<endl;}
			 if (Len==2){CPathline::SetResolution(s_to_i(s[1]));}
			 else       {ImproperFormat(s,l);}
			 break;
			}
		case(610)://-----------------------------------------------------------------------------------
			{//Tracking Precision
			 //string "TrackPrecision", integer resolution
       if (noisy) {cout <<"Tracking Precision"<<endl;}
			 if (Len==2){CParticle::SetPrecision(s_to_i(s[1]));}
			 else				{cout<<"line"<< l <<"is wrong length"<<endl;}
			 break;
			}
		case(611)://-----------------------------------------------------------------------------------
			{//Interpolation Grid file Attached
			 //string "InterpolationGridFile"
				if (noisy) {cout <<"Interpolation Grid File Read"<< endl;}
				if (Len==1){
					ifstream INTGRID;
					INTGRID.open("InterpolationGrid.bbg");
					int cl(0);
					CRectGrid *pGrid;
					pGrid=NULL;
					pGrid=CRectGrid::ReadBBGFile(INTGRID,cl);
					INTGRID.close();
					if (pGrid!=NULL){
						layer[cl]->SetInterpolationGrid(pGrid);	
						eng.interpolate=true;
					}
				}
				else {ImproperFormat(s,l);}
				break;   
		  }
		case(612)://-----------------------------------------------------------------------------------
			{//Interpolation Grid 
			 //string "InterpolationGrid", double north double south double east double west int resolutionX {int resolutionY}
				if (noisy) {cout <<"Interpolation Grid "<< endl;}
				if ((Len==6) || (Len==7)) {
					CRectGrid *pGrid;
					if (Len==7){pGrid=new CRectGrid(s_to_d(s[1]),s_to_d(s[2]),s_to_d(s[3]),s_to_d(s[4]),s_to_i(s[5]),s_to_i(s[6]));}
					else       {pGrid=new CRectGrid(s_to_d(s[1]),s_to_d(s[2]),s_to_d(s[3]),s_to_d(s[4]),s_to_i(s[5]),s_to_i(s[5]));}
					if (pGrid!=NULL){
						layer[cl]->SetInterpolationGrid(pGrid);	
						eng.interpolate=true;
					}
				}
				else {ImproperFormat(s,l);}
		    break; 
			}
		case(613)://------------------------------------------------------------------------------------
		  {	//Particle3D
				//string "Particle3D", double x, double y, double z
				if (noisy) {cout <<"3D Particle"<<endl;}
				if (Len==4) {
				  if (numpart<MAX_PARTICLES){
			      particles[numpart]= new CPathline(numpart, 
																						  aq, 
																							CVector(s_to_d(s[1]),s_to_d(s[2]),s_to_d(s[3])), 
																							TRACK_FORWARD);
					  numpart++;
					}
					else{
						ExitGracefully("Parse: too many particles",BAD_DATA);
					}
				}
				else {ImproperFormat(s,l);}
				break;
		  }
	  case(615)://------------------------------------------------------------------------------------
		  {	//Read Pollock's method grid, use pollocks method
				//string "PollocksMethodOn"
				if (noisy) {cout <<"Reading Tracking Grid (for Pollocks Method)..."<< endl;}
				if (Len==1){
					ifstream TRACKGRID;
					int cl(0);
					CRectGrid *pGrid=NULL;
					CFluxGrid *pFgrid=NULL;
					TRACKGRID.open("concgrid.bbg"); //TMP DEBUG
					pGrid=CRectGrid::ReadBBGFile(TRACKGRID,cl);
					TRACKGRID.close();
					if (pGrid!=NULL){
						pFgrid=new CFluxGrid(pGrid);
						CPathline::SetPollockGrid(pFgrid);
					}
					else {
						ExitGracefully("Parse:Bad Grid for pollock's method",BAD_DATA);
					}
				}
				else {ImproperFormat(s,l);}
				break;   
		  }
		case(701)://-----------------------------------------------------------------------------------
			{	//HeadWater
				//string "Headwater", x, y, inflow
				int junk;
				if (noisy)	{cout<<"Headwater"<<endl;}
				if (Len==4){CFlowNode::AddFlowNode(NULL,s_to_c(s[1],s[2]),FLAT,s_to_d(s[3]),junk);}
				else				{ImproperFormat(s,l);}
				break;
		  }
		case(702)://-----------------------------------------------------------------------------------
			{	//Coastal Aquifer
				//string "Coastal", Sea Level, Saltwater Specific Gravity
				if (noisy)	{cout<<"Coastal Aquifer"<<endl;}
				ExitGracefullyIf(!singlelayer,"Parse::Cannot simulate coastal multilayer aquifer",BAD_DATA);
				if (Len==3){CSingleLayer::SetCoastalInfo(s_to_d(s[1]),s_to_d(s[2]));eng.Grid.saltH20=true; }
				else				{ImproperFormat(s,l);}
				break;
		  }
		case(703)://-----------------------------------------------------------------------------------
			{ //Surface Water Pipeline
				//string "Pipeline", x1, y1, x2, y2
				if (noisy) {cout << "Surface Water Pipeline"<<endl;}
				if (Len==5){
					cmplex ends[2];
					double zeroes[2];
					CQSpec *pPipeline;
					ends[0]=s_to_c(s[1],s[2]);zeroes[0]=0.0;
					ends[1]=s_to_c(s[3],s[4]);zeroes[1]=0.0;
					cast=pPipeline = new CQSpec("pipeline",layer[cl],ends,1,zeroes);
					DynAddElement(cast,GIVEN_ELEM,AllElems[cl],ElemTypes[cl],elemsparsed[cl]);
					cout <<"Added Pipeline"<<endl;
				  //int junk;	cin>>junk;
				}
				else				{ImproperFormat(s,l);}
				break;
			}
		case(801)://-----------------------------------------------------------------------------------
			{	//Mesh Generation Only
				//string "MeshGenerateOnly"
        if (noisy)	{cout<<"Mesh Generate Only"<<endl;}
				eng.meshgenerate=true; 
				eng.solve			  =false;
				eng.grid			  =false;
				eng.track			  =false;
				eng.transport	  =false;	 
				spacetype=USER_SPECIFIED;
				break;
		  }
		case(802)://-----------------------------------------------------------------------------------
			{	//Mesh
			  //string "MeshBoundaries", string filename
				//MaxNodes (or <=zero if default)
				//{x y F} x (NumPoints)
				//&
				//{p1 p2} x (NumBoundarySegs) 
				//&
				//{p1 p2} x (NumInternalSegs) 
				//&
				if (noisy) {cout <<"Generate Mesh"<<endl;}
				eng.meshgenerate=true; 
				if (Len>=1){	
					thisname=strcpy(thisname,blank);
					for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
					pMesh=CTriMesh::DynamicParse(BBD,l,thisname);//builds initial triangulation
				}
				else {ImproperFormat(s,l);}
				break;
		  }
		case(803)://-----------------------------------------------------------------------------------
			{	//Mesh Generation
				//string "MeshDXFOutputOn"
				eng.DXFMeshOn=true; 
				break;
		  }
		case(804)://-----------------------------------------------------------------------------------
			{	//Mesh Generation
				//string "MeshBNAOutputOn"
				eng.BNAMeshOn=true; 
				break;
		  }
		case(805)://-----------------------------------------------------------------------------------
			{	//Courant Mesh Generation
				//string "CourantMeshGenerate" double dt, double minspacing, double maxspacing, double radius
				if (Len==5){
					eng.meshgenerate=true; 
					spacetype=COURANT;
					spaceparam=s_to_d(s[1]);spacemin=s_to_d(s[2]);spacemax=s_to_d(s[3]);spacerad=s_to_d(s[4]);
				}
				break;
		  }
		case(806)://-----------------------------------------------------------------------------------
			{	//Mesh Generation (in addition to solve, grid, track)
				//string "MeshGenerate"
				eng.meshgenerate=true;  
				spacetype=USER_SPECIFIED;
				break;
		  }
		case(901)://-----------------------------------------------------------------------------------
			{//Voronoi
			 if (noisy)	{cout<<"Voronoi Test"<<endl;}
			 if (Len>=1){
				 thisname=strcpy(thisname,blank);
			   for(i=1; i<Len; i++){thisname=strcat(thisname,s[i]);thisname=strcat(thisname,space);}
				 CTriMesh *pMesh;
				 pMesh=CTriMesh::Parse(BBD,l,thisname);
				 pMesh->WriteGeometry(eng.DXFMeshOn,eng.BNAMeshOn);
				 delete pMesh;
				 ExitGracefully("End Voronoi Test",OTHER);
			 }
			 else				{ImproperFormat(s,l);}
		   break;
		  }
	  case(-1)://------------------------------------------------------------------------------------
		  {//comment out
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
				if (!elemdisabled) {cout <<"flow input file: command unrecognized at line "<< l << endl;} 
				break;   
		  }
		case(-10)://------------------------------------------------------------------------------------
		  {//debug mode on
				if (noisy) {cout <<"Debug Mode on" << endl;} 
				eng.debug=true;
				eng.Grid.pot=true;
				eng.Grid.QxQy=true;
				break;   
		  }
		default: {
			  if (noisy) {cout << "##"<<endl;} 
				break;
			}//{cout<< "BAD BBL FILE INPUT, line "<< l+1 <<endl;return false;}
	}//switch(code)
  }//while (tokenize...)
  cout << endl;

  //*******************************************************************************************
	//                       DONE PARSING- CREATE & BUILD SYSTEM
  //*******************************************************************************************
	if (eng.solve || eng.grid || eng.track || eng.transport)
	{
		
		//------------------------------------------------------------------------------------
		//Create Aquifer
		//------------------------------------------------------------------------------------
		aq = new CAquifer("main aquifer");

		//------------------------------------------------------------------------------------
		//check for agreement of level tops & bottoms
		//------------------------------------------------------------------------------------
		for (L=0; L<numlevels-1; L++){
			if ((aqthickness[L]+aqbase[L])!=aqbase[L+1]) {
				cout << "Layer thicknesses and bases do not agree at level "<< L << ". Please repair."<<endl;
				ExitGracefully("Parse:Layer thicknesses and bases do not agree at level",BAD_DATA);
			}
		}
		//Should make sure that reference points for each layer are in same location

		//------------------------------------------------------------------------------------
		//Do Layer/Multilayer stuff
		//------------------------------------------------------------------------------------
    //TMP DEBUG- ML MIGRATION - should do this every time "Next Layer" Command is called or if input layer is finished
		
		// should have new routine for aquifer "AppendLayerSetToAquifer(CSingleLayer or CMultiLayer)", called at end
		for (L=0; L<numlevels; L++){
      
			if (true){//(IsSingleLayer[L]){
				cout << "Building Single Layer #"<<L+1<<" ("<<elemsparsed[L]+transelems[L]<< " elements)"<<endl;

				//add elements to each layer, transient first
				layer[L]->AddToLayer(ffp[L]);
				for (i=0; i<transelems[L];  i++){
					layer[L]->AddToLayer(TransientElems[L][i],TRANSIENT_ELEM ); 
					TransientElems[L][i]->WriteGeometry(BASEMAP);
				}
				for (i=0; i<elemsparsed[L]; i++){
					layer[L]->AddToLayer(AllElems      [L][i],ElemTypes[L][i]); 
					AllElems      [L][i]->WriteGeometry(BASEMAP);
				}
				
				//add property zones to each layer
				for (i=0; i<numpzones[L]; i++)  {layer[L]->AddToLayer(PropZones[L][i]);}
				
				//set black hole 
				layer[L]->SetBlackHole();

				//Set layer values
				layer[L]->SetBackgroundValues(aqbase[L], aqthickness[L],aqconductivity[L],aqporosity[L], L);

				//For Multilevel aquifer,
				//associate each with leaky layers, set leaky layer domain size/location
				if (numlevels>1) {
					cout << "Building Aquiclude #"<<L+1<<" ("<<leakelems[L]<< "elems)"<<endl;
					if (L!=numlevels-1){layer[L]->SetLevelAbove(aquicludes[L+1]);}
	  			for (i=0; i<leakelems[L]; i++)   {aquicludes[L]->AddToLayer(AllLeakElems [L][i]);AllLeakElems[L][i]->WriteGeometry(BASEMAP);}
					for (i=0; i<numlpzones[L]; i++)  {aquicludes[L]->AddToLayer(LeakPropZones[L][i]);}

					layer[L]->SetLevelBelow(aquicludes[L]);
				}
				else {
					aquicludes[L]=NULL;
				}

				aq->AddSingleLayer(layer[L],aquicludes[L]); 
			}
			else{
				cout << "Building MultiLayer #"<<L+1<<" ("<<mlelemsparsed[L]<< " elements)"<<endl;

				//add elements to each layer
				mlayer[L]->AddToLayer(MLffp[L]);
				for (i=0; i<mlelemsparsed[L]; i++){
					mlayer[L]->AddToLayer(AllMLElems[L][i],MLElemTypes[L][i]);
					AllMLElems[L][i]->WriteGeometry(BASEMAP);
				}
				
				//add property zones to each layer
				for (i=0; i<nummlpzones[L]; i++)  {mlayer[L]->AddToLayer(MLPropZones[L][i]);}
				
				//set black hole location 
				mlayer[L]->SetBlackHole();

				//Set multilayer values
				//mlayer[L]->SetBackgroundValues(aqbase[L], aqthickness[L],aqconductivity[L],aqporosity[L], L);

				//For Multilevel aquifer,
				//associate each with leaky layers, set leaky layer domain size/location
				if (numlevels>1) {
					cout << "Building Aquiclude #"<<L+1<<" ("<<leakelems[L]<< "elems)"<<endl;
					if (L!=numlevels-1){layer[L]->SetLevelAbove(aquicludes[L+1]);}
	  			for (i=0; i<leakelems[L]; i++)   {aquicludes[L]->AddToLayer(AllLeakElems [L][i]);AllLeakElems[L][i]->WriteGeometry(BASEMAP);}
					for (i=0; i<numlpzones[L]; i++)  {aquicludes[L]->AddToLayer(LeakPropZones[L][i]);}

					mlayer[L]->SetLevelBelow(aquicludes[L]);
				}
				else {
					aquicludes[L]=NULL;
				}

				aq->AddMultiLayer(mlayer[L],aquicludes[L]);

			}
		}
		//------------------------------------------------------------------------------------
		//Sift through Inhomogeneities, turn off correct sides
		//------------------------------------------------------------------------------------
		CInhom    ::SiftThroughInhoms();
		CBaseInhom::SiftThroughInhoms();

		//------------------------------------------------------------------------------------
		//Do Superblock stuff
		//------------------------------------------------------------------------------------
		if (blocksolve) {
			cout << "Filling Superblocks..." <<endl;
			cout << "----------"<<endl;
			CSuperblock *Master;
			for (L=0; L<numlevels; L++){
				if (true){//IsSingleLayer[L]{
					Master=CSuperblock::FillSuperblocks(layer[L]->GetExtents(),
																							SuperNest,
																							AllElems [L],elemsparsed[L],
																							PropZones[L],numpzones  [L],
																							blockprec,blockimplicit,blocknested);
					layer[L]->SetMasterBlock(Master);
					if (eng.warm){
						cout <<"Updating Superblock Expansions:"<<endl;
						for (i=0; i<elemsparsed[L]; i++){
							WriteEllipsisToScreen(i,elemsparsed[L],50);
							AllElems[L][i]->UpdateBlock(eng.warmtime);
						}
						cout <<endl<<endl;
					}
				}
			}
			cout <<endl<< "                   ...done filling superblocks."<<endl<<endl;
		}

		for (i=0; i<numpart; i++){
			particles[i]->SetAquifer(aq);
		}

		CFlowNode::BuildFlowNetwork();

	} //end if (eng.solve || eng.grid || eng.track || eng.transport)

	if (eng.meshgenerate) {
		ExitGracefullyIf((pMesh==NULL),"Parse:NULL Mesh for generation: must specify mesh bounds",BAD_DATA);
		pMesh->DynamicInitialize(spacetype,spaceparam,spacemin,spacemax,spacerad); 
	}

  BASEMAP.close();
  SOLUTION.close();
	BBD.close();

	//for (i=0; i<MAXINPUTITEMS; i++){delete [] s[i];}
	//delete [] s;
  //delete [] thisname;

	for (L=0; L<numlevels; L++){
		delete [] AllElems      [L]; //just deletes pointers, not elements
		delete [] TransientElems[L];
		delete [] PropZones     [L]; 		
		
		delete [] AllLeakElems  [L];		
		delete [] LeakPropZones [L];

    delete [] AllMLElems    [L];
	  delete [] MLPropZones   [L];
		
		if (elemsparsed[L]>0)	 {delete [] ElemTypes[L];}
		if (mlelemsparsed[L]>0){delete [] MLElemTypes[L];}
	}

  return true;
}
/**********************************************************************
     DynAddElement
-----------------------------------------------------------------------
Dynamically adds element xptr to array pElems of size "size"
also dynamically adds new etype to types[]
**********************************************************************/
bool DynAddElement(CAnalyticElem *xptr, elementtype etype, CAnalyticElem **&pElems,elementtype *&types, int &size){
	CAnalyticElem **tmp=NULL;
  elementtype    *tmp2=NULL;

  if (xptr==NULL){return false;}
	if ((pElems==NULL) && (size>0)) {return false;}
  size=size+1;																								               //increment size
  tmp =new CAnalyticElem *[size];					                                   //allocate memory 
	tmp2=new elementtype    [size];
	if (tmp2==NULL){ExitGracefully("DynAddElement::Out of memory",OUT_OF_MEMORY);}

  for (int i=0; i<(size-1); i++){
		tmp[i]=pElems[i];
		tmp2[i]=types[i];
	}	                                                                         //copy arrays
  tmp [size-1]=xptr;																		                     //add new element
	tmp2[size-1]=etype;                                                        //add new type
  if (size>1){delete [] pElems;delete [] types;}													   //delete old arrays of pointers														
	pElems=tmp;																				                         //redirect pointers
	types =tmp2;
	return true;
}
/**********************************************************************
     DynAddPZone
-----------------------------------------------------------------------
Dynamically adds prop. zone xptr to array pZones of size "size"
**********************************************************************/
bool DynAddPZone  (CPropZone *xptr, CPropZone **&pZones, int &size){
	CPropZone **tmp=NULL;

  if (xptr==NULL){return false;}
	if ((pZones==NULL) && (size>0)) {return false;}
  size=size+1;																								               //increment size
  tmp =new CPropZone *[size];					                                       //allocate memory 
	if (tmp==NULL){ExitGracefully("DynAddPropZone::Out of memory",OUT_OF_MEMORY);}

  for (int i=0; i<(size-1); i++){
		tmp[i]=pZones[i];
	}	                                                                         //copy arrays
  tmp [size-1]=xptr;																		                     //add new zone
  if (size>1){delete [] pZones;}													                   //delete old arrays of pointers														
	pZones=tmp;																				                         //redirect pointers
	return true;

}