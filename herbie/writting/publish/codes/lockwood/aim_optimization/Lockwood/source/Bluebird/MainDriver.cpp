#include "Bluebird.h"
#include "SparseMatrix.h"

/************************************************************
           MAIN DRIVER
************************************************************/

static CAquifer       *Aquifer=NULL;										//pointer to Aquifer
static CChemDomain2D  *TransDomain=NULL;								//pointer to Transport Domain
static CPathline      *Part[MAX_PARTICLES]={NULL};      //array of pointers to particles
static CTriMesh       *Mesh=NULL;                       //pointer to mesh

/************************************************************
 MAIN
	Primary routine
------------------------------------------------------------*/
int main(int argc, char* argv[]){ 

  char           *file=new char[255];	
	int             i;													//counters;

  if (argc>=2)   {
		for (i=1; i<argc-2; i++){     //assumes filename is last
			if (!strcat(argv[i],"-p")){}//pause=true;}
		}
		file =argv[argc-1];cout <<"Input File: "<<file<<endl;
	}
  else           {file="split.dat";                                  }

	EngineCommands  e;													//engine commands from parser

	double          t(0),tstep(ALMOST_INF);		  //time; timestep for transient solution
  int             nlayers(0),numpart(0);			//number of layers; number of particles;
	int             nZB(0);                     //number of zone budget polygons

  clock_t         t1,t2;											//computational time markers

	Aquifer     = NULL;                         //initialize aquifer, transdomain, particles
	TransDomain = NULL;						     
	for (i=0;i<MAX_PARTICLES;i++)     {Part  [i]=NULL;}

  CStringElem::Prepare(); 
	CLineElem::Prepare();
	InitializeRxNLibrary(ExitGracefully);
	randinit(-1);

	ofstream PROGRESS;
	PROGRESS.open("progress.out");
	ofstream SOL;
 
	/*cmplex z=cmplex(2,3);
	cmplex tmp;
	t1=clock();
	for (i=0; i<1000; i++){tmp=2*arctanhz);}
	t2=clock();
	cout <<" ..."<<float(t2-t1)/CLK_TCK << " seconds elapsed (tanh). " << endl;
  t1=clock();
	for (i=0; i<1000; i++){tmp=log((z+1)/(z-1));}
	t2=clock();
	cout <<" ..."<<float(t2-t1)/CLK_TCK << " seconds elapsed (log). " << endl;
	ExitGracefully("Tanh Test",RUNTIME_ERR);*/

  /*TestCationExchange();
	ExitGracefully("Cation Exchange Test",RUNTIME_ERR);*/
	
	/*TestParentDaughter();
  ExitGracefully("Parent Daughter Exchange Test",RUNTIME_ERR);*/

	cout <<"**********************************************************"<<endl;
	cout <<"            Bluebird / Cardinal version 3.00              "<<endl;
	cout <<"----------------------------------------------------------"<<endl;
	cout <<"  an object-oriented simulation tool for simulation of    "<<endl;
	cout <<"  groundwater flow and reactive contaminant transport with"<<endl;
	cout <<"  the analytic element method (AEM).                      "<<endl;
	cout <<"                                                          "<<endl;
	cout <<"  Author: James R. Craig, Ph.D.                           "<<endl; 
	cout <<"**********************************************************"<<endl<<endl;

  //PARSE/CREATE
  //-----------------------------------------------------------------
  if ((true) && (Parse(file,Aquifer,Part,numpart,e,Mesh))){
		
		//-----------------------------------------------------------------
		// GENERATE MESH
		//-----------------------------------------------------------------
		if (e.meshgenerate) 
		{
			if (Mesh!=NULL)
			{
				Mesh->DynamicBuildMesh(true,true);
				Mesh->WriteGeometry(e.DXFMeshOn,e.BNAMeshOn);//writes bna,dxf
				Mesh->WriteOutput2();//writes output .btm file
			}	
		}

		//Parse Transport File
		//---------------------------------------------------------------------
		if (e.transport)
		{
			CardinalParse("cardinal.dat",TransDomain,Aquifer,Mesh);
		}
		if (Aquifer==NULL && (e.solve || e.grid || e.track || e.transport))
		{
			ExitGracefully("Main Bluebird/Cardinal Driver:Aquifer is NULL",BAD_DATA);
		}
		else if (Aquifer!=NULL){
			cout <<"File Parsed- Aquifer Created with " <<Aquifer->GetNumLayers()<<" layers"<<endl;
			cout << "------------------------------------------------------------"<<endl;
		}

		//Parse Observation File or Head.dat
		//---------------------------------------------------------------------
    if (e.obsfileexists)
		{
			if (!ObservationsParse("observations.dat",Aquifer)){
				OldHeadFileParse(Aquifer); //if head.dat does not exist, exits gracefully
			}
		}

		for (int timespassed=0; timespassed<e.nTimes; timespassed++){
			do {
				cout << "**Time: " << t <<"**"<<endl;
				//-----------------------------------------------------------------
				// SOLVE
				//-----------------------------------------------------------------
				if (e.solve)
				{
					t1=clock();
					if (!e.explicitsolve){Aquifer->IterativeSolve   (t,PROGRESS);}
					else                 {Aquifer->IterExplicitSolve(t,PROGRESS);}
					t2=clock();
					cout <<" ..."<<float(t2-t1)/CLK_TCK << " seconds elapsed (Solve). " << endl;
					if (e.writesolvetime){ofstream TIME;TIME.open("solvetime.txt");TIME<<float(t2-t1)/CLK_TCK<<endl;TIME.close();}
				}
				//-----------------------------------------------------------------
				// TRACK (transient)
				//-----------------------------------------------------------------
				if ((e.track) && (e.transient))
				{
					CPathline::TrackAllPathlines(Part,numpart,min(tstep,max(0.0,CPathline::GetDuration()-t)),PROGRESS); 
				}
				tstep=min(e.timestep,e.outtimes[timespassed]-t);
				t+=tstep;
			}while (t<e.outtimes[timespassed]);

			cout <<"Output Time #" <<timespassed+1<<" - t="<<t<<endl;

			if (e.interpolate){
				Aquifer->InterpolateQxQy(t);
			}
			//-----------------------------------------------------------------
			// GRID
			//-----------------------------------------------------------------
			if (e.grid)
			{
				t1=clock();
				GridDriver(Aquifer,e.Grid,t,PROGRESS);
				t2=clock();
				cout <<" ..."<<float(t2-t1)/CLK_TCK << " seconds elapsed (Grid). " << endl;
			} 
			//-----------------------------------------------------------------
			// TRACK (steady state)
			//-----------------------------------------------------------------
			if ((e.track) && (!e.transient))
			{
				if (CPathline::GetPollockGrid()!=NULL){
					CPathline::GetPollockGrid()->Initialize(Aquifer->GetLayer(0));
				}
				CPathline::TrackAllPathlines(Part,numpart,CPathline::GetDuration(),PROGRESS); 
			}
			//-----------------------------------------------------------------
			// WRITE SNAPSHOT SOLUTION (for specific time)
			//-----------------------------------------------------------------
			if ((e.solve) && (e.writesol)) 
			{
				cout << "Writing solution file..."<<endl;
				SOL.open("solution.bbs"); SOL.precision(10);
				Aquifer->WriteItself(SOL,t); 
				SOL.close();
				cout << "      ...Solution file written"<<endl;			
			} 
			//-----------------------------------------------------------------
			// WRITE SNAPSHOT OUTPUT FILES (for specific time)
			//-----------------------------------------------------------------
			if ((e.solve) && (e.writeout))//Avoids writing after gridonly or trackonly (waste)
			{ 
				cout << "Writing output files..."<<endl;
				ClearOutputFiles();
				Aquifer->WriteOutput(t);
				CAnalysisLocation::CalculateAndWriteAll(t);
				cout << "      ...Output files written"<<endl;
			}
			timespassed++;

		} //end for (int timespassed...
		
		//-----------------------------------------------------------------
		//TRACKING OUTPUT
		//-----------------------------------------------------------------	
		if (e.track)
		{
			CPathline::WriteAllPathlines(Part,numpart); 
		}

		//-----------------------------------------------------------------
		//TRANSPORT (non-transient solution only)
		//-----------------------------------------------------------------
		if (e.transport)
		{
			ClearTransportFiles();
			t1=clock();
			ExitGracefullyIf(TransDomain==NULL,"Bluebird/Cardinal Main Driver:Transport Domain is NULL",BAD_DATA);
			TransDomain->Transport(PROGRESS);
			t2=clock();
			cout <<" ..."<<float(t2-t1)/CLK_TCK << " seconds elapsed. " << endl;
		}
		//-----------------------------------------------------------------
		//SET UP SOCKET CONNECTION
		//-----------------------------------------------------------------
		/*if (e.socket){
			SocketConnect(Aquifer);
		}*/
		//-----------------------------------------------------------------
		//PRINT GENERAL INFORMATION
		//----------------------------------------------------------------- 
		if (e.solve || e.grid || e.track || e.transport)
		{
			cout <<endl << "General Solution Information: "<<endl;
			cout << "------------------------------------------------------------"<<endl;
			cout << "Total Elements: "     << CAnalyticElem::GetTotalElements() <<endl;
			//if (CSuperblock::SuperblocksOn()){
			//cout << "Total Segments: "     << Aquifer->GetLayer(0)->GetMasterBlock()->GetPopulation() <<endl;
			//}
			cout << "net discharge: "      << Aquifer->GetNetDischarge(0.0)        << endl;
		}

  }
	else {ExitGracefully("Parser failed-cannot find input file",NO_INPUT_FILE);}

	PROGRESS<<"finished"<<endl;
	PROGRESS.close();
 	ExitGracefully("Successful completion :)",COMPLETED);

	return 0;
};

/******************************************************************************
 GLOBAL DELETE 
	Destroys all objects not local to function which called ExitGracefully 	
	Insures memory deletion
------------------------------------------------------------------------------*/
void GlobalDelete(bool pause){
	CStringElem      ::Destroy(); 
	CInhom           ::Destroy();
	CBaseInhom       ::Destroy();
	CAnalyticElem    ::DestroyAllElements();
	CPropZone        ::DestroyAllPropZones();
	CAnalysisLocation::DestroyAllAnalysisLocations();
	CFlowNode        ::Destroy();
	CPathline        ::Destroy();
	CTimeSeries      ::DestroyAllTimeSeries();
	delete Aquifer;
  delete TransDomain;	
	delete Mesh;
	delete idum;
	for (int i=0; i<MAX_PARTICLES; i++){delete Part[i];}
  ofstream PROGRESS;
	PROGRESS.open("progress.out",ios::app);
	PROGRESS.close();
	  if(false){
	  //if (pause){
	  //if (true){
		int L; cin >> L; 
	}


	exit(0);
}
//******************************************************************************
void   ClearOutputFiles(){
	ofstream EXTRACT;
	EXTRACT.open("extract.csv");
	EXTRACT << "x,y,extract,cum_extract"<<endl;
  EXTRACT.close();

	ofstream ERRORS;
	ERRORS.open("errors.csv");
  ERRORS << "x,y,type,element_type,abs_errors,rel_errors"<<endl;
  ERRORS.close();

	ofstream ZB;
	ZB.open("ZoneBudget.csv");
	//ZB << "ID,influx,outflux,net_extraction"<<endl;
  ZB.close();

	ofstream TRANSECT;
	TRANSECT.open("TransectAnalysis.csv");
  TRANSECT.close();

	/*ofstream TMPINJ; //TMP DEBUG
	TMPINJ.open("tmpInj.txt");
	TMPINJ.close();*/
	//ofstream DISCONT; //TMP DEBUG
	//DISCONT.open("discontinuity.bna");
	//DISCONT.close();

	ofstream OBSOUT;
	OBSOUT.open("obs_errors.csv");
	OBSOUT.close();
}
//******************************************************************************
void ClearTransportFiles(){


	//TRACK.open("tracks.dat");//in CPathline Output  
	//TRACK.close();

	//ofstream CONC; //in CPathline Output  
	//CONC.open("conc.dat");
	//CONC.close();

	/*
	ofstream TRACK;   //TMP DEBUG -for writing streamlines (overwrites particles)
	TRACK.open("tracks.bna");
	TRACK.close();	

	ofstream STREAMOUT;
	STREAMOUT.open("streamlines_on_2D_grid.csv");
	STREAMOUT << "cell,TOF ,x, y, xc, yc, 2D_i, 2Dx, 2Dy, q ,conc[0]"<<endl;
	STREAMOUT.close();

	STREAMOUT.open("streamline 1D grid.csv");
	STREAMOUT << "type,ind, TOFa, TOFb, TOFave,v, Conc[0], q "<<endl;
	STREAMOUT.close();*/
	
}



	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//TEMP DEBUG- superblock testing -	before	for (int timespassed=0;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	/*int dbgCtr=0;
	int thisorder;
  double thisfold;
	ofstream DBGBLOCKS;
	DBGBLOCKS.open("superblocktimes.txt");
	DBGBLOCKS<< "Elements: "<< numelems[0]<<endl;
	DBGBLOCKS<< "Segments: "<< layer[0]->GetMasterBlock()->GetPopulation()<<endl;
	DBGBLOCKS<< "Iterations: " << CSingleLayer::MaxLayerIterations <<endl;
  DBGBLOCKS<< "Grid_res: " << gridprec <<endl;
	DBGBLOCKS<< "system_prec: "<<endl;
	DBGBLOCKS<<"nesting solvetime gridtime"<<endl;
  for (int nest=0; nest<6; nest+=5){
		CSuperblock::SetPrecision(4,thisorder,thisfold);
		cout << "Filling Superblocks...nesting level: "<<nest <<endl;
		cout << "------------------------------------------------------------"<<endl;
	  CSuperblock::FillSuperblocks(layer[0],nest,E,AllElements[0], numelems[0],thisorder,thisfold);
		cout <<endl<< "                   ...done filling superblocks."<<endl<<endl;
    DBGBLOCKS<< nest<<" ";*/
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


			//TEMP DEBUG- superblock testing - before StringElem::Destroy
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*}
	DBGBLOCKS.close();*/
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




	//TMP DEBUG - Sparse Matrix Debug
	/*CSparseMatrix *S;
	S =new CSparseMatrix();
	int j;
	
	S->DynamicInitialize(6);

	double **A;
	double  *B;
	double  *x;
	A=new double *[6];
	B=new double [6];
	x=new double [6];
	for (i=0; i<6; i++){
		B[i]=1.0;
		x[i]=0.2;
		A[i]=new double [6];
		for (j=0; j<6; j++){
			A[i][j]=0.0;
		}
	}

	for (i=0; i<6; i++){
		for (j=-1; j<=1; j++){
			if (((i+j)>=0) && ((i+j)<6)){
				//if (!S->SetAij(1.0,i,i+j)){cout <<"BAD"<<endl;}
				S->DynamicAdd(1.0,i,i+j);
				A[i][i+j]=1.0;
			}
		}
A[i][i]=2.0;
S->DynamicAdd(2.0,i,i);
	}
	for (i=0; i<6; i++){
	  for (j=0; j<6; j++){
			cout << S->GetAij(i,j)<<" ";
		}
		cout <<endl;
	}
	cout <<endl;

	for (i=0; i<6; i++){
	  for (j=0; j<6; j++){
			cout << A[i][j]<<" ";
		}
		cout <<endl;
	}
	cout <<endl;

	int junk(0);
	Gauss(A,B,x,6,junk);
	for (i=0; i<6; i++){cout<<x[i]<<endl;}cout<<endl;

//enum BCGtestparam{RELATIVE_RESIDUAL,RELATIVE_TRANS_RESIDUAL,NORM_ERROR,MAX_ERROR};
	for (i=0; i<6; i++){x[i]=0.2;B[i]=1.0;}
	double tmp_err(0.0);
	BCG(*S,B,x,6,MAX_ERROR,tmp_err);
	for (i=0; i<6; i++){cout<<x[i]<<endl;}cout<<endl;

delete [] B;
delete [] x;
for (i=0; i<6; i++){delete [] A[i];}
delete [] A;

	ExitGracefully("Done",BAD_DATA);*/


/*ofstream BESSEL;
	BESSEL.open("bessel.csv");
	int factorial(1);
	for (int k=0; k<40;k++){
		double sum(0.0);
		for (int m=1; m<=k;m++){
			sum+=1.0/(double)(k);
		}
		if (k!=0){factorial*=k;}
		BESSEL<<k<<","<<((IM*pi)/2.0-log(0.5*Im)-Euler+0.5*sum)*(pow(Im,2.0*k))*(pow(-0.25,1.0*k)/(double)(factorial*factorial)) <<","<<(pow(Im,2*k))*(pow(-0.25,k)/(double)(factorial*factorial)) <<endl;
	}
	BESSEL.close();
ExitGracefully("Bessel coeff Test",RUNTIME_ERR);*/
	/*ofstream BESSEL;
	BESSEL.open("bessel.csv");
	for (double x=0.0; x<5.0;x+=0.1){
		BESSEL<<x<<","<<bessel_I0(x)<<endl;
	}
	BESSEL<<"KKKKK"<<endl;
	for ( x=0.0; x<3.0;x+=0.1){
		BESSEL<<x<<","<<bessel_K0(x)<<endl;
	}
	BESSEL.close();

  ofstream HANTUSH;
	HANTUSH.open("hantush.csv");

	for (double b=-2.0; b<=2.0;b+=0.25){
		HANTUSH<<pow(10.0,b)<<",";
	}HANTUSH<<endl;

	for (double u=-5.0; u<=1.0;u+=1.0){
		HANTUSH<<1.0/pow(10.0,u)<<",";
		for (double b=-2.0; b<=2.0;b+=0.25){
			HANTUSH<<hantush(pow(10.0,u),pow(10.0,b))<<",";
		}
		HANTUSH<<endl;
	}
	HANTUSH.close();*/




