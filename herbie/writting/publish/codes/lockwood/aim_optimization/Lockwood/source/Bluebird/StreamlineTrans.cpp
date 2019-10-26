//StreamlineTrans.cpp
#include "TransportScheme.h"
#include "Streamline.h"
/************************************************************************
                           Constructors
************************************************************************/
CStreamlineTrans::CStreamlineTrans(C2DDomainABC *pDom)
                 :C2DTransportScheme(pDom){
	nstreamlines=0; 
  nCells      =0;
	lineIDs     =NULL;
	linePart    =NULL;
	nLinesInCell=NULL;
	Recharge    =NULL;
	pStreamlines=NULL;
	update_type =UPDATE_DIRECT; //default
}
//------------------------------------------------------------------------
CStreamlineTrans::~CStreamlineTrans(){

	if (globaldebug){cout << "  DESTROYING STREAMLINE TRANSPORT STYLE"<<endl;}

	if (lineIDs !=NULL){for (int k=0;k<nCells; k++){delete [] lineIDs  [k];}}
	if (linePart!=NULL){for (int k=0;k<nCells; k++){delete [] linePart [k];}}
  delete [] lineIDs;     
	delete [] linePart; 
	delete [] nLinesInCell;
  delete [] Recharge; 
	for (int p=0; p<nstreamlines; p++){delete pStreamlines[p];} delete [] pStreamlines;
}
/************************************************************************
                           MANIPULATORS
************************************************************************/
void CStreamlineTrans::SetParameters(const double ns, const str_update_type type){
	if (ns<=0) {ExitGracefully("CStreamlineTrans::SetNumStreamlines: negative number of streamlines ",BAD_DATA);}
	nstreamlines=ns;
	update_type=type;
}
/************************************************************************
                           Initialize
************************************************************************/
void CStreamlineTrans::Initialize(const double         startt,
																	const transport_type ty,
									                Ironclad1DArray      porosity, 
									                Ironclad1DArray      satthick){

	int           j,k,p;
	int						pieces[MAX_CROSSES_OF_CELL];
	double				TOFs  [MAX_CROSSES_OF_CELL];
	int						ntimesincell(0);
	disp_info     junk;

	process_type=ty;

	pLayer   =pDomain  ->GetLayer();
	pConcGrid=pDomain  ->GetGrid(); 
	nspecies =pDomain  ->GetNumSpecies();
	nCells   =pConcGrid->GetNumCells();

	//Allocate memory,initialize----------------------------------------------------------------
	nLinesInCell=new int     [nCells];
  lineIDs     =new int    *[nCells];
	linePart    =new int    *[nCells];
	lineTOF     =new double *[nCells];
	totalTOF    =new double  [nCells];
	Recharge    =new double  [nCells];
	if (Recharge==NULL){
		ExitGracefully("CStreamlineTrans::Initialize; out of memory",OUT_OF_MEMORY);}
	for (k=0;k<nCells; k++){
		lineIDs     [k]=new int    [MAX_LINES_IN_CELL]; //should make dynamic?
		linePart    [k]=new int    [MAX_LINES_IN_CELL];
		lineTOF     [k]=new double [MAX_LINES_IN_CELL];
	  if (lineTOF [k]==NULL){
			ExitGracefully("CStreamlineTrans::Initialize; out of memory(2)",OUT_OF_MEMORY);}
		nLinesInCell[k]=0;
		totalTOF    [k]=0;
		Recharge    [k]=0;
		for (j=0;j<MAX_LINES_IN_CELL; j++){
			lineIDs   [k][j]=0;
			linePart  [k][j]=0;
			lineTOF   [k][j]=0.0;
		} 
	}

  //Create Initial Streamlines--------------------------------------------------------------
	//----------------------------------------------------------------------------------------	
	CStreamline::CreateStreamlines(pLayer,pDomain,pStreamlines,nstreamlines,startt);

  if (pStreamlines==NULL){
		ExitGracefully("CStreamlineTrans::Initialize: unsuccessful streamline creation ",BAD_DATA);}
	

	//Track Streamlines (until captured or out of grid)---------------------------------------
	//----------------------------------------------------------------------------------------
	for(p=0; p<nstreamlines; p++){
		if (pStreamlines[p]==NULL){
			ExitGracefully("CStreamlineTrans::Initialize:NULL Pathline ",BAD_DATA);}
		if (ProgramAborted()){nCells=0;nstreamlines=0;break;}

		cout <<endl<<"Tracking streamline #"<<p;


		pStreamlines[p]->Track(ALMOST_INF, advectiontype,adv_timestep,adv_spacestep,false,junk,false);
		pStreamlines[p]->CleanPath   ();
    pStreamlines[p]->CPathline::WriteOutput(); //TMP DEBUG- writes bna
		pStreamlines[p]->Create1DGrid();
		pStreamlines[p]->WriteOutput ();
		//Link 2D grid to Streamline[p]---------------------------------------------------------
		for (k=0;k<nCells; k++){ //not the fastest scheme

			pStreamlines[p]->Get2Dinfo(k,ntimesincell,pieces,TOFs);

			for (j=0;j<ntimesincell;j++){

				if (nLinesInCell[k]>MAX_LINES_IN_CELL-1){
					ExitGracefully("CStreamlineTrans::Initialize: Exceeded MAX_LINES_IN_CELL",TOO_MANY);}

				linePart    [k][nLinesInCell[k]]=pieces[j];
				lineIDs     [k][nLinesInCell[k]]=p;
				lineTOF     [k][nLinesInCell[k]]=TOFs[j];
				totalTOF    [k]+=TOFs[j];
				nLinesInCell[k]++;

			}
		}
	}

	//Check all cells, if no streamline exists, trace backwards to cell with nearest real streamline
		//----------------------------------------------------------------------------------------
	/*
  CStreamline  *pNewStreamline;
	for (k=nCells-1;k>=0; k--){//ordering matters... flow left to right, nCells to 0, flow right to left, 0 to nCells
		if (nLinesInCell[k]==0){ //and surrounded by nLinesincell>0??
			pNewStreamline=new CStreamline(p,pAquifer,pDomain, pConcGrid->GetCellCenter(k),30); //TMP DEBUG- require some kind of flux!!!
			pNewStreamline->SetDirection(TRACK_BACKWARD);
			DynArrayAppend((void**&)(pStreamlines),(void*)(pNewStreamline),nstreamlines);
			
			p=nstreamlines-1;
			if (pStreamlines[p]==NULL){ExitGracefully("CStreamlineTrans::Initialize:NULL Pathline (2)",BAD_DATA);}
		  cout <<endl<<"Backtracking streamline #"<<p;
			if (ProgramAborted()){nCells=0;nstreamlines=0;break;}
			
			pStreamlines[p]->Track(ALMOST_INF, advectiontype,adv_timestep,adv_spacestep);
			pStreamlines[p]->CleanPath();
			pStreamlines[p]->ReversePath();
			pStreamlines[p]->CPathline::WriteOutput(); //TMP DEBUG- writes bna
			pStreamlines[p]->Create1DGrid();
			pStreamlines[p]->WriteOutput();
	    
			//Link 2D grid to Streamline[p]---------------------------------------------------
			for (int i1=0;i1<nCells; i1++){
				pStreamlines[p]->Get2Dinfo(i1,ntimesincell,pieces,TOFs);
				for (j=0;j<ntimesincell;j++){
					if (nLinesInCell[i1]>MAX_LINES_IN_CELL-1){
						ExitGracefully("CStreamlineTrans::Initialize: Exceeded MAX_LINES_IN_CELL",other);
					}
					linePart    [i1][nLinesInCell[i1]]=pieces[j];
					lineIDs     [i1][nLinesInCell[i1]]=p;
					lineTOF     [i1][nLinesInCell[i1]]=TOFs[j];
					totalTOF    [i1]+=TOFs[j];
					nLinesInCell[i1]++;
				}
			}
		}
	}*/

	//Calculate/Assign Recharge to Streamlines (may have to be during solution)
	//----------------------------------------------------------------------------------------
	for (k=0;k<nCells; k++){  
		Recharge[k]=pLayer->GetLeakage(c3Dto2D(pConcGrid->GetNodeLocation(k)),startt,FROMTOPANDBOTTOM);
		for (j=0; j<nLinesInCell[k]; j++){
			pStreamlines[lineIDs[k][j]]->SetFlux((Recharge[k]*lineTOF[k][j]/totalTOF[k]),linePart[k][j]);
		}
	}
	for(p=0; p<nstreamlines; p++){
		pStreamlines[p]->CleanFluxes();
	}

}
/************************************************************************
                           TRANSPORT
************************************************************************/
void CStreamlineTrans::Transport(const double     t,
																 const double     tstep,
																 Ironclad2DArray  C,
																 Ironclad2DArray  Cend,
																 //Ironclad2DArray  dCdt_sim
																 Writeable2DArray Cnew){

	int		 k,s,p;
	double tmpconc[MAX_SPECIES],sum[MAX_SPECIES];	

	if (sum==NULL){ExitGracefully("CStreamlineTrans::Transport: Out of memory",OUT_OF_MEMORY);}

	for (p=0; p<nstreamlines; p++){
		WriteEllipsisToScreen(p,nstreamlines,20);
		if (ProgramAborted()){return;}
		pStreamlines[p]->UpdateConcentrations(tstep,t,update_type);
	}
	
	for (k=0;k<nCells;k++){
		
		for (s=0;s<nspecies;s++){sum[s]=0.0;}

		//get concentrations from pStreamlines, average by TOF
		for (p=0;p<nLinesInCell[k];p++){
      pStreamlines[lineIDs[k][p]]->GetConcentrations(linePart[k][p],tmpconc);
			for (s=0;s<nspecies;s++){
				sum[s]+=tmpconc[s]*lineTOF[k][p];
			}
		}
		for (s=0;s<nspecies;s++){
			if (nLinesInCell[k]>0){Cnew[k][s]=sum[s]/totalTOF[k];}
			else                  {Cnew[k][s]=C[k][s];           }
		}
	}//end for (i=0..

}
//************************************************************************
//                           WRITE OUTPUT
//************************************************************************
void CStreamlineTrans::WriteOutput(const int outputstep){
	//no output 
}
	//################################################################################
	/*ofstream STREAMTMP; 
	STREAMTMP.open("cellbycell.txt");
	for (j=nY-1;j>=0;j--){
    for (i=0;i<nX; i++){STREAMTMP << nLinesInCell[i][j] <<" ";}STREAMTMP <<endl;
	}
	STREAMTMP <<endl;
	for (j=nY-1;j>=0;j--){
		//for (k=0;k<MAX_LINES_IN_CELL; k++){STREAMTMP << lineIDs[33][j][k] <<" ";}
    for (i=0;i<nX; i++){STREAMTMP << lineIDs[i][j][0] <<" ";}
		STREAMTMP <<endl;
	}
	STREAMTMP <<endl;
	for (j=nY-1;j>=0;j--){
    for (i=0;i<nX; i++){STREAMTMP << linePart[i][j][0] <<" ";}
		//for (k=0;k<MAX_LINES_IN_CELL; k++){STREAMTMP << linePart[33][j][k] <<" ";}		
		STREAMTMP <<endl;
	}
	STREAMTMP <<endl;

	STREAMTMP.close();*/
	//################################################################################
