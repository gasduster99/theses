//CStreamline.cpp
#include "Streamline.h"
#include "MatrixInclude.h"

int iabs(int i);
/********************************************************************
                     CONSTRUCTORS
*********************************************************************/
CStreamline::CStreamline(){}
//---------------------------------------------------------------------
CStreamline::CStreamline(const int            ID,
												 const CAquiferABC   *pAquifer,
												 const C2DDomainABC  *pDom,    //TMP DEBUG - not necc 2D
												 pt3D                 startpt, 
												 double               flux)
						:CPathline(ID,pAquifer,startpt,TRACK_FORWARD){

	pDomain  =pDom; 
	nspecies =pDomain->GetNumSpecies();
	pConcGrid=pDomain->GetGrid();
	if (pConcGrid->GetType()!=CELL_BASED){
		ExitGracefully("The Streamline method cannot be applied to non-cell-based grid or mesh",BAD_DATA);}

	UnevenCells =NULL;//All of these arrays initialized in CStreamline::Create1DGrid
	EvenCells   =NULL;
	nEvenCells  =0;
	nUnevenCells=0;
  maxtimestep =ALMOST_INF;
	qstart      =flux;
	captured_out_of_grid=false;
}
//---------------------------------------------------------------------
CStreamline::~CStreamline(){
	if (globaldebug){cout<<"DESTROYING STREAMLINE"<<endl;}
	int i;

	for (i=0; i<nEvenCells; i++){
		if (EvenCells  !=NULL){delete [] EvenCells[i].C; }
	}
	
	for (i=0; i<nUnevenCells; i++){
		if (UnevenCells!=NULL){delete [] UnevenCells[i].C; }
	}
	delete [] UnevenCells;
	delete [] EvenCells;
}
/*********************************************************************
                     STATIC MEMBERS/FUNCTIONS
**********************************************************************/
CTransect  **CStreamline::StartingFaces =new CTransect*[MAX_STREAMLINE_FACES];
//---------------------------------------------------------------------
int          CStreamline::nstartingfaces=0;
//---------------------------------------------------------------------
void CStreamline::AddStartingFace(CTransect *newface){
	if (nstartingfaces<MAX_STREAMLINE_FACES-1){
		StartingFaces[nstartingfaces]=newface;
		nstartingfaces++;
	}
	else {ExitGracefully("CStreamline::AddStartingFace-Too many streamline injection faces",TOO_MANY);}	
}
//---------------------------------------------------------------------
void CStreamline::Destroy(){//static Destructor
	delete [] StartingFaces;
}
/********************************************************************
                    CREATE STREAMLINES (called once)
*********************************************************************
 input:  Aquifer (pLayer), transport domain(pDom), time (t)
         NULL array of pointers to streamlines (lines) (NULL & unallocated)
				 desired total number of streamlines-just an estimate (nstreamlines)
 output: array of pointers to streamlines distributed along injection fronts (lines)
				 actual number of streamlines created (nstreamlines)
-------------------------------------------------------------------*/
void CStreamline::CreateStreamlines(const CSingleLayerABC *pLayer,
																		const C2DDomainABC		*pDom, 
																		CStreamline					 **&lines, 
																		int										 &nstreamlines, 
																		double								 t){
	int           f;
	double        TotalQ(0),fluxperline;

	//check if this function can operate-------------------------------
  //-----------------------------------------------------------------
	if (nstartingfaces==0){
		ExitGracefully("CStreamline::CreateStreamlines:No streamline injection faces specfied",BAD_DATA);}

	//calculate Total Influx, flux per streamline----------------------
  //-----------------------------------------------------------------
	TotalQ=0.0;

	double *frontFlux;
	frontFlux=new double [nstartingfaces];

	for (f=0; f<nstartingfaces; f++){
		StartingFaces[f]->CalculateStatistics(t);
    //StartingFaces[f]->WriteOutput(t); 
		frontFlux[f]=StartingFaces[f]->GetInflux();
		TotalQ+=frontFlux[f];
	}
  if (TotalQ<=0){
		ExitGracefully("CStreamline::CreateStreamlines:No influx from streamline injection faces",BAD_DATA);}
	
	fluxperline=TotalQ/(double)(nstreamlines);
	
	cout <<endl<<"************STREAMLINE CREATION****************"<<endl;
	cout <<"Total influx of faces:    "<<TotalQ      <<endl;	
	cout <<"Rough number of lines:    "<<nstreamlines<<endl;
	cout <<"Rough flux per streamline:"<<fluxperline <<endl; 

	//go through each face, creating streamlines-----------------------
  //-----------------------------------------------------------------

	        lines  =new CStreamline *[nstreamlines+4]; //TMP DEBUG- 4 should account for array overrun
  cmplex *zs     =new cmplex       [nstreamlines+4];
	double *fluxes =new double       [nstreamlines+4];
	int     nstream;
  int     s(0);

	for (f=0; f<nstartingfaces; f++){

    nstream      =(int)((max(frontFlux[f],0.0) * nstreamlines)/TotalQ);

		StartingFaces[f]->Subdivide(zs,fluxes,nstream,fluxperline);

	  for (int i=0; i<nstream; i++){

			cout <<"...creating streamline #"<<s<< "(front:"<< f <<", flux: "<<fluxes[i] <<")"<<endl;

			ExitGracefully("CStreamline::CreateStreamlines:MUST FIX THIS ISSUE...",RUNTIME_ERR);
			//lines[s]=new CStreamline(s,pLayer,pDom,c2Dto3D(zs[i]),fluxes[i]);//3DREWRITE
			
			if (lines[s]==NULL)   {ExitGracefully("CStreamline::Create: out of memory",OUT_OF_MEMORY);}
			if (s>nstreamlines-1) {ExitGracefully("CStreamline::Create: too many streamlines",RUNTIME_ERR);} 

			s++;
		}
	}
  nstreamlines=s;
	
	delete [] zs;
	delete [] fluxes;
	delete [] frontFlux;
}	

/*********************************************************************
                     ACCESSOR FUNCTIONS
**********************************************************************/
/*********************************************************************
                     GET 2D INFORMATION
*********************************************************************
input:
		index of cell in 2D grid (i2D)
returns:
		ntimesincell         (number of times the streamline entered the cell(visits))
		pieces[ntimesincell] (indices of streamline segment for each visit)
		TOFi  [ntimesincell] (amount time spent in cell per visit)
--------------------------------------------------------------------*/
void CStreamline::Get2Dinfo(const int        i2D, 
														int             &ntimesincell,
														int* const       pieces, 
														Writeable1DArray TOF) const{
	ntimesincell=0;
	for (int i=0; i<nUnevenCells; i++){
		if (ntimesincell>MAX_CROSSES_OF_CELL){
			ExitGracefully("CStreamline::Get2DInfo- streamline crossed cell too many times",TOO_MANY);
		}
		if (UnevenCells[i].ref_ind==i2D){
			pieces[ntimesincell]=i;
			TOF   [ntimesincell]=UnevenCells[i].delTOF;
			ntimesincell++;
		}
	}
}
/********************************************************************
                     GET CONCENTRATIONS
*********************************************************************
 input: Streamline Uneven cell index (i)
output: array of concentrations associated with cell i (tmpconc)
-------------------------------------------------------------------*/
void CStreamline::GetConcentrations(const int i, Writeable1DArray tmpconc) const{
	for (int s=0; s<nspecies; s++){
		tmpconc[s]=UnevenCells[i].C[s];
	}
}
//---------------------------------------------------------------------
double CStreamline::GetMaxTimeStep() const{return maxtimestep;}

/*********************************************************************
                     MODIFIER FUNCTIONS
*********************************************************************/

void CStreamline::SetFlux(const double delQ, const int piece){
	if (piece+1>=nUnevenCells){
		return; //TMP DEBUG 
		ExitGracefully("CStreamline::SetFlux: bad piece specified",RUNTIME_ERR);}
	UnevenCells[piece+1].q=delQ;
}
/*********************************************************************
                     CleanFluxes
*********************************************************************
converts from delta flux to cumulative flux along streamline
--------------------------------------------------------------------*/
void CStreamline::CleanFluxes(){
	for (int s=1; s<nUnevenCells; s++){
		UnevenCells[s].q+=UnevenCells[s-1].q;
	}
}
/*********************************************************************
                     INHERITED PARTICLE FUNCTIONS
**********************************************************************/
bool CStreamline::HasBeenCaptured(){
	captured_out_of_grid=false;
	if (CPathline::HasBeenCaptured()){return true;}
	else if (!pConcGrid->IsInside(lastnode->pt)){ 	//if no longer in grid bounds
		captured_out_of_grid=true;
		return true;
	}
	return false;
}
//--------------------------------------------------------------------
void CStreamline::Capture(const double &endtime){
	//cout << "Captured streamline out of grid @ " << lastnode->pt.x << " "<< lastnode->pt.y<<endl;
	if (!captured){
		if (!captured_out_of_grid){
			CPathline::Capture(endtime);
		}
		else{
			//temporary variables for use of pConcGrid->GetGridBoundaryIntercepts
			pt3D ptint[MAX_CELLS_JUMPED];
			int	 itmp [MAX_CELLS_JUMPED+1];
			int  numjumps(MAX_CELLS_JUMPED);

			pathnode *in,*out;
	
			in=lastnode->last;
			out=lastnode;	
			if ((in==NULL) || (out==NULL)){
				ExitGracefully("CStreamline::CheckIfCaptured: NULL pathnode",RUNTIME_ERR);}
			
			//TMP DEBUG: must handle new revisions to GetMeshBoundary intercepts
			ExitGracefully("CStreamline:Capture: James: must revise code",RUNTIME_ERR);
			pConcGrid->GetMeshBoundaryIntercepts(in->pt,out->pt,ptint,itmp,numjumps); 


			if (numjumps!=0){ //otherwise, streamline started outside grid
				lastnode->t=LinInterp(in->t,out->t,abs(in->pt),abs(out->pt),abs(ptint[numjumps-1]));//interpolate time & distance when it leaves grid
				lastnode->pt=ptint[numjumps-1];
				//cout << "CStreamline::CheckIfCaptured: " <<  ptint[numjumps-1].x <<" "<<  ptint[numjumps-1].y <<" "<< numjumps <<endl;
			}

			lastnode->next=new pathnode;    //create new node

			lastnode->next->pt=lastnode->pt;//initialize new node
			lastnode->next->s =lastnode->s;
			lastnode->next->t=endtime;     
			lastnode->next->last=lastnode;
			lastnode=lastnode->next;
			lastnode->next=PATHEND;
			numsteps++;
			captured=true;
		}
	}

}
/********************************************************************
       CREATE 1-DIMENSIONAL GRID(s) 
*********************************************************************
(called once per streamline)
(1) follows streamline through 2-D/3-D grid, 
			identifies 2D/3D indices and TOF points of intersections
(2) creates 1-D Uneven grid by dividing into sections if each 2D/3D cell
(3) creates 1-D Even grid
--------------------------------------------------------------------*/
void CStreamline::Create1DGrid(){


	pathnode      *node,*a,*b;									 //next pathline node, last node and current node 

	int           *i2D;                          //indices of ALL cells crossed thus far
	double	      *tmpTOF;                       //TOF at ALL cell crossings thus far
  int           tmpi2D[MAX_CELLS_JUMPED+1];    //indices of cells jumped in single particle movement
  pt3D          tmppt [MAX_CELLS_JUMPED];      //location of intercepts withh cell bounds in a single particle movement

	int            cellscrossed;                 //total number of cells crossed thus far
	int            numjumps;	  		             //total number of cells jumped in single particle movement 
	int			       step;                         //number of particle movements thus far

	int            i,s;                          //counter variables
	bool		       endingrid(true);              //true if streamline ends at boundary internal to grid

	int            expected_max(3*(int)(sqrt((double)(pConcGrid->GetNumCells())))); //expected maximum # of entered cells

	tmpTOF			=new double [expected_max+2];		 //initialize arrays
	i2D				  =new int    [expected_max+1]; 

  node        =firstnode;                      //initialize variables
	a           =NULL;
	b           =NULL;
	cellscrossed=0;
	step        =0;
	tmpTOF[0]   =0;
	numjumps    =0;

	if (pConcGrid->GetCellIndex(firstnode->pt,i2D[cellscrossed])){//starting point -> cellscrossed=0 
    
		//follow streamline forward in straight line segments
		//find the TOFs & cell indices where it crosses the concentration grid
		//This is the basis for the 1D uneven concentration Grid
		//------------------------------------------------------------------------------
		while ((node->next!=PATHEND) && 
					 (step<numsteps      ) && 
					 (endingrid          )){
		//for (node=firstnode;node->next!=PATHEND; node=node->next){
			node=node->next;
			a   =node->last;
			b   =node;
			step++;
			
			numjumps =MAX_CELLS_JUMPED;

			//TMP DEBUG: must handle new revisions to GetMeshBoundary intercepts
			ExitGracefully("CStreamline:Create1DGrid: James: must revise code",RUNTIME_ERR);	
			
			endingrid=pConcGrid->GetMeshBoundaryIntercepts(a->pt,b->pt,tmppt,tmpi2D,numjumps);


			if ((cellscrossed+numjumps)>(expected_max-1)){
				ExitGracefully("CStreamline:Create1DGrid: exceeded expected number of grid crossings(a)",RUNTIME_ERR);
			}
      
			for (i=0; i<numjumps; i++){
				if (tmpi2D[i+1]!=OFF){ //accounts for last index if jumps out of grid (only if !endingrid)
					i2D   [cellscrossed+i+1]=tmpi2D[i+1]; //+1 accounts for the fact that the starting ID from the segment has been accounted for in the preivous loop
				}
				tmpTOF[cellscrossed+i+1]=LinInterp(a->t,b->t,abs(a->pt),abs(b->pt),abs(tmppt[i])); 
				//cout << tmpTOF[cellscrossed+i+1]<<" "<<a->t<<" "<<b->t<<" "<<a->z.real()<<" "<<b->z.real()<<" "<<" "<<tmpz[i].real()<<endl;
			}
      cellscrossed+=numjumps;
		
		}//end while-path has been traversed 

		//add final tof if streamline ends within grid----------------------------------
    if (endingrid){
			tmpTOF[cellscrossed+1]=b->last->t;//real final time includes entire duration of capture
      //pConcGrid->GetCellIndex(b->pt,i2D[cellscrossed+1]);
			//if (i2D[cellscrossed+1]==OFF){ExitGracefully("CStreamline:Create1DGrid: should never happen",RUNTIME_ERR);}
			//cout << "endingrid "<<b->last->t<<endl;			
			cellscrossed++;
		}
		
		//Create 1D uneven concentration TOF grid
		//------------------------------------------------------------------------------
		nUnevenCells =cellscrossed;

		if ((tmpTOF[cellscrossed]<0) ||(tmpTOF[cellscrossed]>=ALMOST_INF)){
			ExitGracefully("CStreamline:Create1DGrid: poorly finished streamline",RUNTIME_ERR);}

		UnevenCells  =new strcell [nUnevenCells];
		if (UnevenCells==NULL){
			ExitGracefully("CStreamline:Create1DGrid: out of memory",OUT_OF_MEMORY);}

		for (i=0; i<nUnevenCells; i++){
			UnevenCells[i].TOFa   =tmpTOF[i];
			UnevenCells[i].TOFb   =tmpTOF[i+1];
      UnevenCells[i].delTOF =tmpTOF[i+1]-tmpTOF[i];
			UnevenCells[i].v      =abs(GetLocation(UnevenCells[i].TOFb)-GetLocation(UnevenCells[i].TOFa))/(UnevenCells[i].delTOF); 
			UnevenCells[i].pt     =GetLocation(0.5*(UnevenCells[i].TOFb + UnevenCells[i].TOFa));
			UnevenCells[i].ref_ind=i2D   [i];
			UnevenCells[i].index  =i;
			UnevenCells[i].q      =0.0;//qstart; //TMP DEBUG-should obtain from 2D grid
			UnevenCells[i].C =new double [nspecies];

			if (UnevenCells[i].C==NULL){ExitGracefully("CStreamline:Create1DGrid: out of memory",OUT_OF_MEMORY);}
	
			for (s=0; s<nspecies; s++){
				UnevenCells[i].C[s]   =
					pDomain->GetConcentration(UnevenCells[i].pt,0.0,s);   
					//pDomain->GetConcentration(UnevenCells[i].ref_ind,0.0,s);
			}

		}

		//Create 1D Even TOF grid
		//------------------------------------------------------------------------------
		//should calculate maximum grid spacing: maxspace - peclet number, adjacent cell sizes

		double spacing, currTOF;

		//identify size of evenly spaced grid,
		nEvenCells  =2*nUnevenCells;//arbitrary. This is what Batycky(1999) used		
		spacing    =UnevenCells[nUnevenCells-1].TOFb/nEvenCells; 
		maxtimestep=0.9*spacing;

		if ((spacing<=0.0) || (spacing>1e6)){
			//cout <<UnevenCells[nUnevenCells-1].TOFb<<" "<<spacing<<" "<< nEvenCells<<" "<<nUnevenCells<<endl;
			ExitGracefully("CStreamline:Create1DGrid: negative or huge spacing",RUNTIME_ERR);}

    //allocate memory for grid
		EvenCells=new strcell[nEvenCells];
		if (EvenCells==NULL){
			ExitGracefully("CStreamline:Create1DGrid: out of memory",OUT_OF_MEMORY);}

		//initialize grid cell information		
		currTOF=0.0;
		for (i=0; i<nEvenCells; i++){
      EvenCells[i].delTOF =spacing;
			EvenCells[i].TOFa   =currTOF;
			EvenCells[i].TOFb   =currTOF+spacing;
			currTOF=currTOF+spacing;
			EvenCells[i].pt     =GetLocation(0.5*(EvenCells[i].TOFb + EvenCells[i].TOFa));
			EvenCells[i].v=abs(GetLocation(EvenCells[i].TOFb)-GetLocation(EvenCells[i].TOFa))/(EvenCells[i].delTOF); //velocity over TOF step
			
			EvenCells[i].ref_ind=0;//unused for now
			EvenCells[i].index  =i;
			EvenCells[i].q      =0;//unused for now UnevenCells[i].q;
			EvenCells[i].C =new double [nspecies];
			if (EvenCells[i].C==NULL){ExitGracefully("CStreamline:Create1DGrid: out of memory",OUT_OF_MEMORY);}
			for (s=0; s<nspecies; s++){
				EvenCells[i].C[s]      =0.0;
			}
		}

		TranslateUnevenToEven();
	}
	else {cout <<"Streamline started outside of grid!" << node->pt.x<<" "<<node->pt.y<<endl;}

	delete [] tmpTOF;     
	delete [] i2D; 
}

/*********************************************************************
		UPDATE CONCENTRATIONS 
----------------------------------------------------------------------
Performs advection and dispersion on 1D TOF grid 
*********************************************************************/
void CStreamline::UpdateConcentrations(const double          &tstep, 
																			 const double          &t,
																       const str_update_type  upty){

	int    i;

	if (tstep<=0.0){	
				ExitGracefully("CStreamline:UpdateConcentrations: bad timestep",RUNTIME_ERR);}

	if      (upty==UPDATE_FROM_UNEVEN){ 
		/*only done when other transformations performed on 2D grid (i.e. reaction)
		  get concentrations from 2D/3D grid, translate to Uneven 1D grid then to Even 1D grid
		-----------------------------------------------------------*/
		for (i=0; i<nUnevenCells;i++){
			pDomain->GetConcentrations(UnevenCells[i].pt,t,UnevenCells[i].C);  
		//pDomain->GetConcentrations(UnevenCells[i].ref_ind ,t,UnevenCells[i].C);//Keep here-may need later  
		}
		TranslateUnevenToEven();
	}
	else if (upty==UPDATE_DIRECT){
		for (i=0; i<nEvenCells;i++){
			pDomain->GetConcentrations(EvenCells[i].pt,t,EvenCells[i].C);    
		}
		TranslateEvenToUneven();
	}
	else if (upty==NO_UPDATE){
	}
	
	UpdateImplicit(tstep,t); 

	//Translate from even back to uneven grid
  //--------------------------------------------------
	TranslateEvenToUneven();

}
/*********************************************************************
		UPDATE IMPLICIT 
----------------------------------------------------------------------
	performs advection and dispersion on 1D TOF grid using FD implicit formulation
	ADR in terms of TOF:
	dC	=	-dC	+	a_s*d^2C	+ Nr*(Cs-C)
	dt		 dT		v_s	dT^2		 n*H
*********************************************************************/
void CStreamline::UpdateImplicit(const double &tstep, 
																 const double &t){

	double *e,*f,*g,*b,*sol;
	double *lterm;
	double  Dterm,vterm;
	double  del1,del2;
	double  local_tstep,last_local_tstep,local_t;
	double  poro, satthick;
	int     i,s;
	pt3D    thispt;
	double  alpha,vel,Nr(0);

	if (tstep<=0.0){	
		ExitGracefully("CStreamline:UpdateImplicit: bad timestep",RUNTIME_ERR);}
	
	local_t=0.0;
	local_tstep=min(maxtimestep,tstep-t);
	last_local_tstep=local_tstep;

	lterm =new double [nEvenCells];
	e			=new double [nEvenCells];
	f			=new double [nEvenCells];
	g			=new double [nEvenCells];
	b			=new double [nEvenCells];
	sol		=new double [nEvenCells];
	if (sol==NULL){
		ExitGracefully("CStreamline:UpdateImplicit: out of memory",OUT_OF_MEMORY);}

	//Build Matrix (reusable for this time step)
	//------------------------------------------------------------------------
	for (i=0; i<nEvenCells; i++){
		 
		thispt  =GetLocation((EvenCells[i].TOFa+EvenCells[i].TOFb)/2.0);        //center of TOF grid
		vel     =EvenCells[i].v; 

		alpha   =pDomain->GetDispersivity  (thispt,LONGITUDINAL);
		poro    =pAq->GetPoro              (thispt);
		//satthick=pAq->GetSaturatedThickness(thispt,t);
		//
		satthick=1.0;//pLayer->GetSaturatedThickness(c3Dto2D(thispt),t);//TMP DEBUG
		ExitGracefully("CStreamline::UpdateImplicit: work needed",BAD_DATA);
		Nr      =pAq->GetLeakage           (thispt,t,FROMTOP);

		if (i>0)           {del1=EvenCells[i  ].delTOF+EvenCells[i-1].delTOF;}
		else               {del1=EvenCells[i  ].delTOF;                      }
		if (i<nEvenCells-1){del2=EvenCells[i+1].delTOF+EvenCells[i  ].delTOF;}
		else               {del2=del1;                                       }

		if ((vel<=0.0) || (poro<=0.0)){	
				cout <<endl<< vel<<" "<<poro<<endl;
				ExitGracefully("CStreamline:UpdateImplicit: bad velocity or porosity",RUNTIME_ERR);}
		if ((del1<=0.0) || (del2<=0.0)){	
				cout <<endl<< del1<<" "<<del2<<" "<<EvenCells[i  ].delTOF<<endl;
				ExitGracefully("CStreamline:UpdateImplicit: bad delta TOF",RUNTIME_ERR);}

		lterm[i]=Nr/(poro*satthick);
		if (satthick==0) {lterm[i]=0.0;}
		Dterm=8.0*alpha/((del1+del2)*vel);

		vterm=2.0/del1;

		if      (i==0)           {
			e[i]=(0.0             )*local_tstep;
			g[i]=(Dterm/del2      )*local_tstep;
		}
		else if (i==nEvenCells-1){
			e[i]=(           vterm)*local_tstep;
			g[i]=(0.0             )*local_tstep;
		}
		else                     {
			e[i]=(Dterm/del1+vterm)*local_tstep;
			g[i]=(Dterm/del2      )*local_tstep;
		}				
		f[i] =(-Dterm*(1/del1+1/del2)-vterm+lterm[i]);
		f[i]*=last_local_tstep;
		f[i]-=1.0;
	}


	while (local_t<tstep){
		//solve differential equation on 1D even grid
		//------------------------------------------------------------

		for (i=0; i<nEvenCells; i++){
			//TMP DEBUG- should account for source terms
			//thispt=GetLocation(0.5*(EvenCells[i].TOFa+EvenCells[i].TOFb));  
		  //pDomain->GetSourceConcentration(thispt,t+local_t,s); 

			e[i]/=last_local_tstep;
			g[i]/=last_local_tstep;		
			f[i]+=1.0;
			f[i]/=last_local_tstep;

			f[i]*=local_tstep;
			f[i]-=1.0;
			g[i]*=local_tstep;
			e[i]*=local_tstep;
		}

		//Thomas Algorithm (done for each species):
		//------------------------------------------
		for (s=0; s<nspecies; s++){

			//Recompute b vector for each species:
			for (i=0; i<nEvenCells; i++){
				b[i]=-EvenCells[i].C[s]-lterm[i]*local_tstep*pDomain->GetRechargeConcentration(EvenCells[i].pt,t+local_t,s);  
				//cout << e[i] << " "<<f[i]<<" "<<g[i]<<"|"<<b[i]<<endl;
				//cout << e[i] << " "<<f[i]<<" "<<g[i]<<"|"<<b[i]<<" , "<<local_tstep<<","<<del1<<","<<del2<<","<<maxtimestep<<endl;
				if ((b[i]==exp(2004)) || (b[i]==-exp(2004)) || NotANumber(b[i])){
					ExitGracefully("CStreamline::UpdateConcentrations: poor matrix: singular results",SINGMAT);}
				if ((b[i]>ALMOST_INF) || (b[i]<-ALMOST_INF)){
					ExitGracefully("CStreamline::UpdateConcentrations: poor matrix: singular results",SINGMAT);}
			}
			//cout <<endl;

			//Solve tridiagonal matrix
			ThomasAlgorithm(e,f,g,b,sol,nEvenCells);

			//Copy solution
			for (i=0; i<nEvenCells-1; i++){
				EvenCells[i].C[s]=max(sol[i],0.0);//TMP DEBUG: bastardized version to avoid negative flux
				if ((sol[i]>ALMOST_INF) || (sol[i]<-ALMOST_INF) || NotANumber(sol[i])){
					ExitGracefully("CStreamline::UpdateConcentrations: poor matrix: singular results",SINGMAT);}
			}
		}

		last_local_tstep =local_tstep;
		local_tstep      =min(maxtimestep,tstep-local_t);
		local_t         +=local_tstep; 

	}//end while t<timestep

	delete [] lterm;
	delete [] e ;
	delete [] f ;
	delete [] g ;
	delete [] b ;
	delete [] sol;

}
//*********************************************************************
void CStreamline::TranslateEvenToUneven(){
	//use TOF averaging 
	int i,j,s;
	int jlow(0),jupp(0);
	for (i=0; i<nUnevenCells; i++){
		while ((EvenCells[jlow].TOFa<=UnevenCells[i].TOFa) && (jlow<nEvenCells-1)){jlow++;}//identify range of uneven cells overlapping even cell i
		while ((EvenCells[jlow].TOFa> UnevenCells[i].TOFa) && (jlow>0           )){jlow--;}

		while ((EvenCells[jupp].TOFb>=UnevenCells[i].TOFb) && (jupp>0           )){jupp--;}
		while ((EvenCells[jupp].TOFb< UnevenCells[i].TOFb) && (jupp<nEvenCells-1)){jupp++;}

		if (jlow==jupp){ //only one cell overlaps cell i (complete containment-direct assignment)
			for (s=0; s<nspecies; s++){UnevenCells[i].C[s]=EvenCells[jlow].C[s];}
		}
		else{            //multiple uneven cells correspond to single even cell (requires TOF averaging)
			for (s=0; s<nspecies; s++){
				UnevenCells[i].C[s]=0.0;
			  UnevenCells[i].C[s]+=EvenCells[jlow].C[s]*(EvenCells[jlow].TOFb-UnevenCells[i].TOFa);
				UnevenCells[i].C[s]+=EvenCells[jupp].C[s]*(UnevenCells[i].TOFb-EvenCells[jupp].TOFa);
				for (j=jlow+1;j<jupp;j++){
					UnevenCells[i].C[s]+=EvenCells[j].C[s]*EvenCells[j].delTOF;
				}
				UnevenCells[i].C[s]/=UnevenCells[i].delTOF;
			}
		}
	}
}
//*********************************************************************
void CStreamline::TranslateUnevenToEven(){
	int i,j,s;
	int jlow(0),jupp(0);
	for (i=0; i<nEvenCells; i++){
		while ((UnevenCells[jlow].TOFa<=EvenCells[i].TOFa) && (jlow<nUnevenCells-1)){jlow++;}//identify range of uneven cells overlapping even cell i
		while ((UnevenCells[jlow].TOFa> EvenCells[i].TOFa) && (jlow>0             )){jlow--;}
		
		while ((UnevenCells[jupp].TOFb>=EvenCells[i].TOFb) && (jupp>0             )){jupp--;}		
		while ((UnevenCells[jupp].TOFb< EvenCells[i].TOFb) && (jupp<nUnevenCells-1)){jupp++;}

		if (jlow==jupp){ //only one cell overlaps cell i (complete containment-direct assignment)
			for (s=0; s<nspecies; s++){EvenCells[i].C[s]=UnevenCells[jlow].C[s];}
		}
		else{            //multiple uneven cells correspond to single even cell (requires TOF averaging)
			for (s=0; s<nspecies; s++){
				EvenCells[i].C[s]=0.0;
			  EvenCells[i].C[s]+=UnevenCells[jlow].C[s]*(UnevenCells[jlow].TOFb-EvenCells[i].TOFa);
				EvenCells[i].C[s]+=UnevenCells[jupp].C[s]*(EvenCells[i].TOFb-UnevenCells[jupp].TOFa);
				for (j=jlow+1;j<jupp;j++){
					EvenCells[i].C[s]+=UnevenCells[j].C[s]*UnevenCells[j].delTOF;
				}
				EvenCells[i].C[s]/=EvenCells[i].delTOF;
			}
		}
	}
}
//*********************************************************************
void CStreamline::WriteOutput(){
	pt3D pt,pt2,ptcen;
  int    i;
	ofstream STREAMOUT;
	STREAMOUT.open("streamlines_on_2D_grid.csv",ios::app);

	for (i=0; i<nUnevenCells; i++){
		pt=GetLocation(UnevenCells[i].TOFb);
		pt2=GetLocation((UnevenCells[i].TOFa+UnevenCells[i].TOFb)/2.0);
		ptcen=pConcGrid->GetNodeLocation(UnevenCells[i].ref_ind); 
		STREAMOUT<< i											<<","
			       << (UnevenCells[i].TOFa+UnevenCells[i].TOFb)/2.0 <<","
						 << pt.x              <<","
						 << pt.y              <<","
						 << pt2.x             <<","
						 << pt2.y             <<","
						 << UnevenCells[i].ref_ind<<" ,"
						 << ptcen.x               << " ,"
						 << ptcen.y               << " ,"
						 << UnevenCells[i].q      << " ,";
		if (EvenCells[i].C!=NULL){STREAMOUT<< UnevenCells[i].C[0];	}  
		else                     {STREAMOUT<<"NULL";                }
		STREAMOUT<<endl;
	}
	STREAMOUT <<endl;
	STREAMOUT.close();

	STREAMOUT.open("streamline 1D grid.csv",ios::app);
	for (i=0; i<nUnevenCells; i++){
		STREAMOUT<< "U,";
		STREAMOUT<< UnevenCells[i].index <<",";
		STREAMOUT<< UnevenCells[i].TOFa <<",";
		STREAMOUT<< UnevenCells[i].TOFb <<",";
		STREAMOUT<< UnevenCells[i].v    <<",";
		STREAMOUT<< ((UnevenCells[i].TOFb+ UnevenCells[i].TOFa)*0.5) <<",";
		if (UnevenCells[i].C!=NULL){ STREAMOUT<< UnevenCells[i].C[0] <<",";}
		else                       { STREAMOUT<< "NULL,";                  }
		STREAMOUT<< UnevenCells[i].q <<endl;
	}
	for (i=0; i<nEvenCells; i++){
		STREAMOUT<< "E,";
		STREAMOUT<< EvenCells[i].index <<",";
		STREAMOUT<< EvenCells[i].TOFa <<",";
		STREAMOUT<< EvenCells[i].TOFb <<",";
		STREAMOUT<< EvenCells[i].v    <<",";
		STREAMOUT<< ((EvenCells[i].TOFb+ EvenCells[i].TOFa)*0.5) <<",";
		if (EvenCells[i].C!=NULL){ STREAMOUT<< EvenCells[i].C[0] <<",";}
		else                     { STREAMOUT<< "NULL,";                }
		STREAMOUT<< EvenCells[i].q <<endl;
	}
	STREAMOUT.close();
}
