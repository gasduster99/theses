//EulerianDispersion_Rect.cpp
#ifdef JHUNK 

#include "FDEulerian.h"
//***************************************************************************
C1DFDEulerian::C1DFDEulerian(C2DDomainABC *pDom, 
														 CRectGrid    *pConcGrid, 
														 const double  time_weight, 
														 const double  up_weight)
              :C1DTransportScheme(pDom){

	if ((time_weight>1.0) || (time_weight<0.0)){
		ExitGracefully("C1DFDEulerian: Invalid value for temporal weighting term (<0 or >1)",BAD_DATA);}
	if ((up_weight>1.0)   || (up_weight<0.0)  ){
		ExitGracefully("C1DFDEulerian: Invalid value for upstream weighting term (<0 or >1)",BAD_DATA);}

	pGrid =pConcGrid;
	nNodes=nCells=pConcGrid->GetNumNodes(); //Block centered grid
	
	pLayer=pDomain->GetLayer();

	sinkflux    =NULL;
	sourceflux  =NULL;
  rechargeflux=NULL;
	dirichlet   =NULL;
	
	poro  =NULL;
	thick =NULL;

	AW    =NULL;
	A     =NULL;
	B     =NULL;
	Cs    =NULL;
	dC    =NULL;
	
	w=time_weight; //temporal weighing term- crank nicholson by default
	               //=1.0 for fully implicit
	               //=0 for fully explicit
	upwt=up_weight; //upstream weighting term- 
	                //1.0 for fully upstream
                  //zero for central difference
	print_diagnostics=false;
}
//***************************************************************************
C1DFDEulerian::~C1DFDEulerian(){
	if (globaldebug){cout <<"DESTROYING RECTANGULAR DISPERSION SCHEME "<<endl;}

	delete [] sinkflux;
	delete [] sourceflux;
	delete [] rechargeflux;
	delete [] dirichlet;

	//delete matrices
	delete A;
	delete AW;
	delete [] B;

	delete [] dC;
	delete [] Cs;

}
/***************************************************************************
							PRINT GRID DIAGNOSTICS
****************************************************************************
---------------------------------------------------------------------------*/
void C1DFDEulerian::PrintGridDiagnostics(){print_diagnostics=true;}
/***************************************************************************
							INITIALIZE
****************************************************************************
Creates Vertically Averaged Dispersivity Matrices- dispersivity at cell faces
---------------------------------------------------------------------------*/
void C1DFDEulerian::Initialize(const double         startt,
															 const transport_type ty,
									             Ironclad1DArray      porosity, 
													     Ironclad1DArray      satthick){

	ExitGracefullyIf((pDomain==NULL),"C1DFDEulerian::Initialize: NULL Domain",RUNTIME_ERR);
	ExitGracefullyIf((pGrid  ==NULL),"C1DFDEulerian::Initialize: NULL Grid  ",RUNTIME_ERR);
	
	//rectangular grid parameters
	int			nX    =pGrid->nX;
	int     nNodes=pGrid->nNodes; 
	double *nodeX =pGrid->nodeX;
	double *nodeY =pGrid->nodeY;
	double *gridX =pGrid->gridX;
	double *gridY =pGrid->gridY;
	double PeMax(0.0),PeAvg(0.0);
	double CrMax(0.0),CrAvg(0.0);
	double vloc, dloc,loc_len;

  nspecies=pDomain->GetNumSpecies();

	//local function parameters
	int     i,j,k,i1,j1,ktmp;         //counters
	double  t;                       //time,     
	cmplex  tmpz;                     
	double  dx[3],dy[3],locA[3][3];   
	double  tmpA;
	cmplex  pts[5];
	double  sinkterm,tmpin,tmpout;
	cmplex  z1,z2,z3,z4;

	double  dxavg(0),dyavg(0);
	double  avgarea(0);
	double  upweigh;
  double *peclet;
	double  max_aspect_ratio(0.0);
	cout << "Initializing Eulerian Advection/Dispersion Scheme..."<<endl;

	if (ProgramAborted()){return;}

	process_type  =ty;
	poro          =porosity;
  thick         =satthick;
	t             =startt+0.001;//So as to neglect initial conditions
	sugg_tstep    =ALMOST_INF;
	last_loc_tstep=1.0;


	CVector x_orientation(cos(pGrid->orient),sin(pGrid->orient),0.0);


	//Determine which cells are active---------------------------------
  cout <<"  Reserving Dynamic Memory/Identifying Active Cells..."<<endl;	
	dC          =new double [nCells];//explicit component of change in concentration
	Cs          =new double [nCells];//concentration of current species
	active      =new bool   [nCells];
	dirichlet   =new bool   [nCells];
	sinkflux    =new double [nCells];
	sourceflux  =new double [nCells];
	borderoutflx=new double [nCells];
  rechargeflux=new double [nCells];
	peclet      =new double [nCells];
	ExitGracefullyIf((rechargeflux==NULL),"C1DFDEulerian::Initialize(1a)",OUT_OF_MEMORY);
	for (k=0; k<nCells;k++){
		active   [k]=(satthick[k]>0.0);
		sinkflux [k]=sourceflux[k]=rechargeflux[k]=borderoutflx[k]=0.0; 
		dirichlet[k]=false;
	}

	//identify sink & source fluxes------------------------------------
  cout <<"  Initializing Source/Sink Information...";	
	for (k=0; k<nCells;k++){
		WriteEllipsisToScreen(k,nCells,20);

		pGrid->k_to_ij(k,i,j);
		z1=pGrid->LocalToGlobal(cmplex(gridX[i+1],gridY[j  ]));
		z2=pGrid->LocalToGlobal(cmplex(gridX[i  ],gridY[j  ]));
		z3=pGrid->LocalToGlobal(cmplex(gridX[i+1],gridY[j+1]));
		z4=pGrid->LocalToGlobal(cmplex(gridX[i  ],gridY[j+1]));

		//recharge flux is just a single component of source and sink fluxes
		rechargeflux[k]-=pLayer->GetIntegratedLeakage(z1,z2,z3,t,FROMTOPANDBOTTOM);//TMP DEBUG-should differentiate between bottom and top recharge
		rechargeflux[k]-=pLayer->GetIntegratedLeakage(z3,z2,z4,t,FROMTOPANDBOTTOM);//TMP DEBUG-should differentiate between bottom and top recharge

		pLayer->GetIntegratedBudget (z1,z2,z3,t,tmpin,tmpout);
		sinkflux    [k]+=tmpout;
		sourceflux  [k]+=fabs(tmpin); //should be unnecc

		pLayer->GetIntegratedBudget (z3,z2,z4,t,tmpin,tmpout);
		sinkflux    [k]+=tmpout;
		sourceflux  [k]+=fabs(tmpin); //should be unnecc

	}

	InitializeBorderFlux(t);

	//Alternative Explicit/Implicit/CN Format---------------------------
	//[[A]-[I]/dt+Qsink]*{Cn+1}=[AW]{Cn}-{Cn}/dt-{Qs}^T{Cs}
	A =new CSparseMatrix();
	AW=new CSparseMatrix();
	A ->DynamicInitialize(nCells);
	AW->DynamicInitialize(nCells);

	B=new double [nCells];
	for (k=0;k<nCells;k++){B[k]=0.0;}

  cout <<endl<<"  Initializing Matrices";
	//tmp debug
	//ofstream MB;MB.open("Flow Mass Balance.csv");MB<<"k,x,y,sidein,source,sink,net"<<endl;MB.close();

	
	for (j=0; j<nY; j++){

		if (ProgramAborted()){return;}	

		for (i=0;i<nX;i++){
		
			k =pGrid->ij_to_k(i,j);

			WriteEllipsisToScreen(k,nCells,20);

			pts[4]=c3Dto2D(pGrid->GetNodeLocation(k));

			if (IsActive(k)){

				//get local grid spacing,activity on 9-pt stencil
				for (i1=0; i1<=2; i1++){                 
					if      ((i==0   ) && (i1==0)){dx [i1]=gridX[1   ]-gridX[0     ];}//shouldn't matter
					else if ((i==nX-1) && (i1==2)){dx [i1]=gridX[nX  ]-gridX[nX-1  ];}//shouldn't matter
					else                          {dx [i1]=gridX[i+i1]-gridX[i+i1-1];} 
					if      ((j==0   ) && (i1==0)){dy [i1]=gridY[1   ]-gridY[0     ];}//shouldn't matter
					else if ((j==nY-1) && (i1==2)){dy [i1]=gridY[nY  ]-gridY[nY-1  ];}//shouldn't matter
					else                          {dy [i1]=gridY[j+i1]-gridY[j+i1-1];}
					ExitGracefullyIf((dy [i1]<0),"C1DFDEulerian::Initialize: bad dy",RUNTIME_ERR);
          ExitGracefullyIf((dx [i1]<0),"C1DFDEulerian::Initialize: bad dx",RUNTIME_ERR);
				}
				pts[0]=pGrid->LocalToGlobal(cmplex(nodeX[i  ],gridY[j+1])); //top
				pts[1]=pGrid->LocalToGlobal(cmplex(gridX[i+1],nodeY[j  ])); //right
				pts[2]=pGrid->LocalToGlobal(cmplex(nodeX[i  ],gridY[j  ])); //bottom
				pts[3]=pGrid->LocalToGlobal(cmplex(gridX[i  ],nodeY[j  ])); //left

				//cout <<"dx:"<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<endl;
				//cout <<"dy:"<<dy[0]<<" "<<dy[1]<<" "<<dy[2]<<endl;

				upweigh=upwt;
				//turn off upstream weighing near sinks & sources
				/*if ((sinkflux[k]!=0.0) && (sourceflux[k]!=0.0)){upweigh=0.0;}
				ktmp=pGrid->ij_to_k(i+1,j  ); if ((IsActive(ktmp)) && (sinkflux[ktmp]!=0.0) && (sourceflux[ktmp]!=0.0)){upweigh=0.0;}
				ktmp=pGrid->ij_to_k(i+1,j+1); if ((IsActive(ktmp)) && (sinkflux[ktmp]!=0.0) && (sourceflux[ktmp]!=0.0)){upweigh=0.0;}
				ktmp=pGrid->ij_to_k(i  ,j  ); if ((IsActive(ktmp)) && (sinkflux[ktmp]!=0.0) && (sourceflux[ktmp]!=0.0)){upweigh=0.0;}
				ktmp=pGrid->ij_to_k(i  ,j+1); if ((IsActive(ktmp)) && (sinkflux[ktmp]!=0.0) && (sourceflux[ktmp]!=0.0)){upweigh=0.0;}
        *///upweigh=0.0;

				//obtain influence of dispersion/advection, get local velocity and dispersion coefficients
				GetLocalMatrix(k,locA,dx,dy,pts,0,x_orientation,upweigh,vloc,dloc); //TMP DEBUG-MULTISPECIES

				avgarea+=dx[1]*dy[1];
				loc_len=sqrt(dx[1]*dy[1]);

				//assemble equation for matrix row k
				for (i1=0; i1<=2; i1++){ 

					for (j1=0; j1<=2; j1++){ 
					
						ktmp =pGrid->ij_to_k(i+i1-1,j+j1-1); //ktmp is sparse matrix column index of adjacent cell

						if (IsActive(ktmp)){ 

							tmpA=A ->GetAij(k,ktmp)+(  w)*locA[i1][j1]; //k is row index
							A ->DynamicAdd(tmpA,k,ktmp);                //includes temporal weighing
							
							tmpA=AW->GetAij(k,ktmp)+(1-w)*locA[i1][j1];
							AW->DynamicAdd(tmpA,k,ktmp);

              //cout << "k:"<<k<<" ktmp:"<<ktmp<<"A:"<<AW->GetAij(k,ktmp)<< "->"<<tmpA<<endl;
						}
					}
				}
				
				//Add in influence of time derivative term (independent of temporal weighing)
				A ->AddToAij(-(1.0/last_loc_tstep),k,k);
				
				//Add in influence of time-weighted sink term
				sinkterm=sinkflux[k]/poro[k]/satthick[k]/pGrid->GetCellArea(k); //TMP DEBUG
				A ->AddToAij(-(w    )*sinkterm,k,k);
				AW->AddToAij(-(1.0-w)*sinkterm,k,k);

				double bordterm;
				bordterm=max(borderoutflx[k],0.0)/poro[k]/satthick[k]/pGrid->GetCellArea(k)/2.0; //TMP DEBUG-2.0 is unknown factor- works for any w, any upstream factor, 

				A ->AddToAij(-(w    )*bordterm,k,k);
				AW->AddToAij(-(1.0-w)*bordterm,k,k);


				//Evaluate local Peclet Number 
				PeAvg+=loc_len*fabs(vloc)/dloc;
				upperswap(PeMax,loc_len*fabs(vloc)/dloc);
				peclet[k]=loc_len*fabs(vloc)/dloc;

				//Evaluate local Courant number (assume tstep=1)
				CrAvg+=fabs(vloc)/loc_len;
				upperswap(CrMax,fabs(vloc)/loc_len);

				upperswap(max_aspect_ratio,dx[1]/dy[1]);
				upperswap(max_aspect_ratio,dy[1]/dx[1]);
			}
			else{ //inactive cell
				A ->SetAij(-(1.0/last_loc_tstep),k,k);
				AW->SetAij(0.0,k,k);
			}
		}//for j<nY
	}//for i<nX
	
	avgarea/=nCells;
	CrAvg/=nCells;
	PeAvg/=nCells;

	//Identify dirichlet conditions------------------------------
	cout <<endl<<"Identifying Dirichlet conditions";
	int numDir=0;
	int *dir=new int [nCells];
	for (k=0; k<nCells; k++){
		WriteEllipsisToScreen(k,nCells,20);
		dirichlet[k]=(pDomain->GetDirichletConc(pGrid->GetNodeLocation(k),t,0)!=NO_INFLUENCE);
		if (dirichlet[k]){dir[numDir]=k;numDir++;}
	}
	A ->MakeDirichlet(dir,numDir);
	//AW->MakeDirichlet(dir,numDir); //doesn't seem to matter
	cout <<endl<<"...Number of Dirichlet Cells:"<< numDir<<endl;
	delete [] dir;
  //-----------------------------------------------------------

	//TMP DEBUG
	//A->Print(true,50);
	//cout<<"___________________________________"<<endl;
	//AW->Print(false,20);
	//for (i=0;i<nCells;i++){
	//	cout << AW->GetAij(i,i) <<endl;
	//}
	//A->ScalarMult(0.5);
	//for (i=0;i<nCells;i++){
		//cout << B[i] <<endl;
	//}
	//ExitGracefully("discontinuity print",BAD_DATA);


	cout <<"   -NUMBER OF DIRICHLET NODES:         " << numDir				<<endl;	 
	cout <<"   -AVERAGE CELL AREA:                 " << avgarea				<<endl;
	cout <<"   -AVERAGE REPRESENTATIVE LENGTH:     " << sqrt(avgarea)	<<endl;
  cout <<"   -MAXIMUM GRID PECLET NUMBER:        " << PeMax					<<endl;
  cout <<"   -AVERAGE GRID PECLET NUMBER:        " << PeAvg					<<endl;
  cout <<"   -SUGGESTED TIME STEP (Cr(max)=1):   " << ((process_type==ADV_AND_DISP) ? 1.0/CrMax : sugg_tstep)  <<endl;
  cout <<"   -PERMISSIBLE TIME STEP (Cr(avg)=1): " << ((process_type==ADV_AND_DISP) ? 1.0/CrAvg : sugg_tstep)  <<endl;
	
	//Apply Courant constraints (Cr<=1) 
	if (process_type==ADV_AND_DISP){lowerswap(sugg_tstep,1.0/CrMax);}
	//sugg_tstep=ALMOST_INF; //Overrides suggested time step	
  
	StraightSort(peclet,nCells);

  if (print_diagnostics){
		ofstream DIAG;
		DIAG.open("diagnostics.csv");
		DIAG<<"Number_Of_Cells,"          <<nCells				<<endl;
		DIAG<<"Number_Of_Dirichlet_Cells,"<<numDir				<<endl;
		DIAG<<"Average_Cell_Area,"        <<avgarea				<<endl;
    DIAG<<"Max_Aspect_Ratio,"         <<max_aspect_ratio<<endl;
		DIAG<<"Average_Represent_Length," <<sqrt(avgarea)	<<endl;
		DIAG<<"Max_Peclet_Number,"        <<PeMax					<<endl;
		DIAG<<"Mean_Peclet_Number,"       <<PeAvg					<<endl;
		DIAG<<"90%_Peclet_Number,"        <<peclet[(int)(0.90*nCells)]<<endl;
		DIAG<<"80%_Peclet_Number,"        <<peclet[(int)(0.80*nCells)]<<endl;
		DIAG<<"75%_Peclet_Number,"        <<peclet[(int)(0.75*nCells)]<<endl;
		DIAG<<"50%_Peclet_Number,"        <<peclet[(int)(0.50*nCells)]<<endl;
		DIAG<<"25%_Peclet_Number,"        <<peclet[(int)(0.25*nCells)]<<endl;
    DIAG<<"Courant_Time_Step,"			  <<1.0/CrMax     <<endl;
		DIAG.close();
	}

	delete [] peclet;
	//lowerswap(sugg_tstep,loc_sugg_tstep);
  ExitGracefullyIf((sugg_tstep<=0.0),"C1DFDEulerian::Initialize: negative or zero suggested maximum time step",RUNTIME_ERR);
	
	cout <<"...Done Initializing Finite Difference Transport Scheme..."<<endl;

	//cout <<"Translating to MT3D..."<<endl;
	//TranslateToMT3D();
	//cout <<"...MT3D Translation done."<<endl;
}
/***************************************************************************
							Transport
****************************************************************************
Eulerian Finite Difference Advection/Dispersion on rectangular grid- Implicit/Explicit Formulation
---------------------------------------------------------------------------*/
void C1DFDEulerian::Transport(const double     t, 
												      const double     tstep,
												      Ironclad2DArray  C,
												      Ironclad2DArray  Cend,
												      Writeable2DArray Cnew){
	//--local variables------------
	int		  k,s;
	double  local_t,loc_tstep;
	cmplex  tmpz;
	pt3D	  pt;
	double  tmp_err(0.0),area;
	double  aW(0.5); //advective weight=1 for using post-advection conc, 0 for using no-advection conc.
	double  sourceterm,mult;
	
	if (process_type==ADV_AND_DISP){mult=0.5;}else{mult=1.0;}

	ExitGracefullyIf((nCells<=0),   "C1DFDEulerian::Disperse: Dispersion Scheme Not Initialized(1)",RUNTIME_ERR);
	ExitGracefullyIf((active==NULL),"C1DFDEulerian::Disperse: Dispersion Scheme Not Initialized(2)",RUNTIME_ERR);
  ExitGracefullyIf((tstep <=0.0), "C1DFDEulerian::Disperse: negative or zero time step"          ,RUNTIME_ERR);

  //Copy aqueous concentration array
	for (k=0;k<nNodes;k++){for (s=0; s<nspecies; s++){Cnew[k][s]=C[k][s];}} 
	
	//Correct for source terms (cells sometimes damaged by Lagrangian transport)
	for (k=0;k<nCells;k++){
		for (s=0; s<nspecies; s++){
			double source_t=pDomain->GetDirichletConc(pGrid->GetNodeLocation(k),t,s);
			if (source_t!=NO_INFLUENCE){Cnew[k][s]=source_t;}
		}
	}

	local_t=t;
	
	while (local_t<t+tstep){
	
		loc_tstep=min(min(sugg_tstep,tstep/ceil(tstep/sugg_tstep)),t+tstep-local_t);
		local_t +=loc_tstep;

		ExitGracefullyIf((loc_tstep<=0.0),"C1DFDEulerian::Disperse: negative or zero local time step",RUNTIME_ERR);

		//introduce influence of advective/reactive terms--------------------------

		//(Eulerian fluxes are based upon weighted average (weight=aW) of estimated concentrations 
		// at beginning of local time step and end of local time step (calculated for simultaneous processes) )

		for (k=0;k<nNodes;k++){
			for (s=0; s<nspecies; s++){Cnew[k][s]+=(aW)*(loc_tstep/tstep)*(Cend[k][s]-C[k][s]);}}

		//Implicit/Explicit Formulation--------------------------------------------
	
	  //Adjust implicit matrix [A] for new time step
		//  TMP DEBUG-new species change A and AW matrices (diffusion coeff)
		//  TMP DEBUG-transient flow change A and AW matrices (sinks)
		if (last_loc_tstep!=loc_tstep){
			for (k=0;k<nCells;k++){
				A->AddToAij((1.0/last_loc_tstep)-(1.0/loc_tstep),k,k); 
			}
			last_loc_tstep=loc_tstep;
		}
		

		for (s=0; s<nspecies;s++){

			for (k=0;k<nCells;k++){Cs[k]=Cnew[k][s];dC[k]=B[k]=0.0;} //copy species concentration, reinitialize {B}
			
		  AW->MatVectMult(Cs,dC,nCells,false,false);               //obtains explicit portion of RHS (dC)

			//Reevaluate {B} vector 
			for (k=0;k<nCells;k++){
				pt  =pGrid->GetNodeLocation(k);
				area=pGrid->GetCellArea(k);

				if (!dirichlet[k]){
					sourceterm =pDomain->GetSourceConcentration  (k ,local_t,s)*sourceflux  [k]/area/poro[k]/thick[k];
					sourceterm+=pDomain->GetRechargeConcentration(pt,local_t,s)*rechargeflux[k]/area/poro[k]/thick[k];
					//sourceterm+=pDomain->GetSpecifiedMassFlux    (pt,local_t,s)                *area/thick[k]/poro[k]; //GetSpecifiedMassFlux=Mass flux per unit area
					B[k]=-dC[k]-Cs[k]/loc_tstep-sourceterm; // /R[k] 

					if (borderoutflx[k]<0.0){//and advection is being modeled??? //upgradient source condition
						//B[k]-=(1-w)*borderoutflx[k]*Cs[k]                     /area/poro[k]/thick[k]; ////This causes problems when ambient conditions !=0
						B[k]+=(1-w)*borderoutflx[k]*pDomain->GetAmbientConc(s)/area/poro[k]/thick[k];//
					}
					else if (borderoutflx[k]>0.0){//and advection is being modeled??? //downgradient condition - on LHS
						
					}
				}
				else{
					B[k]=(w  )*(pDomain->GetDirichletConc(pGrid->GetNodeLocation(k),local_t+loc_tstep,s))+
						   (1-w)*(pDomain->GetDirichletConc(pGrid->GetNodeLocation(k),local_t          ,s));
				}
			}

			//solve system of equations (for fully explicit, jumps out nearly immediately because [A]=[I])
			A->BCG(B,Cs,nCells,RELATIVE_RESIDUAL,tmp_err,1e-8);

			//pick up solution
			for (k=0;k<nCells;k++){Cnew[k][s]=Cs[k];}
		}
		
		//add in remainder of parallel advective/reactive term influence for this local time step
		for (k=0;k<nNodes;k++){
			for (s=0; s<nspecies; s++){Cnew[k][s]+=(1.0-aW)*(loc_tstep/tstep)*(Cend[k][s]-C[k][s]);}}
		
	}//end while
}

/***************************************************************************
							GET LOCAL MATRIX
****************************************************************************
returns local coefficient matrix
Optimization analysis:
	currently evaluates 5 points-2 of these have already been evaluated for previous cells
	currently evaluates velocity twice for each boundary (Dxx and Dxy/Dyy and Dyx)
	up to 5+5=10 Sat thick evaluations (evaluated for velocity, too)
	up to 5+5=10 Poro evaluations (evaluated for velocity, too)
	5 velocity calculations
	need new routine CAquifer->GetTransportInfo(&n,&T,&v,&effv,&leak1,&leak2,disp, bool effective)
---------------------------------------------------------------------------*/
void C1DFDEulerian::GetLocalMatrix(const int  k,
																	       double  a[3][3],
															     const double  dx[3], 
															     const cmplex  pts[2],
																   const int     species,   
																	 const double  up_weighting, 
																				 double &vloc, 
																				 double &dloc){



	static double hleft,hright,hcen;
	static double Dleft,Dright,Dxx,Dyy;;
	static double denom;
	static double vxright,vxleft,vx(0);
	static double alpha1,alpha2;
	static cmplex vel[12];
	static double R;
	double t(0.0); //TMP DEBUG
	double Df; //TMP DEBUG

	disp_info     disp;
	int           i1,j1;	

	int i,j;

	pGrid->k_to_ij(k,i,j);

	cmplex z[5];
	cmplex facept[4];

	z[0]=cmplex(pGrid->gridX[i  ],pGrid->gridY[j+1]);
	z[1]=cmplex(pGrid->gridX[i+1],pGrid->gridY[j+1]);
	z[2]=cmplex(pGrid->gridX[i+1],pGrid->gridY[j  ]);
	z[3]=cmplex(pGrid->gridX[i  ],pGrid->gridY[j  ]);
	z[4]=(z[0]+z[1]+z[2]+z[3])/4.0;

  facept[0]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j+1]));
  facept[1]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j+1]));
  facept[2]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j  ]));
  facept[3]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j  ]));

	   pDomain->GetDispersivities(t,species,disp); //TMP DEBUG
	R =pDomain->GetRetardation(c2Dto3D(pts[4]),t,species);  //TMP DEBUG-retardation 
	Df=pDomain->GetDiffusionCoeff(species);
	

//	R=pDomain->GetRetardation(c2Dto3D(z[4]),t,species);  //TMP DEBUG-retardation 

	for (i1=0; i1<=2; i1++){
		for (j1=0; j1<=2; j1++){a[i1][j1]=0.0;}}

	//Currently 1 point integration
	htop  =pLayer->GetPoro(pts[0])*pLayer->GetSaturatedThickness(pts[0],t);
	hright=pLayer->GetPoro(pts[1])*pLayer->GetSaturatedThickness(pts[1],t);
	hbot  =pLayer->GetPoro(pts[2])*pLayer->GetSaturatedThickness(pts[2],t);
	hleft =pLayer->GetPoro(pts[3])*pLayer->GetSaturatedThickness(pts[3],t);
	hcen  =pLayer->GetPoro(pts[4])*pLayer->GetSaturatedThickness(pts[4],t);
	
	/*htop  =pLayer->GetPoro(0.5*(z[0]+z[1]))*pLayer->GetSaturatedThickness(0.5*(z[0]+z[1]),t);
	hright=pLayer->GetPoro(0.5*(z[1]+z[2]))*pLayer->GetSaturatedThickness(0.5*(z[1]+z[2]),t);
	hbot  =pLayer->GetPoro(0.5*(z[2]+z[3]))*pLayer->GetSaturatedThickness(0.5*(z[2]+z[3]),t);
	hleft =pLayer->GetPoro(0.5*(z[3]+z[0]))*pLayer->GetSaturatedThickness(0.5*(z[3]+z[0]),t);
	hcen  =pLayer->GetPoro(z[4])*pLayer->GetSaturatedThickness(pts[4],t);*/


	//for correct mass balance, must use integrated normal velocity
	// integrated tangential velocity is more accurate for dispersion cross-terms
	vel[0]= IM*pLayer->GetFluxThruFace(facept[0],facept[1],t)/dx[1]/htop; 
	vel[1]=    pLayer->GetFluxThruFace(facept[3],facept[0],t)/dy[1]/hright;
	vel[2]= IM*pLayer->GetFluxThruFace(facept[3],facept[2],t)/dx[1]/hbot;
	vel[3]=    pLayer->GetFluxThruFace(facept[2],facept[1],t)/dy[1]/hleft;

	if (htop  <=0.0){vel[0]=0.0;}
	if (hright<=0.0){vel[1]=0.0;}
	if (hbot  <=0.0){vel[2]=0.0;}
	if (hleft <=0.0){vel[3]=0.0;}

	
	vxright=vel[1].real()*hright/hcen/R; //right 
	vxleft =vel[3].real()*hleft /hcen/R; //left  
	vx =0.5*(vxright+vxleft );

	vytop  =vel[0].imag()*htop  /hcen/R; //top    
	vybot  =vel[2].imag()*hbot  /hcen/R; //bottom 
	vy =0.5*(vytop+vybot);

	vloc=vx;

	if (process_type==ADV_AND_DISP)
	{
		
		int dir=0;

		//from -dQxC/dx (eqn 7.54 Zheng & Bennett, 2003)---------------------------------- 
		//--------------------------------------------------------------------------------
		alpha1=alpha2=0.5;
		if (IsActive(pGrid->ij_to_k(i+1,j)) && (vxright!=0.0)  && (dir!=2)){
			alpha1=((up_weighting)*(1.0-(vxright)/fabs(vxright))+(1.0-up_weighting))*0.5;}
		if (IsActive(pGrid->ij_to_k(i-1,j)) && (vxright!=0.0)  && (dir!=2)){
			alpha2=((up_weighting)*(1.0-(vxleft )/fabs(vxleft ))+(1.0-up_weighting))*0.5;}

		denom=0.5*((1.0-alpha1)*dx[0]+dx[1]+alpha1*dx[2]);
		denom=dx[1];
		a[1][1]-=(1.0-alpha1)*vxright/denom;
		a[2][1]-=(    alpha1)*vxright/denom;

		denom=0.5*((1.0-alpha2)*dx[0]+dx[1]+alpha2*dx[2]);
		denom=dx[1];
		a[0][1]+=(1.0-alpha2)*vxleft /denom;
		a[1][1]+=(    alpha2)*vxleft /denom;
	
	}




	//1/nh d/dx(nh Dxx dC/dx) (eqn 7.47 Zheng & Bennett, 2003)---------------------------
	//---------------------------------------------------------------------------------
	Dleft  =hleft *pDomain->CalcDispersionCoeff(c2Dto3D(vel[3]),OR_XX,x_orient)/R;
	Dright =hright*pDomain->CalcDispersionCoeff(c2Dto3D(vel[1]),OR_XX,x_orient)/R;
	Dleft +=hleft *Df/R;
	Dright+=hright*Df/R;

	Dxx=Dright;	

	if (IsActive(k+1)){
    denom=(0.5*(dx[2]+dx[1])*dx[1])*hcen;
		a[2][1]+=Dright/denom;//i+1,j
	  a[1][1]-=Dright/denom;//i  ,j
	}
	if (IsActive(k-1){
		denom=(0.5*(dx[1]+dx[0])*dx[1])*hcen;
		a[0][1]+=Dleft /denom;//i-1,j
		a[1][1]-=Dleft /denom;//i  ,j
	}



	//cout <<"Dleft:"<<Dleft<<" Dright: "<<Dright<<" Dbot:"<<Dbot<<" Dtop:"<<Dtop<<endl;
	/*cout <<"local a[][]:"<<endl;
	cout <<a[0][2]<<" "<<a[1][2]<<" "<<a[2][2]<<endl;
	cout <<a[0][1]<<" "<<a[1][1]<<" "<<a[2][1]<<endl;
	cout <<a[0][0]<<" "<<a[1][0]<<" "<<a[2][0]<<endl;
	ExitGracefully("f",BAD_DATA);*/

	dloc=Dxx

}
//----------------------------------------------------------------------
void C1DFDEulerian::WriteOutput(const int outputstep){
}
//----------------------------------------------------------------------
bool C1DFDEulerian::IsActive(const int k) const {
	return ((k>=0) && (k<nCells));
}
//----------------------------------------------------------------------
double C1DFDEulerian::CalculateSystemMass(const        int  s,
																				  Ironclad2DArray   C,
																				  Ironclad2DArray   Cs,
																					         double  &SorbedMass) const{
	//returns system mass in implicit mass units (usually mg)
	static double AqMass(0.0),pb,area;
	static pt3D   pt;
	int    k;
	
	AqMass=0.0;
	SorbedMass=0.0;
	
	for (k=0; k<nCells; k++){
		pt  =pGrid->GetNodeLocation(k);
		area=pGrid->GetCellArea(k);
		pb  =pDomain  ->GetBulkDryDensity(pt);
			
		AqMass    += C[k][s]*area*poro[k]*thick[k]*pDomain->Vratio();	//vratio converts to L from system units (e.g., m^3)
		SorbedMass+=Cs[k][s]*area*pb     *thick[k]*pDomain->Vratio() ;
	}

	return AqMass+SorbedMass;
}
//----------------------------------------------------------------------
void C1DFDEulerian::InitializeBorderFlux(const double &t){
	cmplex z1,z2;
	int i,j,k;
  //nBorderCells=2*pGrid->nX+2*pGrid->nY-4;
  //border outflux is greater than zero for outflux, negative for influx
	for (k=0; k<nCells;k++){
		borderoutflx[k]=0.0;
	}
	for (i=0; i<pGrid->nX;i++){

		z1=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i          ],pGrid->gridY[0         ]));
		z2=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1        ],pGrid->gridY[0         ]));
		k=pGrid->ij_to_k(i,0);
		borderoutflx[k]+=pLayer->GetFluxThruFace(z1,z2,t).real(); //bottom
		
		z1=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1        ],pGrid->gridY[pGrid->nY-1]));
		z2=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i          ],pGrid->gridY[pGrid->nY-1]));
		k=pGrid->ij_to_k(i,pGrid->nY-1);
		borderoutflx[k]+=pLayer->GetFluxThruFace(z1,z2,t).real(); //top
	}

	for (j=0; j<pGrid->nY;j++){

		z1=pGrid->LocalToGlobal(cmplex(pGrid->gridX[0          ],pGrid->gridY[j+1        ]));
		z2=pGrid->LocalToGlobal(cmplex(pGrid->gridX[0          ],pGrid->gridY[j          ]));
		k=pGrid->ij_to_k(0,j);           
		borderoutflx[k]+=pLayer->GetFluxThruFace(z1,z2,t).real();//left
		
		z1=pGrid->LocalToGlobal(cmplex(pGrid->gridX[pGrid->nX-1],pGrid->gridY[j          ]));
		z2=pGrid->LocalToGlobal(cmplex(pGrid->gridX[pGrid->nX-1],pGrid->gridY[j+1        ]));
		k=pGrid->ij_to_k(pGrid->nX-1,j); 
		borderoutflx[k]+=pLayer->GetFluxThruFace(z1,z2,t).real();//right
		
	}
	/*for (k=0; k<nCells;k++){
		pGrid->k_to_ij(k,i,j);
		if (borderoutflx[k]!=0){
		cout <<"borderoutflx["<<i<<","<<j<<"]:"<<borderoutflx[k]<<" area:"<<pGrid->GetCellArea(k) <<endl;
		}
	}*/
}
//----------------------------------------------------------------------
double C1DFDEulerian::CalculateBorderFlux(const int        s,
																				  Ironclad2DArray  C,
																					Ironclad2DArray  Cend,
																					const double    &t,
																					const double    &tstep) const{
  //for each cell along border, calculate net loss of contaminant [mg/t]
	//Assumes zero dispersive flux along system boundary

	int    k;
	double sum=0.0;


	for (k=0; k<nCells; k++){

		//flux out of system (borderflux>0)
		sum+=0.5*(C[k][s]+Cend[k][s])*max(borderoutflx[k],0.0); 

		//if (fabs(borderoutflx[k])>0){cout<<borderoutflx[k]<<" "<<C[k][s]<<" "<<Cend[k][s]<<" "<<poro[k]<<endl;}
	}
	
	return sum*tstep*pDomain->Vratio();

}
//----------------------------------------------------------------------
double C1DFDEulerian::CalculateSourceGain(const int        s,
																		      Ironclad2DArray  C,
																					Ironclad2DArray  Cend,
															            const double    &t,
																					const double    &tstep) const{
	//combination of source fluxes and dirichlet
	CVector x_orientation(cos(pGrid->orient),sin(pGrid->orient),0.0);

	double sum(0.0),dispsum(0.0),dirconc;
	double Dxx1,Dxx2,Dyy1,Dyy2,dx,dy;
	double h[4];
	cmplex Q[4],z[4];
	pt3D pt;
	int i,j,k,tmpk;
	double advsum(0.0);
	for (k=0; k<nCells; k++){
		pt=pGrid->GetNodeLocation(k);
		//source flux-------------------------------------------------------
		sum+=((1.0-w)*pDomain->GetSourceConcentration(k,t      ,s)+
			        (w)*pDomain->GetSourceConcentration(k,t+tstep,s))*sourceflux[k];

		//recharge source---------------------------------------------------
		if (rechargeflux[k]>0.0){
			sum+=((1.0-w)*pDomain->GetRechargeConcentration(pt,t      ,s)+
			          (w)*pDomain->GetRechargeConcentration(pt,t+tstep,s))*rechargeflux[k];
		}

    //flux into system (borderoutflux<0)-----------------------------------
		sum-=pDomain->GetAmbientConc(s)*min(borderoutflx[k],0.0); 

		//dirichlet flux (should be optimized)------------------------------
		dirconc= (1.0-w)*pDomain->GetDirichletConc(pt,t      ,s)+
			           (w)*pDomain->GetDirichletConc(pt,t+tstep,s);
		
		//fixes problems with non-smooth dirichlet conditions from init conditions
		if ((t==0) && (pDomain->GetDirichletConc(pt,t+tstep,s)==NO_INFLUENCE)){return sum*tstep; }

		if (dirconc!=NO_INFLUENCE){
			//for each side, calculate flux q=d/dn(ds*Q_nC +ds*nhD dC/dn)
			//sum+=(sourceflux[k]-sinkflux[k]+rechargeflux[k])*dirconc;

			pGrid->k_to_ij(k,i,j);
			z[0]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j  ]));
			z[1]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j  ]));
			z[2]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j+1]));
			z[3]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j+1]));

			//the slow part:
			Q[0]= pLayer->GetFluxThruFace(z[1],z[0],t);  h[0]=pLayer->GetSaturatedThickness((z[0]+z[1])/2.0,t)*pLayer->GetPoro((z[0]+z[1])/2.0);
			Q[1]=-pLayer->GetFluxThruFace(z[2],z[1],t);  h[1]=pLayer->GetSaturatedThickness((z[1]+z[2])/2.0,t)*pLayer->GetPoro((z[1]+z[2])/2.0); 			
			Q[2]=-pLayer->GetFluxThruFace(z[3],z[2],t);  h[2]=pLayer->GetSaturatedThickness((z[2]+z[3])/2.0,t)*pLayer->GetPoro((z[2]+z[3])/2.0);
			Q[3]= pLayer->GetFluxThruFace(z[0],z[3],t);  h[3]=pLayer->GetSaturatedThickness((z[3]+z[0])/2.0,t)*pLayer->GetPoro((z[3]+z[0])/2.0);

			Dxx1=pDomain->CalcDispersionCoeff(c2Dto3D(   Q[3]/h[3]/abs(z[3]-z[0])),OR_XX,x_orientation);
			Dxx2=pDomain->CalcDispersionCoeff(c2Dto3D(   Q[1]/h[1]/abs(z[2]-z[1])),OR_XX,x_orientation);
			Dyy1=pDomain->CalcDispersionCoeff(c2Dto3D(IM*Q[0]/h[0]/abs(z[1]-z[0])),OR_YY,x_orientation);
			Dyy2=pDomain->CalcDispersionCoeff(c2Dto3D(IM*Q[2]/h[2]/abs(z[3]-z[2])),OR_YY,x_orientation);

			if (nspecies==1){
				Dxx1+=pDomain->GetDiffusionCoeff(0);
				Dxx2+=pDomain->GetDiffusionCoeff(0);
				Dyy1+=pDomain->GetDiffusionCoeff(0);
				Dyy2+=pDomain->GetDiffusionCoeff(0);
			}
      //cout <<Q[0]<<" "<<Q[1]<<" "<<Q[2]<<" "<<Q[3]<<endl;	
			//cout <<"Dxx:"<<Dxx1<<" "<<Dxx2<<" "<<Dyy1<<" "<<Dyy2<<endl;
			//cout <<"h:"<<h[0]<<" "<<h[1]<<" "<<h[2]<<" "<<h[3]<<endl;
			
			double tmpwt=1.0;
			double C1;

			//x-direction component 1 (left)----------------------------
			tmpk=pGrid->ij_to_k(i-1,j);
			if (tmpk!=OFF){ C1=((1.0-w)*C[tmpk][s]+w*Cend[tmpk][s]);}
			else          { C1=pDomain->GetAmbientConc(s);          }

			advsum-=max(Q[3].real()*C1,0.0); //gain
			advsum-=min(Q[3].real()*dirconc,0.0); 
				
			dx=pGrid->nodeX[i  ]-pGrid->nodeX[i-1];
			dy=pGrid->gridY[j+1]-pGrid->gridY[j  ];
			dispsum+=Dxx1*h[3]*(dirconc-C1)/dx*dy*tmpwt;
			
			//x-direction component 2 (right)-----------------------------
			tmpk=pGrid->ij_to_k(i+1,j);
			if (tmpk!=OFF){ C1=((1.0-w)*C[tmpk][s]+w*Cend[tmpk][s]);}
			else          { C1=pDomain->GetAmbientConc(s);          }

			advsum+=max(Q[1].real()*C1,0.0);
			advsum+=min(Q[1].real()*dirconc,0.0);

			dx=pGrid->nodeX[i+1]-pGrid->nodeX[i  ];
			dy=pGrid->gridY[j+1]-pGrid->gridY[j  ];
			dispsum+=Dxx2*h[1]*(dirconc-C1)/dx*dy*tmpwt;
			
			//y-direction component 1 (bottom)----------------------------
			tmpk=pGrid->ij_to_k(i,j-1);
			if (tmpk!=OFF){ C1=((1.0-w)*C[tmpk][s]+w*Cend[tmpk][s]);}
			else          { C1=pDomain->GetAmbientConc(s);          }

			advsum-=max(Q[0].imag()*C1,0.0);
			advsum-=min(Q[0].imag()*dirconc,0.0);

			dx=pGrid->gridX[i+1]-pGrid->gridX[i  ];
			dy=pGrid->nodeY[j  ]-pGrid->nodeY[j-1];
			dispsum+=Dyy1*h[0]*(dirconc-C1)/dy*dx*tmpwt;
			
			//y-direction component 2 (top)--------------------------------
			tmpk=pGrid->ij_to_k(i,j+1);
			if (tmpk!=OFF){ C1=((1.0-w)*C[tmpk][s]+w*Cend[tmpk][s]);}
			else          { C1=pDomain->GetAmbientConc(s);          }

			advsum+=max(Q[2].imag()*C1,0.0);
			advsum+=min(Q[2].imag()*dirconc,0.0);

			dx=pGrid->gridX[i+1]-pGrid->gridX[i  ];
			dy=pGrid->nodeY[j+1]-pGrid->nodeY[j  ];
			dispsum+=Dyy2*h[2]*(dirconc-C1)/dy*dx*tmpwt;

		//	cout <<"advective  net dirichlet flux : "<<sum<<endl;
			//cout <<"dispersive net dirichlet flux : "<<dispsum<<endl;
		}
	}

	//for each const conc. cell, calculate mass flux from sides
	/*cout <<"advsum:"<<sum<<" , dispsum: "<<dispsum<<" lsum:"<<tmpsum<<" rdum:"<<tmpsum2<<" lneg:"<<tmpsum3<<" rneg:"<<tmpsum4<<endl;
	ofstream MBTMP;
	MBTMP.open("crappy.csv",ios::app);
  MBTMP<<sum<<" , "<<dispsum<<" ,"<<tmpsum<<" ,"<<tmpsum2<<" ,"<<tmpsum3<<" ,"<<tmpsum4<<endl;
	MBTMP.close();*/
	if (Dxx1!=0){return (sum+0.5*advsum+dispsum)*tstep*pDomain->Vratio();}//0.5 factor Still not understood
	else        {return (sum+    advsum+dispsum)*tstep*pDomain->Vratio();}
}
//-----------------------------------------------------------
double C1DFDEulerian::CalculateSinkLoss  (const int        s,
																		      Ironclad2DArray  C,
																					Ironclad2DArray  Cend,
															            const double    &t,
																					const double    &tstep) const{
  //returns [mg/t]
	double sum(0.0);
	int k;
  for (k=0; k<nCells; k++){
		if (sinkflux[k]>0.0){
			//should really use upstream concentrations as estimate
			sum+=0.5*(C[k][s]+Cend[k][s])*sinkflux[k];
		}
		//if (rechargeflux[k]<0.0){
			//sum-=0.5*(C[k][s]+Cend[k][s])*rechargeflux[k]; //recharge is portion of sink flux
		//}
	}
	
	return sum*tstep*pDomain->Vratio();
}
//----------------------------------------------------------------------
double C1DFDEulerian::CalculateDecayLoss(const int        s,
														             Ironclad2DArray  C,
																				 Ironclad2DArray  Cend,
														             const double    &t,
																			   const double    &tstep) const{
	double sum(0.0);
	int k;
	double lambda=pDomain->GetSpeciesDecay(s);
  for (k=0; k<nCells; k++){
		sum+=0.5*(C[k][s]+Cend[k][s])*lambda*pGrid->GetCellArea(k)/poro[k]/thick[k];
	}
	
	return sum*tstep*pDomain->Vratio();
}
/***************************************************************************
							CellByCellMassBalance
****************************************************************************
TEMPORARY TEST
---------------------------------------------------------------------------*/
double C1DFDEulerian::CellByCellMB(const int        s,
																	 Ironclad2DArray  C,
																	 Ironclad2DArray  Cend,
																	 const double     up_weighting,
															     const double    &t,
																	 const double    &tstep){
/*  double MBerr(0.0);
	double Dxx1,Dxx2,Dyy1,Dyy2,dx,dy;
	double h[4];
	cmplex Q[4],z[4];

	int i,j;
	int _B=0;
	int _R=1;
	int _T=2;
	int _L=3;

  ofstream MB;
	MB.open("CellByCellMassBalance.csv",ios::app);
	
	for (k=0; k<pGrid->nCells; k++){

		pGrid->k_to_ij(k,i,j);

		kN =pGrid->ij_to_k(i  ,j+1);
		kS =pGrid->ij_to_k(i  ,j-1);
		kE =pGrid->ij_to_k(i+1,j  );
		kW =pGrid->ij_to_k(i-1,j  );
		kNW=pGrid->ij_to_k(i-1,j+1);
		kSW=pGrid->ij_to_k(i-1,j-1);
		kNE=pGrid->ij_to_k(i+1,j+1);
		kSE=pGrid->ij_to_k(i+1,j-1);

		WXP=dx[i+1]/(dx[i]+dx[i+1]);
    WXM=dx[i]  /(dx[i]+dx[i-1]);
    WYP=dy[j+1]/(dy[j]+dy[j+1]);
    WYM=dy[j  ]/(dy[j]+dy[j-1]);

		z[0]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j  ]));
		z[1]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j  ]));
		z[2]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j+1]));
		z[3]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j+1]));

		//the slow part:
		Q[_B]=pLayer->GetFluxThruFace(z[0],z[1],t);h[_B]=pLayer->GetSaturatedThickness((z[0]+z[1])/2.0,t)*pLayer->GetPoro((z[0]+z[1])/2.0);
		Q[_R]=pLayer->GetFluxThruFace(z[1],z[2],t);h[_R]=pLayer->GetSaturatedThickness((z[1]+z[2])/2.0,t)*pLayer->GetPoro((z[1]+z[2])/2.0); 			
		Q[_T]=pLayer->GetFluxThruFace(z[3],z[2],t);h[_T]=pLayer->GetSaturatedThickness((z[2]+z[3])/2.0,t)*pLayer->GetPoro((z[2]+z[3])/2.0);
		Q[_L]=pLayer->GetFluxThruFace(z[0],z[3],t);h[_L]=pLayer->GetSaturatedThickness((z[3]+z[0])/2.0,t)*pLayer->GetPoro((z[3]+z[0])/2.0);

		//cout <<Q[0]<<" "<<Q[1]<<" "<<Q[2]<<" "<<Q[3]<<endl;	  
		
		Dxx1=pDomain->GetDispersionCoeff(c2Dto3D(    Q[_L]/h[_L]),t,OR_XX,x_orientation);
		Dxx2=pDomain->GetDispersionCoeff(c2Dto3D(    Q[_R]/h[_R]),t,OR_XX,x_orientation);
		Dyy1=pDomain->GetDispersionCoeff(c2Dto3D(-IM*Q[_B]/h[_B]),t,OR_YY,x_orientation);
		Dyy2=pDomain->GetDispersionCoeff(c2Dto3D(-IM*Q[_T]/h[_T]),t,OR_YY,x_orientation);
	
		if (kW!=OFF){ //x-direction component 1 (left)
			ww=dx[i]/(dx[i]+dx[i-1]);
			alpha=0.0;
			if (Q[_L]>0){alpha=(up_weighting)+(1.0-up_weighting)*ww;}

			in-=Q[_L]*(C[kW][s]*alpha+C[kW][s]*(1.0-alpha))*dy[j]*h[_L]/poro[k]*tstep;

			in+=Dxx1*(Cnew[k][s]-Cnew[kW][s]);



			in+=Dxy1*(Cnew[tmpkjp1][s]*(1.0-WXM)
               

		}
		MB << MBerr <<",";
	}




	MB.close();*/
	return 0.0;
}



/***************************************************************************
							NINE POINT STENCIL
****************************************************************************
returns approximations of transverse terms, d/di (D dC/dj)

---------------------------------------------------------------------------*/
double C1DFDEulerian::NinePointStencil(const double      D1, 
																	 const double      D2, 
																	 const orientation dir, 
																	 const double      C     [3][3], 
																	 const bool        active[3][3], 
																	 const double      dx[3],
																	 const double      dy[3]){

	/*
			|			|		
	0,2 | 1,2 |	2,2			   |		|
	____|_____|_____		-C1+--C2+---
			|			|						 |		|	
	0,1	|	1,1	|	2,1			-C3+--C4+---
	____|_____|_____			 |    | 
			|			|
	0,0	|	1,0	|	2,0
			|			|			
	*/


	if (!active[1][1]) {return 0.0;}

	double c1,c2,c3,c4; //weighted concentrations at corners of cell 1,1
	double w1,w2,w3,w4; //local weights

	if (active[0][2]){w1=dx[1]*dy[1];}else{w1=0.0;}
	if (active[1][2]){w2=dx[0]*dy[1];}else{w2=0.0;}
	if (active[0][1]){w3=dx[1]*dy[2];}else{w3=0.0;}
	if (active[1][1]){w4=dx[0]*dy[2];}else{w4=0.0;}

	c1=(w1*C[0][2]+w2*C[1][2]+w3*C[0][1]+w4*C[1][1])/(w1+w2+w3+w4);

	if (active[1][2]){w1=dx[2]*dy[1];}else{w1=0.0;}
	if (active[2][2]){w2=dx[1]*dy[1];}else{w2=0.0;}
	if (active[1][1]){w3=dx[2]*dy[2];}else{w3=0.0;}
	if (active[2][1]){w4=dx[1]*dy[2];}else{w4=0.0;}

	c2=(w1*C[1][2]+w2*C[2][2]+w3*C[1][1]+w4*C[2][1])/(w1+w2+w3+w4);

	if (active[0][1]){w1=dx[1]*dy[0];}else{w1=0.0;}
	if (active[1][1]){w2=dx[0]*dy[0];}else{w2=0.0;}
	if (active[0][0]){w3=dx[1]*dy[1];}else{w3=0.0;}
	if (active[1][0]){w4=dx[0]*dy[1];}else{w4=0.0;}

	c3=(w1*C[0][1]+w2*C[1][1]+w3*C[0][0]+w4*C[1][0])/(w1+w2+w3+w4);

	if (active[1][1]){w1=dx[2]*dy[0];}else{w1=0.0;}
	if (active[2][1]){w2=dx[1]*dy[0];}else{w2=0.0;}
	if (active[1][0]){w3=dx[2]*dy[1];}else{w3=0.0;}
	if (active[2][0]){w4=dx[1]*dy[1];}else{w4=0.0;}

	c4=(w1*C[1][1]+w2*C[2][1]+w3*C[1][0]+w4*C[2][0])/(w1+w2+w3+w4);
	
	switch (dir){
		case(OR_XY):{
			return (D1*(c2-c4)-D2*(c1-c3))/(dx[1]*dy[1]);
			break;
		}
		case(OR_YX):{
			return (D1*(c2-c1)-D2*(c4-c3))/(dx[1]*dy[1]);		
			break;
		}
	}
  return 0.0;
}
#endif