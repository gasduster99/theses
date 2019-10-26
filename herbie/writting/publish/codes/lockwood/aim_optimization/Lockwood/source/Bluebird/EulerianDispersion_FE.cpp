//Eulerian Dispersion_FE.cpp

#include "FEEulerian.h"
#include "MasterInclude.h"

const bool theta_on_left=true;
const bool conservative=true;
const bool turn_off_cross_terms=false;
//***************************************************************************
C2DFEEulerian::C2DFEEulerian(C2DDomainABC *pDom, CFEMesh *pConcGrid, const double time_weight, const double up_weight)
              :C2DTransportScheme(pDom){

	ExitGracefullyIf(((time_weight>1.0) || (time_weight<0.0)),"C2DFDEulerian: Invalid value for temporal weighting term (<0 or >1)",BAD_DATA);
	ExitGracefullyIf(((up_weight>1.0)   || (up_weight<0.0)  ),"C2DFDEulerian: Invalid value for upstream weighting term (<0 or >1)",BAD_DATA);

	pMesh				 =pConcGrid;
	nNodes			 =pConcGrid->GetNumNodes();

	poro				 =NULL;
	thick				 =NULL;
	B						 =NULL;
	sourceflux	 =NULL;
	sinkflux		 =NULL;
	borderflux   =NULL;
	rechargeflux =NULL;
	singularities=NULL;
	wallflux		 =NULL;
	walllength	 =NULL;
	wallthick		 =NULL;
	dirflux			 =NULL;
	tmpC				 =NULL;
	M						 =NULL;
	Df					 =NULL;
	A						 =NULL;
	g_htheta		 =NULL;
	dirichlet		 =NULL;

	w								 =time_weight;
	upwt						 =up_weight;

	use_traditional	 =false;
	LinearParams		 =false; 
	SUPG						 =false;
	lumped					 =true;
	gausstype        =TG_7POINT;
	print_diagnostics=false;
}
//***************************************************************************
C2DFEEulerian::~C2DFEEulerian(){
	if (globaldebug){cout <<"DESTROYING FINITE ELEMENT DISPERSION SCHEME "<<endl;}
	delete [] B;
	delete [] sourceflux;
	delete [] sinkflux;
	delete [] borderflux;
	delete [] rechargeflux;
	delete [] singularities;
	delete [] tmpC;
	delete [] wallflux;
	delete [] walllength;
	delete [] wallthick;
	delete [] dirflux;
	if (g_htheta!=NULL){for (int e=0; e<pMesh->NumElems; e++){delete [] g_htheta[e];}}
	delete [] g_htheta;
	delete M;
	delete Df;
	delete A;
}
/***************************************************************************
							PRINT GRID DIAGNOSTICS
****************************************************************************
---------------------------------------------------------------------------*/
void C2DFEEulerian::PrintMeshDiagnostics(){print_diagnostics=true;}
/***************************************************************************
							Parameters
****************************************************************************
---------------------------------------------------------------------------*/
void C2DFEEulerian::UseSUPG(){SUPG=true;}
//--------------------------------------------------------------------------
void C2DFEEulerian::SetParameters (trigausspoints  gausslevel, bool traditional, bool lump){
  gausstype       =gausslevel;
	use_traditional =traditional;
	lumped          =lump;
	LinearParams    =false; 
	SUPG            =false;
	
	//check quality of gauss integration scheme--------------------------------
	int nG=GetNumTriGaussPts(gausstype);
	cmplex zg[MAXGAUSSPTS];
	double gw[MAXGAUSSPTS];
	double sum(0);
	GetTriGaussInfo(0.0,IM,1.0,gausstype,zg,gw);
	for (int n=0;n<nG;n++){sum+=gw[n];}
	ExitGracefullyIf(fabs(sum-1.0)>REALSMALL,
		"C2DFEEulerian::SetParameters: gauss weighting scheme is improperly coded",RUNTIME_ERR);
}
/***************************************************************************
							INITIALIZE
****************************************************************************
Creates Vertically Averaged FE Matrices
---------------------------------------------------------------------------*/
void C2DFEEulerian::Initialize(const double         startt,
		                           const transport_type ty,
									             Ironclad1DArray      porosity, 
												       Ironclad1DArray      satthick){

	ExitGracefullyIf((pDomain==NULL),     "C2DFEEulerian::Initialize: NULL Domain",RUNTIME_ERR);
	ExitGracefullyIf((pMesh  ==NULL),     "C2DFEEulerian::Initialize: NULL Grid  ",RUNTIME_ERR);
  ExitGracefullyIf((ty==ADVECTION_ONLY),"C2DFEEulerian::Initialize: Cannot simulate pure advection",BAD_DATA);
	

	int    e,i,j,k;
	cmplex tmpz;
	int    nde            [3];
	double elementmatrix  [3][3];
	double diffusionmatrix[3][3];
	double sorptionmatrix [3][3];
	double *peclet;
	double tmpa; 
	double PeMax(0.0),PeAvg(0.0),CrMax(0.0),CrAvg(0.0),VelAvg(0.0);
	double vloc,dloc,loc_len;
	pt3D   pt;
	double t=startt+0.0001; //so as to avoid initial conditions in Dirichlet
  

	//ofstream TMP; //TMP DEBUG
	//TMP.open("element_vel.csv");

	cout << "Initializing Finite Element Scheme..."<<endl;
  if (ProgramAborted()){return;}

	sugg_tstep=ALMOST_INF;
	
	process_type =ty;
	poro         =porosity;
  thick        =satthick;
	nspecies     =pDomain->GetNumSpecies();

	borderflux   =new double  [nNodes];
	sourceflux   =new double  [nNodes];
	sinkflux     =new double  [nNodes];
	rechargeflux =new double  [nNodes];
	singularities=new double  [nNodes];

	wallflux     =new double  [nNodes];
	walllength   =new double  [nNodes];
	wallthick    =new double  [nNodes];
	B            =new double  [nNodes];
	tmpC         =new double  [nNodes];
	dirichlet    =new bool    [nNodes];
	dirflux      =new double  [nNodes];
	g_htheta     =new double *[pMesh->NumElems];
	peclet       =new double  [pMesh->NumElems];
 
	//dankwerts1 =new double  [nNodes];
	//dankwerts2 =new double  [nNodes];	

  ExitGracefullyIf((peclet==NULL),"C2DFEEulerian::Initialize: Out of Memory  ",OUT_OF_MEMORY);

	for (k=0;k<nNodes;k++){
		sinkflux     [k]=sourceflux[k]=borderflux[k]=rechargeflux[k]=B[k]=tmpC[k]=0.0; 
		singularities[k]=0.0;
		dirflux      [k]=0.0; 
		dirichlet    [k]=false;
		wallflux     [k]=0.0;
		walllength   [k]=0.0;
		wallthick    [k]=0.2; //TMP DEBUG-should not be hard coded
	}
	for (e=0;e<pMesh->NumElems; e++){
		g_htheta[e]=new double [MAXGAUSSPTS];
		ExitGracefullyIf((g_htheta[e]==NULL),"C2DFEEulerian::Initialize: Out of Memory(1)  ",OUT_OF_MEMORY);
		for (k=0;k<MAXGAUSSPTS;k++){
			g_htheta[e][k]=1.0;
		}
	}

	cout << "  Processing and Discretizing Fluxes..."<<endl;

	//Figure out Source/Sink Flux Boundary conditions---------------------------------------------
	//Figure out Properties of Discontinuous Boundaries-------------------------------------------
	//Figure out Flux along system Boundaries-----------------------------------------------------
	
	pLayer=pDomain->GetLayer();

	cmplex z1,z2,z3;
	int    n1,n2,n3,s;
	double Q1,Q2,Qn1,Qn2,L,sideflux;
	cmplex zg1,zg2;                  //two gauss points along sides

  //count number of sides that are linked to each node---------------------------
	int *count=new int [nNodes];
	for (k=0;k<nNodes;k++){count[k]=0;}
	for (s=0; s<pMesh->NumSides; s++){
		count[pMesh->side[s].n1]++; 
		count[pMesh->side[s].n2]++; 
	}

	//Recharge/Leakage Source/Sink Fluxes------------------------------------------
	for (e=0; e<pMesh->NumElems; e++){
		n1=pMesh->elem[e].ni;  z1=pMesh->node[n1].z;
		n2=pMesh->elem[e].nj;  z2=pMesh->node[n2].z;	
		n3=pMesh->elem[e].nk;  z3=pMesh->node[n3].z;
		
		//ignores curvature- this is OK
		Q1=pLayer->GetIntegratedLeakage(z1,z3,z2,t,FROMTOPANDBOTTOM);

		//TMP DEBUG- should differentiate for source recharge
		if      (fabs(Q1)<REALSMALL){}//cout <<"no leakage/recharge  "<<Q1<<endl;} 
		else if (Q1> 0.0){
			sinkflux  [n1]+=Q1/3.0; sinkflux  [n2]+=Q1/3.0; sinkflux  [n3]+=Q1/3.0; //cout <<"leakage:  "<<Q1<<endl;
		}
		else             {
			rechargeflux[n1]-=Q1/3.0; rechargeflux[n2]-=Q1/3.0; rechargeflux[n3]-=Q1/3.0; 
			//cout <<"recharge:"<<Q1<<" "<< Q1/pMesh->elem[e].area <<endl;
		} 
	}

	for (s=0; s<pMesh->NumSides; s++){

		n1=pMesh->side[s].n1;  z1=pMesh->node[n1].z;
		n2=pMesh->side[s].n2;  z2=pMesh->node[n2].z;
		
    L=abs(z2-z1);

    //Sink/Source Flux Distribution--------------------------------------------
		if ((pMesh->side[s].status!=DOUBLENODE))
		{

			pLayer->GetFluxDistribution(z1,z2,t,Q1,Q2);

			if ( (fabs(Q1)<REALSMALL) || (fabs(Q2)<REALSMALL) )//point source/sink
			{ 
				if      (fabs(Q1)<REALSMALL){} 
				else if (Q1< 0.0){sinkflux  [n1]-=Q1/count[n1]; }//sink
				else             {sourceflux[n1]+=Q1/count[n1]; }//source
				if      (fabs(Q2)<REALSMALL){} 
				else if (Q2< 0.0){sinkflux  [n2]-=Q2/count[n2]; }//sink
				else             {sourceflux[n2]+=Q2/count[n2]; }//source
			}
			else                                              {//linear/circular/elliptical source/sink
				if      (fabs(Q1)<REALSMALL){} 
				else if (Q1< 0.0){sinkflux  [n1]-=Q1;           }//sink
				else             {sourceflux[n1]+=Q1;           }//source
				if      (fabs(Q2)<REALSMALL){} 
				else if (Q2< 0.0){sinkflux  [n2]-=Q2;           }//sink
				else             {sourceflux[n2]+=Q2;           }//source
			}

		}


    //Leaky Wall Boundaries---------------------------------------------------
		else if ((pMesh->side[s].status==DOUBLENODE))
		{
			//wallflux is positive if entering domain

			zg1=0.21132*(z2-z1)+z1;
      zg2=0.78867*(z2-z1)+z1;			
			Qn1=-(conj(IM*(z2-z1)/L)*conj(pLayer->GetW(zg1,0.0))).real();//normal flux at 2 gauss points
			Qn2=-(conj(IM*(z2-z1)/L)*conj(pLayer->GetW(zg2,0.0))).real();

			if      ((pMesh->side[s].ea==OFF) || (!pMesh->elem[pMesh->side[s].ea].on)) 
			{
					
				/*sideflux=pLayer->GetFluxThruFace(z2,z1,t).real(); //Uses lame approximation to integral  | Ni*Qn dX = | N dX * | Qn dX
				wallflux[n1]-=0.5*sideflux;
				wallflux[n2]-=0.5*sideflux;*/

				wallflux[n1]+=L*0.5*(0.78867*Qn1+0.21132*Qn2);      //Uses Gauss integration of | Ni*Qn dX
				wallflux[n2]+=L*0.5*(0.21132*Qn1+0.78867*Qn2);

				wallflux[n1]=wallflux[n2]=0.0; //TMP DEBUG

				walllength[n1]+=0.5*abs(z1-z2);
				walllength[n2]+=0.5*abs(z1-z2);
				wallthick [n1]=wallthick[n2]=0.2; //TMP DEBUG 
			}
			else if ((pMesh->side[s].eb==OFF) || (!pMesh->elem[pMesh->side[s].eb].on)) 
			{
				
				/*sideflux=pLayer->GetFluxThruFace(z1,z2,t).real();
				wallflux[n1]-=0.5*sideflux;
				wallflux[n2]-=0.5*sideflux;*/
	  		
				wallflux[n1]-=L*0.5*(0.78867*Qn1+0.21132*Qn2);//0.78867 & 0.21132 are values of basis function at gauss points
				wallflux[n2]-=L*0.5*(0.21132*Qn1+0.78867*Qn2);

				wallflux[n1]=wallflux[n2]=0.0; //TMP DEBUG

				walllength[n1]+=0.5*abs(z1-z2);
				walllength[n2]+=0.5*abs(z1-z2);

				wallthick[n1]=wallthick[n2]=0.2;//TMP DEBUG
			}
		}

		//System Boundaries-------------------------------------------------------
		if (pMesh->side[s].status==BOUNDARY){
			if      ((pMesh->side[s].ea==OFF) || (!pMesh->elem[pMesh->side[s].ea].on)) 
			{
				sideflux=pLayer->GetFluxThruFace(z2,z1,t).real();
				borderflux[pMesh->side[s].n1]-=0.5*sideflux;
				borderflux[pMesh->side[s].n2]-=0.5*sideflux;
			}
			else if ((pMesh->side[s].eb==OFF) || (!pMesh->elem[pMesh->side[s].eb].on)) 
			{
				sideflux=pLayer->GetFluxThruFace(z1,z2,t).real();
				borderflux[pMesh->side[s].n1]-=0.5*sideflux;
				borderflux[pMesh->side[s].n2]-=0.5*sideflux;
			}
		}

	}
	//Singularities & Dry regions-------------------------------------------
	for (k=0;k<nNodes;k++){
		if (thick[k]*poro[k]<=0){
			sourceflux[k]=sinkflux[k]=borderflux[k]=wallflux[k]=rechargeflux[k]=0.0;
		}
		singularities[k]=pLayer->GetSingularStrength(pMesh->node[k].z,t).real();
	}

	delete [] count;


	//Create steady-state Global matrix [M] and sorption matrix [A]----------  
	cout << "  Creating Global Matrices...";

	/*ofstream FEMB;
	FEMB.open("fe_water_balance.csv");
  FEMB<<"nd,x,y,vx,vy"<<endl;
	FEMB.close();*/

	ofstream TMP_GAUSS;
	TMP_GAUSS.open("gauss_pts.bna");
	TMP_GAUSS.close();

	double avgarea(0.0);

	M =new CSparseMatrix(); M ->DynamicInitialize(nNodes); 
	Df=new CSparseMatrix(); Df->DynamicInitialize(nNodes);
	A =new CSparseMatrix(); A ->DynamicInitialize(nNodes);

	for (e=0; e<pMesh->NumElems; e++)
	{

		if (ProgramAborted()){return;}

		WriteEllipsisToScreen(e,pMesh->NumElems,20);

		avgarea+=pMesh->elem[e].area;

    loc_len=sqrt(2.0*pMesh->elem[e].area); //local length scale

		nde[0]=pMesh->elem[e].ni;
		nde[1]=pMesh->elem[e].nj;
		nde[2]=pMesh->elem[e].nk;

		GetDispersionMatrix(e,t,elementmatrix, diffusionmatrix,vloc,dloc);
		GetSorptionMatrix  (e,t,sorptionmatrix,0);

		for (i=0;i<3;i++){
			for (j=0;j<3; j++){

				tmpa=M->     GetAij(     nde[i],nde[j])+elementmatrix[i][j];//+lambda*sorptionmatrix[i][j];
					   M-> DynamicAdd(tmpa,nde[i],nde[j]);

				tmpa=Df->    GetAij(     nde[i],nde[j])+diffusionmatrix[i][j];
						 Df->DynamicAdd(tmpa,nde[i],nde[j]);

				tmpa=A->     GetAij(     nde[i],nde[j])+sorptionmatrix[i][j];
					   A-> DynamicAdd(tmpa,nde[i],nde[j]); 
			}
		}

		//Evaluate local Peclet Number 
		PeAvg+=loc_len*fabs(vloc)/dloc;
    upperswap(PeMax,loc_len*fabs(vloc)/dloc);
		peclet[e]=loc_len*fabs(vloc)/dloc;

		//Evaluate local Courant number (assume tstep=1)
		CrAvg+=fabs(vloc)/loc_len;
		upperswap(CrMax,fabs(vloc)/loc_len);
		VelAvg+=fabs(vloc);

	} /* end for (e=0; e<NumElems... */

	PeAvg  /=pMesh->NumElems;
	CrAvg  /=pMesh->NumElems;
  avgarea/=pMesh->NumElems;
  VelAvg /=pMesh->NumElems;

	last_loc_tstep=1.0;
	last_Diff     =0.0; 	

  double Diff,htheta,avehtheta,htheta2;
  int    i_dbl;

	//[M] created with equivalent time step of 1.0 and equivalent diffusion coeff of 0.0
	//([M]=[A]+w(tstep)([D]+D*[Df]+[S]));  [S]=[I]*[Qi/nh]; D*=0.0;
	for (i=0; i<nNodes; i++)
	{

		if (!theta_on_left){htheta=poro[i]*satthick[i];}
		else               {htheta=1.0;                }

		//Process Doublenode Boundaries-----------------------------------
		//i_dbl is index of node "linked to" node i across doublenode boundary
		i_dbl=-1;
		for (int l=0; l<pMesh->NumDoublenodes; l++){
			if (pMesh->dnodes[l].n1==i){i_dbl=pMesh->dnodes[l].n2;break;}
			if (pMesh->dnodes[l].n2==i){i_dbl=pMesh->dnodes[l].n1;break;}
		}
		//-----------------------------------------------------------

		for (j=0; j<nNodes; j++){

			//Process Conventional Advection/Dispersion Terms----------------
			tmpa =A->GetAij(i,j);
			tmpa+=(w)*last_loc_tstep*(M->GetAij(i,j));

			//Process Source & Sink Flux Terms-------------------------------
			if (i==j){
				tmpa+=(w)*last_loc_tstep*(sinkflux    [i]/htheta); //sinkflux=negative
				tmpa+=(w)*last_loc_tstep*(sourceflux  [i]/htheta); //sourceflux=positive
				tmpa+=(w)*last_loc_tstep*(rechargeflux[i]/htheta); //rechargeflux=positive for recharge, neg for leakage
									 
				if ((borderflux[i]>0.0) && (process_type!=DISPERSION_ONLY)){
					tmpa+=(w)*last_loc_tstep*(borderflux[i]/htheta);
				}
			}

			//Process Discontinuous Internal Boundary condition--------------
			// here, Fi= -QnC+ + hnD (C- - C+)/tw

      if (i_dbl!=-1){
				avehtheta=0.25*(poro[i]+poro[i_dbl])*(thick[i]+thick[i_dbl]);
				if (!theta_on_left){htheta2=poro[i_dbl]*satthick[i_dbl];}
				else               {htheta2=1.0;                    }  

				Diff=pDomain->GetDispersivity(0.0,LONGITUDINAL)*fabs(wallflux[i])/walllength[i]/avehtheta+pDomain->GetDiffusionCoeff(0);//TMP DEBUG - species 0
		
				//Advective flux into domain from wall (source flux)
				if (wallflux[i]>=0){//
					if      (j==i_dbl){/*tmpa-=(w)*last_loc_tstep*wallflux[i]/htheta;*/ }//Advective flux only-presently included on RHS
					else if (j==i    ){tmpa+=(w)*last_loc_tstep*wallflux[i]/htheta;}  //Advective flux only-functional (other component on RHS)
				}

				//Dispersive flux into domain from wall (source or sink flux)		
				/*if      (wallflux[i]>=0){
	       //tmpa+=(w)*last_loc_tstep*htheta*Diff*htheta*(Cnew[j][s]-Cnew[tmp][s])/wallthick[j];//Purely Explicit formulation
					if      (j==i_dbl){tmpa-=(w)*last_loc_tstep*htheta2*Diff/wallthick[j]*walllength[j];}
					else if (j==i    ){tmpa+=(w)*last_loc_tstep*htheta *Diff/wallthick[j]*walllength[j];}
				}
				else if (wallflux[i]<0){ //dispersive flux into wall from domain	

					if      (j==i_dbl){tmpa-=(w)*last_loc_tstep*htheta2*Diff/wallthick[j]*walllength[j];}	
					else if (j==i    ){tmpa+=(w)*last_loc_tstep*htheta *Diff/wallthick[j]*walllength[j];}
				}*/
			}

			M-> DynamicAdd(tmpa,i,j);

		}
	} /* end for (i=0; i<NumNodes... */

	//Identify dirichlet conditions----------------------------------------------------
	int     numDir =0;
	int    *dir    =new int    [nNodes];
	double *dirconc=new double [nNodes];
	for (k=0; k<nNodes; k++){
		dirconc  [k]=pDomain->GetDirichletConc(pMesh->GetNodeLocation(k),t,0);
		dirichlet[k]=(dirconc[k]!=NO_INFLUENCE);
		if (dirichlet[k]){dir[numDir]=k; numDir++;}
	}
	M->MakeDirichlet(dir,numDir); //Alternative method (may be better)-ignores rows
	delete [] dir;
	delete [] dirconc;
	/*----------------------------------------------------------------------------------
	   Mesh Diagnostics
	----------------------------------------------------------------------------------*/
	double netsink(0),netsrce(0),netrech(0),netwall(0),netbrdr(0),netbrd2(0);
	int numsing(0);
	for (i=0; i<pMesh->NumNodes; i++){
		netsink+=sinkflux    [i];
		netsrce+=sourceflux  [i];
		netrech+=rechargeflux[i];
		netwall+=fabs(wallflux[i])/2.0;
		netbrdr+=borderflux[i]; 
		netbrd2+=fabs(borderflux[i])/2.0;
		if (singularities[i]!=0.0){numsing++;}
	}
	cout <<endl;
	cout <<"   -LUMPING                            " << ((lumped) ? "ON" : "OFF") <<endl;
	cout <<"   -# GAUSS INTEGRATION PTS:           " << GetNumTriGaussPts(gausstype) <<endl;
	cout <<"   -NUMBER OF DIRICHLET NODES:         " << numDir <<endl;	 
	cout <<"   -AVERAGE ELEMENT AREA:              " << avgarea <<endl;
	cout <<"   -AVERAGE REPRESENTATIVE LENGTH:     " << sqrt(2.0*avgarea)<<endl;
	cout <<"   -NET SINK FLUX:                     " << netsink<<endl;
	cout <<"   -NET SOURCE FLUX:                   " << netsrce<<endl;
	cout <<"   -NUMBER OF SINGULARITIES            " << numsing<<endl;
	cout <<"   -NET RECHARGE FLUX:                 " << netrech<<endl;
	cout <<"   -GROSS DISCONTINUITY FLUX:          " << netwall<<endl;
	cout <<"   -NET FLUX THRU MESH:                " << netbrdr<<endl;
	cout <<"   -GROSS FLUX THRU MESH:              " << netbrd2<<endl;
  cout <<"   -MAXIMUM PECLET NUMBER:             " << PeMax <<endl;
  cout <<"   -AVERAGE PECLET NUMBER:             " << PeAvg <<endl;
  cout <<"   -SUGGESTED TIME STEP (Cr(max)=1):   " << ((process_type==ADV_AND_DISP) ? 1.0/CrMax : sugg_tstep)  <<endl;
  cout <<"   -PERMISSIBLE TIME STEP (Cr(avg)=1): " << ((process_type==ADV_AND_DISP) ? 1.0/CrAvg : sugg_tstep)  <<endl;
	cout <<"   -MAX MATRIX COEFF:                  " << M->GetMaxCoeff()<<endl;
	cout <<"   -AVG DIAGONAL COEFF:                " << M->GetAvgDiagCoeff()<<endl;	
	cout <<"...Done initializing Finite Element Transport Scheme..."<<endl;

	StraightSort(peclet,pMesh->NumElems);

  if (print_diagnostics){
		ofstream DIAG;
		DIAG.open("diagnostics.csv");
		DIAG<<"Number_Of_Elements,"       <<pMesh->NumElems<<endl;
		DIAG<<"Number_Of_Nodes,"          <<nNodes         <<endl;
		DIAG<<"Number_Of_Dirichlet_Nodes,"<<numDir				 <<endl;
		DIAG<<"Average_Cell_Area,"        <<avgarea				 <<endl;
		DIAG<<"Average_Represent_Length," <<sqrt(avgarea)	 <<endl;
		DIAG<<"Max_Peclet_Number,"        <<PeMax					 <<endl;
		DIAG<<"Mean_Peclet_Number,"       <<PeAvg					 <<endl;
		DIAG<<"90%_Peclet_Number,"        <<peclet[(int)(0.90*pMesh->NumElems)]<<endl;
		DIAG<<"80%_Peclet_Number,"        <<peclet[(int)(0.80*pMesh->NumElems)]<<endl;
		DIAG<<"75%_Peclet_Number,"        <<peclet[(int)(0.75*pMesh->NumElems)]<<endl;
		DIAG<<"50%_Peclet_Number,"        <<peclet[(int)(0.50*pMesh->NumElems)]<<endl;
		DIAG<<"25%_Peclet_Number,"        <<peclet[(int)(0.25*pMesh->NumElems)]<<endl;
    DIAG<<"Courant_Time_Step,"			  <<1.0/CrMax      <<endl;
		DIAG.close();
	}

	delete [] peclet;


	//Apply Courant constraints (Cr<=1) (Huyakorn & Pinder [1983], pg. 206)------
	if (process_type==ADV_AND_DISP){lowerswap(sugg_tstep,1.0/CrMax);}

	//Check for a poor matrix (usually peclet number is too high)----------------
	for (i=0; i<nNodes; i++){
		for (j=0; j<nNodes; j++){
			tmpa =M->GetAij(i,j);
			//ExitGracefullyIf((fabs(tmpa)>100000.0),"C2DFEEulerian::Initialize: Ill-conditioned matrix",BAD_DATA);
		}
	}
	//M->Print(false,50);
	//sugg_tstep=5.0; //Overrides suggested time step	
  //TMP.close();//TMP DEBUG
}

/***************************************************************************
							TRANSPORT
****************************************************************************
Eulerian Finite Difference Advection/Dispersion on rectangular grid
---------------------------------------------------------------------------*/
void C2DFEEulerian::Transport(const double     t, 
											        const double     tstep,								   
											        Ironclad2DArray  C,      //Concentration at start of tstep
											        Ironclad2DArray  Cend,   //Concentration from simultaneous processes at end of tstep
											        //Ironclad2DArray  dCdt_sim, //average concentration change from simultaneous processes over time step
															Writeable2DArray Cnew){  //New concentration after advection/dispersion/simultaneous
  
	static double R;
	static double local_t,loc_tstep;
	static double tmpcoeff,tmp_err;
	static cmplex tmpz;
	static pt3D	  pt;
  static double Diff,tmpA; 
	static double htheta;
	int		        i,j,k,s,j_dbl;

	tmp_err=0.0;
	R=1.0;
	
  //Copy aqueous concentration array from start of time step
	for (k=0;k<nNodes;k++){
		for (s=0; s<nspecies; s++){Cnew[k][s]=C[k][s];}} 
	
	local_t=t;

	while (local_t<t+tstep){

		loc_tstep=min(min(sugg_tstep,tstep/ceil(tstep/sugg_tstep)),t+tstep-local_t);

		local_t +=loc_tstep;
		ExitGracefullyIf((loc_tstep<=0.0),
			"C2DFEEulerian::Transport: negative or zero local time step",RUNTIME_ERR);

		//if (ProgramAborted()){return;}

		/*slowly introduce influence of paralell processes (simultaneous advective/reactive transport)
		  dispersive flux is based upon weighted average of concentration at beginning of local time step and end of local time step
		  if tstep=loc_tstep, weighting is just 0.5 (as in Konikow, 1996) */
		for (k=0;k<nNodes;k++){
			for (s=0; s<nspecies; s++){Cnew[k][s]+=0.5*(loc_tstep/tstep)*(Cend[k][s]-C[k][s]);}}

		//Solve for change in concentration due to transport over local time step
		for (s=0; s<nspecies; s++){

			/*Update Master [M] matrix (it changes with local time step and species :( )-----
			  [M]=[A]+dt*w*([D]+[Df]+[S])
       ([D]+[Df]+[S])=([Mlast]-[A])/dtlast/w
			  [M]=[A]+dt/dtlast*([Mlast]-[A])*/

			Diff=pDomain->GetDiffusionCoeff(s);

      //change matrix for impact of new time step============================================
			if (loc_tstep!=last_loc_tstep)
			{ 
				for (k=0; k<M->GetNumEntries(); k++){
					if (!M->K_to_IJ(k,i,j))      {
						ExitGracefully("C2DFEEulerian::Transport: Bad Sparse Matrix modify(1)",RUNTIME_ERR);}	

					  tmpA=A->GetAij(i,j);

						tmpcoeff=tmpA+loc_tstep/last_loc_tstep*(M->GetAij(i,j)-tmpA);

					if (!M->SetAij(tmpcoeff,i,j)){
						ExitGracefully("C2DFEEulerian::Transport: Bad Sparse Matrix modify(2)",RUNTIME_ERR);}
				}

        last_loc_tstep=loc_tstep;
			}

      //change matrix for impact of species-specific Diffusion Coefficient=====================
			if (last_Diff!=Diff)
			{ 
				for (k=0; k<M->GetNumEntries(); k++){
					if (!M->K_to_IJ(k,i,j))      {
						ExitGracefully("C2DFEEulerian::Transport: Bad Sparse Matrix modify(1b)",RUNTIME_ERR);}	

					tmpcoeff=M->GetAij(i,j)+(w)*loc_tstep*(Diff-last_Diff)*Df->GetAij(i,j);

					if (!M->SetAij(tmpcoeff,i,j)){
						ExitGracefully("C2DFEEulerian::Transport: Bad Sparse Matrix modify(2b)",RUNTIME_ERR);}
				}
        last_Diff=Diff;			
			}


			//Update Master B vector (also changes with local time step and species)=================
			//{B}=([A]-(1-w)*dt*([D]+[S])){Cn} + dt*{Fn})
			//where ([D]+[S])=([M]-[A])/w/dt ; Substituting,
			//{B}=([A]-(1-w)/w*([M]-[A])){Cn} + dt*{Fn})
			for (j=0; j<nNodes; j++){
				pt = pMesh->GetNodeLocation(j);
				B[j]=0.0;

				if (dirichlet[j]){
					B[j]=(    w)*pDomain->GetDirichletConc(pt,local_t-loc_tstep,s)+
						   (1.0-w)*pDomain->GetDirichletConc(pt,local_t,s);
				}
				else{

					if (!theta_on_left){htheta=poro[j]*thick[j];}
					else               {htheta=1.0;             }

					//Explicit Advective/Dispersive/Decay Terms-------------------------
					for (i=0; i<nNodes; i++){
            B[j]+=( A->GetAij(j,i)-(1.0-w)/(w)*( M->GetAij(j,i)-A->GetAij(j,i)) )*Cnew[i][s];
					}

					//Source terms (wet sources)-----------------------------------------
					if ((sourceflux[j]>0.0) || (rechargeflux[j]>0)) {//TMP DEBUG
					  //B[j]+=loc_tstep*sourceflux  [j]*((1-w)*pDomain->GetSourceConcentration(j ,local_t-loc_tstep,s)+
						//                                 (  w)*pDomain->GetSourceConcentration(j ,local_t          ,s))/htheta; 
					  B[j]+=loc_tstep*rechargeflux[j]*((1-w)*pDomain->GetRechargeConcentration(pt ,local_t-loc_tstep,s)+
						                                 (  w)*pDomain->GetRechargeConcentration(pt ,local_t          ,s))/htheta; 
						B[j]+=loc_tstep*sourceflux  [j]*(0.0)/htheta;//-Cnew[j][s] (on LHS) //TMP DEBUG -clean water only
					}

					//Source Terms (dry sources)-----------------------------------------
					B[j]  +=loc_tstep*              ((1-w)*pDomain->GetSpecifiedMassFlux  (pt,local_t-loc_tstep,s)+  //TMP DEBUG - this returns mass flux per unit area-need to multiply by area (elem area/ne?)
						                               (  w)*pDomain->GetSpecifiedMassFlux  (pt,local_t          ,s));
					
					//Simultaneous processes (i.e., RxN)---------------------------------
          //B[j]  +=loc_tstep*dCdt_sim[j][s];  
					
					//Sink Terms (wet sinks)---------------------------------------------
					if (sinkflux[j]>0.0){
						B[j]+=loc_tstep*sinkflux[j]*Cnew[j][s]/htheta;//-Cnew[j][s] (on LHS)
					}

					//Decay Terms -------------------------------------------------------
					B[j]+=loc_tstep*pDomain->GetSpeciesDecay(s)/htheta;
					
					//Upstream border term (natural condition- incoming ambient conditions)---------
					if (process_type!=DISPERSION_ONLY){
					  B[j]+=loc_tstep*max(borderflux[j],0.0)*pDomain->GetAmbientConc(s)/htheta;//-Cnew[j][s] (on LHS)
					}

					//Leaky wall terms--------------------------------------------------- 
					j_dbl=-1;
					for (int l=0; l<pMesh->NumDoublenodes; l++){
						if (pMesh->dnodes[l].n1==j){j_dbl=pMesh->dnodes[l].n2;break;}
						if (pMesh->dnodes[l].n2==j){j_dbl=pMesh->dnodes[l].n1;break;}
					}
					if ((wallflux[j]>=0.0) && (j_dbl!=-1)){ //flux into domain
						
						B[j]+=loc_tstep*wallflux[j]*Cnew[j_dbl][s]/htheta;//Advective flux only-other part on LHS
						//B[j]-=loc_tstep*wallflux[j]*Cnew[j][s]/htheta;//now on LHS

					}
					
				}
			}

			//Copy old concentrations to temp vector (good initial guess)------------------------
			for (k=0; k<nNodes; k++){tmpC[k]=Cnew[k][s];}

      //Solve system of equations------------------------------------------------------
			M->BCG(B,tmpC,nNodes,RELATIVE_RESIDUAL,tmp_err, 1e-10);//Very sensitive to tolerance parameter
		  
			//Update concentration values-----------------------------------------------------
			for (k=0; k<nNodes; k++){Cnew[k][s]=tmpC[k];}

		} //end for s=0 to nspecies

    //Add remainder of non-dispersive portion for local tstep (simultaneous process)
		for (k=0;k<nNodes;k++){
			for (s=0; s<nspecies; s++){
			  Cnew[k][s]+=0.5*(loc_tstep/tstep)*(Cend[k][s]-C[k][s]); 
			}
		}
	}/* end while */
}
/***************************************************************************
							WriteOutput
****************************************************************************
---------------------------------------------------------------------------*/															
void C2DFEEulerian::WriteOutput(const int outputstep){
}

/***************************************************************************
							GET DISPERSION MATRIX
****************************************************************************
Calculates Element Dispersion Matrix for element e
Calculates Element Diffusion Matrix for element e
Adapted significantly from Istok, 1989, pg. 93
also returns local representatic velocity (vloc) and dispersion coeff (dloc)
---------------------------------------------------------------------------*/
void C2DFEEulerian::GetDispersionMatrix(const int e, 
																				const double &t, 
																				double De[3][3], 
																				double Df[3][3], 
																				double &vloc, 
																				double &dloc) const{

	static cmplex    z   [3],zcen;
	static double    dNdx[3][MAXGAUSSPTS];
	static double    dNdy[3][MAXGAUSSPTS];
	static disp_info disp;
	static int       i,k,n1,n2,n3;
	static int       nGaussPts,n;
	static cmplex    zg      [MAXGAUSSPTS];
	static double    GaussWts[MAXGAUSSPTS];
	static double    Dxx     [MAXGAUSSPTS]; 
	static double    Dxy     [MAXGAUSSPTS];
	static double    Dyy     [MAXGAUSSPTS];
	static double    Dyx     [MAXGAUSSPTS];
	static cmplex    v       [MAXGAUSSPTS];
	static double    W       [MAXGAUSSPTS][3];
	static double    TxTy    [MAXGAUSSPTS];
	static vector    vel     [MAXGAUSSPTS];
	static double    htheta  [MAXGAUSSPTS];
	static double    dvxdx	 [MAXGAUSSPTS];
	static double    dvydy	 [MAXGAUSSPTS];
	static double    Jdet    [MAXGAUSSPTS];
	static double    alpha[3],xi[3];
  static double    v0;
  static double    Pe;
	static vector    v_tradit; 
	static double    elem_tau;

	bool use_nodal(false);
	static vector    vnode[3];

  for (i=0; i<3; i++){
		for (k=0; k<3; k++){De[k][i]=Df[k][i]=0.0;}}

	//offset sligtly inward for distinct properties for each e,n combo
	zcen=pMesh->elem[e].centroid;

	n1=pMesh->elem[e].ni; z[0] =pMesh->node[n1].z;  z[0] +=0.001*(zcen-z[0]);
	n2=pMesh->elem[e].nj; z[1] =pMesh->node[n2].z;  z[1] +=0.001*(zcen-z[1]);
	n3=pMesh->elem[e].nk; z[2] =pMesh->node[n3].z;  z[2] +=0.001*(zcen-z[2]);

  //Get locations of gauss points------------------------------------------- 
  nGaussPts=GetNumTriGaussPts(gausstype);
	GetTriGaussInfo(cmplex(0.001,0.001),
									cmplex(0.001,0.999),
									cmplex(0.999,0.001),gausstype,zg,GaussWts);//gets zg in local coords
  for (n=0;n<nGaussPts; n++){
		zg[n]=pMesh->elemLocalToGlobal(e,zg[n].real(),zg[n].imag(),1-zg[n].real()-zg[n].imag());
	}
	
	if ((pMesh->side[pMesh->elem[e].si].radius!=0.0) ||  
		  (pMesh->side[pMesh->elem[e].sj].radius!=0.0) ||  
			(pMesh->side[pMesh->elem[e].sk].radius!=0.0) ){ 
		ofstream TMP_GAUSS;
		TMP_GAUSS.open("gauss_pts.bna",ios::app);
		for (n=0;n<nGaussPts; n++){
				TMP_GAUSS << "\" e \",  1" <<endl;
				TMP_GAUSS << zg[n].real()<< " , "    <<zg[n].imag()<<endl;
		}
		TMP_GAUSS.close();
	}

	/*for (n=0;n<nGaussPts; n++){
    pMesh->elemGlobalToLocal(e,zg[n],xi[0],xi[1],xi[2]);
		cout.precision(10);
		cout <<z[0]<<","<<z[1]<<","<<z[2]<<","<<zg[n]<<","<<xi[0]<<","<<xi[1]<<","<<xi[2]<<endl;
	}*/

	//Evaluate upstream weights along element sides---------------------------
	//------------------------------------------------------------------------
	static double up_weight;
	up_weight=upwt;
	if (upwt>0.0){up_weight=1.0;}

	//weighting based on average tangential flux->just need direction, not magnitude
  v0=pLayer->GetFluxThruFace(z[1],z[2],t).imag(); if (v0>0){alpha[0]=up_weight;}else if (v0<=0.0){alpha[0]=0.0;} 
	v0=pLayer->GetFluxThruFace(z[2],z[0],t).imag(); if (v0>0){alpha[1]=up_weight;}else if (v0<=0.0){alpha[1]=0.0;}  
	v0=pLayer->GetFluxThruFace(z[0],z[1],t).imag(); if (v0>0){alpha[2]=up_weight;}else if (v0<=0.0){alpha[2]=0.0;} 

	//Calculate velocity in traditional manner--------------------------------
	//------------------------------------------------------------------------
	if (use_traditional)
	{//velocity obtained from head at element nodes- Only works for linear basis function

		double h[3],K;cmplex dhd;
		K   =pLayer->GetCond(zcen);
    h[0]=pLayer->GetHead(z[0],t); h[1]=pLayer->GetHead(z[1],t); h[2]=pLayer->GetHead(z[2],t);

    dhd= pMesh->GetdN(e,zcen,0)*h[0]+pMesh->GetdN(e,zcen,1)*h[1]+pMesh->GetdN(e,zcen,2)*h[2];

		v_tradit=c2Dto3D(-K*dhd);
	}
	if (use_nodal){
		vnode[0]=pLayer->GetVelocity3D(c2Dto3D(z[0]),t);
		vnode[1]=pLayer->GetVelocity3D(c2Dto3D(z[1]),t);
		vnode[2]=pLayer->GetVelocity3D(c2Dto3D(z[2]),t);
	}

	static cmplex dN(0);
	//Evaluate velocity, dispersion coeff, and weights at gauss points---------
	//------------------------------------------------------------------------
	for (n=0;n<nGaussPts; n++){
		
		pMesh->elemGlobalToLocal(e,zg[n],xi[0],xi[1],xi[2]); 

		//if (pMesh->GetElement(zg[n])!=e){ //NO LONGER VALID WITH CURVATURE???
			//cout << "n:"<<n;ExitGracefully("GetFEDispersionMatrix: bad gauss coordinates",RUNTIME_ERR);}

		for (i=0; i<3; i++){
			dN=pMesh->GetdN(e,zg[n],i);
			dNdx[i][n]=dN.real();
			dNdy[i][n]=dN.imag();
		}

		Jdet[n]=pMesh->GetElemJacobianDet(e,zg[n]);

		if      (use_traditional)
		{																		//Single element parameter (from head at nodes)
			vel   [n]=v_tradit; 
			htheta[n]=pLayer ->GetPoro(zcen)*pLayer ->GetSaturatedThickness(zcen,t);
		}
		else if (use_nodal) 
		{                                  //Nodal parameters: works only for linear basis functions
			vel   [n]=CVector(
				        xi[0]*vnode[0].x+xi[1]*vnode[1].x+xi[2]*vnode[2].x,
								xi[0]*vnode[0].y+xi[1]*vnode[1].y+xi[2]*vnode[2].y,0.0);
	  	htheta[n]=xi[0]*poro[n1]*thick[n1]+
				        xi[1]*poro[n2]*thick[n2]+
								xi[2]*poro[n3]*thick[n3];		
		}
		else                
		{																	 //Continuous parameters 
			vel		[n]=pLayer ->GetVelocity3D(c2Dto3D(zg[n]),t);       
			htheta[n]=pLayer ->GetPoro(zg[n])*pLayer ->GetSaturatedThickness(zg[n],t);
		}

		pDomain->GetDispersivities(c2Dto3D(zg[n]),0,disp); //species doesnt matter-unit diffusion coefficient evaluated


		v          [n]=c3Dto2D(vel[n]);
		g_htheta[e][n]=htheta[n];

    Dxx[n]=pDomain->CalcDispersionCoeff(vel[n],OR_XX,1.0);
    Dxy[n]=pDomain->CalcDispersionCoeff(vel[n],OR_XY,1.0);
    Dyy[n]=pDomain->CalcDispersionCoeff(vel[n],OR_YY,1.0);
    Dyx[n]=pDomain->CalcDispersionCoeff(vel[n],OR_YX,1.0);

		//Subtract out influence of singularities

		if (theta_on_left){
			v    [n]*=htheta[n];//v represents Q if theta is on left
			Dxx  [n]*=htheta[n];
			Dxy  [n]*=htheta[n];
			Dyy  [n]*=htheta[n];
			Dyx  [n]*=htheta[n];
			dvxdx[n]*=htheta[n];
			dvydy[n]*=htheta[n];
		}

		if (turn_off_cross_terms){Dxy[n]=Dyx[n]=0.0;}

		W  [n][0]=xi[0]*(1.0-3.0*(alpha[2]*xi[1]-alpha[1]*xi[2]));
		W  [n][1]=xi[1]*(1.0-3.0*(alpha[0]*xi[2]-alpha[2]*xi[0]));
		W  [n][2]=xi[2]*(1.0-3.0*(alpha[1]*xi[0]-alpha[0]*xi[1]));

		if (SUPG){
			elem_tau=up_weight*0.5*sqrt(pMesh->elem[e].area)/abs(v[n]);//From Cirpka notes, chap 5.2

			for (i=0; i<3; i++){
				W[n][i]=xi[i]+elem_tau*(v[n].real()*dNdx[i][n]+
					                      v[n].imag()*dNdy[i][n]);		
			}
		}	
	}// end for (n=0;n<nGaussPts...

	//If element-averaged parameters are used (averaged obtained via integration)
  if ((false) && (!use_traditional)){ 
		double vxavg(0.0),vyavg(0.0),DxxAvg(0.0),DyyAvg(0.0),DxyAvg(0.0),DyxAvg(0.0);
		for (n=0;n<nGaussPts; n++){
			//should multiply by |J|/2A for curved elements
			vxavg +=GaussWts[n]*v[n].real();  vyavg +=GaussWts[n]*v[n].imag();
			DxxAvg+=GaussWts[n]*Dxx[n]; 			DyyAvg+=GaussWts[n]*Dyy[n];
			DyxAvg+=GaussWts[n]*Dyx[n]; 			DxyAvg+=GaussWts[n]*Dxy[n];
		}
		for (n=0;n<nGaussPts; n++){
			v[n]=cmplex(vxavg,vyavg);
			Dxx[n]=DxxAvg;										Dyy[n]=DyyAvg;
			Dyx[n]=DyxAvg;										Dxy[n]=DxyAvg;
		}
	}

	static double tmp;tmp=0.0;

	//==Assemble Element Matrices De & Df==================================
	for (i=0; i<3; i++){
		for (k=0; k<3; k++)
		{
			for (n=0;n<nGaussPts; n++)
			{
				//---De----------------------------------------------------------

				tmp=0.0;
				tmp+=dNdx[i][n]*(Dxx[n]*dNdx[k][n]+Dxy[n]*dNdy[k][n]);
				tmp+=dNdy[i][n]*(Dyx[n]*dNdx[k][n]+Dyy[n]*dNdy[k][n]);
				if (process_type!=DISPERSION_ONLY){
					tmp+=W[n][k]*(v[n].real()*dNdx[i][n]);
					tmp+=W[n][k]*(v[n].imag()*dNdy[i][n]);
					if (false){tmp+=W[n][k]*xi[i]*(dvxdx[n]+dvydy[n]);}//conservative formulation-NOT CORRECT
				}
				De[k][i]+=GaussWts[n]*Jdet[n]*tmp; 

				//---Df----------------------------------------------------------

				tmp=0.0;
				tmp+=dNdx[k][n]*(dNdx[i][n]); 
				tmp+=dNdy[k][n]*(dNdy[i][n]); 

				if (theta_on_left){tmp*=htheta[n];}

				//defined for unit diffusion coeff 
				Df[k][i]+=GaussWts[n]*Jdet[n]*tmp; 
				
			}
			//handle singular terms here..
			//sing[0]=singular[n1]; //place above
			//sing[1]=singular[n2];
			//sing[2]=singular[n3];
			/*if (fabs(sing[k])>0){
				//calculate thetaj, thetak
				if (k==i){
				  De[k][i]+=sing[k]*dNdx[i][n]*0.5*(cos(thetak)-cos(thetaj));
					De[k][i]-=sing[k]*dNdy[i][n]*0.5*(sin(thetak)-sin(thetaj));
				}
				else{

				}
			}*/
		}
	}
	/*ofstream FEMB;
	FEMB.open("fe_water_balance.csv",ios::app);
	cmplex tmp2(0.0);
	for (i=0; i<3; i++){
		for (n=0;n<nGaussPts; n++){
			tmp2+=W[n][i]*v[n]*GaussWts[n]*(pMesh->elem[e].area);
		}
		if (i==0){FEMB<<pMesh->elem[e].ni<<","<< z[i].real()<<","<< z[i].imag() <<","<<tmp2.real()<<","<<tmp2.imag()<<endl;}
		if (i==1){FEMB<<pMesh->elem[e].nj<<","<< z[i].real()<<","<< z[i].imag() <<","<<tmp2.real()<<","<<tmp2.imag()<<endl;}
		if (i==2){FEMB<<pMesh->elem[e].nk<<","<< z[i].real()<<","<< z[i].imag() <<","<<tmp2.real()<<","<<tmp2.imag()<<endl;}
	}
	FEMB.close();*/

	//==Evaluate local representative vel & D ==================
	vloc=0.0;dloc=0.0;
  for (n=0;n<nGaussPts; n++){
		//should multiply by |J|/2A for curved elements
		if (theta_on_left){
			vloc+=GaussWts[n]*abs(v[n])             /htheta[n];//representative local velocity
			dloc+=GaussWts[n]*max(0.3*Dxx[n],Dyy[n])/htheta[n];//representative local dispersion coeff
		}
		else {
			vloc+=GaussWts[n]*abs(v[n]);                       //representative local velocity
			dloc+=GaussWts[n]*max(0.3*Dxx[n],Dyy[n]);          //representative local dispersion coeff		
		}
	}


}


/***************************************************************************
							GET ELEMENT SORPTION MATRIX
****************************************************************************
Calculates Element Sorption Matrix Ae[][] for element e
Must call Element Dispersion Matrix first (for g_htheta)
---------------------------------------------------------------------------*/
void C2DFEEulerian::GetSorptionMatrix  (const int e, const double &t, double Ae[3][3], const int s) const{

	int    i,k;	

	static int    nGaussPts,n;
	static cmplex zg      [MAXGAUSSPTS];
	static double GaussWts[MAXGAUSSPTS];
	static double Rf      [MAXGAUSSPTS];
	static double hthetag [MAXGAUSSPTS];
	static double Jdet    [MAXGAUSSPTS];
	static double N       [MAXGAUSSPTS][3];
	static double xi[3];
	static cmplex z[3];
	static cmplex zcen;
	static double R_ave,htheta_ave;

	zcen=pMesh->elem[e].centroid;
	z[0] =pMesh->node[pMesh->elem[e].ni].z;  z[0] +=0.01*(zcen-z[0]);
	z[1] =pMesh->node[pMesh->elem[e].nj].z;  z[1] +=0.01*(zcen-z[1]);
	z[2] =pMesh->node[pMesh->elem[e].nk].z;  z[2] +=0.01*(zcen-z[2]);

  //Get locations of gauss points------------------------------------------- 
  nGaussPts=GetNumTriGaussPts(gausstype);
	GetTriGaussInfo(cmplex(0.001,0.001),
									cmplex(0.001,0.999),
									cmplex(0.999,0.001),gausstype,zg,GaussWts);//gets zg in local coords
  for (n=0;n<nGaussPts; n++){
		zg[n]=pMesh->elemLocalToGlobal(e,zg[n].real(),zg[n].imag(),1-zg[n].real()-zg[n].imag());
	}

	//Calculate Rf,theta*h,N,|J| at gauss points-------------------------------
	for (n=0;n<nGaussPts; n++){
		Rf[n] =pDomain->GetRetardation(c2Dto3D(zg[n]),t,s);

		Jdet[n]=pMesh->GetElemJacobianDet(e,zg[n]);

		if (use_traditional){
			hthetag[n] =pLayer->GetPoro(zcen) *pLayer->GetSaturatedThickness(zcen,t);
		}
		else{
			hthetag[n] =g_htheta[e][n];
		}
		if (!theta_on_left){hthetag[n]=1.0;}

		pMesh->elemGlobalToLocal(e,zg[n],xi[0],xi[1],xi[2]);

		N[n][0]=xi[0];
		N[n][1]=xi[1];
		N[n][2]=xi[2];
	}

	static double tmp;
	tmp=0.0;	

	//Consistent formulation-generally less well behaved (poor matrix conditioning) but more accurate
	for (i=0; i<3; i++){
		for (k=0; k<3; k++){
			Ae[k][i]=0.0;
			for (int n=0;n<nGaussPts; n++){
				tmp=N[n][k]*N[n][i]*Rf[n]*hthetag[n];
				//Ae[k][i]+=GaussWts[n]*tmp*2.0*(pMesh->elem[e].area); //Jacobian=2*A
				Ae[k][i]+=GaussWts[n]*Jdet[n]*tmp; //Jacobian=2*A
			}
		}
	}
	//Lumped formulation overwrites consistent formulation 
	if (lumped){ 
		if (true)//uses average value of R & htheta obtained via integration
		{
			R_ave=htheta_ave=0.0;
			for (int n=0;n<nGaussPts; n++){
				//should multiply by |J|/2A
				R_ave     +=GaussWts[n]*Rf[n];
				htheta_ave+=GaussWts[n]*hthetag[n]; 
			}	

			if (!theta_on_left){htheta_ave=1.0;}
	
			for (i=0; i<3; i++){
				for (k=0; k<3; k++){
					if (i==k){Ae[i][k]=R_ave*htheta_ave*2.0*(pMesh->elem[e].area)/3.0;}//Istok, pg. 94-not valid for linear triangles
					else     {Ae[i][k]=0.0;}
				}
			}
		}
		else	//Column-summed lumping: probably better (not yet tested)
		{
			double colsum;
			for (k=0; k<3; k++){
				colsum=0.0;
				for (i=0; i<3; i++){colsum+=Ae[i][k];Ae[i][k]=0.0;}
				Ae[k][k]=colsum;
			}
		}
		
	}
}
/***************************************************************************
							Calculate System Mass
****************************************************************************
---------------------------------------------------------------------------*/
double C2DFEEulerian::CalculateSystemMass(const int s, 
																					Ironclad2DArray C,
																					Ironclad2DArray Cs, 
																					double &SorbedMass) const{
	static double AqMass,area;
	static int    e,n1,n2,n3;

	static int    nGaussPts,n;
	static cmplex zg      [MAXGAUSSPTS];
	static double GaussWts[MAXGAUSSPTS];
	static double Jdet    [MAXGAUSSPTS];
	static double xi[3];
	static double Sg,Cg,pbg,hthetag,hthetacen,aveporo;
	static double t=0;//TMP DEBUG

	AqMass    =0.0;
	SorbedMass=0.0;
	double NegMass=0.0;

	for (e=0; e<pMesh->NumElems; e++){
		if (pMesh->elem[e].on){  //necessary?
			n1   =pMesh->elem[e].ni;
			n2   =pMesh->elem[e].nj;
			n3   =pMesh->elem[e].nk;

			area =pMesh->elem[e].area;

			//Integrated aqueous & sorbed mass (using Gauss Integration)

			//Get locations of gauss points------------------------------------------- 
			nGaussPts=GetNumTriGaussPts(gausstype);
			GetTriGaussInfo(cmplex(0.001,0.001),
											cmplex(0.001,0.999),
											cmplex(0.999,0.001),gausstype,zg,GaussWts);//gets zg in local coords
			for (n=0;n<nGaussPts; n++){
				zg[n]=pMesh->elemLocalToGlobal(e,zg[n].real(),zg[n].imag(),1-zg[n].real()-zg[n].imag());
			}

			//Evaluate |J|, h theta, C, Cs at Gauss points----------------------------
			pbg      =pDomain->GetBulkDryDensity(c2Dto3D(pMesh->elem[e].centroid));
      hthetacen=(thick[n1]*poro[n1]+thick[n2]*poro[n2]+thick[n3]*poro[n3])/3.0;
			aveporo  =(          poro[n1]+          poro[n2]+          poro[n3])/3.0;
			for (n=0;n<nGaussPts; n++){

				Jdet[n]=pMesh->GetElemJacobianDet(e,zg[n]);

				pMesh->elemGlobalToLocal(e,zg[n],xi[0],xi[1],xi[2]);

				if (!lumped){
					if (use_traditional){
						hthetag=(xi[0]*thick[n1]*poro[n1]+xi[1]*thick[n2]*poro[n2]+xi[2]*thick[n3]*poro[n3]);
					}
					else{
					  hthetag=g_htheta[e][n];
					}
				}
				else if (lumped){
				  hthetag=hthetacen;
				}

				Cg =(xi[0]* C[n1][s]+  xi[1]*C[n2][s]+  xi[2]*C[n3][s]); //interpolated concentration

				AqMass    +=GaussWts[n]*Jdet[n]*(Cg*hthetag)/2.0; 

				NegMass+=min(GaussWts[n]*Jdet[n]*(xi[0]*C[n1][s]*hthetag)/2.0,0.0);
				NegMass+=min(GaussWts[n]*Jdet[n]*(xi[1]*C[n2][s]*hthetag)/2.0,0.0);
				NegMass+=min(GaussWts[n]*Jdet[n]*(xi[2]*C[n3][s]*hthetag)/2.0,0.0);

				Sg =(xi[0]*Cs[n1][s]+  xi[1]*Cs[n2][s]+ xi[2]*Cs[n3][s])/aveporo*pbg;

				SorbedMass+=GaussWts[n]*Jdet[n]*(Sg*hthetag)/2.0; 
			}
		}
	}
	//cout <<"AqMass:"<<AqMass<<" Neg Mass:"<<NegMass<<" Pos Mass:"<<AqMass-NegMass<<endl;
	SorbedMass*=pDomain->Vratio();
	AqMass    *=pDomain->Vratio();
	return AqMass+SorbedMass;
}
//---------------------------------------------------------------------------
double C2DFEEulerian::CalculateBorderFlux(const int        s,
														              Ironclad2DArray  C,
																					Ironclad2DArray  Cend,
														              const double    &t,
																				  const double    &tstep) const{
  cmplex z1,z2;
	double sum(0.0);
	for (int k=0; k<pMesh->NumNodes; k++){
		if (!dirichlet[k]){
			//sum-=min(borderflux[k],0.0)*((1-upwt)*Cend[k][s]+upwt*C[k][s]);
			sum-=min(borderflux[k],0.0)*(0.5*Cend[k][s]+0.5*C[k][s]);
		}
	}

	return sum*tstep*pDomain->Vratio(); 
}
//---------------------------------------------------------------------------
double C2DFEEulerian::CalculateSourceGain(const int        s,
														              Ironclad2DArray  C,
																					Ironclad2DArray  Cend,
														              const double    &t,
																				  const double    &tstep) const{
	double sum(0.0);
	pt3D   pt;
	int    i,j,k;
	double Bj;
	double dirconcj,dirconck;
	double tmpsum0(0),tmpsum1(0),tmpsum2(0), tmpsum3(0),tmpsum4(0);

	//Dirichlet Source Fluxes----------------------------------------------------------
  for (j=0; j<nNodes; j++){
    dirflux[j]=0.0;
		if (dirichlet[j]){ //Dirichlet node (E.14)
			
			pt=pMesh->GetNodeLocation(j);
			dirconcj=pDomain->GetDirichletConc(pt,t+0.5*tstep,s);

			//Calculate RHS (if no dirichlet conditions existed)
      Bj=0.0;
			for (i=0; i<nNodes; i++){ //correct RHS
				Bj+=( A->GetAij(j,i)-(1.0-w)/(w)*( M->GetAij(j,i)-A->GetAij(j,i)) )*((1-w)*C[i][s]+(w)*Cend[i][s]); //Full 
			} 
			//cout <<"rhsj:"<<Bj<<endl;
			
			for ( i=0; i<nNodes; i++){
				dirflux[j]+=M->GetAij(j,i)*((1-w)*C[i][s]+(w)*Cend[i][s]);
				tmpsum0+=M->GetAij(j,i)*((1-w)*C[i][s]+(w)*Cend[i][s]);
			}
			dirflux[j]+=(dirconcj-((1-w)*C[j][s]+(w)*Cend[j][s]))-Bj; //Full
			//dirichlet influx/outflux boundaries
			//sum+=dirconc*borderflux[j];
      tmpsum3+=Bj;

			tmpsum4+=(dirconcj-((1-w)*C[j][s]+(w)*Cend[j][s]));
		}
		else{ //Non-Dirichlet node (E.13)
			for (k=0; k<nNodes; k++){
				if (dirichlet[k]){
					pt=pMesh->GetNodeLocation(k);
					dirconck=pDomain->GetDirichletConc(pt,t+0.5*tstep,s);
					dirflux[j]+=(M->GetAij(j,k)*dirconck)-(M->GetAij(k,j)*0.5*(C[k][s]+Cend[k][s]));
					tmpsum1+=M->GetAij(j,k)*dirconck;
					tmpsum2+=(M->GetAij(k,j)*0.5*(C[k][s]+Cend[k][s]));
				}
			}

		}
		//tmpsum4+=dirflux[j];

		sum+=dirflux[j]/last_loc_tstep ;
	}

	//Source fluxes----------------------------------------------------------------
  for (k=0; k<nNodes; k++){
		if (sourceflux  [k]>0.0){sum+=0.0*sourceflux[k];}	//TMP DEBUG-source flux has zero concentration
	}

	//Recharge fluxes--------------------------------------------------------------
  for (k=0; k<nNodes; k++){
		pt=pMesh->GetNodeLocation(k);
    if (rechargeflux[k]>0.0){
			sum+=((1-w)*pDomain->GetRechargeConcentration(pt ,t-tstep,s)+
						(  w)*pDomain->GetRechargeConcentration(pt ,t      ,s))*rechargeflux[k];
			 //cout <<sum<<"|";
		} 
	}

	//Ambient influx----------------------------------------------------------------
	for (k=0; k<nNodes; k++){
		sum+=max(borderflux[k],0.0)*pDomain->GetAmbientConc(s);
	}
	//Discontinuous fluxes (TMP DEBUG)---------------------------------------------
	/*double Diff,avehtheta;
	for (k=0; k<nNodes; k++){
		if (wallflux[k]>0){ //flux into domain
			int k_dbl;
			for (int l=0; l<pMesh->NumDoublenodes; l++){
				if (pMesh->dnodes[l].n1==k){k_dbl=pMesh->dnodes[l].n2;break;}
				if (pMesh->dnodes[l].n2==k){k_dbl=pMesh->dnodes[l].n1;break;}
			}
      //Advective Flux
			sum+=wallflux[k]*0.5*(C[k_dbl][s]+Cend[k_dbl][s]);

			avehtheta=0.25*(poro[k]+poro[k_dbl])*(thick[k]+thick[k_dbl]);

			Diff=pDomain->GetDispersivity(0.0,LONGITUDINAL)*fabs(wallflux[k])/walllength[k]/avehtheta+pDomain->GetDiffusionCoeff(s);
      
			//Dispersive Flux
			//sum+=avehtheta*Diff*0.5*((C[k_dbl][s]-C[k][s])+(Cend[k_dbl][s]-Cend[k][s]))/wallthick[k]*walllength[k];
		}
	}*/
	//cout << "total flux: "<<sum<<endl;
	
	return sum*tstep*pDomain->Vratio();
}
//---------------------------------------------------------------------------
double C2DFEEulerian::CalculateSinkLoss(const int        s,
														            Ironclad2DArray  C,
																				Ironclad2DArray  Cend,
														            const double    &t,
																			  const double    &tstep) const{
	double sum(0.0);
	int k;
  for (k=0; k<nNodes; k++){
		sum+=max(sinkflux[k],0.0)*((1-w)*C[k][s]+(w)*Cend[k][s]);
		//really, this operates on the integrated concentration over the adjacent elements
	}

	//Discontinuous fluxes (TMP DEBUG)---------------------------------------------
	/*double Diff,avehtheta;
  for (k=0; k<nNodes; k++){
    if (wallflux[k]<0){//flux out of domain
			int k_dbl=-1;
			for (int l=0; l<pMesh->NumDoublenodes; l++){
				if (pMesh->dnodes[l].n1==k){k_dbl=pMesh->dnodes[l].n2;break;}
				if (pMesh->dnodes[l].n2==k){k_dbl=pMesh->dnodes[l].n1;break;}
			}
			if (k_dbl!=-1){
				//Advective Flux
				sum-=wallflux[k]*((1-w)*C[k][s]+(w)*Cend[k][s]);

				avehtheta=0.25*(poro[k]+poro[k_dbl])*(thick[k]+thick[k_dbl]);

				Diff=pDomain->GetDispersivity(0.0,LONGITUDINAL)*fabs(wallflux[k])/walllength[k]/avehtheta+pDomain->GetDiffusionCoeff(s);
				
				//Dispersive Flux
			//	sum-=avehtheta*Diff*0.5*((C[k_dbl][s]-C[k][s])+(Cend[k_dbl][s]-Cend[k][s]))/wallthick[k]*walllength[k];
			}
		}
	}*/
	
	return sum*tstep*pDomain->Vratio();
}
//---------------------------------------------------------------------------
double C2DFEEulerian::CalculateDecayLoss(const int        s,
														             Ironclad2DArray  C,
																				 Ironclad2DArray  Cend,
														             const double    &t,
																			   const double    &tstep) const{
	double sum(0.0);
	int e;
	double lambda=pDomain->GetSpeciesDecay(s);
  for (e=0; e<pMesh->NumElems; e++){
		double area=pMesh->elem[e].area;
		sum-=lambda*0.5*0.33333*(C[pMesh->elem[e].ni][s]+Cend[pMesh->elem[e].ni][s])*area/poro[pMesh->elem[e].ni]/thick[pMesh->elem[e].ni];
		sum-=lambda*0.5*0.33333*(C[pMesh->elem[e].nj][s]+Cend[pMesh->elem[e].nj][s])*area/poro[pMesh->elem[e].ni]/thick[pMesh->elem[e].ni];
		sum-=lambda*0.5*0.33333*(C[pMesh->elem[e].nk][s]+Cend[pMesh->elem[e].nk][s])*area/poro[pMesh->elem[e].ni]/thick[pMesh->elem[e].ni];

	}
	
	return sum*tstep*pDomain->Vratio();
}

































	/*
	if (LinearParams) {divisor=pMesh->elem[e].area*pMesh->elem[e].area*8.0;}//Continuous Parameters
	else                  {divisor=1.0/pMesh->elem[e].area;                    }//Constant Parameters	
	if (!LinearParams){ 
		//vel =c2Dto3D(CalculateAverageVelocity(e,t));
		vel =c2Dto3D(pLayer->GetVelocity2D(zcen,t));
		R   =pDomain->GetRetardation(c2Dto3D(zcen),t,s);
		poro=pLayer->GetPoro(zcen);
		h   =pLayer->GetSaturatedThickness(zcen,t); 
	}
	double delta;
	int j;
	for (i=0; i<3; i++){
		for (k=0; k<3; k++){
			De[i][k]=0.0;

			if (LinearParams) {  
			  //Continuous Parameters---------------------------------------------------------
				for (j=0; j<3; j++){ 
					if (i==j){delta=1.0;}
					else     {delta=0.0;}

					vel =c2Dto3D(pLayer ->GetVelocity2D(z[j],t));
          R   =pDomain->GetRetardation(c2Dto3D(z[j]),t,s);
					poro=pLayer ->GetPoro(z[j]);
					h   =pLayer ->GetSaturatedThickness(z[j],t); 

					De[i][k]-=pDomain->GetDispersionCoeff(vel,OR_XX,1.0)/R*dNdx[i][n]*dNdx[j]*dNdx[k]/divisor;
					De[i][k]-=pDomain->GetDispersionCoeff(vel,OR_XY,1.0)/R*dNdx[i][n]*dNdx[j]*dNdy[k]/divisor;
					De[i][k]-=pDomain->GetDispersionCoeff(vel,OR_YY,1.0)/R*dNdy[i]*dNdy[j]*dNdy[k]/divisor;
					De[i][k]-=pDomain->GetDispersionCoeff(vel,OR_YX,1.0)/R*dNdy[i]*dNdy[j]*dNdx[k]/divisor;

					if (fabs(De[i][k])>10000.0){ExitGracefully("C2DFEEulerian::GetDispersionMatrix: Ill-conditioned dispersion matrix element (low retardation or small element area)",BAD_DATA);}
				
					if (process_type!=DISPERSION_ONLY){
            if (pDomain->UseEffectiveParams()){ vel=c2Dto3D(pLayer->GetEffVelocity(z[j],t,disp));}
					 	De[i][k]+=vel.x*dNdx[k]*(1+delta)/24.0/R;
						De[i][k]+=vel.y*dNdy[k]*(1+delta)/24.0/R;				
					} //end !DISPERSION_ONLY
				}//end for (j=0...
			}//end Continuous Params

			else{                    
				//Constant Parameters------------------------------------------------------------
				//Istok, pg. 93
        
				if (abs(vel)>1e8){ExitGracefully("C2DFEEulerian::GetDispersionMatrix: Ill-conditioned dispersion matrix element (bad velocity)",BAD_DATA);}

				De[i][k]+=pDomain->GetDispersionCoeff(vel,OR_XX,1.0)/poro/h/R*dNdx[i][n]*dNdx[k]*pMesh->elem[e].area;
				De[i][k]+=pDomain->GetDispersionCoeff(vel,OR_XY,1.0)/poro/h/R*dNdx[i][n]*dNdy[k]*pMesh->elem[e].area;
				De[i][k]+=pDomain->GetDispersionCoeff(vel,OR_YY,1.0)/poro/h/R*dNdy[i]*dNdy[k]*pMesh->elem[e].area;
				De[i][k]+=pDomain->GetDispersionCoeff(vel,OR_YX,1.0)/poro/h/R*dNdy[i]*dNdx[k]*pMesh->elem[e].area;

				if (fabs(De[i][k])>10000.0){ExitGracefully("C2DFEEulerian::GetDispersionMatrix: Ill-conditioned dispersion matrix element (low retardation or small element area)",BAD_DATA);}
				

				if (process_type!=DISPERSION_ONLY){
					//should use upstream-weighted petrov-galerkin
					//int SS[dW/dx v_x dNdx]
					if (pDomain->UseEffectiveParams()){ 
						pDomain->GetDispersivities(c2Dto3D(z[j]),s,disp);
						vel=c2Dto3D(pLayer->GetEffVelocity(zcen,t,disp));
					}
					De[i][k]+=vel.x*dNdx[k]/6.0/R; 
					De[i][k]+=vel.y*dNdy[k]/6.0/R;
				}
			}
		}
	}*/

	//TEMPORARY TESTING OF EFFECTIVE VELOCITY---------------------------------
	//Non-continuous velocity representation

	/*cmplex vavg(0.0);
	double dxxavg(0.0),dyyavg(0.0),dxyavg(0.0),dyxavg(0.0);
	cmplex vstddev(0.0);

	for (n=0;n<nGaussPts; n++){
	  vavg  +=GaussWts[n]*ve[n];
		dxxavg+=GaussWts[n]*Dxx[n];
		dxyavg+=GaussWts[n]*Dxy[n];
		dyyavg+=GaussWts[n]*Dyy[n];
		dyxavg+=GaussWts[n]*Dyx[n];
	}*/
	/*for (n=0;n<nGaussPts; n++){
    vstddev+=cmplex(pow(ve[n].real()-vavg.real(),2),pow(ve[n].imag()-vavg.imag(),2));

		//cout << "x*-"<<ve[n].real() <<" "<<vavg.real()<<"    "<<100.0*(ve[n].real()-vavg.real())/fabs(vavg.real())<<"%"<<endl;
		//cout << "y*-"<<ve[n].imag() <<" "<<vavg.imag()<<"    "<<100.0*(ve[n].imag()-vavg.imag())/fabs(vavg.imag())<<"%"<<endl;
	  ve  [n]=vavg;
		Dxx [n]=dxxavg;
		Dxy [n]=dxyavg;
		Dyy [n]=dyyavg;
		Dyx [n]=dyxavg;
		TxTy[n]=0.0;
	}*/
//cout<<"elem "<<e<<":  avg vel err:"<<100*sqrt(abs(vstddev))/max(abs(vavg),0.0000001)/nGaussPts<<endl;
	//perform gauss integration of 9 combinations (ii,ij,ik,ji...)------------
	//------------------------------------------------------------------------
		/*//obtain influence of unit diffusion coeff on eff. velocity (twice as much computation for multispecies)
		if (nspecies>1){ 
      disp.D=0.0;
      ve  [n]=pLayer ->GetEffVelocity2D(zg[n],t,disp);
			disp.D=1.0;
			TxTy[n]=(ve[n].real()-pLayer ->GetEffVelocity2D(zg[n],t,disp).real());
		}
		else{
      ve  [n]=pLayer ->GetEffVelocity2D(zg[n],t,disp);
			TxTy[n]=0.0;
		}*/
