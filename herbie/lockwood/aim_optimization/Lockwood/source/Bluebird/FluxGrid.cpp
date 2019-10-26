//FluxGrid.cpp
#include "FluxGrid.h"
/*-------------------------------------------------------------------
   Constructors
-------------------------------------------------------------------*/
CFluxGrid::CFluxGrid(){
	pGrid=NULL;
	vx=NULL;
	vy=NULL;
}
//-------------------------------------------------------------------
CFluxGrid::CFluxGrid(CRectGrid *grid){
	ExitGracefullyIf(grid==NULL,"CFluxGrid::Constructor: NULL input",BAD_DATA);
	pGrid=grid;

	vx = new double *[pGrid->nX+1];
	vy = new double *[pGrid->nX+1];
	ExitGracefullyIf(vy==NULL,"CFluxGrid::Constructor(1)",OUT_OF_MEMORY);
	for (int i=0; i<=pGrid->nX; i++){
		vx[i]=new double [pGrid->nY+1];
		vy[i]=new double [pGrid->nY+1];
		ExitGracefullyIf(vy[i]==NULL,"CFluxGrid::Constructor(2)",OUT_OF_MEMORY);
		for (int j=0; j<=pGrid->nY; j++){
			vx[i][j]=vy[i][j]=0.0;
		}
	}
}
//-------------------------------------------------------------------
CFluxGrid::~CFluxGrid(){
	for (int i=0; i<pGrid->nX; i++){
		delete [] vx[i];
		delete [] vy[i]; 
	}
	delete [] vx;
	delete [] vy;
}
/*******************************************************************************
           INITIALIZE
********************************************************************************

------------------------------------------------------------------------------*/
void  CFluxGrid:: Initialize(const CSingleLayerABC *pLayer){
	cmplex zc[4];
	cmplex midpt[4];
  int i,j;
	double hleft,hright,htop,hbot,dx,dy;
	double t(0.0);
	/*
						vy[i][j+1]
							 /\
							  |
						zc1---zc0
						 |     |
 vx[i][j]--> |     | -->vx[i+1][j]
						 |     |	
						zc2---zc3
							 /\
							  |
						vy[i][j]

	*/

  cout <<" Initializing Flux Grid..."<<endl;
	for (i=0; i<pGrid->nX; i++){
		for (j=0; j<pGrid->nY; j++){
			zc    [0]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j+1]));
			zc    [1]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j+1]));
			zc    [2]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->gridY[j  ]));
			zc    [3]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->gridY[j  ]));

			midpt [0]=pGrid->LocalToGlobal(cmplex(pGrid->nodeX[i  ],pGrid->gridY[j+1]));
			midpt [1]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i+1],pGrid->nodeY[j  ]));
			midpt [2]=pGrid->LocalToGlobal(cmplex(pGrid->nodeX[i  ],pGrid->gridY[j  ]));
			midpt [3]=pGrid->LocalToGlobal(cmplex(pGrid->gridX[i  ],pGrid->nodeY[j  ]));

			//Currently 1 point integration
			htop  =pLayer->GetPoro(midpt[0])*pLayer->GetSaturatedThickness(midpt[0],t);
			hright=pLayer->GetPoro(midpt[1])*pLayer->GetSaturatedThickness(midpt[1],t);
			hbot  =pLayer->GetPoro(midpt[2])*pLayer->GetSaturatedThickness(midpt[2],t);
			hleft =pLayer->GetPoro(midpt[3])*pLayer->GetSaturatedThickness(midpt[3],t);


			dx=abs(zc[2]-zc[1]);
			dy=abs(zc[1]-zc[0]);

			//for correct mass balance, must use integrated normal velocity
			// integrated tangential velocity is more accurate for dispersion cross-terms

			vy[i  ][j  ]= pLayer->GetFluxThruFace(zc[3],zc[2],t).real()/dx/hbot;
			vx[i  ][j  ]= pLayer->GetFluxThruFace(zc[2],zc[1],t).real()/dy/hleft;		
			if (j==pGrid->nY-1){
				vy[i  ][j+1]= pLayer->GetFluxThruFace(zc[0],zc[1],t).real()/dx/htop; 
			}
			if (i==pGrid->nX-1){
				vx[i+1][j  ]= pLayer->GetFluxThruFace(zc[3],zc[0],t).real()/dy/hright;
			}
			//cout << vx[i][j]<<" "<<vy[i][j]<<endl;
		}
	}
  cout <<" ...Done Initializing Flux Grid"<<endl;
}
/*******************************************************************************
           POLLOCK TRACK
********************************************************************************

------------------------------------------------------------------------------*/
pt3D CFluxGrid::PollockTrack (const pt3D &pt, const double &tstep,const bool forward){

	int i,j,k(OFF);
	bool exr,exl,eyt,eyb;
	double tloc,dx,dy,vxp,vyp,Ax,Ay,dtx,dty,dt;
	double mx, my;
	cmplex newpt;
	pGrid->GetCellIndex(pt,k);
	if (k==OFF){return pt;} //out of grid 
	bool captured=false;
  double relx,rely;
	double mult=1.0;
	if (!forward){mult=-1.0;}
	pGrid->k_to_ij(k,i,j);

	newpt =pGrid->GlobalToLocal(c3Dto2D(pt));
	
	tloc=0.0;

	do 
	{

		//cout <<"i: "<<i<<" j: "<<j<<endl;//" vx1: "<<vx[i][j]<<" vx2 "<<vx[i+1][j]<<" vy1 "<<vy[i][j]<<" vy2 "<<vy[i][j+1]<<endl;
		
		dx=pGrid->gridX[i+1]-pGrid->gridX[i];
		dy=pGrid->gridY[j+1]-pGrid->gridY[j];
		
		ExitGracefullyIf(((dx<=0)||(dy<=0.0)),"Bad grid spacing",RUNTIME_ERR);

		Ax=mult*(vx[i+1][j  ]-vx[i][j])/dx;
		Ay=mult*(vy[i  ][j+1]-vy[i][j])/dy;

		//relx=(newpt.real()-pGrid->gridX[i])/dx;
		//rely=(newpt.imag()-pGrid->gridY[j])/dy;
		//cout <<" relx: "<< relx<<" rely: "<<rely<<endl;

		vxp=Ax*(newpt.real()-pGrid->gridX[i])+mult*vx[i][j];
		vyp=Ay*(newpt.imag()-pGrid->gridY[j])+mult*vy[i][j];


		if ((vx[i+1][j  ]*vx[i][j])>0.0){exr=(mult*vx[i][j]>0); exl=!exr;}//same sign- flowthrough
		else if (vx[i][j]>0.0)          {exr=false;             exl= exr;}//sink
		else                            {exr=(vxp>0.0);         exl=!exr;}//source

		if ((vy[i  ][j+1]*vy[i][j])>0.0){eyt=(mult*vy[i][j]>0); eyb=!eyt;}//same sign- flowthrough
		else if (mult*vx[i][j]>0.0)     {eyt=false;             eyb= eyt;}//sink
		else                            {eyt=(vyp>0.0);         eyb=!eyt;}//source
		
		dtx=ALMOST_INF;
		dty=ALMOST_INF;
		if (exr){if (Ax!=0.0){dtx=log(mult*vx[i+1][j]/vxp)/Ax;}else{dtx=(pGrid->gridX[i+1]-newpt.real())/vxp;}}
		if (exl){if (Ax!=0.0){dtx=log(mult*vx[i  ][j]/vxp)/Ax;}else{dtx=(pGrid->gridX[i  ]-newpt.real())/vxp;}}
		if (eyt){if (Ay!=0.0){dty=log(mult*vy[i][j+1]/vyp)/Ay;}else{dty=(pGrid->gridY[j  ]-newpt.imag())/vyp;}}
		if (eyb){if (Ay!=0.0){dty=log(mult*vy[i][j  ]/vyp)/Ay;}else{dty=(pGrid->gridY[j  ]-newpt.imag())/vyp;}}

		dt=min(dtx,min(dty,(tstep-tloc)));

		//cout <<"vxp: "<<vxp<<" vyp: "<<vyp<<" dtx: "<<dtx<<" dty: "<<dty<<" dt: "<<dt<<endl;
		
		//cout <<"Ax: "<<Ax<<" Ay: "<<Ay<<endl;

		//cout <<"x: "<<newpt.real()<<" y: "<<newpt.imag()<<endl;

		if (Ax!=0.0){mx=(vxp*exp(Ax*dt)-vxp)/Ax;}
		else        {mx=vxp*dt;                      }

		if (Ay!=0.0){my=(vyp*exp(Ay*dt)-vyp)/Ay;}
		else        {my=vyp*dt;                      }

		newpt+=cmplex(mx,my);

		//cout <<"nx: "<<newpt.real()<<" ny: "<<newpt.imag()<<endl;
		relx=(newpt.real()-pGrid->gridX[i])/dx;
		rely=(newpt.imag()-pGrid->gridY[j])/dy;

		//cout <<" nrelx: "<< relx<<" nrely: "<<rely<<endl;

		//ExitGracefullyIf( (fabs(mx)>1.1*dx),"CFluxGrid::PollockTrack: bad movement(x)",RUNTIME_ERR);
		//ExitGracefullyIf( (fabs(my)>1.1*dy),"CFluxGrid::PollockTrack: bad movement(y)",RUNTIME_ERR);
		if ((fabs(mx)>1.1*dx) || (fabs(my)>1.1*dy)){captured=true;}

		if (dt==(tstep-tloc)){
		}
		else if (dt==dtx){//exits through right or left
			if (exr){i++;}else{i--;}
		}
		else if (dt==dty){//exits through top or bottom
			if (eyt){j++;}else{j--;}
		}
		else{
			ExitGracefully("CFluxGrid:PollockTrack: bad local tstep",RUNTIME_ERR);
		}
    tloc+=dt;
		
		k=pGrid->ij_to_k(i,j);

	} while ((tloc<tstep) && (k!=OFF) && (!captured));
	return c2Dto3D(pGrid->LocalToGlobal(newpt));
}
/*******************************************************************************
           WRITE GEOMETRY
********************************************************************************
------------------------------------------------------------------------------*/
void   CFluxGrid::WriteGeometry() const{
	pGrid->WriteGeometry();
}
