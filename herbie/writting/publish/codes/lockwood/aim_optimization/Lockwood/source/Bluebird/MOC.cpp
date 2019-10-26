//MOC.cpp
#include "CharacteristicMethods.h"
#include "MassParticle.h"
/************************************************************************
                           Constructors
*************************************************************************/
CMOC::CMOC(C2DDomainABC *pDom)
     :C2DTransportScheme(pDom){
	 pEulerian=NULL;
	 pMP=NULL;
	 nparts=NULL;
	 trackedpoints=NULL;
	 tracked=false;
   nCells=0;
	 distribution=RANDOM_PATTERN;//default
	 nparthigh=4;								 //default
	 npartlow=2;								 //default
}
//------------------------------------------------------------------------
CMOC::~CMOC(){
	if (globaldebug){cout << "  DESTROYING MMOC STYLE"<<endl;}
	if (pMP!=NULL){
		for (int s=0; s<nspecies; s++){delete pMP[s];}
		delete [] pMP;
	}
	if (trackedpoints!=NULL){
		for (int k=0; k<nCells; k++){
			for (int s=0; s<nspecies; s++){
				delete [] trackedpoints[k][s];
			}
			delete [] trackedpoints[k];
		}
		delete [] trackedpoints;
	}
}
/*************************************************************************
                           SetParameters
*************************************************************************/
void CMOC::SetParameters(distribution_type ty, const int nplow, const int nphigh){
	distribution=ty;
	nparthigh=nphigh;
	npartlow=nplow;
}
/*************************************************************************
                           Initialize
*************************************************************************/
void CMOC::Initialize(const double         startt,
											const transport_type ty,
									    Ironclad1DArray      porosity, 
									    Ironclad1DArray      satthick){
	int     k,s; 

	if (ty==DISPERSION_ONLY){
		ExitGracefully("CMMOC::Initialize: cannot simulate dispersion without advection",BAD_DATA);}

	process_type=ty;
	nspecies =pDomain->GetNumSpecies();	
	pConcGrid=pDomain->GetGrid();
	if (pConcGrid->GetType()!=CELL_BASED){
		ExitGracefully("Method of characteristics cannot yet be applied to non-cell-based grid or mesh",BAD_DATA);}
	pAquifer =pDomain->GetAquifer();
  nCells   =pConcGrid->GetNumCells();
	tracked  =false;

	//create particles
	pMP=new CMassParticle *[nspecies];
	for (s=0; s<nspecies; s++){
		pMP[s]=new CMassParticle(pAquifer,pDomain,0.0,TRACK_FORWARD,0.0,s);
	}
	nparts  =new int         [nCells];
	for (k=0; k<nCells; k++){nparts[k]=npartlow;}

	trackedpoints=new pt3D **[nCells];
	for (k=0; k<nCells; k++){
		trackedpoints[k]=new pt3D *[nspecies];
		for (s=0; s<nspecies; s++){
			trackedpoints[k][s]=new pt3D [nparthigh];
			if (trackedpoints[k][s]==NULL){ExitGracefully("CMOC::Initialize(): Out of memory- too many particles",OUT_OF_MEMORY);}
			for (int j=0; j<nparthigh; j++){
				trackedpoints[k][s][j]=0.0;
			}
		}
	}
}
/************************************************************************
                           TRANSPORT
*************************************************************************/
void CMOC::Transport(const double     t,
										 const double     tstep,
										 Ironclad2DArray  C,
										 Ironclad2DArray  Cend,
										 //Ironclad2DArray  dCdt_sim,
										 Writeable2DArray Cnew){
	int     j,k,s;	
	pt3D    *points=new pt3D[nparthigh]; 
	int     newcellID;
	int   **nLanded;

	//initialize concentrations, number of particles in cell (TMP DEBUG-should be member)
	nLanded =new int   *[nCells];
	for (k=0;k<nCells;k++){
		nLanded[k]=new int [nspecies];
		for (s=0;s<nspecies;s++){
			nLanded[k][s]=0;
		  Cnew   [k][s]=0.0;
		}                      
	}

  AllocateParticles(C);

	disp_info disp[MAX_SPECIES];
	for (s=0;s<nspecies;s++){
		pDomain->GetDispersivities(0.0,s,disp[s]); //location doesnt matter
	}

	//pre-identify particle locations,track particles forwards
	if (!tracked){
		cout <<endl<<"Initializing MOC characteristic paths";
		for (k=0;k<nCells;k++){
			WriteEllipsisToScreen(k,nCells,20);
			pConcGrid->DistributePoints(k,RANDOM_PATTERN,nparthigh,points);
			//pConcGrid->DistributePoints(i,RANDOM_PATTERN,npartlow,points);
			for (s=0;s<nspecies;s++){
				//for (k=0;k<nparts[i];k++){
				for (j=0;j<nparthigh;j++){
					pMP[s]->SetLocation(points[j]);
					pMP[s]->Track      (tstep,advectiontype,adv_timestep,adv_spacestep,pDomain->UseEffectiveParams(),disp[s],false); 
					trackedpoints[k][s][j]=pMP[s]->GetLocation();
				}
			}
		}
		cout <<endl;
		tracked=true;
	}

	for (k=0;k<nCells;k++){//for each cell...

		WriteEllipsisToScreen(k,nCells,20);
    //pConcGrid->DistributePoints(k,distribution,nparts[k],points);
    pConcGrid->DistributePoints(k,distribution,nparthigh,points);
		
		for (s=0;s<nspecies;s++){
			for (j=0;j<nparts[k];j++){//for each particle in cell,
			
				/*if (j>npartlow){
					pMP[s]->SetLocation      (points[j]);
					pMP[s]->Track            (tstep,advectiontype,adv_timestep,adv_spacestep); 
				}
				else{*/
					pMP[s]->SetLocation      (trackedpoints[k][s][j]);
				//}
				if (process_type==ADV_AND_DISP){
					pMP[s]->Disperse(tstep,t);
				}
				if (pConcGrid->GetCellIndex(pMP[s]->GetLocation(),newcellID)){
					nLanded[newcellID][s]++;
					//calculate average concentration of particles that have landed in the new cell
					//TMP DEBUG - should be volume weighted based on volume of old cell
					//C[k][s] could be interpolated from original position!?
					Cnew[newcellID][s]=(((double)(nLanded[newcellID][s]-1)*Cnew[newcellID][s])+pDomain->GetConcentration(points[j],t,s))/(double)(nLanded[newcellID][s]);									
					//Cnew[newcellID][s]=(((double)(nLanded[newcellID][s]-1)*Cnew[newcellID][s])+C[k][s])/(double)(nLanded[newcellID][s]);				
				}
			}
		}
	}
	for (k=0; k<nCells; k++){         //for all cells without particles at end of step, use previous concentration
		if (nLanded[k]==0){	
			for (s=0;s<nspecies;s++){
			  Cnew[k][s]=C[k][s];
			}
		}
	}

	/*if (pEulerian!=NULL){
		pEulerian->Transport(t,tstep,C,Cnew(updated),Cnew);
	}*/

	for (k=0; k<nCells; k++){delete [] nLanded[k];}
	delete [] nLanded;
	delete [] points;
}
/*************************************************************************
                           Allocate Particles
**************************************************************************
allocate proper number of particles, based on concentration differences
------------------------------------------------------------------------*/
void CMOC::AllocateParticles(Ironclad2DArray C){

  int     k,s;
	double Cmax[MAX_SPECIES],Cmin[MAX_SPECIES];
	double CmaxLoc,CminLoc;

	for (s=0; s<nspecies; s++){
		//identify max and min for entire grid
		Cmax[s]=-ALMOST_INF;
		Cmin[s]= ALMOST_INF;
		for (k=0; k<nCells;k++){
			upperswap(Cmax[s],C[k][s]);
			lowerswap(Cmin[s],C[k][s]);
		}
	}
	int    neigh[MAX_CELL_NEIGHBORS];
	double junk [MAX_CELL_NEIGHBORS];
  int numneigh(MAX_CELL_NEIGHBORS);

	for (k=0;k<nCells;k++){nparts[k]=npartlow;}

	for (s=0; s<nspecies; s++){	
		for (k=0;k<nCells;k++){
			//identify local max and min
			CmaxLoc=-ALMOST_INF;
			CminLoc= ALMOST_INF;
			numneigh=MAX_CELL_NEIGHBORS;
			pConcGrid->GetNeighbors(k,neigh,junk,numneigh);
			for (k=0;k<numneigh; k++){
				upperswap(CmaxLoc,C[neigh[k]][s]);
				lowerswap(CminLoc,C[neigh[k]][s]);	
			}
			nparts[k]=nparthigh;//TMP DEBUG
			if (((CmaxLoc-CminLoc)/(Cmax[s]-Cmin[s]))>MOC_ALLOC_TOLERANCE){
				nparts[k]=nparthigh;
			}
		}
	}
}
/************************************************************************
                           WRITE OUTPUT
*************************************************************************/
void CMOC::WriteOutput(const int outputstep){
	//particle allocation
	int k;
	ostrstream file;
	char filename[FILENAME_SIZE];
	filename[0]='\0';  
	file << "MOC_alloc-"<<outputstep<<".txt"<<'\0';
  strcpy(filename,file.str());

	if (outputstep>0){
		ofstream MOCOUT;
		MOCOUT.open(filename);
		MOCOUT << "x y n_part"<<endl;
		for (k=0;k<nCells;k++){
			pt3D pt=pConcGrid->GetNodeLocation(k); 
			MOCOUT << pt.x <<" "<<pt.y<<" "<<nparts[k]<<endl;
		}
		MOCOUT.close();
	}
}