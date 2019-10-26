//MMOC.cpp
#include "CharacteristicMethods.h"
#include "MassParticle.h"
#include "FluxGrid.h"
/************************************************************************
                           Constructors
************************************************************************/
CMMOC::CMMOC(C2DDomainABC *pDom)
      :C2DTransportScheme(pDom){
	 pMP=NULL;
   nNodes=0;
	 BacktrackLocations=NULL;
	 tracked=false;
	 pEulerian=NULL;
	 pPollockGrid=NULL;
}
//------------------------------------------------------------------------
CMMOC::~CMMOC(){
	if (globaldebug){cout << "  DESTROYING MMOC STYLE"<<endl;}
	if (pMP!=NULL){
		for (int s=0; s<nspecies; s++){delete pMP[s];}
		delete [] pMP;
	}
	if (BacktrackLocations!=NULL){
		for (int k=0;k<nNodes;k++){
			delete [] BacktrackLocations[k];}}
	delete [] BacktrackLocations;
	delete pPollockGrid;
}
void CMMOC::SetParameters(const int numpercell, CRectGrid *pGrid){
	if (pGrid!=NULL){pPollockGrid=new CFluxGrid(pGrid);}
}
/************************************************************************
                           Initialize
************************************************************************/
void CMMOC::Initialize(const double         startt,
											 const transport_type ty,
									     Ironclad1DArray      porosity, 
									     Ironclad1DArray      satthick){

	if (ty==DISPERSION_ONLY){
		ExitGracefully("CMMOC::Initialize: cannot simulate dispersion without advection",BAD_DATA);}

	//if (pEulerian==NULL){ExitGracefully("CMMOC::Initialize: NULL Eulerian subscheme",BAD_DATA);}

	int     k,s;

	process_type=ty;

	nspecies =pDomain  ->GetNumSpecies(); 
	pConcGrid=pDomain  ->GetGrid();
	pAquifer =pDomain  ->GetAquifer();
	nNodes   =pConcGrid->GetNumNodes();
	tracked  =false;

	//create single particles (one per species)->concentrations are blank, location at (x,y)=(0,0)
	pMP=new CMassParticle *[nspecies];
	for (s=0; s<nspecies; s++){
		pMP[s]	=new CMassParticle(pAquifer,pDomain,0.0,TRACK_BACKWARD,0.0,s);
	}

	BacktrackLocations=new pt3D *[nNodes];
	for (k=0;k<nNodes;k++){
		BacktrackLocations[k]=new pt3D [nspecies];
		for (s=0; s<nspecies; s++){
			BacktrackLocations[k][s]=0.0;
		}
	}
	if (pPollockGrid!=NULL){pPollockGrid->Initialize(pDomain->GetLayer());}
}
/************************************************************************
                           TRANSPORT
************************************************************************/
void CMMOC::Transport(const double     t,
											const double     tstep,
											Ironclad2DArray  C,
                      Ironclad2DArray  Cend,
											//Ironclad2DArray  dCdt_sim,
											Writeable2DArray Cnew){
	int			k,s;

	static disp_info disp[MAX_SPECIES];
	for (s=0;s<nspecies;s++){
		pDomain->GetDispersivities(0.0,s,disp[s]); 
	}

	if (pEulerian==NULL){ExitGracefully("CMMOC::Transport: NULL Eulerian scheme",BAD_DATA);}

	//Pre-track particles for steady state (should also be repeated for non-steady state simulations)
	if (!tracked){//
		cout <<endl<<"Initializing MMOC backward characteristic paths";

		for (k=0;k<nNodes;k++){
			
			if (ProgramAborted()){break;}
			WriteEllipsisToScreen(k,nNodes,20);

			for (s=0;s<nspecies;s++){
				pMP[s]->SetLocation(pConcGrid->GetNodeLocation(k));
				if (pPollockGrid==NULL)
				{
					pMP[s]->Track      (tstep,advectiontype,adv_timestep,adv_spacestep,pDomain->UseEffectiveParams(),disp[s],false); 
				}
				else
				{
					pMP[s]->TrackPollock(tstep,1,pPollockGrid,false,false);
				}
				BacktrackLocations[k][s]=pMP[s]->GetLocation();
			}
		}	
		cout <<endl;
		tracked=true;
	}

	for (k=0;k<nNodes;k++){
		WriteEllipsisToScreen(k,nNodes,20);
		for (s=0;s<nspecies;s++){
			pMP[s]->SetLocation(BacktrackLocations[k][s]);
			if (process_type==ADV_AND_DISP){
				pMP[s]->Disperse(tstep,t);
			}
			//assign particle Concentrations at start of time step to cells
			Cnew[k][s]=pDomain->GetConcentration(pMP[s]->GetLocation(),t,s); 
		}
	}//for (k=0..

	/*if (pEulerian!=NULL){
		pEulerian->Transport(t,tstep,C,Cnew(updated),Cnew);
	}*/
  
}
//************************************************************************
//                           WRITE OUTPUT
//************************************************************************
void CMMOC::WriteOutput(const int outputstep){
	//no output 
}