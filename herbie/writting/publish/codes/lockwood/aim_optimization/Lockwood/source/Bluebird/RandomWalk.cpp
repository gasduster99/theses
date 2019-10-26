//RandomWalk.cpp
#include "TransportScheme.h"
#include "MassParticle.h"
/************************************************************************
                           Constructors
************************************************************************/
CRandomWalk::CRandomWalk(C2DDomainABC *pDom)
            :C2DTransportScheme(pDom){

   nCells=0;
	 nMassParticles=1000; //Default value
	 UseEPVA    =true;
	 pMP        =NULL;
	 CellIndices=NULL;
	 CellVolume =NULL;
	 active     =NULL;
	 pEulerian  =NULL;
}
//------------------------------------------------------------------------
CRandomWalk::~CRandomWalk(){
	if (globaldebug){cout << "  DESTROYING RANDOM WALK STYLE"<<endl;}
	delete [] CellIndices;
	delete [] CellVolume;
	delete [] active;

	if (pMP!=NULL){for (int m=0;m<nMassParticles;m++){delete pMP[m]; }} delete [] pMP;
}
/************************************************************************
                           SetParams
************************************************************************/
void CRandomWalk::SetParameters(const int numparts, bool EPVA){
	nMassParticles=numparts;
	UseEPVA=EPVA;
}
/************************************************************************
                           Initialize
************************************************************************/
void CRandomWalk::Initialize(const double         startt,
														 const transport_type ty,
									           Ironclad1DArray      porosity, 
									           Ironclad1DArray      satthick){

  int    k,m;
  int    nparts, partscreated(0);
	pt3D   pts[20000],pt;
	int    mp;
	double totalmass(0.0),mass;

	if (ty==DISPERSION_ONLY){
		ExitGracefully("CRandomWalk::Initialize: cannot simulate dispersion without advection",BAD_DATA);}
	if (nMassParticles<=0){
		ExitGracefully("CRandomWalk::Initialize: number of particles not set",BAD_DATA);}

	process_type=ty;

	nspecies =pDomain  ->GetNumSpecies();
	pAquifer =pDomain  ->GetAquifer();
	pConcGrid=pDomain  ->GetGrid();
	nCells   =pConcGrid->GetNumCells();

	if (pConcGrid->GetType()!=CELL_BASED){
		ExitGracefully("Random Walk cannot be applied to non-cell-based grid or mesh",BAD_DATA);}
	if (nspecies>1){
		ExitGracefully("Random Walk can only simulate single species transport (for now)",BAD_DATA);
	}
  cin>>m;
	//if (pEulerian==NULL){ExitGracefully("CRandomWalk::Initialize: NULL Eulerian scheme",BAD_DATA);}
//	pEulerian->Initialize(startt,DISPERSION_ONLY,porosity,satthick);

	pMP	       =new CMassParticle *[nMassParticles];
	CellVolume =new double         [nCells];
	CellIndices=new int            [nMassParticles];
	active     =new bool           [nMassParticles];
	if (active==NULL){ExitGracefully("CRandomWalk::Initialize(1)",OUT_OF_MEMORY);} 
	for (k=0; k<nCells; k++){
		CellVolume[k]=satthick[k]*porosity[k]*pConcGrid->GetCellArea(k);
	}
	for (m=0; m<nMassParticles; m++){pMP[m]=NULL;}
  //create particles---------------------------------------------------
  
  //Calculate total mass in domain
	for (k=0; k<nCells;k++){
		pt=pConcGrid->GetNodeLocation(k);
		totalmass+=pDomain->GetConcentration(pt,0.0,0)*CellVolume[k];
	}
	cout <<"total mass in domain: "<<totalmass<<endl;

	//Allocate mass to particles
	for (k=0; k<nCells;k++){

		mass=pDomain->GetConcentration(pConcGrid->GetNodeLocation(k),0.0,0)*CellVolume[k];
		
		if (mass>0){
		  nparts=(int)(mass/totalmass*0.98*nMassParticles);

			cout<<nparts; cin>>m;
			if ((nparts>=20000) || (nparts<=0)){
				ExitGracefully("Random Walk: exceeded max particles per cell",BAD_DATA);}

			pConcGrid->DistributePoints(k,RANDOM_PATTERN,nparts,pts); 

			for (int m=0; m<nparts; m++){
				mp=m+partscreated;
				if (mp>=nMassParticles){
					ExitGracefully("Random Walk: created too many particles",BAD_DATA);}
				
				pMP        [mp]=new CMassParticle(pAquifer,pDomain,pts[m],TRACK_FORWARD,0.0,0);//single species, for now 
				if (pMP[mp]==NULL){ExitGracefully("CRandomWalk::Initialize(2)",OUT_OF_MEMORY);} 
				pMP        [mp]->SetMass(mass/nparts);
				CellIndices[mp]=k;
				active     [mp]=true;

				//cout << "particle created (in cell: "<<k<<")- mass: "<<mass/nparts<<endl;
			}
			partscreated+=nparts;
		}
	}

	cout << partscreated << " out of " << nMassParticles << "created."<<endl;

	//inactive particles
	for (m=0; m<nMassParticles-partscreated; m++){
		mp=m+partscreated;	
		ExitGracefullyIf(mp>=nMassParticles,"Random Walk: created too many particles",BAD_DATA); //Shouldn't be able to happen
		pMP        [mp]=new CMassParticle(pAquifer,pDomain,0.0,TRACK_FORWARD,0.0,0);//single species, for now 
	  if (pMP[m+partscreated]==NULL){ExitGracefully("CRandomWalk::Initialize(3)",OUT_OF_MEMORY);} 
		pMP        [mp]->SetMass(0.0);
		CellIndices[mp]=0;
		active     [mp]=false;
	}	

	if (partscreated<=0){
		ExitGracefully("CRandomWalk::Initialize: created no particles (zero initial mass?}",BAD_DATA);}

  ofstream DBGPARTLOC;
	DBGPARTLOC.open("ParticleLocations.csv");
	DBGPARTLOC.close();

}
/************************************************************************
                           TRANSPORT
************************************************************************/
void CRandomWalk::Transport(const double     t,
														const double     tstep,
														Ironclad2DArray  C,
														Ironclad2DArray  Cend,
														//Ironclad2DArray  dCdt_sim,
														Writeable2DArray Cnew){

  int           k,m,s;



	//Re-assign Mass to particles--------------------------------
	/*for (m=0; m<nMassParticles; m++){
		if (active[m]){
			for (s=0; s<nspecies;s++){
				pMP[m]->SetMass(C[CellIndices[m]][0]*CellVolume[CellIndices[m]]);//one species for now
			}
		}
	}*/

	disp_info disp;
	pDomain->GetDispersivities(0.0,0,disp); //TMP DEBUG

	//cout <<"tracking"<<endl;
	//track all of the particles - advection & dispersion---------------
	for (m=0; m<nMassParticles; m++){
		if (active[m]){
			WriteEllipsisToScreen(m,nMassParticles,20);
			pMP[m]->Track   (tstep,advectiontype,adv_timestep,adv_spacestep,UseEPVA,disp,false); //TMP DEBUG - first boolean should be true for eff.vel
			pMP[m]->Disperse(tstep,t);
		}
	}

	//place on grid, reassign grid C -------

	//associate particles with new grid/mesh cells---------------------------

	//cout <<"re-associating"<<endl;
	for (m=0; m<nMassParticles; m++){
		if (active[m]){
			if (pConcGrid->GetCellIndex(pMP[m]->GetLocation(),k)){ CellIndices[m]=k;                }
			else                                                 { CellIndices[m]=0;active[m]=false;}//not on grid
		}
	}

	//cout <<"re-evaluating"<<endl;
	//re-evaluate C------------------------------------------------------
	for (k=0; k<nCells; k++){
		for (s=0;s<nspecies; s++){
			Cnew[k][s]=0.0;
		}
	}
	for (m=0; m<nMassParticles; m++){
		if ((active[m]) && (CellIndices[m]>=0) && (CellIndices[m]<nCells)){
			Cnew[CellIndices[m]][0]+=pMP[m]->GetMass(); //TMP DEBUG - one species only
		}
	}
  for (k=0; k<nCells; k++){
		for (s=0;s<nspecies; s++){
			Cnew[k][s]/=CellVolume[k];
		}
	}

	//cout <<"done"<<endl;
	//delete particles if too many are in a cell------------------------

	//source zones------------------------------------------------------

}
/************************************************************************
                           MASS BALANCE OPERATIONS
************************************************************************/
double CRandomWalk::CalculateSystemMass(const int        s, 
		                                    Ironclad2DArray  C,
																				Ironclad2DArray  Cs,
																		    double &SorbedMass) const{
	int k,m;
	double TotalMass(0.0);

	if (s==0){
		for (m=0; m<nMassParticles; m++){
			if (active[m]){TotalMass+=pMP[m]->GetMass();}
		}
		for (k=0; k<nCells; k++){
			/*SorbedMass+=
				CellVolume[k]*pDomain->TranslateToSorbed(pConcGrid->GetCellCenter(k), , const int s) const{*/
		}

		return TotalMass;
	}
	else{
		return 0;//	return pEulerian->CalculateSystemMass(s,s2,C,Cs,SorbedMass);
	}
}
//------------------------------------------------------------------------
double CRandomWalk::CalculateBorderFlux(const int        s,
														            Ironclad2DArray  C,
																				Ironclad2DArray  Cend,
														            const double    &t,
																				const double    &tstep) const{
	int m;
	double TotalMass(0.0);

	if (s==0){
		for (m=0; m<nMassParticles; m++){
			if (!active[m]){
				TotalMass+=pMP[m]->GetMass();
				pMP[m]->SetMass(0.0);
			}
		}
		return TotalMass;
	}
	else{
		return 0.0;//pEulerian->CalculateBorderFlux(s,C,Cend,t,tstep);
	}
}
//------------------------------------------------------------------------
double CRandomWalk::CalculateSourceGain(const int        s,
														            Ironclad2DArray  C,
																				Ironclad2DArray  Cend,
														            const double    &t,
																				const double    &tstep) const{
	return 0.0;//pEulerian->CalculateSourceGain(s,C,Cend,t,tstep);
}
//------------------------------------------------------------------------
double CRandomWalk::CalculateSinkLoss(const int        s,
														          Ironclad2DArray  C,
																			Ironclad2DArray  Cend,
														          const double    &t,
																			const double    &tstep) const{
	int m;
	double TotalMass(0.0);

	if (s==0){
		for (m=0; m<nMassParticles; m++){
			if ((active[m]) && (pMP[m]->IsCaptured())){
				TotalMass+=pMP[m]->GetMass();
				pMP[m]->SetMass(0.0);
				active[m]=false;
			}
		}
		return TotalMass;
	}
	else{
		return 0.0;//pEulerian->CalculateSinkLoss(s,C,Cend,t,tstep);
	}
}
//------------------------------------------------------------------------
double CRandomWalk::CalculateDecayLoss(const int        s,
														          Ironclad2DArray  C,
																			Ironclad2DArray  Cend,
														          const double    &t,
																			const double    &tstep) const{
	return 0.0;//pEulerian->CalculateDecayLoss(s,C,Cend,t,tstep);
}
/************************************************************************
                           WRITE OUTPUT
************************************************************************/
void CRandomWalk::WriteOutput(const int outputstep){
	int on;
  ofstream DBGPARTLOC;
	DBGPARTLOC.open("ParticleLocations.csv",ios::app);
	for (int m=0; m<nMassParticles; m++){
		if (active[m]){on=1;}else {on=0;}
		pt3D tmp=pMP[m]->GetLocation();
		DBGPARTLOC<<tmp.x<<","<<tmp.y<<","<<tmp.z<<","<<on<<","<<outputstep<<endl;
	}
	DBGPARTLOC.close();
}