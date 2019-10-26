//2DAnalyticTransport.cpp
//Transport scheme for superimposed analytic solutions

#include "MasterInclude.h"
#include "TransportScheme.h"
#include "AbstractDomain.h" 
#include "AbstractMesh.h"
#include "AnalyticSolutions.h"

/************************************************************************
                           Constructor
************************************************************************/
C2DAnalytic::C2DAnalytic(C2DDomainABC *pDom)
						:C2DTransportScheme(pDom){
	pSolutions=NULL;
	nSolutions=0; 
}
//-----------------------------------------------------------------------
C2DAnalytic::~C2DAnalytic(){
  for (int i=0;i<nSolutions;i++){delete pSolutions[i];}
	delete [] pSolutions;
}
/************************************************************************
                           Initialize
************************************************************************/
void C2DAnalytic::Initialize(const double         startt,
														 const transport_type ty,
									           Ironclad1DArray      porosity, 
									           Ironclad1DArray      satthick){
	
	if (pDomain==NULL)              {ExitGracefully("C2DAnalytic::Initialize: NULL Transport Domain",RUNTIME_ERR);}
	if (pDomain->GetNumSpecies()==0){ExitGracefully("C2DAnalytic::Initialize: No species to transport",BAD_DATA);}
	if (ty!=ADV_AND_DISP)           {ExitGracefully("C2DAnalytic::Initialize: Analytic transport method must simulate both dispersion and advection",BAD_DATA);}

	process_type=ty;

	pLayer   =pDomain->GetLayer();
	pConcGrid=pDomain->GetGrid();
	nspecies =pDomain->GetNumSpecies();

	cout <<"Analytic Transport: initializing (" << nSolutions << " local analytic solutions)..."<<endl;
	for (int i=0;i<nSolutions;i++){
		pSolutions[i]->Initialize(startt);
	}
}
/************************************************************************
                           Transport
************************************************************************/
void C2DAnalytic::Transport(const double     t,
														const double     tstep,
														Ironclad2DArray  C,
														Ironclad2DArray  Cend,
														Writeable2DArray Cnew){

	int i,nNodes;
	nNodes=pConcGrid->GetNumNodes();

	double time=t+tstep;

	//cout <<endl<<"Analytic Transport (" << nSolutions << " solutions)";
	for (i=0; i<nNodes; i++){
    if (ProgramAborted()){return;}
		WriteEllipsisToScreen(i,nNodes,20);

		for (int s=0; s<nspecies; s++){
			Cnew[i][s]=0.0;
			for (int j=0;j<nSolutions;j++){
				Cnew[i][s]+=pSolutions[j]->GetConcentration(c3Dto2D(pConcGrid->GetNodeLocation(i)),s,t+tstep);
			}
		}
	}
}
/************************************************************************
                           AddSolution
************************************************************************/
void C2DAnalytic::AddSolution(C2DAnalyticSol *pSol){
	if (!DynArrayAppend((void**&)pSolutions,(void*)(pSol),nSolutions)){
			ExitGracefully("C2DAnalytic::AddSolution: adding NULL analytic transport solution",BAD_DATA);};
}
void C2DAnalytic::WriteOutput(const int outputstep){
}
