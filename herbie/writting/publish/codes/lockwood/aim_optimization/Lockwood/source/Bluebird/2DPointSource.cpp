//-----------------------------------------------------------------------
//2DPointSource.cpp
//Analytic solution for 2D Point Source at specified point
//Hunt 1978 (adapted from Segol 1994)
//velocity field must be uniform
//-----------------------------------------------------------------------

#include "TransportScheme.h"
#include "AbstractDomain.h" 
#include "AbstractMesh.h"
#include "AnalyticSolutions.h"

/************************************************************************
                           Constructor
************************************************************************/
C2DPointSource::C2DPointSource(C2DDomainABC *pDom)
							 :C2DTransportScheme(pDom){
	zc=0.0;
	Mo=100;
	type=INITIAL_CONCENTRATION;
	angle=0.0;
	v=1.0;
	Dl=0.0;
	Dt=0.0;
	//default values
}
//-----------------------------------------------------------------------
C2DPointSource::~C2DPointSource(){}
//-----------------------------------------------------------------------
void C2DPointSource::SetParameters(const cmplex	     zcen,
										               const double	     Mass,
										               const sourcetype  ty){
	zc=zcen;
	Mo=Mass;
	type=ty;
}
/************************************************************************
                           Initialize
-------------------------------------------------------------------------
All values are initialized as if they are homogeneous 
All values are obtained at the source location
This solution will not reflect the flow field unless it is simply uniform
************************************************************************/
void C2DPointSource::Initialize(const double         startt,
																const transport_type ty,
									              Ironclad1DArray      porosity, 
									              Ironclad1DArray      satthick){

	if (pDomain==NULL){ExitGracefully("C2DPointSource::Initialize: NULL Domain",RUNTIME_ERR);}
	if (pDomain->GetNumSpecies()==0){ExitGracefully("C2DPointSource::Initialize:No species to transport",BAD_DATA);}
	if (ty!=ADV_AND_DISP){ExitGracefully("C2DPatchSource::Initialize:must simulate both dispersion and advection",BAD_DATA);}

	process_type=ty;
	//pAquifer =pDomain->GetAquifer();
	pLayer   =pDomain->GetLayer();
	pConcGrid=pDomain->GetGrid();


	cmplex V=pLayer->GetVelocity2D(zc,startt);
	v    =abs(V);
	angle=arg(V);

	cout <<"Velocity: "<<v<<endl;

	Dl			 =pDomain->GetDispersivity(c2Dto3D(zc),LONGITUDINAL)+pDomain->GetDiffusionCoeff(0);
	Dt			 =pDomain->GetDispersivity(c2Dto3D(zc),TRANSVERSE)  +pDomain->GetDiffusionCoeff(0);
	if (process_type==INITIAL_CONCENTRATION){
		Mo=Mo/pLayer->GetSaturatedThickness(zc,startt);
	}
	if (Dl<=0.0){ExitGracefully("C2DPointSource::Initialize: negative or zero longitudinal dispersion",BAD_DATA);}
	if (Dt<=0.0){ExitGracefully("C2DPointSource::Initialize: negative or zero transverse dispersion",BAD_DATA);}
}
/************************************************************************
                           Transport
************************************************************************/
void C2DPointSource::Transport(const double     t,
															 const double     tstep,
															 Ironclad2DArray  C,
															 Ironclad2DArray  Cend,
															 Writeable2DArray Cnew){
	double x,y,r;
	cmplex z;
	int    i,nNodes;
	double time=t+tstep;
	double poro=pLayer->GetPoro(zc); 

	nNodes=pConcGrid->GetNumNodes();

	for (i=0; i<nNodes; i++){
		if (ProgramAborted()){return;}
		WriteEllipsisToScreen(i,nNodes,20);

		z=(c3Dto2D(pConcGrid->GetNodeLocation(i))-zc)*cmplex(cos(angle),-sin(angle)); 
		x=z.real();
		y=z.imag();

		r=sqrt(x*x+y*y*Dl/Dt);
    if (type==INITIAL_CONCENTRATION){
			Cnew[i][0]=exp(-0.25*(x-v*time)*(x-v*time)/(Dl*time)-0.25*y*y/(Dt*time));
			Cnew[i][0]*=Mo/(4.0*PI*poro*time*sqrt(Dl*Dt));
		}
		else if (type==SPECIFIED_CONC){
			Cnew[i][0]=exp(0.5*x*v/Dl);
			Cnew[i][0]*=Mo/(4.0*PI*poro*sqrt(Dl*Dt));
			Cnew[i][0]*=hantush(0.25*r*r/(Dl*time),0.5*v*r/Dl);
		}
		else{
			Cnew[i][0]=0.0;
		}
	}
}
void C2DPointSource::WriteOutput(const int outputstep){
}

/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************/
//-----------------------------------------------------------------------
//2DPointSource
//Analytic solution for 2D Point Source at specified point
//Hunt 1978 (adapted from Segol 1994)
//velocity field must be uniform in vicinity of solution
//-----------------------------------------------------------------------
/************************************************************************
                           Constructor
************************************************************************/
C2DPtSource::C2DPtSource(C2DDomainABC     *pDom,
												 const cmplex	     zcen,
										     const double	    *MassArray,
										     const sourcetype  ty)
						:C2DAnalyticSol(pDom){
	zc=zcen;
	nspecies =pDomain->GetNumSpecies();
	for (int s=0;s<nspecies;s++){
		Mo[s]=MassArray[s];
	}

	type=ty;
	angle=0.0;
	v=1.0;
	Dl=0.0;
	Dt=0.0;
	//default values
}
//-----------------------------------------------------------------------
C2DPtSource::~C2DPtSource(){}
/************************************************************************
                           Initialize
-------------------------------------------------------------------------
All values are initialized as if they are homogeneous 
All values are obtained at the source location
This solution will not reflect the flow field unless it is simply uniform
************************************************************************/
void C2DPtSource::Initialize(const double startt){

	if (pDomain==NULL)              {ExitGracefully("C2DPtSource::Initialize: NULL Domain",RUNTIME_ERR);}
	if (pDomain->GetNumSpecies()==0){ExitGracefully("C2DPtSource::Initialize:No species to transport",BAD_DATA);}
  nspecies=pDomain->GetNumSpecies();

	starttime=startt;

	pLayer=pDomain->GetLayer();

	cmplex V=pLayer->GetVelocity2D(zc,startt);
	v       =abs(V);
	angle   =arg(V);

  //	cout <<"Velocity: "<<v<<endl;


	Dl			 =pDomain->GetDispersivity(c2Dto3D(zc),LONGITUDINAL)*v+pDomain->GetDiffusionCoeff(0);
	Dt			 =pDomain->GetDispersivity(c2Dto3D(zc),TRANSVERSE  )*v+pDomain->GetDiffusionCoeff(0);
	
	cout<< "** Point source initialization **"<<endl;
	cout<< "Mass (species 0): " <<Mo[0]<<endl;
	cout<< "    x:  "<<zc.real()<<"  y : "<<zc.imag()<<endl;
	cout<< "    vx: "<<V.real() <<"  vy: "<<V.imag()<<endl;
	cout<< "    Dl: "<<Dl <<"  Dt: "<<Dt<<endl;	

	if (type==INITIAL_CONCENTRATION){
		for (int s=0;s<nspecies;s++){
			Mo[s]=Mo[s]/pLayer->GetSaturatedThickness(zc,startt);
		}
	}

	if (Dl<=0.0){ExitGracefully("C2DPtSource::Initialize: negative or zero longitudinal dispersion",BAD_DATA);}
	if (Dt<=0.0){ExitGracefully("C2DPtSource::Initialize: negative or zero transverse dispersion"  ,BAD_DATA);}
}
/************************************************************************
                           Transport
************************************************************************/
double C2DPtSource::GetConcentration(const cmplex &z, const int s, const double &t) const{
	
	double x,y,r;
	double C;
	double poro=pLayer->GetPoro(zc); 

	//transform time and space to local coordinates
  double time=t-starttime;
	cmplex zq=(z-zc)*cmplex(cos(angle),-sin(angle)); 
	x=zq.real();
	y=zq.imag();

	r=sqrt(x*x+y*y*Dl/Dt);

  if      (type==INITIAL_CONCENTRATION){
		C=exp(-0.25*(x-v*time)*(x-v*time)/(Dl*time)-0.25*y*y/(Dt*time));
		C*=Mo[s]/(4.0*PI*poro*time*sqrt(Dl*Dt));
	}
	else if (type==SPECIFIED_CONC){
		C=exp(0.5*x*v/Dl);
		C*=Mo[s]/(4.0*PI*poro*sqrt(Dl*Dt));
		C*=hantush(0.25*r*r/(Dl*time),0.5*v*r/Dl);
	}
	else{
		C=1.0;
	}

	return C/pLayer->GetSaturatedThickness(zc,0.0);
}
double  C2DPtSource::GetMass(const int s, const double &t) const{return Mo[s];}
double  C2DPtSource::GetDecayLoss(const int s, const double &t) const{return 0.0;}

/*****************************************************************************
         PARSE
*****************************************************************************/
C2DPtSource *C2DPtSource::Parse(ifstream &INPUT, sourcetype ty, C2DDomainABC *pDom, int &l){
	//string "2DInitialPtSource" or "2DConstantPtSource" or other such thing
	//double x double y 
	//&
	//{double conc}x(numspecies) ->if C<0, then NO_INFLUENCE
	//&

	C2DPtSource *pSrc;

	int          nspecies(pDom->GetNumSpecies());
  bool         eof(false),done(false);
	int          Len,nlines(0),speccount(0);
	char        *s[MAXINPUTITEMS];	
	cmplex			 thisz;
	double       initconc[MAX_SPECIES];
	
  pSrc=NULL;
	if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++; 
	done=false;  
	do {
		if (Len==2) {
      thisz  =s_to_c(s[0],s[1]);
      done=true;
		}
		else if(Len==0) {                                   
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
    else {cout <<"line"<< l << "is wrong length (point source)"<<endl; return NULL;break;}
	} while (!done);

	done=false;                                           
	speccount=0;
	if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
	do {
		if (ProgramAborted()){return NULL;}
		if ((Len==1) && (speccount<nspecies)){
			initconc[speccount]=s_to_d(s[0]); 
			//cout <<"species "<<speccount<<" added"<<endl;
			speccount++;                                      
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
		else if ((Len==1) && (speccount<nspecies) && (strcmp(s[0],"&"))){
			initconc[speccount]=s_to_d(s[0]); 
			speccount++;                                     
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
		else if ((Len>=1)  && (speccount==nspecies) && (strcmp(s[0],"&"))){ 
			cout <<"Too many species associated with source"<<endl;  return NULL;	
		}
		else if ((Len==1)  && (speccount<nspecies) && (!strcmp(s[0],"&"))){ 
			cout <<"Not enough species associated with source"<<endl;  return NULL;	
		}		
		else if ((Len==1) && (speccount==nspecies) && (!strcmp(s[0],"&")) ) {//hit a &
			pSrc = new C2DPtSource(pDom,
														 thisz,
														 initconc,
														 ty);
			done=true;
		}
		else if (Len==0){                                   
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
	}while ((!done) && (!eof));

  if (eof) {cout <<"Reached end of file during parse of analytic plane source"<<endl;return NULL;}
	else     {return pSrc;}
}