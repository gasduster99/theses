//2DPatchSource.cpp
//Analytic solution for 2D Patch Source at origin

//Cleary and Ungs (from Segol 1994)
//velocity field must be uniform
#include "TransportScheme.h"
#include "AbstractDomain.h" 
#include "AbstractMesh.h"
#include "AnalyticSolutions.h"

/************************************************************************
                           Constructor
************************************************************************/
C2DPatchSource::C2DPatchSource(C2DDomainABC *pDom)
							 :C2DTransportScheme(pDom){
	zc   =0.0;
	a    =10.0;
	alpha=0.0;
	Co   =100.0; //default values
}
//-----------------------------------------------------------------------
C2DPatchSource::~C2DPatchSource(){}
void C2DPatchSource::SetParameters(const cmplex	&z1,
																	 const cmplex &z2,
																	 const double &initconc,
																	 const double &deteriorate){
	zc   =(z1+z2)*0.5;
	a    =abs(z1-z2)*0.5;
	alpha=deteriorate;
	Co   =initconc;
}
/************************************************************************
                           Initialize
************************************************************************/
void C2DPatchSource::Initialize(const double         startt,
																const transport_type ty,
									              Ironclad1DArray      porosity, 
									              Ironclad1DArray      satthick){
	cmplex v;

	if (pDomain==NULL){ExitGracefully("C2DPatchSource::Initialize: NULL Domain",RUNTIME_ERR);}
	if (pDomain->GetNumSpecies()==0){ExitGracefully("C2DPatchSource::Initialize:No species to transport",BAD_DATA);}
	if (ty!=ADV_AND_DISP){ExitGracefully("C2DPatchSource::Initialize:must simulate both dispersion and advection",BAD_DATA);}

	process_type=ty;

	pLayer   =pDomain->GetLayer();
	pConcGrid=pDomain->GetGrid();

	v        =pLayer->GetVelocity2D(zc,startt);
	vx			 =v.real();
	vy			 =v.imag();

	Dx       =pDomain->GetDispersionCoeff(c2Dto3D(zc),startt,OR_XX,1.0)+pDomain->GetDiffusionCoeff(0); 
	Dy       =pDomain->GetDispersionCoeff(c2Dto3D(zc),startt,OR_YY,1.0)+pDomain->GetDiffusionCoeff(0);
	
	R				 =pDomain->GetRetardation    (c2Dto3D(zc),startt,0);

	cout<< "** Patch source initialization **"<<endl;
	cout<< "    vx: "<<vx <<"  vy: "<<vy<<endl;
	cout<< "    Dx: "<<Dx <<"  Dy: "<<Dy<<endl;
	cout<< "    R:  "<<R  <<              endl;	

	if (R<=0.0)								{ExitGracefully("C2DPatchSource::Initialize:negative or zero retardation",BAD_DATA);}
	if (Dx<=0.0)							{ExitGracefully("C2DPatchSource::Initialize:negative or zero longitudinal dispersion",BAD_DATA);}
	if ((Dy<=0.0) && (vy>0.0)){ExitGracefully("C2DPatchSource::Initialize:negative or zero transverse dispersion",BAD_DATA);}
}
/************************************************************************
                           Transport
************************************************************************/
void C2DPatchSource::Transport(const double     t,
															 const double     tstep,
															 Ironclad2DArray  C,
															 Ironclad2DArray  Cend,
															 Writeable2DArray Cnew){
	double integral(0.0),tmp;
	double tau;
	double x,y;
	int i,nNodes;
	int numdiv(DEFAULT_DIVS);
	nNodes=pConcGrid->GetNumNodes();

	double time=t+tstep;

	for (i=0; i<nNodes; i++){
		if (ProgramAborted()){return;}
		WriteEllipsisToScreen(i,nNodes,20);

		x=pConcGrid->GetNodeLocation(i).x-zc.real();
		y=pConcGrid->GetNodeLocation(i).y-zc.imag();

		integral=0.0;
		
		if (x>0.0){
			if (vy!=0.0){ //with vy component--------------------------
				for (tau=(time/R/numdiv)/2.0; tau<(time/R); tau+=(time/R/numdiv)){
					tmp =my_erf(0.5*(a-y)/pow(Dy*tau,0.5)+0.5*vy*pow(tau/Dy,0.5));
					tmp+=my_erf(0.5*(a+y)/pow(Dy*tau,0.5)-0.5*vy*pow(tau/Dy,0.5));
					integral+=(tmp*
						         exp(    -(-alpha*R+(0.25*vx*vx/Dx))*tau    -(0.25*x*x/(Dx*tau))   ) *pow(tau,-1.5)
										 )*(time/R/numdiv);
				}
			}
			else{        //vx component only---------------------------
				for (tau=(time/R/numdiv)/2.0; tau<(time/R); tau+=(time/R/numdiv)){
					tmp =my_erf(0.5*(a-y)/pow(Dy*tau,0.5));
					tmp+=my_erf(0.5*(a+y)/pow(Dy*tau,0.5));
					integral+=(tmp*
						         exp(    -(-alpha*R+(0.25*vx*vx/Dx))*tau    -(0.25*x*x/(Dx*tau))   ) *pow(tau,-1.5)
										 )*(time/R/numdiv);
				}
			}
			integral*=0.25*Co*x/pow(PI*Dx,0.5);
			integral*=exp(0.5*vx*x/Dx-alpha*time);
		}
		else {
			if ((y<a) && (y>-a)){
				integral=Co;
			}
			//integral=0.0;
		}
		//if (integral<0.0){ExitGracefully("C2DPatchSource::Transport : neg conc",RUNTIME_ERR);}
		//if (integral>1.1*Co){ExitGracefully("C2DPatchSource::Transport : bad concentration",RUNTIME_ERR);}
		Cnew[i][0]=integral;
	}
}
void C2DPatchSource::WriteOutput(const int outputstep){
}
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
/************************************************************************
                           Constructor
************************************************************************/
C2DPlaneSource::C2DPlaneSource(C2DDomainABC      *pDom,
		                           const cmplex	     z1a,
								               const cmplex      z2a,
															 const double      deteriorationrate,
							                 const double	    *ConcArray)
							 :C2DAnalyticSol(pDom){
	z1=z1a;
	z2=z2a;
	neg=false;
	alpha=deteriorationrate;
	for (int s=0; s<pDom->GetNumSpecies(); s++){
		Co[s]=ConcArray[s];
	}
}
//-----------------------------------------------------------------------
C2DPlaneSource::~C2DPlaneSource(){}

/************************************************************************
                           Initialize
************************************************************************/
void C2DPlaneSource::Initialize(const double         startt){
	
	if (pDomain==NULL)              {ExitGracefully("C2DPlaneSource::Initialize: NULL Domain",RUNTIME_ERR);}
	if (pDomain->GetNumSpecies()==0){ExitGracefully("C2DPlaneSource::Initialize: No species to transport",BAD_DATA);}

	double vxloc,vyloc;
	double cosa,sina,d;
	cmplex perp_vec;
	cmplex zc=0.5*(z1+z2);
	cmplex v;

	pLayer=pDomain->GetLayer();

	//need to rotate coordinate system
	v        =pLayer->GetVelocity2D(zc,startt);
	vx			 =v.real();
	vy			 =v.imag();

	perp_vec=IM*(z2-z1)/abs(z2-z1);
	d=abs(perp_vec);
	cosa=perp_vec.real()/d;
	sina=perp_vec.imag()/d;
//	cout <<"angle: "<<arg(perp_vec)*180/PI<<endl;
	vxloc=vx*cosa+vy*sina;
	vyloc=vy*cosa-vx*sina;

        CVector tmp = CVector(vx,vy,0.0);
	Dx       =pDomain->CalcDispersionCoeff(tmp,OR_XX,c2Dto3D(IM*(z2-z1)/abs(z2-z1))); 
        CVector t2 = CVector(vx,vy,0.0);
	Dy       =pDomain->CalcDispersionCoeff(t2,OR_YY,c2Dto3D(IM*(z2-z1)/abs(z2-z1)));

  neg=(vxloc<0);
	if (neg){vxloc*=-1.0;vyloc*=-1.0;}

	vx=vxloc;
	vy=vyloc;

	cout<< "** Plane source initialization **"<<endl;
	cout<< "Concentration (species 0): " <<Co[0]<<endl;
	cout<< "    vx: "<<vx <<"  vy: "<<vy<<endl;
	cout<< "    Dx: "<<Dx <<"  Dy: "<<Dy<<endl;	

	if (Dx<=0.0)							{ExitGracefully("C2DPlaneSource::Initialize:negative or zero longitudinal dispersion",BAD_DATA);}
	if ((Dy<=0.0) && (vy>0.0)){ExitGracefully("C2DPlaneSource::Initialize:negative or zero transverse dispersion",BAD_DATA);}
}
/************************************************************************
                           Transport
************************************************************************/
double  C2DPlaneSource::GetConcentration(const cmplex &z,const int s,const double &t) const {
	double integral(0.0),tmp;
	double tau,R;
	double x,y;
	double Df;
	
	int numdiv(DEFAULT_DIVS);

	cmplex zc=0.5*(z1+z2);
  double a= 0.5*abs(z2-z1);

	//need to rotate coordinate system
	double angle=arg(IM*(z2-z1)/abs(z2-z1));
	x=(z.real()-zc.real())*cos(angle)+(z.imag()-zc.imag())*sin(angle);
	y=(z.imag()-zc.imag())*cos(angle)-(z.real()-zc.real())*sin(angle);

	if (neg){x*=-1.0;y*=-1.0;}

	integral=0.0;
	R	=pDomain->GetRetardation (c2Dto3D(zc),starttime,s);
	Df=pDomain->GetDiffusionCoeff(s);
	if (x>=0.0){
		if (vy!=0.0){ //with vy component--------------------------
			for (tau=(t/R/numdiv)/2.0; tau<(t/R); tau+=(t/R/numdiv)){
				tmp =my_erf(0.5*(a-y)/pow((Dy+Df)*tau,0.5)+0.5*vy*pow(tau/(Dy+Df),0.5));
				tmp+=my_erf(0.5*(a+y)/pow((Dy+Df)*tau,0.5)-0.5*vy*pow(tau/(Dy+Df),0.5));
				integral+=(tmp*
						       exp(    -(-alpha*R+(0.25*vx*vx/(Dx+Df)))*tau    -(0.25*x*x/((Dx+Df)*tau))   ) *pow(tau,-1.5)
									)*(t/R/numdiv);
			}
		}
		else{        //vx component only---------------------------
			for (tau=(t/R/numdiv)/2.0; tau<(t/R); tau+=(t/R/numdiv)){
				tmp =my_erf(0.5*(a-y)/pow((Dy+Df)*tau,0.5));
				tmp+=my_erf(0.5*(a+y)/pow((Dy+Df)*tau,0.5));
				integral+=(tmp*
					         exp(    -(-alpha*R+(0.25*vx*vx/(Dx+Df)))*tau    -(0.25*x*x/((Dx+Df)*tau))   ) *pow(tau,-1.5)
									 )*(t/R/numdiv);
			}
		}
		integral*=0.25*Co[s]*x/pow(PI*(Dx+Df),0.5);
		integral*=exp(0.5*vx*x/(Dx+Df)-alpha*t);
	}
	else {
		if ((y<a) && (y>-a)){
			integral=0.0;//Co[s];
		}
		integral=0.0;
	}
	//if (integral<0.0){ExitGracefully("C2DPlaneSource::Transport : neg conc",RUNTIME_ERR);}
	//if (integral>1.1*Co){ExitGracefully("C2DPlaneSource::Transport : bad concentration",RUNTIME_ERR);}
	return integral;
	
}
double C2DPlaneSource::GetDecayLoss(const int s, const double &t) const{return 0.0;}
double C2DPlaneSource::GetMass(const int s, const double &t) const{return 0.0;}//TMP DEBUG

/*****************************************************************************
         PARSE
*****************************************************************************/
C2DPlaneSource *C2DPlaneSource::Parse(ifstream &INPUT, sourcetype ty, C2DDomainABC *pDom, int &l){
	//string "2DPlaneSource" 
	//double x1 double y1 double x2 double y2 double deteriorationrate
	//&
	//{double conc}x(numspecies) ->if C<0, then NO_INFLUENCE
	//&

	C2DPlaneSource *pSrc;

	int          nspecies(pDom->GetNumSpecies());
  bool         eof(false),done(false);
	int          Len,nlines(0),speccount(0);
	char        *s[MAXINPUTITEMS];	
	cmplex			 thisz1,thisz2;
	double       thisdec;
	double       initconc[MAX_SPECIES];
	
  pSrc=NULL;
  if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++; 
	done=false;  
	do {
		if (Len==5) {
      thisz1  =s_to_c(s[0],s[1]);
			thisz2  =s_to_c(s[2],s[3]);
			thisdec =s_to_d(s[4]);                           
      done=true;
		}
		else if(Len==0) {
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
    else {cout <<"line"<< l << "is wrong length (plane source)"<<endl; return NULL; break;}
	} while (!done);

	done=false;                                           
	if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
	speccount=0;
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
			pSrc = new C2DPlaneSource(pDom,
																thisz1,
																thisz2,
																thisdec,
																initconc);
			done=true;
		}
		else if (Len==0){                                   
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
	}while ((!done) && (!eof));

  if (eof) {cout <<"Reached end of file during parse of analytic plane source"<<endl;return NULL;}
	else     {return pSrc;}
}
