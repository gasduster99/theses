#include "MassParticle.h"

//*********************************************************************
CMassParticle::CMassParticle(){}
//---------------------------------------------------------------------
CMassParticle::CMassParticle(const CAquiferABC   *pAquifer, 
														 const C2DDomainABC  *pDom,  //Not neccesary-could be 3D
								             pt3D							    startpt,
														 trackdir             direct,
								             double						    initconc,
														 int                  species_ind):
							     CParticle(0,pAquifer,direct){

	lastpt[0]=lastpt[1]=lastpt[2]=currpt=startpt;
	lastt [0]=lastt [1]=lastt [2]=currt=0.0;
 
	captured=false;
	s=species_ind;
	pDomain =pDom;
	concentration=initconc;

}
//---------------------------------------------------------------------
CMassParticle::~CMassParticle(){
	if (globaldebug){cout<<"DESTROYING MASS PARTICLE"<<endl;}
}
/*********************************************************************
   ACCESSOR FUNCTIONS
*********************************************************************/
bool   				CMassParticle::IsCaptured()			   const {return captured;}
//---------------------------------------------------------------------
double				CMassParticle::GetCurrentT()			 const {return currt;}
//---------------------------------------------------------------------
double				CMassParticle::GetLastT   ()			 const {return lastt[0];}    
//---------------------------------------------------------------------
pt3D  				CMassParticle::GetLocation()       const {return currpt;}
//---------------------------------------------------------------------
const double  CMassParticle::GetConcentration()  const {return concentration;}
//---------------------------------------------------------------------
const double  CMassParticle::GetMass()           const {return concentration;}
/*********************************************************************
   MANIPULATOR FUNCTIONS
*********************************************************************/
void CMassParticle::SetConcentration(const double conc){
	 concentration=conc;
}
//---------------------------------------------------------------------
void CMassParticle::SetMass(const double mass){
	 concentration=mass;
}
//---------------------------------------------------------------------
void CMassParticle::SetLocation(const pt3D &pt){
  lastpt[0]=lastpt[1]=lastpt[2]=currpt=pt;
	lastt [0]=lastt [1]=lastt [2]=currt;
	captured=false;
}
/*********************************************************************
   INHERITED MEMBER FUNCTIONS
*********************************************************************/
bool CMassParticle::HasBeenCaptured(){
	if (captured){return true;}
  if ((lastt[0]!=lastt[1]) &&
			(lastt[1]!=currt) &&
		 ((abs(currpt-lastpt[1])/
			(abs(currpt-lastpt[0])+abs(lastpt[0]-lastpt[1]) ) )<=BAD_TURN) &&
		 ((abs(lastpt[0]-lastpt[2])/
		  (abs(lastpt[0]-lastpt[1])+abs(lastpt[1]-lastpt[2])))<=BAD_TURN)){

		return true;	
	}
	return false;
}
//---------------------------------------------------------------------
void CMassParticle::Capture(const double &endtime){
	currt=endtime;
	lastt[0]=lastt[1]=lastt[2]=endtime;
	currpt=lastpt[2];
	captured=true;
}
//---------------------------------------------------------------------
void CMassParticle::Advect(const vector &movement, const double &tstep){
	static double R;
	lastpt[2]=lastpt[1]; lastpt[1]=lastpt[0]; lastpt[0]=currpt;
	lastt [2]=lastt [1]; lastt [1]=lastt [0]; lastt [0]=currt;
  if (!captured){
	  R=pDomain->GetRetardation(currpt,currt,s); 
		currpt+=(1.0/R)*movement;
	}
  currt+=tstep;
}
//---------------------------------------------------------------------
void CMassParticle::Disperse(const double &tstep, const double &t){

	static double R;
	R=pDomain->GetRetardation (currpt,currt,s);

	currpt.x+=2.0/R*pDomain->GetDispersionCoeff(currpt,t,OR_XX,1.0)*tstep*gaussrand();
	currpt.y+=2.0/R*pDomain->GetDispersionCoeff(currpt,t,OR_YY,1.0)*tstep*gaussrand();
	currpt.z+=2.0/R*pDomain->GetDispersionCoeff(currpt,t,OR_ZZ,1.0)*tstep*gaussrand();
}
//---------------------------------------------------------------------
void CMassParticle::BackTrack(){
	currpt=lastpt[0]; lastpt[0]=lastpt[1]; lastpt[1]=lastpt[2];
	currt =lastt [0]; lastt [0]=lastt [1]; lastt [1]=lastt [2];
}
//*********************************************************************
//   MEMBER FUNCTIONS
//*********************************************************************
void CMassParticle::WriteOutput() const{}

