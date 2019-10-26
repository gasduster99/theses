#include "PointElems.h"

/************************************************************************
                  CPointElem Class
************************************************************************/
/***********************************************************************
					CONSTRUCTOR
***********************************************************************/
CPointElem::CPointElem(){}
//-----------------------------------------------------------------------
CPointElem::CPointElem(char      *Name, 
											 CSingleLayerABC *pLay, 
											 pointtype  PType, 
											 double     ang, 
											 double     str, 
											 cmplex     z):
						CAnalyticElem(Name,pLay){
  Ptype=PType;    
	angle=0.0;     
	Q=str;           
	zp=z;
}
//-----------------------------------------------------------------------
CPointElem::~CPointElem(){
  if (globaldebug){cout <<"   DESTROYING POINTELEM "<<name<<endl;}
}

/***********************************************************************
					ACCESSORS
***********************************************************************/
cmplex CPointElem::GetZ()    const {return zp;}
//-----------------------------------------------------------------------
double CPointElem::GetArea() const {return 0.0;}
//-----------------------------------------------------------------------
bool   CPointElem::HasFlux() const {return (Ptype==WELL);}
/***********************************************************************
					SetStrength
***********************************************************************/
void CPointElem::SetStrength(double strength){Q=strength;}
//-----------------------------------------------------------------------
void CPointElem::UpdateBlock(const double &t) const {
	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,NOT_A_STRING,t);
	}
}
//-----------------------------------------------------------------------
double CPointElem::GetMaxError(const double &t) const	{return 0.0;}
/***********************************************************************
					GEOMETRY MEMBER FUNCTIONS
***********************************************************************/
void CPointElem::WriteGeometry(ofstream &BASEMAP) const              {
	BASEMAP << "\" PointElem \",  1" <<endl;
	BASEMAP << zp.real()<< " , "    <<zp.imag()<<endl;
}
//-----------------------------------------------------------------------
cmplex CPointElem::Centroid() const {return zp;}
//-----------------------------------------------------------------------
bool CPointElem::IsInside(const cmplex &z) const{
  if (zp==z){return true;} else {return false;}
}
//-----------------------------------------------------------------------
bool CPointElem::IsInSquare(const cmplex &zc,const double w) const{
	if   ((zp.real() > (zc.real() + (w/2.0))) ||
        (zp.real() < (zc.real() - (w/2.0))) ||
        (zp.imag() > (zc.imag() + (w/2.0))) ||
        (zp.imag() < (zc.imag() - (w/2.0)))) {return false;} else {return true;}
}
//-----------------------------------------------------------------------
bool CPointElem::IsInCircle(const cmplex &zc,const double r) const{
	if (abs(zp-zc)>r) {return false;} else {return true;}
}
//-----------------------------------------------------------------------
bool CPointElem::PartInCircle    (const cmplex &zc,const double r) const{
	if (abs(zp-zc)>r) {return false;} else {return true;}
}
//-----------------------------------------------------------------------
bool CPointElem::SharesNode(const cmplex &zn) const{
	if (zn==zp){return true;}
	return false;
}
/***********************************************************************
					GetDischargePotential
----------------------------------------------------------------------*/
cmplex CPointElem::GetDischargePotential(const cmplex &z,const double &t) const{
	static cmplex omega;
	if (Ptype==WELL){	
		if (z!=zp){omega=(IOVER2PI*Q)*log((z-zp)              /abs(zp-(pLayer->GetBlackHole())));}
		else      {omega=(IOVER2PI*Q)*log((z+(MOVE_DIST*IM)-zp)/abs(zp-(pLayer->GetBlackHole())));}

		if (CAnalyticElem::BranchcutLocus){
			double argz=arg((z              -zp));
			double arg1=arg(-(zBranchcutLocus-zp));//if negative, away from BCL, else towards
			if      ((arg1>=0) && (argz>0.0) && (argz>arg1)) {omega-=IM*Q;}
			else if ((arg1< 0) && (argz<0.0) && (argz<arg1)) {omega+=IM*Q;}
		}
	}
	else if (Ptype==VORTEX){
		if (z!=zp){omega=(-IM*IOVER2PI*Q)*log((z-zp)              /abs(zp-(pLayer->GetBlackHole())));}
		else      {omega=(-IM*IOVER2PI*Q)*log((z+(MOVE_DIST*IM)-zp)/abs(zp-(pLayer->GetBlackHole())));}
	} 
	else if (Ptype==PTDIPOLE){
		if (z!=zp){omega=(-Q*cmplex(cos(angle),sin(angle))/(2.0*PI*(z-zp)));              }
		else      {omega=(-Q*cmplex(cos(angle),sin(angle))/(2.0*PI*(z+(MOVE_DIST*IM)-zp)));}
	} 
	else {omega=0.0;}

	return omega;
}

/***********************************************************************
					GetW
----------------------------------------------------------------------*/
cmplex CPointElem::GetW(const cmplex &z,const double &t) const{
	if (Ptype==WELL){
		if (z!=zp){return -IOVER2PI*Q/(z-zp);              }
		else{      return -IOVER2PI*Q/(z+(MOVE_DIST*IM)-zp);}
	}
	else if (Ptype==VORTEX){
		if (z!=zp){return IM*IOVER2PI*Q/(z-zp);              }
		else{      return IM*IOVER2PI*Q/(z+(MOVE_DIST*IM-zp));}
	} 
	else if (Ptype==PTDIPOLE){
		if (z!=zp){return (-Q*cmplex(cos(angle),sin(angle))/(4.0*PI*(z-zp)*(z-zp)));              }
		else{      return (-Q*cmplex(cos(angle),sin(angle))/(4.0*PI*(z+(MOVE_DIST*IM))*(z+(MOVE_DIST*IM))));}
	} 
	else {return 0.0;}
}
/***********************************************************************
							 Get Flux Through Face
----------------------------------------------------------------------*/
cmplex CPointElem::GetFluxThruFace(const cmplex &z1, const cmplex &z2, const double &t) const{

	//return GetDischargePotential(z1,t).imag()-GetDischargePotential(z2,t).imag();
  return IM*conj(GetDischargePotential(z1,t)-GetDischargePotential(z2,t)); //complex flux thru face

}
/***********************************************************************
					GetGx
----------------------------------------------------------------------*/
cmplex CPointElem::GetGx(const cmplex &z,const double &t) const{
	if      (Ptype==WELL){
		if (z!=zp){return (IOVER2PI*Q/(z-zp)/(z-zp));        }
		else{      return (IOVER2PI*Q/(z+(MOVE_DIST*IM)-zp)/(z+(MOVE_DIST*IM)-zp));}
	}
	else if (Ptype==VORTEX){
		if (z!=zp){return -IM*IOVER2PI*Q/(z-zp)/(z-zp);              }
		else{      return -IM*IOVER2PI*Q/(z+(MOVE_DIST*IM)-zp)/(z+(MOVE_DIST*IM)-zp);}
	} 
	else if (Ptype==PTDIPOLE){
		//if (z!=zp){return (Q*cmplex(cos(angle),sin(angle))/(4.0*pi*(z-zp)*(z-zp)*(z-zp)));              }
		//else{      return (Q*cmplex(cos(angle),sin(angle))/(4.0*pi*(z+(MOVE_DIST*IM))*(z+(MOVE_DIST*IM))*(z+(MOVE_DIST*Im))));}
		return 0.0; //TMP DEBUG
	} 
	else {return 0.0;}
}
/***********************************************************************
					GetLeakage
----------------------------------------------------------------------*/
double CPointElem::GetLeakage(const cmplex &z, const double &t,const leak_type ltype) const{
	return 0.0;
}
/***********************************************************************
					GetNetDischarge
----------------------------------------------------------------------*/
double CPointElem::GetNetDischarge(const double &t) const{
	if (Ptype==WELL) {return Q;}
	else             {return 0.0;} 
} 
/***********************************************************************
				Get Integrated Budget Function
----------------------------------------------------------------------*/
void CPointElem::GetIntegratedBudget(const cmplex &z1, const cmplex &z2, const cmplex &z3, 
																			const double &t, double &inQ, double &outQ) const{

	inQ=outQ=0.0;
	if(Ptype==WELL){
		if (InTriangle(zp,z1,z2,z3)){
			if (Q>0.0) {outQ=fabs(Q); inQ =0.0;}//discharge well
			else       {inQ =fabs(Q); outQ=0.0;}//injection well
		}
	}
}
/************************************************************************
      Get Flux Distribution
	Gets integrated extraction between points z1 and z2 along line
	Linearly Distributes these fluxes to endpoints Q1 and Q2
	Used for finite element source distribution
************************************************************************/
void CPointElem::GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const{
	Q1=Q2=0.0;
	if(Ptype==WELL){
		if      (abs(z1-zp)/abs(z2-z1)<NEAR_FEATURE){Q1= -Q;Q2=0.0;}
		else if (abs(z2-zp)/abs(z2-z1)<NEAR_FEATURE){Q1=0.0;Q2= -Q;}
	}
}
/************************************************************************
      Get Singular Strength
************************************************************************/
cmplex CPointElem::GetSingularStrength  (const cmplex &z, const double &t) const{
	if      (abs(z-zp)<NEAR_FEATURE){
		if      (Ptype==WELL)  {return cmplex(IOVER2PI*Q,0.0);}
		else if (Ptype==VORTEX){return cmplex(0.0,IOVER2PI*Q);}
	}
	return 0.0;
}
//-----------------------------------------------------------------------
void CPointElem::WriteItself(ofstream &SOL, const double &t) const{
	SOL<< name<<endl;
	SOL<< Q   <<endl;
}
//-----------------------------------------------------------------------
bool CPointElem::ReadItself(ifstream &SOL){
	int Len(0);
	char *s[MAXINPUTITEMS];
	bool warm(true);

	if ( !TokenizeLine(SOL,s,Len)) {} //skip name of element
	if ((!TokenizeLine(SOL,s,Len)) && (Len==1)) {
		if (Ptype!=WELL){	
			Q=s_to_d(s[0]);
		}
	} 
	else {warm=false;}
	return warm;
}

//##################################################################################################
/************************************************************************
*************************************************************************
      CDischargeWell Class
*************************************************************************
************************************************************************/

CDischargeWell::CDischargeWell(char   *Name,
															 CSingleLayerABC *pLay, 
															 cmplex  z,
															 double  prate):
							  CPointElem(Name,pLay,WELL,0.0,prate,z){
	given=true;
}
//***********************************************************************
void   CDischargeWell::SolveItself(double &change, double &objective,const double t){
	if ((pBlock!=NULL) && (pBlock->IsOn())){pBlock->Update(myBlockID,NOT_A_STRING,t);}
  change=0.0;
  objective=0.0;
}
//***********************************************************************
void CDischargeWell::WriteOutput(const double &t) const{}
//------------------------------------------------------------------------
double CDischargeWell::GetMaxError(const double &t) const{return 0.0;}
//------------------------------------------------------------------------
CDischargeWell *CDischargeWell::Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name){
  //string "Well", string name 
  //double x double y double Q {double r (ignored)}
  //&
	CDischargeWell *pWell;
  pWell=NULL;
  bool     eof(false);
	int      Len;
	char   *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Discharge Well"<<endl;}  
																								eof=TokenizeLine(input,s,Len); l++;
	do{
		if ((Len==4) || (Len==3)){
      pLay->UpdateExtents(s_to_c(s[0],s[1]));
 		  pWell=new CDischargeWell(Name, 
					 										 pLay,
														   s_to_c(s[0],s[1]),
		                           s_to_d(s[2]));
			                                          eof=TokenizeLine(input,s,Len); l++;
			}
	  else if (Len==0)                           {eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) {break;}
    else                                       {ImproperFormat(s,l); return NULL;}
	} while ((!strcmp(s[0],"&")) && (!eof));
	
	if (eof){return NULL;}
	else    {return pWell;}
}
//##################################################################################################
/***********************************************************************
***********************************************************************
                  CHeadSpecifiedWell Class
***********************************************************************
**********************************************************************/

/***********************************************************************
												CONSTRUCTORS
***********************************************************************/
CHeadSpecifiedWell::CHeadSpecifiedWell(char			 *Name, 
																			 CSingleLayerABC *pLay, 
																			 cmplex			zw, 
																			 double			Head,  
																			 double			r):
								    CPointElem(Name,pLay,WELL,0.0,0.0,zw){
	
	radius=r;   H=Head;        
}
/***********************************************************************
												ACCESSORS
***********************************************************************/
double CHeadSpecifiedWell::GetRadius() const {return radius;}
//-----------------------------------------------------------------------
double CHeadSpecifiedWell::GetHead()   const {return H;}

/***********************************************************************
												SolveItself
***********************************************************************/
void CHeadSpecifiedWell::SolveItself(double &change, double &objective,const double t){
  double TotalPot, LocalPot; 
  double DesiredPot,Qold,K,T,B;

  K=pLayer->GetCond (zp+radius);
  T=pLayer->GetThick(zp+radius);
	B=pLayer->GetBase (zp+radius);

  //identify the required potential for actual head (confined or unconf)
  TotalPot=pLayer->GetDischargePotential(zp+radius,t).real();
  LocalPot=GetDischargePotential(zp+radius,t).real();
  Qold=Q;

	DesiredPot=ConvertToPotential(H-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());

  Q=((DesiredPot-(TotalPot-LocalPot))*2.0*PI)/
    (log((radius)/abs(zp-(pLayer->GetBlackHole()))));

	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,NOT_A_STRING,t);
	}
  objective=0.0;             //always met exactly
  change=fabs((Q-Qold)/Q);
}
//-----------------------------------------------------------------------
double CHeadSpecifiedWell::GetMaxError(const double &t) const	{
	double B=pLayer->GetBase(zp+radius);
	return (pLayer->GetHead (zp+radius,t)-(H-B))/(H-B);
}
/***********************************************************************
												Read/write
***********************************************************************/
void CHeadSpecifiedWell::WriteOutput(const double &t) const{} 
//------------------------------------------------------------------------
CHeadSpecifiedWell *CHeadSpecifiedWell::Parse(ifstream &input, int &l, 
																							CSingleLayerABC *pLay, char * Name){
	//string "HSWell", string name
  //double x double y double Head double r
  //&  
	CHeadSpecifiedWell *pWell;
  pWell=NULL;
  bool     eof(false);
	int      Len;
	char   *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Head-Specified well"<<endl;}       
																								eof=TokenizeLine(input,s,Len); l++;
	do{
		if (Len==4){
      pLay->UpdateExtents(s_to_c(s[0],s[1]));
 		  pWell=new CHeadSpecifiedWell(Name, 
					 										     pLay,
														       s_to_c(s[0],s[1]),
		                               s_to_d(s[2]),
															     s_to_d(s[3]));
			                                          eof=TokenizeLine(input,s,Len); l++;
			}
		else if (Len==0)                           {eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) {break;}
    else                                       {ImproperFormat(s,l); return NULL;}
	} while ((!strcmp(s[0],"&")) && (!eof));
	
	if (eof){return NULL;}
	else    {return pWell;}
}
//##################################################################################################
//***********************************************************************
//***********************************************************************
//                  CDryWell Class
//***********************************************************************
//***********************************************************************

//***********************************************************************
//												CONSTRUCTORS
//***********************************************************************
CDryWell::CDryWell(char			 *Name, 
									 CSingleLayerABC *pLay, 
									 cmplex		  zw, 
									 double			prate,  
									 double			r):
					CPointElem(Name,pLay,WELL,0.0,prate,zw){
	radius=r;	
	Qpump=prate;
}
//***********************************************************************
//												ACCESSORS
//***********************************************************************
double CDryWell::GetRadius() const{return radius;}       
//-----------------------------------------------------------------------                                                                                                                               
double CDryWell::GetDischarge() const{return Q;}
//***********************************************************************
//												SolveItself
//***********************************************************************
void CDryWell::SolveItself(double &change, double &objective,const double t){
  double TotalPot, LocalPot; 
  double Qold,B;

	B=pLayer->GetBase(zp+radius);


  Qold=Q;
	Q=Qpump;

  //identify the potential at the radius of the well at maximum pumping rate, Qpump
  TotalPot=pLayer->GetDischargePotential(zp+radius,t).real();
  LocalPot=GetDischargePotential(zp+radius,t).real();

	if ((TotalPot)<=0){  
		//if well is pumping more water than exists, well extracts less or not at all (never adds water)
		Q=max((-(TotalPot-LocalPot)*2.0*PI)/log(radius/abs(zp-(pLayer->GetBlackHole()) ) ),0.0);
	}
	else{
		//otherwise, it pumps as much as it wants to
		Q=Qpump;
	}

	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,NOT_A_STRING,t);
	}

  objective=0.0;             //always met exactly
  change=fabs((Q-Qold)/Q);
}
//-----------------------------------------------------------------------
double CDryWell::GetMaxError(const double &t) const	{
	double H=pLayer->GetHead (zp+radius,t);
	if (H>0){return 0.0;}
	else    {
		double TotalPot=pLayer->GetDischargePotential(zp+radius,t).real();
    double LocalPot=GetDischargePotential(zp+radius,t).real();
		double Qtmp=max((-(TotalPot-LocalPot)*2.0*PI)/log(radius/abs(zp-(pLayer->GetBlackHole()) ) ),0.0);
		return (Q-Qtmp)/Qtmp;
	}
}
//***********************************************************************
//												Read/write
//***********************************************************************
void CDryWell::WriteOutput(const double &t) const{} 

//------------------------------------------------------------------------
CDryWell *CDryWell::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){
	//string "DryWell", string name
  //double x double y double Qpump double r
  //&  
	CDryWell *pWell;
  pWell=NULL;
  bool     eof(false);
	int      Len;
	char   *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Drying Well"<<endl;}  eof=TokenizeLine(input,s,Len); l++;
	do{
		if (Len==4){
      pLay->UpdateExtents(s_to_c(s[0],s[1]));
 		  pWell=new CDryWell(Name, 
					 							 pLay,
												 s_to_c(s[0],s[1]),
		                     s_to_d(s[2]),
												 s_to_d(s[3]));
			                                          eof=TokenizeLine(input,s,Len); l++;
		}
	  else if (Len==0){                           eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) {break;}
    else {cout <<"line"<< l << "is wrong length"<<endl;break;}
	} while ((!strcmp(s[0],"&")) && (!eof));
	
	if (eof){return NULL;}
	else    {return pWell;}
}
//##################################################################################################

//************************************************************************
//************************************************************************
//      CVortex Class
//************************************************************************
//************************************************************************
CVortex::CVortex(char			 *Name, 
								 CSingleLayerABC *pLay, 
								 cmplex		  zv, 
								 double			strength):
				 CPointElem(Name,pLay,VORTEX,0.0,strength,zv){
}
//***********************************************************************
double CVortex::GetStrength() const{return Q;}
//***********************************************************************
void   CVortex::SolveItself(double &change, double &objective,const double t){
	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,NOT_A_STRING,t);
	}
  change=objective=0.0;
}
//***********************************************************************
CVortex *CVortex::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){
  //string "Vortex" , string name 
	//double x double y double strength
	//&
	CVortex *pVortex;
  pVortex=NULL;
  bool     eof(false);
	int      Len;
	char   *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "VORTEX"<<endl;}    eof=TokenizeLine(input,s,Len); l++;
	do{
		if (Len==4){
      pLay->UpdateExtents(s_to_c(s[0],s[1]));
 		  pVortex=new CVortex(Name, 
					 								pLay,
													s_to_c(s[0],s[1]),
		                      s_to_d(s[2]));
			                                          eof=TokenizeLine(input,s,Len); l++;
			}
	  else if (Len==0){                           eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) {break;}
    else {cout <<"line"<< l << "is wrong length"<<endl;break;}
	} while ((!strcmp(s[0],"&")) && (!eof));
	
	if (eof){return NULL;}
	else    {return pVortex;}
}
//------------------------------------------------------------------------
void CVortex::WriteOutput(const double &t) const{}
//------------------------------------------------------------------------
double CVortex::GetMaxError(const double &t) const{return 0.0;}
//##################################################################################################

//************************************************************************
//************************************************************************
//      CPtDipole Class
//************************************************************************
//************************************************************************
CPtDipole::CPtDipole(char			 *Name, 
										 CSingleLayerABC *pLay, 
										 cmplex		  z, 
										 double			strength, 
										 double			orientation):
				 CPointElem(Name,pLay,PTDIPOLE,orientation,strength,z){
}
//***********************************************************************
double CPtDipole::GetStrength()    const {return Q;}
//-----------------------------------------------------------------------
double CPtDipole::GetOrientation() const {return angle;}
//***********************************************************************
void   CPtDipole::SolveItself(double &change, double &objective,const double t){
	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,NOT_A_STRING,t);
	}
  change=0.0;
  objective=0.0;
}
//***********************************************************************
CPtDipole *CPtDipole::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){
  //string "Dipole" , string name 
	//double x double y double strength double orientation
	//&
	CPtDipole *pDipole;
  pDipole=NULL;
  bool     eof(false);
	int      Len;
	char    *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Point Dipole"<<endl;}    
																							  eof=TokenizeLine(input,s,Len); l++;
	do{
		if (Len==4){
      pLay->UpdateExtents(s_to_c(s[0],s[1]));
 		  pDipole=new CPtDipole(Name, 
					 								  pLay,
													  s_to_c(s[0],s[1]),
														s_to_d(s[2]),
														s_to_d(s[3]));
			                                          eof=TokenizeLine(input,s,Len); l++;
			}
	  else if (Len==0){                           eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) {break;}
    else {cout <<"line"<< l << "is wrong length"<<endl;break;}
	} while ((!strcmp(s[0],"&")) && (!eof));

	if (eof){return NULL;}
	else    {return pDipole;}
}

//------------------------------------------------------------------------
void CPtDipole::WriteOutput(const double &t) const{}
//------------------------------------------------------------------------
double CPtDipole::GetMaxError(const double &t) const{return 0.0;}

