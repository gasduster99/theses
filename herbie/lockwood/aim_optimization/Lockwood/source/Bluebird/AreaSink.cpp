#include "AreaSink.h"

/************************************************************************
                           CASinkDoublet
************************************************************************/
CASinkDoublet::CASinkDoublet(char									 *Name,
														 const CSingleLayerABC *pLay,
														 const cmplex					 *points,
														 const int							NumOfLines,
														 const int							prec,
													         CAreaSink			 *pASink)
	            :CStringElem(Name,pLay,points,NumOfLines,true,DOUBLET,prec){
	pAreaSink=pASink;
}
//-----------------------------------------------------------------------
void CASinkDoublet::SolveItself(double &change,double &objective,const double t){
  int i,m;
  double  maxchange(0.0),maxobjective(0.0);

  for (i=0; i<NLines; i++){
    
    //solve corresponding doublet string (removes jump in potential at boundary)
    //*********************************************************
    
		//set RHS for system of equations
    for (m=0; m<nlinecontrol; m++){
      rhs[m]=-pAreaSink->GetInteriorPotential(zctrl[i][m],t).real(); //-om.real();
    }
   
    //solve for doublet jump coefficients
		SetConstraints(i,1,0,0,0,UNCONSTRAINED,0.0,0.0);
    GenSolve(JumpCoeff[i],i,ltype,false,1.0,objective,change);
		upperswap(maxchange,change);
		upperswap(maxobjective,objective);

		SetFarFieldCoeff(i);

	} 
	change=maxchange;
	objective=maxobjective;
}
//-----------------------------------------------------------------------
double CASinkDoublet::GetMaxError          (const double &t) const{return 0.0;}//TMP DEBUG
//-----------------------------------------------------------------------
void   CASinkDoublet::WriteOutput          (const double &t) const{}

/************************************************************************
                           CASinkLinesink
************************************************************************/
CASinkLinesink::CASinkLinesink(char									 *Name,
															 const CSingleLayerABC *pLay,
															 const cmplex					 *points,
															 const int							NumOfLines,
															 const int							prec,
													           CAreaSink			 *pASink)
	            :CStringElem(Name,pLay,points,NumOfLines,true,LINESINK,prec){
	pAreaSink=pASink;
}
//-----------------------------------------------------------------------
void CASinkLinesink::SolveItself(double &change,double &objective,const double t){

  double  maxchange(0.0),maxobjective(0.0);
	double  length;
  double  mu,mu1,mu2,muadd,Yc,Xn,Yn,term1,term2;//,r;
  cmplex  Zn,Z,om,junk1;
  int     i,m,n; 

  mu=mu1=mu2=muadd=0.0;
  for (i=0; i<NLines; i++){
    //obtain jump at node points, additive constant
    length=abs(ze[i+1]-ze[i]);
    mu1=mu2; muadd=0;

    Z=(Centroid()-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
    Yc=Z.imag();

		//obtain constraint muadd, mu2 (mu2-mu1)=Net flux from side?
    muadd+=pAreaSink->Relax*pAreaSink->ave_leak*(-pow(length,2)*Yc/8.0);
    mu2+=  pAreaSink->Relax*pAreaSink->ave_leak*(-pow(length,2)*Yc/4.0);

    for (n=0; n<pAreaSink->MQorder; n++){
      Zn=(pAreaSink->zMQbasis[n]-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
      Xn=Zn.real();
      Yn=Zn.imag();
      if (Yn!=0){
				term1=sqrt(pow(-1-Xn,2)+pow(Yn,2));
				term2=sqrt(pow( 1-Xn,2)+pow(Yn,2));
				muadd+=pAreaSink->MQcoeff[n]*Yn*pow(0.5*length,3)/6.0*
						  (  (-1-Xn)*term1 + pow(Yn,2)*log(term1-Xn-1));
				mu2+=  pAreaSink->MQcoeff[n]*Yn*pow(0.5*length,3)/6.0*
							( ((-1-Xn)*term1 + pow(Yn,2)*log(term1-Xn-1))-
								(( 1-Xn)*term2 + pow(Yn,2)*log(term2-Xn+1)));
      }
    }	

    //set RHS for system of equations
    for (m=0; m<nlinecontrol; m++){
      mu=mu1+muadd+pAreaSink->Relax*pAreaSink->ave_leak*(-pow(length,2)*Yc*X[m]/8);
		  for (n=0; n<pAreaSink->MQorder; n++){
				Zn=(pAreaSink->zMQbasis[n]-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
        Xn=Zn.real();
				Yn=Zn.imag();
				term1=sqrt(pow(X[m]-Xn,2)+pow(Yn,2));
				if (Yn!=0){
					mu-=pAreaSink->MQcoeff[n]*Yn*pow(0.5*length,3)/6*
							((X[m]-Xn)*term1 + pow(Yn,2)*log((X[m]-Xn)+term1));
				}
      }
      rhs[m]=mu;
    }

    //solve for linesink jump coefficients
		SetConstraints(i,1,0,0,0,ASSIGNED,mu1,mu2);
    GenSolve(JumpCoeff[i],i,LINESINK,false,1.0,objective, maxchange);
		upperswap(maxobjective,objective);

    SetFarFieldCoeff(i);

	} 
	change=maxchange;
	objective=maxobjective;
}
//-----------------------------------------------------------------------
double CASinkLinesink::GetMaxError          (const double &t) const{return 0.0;}//TMP DEBUG
//-----------------------------------------------------------------------
void   CASinkLinesink::WriteOutput          (const double &t) const{}


/************************************************************************
                           CAreaSink
************************************************************************/
//Constructor
CAreaSink::CAreaSink(char									 *Name, 
										 const CSingleLayerABC *pLay,
										 const cmplex					 *points,
										 const int							NumOfLines,
		                 const double					 *leakage,
										 const cmplex    			 *givenlkpts,
										 const int        			prec,
										 const int							numlkpts)
	         :CAnalyticElem(Name,pLay){
	int m;

  //CAreaSink Initializers
  //allocate memory for dynamic arrays
	nleakctrl=numlkpts; 
  leakctrl=     new double   [nleakctrl];
  zleakctrl=    new cmplex   [nleakctrl]; 

	MQorder  =numlkpts;    //these are the same - want to meet leakage exactly @ points (collocation) 
  zMQbasis=     new cmplex   [MQorder  ];
  MQcoeff=      new double   [MQorder  ]; 

  //copy leakage array and ctrlpt array, initialize MQorder and ave_leak
	ave_leak=0.0;
  for(m=0;m<nleakctrl;m++){
		leakctrl  [m]=leakage[m]; 
	  zleakctrl [m]=givenlkpts[m]; 

	  ave_leak+=leakctrl[m]/double(nleakctrl);

	  MQcoeff   [m]=0.0;
	}

  //must create areal cont. point loc.
	//(this is the given part - control same as basis: different for AreaSink subclasses)
  for(m=0;m<nleakctrl;m++){zMQbasis[m]=zleakctrl[m];}

	//Create boundary elements
  pLinesinkBoundary=new CASinkLinesink(name,pLay,points,NumOfLines,prec,this);
  pDoubletBoundary =new CASinkDoublet  (name,pLay,points,NumOfLines,prec,this);

}
//------------------------------------------------------------------------
CAreaSink::CAreaSink(char									 *Name, 
										 const CSingleLayerABC *pLay,
										 const cmplex					 *points,
										 const int							NumOfLines,
										 const int							prec)
	        :CAnalyticElem(Name,pLay){

  //CAreaSink Initializers
	ave_leak =0.0;   
	nleakctrl=0;
  leakctrl =NULL;
  zleakctrl=NULL;

	MQorder  =0;
  zMQbasis =NULL;
  MQcoeff  =NULL;

	pLinesinkBoundary=new CASinkLinesink(name,pLay,points,NumOfLines,prec,this);
  pDoubletBoundary =new CASinkDoublet (name,pLay,points,NumOfLines,prec,this);
}
//------------------------------------------------------------------------
CAreaSink::~CAreaSink(){
	if (globaldebug){cout <<"   DESTROYING AREASINK "<<name<<endl;}
  delete [] leakctrl;
	delete [] zleakctrl;
  delete [] zMQbasis;
  delete [] MQcoeff;
}
//------------------------------------------------------------------------
void CAreaSink::SetBlockOwner(COwnerABC *BlockPtr, int seg, int IDinBlock){
	pBlock   =BlockPtr;
	myBlockID=IDinBlock;
}
//------------------------------------------------------------------------
void CAreaSink::UpdateBlock(const double &t) const{
	CAnalyticElem    ::UpdateBlock(t);
  pLinesinkBoundary->UpdateBlock(t);
	pDoubletBoundary ->UpdateBlock(t);
}
//-----------------------------------------------------------------------
bool   CAreaSink::HasFlux() const {return true;}
/************************************************************************
                          STATIC MEMBERS
************************************************************************/
int    CAreaSink::RelaxStartIter=0;
int    CAreaSink::RelaxEndIter=0;
double CAreaSink::Relax=1.0;
//-----------------------------------------------------------------------
void   CAreaSink::SetRelaxation(const int start, const int end){
	if (RelaxEndIter>RelaxStartIter){
		RelaxStartIter=start;
		RelaxEndIter=end;
	}
}
/************************************************************************
                          SolveItself
************************************************************************/
void CAreaSink::SolveItself(double &change, double &objective,const double t){
  bool    solved;
  double  maxobjective,maxchange;

	solved=CalcLeakage(t);

  if(!solved){  

		SolveMQCoeff(t);

		maxobjective=maxchange=0.0;

		pLinesinkBoundary->SolveItself(maxchange,maxobjective,t);
		upperswap(change,maxchange);
		upperswap(objective,maxobjective);

		pDoubletBoundary ->SolveItself(maxchange,maxobjective,t);
		upperswap(change,maxchange);
		upperswap(objective,maxobjective);


  } //if (!solved)...

	if ((pBlock!=NULL) && (pBlock->IsOn())){pBlock->Update(myBlockID,NOT_A_STRING,t);}

	objective=change=0.0;
}
//------------------------------------------------------------------------
double CAreaSink::GetMaxError          (const double &t) const{return 0.0;}
//************************************************************************
bool CAreaSink::CalcLeakage(const double t){
	if (CAreaSink::RelaxEndIter>0){ 
		int thisiter=pLayer->GetCurrentIter();
		Relax=(double)(min(thisiter-RelaxStartIter,0))/(RelaxEndIter-RelaxStartIter);
		return true;
	}
	else{
		return false;
	}
}  //does nothing for given leakage
//************************************************************************
void CAreaSink::SolveMQCoeff(const double t){
	//can use general routine
	  MQInterpolate(leakctrl,zleakctrl,nleakctrl,
								  MQcoeff ,zMQbasis ,MQorder,0.0,
								  ave_leak);
}
//************************************************************************
cmplex CAreaSink::GetInteriorPotential (const cmplex &z,const double &t) const{
	double r;	
  int    n;
	cmplex omega(0.0);

	r=abs(z-pLinesinkBoundary->Centroid());
	omega=0.25*Relax*ave_leak*r*r;
	for(n=0; n<MQorder;n++){
		r=abs(z-zMQbasis[n]);
		omega+=MQcoeff[n]*r*r*r/9;
	}
	return omega;
}
//************************************************************************
cmplex CAreaSink::GetDischargePotential(const cmplex &z,const double &t) const{

	cmplex omega(0.0);

	//from leakage
  if (pLinesinkBoundary->IsInside(z)){
		omega+=GetInteriorPotential(z,t);
	}
  omega+=pLinesinkBoundary->GetDischargePotential(z,t);
  omega+=pDoubletBoundary ->GetDischargePotential(z,t);

  return omega;
  
}
//************************************************************************
cmplex CAreaSink::GetW(const cmplex &z,const double &t) const{
	
	cmplex W(0.0);
	double r;

  if (pLinesinkBoundary->IsInside(z)){

    W=-0.5*Relax*ave_leak*conj(z-pLinesinkBoundary->Centroid());
    for(int n=0; n<MQorder;n++){
      r=abs(z-zMQbasis[n]);
      W-=MQcoeff[n]*r*conj(z-zMQbasis[n])/3.0;
		}
	}

  W+=pLinesinkBoundary->GetW(z,t);
  W+=pDoubletBoundary ->GetW(z,t);

  return W;
}
//************************************************************************
cmplex CAreaSink::GetGx(const cmplex &z,const double &t) const{
	cmplex Gx(0.0);
	int    n;

  if (pLinesinkBoundary->IsInside(z)){
		double r;
		Gx=-0.5*Relax*ave_leak; //This works 
    for(n=0; n<MQorder;n++){
      r=abs(z-zMQbasis[n]);
      Gx-=1.0/6.0*MQcoeff[n]*(3.0*r+r*conj(z-zMQbasis[n])/(z-zMQbasis[n])); 
		}
	}

	Gx+=pLinesinkBoundary->GetGx(z,t);
  Gx+=pDoubletBoundary ->GetGx(z,t);

  return Gx;
}
//************************************************************************
double CAreaSink::GetLeakage(const cmplex &z, const double &t,const leak_type ltype) const{ 
	//assumes that recharge is from top, leakage is from bottom
  double leakage(0.0);
  if (pLinesinkBoundary->IsInside(z)){
    if (MQorder<=1){
			leakage=Relax*ave_leak;
		}
    else {
      for(int n=0; n<MQorder;n++){
				leakage+=MQcoeff[n]*abs(z-zMQbasis[n]);
			}
      leakage+=Relax*ave_leak;
		}
		if (ltype==FROMTOPANDBOTTOM){return leakage;       }
		if (ltype==FROMBOTTOM)      {return max(0.0,leakage);} //returns zero if recharge asink
		if (ltype==FROMTOP)         {return min(0.0,leakage);} //returns zero if leakage  asink
	}
  return 0.0;
}  
/************************************************************************
GetDivFluxThruFace
returns integrated normal discharge through face defined by z1 and z2
WARNING- currently cannot handle more than one intersection with Area sink boundary
************************************************************************/
double CAreaSink::GetDivFluxThruFace   (const cmplex &z1, const cmplex &z2, const double &t) const{
	
	//double length(abs(z2-z1));	

  double  Yc,Xn,Yn,term1,term2;
  cmplex  Zn,Z,z1a,z2a,v1,v2,zint;
  int     n;// ,m
  double  Q=0.0;
	bool    in1,in2;
	
	in1=pLinesinkBoundary->IsInside(z1);
	in2=pLinesinkBoundary->IsInside(z2);

	if      (in1 && in2)      {z1a=z1;z2a=z2;}
	else if ((!in1) && (!in2)){return 0.0;} 
	else if (!(in1 && in2))   {
		//find intersection point (only one allowed)
	  int N=pLinesinkBoundary->GetNumSegs();
		for (int i=0; i<N; i++){
			if (i!=0)  {v1=pLinesinkBoundary->GetZArray()[i-1];}
			else       {v1=pLinesinkBoundary->GetZArray()[N-1];}
			v2=pLinesinkBoundary->GetZArray()[i ];

			//Get intersections of side with transect (3 trys per side)
			if (Intersect(v1,v2,z1,z2,zint)==INTERSECTED){
				if (abs(zint-z1)/abs(z2-z1)==0.0){return 0.0;}
				if (abs(zint-z1)/abs(z2-z1)==1.0){z1a=z1;z2a=z2;}//handles disconnect of intersect and isinside
				if (!in1){z1a=zint;z2a=z2;}
				else     {z2a=zint;z1a=z1;} 
				//cout <<"intersect: "<<abs(zint-z1)/abs(z2-z1)<<endl;
				break;
			}
		}
	}

	double length(abs(z2a-z1a));
	if (length==0.0){return 0.0;}

	Z=(Centroid()-0.5*(z1a+z2a))/(0.5*(z2a-z1a));
  Yc=Z.imag();

	//obtain net flux from side analytically
  Q=Relax*ave_leak*(-pow(length,2)*Yc/4.0);

  for (n=0; n<MQorder; n++){
    Zn=(zMQbasis[n]-0.5*(z1a+z2a))/(0.5*(z2a-z1a));
    Xn=Zn.real();
    Yn=Zn.imag();
    if (Yn!=0){
		  term1=sqrt(pow(-1-Xn,2)+pow(Yn,2));
			term2=sqrt(pow( 1-Xn,2)+pow(Yn,2));
			Q+=MQcoeff[n]*Yn*pow(0.5*length,3)/6.0*
				 ( ((-1-Xn)*term1 + pow(Yn,2)*log(term1-Xn-1))-
					 (( 1-Xn)*term2 + pow(Yn,2)*log(term2-Xn+1)));
    }
  }	

	/*cmplex unitY=IM*(z2a-z1a)/abs(z2a-z1a);
	int    numdiv=10;
  cmplex W(0.0);
	double r;
  double Q2(0.0);
	
	//Numerical Integration
	for (int m=0; m<numdiv; m++){
		zint=((double)(m)+0.5)/(double)(numdiv)*(z2a-z1a)+z1a;

		W=-0.5*Relax*ave_leak*conj(zint-pLinesinkBoundary->Centroid());
		for(int n=0; n<MQorder;n++){
			r=abs(zint-zMQbasis[n]);
			W-=MQcoeff[n]*r*conj(zint-zMQbasis[n])/3.0;
		}
		
		Q2-=(conj(unitY)*conj(W)).real()*(length/(double)(numdiv));
	}*/
	//}
	//cout <<"CAreaSink::GetDivFluxThruFace: "<<Q<<" "<<Q2<<endl;

	return Q;
}
/************************************************************************
                           GetFluxThruFace
************************************************************************/
cmplex CAreaSink::GetFluxThruFace   (const cmplex &z1, const cmplex &z2, const double &t) const{
	cmplex Qnet(0.0);
	Qnet+=GetDivFluxThruFace(z1,z2,t); //Normal integrated flux
	//Qnet+=IM*(GetInteriorPotential(z1,t).real()-GetInteriorPotential(z2,t).real()) //tangential integrated flux (only valid inside element)
	Qnet+=pLinesinkBoundary->GetFluxThruFace(z1,z2,t);
  Qnet+=pDoubletBoundary ->GetFluxThruFace(z1,z2,t);
	return Qnet;
}
/************************************************************************
                           GetIntegratedLeakage
Gets integrated leakage [L^3/T] into CLOCKWISE triangle defined by z1,z2,z3
Assumes that largest triangle side is smaller than smallest Area Sink side
Also assumes distance between any two ASink vertices greater than smallest Asink side
*Always?* works on convex polygons that meet these conditions
************************************************************************/
double CAreaSink::GetIntegratedLeakage (const cmplex    &z1, 
																				const cmplex    &z2, 
																				const cmplex    &z3, 
																				const double    &t, 
																				const leak_type  ltype) const{

  int    i;
	int    nInAsink(0),nInTriangle(0);
	cmplex zint1,zint2,zint,v1,v2,v3,zin1(0.0),zin2,zasinkin,ztmp,Z1,Z2;
	bool   in[3], convex;
	double Q;
	double A=TriArea(z1,z2,z3);

	if (A >0.0){ExitGracefully("CAreaSink::GetIntegratedLeakage: points should be clockwise",RUNTIME_ERR);}
  if (A==0.0){return 0.0;}

	if (ave_leak==0){ return 0.0;}
	//assumes that recharge is from top, leakage is from bottom
	if ((ltype==FROMTOP   ) && (ave_leak>0)){return 0.0;}
	if ((ltype==FROMBOTTOM) && (ave_leak<0)){return 0.0;}

	int    N=pLinesinkBoundary->GetNumSegs();

	//Check for triangle points inside area sink
	in[0]=pLinesinkBoundary->IsInside(z1); if (in[0]){nInAsink++;if (zin1!=0.0){/*BAD*/ }else{zin1=z1;}}
	in[1]=pLinesinkBoundary->IsInside(z2); if (in[1]){nInAsink++;if (zin1!=0.0){zin2=z2;}else{zin1=z2;}}
	in[2]=pLinesinkBoundary->IsInside(z3); if (in[2]){nInAsink++;if (zin1!=0.0){zin2=z3;}else{zin1=z3;}}	
	
	//Check for area sink points inside triangle
	for (i=0; i<N; i++){
		if (i!=0)  {v1=pLinesinkBoundary->GetZArray()[i-1];}else{v1=pLinesinkBoundary->GetZArray()[N-1];}
		if (i!=N-1){v3=pLinesinkBoundary->GetZArray()[i+1];}else{v3=pLinesinkBoundary->GetZArray()[0];}
		v2=pLinesinkBoundary->GetZArray()[i ];

		if (InTriangle(v2,z1,z2,z3)){
			//Get intersections of sides with sides (6 trys)
			if      (Intersect(v1,v2,z1,z2,zint)==INTERSECTED){zint1=zint-(zint-v2)*0.001;}
			else if (Intersect(v1,v2,z2,z3,zint)==INTERSECTED){zint1=zint-(zint-v2)*0.001;}
			else if (Intersect(v1,v2,z3,z1,zint)==INTERSECTED){zint1=zint-(zint-v2)*0.001;}

			if      (Intersect(v2,v3,z1,z2,zint)==INTERSECTED){zint2=zint-(zint-v2)*0.001;}
			else if (Intersect(v2,v3,z2,z3,zint)==INTERSECTED){zint2=zint-(zint-v2)*0.001;}
			else if (Intersect(v2,v3,z3,z1,zint)==INTERSECTED){zint2=zint-(zint-v2)*0.001;}

			zasinkin=v2;
			nInTriangle++;
		}
	}

	
	if (nInTriangle==0){ //Triangle does not overlap an area sink corner====================================
		
		if      (nInAsink==0){ //No triangle/asink overlap----------------------------------------------------
			return 0.0;
		}

		else if (nInAsink==3){ //Triangle completely enclosed-------------------------------------------------
			return GetDivFluxThruFace(z1,z2,t)+
				     GetDivFluxThruFace(z2,z3,t)+
				     GetDivFluxThruFace(z3,z1,t);
		}

		else if (nInAsink>0){  //Triangle overlaps single side------------------------------------------------ 
			//cout <<"a    "<<nInAsink<<endl;
			
			zint1=0.0;
			for (i=0; i<N; i++){
				if (i!=0) {v1=pLinesinkBoundary->GetZArray()[i-1];}
				else      {v1=pLinesinkBoundary->GetZArray()[N-1];}
				v2=pLinesinkBoundary->GetZArray()[i ];

					//Get intersections of sides with sides (3 trys per side)
				if (Intersect(v1,v2,z1,z2,zint)==INTERSECTED){
					if  (zint1!=0.0){/*BAD*/                      }
					else            {zint1=zint-(zint-zin1)*0.001;}
				}
				if (Intersect(v1,v2,z2,z3,zint)==INTERSECTED){
					if (zint1!=0.0) {zint2=zint-(zint-zin1)*0.001;}
					else            {zint1=zint-(zint-zin1)*0.001;}
				}
				if (Intersect(v1,v2,z3,z1,zint)==INTERSECTED){
					if (zint1!=0.0) {zint2=zint-(zint-zin1)*0.001;}
					else            {/*BAD*/                      }
				}

				if (zint1!=0.0){break;} //intersection found, cannot intersect another line
			}

			zint1   =zint1   -(zint1   -zint2)*0.001;
			zint2   =zint2   -(zint2   -zint1)*0.001; //pull in so that faces are definitely inside asink

			if (TriArea(zint1,zint2,zin1)>0.0){ztmp=zint1;zint1=zint2;zint2=ztmp;}//makes clockwise
			
			return GetDivFluxThruFace(zint1,zint2,t)+
				     GetDivFluxThruFace(zint2,zin1,t)+
				     GetDivFluxThruFace(zin1, zint1,t);
			
		}
	}

	else if (nInTriangle==1){ //Triangle DOES overlap an area sink corner====================================

		if      (nInAsink==0){  //One subtriangle-------------------------------------------------------------- 
			//cout <<"c    "<<endl;
			zasinkin=zasinkin-(zasinkin-(zint1+zint2)/2.0 )*0.001;
			zint1   =zint1   -(zint1   -zint2)*0.001;
			zint2   =zint2   -(zint2   -zint1)*0.001; //pull in so that faces are definitely inside asink

			if (TriArea(zint1,zint2,zasinkin)>0.0){ztmp=zint1;zint1=zint2;zint2=ztmp;}//makes clockwise

			return  GetDivFluxThruFace(zint1    ,zint2   ,t)+
				      GetDivFluxThruFace(zint2    ,zasinkin,t)+
				      GetDivFluxThruFace(zasinkin ,zint1   ,t);	
		}

		else if (nInAsink==1){ //Two subtriangles with shared side---------------------------------------
			//cout <<"d    "<<endl;
			
			//Check for convex or concave asink corner
			Z1=(zasinkin-0.5*(zint1+zint2))/(0.5*(zint2-zint1));
			Z2=(zin1    -0.5*(zint1+zint2))/(0.5*(zint2-zint1));
			convex=(Z1.imag()*Z2.imag()>0.0);

			zasinkin=zasinkin-(zasinkin-zin1 )*0.001;//pull in so that faces are definitely inside asink
			if (convex){
				if (TriArea(zint1,zint2,zasinkin)>0.0){ztmp=zint1;zint1=zint2;zint2=ztmp;}
				zint1   =zint1   +(zint1   -zint2)*0.001;
				zint2   =zint2   +(zint2   -zint1)*0.001;//push out so that faces are definitely inside asink			
			}
			else       {
				if (TriArea(zint1,zint2,zasinkin)<0.0){ztmp=zint1;zint1=zint2;zint2=ztmp;}
				zint1   =zint1   -(zint1   -zint2)*0.001;
				zint2   =zint2   -(zint2   -zint1)*0.001;//pull in so that faces are definitely inside asink
			}

			//if (TriArea(zint1,zasinkin,zin1)>0.0){cout <<"CAreaSink::GetIntegratedLeakage:Something wrong"<<endl;}
			//if (TriArea(zint2,zin1,zasinkin)>0.0){cout <<"CAreaSink::GetIntegratedLeakage:Something wrong"<<endl;}

			Q= GetDivFluxThruFace(zint1    ,zasinkin,t)+
				 GetDivFluxThruFace(zasinkin ,zin1    ,t)+//*cancels
				 GetDivFluxThruFace(zin1     ,zint1   ,t)+
				 GetDivFluxThruFace(zint2    ,zin1    ,t)+
				 GetDivFluxThruFace(zin1     ,zasinkin,t)+//*cancels
				 GetDivFluxThruFace(zasinkin ,zint2   ,t);

			return Q;
		}

		else if (nInAsink==2){ //Three subtriangles -----------------------------------------------------------

			//cout <<"e    "<<endl;

			//Check for convex or concave asink corner
			Z1=(zasinkin-0.5*(zint1+zint2))/(0.5*(zint2-zint1));
			Z2=(zin1    -0.5*(zint1+zint2))/(0.5*(zint2-zint1));
			
			convex=(Z1.imag()*Z2.imag()>0.0);

			//Adjust points 
			if (TriArea(zin1 ,zin2 ,zasinkin)>0.0){ztmp=zin1;zin1=zin2;zin2=ztmp;}
      zasinkin=zasinkin-(zasinkin-(zin1+zin2)/2.0 )*0.001;
			if (convex){
        if (TriArea(zint1,zint2,zasinkin)>0.0){ztmp=zint1;zint1=zint2;zint2=ztmp;}
				zint1   =zint1   +(zint1   -zint2)*0.001;
				zint2   =zint2   +(zint2   -zint1)*0.001;//push out so that faces are definitely inside asink			
			}
			else       {
        if (TriArea(zint1,zint2,zasinkin)<0.0){ztmp=zint1;zint1=zint2;zint2=ztmp;}
				zint1   =zint1   -(zint1   -zint2)*0.001;
				zint2   =zint2   -(zint2   -zint1)*0.001;//pull in so that faces are definitely inside asink
			}

			//if (TriArea(zin2,zint1,zasinkin)>0.0){cout <<"CAreaSink::GetIntegratedLeakage:Something wrong"<<endl;}
			//if (TriArea(zint2,zin1,zasinkin)>0.0){cout <<"CAreaSink::GetIntegratedLeakage:Something wrong"<<endl;}

			Q= GetDivFluxThruFace(zin1    ,zin2    ,t)+
				 GetDivFluxThruFace(zin2    ,zasinkin,t)+//*cancels
				 GetDivFluxThruFace(zasinkin,zin1    ,t);//*cancels
			Q+=GetDivFluxThruFace(zin2 ,   zint1   ,t)+
				 GetDivFluxThruFace(zint1  , zasinkin,t)+
				 GetDivFluxThruFace(zasinkin,zin2    ,t);//*cancels
			Q+=GetDivFluxThruFace(zint2   ,zin1    ,t)+
				 GetDivFluxThruFace(zin1    ,zasinkin,t)+//*cancels
				 GetDivFluxThruFace(zasinkin,zint2   ,t);			

			return Q;			
		}
	}
	else {//Multiple area sink points inside triangle===============================================
		ExitGracefully("CAreaSink::GetIntegratedLeakage: Area sink geometry too complex",RUNTIME_ERR);
		//Should perform numerical integration
	}

	return 0.0;

	
}
//************************************************************************
void CAreaSink::GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const{
	double Q=GetIntegratedLeakage(z1,z2,z3,t,FROMTOPANDBOTTOM);
	if (Q>0.0) {outQ=fabs(Q); inQ =0.0;}
	else       {inQ =fabs(Q); outQ=0.0;}
}
//************************************************************************
void   CAreaSink::WriteItself(ofstream &SOL, const double &t) const{

	pLinesinkBoundary->WriteItself(SOL,t);
  pDoubletBoundary ->WriteItself(SOL,t);
  
}
//----------------------------------------------------------------------
bool   CAreaSink::ReadItself(ifstream &SOL){
  
	bool warm(true);
	warm=pLinesinkBoundary->ReadItself(SOL);
  warm=pDoubletBoundary ->ReadItself(SOL);

	SolveMQCoeff(0.0);
  
	return warm;
}
//------------------------------------------------------------------------
void CAreaSink::WriteOutput(const double &t) const{}
cmplex CAreaSink::Centroid						 () const                               {return pLinesinkBoundary->Centroid();}
bool   CAreaSink::IsInSquare					 (const cmplex &zc,const double w) const{return pLinesinkBoundary->IsInSquare(zc,w);}
bool   CAreaSink::IsInCircle					 (const cmplex &zc,const double r) const{return pLinesinkBoundary->IsInCircle(zc,r);}
bool   CAreaSink::PartInCircle				 (const cmplex &zc,const double r) const{return pLinesinkBoundary->PartInCircle(zc,r);}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string AreaSink, string name 
	{double x double y}x(numlines+1) 
	& 
	{double x double y double leakage}x(numcontrolpts)
	&[int precision]
----------------------------------------------------------------------*/
CAreaSink *CAreaSink::Parse(ifstream &INPUT, int &l, 
														CSingleLayerABC *pLay, char * Name){

	CAreaSink *pAreaSink=NULL;
  bool       done(false);
	int        Len,thisprec,nlines(0),npoints(0),i;
  cmplex	   stringp[MAXLINES];       
  cmplex	   leakzp [MAXLINES];
  double     leakp  [MAXLINES];
	char      *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Area Sink element"<<endl;}  
	if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
  do {
		if (nlines>=MAXLINES) { ExitGracefully("CAreaSink::Parse- too many line segments in area sink",TOO_MANY);}
		if  (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]); nlines++; 
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
    else if ((Len>=1) && (!strcmp(s[0],"&"))) {
      stringp[nlines]=stringp[0];  nlines++;       	
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
			done=true; 
		}
    else {
			ImproperFormat(s,l); return NULL;
		}
	} while (!done);

	done=false;                                      
	do {
    if (npoints>=MAXLINES) { ExitGracefully("CAreaSink::Parse- too many control points in area sink",TOO_MANY);}
    if ((Len==3) && (strcmp(s[0],"&"))){
      leakzp[npoints]=s_to_c(s[0],s[1]);
			leakp [npoints]=s_to_d(s[2]);
			npoints++;                                   
			if (TokenizeLine(INPUT,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);}
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pAreaSink = new CAreaSink(Name,
				  											pLay,
				  											stringp,
																nlines-1,
				  											leakp,
				   											leakzp,
				  											thisprec,
				  											npoints);
			done=true;
		}

		else {ImproperFormat(s,l); return NULL;}
	} while (!done);

  return pAreaSink;
}
//************************************************************************
void CAreaSink::SetPrecision(const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("SetPrecision::Improper precision level specified",BAD_DATA);}
	switch(Precision){
		case(0): {order=1; fold=1.0; break;}
		case(1): {order=5; fold=1.2; break;}
		case(2): {order=5; fold=1.3; break;}
		case(3): {order=10;fold=1.4; break;}
	  case(4): {order=10;fold=1.4; break;}
	  case(5): {order=12;fold=1.5; break;}
	  case(9): {order=15;fold=1.5; break;}
	}
}


