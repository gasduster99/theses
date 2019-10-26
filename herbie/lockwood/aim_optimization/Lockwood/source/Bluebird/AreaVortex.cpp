#include "AreaVortex.h"
/************************************************************************
                           CAVorDipole
************************************************************************/
CAVorDipole::CAVorDipole(const CSingleLayerABC   *pLay,
										     const cmplex						 *points,
										     const int								NumOfLines,
										     const int								prec,
													     CAreaVortex			 *pAVor)
	            :CStringElem("AVor boundary 1",pLay,points,NumOfLines,true,DIPOLE,prec){
	pAreaVortex=pAVor;
}
//-----------------------------------------------------------------------
void CAVorDipole::SolveItself(double &change,double &objective,const double t){
  int i,m;
  double  maxchange(0.0),maxobjective(0.0);

  for (i=0; i<NLines; i++){
    
    //solve corresponding dipole string (removes jump in stream function at boundary)
    //*********************************************************
    
		//set RHS for system of equations
    for (m=0; m<nlinecontrol; m++){
      rhs[m]=-pAreaVortex->GetInteriorPotential(zctrl[i][m],t).imag(); //-om.real();
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
double CAVorDipole::GetMaxError          (const double &t) const{return 0.0;}//TMP DEBUG
//-----------------------------------------------------------------------
void   CAVorDipole::WriteOutput          (const double &t) const{}

/************************************************************************
                           CAVorLinevortex
************************************************************************/
CAVorLinevortex::CAVorLinevortex(const CSingleLayerABC *pLay,
															   const cmplex					 *points,
															   const int							NumOfLines,
															   const int							prec,
													             CAreaVortex		 *pAVor)
	            :CStringElem("AVor boundary 2",pLay,points,NumOfLines,true,LINEVORTEX,prec){
	pAreaVortex=pAVor;
}
//-----------------------------------------------------------------------
void CAVorLinevortex::SolveItself(double &change,double &objective,const double t){

  double  maxchange(0.0),maxobjective(0.0);
	int i,m;

  double  mu,mu1,mu2,muadd,Yc,Xn,Yn,term1,term2,length,chi1,chi2;//,r;
  cmplex  Zn,Z,om,junk1;
  int     n; 

  mu=mu1=mu2=muadd=0.0;

  for (i=0; i<NLines; i++){
    //obtain jump at node points, additive constant
    length=abs(ze[i+1]-ze[i]);
    mu1=mu2; muadd=0.0;

    Z=(Centroid()-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
    Yc=Z.imag();

    muadd+=pAreaVortex->ave_curl*(-length*length*Yc/8.0);
    mu2+=  pAreaVortex->ave_curl*(-length*length*Yc/4.0);
    
    for (n=0; n<pAreaVortex->MQorder; n++){
      Zn=(pAreaVortex->zMQbasis[n]-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
      Xn=Zn.real();
      Yn=Zn.imag();
			chi1=-1.0-Xn; //eq 37 in asink paper
			chi2= 1.0-Xn; //eq 43 in asink paper
      if (Yn!=0.0){
				term1=sqrt(chi1*chi1+Yn*Yn); 
				term2=sqrt(chi2*chi2+Yn*Yn);
				muadd+=pAreaVortex->MQcoeff[n]*Yn*pow(0.5*length,3)/6.0*
						  (  chi1*term1 + Yn*Yn*log(chi1+term1));             //eq 40 without \mu_j term at end
				mu2+=  pAreaVortex->MQcoeff[n]*Yn*pow(0.5*length,3)/6.0*
							( (chi1*term1 + Yn*Yn*log(chi1+term1))-
								(chi2*term2 + Yn*Yn*log(chi2+term2)));            //eq 40 plus eq 42 : \mu_{j+1}    
      }
    }	

    //set RHS for system of equations
    for (m=0; m<nlinecontrol; m++){
      mu=mu1+muadd+pAreaVortex->ave_curl*(-pow(length,2)*Yc*X[m]/8.0);
		  for (n=0; n<pAreaVortex->MQorder; n++){
				Zn=(pAreaVortex->zMQbasis[n]-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
        Xn=Zn.real();
				Yn=Zn.imag();
				term1=sqrt(pow(X[m]-Xn,2)+Yn*Yn);                     //X[m]-Xn is chi
				if (Yn!=0){
					mu-=pAreaVortex->MQcoeff[n]*Yn*pow(0.5*length,3)/6.0*
							((X[m]-Xn)*term1 + Yn*Yn*log((X[m]-Xn)+term1));     //-eq 40 without \mu_j term at end
				}
      }
      rhs[m]=mu;
    }

    //solve for linesink jump coefficients
		SetConstraints(i,1,0,0,0,ASSIGNED,mu1,mu2);
    GenSolve(JumpCoeff[i],i,ltype,false,1.0,objective,change);
		upperswap(maxobjective,objective);
    upperswap(maxchange,change);

    SetFarFieldCoeff(i);
	}
	//cout << "Net Curl?:"<<-mu2<<endl;
	cout <<"LV net curl (should be zero):"<<GetNetCurl(t)<<endl;

		
	/*double PotOther,Qxother;
	double Kin,Kout;
  for (i=0; i<NLines; i++){		
		cmplex unitX=(ze[i+1]-ze[i])/abs(ze[i+1]-ze[i]);
		
		disabled[i]=true;

		for (m=0; m<nlinecontrol; m++){
			Kin =pLayer->GetCond(zctrl[i][m]+(ze[i+1]-ze[i])*IM*MOVE_DIST);
			Kout=pLayer->GetCond(zctrl[i][m]-(ze[i+1]-ze[i])*IM*MOVE_DIST);

			Qxother=(conj(pLayer->GetW(zctrl[i][m],t))*unitX).real();

      PotOther=pLayer->GetDischargePotential(zctrl[i][m],t).real();
			//rhs[m]=2.0*(Kin-Kout)*max(PotOther,0.0)/(Kout+Kin); //-some function of dk/dx??
			rhs[m]=2.0*(Kin-Kout)*Qxother/(Kout+Kin);
			//cout <<rhs[m]<<endl;
		}

		disabled[i]=false;
		  SetConstraints(i,0,1,0,0,UNCONSTRAINED,0.0,0.0);
		//SetConstraints(i,1,0,0,0,UNCONSTRAINED,0.0,0.0);
		GenSolve(JumpCoeff[i],i,ltype,false,1.0,objective,change);
		upperswap(maxchange   ,change   );
		upperswap(maxobjective,objective);

    SetFarFieldCoeff(i);

	}*/
	
	change=maxchange;
	objective=maxobjective;
}
//-----------------------------------------------------------------------
double CAVorLinevortex::GetMaxError          (const double &t) const{return 0.0;}//TMP DEBUG
//-----------------------------------------------------------------------
void   CAVorLinevortex::WriteOutput          (const double &t) const{}

/************************************************************************
                           CAreaVortex
************************************************************************/
//Constructor
CAreaVortex::CAreaVortex(char									 *Name, 
												 const CSingleLayerABC *pLay,
											   const cmplex					 *points,
												 const int							NumOfLines,
												 const double					 *curl,
												 const cmplex					 *givenlkpts,
												 const int							prec,
												 const int							numlkpts)
						:CAnalyticElem(Name,pLay){
  int m;

  //CAreaVortex Initializers
 
	//allocate memory for dynamic arrays	
	ncurlctrl=numlkpts; 
  curlctrl =new double   [ncurlctrl];
  zcurlctrl=new cmplex   [ncurlctrl];
	ave_curl =0.0; 

	MQorder  =numlkpts;    //these are the same - want to meet curl exactly @ points (collocation) 
  zMQbasis= new cmplex   [MQorder  ];
  MQcoeff=  new double   [MQorder  ]; 

  //copy curl array and ctrlpt array, initialize MQorder and ave_curl
	ave_curl=0.0;
  for(m=0;m<ncurlctrl;m++){
		curlctrl  [m]=curl[m]; 
	  zcurlctrl [m]=givenlkpts[m]; 
	  ave_curl+=curlctrl[m]/double(ncurlctrl);
	  MQcoeff   [m]=0.0;
	}

  //must create areal cont. point loc.(this is the given part - control same as basis: different for AreaVortex subclasses)
  for(m=0;m<ncurlctrl;m++){zMQbasis[m]=zcurlctrl[m];}
 

	//Create boundary elements
  pLinevortexBoundary=new CAVorLinevortex(pLay,points,NumOfLines,prec,this);
  pDipoleBoundary    =new CAVorDipole    (pLay,points,NumOfLines,prec,this);

}
//------------------------------------------------------------------------
CAreaVortex::CAreaVortex(char						 *Name, 
												 const CSingleLayerABC *pLay,
												 const cmplex		 *points,												
												 const int			  NumOfLines, 
												 const int			  prec)
						:CAnalyticElem(Name,pLay){

  //CAreaVortex Initializers
	ave_curl =0.0;   
	ncurlctrl=0;
  curlctrl =NULL;
  zcurlctrl=NULL;

	MQorder  =0;
  zMQbasis =NULL;
  MQcoeff  =NULL;

	//Create boundary elements
  pLinevortexBoundary=new CAVorLinevortex(pLay,points,NumOfLines,prec,this);
  pDipoleBoundary    =new CAVorDipole    (pLay,points,NumOfLines,prec,this);
}
//------------------------------------------------------------------------
CAreaVortex::~CAreaVortex(){
	if (globaldebug){cout <<"   DESTROYING AREA VORTEX "<<name<<endl;}
  delete [] curlctrl;
	delete [] zcurlctrl;

  delete [] MQcoeff;
  delete [] zMQbasis;
}
//------------------------------------------------------------------------
void CAreaVortex::SetBlockOwner(COwnerABC *BlockPtr, int seg, int IDinBlock){
	pBlock   =BlockPtr;
	myBlockID=IDinBlock;
}
//------------------------------------------------------------------------
void CAreaVortex::UpdateBlock(const double &t) const{
	CAnalyticElem::UpdateBlock(t);
	pLinevortexBoundary->UpdateBlock(t);
  pDipoleBoundary    ->UpdateBlock(t);
}
/************************************************************************
                          SolveItself
************************************************************************/
void CAreaVortex::SolveItself(double &change, double &objective,const double t){
  bool    solved(false);
  double  maxobjective(0.0),maxchange(0.0);

	solved=CalcCurl(t);

  if(!solved){  

	  SolveMQCoeff(t);
		
		maxobjective=maxchange=0.0;

		pLinevortexBoundary->SolveItself(maxchange,maxobjective,t);
		upperswap(change,maxchange);
		upperswap(objective,maxobjective);

		pDipoleBoundary->SolveItself(maxchange,maxobjective,t);
		upperswap(change,maxchange);
		upperswap(objective,maxobjective);

  } //if (!solved)...


	if ((pBlock!=NULL) && (pBlock->IsOn())){pBlock->Update(myBlockID,NOT_A_STRING,t);}

  change=maxchange;
  objective=maxobjective;
}
//------------------------------------------------------------------------
double CAreaVortex::GetMaxError          (const double &t) const{return 0.0;}
//************************************************************************
void CAreaVortex::SolveMQCoeff(const double t){
	//solves for Multiquadric coefficients based upon curl @ ctrl pts
  MQInterpolate(curlctrl,zcurlctrl,ncurlctrl,
							  MQcoeff,zMQbasis,MQorder,0.0,
							  ave_curl);
}
//************************************************************************
cmplex CAreaVortex::GetInteriorPotential (const cmplex &z,const double &t) const{
	double r;	
	cmplex omega(0.0);

	r=abs(z-pLinevortexBoundary->Centroid());
	omega=-IM*0.25*ave_curl*r*r;
	for(int n=0; n<MQorder;n++){
		r=abs(z-zMQbasis[n]);
		omega+=-IM*MQcoeff[n]*r*r*r/9.0;
	}
	return omega;
}
//************************************************************************
cmplex CAreaVortex::GetInteriorW (const cmplex &z,const double &t) const{
	double r;	
	cmplex W(0.0);

  W=IM*0.5*ave_curl*conj(z-pLinevortexBoundary->Centroid());
  for(int n=0; n<MQorder;n++){
    r=abs(z-zMQbasis[n]);
    W-=-IM*MQcoeff[n]*r*conj(z-zMQbasis[n])/3.0;
	}
	return W;
}
//************************************************************************
bool CAreaVortex::CalcCurl(const double t){return true;} //for given curl (default), nothing changes 
//************************************************************************
cmplex CAreaVortex::GetDischargePotential(const cmplex &z,const double &t) const{
		
	cmplex omega(0.0);	
	static bool integ_disabled(false);
		
	if (pLinevortexBoundary->IsInside(z)){

		omega=GetInteriorPotential(z,t);
		
		if (!integ_disabled){ //used to avoid infinite loop
			integ_disabled=true;

			//calculate actual_head = head_on_boundary-integrate(Qx/K dx)
      double head;
			/*cmplex zint,dz,startz,unit;
      double dist,Ql,kloc;
			int iter(0);

      startz=cmplex(10.0,60.0); //TMP DEBUG
			dz=(z-startz)/100.0;
			dist=abs(z-startz);
			unit=(z-startz)/abs(z-startz);

			head=pLayer->GetHead(startz,t);
			for (zint=startz;abs(z-zint)>REALSMALL;zint+=dz){
        
				kloc=pLayer->GetCond(zint);

				Ql=(pLayer->GetW(zint,t)*unit).real();  //component of Q in integration direction

				head-=Ql/kloc/min(head,pLayer->GetThick(zint))*abs(dz);

				if (abs(dz)>abs(z-zint)){dz=z-zint;}
				iter++;
			}
			//cout <<"iter:"<<iter<<endl;*/
			head=27;

			omega+=ConvertToPotential(head-pLayer->GetBase(z),pLayer->GetCond(z),pLayer->GetThick(z),pLayer->GetSeaLevel(),pLayer->GetSaltwaterSG());
			omega-=pLayer->GetDischargePotential(z,t).real();//hurt by static representation in Layer

			//omega.real=ConvertToPotential(act_head,K,H,Hsea,alpha)-Phi not inside
			integ_disabled=false;
		}
		//if (ProgramAborted()){ExitGracefully("",BAD_DATA);}
	}

	omega+=pLinevortexBoundary->GetDischargePotential(z,t);
	omega+=pDipoleBoundary    ->GetDischargePotential(z,t);

	//omega+=pLinevortexBoundary->CleanStreamFunct(z,t);

  return omega; 
  
}
//************************************************************************
cmplex CAreaVortex::GetW(const cmplex &z,const double &t) const{

	cmplex W(0.0);
  
	if (pLinevortexBoundary->IsInside(z)){
		W=GetInteriorW(z,t);
	}

	W+=pLinevortexBoundary->GetW(z,t);
	W+=pDipoleBoundary    ->GetW(z,t);
	
	return W; 
}
//************************************************************************
double CAreaVortex::GetCurl(const cmplex &z, const double &t) const{ 
  double curl(0.0);
  if (pLinevortexBoundary->IsInside(z)){
    if (MQorder<=1){
			curl=ave_curl;
		}
    else {
      for(int n=0; n<MQorder;n++){
				curl+=MQcoeff[n]*abs(z-zMQbasis[n]);
			}
      curl+=ave_curl;
		}
		return curl;       
	}
  return 0.0;
} 
//************************************************************************
void   CAreaVortex::WriteItself(ofstream &SOL, const double &t) const{
	pLinevortexBoundary->WriteItself(SOL,t);
	pDipoleBoundary    ->WriteItself(SOL,t);
}
//----------------------------------------------------------------------
bool   CAreaVortex::ReadItself(ifstream &SOL){
	bool warm(true);
	warm=pLinevortexBoundary->ReadItself(SOL);
	warm=pDipoleBoundary    ->ReadItself(SOL);
	return warm;
}
//------------------------------------------------------------------------
void CAreaVortex::WriteOutput(const double &t) const{
  //write errors file
 /* double error(0.0);

  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	
	for (int m=0; m<ncurlcontrol; m++){
	   error=???
		  ERRORS <<zcurlctrl[m].real()<<","
			       <<zcurlctrl[m].imag()<<","
						 <<FLUX_ERROR_TAG     <<","<<element_type<<","<<","<<error<<" "<<rel_error<< endl;
      //TMP DEBUG
		}

	ERRORS.close();*/
}
//************************************************************************
cmplex CAreaVortex::Centroid						 () const                               {return pLinevortexBoundary->Centroid();}
bool   CAreaVortex::IsInSquare					 (const cmplex &zc,const double w) const{return pLinevortexBoundary->IsInSquare(zc,w);}
bool   CAreaVortex::IsInCircle					 (const cmplex &zc,const double r) const{return pLinevortexBoundary->IsInCircle(zc,r);}
bool   CAreaVortex::PartInCircle				 (const cmplex &zc,const double r) const{return pLinevortexBoundary->PartInCircle(zc,r);}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string AreaSink, string name 
	{double x double y}x(numlines+1) 
	& 
	{double x double y double curl}x(numcontrolpts)
	&[int precision]
----------------------------------------------------------------------*/
/*CAreaVortex *CAreaVortex::Parse(ifstream &input, int &l, 
														CSingleLayerABC *pLay, char * Name){

	CAreaVortex *pAreaSink=NULL;
  bool       eof(false),done(false);
	int        Len,thisprec,nlines(0),npoints(0),i;
  cmplex	   stringp[MAXLINES];       
  cmplex	   leakzp [MAXLINES];
  double     leakp  [MAXLINES];
	char      *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Area Sink element"<<endl;}  eof=TokenizeLine(input,s,Len); l++; 
  do {
		if (nlines>=MAXLINES) { ExitGracefully("CAreaVortex::Parse- too many line segments in area sink",TOO_MANY);}
		if       (Len==0){                             eof=TokenizeLine(input,s,Len); l++;}
	  else if  (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]); nlines++; eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len<=3) && (!strcmp(s[0],"&"))) {
      stringp[nlines]=stringp[0]; nlines++; done=true; 
		}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while (!done);

	done=false;                                      eof=TokenizeLine(input,s,Len); l++;
	do {
    if (npoints>=MAXLINES) { ExitGracefully("CAreaVortex::Parse- too many control points in area sink",TOO_MANY);}
		if ((Len==3) && (strcmp(s[0],"&"))){
      leakzp[npoints]=s_to_c(s[0],s[1]);
			leakp [npoints]=s_to_d(s[2]);
			npoints++;                                   eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);}
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pAreaSink = new CAreaVortex(Name,
				  											pLay,
				  											stringp,
				  											leakp,
				   											leakzp,
				  											thisprec,
				  											npoints,
																nlines-1);
			done=true;
		}
    else if (Len==0) {                             eof=TokenizeLine(input,s,Len); l++;}
		else {cout <<"line"<< l << "is wrong length: "<<Len<<endl; break;}
	} while (!done);

  if (eof) {return NULL;}
	else     {return pAreaSink;}
}*/
//************************************************************************
void CAreaVortex::SetPrecision(const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("SetPrecision::Improper precision level specified",BAD_DATA);}
	switch(Precision){
		case(0): {order=1; fold=1.0; break;}
		case(1): {order=5; fold=1.2; break;}
		case(2): {order=5; fold=1.3; break;}
		case(3): {order=12;fold=1.5; break;}
	  case(4): {order=10;fold=1.4; break;}
	  case(5): {order=12;fold=1.5; break;}
	  case(9): {order=15;fold=1.5; break;}
	}
}
//************************************************************************
/*void CAreaVortex::MQuadric(cmplex z, cmplex &om,cmplex &W, double &leak) const{

  double r=abs(z-Centroid());

  leak=ave_curl;
  om=  cmplex(.25*ave_curl*r*r,0.0);
  W = -0.5*ave_curl*conj(z-Centroid());

  for(int n=0; n<MQorder;n++){
    r=abs(z-zMQbasis[n]);

    leak+= MQcoeff[n]*r;
    om+=   MQcoeff[n]*r*r*r/9;
    W+=   -MQcoeff[n]*r*conj(z-zMQbasis[n])/3.0;
  }
}*/
		//TMP DEBUG
		/*if (CAnalyticElem::Cosmetic){
			//numerically integrate stream function from centroid
			double strm(0.0);
			cmplex W;
			cmplex cent=Centroid();
			cmplex zt;
			double x,y;
			for (y=cent.imag(); fabs(y-cent.imag())<fabs(z.imag()-cent.imag()); y+=(z.imag-cent.imag())/500.0){
				zt=cmplex(cent.real(),y);
//zt=z;
				W=-0.5*Relax*ave_curl*conj(zt-cent);
				for(n=0; n<MQorder;n++){
					r=abs(zt-zMQbasis[n]);
					W-=MQcoeff[n]*r*conj(zt-zMQbasis[n])/3.0;
				}

				strm+=(W.real())*(z.imag-cent.imag())/500.0;
			}
			
			for (x=cent.real(); fabs(x-cent.real())<fabs(z.real()-cent.real()); x+=(z.real-cent.real())/500.0){
				zt=cmplex(x,z.imag());

				W=-0.5*Relax*ave_curl*conj(zt-cent);
				for(n=0; n<MQorder;n++){
					r=abs(zt-zMQbasis[n]);
					W-=MQcoeff[n]*r*conj(zt-zMQbasis[n])/3.0;
				}

				strm+=(-W.imag())*(z.real-cent.real())/500.0;
			}

			omega+=IM*strm;
		}
		return omega;//TMP DEBUG */
