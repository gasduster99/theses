#include "StringElem.h"

/************************************************************************
                           CStringElem
************************************************************************	
	  		                   CONSTRUCTOR
************************************************************************/
CStringElem::CStringElem(){NLines=0;}
//------------------------------------------------------------------------
CStringElem::CStringElem(  char									 *Name, 
													 const CSingleLayerABC *pLay, 
													 const cmplex					 *points, 
													 const int							NumOfLines,
													 const bool							closedstring, 
													 const linetype					LType,
													 const int							prec) 
					  :CAnalyticElem(Name,pLay){

	int i,n;
	bool clock;
	double fold;

  //initialize variable values   
	NLines      =NumOfLines;    
	ltype       =LType;
	closed      =closedstring;

	SetPrecision(prec,order,fold);
	nlinecontrol=int(order*fold);

  c[0]=NULL;	   c[1]=NULL;      c[2]=NULL;       c[3]=NULL; 
	
  if (order>MAX_ORDER){ExitGracefully("StringElem::Order too large",TOO_MANY);}
  if (nlinecontrol>MAXLINECONTROL){ExitGracefully("StringElem::Too many control points on element",TOO_MANY);}

	//Create Dynamic Arrays
	X=            new double        [nlinecontrol];
  ze =          new cmplex        [NLines+1];
	FForder=      new int           [NLines];
  myBlockIDs=   new int           [NLines];
	disabled=     new bool          [NLines];
  pBlocks=		  new COwnerABC    *[NLines];

  zctrl=        new cmplex       *[NLines];
  JumpCoeff=    new double       *[NLines];
  FarFldCoeff=  new double       *[NLines];

	bn=           new double       *[NLines];
	dn=           new double       *[NLines];
	en=           new double       *[NLines];
	gn=           new double       *[NLines];
	hn=           new double       *[NLines];
	if (hn==NULL){ExitGracefully("CStringElem:Constructor",OUT_OF_MEMORY);}

  for(i=0;i<NLines;i++){
    zctrl[i]=       new cmplex   [nlinecontrol]; 
		JumpCoeff[i]=   new double   [order+1];
    FarFldCoeff[i]= new double   [MAX_FAR_FIELD];
		bn[i]         = new double   [order+1];
		dn[i]         = new double   [order+1];
		en[i]         = new double   [order+1];
		gn[i]         = new double   [order+1];
		hn[i]         = new double   [order+1];
		if (hn[i]==NULL){ExitGracefully("CStringElem:Constructor",OUT_OF_MEMORY);}
		pBlocks[i]=		  NULL;
  }


  //copy coordinate array
  for(i=0;i<=NLines;i++){ze[i]=points[i];}
  if (closed){
		Area(clock);
    if (clock){for(i=0;i<=NLines;i++){ze[i]=points[NLines-i];}}
	}

  //create matrices of control point location 
  Interpolate(X,nlinecontrol,-1,1);
  for (i=0; i<NLines; i++){
    Interpolate(zctrl[i],nlinecontrol,ze[i],ze[i+1]);
  }

  //initialize all coeff to 0
  for (i=0; i<NLines; i++){
		FForder   [i]=MAX_FAR_FIELD;
		disabled  [i]=false;
    for (n=0; n<=order;     n++)  {   
			JumpCoeff[i][n]=0.0;
			bn[i][n]=0.0;
			dn[i][n]=0.0;
			en[i][n]=0.0;
			gn[i][n]=0.0;
			hn[i][n]=0.0;		
		}
    for (n=0; n<FForder[i]; n++)  { FarFldCoeff[i][n]=0.0;}
	}
}

//-----------------------------------------------------------------------
CStringElem::~CStringElem(){
  if (globaldebug){cout <<"   DESTROYING STRINGELEM "<<name<<endl;}
  delete [] pBlocks; //just deletes pointers to blocks
	delete [] myBlockIDs;
	delete [] X;
	delete [] ze;
	delete [] disabled;
  for(int i=0;i<NLines;i++){
    delete [] zctrl[i];
    delete [] JumpCoeff[i];
    delete [] FarFldCoeff[i];
		delete [] bn[i];
		delete [] dn[i];
		delete [] en[i];
		delete [] gn[i];
		delete [] hn[i];
  }
  delete [] zctrl;
  delete [] JumpCoeff;
  delete [] FarFldCoeff;
	delete [] bn;
	delete [] dn;
	delete [] en;
	delete [] gn;
	delete [] hn;
} 
/***********************************************************************
				ACCESSOR FUNCTIONS
***********************************************************************/
cmplex  CStringElem::GetZ(const int pointno) const {return ze[pointno];}
//-----------------------------------------------------------------------
cmplex *CStringElem::GetZArray()             const {return ze;}
//-----------------------------------------------------------------------
int     CStringElem::GetNumSegs()            const {return NLines;}
//-----------------------------------------------------------------------
bool    CStringElem::IsDisabled(const int i) const {return disabled[i];}
//-----------------------------------------------------------------------
bool    CStringElem::HasFlux()               const {return ltype==LINESINK;}           

/***********************************************************************
				MANIPULATOR FUNCTIONS
***********************************************************************/
void CStringElem::SetBlockOwner(COwnerABC *BlockPtr, int seg,int IDinBlock){
  pBlocks   [seg]=BlockPtr;
  myBlockIDs[seg]=IDinBlock;
}
//-----------------------------------------------------------------------
void CStringElem::SetAsDisabled(int i){disabled[i]=true;}
/************************************************************************
        Global to Local coordinates
-----------------------------------------------------------------------*/
cmplex CStringElem::GlobalToLocal(const cmplex &z, const int i) const{
   return (z-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
}

/************************************************************************
        Get String Discharge Potential
-----------------------------------------------------------------------*/
cmplex CStringElem::GetDischargePotential(const cmplex &z,const double &t) const{
  static cmplex omega;
	omega=0.0;
  for (int i=0; i<NLines; i++){if(!disabled[i]){omega+=GetSegmentPotential(i,z,t);}}

	if ((CAnalyticElem::Cosmetic)  && (ltype==LINESINK  )){omega+=IM*CleanStreamFunct(z,t);}
	else if                           (ltype==LINEVORTEX) {omega+=   CleanStreamFunct(z,t);}

  return omega;
}
/************************************************************************
        Get String Discharge Function
-----------------------------------------------------------------------*/
cmplex CStringElem::GetW(const cmplex &z,const double &t) const{
  static cmplex W;
	W=0.0;
  for (int i=0; i<NLines; i++){if(!disabled[i]){W+=GetSegmentW(i,z,t);}}
  return W;
}
/************************************************************************
        Get String Discharge Derivative Function
-----------------------------------------------------------------------*/
cmplex CStringElem::GetGx(const cmplex &z,const double &t) const{
  cmplex Gx(0.0);
  for (int i=0; i<NLines; i++){if(!disabled[i]){Gx+=GetSegmentGx(i,z,t);}}
  return Gx;
}
/************************************************************************
        Get String Discharge Derivative Function
-----------------------------------------------------------------------*/
cmplex CStringElem::GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const double &t) const{
	cmplex Qflux(0.0);
  for (int i=0; i<NLines; i++){if(!disabled[i]){Qflux+=GetSegFluxThruFace(i,z1,z2,t);}}
  return Qflux;
}

/***********************************************************************
				Get String Leakage Function
----------------------------------------------------------------------*/
double CStringElem::GetLeakage(const cmplex &z, const double &t,const leak_type ltype) const{
  double Leak(0.0);
  for (int i=0; i<NLines; i++){if(!disabled[i]){Leak+=GetSegmentLeakage(i,z,t,ltype);}}
  return 0.0;
}

/***********************************************************************
				Get String Net Discharge Function
----------------------------------------------------------------------*/
double CStringElem::GetNetDischarge(const double &t) const{
	double Discharge(0.0);
	if (ltype==LINESINK){
		for (int i=0; i<NLines; i++){if(!disabled[i]){Discharge+=GetSegmentDischarge(i,t);}}
	}
	return Discharge;
}
//----------------------------------------------------------------------
double CStringElem::GetNetCurl     (const double &t) const{
	double curl(0.0);
	if (ltype==LINEVORTEX){
		for (int i=0; i<NLines; i++){if(!disabled[i]){curl+=GetSegmentCurl(i,t);}}
	}
	return curl;
}
/***********************************************************************
				Get String Integrated Budget Function
----------------------------------------------------------------------*/
void CStringElem::GetIntegratedBudget(const cmplex &z1, const cmplex &z2, const cmplex &z3, 
																			const double &t, double &inQ, double &outQ) const{
  double tmpin,tmpout;
	inQ=outQ=0.0;
	if (ltype==LINESINK){
		for (int i=0; i<NLines; i++){
			if(!disabled[i]){
			  GetSegmentBudget(i,z1,z2,z3,t,tmpin,tmpout);
				inQ+=tmpin;
				outQ+=tmpout;
			}
		}
	}
}
/***********************************************************************
				Get String Flux Distribution Function
----------------------------------------------------------------------*/
void CStringElem::GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const{
	Q1=Q2=0.0;
	if (ltype==LINESINK){
    double tmpQ1,tmpQ2;
		for (int i=0; i<NLines; i++){
			if(!disabled[i]){
				GetSegDistribution(i,z1,z2,t,tmpQ1,tmpQ2);
				Q1+=tmpQ1;
			  Q2+=tmpQ2;
			}
		}
	}
}
/************************************************************************
      Get String Singular Strength
************************************************************************/
cmplex CStringElem::GetSingularStrength  (const cmplex &z, const double &t) const{
	if (ltype==LINESINK)  {
		if (abs(z-ze[0])<NEAR_FEATURE){
			return 0.0;//TMP DEBUG
		}
		else if (abs(z-ze[NLines])<NEAR_FEATURE){
			return 0.0;//TMP DEBUG
		}
	}
	return 0.0;
}
/***********************************************************************
				Get Segment Discharge Potential
----------------------------------------------------------------------*/
cmplex CStringElem::GetSegmentPotential(const int i,const cmplex &z,const double &t) const{
  cmplex  omega(0.0),mult(1.0,0);
	cmplex  Z=GlobalToLocal(z,i);
	double  modulus(abs(Z));

	if (disabled[i]){return 0.0;}

	if ((ltype==DOUBLET) || (ltype==LINEVORTEX)){mult=-IM;}

	//far-field 
  if (modulus>ISFAR){ 
    omega =mult*LaurentRe(Z,(int)(FForder[i]*pow(modulus/ISFAR,FF_POWER)),FarFldCoeff[i]);
		omega+=mult*GetSegmentLog(i,z,t);
  }  
  //near-field
  else if (((abs(Z-E1)>MOVE_DIST) && (abs(Z-E2)>MOVE_DIST))){

		cmplex logterm=log(Zm1oZp1(Z));

		for (int n=0; n<=order; n++){cheby[n]=JumpCoeff[i][n]*logterm+bn[i][n];}

    omega =mult*IOVER2PI*cClenshaw(order,Z,cheby); 
		omega+=mult*GetSegmentLog(i,z,t);
  }
  //singularities
  else {
		if (abs(Z-E1) == abs(Z-E2)){
			return 0.0;
			ExitGracefully("CStringElem::GetPotential: Bad (infinite) Input location",RUNTIME_ERR);}
		else {
			omega=GetSegmentPotential(i,z+MOVE_DIST*abs(ze[i+1]-ze[i]),t);//trouble for inhomogeneities-direction matters
		} 
	}


  return omega;
}
/***********************************************************************
				Get Segment Discharge Function
----------------------------------------------------------------------*/
cmplex CStringElem::GetSegmentW(const int i,const cmplex &z,const double &t) const{

  cmplex  W(0.0),mult(1.0/(PI*(ze[i+1]-ze[i])));  
	cmplex  Z=GlobalToLocal(z,i);
	double  end1(0.0),end2(0.0),modulus(abs(Z)),sum;
  int     n;


	if (disabled[i]){return 0.0;}

  if ((ltype==LINESINK) || (ltype==LINEVORTEX)){
		sum=-1.0;
    for (n=0; n<=order; n++){
			sum*=-1.0;
      end1+=JumpCoeff[i][n]; 
	    end2+=JumpCoeff[i][n]*sum; 
		}
  }
	if ((ltype==DOUBLET) || (ltype==LINEVORTEX)){mult*=-IM;}

  //far-field
  if (modulus>ISFAR){ 
		W+=mult*2.0*PI*LaurentDerRe(Z,(int)(FForder[i]*pow(modulus/ISFAR,FF_POWER)),FarFldCoeff[i]);
    if ((ltype==LINESINK) || (ltype==LINEVORTEX)){
			W-=mult*(end2/(Z+1.0)-end1/(Z-1.0));
		}
  }  
  //near-field
  else if ((abs(Z-E1)>MOVE_DIST) && (abs(Z-E2)>MOVE_DIST)){

		cmplex logterm (log(Zm1oZp1(Z)));
		cmplex logterm2(2.0/(Z+1.0)/(Z-1.0));

		//faster than gZj-> n! simplified to n
		for (n=0; n<=order; n++){
			cheby[n]= JumpCoeff[i][n]*logterm2+dn[i][n]*logterm+en[i][n];
		}

		W-=mult*cClenshaw(order,Z,cheby);  
		if ((ltype==LINESINK  ) || (ltype==LINEVORTEX)){
			W-=mult*(end2/(Z+1.0)-end1/(Z-1.0));
		}		
  }
	
  //singularities
  else {
		if (abs(Z-E1) == abs(Z-E2)){
			return 0.0;
			ExitGracefully("CStringElem::GetW: Bad (infinite) Input location",RUNTIME_ERR);}
		else {W=GetSegmentW(i,z+cmplex(MOVE_DIST*abs(ze[i+1]-ze[i]),0),t);}
	}

  return W;
}
/***********************************************************************
				Get Segment Discharge Derivative Function
----------------------------------------------------------------------*/
cmplex CStringElem::GetSegmentGx(const int i,const cmplex &z,const double &t) const{
  cmplex  Gx(0.0),jump(0.0),mult(1.0/(PI*(ze[i+1]-ze[i])*(ze[i+1]-ze[i])));   
	cmplex  Z=GlobalToLocal(z,i);
	double  end1(0.0),end2(0.0),modulus(abs(Z)),sum;
  int     n;

	if (disabled[i]){return 0.0;}

  if ((ltype==LINESINK) || (ltype==LINEVORTEX)){
		sum=-1.0;
    for (n=0; n<=order; n++){
			sum*=-1.0;
      end1+=JumpCoeff[i][n]; 
	    end2+=JumpCoeff[i][n]*sum; 
		}
  }

	if ((ltype==DOUBLET) || (ltype==LINEVORTEX)){mult=-IM/(PI*(ze[i+1]-ze[i])*(ze[i+1]-ze[i]));}

  if (modulus>ISFAR){ 
    Gx=-mult*4.0*PI*LaurentDxxRe(Z,(int)(FForder[i]*pow(modulus/ISFAR,FF_POWER)),FarFldCoeff[i]);
		if ((ltype==LINESINK) || (ltype==LINEVORTEX)){
			Gx+=2.0*mult*(end2/(Z+1.0)/(Z+1.0)-end1/(Z-1.0)/(Z-1.0));
		}
  }  
	//near-field
  else if ((abs(Z-E1)>MOVE_DIST) && (abs(Z-E2)>MOVE_DIST)){

		cmplex logterm (log(Zm1oZp1(Z)));
		cmplex logterm2(4.0/(Z+1.0)/(Z-1.0));
		cmplex logterm3(-logterm2*Z/(Z+1.0)/(Z-1.0));

		for (n=0; n<=order; n++){
			               cheby[n]=JumpCoeff[i][n]*logterm3;   
			if (n<order)  {cheby[n]+=dn[i][n]*logterm2;}
			if (n<order-1){cheby[n]+=gn[i][n]*logterm;}
			if (n<order-2){cheby[n]+=hn[i][n];}
		}
		jump=cClenshaw(order,Z,cheby);

		Gx =-2.0*mult*jump; 

		if ((ltype==LINESINK  ) || (ltype==LINEVORTEX)){
			Gx+=2.0*mult*(end2/(Z+1.0)/(Z+1.0)-end1/(Z-1.0)/(Z-1.0));
		}		
  }
  //singularities
  else {
		if (abs(Z-E1) == abs(Z-E2)){
			return 0.0;
			ExitGracefully("CStringElem::GetGx: Bad (infinite) Input location",RUNTIME_ERR);}
		else {Gx=GetSegmentGx(i,z+cmplex(MOVE_DIST*abs(ze[i+1]-ze[i]),0),t);}
	}

	return Gx; 
}
/***********************************************************************
				Get Segment Flux Through Face Function
----------------------------------------------------------------------*/
cmplex CStringElem::GetSegFluxThruFace   (const int i, const cmplex &z1,const cmplex &z2, const double &t) const{
  static cmplex zint;
	if (Intersect(z1,z2,ze[i],ze[i+1],zint)!=INTERSECTED){
		//one-point numerical integration

		//cout << -abs(z2-z1)*(conj(IM*(z2-z1)/abs(z2-z1))*conj(GetSegmentW(i,0.5*(z2+z1),t))).real()/(IM*conj(GetSegmentPotential(i,z1,t)-GetSegmentPotential(i,z2,t))).real()<<endl;

		//return -abs(z2-z1)*conj(IM*(z2-z1)/abs(z2-z1))*conj(GetSegmentW(i,0.5*(z2+z1),t));

		return IM*conj(GetSegmentPotential(i,z1,t)-GetSegmentPotential(i,z2,t)); //complex flux thru face

	}
	else{
		return IM*conj(GetSegmentPotential(i,z1,t)-GetSegmentPotential(i,zint-(z2-z1)*REALSMALL,t)-
                   GetSegmentPotential(i,z2,t)+GetSegmentPotential(i,zint+(z2-z1)*REALSMALL,t));
	}
}
/***********************************************************************
				Get Segment Leakage Function
----------------------------------------------------------------------*/
double CStringElem::GetSegmentLeakage (const int i, const cmplex &z, const double &t,const leak_type ltype) const{
	return 0.0;
}
/**********************************************************************
				Get Segment Net Discharge Function
----------------------------------------------------------------------*/
double CStringElem::GetSegmentDischarge(const int i,const double &t) const{
  //SHOULD BE SIMPLIFIED- test with superblocks
  static double Discharge;
	static double sum;
	static double end1;
	static double end2;

	Discharge=sum=end1=end2=0.0;

	if (disabled[i]){return 0.0;}

  if (ltype==LINESINK){
		sum=-1.0;
    for (int n=0; n<=order; n++){
			sum*=-1.0;
			end1+=JumpCoeff[i][n];  
			end2+=JumpCoeff[i][n]*sum;//sum=(-1)^n
		}
		Discharge=end2-end1;
		/*for (int n=0; n<=order; n++){
			sum*=-1.0; 
			Discharge+=JumpCoeff[i][n]*(sum-1.0);//sum=(-1)^n
		}*/
	}
	return Discharge;
}
/**********************************************************************
				Get Segment Net Curl Function
----------------------------------------------------------------------*/
double CStringElem::GetSegmentCurl(const int i,const double &t) const{
  double curl(0.0);
	double sum(0),end1(0),end2(0);
  if (ltype==LINEVORTEX){
		sum=-1.0;
    for (int n=0; n<=order; n++){
			sum*=-1.0;
			end1+=JumpCoeff[i][n];  
			end2+=JumpCoeff[i][n]*sum;//sum=(-1)^n
		}
		curl=end2-end1;
	}
	return curl;
}
/**********************************************************************
				Get Segment Integrated Budget Function
	Gets influx and outflux [L^3/T] added to the specified triangular area from segment i
----------------------------------------------------------------------*/
void CStringElem::GetSegmentBudget(const int i, const cmplex &z1, const cmplex &z2, const cmplex &z3, 
																	 const double &t, double &inQ, double &outQ) const{
	inQ=outQ=0.0;
	if (ltype==LINESINK){
		cmplex zint1,zint2,zint;
		double inQ1,outQ1,inQ2,outQ2;
		int  nint=0;		
		bool in1 =InTriangle(ze[i],z1,z2,z3);
		bool in2 =InTriangle(ze[i+1],z1,z2,z3);		

    //find intersections
		if (Intersect(z1,z2,ze[i+1],ze[i],zint)==INTERSECTED){zint1=zint;nint++;}
		if (Intersect(z2,z3,ze[i+1],ze[i],zint)==INTERSECTED){if (nint==0){zint1=zint;}else{zint2=zint;}nint++;}
		if (Intersect(z3,z1,ze[i+1],ze[i],zint)==INTERSECTED){if (nint==0){zint1=zint;}else{zint2=zint;}nint++;}
		
		if ((!in1) && (!in2) && (nint==0)){return;} //no element/triangle overlap

		double X1(-1.0),X2(1.0); //default- element totally contained
		
		if (nint>0){             //element partially contained
			if      ((nint==1) && (in2)){X1=2.0*abs(zint1-ze[i])/abs(ze[i+1]-ze[i])-1.0;}
			else if ((nint==1) && (in1)){X2=2.0*abs(zint1-ze[i])/abs(ze[i+1]-ze[i])-1.0;}
			else if  (nint==2)          {
				X1=2.0*abs(zint1-ze[i])/abs(ze[i+1]-ze[i])-1.0;
				X2=2.0*abs(zint2-ze[i])/abs(ze[i+1]-ze[i])-1.0;  
				if (X2<X1){double Xtmp=X1;X1=X2;X2=Xtmp;}
			}
		}
		if ((X2<-1.0) || (X2>1.0) || (X1<-1.0) || (X2>1.0)){
			ExitGracefully("CStringElem::GetSegmentBudget: bad intercepts",RUNTIME_ERR);}
		//integrate cumulative extraction from X1 to X2
		GetSegCumExtraction(i,X1,X2,inQ1,outQ1,inQ2,outQ2);
		inQ = inQ1 +  inQ2;
		outQ=outQ1 + outQ2;
	}
}
/**********************************************************************
				Get Segment Flux Distribution 
	Gets integrated extraction between points z1 and z2 along line
	Linearly Distributes these fluxes to endpoints Q1 and Q2
	Used for finite element source distribution
----------------------------------------------------------------------*/
void CStringElem::GetSegDistribution(const int i, const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const{
	Q1=Q2=0.0;
	if (ltype==LINESINK){
		cmplex Z1=GlobalToLocal(z1,i);
		cmplex Z2=GlobalToLocal(z2,i);

		if (fabs(abs(Z1-1.0)+abs(Z1+1.0)-2.0)>0.02){return;} //not near line
    if (fabs(abs(Z2-1.0)+abs(Z2+1.0)-2.0)>0.02){return;} //not near line

		if ((fabs(Z1.imag())<0.02) && (fabs(Z2.imag())<0.02)){

			//cout <<"Z1:"<<Z1.real()<<" Z2:"<<Z2.real()<<endl;
			double inj1(0),ext1(0),inj2(0),ext2(0);
			//double end1,end2;
			GetSegCumExtraction(i,max(min(Z1.real(),1.0),-1.0),max(min(Z2.real(),1.0),-1.0),inj1,ext1,inj2,ext2);
			//end1=GetWJump(i,max(min(Z1.real(),1.0),-1.0),t).imag();
			//end2=GetWJump(i,max(min(Z2.real(),1.0),-1.0),t).imag();
			//cout <<"---"<<endl;
			//Q1=fabs(end1)/(fabs(end1)+fabs(end2))*(inj1+inj2-ext1-ext2); 
			//Q2=fabs(end2)/(fabs(end1)+fabs(end2))*(inj1+inj2-ext1-ext2);
			//cout <<Q1<<" "<<Q2<<" "<<Q1+Q2<<endl;
			Q1=inj1-ext1;
			Q2=inj2-ext2;
			//cout <<Q1<<" "<<Q2<<" "<<Q1+Q2<<endl;

		}
	}
}

/**********************************************************************
				Get Segment Cumulative Extraction
	Gets integrated extraction and injection (both positive values) between points X1 and X2 along line
----------------------------------------------------------------------*/
void CStringElem::GetSegCumExtraction(const int i, const double &X1, const double &X2, double &inj1, double &ext1, double &inj2, double &ext2) const{
  

	if ((X1<-1.0) || (X1>1.0)){ExitGracefully("CStringElem::GetSegCumExtraction: bad extents specified(1)",RUNTIME_ERR);}
	if ((X2<-1.0) || (X2>1.0)){ExitGracefully("CStringElem::GetSegCumExtraction: bad extents specified(2)",RUNTIME_ERR);}
	
	if (ltype==LINESINK){
		//Use numerical integration
		int    numpts=50;
		double Xl,Xr,tmp,tmppos,tmpneg,N1,N2,dX,Xtmp;
		double t=0; //TMP DEBUG

		if (fabs(X1)==1.0){Xl=X1*0.9999;}else{Xl=X1;}
		if (fabs(X2)==1.0){Xr=X2*0.9999;}else{Xr=X2;}
		if (Xr<Xl){tmp=Xr; Xr=Xl; Xl=tmp;}
		dX=(Xr-Xl)/numpts;
		inj1=ext1=inj2=ext2=0.0;
		for (Xtmp=Xl; Xtmp<(Xr-0.5*dX); Xtmp+=dX){
			N1=(Xtmp+dX/2-Xl)/(Xr-Xl); //basis functions (assumed linear)
			N2=1.0-N1;
			IntegrateLine(Xtmp,Xtmp+dX,-GetWJump(i,Xtmp,t).imag(),-GetWJump(i,Xtmp+dX,t).imag(),tmppos,tmpneg); 
			inj1+=tmppos*N1; inj2+=tmppos*N2;
			ext1-=tmpneg*N1; ext2-=tmpneg*N2;
		}
		inj1*=abs(ze[i+1]-ze[i])/2.0; 
		ext1*=abs(ze[i+1]-ze[i])/2.0;
		inj2*=abs(ze[i+1]-ze[i])/2.0;
		ext2*=abs(ze[i+1]-ze[i])/2.0;

	//prefer an analytic solution
	/*double sum(0.0); 
	double tau1=acos(X1);
	double tau2=acos(X2);
  for (int n=0; n<order-1; n++){
		if (sin(n*tau1)!=0.0){
			sum+=dn[i][n]*(log(fabs(sin(n*tau2)))-log(fabs(sin(n*tau1))) )/(double)(n);
		}
	}
	cout <<"GetSegCumExtraction:test:"<<sum*abs(ze[i+1]-ze[i])/2.0<<" "<<inj-ext<<endl;*/
	}
}
/***********************************************************************
 Get Segment Logarithm Potential Function
	Returns the potential Omega=Phi + i Psi from the endpoint singularities of a linesink or line vortex
	Corrects branch cut if desired
----------------------------------------------------------------------*/
cmplex  CStringElem::GetSegmentLog(const int i, const cmplex &z, const double &t) const{
	cmplex omega(0.0);

	if ((ltype==LINESINK) || (ltype==LINEVORTEX)){
		cmplex Z  =GlobalToLocal(z,i);
		cmplex ZBH=GlobalToLocal(pLayer->GetBlackHole(),i);
    cmplex logend1=log((Z-1.0)/abs(ZBH-1.0));
    cmplex logend2=log((Z+1.0)/abs(ZBH+1.0));  

		double sum(-1.0); 
    double end1(0.0),end2(0.0);
    for (int n=0; n<=order; n++){
			sum*=-1.0;
			end1+=JumpCoeff[i][n];  
			end2+=JumpCoeff[i][n]*sum;//sum=(-1)^n
		}

		omega+=IOVER2PI*(-(end1*logend1)+(end2*logend2));

		if ((CAnalyticElem::BranchcutLocus) && (ltype==LINESINK)){
			double argz=arg((z              -ze[i])/(ze[i+1]-ze[i]));;
			double arg1=arg(-(zBranchcutLocus-ze[i])/(ze[i+1]-ze[i]));;//if negative, away from BCL, else towards
			if      ((arg1>=0) && (argz>0.0) && (argz>arg1)) {omega-=IM*(end2-end1);}
			else if ((arg1< 0) && (argz<0.0) && (argz<arg1)) {omega+=IM*(end2-end1);}
		}
  }

	return omega;
}

/***********************************************************************
				Jump Functions
************************************************************************
 Returns the jump in potential DOmega=DPhi + i DPsi across the element
----------------------------------------------------------------------*/
cmplex CStringElem::GetOmJump(int i,const double &X, const double &t) const{
  cmplex jump(0.0);
	cmplex mult(IM);
	if ((ltype==LINESINK) || (ltype==DIPOLE)){mult=cmplex(1.0,0);}
  for(int n=0; n <=order; n++){jump+=mult*JumpCoeff[i][n]*cos(n*acos(X));} //TMP DEBUG- need to change
  return -IM*jump;
}
/* Returns the jump in complex discharge DW=DQX - iDQY across the element
  the imaginary point is the extraction per unit length
-----------------------------------------------------------------------*/
cmplex CStringElem::GetWJump(int i,const double &X, const double &t) const{
	//Can probably optimize as SUM(dn*cos(n*acos(X)))
  double sum, length(abs(ze[i+1]-ze[i])),Xtmp(X);
	cmplex jump(0.0);
  cmplex mult(IM);
	if ((ltype==LINESINK) || (ltype==DIPOLE)){mult=cmplex(1.0,0);}
	if (fabs(X)==1.0){Xtmp=0.9999*X;}
  for(int n=0; n <=order; n++){  
    sum=0.0;
    for(int k=n-1; k>=0; k-=2){
			if (k!=0){sum+=2.0*n*cos(k*acos(Xtmp));}
			else     {sum+=    n*cos(k*acos(Xtmp));}
		} 
		jump+=mult*JumpCoeff[i][n]*2.0/length*sum;
  }
  return -IM*jump;  
}
/*Returns the Cumulative Extraction of the element up to point X
  if direct==DIR_LEFT, integrates from X=-1 to X
  else,                integrates from X= 1 to X
-----------------------------------------------------------------------*/
double CStringElem::GetCumExtract(const int i, const double &X, const leftright direct, const double &t) const{
  double cumext(0);
	int    n;
	double mult(1.0);
	double sum(-1.0);
	if (direct==DIR_LEFT){mult=-1.0;}
  if (ltype==LINESINK){
		for (n=0; n<=order; n++){
			sum*=-1.0;
			cumext+=JumpCoeff[i][n]*(-sum-cos(n*acos(mult*X)));
		}
		return cumext;
	}
	return 0.0;		
}
/***********************************************************************
  Get Cumulative Flux
	Returns cumulative extraction ("base flow") along surface water network
	if this is an extraction element (e.g., linesink) this will return a value 
	if the point z is within NEAR_FEATURE of the element, otherwise, returns zero
--------------------------------------------------------------------------------*/
double CStringElem::GetCumulativeFlux    (const cmplex &z, const double &t) const{
  if (ltype==LINESINK){ //the only type of element with base flow
		double tmpflux=0.0;
		for (int i=0; i<NLines;i++){
			tmpflux=GetSegCumulativeFlux(i,z,t);
			if (tmpflux!=0){return tmpflux;}
		}
	}
	return 0.0;
}
/***********************************************************************
  Get Segment Cumulative Flux
	(VIRTUAL: OVERRIDEN FOR SURFACE WATER ELEMENTS)
--------------------------------------------------------------------------------*/
double CStringElem::GetSegCumulativeFlux    (const int i, const cmplex &z, const double &t) const{
	return 0.0;
}
/***********************************************************************
 CLEAN STREAM FUNCTION
  if Cosmetic=true, insures only single branch cut for a string of linesinks
	if BranchCutLocus=true, insures all branch cuts intersect zBranchCutLocus
-----------------------------------------------------------------------*/
double CStringElem::CleanStreamFunct(const cmplex &z,const double &t) const{

	if (ltype==LINESINK){
		
		double DischargeSum(0),streamadd(0),argz,arg1;

		if ((CAnalyticElem::Cosmetic) && (!CAnalyticElem::BranchcutLocus) ){
			for (int i=NLines-1; i>0; i--){
				DischargeSum-=GetSegmentDischarge(i,t);
				argz=arg((z      -ze[i])/(ze[i+1]-ze[i]));
				arg1=arg((ze[i-1]-ze[i])/(ze[i+1]-ze[i]));
				if      ((arg1>=0) && (argz>0.0) && (argz>arg1)) {streamadd+=DischargeSum;}
				else if ((arg1< 0) && (argz<0.0) && (argz<arg1)) {streamadd-=DischargeSum;}
			}

			DischargeSum-=GetSegmentDischarge(0,t);//last line- horizontally to the left
			argz=arg((z        -ze[0])/(ze[1]-ze[0]));
			arg1=arg((-1.0           )/(ze[1]-ze[0]));
			if      ((arg1>=0) && (argz>0.0) && (argz>arg1)) {streamadd+=DischargeSum;}
			else if ((arg1< 0) && (argz<0.0) && (argz<arg1)) {streamadd-=DischargeSum;}

		}
		return streamadd;
	}
	
	if (ltype==LINEVORTEX){

		double curlsum(0.0),potadd(0.0),argz,arg1;

		for (int i=NLines-1; i>0; i--){
			curlsum-=GetSegmentCurl(i,t);
		  argz=arg((z      -ze[i])/(ze[i+1]-ze[i]));
			arg1=arg((ze[i-1]-ze[i])/(ze[i+1]-ze[i]));
			if      ((arg1>0.0) && (argz>0.0) && (argz> arg1)) {potadd+=curlsum;}
			else if ((arg1<0.0) && (argz<0.0) && (argz< arg1)) {potadd-=curlsum;}
		}
		curlsum-=GetSegmentCurl(0,t);
		argz=arg((z           -ze[0])/(ze[1]-ze[0]));
	  arg1=arg((ze[NLines-1]-ze[0])/(ze[1]-ze[0]));
		if      ((arg1>0.0) && (argz>0.0) && (argz> arg1)) {potadd+=curlsum;}
		else if ((arg1<0.0) && (argz<0.0) && (argz< arg1)) {potadd-=curlsum;}
		return potadd;
	}
	return 0.0;
}
/***********************************************************************
				PURELY VIRTUAL MEMBER FUNCTIONS 
************************************************************************/
void CStringElem::GetSegMatrixBuildInfo (const int i,MatrixInfo &info){
	ExitGracefully("CStringElem::GetSegMatrixBuildInfo: Virtual function. Should never get called.",VIRTUAL_ERROR);}

/***********************************************************************
				MEMBER FUNCTIONS
***********************************************************************/
void CStringElem::UpdateBlock(const double &t) const {
	for (int i=0; i<NLines; i++){
    if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){
			pBlocks[i]->Update(myBlockIDs[i],i,t);
		}
	}
}
//-----------------------------------------------------------------------
void CStringElem::WriteItself(ofstream &SOL, const double &t) const                             {
	// Name
	// * [] [] [] [] [] [] [] [] ... ([] = coeff from a_0 to a_N)
  int i,n;
	SOL << name;
  for (i=0; i<NLines; i++){
		SOL << endl << " * ";
	  for (n=0; n<=order; n++){
			SOL<< JumpCoeff[i][n] <<" ";
		}
	}
  SOL<<endl;
}
//-----------------------------------------------------------------------
bool CStringElem::ReadItself(ifstream &SOL){
	// Name
	// * [] [] [] [] [] [] [] [] ... ([] = coeff from a_0 to a_N) x num segments
	// if old solution has different number of coefficients, then truncation is used

	int Len(0);
	char *s[MAXINPUTITEMS];
	bool warm(true);
	int i,n;

	if (!TokenizeLine(SOL,s,Len)) {} //skip name of element

  for (i=0; i<NLines; i++){
		if (!TokenizeLine(SOL,s,Len)) {

			if ((strcmp(s[0],"*"))){return false;}

			for (n=0; n<=min(order,Len-2); n++){
				JumpCoeff[i][n]=s_to_d(s[n+1]);          
			}
			for (n=min(order,Len-2)+1; n<=order; n++){
				JumpCoeff[i][n]=0.0;          
			}
		}
		else {warm=false; break;}

    SetFarFieldCoeff(i);
	}
	//cout << "warm= " <<warm<<endl;
  return warm;

}
//-----------------------------------------------------------------------
double CStringElem::GetMaxError(const double &t) const	{return 0.0;}
/************************************************************************
*************************************************************************
				EXPLICIT SOLVER FUNCTIONS 
*************************************************************************
************************************************************************/
int  CStringElem::GetSegDegreesOfFreedom(const int i) const {return order+1;}//if simple, return 1
//-----------------------------------------------------------------------
void CStringElem::GetSegUnitInfluences  (const int i, const int n, const cmplex *pts,
																				 const int NumCtrl, double *uPhi, cmplex *uQ, const double t){
  int m,n2;
  double oldJump[MAX_ORDER+1];

	//set pJumpCoeff[n] to 1, others to zero
	for (n2=0; n2<=order; n2++){
		oldJump[n2]=JumpCoeff[i][n2];
		JumpCoeff[i][n2]=0.0;
	}
	JumpCoeff[i][n]=1.0;
	SetFarFieldCoeff(i);

	//obtain unit influences	(a_n=1.0)
	for (m=0;m<NumCtrl;m++){
		uPhi[m]=GetSegmentPotential(i,pts[m],t).real();
		uQ  [m]=GetSegmentW        (i,pts[m],t);
	}

	//set back to previous state
	for (n2=0; n2<=order; n2++){
		JumpCoeff[i][n2]=oldJump[n2];
	}
	SetFarFieldCoeff(i);
}
//-----------------------------------------------------------------------
void CStringElem::SetSegCoeff           (const int i, double *coeff){
	//copies coefficients
	for (int n=0; n<=order; n++){JumpCoeff[i][n]=coeff[n];}
	SetFarFieldCoeff(i);
}
/************************************************************************
*************************************************************************
				Geometry Functions
*************************************************************************
************************************************************************/

/************************************************************************
				Interpolate Control Points
************************************************************************/
void CStringElem::Interpolate(double *Array, int size, double v1, double v2){
  for (int m=0; m<size; m++){
    Array[m]=((v2-v1)*0.5*CONTROLEND*cos(((double)(m+1)-0.5)/(double)(size)*PI)+(v1+v2)/2.0);
	}
}
//------------------------------------------------------------------------
void CStringElem::Interpolate(cmplex *Array, int size, cmplex v1, cmplex v2){
  for (int m=0; m<size; m++){
    Array[m]=((v2-v1)*0.5*CONTROLEND*(double)(cos(((double)(m+1)-.5)/(double)(size)*PI))+(v1+v2)/2.0);
	}
}
//*********************************************************************** 
void CStringElem::WriteGeometry(ofstream &BASEMAP) const              {
	int i;
	if (closed){
		BASEMAP << "\" StringElem \",  "<<(NLines+1)<<endl;	 
		for (i=0; i<=NLines; i++){
			BASEMAP <<ze[i].real()<< " , " <<ze[i].imag()<<endl;
		}
	}
	else{
		BASEMAP << "\" StringElem \",  "<<-(NLines+1)<<endl;	 
		for (i=0; i<=NLines; i++){
			BASEMAP <<ze[i].real()<< " , " <<ze[i].imag()<<endl;
		}
	}
}
//------------------------------------------------------------------------
cmplex CStringElem::Centroid() const{
  double x(0.0),y(0.0);
  for (int i=0; i<NLines; i++){
    x+=(ze[i+1].real()+ze[i].real())/2.0; 
    y+=(ze[i+1].imag()+ze[i].imag())/2.0; 
	}
  return cmplex (x/NLines,y/NLines); 
}
//------------------------------------------------------------------------
bool CStringElem::IsInSquare(const cmplex &zc,const double w) const{
  for (int i=0; i<=NLines; i++){
	if (  (ze[i].real() > zc.real() + (w/2.0)) ||
        (ze[i].real() < zc.real() - (w/2.0)) ||
        (ze[i].imag() > zc.imag() + (w/2.0)) ||
        (ze[i].imag() < zc.imag() - (w/2.0)))     {return false;}
  }
  return true;
}
//------------------------------------------------------------------------
bool CStringElem::IsInCircle(const cmplex &zc,const double r) const{
  for (int i=0; i<=NLines; i++){
		if (abs(ze[i]-zc)>r){return false;}
	}
	return true;
}
//------------------------------------------------------------------------
bool CStringElem::PartInCircle(const cmplex &zc,const double r) const{	for (int i=0; i<NLines;i++){
		if (SegPartInCir(zc,r,i)) {return true;}
	}
  return false;
}
//------------------------------------------------------------------------
bool CStringElem::SharesNode(const cmplex &zn) const{
	for (int i=0; i<=NLines;i++){
		if (zn==ze[i]) {return true;}
	}
	return false;
}
//************************************************************************
cmplex CStringElem::SegCentroid(const int seg) const{
  return cmplex ((ze[seg+1].real()+ze[seg].real())/2.0,
								 (ze[seg+1].imag()+ze[seg].imag())/2.0); 
}
//------------------------------------------------------------------------
bool CStringElem::SegIsInSquare(const cmplex &zc,const double w,const int seg) const{
  if ((ze[seg].real()   > zc.real() + (w/2.0)) ||
      (ze[seg].real()   < zc.real() - (w/2.0)) ||
      (ze[seg].imag()   > zc.imag() + (w/2.0)) ||
      (ze[seg].imag()   < zc.imag() - (w/2.0)) ||
			(ze[seg+1].real() > zc.real() + (w/2.0)) ||
      (ze[seg+1].real() < zc.real() - (w/2.0)) ||
      (ze[seg+1].imag() > zc.imag() + (w/2.0)) ||
      (ze[seg+1].imag() < zc.imag() - (w/2.0)))     {return false;}
  return true;
}
//------------------------------------------------------------------------
bool CStringElem::SegIsInCircle(const cmplex &zc,const double r,const  int seg) const{
	if ((abs(ze[seg]-zc)>r) || (abs(ze[seg+1]-zc)>r)) {return false;}
	else                                              {return true;}
}
//------------------------------------------------------------------------
bool CStringElem::SegPartInCir   (const cmplex &zc, const double r,const int seg) const{
	double a,b,c;
  //solve quadratic
	a = pow(ze[seg+1].real() - ze[seg].real(),2)+ pow(ze[seg+1].imag() - ze[seg].imag(),2) ;

  b = 2.0*((ze[seg+1].real() - ze[seg].real())*(ze[seg].real() - zc.real()) + 
					 (ze[seg+1].imag() - ze[seg].imag())*(ze[seg].imag() - zc.imag())); 
      
	c = pow(     zc.real(),2) + pow(     zc.imag(),2) + 
			pow(ze[seg].real(),2) + pow(ze[seg].imag(),2) - 
			2.0*(zc.real()* ze[seg].real() + 
					 zc.imag()* ze[seg].imag()) - pow(r,2);

	if ((abs(ze[seg  ]-zc)<=r) || 
			(abs(ze[seg+1]-zc)<=r) ||
			(b*b-4.0*a*c>=0.0)){
    return true;
	}
	else {
		return false;
	}
} 
//************************************************************************
bool CStringElem::IsString() const {return true;}
//***********************************************************************
bool CStringElem::IsInside(const cmplex &z) const{
//Should OPTIMIZE!
	double sum(0.0);
	cmplex Z;
  if (closed){
		for (int i=0;i<NLines; i++){
		/*	if ((abs((z.imag()-ze[i].imag())*(ze[i+1].real()-ze[i].real())-(z.real()-ze[i].real())*(ze[i+1].imag()-ze[i].imag()))<REALSMALL) &&
					( (abs(z-ze[i])+abs(z-ze[i+1]))<=abs(ze[i]-ze[i+1]))
					){return true;}*/ //if on boundary
			Z=GlobalToLocal(z,i);
			sum+=log((Z-1.0)/(Z+1.0)).imag();
    }
	  return (fabs(sum)<REALSMALL ? false : true);
	}
	else {return false;}
}
//-------------------------------------------------------------------------
double CStringElem::Area(bool &clockwise) const{ 
  if (closed) {
    //borrowed from Paul Bourke <http://astronomy.swin.edu.au/pbourke/>
    int    j;
    double area(0.0);

    for (int i=0;i<NLines;i++) {
      j = (i+1)%NLines;
      area+= ze[i].real()*ze[j].imag()-ze[i].imag()*ze[j].real();
    }
    area/=2.0;
    clockwise=((area>=0.0) ? false : true);
    return fabs(area);
	}
	else { //used most on closed elements
    clockwise=false; return 0.0;
  }
}
//-------------------------------------------------------------------------
double CStringElem::GetArea      () const{bool junk;return Area(junk);}

/*bool CStringElem::IsInside(const cmplex &z) const{
	double test(0.0);
	int right(0);
	int left(0);
  if (closed){
		for (int i=0;i<NLines; i++){
      test=(z.imag()-ze[i].imag())*(ze[i+1].real()-ze[i].real())-(z.real()-ze[i].real())*(ze[i+1].imag()-ze[i].imag());
			if (test <0) {left++;}
			else if (test>0){right++;}
			else if (test==0){return IsInside(z+MOVE_DIST*IM*abs(ze[i+1]-ze[i]));}
		}
	  return ((left-right%2==0) ? false : true);
	}
	else {return false;}
}*/
/*bool CStringElem::IsInside(const cmplex &z) const{
  if (!closed){return false;}
  else {
    //borrowed from Wm. Randolph Franklin <wrf@ecse.rpi.edu>
    double y(z.imag());
		double x(z.real()); 
		bool   in(false);
    int    i,j;
    for (i=0, j=NLines-1; i<NLines; j=i++) {
      if ((((ze[i].imag()<=y) && (y<ze[j].imag())) ||
	         ((ze[j].imag()<=y) && (y<ze[i].imag()))) &&
         (x<(ze[j].real()-ze[i].real())*(y-ze[i].imag())/
	          (ze[j].imag()-ze[i].imag())+ze[i].real()))
     in=!in;}
    return in;
  }
}*/
	/*
	Old GetSegment Leakage Function
		if ((LJumpCoeff!=NULL) && (LFarFldCoeff!=NULL)){
  
		double  logend1,logend2,leak(0.0),end1(0.0),end2(0.0);
	  cmplex  Z((z-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i])));
		cmplex  jump(0.0);
		double  modulus(abs(Z));

		if (Leakltype==LINESINK){
			cmplex  ZBH((pLayer->GetBlackHole()-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i])));
			logend1=log((Z-1.0)/abs(ZBH-1.0)).real();
			logend2=log((Z+1.0)/abs(ZBH+1.0)).real();
			for (int n=0; n<=order; n++){
				end1+=LJumpCoeff[i][n].real();
				end2+=LJumpCoeff[i][n].real()*pow(-1.0,n);
			}
		}
		//far-field (without truncation)
		if (modulus>ISFAR){ 
			leak+=Outside(Z,(int)(FForder[i]*pow(modulus/ISFAR,FF_POWER)),LFarFldCoeff[i]).real();
			if (Leakltype==LINESINK) {
				leak+=IOVER2PI*(-(end1*logend1)+(end2*logend2));
			}
		}  
		//near-field
		else if ((abs(Z-e1)>MOVE_DIST) && (abs(Z-e2)>MOVE_DIST)){
			cmplex logterm(log(Zm1oZp1(Z)));
			for (int n=0; n<=order; n++){     
				jump+=conj(LJumpCoeff[i][n])*cClenshaw2(n,Z,f[n],logterm);
			}
      leak+=IOVER2PI*jump.real();
		  if (Leakltype==LINESINK){
				leak+=IOVER2PI *(-(end1*logend1)+(end2*logend2));
			}
  
		}
		//singularities
		else {leak=GetSegmentLeakage(i,z+MOVE_DIST*abs(ze[i+1]-ze[i]),c,t);}

		return leak/c;
	}
	else {return 0.0;}*/

/*
//***********************************************************************
//				Get Segment Complex Potential due to Leakage  
//***********************************************************************
cmplex CStringElem::GetSegmentLeakOm(const int i,const cmplex &z,const double &t) const{
	double  logend1,logend2,omega(0.0),modulus;
	double  end1(0.0),end2(0.0),end3(0.0),end4(0.0),end5(0.0),end6(0.0);
	cmplex  Z,ZBH,jump(0.0),jump2(0.0),jump3(0.0);
  cmplex *fZj;          
  int     n,j;
  double L(abs(ze[i]-ze[i+1]));
	if ((LJumpCoeff!=NULL) && (LFarFldCoeff!=NULL)){
    
		fZj=new cmplex [order+1];

		Z=  (           z-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
		ZBH=(pLayer->GetBlackHole()-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
		modulus=abs(Z);
		logend1=log((Z-1.0)/abs(ZBH-1.0)).real();
		logend2=log((Z+1.0)/abs(ZBH+1.0)).real();


		for (n=0; n<=order; n++){
			end1+=h(LJumpCoeff[i],n).real();
			end2+=h(LJumpCoeff[i],n).real()*pow(-1.0,n);
			end3+=h(LJumpCoeff[i],n).imag();
      end4+=h(LJumpCoeff[i],n).imag()*pow(-1.0,n);
			end5+=hh(LJumpCoeff[i],n).imag();
      end6+=hh(LJumpCoeff[i],n).imag()*pow(-1.0,n);  
		}

		//far-field (without truncation)
		if (modulus>ISFAR){ 
			//int nterms=(int)(min(FForder[i],
			//                 MIN_FAR_FIELD+(FForder[i]-MIN_FAR_FIELD)*(ISFAR/modulus)));
			//leak+=Outside(Z,FForder[i],LFarFldCoeff[i]).real();
			//if (Leakltype==LINESINK){leak+=IOVER2PI*(-(end1*logend1)+(end2*logend2));}
		}  
		//near-field
		else if ((abs(Z-e1)>MOVE_DIST) && (abs(Z-e2)>MOVE_DIST)){
			for (n=0; n<=order; n++){
				for (j=0; j<n; j++){fZj[j]=f[n][j];}
        fZj[n]=log(Zm1oZp1(Z));      
        jump +=conj( h(LJumpCoeff[i],n)).imag()*cClenshaw(n,Z,fZj);   //doublet term (without log)
        jump2+=conj( h(LJumpCoeff[i],n)).real()*cClenshaw(n,Z,fZj);   //dipole term 1 (without logs)
				jump3+=conj(hh(LJumpCoeff[i],n)).real()*cClenshaw(n,Z,fZj);   //dipole term 2 (without logs)
			}
			//doublet
			omega=(L*(Z-conj(Z))*IOVER2PI/16.0*(jump-(end1*logend1)+(end2*logend2))).real(); 
      //dipole 
			omega =((Z-conj(Z))*L*IOVER2PI/16.0*(jump2-(end3*logend1)+(end4*logend2))).real();
			omega+=(            L*IOVER2PI/8.0 *(jump2-(end5*logend1)+(end6*logend2))).real();
      omega+=(           -L*IOVER2PI/8.0 *(end3*((Z-1.0)*logend1-(Z-1.0)))).real();
			omega+=(            L*IOVER2PI/8.0 *(end4*((Z+1.0)*logend2-(Z+1.0)))).real();

		  //add multipoles for linesinks 
			//add taylor series
		}
		//singularities
		//else {leak=GetSegmentLeakage(i,z+MOVE_DIST*abs(ze[i+1]-ze[i]),c,t);}

		delete [] fZj;


		return omega;// /c;
	}
	else {return 0.0;}
}*/
//***********************************************************************
//				Get Segment Complex Discharge due to Leakage
//***********************************************************************
/*cmplex CStringElem::GetSegmentLeakW(const int i,const cmplex &z,const double &t) const{
	return 0.0;
}
*/
//***********************************************************************
//				Get Segment Log Centroid
//***********************************************************************
/*cmplex  CStringElem::GetSegmentLogCentroid(const int i) const{
	cmplex centroid((ze[i+1]+ze[i])/2.0);

	if (ltype==LINESINK){
    double  end1(0.0),end2(0.0),sum(0.0);
		int     n; 

		sum=-1.0;
    for (n=0; n<=order; n++){
			sum*=-1.0;
			end1-=JumpCoeff[i][n];  
			end2+=JumpCoeff[i][n]*sum;//sped up this pow
		}
    centroid=(ze[i+1]*end1+ze[i]*end2)/(end1+end2);
		cout <<"CStringElem::GetSegmentLogCentroid:  "<<end1<< " "<<end2 << " "<<"Q:" << end2-end1<< centroid<<endl;
  }

	return centroid;
}*/
//-----------------------------------------------------------------------
/*double CStringElem::GetHeadJump(int i,double X,double t) const{
  double jump(0.0);
  for(int n=0; n <=order; n++){jump+=LJumpCoeff[i][n].imag()*cos(n*acos(X));}
  return jump;
}
//-----------------------------------------------------------------------
double CStringElem::GetHeadGradJump(int i,double X,double t) const{
  double sum, length(abs(ze[i+1]-ze[i]));
	cmplex jump(0);

  for(int n=0; n <=order; n++){  
    sum=0.0;
    for(int k=n-1; k>=0; k-=2){
			if (k!=0){sum+=2.0*n*cos(k*acos(X));}
			else     {sum+=    n*cos(k*acos(X));}
		} 
    jump+=conj(LJumpCoeff[i][n])*2.0/length*sum;
	}
  return jump.real();
}*/