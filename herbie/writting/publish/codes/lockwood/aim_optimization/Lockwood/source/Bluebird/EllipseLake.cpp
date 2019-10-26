//EllLake.cpp
//EllipseInhom.cpp
#include "EllipseElem.h"
/***********************************************************************
                           CEllLake
************************************************************************/
CEllLake::CEllLake(){
 elev=0.0;
}
//----------------------------------------------------------------------
CEllLake::CEllLake(char			       *Name,
														 const CSingleLayerABC *pLay,  
														 const double     head,  
														 const cmplex		  zcen, 
														 const double     MajorAxis,
														 const double		  MinorAxis,
														 const double     angle,
														 const int			  prec):        
							CEllipseElem(Name,pLay,zcen,MajorAxis,MinorAxis,angle,prec){
	elev=head;

	flownode=CFlowNode::AddFlowNode((CAnalyticElem*)(this),zcen,UPSTREAM,0.0,outflow_ID);
}
//----------------------------------------------------------------------
CEllLake::~CEllLake(){}

/***********************************************************************
                           ACCESSOR FUNCTIONS
************************************************************************/
double CEllLake::GetElev() const{return elev;}
//-----------------------------------------------------------------------
bool   CEllLake::HasFlux() const {return true;}
/***********************************************************************
                           SolveItself
***********************************************************************/
void CEllLake::SolveItself(double &change,double &objective,double t){
	double *PotOther;
	int     m,n;
	cmplex  tmp,sum;
	double  AvgPot,DesiredPot,K,B,T,psi;

	change=objective=0.0;

  K=pLayer->GetCond(zc);
	B=pLayer->GetBase(zc);
	T=pLayer->GetThick(zc);


	//Get Potential from other elements along boundary--------------
	PotOther=new double [nellcontrol];
	for (m=0; m<nellcontrol;m++){
		PotOther[m]=pLayer->GetDischargePotential(zctrl[m],t).real()-GetDischargePotential(zctrl[m],t).real();
	}
	
	AvgPot=0.0;
	for (m=0; m<nellcontrol;m++){
		AvgPot+=PotOther[m];
	}	
	AvgPot/=(double)(nellcontrol);
	
  //calculate the desired potential value	
	DesiredPot=ConvertToPotential(elev-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());

  Q=(DesiredPot-AvgPot)/(eta0-GlobalToLocal(pLayer->GetBlackHole()).real());

	InCoeff [0]=(DesiredPot-AvgPot)/2.0;

	//all other coefficients
  for(n=1; n<order; n++){

    sum=0.0; 

	  for(m=0;m<nellcontrol;m++){
			psi =-PI+2.0*PI*(double)(m)/(double)(nellcontrol);
			sum+=PotOther[m]*cmplex(cos(n*psi),sin(n*psi));
		}

    OutCoeff[n]=-2.0*exp(n*eta0)*sum/(double)(nellcontrol); 
		//InCoeff[n]=-2.0*exp(n*eta0)*sum/(double)(nellcontrol); 
		//InCoeff[n]=0.5*exp(n*eta0)*(sum)/(double)(nellcontrol); //Close...
		double areal,aimag;
		areal=sum.real()/(double)(nellcontrol)/(cosh(n*eta0)+sinh(n*eta0));
		aimag=sum.imag()/(double)(nellcontrol)/(sinh(n*eta0)+cosh(n*eta0));

		InCoeff[n]=-2.0*cmplex(areal,aimag);

		//cout << eta0<<" "<<c<< " "<<a<<endl;

		//InCoeff[n]=conj(OutCoeff[n]);//-OutCoeff[n]*exp(2.0*(double)(n)*eta0);
  }


	//evaluate objective
	double error;
	for (m=0; m<nellcontrol; m++){
	  error=PotOther[m]+GetDischargePotential(zctrl[m],t).real()-DesiredPot;
		upperswap(objective,error);
	} 
  //change/=pLayer->GetDeltaPot();

	//update baseflow downstream
	flownode->UpdateFluxes(GetNetDischarge(t),outflow_ID,t);

  delete [] PotOther;
}
/***********************************************************************
                           GetMaxError
***********************************************************************/
double CEllLake::GetMaxError(const double &t) const{
	double h,error(0);
	double B=pLayer->GetBase(zc);
	for (int m=0; m<nellcontrol; m+=max((int)(TEST_ERROR_RATIO*nellcontrol),1)){
		h=pLayer->GetHead(zctrl[m],t);
		upperswap(error,fabs(h-elev)); //absolute head error
		//upperswap(error,fabs(h-elev)/max((elev-B),1.0)); //relative head error
	}
	return error;
}

/***********************************************************************
                           PARSE
************************************************************************
Format:
	string "EllipticalLake"
	double x1 double y1 double elev double A double B double angle (degrees)
  & [int precision]
----------------------------------------------------------------------*/
CEllLake  *CEllLake::Parse(ifstream &input,int &l,CSingleLayerABC *pLay,char * Name){



	CEllLake *pEllInhom;
  bool					 eof(false),done(false);
  double				 thiselev(0);
	cmplex				 thisz;
	double         thisa,thisb,thisang;
	int						 Len,thisprec,nlines(0);
	char					*s[MAXINPUTITEMS];

	pEllInhom=NULL;
	if (parserdebug) {cout << "Elliptical Inhom element"<<endl;}   
																												 eof=TokenizeLine(input,s,Len); l++; 
	do {
		if (Len==6) {
      thisz   =s_to_c(s[0],s[1]);
			thiselev=s_to_d(s[2]);
			thisa   =s_to_d(s[3]);
			thisb   =s_to_d(s[4]);
			thisang =s_to_d(s[5]);                            
			thisang =PI/180*thisang;                           eof=TokenizeLine(input,s,Len); l++; 
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			pLay->UpdateExtents(thisz,thisa);
			if (Len==2){thisprec= s_to_i(s[1]);}
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pEllInhom = new CEllLake(Name,
															 pLay,
															 thiselev,
															 thisz,
															 thisa,
															 thisb,
															 thisang,
															 thisprec); 
			done=true;
		}
		else if(Len==0) {                                    eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while (!done);

  if (eof) {return NULL;}
	else     {return pEllInhom;}
}

/***********************************************************************
                           WriteOutput
***********************************************************************/
void CEllLake::WriteOutput(const double &t) const{
	//no output yet
}