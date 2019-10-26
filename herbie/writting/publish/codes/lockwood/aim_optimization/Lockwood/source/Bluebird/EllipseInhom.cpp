//EllipseInhom.cpp
#include "EllipseElem.h"
/***********************************************************************
                           CEllipseInhom
************************************************************************/
CEllipseInhom::CEllipseInhom(){
 kin=0.0;
}
//----------------------------------------------------------------------
CEllipseInhom::CEllipseInhom(char			       *Name,
														 const CSingleLayerABC *pLay,  
														 const double     cond,  
														 const cmplex		  zcen, 
														 const double     MajorAxis,
														 const double		  MinorAxis,
														 const double     angle,
														 const int			  prec):        
							CEllipseElem(Name,pLay,zcen,MajorAxis,MinorAxis,angle,prec){
	kin=cond;
	Q=0;
	CondZone = new CEllPropZone(hydraulic_conductivity, zc, a,b,angle, cond);
}
//----------------------------------------------------------------------
CEllipseInhom::~CEllipseInhom(){}

/***********************************************************************
                           ACCESSOR FUNCTIONS
************************************************************************/
double CEllipseInhom::GetK() const{return kin;}
CPropZone  *CEllipseInhom::GetCondZone() const{return CondZone;}
/***********************************************************************
                           SolveItself
***********************************************************************/
void CEllipseInhom::SolveItself(double &change,double &objective,double t){
	double *PotOther;
	int     m,n;
	cmplex  sum,tmp;
	double  kout, psi, areal,aimag, AvgPot;

	change=objective=0.0;

  kout=pLayer->GetCond(zctrl[0]+(zctrl[0]-zc)*a*REALSMALL); 


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

	//explicitly calculate coefficients------------------------------
	//1st coefficient

  tmp=AvgPot*(kin-kout)/kout/2.0;
	upperswap(change,abs(InCoeff[0]-tmp));
	InCoeff[0]=tmp;

	//all other coefficients
  for(n=1; n<order; n++){

    sum=0.0; 

	  for(m=0;m<nellcontrol;m++){
			psi =-PI+2.0*PI*(double)(m)/(double)(nellcontrol);
			sum+=PotOther[m]*cmplex(cos(n*psi),sin(n*psi));
		}

		areal=sum.real()*(kin-kout)/(double)(nellcontrol)/(kout*cosh(n*eta0)+kin*sinh(n*eta0));
		aimag=sum.imag()*(kin-kout)/(double)(nellcontrol)/(kout*sinh(n*eta0)+kin*cosh(n*eta0));

		tmp=cmplex(areal,aimag);
		upperswap(change,abs(InCoeff[n]-tmp));
		InCoeff[n]=tmp;          	   

		//cout <<change<<" "<<abs(InCoeff[n]-tmp)<<endl;


		if (2.0*(double)(n)*eta0<MAXEXP) {
			OutCoeff[n]=conj(InCoeff[n])-InCoeff[n]*exp(2.0*(double)(n)*eta0);
		}
		else{
			OutCoeff[n]=0.0; InCoeff [n]=0.0;
		}
	}
	Q=0.0;

	//evaluate objective
	double error;
	for (m=0; m<nellcontrol; m++){
	  error=kin *(PotOther[m]+GetDischargePotential(zctrl[m]+REALSMALL*(zctrl[m]-zc)*a,t).real())-
			    kout*(PotOther[m]+GetDischargePotential(zctrl[m]-REALSMALL*(zctrl[m]-zc)*a,t).real());
		upperswap(objective,error);
	} 
  //change/=pLayer->GetDeltaPot();


  delete [] PotOther;
}
/***********************************************************************
                           GetMaxError
***********************************************************************/
double CEllipseInhom::GetMaxError(const double &t) const{
	double h1,h2,error(0);
	for (int m=0; m<nellcontrol; m+=max((int)(TEST_ERROR_RATIO*nellcontrol),1)){
		h1=pLayer->GetHead(zctrl[m]+REALSMALL*(zctrl[m]-zc),t);
		h2=pLayer->GetHead(zctrl[m]-REALSMALL*(zctrl[m]-zc),t);
		upperswap(error,fabs(h1-h2)); //absolute head error
		//upperswap(error,fabs(h1-h2)/max(h1,1.0)); //relative head error
	}
	return error;
}

/***********************************************************************
                           PARSE
************************************************************************
Format:
	string "EllipticalInhom"
	double x1 double y1 double K double A double B double angle (degrees)
  & [int precision]
----------------------------------------------------------------------*/
CEllipseInhom  *CEllipseInhom::Parse(ifstream &input,int &l,CSingleLayerABC *pLay,char * Name){



	CEllipseInhom *pEllInhom;
  bool					 eof(false),done(false);
  double				 thiscond(0);
	cmplex				 thisz;
	double         thisa,thisb,thisang;
	int						 Len,thisprec;
	char					*s[MAXINPUTITEMS];

	pEllInhom=NULL;
	if (parserdebug) {cout << "Elliptical Inhom element"<<endl;}   
																												 eof=TokenizeLine(input,s,Len); l++; 
	do {
		if (Len==6) {
      thisz   =s_to_c(s[0],s[1]);
			thiscond=s_to_d(s[2]);
			thisa   =s_to_d(s[3]);
			thisb   =s_to_d(s[4]);
			thisang =s_to_d(s[5]);                            
			thisang =PI/180*thisang;                           eof=TokenizeLine(input,s,Len); l++; 
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			pLay->UpdateExtents(thisz,thisa);
			if (Len==2){thisprec= s_to_i(s[1]);}
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pEllInhom = new CEllipseInhom(Name,
																		pLay,
																		thiscond,
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
void CEllipseInhom::WriteOutput(const double &t) const{
	//no output yet
}