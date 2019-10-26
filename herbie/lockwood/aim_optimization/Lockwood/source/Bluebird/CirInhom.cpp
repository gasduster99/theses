#include "Circle.h"

/************************************************************************
*************************************************************************
                           CCirInhom
*************************************************************************
************************************************************************/

/***********************************************************************
							 CONSTRUCTOR
***********************************************************************/
CCirInhom::CCirInhom(char						 *Name, 
										 const CSingleLayerABC *pLay,
										 const cmplex		  zc,
										 const double			rad, 
										 const double			k, 
										 const int				prec) :
           CCircularElem(Name,pLay,zc,rad,prec) {
	cond=k;        
	Q=0.0;
	CondZone = new CCirPropZone(hydraulic_conductivity, zc, rad, cond);
}
//-----------------------------------------------------------------------
CCirInhom::~CCirInhom(){
  if (globaldebug){cout <<"   DESTROYING CIRINHOM "<<name<<endl;}
}
/***********************************************************************
							 ACCESSORS
***********************************************************************/
double      CCirInhom::GetKin() const{return cond;} 
//-----------------------------------------------------------------------
CPropZone  *CCirInhom::GetCondZone() const{return CondZone;} 
/***********************************************************************
							 SolveItself
***********************************************************************/
void     CCirInhom::SolveItself(double &change, double &objective, const double t){
  int     m,n;
  cmplex  sum(0.0),temp;
  double  kout,angle,error,maxobjective(0.0),maxchange(0.0);
	double *Potential;

  Potential= new double [ncircontrol];

	kout=pLayer->GetCond(zctrl[0]+cmplex(MOVE_DIST*R,0.0)); //just to right
  
	//get potential along boundary, identify first coeff---------------------------------
	for (m=0; m<ncircontrol; m++){
		Potential[m]=max((pLayer->GetDischargePotential(zctrl[m],t)-
			                        GetDischargePotential(zctrl[m],t)).real(),0.0); 

		sum+=Potential[m];
	}
	temp=(cond-kout)*sum/((double)(ncircontrol)*kout); 
  OutCoeff[0]=0.0;
	InCoeff [0]=temp;
  
	//identify remaining coeff-----------------------------------------------------------
	for (n=1;n<order;n++){
		sum=0.0;
		for (m=0; m<ncircontrol; m++){
		  angle=2.0*PI*(double)(m)/(double)(ncircontrol);
      sum+=(Potential[m]-OutCoeff[0])*cmplex(cos((double)(n)*angle),-sin((double)(n)*angle));
		}
		temp=(2.0*(cond-kout)*sum/((cond+kout)*(double)(ncircontrol)));
		if (abs(-conj(temp)-OutCoeff[n])>maxchange){maxchange=abs(-conj(temp)-OutCoeff[n]);}
		OutCoeff[n]=-conj(temp);
		InCoeff [n]=temp;
	}

	if ((pBlock!=NULL) && (pBlock->IsOn())){pBlock->Update(myBlockID,NOT_A_STRING,t);}
  
	//evaluate objective

	for (m=0; m<ncircontrol; m++){
	  error=cond*(Potential[m]+GetDischargePotential(zctrl[m]+REALSMALL*(zctrl[m]-zcen),t).real())-
			    kout*(Potential[m]+GetDischargePotential(zctrl[m]-REALSMALL*(zctrl[m]-zcen),t).real());
		upperswap(maxobjective,error);
	} 

  change   =maxchange/pLayer->GetDeltaPot();
	objective=maxobjective;
 //cout << OutCoeff[0] << " "<<R<< " "<<OutCoeff[1]<<endl;

  delete [] Potential; 
}
//-----------------------------------------------------------------------
double CCirInhom::GetMaxError(const double &t) const{
	double h1,h2,error(0);
	for (int m=0; m<ncircontrol; m+=max((int)(TEST_ERROR_RATIO*ncircontrol),1)){
		h1=pLayer->GetHead(zctrl[m]+REALSMALL*(zctrl[m]-zcen),t);
		h2=pLayer->GetHead(zctrl[m]-REALSMALL*(zctrl[m]-zcen),t);
		upperswap(error,fabs(h1-h2)); //absolute head error
		//upperswap(error,fabs(h1-h2)/max(h1,1.0)); //relative head error
	}
	return error;
}
//*******************************************************************************
void CCirInhom::WriteOutput(const double &t) const{

	ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	double headin,headout;
	int m;
	for (m=0; m<ncircontrol; m++){
	  headin =pLayer->GetHead(zctrl[m]+REALSMALL*(zctrl[m]-zcen),t);
		headout=pLayer->GetHead(zctrl[m]-REALSMALL*(zctrl[m]-zcen),t);
		ERRORS << zctrl[m].real() << ","<<zctrl[m].imag() << ","; //coordinate
		ERRORS << HEAD_ERROR_TAG  <<","<<0<<",";
		ERRORS << headin-headout                          << ",";	//absolute head error 
		ERRORS << fabs(headin-headout)/max(headout,1.0)   <<endl;	//relative head error
	} 
	ERRORS.close();

}
/***********************************************************************
                           PARSE
************************************************************************
Format:
	string "CirInhom", string name 
	double x double y double cond double radius
	& [int precision] 
----------------------------------------------------------------------*/
CCirInhom *CCirInhom::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){

	CCirInhom *pCirInhom=NULL;

  double     thisradius,thiscond(0);
	cmplex     thiszcen;
	int        Len,thisprec,nlines(0);
	char     *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Circular Inhom element"<<endl;}  
	
  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if (Len==4) {
    thiszcen=s_to_c(s[0],s[1]);
		thiscond=s_to_d(s[2]);
		thisradius=s_to_d(s[3]);         
	}
  else            {ImproperFormat(s,l); return NULL;}

  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
  if ((Len<=2) && (!strcmp(s[0],"&"))) {
			pLay->UpdateExtents(thiszcen,thisradius);
			if (Len==2){thisprec= s_to_i(s[1]);}
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pCirInhom = new CCirInhom(Name,
														    pLay,
														    thiszcen,
														    thisradius,
														    thiscond,
														    thisprec); 
	}
		
  return pCirInhom;
}

