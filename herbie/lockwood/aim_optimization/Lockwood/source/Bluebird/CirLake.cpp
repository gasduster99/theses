#include "Circle.h"

/************************************************************************
*************************************************************************
                           CCirLake
*************************************************************************
************************************************************************/

/***********************************************************************
							 CONSTRUCTOR
***********************************************************************/
CCirLake::CCirLake     (char						*Name, 
												const CSingleLayerABC *pLay,
												const cmplex		 zc,
												const double		 rad, 
												const double		 head, 
												const int				 prec):
          CCircularElem(Name,pLay,zc,rad,prec) {
	elev=head;        
	Q=0.0;

	flownode=CFlowNode::AddFlowNode((CAnalyticElem*)(this),zc,UPSTREAM,0.0,outflow_ID);
}

/***********************************************************************
							 ACCESSORS
***********************************************************************/
double   CCirLake::GetHead() const {return elev;}
//-----------------------------------------------------------------------
bool     CCirLake::HasFlux() const {return true;}
/***********************************************************************
							 SolveItself
***********************************************************************/
void     CCirLake::SolveItself(double &change, double &objective, const double t){
  int    m,n;
  double AvgPot,oldQ;
	cmplex sum, temp2;
  double Potential[MAXCIRCONTROL],angle,error,error2,maxobjective(0.0),maxchange(0.0),B,K,T,DesiredPot;

  maxobjective=0.0;
  B=pLayer->GetBase (zcen);
	T=pLayer->GetThick(zcen);
	K=pLayer->GetCond (zcen);
	DesiredPot=ConvertToPotential(elev-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());

	//get potential along boundary, identify Q
	oldQ=Q;

	AvgPot=0.0;
	for (m=0; m<ncircontrol; m++){
		Potential[m]=(pLayer->GetDischargePotential(zctrl[m],t)-GetDischargePotential(zctrl[m],t)).real(); 
		AvgPot+=Potential[m]/(double)(ncircontrol);
	}
	upperswap(maxchange,abs(InCoeff[0]-DesiredPot+AvgPot));

	Q          =2.0*PI*(DesiredPot-AvgPot)/log(R/abs(pLayer->GetBlackHole()-zcen));
  OutCoeff[0]=0.0;
	InCoeff [0]=DesiredPot-AvgPot;                          //0.5*Q/pi*log(R/abs(pLayer->GetBlackHole()-zcen));

	upperswap(maxchange,fabs(oldQ-Q));
  
	//identify remaining coeff
	for (n=1;n<order;n++){
		sum=0.0;
		for (m=0; m<ncircontrol; m++){
		  angle=2.0*PI*(double)(m)/(double)(ncircontrol);
      sum+=(DesiredPot-AvgPot-Potential[m])*
				      cmplex(cos((double)(n)*angle),sin((double)(n)*angle));
		}
		temp2=2.0*sum/((double)(ncircontrol));
		upperswap(maxchange,abs(temp2-OutCoeff[n]));
		OutCoeff[n]=temp2;
		InCoeff [n]=conj(temp2);
	}

	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,NOT_A_STRING,t);
	}
  
	//evaluate objective
	for (m=0; m<ncircontrol; m++){
	  error =(Potential[m]+GetDischargePotential(zctrl[m]-(zctrl[m]-zcen)*REALSMALL,t).real())-DesiredPot;
	  error2=(Potential[m]+GetDischargePotential(zctrl[m]+(zctrl[m]-zcen)*REALSMALL,t).real())-DesiredPot;
		upperswap(maxobjective,error);
		upperswap(maxobjective,error2);
	} 

  change   =maxchange;
	objective=maxobjective;

	//update baseflow downstream
	flownode->UpdateFluxes(GetNetDischarge(t),outflow_ID,t);
}
//-----------------------------------------------------------------------
double CCirLake::GetMaxError(const double &t) const{
	double h,error(0);
	double B=pLayer->GetBase(zcen);
	for (int m=0; m<ncircontrol; m+=max((int)(TEST_ERROR_RATIO*ncircontrol),1)){
		h=pLayer->GetHead(zctrl[m],t);
		upperswap(error,fabs(h-elev)); //absolute head error
		//upperswap(error,fabs(h-elev)/max((elev-B),1.0)); //relative head error
	}
	return error;
}
//****************************************************************************
void CCirLake::WriteOutput(const double &t) const{
	double error;

  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);

	double B=pLayer->GetBase (zcen);
	for (int m=0; m<ncircontrol; m++){
		error=pLayer->GetHead(zctrl[m],t)-elev; 
		ERRORS << zctrl[m].real() <<","<<zctrl[m].imag()<<","; //coordinate
		ERRORS << HEAD_ERROR_TAG  <<","<<0<<",";
		ERRORS << error                                 <<","; //absolute head error
		ERRORS << error/max((elev-B),1.0)               <<endl;//relative head error

	} 
	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
	string "CirLake", string name 
	double x double y double elev double radius
	& [int precision] 
----------------------------------------------------------------------*/
CCirLake *CCirLake::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){
  

	CCirLake *pCirLake;
  bool      eof(false),done(false);
  double    thisradius,thishead(0);
	cmplex    thiszcen;
	int       Len,thisprec,nlines(0);
	char    *s[MAXINPUTITEMS];

	if (parserdebug) {
		cout << "Circular Lake element"<<endl;}              eof=TokenizeLine(input,s,Len); l++; 
	do {
		if (Len==4) {
      thiszcen  =s_to_c(s[0],s[1]);
			thishead  =s_to_d(s[2]);
			thisradius=s_to_d(s[3]);                           eof=TokenizeLine(input,s,Len); l++; 
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			pLay->UpdateExtents(thiszcen,thisradius);
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pCirLake = new CCirLake(Name,
														 pLay,
														 thiszcen,
														 thisradius,
														 thishead,
														 thisprec); 
			done=true;
		}
		else if(Len==0) {                                    eof=TokenizeLine(input,s,Len); l++;}
    else            {ImproperFormat(s,l); return NULL;}
	} while (!done);

  if (eof) {return NULL;}
	else     {return pCirLake;}
}
