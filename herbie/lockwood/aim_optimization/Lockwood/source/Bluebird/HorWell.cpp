#include "HorWell.h"

/************************************************************************
                           CHorizontalWell
************************************************************************/
//Constructor
CHorizontalWell::CHorizontalWell(char						 *Name, 
																 const CSingleLayerABC *pLay, 
																 const cmplex		 *points,
																 const int				NumOfLines,
																 const double		  prate,
																// double  cond,
																// double  thickness,
																 const int			  prec)  
            :CStringElem(Name,pLay,points,NumOfLines,false,LINESINK,prec){
   Q=prate;
	 //k=cond;
	 //thick=thickness;
}
//-----------------------------------------------------------------------
CHorizontalWell::~CHorizontalWell(){}
/************************************************************************
                           SolveItself
***********************************************************************/
void CHorizontalWell::SolveItself(double &change, double &objective,double t){
  double  error(0.0),maxobjective(0.0),maxchange(0.0),relax(1.0);
  double  Qxother;
	cmplex  unitX;
  int     i,m;

  
  for (i=0; i<NLines; i++){
    unitX=(ze[i+1]-ze[i])/abs(ze[i+1]-ze[i]);

		disabled[i]=true;

		//Solve Horizontal Well Element
    //******************************************************************
    //set right hand side for system of equations- 
    for (m=0; m<nlinecontrol; m++){
      
      //identify non-segment tangential flux at control points 
      Qxother=((pLayer->GetW(zctrl[i][m],t))*unitX).real();

			rhs[m]=-Qxother;

    }

		disabled[i]=false;

    //solve for jump coefficients
		//SetConstraints(i,0,0,1,0,UNCONSTRAINED,0.0,0.0);
		SetConstraints(i,0,0,1,0,NET_DISCHARGE,Q,0.0);
    GenSolve(JumpCoeff[i],i,LINESINK,false,relax,objective,change);
		upperswap(maxchange   ,change   );
		upperswap(maxobjective,objective);
		SetFarFieldCoeff(i);
		if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){
			pBlocks[i]->Update(myBlockIDs[i],i,t);
		}
  }

  change=maxchange;
  objective=maxobjective;
}
//------------------------------------------------------------------------
double CHorizontalWell::GetMaxError(const double &t) const{
	double error(0.0),AqHead;
	double maxhead(-ALMOST_INF);
	double minhead( ALMOST_INF);
	for (int i=0; i<NLines; i++){
		for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
			AqHead=pLayer->GetHead(zctrl[i][m],t); //Global head
			upperswap(maxhead,AqHead);
			lowerswap(minhead,AqHead);
		}
	}
	//return maxhead-minhead;//absolute head error
	return 0.0;
}
/***********************************************************************
              Read/Write
***********************************************************************/
void CHorizontalWell::WriteOutput(const double &t) const{
  //write errors file
  double error(0.0),maxhead,minhead,AqHead,B;
  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
  
	for (int i=0; i<NLines; i++){
		B=pLayer->GetBase(int(nlinecontrol/2));
		minhead=1e99;maxhead=-1e99;
		for (int m=0; m<nlinecontrol; m++){
			AqHead=pLayer->GetHead(zctrl[i][m],t); //Global head
			upperswap(maxhead,AqHead);
			lowerswap(minhead,AqHead);
		}
		error=maxhead-minhead;
		ERRORS <<zctrl[i][int(nlinecontrol/2)].real()<<","<<zctrl[i][int(nlinecontrol/2)].imag()<<",";
		ERRORS <<HEAD_ERROR_TAG    <<","<<0<<",";
		ERRORS <<error<<","<<error/(0.5*(maxhead+minhead)-B)<<endl;
	}
	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "HWell", string name 
  double Pumping Rate
	{double x double y}x(numlines+1)
  & [int precision] 
----------------------------------------------------------------------*/
CHorizontalWell *CHorizontalWell::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){


 	CHorizontalWell  *pHorWell;
  
	bool     done(false);
  double   thisQ(0);
	int      Len,thisprec,nlines(0),i;
  cmplex	 stringp[MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Horizontal Well "<<endl;}  
	
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if      (Len==1) {thisQ=s_to_d(s[0]);}
  else             {ImproperFormat(s,l); return NULL;}

  done=false; 
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CHorizontalWell::Parse- too many lines in horizontal well",TOO_MANY);}
    if ((Len==2) && (strcmp(s[0],"&"))){
      stringp[nlines]=s_to_c(s[0],s[1]); nlines++;       	
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pHorWell = new CHorizontalWell(Name,
																		 pLay,
																		 stringp,
																		 nlines-1,
																		 thisQ,
																		 thisprec); 
			done=true;
		}

    else             {ImproperFormat(s,l); return NULL;}
	} while (!done);
	
	return pHorWell;
}


