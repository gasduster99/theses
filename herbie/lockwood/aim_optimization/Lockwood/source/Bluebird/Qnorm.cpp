#include "Qnorm.h"

/************************************************************************
                           CQnorm
************************************************************************/
//Constructor
CQnorm::CQnorm(char									 *Name, 
							 const CSingleLayerABC *pLay, 
							 const cmplex					 *points,
							 const int							NumOfLines,
							 const double					 *discharge, 
							 const int							prec)
       :CStringElem(Name,pLay,points,NumOfLines,false,DOUBLET,prec){

  int i;
	
  //allocate memory for dynamic arrays
  pQ=           new double  [NLines+1];
  Qctrl=        new double *[NLines];
  for(i=0;i<NLines;i++){
    Qctrl[i]=   new double  [nlinecontrol];
  }
	
  //copy discharge array,create matrices of desired discharges 
  for (i=0;i<=NLines;i++){pQ[i]=discharge[i];}
  for (i=0;i< NLines;i++){
    Interpolate(Qctrl[i],nlinecontrol,pQ[i],pQ[i+1]);
  }
}
//-----------------------------------------------------------------------
CQnorm::~CQnorm(){
	if (globaldebug){cout <<"   DESTROYING QNORM "<<name<<endl;}
	for(int i=0;i<NLines;i++){
		delete [] Qctrl[i];
	}
	delete [] Qctrl;
	delete [] pQ;
}
/************************************************************************
                           SolveItself
***********************************************************************/
void CQnorm::SolveItself(double &change, double &objective,double t){
  double maxobjective(0.0),maxchange(0.0),K,relax;
  cmplex WOther;
  int    i,m;

  
  for (i=0; i<NLines; i++){

		K=pLayer->GetCond(zctrl[i][int(nlinecontrol/2)]);

		disabled[i]=true;

    //set right hand side for system of equations (the specified normal flux)
    for (m=0; m<nlinecontrol; m++){

      //identify Total,Local normal flux at control points
      WOther=pLayer->GetW(zctrl[i][m],t);

      //set rhs
      rhs[m]=Qctrl[i][m]-(conj(WOther)*abs(ze[i+1]-ze[i])/(ze[i+1]-ze[i])).imag();
    }
		relax=1.0;

		disabled[i]=false;    
		
		//solve for jump coefficients
		SetConstraints(i,0,0,1,0,ZERO_AT_END,0.0,0.0);
    GenSolve(JumpCoeff[i],i,ltype,false,relax,objective,change);
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
double CQnorm::GetMaxError(const double &t) const{
	double error(0.0);

	return error;
}
//***********************************************************************

//***********************************************************************
void CQnorm::WriteOutput(const double &t) const{
  //write errors file
  double error(0);
  /*ofstream ERRORS; //TMP DEBUG
	ERRORS.open("errors.csv", ios::app);
	
	for (int i=0; i<NLines; i++){
		for (int m=0; m<nlinecontrol; m++){
		  ERRORS<<zctrl[i][m].real()<<","<<zctrl[i][m].imag()<<",";
			ERRORS<<FLUX_ERROR_TAG<<","<<0<<",";
			ERRORS<<error<<",";
			ERRORS<<error<<endl;
		}
	}
	ERRORS.close();*/
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "Qnorm", string name 
  {double x double y double Q}x(numlines+1)
  &[int precision]
----------------------------------------------------------------------*/
CQnorm *CQnorm::Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name){

	CQnorm  *pQnorm=NULL;

  bool     done(false);
	int      Len,thisprec,nlines(0),i;
  cmplex	 stringp[MAXLINES];
  double   Qp     [MAXLINES];

	char    *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Q Normal "<<endl;} 
	
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CQnorm::Parse- too many lines in normal discharge element",TOO_MANY);}
		if (Len==3){
      stringp[nlines]=s_to_c(s[0],s[1]);
      Qp     [nlines]=s_to_d(s[2]     ); 
			nlines++;      	
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pQnorm = new CQnorm(Name,
													pLay,
													stringp,
													nlines-1,
													Qp,
													thisprec);
			done=true;                                                                  }
    else           {ImproperFormat(s,l); return NULL;                             }
	} while (!done);

  return pQnorm;
}
//******************************************************************

