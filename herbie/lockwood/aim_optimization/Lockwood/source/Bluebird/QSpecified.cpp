#include "QSpecified.h"

/************************************************************************
                           CQSpec
************************************************************************/
//Constructor
CQSpec::CQSpec(char						 *Name, 
							 const CSingleLayerABC *pLay, 
							 const cmplex		 *points,
							 const int				NumOfLines,
							 const double		 *discharge) 
       :CStringElem(Name,pLay,points,NumOfLines,false,LINESINK,1){
  int i;

	solved=false;
	
  //allocate memory for dynamic arrays
  pQ=           new double    [NLines+1];
  Qctrl=        new double   *[NLines];
	flownodes=    new CFlowNode*[NLines+1];
	outflow_ID=   new int       [NLines+1];
  for(i=0;i<NLines;i++){
    Qctrl[i]=   new double  [nlinecontrol];
  }
	
  //copy discharge array,create matrices of desired discharges 
  for (i=0;i<=NLines;i++){pQ[i]=discharge[i];}
  for (i=0;i< NLines;i++){
    Interpolate(Qctrl[i],nlinecontrol,pQ[i],pQ[i+1]);
  }
	//specify flow nodes
	int junk(0);
  for(i=0;i<NLines;i++){//nodes are upstream of river segments (outflow ID not assigned)
		flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],DOWNSTREAM,0.0,junk);
  }
  for(i=1;i<=NLines; i++){ //nodes are downstream of river segments
	  flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],UPSTREAM,0.0,outflow_ID[i]);
	}
}
//-----------------------------------------------------------------------
CQSpec::~CQSpec(){
	if (globaldebug){cout <<"   DESTROYING EXTRACTION "<<name<<endl;}
	for(int i=0;i<NLines;i++){
		delete [] Qctrl[i];
	}
	delete [] flownodes;
	delete [] outflow_ID;
	delete [] Qctrl;
	delete [] pQ;
}
void CQSpec::SetPrecision (const int Precision,int &order, double &fold){
	order=2; fold=1.0;
}
/************************************************************************
                           SolveItself
***********************************************************************/
void CQSpec::SolveItself(double &change, double &objective,double t){
  double maxobjective(0.0),maxchange(0.0),relax;
  int    i,m;

  if (!solved){
		for (i=0; i<NLines; i++){

			//set right hand side for system of equations (the specified normal flux)
			for (m=0; m<nlinecontrol; m++){
				rhs[m]=Qctrl[i][m];
				//cout << Qctrl[i][m]<<endl;
			}
			relax=1.0;
    
			//solve for jump coefficients
			SetConstraints(i,0,1,0,0,ZERO_AT_END,0.0,0.0);
			GenSolve(JumpCoeff[i],i,ltype,false,relax,objective,change);
			upperswap(maxchange   ,change   );
			upperswap(maxobjective,objective);

			SetFarFieldCoeff(i);
			if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){
				pBlocks[i]->Update(myBlockIDs[i],i,t);
			}
			solved=true;
		}

		change=maxchange;
		objective=maxobjective;
	}
	else {objective=change=0.0;}

	for (i=0; i<NLines; i++){
		//update baseflow downstream
		flownodes[i+1]->UpdateFluxes(flownodes[i]->GetOutflow(t)+GetSegmentDischarge(i,t),outflow_ID[i+1],t);
	}
}
//***********************************************************************
double CQSpec::GetMaxError(const double &t) const{
	return 0.0; //always true!
}

/***********************************************************************
              GetSegCumulativeFlux
***********************************************************************/
double CQSpec::GetSegCumulativeFlux(const int i, const cmplex &z, const double &t) const{
	cmplex Z;
	Z=GlobalToLocal(z,i);
	//cout <<"CRiver::GetSegCumulativeFlux : GOT SOMETHING"<<endl;
	if ((abs(Z.imag())< NEAR_FEATURE) && 
		  (Z.real()>=-1.0-NEAR_FEATURE) && 
			(Z.imag()<= 1.0+NEAR_FEATURE)){ //if this point is close enough to the river segment
	  return flownodes[i]->GetOutflow(t)+GetCumExtract(i,Z.real(),DIR_RIGHT,t); //cumulative baseflow
	}
	else {
		return 0.0;
	}
}
//***********************************************************************
void CQSpec::WriteOutput(const double &t) const{
  //write errors, extraction file
  double error(0),QyJump;
  int i, m;
  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);

	ofstream EXTRACT;
	EXTRACT.open("extract.csv", ios::app);
	for (i=0; i<NLines; i++){
		for (m=0; m<nlinecontrol; m++){
		  
			QyJump=GetWJump(i,X[m],0).imag(); 	
			
			ERRORS<<zctrl[i][m].real()<<","<<zctrl[i][m].imag()<<",";
			ERRORS<<FLUX_ERROR_TAG<<","<<0<<",";
			ERRORS<<(QyJump+Qctrl[i][m])/(max(fabs(Qctrl[i][m]),1.0))<<","; //relative err
			ERRORS<<(QyJump+Qctrl[i][m])<<endl;//absolute err
			
			EXTRACT <<zctrl[i][m].real()<<","<<zctrl[i][m].imag()             <<",";		//coordinate
			EXTRACT <<-QyJump                                                 <<",";		//point extraction
			//EXTRACT <<GetCumExtract(i,X[m],DIR_RIGHT,t)                       <<",";	//cumulative extraction
			EXTRACT <<flownodes[i]->GetOutflow(t)+GetCumExtract(i,X[m],DIR_RIGHT,t)<<endl;	//cumulative extraction			 
		
		}
	}	
	ERRORS.close();
  EXTRACT.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "QSpec", string name 
  {double x double y double Q}x(numlines+1)
  &[int precision]
----------------------------------------------------------------------*/
CQSpec *CQSpec::Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name){
  
	CQSpec  *pQSpec=NULL;

  bool     eof(false),done(false);
	int      Len,nlines(0),i;
  cmplex	 stringp[MAXLINES];
  double   Qp     [MAXLINES];

	char    *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Extraction Specified element"<<endl;}          	
	
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;  
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CQSpec::Parse- too many lines in discharge specified element",TOO_MANY);}
		if (Len==3){
      stringp[nlines]=s_to_c(s[0],s[1]);
      Qp     [nlines]=s_to_d(s[2]     ); nlines++;     	
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len>=1) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			pQSpec = new CQSpec(Name,
													pLay,
													stringp,
													nlines-1,
													Qp);
			done=true;                                                                               }
    else             {ImproperFormat(s,l); return NULL;                                        }
	} while (!done);

  return pQSpec;
}
//******************************************************************

