#include "SurfaceDrain.h"
/************************************************************************
                           CSurfaceDrain
************************************************************************/
//Constructor
CSurfaceDrain::CSurfaceDrain(char						 *Name, 
														 const CSingleLayerABC *pLay, 
														 const cmplex		 *points,
														 const int				NumOfLines,
														 const double		 *heads, 
														 const int				prec) 
       :CStringElem(Name,pLay,points,NumOfLines,false,LINESINK,prec){

  int i;

	//allocate memory for dynamic arrays
  head =        new double        [NLines+1];
  headctrl=     new double       *[NLines];
	flownodes=    new CFlowNode    *[NLines+1];
	outflow_ID=   new int           [NLines+1];
  for(i=0;i<NLines;i++){
    headctrl[i]=     new double   [nlinecontrol];
  }
  //identify streamflow direction based upon 1st & last elevations
	//copy coordinate array, head array,in downstream direction 
  if (heads[0]<heads[NLines]){
		for(i=0;i<=NLines;i++){
			ze   [i]=points[NLines-i];
			head[i]=heads [NLines-i];
		}
	}
	else{
		for(i=0;i<=NLines;i++){
			head[i]=heads [i];
		}
	}
	
	int junk(0);
  for(i=0;i<NLines;i++){//nodes are upstream of river segments (outflow ID not assigned)
		flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],DOWNSTREAM,0.0,junk);
  }
  for(i=1;i<=NLines; i++){ //nodes are downstream of river segments
	  flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],UPSTREAM,0.0,outflow_ID[i]);
	}

  for (i=0; i<NLines; i++){
    Interpolate(headctrl[i],nlinecontrol,head[i],head[i+1]);
		Interpolate(zctrl   [i],nlinecontrol,ze  [i],  ze[i+1]); //if reordering has occured
  }
}
//-----------------------------------------------------------------------
CSurfaceDrain::~CSurfaceDrain(){
	if (globaldebug){cout <<"   DESTROYING SURFACE DRAIN "<<name<<endl;}
	delete [] head;
	delete [] flownodes;
	delete [] outflow_ID;
	for(int i=0;i<NLines;i++){
		delete [] headctrl[i];
	}
	delete [] headctrl;
}
//------STATIC MEMBERS---------------------------------------------------
double CSurfaceDrain::relax=1.0;
//-----------------------------------------------------------------------
void   CSurfaceDrain::SetRelaxation(double relaxation){relax=relaxation;}
/***********************************************************************
                           SolveItself
***********************************************************************/
void CSurfaceDrain::SolveItself(double &change, double &objective,double t){
  double  error(0.0),DesiredPot,maxobjective(0.0),maxchange(0.0),K,B,T,AqHead;
	double *PotOther;
  int     i,m;
	double  baseflow_in;

	PotOther=new double [nlinecontrol];
  
  for (i=0; i<NLines; i++){

		//Solve Surface Drain Element
    //******************************************************************
		K=pLayer->GetCond (zctrl[i][int(nlinecontrol/2)]);
		B=pLayer->GetBase (zctrl[i][int(nlinecontrol/2)]); 
		T=pLayer->GetThick(zctrl[i][int(nlinecontrol/2)]); 			
		
	  baseflow_in=flownodes[i]->GetOutflow(t);

		disabled[i]=true;

    //set right hand side for system of equations- 
    for (m=0; m<nlinecontrol; m++){

      //identify potential from other elements & sides at control points
			PotOther[m]=pLayer->GetDischargePotential(zctrl[i][m],t).real();

      //calculate the desired potential value	(for specified head)
			DesiredPot=ConvertToPotential(headctrl[i][m]-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());

      //calculate neccesary jump in potential
      rhs[m]=min(DesiredPot-PotOther[m],0.0);

    }

		//Gibbs phenomenon correction (via laplacian smoothing)
		double wt=0.5; // 1 - 2 -1
		//double wt=0.333; //  1 - 4- 1
		//double wt=0.25;	 //  1 - 6 -1

		for (int iter=1; iter<15; iter++){
			if (iter%2==0){
				for (m=1; m<nlinecontrol-1; m++){
					if ((rhs[m]*rhs[m-1]*rhs[m+1])==0){ rhs[m]=(1.0-wt)*rhs[m]+wt*0.5*(rhs[m-1]+rhs[m+1]);}
				}
			}
			else{
				for (m=nlinecontrol-2; m>=1; m--){
					if ((rhs[m]*rhs[m-1]*rhs[m+1])==0){rhs[m]=(1.0-wt)*rhs[m]+wt*0.5*(rhs[m-1]+rhs[m+1]);}
				}			
			}
		}

		disabled[i]=false;
    
		//solve for jump coefficients
		SetConstraints(i,0,0,0,1,ZERO_AT_END,0.0,0.0);
    GenSolve(JumpCoeff[i],i,LINESINK,false,relax,objective,change);
		upperswap(maxchange   ,change   );
		upperswap(maxobjective,objective);

		//evaluate error (necc?)
		for (m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
			AqHead=ConvertToHead(PotOther[m]+GetSegmentPotential(i,zctrl[i][m],t).real(),K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());
			if (AqHead+B>headctrl[i][m]){
				error=fabs(AqHead+B-headctrl[i][m])/fabs(headctrl[i][m]-B);}
			else{
				error=0.0;}
			upperswap(maxobjective,error);
		}
		
		SetFarFieldCoeff(i);

		if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){
			pBlocks[i]->Update(myBlockIDs[i],i,t);
		}
		
		//update baseflow downstream
		flownodes[i+1]->UpdateFluxes(baseflow_in+GetSegmentDischarge(i,t),outflow_ID[i+1],t);

  }

	delete [] PotOther;

  change=maxchange;
  objective=maxobjective;
}
//-----------------------------------------------------------------------
double CSurfaceDrain::GetMaxError(const double &t) const{
	static double error,AqHead;
	error=0.0;
	for (int i=0; i<NLines; i++){
		double B=pLayer->GetBase(zctrl[i][(int)(nlinecontrol/2)]);
		for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
			AqHead=pLayer->GetHead(zctrl[i][m],t); //Global head
			//upperswap(error,fabs(max(AqHead-headctrl[i][m],0.0))/(headctrl[i][m]-B));//relative head error
			upperswap(error,fabs(max(AqHead-headctrl[i][m],0.0))); //absolute head error
		}
	}
	return error;
}
/***********************************************************************
              READ/WRITE 
***********************************************************************/
void CSurfaceDrain::WriteOutput(const double &t) const{
	//write extraction file
  //write errors file
  double B;
  double AqHead,QyJump;

	ofstream EXTRACT;
  ofstream ERRORS;
	EXTRACT.open("extract.csv", ios::app);
	ERRORS.open ("errors.csv" , ios::app);
	//cout <<"OUTPUT (should be zero): "<<flownodes[0]->GetOutflow(t)<<endl;
	for (int i=0; i<NLines; i++){
		//cout << "*****"<<flownodes[i]->GetOutflow(t)<<endl;
		B=pLayer->GetBase(zctrl[i][(int)(nlinecontrol/2)]);
		for (int m=0; m<nlinecontrol; m++){
		  
			QyJump=GetWJump(i,X[m],0).imag(); 
			EXTRACT <<zctrl[i][m].real()<<","<<zctrl[i][m].imag()            <<",";		//coordinate
			EXTRACT <<QyJump                                                 <<",";		//point extraction
			EXTRACT <<flownodes[i]->GetOutflow(t)+GetCumExtract(i,X[m],DIR_RIGHT,t)<<endl;	//cumulative extraction
			 
			
			AqHead=pLayer->GetHead(zctrl[i][m],t);
			ERRORS << zctrl[i][m].real()<<","<<zctrl[i][m].imag()   << ",";	//coordinate
			ERRORS << HEAD_ERROR_TAG    <<","<<0<<",";
			ERRORS << fabs(min(AqHead-headctrl[i][m],0.0))/(headctrl[i][m]-B)<< ",";	//relative error	
			ERRORS << min(AqHead-headctrl[i][m],0.0)                         <<endl;	//absolute error
		
		}
	}
  EXTRACT.close();
	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "Head" or "River", string name 
  {double x double y double head}x(numlines+1)
  &[int precision]
----------------------------------------------------------------------*/
CSurfaceDrain *CSurfaceDrain::Parse(ifstream &input,int &l,  CSingleLayerABC *pLay,char * Name){

	CSurfaceDrain  *pRiver=NULL;

  bool     eof(false),done(false);
	int      Len(0),thisprec(3),nlines(0),i;
  cmplex	 stringp[MAXLINES];
  double   headp  [MAXLINES];
	char    *s      [MAXINPUTITEMS];

	if(parserdebug){cout<<"Surface Drain element"<<endl;}
	
	eof=TokenizeLine(input,s,Len); l++;  
  do {
		if (Len==3){
			if (nlines>=MAXLINES) { ExitGracefully("CSurfaceDrain::Parse- too many line segments in surface drain",TOO_MANY);}
      stringp[nlines]=s_to_c(s[0],s[1]);
      headp  [nlines]=s_to_d(s[2]     ); 
			nlines++; 
			eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pRiver = new CSurfaceDrain(Name,
																 pLay,
																 stringp,
																 nlines-1,
																 headp,
																 thisprec);
			done=true;                                                                   }
    else             {ImproperFormat(s,l); return NULL;                            }
	} while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pRiver;}
}
//******************************************************************
void CSurfaceDrain::GetSegMatrixBuildInfo (const int i, MatrixInfo &info){

	//Same as river
	int    m,n;
	double K=pLayer->GetCond (zctrl[i][int(nlinecontrol/2)]);
	double B=pLayer->GetBase (zctrl[i][int(nlinecontrol/2)]); 
  double T=pLayer->GetThick(zctrl[i][int(nlinecontrol/2)]);

	info.phiCoeff=-1.0;
	info.QxCoeff =0.0;
	info.QyCoeff =0.0;
	info.nctrl   =nlinecontrol;

	SetConstraints(i,0,0,0,1,ZERO_AT_END,0.0,0.0);
	BuildUnitMatrix(i,LINESINK,false);
	//BuildSimpleUnitMatrix(i,LINESINK,false);

	cout << name << " Sending Explicit Matrix info"<<endl;

  //set right hand side for system of equations- 
  for (m=0; m<nlinecontrol; m++){
		//copy unit matrix 
		for (n=0;n<=order;n++){
			info.unit[n][m]=unit[m][n]; //note reversal of direction
		}
		//copy control points
		info.zctrl[m]=zctrl[i][m];

    //calculate element donation to RHS
    info.elemrhs[m]=ConvertToPotential(headctrl[i][m]-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());
	}

}
//#################################################################################################################

