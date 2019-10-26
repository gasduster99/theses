#include "River.h"
/************************************************************************
                           CRiver
************************************************************************/
//Constructor
CRiver::CRiver(char						 *Name, 
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
  for (i=0; i<NLines; i++){
    Interpolate(headctrl[i],nlinecontrol,head[i],head[i+1]);
		Interpolate(zctrl   [i],nlinecontrol,ze  [i],  ze[i+1]); //if reordering has occured
  }	

	int junk(0);
  for(i=0;i<NLines;i++){//nodes are upstream of river segments (outflow ID not assigned)
		flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],DOWNSTREAM,0.0,junk);
  }
  for(i=1;i<=NLines; i++){ //nodes are downstream of river segments
	  flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],UPSTREAM,0.0,outflow_ID[i]);
	}


}
//-----------------------------------------------------------------------
CRiver::~CRiver(){
	if (globaldebug){cout <<"   DESTROYING RIVER "<<name<<endl;}
	delete [] head;
	delete [] flownodes;
	delete [] outflow_ID;
	for(int i=0;i<NLines;i++){
		delete [] headctrl[i];
	}
	delete [] headctrl;
}
//------STATIC MEMBERS---------------------------------------------------
double CRiver::relax=1.0;
//-----------------------------------------------------------------------
void   CRiver::SetRelaxation(double relaxation){relax=relaxation;}
/***********************************************************************
                           SolveItself
***********************************************************************/
void CRiver::SolveItself(double &change, double &objective,double t){
  double  error(0.0),DesiredPot,maxobjective(0.0),maxchange(0.0),K,B,T;//,AqHead;
	double *PotOther;
  //cmplex  TotalOmega,LocalOmega;
  int     i,m;
	double  baseflow_in;

	PotOther=new double [nlinecontrol];
  
  for (i=0; i<NLines; i++){

		//Solve River Element
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

			//lowerswap(PotOther[m],0.0); //TMP DEBUG- requires smoothing 
			//upperswap(PotOther[m],2.0*K*headctrl[i][m]*headctrl[i][m]);//TMP DEBUG- requires smoothing 

      //calculate the desired potential value	
			DesiredPot=ConvertToPotential(headctrl[i][m]-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());

      //calculate neccesary jump in potential
      rhs[m]=DesiredPot-PotOther[m];

    }

		disabled[i]=false;
    
		//solve for jump coefficients
		SetConstraints(i,0,0,0,1,ZERO_AT_END,0.0,0.0);
    GenSolve(JumpCoeff[i],i,LINESINK,false,relax,objective,change);
		upperswap(maxchange   ,change   );
		upperswap(maxobjective,objective); //relative to range of rhs-useful

		//evaluate error (necc?)
		/*for (m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
			AqHead=ConvertToHead(PotOther[m]+GetSegmentPotential(i,zctrl[i][m],t).real(),K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());
		  error=fabs(AqHead+B-headctrl[i][m])/fabs(headctrl[i][m]-B); //relative head error
			error=fabs(AqHead+B-headctrl[i][m]); //absolute head error
			upperswap(maxobjective,error);
		}*/
		
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
double CRiver::GetMaxError(const double &t) const{
	static double error,AqHead;
	error=0.0;
	for (int i=0; i<NLines; i++){
		double B=pLayer->GetBase(zctrl[i][(int)(nlinecontrol/2)]);
		for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
			AqHead=pLayer->GetHead(zctrl[i][m],t); //Global head
			//upperswap(error,fabs(AqHead-headctrl[i][m])/(headctrl[i][m]-B));//relative head error
			upperswap(error,fabs(AqHead-headctrl[i][m])); //absolute head error
		}
	}
	return error; //Absolute head error
}	

/***********************************************************************
              GetSegCumulativeFlux
***********************************************************************/
double CRiver::GetSegCumulativeFlux(const int i, const cmplex &z, const double &t) const{
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
/***********************************************************************
              READ/WRITE 
***********************************************************************/
void CRiver::WriteOutput(const double &t) const{
	//write extraction file
  //write errors file
  double B;
  double AqHead,QyJump;
	bool IsLake;
  IsLake=(ze[0]==ze[NLines]);
		
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
			//EXTRACT <<flownodes[i]->GetOutflow(t)+GetCumExtract(i,X[m],DIR_RIGHT,t)<<endl;	//cumulative extraction
			if (!IsLake){
			EXTRACT <<flownodes[i]->GetOutflow(t)+GetCumExtract(i,X[m],DIR_RIGHT,t)<<endl;	//cumulative extraction
			}
			else{
			EXTRACT <<0<<endl;
			}
			
			AqHead=pLayer->GetHead(zctrl[i][m],t);
			ERRORS << zctrl[i][m].real()<<","<<zctrl[i][m].imag()   << ",";	//coordinate
			ERRORS << HEAD_ERROR_TAG    <<","<<0<<",";
			ERRORS << fabs(AqHead-headctrl[i][m])/(headctrl[i][m]-B)<< ",";	//relative error	
			ERRORS << AqHead-headctrl[i][m]                         <<endl;	//absolute error
		
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
CRiver *CRiver::Parse(ifstream &input,int &l,  CSingleLayerABC *pLay,char * Name){

	CRiver  *pRiver=NULL;

  bool     eof(false),done(false);
	int      Len(0),thisprec(3),nlines(0),i;
  cmplex	 stringp[MAXLINES];
  double   headp  [MAXLINES];
	char    *s      [MAXINPUTITEMS];

	if(parserdebug){cout<<"River element"<<endl;}
	
	eof=TokenizeLine(input,s,Len); l++;  
  do {
		if (Len==3){
			if (nlines>=MAXLINES) { ExitGracefully("CRiver::Parse- too many line segments in river",TOO_MANY);}
      stringp[nlines]=s_to_c(s[0],s[1]);
      headp  [nlines]=s_to_d(s[2]     ); 
			nlines++; 
			eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pRiver = new CRiver(Name,
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
void CRiver::GetSegMatrixBuildInfo (const int i, MatrixInfo &info){
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


/************************************************************************
                           CRiverSegment
************************************************************************/

CRiverSegment::CRiverSegment(){headctrl=NULL;}
//-----------------------------------------------------------------------
CRiverSegment::CRiverSegment(char						 *Name,   
								             const CSingleLayerABC *pLay,   
														 const cmplex		  end1,   
														 const cmplex		  end2,   
														 const double		  h1, 
														 const double     h2,   
														 const cmplex     branch,
														 const int				prec):
               CLineElem(Name,pLay,end1,end2,LINESINK,prec){

	head1=h1;
	head2=h2;
	zbranch=branch;
	int junk(0);

	headctrl=new double[nlinecontrol];
	
	//node is upstream of river segments (outflow ID not assigned)
	flownode1=CFlowNode::AddFlowNode((CAnalyticElem*)(this),z1,DOWNSTREAM,0.0,junk);
  
  //node is downstream of river segments
	flownode2=CFlowNode::AddFlowNode((CAnalyticElem*)(this),z2,UPSTREAM,0.0,outflowID);

  Interpolate(headctrl,nlinecontrol,head1,head2);
	Interpolate(zctrl   ,nlinecontrol,z1   ,  z2); //if reordering has occured

}
//-----------------------------------------------------------------------
CRiverSegment::~CRiverSegment(){
  delete [] headctrl;
}
//-----------------------------------------------------------------------
double CRiverSegment::relax=1.0;
//-----------------------------------------------------------------------
void CRiverSegment::SetRelaxation(double relaxation){relax=relaxation;}
//-----------------------------------------------------------------------
void CRiverSegment::SolveItself(double &change, double &objective,const double t){
  double  error(0.0),DesiredPot,maxobjective(0.0),maxchange(0.0),K,B,T,AqHead;
	double *PotOther;
  int     m;
	double  baseflow_in;

	PotOther=new double [nlinecontrol];
  
	//Solve River Element
  //******************************************************************
	K=pLayer->GetCond (zctrl[int(nlinecontrol/2)]);
	B=pLayer->GetBase (zctrl[int(nlinecontrol/2)]); 
	T=pLayer->GetThick(zctrl[int(nlinecontrol/2)]); 			
		
	baseflow_in=flownode1->GetOutflow(t);

	disabled=true;

  //set right hand side for system of equations- 
  for (m=0; m<nlinecontrol; m++){

    //identify potential from other elements & sides at control points
		PotOther[m]=pLayer->GetDischargePotential(zctrl[m],t).real();

    //calculate the desired potential value	
		DesiredPot=ConvertToPotential(headctrl[m]-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());

    //calculate neccesary jump in potential
		//rhs[m]=DesiredPot-PotOther[m];
    rhs[m]=DesiredPot-max(PotOther[m],0.0);//perhaps fix so system never dries out?

  }

	disabled=false;
    
	//solve for jump coefficients
	SetConstraints(0,0,0,1,ZERO_AT_END,0.0,0.0);
  GenSolve(JumpCoeff,LINESINK,false,relax,objective,change);
	upperswap(maxchange   ,change   );
	upperswap(maxobjective,objective);

	//evaluate error (necc?)
	for (m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
		AqHead=ConvertToHead(PotOther[m]+GetDischargePotential(zctrl[m],t).real(),K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());
		error=fabs(AqHead+B-headctrl[m])/fabs(headctrl[m]-B);
		upperswap(maxobjective,error);
	}
		
	SetFarFieldCoeff();

	if ((pBlock!=NULL) && (pBlock->IsOn())){
		pBlock->Update(myBlockID,0,t);
	}
		
	//update baseflow downstream
	flownode2->UpdateFluxes(baseflow_in+GetNetDischarge(t),outflowID,t);

	delete [] PotOther;

  change=maxchange;
  objective=maxobjective;
}
//-----------------------------------------------------------------------
void CRiverSegment::GetMatrixBuildInfo (MatrixInfo &info){
	int    m,n;
	double K=pLayer->GetCond (zctrl[int(nlinecontrol/2)]);
	double B=pLayer->GetBase (zctrl[int(nlinecontrol/2)]); 
  double T=pLayer->GetThick(zctrl[int(nlinecontrol/2)]);

	info.phiCoeff=-1.0;
	info.QxCoeff =0.0;
	info.QyCoeff =0.0;
	info.nctrl   =nlinecontrol;

	SetConstraints(0,0,0,1,ZERO_AT_END,0.0,0.0);
	BuildUnitMatrix(LINESINK,false);
	//BuildSimpleUnitMatrix(i,LINESINK,false);

	cout << name << " Sending Explicit Matrix info"<<endl;

  //set right hand side for system of equations- 
  for (m=0; m<nlinecontrol; m++){
		//copy unit matrix 
		for (n=0;n<=order;n++){
			info.unit[n][m]=unit[m][n]; //note reversal of direction
		}
		//copy control points
		info.zctrl[m]=zctrl[m];

    //calculate element donation to RHS
    info.elemrhs[m]=ConvertToPotential(headctrl[m]-B,K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());
	}
}
//-----------------------------------------------------------------------
cmplex CRiverSegment::GetDischargePotential(const cmplex &z, const double &t) const{
	cmplex omega;
	omega=CLineElem::GetDischargePotential(z,t);

	if (CAnalyticElem::Cosmetic){
		double argz,arg1;
		double discharge=GetNetDischarge(t);//+flownode2->GetBranchFlow(t);
    
		argz=arg((z      -z1)/(z2-z1));
		arg1=arg((zbranch-z1)/(z2-z1));
		
		if      ((arg1>=0.0) && (argz>0.0) && (argz>arg1)) {omega-=IM*discharge;}
		else if ((arg1< 0.0) && (argz<0.0) && (argz<arg1)) {omega+=IM*discharge;}
	}
	return omega;
}
//void CRiverSegment::SetBranchFlow(
//-----------------------------------------------------------------------
double CRiverSegment::GetMaxError(const double &t) const{
	double error(0.0),AqHead;
	double B=pLayer->GetBase(zctrl[(int)(nlinecontrol/2)]);
	for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
		AqHead=pLayer->GetHead(zctrl[m],t);
		upperswap(error,(AqHead-headctrl[m])/(headctrl[m]-B));
	}
	return error;
}
//-----------------------------------------------------------------------
void CRiverSegment::WriteOutput(const double &t) const{
	//write extraction file
  //write errors file
  double B;
  double AqHead,QyJump;


	ofstream EXTRACT;
  ofstream ERRORS;
	EXTRACT.open("extract.csv", ios::app);
	ERRORS.open ("errors.csv" , ios::app);
	//cout <<"OUTPUT (should be zero): "<<flownode1->GetOutflow(t)<<endl;
	//cout << "*****"<<flownode1->GetOutflow(t)<<endl;
	B=pLayer->GetBase(zctrl[(int)(nlinecontrol/2)]);
	for (int m=0; m<nlinecontrol; m++){
		  
		QyJump=GetWJump(X[m],0).imag(); 
		EXTRACT <<zctrl[m].real()<<","<<zctrl[m].imag()            <<",";		//coordinate
		EXTRACT <<QyJump                                           <<",";		//point extraction
		EXTRACT <<flownode1->GetOutflow(t)+GetCumExtract(X[m],DIR_RIGHT,t)<<endl;	//cumulative extraction

			
		AqHead=pLayer->GetHead(zctrl[m],t);
		ERRORS << zctrl[m].real()<<","<<zctrl[m].imag()   << ",";	//coordinate
		ERRORS << HEAD_ERROR_TAG <<","<<0<<",";
		ERRORS << fabs(AqHead-headctrl[m])/(headctrl[m]-B)<< ",";	//relative error	
		ERRORS << AqHead-headctrl[m]                         <<endl;	//absolute error	
	}

  EXTRACT.close();
	ERRORS.close();
}
//-----------------------------------------------------------------------
CRiverSegment **CRiverSegment::Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name, int &NumSegs){
  
	int i;
	CRiverSegment  **pRiverSegs=NULL; 

  NumSegs=0;
  bool     eof(false),done(false);
  int      Len,thisprec,nlines(0);
	cmplex   stringp[MAXLINES];
	double   headp  [MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug){cout << "River"<<endl;}      eof=TokenizeLine(input,s,Len); l++;
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CRiverSegment::Parse- too many lines in inhomogeneity",TOO_MANY);}
		if      (Len==0){                                   eof=TokenizeLine(input,s,Len); l++;}	  
		else if (Len==3){
      stringp[nlines]=s_to_c(s[0],s[1]);       
			headp  [nlines]=s_to_d(s[2]);       nlines++;     eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pRiverSegs=new CRiverSegment *[nlines];
			cout<<"here"<<endl;
			cmplex branch;

			for (i=0; i<nlines-1; i++){
				if (i==0) {branch=stringp[0]-1.0;}
				else      {branch=stringp[i-1];}
				cout<<i<<" "<<stringp[i] << " "<<headp[i] <<" "<<stringp[i+1] << " "<<headp[i+1] <<endl;
				pRiverSegs[i]= new CRiverSegment(Name,
													               pLay,
																				 stringp[i], //should check if downstream
												                 stringp[i+1],
												                 headp[i],
																				 headp[i+1],
																				 branch,
												                 thisprec);
				NumSegs++;
			}			
			//NumSegs=nlines;
			done=true;
		}		
    else            {ImproperFormat(s,l); return NULL;}
	} while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pRiverSegs;}
}

/*		//Solve Associated Divergence Element
    //******************************************************************
		if (CAquifer::Multilayer){
			double trans;
			//set right hand side for headjump coeff
			for (m=0; m<nlinecontrol; m++){
				if ((headctrl[i][m]-B)<T){trans=K*(headctrl[i][m]-B);}
        else                     {trans=K*T;}
				rhs[m]=GetWJump(i,X[m],t).imag()/trans;
			}
	
			//solve for Leakage jump coefficients
			SetConstraints(i,0,1,0,0,ZERO_AT_END);
			//GenSolve(LJumpCoeff,i,LINESINK,false,1.0,true,maxobjective,maxchange);
			if (objective>maxobjective){maxobjective=objective;}

			//solve for Leakage taylor coeff (at bounds of ellipse)
			
		}
		*/