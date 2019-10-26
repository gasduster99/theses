#include "Stage.h"

/************************************************************************
                           CStage
************************************************************************/
CStage::CStage(){}
//------------------------------------------------------------------------
CStage::CStage(char						 *Name, 
							 const CSingleLayerABC *pLay,
							 const cmplex		 *points,
							 const int				NumOfLines,
							 const double		 *heads,
							 const double		 *depths, 
			         const double		  bot_cond, 
							 const double			bot_thick, 
							 const double			bot_width, 
							 const int				prec) 
       :CStringElem(Name,pLay,points,NumOfLines,false,LINESINK,prec){

  int i,m;

	//CStage Initializers
	conductivity=bot_cond;         
	thickness   =bot_thick; 
	width       =bot_width;               
	resistance  =thickness/(conductivity*width); 
  fresh				=true; 

  //allocate memory for dynamic arrays
  c[3]=         new double       [nlinecontrol];   //for Potential coeff
  head =				new double       [NLines+1];
  depth =				new double       [NLines+1];
  headctrl=     new double      *[NLines];         //for const head
  depthctrl=    new double      *[NLines];  
	connected=    new bool        *[NLines];
	flownodes=    new CFlowNode   *[NLines+1];
	outflow_ID=   new int          [NLines+1];
  for(i=0;i<NLines;i++){
    headctrl[i]=     new double  [nlinecontrol];
	  depthctrl[i]=    new double  [nlinecontrol]; 
    connected[i]=    new bool    [nlinecontrol];
  }
  
  //identify streamflow direction based upon 1st & last elevations
	//copy coordinate array, head array,in downstream direction 
  if (heads[0]<heads[NLines]){
		for(i=0;i<=NLines;i++){
			ze   [i]=points[NLines-i];
			head [i]=heads [NLines-i];
			depth[i]=depths[NLines-i];
		}
	}
	else{
		for(i=0;i<=NLines;i++){
			head [i]=heads[i];
			depth[i]=depths[i];
		}
	}
  for(i=0;i<NLines;i++){//nodes are upstream of river segments (outflow ID not assigned)
		flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],DOWNSTREAM,0.0,outflow_ID[i]);
  }
  for(i=1;i<=NLines; i++){ //nodes are downstream of river segments
	  flownodes[i]=CFlowNode::AddFlowNode((CAnalyticElem*)(this),ze[i],UPSTREAM,0.0,outflow_ID[i]);
	}

  for (i=0;i< NLines;i++){
    Interpolate(headctrl [i],nlinecontrol,head [i],head [i+1]);
    Interpolate(depthctrl[i],nlinecontrol,depth[i],depth[i+1]);
		Interpolate(zctrl    [i],nlinecontrol,ze    [i],ze    [i+1]);
	}

  //initialize boundary conditions to connected
  for (i=0; i<NLines; i++){
		for (m=0; m<nlinecontrol; m++){connected[i][m]=true;}
	}	 
}
//-----------------------------------------------------------------------
CStage::~CStage(){
  if (globaldebug){cout <<"   DESTROYING STAGE-SPEC "<<name<<endl;}
  for(int i=0;i<NLines;i++){
    delete [] headctrl[i];
	  delete [] depthctrl[i];
    delete [] connected[i];
  }
  delete [] head;
  delete [] depth;
  delete [] headctrl;       
  delete [] depthctrl;  
	delete [] connected;
}				
//************************************************************************
void CStage::GetDepth(int i, const double inflow, const double t){}	
//************************************************************************
void CStage::SetSegCoeff  (const int i, double *coeff){
	CStringElem::SetSegCoeff(i,coeff);
	fresh=false;
}
/************************************************************************
                           SolveItself	
************************************************************************/
void CStage::SolveItself(double &change, double &objective, double t){
  double          error(0.0),maxobjective(0.0),maxchange(0.0),K,B,T,TotalPot,LocalPot;
	double          QyJump,Aqhead,oldpot;
	double         *PotOther;
  int             i,m;
  double          baseflow_in;

	PotOther=new double [nlinecontrol];

  for (i=0; i<NLines; i++){

	  K=pLayer->GetCond (zctrl[i][int(nlinecontrol/2)]);
		B=pLayer->GetBase (zctrl[i][int(nlinecontrol/2)]);
		T=pLayer->GetThick(zctrl[i][int(nlinecontrol/2)]);

		baseflow_in=flownodes[i]->GetOutflow(t);
    //GetDepth(baseflow_in)  //does nothing for constant stage
		for (m=0; m<nlinecontrol; m++){
 
      //identify current total,local potential, total Qy jump at control points
      TotalPot=pLayer->GetDischargePotential(zctrl[i][m],t).real();	
      LocalPot=        GetSegmentPotential(i,zctrl[i][m],t).real();
   	  PotOther[m]=TotalPot-LocalPot;

      if (fresh)
			{// first guess  -full hydraulic connection
				oldpot=ConvertToPotential(headctrl[i][m]-B,K,T,
				                          pLayer->GetSeaLevel()-B,
																	pLayer->GetSaltwaterSG());
			}  
			else     
			{// other guesses-last calculated potential
				oldpot=TotalPot;
			}                                      
			
			fresh=false;
			//if assumed connected
			if (connected[i][m]){
			  if (oldpot<0){oldpot=MIN_LINEARIZED_POTENTIAL;}
				if (!IsConfined(oldpot,K,T)){
					c[3][m]=pow(2*K*oldpot,-0.5);										//slope of linearization
					rhs[m]=(headctrl[i][m]-B)-                      //local layer head of river
		     					pow(oldpot/(2*K),0.5)-                  //y-int of linearization  
									PotOther[m]*pow(oldpot*2*K,-0.5);}      //pot from everything*slope		
				else{
				}
			}       
	    //if not connected:
	    else {
	      c[3][m]=0.0;
	      rhs[m]=depthctrl[i][m];
			}
		}//m loop

	  //solve for jump coefficients
		SetConstraints(i,0,resistance,0,1,ZERO_AT_END,0.0,0.0);
    GenSolve(JumpCoeff[i],i,ltype,true,1.0,objective,change);
		upperswap(maxchange   ,change   );
		upperswap(maxobjective,objective);
		//cout << change <<endl;

		//evaluate error, save new BC assumptions
		for (m=0; m<nlinecontrol; m++){
			Aqhead=ConvertToHead((GetSegmentPotential(i,zctrl[i][m],t).real()+PotOther[m]),K,T,pLayer->GetSeaLevel()-B,pLayer->GetSaltwaterSG());
      QyJump=GetWJump(i,X[m],t).imag();      
			if (Aqhead>(headctrl[i][m]-depthctrl[i][m]-B)){
				error=fabs(((headctrl[i][m]-B)-Aqhead)    -QyJump*resistance);
				connected[i][m]=true;
			}
			else {
				error=fabs(depthctrl[i][m]-QyJump*resistance);
				connected[i][m]=false;
			}
			//upperswap(objective,error);
		}
		//upperswap(maxobjective,objective);
    

    SetFarFieldCoeff(i);
		if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){
			pBlocks[i]->Update(myBlockIDs[i],i,t);
		}

		//update base flow to next element
		flownodes[i+1]->UpdateFluxes(baseflow_in+GetSegmentDischarge(i,t),outflow_ID[i+1],t);

	}//i loop
  
	delete [] PotOther;

  change=maxchange;
  objective=maxobjective;
}
//-----------------------------------------------------------------------
double CStage::GetMaxError(const double &t) const{
	double error(0.0),AqHead,QyJump;
	for (int i=0; i<NLines; i++){
		double B=pLayer->GetBase(zctrl[i][(int)(nlinecontrol/2)]);
		for (int m=0; m<nlinecontrol; m++){
			AqHead=pLayer->GetHead(zctrl[i][m],t); 
			QyJump=GetWJump(i,X[m],t).imag();   
			if (AqHead>(headctrl[i][m]-depthctrl[i][m])){
				upperswap(error,fabs((headctrl[i][m]-AqHead)-QyJump*resistance)); //absolute head error

				//upperswap(error,(fabs((headctrl[i][m]-AqHead)-QyJump*resistance))/max((headctrl[i][m]-AqHead),1.0)); //relative head error
			}
			else {
				upperswap(error,fabs(depthctrl[i][m]-QyJump*resistance)); //absolute head error
				//upperswap(error,fabs(depthctrl[i][m]-QyJump*resistance)/max(depthctrl[i][m],1.0)); //relative head error
			}
		}
	}
	return error;
}
//******************************************************************
void CStage::WriteOutput(const double &t) const{
//write extraction file
//write errors file
//write cross section file

  double B,K,T;
  double Aqhead,QyJump;
	double error(0);

	ofstream EXTRACT;
	EXTRACT.open("extract.csv", ios::app);
  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	
	for (int i=0; i<NLines; i++){
		B=pLayer->GetBase (zctrl[i][(int)(nlinecontrol/2)]);
		K=pLayer->GetCond (zctrl[i][(int)(nlinecontrol/2)]);
		T=pLayer->GetThick(zctrl[i][(int)(nlinecontrol/2)]);
		for (int m=0; m<nlinecontrol; m++){
		  
			QyJump=GetWJump(i,X[m],t).imag();    
			
			EXTRACT<<zctrl[i][m].real()<<","
					   <<zctrl[i][m].imag()<<","
						 <<QyJump            <<","
						 <<flownodes[i ]->GetOutflow(t)+GetCumExtract(i,X[m],DIR_RIGHT,t)<<endl;

			Aqhead=pLayer->GetHead(zctrl[i][m],t); 
			if (Aqhead>(headctrl[i][m]-depthctrl[i][m])){
				error=fabs((headctrl[i][m]-Aqhead)-QyJump*resistance);
			}
			else {
				error=fabs(depthctrl[i][m]-QyJump*resistance);
			}
		  ERRORS<<zctrl[i][m].real()<<","<<zctrl[i][m].imag()<<",";
			ERRORS << FLUX_ERROR_TAG    <<","<<0<<",";
			ERRORS<<error<<","<<error<<endl; //TMP DEBUG
		
		}
	}
  EXTRACT.close();
	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "Stage", string name  
	double thickness, double width
	double cond
	{double x double y double head double depth}x(numlines+1) 
	& [int precsion]
----------------------------------------------------------------------*/
CStage *CStage::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){
	
	CStage  *pStage=NULL;

  bool     eof(false),done(false);
  double   thisthick(1.0),thiscond(0),thiswidth(1);
	int      Len,thisprec,nlines(0),i;
  cmplex	 stringp[MAXLINES];
  double 	 headp  [MAXLINES];
	double   depthp [MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Stage  element"<<endl;}     eof=TokenizeLine(input,s,Len); l++; 
	do{ 
		if      (Len==2) {
			thisthick=s_to_d(s[0]); 
			thiswidth=s_to_d(s[1]);
			done=true;																				 eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {ImproperFormat(s,l); return NULL;                                     }
	} while ((!done) && (!eof));
	done=false;
	do{ 
		if      (Len==1) {thiscond=s_to_d(s[0]); done=true;  eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {ImproperFormat(s,l); return NULL;                                     }
	}  while ((!done) && (!eof));
	done=false;
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CStage::Parse- too many lines in resistance river",TOO_MANY);}
    if      (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==4) {
      stringp[nlines]=s_to_c(s[0],s[1]); 
			headp  [nlines]=s_to_d(s[2]     );
			depthp [nlines]=s_to_d(s[3]     );nlines++;        eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pStage = new CStage(Name,
													pLay,
													stringp,
													nlines-1,
													headp,
													depthp,
													thiscond,
													thisthick,
													thiswidth,
													thisprec);
			done=true;
		}
    else             {ImproperFormat(s,l); return NULL;                                     }
	} while ((!done) && (!eof));

	if (eof) {return NULL;}
	else     {return pStage;}
}
//******************************************************************
void CStage::WriteXSection(double t) const{
  double K,B,pot,Qjump,generr,lastlength;
  char connect;
  int i,m;
  ofstream STAGE;
  
	STAGE.open("STAGEDEBUG.csv");
  if (STAGE.fail()){cout <<endl<< "DIDN'T OPEN RIVER OUTPUT"<<endl;}
  else             {cout <<endl<< "WRITING STAGEDEBUG.TXT"  <<endl;}

  STAGE << "x y base bottom riv_head aq_head leakage connect? generr" << endl;

  lastlength=0;

  for (i=0;i<NLines; i++){
		B=pLayer->GetBase(zctrl[i][int(nlinecontrol/2)]);
    K=pLayer->GetCond(zctrl[i][int(nlinecontrol/2)]);
    for (m=nlinecontrol-1;m>=0; m--){
      
			pot=pLayer->GetDischargePotential(zctrl[i][m],t).real();
	    
			Qjump=GetWJump(i,X[m],t).imag();
      
			if (sqrt(pot*2/K) >(headctrl[i][m]-depthctrl[i][m])){
	  	  connect='y';
			  generr=resistance*Qjump+pow(pot*2*K,-0.5)*pot-headctrl[i][m]+pow(pot/(2*K),0.5);} 
			else {
				connect='n';
				generr=resistance*Qjump-depthctrl[i][m];
			}

			STAGE << zctrl[i][m].real() << ",";
			STAGE << zctrl[i][m].imag() << ",";
			STAGE << B <<",";
		  STAGE << headctrl[i][m]-depthctrl[i][m] <<",";
	    STAGE << headctrl[i][m] <<",";
	    STAGE << pow(pot*2/K,0.5) <<",";
	    STAGE << Qjump <<",";
			STAGE << connect << ",";
	    STAGE << generr << ",";
	    STAGE << endl;
	  }
    lastlength+=abs(ze[i]-ze[i+1]);
	}
  STAGE.close();
	cout << "...................................done."<<endl;
}
	//--------------------------------------------
	//lasttest=     new double       [NLines];
	//oldtest=      new double       [NLines];
	//relax=        new double       [NLines];
  // divergent=    new bool         [NLines];
	//lastobj=      new double       [NLines];
	//constrain=    new bool         [NLines];
	//--------------------------------------------

	/*	//--------------------------------------------
		lasttest[i]=oldtest[i]=0.0;
		relax[i]=1.0;
		divergent[i]=false;


    lastobj[i]=0.0;
		constrain[i]=false;
		//--------------------------------------------
	}*/
	//-----------------------------------
	//double sum,maxchng,chng,test,t1,t2;
	//cmplex oldcoeff[MAX_ORDER];
	//-----------------------------------
			//check for divergence
			//--------------------------------------------------------
		/*	sum=maxchng=0;
			for (n=0; n<=order; n++){
					sum+=fabs(JumpCoeff[i][n]);
					chng=fabs(oldcoeff[n]-JumpCoeff[i][n]);
					if (chng>maxchng){maxchng=chng;}
			}//i think this (maxchnge/sum) the same as maxchange from gensolve
			test=maxchng/sum;
			t1=test-lasttest[i]; t2=lasttest[i]-oldtest[i];
      if ((t1>t2) && oppsign(t1,t2)){
					relax[i]*=.3;
					cout << "div "<<relax[i];
					for(n=0; n<=order; n++){JumpCoeff[i][n]=oldcoeff[n];}
					divergent[i]=true;
			}
      else if (relax[i]<1) {
					relax[i]*=1.2;
					cout << "up: "<<relax[i];
          oldtest[i] =lasttest[i];
			    lasttest[i]=test;
			}*/
			//relax[i]=1.0;
			//constant overrelaxation of 1.2 and assumption of BC1 works real well for low resistance,
			//and OK for higher resistance, breaks when BC2 should be taking place- how can we get the 
			//element to identify when this is happening, without being thrown off by a nearby source/sink
		//--------------------------------------------------------
		//for(n=0; n<=order; n++){oldcoeff[n]=JumpCoeff[i][n];}
    //--------------------------------------------------------
			//-------------------------------------------------------
				//-----------------------------------------------
	      /*divergent[i]=true;
				if (divergent[i]){  headguess=headctrl[i][m];
			                      oldpot=0.5*K*pow(headguess,2);
			                        divergent[i]=false;}*/
				//-----------------------------------------------
        //-----------------------------------------------
				//empirical- probably not the best methods
				
				//if the assumed potential is real large, estimate using ctrl head 
				//if ((fabs(headctrl[i][m]-headguess))>(2*headctrl[i][m])) 
				//{oldpot=0.5*K*pow(headctrl[i][m],2);cout << "!";}

				//if (resistance<.0001){oldpot=0.5*K*pow(headctrl[i][m],2);}
        //-----------------------------------------------
	  //----------------------------------------------------
      //constrain if objective increases
		//	if (objective>=100){constrain[i]=true;cout <<"bad";}
		//	else {constrain[i]=false;cout<<"good";}
		//	lastobj[i]=objective;
      //--------------------------------------------------------
//***********************************************************************
//***********************************************************************
//***********************************************************************
//***********************************************************************

//***********************************************************************
//                  VARIABLE_STAGE_CLASS
//***********************************************************************
CVariableStage::CVariableStage(){}
//-----------------------------------------------------------------------
CVariableStage::CVariableStage(char			 *Name, 
														   CSingleLayerABC *pLay, 
							 								 cmplex		 *points, 
															 int				NumOfLines,
							 						     double		 *base_elevs,
							 						     double		 *zeros,
							 						     double		 *runoff, 
		           						     double			bot_cond, 
							 						     double			bot_thick, 
							 						     double			bot_width, 
							 						     double			mannings_n, 
							 						     int				prec):
		            CStage(Name,pLay,points,NumOfLines,base_elevs,zeros,bot_cond,bot_thick,
									     bot_width,prec){
  int i;
 
	//depth should be sent as an array of zeros

  //CVariableStage Initializers      
	roughness=mannings_n;   

  //allocate memory for dynamic arrays
  pRunoff = new double  [NLines];
	pSlope  = new double  [NLines];
	pBase   = new double  [NLines+1];
	basectrl= new double *[NLines+1];
  for(i=0;i<NLines;i++){
    basectrl[i]=     new double  [nlinecontrol];
	}

	//identify flow direction, copy runoff appropriately
  if (base_elevs[0]<base_elevs[NLines]){
		for(i=0;i<NLines;i++){pRunoff[i]=runoff[NLines-i];}
	}
	else{
		for(i=0;i<NLines;i++){pRunoff[i]=runoff[i];}
	}

	//create slope array (head,ze already in correct order)
  for(i=0;i<NLines;i++){
		pSlope [i]=(head[i]-head[i+1])/abs(ze[i+1]-ze[i]);
	  Interpolate(basectrl[i],nlinecontrol,head [i],head [i+1]);
	}
	//copy Base array
  for(i=0;i<=NLines;i++){
  	pBase[i]=head[i];
	}

}
//***************************************************************************
void CVariableStage::GetDepth(int i, const double inflow, const double t){
  double denom;
  double *flow;
	int m;

	flow= new double [nlinecontrol];

	denom=(1/roughness)*pow(2,.666666)*pow(pSlope[i],0.5)*pow(width,1.666666);
	flow[nlinecontrol-1]=inflow;
	depthctrl[i][nlinecontrol-1]=pow(inflow/denom,1.66666);

  for (m=nlinecontrol-2; m>=0; m--){
	  flow[m]=flow[m+1]+ 
						fabs(X[m]-X[m+1])*
						(GetWJump(i,(X[m]+X[m+1])/2.0,t).imag()+
						pRunoff[i]*fabs(X[m]-X[m+1]));
		depthctrl[i][m]=pow(flow[m]/denom,1.666666);
		headctrl [i][m]=basectrl[i][m]+depthctrl[i][m];
	}	
   

	delete [] flow; 
}
//***************************************************************************
CVariableStage *CVariableStage::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){
	//string "VarStage", string name  
	//double cond
	//double thickness
	//double width
	//double roughness
	//{double x double y double base double depth double runoff}x(numlines+1) 
	//& [int precsion]
	CVariableStage  *pStage=NULL;

  bool     eof(false),done(false);
  double   thisthick(1.0),thiscond(0),thiswidth(1),thisrough(0.02);
	int      Len,thisprec,nlines(0),i;
  cmplex	 stringp[MAXLINES];
	double   depthp [MAXLINES];
  double 	 headp  [MAXLINES];
	double   runoffp[MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug) {
						  cout << "VarStage  element"<<endl;}        eof=TokenizeLine(input,s,Len); l++; 
	do{ 
		if      (Len==1) {thiscond=s_to_d(s[0]); done=true;  eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line"<< l << "is wrong length"<<endl; break;}
	}  while ((!done) && (!eof));
	done=false;
	do{ 
		if      (Len==1) {thisthick=s_to_d(s[0]); done=true; eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	do{ 
		if      (Len==1) {thiswidth=s_to_d(s[0]); done=true; eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
  done=false;
	do{ 
		if      (Len==1) {thisrough=s_to_d(s[0]); done=true; eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
  done=false;  
	do {
    if (nlines>=MAXLINES) { ExitGracefully("CVarStage::Parse- too many lines in resistance river",TOO_MANY);}
		if (Len==5) {
      stringp[nlines]=s_to_c(s[0],s[1]); 
			headp  [nlines]=s_to_d(s[2]     );
			depthp [nlines]=s_to_d(s[3]     );
			runoffp[nlines]=s_to_d(s[4]     );nlines++;        eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
        pStage = new CVariableStage(Name,
																		pLay,
														        stringp,
																		nlines-1,
																		headp,
																		depthp,
																		runoffp,
																		thiscond,
																		thisthick,
																		thiswidth,
																		thisrough,
																		thisprec);
			  done=true;
			}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));

	if (eof) {return NULL;}
	else     {return pStage;}
}
//***************************************************************************
/*void CVariableStage::SetTheStage(int i, int L, double t){
  double inflow=0;
	double denom,slope;
	int m;
  //get upstream conditions
  //if ((i==0) && (pUpstream[0]!=NULL)){inflow+=pUpstream[0]->Outflow();}
  //if ((i==0) && (pUpstream[1]!=NULL)){inflow+=pUpstream[1]->Outflow();}
	if (i!=0)                          {inflow=flowctrl[i-1][0];}

  slope=fabs((pBase[i]-pBase[i+1])/abs(ze[i+1]-ze[i]));
	denom=(1/roughness)*pow(2,.66666666)*pow(slope,0.5)*pow(width,1.666666);
	cout << inflow << pBase[i]-pBase[i+1] << " " <<slope << " " << denom << " ";
	flowctrl[i][nlinecontrol-1]=inflow;
	headctrl[i][nlinecontrol-1]=basectrl[i][nlinecontrol-1]+pow(inflow/denom,1.666666);
	
	//finite difference non-uniform flow model (based on I.C. in space)
  for (m=nlinecontrol-2; m>=0; m--){
    flowctrl[i][m]=flowctrl[i][m+1]+                    //inflow
				           fabs(X[m]-X[m+1])*                   //delta x
									 (GetWJump(i,X[m],L,t)+               //leakage at section
										pRunoff[i]);                         //runoff addition at section
		if (flowctrl[i][m]<0) {flowctrl[i][m]=0;}
    // translate into river stage
		headctrl[i][m]=basectrl[i][m]+pow(flowctrl[i][m]/denom,1.666666);
	}

}*/
//***********************************************************************
/*void CStage::SolveItself(double &change, double &objective, double t){
  double          error,maxobjective,maxchange,conv[3],K,TotalPot,QyJump,Aqhead,oldpot;
  double          headguess,LocalPot, PotOther[MAXLINECONTROL];
  cmplex TotalOmega,LocalOmega;
  bool            icon[3],con[3];
  int             i,j,n,m;

	//local iteration variables
  //-----------------------------------------
	int lociter=0;
	bool BC[MAXLINES][MAXLINECONTROL];
	double ExLast[MAXLINES];
	double ExOld[MAXLINES];
	double oldcoeff[MAXLINES][MAX_ORDER];
	double relax[MAXLINES];
	double t1,t2,NetEx;
	bool BCchange=false;
	bool divergent=false;
	bool thisbad=false;
	//-----------------------------------------

  double *pconv=conv; bool *picon=icon; bool *pcon=con;

  double T=pLayer->GetThick();
  int    L=pLayer->GetLevel();
  
  maxobjective=objective=maxchange=error=0;
  objective=100;
  picon[0]=1;   picon[1]=0;   picon[2]=0;
  pconv[0]=0.0; pconv[1]=0.0; pconv[2]=0.0;

	for (i=0; i<NLines; i++){relax[i]=1.0;}

  while ((lociter<50) && objective>.0001){
	  lociter++;objective=0;maxobjective=0;divergent=BCchange=false;

    for (i=0; i<NLines; i++){
	    K=pLayer->GetCond(zctrl[i][int(nlinecontrol/2)]);
      objective=0;
      //set LHS coeff, RHS for system of equations- 

      for (m=0; m<nlinecontrol; m++){
 
        //identify current total,local potential, total Qy jump at control points
        TotalPot=pLayer->GetDischargePotential(zctrl[i][m],L,t).real();	
        LocalPot=GetSegPot(i,zctrl[i][m],L,t).real();

   	    PotOther[m]=TotalPot-LocalPot;

	      oldpot=TotalPot;

        headguess=sqrt(2*oldpot/K);
				//if (resistance<.0001){    headguess=headctrl[i][m];}
	      //if (CAnalyticElem::fresh){headguess=headctrl[i][m];
	          //                      oldpot=0.5*K*pow(headguess,2);}

	      //if hydraulically connected:
	      if (headguess>=basectrl[i][m]){     //assumes unconfined only!
			  	if (oldpot<0){oldpot=MIN_LINEARIZED_POTENTIAL;}
					if (!BC[i][m]){BC[i][m]=true; BCchange=true;}
				  //set m-specific BC coefficients
				  c[1][m]=resistance;
          c[3][m]=pow(2*K*oldpot,-0.5);     //slope of linearization

          //set rhs
          rhs[m]=headctrl[i][m]-                          //head of river
		     				 pow(oldpot/(2*K),0.5)-                   //y-int of linearization  
		             (oldpot-LocalPot)*pow(oldpot*2*K,-0.5);  //pot from everything*slope
				}       
	      //if not connected:
	      else {
					if (BC[i][m]){BC[i][m]=false; BCchange=true;}
          c[1][m]=resistance;
	        c[3][m]=0.0;
	        rhs[m]=headctrl[i][m]-basectrl[i][m];
				}

			}//m loop
      //--------------------------------------------------------
			for(n=0; n<=order; n++){oldcoeff[i][n]=JumpCoeff[i][n];}
      //--------------------------------------------------------

	    //solve for jump coefficients
      GenSolve(JumpCoeff,i,pconv,picon,pcon,false,true,true,true,relax[i],change,objective);

      //--------------------------------------------------------
      divergent=false;
        NetEx=GetSegmentDischarge(i);
        t1=NetEx    -ExLast[i];
			  t2=ExLast[i]- ExOld[i];
				if (t1>t2){cout << "+";}
				if (oppsign(t1,t2)) {cout << "*";}
				if (BCchange){cout<<"B";}
				if ((t1>t2) && (oppsign(t1,t2))){ divergent=true;}
			
        //set relaxation
		   	if (lociter>1){
					if (divergent){ //divergent
						relax[i]*=0.3;
            for(n=0; n<=order; n++){JumpCoeff[i][n]=oldcoeff[i][n];}
						cout << " div! ";
					}
					else if (BCchange){               //change in BCs
						relax[i]*=.8;
					}
					else {relax[i]*=1.2;}
				}
				if (!divergent){
				  ExLast[i]=NetEx;
			    ExOld[i]=ExLast[i];
				}
			
      //--------------------------------------------------------
      for (m=0; m<nlinecontrol; m++){
        Aqhead=pow((GetSegPot(LINESINK,i,zctrl[i][m],L,t).real()+PotOther[m])*2/K,0.5);
        QyJump=GetWJump(LINESINK,i,X[m],L,t);      
			  if (Aqhead>basectrl[i][m]){error=fabs((headctrl[i][m]-Aqhead)        -QyJump*resistance);}
			  else                      {error=fabs((headctrl[i][m]-basectrl[i][m])-QyJump*resistance);}
        if (error>objective)      {objective=error;} 
			}

		  if (objective>maxobjective){maxobjective=objective;}
      if (change>maxchange      ){maxchange=change;}

      //Obtain values for far-field coefficients 
      pFarFldCoeff[i][0]=(0,0);
      for (j=0; j<FForder-1; j++){
        pFarFldCoeff[i][j+1]=(0,0);
        for (n=0; n<=order; n++){
		      pFarFldCoeff[i][j+1]+=
		        cmplex(CD(j,n)*JumpCoeff[i][n]/(2.0*pi),0);}}
		}//i loop
		cout << "loc iter: " <<lociter<< " objective: "<<objective<<endl;
	}//lociter loop
  /////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  change=maxchange;
  objective=maxobjective;
}*/