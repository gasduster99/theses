#include "Inhom.h"

/************************************************************************
                           CInhom
************************************************************************/
//Constructor
CInhom::CInhom(){}
//------------------------------------------------------------------------
CInhom::CInhom(char						 *Name, 
							 const CSingleLayerABC *pLay,
							 const cmplex		 *points,
							 const int				NumOfLines,
	             const double		  cond, 
							 const int				prec)
       :CStringElem(Name,pLay,points,NumOfLines,true,DOUBLET,prec){
 
  Kin=cond;  
  pCondZone= new CPolyPropZone(hydraulic_conductivity,points,NumOfLines,cond);

	if (!DynArrayAppend((void**&)(AllInhoms),(void*)(this),NumInhoms)){
		ExitGracefully("CInhom::Constructor: creating NULL element",BAD_DATA);};
}
//***********************************************************************
double     CInhom::GetK() const{return Kin;}
//------------------------------------------------------------------------
CPropZone *CInhom::GetCondZone() const{return pCondZone;}
//***********************************************************************
CInhom   **CInhom::AllInhoms = NULL;
int        CInhom::NumInhoms = 0; 
//------------------------------------------------------------------------
void       CInhom::Destroy(){
	delete [] AllInhoms; //Doesnt actually delete elements, just pointers 
}
/************************************************************************
                           SolveItself
***********************************************************************/
void CInhom::SolveItself(double &change,double &objective,double t){
  double maxchange(0.0),maxobjective(0.0),PotOther,Kout;
	//double error(0.0);
  int    i,m;
  double relax(1.0);

	for (i=0; i<NLines; i++){
		if (!disabled[i]){
		//Solve Inhomogeneity Element
    //******************************************************************

    //identify conductivity outside of inhomogeneity (from center of segment)
		Kout=pLayer->GetCond(zctrl[i][int(nlinecontrol/2)]-(ze[i+1]-ze[i])*IM*MOVE_DIST);

		if (Kout!=Kin){

			disabled[i]=true;//so that potential is from everything but this segment

			//set right hand side for system of equations 
			for (m=0; m<nlinecontrol; m++){
				PotOther=pLayer->GetDischargePotential(zctrl[i][m],t).real();
				//rhs[m]=2.0*(Kin-Kout)*(PotOther)/(Kout+Kin); //old BC - save this
				rhs[m]=2.0*(Kin-Kout)*max(PotOther,0.0)/(Kout+Kin); 
				//cout <<"rhs " << m<<" "<<rhs[m]<<endl;
			}

			/*//Alternative Boundary conditions
			if (Kin>Kout){
				for (m=0; m<nlinecontrol; m++){
					PotOther=pLayer->GetDischargePotential(zctrl[i][m],t).real();
					rhs[m]=(PotOther+0.5*GetOmJump(i,X[m],t).real())*(1.0-(Kout/Kin));  //this is working
				}
			}
			else {
				for (m=0; m<nlinecontrol; m++){
					PotOther=pLayer->GetDischargePotential(zctrl[i][m],t).real();
					rhs[m]=(PotOther-0.5*GetOmJump(i,X[m],t).real())*((Kin/Kout)-1.0); //slower
				}
			}*/

			disabled[i]=false;

			//solve for jump coefficients  
			SetConstraints(i,1,0,0,0,UNCONSTRAINED,0.0,0.0);

			GenSolve(JumpCoeff[i],i,ltype,false,relax,objective,change);
			upperswap(maxchange   ,change   );
			upperswap(maxobjective,objective);
    
			SetFarFieldCoeff(i);

			if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){pBlocks[i]->Update(myBlockIDs[i],i,t);}
		}
		}

  }

  change=maxchange;
  objective=maxobjective;
}
//------------------------------------------------------------------------
double CInhom::GetMaxError(const double &t) const{
	double h1,h2,error(0.0);
	for (int i=0;i<NLines;i++){
		for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
			h1=pLayer->GetHead(zctrl[i][m]-(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
			h2=pLayer->GetHead(zctrl[i][m]+(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
			//upperswap(error,fabs(h1-h2)/(0.5*(h1+h2));//relative head error
			upperswap(error,fabs(h1-h2));//absolute head error
		}
	}
	return error;
}
//************************************************************************
//                           READ/WRITE
//************************************************************************
void CInhom::WriteOutput(const double &t) const{
  //write errors file
  double h1,h2,error(0);
	int i;
  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	
	for (i=0; i<NLines; i++){
		double B=pLayer->GetBase(zctrl[i][(int)(nlinecontrol/2)]);
		for (int m=0; m<nlinecontrol; m++){
			h1=pLayer->GetHead(zctrl[i][m]-(ze[i+1]-ze[i])*IM*MOVE_DIST,t); 
			h2=pLayer->GetHead(zctrl[i][m]+(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
		  ERRORS << zctrl[i][m].real()<<","<<zctrl[i][m].imag()<<","; //coordinate
		  ERRORS << HEAD_ERROR_TAG    <<","<<0<<",";
			ERRORS << h1-h2                                      <<","; //absolute error
			ERRORS <<(h1-h2)/max(h1-B,1.0)                       <<endl;//relative error
		}
	}
	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "Inhom", string name 
	double K 
	{double x double y}x(numlines+1)
	&[int precision]
----------------------------------------------------------------------*/
CInhom *CInhom::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){


	CInhom  *pInhom=NULL;

  bool     done(false);
  double   thiscond(0);
	int      Len,thisprec,nlines(0),i;
	cmplex   stringp[MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug){cout << "Inhomogeneity"<<endl;}      
	
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if      (Len==1) {thiscond=s_to_d(s[0]);            }
  else             {ImproperFormat(s,l); return NULL; }
	
	done=false; 
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CInhom::Parse- too many lines in inhomogeneity",TOO_MANY);}
		if ((Len==2) && (strcmp(s[0],"&"))){
      stringp[nlines]=s_to_c(s[0],s[1]);            nlines++;      
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			stringp[nlines]=stringp[0]; nlines++;
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pInhom= new CInhom(Name,
												 pLay,
												 stringp,
												 nlines-1,
												 thiscond,
												 thisprec);
			done=true;}		

    else            {ImproperFormat(s,l); return NULL;}
	} while (!done);

  return pInhom;
}
//******************************************************************
void CInhom::GetSegMatrixBuildInfo (const int i, MatrixInfo &info){
	int    m,n;

	double Kout=pLayer->GetCond(zctrl[i][int(nlinecontrol/2)]-(ze[i+1]-ze[i])*IM*MOVE_DIST);
	
	info.phiCoeff=2.0*(Kin-Kout)/(Kout+Kin);
	info.QxCoeff =0.0;
	info.QyCoeff =0.0;
	info.nctrl   =nlinecontrol;

	SetConstraints(i,1,0,0,0,UNCONSTRAINED,0.0,0.0);
	BuildUnitMatrix(i,ltype,false);

	cout << name << " Sending Explicit Matrix info"<<endl;
  
  for (m=0; m<nlinecontrol; m++){

		//copy unit matrix (note reversal of direction)
		for (n=0;n<=order;n++){info.unit[n][m]=unit[m][n];}

		//copy control points
		info.zctrl[m]=zctrl[i][m];

    //element donation to RHS
    info.elemrhs[m]=0.0;
	}
}
//************************************************************************
//                          Sift Thru Inhoms
//************************************************************************
void CInhom::SiftThroughInhoms(){
  //Sifts through inhomogeneities and disables one of the shared sides
	int			Points1,Points2;
	int			i,j,k,l;
	double	Cond1,Cond2;
	
  if (NumInhoms<2){return;}
	else {
		cout <<"Sifting through " << NumInhoms <<" Inhomogeneities"<<endl;
		for (i=NumInhoms-1; i>=0; i--){
			Points1=AllInhoms[i]->GetNumSegs()+1;
			Cond1=  AllInhoms[i]->GetK();
			for (j=0; j<i; j++){
				Points2=AllInhoms[j]->GetNumSegs()+1;
				Cond2=  AllInhoms[j]->GetK();		
				for (k=0; k<Points1-1; k++){
				  for (l=0; l<Points2-1; l++){
						if (i==j){cout << "Bad Sift2"<<endl;}
						if (((AllInhoms[i]->GetZ(k)==AllInhoms[j]->GetZ(l))   && (AllInhoms[i]->GetZ(k+1)==AllInhoms[j]->GetZ(l+1))) ||
							  ((AllInhoms[i]->GetZ(k)==AllInhoms[j]->GetZ(l+1)) && (AllInhoms[i]->GetZ(k+1)==AllInhoms[j]->GetZ(l)))){
							if (Cond1>=Cond2) {AllInhoms[i]->SetAsDisabled(k);}
							else						 {AllInhoms[j]->SetAsDisabled(l);}
							if (Cond1==Cond2){AllInhoms[i]->SetAsDisabled(k);
																AllInhoms[j]->SetAsDisabled(l);}
						} 
					}
				}
			}
		}
	}
}
//************************************************************************
//                          Create Box of Inhomogeneities
//************************************************************************
CInhom **CInhom::CreateBoxOfInhoms(CSingleLayerABC *pLay, double xmin, double xmax, double ymin,double ymax, 
																	 int nx, double mincond, double maxcond){
  CInhom **tmp;
	char    *thisname;
  double   thiscond;
	int      j(0);
	cmplex   stringp[5];
	double   xincrem((xmax-xmin)/nx); 
	double   yincrem((ymax-ymin)/nx); 
	tmp = new CInhom *[nx*nx];
  thisname="Cell";
	srand(333);
  pLay->UpdateExtents(cmplex(xmin,ymin));
	pLay->UpdateExtents(cmplex(xmax,ymax));
  
	for   (double xgrid=xmin; xgrid<xmax; xgrid+=xincrem){
		for (double ygrid=ymin; ygrid<ymax; ygrid+=yincrem){

			//assign random conductivity between ranges
			thiscond=log10(mincond)+(log10(maxcond)-log10(mincond))*(double)(rand())/(double)(RAND_MAX);
			thiscond=pow(10,thiscond);

			//create shape
			stringp[0]=cmplex(xgrid,ygrid);
			stringp[1]=cmplex(xgrid+xincrem,ygrid);
			stringp[2]=cmplex(xgrid+xincrem,ygrid+yincrem);
			stringp[3]=cmplex(xgrid,ygrid+yincrem);
			stringp[4]=cmplex(xgrid,ygrid);
			//create element
			tmp[j] =new CInhom(thisname,
												 pLay,
												 stringp,
												 4,
												 thiscond,
												 CAnalyticElem::DefaultPrecision);
			j++;
		}
	}
	return tmp;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

/************************************************************************
                           CInhom
************************************************************************/
//Constructor
CInhomSegment::CInhomSegment(){}
//------------------------------------------------------------------------
CInhomSegment::CInhomSegment(char						 *Name,
								             const CSingleLayerABC  *pLay,
														 const cmplex		  end1,
														 const cmplex		  end2,
														 const double		  cond,
														 const int				prec)
       :CLineElem(Name,pLay,end1,end2,DOUBLET,prec){
 
  //CInhomSegment Initializers
  Kin=cond;  

	//AllInhomSegs[NumInhomSegs]=this;
	//NumInhomSegs++;

	if (!DynArrayAppend((void**&)(AllInhomSegs),(void*)(this),NumInhomSegs)){
		ExitGracefully("CInhomSegment::Constructor: creating NULL element",BAD_DATA);};
}
double CInhomSegment::GetK() const{return Kin;}
//***********************************************************************
CInhomSegment **CInhomSegment::AllInhomSegs = NULL;
int             CInhomSegment::NumInhomSegs = 0; 
//------------------------------------------------------------------------
void       CInhomSegment::Destroy(){
	delete [] AllInhomSegs;
}
/************************************************************************
                           SolveItself
***********************************************************************/
void CInhomSegment::SolveItself(double &change,double &objective,double t){
  double maxchange(0.0),maxobjective(0.0),PotOther,Kout;
  int    m;
  double relax(1.0);

	if (!disabled){
		//Solve Inhomogeneity Segment
    //******************************************************************

    //identify conductivity outside of inhomogeneity (from center of segment)
		Kout=pLayer->GetCond(zctrl[int(nlinecontrol/2)]-(z2-z1)*IM*MOVE_DIST);

		disabled=true;//so that potential is from everything but this segment

    //set right hand side for system of equations 
	  for (m=0; m<nlinecontrol; m++){
      PotOther=pLayer->GetDischargePotential(zctrl[m],t).real();
			//rhs[m]=2.0*(Kin-Kout)*(PotOther)/(Kout+Kin); //old BC - save this
			rhs[m]=2.0*(Kin-Kout)*max(PotOther,0.0)/(Kout+Kin); 
			//cout <<"rhs " << m<<" "<<rhs[m]<<endl;
		}

		/*//Alternative Boundary conditions
		if (Kin>Kout){
      for (m=0; m<nlinecontrol; m++){
        PotOther=pLayer->GetDischargePotential(zctrl[i][m],t).real();
			  rhs[m]=(PotOther+0.5*GetOmJump(i,X[m],t).real())*(1.0-(Kout/Kin));  //this is working
			}
		}
    else {
	    for (m=0; m<nlinecontrol; m++){
				PotOther=pLayer->GetDischargePotential(zctrl[i][m],t).real();
			  rhs[m]=(PotOther-0.5*GetOmJump(i,X[m],t).real())*((Kin/Kout)-1.0); //slower
			}
    }*/

		disabled=false;

		//solve for jump coefficients  
		SetConstraints(1,0,0,0,UNCONSTRAINED,0.0,0.0);

    GenSolve(JumpCoeff,ltype,false,relax,objective,change);
		upperswap(maxchange   ,change   );
		upperswap(maxobjective,objective);
    
    SetFarFieldCoeff();

		if ((pBlock!=NULL) && (pBlock->IsOn())){pBlock->Update(myBlockID,-1,t);}
	
  }

  change=maxchange;
  objective=maxobjective;
}
//------------------------------------------------------------------------
double CInhomSegment::GetMaxError(const double &t) const{
	double h1,h2,error(0.0);
	for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
		h1=pLayer->GetHead(zctrl[m]-(z2-z1)*IM*MOVE_DIST,t);
		h2=pLayer->GetHead(zctrl[m]+(z2-z1)*IM*MOVE_DIST,t);
		upperswap(error,(h1-h2)/max(h1,1.0));
	}
	return error;
}
//************************************************************************
//                           READ/WRITE
//************************************************************************
void CInhomSegment::WriteOutput(const double &t) const{
  //write errors file
  double h1,h2,error(0);
  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	
  double B=pLayer->GetBase(zctrl[(int)(nlinecontrol/2)]);
	for (int m=0; m<nlinecontrol; m++){
		h1=pLayer->GetHead(zctrl[m]-(z2-z1)*IM*MOVE_DIST,t); 
		h2=pLayer->GetHead(zctrl[m]+(z2-z1)*IM*MOVE_DIST,t);
		ERRORS << zctrl[m].real()<<","<<zctrl[m].imag()<<","; //coordinate
		ERRORS << FLUX_ERROR_TAG    <<","<<0<<",";
		ERRORS << h1-h2                                <<","; //absolute error		
		ERRORS <<(h1-h2)/max(h1-B,1.0)                 <<endl; //relative error
	}
	ERRORS.close();
	
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "Inhom", string name 
	double K 
	{double x double y}x(numlines+1)
	&[int precision]
----------------------------------------------------------------------*/
CInhomSegment **CInhomSegment::Parse(ifstream   &input, 
																		 int        &l, 
																		 CSingleLayerABC * pLay, 
																		 char      * Name, 
																		 int        &NumSegs, 
																		 CPropZone *&kzone){

int i;
	CInhomSegment  **pInhomSegs=NULL; 

  NumSegs=0;
	kzone=NULL;
  bool     eof(false),done(false);
  double   thiscond(0);
	int      Len,thisprec,nlines(0);
	cmplex   stringp[MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug){cout << "Inhomogeneity"<<endl;}      eof=TokenizeLine(input,s,Len); l++;
	do{ 
		if      (Len==1) {thiscond=s_to_d(s[0]); done=true; eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                  eof=TokenizeLine(input,s,Len); l++;}
		else             {ImproperFormat(s,l); return NULL;                                    }
	} while ((!done) && (!eof));
	done=false;  
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CInhomSegment::Parse- too many lines in inhomogeneity",TOO_MANY);}
		if      (Len==0){                                   eof=TokenizeLine(input,s,Len); l++;}	  
		else if (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]);nlines++;       eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			stringp[nlines]=stringp[0]; nlines++;
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pInhomSegs=new CInhomSegment *[nlines];
			for (i=0; i<nlines-1; i++){
				pInhomSegs[i]= new CInhomSegment(Name,
													               pLay,
																				 stringp[i], //should check if clockwise
												                 stringp[i+1],
												                 thiscond,
												                 thisprec);
				NumSegs++;
			}

			kzone=new CPolyPropZone(hydraulic_conductivity,
				                      stringp,
															nlines-1,
															thiscond);
			
			done=true;
		}		
    else            {ImproperFormat(s,l); return NULL;}
	} while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pInhomSegs;}
}
//******************************************************************
void CInhomSegment::GetMatrixBuildInfo (MatrixInfo &info){
	int    m,n;

	double Kout=pLayer->GetCond(zctrl[int(nlinecontrol/2)]-(z2-z1)*IM*MOVE_DIST);
	
	info.phiCoeff=2.0*(Kin-Kout)/(Kout+Kin);
	info.QxCoeff =0.0;
	info.QyCoeff =0.0;
	info.nctrl   =nlinecontrol;

	SetConstraints(1,0,0,0,UNCONSTRAINED,0.0,0.0);
	BuildUnitMatrix(ltype,false);

	cout << name << " Sending Explicit Matrix info"<<endl;
  
  for (m=0; m<nlinecontrol; m++){

		//copy unit matrix (note reversal of direction)
		for (n=0;n<=order;n++){info.unit[n][m]=unit[m][n];}

		//copy control points
		info.zctrl[m]=zctrl[m];

    //element donation to RHS
    info.elemrhs[m]=0.0;
	}
}
//************************************************************************
//                          Sift Thru Inhoms
//************************************************************************
void CInhomSegment::SiftThroughInhomSegs(){
  //Sifts through inhomogeneities and disables one of the shared sides
	int			i,j;
	double	Cond1,Cond2;
	
  if (NumInhomSegs<3){return;}
	else {
		cout <<"Sifting through " << NumInhomSegs <<" Inhomogeneity Segments"<<endl;
		for (i=NumInhomSegs-1; i>=0; i--){
			Cond1=  AllInhomSegs[i]->GetK();
			for (j=0; j<i; j++){
				Cond2=  AllInhomSegs[j]->GetK();		
				if (i==j){cout << "Bad Sift2"<<endl;}
				if (((AllInhomSegs[i]->GetZ(0)==AllInhomSegs[j]->GetZ(1))   && 
					   (AllInhomSegs[i]->GetZ(1)==AllInhomSegs[j]->GetZ(0))) ||
						((AllInhomSegs[i]->GetZ(0)==AllInhomSegs[j]->GetZ(1)) && 
						 (AllInhomSegs[i]->GetZ(0)==AllInhomSegs[j]->GetZ(1)))){
					if (Cond1>=Cond2) {AllInhomSegs[i]->SetAsDisabled();}
					else						  {AllInhomSegs[j]->SetAsDisabled();}
					if (Cond1==Cond2) {AllInhomSegs[i]->SetAsDisabled();
														 AllInhomSegs[j]->SetAsDisabled();}
				}
			}
		}
	}
}
//************************************************************************
//                          Create Box of Inhomogeneities
//************************************************************************
CInhomSegment **CInhomSegment::CreateBoxOfInhoms(CSingleLayerABC *pLay, double xmin, double xmax, double ymin,double ymax, 
																	 int nx, double mincond, double maxcond){
  CInhomSegment **tmp;
	char    *thisname;
  double   thiscond;
	int      j(0);
	cmplex   stringp[5];
	double   xincrem((xmax-xmin)/nx); 
	double   yincrem((ymax-ymin)/nx); 
	tmp = new CInhomSegment *[nx*nx];
  thisname="Cell";
	srand(333);
  pLay->UpdateExtents(cmplex(xmin,ymin));
	pLay->UpdateExtents(cmplex(xmax,ymax));
  
	for   (double xgrid=xmin; xgrid<xmax; xgrid+=xincrem){
		for (double ygrid=ymin; ygrid<ymax; ygrid+=yincrem){

			//assign random conductivity between ranges
			thiscond=log10(mincond)+(log10(maxcond)-log10(mincond))*(double)(rand())/(double)(RAND_MAX);
			thiscond=pow(10,thiscond);

			//create shape
			stringp[0]=cmplex(xgrid,ygrid);
			stringp[1]=cmplex(xgrid+xincrem,ygrid);
			stringp[2]=cmplex(xgrid+xincrem,ygrid+yincrem);
			stringp[3]=cmplex(xgrid,ygrid+yincrem);
			stringp[4]=cmplex(xgrid,ygrid);
			//create line elements
			for (int i=0; i<4; i++){
				tmp[j] =new CInhomSegment(thisname,
																	pLay,
																	stringp[i],
																	stringp[i+1],
																	thiscond,
																	CAnalyticElem::DefaultPrecision);
				j++;
			}
		}
	}
	return tmp;
}











































	//relax=1.0-(0.31*(Kin-Kout)/(Kout+Kin)); //good trick sometimes


//******************************************************************
//ANALYTIC SOLUTION METHOD
    //should use analytic method for coeff.(doesn't work for CONTROLEND!=1.0)
		/*double sum=0;
     for (m=0; m<nlinecontrol; m++){
      sum+=rhs[m];
			}
		 cout <<"0 "<<sum/(double)(nlinecontrol) << " "<<JumpCoeff[i][0]<<endl;
		for (n=1; n<=order; n++){
       sum=0;
    for (m=0; m<nlinecontrol; m++){
      sum+=rhs[m]*cos((double)(n)*pi*((double)(m+1)-0.5)/(double)(nlinecontrol));
			}
      cout <<n<<" "<<2.0*sum/(double)(nlinecontrol)<<" "<< JumpCoeff[i][n]<<endl;
		}*/
//******************************************************************

