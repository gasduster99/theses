#include "StringElem.h"
#include "MatrixInclude.h"

/************************************************************************
				BUILD UNIT MATRIX 
*************************************************************************
  Builds generic matrix of unit influences for this segment 
  details of implementation in Jankovic [1997] (PhD Dissertation)
-----------------------------------------------------------------------*/
void CLineElem::BuildUnitMatrix(const linetype LType, 
																const bool     multi_val_c){


	

	int    j,k,m,n;                        //counters
	double sum;                            //temporary variable 
  double length(abs(z2-z1));             //segment length
	double *gZn;                           //row of g matrix (for Chebyshev polynomial coeff)
	
	gZn = new double [order+1];

	//initialize unit influences, gZn
  for(m=0; m<nlinecontrol; m++)   {for(n=0; n<=order; n++){unit[m][n]=0.0;}}
  for(n=0; n<=order; n++) {gZn[n]=0.0;}

  //-----------------------------------------------------------------------------
  // evaluate unit influences for all orders with c1 (jump-specified)
  if(c1!=0){
    for(m=0; m<nlinecontrol; m++){
	    if((multi_val_c) && (c[0]!=NULL)){c1=c[0][m];}
			
      for(n=0; n<=order; n++){ 
				unit[m][n]=c1*cos(n*acos(X[m]));
			}
    }
		if (c1==0){c1++;}
  }
  //-----------------------------------------------------------------------------
  // evaluate unit influences for all orders with c2 (gradient jump specified)
  if(c2!=0){
		for(m=0; m<nlinecontrol; m++){
			if((multi_val_c) && (c[1]!=NULL)){c2=c[1][m];}
			for(n=0; n <=order; n++){  
				sum=0.0;
				for(k=n-1; k>=0; k-=2){
					if (k!=0){sum+=2.0*n*cos(k*acos(X[m]));}
					else     {sum+=    n*cos(k*acos(X[m]));} 
					if(LType==DOUBLET){unit[m][n]-=c2*2.0/length*sum;}
					else              {unit[m][n]+=c2*2.0/length*sum;}
				}
			}
		}
		if (c2==0){c2++;}
  }    
  //-----------------------------------------------------------------------------
  // evaluate unit influences for all orders with c3 (gradient specified)
  if (c3!=0){
    double logterm,two_over;
		for(m=0; m<nlinecontrol; m++){
			two_over=  (2.0/((X[m]-1.0)*(X[m]+1.0)));
      logterm=log(fabs((X[m]-1.0)/(X[m]+1.0)));
      if((multi_val_c) && (c[2]!=NULL)){c3=c[2][m];}
      for(n=0; n <=order; n++){
				for (j=0; j<n; j++)    {gZn[j]=g[n][j];}
			  for (j=n-1; j>=1; j-=2){gZn[j]=2.0*(double)(n)*logterm;}
        if  ((n%2)==1)         {gZn[0]=(double)(n)*logterm;}
			  gZn[n]=two_over;
				unit[m][n]-=(c3*Clenshaw(n,X[m],gZn)/(PI*length));
      }
    }
		if (c3==0){c3++;}
  }
  //-----------------------------------------------------------------------------
  // evaluate unit influences for all orders with c4 (value specified)
  if (c4!=0){
    cmplex ZBH((pLayer->GetBlackHole()-0.5*(z1+z2))/(0.5*(z2-z1)));
    double dist1(abs(ZBH+1.0));     
    double dist2(abs(ZBH-1.0));
    
    for(m=0;m<nlinecontrol;m++){	
			if((multi_val_c) && (c[3]!=NULL)){c4=c[3][m];}
			for(n=0; n<=order; n++){	
	
				//fill fZn, identify influence
        f[n][n]=log(fabs((X[m]-1.0)/(X[m]+1.0)));
				if (LType==DOUBLET){unit[m][n]-=(c4*IOVER2PI*Clenshaw(n,X[m],f[n]));}
				else               {unit[m][n]+=(c4*IOVER2PI*Clenshaw(n,X[m],f[n]));}

				//fix for linesinks
				if (LType==LINESINK){
					unit[m][n]+=c4*IOVER2PI*(- log((1.0-X[m])/dist2)
				                           +(log((1.0+X[m])/dist1)*pow(-1.0,n)));
				}
			}
		}
		if (c4==0){c4++;} //this makes c4 junk, but insures evaluation for multi val c4
  }
	/*for (n=0; n<order;n++){
    for(m=0;m<nlinecontrol;m++){
			cout <<n<<"-"<< unit[m][n] <<endl;
		}
	}*/
	delete [] gZn;
}
/************************************************************************
				GENSOLVE
*************************************************************************
  Generic coefficient solver for a 2D line element
  details of implementation in Jankovic [1997] (PhD Dissertation)
-----------------------------------------------------------------------*/
void CLineElem::GenSolve(double *&JumpCoeff, const linetype LType,
												 const bool multi_val_c, const double relax, 
												 double &objective,double &maxchange){

  double  *aold;                      //old values of jump coeff                        
  double	 change,sum;                //objective indicator data 
  int   	 k,m,n,s;										//counters


	//create arrays to store old jump coeff 
	aold= new double [order+1];

  //A,b,unit & sol are static, so that they dont take up excess memory and 
  //yet can still be passed to subfunctions

  // initialize aold
  for(n=0; n<=order; n++) {aold[n]=0.0;}

	//***********************************************************************
  //build matrix of unit influences

	BuildUnitMatrix(LType,multi_val_c);
  
	//***********************************************************************
  //set system of equations

  // initialize matrices
  for(s=0; s<=MAX_ORDER; s++){ 
    for(n=0; n<=MAX_ORDER; n++){A[s][n]=0.0;} b[s]=0.0;sol[s]=0.0;}
  
  //set equations for coefficients (matrix w/o last 5 rows, last 5 cols)
  for(m=0; m<nlinecontrol; m++){
		for(s=0; s<=order; s++){ 
      for(n=0; n<=order; n++){
				A[s][n]+= unit[m][n]*unit[m][s];
      }
      b[s]+= rhs[m]*unit[m][s];
    }
  }

	//this whole row can be removed!!
	for(n=0;n<=order+4;n++){
    A[order+1][n]=0.0;
  }
  A[order+1][order+1]=1.0;
  b[order+1]=0.0;
  
  //set eq. for constraints
  for (s=0;s<=order;s++){
    A[s][order+2]=0.5*pow(-1.0,s);
    A[s][order+3]=0.5;
    A[s][order+4]=0.5*(1.0+pow(-1.0,s));
	}
  for (n=0;n<=order;n++){
    A[order+2][n]=pow(-1.0,n);
    A[order+3][n]=1.0;
    A[order+4][n]=1.0+pow(-1.0,n);
  }
  b[order+2]=ConstrntVal[0];
  b[order+3]=ConstrntVal[1];
  b[order+4]=ConstrntVal[2];

  //turn off unused constraints
  for (k=0;k<3;k++){
    if (!ConstrntOn[k]){
      for (n=0; n<=order+4;n++){
				A[order+k+2][n]=0.0;
      }
      A[order+k+2][order+k+2]=1.0;
      b[order+k+2]=0.0;
    }
  }

  /*cout <<endl<< "A matrix"<<endl;
  for(s=0; s<=order+4; s++){ cout<<'|';  for(n=0; n<=order+4; n++){ 
	if(fabs(A[s][n])>REALSMALL){cout.width(6);cout<<A[s][n]<<" ,";} 
  else{cout.width(7); cout<<0.0<<" ,";}} cout<<" | "<<b[s]<<endl;}  
  cout << "b vector"<<endl; for(n=0; n<=order+4; n++){cout <<b[n]<<" ";} 
  cout <<endl;*/

  //***********************************************************************
  // solve this system
  if (!Gauss(A,b,sol,(order+5))){ExitGracefully("CLineElem::GenSolve: gauss routine failed",SINGMAT);}
	//if (!SVD(A,b,sol,(order+5)))     {ExitGracefully("CLineElem::GenSolve: SVD routine failed"  ,SINGMAT);}

	//***********************************************************************
  //pick up the solutions and compute max change     
	//------------------------------------------------------------------------
  maxchange=objective=0.0;
	order=order;

  for(n=order; n>=0; n--){ 
		aold[n]=JumpCoeff[n];
    JumpCoeff[n]+=relax*(sol[n]-JumpCoeff[n]);
		change=fabs(JumpCoeff[n]-aold[n])/pLayer->GetDeltaPot();
    if(change>maxchange){maxchange=change;}
  }
	//for(n=0; n<=order; n++){cout <<JumpCoeff[n] << " ";}cout <<endl;

  // check the quality of fit
	//------------------------------------------------------------------------
  for(m=0; m<nlinecontrol; m++){
    sum=0;
		if ((LType==LINESINK) || (LType==DIPOLE)){
			for(n=0; n<=order; n++){sum+=JumpCoeff[n]*unit[m][n];}
		}
		else{
			for(n=0; n<=order; n++){sum+=JumpCoeff[n]*unit[m][n];}
		}
    if(fabs(sum-rhs[m])>objective){      
      objective=fabs(sum-rhs[m]);        
		}
  }

	delete [] aold;
}
/************************************************************************
				Set Constants , Constraints
*************************************************************************
  Sets constraints on and type of  line element solution
  c1,c2,c3,c4 defined in Jankovic [1997] (PhD Dissertation)
  ConstrntOn=true if constraint is on
  ConstrntVal[0]=constrained value of jump function at X=-1 (val1 if CType==ASSIGNED)
  ConstrntVal[1]=constrained value of jump function at X= 1 (val2 if CType==ASSIGNED)
  ConstrntVal[2]=constrained net discharge from element     (val1 if CType==NET_DISCHARGE)
-----------------------------------------------------------------------*/
void CLineElem::SetConstraints(const double tmp1,const double tmp2, const double tmp3, 
											         const double tmp4, const constrainttype CType, const double val1, const double val2){

	c1=tmp1;
	c2=tmp2;
	c3=tmp3;
	c4=tmp4;

	if      (CType==UNCONSTRAINED){
    ConstrntOn[0] =ConstrntOn [1]=ConstrntOn [2]=false;
		ConstrntVal[0]=ConstrntVal[1]=ConstrntVal[2]=0.0;
	}
	else if (CType==ZERO_AT_END){
    ConstrntOn [0]=true; 
		ConstrntOn [1]=ConstrntOn [2]=false;
		ConstrntVal[0]=ConstrntVal[1]=ConstrntVal[2]=0.0;
	}
	else if (CType==ALL_SINGULARITIES){
    /*ConstrntOn [2]=false;
		int n;
		if     (StringID==0)           {ConstrntOn [0]=ConstrntOn [1]=false;}
		else if(StringID==(pString->GetNumLines()-1))    {
			                       ConstrntOn [0]=ConstrntOn [1]=true; 
				  									 ConstrntVal[0]=ConstrntVal[1]=0.0;
			for (n=0;n<=order;n++){ConstrntVal[0]+=JumpCoeff[i-1][n];
														 ConstrntVal[1]+=JumpCoeff[0]  [n]*pow(-1.0,n);
			}
		}
		else                    {ConstrntOn [0]=true; ConstrntOn[1]=false; 
														 ConstrntVal[0]=ConstrntVal[1]=0.0;
			for (n=0;n<=order;n++){ConstrntVal[0]+=JumpCoeff[i-1][n];
			}
		}*/
	}
	else if (CType==ASSIGNED){	
		ConstrntOn [0]=true;
		ConstrntOn [1]=true; 
		ConstrntOn [2]=false;
		ConstrntVal[0]=val1;
		ConstrntVal[1]=val2;
	}
	else if (CType==NET_DISCHARGE){	
		ConstrntOn [0]=false;
		ConstrntOn [1]=false; 
		ConstrntOn [2]=true;
		ConstrntVal[0]=ConstrntVal[1]=0;
		ConstrntVal[2]=val1;
	}
}
/***********************************************************************
				Clenshaw Algorithm
***********************************************************************/
double CLineElem::Clenshaw(const int n,const double &x, Unchangeable1DArray pcoeff) const{
  if (n>1){
    double t((2*x*pcoeff[n])+pcoeff[n-1]);
		double u(pcoeff[n]); 
		double v(0.0);  
    for(int i=n-2;i>=1;i--){v=u; u=t; t=(2.0*x*u)-v+pcoeff[i];}
    return (x*t)-u+pcoeff[0];
  }
  else if (n==1){ return (pcoeff[1]*x)+pcoeff[0];}
  else          { return pcoeff[0];}

}
//-----------------------------------------------------------------------
cmplex CLineElem::cClenshaw(const int n,const cmplex &z, Unchangeable1DArray_z pcoeff) const{
  if (n>1){
    cmplex t((2.0*z*pcoeff[n])+pcoeff[n-1]); 
		cmplex u(pcoeff[n]); 
		cmplex v(0.0);  
    for(int i=n-2;i>=1;i--){
			v=u; 
			u=t; 
			t=(2.0*z*u)-v+pcoeff[i];
		}
    return (z*t)-u+pcoeff[0];
  }
  else if (n==1){ return (pcoeff[1]*z)+pcoeff[0];}
  else          { return pcoeff[0];}
}
//-----------------------------------------------------------------------
cmplex CLineElem::cClenshaw2(const int n,const cmplex &z, Unchangeable1DArray pcoeff, const cmplex &nthterm) const{
	//saves a few cmplex additions, but requires casting anyway
	//saves creation &deletion of fzj & repeat calc of logterm
  if (n>1){
    cmplex t((2.0*z*nthterm)+pcoeff[n-1]); 
		cmplex u(nthterm); 
		cmplex v(0.0);  
    for(int i=n-2;i>=1;i--){
			v=u; 
			u=t; 
			t=(2.0*z*u)-v+pcoeff[i];
		}
    return (z*t)-u+pcoeff[0];
  }
  else if (n==1){ return (nthterm*z)+pcoeff[0];}
  else          { return nthterm;}
}
/***********************************************************************
				Set Far Field Coeff
***********************************************************************/
void CLineElem::SetFarFieldCoeff(){
	MatVectMult(CD,JumpCoeff,FarFldCoeff,MAX_FAR_FIELD,order+1);
	FarFldCoeff[0]=0;
}
//**********************************************************************
void CLineElem::SetPrecision(const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("SetPrecision::Improper precision level specified",BAD_DATA);}
	switch(Precision){
		case(0): {order=0; fold=1.0; break;} //Linear function
		case(1): {order=2; fold=2.0; break;}
		case(2): {order=5; fold=1.5; break;}
		case(3): {order=12;fold=1.5; break;}
	  case(4): {order=20;fold=1.5; break;}
	  case(5): {order=30;fold=1.5; break;}
	  case(9): {order=45;fold=1.5; break;}
	}
}
   //TEMPORARY DEBUG (prints matrix to solve- tests gauss algorithm
  ///////////////////////////////////////////////////////////////
  /*cout <<endl<< "A matrix"<<endl;
  for(s=0; s<=order+4; s++){ 
    cout<<'|';  
    for(n=0; n<=order+4; n++){ 
      if(fabs(A[s][n])>REALSMALL){cout.width(6);cout<<A[s][n]<<" ,";} 
      else{cout.width(7); cout<<0.0<<" ,";}} 
    cout<<" | "<<b[s]<<endl; 
  }  
  cout << "b vector"<<endl; 
  for(n=0; n<=order+4; n++){cout <<b[n]<<" ";} 
  cout <<endl;*/
  ////////////////////////////////////////////////////////////


  //test clenshaw recursion
	/*double fferror(0),tmperr,tmperr2,angle;
	cmplex tmpz,tmze;
	for (int m=0; m<200; m++){
		angle=2.0*pi*(double)(m)/200.0;
    tmpz= 1.001*abs(ze[i]-ze[i+1])*cmplex(cos(angle), sin(angle))+(ze[i]+ze[i+1])/2.0;
		tmze=((tmpz-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i])));
    tmperr=fabs(GetSegmentPotential(i,tmpz,0).real()-Outside(tmze,FForder,FarFldCoeff[i]).real());
    if (tmperr>fferror){fferror=tmperr;}
	}
	for (double r=1.01; r<=5; r+=.01){
		tmperr2=0;
		for (int m=0; m<200; m++){
			angle=2.0*pi*(double)(m)/200.0;
			tmpz= r*abs(ze[i]-ze[i+1])*cmplex(cos(angle), sin(angle))+(ze[i]+ze[i+1])/2.0;
			tmze=((tmpz-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i])));
			tmperr=fabs(GetSegmentPotential(i,tmpz,0).real()-Outside(tmze,FForder,FarFldCoeff[i]).real() );
			if (tmperr>tmperr2){tmperr2=tmperr;}
		}
		//cout <<"max_pot "<< r << " "<<tmperr2<<endl;
    if ((tmperr2>fferror) && (tmperr2>REALSMALL)){cout <<"bad_radius= " <<r << " "<< fferror << " "<<tmperr2<<endl;}	
	  else {cout <<"good radius= " <<r << " "<< fferror << " "<<tmperr2<<endl;}
	}
*/
/*
	//identify desired FForder (in SetFarFldCoeff)
	double Km[10], Tm[10],best[10],angle,err,maxerr;
	cmplex tmpz, tmzem[10];

	for (int m=0; m<10;m++){
		angle=2.0*pi*(double)(m)/10.0;
    tmpz= ISFAR*abs(ze[i]-ze[i+1])*cmplex(cos(angle), sin(angle))+(ze[i]+ze[i+1])/2.0;
		tmzem[m]=((tmpz-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i])));
	  Km[m]=pLayer->GetCond(tmpz);
		Tm[m]=pLayer->GetThick(tmpz);
		best[m]=OutsideRe(tmzem[m],MAX_FAR_FIELD,FarFldCoeff[i]).real();
	}
  //identify error
	for (j=min(FForder[i]+3,MAX_FAR_FIELD); j>MIN_FAR_FIELD; j--){
		maxerr=0;
		for (m=0; m<10;m++){
		  err=ConvertToHead(fabs(OutsideRe(tmzem[m],j,FarFldCoeff[i]).real()-best[m]),Km[m],Tm[m]);
			if (err>maxerr) {maxerr=err;}
		}
		if (maxerr>FF_ERROR){
			FForder[i]=min(j+1,MAX_FAR_FIELD); j=0;
		  //cout <<"FForder Found: " <<FForder[i]<<endl;
		}
	}*/