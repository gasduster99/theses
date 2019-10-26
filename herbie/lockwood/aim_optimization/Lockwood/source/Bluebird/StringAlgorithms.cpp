#include "StringElem.h"
#include "MatrixInclude.h"
/************************************************************************
				BUILD UNIT MATRIX 
*************************************************************************
  Builds generic matrix of unit influences for this segment 
  details of implementation in Jankovic [1997] (PhD Dissertation)
-----------------------------------------------------------------------*/
void CStringElem::BuildUnitMatrix(const int      i,
																	const linetype LType, 
																	const bool     multi_val_c){


	

	int    j,k,m,n;                        //counters
	static double sum,length;              //temporary variable,segment length
	static double gZnstat[MAX_ORDER+1];    //row of g matrix (for Chebyshev polynomial coeff)
	double *gZn=gZnstat;                   //so that matrix may be sent as variable
	
	length=abs(ze[i+1]-ze[i]);

	//initialize unit influences, gZn
  for(m=0; m<nlinecontrol; m++)   {for(n=0; n<=order; n++){unit[m][n]=0.0;}}
  for(n=0; n<=order; n++) {gZn[n]=0.0;}

  //-----------------------------------------------------------------------------
  // evaluate unit influences for all orders with c1 (jump-specified)
  if(c1!=0.0){
    for(m=0; m<nlinecontrol; m++){
	    if((multi_val_c) && (c[0]!=NULL)){c1=c[0][m];}
			
      for(n=0; n<=order; n++){ 
				unit[m][n]=c1*cos(n*acos(X[m]));
			}
			if (LType==LINEVORTEX){
				//unit[m][0]+=c1;//???
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
					if((LType==DOUBLET) || (LType==LINEVORTEX)){unit[m][n]-=c2*2.0/length*sum;}
					else                                       {unit[m][n]+=c2*2.0/length*sum;}
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
				if (LType==LINEVORTEX){
					//unit[m][n]-=c3/(PI*length)*(-(1.0/(1.0-X[m])) -((1.0/(1.0+X[m]))*pow(-1.0,n)));	 //??
				}
				else if (LType==LINESINK){
					unit[m][n]+=c3/(PI*length)*(-(1.0/(1.0-X[m])) -((1.0/(1.0+X[m]))*pow(-1.0,n)));	
				}

      }
    }
	if (c3==0){c3++;}
  }
  //-----------------------------------------------------------------------------
  // evaluate unit influences for all orders with c4 (value specified)
  if (c4!=0){
		cmplex ZBH=GlobalToLocal(pLayer->GetBlackHole(),i);
    double dist1(abs(ZBH+1.0));     
    double dist2(abs(ZBH-1.0));
    
    for(m=0;m<nlinecontrol;m++){	
			if((multi_val_c) && (c[3]!=NULL)){c4=c[3][m];}
			for(n=0; n<=order; n++){	
	
				//fill fZn, identify influence
        f[n][n]=log(fabs((X[m]-1.0)/(X[m]+1.0)));
				if ((LType==DOUBLET)|| (LType==LINEVORTEX)){unit[m][n]-=(c4*IOVER2PI*Clenshaw(n,X[m],f[n]));}
				else                                       {unit[m][n]+=(c4*IOVER2PI*Clenshaw(n,X[m],f[n]));}

				//fix for linesinks
				if      (LType==LINESINK  ){
					unit[m][n]+=c4*IOVER2PI*(- log((1.0-X[m])/dist2)
				                           +(log((1.0+X[m])/dist1)*pow(-1.0,n)));
				}
				else if (LType==LINEVORTEX){
					unit[m][n]-=c4*IOVER2PI*(- log((1.0-X[m])/dist2)
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

}
/************************************************************************
				GENSOLVE
*************************************************************************
  Generic coefficient solver for a 2D line element
  details of implementation in Jankovic [1997] (PhD Dissertation)
-----------------------------------------------------------------------*/
void CStringElem::GenSolve(      double   *&JumpCoeff, 
													 const int        i,
													 const linetype   LType,
													 const bool				multi_val_c, 
													 const double			relax, 
												 //const bool       a1constant,
													       double		 &objective,
													       double		 &maxchange){

  static double  aold[MAX_ORDER+1];   //old values of jump coeff                        
  static double	 change,sum;          //objective indicator data
  int      rank;                      //for gauss subroutine   
  int   	 k,m,n,s;										//counters

  //A,b,unit & sol are static, so that they dont take up excess memory and 
  //yet can still be passed to subfunctions

  // initialize aold
  for(n=0; n<=order; n++) {aold[n]=0.0;}//TMP DEBUG-A1REVISE

	//***********************************************************************
  //build matrix of unit influences

	BuildUnitMatrix(i,LType,multi_val_c);
  
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
	//TMP DEBUG-A1REVISE    
  //keep a1 constant-used for universal flux solution
	//may be numerically challenging (SVD should handle)
	/*for(n=0; n<=order; n++){A[1][n]=0.0;}
	A[1][1]=1.0
	b[1]=JumpCoeff[i][1];*/


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
    A[s][order+4]=0.5*(1.0-pow(-1.0,s));
	}
  for (n=0;n<=order;n++){
    A[order+2][n]=pow(-1.0,n);
    A[order+3][n]=1.0;
    A[order+4][n]=1.0-pow(-1.0,n);
  }
	//constraint(i,3)=(1.d0-(-1)**i)
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

  //if (!Gauss(A,b,sol,(order+5),rank)){ExitGracefully("CStringElem::GenSolve: gauss routine failed",SINGMAT);}

	if (!SVD(A,b,sol,(order+5))){ExitGracefully("CStringElem::GenSolve: SVD routine failed"  ,SINGMAT);} //TMP DEBUG-A1REVISE
  rank=0;

	//***********************************************************************
  //pick up the solutions and compute max change     
	//------------------------------------------------------------------------
  maxchange=objective=0.0;
	order=order;

  for(n=order; n>=0; n--){ 
		//TMP DEBUG-A1REVISE if (n!=1){
		aold[n]=JumpCoeff[n];
    JumpCoeff[n]+=relax*(sol[n]-JumpCoeff[n]);
		change=fabs(JumpCoeff[n]-aold[n])/pLayer->GetDeltaPot()/relax;
    if(change>maxchange){maxchange=change;}
  }

	//for(n=0; n<=order; n++){cout <<JumpCoeff[n] << " ";}cout <<endl;//TMP DEBUG

  // check the quality of fit
	//------------------------------------------------------------------------
	double max_rhs(-ALMOST_INF);
	double min_rhs( ALMOST_INF);	

	for(m=0; m<nlinecontrol; m++){
		upperswap(max_rhs,rhs[m]);
		lowerswap(min_rhs,rhs[m]);
  }
	if ((max_rhs-min_rhs)<=0.0){objective=0.0; return;}

  for(m=0; m<nlinecontrol; m++){
    sum=0.0;
		if ((LType==LINESINK) || (LType==DIPOLE)){
			for(n=0; n<=order; n++){sum+=JumpCoeff[n]*unit[m][n];}
		}
		else{
			for(n=0; n<=order; n++){sum+=JumpCoeff[n]*unit[m][n];}
		}
		upperswap(objective,fabs(sum-rhs[m])/(max_rhs-min_rhs));
  }

}
/************************************************************************
				Set Constants , Constraints
*************************************************************************
  Sets constraints on and type of  line element solution
  c1,c2,c3,c4 defined in Jankovic [1997] (PhD Dissertation)
  ConstrntOn=true if constraint is on
  ConstrntVal[0]=constrained value of jump function at X=-1
  ConstrntVal[1]=constrained value of jump function at X= 1
  ConstrntVal[2]=constrained net discharge from element
-----------------------------------------------------------------------*/
void CStringElem::SetConstraints(const int i,
																 const double tmp1,
																 const double tmp2, 
																 const double tmp3, 
											           const double tmp4, 
																 const constrainttype CType, 
																 const double val1, 
																 const double val2){
  int n;

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
    ConstrntOn [2]=false;
		if     (i==0)           {ConstrntOn [0]=ConstrntOn [1]=false;}
		else if(i==NLines-1)    {ConstrntOn [0]=ConstrntOn [1]=true; 
				  									 ConstrntVal[0]=ConstrntVal[1]=0.0;
			for (n=0;n<=order;n++){ConstrntVal[0]+=JumpCoeff[i-1][n];
														 ConstrntVal[1]+=JumpCoeff[0]  [n]*pow(-1.0,n);
			}
		}
		else                    {ConstrntOn [0]=true; ConstrntOn[1]=false; 
														 ConstrntVal[0]=ConstrntVal[1]=0.0;
			for (n=0;n<=order;n++){ConstrntVal[0]+=JumpCoeff[i-1][n];
			}
		}
	}
	else if (CType==ASSIGNED){	
		ConstrntOn [0]=true;
		ConstrntOn [1]=true; 
		ConstrntOn [2]=false;
		ConstrntVal[0]=val1;
		ConstrntVal[1]=val2;
	}
	else if (CType==NET_DISCHARGE){	
		ConstrntOn [0]=true;
		ConstrntOn [1]=true; 
		ConstrntOn [2]=false;
		ConstrntVal[0]=0.0;
		ConstrntVal[1]=-val1;
		ConstrntVal[2]=-val1;
	}
}
/***********************************************************************
				Clenshaw Algorithm
***********************************************************************/
double CStringElem::Clenshaw(const int n,const double &x, Unchangeable1DArray pcoeff) const{
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
cmplex CStringElem::cClenshaw(const int n,const cmplex &z, Unchangeable1DArray_z pcoeff) const{
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
cmplex CStringElem::cClenshaw2(const int n,const cmplex &z, Unchangeable1DArray pcoeff, const cmplex &nthterm) const{
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
void CStringElem::SetFarFieldCoeff(const int i){

	MatVectMult(CD,JumpCoeff[i],FarFldCoeff[i],MAX_FAR_FIELD,order+1);
	FarFldCoeff[i][0]=0;
	
	MatVectMult(DD,JumpCoeff[i],dn[i],order+1,false); 
  MatVectMult(DD,       dn[i],gn[i],order+1,false);
	
	MatVectMult(B, JumpCoeff[i],bn[i],order+1,false);
	MatVectMult(DD,bn[i],       en[i],order+1,false);
	MatVectMult(DD,en[i],       hn[i],order+1,false);

}
//**********************************************************************
void CStringElem::SetPrecision(const int Precision,int &order, double &fold){
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