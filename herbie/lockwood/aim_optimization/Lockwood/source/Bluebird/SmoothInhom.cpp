//SmoothInhom.cpp

#include "SmoothInhom.h"
#include "TriagonalMesh.h"
#include "MatrixInclude.h"
//********************************************************************	
//Multiquadric constructor
CSmoothInhom::CSmoothInhom (char						*Name, 
						  							const CSingleLayerABC *pLay,						
														const cmplex		*points,
														const int				 NumOfLines,
														const cmplex    *condcontrol,
														const double		*cond, 
														const int        numctrl,
														const int			   prec)
							:CAreaVortex(Name,pLay,points,NumOfLines,prec){

	int m,n;

	pCondZone=new CMQPolyPropZone(hydraulic_conductivity,points,NumOfLines,condcontrol,cond,numctrl,0.0,1.0,0.0); //TMP DEBUG
	
	
	ncurlctrl=50; //TMP DEBUG
	cout <<"ncurlctrl in: "<<ncurlctrl<<endl;	
	CTriMesh::CreatePolygonControlPoints(points,NumOfLines,zcurlctrl,ncurlctrl);	
	cout <<"ncurlctrl out: "<<ncurlctrl<<endl;
	cin>>m;
	if (ncurlctrl==0){ExitGracefully("CSmoothInhom:Constructor: unable to generate control points", BAD_DATA);}

	curlctrl =new double [ncurlctrl];
	ave_curl=0.0;
	for (m=0; m<ncurlctrl; m++){curlctrl [m]=0.0;}

  //ExitGracefully("",BAD_DATA);
	MQorder=20; //TMP DEBUG

	cout <<"MQorder in:  "<< MQorder <<endl;	
	CTriMesh::CreatePolygonControlPoints(points,NumOfLines,zMQbasis,MQorder);	
	cout <<"MQorder out: "<< MQorder <<endl;
  MQorder=min(MQorder,MAX_MQ_ORDER);

	if (MQorder>ncurlctrl){ExitGracefully("CSmoothInhom:Constuctor: unable to overspecify",BAD_DATA);}
	if (MQorder==0){ExitGracefully("CSmoothInhom:Constructor: unable to generate basis points", BAD_DATA);}
	
	MQcoeff =new double [MQorder];
	for (n=0; n<MQorder; n++){MQcoeff[n]=0.0;}

}
//********************************************************************	
//Linear constructor
CSmoothInhom::CSmoothInhom (char									*Name, 
						  							const CSingleLayerABC *pLay,						
														const cmplex					*points,
														const int							 NumOfLines,
														const cmplex					 zcond,
														const double					 cond,
														const cmplex					 kslope,
														const int							 prec)
						 :CAreaVortex(Name,pLay,points,NumOfLines,prec){

	int m,n;

	pCondZone=new CLinPolyPropZone(hydraulic_conductivity,points,NumOfLines,zcond,cond,kslope);
	
	ncurlctrl=15; //TMP DEBUG
	cout <<"ncurlctrl in: "<<ncurlctrl<<endl;	
	CTriMesh::CreatePolygonControlPoints(points,NumOfLines,zcurlctrl,ncurlctrl);	
	cout <<"ncurlctrl out: "<<ncurlctrl<<endl;

	curlctrl =new double [ncurlctrl];
	ave_curl=0.0;
	for (m=0; m<ncurlctrl; m++){
		curlctrl [m]=0.0;
	}

//ExitGracefully("",BAD_DATA);
	MQorder=30; //TMP DEBUG

	cout <<"MQorder in: "<<MQorder<<endl;	
	CTriMesh::CreatePolygonControlPoints(points,NumOfLines,zMQbasis,MQorder);	
	cout <<"MQorder out: "<<MQorder<<endl;
  MQorder=30;//min(MQorder,MAX_MQ_ORDER);

	if (MQorder>ncurlctrl){ExitGracefully("CSmoothInhom:Constuctor: unable to overspecify",BAD_DATA);}
	MQcoeff =new double [MQorder];
	for (n=0; n<MQorder; n++){
		MQcoeff[n]=0.0;
	}

}
//********************************************************************
CSmoothInhom::~CSmoothInhom(){	 
}
//********************************************************************
bool      CSmoothInhom::CalcCurl(const double t){  
	//Calculates curl at control points (neccesary?)
	for (int m=0; m<ncurlctrl; m++){
		curlctrl[m]=GetDesiredCurl(zcurlctrl[m],t);
	}
	return false;
}
//********************************************************************
double CSmoothInhom::GetDesiredCurl (const cmplex &z,const double &t) const{

	cmplex kprime,Qxy;
	double kloc;

	kloc  =pCondZone->GetValue(z);
	kprime=pCondZone->GetSlopes(z);
			
	Qxy   =conj(pLayer->GetW(z,t));

	return 1.0/kloc*(Qxy.real()*kprime.imag()-Qxy.imag()*kprime.real());

}
//********************************************************************
void	CSmoothInhom::SolveMQCoeff(const double t){   
	//solves for Multiquadric curl coefficients based upon k @ ctrl pts
	int m;

	//int n,n1,n2,s,rank;
	//double relax;
	//cout <<"SOLVING MQ COEFF"<<endl;
	if (MQorder<=1){ //piecewise constant- no need- no curl 
    ave_curl=0.0; 
    if (MQorder==1){MQcoeff[0]=0.0;}
  }
  else{  		
		cmplex kprime,Qxy;
		double kloc;
		double rs,rn,Xs,Xn,Ym;
		int    s,n;  
		
    // initialize matrices
		double **A;
		double  *b;
		double  *sol;
		A  =new double *[MQorder+1];
		b  =new double  [MQorder+1];
		sol=new double  [MQorder+1];
    for(n=0; n<=MQorder; n++){
			A[n]=new double [MQorder+1];
      for(int n1=0; n1<=MQorder; n1++){ 
		    A[n][n1]=0.0;
			}
      b  [n]=0.0;
			sol[n]=0.0;
		}

		cmplex zcen=Centroid();
    for (m=0; m<ncurlctrl; m++){

			//get values of k, k'x, k'y, and Qx/Qy non j
			kloc  =pCondZone->GetValue (zcurlctrl[m]);
		  kprime=pCondZone->GetSlopes(zcurlctrl[m]);

			curlctrl[m]=GetDesiredCurl(zcurlctrl[m],t); //needed?

			Qxy   =conj(pLayer->GetW(zcurlctrl[m],t)-GetInteriorW(zcurlctrl[m],t));

			Ym=(1.0/2.0/kloc*(kprime.real()*(zcurlctrl[m]-zcen).real()+
					              kprime.imag()*(zcurlctrl[m]-zcen).imag()))-1.0;	

			//Assemble system of equations (needs to be revised) 
      for (s=0; s<MQorder; s++){

	     	rs=abs(zMQbasis[s]-zcurlctrl[m]);
				Xs=rs*(1.0/3.0/kloc*(kprime.real()*(zcurlctrl[m]-zMQbasis[s]).real()+
					                   kprime.imag()*(zcurlctrl[m]-zMQbasis[s]).imag())-1.0);

				for (n=0; n<MQorder; n++){
					rn=abs(zMQbasis[n]-zcurlctrl[m]);
					Xn=rn*(1.0/kloc/3.0*(kprime.real()*(zcurlctrl[m]-zMQbasis[n]).real()+
															 kprime.imag()*(zcurlctrl[m]-zMQbasis[n]).imag())-1.0);

					A[s][n]+=(rs*rn);
					//A[s][n]+=Xs*Xn; //TMP DEBUG
				}
				A[s][MQorder]+=rs;
				//A[s][MQorder]+=Xs*Ym;;
				b[s]         +=rs*curlctrl[m];
				//b[s]         +=Xs*(-1.0/kloc*(Qxy.imag()*kprime.real()-Qxy.real()*kprime.imag()));
      }
    }

		//Last row
		//Conventional-summation of coefficients equal to one---------------------
		
		for (n=0; n<MQorder; n++){
       A[MQorder][n]=1.0;
		}
    A[MQorder][MQorder]=0.0;
    b[MQorder]         =0.0;
    
    
		//Alternative: Insures net curl along boundary to be zero------------------------------
		cmplex *ze;
		int     Nlines;
		double  length=0;

		Nlines=pLinevortexBoundary->GetNumSegs();
    ze=pLinevortexBoundary->GetZArray();

    /*
		double r1,r2;		
		for (n=0; n<MQorder; n++){
			//int_perim curl ds=0
			A[MQorder][n]=0.0;
			

      A[MQorder][MQorder]=0.0;

			for (int i=0;i<Nlines;i++){

			  length=abs(ze[i+1]-ze[i]);
			  r1=abs(zMQbasis[n]-1.5773502692/2.0*(ze[i+1]-ze[i])+ze[i]);
				r2=abs(zMQbasis[n]-0.4226497308/2.0*(ze[i+1]-ze[i])+ze[i]);
				A[MQorder][n]+=(r1+r2)*length/2.0;

				A[MQorder][MQorder]+=length;
			}
    }*/
		/*double term1,term2,chi1,chi2,Yn,Yc;
		cmplex Zc,Zn;
		for (n=0; n<MQorder; n++){A[MQorder][n]=0.0;} A[MQorder][MQorder]=0.0;

		for (int i=0; i<Nlines; i++){
			//obtain jump at node points, additive constant
			length=abs(ze[i+1]-ze[i]); 
			Zc=(Centroid()-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
			Yc=Zc.imag();

			//netcurl+=  pAreaVortex->ave_curl*(-length*length*Yc/4.0);
			A[MQorder][MQorder]+=(-length*length*Yc/4.0);
    
			for (n=0; n<MQorder; n++){
				Zn=(zMQbasis[n]-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
				Xn=Zn.real();
				Yn=Zn.imag();
				chi1=-1.0-Xn; //eq 37 in asink paper
				chi2= 1.0-Xn; //eq 43 in asink paper
				if (Yn!=0.0){
					term1=sqrt(chi1*chi1+Yn*Yn); 
					term2=sqrt(chi2*chi2+Yn*Yn);
					//netcurl+=  pAreaVortex->MQcoeff[n]*Yn*pow(0.5*length,3)/6.0*
					//			( (chi1*term1 + Yn*Yn*log(chi1+term1))-
					//				(chi2*term2 + Yn*Yn*log(chi2+term2)));            //eq 40 plus eq 42 : \mu_{j+1}  
					A[MQorder][n]+=Yn*pow(length,3)/48.0*
								        ( (chi1*term1 + Yn*Yn*log(chi1+term1))-
									        (chi2*term2 + Yn*Yn*log(chi2+term2)));
				}
			}
		}
    b[MQorder]         =0.0;
		
*/
   /* cout <<endl<< "A matrix"<<endl;
    for(s=0; s<=MQorder; s++){ 
		cout<<'|';  
    for(n=0; n<=MQorder; n++){ 
			if(fabs(A[s][n])>REALSMALL){cout.width(6); cout<<A[s][n]<<" ,";} 
			else                  {cout.width(7); cout<<0.0    <<" ,";}} cout<<" | "<<b[s]<<endl;}  
    cout << "b vector"<<endl; for(n=0; n<=MQorder; n++){cout <<b[n]<<" ";} 
    cout <<endl;*/
    
		if (!SVD(A,b,sol,(MQorder+1))){
			ExitGracefully("CConductInhom::SolveMQCoeff:SVD routine failed",SINGMAT);}

    //pick up solution
		double relax=1.0;//TMP DEBUG
    ave_curl=relax*(sol[MQorder]-ave_curl)+ave_curl;
		for (n=0; n<MQorder; n++){
			MQcoeff[n]=relax*(sol[n]-MQcoeff[n])+MQcoeff[n];
		  //cout <<MQcoeff[n]<<endl;
		}
		for(n=0; n<=MQorder; n++){delete [] A[n];}
		delete [] A;		
		delete [] b;
		delete [] sol;
	}
	
  //Just normal MQ interpolation---------------------------------------------

	double maxerror(0.0);
	maxerror=MQInterpolate(curlctrl,zcurlctrl,ncurlctrl,MQcoeff,zMQbasis,MQorder,0.0,ave_curl);
	cout <<"ERROR: "<< maxerror<<endl;


}
//********************************************************************
CPropZone *CSmoothInhom::GetCondZone() const{return pCondZone;}
//********************************************************************
double CSmoothInhom::GetMaxError   (const double &t) const{
	double error(0),maxerror(0.0);

	for (int m=0; m<ncurlctrl;m++){
		error=GetDesiredCurl(zcurlctrl[m],t)-GetCurl(zcurlctrl[m],t);//curl=length/time
		upperswap(maxerror,error);
	}
	return maxerror;
	
}
//********************************************************************
void CSmoothInhom::WriteOutput   (const double &t) const{
  //write errors file
	double error(0.0),maxcurl(-ALMOST_INF),mincurl(ALMOST_INF);
  int m;

	ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);

	for (m=0; m<ncurlctrl;m++){
		upperswap(maxcurl,curlctrl[m]);
		lowerswap(mincurl,curlctrl[m]);
	}
	for (m=0; m<ncurlctrl;m++){
		error=GetDesiredCurl(zcurlctrl[m],t)-GetCurl(zcurlctrl[m],t);

		ERRORS << zcurlctrl[m].real()<<","<<zcurlctrl[m].imag()<<","; //coordinate
		ERRORS << FLUX_ERROR_TAG    <<","<<0<<",";
		ERRORS << error<< ",";                                        //absolute error 
		ERRORS << min(error/(maxcurl-mincurl),error) << endl;         //relative error
    //ERRORS <<GetDesiredCurl(zcurlctrl[m],t)<<",";
    //ERRORS <<CAreaVortex::GetCurl(zcurlctrl[m],t) <<endl;

	}
	for(int n=0; n<MQorder;n++){cout << "MQcoeff["<<n<<"]: "<<MQcoeff[n]<<endl;} //TMP DEBUG

	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string SmoothInhom, string name 
	{double x double y}x(numlines+1) 
	& 
	{double x double y double k}x(numcontrolpts)
	&[int precision]
----------------------------------------------------------------------*/
CSmoothInhom *CSmoothInhom::Parse(ifstream &input, int &l, 
														CSingleLayerABC *pLay, char * Name){

	CSmoothInhom *pSmoothInhom=NULL;
  bool					eof(false),done(false);
	int						Len,thisprec,nlines(0),npoints(0),i;
  cmplex				stringp[MAXLINES];       
  cmplex				condzp [MAXLINES];
  double				condp  [MAXLINES];
	char				 *s[MAXINPUTITEMS];
	bool          linear(false);

	if (parserdebug) {cout << "Smooth Inhomogeneity element"<<endl;}  eof=TokenizeLine(input,s,Len); l++; 
  do {
		if (nlines>=MAXLINES) { ExitGracefully("CSmoothInhom::Parse- too many line segments on smooth inhomogeneity boundary",TOO_MANY);}
		if       (Len==0){                             eof=TokenizeLine(input,s,Len); l++;}
	  else if  (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]); nlines++; eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len<=3) && (!strcmp(s[0],"&"))) {
      stringp[nlines]=stringp[0];        nlines++; eof=TokenizeLine(input,s,Len); l++;
			done=true;}
    else {
			ImproperFormat(s,l); return NULL;}
	} while (!done);

	done=false;                                     
	do {
    if (npoints>=MAXLINES) { ExitGracefully("CSmoothInhom::Parse- too many control points in smooth inhomogeneity",TOO_MANY);}
    if (Len==0) {                                  eof=TokenizeLine(input,s,Len); l++;}		
		else if ((Len==3) && (strcmp(s[0],"&"))){
      condzp[npoints]=s_to_c(s[0],s[1]);
			condp [npoints]=s_to_d(s[2]);
			npoints++;                                   eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==5) {
      condzp[0]=s_to_c(s[0],s[1]);
      condzp[1]=s_to_c(s[2],s[3]);
			condp [0]=s_to_d(s[4]);
			linear=true;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);}
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			if (!linear){
			pSmoothInhom = new CSmoothInhom(Name,
				  													  pLay,
																			stringp,
                                      nlines-1,
																			condzp,
																			condp,
																			npoints,
																			thisprec);
			}
			else{
		/*	pSmoothInhom = new CSmoothInhom(Name,
				  													  pLay,
																			stringp,
																			condzp[0],
																			condp,
																			condzp[1],
																			thisprec,
																			nlines-1);*/
			}
			done=true;
		}
		else{ImproperFormat(s,l); return NULL;}
	} while ((!done) & (!eof));

  if (eof) {return NULL;}
	else     {return pSmoothInhom;}
}

