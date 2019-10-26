
#include "AreaSink.h"
#include "TriagonalMesh.h"
#include "MatrixInclude.h"
//************************************************************************
//                           CConductInhom
//************************************************************************
//Constructor
CConductInhom::CConductInhom(char            *Name, 
		  											 CAquicludeABC   *pAq,
														 const CSingleLayerABC *pLay,
			  										 const cmplex    *points,
														 const int        NumOfLines,
														 const double     cond,
														 const int        prec,
								  					 const int        Precision)
							:CAreaSink(Name,pLay,points,NumOfLines,prec){


  int m,n;
	conduct=cond;
	pAquiclude=pAq;
	relax=1.0;

	//int ordfact=int(sqrt((double)(NLines)));
	int ordfact=3*NumOfLines;
	switch(Precision){
		case(0): {MQorder=1*ordfact; nleakctrl=(int)(1.1*MQorder); break;} 
		case(1): {MQorder=2*ordfact; nleakctrl=(int)(1.2*MQorder); break;}
		case(2): {MQorder=3*ordfact; nleakctrl=(int)(1.3*MQorder); break;}
		case(3): {MQorder=4*ordfact; nleakctrl=(int)(1.4*MQorder); break;}
	  case(4): {MQorder=5*ordfact; nleakctrl=(int)(1.5*MQorder); break;}
	  case(5): {MQorder=6*ordfact; nleakctrl=(int)(1.6*MQorder); break;}
	  case(9): {MQorder=7*ordfact; nleakctrl=(int)(2.0*MQorder); break;}
	}
	MQorder=min(MQorder,MAX_MQ_ORDER-2);

  nleakctrl=MQorder; //TMP DEBUG - should be removed

  //allocate memory for dynamic arrays 
  leakctrl=     new double  [nleakctrl];
  zleakctrl=    new cmplex  [nleakctrl];

  zMQbasis=     new cmplex  [MQorder];
  MQcoeff=      new double  [MQorder]; 

  //initialize leakage array and ctrlpt array 
  for(m=0;m<nleakctrl;m++){leakctrl[m]=0.0; zleakctrl[m]=points[0];}
	ave_leak=0;

  for(n=0;n<MQorder;n++)  {MQcoeff [n]=0.0; zMQbasis [n]=points[0];}

	//Get control points, basis points from voronoi
	cout <<"nleakctrl in: "<<nleakctrl;	
	CTriMesh::CreatePolygonControlPoints(points,NumOfLines,zleakctrl,nleakctrl);	
	cout <<", nleakctrl out: "<<nleakctrl<<endl;

	cout <<"MQorder in: "<<MQorder;	
	CTriMesh::CreatePolygonControlPoints(points,NumOfLines,zMQbasis,MQorder);	
	cout <<", MQorder out: "<<MQorder<<endl;
  MQorder=min(MQorder,MAX_MQ_ORDER);

  ofstream RES;

	RES.open("conductctrl.bna");
	for (m=0;m<nleakctrl;m++){
		RES << "\" AreaSink \",  "<<1<<endl;	 
		RES <<zleakctrl[m].real()<< " , " <<zleakctrl[m].imag()<<endl;
	}
	RES.close();
	RES.open("conductbasis.bna");
	for (n=0;n<MQorder;n++){
		RES << "\" Res \",  "<<1<<endl;	 
		RES <<zMQbasis[n].real()<< " , " <<zMQbasis[n].imag()<<endl;
	}
	RES.close();

	pConductZone= new CPolyPropZone(hydraulic_conductance,points,NumOfLines,conduct);

}
//------------------------------------------------------------------------
CConductInhom::~CConductInhom(){
	//all arrays deleted in ~CAreaSink
}
//************************************************************************
CPropZone *CConductInhom::GetConductZone() const{return pConductZone;}
//************************************************************************
bool CConductInhom::CalcLeakage(const double t){
	for (int m=0; m<nleakctrl; m++){
	  leakctrl[m]=(pAquiclude->GetDesiredLeakage(zleakctrl[m],t)-
			          (pAquiclude->GetLeakage       (zleakctrl[m],t)-GetLeakage(zleakctrl[m],t,FROMTOPANDBOTTOM)));	
		leakctrl[m]= pAquiclude->GetDesiredLeakage(zleakctrl[m],t); //TMP DEBUG - no influence of other elements
		//cout<<"leakctrl"<<leakctrl[m]<<endl;
	}
	return false;
}
//************************************************************************
void CConductInhom::SolveMQCoeff(const double t){
	int m,n,n1,n2,s;
	double rs,rn,r0,Xs,Xn,X0,lambda;
	double pota,potb,Ka,Kb,Ha,Hb;
	char ind;

	double **A;
	double  *b;
	double  *sol;

	A  =new double *[MQorder+1];
	b  =new double  [MQorder+1];
	sol=new double  [MQorder+1];
	for (n=0; n<=MQorder; n++){A[n]=new double [MQorder+1];}

	CLayerSetABC *Above=pAquiclude->GetLayerAbove();
	CLayerSetABC *Below=pAquiclude->GetLayerBelow();
	Ha=Above->GetThickness(zleakctrl[0]);
  Hb=Below->GetThickness(zleakctrl[0]);			
	Ka=Above->GetCond     (c2Dto3D(zleakctrl[0]));//TMP DEBUG - transfer to MultiLayer
	Kb=Below->GetCond     (c2Dto3D(zleakctrl[0]));//TMP DEBUG - transfer to MultiLayer
	if (MQorder<=1){
    ave_leak=0.0;
    for(m=0; m<nleakctrl; m++){ave_leak+=leakctrl[m]/double(nleakctrl);} //this is actually not the proper formulation
    if (MQorder==1){MQcoeff[0]=0;}
  }
  else{                                      
    // initialize matrices
    for(n1=0; n1<=MAX_MQ_ORDER; n1++){ 
      for(n2=0; n2<=MAX_MQ_ORDER; n2++){A[n1][n2]=0.0;}
      b[n1]=0.0;sol[n1]=0.0;
		}
    for (m=0; m<nleakctrl; m++){
			//TMP DEBUG - transfer to MultiLayer
			pota=0.0;//Above->GetDischargePotential(zleakctrl[m],t).real();
			potb=0.0;//Below->GetDischargePotential(zleakctrl[m],t).real();
			if (IsConfined(pota,Ka,Ha)) {
				if (IsConfined(potb,Kb,Hb))  {lambda=-(1.0/(Ka*Ha)          +1.0/(Kb*Hb));ind='a';}
				else                         {lambda=-(1.0/(Ka*Ha));                      ind='b';}}
			else{
				if (pota>1){
					if (IsConfined(potb,Kb,Hb)){lambda=-(1.0/sqrt(2.0*Ka*pota)+1.0/(Kb*Hb));ind='c';}
					else                       {lambda=-(1.0/sqrt(2.0*Ka*pota));            ind='d';}
				}
				else{
					if (IsConfined(potb,Kb,Hb)){lambda=-(1.0/(Kb*Hb));                      ind='e';}
					else                       {lambda=0;                                   ind='f';}				
				}
			}
			//cout <<"lambda: " <<lambda<<"("<<ind<<")"<<endl;
			r0=abs(zleakctrl[m]-Centroid());
			X0=(1.0-(conduct*lambda*pow(r0,2)/4.0));
      for (s=0; s<MQorder; s++){
	     	rs=abs(zMQbasis[s]-zleakctrl[m]);
				Xs=rs-conduct*lambda*pow(rs,3)/9.0;

				for (n=0; n<MQorder; n++){
					rn=abs(zMQbasis[n]-zleakctrl[m]);
				  Xn=rn-conduct*lambda*pow(rn,3)/9.0;

					//A[s][n]+=(rs*rn);
					A[s][n]+=Xs*Xn;
				}
				//A[s][MQorder]+=rs;
				A[s][MQorder]+=Xs*X0;
				//b[s]         +=rs*leakctrl[m];
				b[s]+=Xs*leakctrl[m];
      }
    }
    for (n=0; n<MQorder; n++){
      A[MQorder][n]=1.0;
    }
    A[MQorder][MQorder]=0;
    b[MQorder]         =0;
    /*cout <<endl<< "A matrix"<<endl;
  for(s=0; s<=MQorder; s++){ 
		cout<<'|';  
		for(n=0; n<=MQorder; n++){ 
			if(fabs(A[s][n])>REALSMALL){cout.width(6); cout<<A[s][n]<<" ,";} 
			else                  {cout.width(7); cout<<0.0    <<" ,";}} cout<<" | "<<b[s]<<endl;}  
  cout << "b vector"<<endl; for(n=0; n<=MQorder; n++){cout <<b[n]<<" ";} 
  cout <<endl;*/
    if (!Gauss(A,b,sol,(MQorder+1))){
			ExitGracefully("CConductInhom::SolveMQCoeff:gauss routine failed",SINGMAT);}

    //pick up solution
    ave_leak=relax*(sol[MQorder]-ave_leak)+ave_leak;
    for (n=0; n<MQorder; n++){MQcoeff[n]=relax*(sol[n]-MQcoeff[n])+MQcoeff[n];}
  }
	for (n=0; n<MQorder; n++){delete [] A[n];}
	delete [] A;
	delete [] b;
	delete [] sol;
}

//-----------------------------------------------------------------------
void CConductInhom::SolveItself(double &change, double &objective,const double t){
	double tmpchng(0),tmpobj(0),error(0),old_ave;
	change=0.0;
	objective=0.0;
	double *oldCoeff=new double [MQorder];
	int n;
	for (n=0;n<MQorder;n++){
		oldCoeff[n]=MQcoeff[n];
	}
	//relax=1;
	old_ave=ave_leak;
	CAreaSink::SolveItself(tmpchng,tmpobj,t);
	upperswap(change   ,tmpchng);
	upperswap(objective,tmpobj );
	cout << "Avg Leakage"<<ave_leak<<endl;

	//calculate errors in leakage representation
	for (int m=0; m<nleakctrl;m++){
		error=leakctrl[m]-GetLeakage(zleakctrl[m],t,FROMTOPANDBOTTOM);
	  upperswap(objective,error);
	}
	//local iterations
	int local_iter=1;
	while (((ave_leak*old_ave<0.0) && (local_iter<10)) || (local_iter<0)){
		//change=0.0;
		//objective=0.0;
		old_ave=ave_leak;
		relax*=0.9;
		CAreaSink::SolveItself(tmpchng,tmpobj,t);
		cout << "    local iter: " <<local_iter<< " "<<ave_leak*old_ave<< " "<<ave_leak<<endl;
		upperswap(change   ,tmpchng);
		upperswap(objective,tmpobj );
		local_iter++;
	}

	for (n=0;n<MQorder;n++){
		tmpchng=fabs(MQcoeff[n]-oldCoeff[n])/max(fabs(ave_leak),1.0); //relative change
		upperswap(change,tmpchng);
	}
	upperswap(change,fabs(ave_leak-old_ave));

	delete [] oldCoeff;
}
//-----------------------------------------------------------------------
double CConductInhom::GetMaxError   (const double &t) const{
	double leakerr(0),error(0);
	for (int m=0; m<nleakctrl;m+=max((int)(TEST_ERROR_RATIO*nleakctrl),1)){
		leakerr=leakctrl[m]-GetLeakage(zleakctrl[m],t,FROMTOPANDBOTTOM);
	  upperswap(error,fabs(leakerr/conduct)); //absolute head error
	}
	return error;
}
//************************************************************************
void CConductInhom::WriteOutput(const double &t) const{
  //write errors file
  double error(0);

  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);

	for (int m=0; m<nleakctrl;m++){
		error=pAquiclude->GetDesiredLeakage(zleakctrl[m],t)-GetLeakage(zleakctrl[m],t,FROMTOPANDBOTTOM);//TMP DEBUG no influence of other elements
		ERRORS << zleakctrl[m].real()<<","<<zleakctrl[m].imag()<<","; //coordinate
		ERRORS << FLUX_ERROR_TAG    <<","<<0<<",";
		ERRORS << error<< ",";                                     //absolute error 
		ERRORS << error/max(leakctrl[m],1.0) << endl;               //relative error
	}
	
	ERRORS.close();
}

//-----------------------------------------------------------------------
CConductInhom *CConductInhom::Parse(ifstream &input,int &l, 
											              CAquicludeABC *pAq,CSingleLayerABC *pLay,char * Name){
	//string "ConductanceInhom", string name 
	//double conduct		
	//{double x double y}x(numlines+1)
	//&[int precision]

	CConductInhom  *pLeakSink=NULL;
  bool        eof(false),done(false);
  double      thiscond(0);
	int         Len,thisprec,nlines(0),i;
	cmplex      stringp[MAXLINES];
	char       *s[MAXINPUTITEMS];

	if ((pAq==NULL) || (pLay==NULL)){
		ExitGracefully("CConductInhom::Parse: Recieved NULL pointer to layer or Aquiclude",BAD_DATA);}
	
	if (parserdebug) {cout << "Conductance Inhom"<<endl;} eof=TokenizeLine(input,s,Len); l++;
	do{ 
		if      (Len==1) {thiscond=s_to_d(s[0]); done=true; eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                  eof=TokenizeLine(input,s,Len); l++;}
		else             {ImproperFormat(s,l); return NULL;                                    }
	} while ((!done) && (!eof));
	done=false;  
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CConductInhom::Parse- too many lines in inhomogeneity",TOO_MANY);}
		if      (Len==0){                                   eof=TokenizeLine(input,s,Len); l++;}	  
		else if (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]);nlines++;       eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			stringp[nlines]=stringp[0]; nlines++;	 
			for (i=0; i<nlines; i++){pAq->UpdateExtents(stringp[i]);}
			if  (Len==2){thisprec= s_to_i(s[1]);                   }
			else        {thisprec= CAnalyticElem::DefaultPrecision;}
			pLeakSink= new CConductInhom(Name,
																	 pAq,
																   pLay,
																	 stringp,
																	 nlines-1,
																   thiscond,
												           thisprec,
															     thisprec);
			done=true;
		}		
    else            {ImproperFormat(s,l); return NULL;                                    }
	} while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pLeakSink;}
}
