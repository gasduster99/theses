#include "Aquitard.h"

//************************************************************************
//************************************************************************
//                           CLeakyLayer
//************************************************************************
//************************************************************************

//***********************************************************************
//													CONSTRUCTOR
//***********************************************************************
CLeakyLayer::CLeakyLayer(){order=0;nLeakCtrl=0;}
//-----------------------------------------------------------------------
CLeakyLayer::CLeakyLayer(CAquicludeABC *pAq,
												 int					 ord, 
												 double				 OS){  

  //CLeakyLayer Initializers
	pAquiclude  =pAq;
  zbl         =0;                 
	length      =1.0;                  
	order       =ord;              
	fold        =OS;    
	
	nLeakCtrl=(int)(order*fold);

	//Array Creation and Initialization
	zctrl=    new cmplex *[nLeakCtrl];
  leakctrl=     new double *[nLeakCtrl];
  B=            new double *[order+1];
	for (int i=0; i<nLeakCtrl; i++){
	  zctrl[i]=    new cmplex [nLeakCtrl];
    leakctrl [i]=    new double [nLeakCtrl];
    for (int j=0; j<nLeakCtrl; j++){
      zctrl[i][j]=0.0;
      leakctrl [i][j]=0.0;
		}
	}
  for (int n=0; n<=order; n++){
		B[n]=new double [order+1];
    for (int n2=0; n2<order; n2++){
			B[n][n2]=0.0;
		}
	}
}
//-----------------------------------------------------------------------
CLeakyLayer::~CLeakyLayer(){
	if (globaldebug){cout <<" DESTROYING LEAKY LAYER "<<endl;}
	 for (int i=0; i<nLeakCtrl; i++){
     delete [] zctrl[i];
     delete [] leakctrl[i];
	 }
	 delete [] zctrl;
	 delete [] leakctrl;

	 for (int n=0; n<=order; n++){
		 delete [] B[n];	 
	 }
	 delete []B;
}
//***********************************************************************
//							 STATIC INITIALIZERS
//***********************************************************************
double**			  CLeakyLayer::SinTerm= NULL;
double**			  CLeakyLayer::Buffer = NULL;
int						  CLeakyLayer::order =80;
double					CLeakyLayer::fold = 2.5;
int							CLeakyLayer::nLeakCtrl = 200;

//***********************************************************************
//							 ASSIGNMENT FUNCTIONS
//***********************************************************************
void						CLeakyLayer::SetBounds (window Extents){
	zbl=cmplex(Extents.w-domain_buffer*(Extents.e-Extents.w),
						 Extents.s-domain_buffer*(Extents.n-Extents.s));

  length=max(((1+(2.0*domain_buffer))*(Extents.e-Extents.w)),
						 ((1+(2.0*domain_buffer))*(Extents.n-Extents.s)));

	for (int i=0; i<nLeakCtrl; i++){
    for (int j=0; j<nLeakCtrl; j++){
      zctrl[i][j]=cmplex(zbl.real()+length*(double)(i)/(double)(nLeakCtrl),
                             zbl.real()+length*(double)(j)/(double)(nLeakCtrl));
		}
	}
	Prepare();
	//cout << "SETTTING BOUNDS"<<endl;
}
//***********************************************************************
//							 GetDischargePotential
//***********************************************************************
cmplex CLeakyLayer::GetDischargePotential(const cmplex &z,const double &t) const{
	int    n,n2;
	double pot(0.0);
	cmplex Z((z-zbl)*PI/length);
	double conductance;

	conductance=pAquiclude->GetConductance(z);
	if (conductance!=0){
		for (n=0; n<=order; n++){
			for (n2=0; n2<=order; n2++){
				if (!((n==0) && (n2==0))){
					pot+=(-1/(double)(n*n+n2*n2))*B[n][n2]*sin((double)(n)*Z.real())*sin((double)(n2)*Z.imag());}
			}
		}
		pot*=pow(length/PI,2);
		//pot+=0.25*B[0][0]*( pow(z.real(),2)+pow(z.imag(),2) ); //a possibility
    pot+=0.25*B[0][0]*( pow(z.real()-(zbl.real()+0.5*length),2)+
												pow(z.imag()-(zbl.imag()+0.5*length),2) );
		//return 0; //TEMP DEBUG
	}
	return cmplex(pot,0);
}

//***********************************************************************
//							 GetW
//***********************************************************************
cmplex CLeakyLayer::GetW(const cmplex &z,const double &t) const{
	int    n,n2;
	double Qx(0.0);
	double Qy(0.0);
	cmplex Z((z-zbl)*PI/length);

	double conductance;
	conductance=pAquiclude->GetConductance(z);
	if (conductance!=0.0){
		for (n=0; n<=order; n++){
			for (n2=0; n2<=order; n2++){
				if (!((n==0) && (n2==0))){
				Qx+=((double)(n) /(double)(n*n+n2*n2))*B[n][n2]*cos((double)(n)*Z.real())*sin((double)(n2)*Z.imag());
				Qy+=((double)(n2)/(double)(n*n+n2*n2))*B[n][n2]*sin((double)(n)*Z.real())*cos((double)(n2)*Z.imag());}
			}
		}
		Qx*=length/PI;
		Qy*=length/PI;
		return 0; //TEMP DEBUG
	}
	return cmplex(Qx,Qy);
}
//***********************************************************************
//							 GetLeakage
//***********************************************************************
double CLeakyLayer::GetLeakage (const cmplex &z, const double &t) const{
	int    n,n2;
	double leakage(0.0);
	cmplex Z((z-zbl)*PI/length);

	double conductance;
	conductance=pAquiclude->GetConductance(z);
	if (conductance!=0.0){
		for (n=0; n<=order; n++){
			for (n2=0; n2<=order; n2++){
			 leakage+=B[n][n2]*sin((double)(n)*Z.real())*sin((double)(n2)*Z.imag());
			}
		}
		leakage+=B[0][0]; //add average leakage
		//cout <<leakage<<endl;
	}
	return leakage;
}
//***********************************************************************
//							 GetNetDischarge
//***********************************************************************
double          CLeakyLayer::GetNetDischarge(const double &t) const{return 0.0;}

//***********************************************************************
//							 SolveItself
//***********************************************************************
void            CLeakyLayer::SolveItself(double &change, double &objective,const double t){
  int							i,j,n,n2;
	double					maxchange(0.0),maxobj(0.0),obj,sum;
	double					temp[MAX_LEAK_ORDER][MAX_LEAK_ORDER];
	double					averageleak;
	double					maxL(-ALMOST_INF);
  double				  minL( ALMOST_INF);
	cmplex          Z;
	double conductance;
	conductance=pAquiclude->GetConductance(zbl);

  cout <<"conductance= "<<conductance<<endl; //TEMP DEBUG
  
	
	if (conductance!=0.0){
	//get leakage at control points
	cout <<"getting leakage values"<<endl;
  for (i=0; i<nLeakCtrl; i++){
    for (j=0; j<nLeakCtrl; j++){
			if (false){//TEMP DEBUG
			//get global pressure head from local head & base elevations
			leakctrl[i][j]=(pAquiclude->GetDesiredLeakage(zctrl[i][j],t)-
										  pAquiclude->GetLeakage(zctrl[i][j],t)); 						 		
			}//TEMP DEBUG

			//leakage from line elements above (negative)
			//leakctrl[i][j]-=pLayerAbove->GetLeakage(zctrl[i][j], t);

			//TEMP DEBUG
			//cmplex temppot;
 			//temppot=pLayerAbove->GetElem(1)->GetDischargePotential(zctrl[i][j],0,t);
			//cmplex Z=(zctrl[i][j]-0.5*(cmplex(25.0,75.0)+cmplex(75.0,25.0)))/(0.5*(cmplex(75.0,25.0)-cmplex(25.0,75.0)));
      //leakctrl[i][j]=(temppot*(Z-conj(Z))/(32.0*pi)).real();


			//leakage from line elements below (positive)     
			/*if (pLayerBelow!=NULL) {
        leakctrl[i][j]+=pLayerBelow->GetLeakage(zctrl[i][j],t);
		  }*/
      if (leakctrl[i][j]>maxL) {maxL=leakctrl[i][j];}
      if (leakctrl[i][j]<minL) {minL=leakctrl[i][j];}
		}
	}

	cout <<"solving for coeff"<<endl; //TEMP DEBUG

  //solve for coefficients
	sum=0;
  for (i=0; i<nLeakCtrl; i++){
		for (j=0; j<nLeakCtrl; j++){
			sum+=leakctrl[i][j];
		}
	}
	averageleak=sum/nLeakCtrl/nLeakCtrl;
	//cout <<averageleak<<endl;

  for (n=0; n<=order; n++){
    for (n2=0; n2<=order; n2++){ 
			sum=0.0;
			for (i=0; i<nLeakCtrl; i++){
				for (j=0; j<nLeakCtrl; j++){
					//Z=(zctrl[i][j]-zbl)*pi/length;
				//	sum+=BufferFunction(zctrl[i][j])*(leakctrl[i][j]-averageleak)*sin((double)(n) *Z.real())*
					//																																	sin((double)(n2)*Z.imag());		
					//sum+=BufferFunction(zctrl[i][j])*(leakctrl[i][j]-averageleak)*SinTerm[n][i]*SinTerm[n2][j];	
					sum+=Buffer[i][j]*(leakctrl[i][j]-averageleak)*SinTerm[n][i]*SinTerm[n2][j];	
					
				}
			}
			//cout <<sum <<endl;
      temp[n][n2]=4.0*sum/nLeakCtrl/nLeakCtrl; //TEMP DEBUG
		}
	}
  temp[0][0]=averageleak;

	cout <<"calculating objective"<<endl; //TEMP DEBUG
	//calculate max change 
  for (n=0; n<=order; n++){
    for (n2=0; n2<=order; n2++){ 
			if (fabs(B[n][n2]-temp[n][n2])>maxchange) {maxchange=fabs(B[n][n2]-temp[n][n2]);}
			B[n][n2]=temp[n][n2];
			//cout <<B[n][n2]<< "  "; //TEMP DEBUG
		}
	//	cout <<endl; //TEMP DEBUG
	}

	//calculate objective function
	//only test non-buffered zone
	int start=(int)(((    domain_buffer)*nLeakCtrl)/(1.0+(2.0*domain_buffer)));
  int end  =(int)(((1.0+domain_buffer)*nLeakCtrl)/(1.0+(2.0*domain_buffer)));
obj=0;
  for (i=start; i<=end; i++){
		for (j=start; j<=end; j++){
			//obj=fabs(leakctrl[i][j]-GetLeakage(zctrl[i][j],t));
      //if (obj>maxobj)  {maxobj=obj;}
      obj+=fabs(leakctrl[i][j]-GetLeakage(zctrl[i][j],t));
		}
	}
	maxobj=obj/pow((double)(end-start+1),2);
  change=   maxchange;
	objective=maxobj;
	cout <<"Leaky Layer maxchange: "<< maxchange << endl;
	cout <<"Leaky Layer objective: "<< maxobj << endl;



	//TEMP DEBUG -------------------------------------------------------------
	cout <<"plotting actual leakage" <<endl;
	ofstream DEBUGLEAKY;
	DEBUGLEAKY.open  ("actualleak.grd");
  DEBUGLEAKY << "DSAA"<<endl;
  DEBUGLEAKY <<nLeakCtrl <<" "<<nLeakCtrl <<endl;
  DEBUGLEAKY <<zctrl[0][0].real()<< " "<<zctrl[nLeakCtrl-1][nLeakCtrl-1].real()<<endl;
	DEBUGLEAKY <<zctrl[0][0].imag()<< " "<<zctrl[nLeakCtrl-1][nLeakCtrl-1].imag()<<endl;
  DEBUGLEAKY <<minL <<" "<<maxL <<endl;
  //DEBUGLEAKY << 0.0 <<" "<< 1.0 <<endl; //for buffer function
	for (j=0; j<nLeakCtrl; j++){
		for (i=0; i<nLeakCtrl; i++){
			DEBUGLEAKY <<leakctrl[i][j]<< " ";
			//DEBUGLEAKY <<Buffer[i][j]<< " ";
			//DEBUGLEAKY <<BufferFunction(zctrl[i][j])<<" ";
		}
		DEBUGLEAKY <<endl;
	}
  DEBUGLEAKY.close();
  
	
	cout <<"plotting calculated leakage" <<endl;
  ofstream DEBUGLEAKY2;	
	double temp3,minL2(ALMOST_INF),maxL2(-ALMOST_INF);
  DEBUGLEAKY2.open  ("calculatedleak.grd");
  DEBUGLEAKY2 << "DSAA"<<endl;
  DEBUGLEAKY2 <<nLeakCtrl <<" "<<nLeakCtrl <<endl;
  DEBUGLEAKY2 <<zctrl[0][0].real()<< " "<<zctrl[nLeakCtrl-1][nLeakCtrl-1].real()<<endl;
	DEBUGLEAKY2 <<zctrl[0][0].imag()<< " "<<zctrl[nLeakCtrl-1][nLeakCtrl-1].imag()<<endl;
	for (i=0; i<nLeakCtrl; i++){
		for (j=0; j<nLeakCtrl; j++){
		  temp3=GetLeakage(zctrl[i][j],0); 
			if (temp3<minL2) {minL2=temp3;}
			if (temp3>maxL2) {maxL2=temp3;}
 		}
	}
	DEBUGLEAKY2 <<minL2 <<" "<<maxL2 <<endl;
	for (j=0; j<nLeakCtrl; j++){
    for (i=0; i<nLeakCtrl; i++){
		  DEBUGLEAKY2 <<GetLeakage(zctrl[i][j],0) <<" ";
		}
		DEBUGLEAKY2 <<endl;
	}
  DEBUGLEAKY2.close();
/*
  cout <<"plotting calculated potential" <<endl;
  ofstream DEBUGLEAKY3;	
	double temp4,minL3(ALMOST_INF),maxL3(-ALMOST_INF);
  DEBUGLEAKY3.open  ("calculatedpot.grd");
  DEBUGLEAKY3 << "DSAA"<<endl;
  DEBUGLEAKY3 <<nLeakCtrl <<" "<<nLeakCtrl <<endl;
  DEBUGLEAKY3 <<zctrl[0][0].real()<< " "<<zctrl[nLeakCtrl-1][nLeakCtrl-1].real()<<endl;
	DEBUGLEAKY3 <<zctrl[0][0].imag()<< " "<<zctrl[nLeakCtrl-1][nLeakCtrl-1].imag()<<endl;
	for (i=0; i<nLeakCtrl; i++){
		for (j=0; j<nLeakCtrl; j++){
		  temp4=GetDischargePotential(zctrl[i][j],0,0.0).real(); 
			if (temp4<minL3) {minL3=temp4;}
			if (temp4>maxL3) {maxL3=temp4;}
 		}
	}
	DEBUGLEAKY3 <<minL3 <<" "<<maxL3 <<endl;
	for (j=0; j<nLeakCtrl; j++){
    for (i=0; i<nLeakCtrl; i++){
		  DEBUGLEAKY3 <<GetDischargePotential(zctrl[i][j],0,0).real() <<" ";
		}
		DEBUGLEAKY3 <<endl;
	}
  DEBUGLEAKY3.close();
*/
	cout <<"ready";
  int temp2;
  cin>> temp2;

	} //end if resistance=-1
	else {change=maxchange=0;}
}


//***********************************************************************
//							 WriteItself
//***********************************************************************
void            CLeakyLayer::WriteItself(ofstream &SOL, const double &t) const{
}
//***********************************************************************
//							 Prepare
//***********************************************************************
void CLeakyLayer::Prepare(){
	double		xtemp,ytemp;
	int				i,j,n;
	double		L1(0.5*PI*(1.0-domain_buffer));
	double		L2(0.5*PI);

	//create sinterm matrix
	SinTerm= new double *[MAX_LEAK_ORDER];
  for(n=0;n<=order;n++){
		SinTerm[n]=new double [nLeakCtrl]; 
		for (i=0; i<nLeakCtrl; i++){
			SinTerm[n][i]=sin((double)(n*i)*PI/(double)(nLeakCtrl));
		}
	}

  Buffer = new double *[MAX_LEAK_CTRL ];
  for(i=0;i<nLeakCtrl;i++){
		xtemp=-(PI/2.0)+(double)(i)*PI/(double)(nLeakCtrl);
		Buffer[i]=new double [nLeakCtrl]; 
		for (j=0; j<nLeakCtrl; j++){
			ytemp=-(PI/2.0)+(double)(j)*PI/(double)(nLeakCtrl);
			Buffer[i][j]=1.0;

			//if in center, do nothing
			if ((fabs(xtemp)<L1) && (fabs(ytemp)<L1)){}
			else{
				//right side
				if      ((xtemp>= L1) && (xtemp<= L2)){Buffer[i][j]*= 0.5*cos(PI*(xtemp-L1)/(L2-L1))+0.5;}
				//left side
				else if ((xtemp<=-L1) && (xtemp>=-L2)){Buffer[i][j]*= 0.5*cos(PI*(xtemp+L1)/(L1-L2))+0.5;}

				if      ((ytemp>= L1) && (ytemp<= L2)){Buffer[i][j]*= 0.5*cos(PI*(ytemp-L1)/(L2-L1))+0.5;}
				//left side
				else if ((ytemp<=-L1) && (ytemp>=-L2)){Buffer[i][j]*= 0.5*cos(PI*(ytemp+L1)/(L1-L2))+0.5;}		
				//else {Buffer[i][j]=2.0;}
			}

		}
	}
}
//***************************************************************************
void CLeakyLayer::Destroy(){
	if (globaldebug){cout <<"DESTROYING LEAKY LAYER STATIC DATA"<<endl;}
	int i,n;
  for(n=0;n<=order;n++){
		delete [] SinTerm[n];
	}
	delete [] SinTerm;
  for(i=0;i<nLeakCtrl;i++){
		delete [] Buffer[i];
	}
	delete [] Buffer;
}
//***************************************************************************
void CLeakyLayer::SetPrecision (const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("SetPrecision::Improper precision level specified",BAD_DATA);}
  switch(Precision){
		case(0): {order=5;  fold=1.0; break;}
		case(1): {order=10; fold=1.5; break;}
		case(2): {order=30; fold=1.8; break;}
		case(3): {order=40; fold=2.0; break;}
	  case(4): {order=60; fold=2.2; break;}
	  case(5): {order=80; fold=2.5; break;}
	  case(9): {order=100;fold=3.0; break;}
	}
}





























//***********************************************************************
//							 Buffercreator
//***********************************************************************
/*double            CLeakyLayer::BufferFunction(cmplex z){

	//X from -length to length
	//Y from -length to length 
	double L1(0.5*length/(1.0+domain_buffer));
	double L2(0.5*length);
  double temp(1.0);
  cmplex Z(z-(zbl+cmplex(length/2.0,length/2.0)));
  cmplex z1,z2;


	//if in center, do nothing
	if ((abs(Z.real())<L1) && (abs(Z.imag())<L1)){ 
			return temp;
	}
	//right side
  if ((Z.real()>=L1) && (Z.real()<=L2)){
			temp*= 0.5*cos(pi*(Z.real()-L1)/(L2-L1))+0.5;
	}
	//left side
  else if ((Z.real()<=-L1) && (Z.real()>=-L2)){
			temp*= 0.5*cos(pi*(Z.real()+L1)/(L1-L2))+0.5;
	}
  if ((Z.imag()>=L1) && (Z.imag()<=L2)){
			temp*= 0.5*cos(pi*(Z.imag()-L1)/(L2-L1))+0.5;
	}
	//left side
  else if ((Z.imag()<=-L1) && (Z.imag()>=-L2)){
			temp*= 0.5*cos(pi*(Z.imag()+L1)/(L1-L2))+0.5;
	}
	return temp;
}*/


