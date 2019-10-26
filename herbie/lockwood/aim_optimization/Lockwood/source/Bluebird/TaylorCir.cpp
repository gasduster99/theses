#include "TaylorCir.h"
//************************************************************************
//                           CTaylorCir
//													 CONSTRUCTORS
//************************************************************************
CTaylorCir::CTaylorCir(){}
//------------------------------------------------------------------------
CTaylorCir::CTaylorCir(CSingleLayerABC *pLay, 
											 cmplex			zcenter, 
											 double			rad, 
											 int				ord, 
											 double			OS){
	int m,n;
	double angle;

	pLayer=pLay;
	zcen=zcenter;
	radius=rad;
	order=ord;
	fold=OS;

	ncircontrol=(int)(order*fold);

  //allocate memory for dynamic arrays (rest allocated in fill routine)
  pTaylorCoeff= new cmplex [order];  
  zctrl=        new cmplex [ncircontrol]; 

	//fill array of control points, initialize all coeff to 0
  for (m=0; m<ncircontrol; m++){                   
		angle=2.0*PI*(double)(m)/(double)(ncircontrol);
    zctrl[m]= radius*cmplex(cos(angle), sin(angle))+zcen;
	}

	for (n=0; n<order; n++){pTaylorCoeff[n]=0.0;} 

}
//------------------------------------------------------------------------
CTaylorCir::~CTaylorCir(){
	delete [] zctrl;
	delete [] pTaylorCoeff;
}
//************************************************************************
//													 STATIC MEMBER FUNCTIONS
//************************************************************************
bool CTaylorCir::Active=false;
//------------------------------------------------------------------------
void CTaylorCir::SetPrecision(const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("SetPrecision::Improper precision level specified",BAD_DATA);}
	switch(Precision){
		case(0): {order=5;  fold=1.1; break;}
		case(1): {order=10; fold=1.4; break;}
		case(2): {order=20; fold=1.5; break;}
		case(3): {order=30; fold=1.5; break;}
	  case(4): {order=40; fold=1.6; break;}
	  case(5): {order=60; fold=2.0; break;}
	  case(9): {order=100;fold=3.0; break;}
	}
}
//------------------------------------------------------------------------
CTaylorCir *CTaylorCir::Parse(ifstream &input, int &l,  
															CSingleLayerABC *pLay, int defaultprecision){
	//string TaylorCir 
	//double xcen double ycen double radius
	//& [precision]
	cmplex			thiszcen(0);
	double			thisrad(1),thisfold(1);
	CTaylorCir *TC;
	bool				done(false),eof(false);
	int					Len,thisorder,thisprec(defaultprecision);
	char			 *s[MAXINPUTITEMS];
	
	TC=NULL;
  if (parserdebug) {cout << "Taylor Series Circle"<<endl;}  eof=TokenizeLine(input,s,Len); l++;
	do{
		if (Len==3){
			thiszcen=s_to_c(s[0],s[1]);
			thisrad=s_to_d(s[2]);
		}
		else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			if (Len==2){thisprec=s_to_i(s[1]);}
			CTaylorCir::SetPrecision(thisprec,thisorder,thisfold);
			TC=new CTaylorCir(pLay,thiszcen,thisrad,thisorder,thisfold);
		}
		else if (Len==0){																				eof=TokenizeLine(input,s,Len); l++;}
		else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	return TC; 
}
//------------------------------------------------------------------------
void CTaylorCir::FillCircle(CAnalyticElem **ElemArray, int numElems){
	for (int i=0; i<numElems; i++){
		if (ElemArray[i]->PartInCircle(zcen,radius)){
		  AddToCir(ElemArray[i]);
		}
	}
}
//************************************************************************
//													 MEMBER FUNCTIONS
//************************************************************************
void CTaylorCir::AddToCir(CAnalyticElem *Elemptr){
  size=size+1;																								//increment size
  CAnalyticElem **ptrArray=new CAnalyticElem *[size];					//allocate memory 
  for (int i=0; i<(size-1); i++){ptrArray[i]=pElemArray[i];}	//copy array
  ptrArray[size-1]=Elemptr;																		//add new pointer
  if (size>1){delete [] pElemArray;}													//delete old array
  pElemArray=&ptrArray[0];																		//redirect pointer
}
//------------------------------------------------------------------------
void CTaylorCir::SetTaylor(const double &t){
	double angle,*Pot;
	cmplex sum;
  int i,m,n;

	Pot=new double[ncircontrol];
	for (m=0; m<ncircontrol; m++){
		Pot[m]=pLayer->GetDischargePotential(zctrl[m],t).real();
		for (i=0; i<size; i++){
			Pot[m]-=pElemArray[i]->GetDischargePotential(zctrl[m],t).real();
		}
    pTaylorCoeff[0]+=Pot[m]/(double)(ncircontrol);
	}
	for (n=1; n<order; n++){
		sum=0.0;
		for (m=0; m<ncircontrol; m++){
			angle=2.0*PI*(double)(m)/(double)(ncircontrol); 
			sum+=Pot[m]*cmplex(cos((double)(n)*angle),-sin((double)(n)*angle)); 
		}
		pTaylorCoeff[n]=2.0*sum/(double)(ncircontrol);
	}
	delete [] Pot;
}
//************************************************************************
cmplex CTaylorCir::GetDischargePotential(const cmplex &z,const double &t) const{
	cmplex omega(0.0);
	for (int i=0; i<size; i++){
		omega+=pElemArray[i]->GetDischargePotential(z,t);
	}
	omega+=Taylor((z-zcen)/radius,order,pTaylorCoeff);
	return omega;
}
//------------------------------------------------------------------------
cmplex CTaylorCir::GetW(const cmplex &z,const double &t) const{
	cmplex omega(0.0);
	for (int i=0; i<size; i++){
		omega+=pElemArray[i]->GetW(z,t);
	}
	omega+=TaylorDer((z-zcen)/radius,order,pTaylorCoeff,radius);
	return omega;
}
//************************************************************************
bool CTaylorCir::IsInside(const cmplex z) const{
  if (abs(z-zcen)<=radius){return true; } 
	else                    {return false;}
}  
//***********************************************************************
//											AddINTToBlock
//************************************************************************
/*void CSuperblock::AddToBlock(CAnalyticElem *Elemptr, int seg){
  sizeINT=sizeINT+1;                                           //increment size
  CAnalyticElem **ptrArray = new CAnalyticElem *[size];  //allocate memory
	int *segArray=new int[size];
  for (int i=0; i<(size-1); i++){
		ptrArray[i]=pElemArrayINT[i];
		segArray[i]=pElemSegmentINT[i];
	}                                                      //copy array
  ptrArray[size-1]=Elemptr;                              //add new pointer
  segArray[size-1]=seg;
	delete [] pElemSegmentINT;
  delete [] pElemArrayINT;                                  //delete old array
	pElemSegment=&segArray[0];
  pElemArray=  &ptrArray[0];                             //redirect pointer
} */