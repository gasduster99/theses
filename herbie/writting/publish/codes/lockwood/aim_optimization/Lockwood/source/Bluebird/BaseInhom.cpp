#include "BaseInhom.h"

/************************************************************************
                           CBaseInhom
************************************************************************/
//Constructor
CBaseInhom::CBaseInhom(char						 *Name, 
											 const CSingleLayerABC *pLay,
											 const cmplex		 *points,
											 const int				NumOfLines,
	                     const double		  base,
											 const double		  thick,
											 const int				prec)
           :CStringElem(Name,pLay,points,NumOfLines,true,DOUBLET,prec){
  int m;
	
	Bin=base;  	   
	Tin=thick;

	pBaseZone=  new CPolyPropZone(base_elevation,points,NumOfLines,base);
	pThickZone= new CPolyPropZone(layer_thickness,points,NumOfLines,thick);
  
	c[2]=         new double       [nlinecontrol];  
	c[0]=         new double       [nlinecontrol];   

	for (m=0;m<nlinecontrol;m++){c[2][m]=0.0;c[0][m]=0.0;}

	//AllInhoms[NumInhoms]=this;
	//NumInhoms++;

	if (!DynArrayAppend((void**&)(AllInhoms),(void*)(this),NumInhoms)){
		 ExitGracefully("CBaseInhom::Constructor: creating NULL element",BAD_DATA);};
}
//-----------------------------------------------------------------------
CBaseInhom::~CBaseInhom(){
	if (globaldebug){cout <<"   DESTROYING BASE/THICKNESS INHOMOGENEITY "<<name<<endl;}
}
//***********************************************************************
double			 CBaseInhom::GetBase()      const{return Bin;}
//-----------------------------------------------------------------------
double			 CBaseInhom::GetThick()     const{return Tin;}
//-----------------------------------------------------------------------
CPropZone*	 CBaseInhom::GetBaseZone () const{return pBaseZone;}
//-----------------------------------------------------------------------
CPropZone*   CBaseInhom::GetThickZone() const{return pThickZone;}
//***********************************************************************
CBaseInhom **CBaseInhom::AllInhoms = NULL; 
int          CBaseInhom::NumInhoms = 0; 
//------------------------------------------------------------------------
void         CBaseInhom::Destroy(){
	delete [] AllInhoms;
}
//***********************************************************************
//                           SolveItself
//***********************************************************************
void CBaseInhom::SolveItself(double &change,double &objective,double t){
  double maxchange(0.0),maxobjective(0.0),error(0.0),Bout,K,Tout;
	double LastJump,delB,delT,LastPhiOut,LastPhiIn,alpha,HeadIn,HeadOut;
	cmplex WOther, UnitY;
	double OutPot;
  int    i,m;
	
  for (i=0; i<NLines; i++){
		if (!disabled[i]){
		UnitY=IM*(ze[i+1]-ze[i])/abs(ze[i+1]-ze[i]);
   	K=   pLayer->GetCond (zctrl[i][int(nlinecontrol/2)]);
		Tout=pLayer->GetThick(zctrl[i][int(nlinecontrol/2)]-(ze[i+1]-ze[i])*IM*MOVE_DIST);
		Bout=pLayer->GetBase (zctrl[i][int(nlinecontrol/2)]-(ze[i+1]-ze[i])*IM*MOVE_DIST);
    delB=Bout-Bin;
		delT=Tout-Tin;


		//cout <<"delB: "<<delB<<" Tin: "<<Tin << " Tout: " << Tout<<endl;

		//Solve Base Inhomogeneity Doublet 
    //******************************************************************
    //set right hand side for system of equations 

		char typ;
    for (m=0; m<nlinecontrol; m++){

      //identify Total,Local potential at control points
			LastJump=GetOmJump(i,X[m],t).real();
			disabled[i]=true;
			OutPot=max(pLayer->GetDischargePotential(zctrl[i][m],t).real(),0.0);
			disabled[i]=false;

			LastPhiOut=OutPot-0.5*LastJump;
			LastPhiIn =OutPot+0.5*LastJump;

			HeadOut=ConvertToHead(LastPhiOut,K,Tout,pLayer->GetSeaLevel()-Bout,pLayer->GetSaltwaterSG());
			HeadIn =ConvertToHead(LastPhiIn ,K,Tin ,pLayer->GetSeaLevel()-Bin ,pLayer->GetSaltwaterSG());

			//if ((LastPhiOut<=0) && (LastPhiIn<=0)){                 //head below base 
			if (false){
				//ExitGracefully("CBaseInhom::Negative potential along base inhomogeneity boundary. Cannot solve",RUNTIME_ERR);}
				c[2][m]=0.0;
				c[0][m]=1.0;
				rhs[m]=0.0; //no jump in potential
				typ='f';
			}
			                 
      //else if (LastPhiOut*LastPhiIn<0.0){                     //if one potential is negative
			//no flow boundary QY=0
		  //else if (true){
			else if (false){
				c[2][m]=1.0;
				c[0][m]=0.0;
			
				//identify normal flux from all other elements at control points
				disabled[i]=true;
				WOther=pLayer->GetW(zctrl[i][m],t);
				disabled[i]=false;
				
				typ='e';
				//cout << TotalW-LocalW <<endl;
				//rhs[m]=-((TotalW-LocalW)*UnitY).real();
				rhs[m]=(WOther*(ze[i+1]-ze[i])/abs(ze[i+1]-ze[i])).imag();
			}
			else{																	                  //continuity of head phi_in=phi_out
				c[2][m]=0.0;
				c[0][m]=1.0;																					//four situations:____________________________
				if (IsConfined(LastPhiOut,K,Tout)){				            //confined outside____________________
					if (IsConfined(LastPhiIn,K,Tin)){			              //confined inside
						alpha=Tin/Tout;                                   typ='a';
						rhs[m]=K*Tin*(delT+2*delB)+2*(alpha-1.0)*OutPot;
					}
					else{																								//unconfined inside
						alpha=sqrt(0.5*K*Tout*Tout/LastPhiIn);            typ='b';
						rhs[m]=K*Tout*(2.0*delB+Tout)-Tout*sqrt(2.0*K*LastPhiIn)+2.0*(1.0-alpha)*OutPot;
					}
				}
				else{																									//unconfined outside__________________
					if (IsConfined(LastPhiIn,K,Tin)){										//confined inside 
						alpha=sqrt(0.5*K*Tin*Tin/LastPhiOut);             typ='c';  
						rhs[m]=K*Tin *(2.0*delB-Tin )+Tin*sqrt(2.0*K*LastPhiOut)+2*(alpha-1.0)*OutPot;
					}
					else{																								//unconfined inside
						alpha=delB*sqrt(K/(8.0*LastPhiOut));	            typ='d';
 						rhs[m]=0.5*K*delB*delB+delB*sqrt(0.5*K*LastPhiOut)+2.0*alpha*OutPot;
					}  
				}
				rhs[m]/=(1+alpha);
			}
			//cout << typ<< " ";
    }}
		//cout <<endl;

		//solve for jump coefficients 
		SetConstraints(i,1,0,1,0,UNCONSTRAINED,0.0,0.0);
    GenSolve(JumpCoeff[i],i,ltype,true,1.0,objective,change);
		upperswap(maxchange,change);
		upperswap(maxobjective,objective);

		//evaluate objective function
    /*for (m=0; m<nlinecontrol; m++){
			//check
			error=pLayer->GetHead(zctrl[i][m]-(ze[i+1]-ze[i])*IM*MOVE_DIST,t)-
			      pLayer->GetHead(zctrl[i][m]+(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
		}*/

    SetFarFieldCoeff(i);
		if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){pBlocks[i]->Update(myBlockIDs[i],i,t);}
  }

  change=maxchange;
  objective=maxobjective;
}
//------------------------------------------------------------------------
double CBaseInhom::GetMaxError          (const double &t) const{
	double h1,h2,error(0.0);
	for (int i=0;i<NLines;i++){
		for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
			h1=pLayer->GetHead(zctrl[i][m]-(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
			h2=pLayer->GetHead(zctrl[i][m]+(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
			if (h1*h2>0){
				//upperswap(error,fabs(h1-h2)/max(h2-Bin,1.0)); //relative error
				upperswap(error,fabs(h1-h2));                 //absolute error
			}
			else{
				//check normal flux
			}
		}
	}
	return error;
}
//************************************************************************
//                           READ/WRITE
//************************************************************************
void CBaseInhom::WriteOutput(const double &t) const{
  //write errors file
	double hin,hout;
	double Bout;
  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	
	for (int i=0; i<NLines; i++){
		Bout  =pLayer->GetBase(zctrl[i][int(nlinecontrol/2)]-(ze[i+1]-ze[i])*IM*MOVE_DIST);
		for (int m=0; m<nlinecontrol; m++){
			hin =pLayer->GetHead(zctrl[i][int(nlinecontrol/2)]-(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
			hout=pLayer->GetHead(zctrl[i][int(nlinecontrol/2)]+(ze[i+1]-ze[i])*IM*MOVE_DIST,t);
		  ERRORS <<zctrl[i][m].real()<<","<<zctrl[i][m].imag()<<","; //coordinate
			ERRORS <<HEAD_ERROR_TAG    <<","<<0<<",";
			ERRORS <<(hin-hout)                                 <<","; //absolute error
			ERRORS <<(hin-hout)/min(hin-Bin,hout-Bout)          <<endl;//relative error
		}
	}
	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "BaseInhom", string name 
	double B double T
	{double x double y}x(numlines+1)
	&[int precision]
----------------------------------------------------------------------*/
CBaseInhom *CBaseInhom::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){


	CBaseInhom  *pInhom=NULL;

  bool     done(false);
  double   thisB(0),thisT(0);
	int      Len,thisprec,nlines(0),i;
	char    *s[MAXINPUTITEMS];
  cmplex	 stringp[MAXLINES];

	if (parserdebug){cout << "Base Inhom Element"<<endl;} 
	
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if      (Len==2) {thisB=s_to_d(s[0]); thisT=s_to_d(s[1]);}
  else             {ImproperFormat(s,l); return NULL;}
	
	done=false;  
  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	do {
    if (nlines>=MAXLINES) { ExitGracefully("CBaseInhom::Parse- too many lines in base inhomogeneity",TOO_MANY);}
	  if (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]);nlines++;       
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			stringp[nlines]=stringp[0]; nlines++;
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if  (Len==2){thisprec= s_to_i(s[1]);}
			else        {thisprec= CAnalyticElem::DefaultPrecision;}
			pInhom= new CBaseInhom(Name,
														 pLay,
														 stringp,
														 nlines-1,
														 thisB,
														 thisT,
														 thisprec);
			done=true;
		}		
    else            {ImproperFormat(s,l); return NULL;}
	} while (!done);

  return pInhom;
}
//************************************************************************
//                          Sift Thru Inhoms
//************************************************************************
void CBaseInhom::SiftThroughInhoms(){
  //Sifts through base/thickness inhomogeneities and disables one of the shared sides
	int			Points1,Points2;
	int			i,j,k,l, count(0);
	double	Base1,Base2,Thick1,Thick2;

  if (NumInhoms<2){return;}
	else {
		cout <<"Sifting through " << NumInhoms << " Base/Thickness Inhomogeneities..."<<endl;
		for (i=NumInhoms-1; i>=0; i--){
			Points1=AllInhoms[i]->GetNumSegs()+1;
			Base1=  AllInhoms[i]->GetBase();
			Thick1=	AllInhoms[i]->GetThick();
			for (j=0; j<i; j++){
				Points2=AllInhoms[j]->GetNumSegs()+1;	
				Base2=  AllInhoms[i]->GetBase();
			  Thick2=	AllInhoms[i]->GetThick();
				for (k=0; k<Points1-1; k++){
				  for (l=0; l<Points2-1; l++){
						if (((AllInhoms[i]->GetZ(k)==AllInhoms[j]->GetZ(l))   && (AllInhoms[i]->GetZ(k+1)==AllInhoms[j]->GetZ(l+1))) ||
							  ((AllInhoms[i]->GetZ(k)==AllInhoms[j]->GetZ(l+1)) && (AllInhoms[i]->GetZ(k+1)==AllInhoms[j]->GetZ(l)))){
							if ((Base1>=Base2) || (Thick1>=Thick2))  {AllInhoms[i]->SetAsDisabled(k);count++;}
							else                                     {AllInhoms[j]->SetAsDisabled(l);count++;}
							//if ((Base1==Base2) && (Thick1==Thick2))  {AllInhoms[i]->SetAsDisabled(k);count++;
							//																				  AllInhoms[j]->SetAsDisabled(l);count++;}
						} 
					}
				}
			}
		}
	}
	cout << "..."<<count << " sides disabled"<<endl;
}