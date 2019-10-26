#include "Leaky.h"

/************************************************************************
                           CLeakyWall
************************************************************************
                           CONSTRUCTOR
************************************************************************/
CLeakyWall::CLeakyWall(char						 *Name, 
							         const CSingleLayerABC *pLay,
											 const cmplex		 *points,
											 const int				NumOfLines,
											 const double			cond,
											 const double			thick, 
											 const int				prec)
          :CStringElem(Name,pLay,points,NumOfLines,false,DOUBLET,prec){

	if (cond<0){ExitGracefully("CLeakyWall::CLeakyWall: negative conductivity not allowed",BAD_DATA);}
	Kwall=cond;    
	if (thick<=0){ExitGracefully("CLeakyWall::CLeakyWall: zero or negative thickness not allowed",BAD_DATA);}
	thickness=thick;
}
//------------------------------------------------------------------------
CLeakyWall::~CLeakyWall(){}
/************************************************************************
                           ACCESSORS
***********************************************************************/
double CLeakyWall::GetWallCond() const{return Kwall;}

/************************************************************************
                           SolveItself
************************************************************************/
void CLeakyWall::SolveItself(double &change, double &objective ,double t){

  double  error(0.0),maxchange(0.0),maxobjective(0.0),K;
  cmplex  WOther,Z,unitY;
  int     i,m;

  for (i=0; i<NLines; i++){
    K=pLayer->GetCond(zctrl[i][int(nlinecontrol/2)]);
    
		//Solve Leaky Wall Element
    //******************************************************************
    unitY=IM*(ze[i+1]-ze[i])/abs(ze[i+1]-ze[i]);

		disabled [i]=true;

    //set right hand side for system of equations 
    for (m=0; m<nlinecontrol; m++){

      //identify Total,Local normal flux at control points
      WOther=pLayer->GetW(zctrl[i][m],t);

      //set rhs
			rhs[m]=-(conj(unitY)*conj(WOther)).real();
    }

		disabled[i]=false;

    //solve for jump coefficients
		SetConstraints(i,Kwall/(K*thickness),0,1,0,UNCONSTRAINED,0.0,0.0);
    GenSolve(JumpCoeff[i],i,ltype,false,1.0,objective,change);
		upperswap(maxchange   ,change   );
		upperswap(maxobjective,objective);

		SetFarFieldCoeff(i);
		if ((pBlocks[i]!=NULL) && (pBlocks[i]->IsOn())){
			pBlocks[i]->Update(myBlockIDs[i],i,t);
		}
  }

  change=maxchange;
  objective=maxobjective;
}
//------------------------------------------------------------------------
double CLeakyWall::GetMaxError(const double &t) const{
	double error(0.0);

	return error;
}
/************************************************************************
                           READ/WRITE
************************************************************************/
void CLeakyWall::WriteOutput(const double &t) const{
  //write errors file
  double error(0);
 /* ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	
	for (int i=0; i<NLines; i++){
		for (int m=0; m<nlinecontrol; m++){
			ERRORS <<zctrl[i][m].real()<<","<<zctrl[i][m].imag()<<",";
		  ERRORS <<FLUX_ERROR_TAG    <<","<<0<<",";
			ERRORS <<error<<",";
			ERRORS <<error<<endl; //TMP DEBUG
		}
	}
	ERRORS.close();*/
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
  string "LeakyWall", string name 
  double thickness
	double Kwall
	{double x double y}x(numlines+1)
  &[int precision] 
----------------------------------------------------------------------*/
CLeakyWall *CLeakyWall::Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name){
  

	CLeakyWall  *pLeaky=NULL;

  bool     done(false);
  double   thiscond(0),thisthick(1.0);
	int      Len,thisprec,nlines(0);
  cmplex	 stringp[MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Leaky wall element"<<endl;} 

  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if      (Len==1) {thisthick=s_to_d(s[0]); done=true;}
  else             {ImproperFormat(s,l); return NULL; }
	
	if (TokenizeLine(input,s,Len)){return NULL;}; l++; 
	if     (Len==1)  {thiscond=s_to_d(s[0]); done=true;}
	else             {ImproperFormat(s,l); return NULL;}

	done=false; 
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
  do {
    if (nlines>=MAXLINES) { ExitGracefully("CLeakyWall::Parse- too many lines in leaky wall",TOO_MANY);}
		if ((Len==2) && (strcmp(s[0],"&"))) {
      stringp[nlines]=s_to_c(s[0],s[1]); nlines++;      
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (int i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pLeaky = new CLeakyWall(Name,
													pLay,
													stringp,
													nlines-1,
													thiscond,
													thisthick,
													thisprec); 
			done=true;                                                                            }
    else             {ImproperFormat(s,l); return NULL;                                     }
	} while (!done);

	return pLeaky;
}
//******************************************************************
  /*  //***************************************************************
		//Associated Divergence Elements
    //***************************************************************
		if (CAquifer::Multilayer){
			double *potplus;
			double *potminus;
			double *transplus;
			double *transminus;
      potplus=   new double [nlinecontrol];
      potminus=  new double [nlinecontrol];
			transplus= new double [nlinecontrol];
      transminus=new double [nlinecontrol];

			
			//Solve for divergence doublet-removes jump in head
      //***************************************************************
			//set right hand side for headjump coeff
			for (m=0; m<nlinecontrol; m++){
			  potplus[m]= pLayer->GetDischargePotential(zctrl[i][m]+(unitY*abs(ze[i+1]-ze[i])*MOVE_DIST),t).real();
				potminus[m]=potplus[m]-GetOmJump(i,X[m],t).real();
				transplus[m]=GetTransmissivity(potplus[m],K,T);
        transminus[m]=GetTransmissivity(potminus[m],K,T);
				//rhs[m]=((potplus[m]/transplus[m])-(potminus[m]/transminus[m]));
        rhs[m]=2.0*((potplus[m]/transplus[m])-(potminus[m]/transminus[m])); //works, but why?
			}

			//solve for Leakage jump coefficients
			SetConstraints(i,1,0,0,0,UNCONSTRAINED);
			//GenSolve(LJumpCoeff,i,doublet,false,1.0,true,maxobjective,maxchange);
      if (error>maxobjective){maxobjective=error;}

			//Solve for divergence linesink-removes normal gradient 
      //***************************************************************
			for (m=0; m<nlinecontrol; m++){
				rhs[m]=Kwall/(K*thickness)*((1.0/transplus[m])-(1.0/transminus[m]))*GetOmJump(i,X[m],t).real();
				//cout <<transplus[m]-transminus[m]<<endl;//cout <<"rhs"<<rhs[m];
			}

			//solve for Leakage gradient jump coefficients
			SetConstraints(i,0,1,0,0,ZERO_AT_END);
			//GenSolve(LJumpCoeff,i,LINESINK,false,1.0,true,maxobjective,maxchange);
      if (error>maxobjective){maxobjective=error;}
			
			delete [] transplus;
			delete [] transminus;
			delete [] potplus;
			delete [] potminus;
		}*/