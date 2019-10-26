#include "Drain.h"

/************************************************************************
                           CDrain
************************************************************************
                           CONSTRUCTOR
************************************************************************/
CDrain::CDrain(char						 *Name, 
							 const CSingleLayerABC *pLay,
							 const cmplex		 *points,
							 const int				NumOfLines,
	             const double		  cond,
							 const double			thick, 
							 const int				prec)
       :CStringElem(Name,pLay,points,NumOfLines,false,DIPOLE,prec){

  Kdrain=cond;    
	if (thick<=0){ExitGracefully("CDrain:CDrain: zero or negative thickness not allowed",BAD_DATA);}
	thickness=thick;
}
/************************************************************************
                           ACCESSORS
************************************************************************/
double CDrain::GetDrainCond() const {return Kdrain;}

/************************************************************************
                           SolveItself
************************************************************************/
void CDrain::SolveItself(double &change, double &objective, const double t){
 
  double maxchange(0.0),maxobjective(0.0),error(0.0),Qxother,K;
  cmplex Z,unitX;
  int    i,m;     

  for (i=0; i<NLines; i++){
    unitX=(ze[i+1]-ze[i])/abs(ze[i+1]-ze[i]);
	  K=pLayer->GetCond(zctrl[i][int(nlinecontrol/2.0)]);

		disabled[i]=true;

		//Solve Drain Element (dipole)
    //******************************************************************
    //set right hand side for system of equations 
    for (m=0; m<nlinecontrol; m++){

      //identify non-segment tangential flux at control points 
      Qxother=((pLayer->GetW(zctrl[i][m],t))*unitX).real();

			rhs[m]=-Qxother;
    }
		
		disabled[i]=false;

    //solve for jump coefficients
		SetConstraints(i,K/(Kdrain*thickness),0,1,0,UNCONSTRAINED,0.0,0.0);
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
double CDrain::GetMaxError(const double &t) const{
	double error(0.0);
	cmplex Qtotal;
	for (int i=0;i<NLines;i++){
		cmplex unitX=abs(ze[i+1]-ze[i])/(ze[i+1]-ze[i]);
	  double K=pLayer->GetCond(zctrl[i][int(nlinecontrol/2.0)]);
		for (int m=0; m<nlinecontrol; m+=max((int)(TEST_ERROR_RATIO*nlinecontrol),1)){
      Qtotal=(conj(pLayer->GetW(zctrl[i][m],t))*unitX);
  		upperswap(error,(Qtotal.real()+K/(Kdrain*thickness)*GetOmJump(i,X[m],t).imag())/max(abs(Qtotal),1.0)); //relative flux error
		}
	}
	return error;
}
//***********************************************************************
void CDrain::WriteOutput(const double &t) const{
  //write errors file
  double error(0),K;
	cmplex Qtotal,unitX;

  ofstream ERRORS;
	ERRORS.open("errors.csv", ios::app);
	
	for (int i=0; i<NLines; i++){
		unitX=abs(ze[i+1]-ze[i])/(ze[i+1]-ze[i]);
	  K=pLayer->GetCond(zctrl[i][int(nlinecontrol/2.0)]);
		for (int m=0; m<nlinecontrol; m++){
      Qtotal=(conj(pLayer->GetW(zctrl[i][m],t))*unitX);
  		error=Qtotal.real()+K/(Kdrain*thickness)*GetOmJump(i,X[m],t).imag();
		  ERRORS <<zctrl[i][m].real()<<","<<zctrl[i][m].imag()<<","; //coordinate
			ERRORS <<FLUX_ERROR_TAG    <<","<<0<<",";
			ERRORS <<error                                      <<","; //absolute error
			ERRORS <<error/max(abs(Qtotal),1.0)                 <<endl;//relative error

		}
	}
	ERRORS.close();
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
	string "Drain", string name 
	double thickness
	double Kdrain
	{double x double y}x(numlines+1)
  & [int precision] 
----------------------------------------------------------------------*/
CDrain *CDrain::Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name){


	CDrain  *pDrain;

  bool     done(false);
  double   thisthick(1.0),thiscond(10000);
	int      Len,thisprec,nlines(0),i;
  cmplex   stringp[MAXLINES];
	char    *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Drain element"<<endl;}      
	
  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if      (Len==1) {thiscond=s_to_d(s[0]); done=true;  }
	else             {ImproperFormat(s,l); return NULL;}

  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if      (Len==1) {thisthick=s_to_d(s[0]);  done=true;}
	else             {ImproperFormat(s,l); return NULL;}

  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	done=false;
  do {
    if (nlines>=MAXLINES) {ExitGracefully("CDrain::Parse- too many lines in drain",TOO_MANY);}
		if  (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]); nlines++;       
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			for (i=0; i<nlines; i++){pLay->UpdateExtents(stringp[i]);}
			if (Len==2){thisprec= s_to_i(s[1]);                   }
			else       {thisprec= CAnalyticElem::DefaultPrecision;}
			pDrain = new CDrain(Name,
													pLay,
													stringp,
													nlines-1,
													thiscond,
													thisthick,
													thisprec); 
			done=true;
		}

    else             {ImproperFormat(s,l); return NULL;}
	} while (!done);

	return pDrain;
}
//***********************************************************************
  	//evaluate error
    /*for (m=0; m<nlinecontrol; m++){
		  error=(conj(GetSegmentW(i,zctrl[i][m],t))*unitX).real()+
						(K/(Kdrain*thickness)*GetOmJump(i,X[m],t).imag())-
						rhs[m];
	    if (error>maxobjective){maxobjective=error;}
		}*/
			/*//evaluate error
			for (m=0; m<nlinecontrol; m++){
				error=GetHeadGradJump(i,X[m],t)-rhs[m];
				//cout <<error<<endl;
	      if (error>maxobjective){maxobjective=error;}
			}*/

		//Solve Corresponding Divergence Element
    //******************************************************************
		/*if (CAquifer::Multilayer){
			double Pot;
			double trans;
			
			//may wish to set extra multilogs to represent leakage wells at ends
			// should be equal to sum(n=0 to N) a_n_drain (+-1)^n
			//set right hand side for headjump coeff
			for (m=0; m<nlinecontrol; m++){
				Pot=pLayer->GetDischargePotential(zctrl[i][m],t).real();
				if (Pot<0.5*K*pow(T,2.0)){trans=K*pow(2*Pot/K,0.5);} //check for 0 cond?
        else {	 								  trans=K*T;}
				rhs[m]=GetWJump(i,X[m],t).imag()/trans;
			}
	
			//solve for Leakage jump coefficients
			SetConstraints(i,0,1,0,0,UNCONSTRAINED);
			//GenSolve(LJumpCoeff,i,ltype,false,1.0,true,objective,change);  //CMPTODBchange
		  if (objective>maxobjective){maxobjective=objective;}

		}*/
