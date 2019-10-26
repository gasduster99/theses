

#include "CationExchange.h"
#include "MatrixInclude.h"

/*************************************************************
 Constructors
**************************************************************
------------------------------------------------------------*/
CCationExchange::CCationExchange():
                 CReactionScheme(true){
	nCations=0;
	pSoilArray=NULL;
	pSpeciesArray=NULL;
	CationInd=NULL;
	valence  =NULL;
	SelCo    =NULL;
	C        =NULL;
	S        =NULL;
	T        =NULL;
	Q				     =0.0;
  AnionConc    =0.0;
	DiActivity   =0.0;
  MonoActivity =0.0;
}
//-----------------------------------------------------------
CCationExchange::CCationExchange(int     NumCations, 
																 int    *indices, 
																 double *selectivities, //mol/L
																 int    *valences,      
																 double  capacity,      //meq/g ~ eq/kg
																 double  AnionicConc,   //mol/L ~ eq/L 
																 int     soiltype_index):
                 CReactionScheme(true){   
	int i,j;
	nCations     =NumCations;
	Q				     =capacity/1000;//eq/kg->eq/g
  AnionConc    =AnionicConc;
	soiltype     =soiltype_index;
	DiActivity   =0.0;
  MonoActivity =0.0;

	pSoilArray   =NULL;
	pSpeciesArray=NULL;

	CationInd=new int    [nCations];
	valence  =new int    [nCations];
	SelCo    =new double [nCations];
	C        =new double [nCations];
	S        =new double [nCations];
	Cneg     =new double [nCations];
	Sneg     =new double [nCations];
	T        =new double [nCations];

  for (i=0; i<nCations; i++){
		CationInd[i]=indices      [i];
		valence  [i]=valences     [i];
    SelCo    [i]=selectivities[i];
	  C        [i]=0.0;
		S        [i]=0.0;
		T        [i]=0.0;
		Cneg     [i]=0.0;
		Sneg     [i]=0.0;
	}

  //Reserve memory for Jacobian matrix 
	J  =new double   *[nCations*2];
	B  =new double    [nCations*2];
	dC =new double    [nCations*2];
	if (dC==NULL){RxNExitGracefully("CCationExchange::React",OUT_OF_MEMORY);}
	for (i=0; i<nCations*2; i++){
		J[i]=new double [nCations*2];	
		if (J[i]==NULL){RxNExitGracefully("CCationExchange::React",OUT_OF_MEMORY);}
		for (j=0; j<nCations*2; j++){
		  J[i][j]=0.0;
		}
		B[i]=0.0;
		dC[i]=0.0;
	}
	if (Q<0.){
		RxNExitGracefully("CCationExchange::Constructor: Need positive cation exchange capacity",BAD_DATA);}
	if (AnionConc<0.0){
		RxNExitGracefully("CCationExchange::Constructor: Need positive anion concentration",BAD_DATA);}
}
//-----------------------------------------------------------
CCationExchange::~CCationExchange(){
	delete [] CationInd;
	delete [] valence;
	delete [] SelCo;
	delete [] C;
	delete [] S;
	delete [] T;
	delete [] Cneg;
	delete [] Sneg;

	if (J!=NULL){for (int i=0;i<nCations*2; i++){delete [] J[i];}}
	delete [] J;
	delete [] B;
	delete [] dC;
}
/***************************************************************
 Initialize
**************************************************************
Initialize Reaction
------------------------------------------------------------*/
void CCationExchange::Initialize(){

	if (nCations>nspecies){
		RxNExitGracefully("CCationExchange::Constructor: More cations than existing species",BAD_DATA);}
	if (pSpeciesArray==NULL){
		RxNExitGracefully("CCationExchange::Constructor: Species Array not set",BAD_DATA);}
	if (pSoilArray==NULL){
		RxNExitGracefully("CCationExchange::Constructor: Soil Array not set",BAD_DATA);}
	//Quality check
  for (int i=0; i<nCations; i++){
		if (CationInd[i]>nspecies){
			RxNExitGracefully("CCationExchange::Constructor: Invalid cation species index",BAD_DATA);}
		if ((valence[i]!=1) && (valence[i]!=2)) {
			RxNExitGracefully("CCationExchange::Constructor: only monovalent and divalent species allowed",BAD_DATA);}
		if (pSpeciesArray[CationInd[i]]->GetMolecularWeight()==0){
			RxNExitGracefully("CCationExchange::Constructor: molecular weights not initialized",BAD_DATA);}
	}
	//Make sure initial sorbed conditions match cation exchange capacity 


	if (valence[0]!=1){RxNExitGracefully("CCationExchange::Constructor: divalent reference species not allowed",BAD_DATA);}
}
/*************************************************************
 Calculate Activity Coefficients
**************************************************************
from actcoe soubroutine in Mouser
Concentration, C, and anion concentration, AnionConc, must be in units of mol/L !
------------------------------------------------------------*/
void CCationExchange::CalculateActivityCoeff(bool estimateonly){

	static double CationConc,tmp;    //Ionic strength due to cations
	int    s;                    //species counter

	CationConc=0.0;
	if (estimateonly){
		CationConc=AnionConc;
	}
	else{
		for (s=0; s<nCations; s++){
			CationConc+=0.5*C[s]*valence[s]*valence[s];  //[mol/L]
		}
	}

  //ionic strength in moles per kilogram solvent (water)=[mol/L]
	IonicStrength=AnionConc+CationConc; 

	//Davies Equation: calculates the log monovalent and divalent activity coeffs. ()
	tmp=-0.5*(pow(IonicStrength,0.5)/(1+pow(IonicStrength,0.5))-0.2*IonicStrength);

	MonoActivity=pow(10,tmp);
	DiActivity  =pow(10,4.0*tmp);

}
/***************************************************************
 GetEquilibriumAqueous
**************************************************************
Calculates Equilibrium Aqueous concentrations given sorbed concentrations and 
------------------------------------------------------------*/
void CCationExchange::GetEquilibriumAqueous(Ironclad1DArray  Cs,
		                                        Writeable1DArray Caq_new, 
																						const double    &Caq_ref,
		 						                            const int        soil_index,
																						const double    &n){
  int    i,iter;
	double molwt,lastionstr,sorbed(0.0);	
	bool done;

  if (soil_index!=soiltype){return;} //if not in proper soiltype, does not calculate

	double pb=pSoilArray[soil_index]->GetBulkDryDensity()*1000; //convert density into g/L 
	
	//Translate model species to local species with molar concentrations (mol/L)
	//Calculate total concentrations

	

	for (i=0; i<nCations; i++){
		molwt=pSpeciesArray[CationInd[i]]->GetMolecularWeight();
    
    if (i==0){C[i]=Caq_ref/(molwt*1000);}//mg/L  * mol/mg =mol/L
		else     {C[i]=0.0;                 }

		S[i]=Cs [CationInd[i]]/(molwt*1000)/1000; //mg/kg * mol/mg * kg/g =mol/g

		sorbed+=valence[i]*S[i];
	}

	cout <<"Sorbed: "<<sorbed<<" Q:"<<Q<<" "<<pb<<" "<<C[0]<<endl;
	if (fabs((sorbed-Q))>0.0001*Q){
		RxNExitGracefully("CCationExchange::GetEquilibriumAqueous: not enough sorbed material to meet capacity",BAD_DATA);}

  CalculateActivityCoeff(true); //estimates as being equal to anion conc
  done =false;
  lastionstr=IonicStrength;
  iter=0;
	while (!done){
		cout <<"iter: "<<iter<<" IonicStr:"<<IonicStrength<<endl;

		
		for (i=1; i<nCations; i++){
      cout<<C[i]<<" "<<S[i]<<endl;
			if (valence[i]==1){
				C[i]=C[0]*S[i]/(S[0]*SelCo[i]);
			}
			else if (valence[i]==2){
				C[i]=2.0*Q*MonoActivity*MonoActivity*C[0]*C[0]*S[i]/(DiActivity*S[0]*S[0]*SelCo[i]);
			}
		}
		
    CalculateActivityCoeff(false);
		if (fabs(lastionstr-IonicStrength)/IonicStrength<CAT_EX_TOLERANCE){done=true;}
    lastionstr=IonicStrength;
		iter++;
		if (iter>400){done=true;}
	}

	//Convert back to mass concentration (currently not weighted) 
	for (i=0; i<nCations; i++){
    Caq_new[CationInd[i]]=C[i]*(pSpeciesArray[CationInd[i]]->GetMolecularWeight()*1000);       //mol/L *mg/mol=mg/L
		//Cs_new [CationInd[i]]=S[i]*(pSpeciesArray[CationInd[i]]->GetMolecularWeight()*1000)*1000;  //mol/g *mg/mol*g/kg=mg/kg
	}

}
/***************************************************************
 React
**************************************************************
combination of newton and Exchange routines in mouser
------------------------------------------------------------*/
void CCationExchange::React(Ironclad1DArray  Caq,
														Ironclad1DArray  Caq_end,
														Writeable1DArray Caq_new,
														Ironclad1DArray  Cs,
														Ironclad1DArray  Cs_end,
														Writeable1DArray Cs_new,	
														const int        soil_index,
														const double    &n,  
														const double    &tstep){
  int    i,k;	
	int    iter(0);
	bool   done(false);
	static double rel_error, relax, decay,molwt;
  static double Sorbed_to_sites,diff, change,available,pb;

	for (i=0;i<nspecies;i++){
		Cs_new [i]=Cs_end[i];
		Caq_new[i]=Caq_end[i];
	}
  if (soil_index!=soiltype){return;} //if not in proper soiltype, does not calculate

	pb=pSoilArray[soil_index]->GetBulkDryDensity()*1000; //convert density into g/L 
	
	//To avoid problems with reaction solution, set negative aqueous concentrations to zero
  //Save negative values to preserve mass balance in transport solution
	//Translate model species to local species with molar concentrations (mol/L)
	//Calculate total concentrations
	//Currently not weighted for start and end concentrations
  available=0.0;
	for (i=0; i<nCations; i++){

		decay=pSpeciesArray[CationInd[i]]->GetDecayRate();
		molwt=pSpeciesArray[CationInd[i]]->GetMolecularWeight();

    C[i]=Caq_end[CationInd[i]]*exp(-decay*tstep)/(molwt*1000);      //mg/L  * mol/mg =mol/L
		S[i]=Cs_end [CationInd[i]]*exp(-decay*tstep)/(molwt*1000)/1000; //mg/kg * mol/mg * kg/g =mol/g
		
		Cneg[i]=min(C[i],1e-20);//zero concentrations may cause problems with jacobian
		Sneg[i]=min(S[i],1e-20);
		C[i]-=Cneg[i];
		S[i]-=Sneg[i];

		T[i]=C[i]+S[i]*pb/n; //in mol/L

		available+=valence[i]*T[i]*n/pb;  //total available equivalents per gram for sorption	
	}

	if (available<Q)
	{
		//RxNExitGracefully("CCationExchange::React: Not enough cation mass to sorb to all surface sites",RUNTIME_ERR);
		//optional-assume all mass is instantaneously sorbed
		for (i=0; i<nCations;i++){ 
		  change=C[i]*n/pb;
			C[i]-=change*pb/n;	
			S[i]+=change;

			if (fabs(C[i]+S[i]*pb/n-T[i])>REALSMALL){
				RxNExitGracefully("CCationExchange::React: Internal mass balance error: implementation mistake(1)",RUNTIME_ERR);}
		}	
	}
	else{

	//cout <<"r";
	/*cout <<"Q: "<<Q<<endl;
	for (i=0; i<nCations; i++){
		cout <<"C["<<i<<"]: "<<C[i]<<"      S["<<i<<"]: "<<S[i]<<"      K["<<i<<"]: "<<SelCo[i]<<endl;
	}*/

	//-----------------------------------------------------------------------------
	/*Challenging:
	 If soil is not filled to capacity or overfilled (sum of v*S != Q), then problem is numerically unsolvable
	 To correct for this problem, if there is too little sorbed mass,
	   (a) assume that lower K species instantly equilibrate with soil until aqueous mass is zero
	   (b) if soil still not filled to capacity, then move to next species
	   (c) if soil *still* not filled to capcity after all aqueous mass is removed, then kill the simulation
	 To correct for this problem, if there is too much sorbed mass,
	   (a) assume that lower K species instantly equilibrate with soil until sorbed mass is zero
	   (b) if soil still filled over capacity, then move to next species
  */

	Sorbed_to_sites=0.0;
	for (i=0; i<nCations;i++){
		Sorbed_to_sites+=valence[i]*S[i];
	}
	diff=Q-Sorbed_to_sites;

	if (fabs(diff*pb/n)>REALSMALL)//Should use some kind of tolerance (diff/Q<~0.0001?)
	{
		if (diff>0.0)//Too little sorbed mass (diff>0)
		{ 
			for (i=0; i<nCations;i++)//TMP DEBUG: for now, assumes species are ordered by selectivity coeff
			{ 
				change=min(C[i]*n/pb,diff/valence[i]);
				C[i]-=change*pb/n;	
				S[i]+=change;
				diff-=change*valence[i];

				if (fabs(C[i]+S[i]*pb/n-T[i])>REALSMALL){
					RxNExitGracefully("CCationExchange::React: Internal mass balance error: implementation mistake(1)",RUNTIME_ERR);}
			}		
		}
		else if (diff<0.0)//Too much sorbed mass (diff<0)
		{ 
			for (i=0; i<nCations;i++){//TMP DEBUG: for now, assumes species are ordered by selectivity coeff
				
				change=min(S[i],-diff/valence[i]);
				S[i]-=change;
				C[i]+=change*pb/n;
				diff+=change*valence[i];

				if (fabs(C[i]+S[i]*pb/n-T[i])>REALSMALL){
					RxNExitGracefully("CCationExchange::React: Internal mass balance error: implementation mistake(2)",RUNTIME_ERR);}
			}			
		}
		else //WAY too much sorbed mass (more than 50% over capcity)
		{ 
			RxNExitGracefully("CCationExchange::React: Sorbed well past capacity!!!",RUNTIME_ERR);
		}
		for (i=0; i<nCations;i++)//zero concentrations may cause problems with jacobian
		{ 
			upperswap(C[i],1e-20); 
			upperswap(S[i],1e-20);
		}	
	}
	//-----------------------------------------------------------------------------

	while ((iter<CAT_EX_MAXITER) && (!done))
	{

		//Calculate Activity Coefficients 
		CalculateActivityCoeff(false);
		
		//Build System of equations
		BuildJacobian(J,B,pb,n);

		//Solve System of Equations
		if (!Gauss(J,B,dC,nCations*2))
		{
			/*cout << "Available: "<<available<<endl;
			cout << "Q:         "<<Q<< endl;
      Sorbed_to_sites=0.0;
		  for (i=0;i<nCations; i++){
			  cout << "C["<<i<<"]: "<<C[i]<<" S["<<i<<"]: "<<S[i]<<endl;
		    Sorbed_to_sites+=valence[i]*S[i]; //pb/n used to convert error to mol/L
			
			}
			cout << "Sorbed to sites:         "<<Sorbed_to_sites<< endl;
		  for (i=0; i<nCations; i++){
				cout << "K["<<i<<"    :"<<SelCo    [i]<<endl;
			}*/
			RxNExitGracefully("CCationExchange::React: Unable to solve system of equations",RUNTIME_ERR);
		}

		//Pick up solution
		relax=1.0;
		for (i=0;i<nCations; i++){
			if ((C[i]+dC[i]         )<0.0){lowerswap(relax,fabs(0.95*C[i]/dC[i         ]));}
      if ((S[i]+dC[i+nCations])<0.0){lowerswap(relax,fabs(0.95*S[i]/dC[i+nCations]));}
		}

		for (i=0;i<nCations; i++){
			C[i]+=relax*dC[i]; 
			S[i]+=relax*dC[i+nCations];
			//cout << "C["<<i<<"]: "<<C[i]<<" S["<<i<<"]: "<<S[i]<<endl;
		}

    //Evaluate relative change in concentration
    rel_error=0.0;
	  for (i=0;i<nCations; i++){
			k=i+nCations;
			rel_error+=fabs(dC[i]/(relax*C[i]));
			rel_error+=fabs(dC[k]/(relax*S[i]));
		}		
	
    if (rel_error<=CAT_EX_TOLERANCE){done=true;}

		iter++;

    //cout <<"iter: "<<iter<<" rel_error:"<<rel_error<<" relax: "<<relax<<endl;
	}
	}//if available <Q
	//cout <<"Error: "<<Validate(pb,n)<<endl;

	//Add back any negative aqueous concentrations to preserve mass balance
	//Convert back to mass concentration (currently not weighted) 
	for (i=0; i<nCations; i++){
		C[i]+=Cneg[i];
		S[i]+=Sneg[i];
    Caq_new[CationInd[i]]=C[i]*(pSpeciesArray[CationInd[i]]->GetMolecularWeight()*1000);       //mol/L *mg/mol=mg/L
		Cs_new [CationInd[i]]=S[i]*(pSpeciesArray[CationInd[i]]->GetMolecularWeight()*1000)*1000;  //mol/g *mg/mol*g/kg=mg/kg
	}




}
/***************************************************************
   BUILD JACOBIAN
**************************************************************
 Builds Jacobian Matrix and RHS for system of equations
 Based upon a generic version of MOUSER's eqns subroutine
------------------------------------------------------------*/
void CCationExchange::BuildJacobian(Writeable2DArray  A, 
																		Writeable1DArray  B, 
																		const double     &pb, //kg/L
																		const double     &n){ //-

	int    i,j,k,s;
	
	static double d_mono_mono; //first derivative of a MONOvalent act. coeff w.r.t. a MONOvalent ion
	static double d_mono_di;   //first derivative of a MONOvalent act. coeff w.r.t. a DIvalent   ion
	static double d_di_mono;   //first derivative of a DIvalent   act. coeff w.r.t. a MONOvalent ion
	static double d_di_di;     //first derivative of a DIvalent   act. coeff w.r.t. a DIvalent   ion
	static double d_sq_mono;   //first derivative of the square of a monovalent act. coeff w.r.t. a MONOvalent ion
	static double d_sq_di;     //first derivative of the square of a monovalent act. coeff w.r.t. a DIvalent ion

  static double SRA  = pow(IonicStrength,0.5);
  static double SRAM = SRA*(1.0+SRA)*(1.0+SRA);

	d_mono_mono = MonoActivity             *LN_10*(SRAM-5.0)/(40.0*SRAM);
  d_mono_di   = MonoActivity             *LN_10*(SRAM-5.0)/(10.0*SRAM);
  d_di_mono   = DiActivity               *LN_10*(SRAM-5.0)/(10.0*SRAM);
  d_di_di     = DiActivity               *LN_10*(SRAM-5.0)/( 2.5*SRAM);
  d_sq_mono   = MonoActivity*MonoActivity*LN_10*(SRAM-5.0)/(20.0*SRAM);
  d_sq_di     = MonoActivity*MonoActivity*LN_10*(SRAM-5.0)/( 5.0*SRAM);

	//Equilibrium Equations************************************************************
	//one for each cation except reference cation
	// (i=0 to N-2)
	for (i=0; i<nCations-1; i++){  
    s=i+1; //Equilibrium equation i deals with aqueous species i+1 (there is no equilibrium equation for reference species)
		
		//Derivative with regard to aqueous reference species-------------
		j=0;
		if      (valence[s]==1){ 
			A[i][j]=d_mono_mono*C[s]-(      S[s]/(S[0]     *SelCo[s]))*(d_mono_mono*C[0]+MonoActivity); //pg 90 Bandilla
		}
		else if (valence[s]==2){
			A[i][j]=d_di_mono  *C[s]-(2.0*Q*S[s]/(S[0]*S[0]*SelCo[s]))*(d_sq_mono*C[0]*C[0]+2.0*MonoActivity*MonoActivity*C[0]);//pg 91 Bandilla
		}

		//Derivatives with regard to aqueous species (dEQ_i/dCj)------
		for (j=1; j<nCations; j++){
			k=j;
			if      (valence[s]==1){ 
				if      (valence[k]==1){A[i][j]=d_mono_mono*(C[s]-C[0]*S[s]/(S[0]*SelCo[s]));}//monovalent/monovalent
				else if (valence[k]==2){A[i][j]=d_mono_di  *(C[s]-C[0]*S[s]/(S[0]*SelCo[s]));}//monovalent/divalent
			}
			else if (valence[s]==2){
				if      (valence[k]==1){A[i][j]=d_di_mono*C[s]-2.0*Q*d_sq_mono*C[0]*C[0]*S[s]/(S[0]*S[0]*SelCo[s]);}//divalent/monovalent
				else if (valence[k]==2){A[i][j]=d_di_di  *C[s]-2.0*Q*d_sq_di  *C[0]*C[0]*S[s]/(S[0]*S[0]*SelCo[s]);}//divalent/divalent
			}
			if (s==k){
				if      (valence[k]==1){A[i][j]+=MonoActivity;}
				else if (valence[k]==2){A[i][j]+=DiActivity;  }
			}
		}

		//Derivative with regard to sorbed reference species (dEQ_i/dSo)--
    j=nCations;
		k=0;
		if      (valence[s]==1){ 
			A[i][j]=MonoActivity*C[0]*S[s]/(S[0]*S[0]*SelCo[s]); //pg 90 Bandilla thesis 
		}
		else if (valence[s]==2){
			A[i][j]=4.0*Q*MonoActivity*MonoActivity*C[0]*C[0]*S[s]/(S[0]*S[0]*S[0]*SelCo[s]); //pg 91 Bandilla thesis 
		}

    //Derivatives with regard to sorbed species (dEQ_i/dSk)-------
		for (j=nCations+1; j<2*nCations; j++){
			k=j-nCations; 
			if (s==k){
				if      (valence[s]==1){ 
					A[i][j]=-MonoActivity*C[0]/(S[0]*SelCo[s]); //pg 90 Bandilla thesis 
				}
				else if (valence[s]==2){
					A[i][j]=-2.0*Q*MonoActivity*MonoActivity*C[0]*C[0]/(S[0]*S[0]*SelCo[s]);//pg 91 Bandilla thesis 
				}
			}
			else{
				A[i][j]=0.0;
			}
		}
		//RHS---------------------------------------------------------
    if      (valence[s]==1){
			B[i]=(MonoActivity*C[0]*S[s]/(S[0]*SelCo[s])-MonoActivity*C[s]); 
		}
		else if (valence[s]==2){
		  B[i]=(2.0*Q*MonoActivity*MonoActivity*C[0]*C[0]*S[s]/(S[0]*S[0]*SelCo[s])-DiActivity*C[s]);
		}
	}

	//Mass Balance Equations***********************************************************
	//One for each cation, including reference Cation
	// (i=N-1 to N=2N-2)
	for (i=nCations-1; i<2*nCations-1; i++){
		k=i-nCations+1;
		for (j=0; j<nCations; j++){
			if (j==k)         {A[i][j]=1.0;}
			else              {A[i][j]=0.0;}
		}
		for (j=nCations; j<2*nCations; j++){
			if (j==i+1)       {A[i][j]=pb/n;}
			else              {A[i][j]=0.0;}		 
	  }		 
	  B[i]=T[k]-C[k]-S[k]*pb/n;
	}

	//Sorbed Mass Balance Equation*****************************************************
	// last equation (i=2N-1)
	i=2*nCations-1;
 	for (j=0; j<nCations; j++){
	  A[i][j]=0.0;
	}
	for (j=nCations; j<2*nCations; j++){
		k=j-nCations;
		A[i][j]=valence[k];
	}
	B[i]=Q;
	for (k=0;k<nCations; k++){
	  B[i]-=valence[k]*S[k];   
	}

  //Print Matrix (TMP DEBUG)*******************************************************
	/*cout <<endl<< "Jacobian matrix"<<endl;
  for(s=0; s<2*nCations; s++){ 
		cout<<'|';  
		for(i=0; i<2*nCations; i++){ 
	     if(fabs(A[s][i])>REALSMALL){cout.width(6);cout<<A[s][i]<<" ,";} 
       else                       {cout.width(7);cout<<0.0    <<" ,";}
		} 
	  cout<<" | "<<B[s]<<endl;
	}  
  cout << "B vector"<<endl; for(i=0; i<2*nCations; i++){cout <<B[i]<<" ";} 
  cout <<endl;*/
}
/***************************************************************
   VALIDATE
**************************************************************
 Check Equilibrium & Mass Balance Equations
 returns error in units of aqueous molar concentrations

at this point C=mol/L and S=mol/g and pb=kg/L and Q=eq/g
------------------------------------------------------------*/
double CCationExchange::Validate(const double &pb, const double &n) const{
	int    i,s;
	double error(0.0); //units of aqueous molar concentration [mol/L]
	double MBerror;

  //Check Equilibrium Equations--------------------------
	//one for each cation except reference cation
	for (i=0; i<nCations-1; i++){
    s=i+1; //Equilibrium equation i deals with species i+1 (there is no equilibrium equation for reference species)
		
		if (valence[s]==1){
			MBerror=MonoActivity*(C[s]-C[0]*S[s]/(S[0]*SelCo[s]));
			error  +=MBerror;

			cout <<"EQ error "<<s <<": " <<MBerror<<endl;
		}

		else if (valence[s]==2){
			MBerror=C[s]-2.0*Q*MonoActivity*MonoActivity*C[0]*C[0]*S[s]/(DiActivity*S[0]*S[0]*SelCo[s]);
      error +=MBerror;

			cout <<"EQ error "<<s <<": " <<MBerror<<endl;
		}
	}
  //Check Mass Balance Equations-------------------------
	for (i=0; i<nCations;i++){
		MBerror=C[i]+S[i]*pb/n-T[i];
		error+=MBerror;

    cout <<"MB error "<<i+1 <<": " <<MBerror<<endl;
	}

	//Check Capcity Mass Balance----------------------------
	MBerror=-Q*pb/n;
	for (i=0; i<nCations;i++){
		MBerror+=valence[i]*S[i]*pb/n; //pb/n used to convert error to mol/L
   cout <<"Q/Qn"<< Q<<" "<<(MBerror+Q*pb/n)*n/pb<<endl;
	}
  
  error+=MBerror;

	cout <<"MB error (Q): " <<MBerror<<endl;

	return error;
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
	string "CationExchangeRxN"
	double capacity double anionconc int soiltype
	{int index int valence double selectivity}x(NumCations)
  & 
	First selectivity selco[0] is not used, same with valence 
----------------------------------------------------------------------*/
CCationExchange *CCationExchange::Parse(ifstream &input, int &l){


	CCationExchange  *pReaction=NULL;
  int      soiltype_ind;
  bool     done(false);
  double   thisQ(0),thisAA(0);
	int      Len,N(0);
  int      indices[MAX_SPECIES];
	int      valence[MAX_SPECIES];
  double   coeff  [MAX_SPECIES];
	char    *s      [MAXINPUTITEMS];

	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if (Len==3) {thisQ =s_to_d(s[0]); thisAA=s_to_d(s[1]);soiltype_ind=s_to_i(s[2]);}
	else        {ImproperFormat(s,l); return NULL;}

  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	done=false;
  do {
    if (N>=MAX_SPECIES) { RxNExitGracefully("CCationExchange::Parse- too many species in reaction data",TOO_MANY);}
		if  (Len==3){
      indices[N]=s_to_i(s[0]); 
      valence[N]=s_to_i(s[1]);
			coeff  [N]=s_to_d(s[2]);
			N++;                                               
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			pReaction = new CCationExchange(N,indices,coeff,valence,thisQ,thisAA,soiltype_ind); 
			done=true;
		}
    else             {ImproperFormat(s,l); return NULL;}
	} while (!done);

	return pReaction;
}
/*************************************************************
Test Module
**************************************************************
------------------------------------------------------------*/
void TestCationExchange(){
	CCationExchange *pRxN;

	CSpecies   *pSArray[6];
	CSoilType  *pSoilarray[1];
	int    N=6; //TMP DEBUG - no strontium 90
	int    s;
	int    ind[6],v[6];
	double select[6];        //mol/L (divalent) or no units (monovalent)
	double cap    =1.143895; //eq/kg= 1.144 meq/g CORRECT!!!
	double Aconc  =0.003504; //0.003504 eq/L ~ mol/L CORRECT!!!

	double Cold[6];
	double Cnew[6];
	double Sold[6];
	double Snew[6];
	double Cinit[6];
	double Sinit[6];

	for (s=0;s<N;s++){
		ind [s]=s;
		Cnew[s]=0.0;
		Snew[s]=0.0;
	}  
	
	cout <<"Setting Test Variables..."<<endl;

	double pb=0.8282;//kg/L
	double n=0.6056;		

	//selectivity (dimless or mol/L), Cold (mg/L), Cnew (mg/L)
	// mol/L (for v=2), dimless (for v=1)                          g/mol      mg/L                mg/kg 
	select[0]=0;    v[0]=1; pSArray[0]=new CSpecies("Na"  ,0.0,23  ,0.0); Cold[0]=216;        Sold[0]= 6000;   
  select[1]=28;   v[1]=1; pSArray[1]=new CSpecies("K"   ,0.0,39  ,0.0); Cold[1]=3.5;        Sold[1]=25000;
  select[2]=0.30; v[2]=2; pSArray[2]=new CSpecies("Mg"  ,0.0,24  ,0.0); Cold[2]=15;         Sold[2]=  120;
  select[3]=0.35; v[3]=2; pSArray[3]=new CSpecies("Ca"  ,0.0,40  ,0.0); Cold[3]=100;        Sold[3]= 4600;
  select[4]=4.80; v[4]=2; pSArray[4]=new CSpecies("Sr"  ,0.0,87.6,0.0); Cold[4]=0.98;       Sold[4]=   87.6;
  select[5]=4.80; v[5]=2; pSArray[5]=new CSpecies("Sr90",0.0,90  ,0.0); Cold[5]=0.00000537; Sold[5]=0.0;    //decay=6.6e-5/day

	pSoilarray[0]=new CSoilType("Zeolite", pb);

	//if sodium the only initial sorbed cation present
	/*Sold[0]= 1.143895*1000*23.0;//eq/kg* mg/g * g/mol = mg/kg ~ 26 grams per kilogram (saltwater?)
	for (s=1;s<N;s++){
		Sold[s]=0.0;
	}  */
	

	for (s=0;s<N;s++){
		Cinit[s]=Cold[s];
		Sinit[s]=Sold[s];
	}  

	cout <<"Constructing..."<<endl;
	pRxN=new CCationExchange(N,ind,select,v,cap,Aconc,0);

	cout <<"Setting Soiltype Array..."<<endl;
	pRxN->SetSoilArray(pSoilarray,1);

	cout <<"Setting Species Array..."<<endl;
	pRxN->SetSpeciesArray(pSArray,N);
	
	cout <<"Initializing..."<<endl;
  pRxN->Initialize();

	cout <<"Reacting..."<<endl;
	clock_t t1,t2;
	t1=clock();
	for (int i=0;i<1; i++){		
		pRxN->React(Cold,Cold,Cnew,Sold,Sold,Snew,0,n,1.0);
	}
	for (s=0;s<N;s++){
		cout <<" Species: "<<pSArray[s]->GetName() <<endl;
		cout <<"    aqueous[mg/L] : "<<Cinit[s]               <<"-------->"<<Cnew[s]             <<endl;
		cout <<"    sorbed [mg/L] : "<<Sinit[s]*pb/n          <<"-------->"<<Snew[s]*pb/n        <<endl;
		cout <<"     sorbed[mg/kg]: "<<Sinit[s]               <<"-------->"<<Snew[s]             <<endl;
		cout <<"    total  [mg/L] : "<<Cinit[s] +Sinit[s]*pb/n<<"-------->"<<Cnew[s]+Snew[s]*pb/n<<endl;
		cout <<"           %sorbed: "<<Sinit[s]*pb/n /(Cinit[s] +Sinit[s]*pb/n)*100<<"-------->";
		cout                         <<Snew [s]*pb/n /(Cnew [s] +Snew [s]*pb/n)*100 <<endl;
    cout <<"    Selco           "<<select[s]<<endl;
	}

	t2=clock();
	cout <<" ..."<<float(t2-t1)/CLK_TCK << " seconds elapsed (Cation Exchange). " << endl;
  //For 4 species, 5 iterations, 100,000 runs: ~12.1 seconds
  //For 5 species, 5 iterations, 100,000 runs: ~16.7 seconds
  //For 6 species, 5 iterations, 100,000 runs: ~22.1 seconds

	for (s=0;s<N;s++){Snew[s]=Sold[s];}Cnew[0]=Cold[0];

	pRxN->GetEquilibriumAqueous(Snew,Cnew,Cnew[0],0,n);

	cout<<"After Equilibrium Calulation:"<<endl;

	/*for (s=0;s<N;s++){
		
		cout <<" Species: "<<pSArray[s]->GetName() <<endl;
		cout <<"    aqueous[mg/L] : "<<Cinit[s]               <<"-------->"<<Cnew[s]             <<endl;
		cout <<"    sorbed [mg/L] : "<<Sinit[s]*pb/n          <<"-------->"<<Snew[s]*pb/n        <<endl;
		cout <<"     sorbed[mg/kg]: "<<Sinit[s]               <<"-------->"<<Snew[s]             <<endl;
		cout <<"    total  [mg/L] : "<<Cinit[s] +Sinit[s]*pb/n<<"-------->"<<Cnew[s]+Snew[s]*pb/n<<endl;
		cout <<"           %sorbed: "<<Sinit[s]*pb/n /(Cinit[s] +Sinit[s]*pb/n)*100<<"-------->";
		cout                         <<Snew [s]*pb/n /(Cnew [s] +Snew [s]*pb/n)*100 <<endl;
    cout <<"    Selco           "<<select[s]<<endl;
		
	}*/
	for (s=0;s<N;s++){
		delete pSArray[s];
	}

	delete pRxN;
	
}