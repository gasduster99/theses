

#include "ReactionScheme.h"
#include "ParentDaughter.h"
/*************************************************************
 Constructors
**************************************************************
------------------------------------------------------------*/
CParentDaughter::CParentDaughter(){
	NS=0;
	Indices  =NULL;
	K        =NULL;
	y        =NULL;
	C        =NULL;
	S        =NULL;
	Cneg     =NULL;
	Sneg     =NULL;
}
//-----------------------------------------------------------
CParentDaughter::CParentDaughter(int     NumSpecies, 
																 int    *indices, 
																 double *decay_coeff,
																 double *yield_coeff, 
																 int     soiltype_index):
                 CReactionScheme(false){ 

	//must be input so that species[indices[i]] is parent of species[indices[i+1]]
	//species[indices[NumSpecies-1]] (i.e., Vinyl Chloride) has no children
	int i;

  NS       =NumSpecies;
	soiltype =soiltype_index;
	sorption_reaction=false;
	Indices  =new int    [NS];
	K        =new double [NS];
	y        =new double [NS];
	C        =new double [NS];
	S        =new double [NS];
	Cneg     =new double [NS];
	Sneg     =new double [NS];
	a        =new double [NS];
  if (a==NULL){RxNExitGracefully("CParentDaughter:Constructor",OUT_OF_MEMORY);}

  for (i=0; i<NS; i++){
		Indices  [i]=indices      [i];
    K        [i]=decay_coeff  [i];
    y        [i]=yield_coeff  [i];
		C        [i]=0.0;
		S        [i]=0.0;
		Cneg     [i]=0.0;
		Sneg     [i]=0.0;
		a				 [i]=0.0;
	}
}
//-----------------------------------------------------------
CParentDaughter::~CParentDaughter(){
	delete [] Indices;
	delete [] y;
	delete [] K;
	delete [] a;
	delete [] C;
	delete [] S;
}
/***************************************************************
 Initialize
**************************************************************
Initialize Reaction
------------------------------------------------------------*/
void CParentDaughter::Initialize(){

	if (NS>nspecies){
		RxNExitGracefully("CParentDaughter::Initialize: More parent-daughter products than existing species",BAD_DATA);}
	if (pSpeciesArray==NULL){
		RxNExitGracefully("CParentDaughter::Initialize: Species Array not set",BAD_DATA);}
	//Quality check
  for (int i=0; i<NS; i++){
		if ((Indices[i]>nspecies)|| (Indices[i]<0)){
			RxNExitGracefully("CParentDaughter::Initialize: Invalid species index",BAD_DATA);}
		if ((y[i]>1.0)|| (y[i]<0)){
			RxNExitGracefully("CParentDaughter::Initialize: Invalid parent daughter decay yield ratio",BAD_DATA);}
		if (K[i]<0.0){
			RxNExitGracefully("CParentDaughter::Initialize: Parent daughter decay rates must be positive",BAD_DATA);}
		for (int j=0; j<i; j++){//corrects for equal decay values (division by zero)
			if (K[i]==K[j]){
				K[i]+=K[i]/100000;
			}
		}
	}
}
/***************************************************************
 React
**************************************************************
combination of newton and Exchange routines in mouser
------------------------------------------------------------*/
void CParentDaughter::React(Ironclad1DArray  Caq,
														Ironclad1DArray  Caq_end,
														Writeable1DArray Caq_new,
														Ironclad1DArray  Cs,
														Ironclad1DArray  Cs_end,
														Writeable1DArray Cs_new,	
														const int        soil,
														const double    &n,  
														const double    &tstep){
  int    i,j,k,s;
	static double prod,sum;

  //To avoid problems with reaction solution, set negative aqueous concentrations to zero
  //Save negative values to preserve mass balance in transport solution
  for (s=0;s<nspecies;s++){
		Caq_new[s]=Caq_end[s];
		Cs_new [s]=Cs_end[s];
	}
	//cout <<"C: ";
	for (i=0;i<NS;i++){
		Cneg[i]=min(Caq_end[Indices[i]],0.0);
		Sneg[i]=min(Cs_end [Indices[i]],0.0);
		C[i]   =max(Caq_end[Indices[i]],0.0);
		S[i]   =max(Cs_end [Indices[i]],0.0);
		a[i]=0.0;
		//cout <<C[i]<<",    ";
	}
	//cout <<endl;

	//forward transformation to surrogate concentration a

	//Transform----------------------------
	for (i=0;i<NS;i++){
		sum=0.0;
		for (j=0;j<i; j++){
			prod=1.0;
			for (k=j;k<i; k++){
				prod*=(y[k]*K[k])/(K[k]-K[i]);
			}
			sum+=prod*C[j];
		}
		a[i]=sum+C[i];
	}

	//react--------------------------------
	for (i=0;i<NS;i++){
		a[i]*=exp(-K[i]*tstep);
	}

	//cout <<"Cend: ";	
	//Transform back-----------------------
	for (i=0;i<NS;i++){
		sum=0.0;
		for (j=0;j<i; j++){
			prod=1.0;
			for (k=j;k<i; k++){
				prod*=(y[k]*K[k])/(K[k]-K[i]);
			}
			sum+=prod*C[j];
		}
		C[i]=a[i]-sum;
		//cout <<C[i]<<",    ";
	}
	//cout <<endl;	

	//Add back any negative aqueous concentrations to preserve mass balance
	for (i=0;i<NS;i++){
		Caq_new[Indices[i]]=C[i]+Cneg[i];
		Cs_new [Indices[i]]=S[i]+Sneg[i];
	}		

}
/***************************************************************
   VALIDATE
**************************************************************
 Check Equilibrium & Mass Balance Equations
 returns error in units of aqueous molar concentrations
------------------------------------------------------------*/
double CParentDaughter::Validate(const double &pb, const double &n) const{
	double error(0.0); //units of aqueous molar concentration

	return error;

}

/***********************************************************************
                           PARSE
************************************************************************
Format:
	string "ParentDaughterRxN"
	{int soiltype}
	{int index double decay_coeff double yield_coeff}x(NumSpecies)
  & 
----------------------------------------------------------------------*/
CParentDaughter *CParentDaughter::Parse(ifstream &input, int &l){


	CParentDaughter  *pReaction=NULL;

  bool     done(false);
  double   thisQ(0),thisAA(0);
	int      Len,N(0);
	int      soilindex(0);
  int      indices[MAX_SPECIES];
  double   decay  [MAX_SPECIES];
	double   yld    [MAX_SPECIES];
	char    *s      [MAXINPUTITEMS];

	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	if (Len==1)      {soilindex=s_to_i(s[0]);}
	else             {ImproperFormat(s,l); return NULL;}

  if (TokenizeLine(input,s,Len)){return NULL;}; l++;
	done=false;
  do {
    if (N>=MAX_SPECIES) { RxNExitGracefully("CParentDaughter::Parse- too many species in reaction data",TOO_MANY);}
		if  (Len==3){
      indices[N]=s_to_i(s[0]); 
			decay  [N]=s_to_d(s[1]);
      yld    [N]=s_to_d(s[2]);
			N++;                                               
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			pReaction = new CParentDaughter(N,indices,decay,yld,soilindex); 
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
void TestParentDaughter(){
	CParentDaughter *pRxN;
	CSoilType  *pSoilarray[1];
	CSpecies *pSArray[4];
	int    N=4; 
	int    s;
	int    ind[4];
	double decay_coeff[4],frac[4];

	double Cold[4];
	double Cnew[4];
	double Sold[4];
	double Snew[4];
	double Cinit[4];
	double Sinit[4];

	for (s=0;s<N;s++){
		ind  [s]=s;
		Cnew [s]=0.0;
		Snew [s]=0.0;
	}  
	
	cout <<"Setting Test Variables..."<<endl;

	double pb=0.8;//kg/L
	double n=0.3;		

	//selectivity (mol/L), Cold (mg/L), Cnew (mg/L)
	// d^-1                                                     g/mol          mg/L                mg/kg 
	decay_coeff[0]=3.61/24; frac[0]=0.40; pSArray[0]=new CSpecies("TCE" ,0.0,0.0,0.0); Cold[0]=10000;   Sold[0]=0.0;   
  decay_coeff[1]=5.73/24; frac[1]=0.02; pSArray[1]=new CSpecies("PCE" ,0.0,0.0,0.0); Cold[1]=1000;		Sold[1]=0.0;
  decay_coeff[2]=2.97/24; frac[2]=0.01; pSArray[2]=new CSpecies("DCE" ,0.0,0.0,0.0); Cold[2]=100;			Sold[2]=0.0;
  decay_coeff[3]=3.61/24; frac[3]=0.00; pSArray[3]=new CSpecies("VC"  ,0.0,0.0,0.0); Cold[3]=10;			Sold[3]=0.0;

	pSoilarray[0]=new CSoilType("Zero-valent Iron", pb);

	for (s=0;s<N;s++){
		Cinit[s]=Cold[s];
		Sinit[s]=Sold[s];
	}  


	//if sorbed K concentration is too low, negative concentrations resulted
	cout <<"Constructing..."<<endl;
	pRxN=new CParentDaughter(N,ind,decay_coeff,frac,0);

	cout <<"Setting Species Array..."<<endl;
	pRxN->SetSpeciesArray(pSArray,N);
	
	cout <<"Setting Soiltype Array..."<<endl;
	pRxN->SetSoilArray(pSoilarray,1);

	cout <<"Initializing..."<<endl;
  pRxN->Initialize();

	cout <<"Reacting..."<<endl;
	ofstream PDDECAY;
	PDDECAY.open("PDdecay.csv");
	PDDECAY<<"t,TCE,PCE,DCE,VC"<<endl;
	for (double t=0.0; t<60; t+=0.25){
		pRxN->React(Cold,Cold,Cnew,Sold,Sold,Snew,0,n,0.25);
		PDDECAY<<t<<","<<Cnew[0]<<","<<Cnew[1]<<","<<Cnew[2]<<","<<Cnew[3]<<endl;
		for (int i=0;i<4;i++){Cold[i]=Cnew[i]; Sold[i]=Snew[i];}
	}
	PDDECAY.close();

	for (s=0;s<N;s++){
		cout <<" for one liter:"<<endl;
		cout <<"         aqueous[mg/L]: "<<Cnew[s]              <<endl;
		cout <<"          sorbed[mg/L]: "<<Snew[s]*pb/n         <<endl;
		cout <<"  initial total[mg/kg]: "<<Cinit[s] +Sinit[s]*pb/n<<endl;
		cout <<"           total[mg/L]: "<<Cnew[s] +Snew[s]*pb/n<<endl;
		cout <<"         sorbed[mg/kg]: "<<Snew[s]              <<endl;
		cout <<"               %sorbed: "<<Snew[s]*pb/n /(Cnew[s] +Snew[s]*pb/n)*100 <<endl;
		delete pSArray[s];
	}

	delete pRxN;
	
}