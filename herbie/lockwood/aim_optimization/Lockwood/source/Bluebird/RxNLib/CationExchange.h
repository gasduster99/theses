//ReactionScheme.h
#ifndef CATEXCHANGE_H
#define CATEXCHANGE_H

#include "ReactionScheme.h"

const double CAT_EX_TOLERANCE=1e-6; //From Mouser
const int    CAT_EX_MAXITER  =100;  //From Mouser

/**************************************************************
	Class CCationExchange
	Data Abstraction for Cation Exchange module
**************************************************************/
class CCationExchange:public CReactionScheme{
 private:/*--------------------------------------------------*/
	int     *CationInd;        //array of indices [nCations] of cation species->reference species must hold first place
	int     *valence;          //array of valences [nCations] of cation species->reference species must hold first place
	int      nCations;         //number of cations

	int      soiltype;         //global index of applicable soil

	double   Q;                //cation exchange capacity [eq/g]

	double  *SelCo;            //Array of selectivity coefficients (against reference species)[mol/L]
	                           //[1..nCations] (first place-reference not used)

	double   AnionConc;        //concentration of ions in solution [eq/L]
	double   MonoActivity;     //monovalent activity coefficient [-]
  double   DiActivity;       //divalent activity coefficient [-]
	double   IonicStrength;    //Ionic strength of aqueous solution [mol/L]

	double  *C;                //molar Caq [nCations]->reference species must hold first place
	double  *S;                //molar Csorbed [nCations]->reference species must hold first place
	double  *T;                //total mass of each Cation Species [nCations]->reference species must hold first place
  double  *Cneg;
	double  *Sneg;

  double **J;                //Jacobian Matrix 
	double  *B;                //Right Hand Side Vector
	double  *dC;

	void   CalculateActivityCoeff(bool estimateonly);
	void                Translate(Ironclad1DArray Caq, Ironclad1DArray  Cs, const double &pb, const double &n);
	void            BuildJacobian(Writeable2DArray A,  Writeable1DArray B, const double &pb, const double &n);
	double               Validate(const double &pb, const double &n) const;

 public:/*---------------------------------------------------*/
	//Constructors
	CCationExchange();
  CCationExchange::CCationExchange(int     NumCations, 
																   int    *indices, 
																   double *selectivities, //[mol/L]
																   int    *valences,
																   double  capacity,      //[meq/g] ~ [eq/kg]
																   double  AnionicConc,   //[mol/L] ~ [eq/L]
																	 int     soiltype_index); //index of appropriate soil
 ~CCationExchange();

  static CCationExchange *CCationExchange::Parse(ifstream &input, int &l);
  
	void Initialize();
	void React(Ironclad1DArray  Caq,      //all concentrations in [mg/L], for now
						 Ironclad1DArray  Caq_end,
						 Writeable1DArray Caq_new,
		         Ironclad1DArray  Cs,       //all sorbed concentrations in [mg/kg], for now
						 Ironclad1DArray  Cs_end,
						 Writeable1DArray Cs_new,	
						 const int        soil,
						 const double    &n,  
						 const double    &tstep);   //local time coordinates

	void GetEquilibriumAqueous(Ironclad1DArray  Cs,
		                         Writeable1DArray Caq_new,
														 const double    &Caq_ref,
		 						             const int        soil,
						                 const double    &n); 

};

void TestCationExchange();

#endif