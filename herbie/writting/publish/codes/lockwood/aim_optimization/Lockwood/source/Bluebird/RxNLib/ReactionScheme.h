//ReactionScheme.h
#ifndef RXNSCHEME_H
#define RXNSCHEME_H

//#include "CardinalInclude.h"
//#include "MasterInclude.h"
#include "RxNLibraryInclude.h"
#include "SpeciesAndSoils.h"

/**************************************************************
	Class CReactionScheme
	Data Abstraction for all batch reaction modules
**************************************************************/
class CReactionScheme{
 protected:/*------------------------------------------------*/

  CSpecies  **pSpeciesArray; //array of pointers to species
	CSoilType **pSoilArray;    //array of pointers to soils
  
	int         nsoils;
	int         nspecies;

	double      localtimestep;
	bool        sorption_reaction;

 public:/*---------------------------------------------------*/
	CReactionScheme(){
		pSpeciesArray=NULL;
		localtimestep=ALMOST_INF;
		sorption_reaction=false;
	}
	CReactionScheme(const bool sorptionRxn){
		pSpeciesArray=NULL;
		localtimestep=ALMOST_INF;
		sorption_reaction=sorptionRxn;
	}

	virtual ~CReactionScheme(){
		 if (globaldebug){cout<<"DESTROYING ABSTRACT REACTION SCHEME"<<endl;}
		 delete [] pSpeciesArray;
		 delete [] pSoilArray;
	 }

	void SetSpeciesArray(CSpecies **spec, const int nspec){
		if (spec==NULL){
			RxNExitGracefully("CReactionScheme::SetSpeciesArray: NULL species",RUNTIME_ERR);}
		if ((nspec<0) || (nspec>MAX_SPECIES)){
			RxNExitGracefully("CReactionScheme::SetSpeciesArray: invalid number of  species",RUNTIME_ERR);}
		nspecies=nspec;
		pSpeciesArray=new CSpecies *[nspecies];
		for (int s=0; s<nspecies; s++){
			pSpeciesArray[s]=spec[s];//just holds pointer to species
		}
	}
	void SetSoilArray(CSoilType **soil, const int nsoil){
		if (soil==NULL){
			RxNExitGracefully("CReactionScheme::SetSoilArray: NULL species",RUNTIME_ERR);}
		nsoils=nsoil;
		pSoilArray=new CSoilType *[nsoils];
		for (int s=0; s<nsoils; s++){
			pSoilArray[s]=soil[s]; //just holds pointer to soils
		}
	}
	virtual bool IsSorptionReaction(){return sorption_reaction;}

  virtual void Initialize()=0;
  virtual void React(Ironclad1DArray  Caq,
										 Ironclad1DArray  Caq_end,
                     Writeable1DArray Caq_new,
										 Ironclad1DArray  Cs,
										 Ironclad1DArray  Cs_end,
                     Writeable1DArray Cs_new,										 
										 const int        soil,
										 const double    &n,  
										 const double    &tstep)=0;
	
	virtual double GetRetardation(const double  pb,
																const double  n, 
																const double *conc, 
																const int			s){return 1.0;}

};
/**************************************************************
	Class CSimpleDecay
	Data Abstraction for Simple Decay module
**************************************************************/
class CSimpleDecay:public CReactionScheme{
 private:/*--------------------------------------------------*/
  double *DecayRates;

 public:/*---------------------------------------------------*/
	CSimpleDecay();
	CSimpleDecay(double *decayrates);
 ~CSimpleDecay();

	void Initialize(){}
	void React(Ironclad1DArray  Caq,
						 Ironclad1DArray  Caq_end,
						 Writeable1DArray Caq_new,
		         Ironclad1DArray  Cs,
						 Ironclad1DArray  Cs_end,
						 Writeable1DArray Cs_new,	
						 const int        soil,
						 const double    &n,  
						 const double    &tstep);
};

#endif