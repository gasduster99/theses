//SpeciesAndSoils.h (includes soil & Species Data)
#ifndef SPECIES_SOIL_H
#define SPECIES_SOIL_H

#include "RxNLibraryInclude.h"

const int     MAX_SPECIES      =7;        // maximum number of species

/****************************************************
 *  Class CSpecies
 *  Layer Chemical Species Data Abstraction
 ***************************************************/
class CSpecies{
 private:/*----------------------------------------*/
  static int  TotalSpecies;

	//Member Variables
  char       *name;							//Species Name
	int         index;            
	double      Diff;             //species-specific diffusion coeff
	double      molweight;        //molecular weight of species
	double      decay;            //decay rate (1/T)

 public:/*-----------------------------------------*/
	//Constructors
	CSpecies();
	CSpecies(char *Name, double D,double mol, double decay_rate);
 ~CSpecies();

	//static Accessors
//  static int GetNumSpecies ();

	//Accessor functions
	char   *GetName           () const;
	double  GetDiffusionCoeff () const;
  double  GetMolecularWeight() const;
	double  GetDecayRate      () const;
};
/****************************************************
 *  Class CSoilType
 *  Layer Chemical Soil Property Data Abstraction
 ***************************************************/
class CSoilType{
 private:/*----------------------------------------*/
  static int  TotalSoils;

	//Member Variables
  char       *name;							//Soiltype Name
	double      pb;								//Bulk Dry Density [kg/L]
	int         index;            //soil index

 public:/*-----------------------------------------*/
	//Constructors
	CSoilType();
	CSoilType(char *Name, const double density);//[kg/L]
 ~CSoilType();

	//static Accessors
  static int GetNumSoils ();

	//Accessor functions
	char   *GetName       () const;
	double  GetBulkDryDensity() const;
};
//================================================================================================
//================================================================================================


#endif