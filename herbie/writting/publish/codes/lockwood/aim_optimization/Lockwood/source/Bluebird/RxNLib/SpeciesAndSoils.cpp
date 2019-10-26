#include "SpeciesAndSoils.h"

/************************************************************************
           SPECIES CLASS
************************************************************************/
CSpecies::CSpecies(){}//should not be called
//------------------------------------------------------------------------------
CSpecies::CSpecies(char *Name, double D, double molwt, double decay_rate){
  name     =Name;
	Diff     =D;
	molweight=molwt;
	decay    =decay_rate;
	index    =TotalSpecies;
	TotalSpecies++;
}
//------------------------------------------------------------------------------
CSpecies::~CSpecies(){}
//------------------------------------------------------------------------------
int    CSpecies::TotalSpecies=0;
//------------------------------------------------------------------------------
//int    CSpecies::GetNumSpecies     ()      {return TotalSpecies;}
//------------------------------------------------------------------------------
char*  CSpecies::GetName           () const{return name;}
//------------------------------------------------------------------------------
double CSpecies::GetDiffusionCoeff () const{return Diff;};
//------------------------------------------------------------------------------
double CSpecies::GetMolecularWeight() const{return molweight;};
//------------------------------------------------------------------------------
double CSpecies::GetDecayRate      () const{return decay;}

/************************************************************************
           SOILTYPE CLASS
************************************************************************/
CSoilType::CSoilType(){}//should not be called
//------------------------------------------------------------------------------
CSoilType::CSoilType(char *Name, double density){
  name=Name;
	//if (density<=0){ExitGracefully("CSoilType:Bad bulk dry density",BAD_DATA);}
	pb=density;
	index=TotalSoils;
	TotalSoils++;
}
//------------------------------------------------------------------------------
CSoilType::~CSoilType(){}
//------------------------------------------------------------------------------
int    CSoilType::TotalSoils=0;
//------------------------------------------------------------------------------
int    CSoilType::GetNumSoils      ()      {return TotalSoils;}
//------------------------------------------------------------------------------
char*  CSoilType::GetName          () const{return name;}
//------------------------------------------------------------------------------
double CSoilType::GetBulkDryDensity() const{return pb;}
