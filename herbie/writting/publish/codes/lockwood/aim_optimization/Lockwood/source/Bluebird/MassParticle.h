//MassParticle.h
#ifndef MASSPARTICLE_H
#define MASSPARTICLE_H

#include "MasterInclude.h"
#include "BluebirdLibrary.h"
#include "Particle.h"
#include "AbstractLayer.h"
#include "AbstractDomain.h"  

class CParticle;
class CPathline;
/****************************************************
 *  Class CMassParticle
 *  Contaminant Particle Tracking Data Abstraction
 *  Strong for random walk
 ***************************************************/
class CMassParticle:public CParticle{
 private:/*----------------------------------------------------------*/
	//Member Variables
	double  concentration;          //array of current concentrations or masses for different species    
	int     s;							        //index of species

	const   C2DDomainABC *pDomain;	//pointer to chemical domain //TMP DEBUG- should be 3D

	pt3D    lastpt[3];							//previous placement of particle
	double  lastt[3];								//previous times of particle
	pt3D    currpt;                 //current location
	double  currt;                  //current time
	bool    captured;               //true if captured, false otherwise
  
	//Member Functions (inherited from particle for tracking)															
	void   Advect          (const vector &movement, const double &tstep); //includes retardation factor (single species)
	void   Capture         (const double &endtime);
	bool	 HasBeenCaptured ();																											  
  void	 BackTrack       ();       

 public:/*-----------------------------------------------------------*/
	CMassParticle();
	CMassParticle(const CAquiferABC   *pAquifer,
								const C2DDomainABC  *pDom, //TMP DEBUG- should be 3D
								pt3D						     startpt, 
								trackdir					   direct, 
								double						   initconc,
								int                  species_ind);
 ~CMassParticle();

	//Accessor Functions (inherited from particle)
  bool   IsCaptured          () const;
	pt3D	 GetLocation         () const;                              
	double GetCurrentT         () const;  
	double GetLastT            () const;

	//Accessor Functions
	const double GetConcentration() const;
	const double GetMass() const;

  //Assignment Functions
  void	 SetMass             (const double mass);
  void	 SetConcentration    (const double conc);
	void	 SetLocation         (const pt3D &pt);
	
	//member functions 
  void	 Disperse            (const double &tstep, const double &t);
	void	 WriteOutput         () const;

	void   CMassParticle::WriteOutput			 (const double &t) const{}
  void	 CMassParticle::CleanPath        (){};

};
#endif