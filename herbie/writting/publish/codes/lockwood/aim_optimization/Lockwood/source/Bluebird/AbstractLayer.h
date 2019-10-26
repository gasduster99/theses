//AbstractLayer.h

#ifndef ABSTRACTLAYER_H
#define ABSTRACTLAYER_H

#include "BluebirdLibrary.h"
/****************************************************
 *  Class CLayerSetABC
 *  Layer Set Data Abstraction: can be single layer (Girinskii Potential) or multilayer (i.e., Bakker & Strack)
 *  for minimal knowledge access from transport schemes or particles
 ***************************************************/
class CLayerSetABC{  
 private:/*----------------------------------------------------------*/

 public:/*-----------------------------------------------------------*/
  //Constructors:
	CLayerSetABC(){}
	virtual ~CLayerSetABC(){if (globaldebug){cout <<" DESTROYING ABSTRACT LAYER SET "<<endl;}}

  //Accessor functions
  virtual double					GetCond							 (const pt3D   &pt) const=0;
	virtual double          GetPoro							 (const cmplex &z ) const=0;

  virtual double					GetBottom						 (const cmplex &z ) const=0;
  virtual double					GetThickness 				 (const cmplex &z ) const=0;
	virtual anisotropy      GetAnisotropy        (const pt3D   &pt) const=0;

  //Manipulator functions (horizontal extents are the only thing an element may change)
  virtual void            UpdateExtents        (const cmplex z)=0;												
  virtual void            UpdateExtents        (const cmplex z, const double r)=0;
  
  //member functions
  virtual double          GetHead              (const pt3D   &pt, const double &t) const=0;
  //virtual double          GetSaturatedThickness(const cmplex &z, const double &t) const=0;
	//virtual void            GetHeadAndPotential  (const cmplex &z, cmplex &omega,
//		                                            double &head,    const double &t) const=0;
	//virtual cmplex          GetDischargePotential(const cmplex &z, const double &t) const=0;
  //virtual cmplex          GetW                 (const cmplex &z, const double &t) const=0;
	//virtual cmplex          GetGx                (const cmplex &z, const double &t) const=0;
  //virtual double          GetCurl              (const cmplex &z, const double &t) const=0;
	virtual double          GetLeakage           (const cmplex &z, const double &t,
		                                            const leak_type ltype) const=0;
	virtual double          GetNetDischarge      (const double &t)                  const=0;
  virtual cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2,const double &t) const=0;
	virtual double          GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, 
		                                            const cmplex &z3, const double &t,
																								const leak_type ltype) const=0;
	virtual void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, 
		                                            const cmplex &z3, const double &t, 
																								double &inQ, double &outQ) const=0;
  //virtual void            GetFluxDistribution  (const cmplex &z1, const cmplex &z2, 
	//	                                            const double &t, double &Q1, 
	//																							double &Q2) const=0;
  virtual vector          GetVelocity          (const pt3D &pt,  const double &t) const=0;                
	virtual vector          GetEffVelocity       (const pt3D &pt,  const double &t,const disp_info &disp) const=0;

	virtual double          GetBaseflow          (const cmplex &z, const double t) const=0;

};
/****************************************************
 *  Class CSingleLayerABC
 *  Single Layer Data Abstraction
 *  for minimal knowledge access from elements or particles
 ***************************************************/
class CSingleLayerABC{  
 private:/*----------------------------------------------------------*/

 public:/*-----------------------------------------------------------*/
  //Constructors:
	CSingleLayerABC(){}
	virtual ~CSingleLayerABC(){if (globaldebug){cout <<" DESTROYING ABSTRACT SINGLELAYER "<<endl;}}

  //Accessor functions
  virtual double					GetCond							 (const cmplex &z) const=0;
  virtual double					GetBase							 (const cmplex &z) const=0;
  virtual double					GetThick 						 (const cmplex &z) const=0;
	virtual double          GetPoro							 (const cmplex &z) const=0;
	virtual anisotropy      GetAnisotropy        (const cmplex &z) const=0;
	virtual double          GetDeltaPot  				 ()      const=0;														
  virtual cmplex          GetBlackHole      	 ()      const=0;
	virtual double          GetSeaLevel          ()      const=0;
	virtual double          GetSaltwaterSG       ()      const=0;
	virtual int             GetCurrentIter       ()      const=0;

  //Manipulator functions (extents are the only thing an element may change)
  virtual void            UpdateExtents        (const cmplex z)=0;												
  virtual void            UpdateExtents        (const cmplex z, const double r)=0;
  
  //member functions
  virtual double          GetHead              (const cmplex &z, const double &t) const=0;
  virtual double          GetSaturatedThickness(const cmplex &z, const double &t) const=0;
	virtual void            GetHeadAndPotential  (const cmplex &z, cmplex &omega,
		                                            double &head,    const double &t) const=0;
	virtual cmplex          GetDischargePotential(const cmplex &z, const double &t) const=0;
  virtual cmplex          GetW                 (const cmplex &z, const double &t) const=0;
	virtual cmplex          GetGx                (const cmplex &z, const double &t) const=0;
  virtual double          GetCurl              (const cmplex &z, const double &t) const=0;
	virtual double          GetLeakage           (const cmplex &z, const double &t,
		                                            const leak_type ltype) const=0;
	virtual double          GetNetDischarge      (const double &t)                  const=0;
  virtual cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2,const double &t) const=0;
	virtual double          GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, 
		                                            const cmplex &z3, const double &t,
																								const leak_type ltype) const=0;
	virtual void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, 
		                                            const cmplex &z3, const double &t, 
																								double &inQ, double &outQ) const=0;
  virtual void            GetFluxDistribution  (const cmplex &z1, const cmplex &z2, 
		                                            const double &t, double &Q1, 
																								double &Q2) const=0;
	virtual cmplex          GetSingularStrength  (const cmplex &z, const double &t) const=0;
  virtual cmplex          GetVelocity2D        (const cmplex &z, const double &t) const=0;
	virtual vector          GetVelocity3D        (const pt3D &pt,  const double &t) const=0;

	virtual cmplex          GetEffVelocity2D     (const cmplex &z, const double &t, 
																								const disp_info &disp) const=0;                  
	virtual vector          GetEffVelocity3D     (const pt3D &pt,  const double &t,   
		                                            const disp_info &disp) const=0;

	virtual double          GetBaseflow          (const cmplex &z, const double t) const=0;

};

/****************************************************
 *  Class CAquicludeABC
 *  Aquiclude Data Abstraction
 *  for minimal knowledge access from elements or particles
 ***************************************************/
class CAquicludeABC{  
 private:/*----------------------------------------------------------*/

 public:/*-----------------------------------------------------------*/
  //Constructors:
	CAquicludeABC(){}
	virtual ~CAquicludeABC(){if (globaldebug){cout <<" DESTROYING ABSTRACT AQUICLUDE "<<endl;}}

  //Accessor functions
  virtual double					GetConductance			 (const cmplex &z) const=0;
  virtual double					GetElevation				 (const cmplex &z) const=0;

	virtual double          GetDeltaLeakage  		 () const=0;

	virtual CLayerSetABC   *GetLayerAbove				 () const=0;
	virtual CLayerSetABC   *GetLayerBelow				 () const=0;
  
	//Manipulator functions (extents are the only thing an element may change)
  virtual void            UpdateExtents        (const cmplex &z)=0;
  virtual void            UpdateExtents        (const cmplex &z, const double &r)=0;
 
  //member functions
	virtual cmplex          GetDischargePotential(dir direct, const cmplex &z,const double &t) const=0;
  virtual cmplex          GetW                 (dir direct, const cmplex &z,const double &t) const=0;
  virtual double          GetLeakage           (const cmplex &z,const double &t) const=0;
	virtual double          GetDesiredLeakage    (const cmplex &z,const double &t) const=0;
  virtual double          GetNetFlux           (const double &t) const=0;
 
};

/****************************************************
 *  Class COwnerABC
 *  Owner Data Abstraction
 *  for minimal knowledge access from elements 
 ***************************************************/
class COwnerABC{
 private:/*----------------------------------------------------------*/

 public:/*----------------------------------------------------------*/
	 virtual bool          IsOn                 ()	const=0;
	 virtual void          Update               (int OwnerID, int segment, double t)=0;
};
#endif