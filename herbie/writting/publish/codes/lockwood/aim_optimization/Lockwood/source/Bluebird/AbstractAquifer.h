//AbstractAquifer.h

#ifndef ABSTRACTAQUIFER_H
#define ABSTRACTAQUIFER_H

#include "BluebirdLibrary.h"

/****************************************************
 *  Class CAquiferABC
 *  3D Aquifer Data Abstraction
 *  for minimal knowledge access from particles and transport domains
 ***************************************************/
class CAquiferABC{  
 private:/*----------------------------------------------------------*/

 public:/*-----------------------------------------------------------*/
  //Constructors:
	CAquiferABC(){}
	virtual ~CAquiferABC(){}

  //Accessor functions
	virtual double GetCond              (const pt3D &pt) const=0;
	virtual double GetPoro              (const pt3D &pt) const=0;
	virtual double GetHead              (const pt3D &pt,const double &t) const=0;                       
	virtual vector GetVelocity2D        (const pt3D &pt,const double &t) const=0;                       
	virtual vector GetVelocity3D        (const pt3D &pt,const double &t) const=0;                       
	virtual vector GetEffVelocity2D     (const pt3D &pt,const double &t, const disp_info &disp) const=0;
	virtual vector GetEffVelocity3D     (const pt3D &pt,const double &t, const disp_info &disp) const=0;
	
  virtual double GetLeakage           (const pt3D &pt,const double &t,const leak_type ltype) const=0; 
  virtual double GetBaseflow          (const pt3D &pt, const double t) const=0;
};

#endif