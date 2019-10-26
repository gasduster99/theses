#ifndef MLPROPZONE_H
#define MLPROPZONE_H

#include "MasterInclude.h"
#include "PropertyZone.h"


/****************************************************
 *  Class CmlPropZone
 *  Polygonal & Circular Multilayer Property zone Data Abstraction
 *  used as area of constant (double) property value
 ***************************************************/

class CmlPropZone {
 protected:/*----------------------------------------------------------*/

	static CmlPropZone  *pAllPropZones[MAX_PROPZONES];
	static int	         TotalPropZones;

	//Member Variables
	double           values[MAX_MLAYERS];//parameter values in zone	
	proptype         ptype;              //property type
	double           area;               //area of zone

 public:/*-------------------------------------------------------------*/
	 //Constructors
	 CmlPropZone(const proptype  typ,
		           const double   *val);
   virtual ~CmlPropZone();

   static void							DestroyAllMLPropZones();
	 static Ironclad1DArray		NestedSift      (const cmplex			 &z, 
																						 CmlPropZone			**pZones, 
																						 const int					nzones, 
																						 Ironclad1DArray		back_val);

	 //accessor functions
	 virtual Ironclad1DArray   GetValues       (const cmplex &z) const;
	 Ironclad1DArray					 GetValues       () const;
	 proptype									 GetType         () const;
   double										 GetArea         () const;

	 //Virtual Member Functions
	 virtual bool   IsInside  (const cmplex &z) const;
   virtual bool   IsInSquare(const cmplex &zc,const double w) const;
	
};

#endif