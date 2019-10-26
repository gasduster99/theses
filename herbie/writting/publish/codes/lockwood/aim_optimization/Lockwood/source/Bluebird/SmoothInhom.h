//SmoothInhom.h
#ifndef SMOOTHINHOM_H
#define SMOOTHINHOM_H

#include "AreaVortex.h"

/************************************************************************
 *  Class CSmoothInhom
 *  Analytic smooth conductivity area vortex Element Data Abstraction
 *  parents: CAnalyticElem, CStringElem, CAreaVortex 
 ***********************************************************************/

class CSmoothInhom: public CAreaVortex {
 private:/*----------------------------------------------------------*/

	CPolyPropZone     *pCondZone;         //conductivity property zone

	bool      CalcCurl(const double t);   //Calculates curl at control points

  void	SolveMQCoeff(const double t);   //solves for Multiquadric coefficients based upon curl @ ctrl pts

	double GetDesiredCurl(const cmplex &z,const double &t) const;

 public:/*-----------------------------------------------------------*/
	//Constructor		
  CSmoothInhom(char									 *Name,            //Name of smooth inhomogeneity
							 const CSingleLayerABC *pLay,						 //layer of element 
							 const cmplex					 *points,          //vertices of polygon (first/last counted twice)
               const int							NumOfLines,      //number of polygon sides
							 const cmplex					 *condcontrol,     //numctrl locations of specified conductivity 
							 const double					 *cond,            //conductivity at numctrl given locations
							 const int							numctrl,         //number of specified conductivity locations
							 const int							prec);           //precision of element

  CSmoothInhom(char									 *Name,            //name of smooth inhomogeneity
		           const CSingleLayerABC *pLay,            //layer of element
							 const cmplex					 *points,          //vertices of polygon (first/last counted twice)
							 const int							NumOfLines,      //number of polygon sides
							 const cmplex						zcond,           //reference location of specified conductivity
							 const double						cond,            //conductivity at this location
							 const cmplex						kslope,          //slope of linearly varying conductivity function
							 const int							prec);	         //precision of element
 ~CSmoothInhom();

  static CSmoothInhom *Parse(ifstream  &input,int &l, CSingleLayerABC *pLay, char *Name);

	//Accessor Functions
 	CPropZone      *GetCondZone() const;

	//Member Functions (Inherited from CAnalyticElem via StringElem via AreaSink):
	double          GetMaxError   (const double &t) const;
	void            WriteOutput   (const double &t) const;

};

#endif