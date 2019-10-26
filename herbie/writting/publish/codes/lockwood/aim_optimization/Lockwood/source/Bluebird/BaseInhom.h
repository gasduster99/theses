#ifndef BASEINHOM_H
#define BASEINHOM_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "PropertyZone.h"

/************************************************************************
 *  Class CBaseInhom
 *  Analytic Base and/or Thickness Inhomogeneity doublet Element Data Abstraction
 *  parents: CAnalyticElem, CStringElem
 ***********************************************************************/

class CBaseInhom: public CStringElem {
 private:/*----------------------------------------------------------*/
  double              Bin;          //base within inhomogeneiety
	double              Tin;          //thickness within inhomogeneity

	CPolyPropZone      *pBaseZone;    //pointer to base property zone
	CPolyPropZone      *pThickZone;   //pointer to thickness property zone

	static CBaseInhom **AllInhoms;    //array of pointers to all base inhoms
	static int          NumInhoms;    //number of base inhoms created

 public:/*-----------------------------------------------------------*/
  //Constructor
  CBaseInhom(char									 *Name,         //name of element
						 const CSingleLayerABC *pLay,         //layer of element
						 const cmplex					 *points,       //vertices of element polygon (first/last counted twice)
						 const int							NumOfLines,   //number of polygon sides
						 const double						base,         //base within polygon
	           const double						thick,        //thickness within polygon
						 const int							prec);        //precision of element
 ~CBaseInhom();

  static CBaseInhom* Parse (ifstream &input,int &l,CSingleLayerABC *pLay, char * Name);

  static void        SiftThroughInhoms();
	static void               Destroy   ();

  //Accessor Functions:
	double          GetBase     () const;
	double          GetThick    () const;
	CPropZone*      GetBaseZone () const;
	CPropZone*      GetThickZone() const;

  //Member Functions (Inherited from CAnalyticElem via CStringElem
  void            SolveItself(double &change, double &objective,const double t);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};

#endif