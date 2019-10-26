#ifndef DRAIN_H
#define DRAIN_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"

/***********************************************************************
  Class CDrain
  Analytic Drainy Crack Dipole Element Data Abstraction
------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
***********************************************************************/

class CDrain: public CStringElem {
 private:/*----------------------------------------------------------*/
  double           Kdrain;      //conductivity of drain
	double           thickness;   //drain thickness
				
 public:/*-----------------------------------------------------------*/
  //Constructor
  CDrain(char									 *Name,          //name of element
				 const CSingleLayerABC *pLay,          //layer of element
				 const cmplex					 *points,        //vertices of polyline
				 const int							NumOfLines,    //number of line segments in polyline
				 const double						cond,          //conductivity of drain
				 const double						thick,         //thickness of drain
	       const int							prec);         //precision of element

  static CDrain *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Accessor Functions:
  double          GetDrainCond() const;

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &objective, const double t);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};


//#include "AreaSink.h"
/*
class CReservoir: public CStringElem {
 private:
  CDrain     *drain;
	CAreaSink  *asink;
				
 public:
  //Constructor
  CReservoir(char * Name, CSingleLayerABC *pLay,cmplex *points,double extract,int ord, double OS,int NumOfLines);

  static CReservoir *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Accessor Functions:
  double          GetExtraction() const;

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &objective, const double t);
  void            WriteItself(ofstream &SOL) const;
	void            WriteOutput() const;
};*/

#endif