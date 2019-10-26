#ifndef HORWELL_H
#define HORWELL_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"

/************************************************************************
  Class CHorizontalWell
  Analytic Horizontal Well linesink Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
************************************************************************/

class CHorizontalWell: public CStringElem {
 private:/*----------------------------------------------------------*/
  double           Q;     //pumping rate of well
	//double         k;     //conductivity around well (not used)
	//double         thick; //thickness of layer around well (not used)

 public:/*-----------------------------------------------------------*/
  //Constructor
  CHorizontalWell(char									*Name,        //name of element
								  const CSingleLayerABC *pLay,        //layer of element
									const cmplex					*points,      //vertices of element polyline
									const int							 NumOfLines,  //number of line segments in polyline
									const double					 prate,       //pumping rate of well
									const int							 prec);       //precision of element 
 ~CHorizontalWell();

  static CHorizontalWell *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Accessor Functions:
  double         GetPumpingRate() const;

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &maxchange, const double t);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};

#endif