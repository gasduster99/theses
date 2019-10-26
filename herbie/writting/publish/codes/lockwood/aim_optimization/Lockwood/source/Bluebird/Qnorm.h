#ifndef QNORM_H
#define QNORM_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"

/************************************************************************
 *  Class CQnorm
 *  Analytic Normal Discharge Specified doublet Element Data Abstraction
 *  parents: CAnalyticElem, CStringElem
 ***********************************************************************/

class CQnorm : public CStringElem {
 private:/*----------------------------------------------------------*/
  double           *pQ;        //pointer to array of unit Discharge at endpoints
  double          **Qctrl;     //array of unit discharges at all control points

 public:/*-----------------------------------------------------------*/
  //Constructor
  CQnorm(char									 *Name,            //name of element
		     const CSingleLayerABC *pLay,            //layer of element
				 const cmplex					 *points,          //vertices of polyline
				 const int							NumOfLines,      //number of polyline segments
				 const double					 *discharge,       //specified discharge at element vertices
	       const int							prec);           //precision of element
	~CQnorm();

	static CQnorm *Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &objective,const double t);
	double          GetMaxError(const double &t) const;
  void            WriteOutput(const double &t) const;
};

#endif