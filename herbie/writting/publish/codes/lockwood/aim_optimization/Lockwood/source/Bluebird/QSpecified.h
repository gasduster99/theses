//QSpecified.h
#ifndef QSPEC_H
#define QSPEC_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "Flownode.h"

/************************************************************************
  Class CQSpec
  Analytic Discharge Specified doublet Element Data Abstraction
  parents: CAnalyticElem, CStringElem
************************************************************************/

class CQSpec : public CStringElem {
 private:/*----------------------------------------------------------*/
  double           *pQ;        //pointer to array of unit Discharge at endpoints
  double          **Qctrl;     //array of unit discharges at all control points

	CFlowNode       **flownodes;    //array of pointers to surface flow nodes
	int              *outflow_ID;   //array of surface flow IDs
	flowdir           direction;    //direction of surface flow

	bool              solved;

  double            GetSegCumulativeFlux(const int i, const cmplex &z, const double &t) const;		
		
 public:/*-----------------------------------------------------------*/
  //Constructor
  CQSpec(char									 *Name,          //name of element
				 const CSingleLayerABC *pLay,          //layer of element
				 const cmplex					 *points,        //vertices of element polyline
				 const int							NumOfLines,    //number of lines
				 const double					 *discharge);    //specified discharge at element vertices
 ~CQSpec();

	static CQSpec *Parse        (ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);
	void           SetPrecision (const int Precision,int &order, double &fold);

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &objective,const double t);
  void            WriteOutput(const double &t) const;
	double          GetMaxError(const double &t) const;


};

#endif
