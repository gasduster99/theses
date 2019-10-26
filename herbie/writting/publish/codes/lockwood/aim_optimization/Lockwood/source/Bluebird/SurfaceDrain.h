
#ifndef SURFACEDRAIN_H
#define SURFACEDRAIN_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "Flownode.h"

/************************************************************************
  Class CSurfaceDrain
  Analytic Head-Specified Extraction Only Linesink Element Data Abstraction
	For modeling of seepage faces, etc.
  parents: CAnalyticElem, CStringElem
************************************************************************/
class CSurfaceDrain : public CStringElem {
 protected:/*----------------------------------------------------------*/
	static double     relax;        //relaxation coefficient
  
	double           *head;         //array of Heads at endpoints
  double          **headctrl;     //array of heads at all control points

	CFlowNode       **flownodes;    //array of pointers to surface flow nodes
	int              *outflow_ID;   //array of surface flow IDs
	flowdir           direction;    //direction of surface flow

 public:/*-------------------------------------------------------------*/
  //Constructor
  CSurfaceDrain(char									*Name,          //name of element
							  const CSingleLayerABC *pLay,          //layer of element
								const cmplex					*points,        //vertices of surface drain polyline
				        const int							 NumOfLines,    //number of polyline segments
				        const double					*heads,         //specified head at vertices
	              const int							 prec);         //precision of element
 ~CSurfaceDrain();

	static void    SetRelaxation(double relaxation);

  static CSurfaceDrain *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &objective,const double t);
  void GetSegMatrixBuildInfo (const int i,MatrixInfo &info);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};


#endif
