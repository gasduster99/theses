
#ifndef RIVER_H
#define RIVER_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "Flownode.h"

/************************************************************************
  Class CRiver
  Analytic Head-Specified Linesink Element Data Abstraction
  parents: CAnalyticElem, CStringElem
************************************************************************/

class CRiver : public CStringElem {
 protected:/*----------------------------------------------------------*/
	static double     relax;        //relaxation coefficient
  
	double           *head;         //array of Heads at endpoints
  double          **headctrl;     //array of heads at all control points

	CFlowNode       **flownodes;    //array of pointers to surface flow nodes
	int              *outflow_ID;   //array of surface flow IDs
	flowdir           direction;    //direction of surface flow

	double GetSegCumulativeFlux(const int i, const cmplex &z, const double &t) const;

 public:/*-------------------------------------------------------------*/
  //Constructor
  CRiver(char									 *Name,          //name of element
		     const CSingleLayerABC *pLay,          //layer of element
				 const cmplex					 *points,        //vertices of river polyline
				 const int							NumOfLines,    //number of polyline segments
				 const double					 *heads,         //specified head at vertices
	       const int							prec);         //precision of element
 ~CRiver();

	static void    SetRelaxation(double relaxation);

  static CRiver *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &objective,const double t);
  void GetSegMatrixBuildInfo (const int i,MatrixInfo &info);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};


/***********************************************************************
  Class CRiverSegment
  Analytic River Linesink Element Data Abstraction
------------------------------------------------------------------------
  parents: CAnalyticElem, CLineElem

	Issues with this representation
	  -Cosmetic Branch cuts require string communication
		-Should be done through CFlownode
***********************************************************************/
class CRiverSegment: public CLineElem {
 protected:/*----------------------------------------------------------*/
	
  static double relax;

	     
	double     head1,head2;
	double    *headctrl;
	
	CFlowNode  *flownode1,*flownode2;
	int        outflowID;
	flowdir    direction;

	cmplex     zbranch;

 public:/*-----------------------------------------------------------*/
	CRiverSegment();
	CRiverSegment(char									 *Name,   //name of element
								const CSingleLayerABC  *pLay,   //layer of element
								const cmplex						end1,   //first vertex of inhom segment
								const cmplex						end2,   //second vertex of inhom segment
								const double						h1,     //head at end1
								const double						h2,     //head at end2
								const cmplex						zbranch,//location of branch cut target direction
								const int								prec);  //precision of element
	~CRiverSegment();

	static void     SetRelaxation(double relaxation);

  static CRiverSegment **Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name, int &NumSegs);


  //Member Functions (Inherited from CAnalyticElem via CLineElem)
  cmplex          GetDischargePotential(const cmplex &z, const double &t) const;

  void            SolveItself(double &change, double &objective,const double t);
	void            GetMatrixBuildInfo (MatrixInfo &info);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};

#endif
