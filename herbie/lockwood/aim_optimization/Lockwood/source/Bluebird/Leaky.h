#ifndef LEAKYWALL_H
#define LEAKYWALL_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"

/************************************************************************
  Class CLeakyWall
  Analytic Leaky Wall doublet Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
************************************************************************/

class CLeakyWall: public CStringElem {
 private:/*----------------------------------------------------------*/
  double           Kwall;       //conductivity of wall
	double           thickness;   //wall thickness

 public:/*-----------------------------------------------------------*/
  //Constructor
  CLeakyWall(char									 *Name,          //name of element
		         const CSingleLayerABC *pLay,          //layer of element
				     const cmplex					 *points,        //vertices of polyline
				     const int							NumOfLines,    //number of polyline segments 
				     const double						cond,          //conductivity of leaky wall
				     const double						thick,         //thickness of leaky wall
				     const int							prec);         //precision of element
	~CLeakyWall();

	static CLeakyWall *Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);

  //Accessor Functions:
  double          GetWallCond() const;

  //Member Functions (Inherited from CAnalyticElem via CStringElem):
  void            SolveItself(double &change, double &maxchange, const double t);
	double          GetMaxError(const double &t) const;
  void            WriteOutput(const double &t) const;
};

#endif