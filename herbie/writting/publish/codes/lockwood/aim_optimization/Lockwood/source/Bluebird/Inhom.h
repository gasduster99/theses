#ifndef INHOM_H
#define INHOM_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "PropertyZone.h"

/***********************************************************************
  Class CInhom
  Analytic Inhomogeneity doublet Element Data Abstraction
------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
***********************************************************************/

class CInhom: public CStringElem {
 private:/*----------------------------------------------------------*/
  double           Kin;         //conductivity of inhomogeneity
  bool             adjacent;	  //true if touching another inhom.

	CPolyPropZone   *pCondZone;   //pointer to conductivty property zone

	static CInhom  **AllInhoms;   //array of pointers to all inhoms
	static int       NumInhoms;   //number of inhomogeneities created

 public:/*-----------------------------------------------------------*/
  //Constructor
	CInhom();
  CInhom(char									 *Name,          //name of element
				 const CSingleLayerABC *pLay,          //layer of element
				 const cmplex					 *points,        //vertices of element polygon (first/last counted twice)
				 const int							NumOfLines,    //number of polygon sides
				 const double						cond,          //conductivity of polygon
	       const int							prec);         //precision of element

  static CInhom  *Parse            (ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);
	static CInhom **CreateBoxOfInhoms(CSingleLayerABC *pLay, 
		                                double					 xmin,    double xmax, 
																		double					 ymin,    double ymax, 
		                                int							 nx, 
																		double					 mincond, double maxcond);
  static void     SiftThroughInhoms();
	static void               Destroy();

  //Accessor Functions:
	double          GetK() const;
	CPropZone      *GetCondZone() const;

	//Manipulator Functions
	void            SetAsAdjacent();

  //Member Functions (Inherited from CAnalyticElem via CStringElem)
  void            SolveItself(double &change, double &objective,const double t);
	void            GetSegMatrixBuildInfo (const int i,MatrixInfo &info);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;

	//unique member functions
	void            SolveVertices(double &change, double &objective,const double t);
};

/***********************************************************************
  Class CInhomSegment
  Analytic Inhomogeneity doublet Element Data Abstraction
------------------------------------------------------------------------
  parents: CAnalyticElem, CLineElem
***********************************************************************/
class CInhomSegment: public CLineElem {
 private:/*----------------------------------------------------------*/
	double Kin;

	static CInhomSegment  **AllInhomSegs;   //array of pointers to all inhoms
	static int              NumInhomSegs;   //number of inhomogeneities created


 public:/*-----------------------------------------------------------*/
	CInhomSegment();
	CInhomSegment(char									 *Name,   //name of element
								const CSingleLayerABC  *pLay,   //layer of element
								const cmplex						end1,   //first vertex of inhom segment
								const cmplex						end2,   //second vertex of inhom segment
								const double						cond,   //conductivity of element
								const int								prec);  //precision of element

	static void     SiftThroughInhomSegs();
	static void                  Destroy();

  static CInhomSegment **Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name, int &NumSegs, CPropZone *&kzone);

	static CInhomSegment **CreateBoxOfInhoms(CSingleLayerABC *pLay, 
		                                 double		 xmin,    double xmax, 
																		 double		 ymin,    double ymax, 
		                                 int				 nx, 
																		 double		 mincond, double maxcond);

	double          GetK() const;

  //Member Functions (Inherited from CAnalyticElem via CLineElem)
  void            SolveItself(double &change, double &objective,const double t);
	void            GetMatrixBuildInfo (MatrixInfo &info);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};

#endif