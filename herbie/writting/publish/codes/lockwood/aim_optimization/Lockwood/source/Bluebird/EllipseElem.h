//EllipseElem.h
#ifndef ELLIPTICAL_H
#define ELLIPTICAL_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "PropertyZone.h"
#include "Flownode.h"

/*********************************************************************
* Class CEllipseElem
*	Analytic Elliptical Element Data Abstraction
*	parents: CAnalyticElem
********************************************************************/

class CEllipseElem : public CAnalyticElem {
 protected:/*----------------------------------------------------------*/
	//Member Variables
  cmplex					*OutCoeff;    //pointer to array of Laurent Series  coeff. 
	cmplex          *InCoeff;     //pointer to array of Taylor Series Coeff

  cmplex					 zc;          //center of ellipse
  double           a;           //1/2 major axis
	double           b;           //1/2 minor axis
	double           c;           //Focal distance
	double           angle;       //angle of ellipse (in radians)
	double           Q;           //net discharge from ellipse
	double           eta0;        //"radius" of ellipse in local coords

  cmplex					*zctrl;       //array of control pts in global coords
  int              nellcontrol; //number of control pts on ellipse

  int              order;       //order of Laurent representation
  
  //Protected Member Functions (available to all Ellipses)
	cmplex GlobalToLocal(const cmplex &z) const;      //returns zeta = eta + i psi
  cmplex LocalToGlobal(const cmplex &zeta) const; 

 public:/*-------------------------------------------------------------*/
  //Constructor:
  CEllipseElem(); 
	CEllipseElem( char									*Name,         //name of element
								const CSingleLayerABC *pLay,         //layer of element
								const cmplex					 zcen,         //center of ellipse
								const double					 MajorAxis,    //major axis of ellipse
								const double					 MinorAxis,    //minor axis of ellipse
								const double					 angle,				 //angle (in radians)
								const int							 prec);        //precision of element
 ~CEllipseElem();

 	void            SetPrecision         (const int Precision,int &order, double &fold);

  //Accessor Functions
  cmplex  				GetCenter            () const;
  void            GetAxes              (double &Major, double &Minor) const;

  //Member Functions (inherited from AnalyticElem):

	bool            HasFlux              () const;

  cmplex					GetDischargePotential(const cmplex &z, const double &t) const;
  cmplex					GetW                 (const cmplex &z, const double &t) const;
	cmplex          GetGx                (const cmplex &z, const double &t) const;
  cmplex          GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const double &t) const; 
	double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;
  double          GetNetDischarge      (const double &t) const;

  void            SolveItself	         (double &change, double &maxchange,const double t)=0;
	double          GetMaxError          (const double &t) const=0;
	void            WriteOutput          (const double &t) const=0;

  void            WriteItself          (ofstream &SOL, const double &t) const;
	bool            ReadItself           (ifstream &SOL);

	void            UpdateBlock          (const double &t) const;


  cmplex					Centroid             () const;
  bool            IsInside						 (const cmplex &z) const; 
  bool            IsInSquare					 (const cmplex &zc,const double w) const;
  bool            IsInCircle           (const cmplex &zc,const double r) const;
  bool            PartInCircle         (const cmplex &zc,const double r) const;
  double          GetArea              () const;
  bool            SharesNode           (const cmplex &zn) const; 
	void            WriteGeometry        (ofstream &BASEMAP) const;

};
/*********************************************************************
*	Class CEllipseElem
*	Analytic Elliptical Inhomogeneity Data Abstraction
* parents: CAnalyticElem, CEllipseElem
********************************************************************/
class CEllipseInhom : public CEllipseElem {
 private:
	double kin;
  CEllPropZone *CondZone;

 public:
	CEllipseInhom();
 ~CEllipseInhom();
	CEllipseInhom(char									*Name,         //name of element
								const CSingleLayerABC *pLay,         //layer of element
							  const double					 cond,         //conductivity of elliptical zone
								const cmplex					 zcen,         //center of ellipse
								const double					 MajorAxis,    //major axis of ellipse
								const double					 MinorAxis,    //minor axis of ellipse
								const double					 angle,				 //angle (in radians)
								const int							 prec);        //precision of element

  static CEllipseInhom  *Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);
	
	//Accessor Functions
	double GetK() const;
	CPropZone  *GetCondZone() const;

	//Member Functions (Inherited from CAnalyticElem via CEllipseElem)
  void            SolveItself(double &change, double &objective,const double t);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};
/*********************************************************************
*	Class CEllLake
*	Analytic Elliptical Head-specified Lake Data Abstraction
* parents: CAnalyticElem, CEllipseElem
********************************************************************/
class CEllLake : public CEllipseElem {
 private:

	double     elev;
	CFlowNode *flownode;
	int        outflow_ID;

 public:
	CEllLake();
 ~CEllLake();
	CEllLake(char									 *Name,        //name of element
					 const CSingleLayerABC *pLay,        //layer of element
					 const double						head,        //specified head
					 const cmplex						zcen,        //center of ellipse
					 const double						MajorAxis,   //major axis of ellipse
					 const double						MinorAxis,   //minor axis of ellipse
					 const double						angle,			 //angle (in radians)
					 const int							prec);       //precision of element

  static CEllLake  *Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);
	
	//Accessor Functions
	double GetElev() const;
	bool            HasFlux() const;

	//Member Functions (Inherited from CAnalyticElem via CEllipseElem)
  void            SolveItself(double &change, double &objective,const double t);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};
/*********************************************************************
*	Class CEllReservoir
*	Analytic Elliptical Discharge-specified Lake Data Abstraction
* parents: CAnalyticElem, CEllipseElem
********************************************************************/
class CEllReservoir : public CEllipseElem {
 private:

	double     elev;
	CFlowNode *flownode;
	int        outflow_ID;

 public:
	CEllReservoir();
 ~CEllReservoir();
	CEllReservoir(char									 *Name,        //name of element
								const CSingleLayerABC  *pLay,        //layer of element
								const double						head,        //specified head
								const cmplex						zcen,        //center of ellipse
								const double						MajorAxis,   //major axis of ellipse
								const double						MinorAxis,   //minor axis of ellipse
								const double						angle,			 //angle (in radians)
								const int								prec);       //precision of element

  static CEllReservoir  *Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);
	
	//Accessor Functions
	double GetFlux() const;

	//Member Functions (Inherited from CAnalyticElem via CEllipseElem)
  void            SolveItself(double &change, double &objective,const double t);
	double          GetMaxError(const double &t) const;
	void            WriteOutput(const double &t) const;
};
#endif