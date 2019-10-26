#ifndef CIRCULAR_H
#define CIRCULAR_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "PropertyZone.h"
#include "Flownode.h"


const int    MAXCIRCONTROL=   500;      // maximum number of control points

/*********************************************************************
  Class CCircularElem
  Circular Analytic Element Data Abstraction
    parents: CAnalyticElem
*********************************************************************/

class CCircularElem : public CAnalyticElem {
 protected:/*----------------------------------------------------------*/
	//Member Variables
  cmplex					*OutCoeff;		//pointer to array of Laurent Series  coeff. 
	cmplex          *InCoeff;			//pointer to array of Taylor Series Coeff
	double           Q;           //net discharge from circle
  int              order;       //order of Laurent representation

  cmplex					 zcen;        //center of circle
  double           R;           //radius of circle

  cmplex					*zctrl;       //pointer to dynamic array of control pts
  //                            //in standard domain coords[lines][m]
  int              ncircontrol; //number of control pts per circle

 public:/*-------------------------------------------------------------*/
  //Constructor:
  CCircularElem(); 
	CCircularElem(char									*Name,   //name of element
								const CSingleLayerABC *pLay,   //layer of element
								const cmplex					 zc,     //center of circle
								const double					 rad,    //radius of circle
								const int							 prec);  //precision of element
 ~CCircularElem();

 	void            SetPrecision         (const int Precision,int &order, double &fold);


  //Accessor Functions
  cmplex					GetZcen              () const;
  double          GetRadius            () const;

  //Member Functions (inherited from AnalyticElem):	
	bool            HasFlux              () const;

  cmplex					GetDischargePotential(const cmplex &z, const double &t) const;
  cmplex					GetW                 (const cmplex &z, const double &t) const;
	cmplex          GetGx                (const cmplex &z, const double &t) const;
  cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const double &t) const; 
	double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;
  double          GetNetDischarge      (const double &t) const;

  void            SolveItself	         (double &change, double &maxchange,const double t)=0;
	double          GetMaxError          (const double &t) const=0;
	void            WriteOutput          (const double &t) const=0;

	bool            ReadItself           (ifstream &SOL);
  void            WriteItself          (ofstream &SOL, const double &t) const;

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
 *  Class CCirInhom
 *  Analytic Circular Inhomogeneity Element Data Abstraction
 *  parents: CAnalyticElem, CCircularElem
 ********************************************************************/
class CCirInhom : public CCircularElem {

 private:/*----------------------------------------------------------*/
  double        cond;
  CCirPropZone *CondZone;

 public:/*-----------------------------------------------------------*/
  //Constructor
  CCirInhom(char									*Name,    //name of element
		        const CSingleLayerABC *pLay,    //layer of element
						const cmplex					 zc,      //center of circle
						const double					 rad,     //radius of circle
						const double					 k,       //conductivity of circular zone
						const int							 prec);   //precision of element
  ~CCirInhom();

	//static functions
  static CCirInhom *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Accessor Functions:
  double          GetKin() const;
	CPropZone      *GetCondZone() const;

  //Member Functions (Inherited from CAnalyticElem):
  void            SolveItself(double &change, double &objective,const double t);
  void            WriteOutput(const double &t) const;
	double          GetMaxError(const double &t) const;
};

/*********************************************************************
 *  Class CCirLake
 *  Analytic Circular Head-Specified Element Data Abstraction
 *  parents: CAnalyticElem, CCircularElem
 ********************************************************************/
class CCirLake : public CCircularElem {

private:/*----------------------------------------------------------*/
  double     elev;
	
	CFlowNode *flownode;
	int        outflow_ID;

public:/*-----------------------------------------------------------*/
  //Constructor
  CCirLake(char									 *Name,       //name of element
		       const CSingleLayerABC *pLay,       //layer of element
					 const cmplex						zc,         //center of circle
					 const double						rad,        //radius of circle
					 const double						head,       //specified head on circular "lake" boundary
					 const int							prec);      //precision of element

	//static functions
  static CCirLake *Parse(ifstream &input, int &l, CSingleLayerABC *pLay, char * Name);

  //Accessor Functions:
  double          GetHead() const;
	bool            HasFlux() const;

  //Member Functions (Inherited from CAnalyticElem):
  void            SolveItself(double &change, double &objective,const double t);
  void            WriteOutput(const double &t) const;
	double          GetMaxError(const double &t) const;

};
#endif

