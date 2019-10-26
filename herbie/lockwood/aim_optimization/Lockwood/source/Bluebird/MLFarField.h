#ifndef MLFARFIELD_H
#define MLFARFIELD_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "MLAnalyticElem.h"

/****************************************************
  Class CmlFarfield
  Analytic Multi-layer Farfield Representation Data Abstraction
  parents: CmlAnalyticElem
****************************************************/

class CmlFarField : public CmlAnalyticElem {
 private:/*----------------------------------------------------------*/
  cmplex					zref;              //z-coord of reference point
  double          Constant;          //potential at reference point
  cmplex					Qo;                //uniform flow rate
  double          alpha;             //flow angle 
  bool            RefPoint;          //true if reference point specified,
  //                                   false if net extraction specified.
  double          RefHead;           //head at reference point
  double          NetExt;            //specified net extraction

	double          OldC[10];
	double          OldQ[10];
  
 public:/*-----------------------------------------------------------*/
  //Constructors  
  CmlFarField(CMultiLayerABC *pLay,            //constructor for reference point type
							cmplex					Gradient, 
							cmplex					Zref, 
							double					refhead);   
 ~CmlFarField();

  //Accessors
  cmplex						GetZref              () const;  
  cmplex						GetUnifFlow          () const;                  
  bool							IsReferencePt        () const;

  //Manipulator Functions
  void							SetZref              (cmplex z);

  //Member functions (inherited from CAnalyticElem):
	double						GetCompPotential    (const cmplex &z, const double &t) const;
	Ironclad1DArray		GetPotentialVect		(const cmplex &z, const double &t) const;
	Ironclad1DArray		GetHelmPotentialVect(const cmplex &z, const double &t) const;       
  Ironclad1DArray_z	GetQxQyVect					(const cmplex &z, const double &t) const; 
	Ironclad1DArray_z	GetGxVect						(const cmplex &z, const double &t) const;
	cmplex						GetFluxThruFace			(const cmplex &z1,const cmplex &z2,const int lev,  const double &t) const;
  Ironclad1DArray		GetLeakageVect			(const cmplex &z, const double &t,const leak_type ltype) const;
  double            GetNetDischarge     (const int    lev,const double &t) const;
	double						GetNetDischarge			(const double &t) const; 

  void							SolveItself          (double &change,double &objective,const double t);                                                       
  void							WriteItself          (ofstream &SOL, const double &t) const;
  bool							ReadItself           (ifstream &SOL);
	void							WriteOutput          (const double &t) const;
	void							UpdateBlock          (const double &t) const{}
	double						GetMaxError          (const double &t) const;

	//virtual member functions- explicit solver
	/*int								GetDegreesOfFreedom  () const;
	void							GetMatrixBuildInfo   (MatrixInfo &info);
	void							GetUnitInfluences    (const int n,const cmplex *pts,const int NumCtrl,double *uPhi, cmplex *uQ, const double t);
	void							SetCoeff             (double *coeff);*/

	//Geometric functions
  cmplex						Centroid             () const;
	bool							IsInSquare           (const cmplex &zc,const double w) const;
  bool							IsInCircle           (const cmplex &zc,const double r) const;
  bool							PartInCircle         (const cmplex &zc,const double r) const;
	double						GetArea              () const;  
	bool							SharesNode           (const cmplex &zn) const; 
};

#endif