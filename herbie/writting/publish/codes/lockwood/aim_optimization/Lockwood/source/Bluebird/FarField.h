#ifndef FARFIELD_H
#define FARFIELD_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"

/****************************************************
  Class CFarfield
  Analytic Farfield Representation Data Abstraction
  parents: CAnalyticElem
****************************************************/

class CFarField : public CAnalyticElem {
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
  CFarField();
  CFarField(CSingleLayerABC *pLay,            //constuctor for net extraction type
		        cmplex					 UnifFlow, 
						double		 			 NetExt);       
  CFarField(CSingleLayerABC *pLay,            //constructor for reference point type
		        cmplex					 UnifFlow, 
						cmplex					 Zref, 
						double					 refhead);   
 ~CFarField();

  //Accessors
  cmplex					GetZref              () const;  
  cmplex					GetUnifFlow          () const;                  
  bool            IsReferencePt        () const;

  //Manipulator Functions
  void            SetZref              (cmplex z);

  //Member functions (inherited from CAnalyticElem):
  cmplex					GetDischargePotential(const cmplex &z, const double &t) const;           
  cmplex					GetW                 (const cmplex &z, const double &t) const; 
	cmplex          GetGx                (const cmplex &z, const double &t) const;
	cmplex          GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const double &t) const;
  double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;
  double          GetNetDischarge      (const double &t) const; 

  void            SolveItself          (double &change,double &objective,const double t);                                                       
  void            WriteItself          (ofstream &SOL, const double &t) const;
  bool            ReadItself           (ifstream &SOL);
	void            WriteOutput          (const double &t) const;
	void            UpdateBlock          (const double &t) const{}
	double          GetMaxError          (const double &t) const;

	//virtual member functions- explicit solver
	int             GetDegreesOfFreedom  () const;
	void            GetMatrixBuildInfo   (MatrixInfo &info);
	void            GetUnitInfluences    (const int n,const cmplex *pts,const int NumCtrl,double *uPhi, cmplex *uQ, const double t);
	void            SetCoeff             (double *coeff);

  cmplex					Centroid             () const;
	bool            IsInSquare           (const cmplex &zc,const double w) const;
  bool            IsInCircle           (const cmplex &zc,const double r) const;
  bool            PartInCircle         (const cmplex &zc,const double r) const;
	double          GetArea              () const;  
	bool            SharesNode           (const cmplex &zn) const; 
};

#endif