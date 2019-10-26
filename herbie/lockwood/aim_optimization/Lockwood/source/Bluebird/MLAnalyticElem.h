#ifndef ML_ANALYTICELEM_H
#define ML_ANALYTICELEM_H

#include "BluebirdLibrary.h"
#include "AbstractMultiLayer.h"
#include "AbstractLayer.h"  //contains owner info

//information packet sent for explicit matrix construction
/*---------------------------------------------------------------------*/
/*struct MatrixInfo{
	int    nctrl;                           //number of control points
	cmplex zctrl  [MAXCONTROLPTS];          //locations of control points
	double elemrhs[MAXCONTROLPTS];          //element rhs at control points
	double unit   [MAXDOF][MAXCONTROLPTS];  //unit influence array at control points
	double phiCoeff;                        //
	double QxCoeff;
	double QyCoeff;
};*/

/***********************************************************************
  Class CMLAnalyticElem
  Analytic Multilayer Element Data Abstraction
------------------------------------------------------------------------
Generic parent class for all MultiLayer analytic elements. Each element must:
	1) Be able to solve for itself
	2) Be associated with potential function (comprehensive and helmholtz)
	    W, Gx and discharge function
***********************************************************************/
class CmlAnalyticElem { 
 protected:/*----------------------------------------------------------*/

	//Static Variables---------------------------------------------------
  static int								TotalElems;				///holds the number of elements created
	static CmlAnalyticElem	**pAllElements;			///an array of pointers to every element created
	static int								DefaultPrecision;	///the default precision for all elements

	//Member Variables---------------------------------------------------
  int												elemID;						///identification num of element
  char										 *name;							///name of element

  int												nSubElems;				///number of elements in container
  CmlAnalyticElem					**pSubElemArray;		///pointer to an array of pointers to Analytic Elements 
																							///(the subelements of an aggregate) 

  const CMultiLayerABC		 *pLayer;						///Pointer to Layer container (to get W and Omega)

  //COwnerABC							 *pOwner;						///Pointer to superblock container (to update laurent series)
	//int											myOwnerID;				///Element's Superblock ID number 

	bool											given;						///true if behavior is given (does not need to be solved for)

 public:/*-------------------------------------------------------------*/

  //Static functions----------------------------------------------------
  static int								GetTotalElements     ();
	static void								DestroyAllMLElements ();
  static void								SetDefaultPrecision  (const int prec);

	//Member functions----------------------------------------------------

  //Constructors: 
  CmlAnalyticElem					 (char									 *Name, 
														const CMultiLayerABC   *pLay);
  CmlAnalyticElem					 (			CmlAnalyticElem **pElems, 
														const CMultiLayerABC   *pLay,
														int											numelems);
  virtual ~CmlAnalyticElem();
  
  //Accessor functions
  int												GetID                () const;        
	char										 *GetName              () const;
  CmlAnalyticElem					 *GetAllSubElems       () const;
  int												GetNumSubElems       () const;
	bool											IsGiven              () const;
	virtual bool							HasFlux              () const;

  //Assignment functions
  void											AddToContainer       (CmlAnalyticElem *pElem); 
  //virtual void							SetOwner			       (COwnerABC *pOwner, int seg, int IDinOwner);
 
  //Virtual member functions- Modeling Functionality
	virtual double						GetCompPotential     (const cmplex &z, const double &t) const;
	virtual Ironclad1DArray		GetPotentialVect		 (const cmplex &z, const double &t) const;
	virtual Ironclad1DArray		GetHelmPotentialVect (const cmplex &z, const double &t) const;
  virtual Ironclad1DArray_z GetQxQyVect          (const cmplex &z, const double &t) const;
	virtual Ironclad1DArray_z GetGxVect            (const cmplex &z, const double &t) const;

	virtual Ironclad1DArray		GetLeakageVect       (const cmplex &z, const double &t, const leak_type ltype) const;  
	virtual Ironclad1DArray		GetCurlVect          (const cmplex &z, const double &t) const; 
	virtual double            GetNetDischarge      (const int    lev,const double &t) const;
  virtual double            GetNetDischarge      (const double &t) const;
  
	virtual cmplex            GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const int lev, const double &t) const; 
  virtual double            GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, const cmplex &z3, const int lev, const double &t, const leak_type ltype) const; 
	virtual void              GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const int lev, const double &t, double &inQ, double &outQ) const; 
  virtual void              GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const int lev, const double &t, double &Q1, double &Q2) const;
  virtual double            GetCumulativeFlux    (const cmplex &z , const double &t) const;

  virtual void							SolveItself          (double &change, double &objective,const double t);
	virtual void							WriteItself          (ofstream &SOL, const double &t) const=0;
	virtual bool							ReadItself           (ifstream &SOL);
  virtual void							WriteOutput          (const double &t) const;
//  virtual void							UpdateOwner          (const double &t) const;
	virtual double						GetMaxError          (const double &t) const;

	//Virtual member functions- Explicit solvers
/*virtual int								GetDegreesOfFreedom    () const;
	virtual void							GetMatrixBuildInfo     (MatrixInfo &info);
	virtual void							GetUnitInfluences      (const int n,const cmplex *pts,const int NumCtrl,double *uPhi, cmplex *uQ, const double t);
	virtual void							SetCoeff               (double *coeff);
	
	virtual void							GetFluxMatrixBuildInfo (MatrixInfo &info);
	virtual void							GetFluxUnitInfluences  (const cmplex *pts,const int NumCtrl,double *uPhi, const double t);
	virtual void							SetFluxCoeff           (double coeff);*/

  //virtual member functions- Geometric Functionality
  virtual cmplex          Centroid             () const;                   
  virtual bool            IsInside             (const cmplex &z) const;   
  virtual bool            IsInSquare           (const cmplex &zc,const double w) const;
  virtual bool            IsInCircle           (const cmplex &zc,const double r) const;
	virtual bool            PartInCircle         (const cmplex &zc,const double r) const;
  virtual bool            SharesNode           (const cmplex &zn) const; 
	virtual void            WriteGeometry        (ofstream &BASEMAP) const;

  virtual bool            IsDisabled					 (const int i) const {
		ExitGracefully("CmlAnalyticElem::IsDisabled",VIRTUAL_ERROR);return false;}
	
};
#endif
