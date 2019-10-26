//AbstractMultiLayer.h

#ifndef MULTILAYER_H
#define MULTILAYER_H

#include "BluebirdLibrary.h"
#include "AbstractMultiLayer.h"
#include "MLAnalyticElem.h"
#include "MLFarField.h"
#include "MLPropZone.h"

//class CmlFarField; //TMP DEBUG
//class CmlPropZone; //TMP DEBUG

enum mlelementtype{REGULAR_MLELEM,GIVEN_MLELEM,HEAD_SPECIFIED_MLELEM,TRANSIENT_MLELEM};

/****************************************************
 *  Class CMultiLayer
 *  Multi-Layer Data Abstraction
 *  container for multilayer elements
 ***************************************************/
class CMultiLayer:public CMultiLayerABC{  
 private:/*----------------------------------------------------------*/
 private:/*----------------------------------------------------------*/
	//Static Variables: All for iterative solve routine----------------

	static bool				solved;              //true if iteration complete   
	static int				iter;                //current iteration

	static int				MinMLayerIterations;  //minimum local solve iterations
	static int				MaxMLayerIterations;  //maximum local solve iterations
  static double			MLayerTolerance;      //maximum acceptable local tolerance

	static int				WriteInterval;       //solution written every this many iterations //should be elsewhere

  static cmplex			black_hole;          //location where log terms go to zero

	//Member Variables------------------------------------------------ 

	int               nLayers;

  CmlAnalyticElem **pElemArray;          //Array of pointers to all Analytic Elements in layer
  int               nElems;              //total number of elements in Layer
	CmlAnalyticElem **pHSElems;            //Array of pointers to Head-Specified Analytic Elements
	int               nHSElems;            //number of head-specified elements
	CmlAnalyticElem **pGivenElems;         //Array of pointers to Given Elements
	int               nGivenElems;         //Number of given Elements
	CmlFarField      *pFarField;           //Pointer to far field element

  double           *conductivity;        //background conductivity
  double           *thickness;           //background aquifer thickness
  double           *porosity;            //background porosity
	double           *aqtardporo;          //background aquitard porosity
	double           *aqtardcond;          //background aquitard K
	double           *aqtardthick;         //background aquitard thickness
	double           *conductance;         //background aquitard conductance (K*t)

	double            top_elev;            //elevation of multilayer top
	double            base_elev;           //elevation of multilayer base 
	double            total_transmissivity;//total transmissivity of system

	int								nCondZones;					 //number of conductivity inhomogeneities
	int								nThickZones;				 //number of thickness inhomogeneities
	int								nPoroZones;          //number of porosity inhomogeneities
	int								nAqCondZones;				 //number of aquitard conductivity inhomogeneities
	int								nAqThickZones;			 //number of aquitard thickness inhomogeneities
	int								nAqPoroZones;        //number of aquitard porosity inhomogeneities
	CmlPropZone			**pCondZoneArray;      //Array of pointers to conductivity inhomogeneity zones
	CmlPropZone			**pThickZoneArray;     //Array of pointers to conductivity inhomogeneity zones
	CmlPropZone     **pPoroZoneArray;      //Array of pointers to porosity inhomogeneity zones
	CmlPropZone			**pAqCondZoneArray;    //Array of pointers to conductivity inhomogeneity zones
	CmlPropZone			**pAqThickZoneArray;   //Array of pointers to conductivity inhomogeneity zones
	CmlPropZone     **pAqPoroZoneArray;    //Array of pointers to porosity inhomogeneity zones

  CAquicludeABC		 *pLayerAbove;         //pointer to aquiclude above	
  CAquicludeABC		 *pLayerBeneath;       //pointer to aquiclude below 

	double						DeltaPot;						 //Range of comprehensive potential values in layer domain (Potmax-Potmin)
  window						Extents;						 //Extents of layer domain

	int								GetLevel(const cmplex &z, const double &elev, bool &inaquitard) const;

	void							CalculateDeltaPotential(const double t);
  
  double						ConvertToVertVelocity (const double &z,      const double &K,
																					 const double &T,      const double &B, 
																					 const double &n,      const double &head,
																					 const double &pot,    const cmplex &W,
																					 const double &leaktop,const double &leakbot) const;


 public:/*-----------------------------------------------------------*/
  //Constructors:
	CMultiLayer();
	~CMultiLayer();

	//Static Manipulator Functions:
  static void				SetSolveData         (const int miniter, const int maxiter, const double tolerance);	 
	
	void							AddToLayer           (CmlAnalyticElem *Elemptr, mlelementtype type); 
	void							AddToLayer           (CmlFarField     *FFptr);
	void							AddToLayer           (CmlPropZone     *Propptr);

	//Manipulator Functions (visible only to aquifer)
	void              SetNumLayers         (int N);
	void							SetBackgroundValues  (double					top,
																					Ironclad1DArray T, 
																					Ironclad1DArray K, 
																					Ironclad1DArray N, 
																					Ironclad1DArray t, 
																					Ironclad1DArray k,
																					Ironclad1DArray n);

  void							SetLevelAbove        (CAquicludeABC *above);
  void							SetLevelBelow        (CAquicludeABC *below);
	void							SetBlackHole         ();       //Can be called from any layer, but only affects static data
  
	//Accessor functions (virtual->inherited from CMultiLayerABC)
	int								GetNumLayers         () const;

	double						GetBase							 (const cmplex &z) const;
  double						GetTotalThickness 	 (const cmplex &z) const;

  Ironclad1DArray		GetCondVect					 (const cmplex &z) const;
  Ironclad1DArray		GetThicknessVect  	 (const cmplex &z) const;
	Ironclad1DArray		GetPoroVect					 (const cmplex &z) const;

  Ironclad1DArray		GetAqCondVect				 (const cmplex &z) const;
  Ironclad1DArray		GetAqThicknessVect   (const cmplex &z) const;
	Ironclad1DArray		GetAqPoroVect        (const cmplex &z) const;

	Ironclad1DArray		GetBottomVect				 (const cmplex &z) const;

	Ironclad1DArray		GetConductVect       (const cmplex &z) const;
	Ironclad1DArray		GetNormTransmissVect (const cmplex &z) const;
  Ironclad1DArray		GetSystemEigenvector (const cmplex &z) const;

	anisotropy			  GetAnisotropy        (const cmplex &z, const int lev) const;
	

	//Special Accessors (only really useful to elements)
	double						GetDeltaPotential  	 ()      const;														
  cmplex						GetBlackHole      	 ()      const;
	int								GetCurrentIter       ()      const;

  //Manipulator functions (extents are the only thing an element may change)
  void							UpdateExtents        (const cmplex z);												
  void							UpdateExtents        (const cmplex z, const double r);

  
	int      					GetNumElems  				 ()      const; 
	CmlAnalyticElem	 *GetElem  						 (int i) const;
	window						GetExtents					 ()      const;
	double						GetDeltaPot  				 ()      const;

	double            GetHead            (const pt3D &pt, const double &t) const;

  //member functions
  Ironclad1DArray		GetHeadVect        (const cmplex &z, const double &t) const;
	Ironclad1DArray   GetPotentialVect   (const cmplex &z, const double &t) const;
  Ironclad1DArray_z GetQxQyVect        (const cmplex &z, const double &t) const;
	Ironclad1DArray_z GetGxVect          (const cmplex &z, const double &t) const;
  Ironclad1DArray   GetCurlVect        (const cmplex &z, const double &t) const;
	Ironclad1DArray   GetLeakageVect     (const cmplex &z, const double &t,
		                                            const leak_type ltype) const;
	double            GetCompPotential   (const cmplex &z, const double &t) const;

	double						GetNetDischarge    (const int lev, const double &t) const;
  double						GetNetDischarge    (const double &t)								const;

	//I can stub out these functions and later add transport capabilities
  cmplex						GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const int lev,const double &t) const;
	double						GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, 
																					const cmplex &z3, const int lev, const double &t,
																					const leak_type ltype) const;
	void							GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, 
																					const cmplex &z3, const int lev, const double &t, 
																					double &inQ, double &outQ) const;
  void							GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const int lev, 
																					const double &t, double &Q1, 
																					double &Q2) const;
	/*cmplex						GetEffVelocity2D     (const cmplex &z, const int lev, const double &t, 
																					const disp_info &disp) const;                  
	vector						GetEffVelocity3D     (const pt3D &pt,  const double &t,   
																					const disp_info &disp) const;*/
  
	//must be able to calculate these for particle tracking
  cmplex						GetVelocity2D        (const cmplex &z, const int lev, const double &t) const;
	vector						GetVelocity3D        (const pt3D &pt,  const double &t) const;

//	double          GetBaseflow          (const cmplex &z, const double t) const;

	void            IterativeSolve       (double &maxchange, double &maxobj, const double &t, ofstream &PROGRESS);
 // void            ExplicitSolve        (ofstream &PROGRESS, const double &t); 
//	void            SolveFluxMatrix      (double &maxchange, double &maxobj, const double &t);

  void            WriteItself          (ofstream &SOL, const double &t) const;
  void            WriteOutput          (const double &t) const;

};
#endif