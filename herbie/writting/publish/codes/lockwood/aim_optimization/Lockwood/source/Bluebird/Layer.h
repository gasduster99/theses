#ifndef LAYER_H
#define LAYER_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "FarField.h"
#include "Superblock.h"
#include "PropertyZone.h"
#include "TaylorCir.h"
#include "Grid.h"
#include "AbstractMesh.h"

/****************************************************
  Class CSingleLayer
  Single Layer Element Container Data Abstraction
	Can represent confined and unconfined conditions using Girinskii Potential
  parents: CSingleLayerABC, CLayerSetABC (Abstract Base Class)
****************************************************/

enum elementtype{REGULAR_ELEM,GIVEN_ELEM,HEAD_SPECIFIED_ELEM,TRANSIENT_ELEM};

class CSingleLayer : public CSingleLayerABC{  

 private:/*----------------------------------------------------------*/
	//Static Variables: All for iterative solve routine----------------

	static bool     solved;              //true if iteration complete   
	static int      iter;                //current iteration

	static bool     directflux;          //true if head-specified elems solved directly

	static bool     GridWhileSolve;      //Grid while solve on
	static window   GWSW;								 //Grid while solve window
	static int      GWSresolution;			 //Grid while solve gridding resolution
	
	static bool     FullAnalyze;

	static int      MinLayerIterations;  //minimum local solve iterations
	static int      MaxLayerIterations;  //maximum local solve iterations
  static double   LayerTolerance;      //maximum acceptable local tolerance

	static bool     FFeverytime;				 //if farfield is solved after every element	
	static int      WriteInterval;       //solution written every this many iterations //should be elsewhere

	static double   sealevel;						 //elevation of saltwater 
	static double   saltwaterSG;         //ratio of saltwater density to freshwater density ~1.035

  static cmplex   black_hole;          //location where log terms go to zero

	//Member Variables------------------------------------------------ 

  CAnalyticElem **pElemArray;          //Array of pointers to all Analytic Elements in layer
  int             nElems;              //total number of elements in Layer
	CAnalyticElem **pHSElems;            //Array of pointers to Head-Specified Analytic Elements
	int             nHSElems;            //number of head-specified elements
	CAnalyticElem **pGivenElems;         //Array of pointers to Given Elements
	int             nGivenElems;         //Number of given Elements
	CFarField      *pFarField;           //Pointer to far field element

  CSuperblock    *pMasterBlock;        //pointer to master superblock

  double          conductivity;        //background conductivity
  double          base_elev;           //background base elevation
  double          thickness;           //background aquifer thickness
  double          porosity;            //background porosity

	int             nCondZones;					 //number of conductivity inhomogeneities
	int             nBaseZones;					 //number of base inhomogeneities
	int             nThickZones;         //number of thickness inhomogeneities
	int             nPoroZones;          //number of porosity inhomogeneities
	CPropZone			**pCondZoneArray;      //Array of pointers to conductivity inhomogeneity zones
	CPropZone			**pBaseZoneArray;      //Array of pointers to base inhomogeneity zones
	CPropZone			**pThickZoneArray;     //Array of pointers to thickness inhomogeneity zones
	CPropZone     **pPoroZoneArray;      //Array of pointers to porosity inhomogeneity zones

  int             level;               //vertical location within aquifer system
  CAquicludeABC  *pLayerAbove;         //pointer to leaky layer above	
  CAquicludeABC  *pLayerBeneath;       //pointer to leaky layer below 

	int             nTaylorCirs;				 //number of taylor series circles
	CTaylorCir		**pTaylorCirArray;		 //Array of pointers to taylor series circles

	CMesh          *pInterpGrid;         //Grid for interpolation
	bool            InterpolationOn;     //True if interpolation is on
	double         *QxInterp;            //Interpolated Qx
	double         *QyInterp;            //Interpolated Qy

	double        **E;                   //Direct flux matrix (LHS)
	double        **B;                   //Direct flux matrix (RHS) 

	double          DeltaPot;            //Range of potential values in layer domain (Potmax-Potmin)
  window          Extents;						 //Extents of layer domain

  void            StartAnalysis         (ofstream &ANALYSIS);
	void            PrintAnalysis         (ofstream &ANALYSIS,  const int iter, const int i,
																				 const double change, const double &t);

	void            IdentifyDeltaPotential(const double t);
  
	void            Grid                  (window W, int resolution, int iter, int i, double t) const;
  double          ConvertToVertVelocity (const double &z,      const double &K,
		                                     const double &T,      const double &B, 
																				 const double &n,      const double &head,
																		     const double &pot,    const cmplex &W,
																				 const double &leaktop,const double &leakbot) const;
	cmplex          ConvertToEffective    (const double &pot,    const double &head, 
													               const double &K,      const double &T,  
																				 const double &n,      const cmplex &nprime, 
													               const cmplex &Gx,     const cmplex &Gy, 
																				 const cmplex &W,      const disp_info &disp) const;
	void            BuildDirectFluxMatrix();
	void            SolveDirectFluxMatrix();

 public:/*----------------------------------------------------------*/

  //Static Variables
  static bool     fresh;               //true if first iteration

	//Static Manipulator Functions:
  static void     SetSolveData         (const int miniter, const int maxiter, const double tolerance);	 
  static void     SetFFEveryTime       ();
	static void     SetGridWhileSolve    (const double n, const double e, const double w, const double s, const int res);
	static void     SetWriteInterval     (const int interval);
	static void     SetFullAnalysis      ();
  static void			SetCoastalInfo       (double sea_elev, double brineSG);
  
	void            SetBlackHole         ();       //Can be called from any layer, but only affects static data
  cmplex          GetBlackHole      	 () const; //Can be called from any layer, but only requests static data

	//Constructors:
  CSingleLayer();
 ~CSingleLayer();

  //Accessor functions
  double					GetCond							 (const cmplex &z) const;
  double					GetBase							 (const cmplex &z) const;
  double					GetThick 						 (const cmplex &z) const;
	double          GetPoro							 (const cmplex &z) const;
	anisotropy      GetAnisotropy        (const cmplex &z) const;

	int      			  GetNumElems  				 ()      const; 
	CAnalyticElem  *GetElem  						 (int i) const;
  CSuperblock    *GetMasterBlock			 ()      const;
	window          GetExtents					 ()      const;
	double          GetDeltaPot  				 ()      const;

	double          GetSeaLevel          ()      const;
	double          GetSaltwaterSG       ()      const;
  
  //Manipulator functions
  void            SetBackgroundValues  (double B, double T, double K, double n, int L);
  void            SetMasterBlock       (CSuperblock *Master);
  void						SetLevelAbove        (CAquicludeABC *above);
  void						SetLevelBelow        (CAquicludeABC *below);
	void            SetInterpolationGrid (CMesh *pGrid);

	//Assignment Functions
	void            AddToLayer           (CAnalyticElem *Elemptr, elementtype type); 
	void            AddToLayer           (CFarField     *FFptr);
	void            AddToLayer           (CPropZone     *Propptr);
	void            AddToLayer           (CTaylorCir    *Cirptr);
  
	//Manipulator functions (inherited from ABC)
  void            UpdateExtents        (const cmplex z);
  void            UpdateExtents        (const cmplex z, const double r);

  //Member functions (inherited from abstract base class)
  double          GetHead              (const cmplex &z, const double &t) const;
  double          GetSaturatedThickness(const cmplex &z, const double &t) const;
	void            GetHeadAndPotential  (const cmplex &z,cmplex &omega,double &head,const double &t) const;

  cmplex          GetDischargePotential(const cmplex &z, const double &t) const;
  cmplex          GetW                 (const cmplex &z, const double &t) const;
	cmplex          GetGx                (const cmplex &z, const double &t) const;
	double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;
	double          GetCurl              (const cmplex &z, const double &t) const;
  double          GetNetDischarge      (const double &t) const;
  cmplex          GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const double &t) const;
  double          GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, const leak_type ltype) const; 
	void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const; 
  void            GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const;
	cmplex          GetSingularStrength  (const cmplex &z, const double &t) const;

  cmplex          GetVelocity2D        (const cmplex &z, const double &t) const;
	vector          GetVelocity3D        (const pt3D &pt,  const double &t) const;
  cmplex          GetEffVelocity2D     (const cmplex &z, const double &t, const disp_info &disp) const;	
	vector          GetEffVelocity3D     (const pt3D &pt,  const double &t, const disp_info &disp) const;

	double          GetBaseflow          (const cmplex &z, const double t) const;   		

	//Member functions (invisible to elements)
	int             GetCurrentIter       () const;

  void            IterativeSolve       (double &maxchange, double &maxobj, const double &t, ofstream &PROGRESS);
  void            ExplicitSolve        (ofstream &PROGRESS, const double &t); 
	void            SolveFluxMatrix      (double &maxchange, double &maxobj, const double &t);

	void            InterpolateQxQy      (const double &t);

  void            WriteItself          (ofstream &SOL, const double &t) const;
  void            WriteOutput          (const double &t) const;

};


#endif