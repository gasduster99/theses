#ifndef ANALYTICELEM_H
#define ANALYTICELEM_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"

//information packet sent for explicit matrix construction
/*---------------------------------------------------------------------*/
struct MatrixInfo{
	int    nctrl;                           //number of control points
	cmplex zctrl  [MAXCONTROLPTS];          //locations of control points
	double elemrhs[MAXCONTROLPTS];          //element rhs at control points
	double unit   [MAXDOF][MAXCONTROLPTS];  //unit influence array at control points
	double phiCoeff;                        //
	double QxCoeff;
	double QyCoeff;
};

/***********************************************************************
  Class CAnalyticElem
  Analytic Element Data Abstraction
------------------------------------------------------------------------
Generic parent class for all analytic elements. Each element must:
	1) Be able to solve for itself
	2) Be associated with potential function, W, Gx and discharge function
***********************************************************************/
class CAnalyticElem { 
 protected:/*----------------------------------------------------------*/

	//Static Variables---------------------------------------------------
  static int             TotalElems;				///holds the number of elements created
	static CAnalyticElem **pAllElements;			///an array of pointers to every element created
	static int             DefaultPrecision;	///the default precision for all elements
	static cmplex          zBranchcutLocus;   ///all branch cuts go through this point when "BranchcutLocus" is on

	static bool            Cosmetic;					///true if results should be cosmetic (i.e. nice branch cuts, etc)
	static bool            BranchcutLocus;    ///true if branchcuts forced to go through zBranchcutLocus (for quick flux estimates)

	//Member Variables---------------------------------------------------
  int             elemID;										///identification num of element
  char           *name;											///name of element

  int             nSubElems;								///number of elements in container
  CAnalyticElem **pSubElemArray;						///pointer to an array of pointers to Analytic Elements 
																						///(the subelements of an aggregate) 

  const CSingleLayerABC *pLayer;						///Pointer to Layer container (to get W and Omega)

  COwnerABC             *pBlock;						///Pointer to superblock container (to update laurent series)
	int										 myBlockID;					///Element's Superblock ID number 

	bool									 given;							///true if behavior is given (does not need to be solved for)

 public:/*-------------------------------------------------------------*/

  //Static functions----------------------------------------------------
  static int              GetTotalElements     ();
	static void							DestroyAllElements   ();
  static void             SetDefaultPrecision  (const int prec);
	static void             SetCosmetic          (const bool on);
	static void             SetBranchCutLocus    (const bool on, const cmplex &z);

	//Member functions----------------------------------------------------

  //Constructors: 
  CAnalyticElem();
  CAnalyticElem(char									 *Name, 
								const CSingleLayerABC  *pLay);
  CAnalyticElem(			CAnalyticElem   **p, 
		            const CSingleLayerABC  *p2,
		            int											s);
  virtual ~CAnalyticElem();
  
  //Accessor functions
  int                     GetID                () const;        
	char                   *GetName              () const;
  CAnalyticElem          *GetAllSubElems       () const;
  int                     GetNumSubElems       () const;
	bool										IsGiven              () const;
	virtual bool            HasFlux              () const;

  //Assignment functions
  void                    AddToContainer       (CAnalyticElem *Elemptr); 
  virtual void            SetBlockOwner        (COwnerABC *BlockPtr, int seg, int IDinBlock);
 
  //Virtual member functions- Modeling Functionality
  virtual cmplex          GetDischargePotential(const cmplex &z, const double &t) const;
  virtual cmplex          GetW                 (const cmplex &z, const double &t) const;
	virtual cmplex          GetGx                (const cmplex &z, const double &t) const;
	virtual double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;  
	virtual double          GetCurl              (const cmplex &z, const double &t) const;  
  virtual double          GetNetDischarge      (const double &t) const;
  virtual cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const double &t) const; 
  virtual double          GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, const leak_type ltype) const; 
	virtual void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const; 
  virtual void            GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const;
  virtual double          GetCumulativeFlux    (const cmplex &z, const double &t) const;
	virtual cmplex          GetSingularStrength  (const cmplex &z, const double &t) const;

  virtual void            SolveItself          (double &change, double &objective,const double t);
	virtual void            WriteItself          (ofstream &SOL, const double &t) const=0;
	virtual bool            ReadItself           (ifstream &SOL);
  virtual void            WriteOutput          (const double &t) const;
  virtual void            UpdateBlock          (const double &t) const;
	virtual double          GetMaxError          (const double &t) const;

	//Virtual member functions- Explicit solvers
	virtual int             GetDegreesOfFreedom    () const;
	virtual void            GetMatrixBuildInfo     (MatrixInfo &info);
	virtual void            GetUnitInfluences      (const int n,const cmplex *pts,const int NumCtrl,double *uPhi, cmplex *uQ, const double t);
	virtual void            SetCoeff               (double *coeff);
	
	virtual void            GetFluxMatrixBuildInfo (MatrixInfo &info);
	virtual void            GetFluxUnitInfluences  (const cmplex *pts,const int NumCtrl,double *uPhi, const double t);
	virtual void            SetFluxCoeff           (double coeff);

  //virtual member functions- Geometric Functionality
  virtual cmplex          Centroid             () const;                   
  virtual bool            IsInside             (const cmplex &z) const;   
  virtual bool            IsInSquare           (const cmplex &zc,const double w) const;
  virtual bool            IsInCircle           (const cmplex &zc,const double r) const;
	virtual bool            PartInCircle         (const cmplex &zc,const double r) const;
  virtual bool            SharesNode           (const cmplex &zn) const; 
	virtual void            WriteGeometry        (ofstream &BASEMAP) const;

  //Cast Cheating Functions (necessary evils)
  //oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  //for Layer/inhomogeneity-
  virtual bool            CAnalyticElem::IsString()        const {return false;}
  
  //for superblocks/strings (so pointer to CAnalyticElem can access string segment functions)
  virtual cmplex          CAnalyticElem::SegCentroid        (const int seg)   const                  {ExitGracefully("CAnalyticElem::GetNumSegs",VIRTUAL_ERROR); return 0;}
	virtual int             CAnalyticElem::GetNumSegs         (             )   const                  {ExitGracefully("CAnalyticElem::GetNumSegs",VIRTUAL_ERROR); return 0;}
  virtual bool            CAnalyticElem::SegIsInSquare      (const cmplex &zc,const double w,
		                                                         const int seg)   const                  {ExitGracefully("CAnalyticElem::SegIsInSquare",VIRTUAL_ERROR); return false;}
  virtual bool            CAnalyticElem::SegIsInCircle      (const cmplex &zc,const double r,
																											       const int seg)   const									 {ExitGracefully("CAnalyticElem::SegIsInCircle",VIRTUAL_ERROR); return false;}
  virtual bool            CAnalyticElem::SegPartInCir       (const cmplex &zc,const double r,
																											       const int seg)   const									 {ExitGracefully("CAnalyticElem::SegPartInCir",VIRTUAL_ERROR);  return false;}
	virtual cmplex          CAnalyticElem::GetSegmentPotential(const int i, const cmplex &z,
																													   const double &t) const                  {ExitGracefully("CAnalyticElem::GetSegmentPotential",VIRTUAL_ERROR); return 0.0;}
  virtual cmplex          CAnalyticElem::GetSegmentW		    (const int i, const cmplex &z,
																												     const double &t) const								   {ExitGracefully("CAnalyticElem::GetSegmentW",VIRTUAL_ERROR); return 0.0;}
  virtual cmplex          CAnalyticElem::GetSegmentGx		    (const int i, const cmplex &z,
																												     const double &t) const								   {ExitGracefully("CAnalyticElem::GetSegmentGx",VIRTUAL_ERROR); return 0.0;}
	virtual double          CAnalyticElem::GetSegmentLeakage  (const int i, const cmplex &z,
																												     const double &t,const leak_type ltype) const{
																																																			ExitGracefully("CAnalyticElem::GetSegmentLeakage",VIRTUAL_ERROR); return 0.0;}
	virtual double          CAnalyticElem::GetSegmentDischarge(const int i,const double &t) const      {ExitGracefully("CAnalyticElem::GetSegmentDischarge",VIRTUAL_ERROR); return 0.0;}																
  virtual bool            CAnalyticElem::IsDisabled         (const int i)        const               {ExitGracefully("CAnalyticElem::IsDisabled",VIRTUAL_ERROR);return false;}
	
	//Virtual member functions- explicit solver
	virtual int             CAnalyticElem::GetSegDegreesOfFreedom(const int i) const                   {ExitGracefully("CAnalyticElem::GetSegDegreesOfFreedom",VIRTUAL_ERROR); return 0;}		
	virtual void            CAnalyticElem::GetSegMatrixBuildInfo (const int i, MatrixInfo &info)       {ExitGracefully("CAnalyticElem::GetSegMatrixBuildInfo",VIRTUAL_ERROR);}
	virtual void            CAnalyticElem::GetSegUnitInfluences  (const int i, const int n,
																																const cmplex *pts,const int NumCtrl,
																																double *uPhi, cmplex *uQ, 
																																const double t)											 {ExitGracefully("CAnalyticElem::GetSegUnitInfluences",VIRTUAL_ERROR);}
	virtual void            CAnalyticElem::SetSegCoeff           (const int i, double *coeff)          {ExitGracefully("CAnalyticElem::SetSegCoeff",VIRTUAL_ERROR);}
	virtual void            CAnalyticElem::GetSegFluxMatrixBuildInfo(MatrixInfo &info)                 {ExitGracefully("CAnalyticElem::GetSegFluxMatrixBuildInfo",VIRTUAL_ERROR);}
	virtual void            CAnalyticElem::GetSegFluxUnitInfluences (const cmplex *pts,
		                                                               const int NumCtrl,
																																	 double *uPhi, const double t)     {ExitGracefully("CAnalyticElem::GetSegUnitInfluences",VIRTUAL_ERROR);}
	virtual void            CAnalyticElem::SetSegFluxCoeff       (const int i, double coeff)           {ExitGracefully("CAnalyticElem::SetSegFluxCoeff",VIRTUAL_ERROR);}
	//oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
};
#endif
