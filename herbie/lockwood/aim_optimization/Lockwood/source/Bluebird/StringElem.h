#ifndef STRINGELEM_H
#define STRINGELEM_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"



/********************************************************************
  Class CStringElem
  Analytic String Element Data Abstraction
  parents: CAnalyticElem
********************************************************************/

const double ISFAR=            1.2;     // |Z| at which far-field estimates are used

const int    MAX_FAR_FIELD=   70;       // maximum number of terms in f.f. expansion
const int    MAX_ORDER=       45;       // maximum order of chebyshev series
const int    MAXLINECONTROL= 200;       // maximum number of control points

const double CONTROLEND=       0.98;    // control point location parameter
const double FF_POWER=        -0.5;     // truncation reduction power

const int    MIN_FAR_FIELD=    3;       // minimum number of terms in f.f. expansion (not used)
const double FF_ERROR=         1e-6;    // maximum error allowed in far field truncation (not used) 
                     
class CStringElem : public CAnalyticElem {
 
 protected:/*----------------------------------------------------------*/

	//Enumerated types----------------------------------------------------------------------------------
	enum linetype       {LINESINK, DIPOLE, DOUBLET, LINEVORTEX};
  enum constrainttype {ALL_SINGULARITIES,ZERO_AT_END,UNCONSTRAINED,NET_DISCHARGE,ASSIGNED}; 

	//Static Variables----------------------------------------------------------------------------------
  double   *c[4];				        	//generic boundary condition coeff (non-const)	//TEMP DEBUG should make static:
  static double    c1,c2,c3,c4;		//generic boundary condition coefficients (const)
	static bool      ConstrntOn[3]; //constraints on 
	static double    ConstrntVal[3];//constraint values

	static double  **CD;            //unchangable coefficient arrays            
  static double  **f;
  static double  **g; 

	static double  **DD;   
	static double  **B;

  static double   *rhs;						//globally available modifyable temp arrays
  static double   *sol;						//--A,b,& sol are for Gauss Matrix solution
  static double  **A;							//--unit and rhs are for Gensolve subroutine
  static double   *b;	
  static double  **unit;
	static cmplex   *cheby;         //temporary chebyshev coefficients 

	//Member Variables-----------------------------------------------------------------------------------
  double         **JumpCoeff;     //2-D array of jump coefficients        
	double         **FarFldCoeff;   //2-D array of far field coefficients   

	double         **bn;            //correction polynomial coeff [i][n]
	double         **en;            //1st derivative of correction coeff
	double         **hn;            //2nd derivative of correction coeff
	double         **dn;            //1st derivative of jump coeff
	double         **gn;            //2nd derivative of jump coeff

  int              order;					//order of chebyshev series representation   
	int             *FForder;       //order of far field expansion    

  cmplex          *ze;						//array of polyline endpoints (global coords) 
  cmplex				 **zctrl;					//2-d array of control points (global coords) 
  double          *X;							//array of control points (local coords) 

  int              NLines;				//number of line segments in string 
  int              nlinecontrol;	//number of control pts per line 
  
	bool            *disabled;      //if line segment is disabled, turned to true 

  linetype         ltype;					//type of line element(doublet, dipole, linesink, LS/doublet, Dipole/Doublet)
  bool             closed;				//defines closed or open string 

  COwnerABC      **pBlocks;				//pointer array of ptrs to block owners of line segments
  int             *myBlockIDs;		//array of segment block IDs
  
  //Protected Member Functions -----------------------------------------------------------------------
  cmplex      GlobalToLocal(const cmplex &z, const int i) const;
	
	void       SetConstraints(const int i, const double tmp1, const double tmp2,                  
		                        const double tmp3, const double tmp4, const constrainttype CType, 
														const double val1, const double val2);

  void     SetFarFieldCoeff(const int i);                                                       

	void      BuildUnitMatrix(const int i,const linetype LType, const bool multi_val_c); 
  void             GenSolve(double *&JumpCoeff,const int i,                                     
														const linetype LType, const bool multi_val_c,const double relax, 
														double &objective, double &maxchange);
  
	double           Clenshaw(const int n,const double &x, Unchangeable1DArray   pcoeff) const;
  cmplex					cClenshaw(const int n,const cmplex &z, Unchangeable1DArray_z pcoeff) const;
  cmplex				 cClenshaw2(const int n,const cmplex &z, Unchangeable1DArray   pcoeff, const cmplex &nthterm) const; 
	
	void          Interpolate(double *Array, const int size, const double  v1, const double v2);
  void          Interpolate(cmplex *Array, const int size, const cmplex  v1, const cmplex v2);

	cmplex          GetOmJump(const int i,const double &X,  const double &t) const;                      
  cmplex           GetWJump(const int i,const double &X,  const double &t) const;                      
	cmplex      GetSegmentLog(const int i,const cmplex &z,  const double &t) const;                    
  double      GetCumExtract(const int i,const double &X,  const leftright dir, const double &t) const;
	void  GetSegCumExtraction(const int i,const double &X1, const double &X2, double &pos1, double &neg1, double &pos2, double &neg2) const;
	
	virtual double GetSegCumulativeFlux(const int i, const cmplex &z, const double &t) const;

  double   CleanStreamFunct(const cmplex &z,const double &t) const;         

  bool             IsInside(const cmplex &z) const;  //may want to bring these down to Asink level 
  double               Area(bool &clockwise) const;  //may want to bring these down to Asink level 
  double            GetArea() const;                 //may want to bring these down to Asink level 

 public:/*----------------------------------------------------------*/
  //--------------------------------------------------------------------------------------------------
  //Constructor:
  CStringElem              ();
  CStringElem              (char									*Name,           //name of element
														const CSingleLayerABC *pLay,           //layer of element
														const cmplex					*points,         //vertices of polyline/polygon (first/last of polygon counted twice)
														const int							 NumOfLines,     //number of shape line segments
														const bool						 closedstring,   //true if polygon, false if polyline
														const linetype				 LType,          //linetype (linesink/doublet/dipole)
		                        const int							 prec);          //precision of element
 ~CStringElem              (); 

  //Static Array Creation (neccesary for solution)
  static void     Prepare();
  static void     Destroy();

  void            SetPrecision         (const int Precision,int &order, double &fold);

  //Accessor Functions
  cmplex					GetZ                 (const int pointno) const;
  cmplex				 *GetZArray            () const;
  int             GetNumSegs           () const;  
	int             GetSegDegreesOfFreedom(const int i) const;
	bool            IsDisabled           (const int i) const;

	//Manipulator Functions
  void            SetBlockOwner        (COwnerABC *BlockPtr, int seg, int IDinBlock);
  void            SetAsDisabled        (int i);
	virtual void    SetSegCoeff          (const int i, double *coeff);

  //Inherited Member Functions

	bool            HasFlux              () const;

  cmplex					GetDischargePotential(const cmplex &z, const double &t) const;
  cmplex					GetW                 (const cmplex &z, const double &t) const;
	cmplex          GetGx                (const cmplex &z, const double &t) const;
  cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const double &t) const; 
  double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;
  double          GetNetDischarge      (const double &t) const;
	double          GetNetCurl           (const double &t) const;
	void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const;
  void            GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const;
  double          GetCumulativeFlux    (const cmplex &z, const double &t) const;
	cmplex					GetSingularStrength  (const cmplex &z, const double &t) const;

  cmplex					GetSegmentPotential  (const int i, const cmplex &z, const double &t) const; 
  cmplex					GetSegmentW          (const int i, const cmplex &z, const double &t) const;
	cmplex          GetSegmentGx         (const int i, const cmplex &z, const double &t) const;
  cmplex          GetSegFluxThruFace   (const int i, const cmplex &z1,const cmplex &z2, const double &t) const;
  double          GetSegmentLeakage    (const int i, const cmplex &z, const double &t,const leak_type ltype) const;
  double          GetSegmentDischarge  (const int i, const double &t) const;
	double          GetSegmentCurl       (const int i, const double &t) const;
	void            GetSegmentBudget     (const int i, const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const;
  void            GetSegDistribution   (const int i, const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const;

  void            SolveItself          (double &change, double &maxchange,const double t)=0;
	double          GetMaxError          (const double &t) const=0; 
	void            WriteOutput          (const double &t) const=0;
	
  void            WriteItself          (ofstream &SOL, const double &t) const; 
	bool            ReadItself           (ifstream &SOL);

	void            UpdateBlock          (const double &t) const;

	void            GetSegMatrixBuildInfo(const int i, MatrixInfo &info);
	void            GetSegUnitInfluences (const int i, const int n, const cmplex *pts,
		                                    const int NumCtrl, double *uPhi, cmplex *uQ, const double t);            
	//Geometric Member Functions
  bool            IsString						 () const;                                              

  cmplex					Centroid						 () const;
  bool            IsInSquare					 (const cmplex &zc,const double w) const;
  bool            IsInCircle					 (const cmplex &zc,const double r) const;
  bool            PartInCircle				 (const cmplex &zc,const double r) const;
	bool            SharesNode					 (const cmplex &zn) const; 

  cmplex					SegCentroid					 (const int seg) const;                                  
  bool            SegIsInCircle				 (const cmplex &zc,const double r,const int seg) const;
  bool            SegIsInSquare				 (const cmplex &zc,const double w,const int seg) const;  
  bool            SegPartInCir				 (const cmplex &zc,const double r,const int seg) const; 

	void            WriteGeometry        (ofstream &BASEMAP) const;
};
/********************************************************************
  Class CLineElem
  Analytic Line Element Data Abstraction
  parents: CAnalyticElem
********************************************************************/
class CLineElem : public CAnalyticElem {
 
 protected:/*----------------------------------------------------------*/

	//Enumerated types----------------------------------------------------------------------------------
	enum linetype       {LINESINK, DIPOLE, DOUBLET, LINEVORTEX};
  enum constrainttype {ALL_SINGULARITIES,ZERO_AT_END,UNCONSTRAINED,NET_DISCHARGE,ASSIGNED}; 

	//Static Variables----------------------------------------------------------------------------------
  double   *c[4];				        	//generic boundary condition coeff (non-const)	//TEMP DEBUG should make static:
  static double    c1,c2,c3,c4;		//generic boundary condition coefficients (const)
	static bool      ConstrntOn[3]; //constraints on 
	static double    ConstrntVal[3];//constraint values

	static double  **CD;            //unchangable coefficient arrays            
  static double  **f;
  static double  **g; 

  static double   *rhs;						//globally available modifyable temp arrays
  static double   *sol;						//--A,b,& sol are for Gauss Matrix solution
  static double  **A;							//--unit and rhs are for Gensolve subroutine
  static double   *b;	
  static double  **unit;
	static cmplex   *cheby;         //temporary chebyshev coefficients 

	//Member Variables-----------------------------------------------------------------------------------
	double          *JumpCoeff;     //array of jump coefficients 
  int              order;					//order of chebyshev series representation 	

//	double          *bn;            //correction polynomial coeff
//	double          *dn;            //1st derivative of jump coeff
//	double          *en;            //1st derivative of correction coeff
//	double          *gn;            //2nd derivative of jump coeff
//	double          *hn;            //2nd derivative of correction coeff

	double          *FarFldCoeff;   //array of far field coefficients   
	int              FForder;       //order of far field expansion    
	
  int              nlinecontrol;	//number of control pts per line 

  cmplex           z1,z2;					//line endpoints (global coords) 
  cmplex				  *zctrl;					//array of control points (global coords) 
  double          *X;							//array of control points (local coords) 
  
	bool             disabled;      //if line segment is disabled, turned to true 

  linetype         ltype;					//type of line element(doublet, dipole, linesink, LS/doublet, Dipole/Doublet)

  
  //Protected Member Functions -----------------------------------------------------------------------
  void       SetConstraints(const double tmp1, const double tmp2,                  
		                        const double tmp3, const double tmp4, const constrainttype CType, 
														const double val1, const double val2);
  void     SetFarFieldCoeff();                                                       

	void      BuildUnitMatrix(const linetype LType, const bool multi_val_c); 
  void             GenSolve(double *&JumpCoeff,                                     
														const linetype LType, const bool multi_val_c,const double relax, 
														double &objective, double &maxchange);
  
	double           Clenshaw(const int n,const double &x, Unchangeable1DArray   pcoeff) const;
  cmplex					cClenshaw(const int n,const cmplex &z, Unchangeable1DArray_z pcoeff) const;
  cmplex				 cClenshaw2(const int n,const cmplex &z, Unchangeable1DArray   pcoeff, const cmplex &nthterm) const; 
	
	void          Interpolate(double *Array, const int size, const double  v1, const double v2);
  void          Interpolate(cmplex *Array, const int size, const cmplex  v1, const cmplex v2);

	cmplex          GetOmJump(const double X,const double t) const;                      
  cmplex           GetWJump(const double X,const double t) const;                      
	cmplex             GetLog(const cmplex &z,const double &t) const;                    
  double      GetCumExtract(const double X, const leftright dir, const double t) const;

 public:/*----------------------------------------------------------*/
  //--------------------------------------------------------------------------------------------------
  //Constructor:
  CLineElem              ();
  CLineElem              (char									*Name,
											    const CSingleLayerABC *pLay, 
													const cmplex					 end1, 
													const cmplex					 end2,
													const linetype				 LType,
													const int							 prec);
 ~CLineElem              (); 

  //Static Array Creation (neccesary for solution)
  static void     Prepare();
  static void     Destroy();

  void            SetPrecision         (const int Precision,int &order, double &fold);

  //Accessor Functions
  bool            IsDisabled           () const;
  cmplex					GetZ                 (const int pointno) const;
	bool            HasFlux              () const;

	//Manipulator Functions
  void            SetAsDisabled        ();
	virtual void    SetCoeff             (double *coeff);

  //Inherited Member Functions
  cmplex					GetDischargePotential(const cmplex &z, const double &t) const;
  cmplex					GetW                 (const cmplex &z, const double &t) const;
	cmplex          GetGx                (const cmplex &z, const double &t) const;
  cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const double &t) const; 
  double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;
  double          GetNetDischarge      (const double &t) const;
	void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const; 

  void            SolveItself          (double &change, double &maxchange,const double t)=0;
 	double          GetMaxError          (const double &t) const=0; 
	void            WriteOutput          (const double &t) const=0;

  void            WriteItself          (ofstream &SOL, const double &t) const; 
	bool            ReadItself           (ifstream &SOL);

	int             GetDegreesOfFreedom  () const;
	void            GetMatrixBuildInfo   (MatrixInfo &info);
	void            GetUnitInfluences    (const int n, const cmplex *pts,
		                                    const int NumCtrl, double *uPhi, cmplex *uQ, const double t);            

	//Geometric Member Functions
  cmplex					Centroid						 () const;
  bool            IsInSquare					 (const cmplex &zc,const double w) const;
  bool            IsInCircle					 (const cmplex &zc,const double r) const;
  bool            PartInCircle				 (const cmplex &zc,const double r) const;
	bool            SharesNode					 (const cmplex &zn) const; 

	void            WriteGeometry        (ofstream &BASEMAP) const;
};

#endif

 
	//Unused Divergence Element Functions/Data
	//cmplex         **LJumpCoeff;		//pointer to 2-d array of leakage jump coeff
	//cmplex				 **LFarFldCoeff;	//pointer to 2-d array of Leakage Far Field coeff (perhaps unnecc.)
	//cmplex				 **LTaylorCoeff;  //pointer to 2-d array of Leakage Taylor Coeff
																	//should make above real only, with 2 pairs, or an array of 2
  //linetype         Leakltype;			//type of divergence line element (doublet, dipole, linesink, LS/doublet, Dipole/Doublet)	
	//double        GetHeadJump(const int i,const double X,const double t) const;
  //double    GetHeadGradJump(const int i,const double X,const double t) const;
	//cmplex   GetSegmentLeakOm(const int i,const cmplex &z,const double &t) const;
	//cmplex    GetSegmentLeakW(const int i,const cmplex &z,const double &t) const;
	//cmplex                  h(const cmplex *pCoeff, const int n) const;
	//cmplex                 hh(const cmplex *pCoeff, const int n) const;