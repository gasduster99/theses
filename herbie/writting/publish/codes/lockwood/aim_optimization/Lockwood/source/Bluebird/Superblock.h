#ifndef SUPERBLOCK_H
#define SUPERBLOCK_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "PropertyZone.h"

const double SB_FF_POWER=      -0.5;     // truncation reduction power
const int    MAXBLOCKCONTROL= 500;      // maximum # of block control pts

/************************************************************************
  Class CSuperblock
  Superblock Element Container Data Abstraction
************************************************************************/

class CSuperblock : public COwnerABC{

 private: /*----------------------------------------------------------*/
  static cmplex   **f;
  static cmplex  ***G;
  static int        nblockcontrol;	 //number of control points per block
  static int        order;					 //number of laurent coeff per block
	static bool       Explicit;				 //method for solving parent coeff (false means Randall's trick)
	static bool       Nested;			     //true if using nesting (default, false for "classical" SBs)
  static bool       On;							 //true if superblocks are enabled
	static double     BlockRadius;     //radius of block (as compared to blockdiagonal/2)
  static double     BlockBuffer;     //relative radius of elements within expansion

  CSuperblock      *pParent;				 //pointer to parent block
  CSuperblock      *pChild[4];			 //pointers to children 

  cmplex           *LaurentCoeff;		 //array of laurent coeff for outer expansion 
  double            Q;							 //net discharge from block

  cmplex            zcen;						 //center of superblock
  double            width;					 //width of block square
  double            radius;					 //radius of expansion circle  

  cmplex           *zctrl;					 //array of control pts (global coords)  

  int				        size;						 //number of contained elements (or segments)
  CAnalyticElem   **pElemArray;			 //an array of pointers to Analytic Elements 
  int              *ElemSegment;		 //an array of segment indices
	int               population;			 //number of all contained elements (including those within children)

	int								nCondZones;			 //number of conductivity inhomogeneities
	int								nBaseZones;			 //number of base inhomogeneities
	int								nThickZones;		 //number of thickness inhomogeneities
	int								nPoroZones;			 //number of porosity inhomogeneities
	CPropZone				**pCondZoneArray;  //Array of pointers to conductivity inhomogeneity zones
	CPropZone				**pBaseZoneArray;  //Array of pointers to base inhomogeneity zones
	CPropZone				**pThickZoneArray; //Array of pointers to thickness inhomogeneity zones
	CPropZone				**pPoroZoneArray;  //Array of pointers to porosity inhomogeneity zones

  double          **ElemRHS;				 //[i][m] contains most recent rhs from each element in block
  double           *ElemQ;					 //contains most recent net discharge from each element 

  int               nestlev;				 //level of nesting of this block (0-master)
  int			          quadrant;				 //id of block at nesting level (0-ne, 1-nw, 2-sw, 3-se)
  bool              master;					 //indicator if level 0 block
  bool              leaf;						 //indicator if childless block

  static void       Prepare           ();

  //Element Assignment Functions (for fill superblocks routine, mostly)
  int               AssignToBlock     (CAnalyticElem *Element,const int seg);
  void              AssignToBlock     (CPropZone *PZ);	
  void              AddToBlock        (CAnalyticElem *Elemptr,const int seg);
  void              AddToBlock        (CPropZone *PZ);
  void              AddChild          (CSuperblock *Child, const int index);
  void              AddParent         (CSuperblock *Parent,const int quad);
	void              SetAsLeaf         ();
  void              InitializeBlock   ();

  void              UpdateParent      (const cmplex *CoeffIncrement,const double QIncrement, const int quad);

 public:/*----------------------------------------------------------*/
  //Constructors:
  CSuperblock();
 ~CSuperblock();
  CSuperblock(const cmplex zc, const double Width, const int nestlevel, const int quad);
  

  //Static Member Functions:
  static CSuperblock*
		               FillSuperblocks      (window          Extents  , const int nestlevels, 
	                                       CAnalyticElem **ElemArray, const int numElems,
																				 CPropZone     **PropArray, const int numprops,
																				 const int prec, const bool implicit, const bool nested);
	static void      SetPrecision         (const int Precision,int &order, double &fold);
  static void      Destroy();
	static bool      SuperblocksOn        ();

	//Accessor Functions:
  double           GetWidth             () const;
  int							 GetQuadrant          () const;
  double           GetRadius            () const;
	cmplex           GetCenter            () const;
	int              GetLevel             () const;
  bool             IsLeaf               () const;
  int							 GetPopulation        () const;

	//Inherited from COwnerABC Class (the only visible functions to elements
	bool             IsOn                 () const;
  void             Update               (int BlockID, int segment, double t);

	//Member Functions
	void             SetLaurent           (double &change, double &objective,double t);
  int	             CalcPopulation       ();

  cmplex           GetDischargePotential(const cmplex &z,const double &t) const;
  cmplex           GetW                 (const cmplex &z,const double &t) const;
	cmplex           GetGx                (const cmplex &z,const double &t) const;
	double           GetLeakage           (const cmplex &z,const double &t,const leak_type ltype) const; 
  double           GetNetDischarge      (const double &t) const;

  double		       GetCond              (const cmplex &z) const;
  double		       GetBase              (const cmplex &z) const;
  double		       GetThick             (const cmplex &z) const;
  double		       GetPoro              (const cmplex &z) const;

  cmplex           Centroid             () const;                   
  bool             IsInside             (const cmplex &z) const;   
};

#endif