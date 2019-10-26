//TriagonalMesh.h
#ifndef TRIMESH_H
#define TRIMESH_H

#include "MasterInclude.h"
#include "AbstractMesh.h"

//************************************************************************
const int    MAX_NODES      =15000;
const int    MAX_DYN_PTS		=500;
const int    MAX_VORONOI_PTS=200;
const bool   DOUBLENODETEST =true;

enum elem_state  {ACTIVE,DONE,WAITING};
enum chain_type  {CLOSED,OPEN,INSIDE};
enum status_type {BOUNDARY, INTERIOR, IS_OFF,DOUBLENODE,DNODE_TEMP};
enum dnode_type  {DNODE_JUMP,DNODE_FLUX,DNODE_BOUNDARY};
enum spc_type    {INTERP, COURANT,USER_SPECIFIED};
//************************************************************************

/******************************************************************
  Class CTriMesh
  Triangular Mesh Data Abstraction
******************************************************************/
class CTriMesh:public CMesh{
 protected:/*----------------------------------------------------*/

	//Static Functions-------------------------------------------------------
	static status_type TranslateMark(const int mark);

	//Protected structures---------------------------------------------------
	struct elemstruct{
		int ni, nj, nk;       //node indices
		int ei, ej, ek;       //adjacent element indices
		int si, sj, sk;       //adjacent side indices
		cmplex		 zv;        //center of circumscribed circle							
		double		 R;         //radius of circumscribed circle
		cmplex		 zin;				//center of inscribed circle (never really used)
		double		 r;					//radius of inscribed circle 
		bool       on;
		double     area;      //area of triangle
		cmplex     centroid;  //centroid of triangle
		elem_state state;     //Done, Active or Waiting 
	};

	struct sidestruct{
		int         ea, eb;   //left and right element 
		int         na, nb;		//left and right node
		int         n1, n2;		//start and end node 
		int         A, B;     //coefficients of half plane assoc with side determined by Ax+By=1.0
		status_type status;   //on boundary, interior, turned off, or doublenode
		double      radius;   //radius of curved sides (=0 for non-curved)
	};
	
	struct nodestruct{  
		cmplex      z;        //location
		double      F;        //Spacing function value
		status_type status;   //on boundary, interior, turned off, or doublenode
		cmplex      inv;      //doublenode normal vector (direction of split)
	};

	struct segstruct{
		int         p0,p1;    //start and end point indices
		int         N;        //number of nodes along boundary
		status_type status;   //on boundary, interior or turned off
		double      radius;   //radius of curved segment (=0 for non-curved)
	};
	
	struct chainstruct{
		int        s0, s1;    //start segment
		chain_type type;      //type
	};

	struct dblestruct {
		int        n1, n2;    //colocated nodes
		double     width;     //"distance" between nodes
		int        nsegs;     //number of interstitial nodes
		dnode_type type;      //type
	};
	
	struct spcstruct{
		spc_type  type;
		double    time_step;
		double    min_spacing;
		double    max_spacing;
		double    radius;
	};

	//Member Variables-------------------------------------------------------
	char        *name;      //Name of Mesh

  elemstruct  *elem;			//array of elements
	int					 NumElems;	//number of elements
  sidestruct  *side;			//array of sides
	int					 NumSides;	//number of Sides
	nodestruct  *node;			//array of nodes
	int					 NumNodes;	//number of nodes
	dblestruct	*dnodes;    //array of doublenode pairs
	int          NumDoublenodes; //number of doublenode pairs
	int          nIntNodes; //number of internal nodes (NumDoublenodes x num internal segments per dblenode)

	int					 MaxNodes;	//Maximum Number of nodes

	int					 ugly;			//index of current ugly node (dynamic creation only)      
	double       minSpacing;//minimum spacing between nodes

	cmplex      *sumz;       //sum of adjacent locations [NumNodes]
	int         *Nne;        //number of node neighbors which are interior nodes [NumNodes]
	double      *sumF;       //sum of adjacent weighting functions [NumNodes]
  spcstruct    spacinginfo;//dynamic creation information

	//Function Declarations (protected)-mostly for creation of mesh----------

	virtual void Process   (){return;} //ExitGracefully("CTriMesh::Process",VIRTUAL_ERROR);   

  void	 Initialize      ();
	void	 CleanIndices    ();
	void	 Renumerate      ();

	double CalculateSpacingFunction(const int e, const cmplex &z) const;
	void   CalculateCircles				 (const int e);
	void	 SwapDiagonal            (const int s);

	int		 InsertNode							 (const cmplex z,  
																	const status_type status,
																	const cmplex &invector); 
	void	 UpdateNewSideStatus     (const int         adjacent_n, 
																	const status_type adjacent_status,
																  const double      radius);
	void	 DynamicClassify				 ();
	bool	 DynamicNewNode					 ();
	int		 DynamicGetDoneSide      (const int ug, int &n_across) const;
	void	 DynamicSmooth					 (bool relax, bool mild);
	void	 DynamicFix              (); 

  void	 DoubleNodeFilter        ();
	void	 DoubleNodeFilter2       ();//Not operational
  int		 DoubleNodeBreak         (const int n, const int s1,const int s2);
	bool	 WingInternal            (const cmplex z,const int s1, const int dn, const int s2) const;

	int		 GetElement							 (const cmplex &z) const;
	void	 GetInternalNodes				 (cmplex *&zpoints,int &ncontrol) const;

 public:/*----------------------------------------------------*/
  
	//Static Functions--------------------------------------------------------
	static CTriMesh *DynamicParse  (ifstream &input, int &l,char *Name);
	static CTriMesh *Parse         (ifstream &input, int &l,char *Name);

	static bool CreatePolygonControlPoints(const cmplex *ze,
		                                     const int     NLines,
																				 cmplex      *&zpoints,
																				 int          &ncontrol);

	//Constructors------------------------------------------------------------
	CTriMesh(char              *Name);
	CTriMesh(char              *Name,
		       const nodestruct  *points,		       
					 const int          nP,
					 const segstruct   *segments,
	         const int          nS,
					 const int          maximum_nodes,
					 const spcstruct    spacing_package);

	CTriMesh(char              *Name,
		       const cmplex      *pts,
					 const int          numpts);
 ~CTriMesh();

	//Public functions---------------------------------------------------------

  void   DynamicInitialize       (spc_type  type, double param, double minx, double maxx, double radius);
 	void	 DynamicBuildMesh				 (bool relax_on,bool smooth_on);
	void   WriteOutput();
	void   WriteOutput2();

	//Functions inherited from CMesh (not virtual)
	bool   IsInside                 (const pt3D &pt) const;   
	void   WriteGeometry            (bool DXFon, bool BNAon) const;	

	//Functions inherited from CMesh (all virtual-definition of cell/cell bounds/neighbors differs)
	pt3D   GetNodeLocation					(const int k) const                   {ExitGracefully("CTriMesh::GetNodeLocation",VIRTUAL_ERROR);return 0;}          
 
	double GetCellArea              (const int k) const                   {ExitGracefully("CTriMesh::GetCellArea",VIRTUAL_ERROR);return 0;}
  bool   GetCellIndex             (const pt3D &pt, int &k) const        {ExitGracefully("CTriMesh::GetCellIndex",VIRTUAL_ERROR);return false;}                                           
  bool   IsInCell                 (const pt3D &pt, const int k) const   {ExitGracefully("CTriMesh::IsInCell",VIRTUAL_ERROR);return false;}                                   
	bool   GetMeshBoundaryIntercepts(const pt3D &pt1,const pt3D &pt2, 
		                               pt3D *ptint, int *k, int &nint) const{ExitGracefully("CTriMesh::GetMeshBoundaryIntercepts",VIRTUAL_ERROR);return false;}   
  double InterpolateValue         (const pt3D &pt, 
		                               Ironclad1DArray values) const        {ExitGracefully("CTriMesh::InterpolateValue",VIRTUAL_ERROR);return 0;}
  void   InterpolateValues        (const pt3D &pt, 
		                               Ironclad2DArray values, 
																	 Writeable1DArray intval, 
																	 const int nval) const                {ExitGracefully("CTriMesh::InterpolateValues",VIRTUAL_ERROR);return;}
	double InterpolateDerivative    (const pt3D &pt, 
																	 Ironclad1DArray values, 
																	 const cmplex &dir) const             {ExitGracefully("CTriMesh::InterpolateDerivatives",VIRTUAL_ERROR);return 0.0;}
	void   GetNeighbors             (const int k, 
		                               int *neighbors, 
																	 Writeable1DArray dist, 
																	 int &numneighbors) const             {ExitGracefully("CTriMesh::GetNeighbors",VIRTUAL_ERROR);return;}
  int    GetNumCellSides          (const int k) const                   {ExitGracefully("CTriMesh::GetNumCellSides",VIRTUAL_ERROR);return 0;}
	void   GetCellSideVertices      (const int k, const int j, cmplex &z1, cmplex &z2) const{ExitGracefully("CTriMesh::GetCellSideVertices",VIRTUAL_ERROR);return;}
  void   DistributePoints         (const int k, distribution_type ty,const int npts, pt3D *pts) const{ExitGracefully("CTriMesh::DistributePoints",VIRTUAL_ERROR);return;}


	bool   ReadItself (ifstream &NODE, ifstream &ELEM, ifstream &SIDE, ifstream &DNOD);
	bool   ReadItself2(ifstream &BTM);

}; //end class CTriMesh

/******************************************************************
  Class CFEMesh
  Finite Element Mesh Data Abstraction
	Uses only *active* elements from TriMesh
******************************************************************/

class C2DFEEulerian;

class CFEMesh:public CTriMesh{
  friend class C2DFEEulerian;

 protected:/*----------------------------------------------------*/
	
	void              Process();

  double						GetN   (const int e, const cmplex &z, const int i) const;
	cmplex            GetdN  (const int e, const cmplex &z, const int i) const;

	cmplex  elemLocalToGlobal(const int e, const double &xii, const double &xij, const double &xik) const;
  void    elemGlobalToLocal(const int e, const cmplex &z,double &xii,double &xij, double &xik) const;

	void   GetElemJacobian   (const int e, const double &xii, const double &xij, const double &xik,double J[2][2],double Jinv[2][2], double &Jdet)const;
  double GetElemJacobianDet(const int e, const cmplex &z) const;

	void          Get_d_coeff(const cmplex Z[3], const double R, double d[3]) const;
  double GetRadiusAndRelativeIndices(const int e, int ind[3]) const;

 public:/*----------------------------------------------------*/
  
	//Constructors
	CFEMesh(char              *Name);
	CFEMesh(char              *Name,
		      const nodestruct  *points,
					const int          nP,
					const segstruct   *segments,
					const int          nS,
					const int          maximum_nodes,
					const spcstruct    spacing_package);
	CFEMesh(char              *Name,
		      const cmplex      *pts,
					const int          numpts);
 ~CFEMesh();

	//public functions
	void WriteOutput();

	//functions inherited from CMesh 
	pt3D   GetNodeLocation          (const int k) const;
	double InterpolateValue         (const pt3D &pt, Ironclad1DArray values) const;
  void   InterpolateValues        (const pt3D &pt, Ironclad2DArray values, Writeable1DArray intval, const int nval) const;
	double InterpolateDerivative    (const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const;
	void   GetNeighbors             (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const;
	
	double FullIntegrate            (Ironclad1DArray values, const cmplex z1, const cmplex z2, const cmplex z3) const;

	double IntegrateFunction        (const int e, const trigausspoints gausspts, double (*F)(const cmplex &z)) const;



	//Cell functions (should not be used) - should add VIRTUAL ERRORS
	double GetCellArea              (const int k) const{return 0.0;}
	int    GetNumCellSides          (const int k) const{return 0;} 
	void   GetCellSideVertices      (const int k, const int j, cmplex &z1, cmplex &z2) const{return;} 
  bool   IsInCell                 (const pt3D &pt, const int k) const{return false;}
  bool   GetCellIndex             (const pt3D &pt, int &k) const{return OFF;}  
	bool   GetMeshBoundaryIntercepts(const pt3D &pt1,const pt3D &pt2, pt3D *ptint, int *k, int &nint) const{return false;}
  void   DistributePoints         (const int k, distribution_type ty,const int npts, pt3D *pts) const{return;}
}; //end class CFEMesh

/******************************************************************
  Class CVoronoiMesh
  Voronoi Mesh Data Abstraction
******************************************************************/
class CVoronoiMesh:public CTriMesh{
 protected:/*----------------------------------------------------*/
	struct vorstruct{
		cmplex  z;           //center node (not neccesarily a node)
		cmplex *zv;          //vertices of voronoi region 
		int    *e;           //elements in voronoi region
		int    *n;           //adjacent nodes
		int     nVert;       //number of vertices
		double  area;        //area of cell 
	};

	vorstruct   *voronoi;   //array of voronoi cells (postprocessing: one per node)

	void Process   ();

	void CalcVoronoi             (vorstruct &vor,const cmplex &zp, const int n);

 public:/*----------------------------------------------------*/
  
	//Constructors
	CVoronoiMesh(char              *Name,
				       const nodestruct  *points,
							 const int          nP,
							 const segstruct   *segments,
							 const int          nS,
							 const int          maximum_nodes,
					     const spcstruct    spacing_package);
	CVoronoiMesh(char              *Name,
		           const cmplex      *pts,
					     const int          numpts);
 ~CVoronoiMesh();

	//public functions
	void WriteOutput();

	//functions inherited from CMesh 
	pt3D   GetNodeLocation					(const int k) const;            
  double GetCellArea              (const int k) const;
  bool   GetCellIndex             (const pt3D &pt, int &k) const;                                           
	bool   IsInCell                 (const pt3D &pt, const int k) const;                                     
	bool   GetMeshBoundaryIntercepts(const pt3D &pt1,const pt3D &pt2, pt3D *ptint, int *k, int &nint) const;
  double InterpolateValue         (const pt3D &pt, Ironclad1DArray values) const;
  void   InterpolateValues        (const pt3D &pt, Ironclad2DArray values, Writeable1DArray intval, const int nval) const;
	void   GetNeighbors             (const int k, int *neighbors, double *dist, int &numneighbors) const;
	
	void   WriteGeometry            () const;

}; //end class CVoronoiMesh

#endif
