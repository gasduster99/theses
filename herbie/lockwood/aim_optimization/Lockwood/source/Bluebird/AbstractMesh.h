//Abstract Mesh Header

#ifndef ABSTRACTMESH_H
#define ABSTRACTMESH_H

#include "MasterInclude.h"

const int OFF               =-1;       // element/node/side is switched off 
const int MAX_CELL_NEIGHBORS=12;       //maxmimum number of neighbors for a grid cell

/****************************************************
  Class CMesh
  Purely virtual Grid or Mesh Data Abstraction
****************************************************/
enum distribution_type  {FIXED_PATTERN,RANDOM_PATTERN};
enum mesh_rep           {NODE_BASED,CELL_BASED};

class CMesh {

 protected:/*----------------------------------------------------------*/

	bool     is3D;      //true if is 3D
	window   BBox;		  //bounding box (global coords) around mesh
	polygon  bounds;    //boundary polygon
	int      nCells;    //number of cells (for cell based methods) 
	int      nNodes;    //number of nodes (degrees of Freedom; valid for node & cell based meshes)
	mesh_rep meshtype;  //type of representation (cell-based or node-based)

 public:/*-----------------------------------------------------------*/

	//Constructors 
	CMesh::CMesh(){nCells=0;nNodes=0;is3D=true;bounds.nsides=0; bounds.vertices=NULL;meshtype=CELL_BASED;}

	virtual        CMesh::~CMesh(){delete [] bounds.vertices;}

	//Accessor Functions
	window         CMesh::GetBoundingBox           () const  {return BBox;}
	polygon        CMesh::GetBoundary              () const  {return bounds;}
  mesh_rep       CMesh::GetType                  () const  {return meshtype;}
	int            CMesh::GetNumNodes              () const  {return nNodes;}

  //for node-based meshes

	//For cell-based meshes (Cells must always be convex polygons):
  int            CMesh::GetNumCells              () const  {if (meshtype==CELL_BASED){return nCells;}else {ExitGracefully("CMesh:Requested number of Cells for non-cell based method",RUNTIME_ERR);return 0;}}
  virtual double CMesh::GetCellArea              (const int k) const=0; 
	virtual	pt3D	 CMesh::GetNodeLocation          (const int k) const=0; 
  virtual int    CMesh::GetNumCellSides          (const int k) const=0; 
	virtual void   CMesh::GetCellSideVertices      (const int k, const int j, cmplex &z1, cmplex &z2) const=0; 
  virtual bool   CMesh::IsInCell                 (const pt3D &pt, const int k) const=0;
  virtual bool   CMesh::GetCellIndex             (const pt3D &pt, int &k) const=0; 
  virtual void   CMesh::DistributePoints         (const int k, distribution_type ty,const int npts, pt3D *pts) const=0;

	//inquiry functions
	virtual bool   CMesh::IsInside                 (const pt3D &pt) const=0;
  virtual double CMesh::InterpolateValue         (const pt3D &pt, Ironclad1DArray values) const=0;
  virtual void   CMesh::InterpolateValues        (const pt3D &pt, Ironclad2DArray values, Writeable1DArray intval, const int nval) const=0;
	virtual double CMesh::InterpolateDerivative    (const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const=0;
	virtual void   CMesh::GetNeighbors             (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const=0;
	virtual bool   CMesh::GetMeshBoundaryIntercepts(const pt3D &pt1, const pt3D &pt2, pt3D *ptint, int *k, int &nint) const=0;
	
	//virtual void   CMesh::WriteGeometry          () const=0;

};

#endif