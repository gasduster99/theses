//SubGrid.h

#ifndef SUBGRID_H
#define SUBGRID_H

#include "MasterInclude.h"
#include "AbstractMesh.h"

/****************************************************
  Class CSubGrid
  SubGrid Data Abstraction
****************************************************/

class CSubGrid:public CMesh {

 private:/*----------------------------------------------------------*/
	CMesh *pMasterGrid; //

	int     *li;          //contains correspinding master grid indices (indexed by subgrid indices)       
	int     *mi;          //contains corresponding local grid indices (indexed by master grid indices)
	int      nMasterCells;

 public:/*-----------------------------------------------------------*/
	//Constructors
	CSubGrid();
  CSubGrid(CMesh     *master,
		       const int *indices,
					 const int  numcells);
	~CSubGrid();

	//Manipulator Functions
	void	 SetCellHeight            (const int k, const double &h);

	//Accessor Functions
	window GetBoundingBox           () const;
  pt3D   GetNodeLocation					(const int k) const;
  double GetCellArea              (const int k) const; 
	double GetCellVolume						(const int k) const;

	//Inquiry functions
  bool   GetCellIndex             (const pt3D &pt, int &k) const; 
	bool   IsInside                 (const pt3D &pt) const;
  bool   IsInCell                 (const pt3D &pt, const int k) const;
	bool   GetMeshBoundaryIntercepts(const pt3D &pt1,const pt3D &pt2, pt3D *ptint, int *k, int &nint) const;
  double InterpolateValue         (const pt3D &pt, Ironclad1DArray values) const;
  void   InterpolateValues        (const pt3D &pt, Ironclad2DArray values, Writeable1DArray intval, const int nval) const;
  void   GetNeighbors             (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const;
  void   DistributePoints         (const int k, distribution_type ty,const int npts, pt3D *pts) const;

	void   WriteGeometry() const;

};

#endif 