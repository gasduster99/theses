//RectGrid.h

#ifndef RECTGRID_H
#define RECTGRID_H

#include "MasterInclude.h"
#include "AbstractMesh.h"

class CFluxGrid;
class C2DFDEulerian;
class CChemDomain2D; //TMP DEBUG

/****************************************************
  Class CRectGrid
  Rectangular Finite Difference Grid Data Abstraction
  Arbitrary orientation

Indexing goes from bottom left corner:
	/*       |          |
	 i-1,j+1 | i,j+1    | i+1,j+1  nodeY[j+1]
	_________|__________|_________ gridY[j+1]
			     |          |          
	 i-1,j   | i,j      | i+1,j    nodeY[j]
	_________|__________|_________ gridY[j]
			     |          |
	 i-1,j-1 | i,j-1    | i+1,j-1  nodeY[j-1]
			     |          |          gridY[j-1]
		 |     |    |     |    |
		 |     |    |     |    |
nodeX[i-1] | nodeX[i] | nodeX[i+1]
           |          |
      gridX[i]    gridX[i+1]
*/
/****************************************************/

class CRectGrid:public CMesh {
	friend class CFluxGrid;
	friend class C2DFDEulerian;
  friend class CRectGridp3D;

 private:/*----------------------------------------------------------*/
	double  orient; //angle of orientation
 	cmplex  zg;     //bottom left of grid   

	double *gridX;	//grid lines  (in x-dir) (0 to nX)   //all local coordinates (i.e. gridX[0], gridY[0]=0)
	double *gridY;	//grid lines  (in y-dir) (0 to nY)
	double *nodeX;	//node centers(in x-dir) (0 to nX-1)
	double *nodeY;	//node centers(in y-dir) (0 to nY-1)

	int     nX;			//number of cells in X-dir
	int     nY;			//number of cells in Y-dir

	double  Ymax;		//maximum Y value (in local coords)
	double  Xmax;		//maximum X value (in local coords)

	int     ij_to_k			 (const int i, const int j) const;      //translates from local(i,j) to global(k) indices
	void    k_to_ij			 (const int k, int &i, int &j) const;   //translates from global(k) to local(i,j) indices

	cmplex  GlobalToLocal(const cmplex &z) const;               //translates from global(z) to local(Z) coordinates
	cmplex  LocalToGlobal(const cmplex &Z) const;               //translates from local(Z) to global(z) coordinates

	void    GetInterpolationWeights(const cmplex &Z, double &wx, double &wy, int &i, int &j) const;

 public:/*-----------------------------------------------------------*/
	//Constructors
  CRectGrid(const double north, 
		        const double south,
						const double east, 
						const double west,
						const int    resX,
						const int    resY);
	CRectGrid(const cmplex  zpivot, 
		        const double  angle, 
						const double *X, 
						const double *Y, 	
						const int     numcellsx, 
						const int     numcellsy);
	~CRectGrid();

	//Accessor Functions
	window GetBoundingBox           () const;
  pt3D   GetNodeLocation					(const int k) const;

  double GetCellArea              (const int k) const;
	int    GetNumCellSides          (const int k) const;
	void   GetCellSideVertices      (const int k, const int j, cmplex &z1, cmplex &z2) const;

	//Inquiry functions
  bool   GetCellIndex             (const pt3D &pt, int &k) const; 
	bool   IsInside                 (const pt3D &pt) const;
  bool   IsInCell                 (const pt3D &pt, const int k) const;
	bool   GetMeshBoundaryIntercepts(const pt3D &pt1,const pt3D &pt2, pt3D *ptint, int *k, int &nint) const;
  double InterpolateValue         (const pt3D &pt, Ironclad1DArray values) const;
  void   InterpolateValues        (const pt3D &pt, Ironclad2DArray values, Writeable1DArray intval, const int nval) const;
	double InterpolateDerivative    (const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const;
  void   GetNeighbors             (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const;
  void   DistributePoints         (const int k, distribution_type ty,const int npts, pt3D *pts) const;
	
	static CRectGrid *ReadBBGFile   (ifstream &input,int &l);
	static CRectGrid *ParseSpacing  (ifstream &input,int &l);

	void   WriteGeometry() const;

  void   WriteBBGFile() const;

};
//*********************************************************************
//*********************************************************************
class CRectGridp3D:public CMesh {
 private:/*----------------------------------------------------------*/

	const CRectGrid *p2DGrid;

	double  *gridH;
	double  *nodeH;    //H goes from 0 to 1 (grid surfaces divided at constant H=z/h)
	int     nH;

	double  *top;        //saturated thickness for all k2D=(i,j) locations

	int     k2Dh_to_k			 (const int k2D, const int h) const; //translates from local(i,j,h) to global(k) indices
	void    k_to_k2Dh			 (const int k, int &k2D, int &h) const;   //translates from global(k) to local(i,j,h) indices

	void    GetInterpolationWeights(const pt3D &pt, 
		                              double &wx, double &wy, double &wh, 
																	int    &i,     int &j,     int &h) const;

 public:/*-----------------------------------------------------------*/
	//Constructors
	CRectGridp3D(const CRectGrid *pGrid2D,
							 const double    *H,
							 const int        numcellsh);
	~CRectGridp3D();

	//Accessor Functions
	window GetBoundingBox           () const;
  pt3D   GetNodeLocation					(const int k) const;

  double GetCellArea              (const int k) const;
	int    GetNumCellSides          (const int k) const;
	void   GetCellSideVertices      (const int k, const int j, cmplex &z1, cmplex &z2) const;

	//Inquiry functions
  bool   GetCellIndex             (const pt3D &pt, int &k) const; 
	bool   IsInside                 (const pt3D &pt) const;
  bool   IsInCell                 (const pt3D &pt, const int k) const;
	bool   GetMeshBoundaryIntercepts(const pt3D &pt1,const pt3D &pt2, pt3D *ptint, int *k, int &nint) const;
  double InterpolateValue         (const pt3D &pt, Ironclad1DArray values) const;
  void   InterpolateValues        (const pt3D &pt, Ironclad2DArray values, Writeable1DArray intval, const int nval) const;
	double InterpolateDerivative    (const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const;
  void   GetNeighbors             (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const;
  void   DistributePoints         (const int k, distribution_type ty,const int npts, pt3D *pts) const;
	
	static CRectGridp3D *ReadBBGFile   (ifstream &input,int &l);

	void   WriteGeometry() const;
  void   WriteBBGFile() const;
};

#endif 