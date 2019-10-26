//One-dimensional grid Header

#ifndef ONEDGRID_H
#define ONEDGRID_H

#include "MasterInclude.h"
#include "AbstractMesh.h" //for definition of OFF
/****************************************************
  Class C1DFDGrid
  1-dimensional Grid/mesh abstraction
****************************************************/

class C1DFDGrid {

 protected:/*----------------------------------------------------------*/

	double   minX;
	double  *nodeX;     //nodal coordinates 
	double  *gridX;     //cell/element center coordinates
	int      nCells;    //number of cells (for cell based methods) 
	int      nX;        //number of divisions
	

 public:/*-----------------------------------------------------------*/

	//Constructors 
	C1DFDGrid::C1DFDGrid(){nCells=0;nX=0;minX=0.0;nodeX=NULL; gridX=NULL;}
	C1DFDGrid::C1DFDGrid(double *X,int nX);

	C1DFDGrid::~C1DFDGrid(){delete [] nodeX; delete [] gridX;}

	//Accessor Functions
  int    C1DFDGrid::GetNumCells              () const  {return nCells;}

  double C1DFDGrid::GetCellLength            (const int k) const;
	double C1DFDGrid::GetNodeLocation          (const int k) const;                                             
	bool   C1DFDGrid::IsInCell                 (const double &x, const int k) const;
  bool   C1DFDGrid::GetCellIndex             (const double &x, int &k) const;

	//inquiry functions
	bool   C1DFDGrid::IsInside                 (const double &x) const;

	double C1DFDGrid::InterpolateValue         (const double &x, Ironclad1DArray values) const;               
  void   C1DFDGrid::InterpolateValues        (const double &x, Ironclad2DArray values, Writeable1DArray intval, const int nval) const;
	//double C1DFDGrid::InterpolateDerivative    (const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const;
	
	void   C1DFDGrid::WriteGeometry            ();

};

#endif