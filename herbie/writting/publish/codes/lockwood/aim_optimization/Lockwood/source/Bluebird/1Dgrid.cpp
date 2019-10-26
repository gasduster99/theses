#include "1Dgrid.h"

/*********************************************************************
                  GENERIC 1D Finite Difference Grid
**********************************************************************
                           CONSTRUCTORS
**********************************************************************/
C1DFDGrid::C1DFDGrid(double *X,int numcells){

	int k;

	nX    =numcells+1;
	nCells=numcells;
	
	gridX=new double [numcells+1];
	nodeX=new double [numcells  ];

	minX=X[0];
	gridX[0]=0.0;
	for (k=1; k<nX; k++){
		gridX[k]=X[k]-minX;
		nodeX[k]=0.5*(gridX[k]+gridX[k-1]);
	}
}
/*********************************************************************
                           ACCESSORS
**********************************************************************/
double C1DFDGrid::GetCellLength  (const int k) const{return gridX[k+1]-gridX[k];}
//--------------------------------------------------------------------
double C1DFDGrid::GetNodeLocation(const int k) const{return nodeX[k]+minX;}//cell center
//--------------------------------------------------------------------
bool   C1DFDGrid::IsInCell       (const double &x, const int k) const{
	return ((x-minX>=gridX[k]) && (x-minX>=gridX[k+1]));
}
//--------------------------------------------------------------------
bool   C1DFDGrid::GetCellIndex             (const double &x, int &k) const{


	if (!IsInside(x)){k=OFF; return false;}

	k=(int)((x-minX)/gridX[nX]);
	while ((k<nX-1) && (x-minX>gridX[k+1])){k++;}
	while ((k>0   ) && (x-minX<gridX[k  ])){k--;}

	return true;
}
//--------------------------------------------------------------------
bool   C1DFDGrid::IsInside(const double &x) const{
 return ((x>=minX) && (x-minX<gridX[nCells+1]));
}
/*********************************************************************
                           MEMBER FUNCTIONS
**********************************************************************/
double C1DFDGrid::InterpolateValue  (const double &x, Ironclad1DArray values) const{
	//finite difference- interpolate from cell centers
	double frac;
	int k;
	GetCellIndex(x,k);
	if (x==nodeX[k])
	{
		return values[k];
	}
	else if (x>nodeX[k])
	{
		if (k==nX-1){return values[nX-1];}
		frac=(x-nodeX[k])/(nodeX[k+1]-nodeX[k]);
		return values[k]+frac*(values[k+1]-values[k]);
	}
	else
	{
		if (k==0){return values[0];}
		frac=(x-nodeX[k-1])/(nodeX[k]-nodeX[k-1]);
		return values[k-1]+frac*(values[k]-values[k-1]);
	}
}
//--------------------------------------------------------------------
void C1DFDGrid::InterpolateValues(const double &x, Ironclad2DArray values, Writeable1DArray intval, const int nval) const{
	//finite difference- interpolate from cell centers
	double frac;
	int i,k;
	GetCellIndex(x,k);
	if (x==nodeX[k])
	{
		for (i=0; i<nval; i++){intval[i]=values[k][i];}
	}
	else if (x>nodeX[k])
	{
		if (k==nX-1){
			for (i=0; i<nval; i++){intval[i]=values[nX-1][i];}
		}
		frac=(x-nodeX[k])/(nodeX[k+1]-nodeX[k]);
		for (i=0; i<nval; i++){
			intval[i]=values[k][i]+frac*(values[k+1][i]-values[k][i]);
		}
	}
	else
	{
		if (k==0){
			for (i=0; i<nval; i++){intval[i]=values[0][i];}
		}
		frac=(x-nodeX[k-1])/(nodeX[k]-nodeX[k-1]);
		for (i=0; i<nval; i++){
			intval[i]=values[k-1][i]+frac*(values[k][i]-values[k-1][i]);
		}
	}
}
//--------------------------------------------------------------------
void   C1DFDGrid::WriteGeometry            (){
}
