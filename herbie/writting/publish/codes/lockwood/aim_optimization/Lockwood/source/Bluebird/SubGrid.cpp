//SubGrid.cpp
//RectGrid.cpp
#include "SubGrid.h"

/*******************************************************************************
           CONSTRUCTORS
*******************************************************************************/
CSubGrid::CSubGrid(){
	pMasterGrid=NULL;
	li=NULL;
	mi=NULL;
	BBox.n=0;   BBox.s=0;
	BBox.e=0;   BBox.w=0;
	nCells=0;   nNodes=0;
	nMasterCells=0;
}
//---------------------------------------------------------------------------
CSubGrid::CSubGrid(CMesh *master,
									 const int   *indices,
									 const int    numcells){
	pMasterGrid=master;
	nMasterCells=pMasterGrid->GetNumNodes();
	BBox        =pMasterGrid->GetBoundingBox();
	nNodes      =numcells;
	meshtype    =pMasterGrid->GetType();
	if (meshtype==CELL_BASED){nCells=nNodes;}
	else {ExitGracefully("CSubGrid constructor: TMP DEBUG",RUNTIME_ERR);}

	li=new int [nMasterCells];
	mi=new int [nCells];
	
	int count(0);
	for (int i=0;i<nMasterCells;i++){
		li[i]=OFF;
		if ((indices[count])==i){
			mi[count]=i;
			li[i]=count;
			count++;
		}
	}
}

//---------------------------------------------------------------------------
CSubGrid::~CSubGrid(){
	if (globaldebug){cout<<"  DESTROYING SUBGRID"<<endl;}
	delete [] li;
	delete [] mi;
}
/*******************************************************************************
           ACCESSOR FUNCTIONS 
*******************************************************************************/
pt3D CSubGrid::GetNodeLocation(const int k) const{
	if ((k<0) || (k>nCells)){ExitGracefully("CSubGrid::GetNodeLocation: Bad Cell index",RUNTIME_ERR);}
	return pMasterGrid->GetNodeLocation(mi[k]);
}
//---------------------------------------------------------------------------
double CSubGrid::GetCellArea(const int k) const{
	if ((k<0) || (k>nCells)){ExitGracefully("CSubGrid::GetCellArea: Bad Cell index",RUNTIME_ERR);}
	return pMasterGrid->GetCellArea(mi[k]);
}
//---------------------------------------------------------------------------
window CSubGrid::GetBoundingBox()      const{return BBox;}

/*******************************************************************************
           INQUIRY FUNCTIONS
*******************************************************************************/
/****************************************************************************
			Get Cell Index
-----------------------------------------------------------------------------
 output: k (OFF if outside, local cell index otherwise)
				 bool (true if inside)
****************************************************************************/
bool CSubGrid::GetCellIndex(const pt3D &pt, int &k) const{ 
	int masterind(OFF);
	if (pMasterGrid->GetCellIndex(pt,masterind)){ 
		if (masterind!=OFF){
			for (int i=0; i<nCells; i++){
				if (mi[i]==masterind){
					k=i;
					return true;
				}
			}
		}
	}
	k=OFF;
	return false;
}
//---------------------------------------------------------------------------
bool CSubGrid::IsInside(const pt3D &pt) const{
	for (int k=0; k<nCells; k++){
		if (IsInCell(pt,k)){return true;}
	}
	return false;
}
//---------------------------------------------------------------------------
bool CSubGrid::IsInCell(const pt3D &pt, const int k) const{
	if ((k<0) || (k>nCells)){ExitGracefully("CSubGrid::IsInCell: Bad Cell index",RUNTIME_ERR);}
	return pMasterGrid->IsInCell(pt,mi[k]);
}
/***************************************************************************
                  INTERPOLATE VALUE(S) 
****************************************************************************
input:   z (global coordinates) 
				 nCell-sized array(s) of values to be interpolated, indexed by k (i.e. concentration)
				 InterpolateValues has additional parameter of the number of arrays (nval)
output:  interpolated value(s) ((*intval) for InterpolateValues function)
---------------------------------------------------------------------------*/
double CSubGrid::InterpolateValue(const pt3D &pt,Ironclad1DArray values) const{
	//weighted average of closest values
	if (IsInside(pt)){
		//use getneighbors function, interpolate based upon inverse distance
		int    neigh[MAX_CELL_NEIGHBORS];
		double dist [MAX_CELL_NEIGHBORS];
		int    nn(MAX_CELL_NEIGHBORS),k;
		double sum1(0.0),sum2(0.0),power(2); 
		pt3D   cpt;

		GetCellIndex(pt,k);
		cpt=GetNodeLocation(k);
		GetNeighbors(k,neigh,dist,nn);

		for (int i=0;i<nn;i++){
			sum1+=values[mi[neigh[i]]]/pow(dist[i],power);
			sum2+=1.0/pow(dist[i],power);
		}
		sum1+=values[mi[k]]/pow(abs(cpt-pt),power);
		sum2+=1.0/pow(abs(cpt-pt),power);		

		return sum1/sum2;
	}
	else { //not on grid, cannot interpolate
		return 0.0;
	}
}
//---------------------------------------------------------------------------
void CSubGrid::InterpolateValues(const pt3D &pt, Ironclad2DArray values, Writeable1DArray intval, const int nval) const{
	//average of closest 4 values
	int n;
	if ((values==NULL) || (intval==NULL)){ExitGracefully("CSubGrid::InterpolateValues: cannot interpolate on NULL value array",RUNTIME_ERR);}
	if (IsInside(pt)){
		//use getneighbors function, interpolate based upon inverse distance
		int    neigh[MAX_CELL_NEIGHBORS];
		double dist [MAX_CELL_NEIGHBORS];
		int    nn(MAX_CELL_NEIGHBORS),k;
		double sum1(0.0),sum2(0.0),power(2); 
		pt3D   cpt;

		GetCellIndex(pt,k);
		cpt=GetNodeLocation(k);
		GetNeighbors(k,neigh,dist,nn);
		for (int i=0;i<nn;i++){
			sum2+=1.0/pow(dist[i],power);
		}
		sum2+=1.0/pow(abs(cpt-pt),power);

		for (int n=0;n<nval;n++){
			for (int i=0;i<nn;i++){
				sum1+=values[mi[neigh[i]]][n]/pow(dist[i],power);
			}
			sum1+=values[mi[k]][n]/pow(abs(cpt-pt),power);
			intval[n]=sum1/sum2;
		}
	}
	else { //not on grid, cannot interpolate
    for (n=0;n<nval;n++){intval[n]=0.0;}
	}
}
/***************************************************************************
                  GET NEIGHBORS
****************************************************************************
input:   cell index (k)
				 empty integer array to store neighbor values (pre-allocated)
				 empty double array to store distance to neighbors (pre-allocated)
				 maximum number of neighbors (numneighbors)
output:  neighbor indices (neighbors), numneighbors (maximum 4 for grid, 10 for mesh)
--------------------------------------------------------------------------*/
void CSubGrid::GetNeighbors  (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const{
	int    *locneigh;
	double *locdist;
	int     i,count(0);
	pMasterGrid->GetNeighbors(mi[k],neighbors,dist,numneighbors);
	locneigh=new int    [numneighbors];
	locdist =new double [numneighbors];

	for (i=0; i<numneighbors;i++){ //sift through master neighbors, include only those in subgrid
		if (li[neighbors[i]]!=OFF){ //if neighbor index part of subgrid
			locneigh[count]=li[neighbors[i]];
			locdist [count]=locdist[i];
			count++;
		}	
	}
	numneighbors=count;
	for (i=0;i<numneighbors;i++){
		neighbors[i]=locneigh[i];
		dist[i]     =locdist[i];
	}
	delete [] locneigh;
	delete [] locdist;
}
/****************************************************************************
                  DISTRIBUTE POINTS 
*****************************************************************************
input:   cell index (k)
				 distribution type (ty) :random or fixed pattern
				 empty cmplex array (pts) to store (npt) point locations
output:  cmplex array of locations (pts)
---------------------------------------------------------------------------*/
void CSubGrid::DistributePoints(const int k, distribution_type ty,const int npts, pt3D *pts) const{
	pMasterGrid->DistributePoints(k,ty,npts,pts);

}
/****************************************************************************
                  Get Mesh Boundary Intercepts
***************************************************************************
input:   endpoints pt1 and pt2 of a line
         max array size of ptint,i,and j: nint
output:  intercepts of line and grid faces (zint), sorted from z1 to z2
         i and j indices of all cells crossed (nint+1)
         nint, the total number of intersections
returns: true if both points are in grid
---------------------------------------------------------------------------*/
bool CSubGrid::GetMeshBoundaryIntercepts(const pt3D &pt1, const pt3D &pt2, pt3D *ptint, int *k, int &nint) const{


	//TOUGH COOKIE

	return false;
}
/***************************************************************************
                  WRITE GEOMETRY
****************************************************************************
 Writes grid to Atlas BNA file, rectgrid.bna
--------------------------------------------------------------------------*/
void CSubGrid::WriteGeometry() const{

}
