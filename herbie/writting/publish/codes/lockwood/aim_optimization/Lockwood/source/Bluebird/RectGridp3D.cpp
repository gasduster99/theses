//RectGrid.cpp
#include "RectGrid.h"

/*******************************************************************************
           CONSTRUCTORS
*******************************************************************************/
CRectGridp3D::CRectGridp3D(const CRectGrid *p2D,
													 const double    *H,
													 const int        numcellsh){
	int i;

	p2DGrid=p2D;
	nH=numcellsh;
	nCells=nNodes=p2DGrid->GetNumCells()*nH;
	top   =new double [p2DGrid->GetNumCells()];
	gridH	=new double [nH+1];
	nodeH	=new double [nH];
	if (nodeH==NULL){ExitGracefully("CRectGridp3D::CRectGridp3D: out of memory", OUT_OF_MEMORY);}

	for (i=0; i<=nH; i++){gridH[i]=H[i];}
	for (i=0; i< nH; i++){nodeH[i]=(gridH[i]+gridH[i+1])/2.0;}
  for (i=0; i< p2DGrid->GetNumCells(); i++){
		top[i]=10.0; //default grid thickness 
	}
	is3D=true;
	meshtype=CELL_BASED;
}
//---------------------------------------------------------------------------
CRectGridp3D::~CRectGridp3D(){
	if (globaldebug){cout<<"  DESTROYING PSEUDO-3D RECTGRID"<<endl;}
	delete p2DGrid;
	delete [] gridH;
	delete [] nodeH;
	delete [] top;

}
/*******************************************************************************
           PRIVATE CONVERSION FUNCTIONS 
//SHOULD BE OPTIMIZED- could be inlined
*******************************************************************************/
int CRectGridp3D::k2Dh_to_k			 (const int k2D, const int h) const{
 //translates from local(i,j,h) to global(k) indices
	if ((k2D==OFF) ||(k2D>=p2DGrid->nCells)  || (k2D<0) || 
			(h  ==OFF) ||(h  >=nH             )  || (h  <0)) {return OFF;}
	//if ((j*nX+i>=nCells) || (j*nX+i<0)){ExitGracefully("CRectGridp3D::ij_to_k: bad input", RUNTIME_ERR);}
  
	return h*(p2DGrid->nCells)+k2D;
}
//------------------------------------------------------------------------------
void CRectGridp3D::k_to_k2Dh			 (const int k, int &k2D, int &h) const{
	//translates from global(k) to local(i,j,h) indices
  if ((k==OFF) || (k>=nCells)){k2D=OFF;h=OFF;return;}
  k2D=k%(p2DGrid->nCells);
	h=(k-k2D)/(p2DGrid->nCells);
}
/*******************************************************************************
           ACCESSOR FUNCTIONS 
*******************************************************************************/
//TMP DEBUG: presently only works with base of zero
pt3D CRectGridp3D::GetNodeLocation(const int k) const{
	ExitGracefullyIf((k<0) || (k>nCells),"CRectGridp3D::GetNodeLocation: Bad Cell index",RUNTIME_ERR);
	int k2D,h;
	pt3D pt;
	k_to_k2Dh(k,k2D,h);
  pt=p2DGrid->GetNodeLocation(k2D); //gets x,y location
  pt.z=nodeH[h]*top[k2D];             //sets z location
	return pt; 
}
double CRectGridp3D::GetCellArea        (const int k) const {return p2DGrid->GetCellArea(k);}
//---------------------------------------------------------------------------
int    CRectGridp3D::GetNumCellSides    (const int k) const {return p2DGrid->GetNumCellSides(k);}
//---------------------------------------------------------------------------
window CRectGridp3D::GetBoundingBox     ()            const {return p2DGrid->BBox;}
//---------------------------------------------------------------------------
void   CRectGridp3D::GetCellSideVertices(const int k, const int s, cmplex &z1, cmplex &z2) const{p2DGrid->GetCellSideVertices(k,s,z1,z2);}
/*******************************************************************************
           INQUIRY FUNCTIONS
*******************************************************************************/
/****************************************************************************
			Get Cell Index
-----------------------------------------------------------------------------
 output: k (OFF if outside, local cell index otherwise)
				 bool (true if inside)
****************************************************************************/
bool CRectGridp3D::GetCellIndex(const pt3D &pt, int &k) const{ 
//TMP DEBUG: presently only works with base of zero
	cmplex Z;
	int h,k2D;

  p2DGrid->GetCellIndex(pt,k2D);
  if (k2D==OFF){return false;} //outside horizontal extents of grid

  Z=p2DGrid->GlobalToLocal(c3Dto2D(pt));

	h=(int)(nH*pt.z/top[k2D]);
	while ((h<nH-1) && (pt.z/top[k2D]>gridH[h+1])){h++;}
	while ((h>0   ) && (pt.z/top[k2D]<gridH[h  ])){h--;}
	
	k=k2Dh_to_k(k2D,h);

	return true;
}
//---------------------------------------------------------------------------
bool CRectGridp3D::IsInside(const pt3D &pt) const{
//TMP DEBUG: presently only works with base of zero
	int k2D;
	if (p2DGrid->IsInside(pt)){
		p2DGrid->GetCellIndex(pt,k2D);
		if ((pt.z<0) || (pt.z>top[k2D])){return false;}
		return true;
	}
	else{
		return false;
	}
}
//---------------------------------------------------------------------------
bool CRectGridp3D::IsInCell(const pt3D &pt, const int k) const{
	//TMP DEBUG: presently only works with base of zero
	int k2D,h;
	//if ((k<0) || (k>nCells)){ExitGracefully("CRectGridp3D::IsInCell: Bad Cell index",RUNTIME_ERR);}
	
	k_to_k2Dh(k,k2D,h);

	if (!p2DGrid->IsInCell(pt,k2D)){return false;}
  if (pt.z<=gridH[h  ]*top[k2D]){return false;}
  if (pt.z> gridH[h+1]*top[k2D]){return false;}
  return true;
}

/***************************************************************************
                  GET INTERPOLATION WEIGHTS 
****************************************************************************
input:   Z (local coordinate)
output:  weights in X and Y direction (wx and wy) for bilinear interpolation between nodes
				 index of column to left of point (i) and row beneath point (j) 
--------------------------------------------------------------------------*/
void CRectGridp3D::GetInterpolationWeights(const pt3D &pt, 
																					 double &wx, double &wy, double &wh,
																					 int &i, int &j, int &h) const{
	int k2D;
	double Hloc;
	cmplex Z=p2DGrid->GlobalToLocal(c3Dto2D(pt));
	p2DGrid->GetInterpolationWeights(Z,wx,wy,i,j);

	k2D=p2DGrid->ij_to_k(i,j);

	Hloc=pt.z/top[k2D];
	//---------------z-direction interpolation-------------------------------------------
	if ((Hloc>nodeH[0]) && (Hloc<nodeH[nH-1])){
		h=((int)(nH*(Hloc-nodeH[0])/(nodeH[nH-1]-nodeH[0])));			 //guess @ y-index		
    
		if (h<0)    {h=0;}
		if (h>=nH-1){h=nH-2;} //correction for problem with inequality above for certain nums ((Hloc==nodeH[nH-1]))
    
		while ((h<nH-2) && (Hloc>nodeH[h+1])){h++;}                //guess higher
		while ((h>0   ) && (Hloc<nodeH[h  ])){h--;}								 //guess lower

		wh=(1.0-((Hloc-nodeH[h])/(nodeH[h+1]-nodeH[h])));					 //done: assign weight
	}
	else {//boundary node, must interpolate differently
		if      (Hloc<=nodeH[0   ]){h=0   ; wh=1.0;}
		else if (Hloc>=nodeH[nH-1]){h=nH-2; wh=0.0;} 
		else                           {ExitGracefully("CRectGridp3D::GetInterpolationWeights: bad H",RUNTIME_ERR);}
	}
	if (nH==1){h=0;wh=1.0;}

	if    ((k2D<0) || ((k2D>1)) || (h<0) || ((nH>1) && (h>=nH-1))){
		cout <<"k2D: "<<k2D<<" h:"<<h<<endl;
		ExitGracefully("CRectGridp3D::GetInterpolationWeights: bad vertical index",RUNTIME_ERR);}

	if ((wh>1)||(wh<0)){
		cout <<"wh: "<<wh<<endl;
		ExitGracefully("CRectGridp3D::GetInterpolationWeights: bad weights",RUNTIME_ERR);}
}
/***************************************************************************
                  INTERPOLATE VALUE(S) 
****************************************************************************
input:   z (global coordinates) 
				 nCell-sized array(s) of values to be interpolated, indexed by k (i.e. concentration)
				 InterpolateValues has additional parameter of the number of arrays (nval)
output:  interpolated value(s) ((*intval) for InterpolateValues function)
---------------------------------------------------------------------------*/
double CRectGridp3D::InterpolateValue(const pt3D &pt,Ironclad1DArray values) const{
	//weighted average of closest 4 values
	static double v1,v2,v3,v4;
	static double weightX,weightY,weightH,vfix;
  //TMP DEBUG-piecewise interpolation
	int k;
	GetCellIndex(pt,k); 
  if (k!=OFF){return values[k];}
	else       {return 0.0;}
	
	/*if (IsInside(pt)){

		int    i,j,h;

		
		GetInterpolationWeights(pt,weightX,weightY,weightH,i,j,h);  //i and j always good?

		vfix=values[ij_to_k(i,j)];

		if (p2DGrid->ij_to_k(i  ,j  )==OFF) {v1=0.0;}else{v1=values[ijh_to_k(i  ,j  ,h  )];}
		if (p2DGrid->ij_to_k(i+1,j  )==OFF) {v2=0.0;}else{v2=values[ijh_to_k(i+1,j  ,h  )];}
		if (p2DGrid->ij_to_k(i  ,j+1)==OFF) {v3=0.0;}else{v3=values[ijh_to_k(i  ,j+1,h  )];}
		if (p2DGrid->ij_to_k(i+1,j+1)==OFF) {v4=0.0;}else{v4=values[ijh_to_k(i+1,j+1,h  )];}

		//weighted average of nodal values [i,j], [i,j+1],[i+1,j], [i+1.j+1] 
		//using bilinear interpolation
		return (weightY  )*((weightX)  *(v1)+(1-weightX)*(v2))+
					 (1-weightY)*((weightX)  *(v3)+(1-weightX)*(v4));
	}
	else { //not on grid, cannot interpolate
		return 0.0;
	}*/
}
//---------------------------------------------------------------------------
void CRectGridp3D::InterpolateValues(const pt3D &pt,
																	Ironclad2DArray values, 
																	Writeable1DArray intval, 
																	const int nval) const{
	//average of closest 4 values
	int n;
	static double v1,v2,v3,v4;
	static double weightX,weightY,vfix;
	//TMP DEBUG 
	int k;
	GetCellIndex(pt,k); 
  if (k!=OFF){for (n=0;n<nval;n++){intval[n]=values[k][n];}}
	else       {for (n=0;n<nval;n++){intval[n]=0.0;}}

	/*if ((values==NULL) || (intval==NULL)){ExitGracefully("CRectGridp3D::InterpolateValues: cannot interpolate on NULL value array",RUNTIME_ERR);}
	if (IsInside(pt)){

		int    i,j;
		cmplex Z=GlobalToLocal(c3Dto2D(pt));	
		
		GetInterpolationWeights(Z,weightX,weightY,i,j);

		
    for (n=0;n<nval;n++){
			vfix=values[ij_to_k(i,j)][n];
		//weighted average of nodal values [i,j,k], [i,j+1,k],[i+1,j,k], [i+1.j+1,k] 
			if (ij_to_k(i  ,j  )==OFF) {v1=vfix;}else{v1=values[ij_to_k(i  ,j  )][n];}
			if (ij_to_k(i+1,j  )==OFF) {v2=vfix;}else{v2=values[ij_to_k(i+1,j  )][n];}
			if (ij_to_k(i  ,j+1)==OFF) {v3=vfix;}else{v3=values[ij_to_k(i  ,j+1)][n];}
			if (ij_to_k(i+1,j+1)==OFF) {v4=vfix;}else{v4=values[ij_to_k(i+1,j+1)][n];}

			intval[n]= (weightY  )*((weightX)  *(v1)+(1-weightX)*(v2))+
								 (1-weightY)*((weightX)  *(v3)+(1-weightX)*(v4));
		}
	}
	else { //not on grid, cannot interpolate
    for (n=0;n<nval;n++){intval[n]=0.0;}
	}*/
}
//---------------------------------------------------------------------------
double CRectGridp3D::InterpolateDerivative(const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const{
  ExitGracefully("CRectGridp3D::InterpolateDerivative: STUB Function",BAD_DATA);
	return 0.0;
}
/***************************************************************************
                  GET NEIGHBORS
****************************************************************************
input:   cell index (k)
				 empty integer array to store neighbor values (pre-allocated)
				 empty double array to store distance to neighbors (pre-allocated)
				 maximum number of neighbors (numneighbors)
output:  neighbor indices (neighbors), numneighbors (maximum 4 for grid)
--------------------------------------------------------------------------*/
void CRectGridp3D::GetNeighbors  (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const{
  ExitGracefully("CRectGridp3D::GetNeighbors: STUB Function",BAD_DATA);
	return;

	/*int i, j, count(0);
	k_to_ij(k,i,j);
	if (numneighbors<4){ExitGracefully("CRectGridp3D::GetNeighbors: not enough space allocated", RUNTIME_ERR);}
	if (i>0   ){neighbors[count]=ij_to_k(i-1,j  );dist[count]=fabs(nodeX[i-1]-nodeX[i]);count++;}
	if (i<nX-1){neighbors[count]=ij_to_k(i+1,j  );dist[count]=fabs(nodeX[i+1]-nodeX[i]);count++;}
	if (j>0   ){neighbors[count]=ij_to_k(i  ,j-1);dist[count]=fabs(nodeY[j-1]-nodeY[j]);count++;}
	if (j>nY-1){neighbors[count]=ij_to_k(i  ,j+1);dist[count]=fabs(nodeY[j+1]-nodeY[j]);count++;}
	numneighbors=count;*/
}
/****************************************************************************
                  DISTRIBUTE POINTS 
*****************************************************************************
input:   cell index (k)
				 distribution type (ty) :random or fixed pattern
				 empty cmplex array (pts) to store (npt) point locations
output:  cmplex array of locations (pts)
---------------------------------------------------------------------------*/
void CRectGridp3D::DistributePoints(const int k, distribution_type ty,const int npts, pt3D *pts) const{
	int k2D,h;
	k_to_k2Dh(k,k2D,h);
	
	p2DGrid->DistributePoints(k2D,ty,npts,pts);

	//if ((i==OFF) || (j==OFF)){ExitGracefully("CRectGridp3D::DistributePoints: Bad i or j", RUNTIME_ERR);}
	for (int n=0; n<npts; n++){
		if (ty==RANDOM_PATTERN){
			pts[n].z=quickrandom(gridH[h]*top[k2D],gridH[h+1]*top[k2D]);
		}
		else{
			ExitGracefully("CRectGridp3D::DistributePoints: cannot distrubute in non-random pattern yet",BAD_DATA);
		} //end if random
	}//end for i...
}
/****************************************************************************
                  Get Mesh Boundary Intercepts
***************************************************************************
input:   endpoints pt1 and pt2 of a line
         max array size of ptint,i,and j: ncellscrossed+1
output:  intercepts of grid faces (zint[ncellscrossed+1]), sorted from z1 to z2 (including z1 and z2 if inside grid)
         k indices of all cells crossed k[ncellscrossed]
         nint, the total number of intersections
returns: true if both points are in grid
---------------------------------------------------------------------------*/
bool CRectGridp3D::GetMeshBoundaryIntercepts(const pt3D &pt1, const pt3D &pt2, pt3D *ptint, int *k, int &ncellscrossed) const{

  ExitGracefully("CRectGridp3D::GetMeshBoundaryIntercepts: STUB Function",BAD_DATA);
	return false;
}
/***************************************************************************
                  WRITE GEOMETRY
****************************************************************************
 Writes grid to Atlas BNA file, rectgrid.bna
--------------------------------------------------------------------------*/
void CRectGridp3D::WriteGeometry() const{
	p2DGrid->WriteGeometry();
}
/***************************************************************************
                  WRITE RECTGRID TO BBG FILE
****************************************************************************
 Writes grid to grid.bbg (BlueBird Grid) file 
--------------------------------------------------------------------------*/
void CRectGridp3D::WriteBBGFile() const{
	//p2DGrid->WriteBBGFile();
	int i;
  ofstream BBG;
	BBG.open("grid.bbg");
  BBG<<"FiniteDifferenceGrid3D BBGenerated"<<endl;
  BBG<<p2DGrid->zg.real()<<" "<<p2DGrid->zg.imag()<<endl;
	BBG<<p2DGrid->orient<<endl;
	BBG<<p2DGrid->nX<<" "<<p2DGrid->nY<<" "<<nH<<endl;
	for (i=0; i<=p2DGrid->nX; i++){BBG<<p2DGrid->gridX[i]<<endl;}
	BBG<<"&"<<endl;
	for (i=0; i<=p2DGrid->nY; i++){BBG<<p2DGrid->gridY[i]<<endl;}
	BBG<<"&"<<endl;
	for (i=0; i<=nH; i++){BBG<<gridH[i]<<endl;}
	BBG<<"&"<<endl;
	BBG.close();
};

/***************************************************************************
                  PARSE RECTGRID
****************************************************************************
 Reads grid from *.bbg (BlueBird Grid) file 
--------------------------------------------------------------------------*/
CRectGridp3D *CRectGridp3D::ReadBBGFile(ifstream &input,int &l){
	//string "FiniteDifferenceGrid3D", string label (unused) 
	//double x_bottom_left y_bottom_left
	//double pivot_angle
	//int nX (number of columns) int nY (number of rows) int nH
	//{double X }x(nX+1)
  //& 
	//{double y }x(nY+1)
  //& 
	//{double H }x(nH+1)
	//&

	CRectGrid    *pRectGrid;
  CRectGridp3D *p3DGrid;

  bool     eof(false),done(false);
	double   xp,yp,angle;
	int      Len,nX,nY,nH,i;
  double  *X=NULL;
	double  *Y=NULL; 
	double  *H=NULL;
	char    *s[MAXINPUTITEMS];
	pRectGrid=NULL;

	if (parserdebug) {cout << "Rectangular Grid (3D)"<<endl;}

	//starts at line #2
	eof=TokenizeLine(input,s,Len); l++; 
	eof=TokenizeLine(input,s,Len); l++; 
	do{
    if      (Len==2) {xp=s_to_d(s[0]); 
											yp=s_to_d(s[1]); 
											done=true;												 eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line "<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
  do{ 
		if      (Len==1) {angle=s_to_d(s[0]);     done=true; eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line "<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	do{ 
    if      (Len==3) {nX=s_to_i(s[0]); 
											nY=s_to_i(s[1]); 
											nH=s_to_i(s[2]);
											done=true;												 
											X=new double [nX+1];
											Y=new double [nY+1];
											H=new double [nH+1];
											ExitGracefullyIf(H==NULL,"CRectGridp3D::Parse",OUT_OF_MEMORY);
																												 eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==2) {nX=s_to_i(s[0]); 
											nY=s_to_i(s[1]); 
											nH=1;
											done=true;												 
											X=new double [nX+1];
											Y=new double [nY+1];
											H=new double [2]; H[0]=0.0; H[1]=1.0;
											ExitGracefullyIf(H==NULL,"CRectGridp3D::Parse",OUT_OF_MEMORY);
																												 eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line "<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	i=0;
  do {
		if ((Len==1) && (strcmp(s[0],"&"))){
			if (i>nX)  {ExitGracefully("CRectGridp3D::Parse- too many columns in grid",TOO_MANY);}
			X[i]=s_to_d(s[0]); i++;                            eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			if (i==nX+1){done=true;}
			else        {ExitGracefully("CRectGridp3D::Parse- not enough columns in grid",BAD_DATA);}  
		}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	i=0;                                                   eof=TokenizeLine(input,s,Len); l++;
  do {
		if ((Len==1) && (strcmp(s[0],"&"))){
			if (i>nY)  {ExitGracefully("CRectGridp3D::Parse- too many columns in grid",TOO_MANY);}
			Y[i]=s_to_d(s[0]); i++;                            eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			if (i==nY+1){done=true;}
			else        {ExitGracefully("CRectGridp3D::Parse- not enough rows in grid",BAD_DATA);}  
		}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	i=0;  
	do {
		if ((Len==1) && (strcmp(s[0],"&"))){
			if (i>nH)  {ExitGracefully("CRectGridp3D::Parse- too many layers in grid",TOO_MANY);}
			H[i]=s_to_d(s[0]); i++;                            eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			if (i==nH+1){done=true;}
			else        {ExitGracefully("CRectGridp3D::Parse- not enough layers in grid",BAD_DATA);}  
		}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));

	pRectGrid = new CRectGrid(cmplex(xp,yp),angle,X,Y,nX,nY); 
	
	if (pRectGrid==NULL){return NULL;}

  p3DGrid   = new CRectGridp3D(pRectGrid,H,nH);

	delete [] X;
	delete [] Y;
	delete [] H;

	return p3DGrid;
}
