//RectGrid.cpp
#include "RectGrid.h"

/*******************************************************************************
           CONSTRUCTORS
*******************************************************************************/
CRectGrid::CRectGrid(const double north, 
										 const double south, 
										 const double east, 
										 const double west,
										 const int    resX, 
										 const int    resY){
	//cartesian grid oriented normally
	int i,j;
	Xmax=east-west;
	Ymax=north-south;
	nX=resX;
	nY=resY;
	nCells=nNodes=nX*nY;
	zg=cmplex(west,south);
	orient=0.0;
	gridX =new double [nX+1];
	gridY =new double [nY+1];
	nodeX =new double [nX];
	nodeY =new double [nY];
	if (nodeY==NULL){ExitGracefully("CRectGrid::CRectGrid: out of memory", OUT_OF_MEMORY);}
	for (i=0; i<=nX; i++){gridX[i]=Xmax*(double)(i)/(double)(nX);}
	for (i=0; i< nX; i++){nodeX[i]=(gridX[i]+gridX[i+1])/2.0;}
	for (j=0; j<=nY; j++){gridY[j]=Ymax*(double)(j)/(double)(nY);}
	for (j=0; j< nY; j++){nodeY[j]=(gridY[j]+gridY[j+1])/2.0;}
	BBox.n=north;
	BBox.s=south;
	BBox.e=east;
	BBox.w=west;
	bounds.nsides=4;
	bounds.vertices = new cmplex [5];
	bounds.vertices[0]=cmplex(east,north);//counterclockwise-mild inclusion into middle
	bounds.vertices[1]=cmplex(west,north);
	bounds.vertices[2]=cmplex(west,south);
	bounds.vertices[3]=cmplex(east,south);
	bounds.vertices[4]=cmplex(east,north);
	is3D=false;
	meshtype=CELL_BASED;
}
//---------------------------------------------------------------------------
CRectGrid::CRectGrid(const cmplex  zpivot, 
										 const double  angle,
										 const double *X, 
										 const double *Y,
										 const int     numcellsx,
										 const int     numcellsy){
	int i,j;

	nX=numcellsx;
	nY=numcellsy;
	nCells=nNodes=nX*nY;
	zg=zpivot;
	orient=angle;
	gridX	=new double [nX+1];
	gridY	=new double [nY+1];
	nodeX	=new double [nX];
	nodeY	=new double [nY];
	if (nodeY==NULL){ExitGracefully("CRectGrid::CRectGrid: out of memory", OUT_OF_MEMORY);}

	for (i=0; i<=nX; i++){gridX[i]=X[i];}
	for (i=0; i< nX; i++){nodeX[i]=(gridX[i]+gridX[i+1])/2.0;}
	for (j=0; j<=nY; j++){gridY[j]=Y[j];}
	for (j=0; j< nY; j++){nodeY[j]=(gridY[j]+gridY[j+1])/2.0;}

	Xmax=gridX[nX];
	Ymax=gridY[nY];

	double ymax(-ALMOST_INF),ymin(ALMOST_INF);             //define bounding box
	double xmax(-ALMOST_INF),xmin(ALMOST_INF);

	cmplex z1(LocalToGlobal(cmplex(0.0 ,0.0 )));
	cmplex z2(LocalToGlobal(cmplex(Xmax,0.0 )));
	cmplex z3(LocalToGlobal(cmplex(Xmax,Ymax)));
	cmplex z4(LocalToGlobal(cmplex(0.0 ,Ymax)));
	upperswap(xmax,z1.real());	lowerswap(xmin,z1.real());
	upperswap(xmax,z2.real());	lowerswap(xmin,z2.real());
	upperswap(xmax,z3.real());	lowerswap(xmin,z3.real());
	upperswap(xmax,z4.real());	lowerswap(xmin,z4.real());

	upperswap(ymax,z1.imag());	lowerswap(ymin,z1.imag());
	upperswap(ymax,z2.imag());	lowerswap(ymin,z2.imag());
	upperswap(ymax,z3.imag());	lowerswap(ymin,z3.imag());
	upperswap(ymax,z4.imag());	lowerswap(ymin,z4.imag());
	BBox.n=ymax;								BBox.s=ymin; 
	BBox.e=xmax;								BBox.w=xmin;
	bounds.nsides=4;
	bounds.vertices = new cmplex [5];
	bounds.vertices[0]=z4;//counterclockwise
	bounds.vertices[1]=z3;
	bounds.vertices[2]=z2;
	bounds.vertices[3]=z1;
	bounds.vertices[4]=z4;
	is3D=false;
	//cout <<"CRECTGRID:"<< nX <<" "<<nY<<endl;
	meshtype=CELL_BASED;
}
//---------------------------------------------------------------------------
CRectGrid::~CRectGrid(){
	if (globaldebug){cout<<"  DESTROYING RECTGRID"<<endl;}
	delete [] gridX;
	delete [] gridY;
	delete [] nodeX;
	delete [] nodeY;
}
/*******************************************************************************
           PRIVATE CONVERSION FUNCTIONS 
//SHOULD BE OPTIMIZED- could be inlined
*******************************************************************************/
int  CRectGrid::ij_to_k(const int i, const int j) const{
	if ((i==OFF) || (j==OFF) ||(j>=nY) || (j<0) || (i>=nX) || (i<0)) {return OFF;}
	//if ((j*nX+i>=nCells) || (j*nX+i<0)){ExitGracefully("CRectGrid::ij_to_k: bad input", RUNTIME_ERR);}

	return j*(nX)+i;
}
//---------------------------------------------------------------------------
void CRectGrid::k_to_ij(const int k, int &i, int &j) const{
	if ((k==OFF) || (k>=nCells)){i=OFF;j=OFF;return;}
	i=k%(nX);
	j=(k-i)/(nX);
}
//---------------------------------------------------------------------------
cmplex CRectGrid::GlobalToLocal(const cmplex &z) const{
	return (z-zg)*cmplex(cos(orient),-sin(orient));
}
//---------------------------------------------------------------------------
cmplex CRectGrid::LocalToGlobal(const cmplex &Z) const{
	//return Z/(cmplex(cos(orient),-sin(orient)))+zg;
	return Z*cmplex(cos(-orient),-sin(-orient))+zg;
}

/*******************************************************************************
           ACCESSOR FUNCTIONS 
*******************************************************************************/
pt3D CRectGrid::GetNodeLocation(const int k) const{
	if ((k<0) || (k>nCells)){ExitGracefully("CRectGrid::GetNodeLocation: Bad Cell index",RUNTIME_ERR);}
	int i,j;
	k_to_ij(k,i,j);
	cmplex z=LocalToGlobal(cmplex(nodeX[i],nodeY[j]));
	return pt3D(z.real(),z.imag(),0.0); //may wish to ignore height for now??
}
//---------------------------------------------------------------------------
double CRectGrid::GetCellArea(const int k) const{
	if ((k<0) || (k>nCells)){ExitGracefully("CRectGrid::GetCellArea: Bad Cell index",RUNTIME_ERR);}
	int i,j;
	k_to_ij(k,i,j);
	return (gridX[i+1]-gridX[i])*(gridY[j+1]-gridY[j]);
}
//---------------------------------------------------------------------------
int    CRectGrid::GetNumCellSides          (const int k) const {
	if ((k<0) || (k>nCells)){ExitGracefully("CRectGrid::GetNumCellSides: Bad Cell index",RUNTIME_ERR);}
	return 4;
}
//---------------------------------------------------------------------------
void   CRectGrid::GetCellSideVertices      (const int k, const int s, cmplex &z1, cmplex &z2) const{
	if ((k<0) || (k>nCells)){ExitGracefully("CRectGrid::GetCellSideVertices: Bad Cell index",RUNTIME_ERR);}    

	int i,j;
	k_to_ij(k,i,j);
	if      (s==0){z1=LocalToGlobal(cmplex(gridX[i  ],gridY[j+1])); z2=LocalToGlobal(cmplex(gridX[i+1],gridY[j+1]));}
	else if (s==1){z1=LocalToGlobal(cmplex(gridX[i+1],gridY[j+1])); z2=LocalToGlobal(cmplex(gridX[i+1],gridY[j  ]));}
	else if (s==2){z1=LocalToGlobal(cmplex(gridX[i+1],gridY[j  ])); z2=LocalToGlobal(cmplex(gridX[i  ],gridY[j  ]));}
	else if (s==3){z1=LocalToGlobal(cmplex(gridX[i  ],gridY[j  ])); z2=LocalToGlobal(cmplex(gridX[i  ],gridY[j+1]));}
	else          {ExitGracefully("CRectGrid::GetCellSideVertices: Bad Cell side index",RUNTIME_ERR);}
}
//---------------------------------------------------------------------------
window CRectGrid::GetBoundingBox()      const{return BBox;}

/*******************************************************************************
           INQUIRY FUNCTIONS
*******************************************************************************/
/****************************************************************************
			Get Cell Index
-----------------------------------------------------------------------------
 output: k (OFF if outside, local cell index otherwise)
				 bool (true if inside)
****************************************************************************/
bool CRectGrid::GetCellIndex(const pt3D &pt, int &k) const{ 

	cmplex Z;
	int i,j;

	if (!IsInside(pt)){k=OFF; return false;}

  Z=GlobalToLocal(c3Dto2D(pt));

	i=(int)(nX*Z.real()/Xmax);
	while ((i<nX-1) && (Z.real()>gridX[i+1])){i++;}
	while ((i>0   ) && (Z.real()<gridX[i  ])){i--;}
	
	j=(int)(nY*Z.imag()/Ymax);
	while ((j<nY-1) && (Z.imag()>gridY[j+1])){j++;}
	while ((j>0   ) && (Z.imag()<gridY[j  ])){j--;}

	k=ij_to_k(i,j);

	//cout << "x: "<<z.real() << " y: "<< z.imag() << " i: "<<i<< " j: " <<j<<endl;
	return true;
}
//---------------------------------------------------------------------------
bool CRectGrid::IsInside(const pt3D &pt) const{
	cmplex z=c3Dto2D(pt);
	cmplex Z=GlobalToLocal(z);
	if (Z.real()<0.0 ){return false;}
	if (Z.real()>Xmax){return false;}
	if (Z.imag()<0.0 ){return false;}
  if (Z.imag()>Ymax){return false;}
	return true;
}
//---------------------------------------------------------------------------
bool CRectGrid::IsInCell(const pt3D &pt, const int k) const{
	int i,j;
	//if ((k<0) || (k>nCells)){ExitGracefully("CRectGrid::IsInCell: Bad Cell index",RUNTIME_ERR);}
	
	k_to_ij(k,i,j);

	cmplex Z=GlobalToLocal(c3Dto2D(pt));

	if(Z.real()<gridX[i  ]) {return false;}
	if(Z.real()>gridX[i+1]) {return false;}
	if(Z.imag()<gridY[j  ]) {return false;}
	if(Z.imag()>gridY[j+1]) {return false;}
	return true;
}

/***************************************************************************
                  GET INTERPOLATION WEIGHTS 
****************************************************************************
input:   Z (local coordinate)
output:  weights in X and Y direction (wx and wy) for bilinear interpolation between nodes
				 index of column to left of point (i) and row beneath point (j) 
--------------------------------------------------------------------------*/
void CRectGrid::GetInterpolationWeights(const cmplex &Z, double &wx, double &wy, int &i, int &j) const{

	//---------------X-direction interpolation-------------------------------------------
	if ((Z.real()>nodeX[0]) && (Z.real()<nodeX[nX-1])) {
		i=((int)(nX*(Z.real()-nodeX[0])/(nodeX[nX-1]-nodeX[0])));      //guess @ x-index
		if (i<0)    {i=0;}
    if (i>=nX-1){i=nX-2;} 
		while ((i<nX-2) && (Z.real()>nodeX[i+1])){i++;}                //guess higher
		while ((i>0   ) && (Z.real()<nodeX[i  ])){i--;}                //guess lower
		wx=(1.0-(Z.real()-nodeX[i])/(nodeX[i+1]-nodeX[i]));						 //done: assign weight

	}
	else{//boundary node, must interpolate differently
		if      (Z.real()<=nodeX[0   ]){i=0   ; wx=1.0;}
		else if (Z.real()>=nodeX[nX-1]){i=nX-2; wx=0.0;} 
		else                           {ExitGracefully("CRectGrid::GetInterpolationWeights: bad X",RUNTIME_ERR);}
	}
	if (nX==1){i=0;wx=1.0;}

	//---------------Y-direction interpolation-------------------------------------------
	if ((Z.imag()>nodeY[0]) && (Z.imag()<nodeY[nY-1])){
		j=((int)(nY*(Z.imag()-nodeY[0])/(nodeY[nY-1]-nodeY[0])));			 //guess @ y-index		
    
		if (j<0)    {j=0;}
		if (j>=nY-1){j=nY-2;} //correction for problem with inequality above for certain nums ((Z.imag()==nodeY[nY-1]))
    
		while ((j<nY-2) && (Z.imag()>nodeY[j+1])){j++;}                //guess higher
		while ((j>0   ) && (Z.imag()<nodeY[j  ])){j--;}								 //guess lower
		wy=(1.0-((Z.imag()-nodeY[j])/(nodeY[j+1]-nodeY[j])));					 //done: assign weight
	}
	else {//boundary node, must interpolate differently
		if      (Z.imag()<=nodeY[0   ]){j=0   ; wy=1.0;}
		else if (Z.imag()>=nodeY[nY-1]){j=nY-2; wy=0.0;} 
		else                           {ExitGracefully("CRectGrid::GetInterpolationWeights: bad Y",RUNTIME_ERR);}
	}
	if (nY==1){j=0;wy=1.0;}

	if    ((i<0) || ((nX>1) && (i>=nX-1)) || (j<0) || ((nY>1) && (j>=nY-1))){
		cout <<"i: "<<i<<" j:"<<j<<endl;
		ExitGracefully("CRectGrid::GetInterpolationWeights: bad index",RUNTIME_ERR);}

	if ((wy>1)||(wy<0)||(wx>1)||(wx<0)){
		cout <<"wy: "<<wy<<" wx:"<<wx<<endl;
		ExitGracefully("CRectGrid::GetInterpolationWeights: bad weights",RUNTIME_ERR);}
}
/***************************************************************************
                  INTERPOLATE VALUE(S) 
****************************************************************************
input:   z (global coordinates) 
				 nCell-sized array(s) of values to be interpolated, indexed by k (i.e. concentration)
				 InterpolateValues has additional parameter of the number of arrays (nval)
output:  interpolated value(s) ((*intval) for InterpolateValues function)
---------------------------------------------------------------------------*/
double CRectGrid::InterpolateValue(const pt3D &pt,Ironclad1DArray values) const{
	//weighted average of closest 4 values
	static double v1,v2,v3,v4;
	static double weightX,weightY,vfix;

	if (IsInside(pt)){
		cmplex Z=GlobalToLocal(c3Dto2D(pt));

		int    i,j;

		GetInterpolationWeights(Z,weightX,weightY,i,j);  //i and j always good?

		vfix=values[ij_to_k(i,j)];

		if (ij_to_k(i  ,j  )==OFF) {v1=0.0;}else{v1=values[ij_to_k(i  ,j  )];}
		if (ij_to_k(i+1,j  )==OFF) {v2=0.0;}else{v2=values[ij_to_k(i+1,j  )];}
		if (ij_to_k(i  ,j+1)==OFF) {v3=0.0;}else{v3=values[ij_to_k(i  ,j+1)];}
		if (ij_to_k(i+1,j+1)==OFF) {v4=0.0;}else{v4=values[ij_to_k(i+1,j+1)];}

		//weighted average of nodal values [i,j], [i,j+1],[i+1,j], [i+1.j+1] 
		//using bilinear interpolation
		return (weightY  )*((weightX)  *(v1)+(1-weightX)*(v2))+
					 (1-weightY)*((weightX)  *(v3)+(1-weightX)*(v4));
	}
	else { //not on grid, cannot interpolate
		return 0.0;
	}
}
//---------------------------------------------------------------------------
void CRectGrid::InterpolateValues(const pt3D &pt,
																	Ironclad2DArray values, 
																	Writeable1DArray intval, 
																	const int nval) const{
	//average of closest 4 values
	int n;
	static double v1,v2,v3,v4;
	static double weightX,weightY,vfix;

	if ((values==NULL) || (intval==NULL)){ExitGracefully("CRectGrid::InterpolateValues: cannot interpolate on NULL value array",RUNTIME_ERR);}
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
	}
}
//---------------------------------------------------------------------------
double CRectGrid::InterpolateDerivative(const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const{
	//d/dX is average change in X direction
	//d/dX is average change in Y direction

	if (dir==0.0){return 0.0;}

	cmplex locdir=dir*cmplex(cos(orient),-sin(orient))/dir;//unit direction in local coords

	if (IsInside(pt)){
		cmplex Z=GlobalToLocal(c3Dto2D(pt));
		
		double weightX,weightY,dx,dy;
		int    i,j;

		GetInterpolationWeights(Z,weightX,weightY,i,j);

		dx=   weightY *(values[ij_to_k(i+1,j  )]-values[ij_to_k(i,j  )])+
			 (1-weightY)*(values[ij_to_k(i+1,j+1)]-values[ij_to_k(i,j+1)]);
		dx/=(gridX[i+1]-gridX[i]);

		dy=   weightX *(values[ij_to_k(i  ,j+1)]-values[ij_to_k(i  ,j  )])+
			 (1-weightX)*(values[ij_to_k(i+1,j+1)]-values[ij_to_k(i+1,j  )]);
		dy/=(gridY[j+1]-gridY[j]);

		return sqrt(locdir.real()*dx*locdir.real()*dx+locdir.imag()*dy*locdir.imag()*dy);
	}
	else { //not on grid, cannot interpolate
		return 0.0;
	}	
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
void CRectGrid::GetNeighbors  (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const{
	int i, j, count(0);
	k_to_ij(k,i,j);
	if (numneighbors<4){ExitGracefully("CRectGrid::GetNeighbors: not enough space allocated", RUNTIME_ERR);}
	if (i>0   ){neighbors[count]=ij_to_k(i-1,j  );dist[count]=fabs(nodeX[i-1]-nodeX[i]);count++;}
	if (i<nX-1){neighbors[count]=ij_to_k(i+1,j  );dist[count]=fabs(nodeX[i+1]-nodeX[i]);count++;}
	if (j>0   ){neighbors[count]=ij_to_k(i  ,j-1);dist[count]=fabs(nodeY[j-1]-nodeY[j]);count++;}
	if (j>nY-1){neighbors[count]=ij_to_k(i  ,j+1);dist[count]=fabs(nodeY[j+1]-nodeY[j]);count++;}
	numneighbors=count;
}
/****************************************************************************
                  DISTRIBUTE POINTS 
*****************************************************************************
input:   cell index (k)
				 distribution type (ty) :random or fixed pattern
				 empty cmplex array (pts) to store (npt) point locations
output:  cmplex array of locations (pts)
---------------------------------------------------------------------------*/
void CRectGrid::DistributePoints(const int k, distribution_type ty,const int npts, pt3D *pts) const{
	int i,j;
	k_to_ij(k,i,j);
	//if ((i==OFF) || (j==OFF)){ExitGracefully("CRectGrid::DistributePoints: Bad i or j", RUNTIME_ERR);}
	for (int n=0; n<npts; n++){
		if (ty==RANDOM_PATTERN){
			pts[n]=c2Dto3D(LocalToGlobal(cmplex(quickrandom(gridX[i],gridX[i+1]),
																					quickrandom(gridY[j],gridY[j+1]))));
		}
		else{
			//placed like dots on dice...
			if      (npts==1){
				pts[0]=c2Dto3D(LocalToGlobal(cmplex(                nodeX[i],                 nodeY[j]  )));
			}
			else if (npts==2){
				pts[0]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[1]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));			
			}
			else if (npts==3){
				pts[0]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[1]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));	
				pts[2]=c2Dto3D(LocalToGlobal(cmplex(                nodeX[i],                 nodeY[j]  )));			
			}
			else if (npts==4){
				pts[0]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[1]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));	
				pts[2]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[3]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));			
			}
			else if (npts==5){
				pts[0]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[1]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));	
				pts[2]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[3]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));	
				pts[4]=c2Dto3D(LocalToGlobal(cmplex(                nodeX[i],                 nodeY[j] )));				
			}
			else if (npts==6){
				pts[0]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[1]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));	
				pts[2]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),0.5*(gridY[j  ]+nodeY[j]) )));
        pts[3]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),0.5*(gridY[j+1]+nodeY[j]) )));	
				pts[4]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i  ]+nodeX[i]),                nodeY[j]  )));			
				pts[5]=c2Dto3D(LocalToGlobal(cmplex(0.5*(gridX[i+1]+nodeX[i]),                nodeY[j]  )));
			}
			else {
				//TMP DEBUG- not yet coded
				ExitGracefully("CRectGrid::DistributePoints: cannot evenly distribute more than 6 points",RUNTIME_ERR);
			} //end if npts==1...
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
bool CRectGrid::GetMeshBoundaryIntercepts(const pt3D &pt1, const pt3D &pt2, pt3D *ptint, int *k, int &ncellscrossed) const{


	int     i1,j1,i2,j2,k1,k2;
	int     ijump,jjump,Xjumpdir(1),Yjumpdir(1),itmp,jtmp;
	int     n,add,max_crossings(ncellscrossed-1);
	bool    insidegrid1,insidegrid2;
	double  XJump,YJump;
	double *pctX,*pctY;       //range from zero to 1
	cmplex  Z1,Z2,z1,z2,*zint;
  int    *i,*j;

	//preprocessing of data--------------------------------------------------
	//-----------------------------------------------------------------------
	insidegrid1=IsInside(pt1);
	insidegrid2=IsInside(pt2);
	z1=c3Dto2D(pt1);
	z2=c3Dto2D(pt2);

	Z1=GlobalToLocal(z1);
  Z2=GlobalToLocal(z2);	
	
	cmplex junk;
	if ((!insidegrid1) && (!insidegrid2) && 
		 (Intersect(Z1,Z2,cmplex(gridX[0 ],gridY[0 ]),cmplex(gridX[0 ],gridY[nY]),junk)!=INTERSECTED) &&
		 (Intersect(Z1,Z2,cmplex(gridX[0 ],gridY[nY]),cmplex(gridX[nX],gridY[nY]),junk)!=INTERSECTED) &&
     (Intersect(Z1,Z2,cmplex(gridX[nX],gridY[nY]),cmplex(gridX[nX],gridY[0 ]),junk)!=INTERSECTED) &&
     (Intersect(Z1,Z2,cmplex(gridX[nX],gridY[0 ]),cmplex(gridX[0 ],gridY[0 ]),junk)!=INTERSECTED)
		 ){
    ncellscrossed=0;
		return false;
	}
	
	if (!insidegrid1){
		//interpolate z1 to be entrance location
		if      (Z1.real()<gridX[0 ]){
			Z1=LinInterp(0.0,1.0,Z1.real(),Z2.real(),gridX[0 ]+REALSMALL*(gridX[nX]-gridX[0]))*(Z2-Z1)+Z1;
		}
		else if (Z1.real()>gridX[nX]){
			Z1=LinInterp(0.0,1.0,Z1.real(),Z2.real(),gridX[nX]-REALSMALL*(gridX[nX]-gridX[0]))*(Z2-Z1)+Z1;
		}
		if      (Z1.imag()<gridY[0 ]){
			Z1=LinInterp(0.0,1.0,Z1.imag(),Z2.imag(),gridY[0 ]+REALSMALL*(gridY[nY]-gridY[0]))*(Z2-Z1)+Z1;
		}
		else if (Z1.imag()>gridY[nY]){
			Z1=LinInterp(0.0,1.0,Z1.imag(),Z2.imag(),gridY[nY]-REALSMALL*(gridY[nY]-gridY[0]))*(Z2-Z1)+Z1;
		}
		z1=LocalToGlobal(Z1);
	}
	if (!insidegrid2){
		//interpolate z2 to be slightly inside exit location
		
		if      (Z2.real()<gridX[0 ]){
			Z2=LinInterp(0.0,1.0,Z1.real(),Z2.real(),gridX[0 ]+REALSMALL*(gridX[nX]-gridX[0]))*(Z2-Z1)+Z1;
		}
		else if (Z2.real()>gridX[nX]){
			Z2=LinInterp(0.0,1.0,Z1.real(),Z2.real(),gridX[nX]-REALSMALL*(gridX[nX]-gridX[0]))*(Z2-Z1)+Z1;
		}
		if      (Z2.imag()<gridY[0 ]){
			Z2=LinInterp(0.0,1.0,Z1.imag(),Z2.imag(),gridY[0 ]+REALSMALL*(gridY[nY]-gridY[0]))*(Z2-Z1)+Z1;
		}
		else if (Z2.imag()>gridY[nY]){
			Z2=LinInterp(0.0,1.0,Z1.imag(),Z2.imag(),gridY[nY]-REALSMALL*(gridY[nY]-gridY[0]))*(Z2-Z1)+Z1;
		}
		z2=LocalToGlobal(Z2);	
	}
	//if ((!insidegrid2) && (!insidegrid1)){nint=0; return false;}//TMP DEBUG-should be able to handle this condition

	Z1=GlobalToLocal(z1); 
	Z2=GlobalToLocal(z2);

	GetCellIndex(c2Dto3D(z1),k1);//operate on revised location-should always be inside
	GetCellIndex(c2Dto3D(z2),k2);

	if ((k1==OFF) || (k2==OFF) || (k1>nCells) || (k2>nCells)){
		ExitGracefully("CRectGrid::GetMeshBoundaryIntercepts: Improper indices returned(1)",RUNTIME_ERR);}

	k_to_ij(k1,i1,j1);
	k_to_ij(k2,i2,j2);

	if ((i1==OFF) || (j1==OFF) || (i2==OFF) || (j2==OFF)){
		cout << i1 <<" "<<j1<<" "<<k1<<" "<< i2 <<" "<<j2<<" "<<k2<<" "<<endl;
		ExitGracefully("CRectGrid::GetMeshBoundaryIntercepts: Improper indices returned(2)",RUNTIME_ERR);}

  //initialize X vars------------------------------------------------------   
	XJump   =(Z2-Z1).real(); 	
	ijump   =(i2-i1);                 
	if (ijump!=0){Xjumpdir=ijump/iabs(ijump);}else{Xjumpdir=1;}          

	//initialize Y Vars------------------------------------------------------
	YJump   =(Z2-Z1).imag();	
	jjump   =(j2-j1);
	if (jjump!=0){Yjumpdir=jjump/iabs(jjump);}else{Yjumpdir=1;} 
	
	int nint=iabs(ijump)+iabs(jjump);

	if (nint>=max_crossings){
		ExitGracefully("CRectGrid::GetMeshBoundaryIntercepts: Exceeded allowable intersections",RUNTIME_ERR);}
	
	//Create storage arrays--------------------------------------------------
	if (nint>0){
		pctX=new double [nint];
		pctY=new double [nint];
	}
	else {
		pctX=NULL;
		pctY=NULL;
	}
	zint=new cmplex [nint+2]; 
  if (zint==NULL){ExitGracefully("CRectGrid::GetMeshBoundaryIntercepts: Out of memory",OUT_OF_MEMORY);}

	//process interface jumps in X-direction---------------------------------
	//-----------------------------------------------------------------------
	if (Xjumpdir<0){add=0;}else{add=1;}							                 
  for (n=0; n<iabs(ijump); n++){
		pctX[n]=LinInterp(0.0,1.0,Z1.real(),Z2.real(),gridX[i1+Xjumpdir*n+add]);//cout <<"pctX["<<n<<"]: "<<pctX[n]<<endl;
	}
	//process interface jumps in Y-direction---------------------------------
	//----------------------------------------------------------------------- 
	if (Yjumpdir<0){add=0;}else{add=1;}							                 									     
	for (n=0; n<iabs(jjump); n++){
		pctY[n]=LinInterp(0.0,1.0,Z1.imag(),Z2.imag(),gridY[j1+Yjumpdir*n+add]);//cout <<"pctY["<<n<<"]: "<<pctY[n]<<endl;
	}
	
	//sort intersections, indices from z1 to z2------------------------------
	//-----------------------------------------------------------------------
	itmp=0;
	jtmp=0;	
	i=new int[nint+2];
	j=new int[nint+2];
	i[0]=i1;  j[0]=j1;  k[0]=k1; zint[0]=z1;//initial index same as starting cell 
	
	for (n=0; n<nint; n++){
		if			((jtmp>=iabs(jjump)) ||           //j jumps done or
			       ((itmp<iabs(ijump)) &&           
						  (jtmp<iabs(jjump)) &&           //next j is sooner than next i
						  (pctX[itmp]<pctY[jtmp]))){															 //next X crossing happened before next Y crossing
			zint[n+1]=pctX[itmp]*(Z2-Z1)+Z1;                                    //Zint interpolated
      i   [n+1]=i[n]+Xjumpdir;                                            //i increased/decreased by 1
			j   [n+1]=j[n];                                                     //j stays the same
			itmp++;                                                             //current i moves forward
			//cout <<"a";
		}
		else if ((itmp>=iabs(ijump)) ||           //i jumps done or
						 ((itmp<iabs(ijump)) &&           
						  (jtmp<iabs(jjump)) &&           //next i is sooner than next j
			        (pctX[itmp]>pctY[jtmp]))){                                //next Y crossing happened before next X crossing
			zint[n+1]=pctY[jtmp]*(Z2-Z1)+Z1;																		//Zint interpolated
			i   [n+1]=i[n];																											//i stays the same
      j   [n+1]=j[n]+Yjumpdir;																						//j increased/decreased by 1
			jtmp++;																															//current j moves forward
			//cout <<"b";
		}
		else if ((iabs(ijump)>itmp) && 
						 (iabs(jjump)>jtmp) &&
						 (fabs(pctX[itmp]-pctY[jtmp])<REALSMALL)){                  //X and Y crossing happen at same time (rare)
			zint[n+1]=pctY[jtmp]*(Z2-Z1)+Z1;																		//Zint interpolated
			i   [n+1]=i[n]+Xjumpdir;																						//i increased/decreased by 1
      j   [n+1]=j[n]+Yjumpdir;																						//j increased/decreased by 1
			jtmp++;																															//current j and i move forward
			itmp++;																															
      nint--;                                                             //number of intersections decreased by one (this one doubly counted)
			//cout <<"c"<<pctX[itmp]<<"|"<<pctY[jtmp]<<"|";
		}
		else{
      //cout <<endl<<"d"<<pctX[itmp]<<"|"<<pctY[jtmp]<<"|"<<itmp<<"|"<<jtmp<<"|"<<iabs(ijump)<<"|"<<iabs(jjump)<<"|"<<endl;
      ExitGracefully("CRectGrid::BoundaryIntercepts: Should never happen...",BAD_DATA);
		}
		zint[n+1]=LocalToGlobal(zint[n+1]);																//convert from local to global coords
		k   [n+1]=ij_to_k(i[n+1],j[n+1]);                                  //translate into general index, k
		if (k[n+1]==OFF){ExitGracefully("CRectGrid::BoundaryIntercepts: BadIndex(3)",BAD_DATA);}
	}
	
	//handle last point---------------------------------
	zint[nint+1] =z2; 
	k   [nint]   =k2; 
	ncellscrossed=nint+1;

	//translate to 3D-----------------------------------
	for (n=0; n<nint+2; n++){
		ptint[n]=c2Dto3D(zint[n]);
	}

	//Clean up----------------
	delete [] pctX;
	delete [] pctY;
	delete [] zint;
	delete [] i;
	delete [] j;
  
	return true;

}
/***************************************************************************
                  WRITE GEOMETRY
****************************************************************************
 Writes grid to Atlas BNA file, rectgrid.bna
--------------------------------------------------------------------------*/
void CRectGrid::WriteGeometry() const{
	//create output
	int i,j;
	cmplex znode;
	cmplex z1;
	cmplex z2;
	ofstream RECTGRID;
	RECTGRID.open("rectgrid.bna");
	if ((gridX==NULL) || (gridY==NULL)){return;}
	for (i=0; i<=nX; i++){
		z1=LocalToGlobal(cmplex(gridX[i],0.0));
		z2=LocalToGlobal(cmplex(gridX[i],Ymax));
    RECTGRID << "\" gridline \", -2" <<endl;
    RECTGRID << z1.real() << " " << z1.imag() <<endl;
		RECTGRID << z2.real() << " " << z2.imag() <<endl;
		/*if (i!=nX){
			for (j=0; j<nY; j++){
				znode=LocalToGlobal(cmplex(nodeX[i],nodeY[j]));
				RECTGRID << "\" cen \", 1" <<endl;
				RECTGRID << znode.real() << " " << znode.imag() <<endl;
			}
		}*/
	}
	for (j=0; j<=nY; j++){
		z1=LocalToGlobal(cmplex(0.0 ,gridY[j]));
		z2=LocalToGlobal(cmplex(Xmax,gridY[j]));
    RECTGRID << "\" gridline \", -2" <<endl;
    RECTGRID << z1.real() << " " << z1.imag() <<endl;
		RECTGRID << z2.real() << " " << z2.imag() <<endl;
	}
  RECTGRID.close();
}
/***************************************************************************
                  WRITE RECTGRID TO BBG FILE
****************************************************************************
 Writes grid to grid.bbg (BlueBird Grid) file 
--------------------------------------------------------------------------*/
void CRectGrid::WriteBBGFile() const{
  ofstream BBG;
	int i;
	BBG.open("grid.bbg");
  BBG<<"FiniteDifferenceGrid BBGenerated"<<endl;
  BBG<<zg.real()<<" "<<zg.imag()<<endl;
	BBG<<orient<<endl;
	BBG<<nX<<" "<<nY<<endl;
	for (i=0; i<=nX; i++){BBG<<gridX[i]<<endl;}
	BBG<<"&"<<endl;
	for (i=0; i<=nY; i++){BBG<<gridY[i]<<endl;}
	BBG<<"&"<<endl;
	BBG.close();
	
};

/***************************************************************************
                  PARSE RECTGRID
****************************************************************************
 Reads grid from *.bbg (BlueBird Grid) file 
--------------------------------------------------------------------------*/
CRectGrid *CRectGrid::ReadBBGFile(ifstream &input,int &l){
	//string "FiniteDifferenceGrid", string label (unused) 
	//double x_bottom_left y_bottom_left
	//double pivot_angle
	//double nX (number of columns) double nY (number of rows)
	//{double X }x(nX+1)
  //& 
	//{double y }x(nY+1)
  //& 

	CRectGrid  *pRectGrid;

  bool     eof(false),done(false);
	double   xp,yp,angle;
	int      Len,nX,nY,i;
  double  *X=NULL;
	double  *Y=NULL;
	char    *s[MAXINPUTITEMS];
	pRectGrid=NULL;

	if (parserdebug) {cout << "Rectangular Grid"<<endl;}

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
    if      (Len==2) {nX=s_to_i(s[0]); 
											nY=s_to_i(s[1]); 
											done=true;												 
											X=new double [nX+1];
											Y=new double [nY+1];
											ExitGracefullyIf(Y==NULL,"CRectGrid::Parse",OUT_OF_MEMORY);
																												 eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line "<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	i=0;
  do {
		if ((Len==1) && (strcmp(s[0],"&"))){
			if (i>nX)  {ExitGracefully("CRectGrid::Parse- too many columns in grid",TOO_MANY);}
			X[i]=s_to_d(s[0]); i++;                            eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			if (i==nX+1){done=true;}
			else        {ExitGracefully("CRectGrid::Parse- not enough columns in grid",BAD_DATA);}  
		}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	i=0;                                                   eof=TokenizeLine(input,s,Len); l++;
  do {
		if ((Len==1) && (strcmp(s[0],"&"))){
			if (i>nY)  {cout <<i<<" "<<nY<<endl;
				ExitGracefully("CRectGrid::Parse- too many rows in grid",TOO_MANY);}
			Y[i]=s_to_d(s[0]); i++;                            eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			if  (i==nY+1) {
				pRectGrid = new CRectGrid(cmplex(xp,yp),
																	angle,
																	X,
																	Y,	
																	nX,
																	nY); 
				done=true;
			}
			else{cout <<i<<" "<<nY<<endl; ExitGracefully("CRectGrid::Parse- not enough rows in grid",BAD_DATA);}  
		}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	delete [] X;
	delete [] Y;
	if (eof) {return NULL;}
	else     {return pRectGrid;}
}
/***************************************************************************
                  PARSE RECTGRID
****************************************************************************
--------------------------------------------------------------------------*/
CRectGrid *CRectGrid::ParseSpacing(ifstream &input,int &l){
	//string "GridSpacingFunction", string label (unused) 
	//double x_bottom_left y_bottom_left
	//double pivot_angle
	//double x-length, double y-length
	//double minspacing, double maxspacing
	//double nX (number of columns) double nY (number of rows)
	//{double X double Fx}xNX
  //& 
	//{double Y double Fy}xNY
  //& 

	CRectGrid  *pRectGrid;

  bool     eof(false),done(false);
	double   xpt,ypt,angle;
	int      Len,nX,nY,NX,NY,i;
	double   Xlength,Ylength;
	double   minspacing, maxspacing;
  double  *X=NULL;
	double  *Y=NULL;
	double  *xp=NULL;double *fx=NULL;
	double  *yp=NULL;double *fy=NULL;
	char    *s[MAXINPUTITEMS];
	pRectGrid=NULL;

	if (parserdebug) {cout << "Rectangular Grid"<<endl;}

	//starts at line #2
	//eof=TokenizeLine(input,s,Len); l++; 
	eof=TokenizeLine(input,s,Len); l++; 
	do{
    if      (Len==2) {xpt=s_to_d(s[0]); 
											ypt=s_to_d(s[1]); 
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
    if      (Len==2) {Xlength=s_to_d(s[0]); 
											Ylength=s_to_d(s[1]); 
											done=true;												 eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line "<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	do{
    if      (Len==2) {minspacing=s_to_d(s[0]); 
											maxspacing=s_to_d(s[1]); 
											ExitGracefullyIf(minspacing<=0,"CRectGrid:Parse: negative or zero minimum grid spacing",BAD_DATA);
											done=true;												 eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line "<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	do{ 
    if      (Len==2) {NX=s_to_i(s[0]); 
											NY=s_to_i(s[1]); 
											done=true;												 
											xp=new double [NX];
											fx=new double [NX];
											yp=new double [NY];
											fy=new double [NY];
											ExitGracefullyIf(fy==NULL,"CRectGrid::ParseSpacing",OUT_OF_MEMORY);
																												 eof=TokenizeLine(input,s,Len); l++;}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line "<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	i=0;
  do {
		if (Len==2) {
			if (i>NX) {ExitGracefully("CRectGrid::ParseSpacing- too many points in polynomial spacing function(x)",TOO_MANY);}
			xp[i]=s_to_d(s[0]); fx[i]=s_to_d(s[1]);i++;        eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			if (i==NX){done=true;}
			else      {ExitGracefully("CRectGrid::ParseSpacing- not enough points in polynomial spacing function(x)",BAD_DATA);}  
		}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else        {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;
	i=0;                                                   eof=TokenizeLine(input,s,Len); l++;
  do {
		if (Len==2){
			if (i>NY)  {
				ExitGracefully("CRectGrid::ParseSpacing- too many points in polynomial spacing function(y)",TOO_MANY);}
			yp[i]=s_to_d(s[0]); fy[i]=s_to_d(s[1]);i++;        eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&"))) {
			if  (i==NY) {done=true;}
			else{ExitGracefully("CRectGrid::ParseSpacing- not enough points in polynomial spacing function(y)",BAD_DATA);}  
		}
		else if (Len==0) {                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
  double *ax=new double [NX];
	double *ay=new double [NY];
  X=new double [(int)(Xlength/minspacing)+1];
  Y=new double [(int)(Ylength/minspacing)+1];
	if (Y==NULL){ExitGracefully("CRectGrid::ParseSpacing: outofmem",OUT_OF_MEMORY);}
	PolyFit(xp,fx,ax,NX);
	PolyFit(yp,fy,ay,NY);
	X[0]=0.0;
	Y[0]=0.0;
	nX=nY=0;
	ofstream SPC;
	SPC.open("spacedebug.csv");
  //prints polynomial spacing function
	for (double tmpx=0; tmpx<Xlength; tmpx+=Xlength/100){
		SPC<<tmpx<<","<<PolyEval(tmpx,ax,NX)<<endl;
	}
	for (double tmpy=0; tmpy<Ylength; tmpy+=Ylength/100){
		SPC<<tmpy<<","<<PolyEval(tmpy,ay,NY)<<endl;
	}
	while (X[nX]<Xlength){
		nX++; 
		X[nX]=X[nX-1]+min(max(PolyEval(     X[nX-1]         ,ax,NX),minspacing),maxspacing);
		X[nX]=X[nX-1]+min(max(PolyEval(0.5*(X[nX-1]+X[nX  ]),ax,NX),minspacing),maxspacing);//predictor-corrector
		
		ExitGracefullyIf(nX>=(int)(Xlength/minspacing),"CRectGrid::ParseSpacing: array overflow(x)",RUNTIME_ERR);
    //SPC<<X[nX]<<","<<PolyEval(X[nX-1]*Xlength,ax,NX)/Xlength<<","<<minspacing/Xlength<<","<< maxspacing/Xlength<<endl;
	} 
	while (Y[nY]<Ylength){
		nY++; 
		Y[nY]=Y[nY-1]+min(max(PolyEval(     Y[nY-1]         ,ay,NY),minspacing),maxspacing);
		Y[nY]=Y[nY-1]+min(max(PolyEval(0.5*(Y[nY-1]+Y[nY  ]),ay,NY),minspacing),maxspacing);//predictor-corrector
		
		ExitGracefullyIf(nY>=(int)(Ylength/minspacing),"CRectGrid::ParseSpacing: array overflow(y)",RUNTIME_ERR);
	  //SPC<<Y[nY]<<","<<PolyEval(Y[nY-1]*Ylength,ax,NY)/Ylength<<endl;
	} 

	SPC.close();

	for (i=0; i<=nX; i++){X[i]/=(X[nX]/Xlength);}
	for (i=0; i<=nY; i++){Y[i]/=(Y[nY]/Ylength);}
	cout << "X[nX]:"<<X[nX]<<" "<<Xlength<<endl;
	cout << "Y[nY]:"<<Y[nY]<<" "<<Ylength<<endl;

	pRectGrid = new CRectGrid(cmplex(xpt,ypt),angle,X,Y,nX,nY);
	//pRectGrid->WriteGeometry();//TMP DEBUG
	pRectGrid->WriteBBGFile();//writes grid.bbg

  delete [] ax;
  delete [] ay;
	delete [] xp;
	delete [] yp;
	delete [] fx;
	delete [] fy;
	delete [] X;
	delete [] Y;
	//ExitGracefully("EXITING FROM CRECTGRID:PARSESPACING",OTHER);//TMP DEBUG 
	if (eof) {return NULL;}
	else     {return pRectGrid;}
}