//VoronoiMesh.cpp
#include "TriagonalMesh.h"
/*******************************************************************************
           CONSTRUCTORS
********************************************************************************
------------------------------------------------------------------------------*/
CVoronoiMesh::CVoronoiMesh(char         *Name,
													 const cmplex *pts, 
													 const int     numpts):
							CTriMesh    (Name,pts,numpts){
	voronoi=NULL;
	meshtype=CELL_BASED;
}
//------------------------------------------------------------------------------
CVoronoiMesh::CVoronoiMesh(char              *Name,
													 const nodestruct  *points,
													 const int          nP,
													 const segstruct   *segments,
													 //const int         *segp0,
													 //const int         *segp1,
												   //const status_type *segstatus,
													 const int          nS,
													 const int          maximum_nodes,
													 const spcstruct    spacing_package):
							CTriMesh    (Name,points,nP,segments,nS,maximum_nodes,spacing_package){
//							CTriMesh    (Name,points,nP,segp0,segp1,
//												   segstatus,nS,maximum_nodes,spacing_package){
	voronoi=NULL;
	meshtype=CELL_BASED;
}
//------------------------------------------------------------------------------
CVoronoiMesh::~CVoronoiMesh(){
	for (int v=0; v<MaxNodes;v++){
		if (voronoi!=NULL){
			if (voronoi[v].n!=NULL){delete [] voronoi[v].n;}
			if (voronoi[v].e!=NULL){delete [] voronoi[v].e;}
			if (voronoi[v].zv!=NULL){delete [] voronoi[v].zv;}
		}
	}
	delete [] voronoi;
	meshtype=CELL_BASED;
}
//***********************************************************************
void CVoronoiMesh::Process(){
	nNodes=nCells=NumNodes;
	
	int k;
	voronoi=new vorstruct [nCells];

	//---Identifying Element Centroids--------------------------------
	for (k=0;k<NumElems;k++){ 
		elem[k].centroid=node[elem[k].ni].z+(node[elem[k].nj].z+node[elem[k].nk].z-2.0*node[elem[k].ni].z)/3.0;
	}
	//---Identifying Element Areas--------------------------------
	for (k=0;k<NumElems;k++){ 
		elem[k].area=TriArea(node[elem[k].ni].z,node[elem[k].nj].z,node[elem[k].nk].z);
	}

	for(k=0; k<nCells; k++){
		voronoi[k].area =0.0;           
		voronoi[k].n    =NULL; 
		voronoi[k].e    =NULL;
		voronoi[k].nVert=0;
		voronoi[k].z    =0.0;
		voronoi[k].zv   =NULL; 

		CalcVoronoi(voronoi[k],node[k].z,k); 
	
	}
}
/*************************************************************************
									Inherited Abstract Mesh Functions
**************************************************************************
-------------------------------------------------------------------------*/
pt3D CVoronoiMesh::GetNodeLocation(const int k) const{
	//takes global index (not local)
	return c2Dto3D(node[k].z);//TMP DEBUG 
}     
//---------------------------------------------------------------------------
double CVoronoiMesh::GetCellArea(const int k) const{
	if ((k<0) || (k>nCells)){ExitGracefully("CVoronoiMesh::GetCellVolume: Bad Cell index",RUNTIME_ERR);}
	return voronoi[k].area;
}          
//--------------------------------------------------------------------------                                       
bool   CVoronoiMesh::GetCellIndex             (const pt3D &pt, int &k) const{
	//NEEDS OPTIMIZATION
	if (!IsInside(pt)){k=OFF;return false;}
	double mindist;
	cmplex z=c3Dto2D(pt);
	mindist=ALMOST_INF;
	for (int v=0;v<NumNodes; v++){      //searches through all- needs optimization
		if (abs(voronoi[v].z-z)<mindist){
			mindist=abs(voronoi[v].z-z);
			k=v;
		}
	}
	return true;
}    
//--------------------------------------------------------------------------                             
bool   CVoronoiMesh::IsInCell                 (const pt3D &pt, const int k) const{
	//k is global index, NOT local v
	if ((k<0) || (k>=NumNodes)){return false;}
	cmplex z(c3Dto2D(pt));

	double sum(0.0);
	for (int i=0; i<voronoi[k].nVert-1; i++){
		sum+=TriArea(z, voronoi[k].zv[i] , voronoi[k].zv[i+1]);
	}
	sum+=TriArea(z, voronoi[k].zv[voronoi[k].nVert-1] , voronoi[k].zv[0]);
	if (sum<0){return false;}
	else      {return true;}
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
bool   CVoronoiMesh::GetMeshBoundaryIntercepts(const pt3D &pt1, const pt3D &pt2, pt3D *ptint, int *k, int &nint) const{

	return false; //TMP DEBUG -must fill out
}
//--------------------------------------------------------------------------
double CVoronoiMesh::InterpolateValue(const pt3D &pt, Ironclad1DArray values) const{

	//linear interpolation on triangle
	if (IsInside(pt)){
		int e=GetElement(c3Dto2D(pt));
		//find voronoi indices of element corners, interpolate 
		//TMP DEBUG
		return 0.0;
	}
	else { //not on mesh, cannot interpolate
		return 0.0;
	}
}                                 
//--------------------------------------------------------------------------
void CVoronoiMesh::InterpolateValues(const pt3D       &pt, 
																		 Ironclad2DArray   values, 
																		 Writeable1DArray  intval, 
																		 const int         nval) const{

	//linear interpolation on triangle
	int n;
	if (IsInside(pt)){
		int e=GetElement(c3Dto2D(pt));
		//find voronoi indices of element corners, interpolate 
		//TMP DEBUG
	}
	else { //not on mesh, cannot interpolate
		for (n=0;n<nval;n++){intval[n]= 0.0;}
	}
} 
/***************************************************************************
                  GET NEIGHBORS
****************************************************************************
input:   voronoi cell index (k)
				 empty integer array to store neighbor values (pre-allocated)
				 empty double array to store distance to neighbors (pre-allocated)
				 maximum number of neighbors (numneighbors)
output:  neighbor indices (neighbors), numneighbors 
--------------------------------------------------------------------------*/
void CVoronoiMesh::GetNeighbors  (const int k, int *neighbors, double *dist, int &numneighbors) const{
	if (numneighbors<voronoi[k].nVert){ExitGracefully("CVoronoiMesh::GetNeighbors: not enough space allocated", RUNTIME_ERR);}
	numneighbors=voronoi[k].nVert;
	for (int i=0; i<voronoi[k].nVert; i++){
		neighbors[i]=voronoi[k].n[i];
		dist[i]=abs(voronoi[k].z-voronoi[voronoi[k].n[i]].z);
	}
}
/**************************************************************************
														Calculate Voronoi
***************************************************************************
  This function calculates voronoi region for node n
--------------------------------------------------------------------------*/
void CVoronoiMesh::CalcVoronoi(vorstruct &vor, 
															 const cmplex &zp, 
															 const int n){

	int sv[100]; //TMP DEBUG- should be hardcoded as max_vor_neigh
	int i,s;
	int ea,eb;
	cout <<"**************************"<<endl;
 
	vor.z=zp;

	if (n==OFF){
	  //find all circumcenters this node belongs to...
	}
	else{
		//find number of sides/vertices of voronoi cell
		vor.nVert=0; 
		for (s=0; s<NumSides; s++){
			if (side[s].status!=IS_OFF){
				if ((n==side[s].n1) && (side[s].eb!=OFF)){
					sv[vor.nVert]=s;
					vor.nVert++;
				}
				else if ((n==side[s].n2)&& (side[s].ea!=OFF)){
					sv[vor.nVert]=s;
					vor.nVert++;
				}
			}
		}
		cout <<"z:"<<vor.z<<endl;
		cout <<"nV:"<<vor.nVert<<endl;

		//dimension arrays
		vor.zv=new cmplex [vor.nVert];
		vor.n =new int    [vor.nVert];
		vor.e =new int    [vor.nVert];
    
		if (vor.e==NULL){ExitGracefully("CalcVoronoi: Out of Memory",OUT_OF_MEMORY);}

		//fill arrays 
		vor.area=0;

		for (i=0; i<vor.nVert; i++){
			if ((sv[i]<0) || (sv[i]>NumSides)){ExitGracefully("CalcVoronoi: bad side",RUNTIME_ERR);}
			
			if      (n==side[sv[i]].n1){
				eb=side[sv[i]].eb;
				vor.zv[i]=elem[eb].zv;
				vor.n [i]=side[sv[i]].n2;
				vor.e [i]=eb;
				if(eb==OFF)   {cout << "shouldnt happen"<<endl;vor.zv[i]=0.5*(node[side[sv[i]].n1].z+node[side[sv[i]].n2].z);}
				//else          {vor.zv[i]=vor.z;      }
			}
			else if (n==side[sv[i]].n2){
				ea=side[sv[i]].ea;
				vor.zv[i]=elem[ea].zv;
				vor.n [i]=side[sv[i]].n1;
				vor.e [i]=ea;
				if(ea==OFF){cout << "shouldnt happen"<<endl; vor.zv[i]=0.5*(node[side[sv[i]].n1].z+node[side[sv[i]].n2].z);}
				//else          {vor.zv[i]=vor.z;}		
		
			}
		}

		cout <<"A:"<<vor.area<<endl;
	}
	//must sort vertices in clockwise fashion (straight insertion);
	double a;
	int j;
	cmplex tmpz;
	int  tmpe,tmpn;
	for (j=1; j<vor.nVert; j++){
		a=arg(vor.zv[j]- vor.z);
		tmpz=vor.zv[j];
		tmpn=vor.n[j];
		tmpe=vor.e[j];
		i=j;
		while (i>0 && (arg(vor.zv[i-1]- vor.z)>a)){
			vor.zv[i]=vor.zv[i-1];
			vor.n[i]=vor.n[i-1];
			vor.e[i]=vor.e[i-1];
			i--;
		}
		vor.zv[i]=tmpz;
		vor.n[i]=tmpn;
		vor.e[i]=tmpe;
	}

/*cmplex isect;
int tmps;
	//if one node inside delunay & one outside, move outside to intersection of side and connector
	for (i=0; i<vor.nVert-1; i++){
		tmpe=GetElement(vor.zv[i]);
    if ((tmpe==OFF) || (elem[tmpe].on==false)){//node outside
			tmpe=GetElement(vor.zv[i+1]);

			if ((tmpe!=OFF) && (elem[tmpe].on==true)){ //node inside

				if (((elem[vor.e[i]].si==elem[vor.e[i+1]].si) ||
					   (elem[vor.e[i]].si==elem[vor.e[i+1]].sj) ||
						 (elem[vor.e[i]].si==elem[vor.e[i+1]].sk)) &&
						(elem[vor.e[i]].si!=OFF))                {tmps=elem[vor.e[i]].si;}
				if (((elem[vor.e[i]].sj==elem[vor.e[i+1]].si) ||
					   (elem[vor.e[i]].sj==elem[vor.e[i+1]].sj) ||
						 (elem[vor.e[i]].sj==elem[vor.e[i+1]].sk)) &&
					 (elem[vor.e[i]].sj!=OFF))                 {tmps=elem[vor.e[i]].sj;}
				if (((elem[vor.e[i]].sk==elem[vor.e[i+1]].si) ||
					   (elem[vor.e[i]].sk==elem[vor.e[i+1]].sj) ||
						 (elem[vor.e[i]].sk==elem[vor.e[i+1]].sk)) &&
					 (elem[vor.e[i]].sk!=OFF) )                {tmps=elem[vor.e[i]].sk;}

				if (Intersect(node[side[tmps].n1].z,  //side 
							        node[side[tmps].n2].z,
							        vor.zv[i],               //vertex locations
											vor.zv[i+1],
											isect                   )==INTERSECTED){
							vor.zv[i]=isect;
				}				
			}
		}
	}*/

	vor.area=0;
	for (i=0; i<vor.nVert; i++){
		if (i>1)         {vor.area+=TriArea(vor.z,vor.zv[i],vor.zv[i-1]);}
		if (i==vor.nVert){vor.area+=TriArea(vor.z,vor.zv[i],vor.zv[0]);}
	}
}
void  CVoronoiMesh::WriteGeometry            () const{
		/*ofstream TMP;
	TMP.open("numnodes.dat");
	TMP << "x y N" <<endl;
	for (int v=3; v<NumNodes; v++){
		cout <<voronoi[v].area<<endl;*/
		/*VOR << "\" voronoi \",  "<<-voronoi[v].nVert<<endl;
		for (int i=0;i<voronoi[v].nVert; i++){ 
			VOR <<voronoi[v].zv[i].real()<< " , " <<voronoi[v].zv[i].imag()<<endl;
		}*/
	/*	TMP << voronoi[v].z.real() << " "<<voronoi[v].z.imag()<< " "<<voronoi[v].nVert <<endl;
		VOR << "\" voronoi \",  "<<voronoi[v].nVert+1<<endl;
		for (int i=0;i<voronoi[v].nVert; i++){ 
			VOR <<voronoi[v].zv[i].real()<< " , " <<voronoi[v].zv[i].imag()<<endl;
		}
		VOR <<voronoi[v].zv[0].real()<< " , " <<voronoi[v].zv[0].imag()<<endl;
	
	}
	TMP.close();*/
	//TMP DEBUG- find voronoi cell for single point (node 14)
	/*ofstream TMP;
	TMP.open("temp_vor.bna");
	int n=14;
	cout << node[n].Nne <<endl;
	for (s=0; s<NumSides; s++){
		if (side[s].status !=IS_OFF){
			if ((side[s].n1==n) || (side[s].n2==n)){
				ea=side[s].ea;
				eb=side[s].eb;

				if(ea!=OFF){za=elem[ea].zv;}8
				else       {za=0.5*(node[side[s].n1].z+node[side[s].n2].z);}
 
				if(eb!=OFF){zb=elem[eb].zv;}
				else       {zb=0.5*(node[side[s].n1].z+node[side[s].n2].z);}
				TMP << "\" voronoi \",  "<<-2<<endl;
				TMP <<za.real()<< " , " <<za.imag()<<endl;
				TMP <<zb.real()<< " , " <<zb.imag()<<endl;
			}
		}
	}
	TMP.close();
	*/
}