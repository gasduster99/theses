//FEMesh.cpp

//DelunayMesh.cpp
#include "TriagonalMesh.h"
/*******************************************************************************
           CONSTRUCTORS
********************************************************************************
------------------------------------------------------------------------------*/
CFEMesh::CFEMesh(char         *Name):
				 CTriMesh(Name){
	meshtype=NODE_BASED;
}
//------------------------------------------------------------------------------
CFEMesh::CFEMesh(char         *Name,
								 const cmplex *pts, 
								 const int     numpts):
				 CTriMesh    (Name,pts,numpts){
	meshtype=NODE_BASED;
}
//------------------------------------------------------------------------------
CFEMesh::CFEMesh(char              *Name,
								 const nodestruct  *points, 
								 const int          nP,
								 const segstruct   *segments, 
								 const int          nS,
								 const int          maximum_nodes,
					       const spcstruct    spacing_package):
				 CTriMesh(Name,points,nP,segments,nS,maximum_nodes,spacing_package){
	meshtype=NODE_BASED;
}
//------------------------------------------------------------------------------
CFEMesh::~CFEMesh(){
}
/*******************************************************************************
           PROCESS
********************************************************************************
	Adds information to make this a finite element mesh
	Adds centroid information for elements
	Adds area information for elements
	Calculates A & B for sides
------------------------------------------------------------------------------*/
void CFEMesh::Process(){
	int k;
	double x1,y1,x2,y2;
 
	nNodes=NumNodes; //node-based method->no cells
	nCells=0;

	if (DOUBLENODETEST) {
		nIntNodes=0;
		for (k=0;k<NumDoublenodes;k++){
			nIntNodes+=dnodes[k].nsegs;
		}
		nNodes=NumNodes+nIntNodes;
	}
	
	//---Identifying Element Centroids--------------------------------
	for (k=0;k<NumElems;k++){ 
		elem[k].centroid=node[elem[k].ni].z+(node[elem[k].nj].z+node[elem[k].nk].z-2.0*node[elem[k].ni].z)/3.0;
	}
	//---Identifying Element Areas--------------------------------
	for (k=0;k<NumElems;k++){ 
		elem[k].area=TriArea(node[elem[k].ni].z,node[elem[k].nj].z,node[elem[k].nk].z);
	}

	for (int s=0; s<NumSides; s++){
		x1=node[side[s].n1].z.real();y1=node[side[s].n1].z.imag();
		x2=node[side[s].n2].z.real();y2=node[side[s].n2].z.imag();
		side[s].A=(y2-y1)/(x1*y2-x2*y1);		
		side[s].B=(x1-x2)/(x1*y2-x2*y1);
	}
}
/*************************************************************************
									Inherited Abstract Mesh Functions
**************************************************************************
-------------------------------------------------------------------------*/  
//--------------------------------------------------------------------------
pt3D CFEMesh::GetNodeLocation (const int k) const{
	if (DOUBLENODETEST) {
		int i,count(0);
		if      (k<NumNodes                   ){ return c2Dto3D(node[k].z);}//regular node
		else if (k<NumNodes+nIntNodes){ 
			for (i=0; i<NumDoublenodes;i++){
        if (k<count+dnodes[i].nsegs){return c2Dto3D(node[dnodes[i].n1].z);}
				count+=dnodes[i].nsegs;
			}
			ExitGracefully("CFEMesh::GetNodeLocation: doublenode index error (bad count)",RUNTIME_ERR);
		}//interstitial node
		else {//bad index
			ExitGracefully("CFEMesh::GetNodeLocation: doublenode index error",RUNTIME_ERR);
		}
	}
	return c2Dto3D(node[k].z);
}
//--------------------------------------------------------------------------
double CFEMesh::InterpolateValue(const pt3D &pt, Ironclad1DArray values) const{
	static int    e;                //declared static for speed
	static double Ni,Nj,Nk;
	static cmplex z;
	z=c3Dto2D(pt);

	//linear interpolation between nodes
  e=GetElement(z);

	if (e!=OFF){
		
		elemGlobalToLocal(e,z,Ni,Nj,Nk);

		return Ni*values[elem[e].ni]+
					 Nj*values[elem[e].nj]+
					 Nk*values[elem[e].nk];

	}
	return 0.0;
}                                 
//--------------------------------------------------------------------------
void CFEMesh::InterpolateValues(const pt3D       &pt, 
																Ironclad2DArray   values, 
																Writeable1DArray  intval, 
																const int         nval) const{
	static int    e;
	static double Ni,Nj,Nk;              //declared static for speed
	static cmplex z;
	z=c3Dto2D(pt);

	//linear interpolation between nodes
  e=GetElement(z);

	if (e!=OFF){
    
		elemGlobalToLocal(e,z,Ni,Nj,Nk);

	  for (int i=0;i<nval; i++){
			intval[i]=Ni*values[elem[e].ni][i]+
				        Nj*values[elem[e].nj][i]+
								Nk*values[elem[e].nk][i];
		}
		return;
	}
	for (int i=0;i<nval; i++){intval[i]=0.0;}
}
//-------------------------------------------------------------------------- 
double CFEMesh::InterpolateDerivative    (const pt3D &pt, Ironclad1DArray values, const cmplex &dir) const{
	//TMP DEBUG- STUB FUNCTION. MUST FILL
	return 0.0;
}
//--------------------------------------------------------------------------
double CFEMesh::GetN(const int e, const cmplex &z, const int i) const{
	//Basis function equal to local coordinate
	static double xii,xij,xik;
	elemGlobalToLocal(e, z, xii,xij, xik);
	
	if      (i==0){return xii;}
	else if (i==1){return xij;}
	else if (i==2){return xik;}
	else          {ExitGracefully("CFEMesh::GetN: Bad basis function request",RUNTIME_ERR);return 0.0;}
}
//--------------------------------------------------------------------------
cmplex CFEMesh::GetdN(const int e, const cmplex &z, const int i) const{
	//returns dNdx + i dNdy
	//slope independent of z with linear basis function

  static double d[3];
	static int    ind[3];
	static cmplex dN;
	static double R,junk,chi,eta;

  static cmplex zt[3],Z[3];

	zt[0] =node[elem[e].ni].z;  
	zt[1] =node[elem[e].nj].z;  
	zt[2] =node[elem[e].nk].z; 
	
	R=GetRadiusAndRelativeIndices(e,ind);

	if (R==0)
	{
		if      (i==0){dN= conj(zt[2]-zt[1]);}
		else if (i==1){dN= conj(zt[0]-zt[2]);}
		else if (i==2){dN= conj(zt[1]-zt[0]);}
		else          {
			ExitGracefully("CFEMesh::GetdNdx: Bad basis function request",RUNTIME_ERR);return 0.0;}
	}
	else
	{
    for (int i1=0;i1<3;i1++){Z[i1]=zt[ind[i1]];}

		Get_d_coeff(Z,R,d);

		//cout <<d[0]<<" "<<d[1]<<" "<<d[2]<<endl;
		//d[0]=0;d[1]=0;d[2]=0;
		TriGlobalToLocal(Z[0],Z[1],Z[2],z,junk,chi,eta);

		if      (i==ind[0])   //dNi/dz
		{ 
			dN= (conj(Z[1]-Z[0])*(d[0]*eta-1.0)+conj(Z[0]-Z[2])*(d[0]*chi-1.0));
		}
		else if (i==ind[1])   //dNj/dz
		{
			dN= (conj(Z[1]-Z[0])*(d[1]*eta    )+conj(Z[0]-Z[2])*(d[1]*chi+1.0));
		}
		else if (i==ind[2])   //dNk/dz
		{
			dN= (conj(Z[1]-Z[0])*(d[2]*eta+1.0)+conj(Z[0]-Z[2])*(d[2]*chi    ));
		}
		else          {
			ExitGracefully("CFEMesh::GetdNdx: Bad basis function request(2)",RUNTIME_ERR);return 0.0;}
	}

	return cmplex(dN.imag(),dN.real())/(2.0*elem[e].area);
}
/***************************************************************************
                  GET NEIGHBORS
****************************************************************************
input:   cell index (k)
				 empty integer array to store neighbor values (pre-allocated)
				 empty double array to store distance to neighbors (pre-allocated)
				 maximum number of neighbors (numneighbors)
output:  neighbor indices (neighbors), numneighbors (maximum 3 for trimesh)
--------------------------------------------------------------------------*/
void CFEMesh::GetNeighbors  (const int k, int *neighbors, Writeable1DArray dist, int &numneighbors) const{
	int count(0);

	if (numneighbors<8){ExitGracefully("CFEMesh::GetNeighbors: not enough space allocated", RUNTIME_ERR);}

	numneighbors=count;
}
/***************************************************************************
                  GAUSS INTEGRATION
****************************************************************************

--------------------------------------------------------------------------*/
double CFEMesh::IntegrateFunction(const int e, 
																  const trigausspoints gausspts, 
																  double (*F)(const cmplex &z)) const{

  return 0.0;
}
/***************************************************************************
                 COORDINATE TRANSFORMATION
****************************************************************************
	Returns Radius of sides and indices relative to fulcrum node 
	(e.g., if curved side is j, then index zero is j, index one is k, index two is i)  
--------------------------------------------------------------------------*/
double CFEMesh::GetRadiusAndRelativeIndices(const int e, int ind[3]) const{
  static double R;
	R=0;ind[0]=0;ind[1]=1;ind[2]=2;
	if      (side[elem[e].si].radius!=0)
	{
		R     =side[elem[e].si].radius;   
		
		if (side[elem[e].si].n1!=elem[e].nj){R*=-1.0;}

		ind[0]=0;ind[1]=1;ind[2]=2;
	}
	else if (side[elem[e].sj].radius!=0)
	{
		R     =side[elem[e].sj].radius;   

		if (side[elem[e].sj].n1!=elem[e].nk){R*=-1.0;}

		ind[0]=1;ind[1]=2;ind[2]=0;
	}
	else if (side[elem[e].sk].radius!=0)
	{
		R     =side[elem[e].sk].radius;   
		
		if (side[elem[e].sk].n1!=elem[e].ni){R*=-1.0;}

		ind[0]=2;ind[1]=0;ind[2]=1;	
	}
	return R;
}
/***************************************************************************
                 COORDINATE TRANSFORMATION
****************************************************************************
	Converts back and forth to barycentric coordinates
--------------------------------------------------------------------------*/
cmplex CFEMesh::elemLocalToGlobal(const int e, const double &xii, const double &xij, const double &xik) const{

	static cmplex z;
  static cmplex zt[3];
  static int    ind[3];
	static double Rprime,R,L;
  static cmplex bci,bci2;

	zt[0] =node[elem[e].ni].z;  zt[0] +=0.001*(elem[e].centroid-zt[0]);
	zt[1] =node[elem[e].nj].z;  zt[1] +=0.001*(elem[e].centroid-zt[1]);
	zt[2] =node[elem[e].nk].z;  zt[2] +=0.001*(elem[e].centroid-zt[2]);

	z= xii*zt[0]+xij*zt[1]+xik*zt[2];

  R=GetRadiusAndRelativeIndices(e,ind);

	if (R!=0.0){
		L			=abs(zt[ind[1]]-zt[ind[2]]);
		Rprime=R/L*(1.0-sqrt(1.0-L*L/(4.0*R*R)));

		bci=(zt[ind[2]]-zt[ind[1]]);
		bci2=cmplex(bci.imag(),-bci.real());

		if			(ind[0]==0){//node i is fulcrum
			z-=xij*xik*4.0*Rprime*bci2;
		}
		else if (ind[0]==1){//node j is fulcrum
			z-=xii*xik*4.0*Rprime*bci2;
		}
		else if (ind[0]==2){//node k is fulcrum
			z-=xii*xij*4.0*Rprime*bci2;
		}
	}

	return z;
}
//--------------------------------------------------------------------------
void CFEMesh::elemGlobalToLocal(const int e, const cmplex &z, double &xii,double &xij, double &xik) const{

  static cmplex zt[3],Z[3];
  static int    ind[3];
	static double d[3];
	static double R,chieta;

	zt[0] =node[elem[e].ni].z;  
	zt[1] =node[elem[e].nj].z;  
	zt[2] =node[elem[e].nk].z;  
  
	TriGlobalToLocal(zt[0],zt[1],zt[2],z,xii,xij,xik);
  
  R=GetRadiusAndRelativeIndices(e,ind);

	if (R!=0){
    for (int i=0;i<3;i++){Z[i]=zt[ind[i]];}

		Get_d_coeff(Z,R,d);

		chieta =xij*xik;
																//desired  ind    ind[ind]
		xii+=d[ind[ind[0]]]*chieta; //0 2 1    0 1 2  0 2 1
		xij+=d[ind[ind[1]]]*chieta; //1 0 2    1 2 0  1 0 2
		xik+=d[ind[ind[2]]]*chieta; //2 1 0    2 0 1  2 1 0
	}

	return;
}
//--------------------------------------------------------------------------
void CFEMesh::Get_d_coeff(const cmplex Z[3], const double R, double d[3]) const{
	//Z[] are pre-transformed such that Z[0] is the fulcrum node (side [0] is curved)
	static double LoverR,Rprime,chi4,eta4,chieta4,junk;
	static cmplex z4;

	//Calculate di,dj,dk (need zt,R)
	LoverR=abs(Z[1]-Z[2])/R;
	Rprime=(1.0-sqrt(1.0-0.25*LoverR*LoverR))/LoverR;

	z4=0.5*(Z[1]+Z[2])+Rprime*conj(Z[1]-Z[2]);

	TriGlobalToLocal(Z[0],Z[1],Z[2],z4,junk,chi4,eta4);
	//cout <<"L:"<<L<<" R':"<<Rprime<<" chi4:"<<chi4<<" "<<eta4<<endl;

	chieta4 =(chi4*eta4);
	d[0]=(chi4+eta4-1.0)/chieta4;
	d[1]=(0.5-eta4     )/chieta4;
	d[2]=(0.5-chi4     )/chieta4;
}
//--------------------------------------------------------------------------
double CFEMesh::GetElemJacobianDet(const int e, 
																	 const cmplex &z) const{
  static double Jdet;
  static double R;
	static int    ind[3];

  R=GetRadiusAndRelativeIndices(e,ind);

	Jdet=2.0*elem[e].area;

	if (R!=0){
		static cmplex zt[3],Z[3];
		static double xi[3],d[3],bi,ci,bk,ck;
		static double junk,chi,eta,LoverR,Rprime;
		static cmplex bci,bki;

		zt[0] =node[elem[e].ni].z;  
		zt[1] =node[elem[e].nj].z;  
		zt[2] =node[elem[e].nk].z;  

    for (int i=0;i<3;i++){Z[i]=zt[ind[i]];}

		TriGlobalToLocal(Z[0],Z[1],Z[2],z,junk,chi,eta);		

		LoverR=abs(Z[1]-Z[2])/R;
		Rprime=(1.0-sqrt(1.0-0.25*LoverR*LoverR))/LoverR;

		//Calculate di,dj,dk 
		Get_d_coeff(Z,R,d);

		//Calculate bi,bk,ci,ck
		bci=conj(Z[2]-Z[1]); //=ci+i bi
		bki=conj(Z[1]-Z[0]); //=ck+i bk
		bi=bci.imag();
		ci=bci.real();
		bk=bki.imag();
		ck=bki.real();
		Jdet-=2.0*Rprime*((bi*bk+ci*ck)*(eta*(1+d[1]*chi*eta)));//+
				  //        (ci*cj+bi*bj)*(chi*(1+d[2]*chi*eta)));
	}
	return Jdet;
}
//--------------------------------------------------------------------------
void CFEMesh::GetElemJacobian (const int e, 
															 const double &xii, const double &xij, const double &xik, 
															 double J[2][2],
															 double Jinv[2][2],
															 double &Jdet)const{
	//TMP DEBUG: CURVATURE REVISION
	cmplex zi=node[elem[e].ni].z;
	cmplex zj=node[elem[e].nj].z;
	cmplex zk=node[elem[e].nk].z;



	Jdet=2.0*elem[e].area;//J[0][0]*J[1][1]-J[0][1]*J[1][0];

	J[0][0]=(zj-zi).real();//dx/dLi
	J[0][1]=(zk-zi).real();//dx/dLj
	J[1][0]=(zj-zi).imag();//dy/dLi
	J[1][1]=(zk-zi).imag();//dy/dLj


	Jinv[0][0]= J[1][1]/Jdet;
	Jinv[0][1]=-J[0][1]/Jdet;
	Jinv[1][0]=-J[1][0]/Jdet;
	Jinv[1][1]= J[0][0]/Jdet;
}
