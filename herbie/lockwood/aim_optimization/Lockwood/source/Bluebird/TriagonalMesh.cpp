//TriMesh class
#include "TriagonalMesh.h"
/*******************************************************************************
					 CONSTRUCTORS
********************************************************************************
Constructor #2- builds mesh/voronoi diagram based upon set of internal points
------------------------------------------------------------------------------*/
CTriMesh::CTriMesh(char  *Name,
									 const cmplex *pts, 
									 const int numpts){

	int e,i,s,n;

	nCells=0;
	nNodes=0;
	if (Name!=NULL){
		name = new char[strlen(Name)+1];
		name=strcpy(name,Name);
	}
	is3D=false;

	spacinginfo.type=USER_SPECIFIED;
	spacinginfo.min_spacing=0;
	spacinginfo.max_spacing=0;

	MaxNodes=numpts+3;
	ugly	 =OFF;		//not used
	elem	 =NULL; 	//created in initialize
	side	 =NULL; 	//created in initialize
	node	 =NULL; 	//created in initialize
	NumElems=0;
	NumSides=0;
	NumNodes=0;

	dnodes =NULL; 	//created in doublenodefilter
	NumDoublenodes=0;

	sumz	 =NULL;  //created in initialize
	sumF	 =NULL;  //created in initialize
	Nne 	 =NULL;  //created in initialize

	//initialize extents
	BBox.e=-ALMOST_INF; BBox.w=ALMOST_INF;
	BBox.n=-ALMOST_INF; BBox.s=ALMOST_INF;
	for (i=0; i<numpts; i++){
		upperswap(BBox.e,pts[i].real());
		upperswap(BBox.n,pts[i].imag());
		lowerswap(BBox.w,pts[i].real());
		lowerswap(BBox.s,pts[i].imag());
	}

	Initialize();

	cout <<" Assembling Vornonoi Mesh..."<<endl;

	for (i=0; i<numpts; i++){
		InsertNode(pts[i],INTERIOR,0.0);
	}

	//--Elimination of sides & nodes associated w/ initial element -------
	for(e=0; e<NumElems; e++){
		if(min(min(elem[e].ni, elem[e].nj), elem[e].nk )<3){//checks if associated with initial element: i,j,k<3)
			elem[e].on=false; 	
		}
	}
	for(s=0; s<3; s++){side[s].status=IS_OFF;}
	for(n=0; n<3; n++){node[n].status=IS_OFF;}

	CleanIndices();

	Renumerate();

	Process();

	cout <<"...Triangular Mesh Assembled with "<<NumNodes << " nodes."<<endl;

	//TMP DEBUG - must create boundary polygon!

	WriteGeometry(true,true);

	ExitGracefully("trimesh",BAD_DATA);

}
/**************************************************************************
Constructor #3- creates empty mesh (later read from file)
-------------------------------------------------------------------------*/
CTriMesh::CTriMesh(char *Name){
	
	nCells=0;
	nNodes=0;
	if (Name!=NULL){
		name = new char[strlen(Name)+1];
		name=strcpy(name,Name);
	}
	is3D=false;

	spacinginfo.type=USER_SPECIFIED;
	spacinginfo.min_spacing=0;
	spacinginfo.max_spacing=0;

	MaxNodes=0.0;
	ugly	 =OFF;		//not used
	elem	 =NULL; 	//created in initialize
	side	 =NULL; 	//created in initialize
	node	 =NULL; 	//created in initialize
	dnodes =NULL;
	NumElems=0;
	NumSides=0;
	NumNodes=0;
	NumDoublenodes=0;

	sumz	 =NULL;  //created in initialize
	sumF	 =NULL;  //created in initialize
	Nne 	 =NULL;  //created in initialize

	//initialize extents
	BBox.e=-ALMOST_INF; BBox.w=ALMOST_INF;
	BBox.n=-ALMOST_INF; BBox.s=ALMOST_INF;

}

//**************************************************************************
CTriMesh::~CTriMesh(){
	if (globaldebug){cout << "  DESTROYING TRIANGULAR MESH"<<endl;}
	delete [] elem;
	delete [] side;
	delete [] node;

	delete [] sumz;
	delete [] sumF;
	delete [] Nne;
}

/*************************************************************************
													Initialize
**************************************************************************
	initializes all nodes, elements and sides
	called once from constructor
	accurate mesh extents (BBox) required
-------------------------------------------------------------------------*/
void CTriMesh::Initialize(){

	elem	 =new elemstruct[MaxNodes*2];//shouldnt this be -2?
	side	 =new sidestruct[MaxNodes*3];//shouldnt this be +3?
	node	 =new nodestruct[MaxNodes];
	sumz	 =new cmplex		[MaxNodes];  //should make members of class for speed
	Nne 	 =new int 			[MaxNodes];
	sumF	 =new double		[MaxNodes];
	if (sumF==NULL){cout <<"out of memory"<<endl;}

	for (int e=0; e<MaxNodes*2; e++){
		elem[e].R=0;														 
		elem[e].r=0;														 
		elem[e].ei=OFF;elem[e].ej=OFF;elem[e].ek=OFF;
		elem[e].si=OFF;elem[e].sj=OFF;elem[e].sk=OFF;
		elem[e].ni =OFF;elem[e].nj =OFF;elem[e].nk =OFF;
		elem[e].on=false; 																										
		elem[e].zin =0;
		elem[e].zv =0;
	}
	for (int s=0;s<MaxNodes*3; s++){
		side[s].na=OFF; side[s].nb=OFF; 			
		side[s].n1=OFF; side[s].n2=OFF;
		side[s].ea=OFF; side[s].eb=OFF; 
		side[s].A =0.0; side[s].B =0.0;
		side[s].status=IS_OFF; 		
		side[s].radius=0.0;
	}
	for(int n=0; n<MaxNodes; n++){
		node[n].F =1.0; 				
		node[n].status	=IS_OFF;
		node[n].z =0;
		node[n].inv=0.0;
	}

	//initialize domain with a single dummy element---------------------------
	//------------------------------------------------------------------------

	cmplex zcent=cmplex( 0.5*(BBox.e+BBox.w),0.5*(BBox.n+BBox.s));
	double width=max((BBox.e-BBox.w),(BBox.n-BBox.s));
 
	//initialize domain with a single element----------------------------
	//first three nodes set up as triangle around bounding square

	NumNodes = 3;

	node[2].z = cmplex(zcent.real(),						 zcent.imag() + 2.8*width); 
	node[0].z = cmplex(zcent.real() - 2.0*width, zcent.imag() - 1.4*width); 
	node[1].z = cmplex(zcent.real() + 2.0*width, zcent.imag() - 1.4*width); 

	NumElems=1;
	elem[0].ni=0; 	elem[0].nj=1; 	elem[0].nk=2;
	elem[0].ei=OFF; elem[0].ej=OFF; elem[0].ek=OFF;
	elem[0].si=1; 	elem[0].sj=2; 	elem[0].sk=0;
 
	NumSides=3;
	side[0].n1=0; 	side[0].n2=1; 	side[0].na=2; 	side[0].nb=OFF; 
	side[1].n1=1; 	side[1].n2=2; 	side[1].na=0; 	side[1].nb=OFF; 
	side[2].n1=0; 	side[2].n2=2; 	side[2].na=OFF; side[2].nb=1;  
	side[0].ea= 0;	side[0].eb=OFF;
	side[1].ea= 0;	side[1].eb=OFF; 
	side[2].ea=OFF; side[2].eb= 0;

	NumDoublenodes=0;

	CalculateCircles(0);
}
/*************************************************************************
													CleanIndices
**************************************************************************
	cleans up all references to bad (or OFF) elements
	called once from constructor
-------------------------------------------------------------------------*/
void CTriMesh::CleanIndices(){
	int e,s;
	bool ea_on, eb_on;

	//--Correction of neighbor indices (i.e. elem[e].ei)--------------------
	for(e=0; e<NumElems; e++){
		if (elem[e].ei!=OFF){
			if(elem[elem[e].ei].on==false) {elem[e].ei=OFF;}}
		if (elem[e].ej!=OFF){
			if(elem[elem[e].ej].on==false) {elem[e].ej=OFF;}}
		if (elem[e].ek!=OFF){
			if(elem[elem[e].ek].on==false) {elem[e].ek=OFF;}}
	}

	//--Correction of neighbor indices (i.e. side[s].ea)--------------------
	for(s=3; s<NumSides; s++){
		ea_on=elem[side[s].ea].on;	
		eb_on=elem[side[s].eb].on;
		if ((!ea_on) && (!eb_on)){side[s].status=IS_OFF;}
		else{
			if			(!ea_on){side[s].ea=OFF; side[s].na=OFF;}
			else if (!eb_on){side[s].eb=OFF; side[s].nb=OFF;}
		}
	}
}
/*************************************************************************
													CreatePolygonControlPoints
**************************************************************************
	given a polygon, creates a mesh. The interior node locations of the mesh are returned.
-------------------------------------------------------------------------*/
bool CTriMesh::CreatePolygonControlPoints(const cmplex *ze,
																					const int 		NLines,
																					cmplex			*&zpoints,
																					int 				 &ncontrol){


	int          i,j;
	double       area(0.0),perim(0.0),maxside(0.0);
	double			 spacing;
	bool				 worked(false);

	CTriMesh    *pMesh   =NULL;
	nodestruct	*points  =new nodestruct[NLines];
	segstruct 	*segments=new segstruct [NLines];
	if (segments==NULL){ExitGracefully("CTriMesh::CreatePolygonControlPoints",OUT_OF_MEMORY);}

	if (NLines<=3){ExitGracefully("CTriMesh::CreatePolygonControlPoints: not enough polygon sides",BAD_DATA);}


	for (i=0;i<NLines;i++) {
		j = (i+1)%NLines; 
		area+= 0.5*(ze[i].real()*ze[j].imag()-ze[i].imag()*ze[j].real());
		if (i!=NLines-1){perim+=abs(ze[i]-ze[i+1]);}
		upperswap(maxside,abs(ze[i]-ze[i+1]));
	}

	perim+=abs(ze[0]-ze[NLines-1]);
	upperswap(maxside,abs(ze[0]-ze[NLines-1]));

	//cout <<"Area: "<<area<<" "<<perim<<endl;
	// empirical spacing formula
	spacing=1.7/(pow(sqrt((double)ncontrol)+1.0,2.0))*perim*area/pow(maxside,2.0);
	spacing=5;//TMP DEBUG

	for (i=0; i<NLines; i++){
		points  [i].z=ze[i];
		points  [i].F=spacing;
		points  [i].status =BOUNDARY; 
		segments[i].p0=i;
		segments[i].p1=i+1;
		segments[i].status =BOUNDARY;
		segments[i].radius =0;
	}
	segments[NLines-1].p1=0;

	spcstruct spacing_package;
	spacing_package.type = INTERP;

	pMesh=new CTriMesh("temp",points,NLines,segments,NLines,80*ncontrol,spacing_package);

	if (pMesh!=NULL){
		pMesh->DynamicBuildMesh(true,true);
		pMesh->GetInternalNodes(zpoints,ncontrol);
		worked=true;

		ofstream TMP;
		TMP.open("controlpoints.bna");
		for (i=0; i<ncontrol;i++){
			TMP << "\" ctrl \",  1" <<endl;
			TMP << zpoints[i].real()<< " , "<<zpoints[i].imag()<<endl;
		}
		TMP.close();
	}
	else{ExitGracefully("CreatePolygonControlPoints:unable to generate",RUNTIME_ERR);}

	pMesh->WriteGeometry(false,true);//TMP DEBUG

	delete [] points;
	delete [] segments;
	delete pMesh;

	return worked;
}
//****************************************************************************
void CTriMesh::GetInternalNodes(cmplex *&zpoints,int &ncontrol) const{
	bool *nodesaccounted=new bool [NumNodes];
	zpoints=new cmplex [NumNodes-3];
	ncontrol=0;
	for (int n=0; n<NumNodes-3;n++){
		nodesaccounted[n]=false;
	}
	nodesaccounted[0]=true;
	nodesaccounted[1]=true;
	nodesaccounted[2]=true;
	for(int s=0; s<NumSides; s++){
		if(side[s].status==INTERIOR){
			if ((!nodesaccounted[side[s].n1]) && (node[side[s].n1].status==INTERIOR)){
				zpoints[ncontrol]=node[side[s].n1].z;
				nodesaccounted[side[s].n1]=true;
				ncontrol++;
			} 
			if ((!nodesaccounted[side[s].n2]) && (node[side[s].n2].status==INTERIOR)){
				zpoints[ncontrol]=node[side[s].n2].z;
				nodesaccounted[side[s].n2]=true;
				ncontrol++;
			} 
		}
	}
	delete [] nodesaccounted;
}
/*************************************************************************
													Renumerate
**************************************************************************
	cleans and renumerates all nodes, elements and sides
	used after all points have been inserted
	banded matrix now possible for FE simulation
	called once from constructor
-------------------------------------------------------------------------*/
void CTriMesh::Renumerate(){

	int n, s, e, c, d, i, j, k;
	int new_elem(0), new_node(0), new_side(0);
	int next_e, next_s, lowest;

	
	int *newnodeindex=new int [NumNodes]; 						
	for(n=0; n<NumNodes; n++) {newnodeindex[n]=OFF;}
	int *newelemindex=new int [NumElems];
	for(e=0; e<NumElems; e++) {newelemindex[e]=OFF;}
	int *newsideindex=new int [NumSides];
	for(s=0; s<NumSides; s++) {newsideindex[s]=OFF;}
	
	//Searching for the first element which is on 
	e=0; while (!elem[e].on){e++;}
	
	//Assigning numbers 0 and 1 to the nodes of first element  
	newnodeindex[elem[e].ni]=new_node; new_node++;
	newnodeindex[elem[e].nj]=new_node; new_node++;

	//int i1,j1,k1;
	//---Renumeration of nodes-----------------------------------------
	//Problems with doublenode boundaries
	do{
	 
		lowest = NumNodes+NumNodes;
		next_e = OFF;
		//int tmp;
		for(e=0; e<NumElems; e++){
			if((elem[e].on) && (newelemindex[e]==OFF)){
				
				i=newnodeindex[elem[e].ni];
				j=newnodeindex[elem[e].nj];
				k=newnodeindex[elem[e].nk];

				if( i+j+k+2 == abs(i) + abs(j) + abs(k) ){ //if only one is off, 
					if( (i==OFF) && (j+k) < lowest) {next_e=e; lowest=j+k;}
					if( (j==OFF) && (i+k) < lowest) {next_e=e; lowest=i+k;}
					if( (k==OFF) && (i+j) < lowest) {next_e=e; lowest=i+j;}
				}
				
			}
		}

		int dni,dnj,dnk;
		if(next_e!=OFF){
			i=newnodeindex[elem[next_e].ni];
			j=newnodeindex[elem[next_e].nj];
			k=newnodeindex[elem[next_e].nk];

			dni=OFF;dnj=OFF;dnk=OFF;
			for (int k2=0; k2<NumDoublenodes;k2++){
				if ((i==OFF) && (dnodes[k2].n1==elem[next_e].ni)){dni=dnodes[k2].n2;}
				if ((j==OFF) && (dnodes[k2].n1==elem[next_e].nj)){dnj=dnodes[k2].n2;}
				if ((k==OFF) && (dnodes[k2].n1==elem[next_e].nk)){dnk=dnodes[k2].n2;}
			}

			//Assign a new number to the node  
			if(i ==OFF) {
				newnodeindex[elem[next_e].ni] = new_node; new_node++;
				if ((dni!=OFF) && (newnodeindex[dni]==OFF)) {newnodeindex[dni] = new_node; new_node++;}
			}
			if(j ==OFF) {
				newnodeindex[elem[next_e].nj] = new_node; new_node++;
				if ((dnj!=OFF) && (newnodeindex[dnj]==OFF)) {newnodeindex[dnj] = new_node; new_node++;}
			}
			if(k ==OFF) {
				newnodeindex[elem[next_e].nk] = new_node; new_node++;
				if ((dnk!=OFF) && (newnodeindex[dnk]==OFF)) {newnodeindex[dnk] = new_node; new_node++;}
			}

		}
		if (ProgramAborted()){return;}

	} while(next_e != OFF);


	//---Renumeration of elements-------------------------------------
	do {
		lowest = NumNodes+NumNodes+NumNodes;
		next_e = OFF;

		for(e=0; e<NumElems; e++){
			if(elem[e].on && newelemindex[e]==OFF){ 
				i=newnodeindex[elem[e].ni];
				j=newnodeindex[elem[e].nj];
				k=newnodeindex[elem[e].nk];
				
				if( (i+j+k) < lowest ){next_e=e; lowest=i+j+k;}
			}
		}

		if(next_e!=OFF){
			newelemindex[next_e]=new_elem; new_elem++;
		}
	} while(next_e != OFF);

	//cout << "new_elem:" <<new_elem<<endl;

	//cout <<"Renumerate sides..."<<endl;


	//---Renumeration of sides----------------------------------------
	do {
		lowest = NumNodes+NumNodes;
		next_s = OFF;

		for(s=0; s<NumSides; s++){
			if(side[s].status!=IS_OFF && newsideindex[s]==OFF){
				c=newnodeindex[side[s].n1];
				d=newnodeindex[side[s].n2];
				if( (c+d) < lowest){next_s=s; lowest=c+d;}
			}
		}
		if(next_s!=OFF){
			newsideindex[next_s]=new_side;	new_side++;
		}  
	} while(next_s != OFF);
	//cout << "new_side:" <<new_side<<endl;

	//---Update doublenode array--------------------------------------
	for (n=0; n<NumDoublenodes; n++){
		dnodes[n].n1=newnodeindex[dnodes[n].n1];
		dnodes[n].n2=newnodeindex[dnodes[n].n2];
	}

	//Totally recreate element, side, and node arrays with new indices
	//----------------------------------------------------------------

	nodestruct * tmpnode=new nodestruct[new_node];
	elemstruct * tmpelem=new elemstruct[new_elem];
	sidestruct * tmpside=new sidestruct[new_side];

	//copy nodes in order
	int nextnode=0;
	do {
		for (n=0; n<NumNodes; n++){
			if (newnodeindex[n]==nextnode){
				tmpnode[nextnode].z 		=node[n].z;
				tmpnode[nextnode].F 		=node[n].F;
				tmpnode[nextnode].status=node[n].status;
				tmpnode[nextnode].inv 	=node[n].inv;
				nextnode++;
			}
		}
	} while (nextnode<new_node);

	//copy elements in order
	int nextelem=0;
	do {
		for (e=0; e<NumElems; e++){
			if (newelemindex[e]==nextelem){
				if (elem[e].ei==OFF){tmpelem[nextelem].ei=OFF;}
				else								{tmpelem[nextelem].ei=newelemindex[elem[e].ei];}
				if (elem[e].ej==OFF){tmpelem[nextelem].ej=OFF;}
				else								{tmpelem[nextelem].ej=newelemindex[elem[e].ej];}
				if (elem[e].ek==OFF){tmpelem[nextelem].ek=OFF;}
				else								{tmpelem[nextelem].ek=newelemindex[elem[e].ek];}
				tmpelem[nextelem].si=newsideindex[elem[e].si];		
				tmpelem[nextelem].sj=newsideindex[elem[e].sj];
				tmpelem[nextelem].sk=newsideindex[elem[e].sk];
				tmpelem[nextelem].ni=newnodeindex[elem[e].ni];
				tmpelem[nextelem].nj=newnodeindex[elem[e].nj];
				tmpelem[nextelem].nk=newnodeindex[elem[e].nk];
				tmpelem[nextelem].on=elem[e].on;	//should always be true!!!
				tmpelem[nextelem].R =elem[e].R;
				tmpelem[nextelem].r =elem[e].r; 
				tmpelem[nextelem].zv=elem[e].zv;
				tmpelem[nextelem].zin  =elem[e].zin;	
				tmpelem[nextelem].state=elem[e].state;
				nextelem++;
			}
		}
	} while (nextelem<new_elem);

	//copy sides in order
	int nextside=0;
	do {
		for (s=0; s<NumSides; s++){
			if (newsideindex[s]==nextside){
				if (side[s].ea==OFF) {tmpside[nextside].ea=OFF;}
				else								 {tmpside[nextside].ea=newelemindex[side[s].ea];}
				if (side[s].eb==OFF) {tmpside[nextside].eb=OFF;}
				else								 {tmpside[nextside].eb=newelemindex[side[s].eb];}
				tmpside[nextside].n1=newnodeindex[side[s].n1];
				tmpside[nextside].n2=newnodeindex[side[s].n2];
				if (side[s].na==OFF) {tmpside[nextside].na=OFF;}
				else								 {tmpside[nextside].na=newnodeindex[side[s].na];}
				if (side[s].nb==OFF) {tmpside[nextside].nb=OFF;}
				else								 {tmpside[nextside].nb=newnodeindex[side[s].nb];}
				tmpside[nextside].status =side[s].status;
				tmpside[nextside].radius =side[s].radius;
				nextside++;
			}
		}
		//cout << nextside<<" "<<new_side-1<<endl;
	} while (nextside<new_side);

	delete [] node;
	delete [] elem;
	delete [] side;

	node=tmpnode;
	elem=tmpelem;
	side=tmpside;
	
	NumNodes=new_node;
	NumElems=new_elem;
	NumSides=new_side;

	delete [] newnodeindex; 						
	delete [] newelemindex;
	delete [] newsideindex;
}
/**************************************************************************
	CalculateSpacingFunction:
***************************************************************************
		input: element index (e)
					 node index (n)
	This function calculates the value of the spacing function in  
	a new node 'n' which is inserted in element 'e' by a linear 	 
	approximation from the values of the spacing function in the	 
	elements nodes.  
	->right now very affected by ordering of node placement (only with smoothing)
-------------------------------------------------------------------------*/
double CTriMesh::CalculateSpacingFunction(const int e, const cmplex &z) const{

	cmplex dzji,dzki,dz;
	double det, a, b;

	if ((e<0) || (e>NumElems)){cout << "Spacing function: Bad element"<<endl;return 0.0;}
	
	dzji = node[elem[e].nj].z - node[elem[e].ni].z;
	dzki = node[elem[e].nk].z - node[elem[e].ni].z;
	dz	 =									z - node[elem[e].ni].z;

	det = (dzji.real()*dzki.imag()-dzki.real()*dzji.imag());

	a = ( dzki.imag()*dz.real() - dzki.real()*dz.imag())/det;
	b = (-dzji.imag()*dz.real() + dzji.real()*dz.imag())/det;

	return			node[elem[e].ni].F + 
							a*(node[elem[e].nj].F - node[elem[e].ni].F) +
							b*(node[elem[e].nk].F - node[elem[e].ni].F);

			//if (spacing_package.type=COURANT){
			//}
}
/**************************************************************************
	CalculateCircles:
***************************************************************************
		This function calculates radii of inscribed and
		circumscribed circles for a given element (e)
--------------------------------------------------------------------------*/
void CTriMesh::CalculateCircles(int e){ 

	double x, y;
	double xi, yi, xj, yj, xk, yk; //node locations (x,y)
	double xij, yij, xjk, yjk;		 //side midpoints (x,y)
	double num, den;
	double Li, Lj, Lk;					//length of adjacent sides
	double Perimeter; 					//perimiter of element

	if ((e<0) || (e>NumElems-1)){cout << "CalculateCircles: bad e"<<endl;return;}

	xi=node[elem[e].ni].z.real(); yi=node[elem[e].ni].z.imag();
	xj=node[elem[e].nj].z.real(); yj=node[elem[e].nj].z.imag();
	xk=node[elem[e].nk].z.real(); yk=node[elem[e].nk].z.imag();
	 
	xij=0.5*(xi+xj); yij=0.5*(yi+yj);
	xjk=0.5*(xj+xk); yjk=0.5*(yj+yk);

	//circumscribed circle
	num = (xij-xjk)*(xj-xi) + (yij-yjk)*(yj-yi);
	den = (xj -xi) *(yk-yj) - (xk -xj) *(yj-yi);

	if(den>0){
		x=xjk + num/den*(yk-yj);
		y=yjk - num/den*(xk-xj);
		elem[e].zv = cmplex(x,y); 
		elem[e].R  = sqrt( (xi-x)*(xi-x) + (yi-y)*(yi-y) );
	}

	//inscribed circle
	Li=abs(node[side[elem[e].si].n1].z-node[side[elem[e].si].n2].z);
	Lj=abs(node[side[elem[e].sj].n1].z-node[side[elem[e].sj].n2].z);
	Lk=abs(node[side[elem[e].sk].n1].z-node[side[elem[e].sk].n2].z);

	Perimeter =Li+Lj+Lk;

	elem[e].zin =cmplex( ( xi*Li + xj*Lj + xk*Lk ),( yi*Li + yj*Lj + yk*Lk )) / Perimeter;

	elem[e].r 	= xi*(yj-yk) - xj*(yi-yk) + xk*(yi-yj) / Perimeter;
}


/**************************************************************************
	SwapDiagonal
***************************************************************************
		input: side index
		swaps diagonal of two adjacent elements from long diagonal of 4-sided elem to short diagonal
		rearranges adjacent elements appropriately
-------------------------------------------------------------------------*/
void CTriMesh::SwapDiagonal(const int s){ 
	int 	 na, nb;
	int 	 n1, n2;
	int 	 ea,	eb;
	int 	 eac, ead, ebc, ebd;
	int 	 sad, sac, sbc, sbd;
 
/*-------------------------------------------------------- 
				 |\ 														 
				 |	\ 													 
				 |		\ 												
				 | eac	\ 											 
			 na|________\n2_________		 ___sac__
				/|				/|				/ 		| 			/|				
			/  | ea 	/  | ebd	/ 			| 		/  |	 
		/ 	 |		/ 	 |		/ 		 sad| 	s 	 |sbd 	 
	/  ead |	/  eb  |	/ 					|  /		 |	
/________|/________|/ 						|/_______|
				n1\ 			 |nb								sbc 	
						\  ebc |									 
							\ 	 |												 
								\  |														
									\|															
---------------------------------------------------------*/

	if ((s<0) || (s>NumSides)){cout<< "SwapDiagonal: Bad side specified"<<endl;}

	ea=side[s].ea; 
	eb=side[s].eb;

	if (ea==OFF){cout << "SwapDiagonal::Bad Side (ea)"<<endl;}
	if (eb==OFF){cout << "SwapDiagonal::Bad Side (eb)"<<endl;}

	//--gather information--------------------------------------------
	na=side[s].na; nb=side[s].nb; 
	n1=side[s].n1; n2=side[s].n2;

	if(elem[ea].ei==eb) {ead=elem[ea].ej; eac=elem[ea].ek; sad=elem[ea].sj; sac=elem[ea].sk;}
	if(elem[ea].ej==eb) {ead=elem[ea].ek; eac=elem[ea].ei; sad=elem[ea].sk; sac=elem[ea].si;} 	
	if(elem[ea].ek==eb) {ead=elem[ea].ei; eac=elem[ea].ej; sad=elem[ea].si; sac=elem[ea].sj;}

	if(elem[eb].ei==ea) {ebc=elem[eb].ej; ebd=elem[eb].ek; sbc=elem[eb].sj; sbd=elem[eb].sk;}
	if(elem[eb].ej==ea) {ebc=elem[eb].ek; ebd=elem[eb].ei; sbc=elem[eb].sk; sbd=elem[eb].si;}
	if(elem[eb].ek==ea) {ebc=elem[eb].ei; ebd=elem[eb].ej; sbc=elem[eb].si; sbd=elem[eb].sj;}

	//--modify elements ea and eb-------------------------------------
	elem[ea].ni =na;	elem[ea].nj =nb;	elem[ea].nk =n2;	 
	elem[ea].ei=ebd; elem[ea].ej=ead; elem[ea].ek=eb;  
	elem[ea].si=sbd; elem[ea].sj=sad; elem[ea].sk=s;	
	
	elem[eb].ni =na;	elem[eb].nj =n1;	elem[eb].nk =nb;	
	elem[eb].ei=ebc; elem[eb].ej=ea;	elem[eb].ek=eac;	
	elem[eb].si=sbc; elem[eb].sj=s; 	elem[eb].sk=sac;	

	if(eac!=OFF){
		if(elem[eac].ei==ea){elem[eac].ei=eb;}
		if(elem[eac].ej==ea){elem[eac].ej=eb;}
		if(elem[eac].ek==ea){elem[eac].ek=eb;}
	}
	if(ebd!=OFF){
		if(elem[ebd].ei==eb){elem[ebd].ei=ea;}
		if(elem[ebd].ej==eb){elem[ebd].ej=ea;}
		if(elem[ebd].ek==eb){elem[ebd].ek=ea;} 
	}
	
	//modify all sides-----------------------------------------------
	if(side[sad].ea==ea) {side[sad].na=nb;}
	if(side[sad].eb==ea) {side[sad].nb=nb;}

	if(side[sbc].ea==eb) {side[sbc].na=na;}
	if(side[sbc].eb==eb) {side[sbc].nb=na;}

	if(side[sbd].ea==eb) {side[sbd].ea=ea; side[sbd].na=na;}
	if(side[sbd].eb==eb) {side[sbd].eb=ea; side[sbd].nb=na;}
 
	if(na<nb){
		side[s].n1=na;	side[s].na=n2;
		side[s].n2=nb;	side[s].nb=n1;
		side[s].ea=ea;	side[s].eb=eb;
	}
	else {
		side[s].n1=nb;	side[s].na=n1;	 
		side[s].n2=na;	side[s].nb=n2;
		side[s].ea=eb;	side[s].eb=ea;
	}

	if(side[sac].ea==ea) {side[sac].ea=eb; side[sac].na=nb;}
	if(side[sac].eb==ea) {side[sac].eb=eb; side[sac].nb=nb;}

	if(side[sad].ea==ea) {side[sad].na=nb;}
	if(side[sad].eb==ea) {side[sad].nb=nb;}

	if(side[sbc].ea==eb) {side[sbc].na=na;}
	if(side[sbc].eb==eb) {side[sbc].nb=na;}

	if(side[sbd].ea==eb) {side[sbd].ea=ea; side[sbd].na=na;}
	if(side[sbd].eb==eb) {side[sbd].eb=ea; side[sbd].nb=na;}

	//recalculate circles of elements--------------------------------
	CalculateCircles(ea);
	CalculateCircles(eb);
}
/**************************************************************************
	InsertNode
***************************************************************************
	input: node location (x,y) 
				 boundary condition info (status) 	 
	This function inserts a new node into the mesh	
	returns the element in which the node was located (or off, if point was outside any element)
	Calls To :
		CalculateSpacingFunction(e, n)
		CalculateCircles(e)
		SwapDiagonal(s)
-------------------------------------------------------------------------*/
int CTriMesh::InsertNode(const cmplex 		 z, 
												 const status_type status,
												 const cmplex &invector){

	int 	 n,ni,nj,nk;
	int 	 e,ei,ej,ek;
	int 	 s,si,sj,sk;
	bool	 swap;


	NumNodes++; 												 // one new node 
	n=NumNodes-1; 											 // current node=n

	if (NumNodes>MaxNodes){ExitGracefully("CTriMesh::InsertNode: Exceeded max nodes",TOO_MANY);} 

	node[n].z=z;
	node[n].status=status;
	if (invector!=0.0){
		node[n].inv=invector/abs(invector);
		//cout << "IN VECTOR: "<< node[n].inv<<endl; 
	}
	else{node[n].inv=0.0;}


	// find the existing element which contains new node 
	e = GetElement(z);

	if ((e<0) || (e>NumElems-1)){cout << "CTriMesh::InsertNode: insertion point outside domain"<<endl; return OFF;}

	// interpolate the spacing function in the new node 

	node[n].F=CalculateSpacingFunction(e, z);

	ni=elem[e].ni; nj=elem[e].nj; nk =elem[e].nk;
	ei=elem[e].ei; ej=elem[e].ej; ek=elem[e].ek; 
	si=elem[e].si; sj=elem[e].sj; sk=elem[e].sk; 
 
	NumElems+=2;
	NumSides+=3;

	//--	new elements	-----------------------------------------

	elem[NumElems-2].ni=n;	elem[NumElems-2].nj=nk; 				elem[NumElems-2].nk=ni;
	elem[NumElems-1].ni=n;	elem[NumElems-1].nj=ni; 				elem[NumElems-1].nk=nj; 
 
	elem[NumElems-2].ei=ej; elem[NumElems-2].ej=NumElems-1; elem[NumElems-2].ek=e;
	elem[NumElems-1].ei=ek; elem[NumElems-1].ej=e;					elem[NumElems-1].ek=NumElems-2;
 
	elem[NumElems-2].si=sj; elem[NumElems-2].sj=NumSides-2; elem[NumElems-2].sk=NumSides-3;
	elem[NumElems-1].si=sk; elem[NumElems-1].sj=NumSides-1; elem[NumElems-1].sk=NumSides-2;

	elem[NumElems-1].on=true; 				
	elem[NumElems-2].on=true;
	
	//-- new sides ----------------------------------------------  
	
	side[NumSides-3].n1=nk; 				side[NumSides-3].n2 =n; 						// c-d 
	side[NumSides-3].na=nj; 				side[NumSides-3].nb =ni;						// a-b 
	side[NumSides-3].ea=e;					side[NumSides-3].eb =NumElems-2;
	side[NumSides-3].status=INTERIOR; 
	side[NumSides-3].radius=0.0;
 
	side[NumSides-2].n1=ni; 				side[NumSides-2].n2 =n; 						// c-d 
	side[NumSides-2].na=nk; 				side[NumSides-2].nb =nj;						// a-b 
	side[NumSides-2].ea=NumElems-2; side[NumSides-2].eb =NumElems-1;
	side[NumSides-2].status=INTERIOR; 
	side[NumSides-2].radius=0.0;
	
	side[NumSides-1].n1 =nj;				side[NumSides-1].n2 =n; 						// c-d 
	side[NumSides-1].na =ni;				side[NumSides-1].nb =nk;						// a-b 
	side[NumSides-1].ea=NumElems-1; side[NumSides-1].eb =e; 			
	side[NumSides-1].status=INTERIOR; 
	side[NumSides-1].radius=0.0;

	elem[e].ni = n;
	elem[e].ej = NumElems-2; elem[e].ek = NumElems-1;
	elem[e].sj = NumSides-3; elem[e].sk = NumSides-1;

	if(side[si].na==ni) {side[si].na=n; side[si].ea=e;}
	if(side[si].nb==ni) {side[si].nb=n; side[si].eb=e;}
 
	if(side[sj].na==nj) {side[sj].na=n; side[sj].ea=NumElems-2;}
	if(side[sj].nb==nj) {side[sj].nb=n; side[sj].eb=NumElems-2;}
 
	if(side[sk].na==nk) {side[sk].na=n; side[sk].ea=NumElems-1;} 
	if(side[sk].nb==nk) {side[sk].nb=n; side[sk].eb=NumElems-1;} 

	//if ((ej<-1) || (ej>NumElems-1) || (ek<-1)|| (ek>NumElems-1)){cout << "CTriMesh::InsertNode: Bad Elem"<<endl;}

	if(ej!=OFF){
		if(elem[ej].ei==e) {elem[ej].ei=NumElems-2;}
		if(elem[ej].ej==e) {elem[ej].ej=NumElems-2;}
		if(elem[ej].ek==e) {elem[ej].ek=NumElems-2;}
	}
	if(ek!=OFF){
		if(elem[ek].ei==e) {elem[ek].ei=NumElems-1;}
		if(elem[ek].ej==e) {elem[ek].ej=NumElems-1;}
		if(elem[ek].ek==e) {elem[ek].ek=NumElems-1;}
	}

	// Find circumcenters for two new elements, and for the one who's segment has changed
	CalculateCircles(e);
	CalculateCircles(NumElems-2);
	CalculateCircles(NumElems-1);
	
	//--Bowyer Algorithm--------------------------------------------------------
	//sifts through all sides associated with newly inserted node n
	//swap diagonals if node is in adjacent element's circle
	//iterates until all swaps have been performed
	//swaps are not performed on non-interior segments

	int count(0);
	do{ 
		swap=false;
		for(s=0; s<NumSides; s++){
			if(side[s].status==INTERIOR){
				if		 (side[s].na==n){
					e=side[s].eb; 
					if(e!=OFF){
						if( abs(elem[e].zv-node[n].z) < elem[e].R ){
							SwapDiagonal(s); 
							swap=true;
						}
					}
				}
				else if(side[s].nb==n){
					e=side[s].ea; 
					if(e!=OFF){
						if( abs(elem[e].zv-node[n].z) < elem[e].R ){
							SwapDiagonal(s); 
							swap=true;
						}
					}
				}
			}
		}
		count++;
	} while((swap) && (count<1000));
	
	//check all side lengths to make sure min spacing is not too small
	for(s=0; s<NumSides; s++){
		if (side[s].status!=IS_OFF){
			if (abs(node[side[s].n1].z-node[side[s].n2].z)<minSpacing){
				cout << "Bad side length upon insertion of new node"<<endl;
				e=OFF; //this indicates bad insertion
			}
		}
	}

	return e;
}
/**************************************************************************
	UpdateNewSideStatus
***************************************************************************
 Identify sides associated with newly placed node (node NumNodes-1)
 if side is between prev_n and the new n or between next_n and new n, 
 then the appropriate status is assigned
 only important to sides composed of "specified" nodes (i.e. not interior nodes)
-------------------------------------------------------------------------*/ 
void CTriMesh::UpdateNewSideStatus(const int				 adjacent_n, 
																	 const status_type adjacent_status,
																	 const double      radius){

	for(int s=3; s<NumSides; s++)
	{
		if      (side[s].n1==adjacent_n && side[s].n2==NumNodes-1){
			side[s].status= adjacent_status;
			side[s].radius= radius;
		} 
		else if (side[s].n2==adjacent_n && side[s].n1==NumNodes-1){
			side[s].status= adjacent_status;
			side[s].radius=-radius;
		}
	}

}
/**************************************************************************
	GetElement
***************************************************************************
 Identifies element which contains the point z, returns element index (or OFF if no element)
 optimized by assuming that the element of the point last queried is close 
 to the element of the newly requested point
 -Especially robust after renumeration
 -may improve by testing "is inside" with half planes
-------------------------------------------------------------------------*/ 
int CTriMesh::GetElement(const cmplex &z) const{
  static double xii,xij,xik;
	static int laste(0);
	int 			 e;
	//double	 mult;

	if (laste>=NumElems){laste=0;}//updates static member

	for (int i=0; i<NumElems; i++){

		if			(i==0)		{e=laste;}
		else if ((i%2)==1){e=laste-(i-1)/2-1;}
		else							{e=laste+(i  )/2	;}

		if			(e<0				 ){e+=NumElems;}
		else if (e>NumElems-1){e-=NumElems;}
		
		//Area test
		if (TriArea(z, node[elem[e].ni].z, node[elem[e].nj].z) >= 0.0 && 
				TriArea(z, node[elem[e].nj].z, node[elem[e].nk].z) >= 0.0 &&
				TriArea(z, node[elem[e].nk].z, node[elem[e].ni].z) >= 0.0 ) {
			laste=e;
			return e;
		}
		//Barycentric coordinate Test - slower
		/*TriGlobalToLocal(node[elem[e].ni].z, node[elem[e].nj].z, node[elem[e].nk].z,z,xii,xij,xik);
		if ((xii>=0.0) && (xij>=0) && (xik>=0)){
			laste=e;
			return e;
		}*/
	}
	return OFF;
} 		

//--------------------------------------------------------------------------																			 
bool	 CTriMesh::IsInside(const pt3D &pt) const{
	//TMP DEBUG MUST OPTIMIZE: should be a simple in polygon check
	if (GetElement(c3Dto2D(pt))==OFF){return false;}
	else														 {return true; }
} 	
//*************************************************************************
CTriMesh *CTriMesh::Parse(ifstream &input, int &l,char *Name){
	//TriMesh 
	//{x y} x (NumNodes) 
	//&
	CTriMesh*pMesh=NULL;

	bool		 eof(false),done(false);
	int 		 Len,npts(0);
	cmplex	 pts[MAX_VORONOI_PTS];
	char		*s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Voronoi Mesh"<<endl;}			 eof=TokenizeLine(input,s,Len); l++;	 
	do {
		if (npts>=MAX_VORONOI_PTS) { ExitGracefully("CTriMesh::Parse- too many points specified",TOO_MANY);}
		if (Len==2) {
			pts[npts]=s_to_c(s[0],s[1]); 
			npts++; 																					 eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&"))) {
			pMesh = new CTriMesh(Name,pts,npts);
			done=true;
		}
		else if(Len==0)  {																	 eof=TokenizeLine(input,s,Len); l++;}
		else	{ImproperFormat(s,l); return NULL;}
	} while ((!done) && (!eof));

	if (eof) {return NULL;}
	else		 {return pMesh;}
}
//*********************************************************************
void CTriMesh::WriteGeometry(bool DXFon, bool BNAon) const{

	int 	 s;
	double xc, yc, xd, yd;
	cmplex za, zb;
	ofstream DEL;
	ofstream VOR;
	ofstream DXF;

	if (BNAon)
	{
		DEL.open("delunay.bna");
		if (DEL.fail()){		
			cout << "delunay.bna cannot be opened for output"<<endl;
			return;
		}
		VOR.open("voronoi.bna");
		if (VOR.fail()){		
			cout << "voronoi.bna cannot be opened for output"<<endl;
			return;
		}
	}
	if (DXFon)
	{
		DXF.open("mesh.dxf");
		if (DXF.fail()){		
			cout << "mesh.dxf cannot be opened for output"<<endl;
			return;
		}
	}

	//DXF Header
	if (DXFon)
	{
		DXF<<"0"<<endl<<"SECTION" <<endl;
		DXF<<"2"<<endl<<"ENTITIES"<<endl;
	}

	//Draw boundary  
	if (DXFon)
	{
		for(s=0; s<NumSides; s++){
			if(side[s].status==BOUNDARY) 
			{
			 xc=node[side[s].n1].z.real(); yc=node[side[s].n1].z.imag();
			 xd=node[side[s].n2].z.real(); yd=node[side[s].n2].z.imag();
			 line_dxf(DXF,xc, yc, 0, xd, yd, 0, "boundary",DXF_WHITE);
			}
			else if(side[s].status==DOUBLENODE) 
			{
			 xc=node[side[s].n1].z.real(); yc=node[side[s].n1].z.imag();
			 xd=node[side[s].n2].z.real(); yd=node[side[s].n2].z.imag();
			 line_dxf(DXF,xc, yc, 0, xd, yd, 0, "doublenode",DXF_ORANGE);
			}
		}
	}
	//Draw Delaunay  
	if (BNAon)
	{
		for(s=0; s<NumSides; s++){
			if(side[s].status!=IS_OFF) 
			{
			xc=node[side[s].n1].z.real(); yc=node[side[s].n1].z.imag();
			xd=node[side[s].n2].z.real(); yd=node[side[s].n2].z.imag();
			DEL << "\" delunay \",  "<<-2<<endl<<xc<< " , " <<yc<<endl<<xd<< " , " <<yd<<endl;
			}
		}
	}
	if (DXFon)
	{
		for(s=0; s<NumSides; s++){
			if(side[s].status!=IS_OFF) 
			{
			 xc=node[side[s].n1].z.real(); yc=node[side[s].n1].z.imag();
			 xd=node[side[s].n2].z.real(); yd=node[side[s].n2].z.imag();
			 line_dxf(DXF,xc, yc, 0, xd, yd, 0, "delaunay",DXF_BLUE);
			}
		}
		//Draw Ugly element
		if (ugly!=OFF){
			xc=node[elem[ugly].ni].z.real();yc=node[elem[ugly].ni].z.imag();
			xd=node[elem[ugly].nj].z.real();yd=node[elem[ugly].nj].z.imag();
			line_dxf(DXF,xc, yc, 0, xd, yd, 0, "ugly",DXF_GREEN);
			xc=node[elem[ugly].nk].z.real();yc=node[elem[ugly].nk].z.imag();
			line_dxf(DXF,xc, yc, 0, xd, yd, 0, "ugly",DXF_GREEN);
			xd=node[elem[ugly].ni].z.real();yd=node[elem[ugly].ni].z.imag();
			line_dxf(DXF,xc, yc, 0, xd, yd, 0, "ugly",DXF_GREEN);
		}
	}
	//Draw Voronoi	
	if ((DXFon) || (BNAon))
	{
		int ea, eb;
		for(s=0; s<NumSides; s++){
			if(side[s].status!=IS_OFF){
				ea=side[s].ea;
				eb=side[s].eb;
				if(ea!=OFF){za=elem[ea].zv;}
				else			 {za=0.5*(node[side[s].n1].z+node[side[s].n2].z);} //for non-obtuse elements
 
				if(eb!=OFF){zb=elem[eb].zv;}
				else			 {zb=0.5*(node[side[s].n1].z+node[side[s].n2].z);} //for non-obtuse elements

				if (BNAon){
					VOR << "\" voronoi \",  "<<-2<<endl;
					VOR <<za.real()<< " , " <<za.imag()<<endl;
					VOR <<zb.real()<< " , " <<zb.imag()<<endl;
				}
				if (DXFon){line_dxf(DXF,za.real(), za.imag(), 0, zb.real(), zb.imag(), 0, "voronoi",DXF_VIOLET);}
			}
		}
	}
	//Draw Delunay Centers
	/*for (int e=0; e<NumElems; e++){
		if (elem[e].on){
			xc=elem[e].cen.real();
			yc=elem[e].cen.imag();
			DEL << "\" center \",  "<<1<<endl<<xc<< " , " <<yc<<endl;
		}
	}*/

	if (DXFon)
	{
		DXF<<"0"		 <<endl<<"ENDSEC"<<endl;
		DXF<<"0"		 <<endl<<"EOF"	 <<endl;
		DXF.close();
	}
	if (BNAon){
		DEL.close();
		VOR.close();
	}
}

/**************************************************************************
	Translates Boundary marker from input file to enumerated status
-------------------------------------------------------------------------*/
status_type CTriMesh::TranslateMark(const int mark){
	if			(mark==1){return BOUNDARY;}
	else if (mark==2){return DOUBLENODE;}//double node boundary
	else if (mark<=0){return INTERIOR;}
	else						 {return IS_OFF;} 
}
/**************************************************************************
	ReadItself (OLD)
***************************************************************************
	Reads from Bluebird-created generated mesh files (.n,.e,.s,.d) created in CTriMesh::WriteOutput
-------------------------------------------------------------------------*/
bool CTriMesh::ReadItself(ifstream &NODE, ifstream &ELEM, ifstream &SIDE, ifstream &DNOD){
	
	bool		 eof(false),done(false);
	int 		 Len,npts(0);
	char		*s[MAXINPUTITEMS];
	int 		 e(0),n(0),i(0);
	
	cout <<"Reading Triangular Mesh..."<<endl;
	NumDoublenodes=0;
	
	eof=TokenizeLine(NODE,s,Len);
	if (eof){ExitGracefully("CTriMesh:Unable To Read Nodes file",BAD_DATA);}
	if (Len==1){NumNodes=s_to_i(s[0]);}
	
	eof=TokenizeLine(ELEM,s,Len);
	if (eof){ExitGracefully("CTriMesh:Unable To Read Elems file",BAD_DATA);}
	if (Len==1){NumElems=s_to_i(s[0]);}

	eof=TokenizeLine(SIDE,s,Len);
	if (eof){ExitGracefully("CTriMesh:Unable To Read Sides file",BAD_DATA);}
	if (Len==1){NumSides=s_to_i(s[0]);}

	eof=TokenizeLine(DNOD,s,Len);
	if (eof){NumDoublenodes=0;}
	else if (Len==1){NumDoublenodes=s_to_i(s[0]);}

	this->node	=new nodestruct [NumNodes];
	this->elem	=new elemstruct [NumElems]; 
	this->side	=new sidestruct [NumSides];
	this->dnodes=new dblestruct [NumDoublenodes];

	this->MaxNodes=NumNodes;
	
	while ((!eof) && (n<this->NumNodes)){
		eof=TokenizeLine(NODE,s,Len);
		if ((eof) || ((Len!=6) && (Len!=4))) {ExitGracefully("CTriMesh:Nodes File incorrectly formatted",BAD_DATA);}
		this->node[n].z =s_to_c(s[1],s[2]);
		this->node[n].status =TranslateMark(s_to_i(s[3]));
		this->node[n].F=0.0; //doesn't matter
		if (Len==6){
			this->node[n].inv = s_to_c(s[4],s[5]);
		}
		else{
			this->node[n].inv = 0.0;
		}
		n++;
	}

	while ((!eof) && (e<NumElems)){
		eof=TokenizeLine(ELEM,s,Len);
		if ((eof) || (Len!=13)) {ExitGracefully("CTriMesh:Elems File incorrectly formatted",BAD_DATA);}

		if (s_to_i(s[1])>=NumNodes){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])>=NumNodes){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[3])>=NumNodes){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[1])<0) 			 {ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])<0) 			 {ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[3])<0) 			 {ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		this->elem[e].ni=s_to_i(s[1]);
		this->elem[e].nj=s_to_i(s[2]);
		this->elem[e].nk=s_to_i(s[3]);

		if (s_to_i(s[4])>=NumElems){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad elem index)",BAD_DATA);}
		if (s_to_i(s[5])>=NumElems){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad elem index)",BAD_DATA);}
		if (s_to_i(s[6])>=NumElems){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad elem index)",BAD_DATA);}
		this->elem[e].ei=s_to_i(s[4]);
		this->elem[e].ej=s_to_i(s[5]);
		this->elem[e].ek=s_to_i(s[6]);

		if (s_to_i(s[7])>=NumSides){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad side index)",BAD_DATA);}
		if (s_to_i(s[8])>=NumSides){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad side index)",BAD_DATA);}
		if (s_to_i(s[9])>=NumSides){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad side index)",BAD_DATA);}
		this->elem[e].si=s_to_i(s[7]);
		this->elem[e].sj=s_to_i(s[8]);
		this->elem[e].sk=s_to_i(s[9]);

		//this->elem[e].zv =s_to_c(s[10],s[11]); //optional->use calculate circles
		
		this->elem[e].on=true;

		//Calculated in Process
		this->elem[e].area=TriArea(node[elem[e].ni].z,node[elem[e].nj].z,node[elem[e].nk].z);
		this->elem[e].centroid=0.0;
		
		//Calculated in calc circles
		this->elem[e].zin=0;
		this->elem[e].zv=0;
		this->elem[e].R =0;
		this->elem[e].r=0.0; 
		this->elem[e].state=DONE ;//Doesn't matter

		e++;
	}

	while ((!eof) && (i<NumSides)){
		eof=TokenizeLine(SIDE,s,Len);
		if ((eof) || (Len!=8)) {ExitGracefully("CTriMesh:Sides File incorrectly formatted",BAD_DATA);}

		if (s_to_i(s[1])>=NumNodes){ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])>=NumNodes){ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[3])>=NumNodes){ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[4])>=NumNodes){ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[1])<0) 			 {ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])<0) 			 {ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad node index)",BAD_DATA);}
		this->side[i].n1=s_to_i(s[1]);
		this->side[i].n2=s_to_i(s[2]);
		this->side[i].na=s_to_i(s[3]);
		this->side[i].nb=s_to_i(s[4]);

		if (s_to_i(s[5])>=NumElems){ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad elem index)",BAD_DATA);}
		if (s_to_i(s[6])>=NumElems){ExitGracefully("CTriMesh:Sides File incorrectly formatted (bad elem index)",BAD_DATA);}
		this->side[i].ea=s_to_i(s[5]);
		this->side[i].eb=s_to_i(s[6]);

		this->side[i].status =TranslateMark(s_to_i(s[7]));
		this->side[i].radius=0;

		i++;
	}

	n=0;
	while ((!eof) && (n<this->NumDoublenodes)){
		eof=TokenizeLine(DNOD,s,Len);
		if ((eof) || (Len!=4)) {ExitGracefully("CTriMesh:Doublenodes File incorrectly formatted",BAD_DATA);}
		this->dnodes[n].n1=s_to_i(s[1]);
		this->dnodes[n].n2=s_to_i(s[2]);
		this->dnodes[n].type=DNODE_JUMP;//TranslateDnode(s_to_i(s[2]));//TMP DEBUG
		n++;
	}

	for (e=0; e<NumElems;e++){
		CalculateCircles(e);
	}
	for (n=0; n<NumNodes;n++){
		upperswap(BBox.e,this->node[n].z.real());
		upperswap(BBox.n,this->node[n].z.imag());
		lowerswap(BBox.w,this->node[n].z.real());
		lowerswap(BBox.s,this->node[n].z.imag());
	}
	Process();
	
	cout <<"...Done Reading Triangular Mesh."<<endl;

	return true; 
}
/**************************************************************************
	ReadItself2 (NEW)
***************************************************************************
	Reads from Bluebird-created generated mesh files (.btm) created in CTriMesh::WriteOutput2
-------------------------------------------------------------------------*/
bool CTriMesh::ReadItself2(ifstream &BTM){
	
	bool		 eof(false),done(false);
	int 		 Len,npts(0);
	char		*s[MAXINPUTITEMS];
	int 		 e(0),n(0),i(0);
	
	cout <<"Reading Triangular Mesh..."<<endl;
	NumDoublenodes=0;
	
	eof=TokenizeLine(BTM,s,Len);
	if (eof){ExitGracefully("CTriMesh:Unable To Read Triangular mesh file(1)",BAD_DATA);}
	if (Len==1){NumNodes=s_to_i(s[0]);}
	else {ExitGracefully("CTriMesh:BTM File incorrectly formatted(1)",BAD_DATA);}
	
	eof=TokenizeLine(BTM,s,Len);
	if (eof){ExitGracefully("CTriMesh:Unable To Read Triangular mesh file(2)",BAD_DATA);}
	if (Len==1){NumElems=s_to_i(s[0]);}
	else {ExitGracefully("CTriMesh:BTM File incorrectly formatted(2)",BAD_DATA);}

	eof=TokenizeLine(BTM,s,Len);
	if (eof){ExitGracefully("CTriMesh:Unable To Read Triangular mesh file(3)",BAD_DATA);}
	if (Len==1){NumSides=s_to_i(s[0]);}
	else {ExitGracefully("CTriMesh:BTM File incorrectly formatted(3)",BAD_DATA);}

	eof=TokenizeLine(BTM,s,Len);
	if (eof){ExitGracefully("CTriMesh:Unable To Read Triangular mesh file(4)",BAD_DATA);}
	if (Len==1){NumDoublenodes=s_to_i(s[0]);}
	else {ExitGracefully("CTriMesh:BTM File incorrectly formatted(4)",BAD_DATA);}

	this->node	=new nodestruct [NumNodes];
	this->elem	=new elemstruct [NumElems]; 
	this->side	=new sidestruct [NumSides];
	this->dnodes=new dblestruct [NumDoublenodes];

	this->MaxNodes=NumNodes;
	
	while ((!eof) && (n<this->NumNodes)){
		eof=TokenizeLine(BTM,s,Len);
		if ((eof) || ((Len!=6) && (Len!=4))) {
			ExitGracefully("CTriMesh:Nodes File incorrectly formatted",BAD_DATA);}
		this->node[n].z =s_to_c(s[1],s[2]);
		this->node[n].status =TranslateMark(s_to_i(s[3]));
		this->node[n].F=0.0; //doesn't matter
		if (Len==6){
			this->node[n].inv = s_to_c(s[4],s[5]);
		}
		else{
			this->node[n].inv = 0.0;
		}
		n++;
	}

	while ((!eof) && (e<NumElems)){
		eof=TokenizeLine(BTM,s,Len);
		if ((eof) || (Len!=13)) {ExitGracefully("CTriMesh:Elems File incorrectly formatted",BAD_DATA);}

		if (s_to_i(s[1])>=NumNodes){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])>=NumNodes){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[3])>=NumNodes){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[1])<0) 			 {ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])<0) 			 {ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[3])<0) 			 {ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad node index)",BAD_DATA);}
		this->elem[e].ni=s_to_i(s[1]);
		this->elem[e].nj=s_to_i(s[2]);
		this->elem[e].nk=s_to_i(s[3]);

		if (s_to_i(s[4])>=NumElems){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad elem index)",BAD_DATA);}
		if (s_to_i(s[5])>=NumElems){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad elem index)",BAD_DATA);}
		if (s_to_i(s[6])>=NumElems){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad elem index)",BAD_DATA);}
		this->elem[e].ei=s_to_i(s[4]);
		this->elem[e].ej=s_to_i(s[5]);
		this->elem[e].ek=s_to_i(s[6]);

		if (s_to_i(s[7])>=NumSides){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad side index)",BAD_DATA);}
		if (s_to_i(s[8])>=NumSides){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad side index)",BAD_DATA);}
		if (s_to_i(s[9])>=NumSides){ExitGracefully("CTriMesh:Elems File incorrectly formatted (bad side index)",BAD_DATA);}
		this->elem[e].si=s_to_i(s[7]);
		this->elem[e].sj=s_to_i(s[8]);
		this->elem[e].sk=s_to_i(s[9]);

		//this->elem[e].zv =s_to_c(s[10],s[11]); //optional->use calculate circles
		
		this->elem[e].on=true;

		//Calculated in Process
		this->elem[e].area=TriArea(node[elem[e].ni].z,node[elem[e].nj].z,node[elem[e].nk].z);
		this->elem[e].centroid=0.0;
		
		//Calculated in calc circles
		this->elem[e].zin=0;
		this->elem[e].zv=0;
		this->elem[e].R =0;
		this->elem[e].r=0.0; 
		this->elem[e].state=DONE ;//Doesn't matter

		e++;
	}

	while ((!eof) && (i<NumSides)){
		eof=TokenizeLine(BTM,s,Len);
		if ((eof) || ((Len!=8) && (Len!=9))) {ExitGracefully("CTriMesh:Sides File incorrectly formatted",BAD_DATA);}

		if (s_to_i(s[1])>=NumNodes){ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])>=NumNodes){ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[3])>=NumNodes){ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[4])>=NumNodes){ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[1])<0) 			 {ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad node index)",BAD_DATA);}
		if (s_to_i(s[2])<0) 			 {ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad node index)",BAD_DATA);}
		this->side[i].n1=s_to_i(s[1]);
		this->side[i].n2=s_to_i(s[2]);
		this->side[i].na=s_to_i(s[3]);
		this->side[i].nb=s_to_i(s[4]);

		if (s_to_i(s[5])>=NumElems){ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad elem index)",BAD_DATA);}
		if (s_to_i(s[6])>=NumElems){ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad elem index)",BAD_DATA);}
		this->side[i].ea=s_to_i(s[5]);
		this->side[i].eb=s_to_i(s[6]);

		this->side[i].status =TranslateMark(s_to_i(s[7]));

		this->side[i].radius=0.0;
		if (Len==9){
		  this->side[i].radius=s_to_d(s[8]);
			if (side[i].radius!=0.0){
				if (fabs(side[i].radius)<0.5*abs(node[side[i].n1].z-node[side[i].n2].z) ){
					ExitGracefully("CTriMesh:BTM File (sides) incorrectly formatted (bad radius)",BAD_DATA);}
			}
		}
		i++;
	}

	n=0;
	while ((!eof) && (n<this->NumDoublenodes)){
		eof=TokenizeLine(BTM,s,Len);
		if ((eof) || (Len!=4)) {ExitGracefully("CTriMesh:Doublenodes File incorrectly formatted",BAD_DATA);}
		this->dnodes[n].n1=s_to_i(s[1]);
		this->dnodes[n].n2=s_to_i(s[2]);
		this->dnodes[n].type=DNODE_JUMP;//TranslateDnode(s_to_i(s[2]));//TMP DEBUG
		n++;
	}

	for (e=0; e<NumElems;e++){
		CalculateCircles(e);
	}
	for (n=0; n<NumNodes;n++){
		upperswap(BBox.e,this->node[n].z.real());
		upperswap(BBox.n,this->node[n].z.imag());
		lowerswap(BBox.w,this->node[n].z.real());
		lowerswap(BBox.s,this->node[n].z.imag());
	}
	Process();
	
	cout <<"...Done Reading Triangular Mesh."<<endl;

	return true; 
}
/**************************************************************************
	WriteOutput (OLD)
***************************************************************************
	Writes mesh to 
	 .n (nodes file) Format:	
			[NumNodes]
			{node: x y mark}xNumNodes

	 .e (nodes file) Format:
			[NumElems]
			{[elem]: [ni] [nj] [nk] [ei] [ej] [ek] [si] [sj] [sk] [xv] [yv] [0]} x NumElems 

	 .s (sides file) Format:
			[NumSides]
			{[side]: [n1] [n2] [na] [nb] [ea] [eb] [mark] {optional radius}}xNumSides

	 .d (doublenodes file) format: (not written if NumDoublenodes=0)
			[NumDoublenodes]
			{[dnode]: [n1] [n2] [type(mark)]}xNumDoublenodes
-------------------------------------------------------------------------*/
void CTriMesh::WriteOutput(){

	int mark;
	char filename[100];
	int n,s,e;
	cout<<"Mesh name: " <<name<<"| "<<strlen(name)<<endl;

	//	Node data-------------------------------------------------- 
	sprintf(filename,"ConcMesh.n.%s",name);
	ofstream NODES;
	NODES.open(filename);
	if (NODES.fail()){cout << "Unable to open nodes output file"<<endl;return;}  
	
	cout <<"Writing nodes output file"<<endl;
	NODES<<NumNodes<<endl;
	for(n=0; n<NumNodes; n++){
		if			(node[n].status==INTERIOR  ){mark=0;}
		else if (node[n].status==BOUNDARY  ){mark=1;}
		else if (node[n].status==DOUBLENODE){mark=2;}
		else															{mark=OFF;}
		NODES << n									 <<":  ";
		NODES << node[n].z.real() 	 <<" "<< node[n].z.imag() 	 << " ";
		NODES << mark 							 <<" ";
		NODES << node[n].inv.real()  <<" "<< node[n].inv.imag()  << endl;
	}
	NODES<<"----------------------------------------------------------"<<endl;
	NODES<<"   n:  x      y     mark    lagx   lagy"<<endl;

	NODES.close();

	// -- Element data	-----------------------------------------------
	sprintf(filename,"%s.e"  ,name);
	ofstream ELEMS;
	ELEMS.open("ConcMesh.e");
	if (ELEMS.fail()){cout << "Unable to open elements output file"<<endl;return;}	

	cout <<"Writing element output file"<<endl;
	ELEMS<<NumElems<<endl;
	for(e=0; e<NumElems; e++){
		ELEMS << e								 <<": ";
		ELEMS << elem[e].ni 			 <<" "<< elem[e].nj 			 <<" "<< elem[e].nk <<" ";
		ELEMS << elem[e].ei 			 <<" "<< elem[e].ej 			 <<" "<< elem[e].ek <<" ";
		ELEMS << elem[e].si 			 <<" "<< elem[e].sj 			 <<" "<< elem[e].sk <<" ";
		ELEMS << elem[e].zv.real() <<" "<< elem[e].zv.imag() <<" "<< "0"				<<endl;
	}
	ELEMS<<"---------------------------------------------------";
	ELEMS<<"-------------------------------------------------------"<<endl;
	ELEMS<<"   e:   ni,   nj,   nk,   ei,  ej,  ek,   si,  sj,  sk";
	ELEMS<<"   xV,                    yV                       sign"<<endl;

	ELEMS.close();

	//--	Side data  --------------------------------------------------------
	sprintf(filename,"%s.s"  ,name);
	ofstream SIDES;
	SIDES.open("ConcMesh.s");
	if (SIDES.fail()){cout << "Unable to open sides output file"<<endl;return;}  

	cout <<"Writing sides output file"<<endl;
	SIDES << NumSides <<endl;
	for(s=0; s<NumSides; s++){
		if			(side[s].status==INTERIOR  ){mark=0;}
		else if (side[s].status==BOUNDARY  ){mark=1;}
		else if (side[s].status==DOUBLENODE){mark=2;}
		else																{mark=OFF;}
		SIDES << s					 << ":  ";
		SIDES << side[s].n1  << " " << side[s].n2 	<< " ";
		SIDES << side[s].na  << " " << side[s].nb 	<< " ";
		SIDES << side[s].ea  << " " << side[s].eb  << " ";
		SIDES << mark;
		if (side[s].radius!=0.0){SIDES<<" "<<side[s].radius;}
		SIDES << endl;
	}
	SIDES << "   s:    n1    n2    na   nb   ea   eb   mark"<<endl;
	SIDES << "--------------------------------"<<endl;

	SIDES.close();

	//	Doublenode data-------------------------------------------------- 
	//if (NumDoublenodes>0){
	sprintf(filename,"ConcMesh.d.%s",name);
	ofstream DNOD;
	DNOD.open(filename);
	if (DNOD.fail()){cout << "Unable to open doublenodes output file"<<endl;return;}	
	
	cout <<"Writing doublenodes output file"<<endl;
	DNOD<<NumDoublenodes<<endl;
	for(n=0; n<NumDoublenodes; n++){
		if			(dnodes[n].type==DNODE_JUMP 	 ){mark=0;}
		else if (dnodes[n].type==DNODE_FLUX 	 ){mark=1;}
		else if (dnodes[n].type==DNODE_BOUNDARY){mark=2;}
		else																		{mark=OFF;}
		DNOD << n 									<<":  ";
		DNOD << dnodes[n].n1				<<" ";
		DNOD << dnodes[n].n2		<< " ";
		DNOD << mark								<< endl;
	}
	DNOD<<"----------------------------------------------------------"<<endl;
	DNOD<<"   dn:  n1                      n2                     mark"<<endl;

	DNOD.close();
	//}
}
/**************************************************************************
	WriteOutput2
***************************************************************************
	Writes mesh to 
	 .btm (bluebird triangular mesh) Format:	
			[NumNodes]
			[NumElems]
			[NumSides]
			[NumDoublenodes]
			{[node]:	[x]  [y]	[mark] [lagx] [lagy]}xNumNodes
			{[elem]:	[ni] [nj] [nk] [ei] [ej] [ek] [si] [sj] [sk] [xv] [yv] [0]} x NumElems 
			{[side]:	[n1] [n2] [na] [nb] [ea] [eb] [mark]} x NumSides
			{[dnode]: [n1] [n2] [type(mark)]} x NumDoublenodes
-------------------------------------------------------------------------*/
void CTriMesh::WriteOutput2(){

	int mark;
	char filename[100];
	int n,s,e;
	//cout<<"Mesh name: " <<name<<"| "<<strlen(name)<<endl;

	//	Node data-------------------------------------------------- 
	sprintf(filename,"ConcMesh.btm.%s",name);
	ofstream BTM;
	BTM.open(filename);
	BTM.precision(16);
	if (BTM.fail()){cout << "Unable to create/open triangular mesh output file"<<endl;return;}	
	
	cout <<"Writing triangular mesh output file"<<endl;
	
	BTM << NumNodes 			<<endl;
	BTM << NumElems 			<<endl;
	BTM << NumSides 			<<endl;
	BTM << NumDoublenodes <<endl;

	cout <<"Writing node data"<<endl;
	for(n=0; n<NumNodes; n++){
		if			(node[n].status==INTERIOR  ){mark=0;}
		else if (node[n].status==BOUNDARY  ){mark=1;}
		else if (node[n].status==DOUBLENODE){mark=2;}
		else															{mark=OFF;}
		BTM << n									 <<":  ";
		BTM << node[n].z.real() 	 <<" "<< node[n].z.imag() 	 << " ";
		BTM << mark 							 <<" ";
		BTM << node[n].inv.real()  <<" "<< node[n].inv.imag()  << endl;
	}

	// -- Element data	-----------------------------------------------
	cout <<"Writing element data"<<endl;

	for(e=0; e<NumElems; e++){
		BTM << e								 <<": ";
		BTM << elem[e].ni 			 <<" "<< elem[e].nj 			 <<" "<< elem[e].nk <<" ";
		BTM << elem[e].ei 			 <<" "<< elem[e].ej 			 <<" "<< elem[e].ek <<" ";
		BTM << elem[e].si 			 <<" "<< elem[e].sj 			 <<" "<< elem[e].sk <<" ";
		BTM << elem[e].zv.real() <<" "<< elem[e].zv.imag() <<" "<< "0"				<<endl;
	}


	//--	Side data  --------------------------------------------------------
	cout <<"Writing sides data"<<endl;

	for(s=0; s<NumSides; s++){
		if			(side[s].status==INTERIOR  ){mark=0;}
		else if (side[s].status==BOUNDARY  ){mark=1;}
		else if (side[s].status==DOUBLENODE){mark=2;}
		else																{mark=OFF;}
		BTM << s					 << ":  ";
		BTM << side[s].n1  << " " << side[s].n2 	<< " ";
		BTM << side[s].na  << " " << side[s].nb 	<< " ";
		BTM << side[s].ea  << " " << side[s].eb  << " ";
		BTM << mark;
		if (side[s].radius!=0.0){
		BTM<<" "<<side[s].radius;
		}
		BTM << endl;
	}

	//	Doublenode data-------------------------------------------------- 
	cout <<"Writing doublenodes output file"<<endl;

	for(n=0; n<NumDoublenodes; n++){
		if			(dnodes[n].type==DNODE_JUMP 	 ){mark=0;}
		else if (dnodes[n].type==DNODE_FLUX 	 ){mark=1;}
		else if (dnodes[n].type==DNODE_BOUNDARY){mark=2;}
		else																		{mark=OFF;}
		BTM << n									 <<":  ";
		BTM << dnodes[n].n1 			 <<" ";
		BTM << dnodes[n].n2 	 << " ";
		BTM << mark 							 << endl;
	}

	BTM.close();

}
