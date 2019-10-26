//DynamicTriMesh.cpp
#include "TriagonalMesh.h"
/*******************************************************************************
           CONSTRUCTORS
********************************************************************************
Constructor #1- builds mesh based upon boundary spacing
OUTER BOUNDARY MUST BE LABELED AS CHAIN 0
------------------------------------------------------------------------------*/
CTriMesh::CTriMesh(char              *Name,            //Name of mesh (&output file)
									 const nodestruct  *points,          //initial points (w/ spacing function, status)
									 const int          nP,              //number of points
                   const segstruct   *segments,        //boundary & internal segment info
              		 const int          nS,              //number of segments
									 const int          maximum_nodes,   //maximum number of nodes
					         const spcstruct    spacing_package){ 

	int					 n,s,c,m;															 //counters
	cmplex			 z0;																	 //first point of segment
	cmplex			 zN;																	 //second point of segment
	cmplex			 L;																		 //segment as vector (L=(zN-z0))
	int          NumSegs;															 //Number of segments
	int         *segsperpoint;												 //segments associated with point
	nodestruct  *point;																 //array of points
	int         *ptnode;                               //node associated with point
	segstruct   *segment=NULL;                         //array of segments
	int         *seg_chainind=NULL;                    //indices of chains in which segment belongs
	chainstruct *chain=NULL;	                         //array of chains (dynamic creation only)
	int					 NumChains;	                           //Number of Chains (dynamic creation only)
	int          node_chainind[MAX_DYN_PTS];
	double       dLm, ddL, L_tot;
	int          M;
  int					 NumPoints;	                           //Number of Points (dynamic creation only) 

	if (Name!=NULL){
    name = new char[strlen(Name)+1];
    name=strcpy(name,Name);
	}

	minSpacing=ALMOST_INF;
	nCells=0;
	nNodes=0;
	BBox.e=-ALMOST_INF; BBox.w=ALMOST_INF;
	BBox.n=-ALMOST_INF; BBox.s=ALMOST_INF;
	is3D=false;

	spacinginfo.type=USER_SPECIFIED;
	spacinginfo.min_spacing=0;
	spacinginfo.max_spacing=0;

	MaxNodes=maximum_nodes;
	ugly   =OFF;
	elem   =NULL;	//created in initialize
	side   =NULL;	//created in initialize
	node   =NULL;	//created in initialize
	sumz   =NULL; //created in initialize
	sumF   =NULL; //created in initialize
  Nne    =NULL; //created in initialize

	NumPoints   =nP;
	NumSegs = nS;
  if (NumPoints<=0){cout << "CTriMesh::Constructor(Dynamic):Invalid number of points specified"<<endl;}
	if (NumSegs<=0)  {cout << "CTriMesh::Constructor(Dynamic):Invalid number of segments specified"<<endl;}

	segsperpoint=new int         [NumPoints+1];
	point       =new nodestruct  [NumPoints+1];
	ptnode      =new int         [NumPoints+1];
	seg_chainind=new int         [NumSegs  +1];
	segment     =new segstruct   [NumSegs  +1];
	chain       =new chainstruct [NumSegs  +1]; 
	ExitGracefullyIf (chain==NULL,"CTriMesh::Constructor(Dynamic): Out of memory",OUT_OF_MEMORY);

	//copy points, identify domain extents--------------------------------------
	for(n=0; n<NumPoints; n++)
	{
	  point[n].z        =points[n].z;
	  point[n].F        =points[n].F;
		point[n].status   =points[n].status;

		if (spacing_package.type==COURANT){
			point[n].F=point[n].F; //COURANTSPACING UPDATE
		}
		lowerswap(minSpacing,point[n].F);
		ExitGracefullyIf(point[n].F<=0.0,"CTriMesh::Constructor(Dynamic): Bad (negative or zero) spacing function specified",BAD_DATA);
	  
		ptnode[n]=OFF;                          //indicates node associated with point (initialized to none)
		upperswap(BBox.e,point[n].z.real());
		upperswap(BBox.n,point[n].z.imag());
		lowerswap(BBox.w,point[n].z.real());
		lowerswap(BBox.s,point[n].z.imag());
		segsperpoint[n]=0;
	}
	minSpacing/=10.0;

	//copy segments to array of segment structures------------------------------
	for(s=0; s<NumSegs; s++)
	{
	 	segment[s].p0    =segments[s].p0;
		segment[s].p1    =segments[s].p1;
	  segment[s].status=segments[s].status;
		segment[s].radius=segments[s].radius;
		if ((segment[s].p0<0) || (segment[s].p0>NumPoints) || 
			  (segment[s].p1<0) || (segment[s].p1>NumPoints)) {
			ExitGracefully("CTriMesh::Constructor(Dynamic): bad segment point index specified (<0 or >Numpoints)",BAD_DATA);}
		
		if ((segment[s].radius!=0) && (fabs(segment[s].radius)<0.5*abs(point[segment[s].p0].z-point[segment[s].p1].z))){
			ExitGracefully("CTriMesh::Constructor(Dynamic): bad segment radius specified: radius less than 1/2 length of segment",BAD_DATA);}
	}
	segment[NumSegs].p0=OFF;
	segment[NumSegs].p1=OFF;

	//count the chains----------------------------------------------------------
	NumChains=0;
	chain[NumChains].s0=0;
	for(s=0; s<NumSegs; s++)
	{
		segsperpoint[segment[s].p0]++;
		segsperpoint[segment[s].p1]++;

    seg_chainind[s]=NumChains;

		if(segment[s].p1!=segment[s+1].p0){//chain ended, start a new one
			chain[NumChains].s1=s;
			NumChains++;
			chain[NumChains].s0=s+1;
		}
  }

	//adjust spacing functions on each small segment-----------------------------
	double length;
	for(s=0; s<NumSegs; s++)
	{
		z0=point[segment[s].p0].z; 
		zN=point[segment[s].p1].z; 

		length=abs(zN-z0);
		//CORRECT FOR CURVED SEGMENTS
		//if (segment[s].radius!=0){length=2.0*segment[s].radius*asin(abs(zN-z0)/2.0/segment[s].radius);}

		if((point[segment[s].p0].F+point[segment[s].p1].F > length ) &&
			 (segment[s].p0 != segment[s].p1))
		{
			point[segment[s].p0].F = min(point[segment[s].p0].F,length);
			point[segment[s].p1].F = min(point[segment[s].p1].F,length);
		}
  }

  //set the number of nodes on each segment------------------------------------
	for(s=0; s<NumSegs; s++)
	{
		z0=point[segment[s].p0].z; 
		zN=point[segment[s].p1].z;

		length=abs(zN-z0);
		//CORRECT FOR CURVED SEGMENTS
		//if (segment[s].radius!=0){length=2.0*segment[s].radius*asin(abs(zN-z0)/2.0/segment[s].radius);}
		
		//segment with many intermediate points
		if(point[segment[s].p1].F+point[segment[s].p0].F<=length)
		{
			dLm=0.5*(point[segment[s].p0].F+point[segment[s].p1].F);
			segment[s].N=(int)(ceil(length/dLm));
		}  
		else//single point "segment" or segment less than spacing function
		{ 
			segment[s].N=1;
		}
	}

	//Identify chain type----------------------------------------------------------
	for(n=0; n<NumChains; n++)
	{
		if( segment[chain[n].s0].p0 == segment[chain[n].s1].p1 ){chain[n].type=CLOSED;}
		else                                                    {chain[n].type=OPEN;}
		if( (segsperpoint[segment[chain[n].s0].p0]==1) &&
				(segsperpoint[segment[chain[n].s1].p1]==1) )        {chain[n].type=INSIDE;}
	}

	cout <<"Assembling Unrefined Mesh..."<<endl;
	
	//create, initialize arrays, start domain as triangle around bounding box
	Initialize();

	node_chainind[0]=0;
  node_chainind[1]=0;
  node_chainind[2]=0;

	//Inserting------------------------------------------------------------------
	//---------------------------------------------------------------------------
	cmplex unit_in(0.0),last_unit_in(0.0),unit_vec;
	int    firstsegpt,lastsegpt;

	for(c=0; c<NumChains; c++)											//for each chain
	{                    
		unit_in=0.0;
		if (chain[c].type ==CLOSED)
		{
			firstsegpt=segment[chain[c].s1].p0;
			lastsegpt =segment[chain[c].s1].p1;
			L=point[lastsegpt ].z-point[firstsegpt].z;
			length=abs(zN-z0);
			//CORRECT FOR CURVED SEGMENTS
		  //if (segment[s].radius!=0){length=2.0*segment[s].radius*asin(abs(zN-z0)/2.0/segment[s].radius);}
			unit_in=IM*L/length;
			unit_vec=L/length;
		}
		for(s=chain[c].s0; s<=chain[c].s1; s++)       //for each segment in chain
		{     

			firstsegpt=segment[s].p0;
			lastsegpt =segment[s].p1;

			z0=point[firstsegpt].z; 
			zN=point[lastsegpt ].z;

			L=zN-z0;
			length=abs(zN-z0);
			//CORRECT FOR CURVED SEGMENTS
		  //if (segment[s].radius!=0){length=2.0*segment[s].radius*asin(abs(zN-z0)/2.0/segment[s].radius);}

			dLm=length/segment[s].N;

			last_unit_in=unit_in;
			unit_in =IM*L/length;
			unit_vec=   L/length;
			
			if      ((dLm==0.0) && (firstsegpt==lastsegpt)) //zero length side- single point "segment"
			{
				InsertNode(z0,point[firstsegpt].status,0.0);

				node         [NumNodes-1].F	=point[firstsegpt].F;
				node_chainind[NumNodes-1]   =seg_chainind[s]; // NumNodes-1 is index of inserted node 
				ptnode[firstsegpt]          =NumNodes-1;      //associate point with node
			}
			else if (dLm==0.0)
			{
				ExitGracefully("CTriMesh Constructor: zero length side specified",BAD_DATA);
			}
			else
			{
				//if spacing functions of segment endpoints are less than the length of the segment
				//delta length is modified
				if(point[firstsegpt].F + point[lastsegpt].F <= length)
				{ 
					if(point[firstsegpt].F > point[lastsegpt].F)
					{
						M=-segment[s].N/2; 
						ddL=(point[lastsegpt].F-dLm)/M;
					}
					else
					{
						M=+segment[s].N/2; 
						ddL=(dLm-point[firstsegpt].F)/M;
					}
				}

				//--first point of segment------------------------------------------------
				if(ptnode[firstsegpt]==OFF)                 // if first point in segment not already associated with node
				{             

					if     (s==chain[c].s0)                   // first segment in the chain-just insert 
					{                      
						InsertNode(z0,point[firstsegpt].status,(unit_in+last_unit_in)/2.0);
					}
					else if(s==chain[c].s1 && segment[s].N==1)//last segment in the chain and only point in segment
					{   
						InsertNode(z0,point[firstsegpt].status,unit_in);//must revise vector
						UpdateNewSideStatus(NumNodes-2                     ,segment[s-1].status,segment[s-1].radius);
						UpdateNewSideStatus(ptnode[segment[chain[c].s0].p0],segment[s  ].status,-segment[s  ].radius);
					}
					else                                     //all other segment starting points-insert and adjust side between last node and this one
					{                                     
						InsertNode(z0,point[firstsegpt].status,(unit_in+last_unit_in)/2.0);
						UpdateNewSideStatus(NumNodes-2,segment[s-1].status,segment[s-1].radius);
					}

					node         [NumNodes-1].F	=point[firstsegpt].F;
					node_chainind[NumNodes-1]   =seg_chainind[s]; // NumNodes-1 is index of inserted node 
					ptnode[firstsegpt]          =NumNodes-1;      //associate point with node

				}//end if( point[firstsegpt].new_numb == OFF ){

				//--- middle points in this segment-------------------------------------------
				L_tot=0;
				if (point[firstsegpt].F + point[lastsegpt].F <= length)
				{
					for(m=1; m<segment[s].N; m++){ //sifts through points in segment
					
						L_tot+=(dLm-M*ddL);
						if(point[firstsegpt].F > point[lastsegpt].F){
							M++; if(M==0 && segment[s].N%2==0) {M++;}
						}
						else{
							M--; if(M==0 && segment[s].N%2==0) {M--;}
						}

						if(s==chain[c].s1 && m==(abs(segment[s].N)-1)) //if first segment in chain and last node in middle	
						{ 																										
							InsertNode(z0+unit_vec*L_tot,segment[s].status,unit_in);
							UpdateNewSideStatus(NumNodes-2       ,segment[s].status,segment[s].radius);
							UpdateNewSideStatus(ptnode[lastsegpt],segment[s].status,segment[s].radius);
						}
						else if(m==1)												       		//if first node in middle
						{			
							InsertNode(z0+unit_vec*L_tot,segment[s].status,unit_in);
							UpdateNewSideStatus(ptnode[firstsegpt],segment[s].status,segment[s].radius);
						}
						else				  																//all other middle nodes and first node if not first segment
						{		
							InsertNode(z0+unit_vec*L_tot,segment[s].status,unit_in);
							UpdateNewSideStatus(NumNodes-2           ,segment[s].status,segment[s].radius);
						}

						node         [NumNodes-1].F	=0.5*(node[NumNodes-2].F + (dLm-M*ddL));
						node_chainind[NumNodes-1]   =seg_chainind[s]; 
					}// end for m...
				}//end if (point[firstsegpt].F

				//--- last point of segment (just for the open interior chains)-----------------------
				if((ptnode[lastsegpt] == OFF) && 
					 (s==chain[c].s1))//if node doesnt exist yet and this is the last segment in the chain
				{
					InsertNode(zN,point[lastsegpt].status,unit_in);

					UpdateNewSideStatus(NumNodes-2           ,segment[s].status,segment[s].radius);

					node         [NumNodes-1].F	=point[lastsegpt].F;
					node_chainind[NumNodes-1]   =seg_chainind[s]; 
				}
				//else, must revise node.inv for first point

			}//end if (dLm==0.0...
		}//end for (s=0...
  }//end for (c=0...

	cout <<"...Unrefined Mesh Assembled with "<<NumNodes-3 << " nodes."<<endl;

	int e, a, b;

	//  Negative area check for elimination of elements ---------------------------------- 
	for(e=0; e<NumElems; e++){
		
		if((node_chainind[elem[e].ni]==node_chainind[elem[e].nj]) &&
       (node_chainind[elem[e].nj]==node_chainind[elem[e].nk]) &&
       (chain[node_chainind[elem[e].ni]].type==CLOSED) )
		{ 

			a = min( min(elem[e].ni, elem[e].nj), elem[e].nk );//minimum node index
			c = max( max(elem[e].ni, elem[e].nj), elem[e].nk );//maximum node index
			b = elem[e].ni+elem[e].nj+elem[e].nk - a - c;      //middle  node index
     
			if ((b<0) || (b>NumNodes-1) ){cout << "CleanMesh: bad node(2)"<<endl;}
			
			if(a<3){
				elem[e].on=false;                                      //checks if associated with initial element )i,j,k<3)
			}
			else if(TriArea(node[a].z, node[b].z, node[c].z) <=0.0){ //checks for negative area
				elem[e].on=false;                                      //extremely important step!!
			}
		}
	}

	//--Elimination of sides & nodes---------------------------
	for(s=0; s< 3; s++){side[s].status=IS_OFF;}
	for(n=0; n< 3; n++){node[n].status=IS_OFF;}
		
	CleanIndices();

	//cout <<"...Unrefined Mesh Cleaned"<<endl;

  delete [] segment;
  delete [] segsperpoint;
	delete [] seg_chainind;
	delete [] point;
	delete [] chain;

}
/**************************************************************************
														Dynamic Initialize
***************************************************************************
  Sets parameters for non-user specified mesh generation techniques
--------------------------------------------------------------------------*/
void   CTriMesh::DynamicInitialize(spc_type  type, 
																	 double param,
																	 double minx, 
																	 double maxx, 
																	 double radius){
	spacinginfo.type=type;
	spacinginfo.min_spacing=minx;
	spacinginfo.max_spacing=maxx;
	spacinginfo.time_step=param;
	spacinginfo.radius=radius;
}
/**************************************************************************
														Dynamic Classify
***************************************************************************
  This function defines the strategy for insertion of new nodes    
  This function searches through all elements every time.  
  Some optimisation will definitely be needed                        
	Called only from DynamicBuildMesh
--------------------------------------------------------------------------*/
void CTriMesh::DynamicClassify(){
	int		 ea,  eb;
	int    ei,	ej,	 ek;
	int    si,	sj,	 sk;
	int		 eac, ead, ebc, ebd;
	int		 e,s;
  int    uglytemp;
	double ratio(-ALMOST_INF), F;
	int n;
	int lastugly=ugly;

	//cout <<NumElems<<" "<<NumSides<<" "<<NumNodes<<endl;
	for(e=0; e<NumElems; e++){
		if(elem[e].on==true){
			ei=elem[e].ei; 
			ej=elem[e].ej; 
			ek=elem[e].ek;

			F=(node[elem[e].ni].F + node[elem[e].nj].F + node[elem[e].nk].F)/3.0;

			//if (spacing_package.type=COURANT){
			//}
			
			//set all elements to "waiting" state
			elem[e].state=WAITING;

			//if element radius is less than 0.7 times the average spacing function, it is "done"
		  //(R=0.577*F is ideal triangle)
			if(elem[e].R < 0.7*F){elem[e].state=DONE;} // 0.0866; 0.07 

			// if all elements surrounding an element are done, the element is "done"
			if(ei!=OFF && ej!=OFF && ek!=OFF){
				if(elem[ei].state==DONE && 
					 elem[ej].state==DONE && 
					 elem[ek].state==DONE){
					elem[e].state=DONE;
				}
			}
    }
	}

	//--Diamond check------------------------------------------------------
  //if all four sides of an "element diamond" are off or done, the element is done
	for(s=0; s<NumSides; s++){
		if(side[s].status!=IS_OFF){
			ea=side[s].ea;
			eb=side[s].eb;

			if (ea==OFF)          {ead=OFF;         eac=OFF;        } 
			else{
				if(elem[ea].ei==eb) {ead=elem[ea].ej; eac=elem[ea].ek;}
				if(elem[ea].ej==eb) {ead=elem[ea].ek; eac=elem[ea].ei;}   
				if(elem[ea].ek==eb) {ead=elem[ea].ei; eac=elem[ea].ej;}
			}
			if (eb==OFF)          {ebc=OFF;         ebd=OFF;        }
			else{
				if(elem[eb].ei==ea) {ebc=elem[eb].ej; ebd=elem[eb].ek;}
				if(elem[eb].ej==ea) {ebc=elem[eb].ek; ebd=elem[eb].ei;}
				if(elem[eb].ek==ea) {ebc=elem[eb].ei; ebd=elem[eb].ej;}
			}
			
			if( (eac==OFF || elem[eac].state==DONE) &&
					(ebc==OFF || elem[ebc].state==DONE) &&
					(ead==OFF || elem[ead].state==DONE) &&
					(ebd==OFF || elem[ebd].state==DONE) ){
				if (ea!=OFF){elem[ea].state=DONE;}
				if (eb!=OFF){elem[eb].state=DONE;}
			}
		}
	}

	//--Find Ugly element-------------------------------------------------
	//search through the un-"done" boundary elements
	ugly=OFF;
	uglytemp=OFF;
	for(e=0; e<NumElems; e++){
		if(elem[e].on==true && elem[e].state!=DONE){
			
			si=elem[e].si; 
			sj=elem[e].sj; 
			sk=elem[e].sk;

			//any boundary/doublenode element not done is instantly active
			if(side[si].status!=INTERIOR){elem[e].state=ACTIVE;}
			if(side[sj].status!=INTERIOR){elem[e].state=ACTIVE;}
			if(side[sk].status!=INTERIOR){elem[e].state=ACTIVE;}
  
			//if it is active, and it has a real bad radii, it is the ugly element
			//TMP DEBUG ! Must check if element circumcenter is inside element!!!!!!!!
			if(elem[e].state==ACTIVE && elem[e].R/elem[e].r > ratio){
				uglytemp=e;
				s=DynamicGetDoneSide(uglytemp,n);
				//HOW COULD THIS POSSIBLY HAPPEN!?!?!?
				ExitGracefullyIf((side[s].n1>NumNodes) || (side[s].n2>NumNodes) || (s>NumSides),"Dynamic Classify:Bad node in side",RUNTIME_ERR);
				
				if(TriArea(node[side[s].n1].z, node[side[s].n2].z, elem[uglytemp].zv) *
					 TriArea(node[side[s].n1].z, node[side[s].n2].z, node[n].z ) > 0.0 ){
					upperswap(ratio, elem[e].R/elem[e].r);
					ugly=uglytemp;
				}
			}
    }
	}


	// if ugly element is not found on the boundary 
	// ugly element can be any element next to a "done" element
	if(ugly==OFF){
		for(e=0; e<NumElems; e++){
			if(elem[e].on==true && elem[e].state!=DONE){
				ei=elem[e].ei; 
				ej=elem[e].ej; 
				ek=elem[e].ek;

				if(ei!=OFF){if(elem[ei].state==DONE){ elem[e].state=ACTIVE;}}
				if(ej!=OFF){if(elem[ej].state==DONE){ elem[e].state=ACTIVE;}}
				if(ek!=OFF){if(elem[ek].state==DONE){ elem[e].state=ACTIVE;}}
  
				//if it is active, and it has a real bad radii, it is the ugly element
				if((elem[e].state==ACTIVE) && 
					 ((elem[e].R/elem[e].r) > ratio)){ 
					
					uglytemp=e;

					s=DynamicGetDoneSide(uglytemp,n);

					if(TriArea(node[side[s].n1].z, node[side[s].n2].z, elem[uglytemp].zv) *
						 TriArea(node[side[s].n1].z, node[side[s].n2].z, node[n].z ) > 0.0 ){
						upperswap(ratio, elem[e].R/elem[e].r);
						ugly=uglytemp;
					}
				}
			}
		}
	}
	//if (GetElement(elem[ugly].zv)!=ugly){cout <<"CTriMesh::Dynamic Classify: bad circumcenter"<<endl;}
	//if (ugly==OFF){cout <<"CTriMesh::Dynamic Classify: Ugly node not found"<<endl;}
}
/**************************************************************************
									Dynamic Get Done Side
***************************************************************************
  This function returns the side next to the ugly element (ug) which is done
	it also returns n_across, which is the node across from this side
-------------------------------------------------------------------------*/
int CTriMesh::DynamicGetDoneSide(const int ug, int &n_across) const{
	int s(OFF),n(OFF);

	// find any side of ugly element next to "Done" element
	if (elem[ug].ei!=OFF){
		if(elem[elem[ug].ei].state==DONE)				{s=elem[ug].si; n=elem[ug].ni;}}
  if (elem[ug].ej!=OFF){
		if(elem[elem[ug].ej].state==DONE)				{s=elem[ug].sj; n=elem[ug].nj;}}
  if (elem[ug].ej!=OFF){
		if(elem[elem[ug].ek].state==DONE)				{s=elem[ug].sk; n=elem[ug].nk;}}

	//boundary elements override this 
	if ((side[elem[ug].si].status==BOUNDARY) ||
	    (side[elem[ug].si].status==DOUBLENODE)){s=elem[ug].si; n=elem[ug].ni;}
	if ((side[elem[ug].sj].status==BOUNDARY) ||
	    (side[elem[ug].sj].status==DOUBLENODE)){s=elem[ug].sj; n=elem[ug].nj;}
	if ((side[elem[ug].sk].status==BOUNDARY) || 
	    (side[elem[ug].sk].status==DOUBLENODE)){s=elem[ug].sk; n=elem[ug].nk;}
	n_across=n;
	return s;
}
/**************************************************************************
									Dynamic New Node
***************************************************************************
  This function determines the position of the inserted node.  
	used only for dynamic creation of mesh with boundary spacing:
		not used for simple delunay triangulation of points
	Called only from DynamicBuildMesh
-------------------------------------------------------------------------*/
bool CTriMesh::DynamicNewNode(){

	int    s(OFF);int e;
	int    n;
	cmplex zmid;      //midpoint of side
	double s_length;  //length of side
	cmplex q;         //direction from midpoint to insert node
	double rhoM;
	double rho_M;
	double delta;     //distance from midpoint to insert node 

	cmplex uglyvertex;

	if (ugly==OFF){cout<< "DynamicNewNode:No bad elements to repair"<<endl;return false;}

	//  elements which are near external or internal boundaries will come into play first.                                                                                                                                      
	//  However, some attention has to be paid for the case when two accepted  
	//  elements surround the ugly one                                                                                                                                               

	// find any side of ugly element next to "Done" or boundary element
	s=DynamicGetDoneSide(ugly,n);

	if (s==OFF) {cout<< "DynamicNewNode:Element has side turned off or no adjacent elements"<<endl;return false;}

	s_length   = 0.5*abs(node[side[s].n1].z - node[side[s].n2].z);
	zmid       = 0.5*   (node[side[s].n1].z + node[side[s].n2].z);	
	uglyvertex = elem[ugly].zv;

	//place node in direction of ugly vertex
	q=uglyvertex-zmid;

	//obtain distance from midpoint to insert node (delta)
	rhoM = 0.577 *  0.5*(node[side[s].n1].F + node[side[s].n2].F);

	rho_M = min( max( rhoM, s_length), 0.5*(s_length*s_length+abs(q)*abs(q))/abs(q) );

	if(rho_M < s_length){delta=rho_M;                                    }
	else                {delta=rho_M+sqrt(rho_M*rho_M-s_length*s_length);} 

	// The following line checks whether the new point falls outside the domain.                                                
	// if the ugly element circumcenter is on the same side of the side as the node, 
	// then insert at circumcenter 
	if (abs(delta*q/abs(q))<minSpacing){cout <<"minSpacing beat"<<endl;}
	if(TriArea(node[side[s].n1].z, node[side[s].n2].z, uglyvertex) *
		 TriArea(node[side[s].n1].z, node[side[s].n2].z, node[n].z ) > 0.0 ){
		
		e=InsertNode(zmid+delta*q/abs(q),INTERIOR,0.0);
		//cout <<"ugly:"<< ugly <<" "<< delta<<" "<<s<<endl;
    //if ((ugly==1252) && (delta<21.5)){return false;}
		if (e!=OFF){return true;}
		else       {return false;}
	}
	else {
		//this happens when flat side of obtuse triangle is side of interest
		//could place node @ midpoint?
		cout<< "DynamicNewNode: poorly oriented circumcenter @ x="<<node[n].z.real() <<" y="<<node[n].z.imag()<<endl;
		//e=InsertNode((zmid-(zmid-node[n].z)/2.0),INTERIOR);
		/*int e;
		e=InsertNode(zmid,INTERIOR);
		if (e!=OFF){return true;}
		else       {return false;}
*/
		return false;
	}
/*
	else{ //move opposing node closer to c 
		node[n].z = zmid - delta*q/abs(q);
		node[n].mark=6;   
		for(e=0; e<NumElems; e++) {
			if(elem[e].ni==n || elem[e].nj==n || elem[e].nk==n){
				CalculateCircles(e);
			}
		}
	}
*/

}
/**************************************************************************
									Dynamic Smooth
***************************************************************************
	Cleans up mesh after everthing has been added
	Should be used only for creation of dynamic mesh with boundary spacing:
		not used for simple delunay triangulation of points
-------------------------------------------------------------------------*/
void CTriMesh::DynamicSmooth(bool relax, bool mild){
	int it, s, n, e;
 	int numiter;
	//double *sumW;
	//sumW =new double [NumNodes];
	for (n=0;n<NumNodes;n++){Nne[n]=0;}

	//counts node neighbors
	for(s=0; s<NumSides; s++){
		//if (side[s].status==INTERIOR){ //necc?
    if (side[s].status!=IS_OFF){ 
			Nne[side[s].n1]++;
			Nne[side[s].n2]++;
		}
	}

	//relax mesh
  if (relax) {
	int T, E;
	for(T=6; T>=3; T--){
		for(s=0; s<NumSides; s++){
			if(side[s].status==INTERIOR){
				if ((node[side[s].na].status==INTERIOR) &&
						(node[side[s].nb].status==INTERIOR) &&
						(node[side[s].n1].status==INTERIOR) &&
						(node[side[s].n2].status==INTERIOR) ){
					E =Nne[side[s].n1] + Nne[side[s].n2] - Nne[side[s].na] - Nne[side[s].nb];
					if(E==T) { 
						Nne[side[s].na]++; 
						Nne[side[s].nb]++; 
						Nne[side[s].n1]--; 
						Nne[side[s].n2]--;  
						SwapDiagonal(s);
					}
				}	
			}
		}
	}
	}

	if (mild){numiter=3;}
	else     {numiter=2000;}
	//double L,w,Favg;
  cmplex newz;
	int newe;
	//Smooth mesh
	//should move towards nodes with lower spacing functions
	for(n=0; n<NumNodes; n++){sumz[n]=0.0;sumF[n]=0.0;}// sumW[n]=0.0;}
	for(it=0; it<numiter; it++){
    if (ProgramAborted()){return;}
		for(s=0; s<NumSides; s++){
			//if(side[s].status==INTERIOR){
				//Favg=0.5*(node[side[s].n1].F+node[side[s].n2].F);
				//L=abs(node[side[s].n2].z-node[side[s].n1].z);
				//w=1.0+max(min((L-Favg)/(Favg),0.2),-0.2); //=1 for perfection, >1 <1.2 for error 
				//w=1.0;
				sumF[side[s].n1]+=node[side[s].n2].F; //this is the hard part- re-evaluating the spacing function
				sumF[side[s].n2]+=node[side[s].n1].F;
				sumz[side[s].n1]+=node[side[s].n2].z; //better-handles boundaries
				sumz[side[s].n2]+=node[side[s].n1].z;
				//sumW[side[s].n1]+=w;
				//sumW[side[s].n2]+=w;
			//}
		}
		for(n=0; n<NumNodes; n++){
			if(node[n].status==INTERIOR){
				newz=sumz[n]/(double)(Nne[n]);
				if (mild){ //This is really slow
					newe=GetElement(newz);
					node[n].F=CalculateSpacingFunction(newe,newz);
				}
				else{
					node[n].F=sumF[n]/(double)(Nne[n]);//doesn;t really matter
				}
        node[n].z=newz;
				//node[n].z=sumz[n]/sumW[n];
				//node[n].F=sumF[n]/sumW[n];
			}
      sumz[n]=0.0;
			sumF[n]=0.0;
			//sumW[n]=0.0;
		}
	}

  //correct circumference changes from node movement
	for(e=0; e<NumElems; e++){
		if(elem[e].on==true){CalculateCircles(e);} 
	}
  //delete [] sumW;
}
/**************************************************************************
									Dynamic Fix
***************************************************************************
	Sifts through and looks for really short sides
	Removes bad side and fixes references
	Should be used only for creation of dynamic mesh with boundary spacing:
		not used for simple delunay triangulation of points

         |\											         
				 |  \														 
				 |    \													
				 | ea2  \												 
			 na|________\n2_________		 ___sa2__
		    /|			  /|        /		  |			  /|        
		  /	 | ea	  /  | eb1  /				|   	/  |   
	  /		 |	  /    |    /			 sa1|	  s    |sb1    
	/  ea1 |  /  eb  |  / 	        |  /     |  
/________|/________|/	            |/_______|
				n1\        |nb				        sb2   
					  \  eb2 | 					         
						  \    |            						 
							  \  | 							              
								  \|							            	  
delete s, sa2, sb2
delete ea, eb
delete n2
-------------------------------------------------------------------------*/
void CTriMesh::DynamicFix(){
	//double maxL;
	int n1,n2,na,nb;
	int ea,eb,ea1,ea2,eb1,eb2;
	int s,sa1,sa2,sb1,sb2;

	for(s=0; s<NumSides; s++){

		n1=side[s].n1;
		n2=side[s].n2;
		na=side[s].na;
		nb=side[s].nb;

		//maxL=0.25*max(node[n1].F,node[n2].F);

		if(side[s].status==INTERIOR){
			if (abs(node[n1].z-node[n2].z)<minSpacing){

				cout <<"Bad Side (side "<<s<<" being removed: length:"<<abs(node[n1].z-node[n2].z)<<" min spacing: "<<minSpacing<<endl;
				//return;
				//move nodes
				cmplex newz;
				if      ((node[n1].status!=INTERIOR) && (node[n2].status!=INTERIOR)){
					return; //cant move either node
				}
				else if ((node[n1].status==INTERIOR) && (node[n2].status!=INTERIOR)){
					newz=node[n2].z;
				}
				else if ((node[n2].status==INTERIOR) && (node[n1].status!=INTERIOR)){
					newz=node[n1].z;
				}
				else{
					newz=0.5*(node[n1].z+node[n2].z);
				}
				
				node[n1].z=newz;
				node[n2].z=newz;
		
				//rereference element sides, nodes, & adjacent elements
				//rereference side adjacent elements and nodes
				ea=side[s].ea;
				eb=side[s].eb;

				if (ea!=OFF){//==========================================================================
					//get local information
					if      (elem[ea].ni==na){ea1=elem[ea].ek;ea2=elem[ea].ej;sa1=elem[ea].sk; sa2=elem[ea].sj;}
					else if (elem[ea].nj==na){ea1=elem[ea].ei;ea2=elem[ea].ek;sa1=elem[ea].si; sa2=elem[ea].sk;}
					else if (elem[ea].nk==na){ea1=elem[ea].ej;ea2=elem[ea].ei;sa1=elem[ea].sj; sa2=elem[ea].si;}
					else {ExitGracefully("TriMesh: DynamicFix: bad reference(1)",RUNTIME_ERR);}
					
					cout << "ea1:"<<ea1<<" ea2:"<<ea2<<" sa1:"<<sa1<<" sa2:"<<sa2<<" n1:"<<n1<<" n2:"<<n2<<endl;
					//rereference element adjacent elems, sides, and nodes---------------------------------
					//first element
					if (ea2!=OFF){
						if      (elem[ea2].si==sa2){
							elem[ea2].ei=ea1;
							elem[ea2].si=sa1;
							if      (elem[ea2].nj==n2){elem[ea2].nj=n1;}
							else if (elem[ea2].nk==n2){elem[ea2].nk=n1;}
							else    {ExitGracefully("TriMesh: DynamicFix: bad reference(2)",RUNTIME_ERR);}
						}
						else if (elem[ea2].sj==sa2){
							elem[ea2].ej=ea1;
							elem[ea2].sj=sa1;
							if      (elem[ea2].nk==n2){elem[ea2].nk=n1;}
							else if (elem[ea2].ni==n2){elem[ea2].ni=n1;}
							else    {ExitGracefully("TriMesh: DynamicFix: bad reference(3)",RUNTIME_ERR);}
						}
						else if (elem[ea2].sk==sa2){
							elem[ea2].ek=ea1;
							elem[ea2].sk=sa1;
							if      (elem[ea2].ni==n2){elem[ea2].ni=n1;}
							else if (elem[ea2].nj==n2){elem[ea2].nj=n1;}
							else    {
								cout <<" ni:"<< elem[ea2].ni<< " nj:"<<elem[ea2].nj<<" nk:"<<elem[ea2].nk<<endl;
								ExitGracefully("TriMesh: DynamicFix: bad reference(4)",RUNTIME_ERR);}
						}
						else      {ExitGracefully("TriMesh: DynamicFix: bad reference(5)",RUNTIME_ERR);}
						CalculateCircles(ea2);
					}
					//second element
					if (ea1!=OFF){
						if      (elem[ea1].si==sa1){elem[ea1].ei=ea2;}
						else if (elem[ea1].sj==sa1){elem[ea1].ej=ea2;}
						else if (elem[ea1].sk==sa1){elem[ea1].ek=ea2;}
						else      {ExitGracefully("TriMesh: DynamicFix: bad reference(6)",RUNTIME_ERR);}		
						CalculateCircles(ea1);
					}

					//rereference side nodes and elements--------------------------------------------------
          if      (side[sa1].ea==ea){
						side[sa1].ea=ea2;
						if (ea2==OFF){
							side[sa1].na=OFF;
						}
						else {
							if      (elem[ea2].si=sa1){side[sa1].na=elem[ea2].ni;}
							else if (elem[ea2].sj=sa1){side[sa1].na=elem[ea2].nj;}
							else if (elem[ea2].sk=sa1){side[sa1].na=elem[ea2].nj;}
							else      {ExitGracefully("TriMesh: DynamicFix: bad reference(7)",RUNTIME_ERR);}	
						}
					}
					else if (side[sa1].eb==ea){
						side[sa1].eb=ea2;
						if (ea2==OFF){
							side[sa1].nb=OFF;
						}
						else {
							if      (elem[ea2].si=sa1){side[sa1].nb=elem[ea2].ni;}
							else if (elem[ea2].sj=sa1){side[sa1].nb=elem[ea2].nj;}
							else if (elem[ea2].sk=sa1){side[sa1].nb=elem[ea2].nj;}
							else      {ExitGracefully("TriMesh: DynamicFix: bad reference(8)",RUNTIME_ERR);}	
						}
				
					}
					else        {ExitGracefully("TriMesh: DynamicFix: bad reference(9)",RUNTIME_ERR);}	

					//delete ea,sa2
					elem[ea].on=false;
					side[sa2].status=IS_OFF;
				}


				if (eb!=OFF){//==========================================================================
					//get local information
					if      (elem[eb].ni==nb){eb1=elem[eb].ej;eb2=elem[eb].ek;sb1=elem[eb].sj; sb2=elem[eb].sk;}
					else if (elem[eb].nj==nb){eb1=elem[eb].ek;eb2=elem[eb].ei;sb1=elem[eb].sk; sb2=elem[eb].si;}
					else if (elem[eb].nk==nb){eb1=elem[eb].ei;eb2=elem[eb].ej;sb1=elem[eb].si; sb2=elem[eb].sj;}
					else {ExitGracefully("TriMesh: DynamicFix: bad reference(1b)",RUNTIME_ERR);}
					
					//rereference element adjacent elems, sides, and nodes---------------------------------
					//first element
					if (eb2!=OFF){
						if      (elem[eb2].si==sb2){
							elem[eb2].ei=eb1;
							elem[eb2].si=sb1;
							if      (elem[eb2].nj==n2){elem[eb2].nj=n1;}
							else if (elem[eb2].nk==n2){elem[eb2].nk=n1;}
							else    {ExitGracefully("TriMesh: DynamicFix: bad reference(2b)",RUNTIME_ERR);}
						}
						else if (elem[eb2].sj==sb2){
							elem[eb2].ej=eb1;
							elem[eb2].sj=sb1;
							if      (elem[eb2].nk==n2){elem[eb2].nk=n1;}
							else if (elem[eb2].ni==n2){elem[eb2].ni=n1;}
							else    {ExitGracefully("TriMesh: DynamicFix: bad reference(3b)",RUNTIME_ERR);}
						}
						else if (elem[eb2].sk==sb2){
							elem[eb2].ek=eb1;
							elem[eb2].sk=sb1;
							if      (elem[eb2].ni==n2){elem[eb2].ni=n1;}
							else if (elem[eb2].nj==n2){elem[eb2].nj=n1;}
							else    {ExitGracefully("TriMesh: DynamicFix: bad reference(4b)",RUNTIME_ERR);}
						}
						else      {ExitGracefully("TriMesh: DynamicFix: bad reference(5b)",RUNTIME_ERR);}
						CalculateCircles(eb2);
					}
					//second element
					if (eb1!=OFF){
						if      (elem[eb1].si==sb1){elem[eb1].ei=eb2;}
						else if (elem[eb1].sj==sb1){elem[eb1].ej=eb2;}
						else if (elem[eb1].sk==sb1){elem[eb1].ek=eb2;}
						else      {ExitGracefully("TriMesh: DynamicFix: bad reference(6b)",RUNTIME_ERR);}	
						CalculateCircles(eb1);
					}

					//rereference side nodes and elements--------------------------------------------------
          if      (side[sb1].ea==eb){
						side[sb1].ea=ea2;
						if (eb2==OFF){
							side[sb1].na=OFF;
						}
						else {
							if      (elem[eb2].si=sb1){side[sb1].na=elem[eb2].ni;}
							else if (elem[eb2].sj=sb1){side[sb1].na=elem[eb2].nj;}
							else if (elem[eb2].sk=sb1){side[sb1].na=elem[eb2].nj;}
							else      {ExitGracefully("TriMesh: DynamicFix: bad reference(7b)",RUNTIME_ERR);}	
						}
					}
					else if (side[sb1].eb==eb){
						side[sb1].eb=eb2;
						if (eb2==OFF){
							side[sb1].nb=OFF;
						}
						else {
							if      (elem[eb2].si=sb1){side[sb1].nb=elem[eb2].ni;}
							else if (elem[eb2].sj=sb1){side[sb1].nb=elem[eb2].nj;}
							else if (elem[eb2].sk=sb1){side[sb1].nb=elem[eb2].nj;}
							else      {ExitGracefully("TriMesh: DynamicFix: bad reference(8b)",RUNTIME_ERR);}	
						}
				
					}
					else        {ExitGracefully("TriMesh: DynamicFix: bad reference(9b)",RUNTIME_ERR);}	
					
					//delete eb,sb2
					elem[eb].on=false;
					side[sb2].status=IS_OFF;

				}

				//"delete" n2, s
				side[s].status=IS_OFF;
				node[n2].status=IS_OFF;
			
			}//end if (abs(node[n1].z-node[n2].z)<maxL){
    } //end if (side[s].status==INTERIOR)
	} //end for (s=0...
}
/*************************************************************************
												Double Node filter 2
**************************************************************************
Filters through dynamically generated mesh
Breaks adjacency of elements separated by doublenode boundary
Creates new nodes at doublenode locations
Does not work for complex intersecting doublenode boundaries (acute angles)
-------------------------------------------------------------------------*/
void CTriMesh::DoubleNodeFilter2        (){
	int n,s,news,newn;
	int n1,n2;
	int NewNumSides=NumSides;
	int NewNumNodes=NumNodes;
	int *pairs;
	cmplex move;
	pairs=new int [2*NumNodes];
	for (n=0; n<2*NumNodes;n++){pairs[n]=n;}

	int AdjDnodes(0);

	for (n=0; n<NumNodes;n++){
		if (node[n].status==DOUBLENODE){
			AdjDnodes=0;
			n1=n2=OFF;
			for (s=0;s<NumSides;s++){
				if (side[s].status==DOUBLENODE){
					if ((side[s].n1==n) && (node[side[s].n2].status==DOUBLENODE)){
						AdjDnodes++; 
						if (n1==OFF){n1=side[s].n2;}	else        {n2=side[s].n2;}
					}
					if ((side[s].n2==n) && (node[side[s].n1].status==DOUBLENODE)){
						AdjDnodes++; 
						if (n2==OFF){n2=side[s].n1;}  else        {	n1=side[s].n1;}
					}
				}
			}
			if (AdjDnodes==2){ //if adjacent doublenodes is two 
				if (NewNumNodes>=MaxNodes){ExitGracefully("CTriMesh::DoubleNodeFilter: too many new nodes",BAD_DATA);}
				newn=NewNumNodes;
				if ((n1==OFF) || (n2==OFF)){ExitGracefully("CTriMesh::DoubleNodeFilter: bad adjacency",RUNTIME_ERR);}

				move=node[n].inv*0.01*abs(node[n2].z-node[n1].z);
				
				node[newn].z=node[n].z-move;
				node[n].z   =node[n].z+move;
				node[newn].F=node[n].F; //Not really neccesary
				node[newn].status =DOUBLENODE;
				if (pairs[n]!=n){ExitGracefully("CTriMesh::DoubleNodeFilter: : cannot handle intersection of 3 or more doublenode sides (2)",BAD_DATA);}
				pairs[n]=newn;
				//cout<<"pair: "<<n<<":"<<pairs[n]<<endl;
				NewNumNodes++;
			}
			else if (AdjDnodes>2){
				ExitGracefully("CTriMesh::DoubleNodeFilter: cannot handle intersection of 3 or more doublenode sides",BAD_DATA);
			}
		}
	}

	NumDoublenodes=(NewNumNodes-NumNodes);
	dnodes=new dblestruct[NumDoublenodes];
	int ntmp=0;
	for (n=0; n<NumNodes; n++){
		if (pairs[n]!=n){
			dnodes[ntmp].n1=n;
			dnodes[ntmp].n2=pairs[n];
			dnodes[ntmp].type=DNODE_JUMP; //TMP DEBUG- should have other types here
			ntmp++;
		}
	}
	if (ntmp!=NumDoublenodes){ExitGracefully("CTriMesh::DoubleNodeFilter: problem!",RUNTIME_ERR);}

	for (s=0; s<NumSides;s++){
		if (side[s].status==DOUBLENODE){
			if (NewNumSides>=MaxNodes*3){ExitGracefully("CTriMesh::DoubleNodeFilter: too many new sides",BAD_DATA);}
			
			//create new side
			news=NewNumSides;
			side[news].n1    =pairs[side[s].n1];
			side[news].n2    =pairs[side[s].n2];
			side[news].nb    =side[s].nb;
			side[news].eb    =side[s].eb;
			side[news].na    =OFF;
			side[news].ea    =OFF;
			side[news].status=DOUBLENODE;
			
			//change elements, redirect element nodes
			int ea=side[s].ea;
			int eb=side[s].eb;

			if ((ea!=OFF) && (eb!=OFF)){

				if      (elem[eb].ei==ea){
					elem[eb].ei=OFF;
					elem[eb].si=news;

					//elem[eb].ni=pairs[elem[eb].ni];
					elem[eb].nj=pairs[elem[eb].nj];
					elem[eb].nk=pairs[elem[eb].nk];

					if (side[elem[eb].sj].n1==side[s].n1){side[elem[eb].sj].n1=pairs[side[s].n1];}
					if (side[elem[eb].sj].n1==side[s].n2){side[elem[eb].sj].n1=pairs[side[s].n2];}
					if (side[elem[eb].sj].n2==side[s].n1){side[elem[eb].sj].n2=pairs[side[s].n1];}
					if (side[elem[eb].sj].n2==side[s].n2){side[elem[eb].sj].n2=pairs[side[s].n2];}

					if (side[elem[eb].sk].n1==side[s].n1){side[elem[eb].sk].n1=pairs[side[s].n1];}
					if (side[elem[eb].sk].n1==side[s].n2){side[elem[eb].sk].n1=pairs[side[s].n2];}
					if (side[elem[eb].sk].n2==side[s].n1){side[elem[eb].sk].n2=pairs[side[s].n1];}
					if (side[elem[eb].sk].n2==side[s].n2){side[elem[eb].sk].n2=pairs[side[s].n2];}
				}
				else if (elem[eb].ej==ea){
					elem[eb].ej=OFF;
					elem[eb].sj=news;

					elem[eb].ni=pairs[elem[eb].ni];
					//elem[eb].nj=pairs[elem[eb].nj];
					elem[eb].nk=pairs[elem[eb].nk];

					if (side[elem[eb].sk].n1==side[s].n1){side[elem[eb].sk].n1=pairs[side[s].n1];}
					if (side[elem[eb].sk].n1==side[s].n2){side[elem[eb].sk].n1=pairs[side[s].n2];}
					if (side[elem[eb].sk].n2==side[s].n1){side[elem[eb].sk].n2=pairs[side[s].n1];}
					if (side[elem[eb].sk].n2==side[s].n2){side[elem[eb].sk].n2=pairs[side[s].n2];}

					if (side[elem[eb].si].n1==side[s].n1){side[elem[eb].si].n1=pairs[side[s].n1];}
					if (side[elem[eb].si].n1==side[s].n2){side[elem[eb].si].n1=pairs[side[s].n2];}
					if (side[elem[eb].si].n2==side[s].n1){side[elem[eb].si].n2=pairs[side[s].n1];}
					if (side[elem[eb].si].n2==side[s].n2){side[elem[eb].si].n2=pairs[side[s].n2];}
				}
				else if (elem[eb].ek==ea){
					elem[eb].ek=OFF;
					elem[eb].sk=news;

					elem[eb].ni=pairs[elem[eb].ni];
					elem[eb].nj=pairs[elem[eb].nj];
					//elem[eb].nk=pairs[elem[eb].nk];					

					if (side[elem[eb].si].n1==side[s].n1){side[elem[eb].si].n1=pairs[side[s].n1];}
					if (side[elem[eb].si].n1==side[s].n2){side[elem[eb].si].n1=pairs[side[s].n2];}
					if (side[elem[eb].si].n2==side[s].n1){side[elem[eb].si].n2=pairs[side[s].n1];}
					if (side[elem[eb].si].n2==side[s].n2){side[elem[eb].si].n2=pairs[side[s].n2];}

					if (side[elem[eb].sj].n1==side[s].n1){side[elem[eb].sj].n1=pairs[side[s].n1];}
					if (side[elem[eb].sj].n1==side[s].n2){side[elem[eb].sj].n1=pairs[side[s].n2];}
					if (side[elem[eb].sj].n2==side[s].n1){side[elem[eb].sj].n2=pairs[side[s].n1];}
					if (side[elem[eb].sj].n2==side[s].n2){side[elem[eb].sj].n2=pairs[side[s].n2];}
				}
			}

			//fix node for element not adjacent to side
			//Still problematic - somehow turns sides along doublenode corners off
			for (int e=0; e<NumElems; e++){
				if ((elem[e].ni==side[s].n1) &&
					  (TriArea(node[side[s].n1].z,node[side[s].n2].z,elem[e].zv)<=0.0)){
					//cout << "Corrected added node (i) elem:"<<e<<" node: "<<elem[e].ni<<endl;
			//		elem[e].ni=pairs[elem[e].ni];
				}
				if ((elem[e].nj==side[s].n1)  &&
					  (TriArea(node[side[s].n1].z,node[side[s].n2].z,elem[e].zv)<=0.0)){
					//cout << "Corrected added node (j) elem:"<<e<<" node: "<<elem[e].nj<<endl;
			//		elem[e].nj=pairs[elem[e].nj];
				}
				if ((elem[e].nk==side[s].n1)  &&
					  (TriArea(node[side[s].n1].z,node[side[s].n2].z,elem[e].zv)<=0.0)){
					//cout << "Corrected added node (k) elem:"<<e<<" node: "<<elem[e].nk<<endl;
					
				//	elem[e].nk=pairs[elem[e].nk];
					//if ((elem[e].nk<0) || (elem[e].nk>=NewNumNodes)){ExitGracefully("Doublenodefilter: Bad elem node",RUNTIME_ERR);}
				}
			}

      side[s].eb=OFF;
			side[s].nb=OFF;

			NewNumSides++;
		}
	}

	//adjust elements based upon correct nodes associated with correct sides
	for(int e=0; e<NumElems; e++){
		int si=elem[e].si;
		int sj=elem[e].sj;
		int sk=elem[e].sk;
		int n1,n2;
		if (si!=OFF){
			n1=side[si].n1;
			n2=side[si].n2;
			if      (abs(node[n1].z-node[elem[e].nj].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nj=n1;}
			else if (abs(node[n2].z-node[elem[e].nj].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nj=n2;}
			if      (abs(node[n1].z-node[elem[e].nk].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nk=n1;}
			else if (abs(node[n2].z-node[elem[e].nk].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nk=n2;}
			//if      (abs(node[n1].z-node[elem[e].ni].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].ni=n1;}
			//else if (abs(node[n2].z-node[elem[e].ni].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].ni=n2;}
		}
		if (sj!=OFF){
			n1=side[sj].n1;
			n2=side[sj].n2;
			if      (abs(node[n1].z-node[elem[e].ni].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].ni=n1;}
			else if (abs(node[n2].z-node[elem[e].ni].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].ni=n2;}
			if      (abs(node[n1].z-node[elem[e].nk].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nk=n1;}
			else if (abs(node[n2].z-node[elem[e].nk].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nk=n2;}
			//if      (abs(node[n1].z-node[elem[e].nj].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nj=n1;}
			//else if (abs(node[n2].z-node[elem[e].nj].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nj=n2;}
		}
		if (sk!=OFF){
			n1=side[sk].n1;
			n2=side[sk].n2;
			if      (abs(node[n1].z-node[elem[e].ni].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].ni=n1;}
			else if (abs(node[n2].z-node[elem[e].ni].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].ni=n2;}
			if      (abs(node[n1].z-node[elem[e].nj].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nj=n1;}
			else if (abs(node[n2].z-node[elem[e].nj].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nj=n2;}
			//if      (abs(node[n1].z-node[elem[e].nk].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nk=n1;}
			//else if (abs(node[n2].z-node[elem[e].nk].z)<0.5*abs(node[n1].z-node[n2].z)){elem[e].nk=n2;}
		}
	}

	if (NumDoublenodes>0){
		
	  cout <<"..."<<NewNumSides-NumSides<<" new sides created along doublenode boundaries"<<endl;
		cout <<"..."<<NumDoublenodes      <<" new doubly valued nodes created"<<endl;
	}
	NumSides=NewNumSides;
	NumNodes=NewNumNodes;

	/*for(int e=0; e<NumElems; e++){
		cout <<"dblnode2 " << e << " "<<elem[e].ni<<" " <<elem[e].nj<<" "<<elem[e].nk<<endl;
		if ((elem[e].ni<0) || (elem[e].ni>=NumNodes)){ExitGracefully("Doublenodefilter: Bad elem node i",RUNTIME_ERR);}
		if ((elem[e].nj<0) || (elem[e].nj>=NumNodes)){ExitGracefully("Doublenodefilter: Bad elem node j",RUNTIME_ERR);}
		if ((elem[e].nk<0) || (elem[e].nk>=NumNodes)){ExitGracefully("Doublenodefilter: Bad elem node k",RUNTIME_ERR);}	
	}*/

	delete [] pairs;
}
/*************************************************************************
												Double Node filter
**************************************************************************
Filters through dynamically generated mesh
Breaks adjacency of elements separated by doublenode boundary
Creates new nodes at doublenode locations
-------------------------------------------------------------------------*/
void CTriMesh::DoubleNodeFilter(){
  int s,s2,n1,n2;
  int dn;
	bool gotone;
	
	int estcount=0;
	for (int n=0; n<NumNodes; n++){
		if (node[n].status==DOUBLENODE){estcount++;} //may be a bit high
	}
	dnodes=new dblestruct[estcount];

	NumDoublenodes=0;

	for (s=0; s<NumSides; s++){
		if (side[s].status==DOUBLENODE){
			n1=side[s].n1;
			n2=side[s].n2;
			for (s2=0;s2<NumSides;s2++){
				if ((side[s2].status==DOUBLENODE) && (s!=s2)) {
					gotone=false;
					
					if      (((n1==side[s2].n1) || (n1==side[s2].n2)) && (node[n1].status==DOUBLENODE)){dn=n1;gotone=true;}
					else if (((n2==side[s2].n1) || (n2==side[s2].n2)) && (node[n2].status==DOUBLENODE)){dn=n2;gotone=true;}
					
					if (gotone){
						if      (WingInternal(node[dn].z+node[dn].inv,s,dn,s2)){
							DoubleNodeBreak(s,dn,s2);
							
							dnodes[NumDoublenodes].n1=dn;
							dnodes[NumDoublenodes].n2=NumNodes-1;
							dnodes[NumDoublenodes].type=DNODE_BOUNDARY; //TMP DEBUG
							NumDoublenodes++;
							break;
						}
						else if (WingInternal(node[dn].z+node[dn].inv,s2,dn,s)){
							DoubleNodeBreak(s2,dn,s);
							
							dnodes[NumDoublenodes].n1=dn;
							dnodes[NumDoublenodes].n2=NumNodes-1;
							dnodes[NumDoublenodes].type=DNODE_BOUNDARY; //TMP DEBUG
							NumDoublenodes++;
							break;
						}
						else{
							//actually, wing is between previously split sides
							//ExitGracefully("CTriMesh::DoubleNodeFilter2: wing internal not working",RUNTIME_ERR);
						}
					}
				}
			}
		}
	}

	for (s=0; s<NumSides; s++){
		if (side[s].status==DNODE_TEMP){
			side[s].status=DOUBLENODE;
		}
	}
}
/*************************************************************************
												Double Node Break
**************************************************************************
Breaks mesh at node n along sides s1 and s2
      /|\               //\\
 out/  s2 \  in       /s2s2a\
	/    |    \        /   ||   \
 <-----n----->  ->  <----n|---->
  \    |    /        \   ||   /
	  \  s1 /            \s1s1a/
		  \|/               \\//
creates one new node
creates zero, one, or two new sides
returns number of sides created
-------------------------------------------------------------------------*/
int CTriMesh::DoubleNodeBreak( const int s1,const int n,const int s2){

	bool fix1,fix2; 
	bool b1,b2; //true if side s1 and s2, respectively are "backwards"
	int s1a,s2a,ea,eb,nstart,nend;
	int newn;

	if (NumNodes>=MaxNodes-1)  {ExitGracefully("CTriMesh::DoubleNodeBreak2: Exceeded max nodes",RUNTIME_ERR);}
  if (NumSides>=MaxNodes*3-2){ExitGracefully("CTriMesh::DoubleNodeBreak2: Exceeded max sides",RUNTIME_ERR);}
  if ((side[s1].status!=DOUBLENODE) || 
		  (side[s2].status!=DOUBLENODE) || 
		  (node[ n].status!=DOUBLENODE)){ExitGracefully("CTriMesh::DoubleNodeBreak2: bad input",RUNTIME_ERR);}

	fix1=(!((side[s1].ea==OFF) || (side[s1].eb==OFF))); //if side has already been split, then no splitting
	fix2=(!((side[s2].ea==OFF) || (side[s2].eb==OFF)));

	//create new node
	newn=NumNodes;
	node[newn].z     =node[n].z;
	node[newn].F     =node[n].F;
	node[newn].status=BOUNDARY; //ensures that nodes arent re-broken
	node[newn].inv   =node[n].inv;
	node[n   ].status=BOUNDARY;
	NumNodes=NumNodes+1;

	cout <<"Breaking double node! : "<<fix1<<" "<<fix2<<" "<<n<<" "<<newn<<endl;

	s1a=s1;
	s2a=s2;
	//create 2 new sides-direct copy 
	if (fix1){
		s1a         =NumSides; //newest index
		side[s1a].ea=side[s1].ea;		
		side[s1a].eb=side[s1].eb;		
		side[s1a].na=side[s1].na;		
		side[s1a].nb=side[s1].nb;		
		side[s1a].n1=side[s1].n1;		
		side[s1a].n2=side[s1].n2;		
		side[s1a].A =side[s1].A;		
		side[s1a].B =side[s1].B;	
		side[s1a].radius=-side[s1].radius;		
		side[s1a].status=DOUBLENODE;	
		side[s1 ].status=DNODE_TEMP;	
		NumSides=NumSides+1;
	}
	if (fix2){	
		s2a=NumSides;        //newest index
		side[s2a].ea=side[s2].ea;
		side[s2a].eb=side[s2].eb;
		side[s2a].na=side[s2].na;
		side[s2a].nb=side[s2].nb;
		side[s2a].n1=side[s2].n1;
		side[s2a].n2=side[s2].n2;
		side[s2a].A=side[s2].A;
		side[s2a].B=side[s2].B;
		side[s2a].radius=-side[s2].radius;	
		side[s2a].status=DOUBLENODE;
    side[s2 ].status=DNODE_TEMP;

		NumSides=NumSides+1;
	}
	
	if (side[s1].n1==n){nstart=side[s1].n2;b1=true;}
	else               {nstart=side[s1].n1;b1=false;}
	if (side[s2].n1==n){nend  =side[s2].n2;b2=false;}
	else               {nend  =side[s2].n1;b2=true;}

	//adjust side references
	//----First Side---------------------------------------------------------
	if (fix1){ //first side is not already "split"
		if (!b1){ //side 1 is not backwards (dn==n2)
			cout <<"repair 1a (not backwards)"<<endl;

			ea=side[s1].ea;//ea is non-split side
			eb=side[s1].eb;//eb is split side
   
			if      (elem[ea].si==s1){elem[ea].ei=OFF;}
			else if (elem[ea].sj==s1){elem[ea].ej=OFF;}
			else if (elem[ea].sk==s1){elem[ea].ek=OFF;}

			if      (elem[eb].si==s1){elem[eb].ei=OFF;elem[eb].si=s1a;elem[eb].nj=newn;} 
			else if (elem[eb].sj==s1){elem[eb].ej=OFF;elem[eb].sj=s1a;elem[eb].nk=newn;} 
			else if (elem[eb].sk==s1){elem[eb].ek=OFF;elem[eb].sk=s1a;elem[eb].ni=newn;} 
		
			side[s1 ].eb=OFF;
			side[s1 ].nb=OFF;
			side[s1a].ea=OFF;
			side[s1a].na=OFF;
			side[s1a].n2=newn;
		}
		else if (b1){
			cout <<"repair 1b (backwards)"<<endl;

			ea=side[s1].eb;//ea is non-split side
			eb=side[s1].ea;//eb is split side
   
			if      (elem[ea].si==s1){elem[ea].ei=OFF;}
			else if (elem[ea].sj==s1){elem[ea].ej=OFF;}
			else if (elem[ea].sk==s1){elem[ea].ek=OFF;}

			if      (elem[eb].si==s1){elem[eb].ei=OFF;elem[eb].si=s1a;elem[eb].nj=newn;} 
			else if (elem[eb].sj==s1){elem[eb].ej=OFF;elem[eb].sj=s1a;elem[eb].nk=newn;} 
			else if (elem[eb].sk==s1){elem[eb].ek=OFF;elem[eb].sk=s1a;elem[eb].ni=newn;} 
			
			side[s1 ].ea=OFF;
			side[s1 ].na=OFF;
			side[s1a].eb=OFF;
			side[s1a].nb=OFF;
			side[s1a].n1=newn;

		}
	}

	//----Second Side---------------------------------------------------------
	if (fix2){ //first side is not already "split"
		if (!b2){ //side 2 is not backwards (dn==n1)
			cout <<"repair 2a (not backwards)"<<endl;

			ea=side[s2].ea;//ea is non-split side
			eb=side[s2].eb;//eb is split side
   
			if      (elem[ea].si==s2){elem[ea].ei=OFF;}
			else if (elem[ea].sj==s2){elem[ea].ej=OFF;}
			else if (elem[ea].sk==s2){elem[ea].ek=OFF;}

			if      (elem[eb].si==s2){elem[eb].ei=OFF;elem[eb].si=s2a;elem[eb].nk=newn;} 
			else if (elem[eb].sj==s2){elem[eb].ej=OFF;elem[eb].sj=s2a;elem[eb].ni=newn;} 
			else if (elem[eb].sk==s2){elem[eb].ek=OFF;elem[eb].sk=s2a;elem[eb].nj=newn;} 
		
			side[s2 ].eb=OFF;
			side[s2 ].nb=OFF;
			side[s2a].ea=OFF;
			side[s2a].na=OFF;
			side[s2a].n1=newn;
		}
		else if (b2){
			cout <<"repair 2b (backwards)"<<endl;

			ea=side[s2].eb; //ea is non-split side
			eb=side[s2].ea; //eb is split side
   
			if      (elem[ea].si==s2){elem[ea].ei=OFF;}
			else if (elem[ea].sj==s2){elem[ea].ej=OFF;}
			else if (elem[ea].sk==s2){elem[ea].ek=OFF;}

			if      (elem[eb].si==s2){elem[eb].ei=OFF;elem[eb].si=s2a;elem[eb].nk=newn;} 
			else if (elem[eb].sj==s2){elem[eb].ej=OFF;elem[eb].sj=s2a;elem[eb].ni=newn;} 
			else if (elem[eb].sk==s2){elem[eb].ek=OFF;elem[eb].sk=s2a;elem[eb].nj=newn;}
			
			side[s2 ].ea=OFF;
			side[s2 ].na=OFF;
			side[s2a].eb=OFF;
			side[s2a].nb=OFF;
			side[s2a].n2=newn;

		}
	}
	//adjust element nodal references, side na and sb references
	int si,sj,sk;

	for (int e=0; e<NumElems; e++){

    //Detect all elements "linked" to doublenode (new or old)
		if ((elem[e].ni==n   ) || (elem[e].nj==n   ) || (elem[e].nk==n   ) ||
			  (elem[e].ni==newn) || (elem[e].nj==newn) || (elem[e].nk==newn)){ 
			
			si=elem[e].si;
			sj=elem[e].sj;
			sk=elem[e].sk;
			//check if element is on "adjusted side" of new node
			if (WingInternal(elem[e].zin,s1,n,s2)){
				//if so, modify side and element nodal references  accordingly
				if ((elem[e].ni==n) || (elem[e].ni==newn)){//cout <<"adji: "<<e<<endl;
					if      (side[si].na==n){side[si].na=newn;}
					else if (side[si].nb==n){side[si].nb=newn;}

					if ((sj!=s1) && (sj!=s2)){
						if      (side[sj].n1==n){side[sj].n1=newn;}
						else if (side[sj].n2==n){side[sj].n2=newn;}
					}
					if ((sk!=s1) && (sk!=s2)){
						if      (side[sk].n1==n){side[sk].n1=newn;}
						else if (side[sk].n2==n){side[sk].n2=newn;}
					}
					elem[e].ni=newn;
				}
				else if ((elem[e].nj==n) || (elem[e].nj==newn)){//cout <<"adjj: "<<e<<endl;
					if      (side[sj].na==n){side[sj].na=newn;}
					else if (side[sj].nb==n){side[sj].nb=newn;}

					if ((si!=s1) && (si!=s2)){
						if      (side[si].n1==n){side[si].n1=newn;}
						else if (side[si].n2==n){side[si].n2=newn;}
					}
					if ((sk!=s1) && (sk!=s2)){
						if      (side[sk].n1==n){side[sk].n1=newn;}
						else if (side[sk].n2==n){side[sk].n2=newn;}
					}
					elem[e].nj=newn;
				}
				else if ((elem[e].nk==n) || (elem[e].nk==newn)){//cout <<"adjk: "<<e<<endl;
					if      (side[sk].na==n){side[sk].na=newn;}
					else if (side[sk].nb==n){side[sk].nb=newn;}

					if ((si!=s1) && (si!=s2)){
						if      (side[si].n1==n){side[si].n1=newn;}
						else if (side[si].n2==n){side[si].n2=newn;}
					}
					if ((sj!=s1) && (sj!=s2)){
						if      (side[sj].n1==n){side[sj].n1=newn;}
						else if (side[sj].n2==n){side[sj].n2=newn;}
					}
					elem[e].nk=newn;
				}
				else{
					//cout <<"Weirdness"<<endl;
				}
			}
		}
	}

	if (!fix1){
		//merely move side reference
		if (!b2) {
			side[s2].n1=newn;
		}
		else{
			side[s2].n2=newn;
		}
	}
	if (!fix2){
		//merely move side reference
		if (!b1) {
			side[s1].n2=newn;
		}
		else{
			side[s1].n1=newn;
		}
	}
	//slightly move node
	cmplex move=node[n].inv*0.005*abs(node[nend].z-node[nstart].z);

	node[n].z   -=move;
	node[newn].z+=move;

	int count(0);
  if (fix1){count++;}
	if (fix2){count++;}

	return count;

}
/*************************************************************************
												WingInternal
**************************************************************************
	z+			 /n2
False	  s2/	
			dn /
			  |    z +
		  s1|  TRUE
		    |n1
Works based upon Right hand rule
-------------------------------------------------------------------------*/
bool CTriMesh::WingInternal(const cmplex z,const int s1, const int dn, const int s2) const{

	int n1,n2;
	
	if (s1==s2){
		     ExitGracefully("CTriMesh::WingInternal: bad wing(1)",RUNTIME_ERR);} //bad input

	if      (side[s1].n1==dn){n1=side[s1].n2;}
	else if (side[s1].n2==dn){n1=side[s1].n1;}
	else  {ExitGracefully("CTriMesh::WingInternal: bad wing(2)",RUNTIME_ERR);} //bad input

	if      (side[s2].n1==dn){n2=side[s2].n2;}
	else if (side[s2].n2==dn){n2=side[s2].n1;}
	else  {ExitGracefully("CTriMesh::WingInternal: bad wing(3)",RUNTIME_ERR);} //bad input

	//cout <<z<<" s1:"<<s1<<" s2:"<<s2<<" dn:"<<dn<<" n1:"<<n1<<" n2:"<<n2<<endl;
//	cout <<"side[s1].n1:"<<side[s1].n1<<" side[s2].n1:"<<side[s2].n1<<endl;
//	cout <<"side[s1].n2:"<<side[s1].n2<<" side[s2].n2:"<<side[s2].n2<<endl;

	if (abs(node[n1].z-node[dn].z)==0.0){ExitGracefully("CTriMesh::WingInternal: bad side",RUNTIME_ERR);}
	if (abs(node[n2].z-node[dn].z)==0.0){ExitGracefully("CTriMesh::WingInternal: bad side",RUNTIME_ERR);}

  cmplex Z1=(z-0.5*(node[dn].z+node[n1].z))/(0.5*(node[n1].z-node[dn].z));
	cmplex Z2=(z-0.5*(node[n2].z+node[dn].z))/(0.5*(node[dn].z-node[n2].z));

	if ((Z1.imag()==0.0) || (Z2.imag()==0.0)){
		cout<<Z1<<" "<<Z2<<endl;
		ExitGracefully("CTriMesh::WingInternal: along wing",RUNTIME_ERR);} //bad input


	return ((Z1.imag()>0.0) && (Z2.imag()>0.0));
}
/*************************************************************************
												Dynamic Build Mesh
**************************************************************************
-------------------------------------------------------------------------*/
void CTriMesh::DynamicBuildMesh(bool relax_on, bool smooth_on){

	int count(0);
	bool success(false);
	cout << "Dynamically building triangular mesh...";

	DynamicClassify();
	do{
		if (!DynamicNewNode())    {cout <<"unable to insert new node ("  <<count<<" iterations performed)"<<endl;break;}
		//DynamicFix();
		//DynamicSmooth(false,true);

		if ((NumNodes>=500) && (NumNodes%500==0)) {cout<<NumNodes<<" nodes (of maximum "<<MaxNodes<<")"<<endl;}

		DynamicClassify();
		count++;

		if (ProgramAborted())     {cout <<"stopped by user ("            <<count<<" iterations performed)"<<endl;break;}
		if (NumNodes==MaxNodes-1) {cout <<"exceeded maximum nodes ("     <<count<<" iterations performed)"<<endl;success=true;break;}
		if (ugly==OFF)            {cout <<"success! ("                   <<count<<" iterations performed)"<<endl;success=true;}
	} while (ugly!=OFF);
	int s;
	int curvecount=0;for (s=0;s<NumSides;s++){if (fabs(side[s].radius)>0){curvecount++;}}	cout <<"curvecount"<<curvecount<<endl;

	if (NumNodes==MaxNodes-1){
    cout << "Performing excessive smoothing (this may take a while)";
		//smooth until fixed??
		for (int i=0; i<50; i++){
			cout <<".";
			DynamicSmooth(relax_on,false);
		}
		cout<<endl;
	}

	cout <<"  Improving the grid quality..."<<endl;
	if ((success) && (smooth_on)){DynamicSmooth(relax_on,false);}

	cout <<"  Creating doublenodes..."<<endl;
	if (success){DoubleNodeFilter();}

	cout <<"  Renumerating nodes, elements and sides..."<<endl;
	if (success){Renumerate();}
  
	curvecount=0;for ( s=0;s<NumSides;s++){if (fabs(side[s].radius)>0){curvecount++;}}	cout <<"curvecount"<<curvecount<<endl;

	//Renumerate();//TMP DEBUG- for intermediate images
	
	Process();
	
	if (success){
		cout << "...Triangular Mesh Created:"<<endl;
		cout << "     # of Elems: " <<NumElems<<endl;
		cout << "     # of Sides: " <<NumSides<<endl;
		cout << "     # of Nodes: " <<NumNodes<<endl<<endl; 
	}
}

//*************************************************************************
CTriMesh *CTriMesh::DynamicParse(ifstream &input, int &l,char *Name){
	/*
	TriMesh 
	MaxNodes (or <=zero if default)
	{x y spacing mark} x (NumPoints+NumBoundarySegs) 
	     -mark=0 for interior (e.g. coarsening), 1 for real boundary, 
	&
	{p1 p2 mark {optional radius}} x (NumInternalSegs) 
	     -p1=index of 1st point, p2=index of endpt 
			 -mark=0 for interior, 1 for real boundary, 2 for double noded boundary
			 -radius is positive if arc center is right of seg p1->p2
	&
	*/
	CTriMesh   *pMesh=NULL;
  bool				eof(false),done(false);
	int					npoints(0), nsegs(0);
	int					maxnodes;
	int					Len;
	char       *s[MAXINPUTITEMS];
	
	nodestruct  points  [MAX_DYN_PTS];
	segstruct   segments[MAX_DYN_PTS];


	if (parserdebug) {cout << "Triangular Mesh"<<endl;}  eof=TokenizeLine(input,s,Len); l++; 
	while ((Len==0)&& (!eof)){													 eof=TokenizeLine(input,s,Len); l++;}
	if (Len==1){                               
		maxnodes=s_to_i(s[0]);
		if ((maxnodes<=0) || (maxnodes>MAX_NODES)){maxnodes=MAX_NODES;}
	}
	else       {ImproperFormat(s,l); return NULL;}

	done=false; 																				 eof=TokenizeLine(input,s,Len); l++;
  do {
		if (npoints>=MAX_DYN_PTS) { ExitGracefully("CTriMesh::Parse- too many points associated with mesh",TOO_MANY);}
	  if  (Len==4){
			points[npoints].z=s_to_c(s[0],s[1]);
			points[npoints].F=s_to_d(s[2]);
			points[npoints].status=TranslateMark(s_to_i(s[3])); 
			npoints++;	                                     eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len==1) && (!strcmp(s[0],"&")))    {
			done=true;                                       eof=TokenizeLine(input,s,Len); l++;}
		else if       (Len==0){														 eof=TokenizeLine(input,s,Len); l++;}
    else                  {ImproperFormat(s,l); return NULL;}
	} while ((!done) && (!eof));

	done=false;																					 
  do {
		if (nsegs>=MAX_DYN_PTS) { ExitGracefully("CTriMesh::Parse- too many segments associated with mesh",TOO_MANY);}
	  if ((Len==3) || (Len==4)){
			//checks for repeated segments
			segments[nsegs].p0=s_to_i(s[0]);
			segments[nsegs].p1=s_to_i(s[1]); 
			segments[nsegs].status=TranslateMark(s_to_i(s[2])); 
			segments[nsegs].radius=0;
			if (Len==4){segments[nsegs].radius=s_to_d(s[3]);}
			
			if (segments[nsegs].p0>=npoints || segments[nsegs].p1>=npoints){ 
				 ExitGracefully("CTriMesh::Parse- segment point index greater than number of points ",BAD_DATA);}
      bool repeated=false;
			for (int i=0;i<nsegs; i++){
				if (((segments[i].p0==segments[nsegs].p0) &&
					   (segments[i].p1==segments[nsegs].p1)) || 
						((segments[i].p0==segments[nsegs].p1) &&
					   (segments[i].p1==segments[nsegs].p0))){
					repeated=true;
					//cout <<"CTriMesh::Parse: REPEATED SEGMENT!"<<endl;
				}
			}
			if (!repeated){nsegs++;}	                       eof=TokenizeLine(input,s,Len); l++;
		}
    else if ((Len==1) && (!strcmp(s[0],"&")))         {
			spcstruct spacing_package;
      spacing_package.type = INTERP;
			pMesh=new CTriMesh(Name,points,npoints,segments,nsegs,maxnodes,spacing_package);
			done=true; 
		}
		else if       (Len==0){														 eof=TokenizeLine(input,s,Len); l++;}
    else                  {ImproperFormat(s,l); return NULL;}
	} while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pMesh;}
}

/*double CalcCourantSpacingFunction(const cmplex &z) const{
	
	cmplex v;
	double t(0.0);

	v=CSingleLayer::GetVelocity2D(z,t);
}*/
