#include "Superblock.h"

//***********************************************************************
//											FillSuperBlocks (called only once)
//***********************************************************************
CSuperblock *CSuperblock::FillSuperblocks(window          Extents,  const int    nestlevels,  
																					CAnalyticElem **ElemArray,const int    numElems,
																					CPropZone     **PropArray,const int    numprops, 
																					const int       prec,     const bool   implicit, 
																					const bool      nested){
	int             i,j;
	int             count(0);    //holds count of blocks 
	int             pastcount(0);//holds index in Block[] of first block on this level
	int             startind(0); //holds index in Block[] of first block on last level
	int             row(0);      //row on current level
	int             col(0);      //column on current level
	int             levelstats[10];//holds counters of added element segments
	int             addedtolevel;
	cmplex          zctmp;       //center location of current block
	double          wtmp;        //width of blocks at current level
	double          totwidth;    //total width of extents
  cmplex          topleft;     //top left of domain extents
	CSuperblock    *MasterBlock; //pointer to master block
	double          fold;

	for (i=0; i<=nestlevels; i++){levelstats[i]=0;}

	CSuperblock::Nested  = nested;
	CSuperblock::Explicit=!implicit;
	SetPrecision(prec,CSuperblock::order,fold);
	CSuperblock::nblockcontrol=(int)(fold*order);

	//check for any bad input
	if ((nestlevels<0)     || (nestlevels>5) || 
			(ElemArray==NULL)  || (numElems<0)   || 
			((PropArray==NULL)  && (numprops>0)) || 
			(order<0)          || (fold<1.0)     ||
	    (Extents.n<=Extents.s)               ||
			(Extents.e<=Extents.w)) {
		cout << "nesting:  " <<nestlevels<<"|array: "<<ElemArray<<"|numelems: "<<numElems<<"|order:"<<order<<"|fold:"<<fold<<endl;
		cout << "Extents   " <<Extents.n <<"|"<<Extents.s<<"|"<<Extents.e<<"|"<<Extents.w<<endl;
		cout << "Proparray:" <<PropArray <<"|"<<numprops <<endl;
		ExitGracefully("CSuperblock::FillSuperblocks Bad input",BAD_DATA);
	}

	int totalblocks=0;
	for (i=0;i<nestlevels; i++){
		totalblocks+=ipow(2,2*i);
	}

	CSuperblock **Block =new CSuperblock * [totalblocks];


  ofstream BLOCKS;
  BLOCKS.open("Superblocks.bna");

	//take domain extents, widen slightly (as safety precaution), set domain width & position
  Extents.n+=BIGGERIZE*(Extents.n-Extents.s);
  Extents.s-=BIGGERIZE*(Extents.n-Extents.s)/(1+BIGGERIZE);
  Extents.e+=BIGGERIZE*(Extents.e-Extents.w);
  Extents.w-=BIGGERIZE*(Extents.e-Extents.w)/(1+BIGGERIZE);
  totwidth=max(Extents.n-Extents.s,Extents.e-Extents.w);
  topleft=cmplex (Extents.w, Extents.n);

	//prepare matrices for implicit solution
  Prepare();

	//create block hierarchy
  //--------------------------------------------------------------------------------
	if (globaldebug){cout<<"Creating Nested Superblock Hierarchy"<<endl;}
	for (i=nestlevels; i>=0; i--){																									//for each level,going up from leaf level
		
		wtmp=totwidth/(ipow(2,i));																											//set width of blocks at this level

		for (j=0; j< ipow(4,i); j++){																									//for each block in level, (j counts blocks)
																																									//set zcenter
			zctmp=topleft+cmplex(wtmp/2.0, -wtmp/2.0)+
			      wtmp*cmplex((double)((j)%ipow(2,i)),-(double)(j/ipow(2,i)) );

			Block[count]=new CSuperblock(zctmp,wtmp,i,count-pastcount);                //create block
			if (Block[count]==NULL){ExitGracefully("CSuperblock::FillSuperblocks: out of memory",OUT_OF_MEMORY);}

			if (i!=nestlevels){                                                        //if not leaf, add children from level below
        Block[count]->AddChild(Block[startind+(col+1)+( row   *ipow(2,i+1))],0); //ne child
				Block[count]->AddChild(Block[startind+ col   +( row   *ipow(2,i+1))],1); //nw	child
        Block[count]->AddChild(Block[startind+ col   +((row+1)*ipow(2,i+1))],2); //sw child
        Block[count]->AddChild(Block[startind+(col+1)+((row+1)*ipow(2,i+1))],3); //se child
		    col+=2;
				if ((j+1)%(ipow(2,i))==0){row+=2;col=0;}
			}
			else if (i==nestlevels){Block[count]->SetAsLeaf();}
			count++;
		}
		startind=pastcount;
    pastcount=count;
		row=0;
		col=0;
	}

	MasterBlock=Block[count-1];

	//fill blocks with elements 
	//--------------------------------------------------------------------------------
	if (globaldebug){cout<<"Filling Superblocks"<<endl;}
	for (i=0; i<numElems; i++){
		if (!ElemArray[i]->IsString()){                //put entire element in block
			addedtolevel=MasterBlock->AssignToBlock(ElemArray[i],NOT_A_STRING);
			levelstats[addedtolevel]++;
		}
		else {                                         //put only enabled line segments in blocks
			for (j=0;j<ElemArray[i]->GetNumSegs(); j++){
				if (!ElemArray[i]->IsDisabled(j)){
					addedtolevel=MasterBlock->AssignToBlock(ElemArray[i],j);
					levelstats[addedtolevel]++;
				}
			}
		}
		WriteEllipsisToScreen(i,numElems,10);
	}
  for (i=0; i<numprops; i++){
    MasterBlock->AssignToBlock(PropArray[i]);
	}

	//Add master block to layer, turn on block solve
	CSuperblock::On=true;
	MasterBlock->CalcPopulation(); //recursively calculate populations

	//create dynamic arrays of rhs for each block, write to output
  for (j=0; j<count; j++){
		Block[j]->InitializeBlock();
		if ((Block[j]->GetPopulation()!=0) && (Nested || Block[j]->IsLeaf())) {
			zctmp=Block[j]->GetCenter();
			wtmp=Block[j]->GetWidth();
			i=Block[j]->GetLevel();
			

			BLOCKS << "\" Level "<<i<<" block "<<j<<" \", 2"<<endl;											//write geometry to output						
			BLOCKS << zctmp.real()<<" "<< zctmp.imag()<<endl;
			BLOCKS << BlockRadius*sqrt(2.0*pow(wtmp/2.0,2.0))<<" "<< 0<< endl;

			BLOCKS << "\" Level "<<i<<" block "<<j<<" \", -5"<<endl; 
			BLOCKS << zctmp.real()+wtmp/2<<" "<< zctmp.imag()+wtmp/2<<endl;
			BLOCKS << zctmp.real()+wtmp/2<<" "<< zctmp.imag()-wtmp/2<<endl;
			BLOCKS << zctmp.real()-wtmp/2<<" "<< zctmp.imag()-wtmp/2<<endl;
			BLOCKS << zctmp.real()-wtmp/2<<" "<< zctmp.imag()+wtmp/2<<endl;
			BLOCKS << zctmp.real()+wtmp/2<<" "<< zctmp.imag()+wtmp/2<<endl;
		}
	}



	BLOCKS.close();

	cout << MasterBlock->GetPopulation() <<" elements/element segments added to blocks" <<endl<<endl; 
	for (i=0; i<=nestlevels; i++){
		cout << levelstats[i] << " element/element segments at level " << i << endl;
	}
	delete [] Block;

	return MasterBlock;
}

//***********************************************************************
//											AssignToBlock (recursive) 
//***********************************************************************
int CSuperblock::AssignToBlock(CAnalyticElem *Element,const int seg){
	//ownership is based on the element (or string element segment) 
	//	1) being unable to fit in child 
	//	2) being located entirely within superblock radius
	//	3) having a centroid within the square block bounds

  bool fits_in_child=false;
	int k;
	int addedtolevel;
  if (!leaf){
	  for (k=0; k<4; k++){
			if (seg==NOT_A_STRING){  //for non-string
			  if (!fits_in_child && 
						Element->IsInCircle(pChild[k]->Centroid(),BlockBuffer*pChild[k]->GetRadius()) &&
						InSquare(Element->Centroid(),pChild[k]->Centroid(),pChild[k]->GetWidth())){
				  fits_in_child=true;
				  addedtolevel=pChild[k]->AssignToBlock(Element,seg);
				}
			}
			else{                    //for string segment
 		  	if (!fits_in_child && 
						Element->SegIsInCircle(pChild[k]->Centroid(),BlockBuffer*pChild[k]->GetRadius(),seg) &&
						InSquare(Element->SegCentroid(seg),pChild[k]->Centroid(),pChild[k]->GetWidth())){
			  	fits_in_child=true;
			  	addedtolevel=pChild[k]->AssignToBlock(Element,seg);
				}
			}
		}
	}

	if (leaf || !fits_in_child) {
		this->AddToBlock(Element,seg);
		Element->SetBlockOwner(this,seg,size-1);
		//cout << "  level: "     << nestlev           <<" quad: "  << quadrant 
	  //		 <<" Added element "<<Element->GetName() <<" Segment "<< seg      << endl;
		return nestlev;
	}
	else {return addedtolevel;}
}

//***********************************************************************
//											AssignToBlock (property zone) (recursive) 
//***********************************************************************
void CSuperblock::AssignToBlock(CPropZone *PZ){
	//inhom ownership is based on the property   
	//	1) being unable to fit in child 
	//	2) being located entirely within superblock square

  bool fits_in_child=false;
	int k;
  if (!leaf){
	  for (k=0; k<4; k++){
			if (!fits_in_child && 
				PZ->IsInSquare(pChild[k]->Centroid(),pChild[k]->GetWidth())){  //james' method
				fits_in_child=true;
				pChild[k]->AssignToBlock(PZ);
			}
		}
	}

	if (leaf || !fits_in_child) {
		this->AddToBlock(PZ);
		//cout << "  level: "     << nestlev           <<" quad: "  << quadrant 
		//		 <<" Added Prop Zone with v= "<<PZ->GetValue()     << endl;
	}
}
//***********************************************************************
//											AddToBlock
//************************************************************************
void CSuperblock::AddToBlock(CAnalyticElem *Elemptr, const int seg){
  CAnalyticElem **ptrArray; 
	int            *segArray;
  size=size+1;                                           //increment size
  ptrArray = new CAnalyticElem *[size];									 //allocate memory
	segArray = new int            [size];
  for (int i=0; i<(size-1); i++){
		ptrArray[i]=pElemArray  [i];
		segArray[i]=ElemSegment [i];
	}                                                      //copy array
  ptrArray[size-1]=Elemptr;                              //add new pointer
  segArray[size-1]=seg;
	delete [] ElemSegment;
  delete [] pElemArray;                                  //delete old array
	ElemSegment=segArray;
  pElemArray=  ptrArray;                                 //redirect pointer
} 

//***********************************************************************
//											AddToBlock (Property Zone)
//************************************************************************
void CSuperblock::AddToBlock(CPropZone *PZ){
	CPropZone **ptrArray=NULL;
	int    & thissize=nCondZones;
	void **& pZoneArray=(void**&)(pCondZoneArray);

  if (PZ==NULL){
		ExitGracefully("CSuperblock::AddToBlock- adding NULL property zone",BAD_DATA);
	}
  if (PZ->GetType()==hydraulic_conductivity){
	  thissize=nCondZones;
	  pZoneArray=(void**&)(pCondZoneArray);
	}
	else if (PZ->GetType()==base_elevation){
	  thissize=nBaseZones;
	  pZoneArray=(void**&)(pBaseZoneArray);
	}
	else if (PZ->GetType()==poro){
	  thissize=nPoroZones;
	  pZoneArray=(void**&)(pPoroZoneArray);
	}
	else if (PZ->GetType()==layer_thickness){
	  thissize=nThickZones;
	  pZoneArray=(void**&)(pThickZoneArray);
	}
	else {ExitGracefully("CSuperblock::Unused property type in superblock",OTHER);}

	if (!DynArrayAppend(pZoneArray,(void*)(PZ),thissize)){
		ExitGracefully("CSuperblock::AddToLayer(PZ): adding NULL property zone",BAD_DATA);};
} 
//***********************************************************************
//											InitializeBlock
//***********************************************************************
void CSuperblock::InitializeBlock(){
  int i,m;
	//create dynamic arrays of rhs
  ElemRHS=new double *[size];
	ElemQ  =new double  [size];
  Q=0.0;
	for (i=0; i<size; i++){
    ElemRHS[i]=new double [nblockcontrol];
		ElemQ[i]=0.0;
    for (m=0; m<nblockcontrol; m++){ElemRHS[i][m]=0.0;}
	}
}
//***********************************************************************
//											AddChild,AddParent, SetAsLeaf
//***********************************************************************
void CSuperblock::AddChild (CSuperblock *Child,int index)  {pChild[index]=Child;Child->AddParent(this,index);}
//-----------------------------------------------------------------------
void CSuperblock::AddParent(CSuperblock *Parent, int quad) {if (!master) {pParent=Parent; quadrant=quad;}}
//-----------------------------------------------------------------------
void CSuperblock::SetAsLeaf()                              {leaf=true;}

//***********************************************************************
//				              STATIC MEMBER INITIALIZATION
//***********************************************************************
cmplex **  CSuperblock::f=NULL;
cmplex *** CSuperblock::G=NULL;
//-----------------------------------------------------------------------  
bool       CSuperblock::On=false;
//-----------------------------------------------------------------------  
bool       CSuperblock::Explicit=false;
//----------------------------------------------------------------------- 
bool       CSuperblock::Nested=true;
//-----------------------------------------------------------------------  
int        CSuperblock::nblockcontrol=0;
//-----------------------------------------------------------------------
int        CSuperblock::order=0; 
//-----------------------------------------------------------------------
double     CSuperblock::BlockRadius=1.1;
//-----------------------------------------------------------------------
double     CSuperblock::BlockBuffer=0.85;   
//***********************************************************************
//											Prepare (called once)
//***********************************************************************
void CSuperblock::Prepare(){
  double   angle,childangle;
	cmplex   zctrl,zchild,sum,temp;
	int      k,m,n,n2;
  cmplex **C;
	cmplex  *d[4];
	cmplex **E[4];

	if (globaldebug){cout<<"Preparing Superblock Matrices order:"<<order<< " nblockcontrol:"<<nblockcontrol<<endl;}
	
	if ((iabs(nblockcontrol)>300) || (iabs(order) > 100)){
		ExitGracefully("CSuperblock::Prepare-Matrices too large",TOO_MANY);}	
	if ((nblockcontrol<=0       ) || (order<=0         )){
		ExitGracefully("CSuperblock::Prepare-Bad Input"         ,BAD_DATA);}

	//create, initialize conjugate C matrix
  //------------------------------------------------
	if (globaldebug){cout<<"	C Matrix"<<endl;}
	C=new cmplex *[order+1]; //C is actually conj(C) from writeup
	for (n=0;n<=order; n++){
		C[n]=new cmplex [nblockcontrol];
	}

	for (n=0;n<=order; n++){
		for (m=0; m<nblockcontrol; m++){
      if (n==0){C[n][m]=1.0;}
			else     {
			  angle=-2.0*PI*(double)(m)/(double)(nblockcontrol); //actually -angle
			  C[n][m]=conj(cmplex(cos((double)(n)*angle),sin(double(n)*angle)));
			}
		}
	}

	//Create, initialize d vector,E matrix
  //------------------------------------------------
	if (globaldebug){cout<<"	d vector, E Matrix"<<endl;}
	for (k=0;k<4;k++){
		d[k]=new cmplex  [nblockcontrol];
		E[k]=new cmplex *[order+1];
		for (n=0;n<=order; n++){E[k][n]=new cmplex [nblockcontrol];}
	}

	//based upon a block with diagonal=2, R=BlockRadius
	for (k=0; k<4; k++){
		childangle=0.5*PI*(double)(k)+0.25*PI;    //(0(ne):pi/4,1(nw):3pi/4,2(sw):5pi/4,s(se):7pi/4) childradius=0.5
		zchild=0.5*cmplex(cos(childangle),sin(childangle)); 
    for (m=0; m<nblockcontrol; m++){
			angle=2.0*PI*(double)(m)/(double)(nblockcontrol);
			zctrl=BlockRadius*cmplex(cos(angle), sin(angle)); //controlpts of parent
			d[k][m]=(0.5/PI)*log((2.0*(zctrl-zchild)/zctrl));
			E[k][0][m]=cmplex(1.0,0.0);
			temp=(zctrl-zchild)/(0.5*BlockRadius);
			for (n=1;n<=order; n++){E[k][n][m]=E[k][n-1][m]/temp;} //should reverse order [m],then m
		}
	}

	//Create, initialize f vector, G matrix
  //------------------------------------------------	
	if (globaldebug){cout<<"	f Vector G Matrix"<<endl;}
	f=new cmplex *[4];
  G=new cmplex**[4];
	for (k=0;k<4;k++){
		f[k]=new cmplex  [order+1];
	  G[k]=new cmplex *[order+1];
		for (n=0;n<=order; n++){G[k][n]=new cmplex [order+1];}
	}
  for (k=0;k<4;k++){
    for (n=0; n<=order; n++){
			f[k][n]=0.0;
			//matrix vector multiply (4x)
      for (m=0;m<nblockcontrol;m++){
				f[k][n]+=C[n][m]*d[k][m]/(double)(nblockcontrol); 
			}
      //matrix matrix multiply (4x)
			for (n2=0; n2<=order;n2++){		
        sum=cmplex(0.0,0.0);
				for (m=0; m<nblockcontrol; m++){
					sum+=C[n][m]*E[k][n2][m];
				}																						
				G[k][n][n2]=sum/(double)(nblockcontrol);
			}
		}
	}


/*
	ofstream DEBUG;
	DEBUG.open("debugsuperblockmatrix.txt");
  DEBUG << "C_Matrix(real)------------------------------------------------------------"<<endl<<endl;
	for (n=0; n<=order; n++){for (m=0; m<nblockcontrol; m++){DEBUG << C[n][m].real()<< "    ";}DEBUG<<endl;}
  DEBUG << "C_Matrix(imag)------------------------------------------------------------"<<endl<<endl;
	for (n=0; n<=order; n++){for (m=0; m<nblockcontrol; m++){DEBUG << C[n][m].imag()<< "    ";}DEBUG<<endl;}

	DEBUG << "d_vectors(real)-------------------------------------------------------"<<endl<<endl;
  for (m=0; m<nblockcontrol; m++){DEBUG << d[0][m].real()<< "    ";}DEBUG <<endl;
	for (m=0; m<nblockcontrol; m++){DEBUG << d[1][m].real()<< "    ";}DEBUG <<endl;
	for (m=0; m<nblockcontrol; m++){DEBUG << d[2][m].real()<< "    ";}DEBUG <<endl;
  for (m=0; m<nblockcontrol; m++){DEBUG << d[3][m].real()<< "    ";}DEBUG <<endl;
  DEBUG<<endl<<endl;
	DEBUG << "d_vectors(imag)-------------------------------------------------------"<<endl<<endl;
  for (m=0; m<nblockcontrol; m++){DEBUG << d[0][m].imag()<< "    ";}DEBUG <<endl;
	for (m=0; m<nblockcontrol; m++){DEBUG << d[1][m].imag()<< "    ";}DEBUG <<endl;
	for (m=0; m<nblockcontrol; m++){DEBUG << d[2][m].imag()<< "    ";}DEBUG <<endl;
  for (m=0; m<nblockcontrol; m++){DEBUG << d[3][m].imag()<< "    ";}DEBUG <<endl;
  DEBUG<<endl<<endl;
	DEBUG << "E_Matrix(k=0)(real)------------------------------------------------------------"<<endl<<endl;
  for (n=0; n<=order; n++){for (m=0; m<nblockcontrol; m++){DEBUG << E[0][n][m].real()<< "    ";}DEBUG<<endl;}
	DEBUG << "E_Matrix(k=0)(imag)------------------------------------------------------------"<<endl<<endl;
  for (n=0; n<=order; n++){for (m=0; m<nblockcontrol; m++){DEBUG << E[0][n][m].imag()<< "    ";}DEBUG<<endl;}
  DEBUG << "f_vectors(real)-------------------------------------------------------"<<endl<<endl;
  for (n=0; n<=order; n++){DEBUG << f[0][n].real()<< "    ";}DEBUG <<endl;
	for (n=0; n<=order; n++){DEBUG << f[1][n].real()<< "    ";}DEBUG <<endl;
	for (n=0; n<=order; n++){DEBUG << f[2][n].real()<< "    ";}DEBUG <<endl;
	for (n=0; n<=order; n++){DEBUG << f[3][n].real()<< "    ";}DEBUG <<endl;
  DEBUG<<endl<<endl;
  DEBUG << "f_vectors(imag)-------------------------------------------------------"<<endl<<endl;
  for (n=0; n<=order; n++){DEBUG << f[0][n].imag()<< "    ";}DEBUG <<endl;
	for (n=0; n<=order; n++){DEBUG << f[1][n].imag()<< "    ";}DEBUG <<endl;
	for (n=0; n<=order; n++){DEBUG << f[2][n].imag()<< "    ";}DEBUG <<endl;
	for (n=0; n<=order; n++){DEBUG << f[3][n].imag()<< "    ";}DEBUG <<endl;
  DEBUG<<endl<<endl;
  DEBUG << "G_Matrix(k=0)(real & imag)------------------------------------------------------------"<<endl<<endl;
  for (n=0; n<=order; n++){for (n2=0; n2<=order;n2++){	DEBUG << G[0][n][n2].real()<< "    "<<G[0][n][n2].imag()<< "    ";;}DEBUG<<endl;}
  DEBUG.close();
*/
	if (globaldebug){cout<<"	deleting d,e,& C"<<endl;}
	//delete matrices
  for (k=0; k<4; k++){
		for (n=0;n<=order; n++){delete [] E[k][n];} 
		delete [] d[k];
		delete [] E[k];
	}

  for (n=0; n<=order; n++){delete [] C[n];} 
	delete [] C;
	if (globaldebug){cout<<"	done preparing matrices"<<endl;}
}
//***********************************************************************
void CSuperblock::Destroy(){
	if (globaldebug){cout <<"DESTROYING SUPERBLOCK MATRICES"<<endl;}
	for (int k=0;k<4;k++){
		for (int n=0;n<=order; n++){
			delete [] G[k][n];
		}
		delete [] f[k];
	  delete [] G[k];

	}	
	delete [] f;
  delete [] G;
}
//***********************************************************************
void CSuperblock::SetPrecision(const int Precision,int &order, double &fold){
	if ((Precision<0) || ((Precision>5) && (Precision!=9))){
		ExitGracefully("SetPrecision::Improper precision level specified",BAD_DATA);}
	switch(Precision){
		case(0): {order=5;  fold=1.0; break;}
		case(1): {order=10; fold=1.4; break;}
		case(2): {order=20; fold=1.5; break;}
		case(3): {order=30; fold=1.5; break;}
	  case(4): {order=40; fold=1.6; break;}
	  case(5): {order=60; fold=2.0; break;}
	  case(9): {order=100;fold=3.0; break;}
	}
}