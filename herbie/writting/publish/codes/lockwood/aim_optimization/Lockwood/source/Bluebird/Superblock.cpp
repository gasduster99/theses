#include "Superblock.h"

//************************************************************************
//                           CSuperblock
//													 CONSTRUCTORS
//************************************************************************
CSuperblock::CSuperblock(){
	order=0; nblockcontrol=0;size=0;
	pChild[0]=NULL; pChild[1]=NULL;  pChild[2]=NULL; pChild[3]=NULL;}
//------------------------------------------------------------------------
CSuperblock::CSuperblock(const cmplex zc,const double Width, 
						             const int nestlevel,const int quad){
  int    m,n;
	double angle;	
  
	if ((order<0) || (nblockcontrol<0) || (nblockcontrol>MAXBLOCKCONTROL)){
		ExitGracefully("CSuperblock::Constructor- bad input",BAD_DATA);
	}
  
	zcen=    zc;         
	width=   Width;    
	nestlev= nestlevel;
  Q=       0.0;
	leaf=    false;      
	master=  false;   
	quadrant=quad;

	pElemArray     =NULL; size=0; population=0;
	pCondZoneArray =NULL; nCondZones=0;
	pBaseZoneArray =NULL; nBaseZones=0; 
	pThickZoneArray=NULL; nThickZones=0;
	pPoroZoneArray =NULL; nPoroZones=0;

  ElemQ         =NULL;         
	ElemRHS       =NULL;
	ElemSegment   =NULL;
	pParent       =NULL;
	for (int k=0; k<4; k++){pChild[k]=NULL;}

  radius=BlockRadius*width/sqrt(2.0); 

  //allocate memory for dynamic arrays (rest allocated in initialize)
  LaurentCoeff= new cmplex [order];  
  zctrl=        new cmplex [nblockcontrol]; 

  //fill array of control points, initialize all coeff to 0
  for (m=0; m<nblockcontrol; m++){                   
		angle=2.0*PI*(double)(m)/(double)(nblockcontrol);
    zctrl[m]= radius*cmplex(cos(angle), sin(angle))+zcen;
	}

	for (n=0; n<order; n++){LaurentCoeff[n]=0.0;} 
  if  (nestlev==0)       {master=true;}
}
//------------------------------------------------------------------------
CSuperblock::~CSuperblock(){
	int j;
	if (globaldebug){
		for(j=0; j<nestlev; j++){cout<<" ";}
		cout <<"  DESTROYING SUPERBLOCK LEV:"<<nestlev<<" QUAD:"<<quadrant<<endl;
	}

  //recursive destruction of all child blocks -only called once on master from layer
	for (int k=0; k<4; k++){
		if (pChild[k]!=NULL){delete pChild[k];}
	}
	if (nblockcontrol>0){
		delete [] pElemArray;       //delete pointers (elements deleted by layer)
		delete [] pCondZoneArray;   //delete pointers (prop zones deleted by layer)
		delete [] pBaseZoneArray;   //delete pointers (prop zones deleted by layer) 
	  delete [] pThickZoneArray;  //delete pointers (prop zones deleted by layer) 
	  delete [] pPoroZoneArray;   //delete pointers (prop zones deleted by layer) 
  	delete [] zctrl;
	  delete [] LaurentCoeff;
 		for (int i=0; i<size; i++){
			delete [] ElemRHS[i];
		}
		delete [] ElemRHS;
		delete [] ElemQ;
		delete [] ElemSegment;	
	}
}        
//***********************************************************************
//											ACCESSOR FUNCTIONS
//***********************************************************************
cmplex           CSuperblock::Centroid() const           {return zcen;}
//-----------------------------------------------------------------------
int			         CSuperblock::GetQuadrant() const		     {return quadrant;}
//-----------------------------------------------------------------------
double           CSuperblock::GetRadius() const          {return radius;}
//-----------------------------------------------------------------------
double           CSuperblock::GetWidth() const           {return width;}
//-----------------------------------------------------------------------
cmplex           CSuperblock::GetCenter() const          {return zcen;}
//-----------------------------------------------------------------------
int              CSuperblock::GetLevel() const					 {return nestlev;}
//-----------------------------------------------------------------------
bool             CSuperblock::IsLeaf() const             {return leaf;}
//-----------------------------------------------------------------------
int 	           CSuperblock::GetPopulation() const      {return population;}
//-----------------------------------------------------------------------
bool             CSuperblock::IsOn() const               {return On;}
//-----------------------------------------------------------------------
bool             CSuperblock::SuperblocksOn()            {return On;}
//-----------------------------------------------------------------------
int	             CSuperblock::CalcPopulation()           { //recursive
	int k,pop=0.0;
	if (!leaf){for (k=0;k<4; k++){pop+=pChild[k]->CalcPopulation();}}
	population=pop + size;
  return population;
}
/***********************************************************************
											Update 
************************************************************************
equivalent to solveitself for superblocks
an element notifies the superblock that its coefficients have changed by calling this function
The superblock updates its Laurent series and net flux accordingly
The superblock then requests its parent block to change
-----------------------------------------------------------------------*/
void CSuperblock::Update(const int BlockID,const int segment,const double t){
  int              m,n;
	double           angle,IncQ;
	cmplex           Omega,sum,Z;
  cmplex          *IncLaurent,*b;
  
	// if superblocks are not nested and this is not a real superblock, then do nothing
  //------------------------------------------------------------------------------
	if (!Nested && !leaf) {return;}

  // if the block is empty, set coefficients to zero and do nothing else
  //------------------------------------------------------------------------------
	if (population==0) {for (n=0; n<order; n++){LaurentCoeff[n]=0.0;} Q=0.0; return;}

	IncLaurent=new cmplex[order];
	b=				 new cmplex[order];

  //save old laurent, Q, centroid for incrementing parent
  //------------------------------------------------------------------------------
  for (n=0; n<order; n++){IncLaurent[n]=LaurentCoeff[n];}   IncQ=Q; 

	//subtract old influence of element from Q & Laurent Series coefficients 
  //------------------------------------------------------------------------------
	//discharge
	Q-=ElemQ[BlockID];

	//1st coeff
	for (m=0; m<nblockcontrol; m++){
    LaurentCoeff[0]-=ElemRHS[BlockID][m]/(double)(nblockcontrol);
	}

	//the rest
	for (n=1; n<order; n++){
		sum=0.0;
		for (m=0; m<nblockcontrol; m++){
			angle=2.0*PI*(double)(m)/(double)(nblockcontrol);  //should calc only once??
			sum+=ElemRHS[BlockID][m]*cmplex(cos((double)(n)*angle),sin((double)(n)*angle)); 
    }
	  LaurentCoeff[n]-=2.0*sum/(double)(nblockcontrol);
	}

	// obtain new potential & Discharge from element inside the block at outer radius (update ElemRHS & ElemQ & ElemCentroid)
  //------------------------------------------------------------------------------
	for (m=0; m<nblockcontrol; m++){
		if (segment<0) {ElemRHS[BlockID][m]=(pElemArray[BlockID]->GetDischargePotential      (zctrl[m],t)).real();}
		else           {ElemRHS[BlockID][m]=(pElemArray[BlockID]->GetSegmentPotential(segment,zctrl[m],t)).real();}	
	}
  if (segment<0)	 {ElemQ       [BlockID]=pElemArray[BlockID]->GetNetDischarge(t);           }
  else             {ElemQ       [BlockID]=pElemArray[BlockID]->GetSegmentDischarge(segment,t);}


  // get the new Q, Laurent Series coefficients for this block
  //------------------------------------------------------------------------------
	//discharge
  Q+=ElemQ[BlockID];

	//1st Coeff
	for (m=0; m<nblockcontrol; m++){
    LaurentCoeff[0]+=ElemRHS[BlockID][m]/(double)(nblockcontrol); 
	}

	//the rest
	for (n=1; n<order; n++){
		sum=0.0;
		for (m=0; m<nblockcontrol; m++){
			angle=2.0*PI*(double)(m)/(double)(nblockcontrol); 
			sum+=ElemRHS[BlockID][m]*cmplex(cos((double)(n)*angle),sin((double)(n)*angle)); 
		}
		LaurentCoeff[n]+=2.0*sum/(double)(nblockcontrol);
	}

	//save incremental change for parent (inc=new-old)
	//------------------------------------------------------------------------------
	for (n=0; n<order; n++){IncLaurent[n]=(LaurentCoeff[n]-IncLaurent[n]);}
	IncQ=(Q-IncQ); 
 
	//Update the parent block(s)
	//------------------------------------------------------------------------------
	if(!master){
		if (!Explicit){              pParent->UpdateParent(IncLaurent,IncQ,quadrant);}
		else          {double ch,obj;pParent->SetLaurent(ch,obj,t);                  }
	}

	// compute the objective function
	//------------------------------------------------------------------------------
	/*double objective=0.0;
  double change=0.0;
  double error,pot; 
  for (m=0; m<nblockcontrol; m++){
    pot=0.0;
    Z=(zctrl[m]-zcen)/(radius);
		for (int i=0;i<size; i++){pot+=ElemRHS[i][m];}
		if(!leaf){
			for (int k=0; k<4; k++){
			  pot+=(pChild[k]->GetDischargePotential(zctrl[m],t)).real();
			}
		}
		error = pot-(Laurent(Z,order,LaurentCoeff).real()) - 0.5*Q/pi*log(Z).real();	
		//cout <<error<<endl;
		if (fabs(error)>change)    {change   =fabs(error);}
    if (fabs(error)>objective) {objective=fabs(error);}	
	}
	cout <<"chi_obj:"<<objective <<endl;*/

	/*ofstream STREAMMATCH;
	STREAMMATCH.open("streammatch.csv");
	for (n=0; n<order; n++){
		STREAMMATCH<<LaurentCoeff[n].real()<<","<<LaurentCoeff[n].imag()<<endl;
	}
	STREAMMATCH.close();*/

	delete [] IncLaurent;
	delete [] b;
}
/************************************************************************
											UpdateParent
*************************************************************************
Update parent block's Laurent expansion based upon change in child
Based upon Randal's trick
-----------------------------------------------------------------------*/
void CSuperblock::UpdateParent(const cmplex *CoeffIncrement,const double QIncrement, const int quad){
	  int     n;
		cmplex  sum;
		cmplex *b= new cmplex[order];

		//randall's method
		for (n=0; n<order; n++){
			b[n]=0.0;
      for (int n2=0; n2<order; n2++){b[n]+=G[quad][n][n2]*CoeffIncrement[n2];}	 //G * a //correct
			b[n]+=f[quad][n]*QIncrement;																						   //b=f*Q + G*a	
		}                                        
    //b[0]+=QIncrement/pi*log(2.0);																								 //b0=f*Q +G*a +Q/pi*ln(2)
		
		//update laurent coeff, netQ
	  for (n=0; n<order; n++){LaurentCoeff[n]+=b[n];} 		Q+=QIncrement;

		//send incremental change up to parent
    if (!master){	pParent->UpdateParent(b,QIncrement,quadrant);}

		//Evaluate objective function-------------------------------------------------------
		/*double objective=0.0;
		double change=0.0;
		double error,pot,out; 
		int m;
		double  t(0.0);
		for ( m=0; m<nblockcontrol; m++){
			cmplex Z=(zctrl[m]-zcen)/(radius);
			pot=(GetDischargePotential(zctrl[m]-(0.000001*Z),t)).real(); //slightly inside
			out=Laurent(Z,order,LaurentCoeff).real()+0.5*Q/pi*log(Z).real();
			error = pot- out;	
			if (fabs(error)>change)    {change   =fabs(error);}
			if (fabs(error)>objective) {objective=fabs(error);}	
		}
		cout << "par_obj "<<objective<<endl;*/

    delete [] b; 
 
}
/************************************************************************
											SetLaurent
*************************************************************************
Explicitly calculates Laurent Series coefficients based on update of child
-----------------------------------------------------------------------*/
void CSuperblock::SetLaurent(double &change, double &objective,double t){
	int    i,k,m,n;
	static double angle,Potential[MAXBLOCKCONTROL];
	static cmplex Omega,sum,Z;

  // if the block is empty, set coefficients to zero and do nothing else (should never be used)
	if (population==0) {for (n=0; n<order; n++){LaurentCoeff[n]=0.0;}Q=0; return;}

  //update net discharge from block
	Q=0.0;
	for (i=0; i<size; i++)  {Q+=ElemQ[i];}
	if (!leaf) {
		for (k=0; k<4; k++)   {Q+=pChild[k]->GetNetDischarge(t);}
	}
	// obtain potential from all elems inside the block at outer radius
	for (m=0; m<nblockcontrol; m++){
		Potential[m]=0.0;
		for (i=0; i<size; i++){Potential[m]+=ElemRHS[i][m];}
		if (!leaf) {
		  for (k=0; k<4; k++) {Potential[m]+=pChild[k]->GetDischargePotential(zctrl[m],t).real();}
		}
	}

  // get the Laurent Series coefficients for this block
	LaurentCoeff[0]=0.0;
	for (m=0; m<nblockcontrol; m++){
		LaurentCoeff[0]+=cmplex(Potential[m]/(double)(nblockcontrol),0.0);
	}

	for (n=1; n<order; n++){
		sum=0.0;
		for (m=0; m<nblockcontrol; m++){
			angle=2.0*PI*(double)(m)/(double)(nblockcontrol); 
			sum+=Potential[m]*cmplex(cos((double)(n)*angle),sin((double)(n)*angle));
		}
		
	  LaurentCoeff[n]=2.0*sum/(double)(nblockcontrol);
	}

	if (!master){pParent->SetLaurent(change,objective,t);}


   //=======================================================================
   // compute the objective function
  /*objective=change=0.0;
	double error;
  for (m=0; m<nblockcontrol; m++){
    Z=(zctrl[m]-zcen)/(radius);
    error = Potential[m]-Laurent(Z,order,LaurentCoeff).real() -0.5*Q/pi*log(Z).real();
	  if (fabs(error)>change)    {change   =fabs(error);}
    if (fabs(error)>objective) {objective=fabs(error);}
	}
	cout <<"par_obj:"<< objective <<"   Q: "<<Q<<endl; */
}
//***********************************************************************
//											GetDischargePotential
//***********************************************************************

cmplex CSuperblock::GetDischargePotential(const cmplex &z,const double &t) const{

	cmplex Omega(0.0);
	cmplex Z((z-zcen)/(radius));

  if ((abs(Z)>=1.0) && (Nested || leaf)){			          //outside block radius- Laurent	+ Well
		Omega=Laurent(Z,(int)(order*pow(abs(Z),SB_FF_POWER)),LaurentCoeff)+0.5*Q/PI*log(Z); //laurent
		Omega=Laurent(Z,order,LaurentCoeff)+0.5*Q/PI*log(Z); //laurent
	}
  else {
		if (!leaf){                                         //get potential from children
			for (int k=0; k<4; k++) {Omega+=    pChild[k]->GetDischargePotential(z,t);              }
		}                                                   //get potential from local elements
    for (int i=0; i<size; i++){
			if (ElemSegment[i]<0)   {Omega+=pElemArray[i]->GetDischargePotential(z,t);              }
			else                    {Omega+=pElemArray[i]->GetSegmentPotential(ElemSegment[i],z,t);}
		} 
	}

  return Omega;
}

/***********************************************************************
											GetW
************************************************************************
----------------------------------------------------------------------*/
cmplex CSuperblock::GetW(const cmplex &z,const double &t) const{

  cmplex W(0.0);
	cmplex Z((z-zcen)/radius);

  if ((abs(Z)>=1.0) && (Nested || leaf)){						    //perform laurent expansion
		W=LaurentDer(Z,(int)(order*pow(abs(Z),SB_FF_POWER)),LaurentCoeff,radius)-0.5*Q/(PI*(z-zcen));
	}
  else {
    if (!leaf){																					//get complex W from children
		  for (int k=0; k<4; k++){W+=    pChild[k]->GetW(z,t);                       }
		}																										//get complex W from local elements
		for (int i=0; i<size; i++){
			if (ElemSegment[i]<0)  {W+=pElemArray[i]->GetW(z,t);                       }
			else                   {W+=pElemArray[i]->GetSegmentW(ElemSegment[i],z,t);}
		}
	}

  return W;
}
/***********************************************************************
											GetGx
************************************************************************
----------------------------------------------------------------------*/
cmplex CSuperblock::GetGx(const cmplex &z,const double &t) const{

  cmplex Gx(0.0);
	cmplex Z((z-zcen)/radius);

  if ((abs(Z)>=1.0) && (Nested || leaf)){						    //perform laurent expansion
		Gx=-LaurentDxx(Z,(int)(order*pow(abs(Z),SB_FF_POWER)),LaurentCoeff,radius)+0.5*Q/(PI*(z-zcen)/(z-zcen));
	}
  else {
    if (!leaf){																					//get Gx from children
		  for (int k=0; k<4; k++){Gx+=    pChild[k]->GetGx(z,t);                       }
		}																										//get Gx from local elements
		for (int i=0; i<size; i++){
			if (ElemSegment[i]<0)  {Gx+=pElemArray[i]->GetGx(z,t);                       }
			else                   {Gx+=pElemArray[i]->GetSegmentGx(ElemSegment[i],z,t);}
		}
	}

  return Gx;
}
/***********************************************************************
											GetLeakage
************************************************************************
----------------------------------------------------------------------*/
double CSuperblock::GetLeakage (const cmplex &z, const double &t,const leak_type ltype) const{
  double Leak(0.0);
	cmplex Z((z-zcen)/radius);
  if (abs(Z)>1.0){																			//
    Leak=0.0;
	}
  else {
    if (!leaf){																					//get leakage from children
		  for (int k=0; k<4; k++){Leak+=  pChild[k]->GetLeakage(z,t,ltype);}
		}																										//get leakage from local elements
		for (int i=0; i<size; i++){
			if (ElemSegment[i]<0){ Leak+=pElemArray[i]->GetLeakage(z,t,ltype);										  }
			else                 { Leak+=pElemArray[i]->GetSegmentLeakage(ElemSegment[i],z,t,ltype);}
		}
	}

  return Leak;

}
/************************************************************************
											GetNetDischarge
************************************************************************/
double CSuperblock::GetNetDischarge(const double &t) const{return Q;}
/************************************************************************
											GetCond
************************************************************************/
double CSuperblock::GetCond(const cmplex &z) const{
	double Kloc(NO_VALUE);	
  if (!leaf){                             //get conductivity from children
		for (int k=0; k<4; k++){ 
			if (pChild[k]->IsInside(z)){
			  Kloc=pChild[k]->GetCond(z);
			  if (Kloc !=NO_VALUE) {return Kloc;}    //smallest nested inhoms will be in children 
			}
		}
	}                                       //get conductvity from members
	return CPropZone::NestedSift(z,pCondZoneArray,nCondZones,NO_VALUE);
}
/************************************************************************
											GetBase
************************************************************************/
double CSuperblock::GetBase(const cmplex &z) const{
	double Bloc(NO_VALUE);	
  if (!leaf){                             //get base elev from children
		for (int k=0; k<4; k++){ 
			Bloc=pChild[k]->GetBase(z);
			if (Bloc != NO_VALUE) {return Bloc;} 
		}
	}                                       //get base elev from members
	return CPropZone::NestedSift(z,pBaseZoneArray,nBaseZones,NO_VALUE);
}
/************************************************************************
											GetThick
************************************************************************/
double CSuperblock::GetThick(const cmplex &z) const{
	double Tloc(NO_VALUE);	
  if (!leaf){                             //get base elev from children
		for (int k=0; k<4; k++){ 
			Tloc=pChild[k]->GetThick(z);
			if (Tloc != NO_VALUE) {return Tloc;} 
		}
	}                                       //get base elev from members
	return CPropZone::NestedSift(z,pThickZoneArray,nThickZones,NO_VALUE);
}
/************************************************************************
											GetPoro
************************************************************************/
double CSuperblock::GetPoro(const cmplex &z) const{
	double nloc(NO_VALUE);	
  if (!leaf){                             //get base elev from children
		for (int k=0; k<4; k++){ 
			nloc=pChild[k]->GetPoro(z);
			if (nloc != NO_VALUE) {return nloc;} 
		}
	}                                       //get base elev from members
	return CPropZone::NestedSift(z,pPoroZoneArray,nPoroZones,NO_VALUE);
}
/************************************************************************
											GEOMETRIC FUNCTIONALITY
************************************************************************/
bool CSuperblock::IsInside(const cmplex &z) const{
	if ((z.real()>(zcen.real()-width/2.0))&&
		  (z.real()<(zcen.real()+width/2.0))&&
		  (z.imag()>(zcen.imag()-width/2.0))&&
			(z.imag()<(zcen.imag()+width/2.0))) {return true;}
	else                                    {return false;}
}


//void CSuperblock 
/***********************************************************************
											SetTaylor
************************************************************************/
/*void CSuperblock::SetTaylor(const double &t){
	int     i,n,m;
	double *Potential, angle;
	cmplex  Omega,sum;

	Potential =new double [nblockcontrol];

  if (leaf){
	  // obtain potential from all elems outside the block at outer radius
		for (m=0; m<nblockcontrol; m++){
			Omega=pLayer->GetDischargePotential(zctrl[m],t);  //Block cannot have knowledge of layer
			for (i=0; i<sizeINT; i++){
				Omega-=pElemArrayINT[i]->GetDischargePotential(zctrl[m],t);
			}
		}
		Potential[m]=Omega.real();
		pTaylorCoeff[0]+=cmplex(Potential[m]/(double)(nblockcontrol),0.0);
    
		// get the Taylor Series coefficients for this block
		for (n=1; n<order; n++){
			sum=cmplex(0.0,0.0);
			for (m=0; m<nblockcontrol; m++){
				angle=2.0*pi*(double)(m)/(double)(nblockcontrol); 
				sum+=cmplex( Potential[m]*cos((double)(n)*angle),
					          -Potential[m]*sin((double)(n)*angle)); 
			}
			pTaylorCoeff[n]=-conj(2.0*sum/(double)(nblockcontrol));
		}
	}

	delete [] Potential;
}*/
//TEMP DEBUG-Taylor series expansion for solved mode
/*if (!solved){
	if (!leaf){
		for (int k=0; k<4; k++){
			if (pChild->IsInside(z)) { return pChild[k]->GetDischargePotential(z,t);}
		}
	}
	else {
		for (int i=0; i<sizeINT; i++){
		  if (ElemSegmentINT[i]<0)  {Omega+=pElemArrayINT[i]->GetDischargePotential(z,t);                 }
			else                       {Omega+=pElemArrayINT[i]->GetSegmentPotential(ElemSegmentINT[i],z,t);}
		}
	}
	Omega+=Inside(Z,order,pTaylorCoeff);
	return Omega;
}*/
//TEMP DEBUG-Taylor series expansion for solved mode
/*if (!solved){
	if (!leaf){
		for (int k=0; k<4; k++){
			if (pChild->IsInside(z)) { return pChild->GetW(z,t);}
		}
	}
	else {
		for (int i=0; i<sizeINT; i++){
		  if (ElemSegmentINT[i]<0)  {Omega+=pElemArrayINT[i]->GetW(z,t);                 }
			else                       {Omega+=pElemArrayINT[i]->GetSegmentW(ElemSegmentINT[i],z,t);}
		}
	}
	Omega+=InsideW(Z,order,pTaylorCoeff);
	return Omega;
}*/
		//add influence of well at new location
		/*for (n=1; n<order; n++){
			sum=0.0;
			for (m=0; m<nblockcontrol; m++){
				angle=2.0*pi*(double)(m)/(double)(nblockcontrol);  //should calc only once??
				wellinfluence[m]=0.5*Q/pi*log((zctrl[m]-zcentroid)/abs(pLayer->GetBlackHole()-zcentroid)).real(); //Block cannot have knowledge of layer
				sum+=(wellinfluence[m])*cmplex(cos((double)(n)*angle),sin((double)(n)*angle)); 
			}
			LaurentCoeff[n]+=-conj(2.0*sum/(double)(nblockcontrol));
		}*/