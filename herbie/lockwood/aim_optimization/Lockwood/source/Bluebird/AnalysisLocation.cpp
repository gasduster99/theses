//AnalysisLocation.cpp

#include "AnalysisLocation.h"
/*******************************************************************************
			              CANALYSISLOCATION CLASS:
**********************************************************************
                     CONSTRUCTORS
*********************************************************************/
CAnalysisLocation::CAnalysisLocation(){
	if (!DynArrayAppend((void**&)(pLocations),(void*)(this),nLocations)){
		ExitGracefully("CAnalysisLocation::Constructor: creating NULL analysis location",BAD_DATA);};
}
//--------------------------------------------------------------------
CAnalysisLocation::~CAnalysisLocation(){
}
//--------------------------------------------------------------------
int                 CAnalysisLocation::nLocations=0;
CAnalysisLocation **CAnalysisLocation::pLocations=NULL; 
//--------------------------------------------------------------------
void CAnalysisLocation::CalculateAndWriteAll(const double &t){
	cout <<"Calculating/writing "<<nLocations<<" output locations..."<<endl;
	for (int i=0; i<nLocations; i++){
		pLocations[i]->CalculateStatistics(t);
		pLocations[i]->WriteOutput(t);
	}
}
//--------------------------------------------------------------------
void CAnalysisLocation::DestroyAllAnalysisLocations(){
	for (int i=0; i<nLocations; i++){delete pLocations[i];} delete [] pLocations;
	nLocations=0;
}

/*******************************************************************************
			              TRANSECT CLASS:
**********************************************************************
                     CONSTRUCTORS
*********************************************************************/
CTransect::CTransect():CAnalysisLocation(){			
	z1=z2=0.0;
	influx=0.0;
	outflux=0.0;
	nx=0;
	pLayer=NULL;
	x =NULL;
	Qn=NULL;
	Qt=NULL;
	ID=0;
}
//---------------------------------------------------------------------
CTransect::CTransect(const int              my_ID,
										 const CSingleLayerABC *pLay,				//Layer
										 const cmplex						zstart,			//starting point of line
										 const cmplex						zend,				//ending point of line
										 const int							numdiv) 		//number of integration divisions
										:CAnalysisLocation(){
	double L;

	ID=my_ID;
	ExitGracefullyIf(pLay==NULL,"CTransect: NULL layer",BAD_DATA);
  ExitGracefullyIf(numdiv<=0,"CTransect: negative or zero control points",BAD_DATA);
	
	pLayer=pLay;
	z1=zstart; z2=zend;
	nx=numdiv+1;
	influx=outflux=0.0;
	L=abs(z2-z1);	

	x= new double [nx+1];
	Qn=new double [nx+1];
	Qt=new double [nx+1];
	ExitGracefullyIf(Qt==NULL,"CTransect: constructor: out of memory",OUT_OF_MEMORY);

	//initialize divisions
	x[0]=0;
	for (int i=1; i<=nx; i++){ //may be distributed differently, if desired
		x[i]=x[i-1]+(L/nx); Qn[i]=0; Qt[i]=0;
	}
}
//---------------------------------------------------------------------
CTransect::~CTransect(){
	delete [] x;
	delete [] Qn;
	delete [] Qt;
}
//*********************************************************************
//                     ACCESSORS
//*********************************************************************
double CTransect::GetInflux() const{return influx;}
//---------------------------------------------------------------------
double CTransect::GetOutflux() const{return outflux;}
//---------------------------------------------------------------------
double CTransect::GetTotalFlux() const{return influx-outflux;}//outflux is negative

/*********************************************************************
                     MEMBER FUNCTIONS
*********************************************************************/

void CTransect::CalculateStatistics(const double &t){

	//right hand rule- if flow moving from left of line to right of line, flux "in"
	cmplex  W;
	int     i; 
	double  tmpin,tmpout;

	for (i=0; i<=nx; i++){
		W=pLayer->GetW((x[i]-x[0])/(x[nx]-x[0])*(z2-z1)+z1,t);
		Qn[i]=(conj(W)* IM*conj(z2-z1)/abs(z2-z1)).real(); //normal flux 
		Qt[i]=(conj(W)*    conj(z2-z1)/abs(z2-z1)).real(); //tangential flux 
	}
	//integrate normal flux alont front
	influx=outflux=0.0;
	for (i=0; i<nx; i++){
		IntegrateLine(x[i],x[i+1],Qn[i],Qn[i+1],tmpin,tmpout);
		influx +=tmpin;
		outflux-=tmpout;
		//cout <<"tmpin ("<< i <<") x: "<<x[i]<< " "<<x[i+1]<<" Qn:"<<Qn[i] <<" "<<Qn[i+1]<<" "<<tmpin<<endl;
	}
  //cout <<"Transect: i " << influx << "  o:"<<outflux<<"   a "<<pLayer->GetFluxThruFace(z1,z2,0.0).real();
}
//--------------------------------------------------------------------
/*********************************************************************
									Flux Through Segment
**********************************************************************
calculates the positive area (influx) and negative area (outflux) under 
the normal flux function (Qn) between the points z1, and z2, which MUST lie along the transect
returns two POSITIVE VALUES tmpinflux & tmpoutflux
--------------------------------------------------------------------*/
void CTransect::GetFluxThroughSegment(const cmplex &zint1, 
																	    const cmplex &zint2,
																	          double &tmpinflux, 
																	          double &tmpoutflux){
	

	double xint1,xint2,Qint1,Qint2,tmpin,tmpout;

	tmpinflux=0.0;
	tmpoutflux=0.0;

	if (fabs(abs(zint1-z1)+abs(zint1-z2)-abs(z1-z2))>0.005*abs(z1-z2))
	{
		cout <<fabs(abs(zint1-z1)+abs(zint1-z2)-abs(z1-z2))<<endl;
		cout <<"pt1 "<< z1 <<" "<<z2 <<" "<< zint1<<" "<<zint2<<endl;
		ExitGracefully("CTransect::FluxThroughSegment: point not on transect",RUNTIME_ERR);
	}

	if (fabs(abs(zint2-z1)+abs(zint2-z2)-abs(z1-z2))>0.005*abs(z1-z2))
	{
		cout <<fabs(abs(zint2-z1)+abs(zint2-z2)-abs(z1-z2))<<endl;
		cout <<"pt2 "<< z1 <<" "<<z2 <<" "<< zint1<<" "<<zint2<<endl;
		ExitGracefully("CTransect::FluxThroughSegment: point not on transect",RUNTIME_ERR);
	}

	xint1=abs(zint1-z1);
	xint2=abs(zint2-z1);

	if      (xint1> xint2){double tmpx=xint1; xint1=xint2; xint2=tmpx;} //switch points
	else if (xint1==xint2){return;}

	//cout <<"transect %:"<<(xint2-xint1)/(x[nx])<<endl;

	for (int i=0; i<nx; i++){
		if ((xint1>x[i+1]) || (xint2<x[i])){}
		else if ((xint1>x[i]) && (xint1<x[i+1]))
		{
			if ((xint2>x[i]) && (xint2<x[i+1]))
			{  //special case distance between points smaller than dx
				Qint1=Qn[i]+(xint1-x[i])/(x[i+1]-x[i])*(Qn[i+1]-Qn[i]);
				Qint2=Qn[i]+(xint2-x[i])/(x[i+1]-x[i])*(Qn[i+1]-Qn[i]);
				IntegrateLine(xint1,xint2,Qint1,Qint2,tmpin,tmpout);
				tmpinflux +=tmpin;
				tmpoutflux-=tmpout;
			}
			else
			{
				Qint1=Qn[i]+(xint1-x[i])/(x[i+1]-x[i])*(Qn[i+1]-Qn[i]);
				IntegrateLine(xint1,x[i+1],Qint1,Qn[i+1],tmpin,tmpout);
				tmpinflux +=tmpin;
				tmpoutflux-=tmpout;			
			}
		}
		else if ((xint2>x[i]) && (xint2<x[i+1]))
		{
			Qint2=Qn[i]+(xint2-x[i])/(x[i+1]-x[i])*(Qn[i+1]-Qn[i]);
			IntegrateLine(x[i],xint2,Qn[i],Qint2,tmpin,tmpout);
			tmpinflux +=tmpin;
			tmpoutflux-=tmpout;
		}
		else //segment fully contained by xint1 and xint2
		{ 
			IntegrateLine(x[i],x[i+1],Qn[i],Qn[i+1],tmpin,tmpout);
			tmpinflux +=tmpin;
			tmpoutflux-=tmpout;
		}
	}
	ExitGracefullyIf((tmpinflux <0),"CTransect::GetFluxThroughSegment: bad flux sign(1)",RUNTIME_ERR);//TMP DEBUG
	ExitGracefullyIf((tmpoutflux<0),"CTransect::GetFluxThroughSegment: bad flux sign(2)",RUNTIME_ERR);//TMP DEBUG
	//cout <<" transect in: "<<tmpinflux<<" transect out: "<<tmpoutflux<<" net: "<<tmpinflux-tmpoutflux<<endl;
}

//--------------------------------------------------------------------
double CTransect::CalculateMassFlux(const double &t, double (*conc)(const cmplex &z, const double &t)){
	//calculates total advective mass flux (QC) through transect
	double massflux(0.0);
	double poro;
	cmplex z;

	for (int i=0; i<nx; i++)
	{
		z=(x[i]-x[0])/(x[nx]-x[0])*(z2-z1)+z1;
		poro=pLayer->GetPoro(z);
		massflux+= Qn[i]*fabs(x[i+1]-x[i])*conc(z,t)/poro; //Flux per unit width * unit width * mass per unit volume =mass per unit time
	}
	return massflux;
}
/*********************************************************************
									SUBDIVIDE
**********************************************************************
identifies the positive area under the normal flux function (Qn)
divides this area into equal subareas with flux==fluxperdiv (end-of-line subareas may be less than this)
returns centroids of each equal subarea

input:  fluxperdiv (how to subdivide)
				arrays zc[] and fluxes[] are precreated with size ns(as input)
output: ns (number of divisions), zc (subdivision centroid locations) and fluxes (flux of that subdivision)
--------------------------------------------------------------------*/
void CTransect::Subdivide(cmplex *zc,Writeable1DArray fluxes, int &ns, const double fluxperdiv){
	
	//DOES NOT WORK IF END OF LINE IS NEGATIVE!

	int    count(0), i(0);
	double Qacc(0.0),Qrem(0);				//accumulated flux, flux remaining in segment
	double xs(0),xc(0);							//x location of streamline, x of center of segment
	double xa(x [0]),xb(x [1]);
	double qa(Qn[0]),qb(Qn[1]);
  int    maxtubes=ns;

	//TMP DEBUG-----------------------------------
	//int track(0);
	//double breaks[1000]; //for visual output
	//double xstr  [1000]; 
	//int    bcount(0);    
	//ofstream TMPINJ;
	//TMPINJ.open("tmpInj.txt",ios::app);
	//TMPINJ << "x Qn"<<endl; for (i=0; i<=nx;i++){TMPINJ<<x[i]<<" " <<Qn[i]<<endl;}
	//TMPINJ << "step i xa xb qa qb qrem qacc qi qi+1"<<endl;
	//TMP DEBUG-----------------------------------

	if (fluxperdiv<0) {ExitGracefully("Transect::Subdivide-negative flux per division",OTHER);}

	i=0;
	while ((i<nx) && (count<maxtubes)){
	  if (count>=ns) {ExitGracefully("Transect::Subdivide- Too many subdivisions",OTHER);}
		//if (bcount>=999){ExitGracefully("Transect::Subdivide- bcount exceeds 1000",other);}//TMP DEBUG
		//TMPINJ <<track<<" "<< i<<" "<<xa<<" "<<xb<<" "<<qa<<" "<<qb<<" "<<Qrem<< " "<<Qacc<<" "<<Qn[i]<<" "<<Qn[i+1]<<endl;track++;//TMP DEBUG

		//special cases:
		if ((qa<=0) && (qb<=0) && (i!=nx-1)){						//no effect, negative flux
			i++; xa=x[i];xb=x[i+1];qa=Qn[i];qb=Qn[i+1];
		}
		if ((qa<0) && (qb>0)){													//start new streamtube
			xa=-qa*(xb-xa)/(qb-qa)+xa;
			xs=0.66666*(xb-xa)+xa;
			Qacc=0;
			//breaks[bcount]=xa; bcount++;//TMP DEBUG
		}
		else if ((qa>0) && (qb<0)){											//end of the line- finish off streamtube
			xb=qa*(xb-xa)/(qa-qb)+xa;
			qb=0;
		}

		Qrem=(qa+qb)*(xb-xa)/2.0;												//calculate remaining flux in segment

		if ((Qrem+Qacc)>fluxperdiv){										//new tube done in this segment
			xb=(fluxperdiv-Qacc)/(Qrem)*(xb-xa)+xa;
			qb=(fluxperdiv-Qacc)/(Qrem)*(qb-qa)+qa;
			if (qa>=qb){xc=( (qb*(xb+xa))+( (qa-qb)*(0.33333*(xb-xa)+xa) ) )/(qa+qb);}
			else       {xc=( (qa*(xb+xa))+( (qb-qa)*(0.66667*(xb-xa)+xa) ) )/(qa+qb);}
			xs=((xs*Qacc)+(xc*(fluxperdiv-Qacc)))/fluxperdiv;
			zc    [count]=(xs-x[0])/(x[nx]-x[0])*(z2-z1)+z1;			//add streamline @ xs
			fluxes[count]=fluxperdiv;
			count++;
			Qacc=0; xa=xb; xb=x[i+1]; qa=qb; qb=Qn[i+1];
			//xstr  [count-1]=xs; breaks[bcount]=xa; bcount++;//TMP DEBUG
		}
		else if (((qb<=0) && (qa>0)) || ((i==nx-1) && ((Qacc+Qrem)>0)) ){			//finish off tube
			if (qa>=qb){xc=( (qb*(xb+xa))+( (qa-qb)*(0.33333*(xb-xa)+xa) ) )/(qa+qb);}
			else       {xc=( (qa*(xb+xa))+( (qb-qa)*(0.66667*(xb-xa)+xa) ) )/(qa+qb);}
			xs=((xs*Qacc)+(xc*Qrem))/(Qacc+Qrem);
			zc    [count]=(xs-x[0])/(x[nx]-x[0])*(z2-z1)+z1;			//add streamline @ xs
			fluxes[count]=Qacc+Qrem;
			count++;
			Qacc=0; xa=xb; xb=x[i+1]; qa=qb; qb=Qn[i+1];
			if (i==nx-1){i++;}
			//xstr[count-1]=xs; breaks[bcount]=xa; bcount++;//TMP DEBUG
		}
		else {                                          //keep going to next division
			if (qa>=qb){xc=( (qb*(xb+xa))+( (qa-qb)*(0.33333*(xb-xa)+xa) ) )/(qa+qb);}
			else       {xc=( (qa*(xb+xa))+( (qb-qa)*(0.66667*(xb-xa)+xa) ) )/(qa+qb);}
			xs=((xs*Qacc)+(xc*Qrem))/(Qacc+Qrem);
			Qacc+=Qrem;
			i++; xa=x[i];xb=x[i+1];qa=Qn[i];qb=Qn[i+1];
		}
	} 
	//if (count==maxtubes) {ExitGracefully("Transect::PlaceStreamlines- Too many streamtubes",other);}

	ns=count;

	//TMP DEBUG ##########################################
	/*TMPINJ<<"count xstr flux break "<<endl;
	for (count=0;count<ns; count++){
		TMPINJ<<count<<" "<<xstr[count]<<" " << fluxes[count]<<" "<< breaks[count] << endl;
	}
	for (count=ns;count<bcount; count++){
		TMPINJ<<count<<" 0 0 "<< breaks[count] << endl;
	}
	TMPINJ.close();*/
	//TMP DEBUG ##########################################

}

//---------------------------------------------------------------------
void CTransect::WriteOutput(const double &t) const{
	double head;
	cmplex z;
	int i;

	ofstream TRANSECT;
	TRANSECT.open("TransectAnalysis.csv", ios::app);
	if (TRANSECT.fail()){cout <<"Can't open file TransectAnalysis.csv"<<endl;return;}
	double maxQ(0);
	for (i=0; i<=nx; i++){
		upperswap(maxQ,fabs(Qn[i]));
	}
	if (ID>=0){ //otherwise, part of Zone budget
		//ptID x y head Qnorm Qtan
		TRANSECT << "Transect,"<<ID<<","<<z1.real()<<","<<z1.imag()<<","<<z2.real()<<","<<z2.imag()<<","<<nx+1<<","<<maxQ <<","<<influx-outflux <<endl;
		for (i=0; i<=nx; i++){
			z=(x[i]-x[0])/(x[nx]-x[0])*(z2-z1)+z1;
			head=pLayer->GetHead(z,t);
			TRANSECT <<i       <<",";
			TRANSECT <<z.real()<<","<<z.imag()<<",";
			TRANSECT <<head    <<",";
			TRANSECT <<  Qn[i] <<","<<Qt[i]   <<endl;
		}
	}
}
/*******************************************************************************
			              CIRTRANSECT CLASS:
**********************************************************************
                     CONSTRUCTORS
*********************************************************************/
CCirTransect::CCirTransect():CAnalysisLocation(){			
	zcen=0.0;R=1.0;influx=outflux=budget=0.0;
	ncontrol=0;
	pLayer=NULL;
	zctrl =NULL;
	Qn=NULL;
	Qt=NULL;
	ID=0;
}
//---------------------------------------------------------------------
CCirTransect::CCirTransect(const int							my_ID,
													 const CSingleLayerABC *pLay,				//Layer
													 const cmplex						zc,			    //center of circle
													 const double						radius,			//radius of circle
													 const int							numdiv) 		//number of integration divisions
													:CAnalysisLocation(){
	ID      =my_ID;
	pLayer  =pLay;
	zcen    =zc;
	R       =radius;
	ncontrol=numdiv;
	ExitGracefullyIf(ncontrol<=0,"CCirTransect: negative or zero control points",BAD_DATA);

	zctrl=new cmplex [ncontrol];
	Qn   =new double [ncontrol];
	Qt   =new double [ncontrol];
	ExitGracefullyIf(Qt==NULL,"CCircular Transect: constructor: out of memory",OUT_OF_MEMORY);
	influx=outflux=budget=0.0;

	//initialize divisions
	double angle;
	for (int m=0; m<ncontrol; m++){
    angle=2.0*PI*(double)(m)/(double)(ncontrol);
    zctrl[m]= R*cmplex(cos(angle), sin(angle))+zcen; 
	}
}
//---------------------------------------------------------------------
CCirTransect::~CCirTransect(){
	delete [] zctrl;
	delete [] Qn;
	delete [] Qt;
}
/*********************************************************************
                     MEMBER FUNCTIONS
*********************************************************************/
void CCirTransect::CalculateStatistics(const double &t){

	//if flow moving out of circle, flux positive
	cmplex W;
	int m; 
	for (m=0; m<ncontrol; m++){
		W=pLayer->GetW(zctrl[m],t);
		Qn[m]=(conj(W)*    conj(zctrl[m]-zcen)/R).real(); //normal flux 
		Qt[m]=(conj(W)*-IM*conj(zctrl[m]-zcen)/R).real();//tangential flux 
	}
	//integrate normal flux alont front
	influx=outflux=0;
	double Qnext;
	cmplex znext;
	for (m=0; m<ncontrol;m++){
		if (m==ncontrol-1){Qnext=Qn[0];  znext=zctrl[0];}
		else              {Qnext=Qn[m+1];znext=zctrl[m+1];}
		if ((Qn[m]>=0) && (Qnext>=0)){
			influx +=(Qn[m]+Qnext)*abs(znext-zctrl[m])/2.0;  //use trapezoidal rule
		}
		else if ((Qn[m]<=0) && (Qnext<=0)){
		  outflux+=(Qn[m]+Qnext)*abs(znext-zctrl[m])/2.0;	//use trapezoidal rule
		}
		else{
			//find intercept-use area of triangles
			cmplex zint=-Qn[m]*(znext-zctrl[m])/(Qnext-Qn[m])+zctrl[m];
			if (Qn[m]>Qnext){
				influx +=Qn[m]*abs(zint -zctrl[m])/2.0;
				outflux+=Qnext*abs(znext-zint)/2.0;
			}
			else{
				influx +=Qnext*abs(znext-zint)/2.0;
				outflux+=Qn[m]*abs(zint -zctrl[m])/2.0;
			}
		}
	}
}//---------------------------------------------------------------------
void CCirTransect::WriteOutput(const double &t) const{
	double head;
	int    m;
	double maxQ(0);

	ofstream TRANSECT;
	TRANSECT.open("CircularBudget.csv", ios::app);
	if (TRANSECT.fail()){cout <<"Can't open file CircularBudget.csv"<<endl;return;}

	for (m=0; m<ncontrol; m++){
		upperswap(maxQ,fabs(Qn[m]));
	}
	TRANSECT << "Transect_"<<ID<<","<<zcen.real()<<","<<zcen.imag()<<","<<R<<","<<ncontrol<<","<<maxQ <<","<<budget<<endl;
	//ptID x y head Qnorm Qtan
	for (m=0; m<ncontrol; m++){
		head=pLayer->GetHead(zctrl[m],t);
		TRANSECT <<m              <<",";
		TRANSECT <<zctrl[m].real()<<","<<zctrl[m].imag()<<",";
		TRANSECT <<head           <<",";
		TRANSECT <<Qn[m]          <<","<<Qt[m]          <<endl;
	}
}
/*******************************************************************************
			              ZONE BUDGET CLASS:
**********************************************************************
                     CONSTRUCTORS
*********************************************************************/
//---------------------------------------------------------------------
CZoneBudget::CZoneBudget():CAnalysisLocation(){
	pSides=NULL;
	nsides=0;
	budget=0;
	influx=0;
	outflux=0;
	ID=0;
}
//---------------------------------------------------------------------
CZoneBudget::CZoneBudget(const CSingleLayerABC *pLay,
												 const cmplex					 *zp,     //first/last point counted twice
												 const int							nlines, //number of sides 
												 const int							numdiv)
												:CAnalysisLocation(){
	int i;
	nZB++;
	ID=nZB-1;
	nsides=nlines;
	budget=0;
	influx=0;
	outflux=0;
	pSides=new CTransect*[nsides];

	for(i=0;i<nsides;i++){
		pSides[i]=new CTransect(-1,pLay,zp[i],zp[i+1],numdiv);
	}
	
}
//---------------------------------------------------------------------
CZoneBudget::~CZoneBudget(){
	/*if (pSides!=NULL){
		for (int i=0; i<nsides; i++){delete pSides[i];}
	}*/
	delete [] pSides;
}
//---------------------------------------------------------------------
int CZoneBudget::nZB=0;
//*********************************************************************
void CZoneBudget::CalculateStatistics(const double &t){
	budget=influx=outflux=0.0;
	for (int i=0; i<nsides; i++){
		pSides[i]->CalculateStatistics(t);
		budget +=pSides[i]->GetTotalFlux();
		influx +=pSides[i]->GetInflux();
		outflux+=pSides[i]->GetOutflux();
	}
}
//*********************************************************************
double CZoneBudget::CalculateMassBudget(const double &t, double (*conc)(const cmplex &z, const double &t)){
	double sum(0.0);
	for (int i=0; i<nsides; i++){
		sum-=pSides[i]->CalculateMassFlux(t,conc); //negative because mass is lost through transects oriented outwards
	}
	return sum;
}
//*********************************************************************
void CZoneBudget::WriteOutput(const double &t) const{
	ofstream ZB;
	ZB.open("ZoneBudget.csv",ios::app);
	ZB<<ID<<","<<influx<<","<<outflux<<","<<-budget<<endl; //budget is net outflux in output file
	ZB.close();
}

//*********************************************************************
CZoneBudget *CZoneBudget::Parse(ifstream &input, int &l,CSingleLayerABC *pLay){
  //string "ZoneBudget" 
	//{double x double y}x(numlines+1)
	//&

	CZoneBudget  *pZB=NULL;
  bool          eof(false),done(false);
	int           Len,nlines(0);
	cmplex        stringp[MAXLINES];
	char         *s[MAXINPUTITEMS];

	if (parserdebug){cout << "Zone Budget"<<endl;}        eof=TokenizeLine(input,s,Len); l++; 
  do {
    ExitGracefullyIf(nlines>=MAXLINES,"CZoneBudget::Parse- too many lines in zone budget polygon",TOO_MANY);
	  if (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]);nlines++;       eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			stringp[nlines]=stringp[0]; 
			pZB= new CZoneBudget(pLay,
												   stringp,
												   nlines,
													 DEFAULT_DIVS);
			done=true;}		
		else if (Len==0){                                  eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pZB;}
}



