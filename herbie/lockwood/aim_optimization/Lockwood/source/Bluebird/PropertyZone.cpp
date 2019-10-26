#include "PropertyZone.h"
//************************************************************************
//                           CPropZone Constructors
//************************************************************************
CPropZone::CPropZone(){
	ptype=notype;
	value=0;
	TotalPropZones++;
	pAllPropZones[TotalPropZones-1]=this;
}
//------------------------------------------------------------------------
CPropZone::CPropZone(const proptype typ,const double val){
	ptype  =typ;
	value  =val;
	pAllPropZones[TotalPropZones]=this;
	TotalPropZones++;
}
//------------------------------------------------------------------------
CPropZone::~CPropZone(){
	if (globaldebug){cout<<"  DESTROYING PROPERTY ZONE"<<endl;}
}
//************************************************************************
//									ACCESSORS
//************************************************************************
double   CPropZone::GetValue()                const{return value;}
//------------------------------------------------------------------------
proptype CPropZone::GetType ()                const{return ptype;}
//------------------------------------------------------------------------
double   CPropZone::GetArea ()                const{return area;}
//------------------------------------------------------------------------
double   CPropZone::GetValue(const cmplex &z) const{
	if (IsInside(z)){return value;   }
	else            {return NO_VALUE;}
}
//------------------------------------------------------------------------
cmplex   CPropZone::GetSlopes(const cmplex &z) const{
	return NO_VALUE;//for piecewise continuous, overriden by some children
}

//************************************************************************
bool CPropZone::IsInside  (const cmplex &z                ) const {cout <<"PZVirBug1"<<endl;return false; }
//------------------------------------------------------------------------
bool CPropZone::IsInSquare(const cmplex &zc,const double w) const {cout <<"PZVirBug2"<<endl;return false; }
/*************************************************************************
							STATIC MEMBER FUNCTIONS
*************************************************************************/
int         CPropZone::TotalPropZones=0;
//------------------------------------------------------------------------
CPropZone  *CPropZone::pAllPropZones[];
//------------------------------------------------------------------------
void CPropZone::DestroyAllPropZones(){
	if (globaldebug){cout <<"DESTROYING ALL PROPERTY ZONES"<<endl;}
	for (int i=0; i<TotalPropZones; i++){delete pAllPropZones[i];}
}

/************************************************************************
					Nested Sift
------------------------------------------------------------------------
Searches through an array (size=nzones) of pointers to abstract property zones (pZones)
returns the value of the zone in which z is located 
returns back_val if z is outside of all zones
MUST OPTIMIZE -this would be optimized if zonearray was sorted by size
************************************************************************/
double CPropZone::NestedSift(const cmplex  &z, 
														 CPropZone    **pZones, 
														 const int      nzones, 
														 const double   back_val){
	if (nzones==0){return back_val;}   //quick out

	static double vlocal;
	static double smallArea;
	static double thisval;
	static double tmpArea;	
	vlocal=thisval=NO_VALUE;
	smallArea=ALMOST_INF;
	for (int i=0; i<nzones; i++){      //sift thru zones
    vlocal=pZones[i]->GetValue(z);   //get value due to zone at z
		if (vlocal!=NO_VALUE){                
      tmpArea=pZones[i]->GetArea();  //get area of zone
			if (tmpArea<smallArea){        //if zone is smaller than smallest zone 
				thisval  =vlocal;            //then this is true value (accounts for nesting)
        smallArea=tmpArea;                 
			} 
		}
	}
	if (thisval==NO_VALUE){return back_val;}       //if no value is found, background is used
	else                  {return thisval; }
}
//--------------------------------------------------------------------
cmplex CPropZone::NestedSiftSlopes(const cmplex &z, CPropZone **pZones, const int nzones, const cmplex back_val){
	if (nzones==0){return back_val;}   //quick out

	static cmplex vlocal;
	static double smallArea;
	static cmplex thisval;
	static double tmpArea;	
	vlocal=thisval=NO_VALUE;
	smallArea=ALMOST_INF;
	for (int i=0; i<nzones; i++){      //sift thru zones
    vlocal=pZones[i]->GetValue(z);   //get value due to zone at z
		if (vlocal!=NO_VALUE){                
      tmpArea=pZones[i]->GetArea();  //get area of zone
			if (tmpArea<smallArea){        //if zone is smaller than smallest zone 
				thisval  =vlocal;            //then this is true value (accounts for nesting)
        smallArea=tmpArea;                 
			} 
		}
	}
	if (thisval==NO_VALUE){return back_val;}       //if no value is found, background is used
	else                  {return thisval; }
}
/************************************************************************##################################
					CPolyPropZone
************************************************************************/
CPolyPropZone::CPolyPropZone(){}
//------------------------------------------------------------------------
CPolyPropZone::CPolyPropZone(const proptype typ, 
														 const cmplex *points, 
														 const int NumOfLines, 
														 const double val):
               CPropZone(typ,val){
	int i;
	bool clock(false);
 
	NLines   =NumOfLines;
	ze=new cmplex [NLines+1];
  boundingbox.e=-ALMOST_INF;
  boundingbox.w= ALMOST_INF;
	boundingbox.n=-ALMOST_INF;
	boundingbox.s= ALMOST_INF;
	for(i=0;i<=NLines;i++){
		ze[i]=points[i];
		if (ze[i].real()>boundingbox.e){boundingbox.e=ze[i].real();}
		if (ze[i].real()<boundingbox.w){boundingbox.w=ze[i].real();}
		if (ze[i].imag()>boundingbox.n){boundingbox.n=ze[i].imag();}
		if (ze[i].imag()<boundingbox.s){boundingbox.s=ze[i].imag();}
	}
	//size up bounding box
	boundingbox.n+=BIGGERIZE*(boundingbox.n-boundingbox.s);
  boundingbox.s-=BIGGERIZE*(boundingbox.n-boundingbox.s)/(1.0+BIGGERIZE);
  boundingbox.e+=BIGGERIZE*(boundingbox.e-boundingbox.w);
  boundingbox.w-=BIGGERIZE*(boundingbox.e-boundingbox.w)/(1.0+BIGGERIZE);
	area=Area(clock);
  if (clock){for(i=0;i<=NLines;i++){ze[i]=points[NLines-i];}}
}
//------------------------------------------------------------------------
CPolyPropZone::~CPolyPropZone(){
	if (globaldebug){cout <<"    DESTROYING POLYPROPZONE "<<endl;}
	delete [] ze;
}
//************************************************************************
double CPolyPropZone::GetValue (const cmplex &z) const{return CPropZone::GetValue(z);}
cmplex CPolyPropZone::GetSlopes(const cmplex &z) const{return CPropZone::GetSlopes(z);}
//************************************************************************
bool CPolyPropZone::IsInside(const cmplex &z) const{
	//MUST OPTIMIZE
	double y(z.imag());
	double x(z.real()); 

	//if not in bounding box, return false
  if (x>=boundingbox.e){return false;}
	if (x<=boundingbox.w){return false;}
	if (y>=boundingbox.n){return false;}
	if (y<=boundingbox.s){return false;}

	int  i;
  //check if close to boundary of zone
	bool close(false);
	for (i=0; i<NLines; i++) {
		if (((abs(z-ze[i])+abs(z-ze[i+1]))/(abs(ze[i+1]-ze[i])))<CLOSE_TO_LINE){close=true;break;}
	}

  if (!close){
	  //if not close, use faster algorithm
	  //borrowed from Wm. Randolph Franklin <wrf@ecse.rpi.edu>
		bool   in(false);
    int    j;
    for (i=0, j=NLines-1; i<NLines; j=i++) {
      if ((((ze[i].imag()<=y) && (y<ze[j].imag())) ||
	         ((ze[j].imag()<=y) && (y<ze[i].imag()))) &&
         (x<(ze[j].real()-ze[i].real())*(y-ze[i].imag())/
				 (ze[j].imag()-ze[i].imag())+ze[i].real())){
				in=!in;
			}
		}
    return in;
	}
	else{
		//if close, use complex algortithm
		cmplex Z;
		double sum(0);
		for (int i=0;i<NLines; i++){
			Z=(z-0.5*(ze[i]+ze[i+1]))/(0.5*(ze[i+1]-ze[i]));
			sum+=log(Zm1oZp1(Z)).imag();
		}
		return (fabs(sum)<REALSMALL ? false : true);
	}
}
//************************************************************************
double CPolyPropZone::Area(bool &clockwise) const{
  //borrowed from Paul Bourke <http://astronomy.swin.edu.au/pbourke/>
  int    j;
  double area(0.0);

  for (int i=0;i<NLines;i++) {
    j = (i+1)%NLines;
    area+= ze[i].real()*ze[j].imag()-ze[i].imag()*ze[j].real();
  }
  area/=2.0;
  clockwise=((area>=0.0) ? false : true);
  return fabs(area);
}
//************************************************************************
window  *CPolyPropZone::GetBoundingBox() const{
  window *tmp;
	tmp=new window;
	*tmp=boundingbox;
	return tmp;
}
//************************************************************************
bool CPolyPropZone::IsInSquare(const cmplex &zc,const double w) const{
	for (int i=0; i<=NLines; i++){
	if (  (ze[i].real() > zc.real() + (w/2.0)) ||
        (ze[i].real() < zc.real() - (w/2.0)) ||
        (ze[i].imag() > zc.imag() + (w/2.0)) ||
        (ze[i].imag() < zc.imag() - (w/2.0)))     {return false;}
  }
  return true;
}
//************************************************************************
CPropZone *CPolyPropZone::Parse(ifstream &input, int &l,proptype typ){
	//string "PropZone", string name 
	//double value
	//{double x double y}x(numlines+1)
	//&[int precision]

	CPropZone  *pZone;
  bool     eof(false),done(false);
  double   thisval(0);
	int      Len,nlines(0);
  cmplex	 stringp[MAXPZLINES];
	char    *s[MAXINPUTITEMS];

  pZone=NULL;
	if (parserdebug) {cout << "Property Zone"<<endl;}    eof=TokenizeLine(input,s,Len); l++;
	do{ 
		if      (Len==1) {thisval=s_to_d(s[0]); done=true; eof=TokenizeLine(input,s,Len); l++;}
    else if (Len==0) {                                 eof=TokenizeLine(input,s,Len); l++;}
		else             {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;  
  do {
	  if (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]);nlines++;       eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			if (nlines>MAXPZLINES){ExitGracefully("PolyPropZone::Parse: too many lines specified for polygon",BAD_DATA);}
			stringp[nlines]=stringp[0]; nlines++;
			pZone= new CPolyPropZone(typ,
															 stringp,
												       nlines-1,
															 thisval);
			done=true;
		}		
		else if (Len==0){                                  eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pZone;}
}
//************************************************************************
void CPolyPropZone::GetPoly(cmplex *points, int &nlines){
	points=new cmplex[NLines+1];
  for (int i=0; i<=NLines; i++){
	  points[i]=ze[i];
	}
	nlines=NLines;
}

//************************************************************************##################################
//
//************************************************************************
CMQPolyPropZone::CMQPolyPropZone(){
	values=NULL;
	Zctrl=NULL;
	nbasispts=0;
	anisotropy=1.0;
	anis_angle=0.0;
	zref=0.0;
}
//************************************************************************
CMQPolyPropZone::CMQPolyPropZone(const proptype  typ,
																 const cmplex   *vertices,
																 const int       NumOfLines,
																 const cmplex   *basispts,
																 const double   *specvalues,
																 const int       nbasis,
																 const double    R,
																 const double    anis,
																 const double    anis_ang):
                 CPolyPropZone(typ,vertices,NumOfLines,NO_VALUE){
	nbasispts=nbasis;
	radius=R;
	values =new double [nbasispts];
	Zctrl  =new cmplex [nbasispts];
	MQcoeff=new double [nbasispts];
	anisotropy=anis;
	anis_angle=anis_ang;
	zref=basispts[0];

	
	for(int n=0; n<nbasispts;n++){
		values[n]=specvalues[n];
		Zctrl [n]=cmplex((basispts[n]-zref).real(),(basispts[n]-zref).imag()*anisotropy)+zref;//TMP DEBUG - must include orientation
		value+=values[n]/(double)(nbasispts);//average value
	}
	cout<<"Calculating Multiquadric coefficients"<<endl;
	MQInterpolate(values,Zctrl,nbasispts,MQcoeff,Zctrl,nbasispts,radius,value);
	cout<<" ...Multiquadric coefficients calculated successfully"<<endl;
}
//************************************************************************
CMQPolyPropZone::~CMQPolyPropZone(){
	if (globaldebug){cout <<"  DESTROYING MQ POLYPROPZONE "<<endl;}
	delete [] values;
	delete [] Zctrl;
	delete [] MQcoeff;
}
//************************************************************************
double CMQPolyPropZone::GetValue(const cmplex &z) const{
	double v(0.0);
	cmplex Z;

	if (IsInside(z)){
		if (nbasispts==1){
			return MQcoeff[0]+value; 
		}
    else {
			
		//	Z=(z-zref)*anisotropy+zref;//TMP DEBUG - must include orientation
		  Z=cmplex((z-zref).real(),(z-zref).imag()*anisotropy)+zref;
      for(int n=0; n<nbasispts;n++){
				v+=MQcoeff[n]*sqrt(pow(abs(Z-Zctrl[n]),2)+radius*radius);
			}
      v+=value; //add average value
		}
		return max(0.0,v);
	}
	return NO_VALUE;
}
//-------------------------------------------------------------------
cmplex CMQPolyPropZone::GetSlopes(const cmplex &z) const{
	cmplex v(0.0);
	if (IsInside(z)){
		if (nbasispts<=1){
			return 0.0; //no slope 
		}
    else {
      for(int n=0; n<nbasispts;n++){
				return -cmplex((GetValue(z)-GetValue(z+   0.01))/0.01, //TMP DEBUG
					             (GetValue(z)-GetValue(z+IM*0.01))/0.01);
			//	if (abs(z-zctrl[n])>0.0){v+=MQcoeff[n]*(z-zctrl[n])/abs(z-zctrl[n]);} //d/dx=an*(x-xn)/r , d/dy=an*(y-yn)/r
			}
		}
	}
	return NO_VALUE;//outside
}
//************************************************************************
CMQPolyPropZone *CMQPolyPropZone::Parse(ifstream &input, int &l,proptype typ){
	//string "MQPropZone", string name 
	//{double x double y}x(numlines+1)
	//&
	//{double x double y}x(numbasis+1)
	//&

	CMQPolyPropZone *pZone=NULL;
  bool       done(false);
	int        Len,nlines(0),npoints(0);
  cmplex	   stringp[MAXPZLINES];       
  cmplex	   leakzp [MAXPZPTS];
  double     leakp  [MAXPZPTS];
	char      *s[MAXINPUTITEMS];

	if (parserdebug) {cout << "Multiquadric Property Zone"<<endl;}  
	
	if (TokenizeLine(input,s,Len)){return NULL;}; l++;
  do {
		if (nlines>=MAXPZLINES) { ExitGracefully("CMQPolyPropZone::Parse- too many line segments in MQ property zone",TOO_MANY);}
		if  (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]); nlines++; 
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=3) && (!strcmp(s[0],"&"))) {
      stringp[nlines]=stringp[0];  nlines++;       	
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
			done=true; 
		}
    else {
			ImproperFormat(s,l); return NULL;
		}
	} while (!done);

	done=false;                                      
	do {
    if (npoints>=MAXPZPTS) { ExitGracefully("CMQPolyPropZone::Parse- too many control points in MQ Propzone",TOO_MANY);}
    if ((Len==3) && (strcmp(s[0],"&"))){
      leakzp[npoints]=s_to_c(s[0],s[1]);
			leakp [npoints]=s_to_d(s[2]);
			npoints++;                                   
			if (TokenizeLine(input,s,Len)){return NULL;}; l++;
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			pZone     = new CMQPolyPropZone(typ,
				                             stringp,
																		 nlines-1,
																		 leakzp,
																		 leakp,
																		 npoints,0.0,1.5,0.0);//zero, for now
			done=true;
		}

		else {ImproperFormat(s,l); return NULL;}
	} while (!done);

  return pZone;
}
//************************************************************************##################################
//         LINPOLYPROPZONE
//************************************************************************
CLinPolyPropZone::CLinPolyPropZone(){}
CLinPolyPropZone::CLinPolyPropZone(const proptype  typ,
																	 const cmplex   *vertices,
																	 const int       NumOfLines,
																	 const cmplex    center,
																	 const double    valcen,
																	 const cmplex    slope_dir):
                  CPolyPropZone(typ,vertices,NumOfLines,NO_VALUE){
	zcen=center;
	center_value=valcen;
	slope=slope_dir;
}
CLinPolyPropZone::~CLinPolyPropZone(){}
//------------------------------------------------------------------------
//Static Member Functions
CLinPolyPropZone *CLinPolyPropZone::Parse(ifstream &input, int &l,proptype typ){
	return NULL; //TMP DEBUG- must fill out 
}
//------------------------------------------------------------------------
//Accessor Functions (inherited from CPropzone)
double CLinPolyPropZone::GetValue      (const cmplex &z) const{
	return slope.real()*(z-zcen).real()+
		     slope.imag()*(z-zcen).imag()+center_value;
}
//------------------------------------------------------------------------
cmplex CLinPolyPropZone::GetSlopes     (const cmplex &z) const{
	return slope;
}
//************************************************************************##################################
//
//************************************************************************
CCirPropZone::	 CCirPropZone(const proptype typ, 
															const cmplex	 zc, 
															const double	 rad, 
															const double	 val):
                 CPropZone(typ,val){
	zcen=zc;
	radius=rad;
	area=PI*radius*radius;
}
//------------------------------------------------------------------------
CCirPropZone::~CCirPropZone(){
	if (globaldebug){cout <<"    DESTROYING CIRPROPZONE "<<endl;}
}
//************************************************************************
bool CCirPropZone::IsInside(const cmplex &z) const{
	return ((abs(z-zcen)<=radius) ? true : false);
}

//************************************************************************
bool CCirPropZone::IsInSquare(const cmplex &zc,const double w) const{
  if ((zcen.real()-radius > zc.real()-w/2) &&
      (zcen.real()+radius < zc.real()+w/2) &&
      (zcen.imag()-radius > zc.imag()-w/2) &&
      (zcen.imag()+radius < zc.imag()+w/2)) {return true; }
	else                                      {return false;}
}

//************************************************************************##################################
//
//************************************************************************
CEllPropZone::	 CEllPropZone(const proptype typ, 
															const cmplex	 zcen, 
															const double	 major, 
															const double	 minor, 
															const double	 orient, 
															const double	 val):
                 CPropZone(typ,val){
	zc=zcen;
	a=major;
	b=minor;
	angle=orient;
	area=PI*a*b;
}
//------------------------------------------------------------------------
CEllPropZone::~CEllPropZone(){
	if (globaldebug){cout <<"    DESTROYING ELLPROPZONE "<<endl;}
}
//************************************************************************
bool CEllPropZone::IsInside(const cmplex &z) const{
	cmplex Z=(z-zc)*exp(-IM*angle);
	return ((Z.real()*Z.real()/(a*a) + Z.imag()*Z.imag()/(b*b)) <= 1);
}

//************************************************************************
bool CEllPropZone::IsInSquare(const cmplex &zc,const double w) const{
  if ((zc.real()-a > zc.real()-w/2) &&
      (zc.real()+a < zc.real()+w/2) &&
      (zc.imag()-a > zc.imag()-w/2) &&
      (zc.imag()+a < zc.imag()+w/2)) {return true; }
	else                                 {return false;}
}
/***********************************************************************
                           PARSE
************************************************************************
Format:
	string Tag [string name (ignored)]
	double x1 double y1 double value double A double B double angle (degrees)
  & [int precision]
----------------------------------------------------------------------*/
CEllPropZone  *CEllPropZone::Parse(ifstream &input,int &l,proptype type){

	CEllPropZone *pEllPZone;
  bool					eof(false),done(false);
  double				thisval(0);
	cmplex				thisz;
	double        thisa,thisb,thisang;
	int						Len;
	char				 *s[MAXINPUTITEMS];

	pEllPZone=NULL;
	if (parserdebug) {cout << "Elliptical Inhom element"<<endl;}   
																												 eof=TokenizeLine(input,s,Len); l++; 
	do {
		if (Len==6) {
      thisz   =s_to_c(s[0],s[1]);
			thisval =s_to_d(s[2]);
			thisa   =s_to_d(s[3]);
			thisb   =s_to_d(s[4]);
			thisang =s_to_d(s[5]);                            
			thisang =PI/180.0*thisang;                         eof=TokenizeLine(input,s,Len); l++; 
		}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			pEllPZone = new CEllPropZone(type,
																	 thisz,
																	 thisa, 
																	 thisb,
																	 thisang,
																	 thisval); 
			done=true;
		}
		else if(Len==0) {                                    eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while (!done);

  if (eof) {return NULL;}
	else     {return pEllPZone;}
}
