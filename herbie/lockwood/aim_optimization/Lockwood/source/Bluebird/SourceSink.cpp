//#include "SpeciesAndSoils.h"
#include "RxNLib.h"
#include "SourceZone.h"

/************************************************************************
           SOURCE ZONE CLASS
************************************************************************/
//CONSTRUCTORS
CAreaSource::CAreaSource(sourcetype         ty,
												 double						 *conc,
												 double            *sorb,
												 const double       startt,
												 const double       endt,
												 const int          NumSpecies){
	//if type=constant, initial, or recharge, conc is the concentration 
	//if type=flux, then conc is the concentration flux C/t
  ConcArray=NULL;
	SorbArray=NULL;
	ConcSeries=NULL;
	type     =ty;
	starttime=startt;
	endtime  =endt;
	if (endtime<starttime){ExitGracefully("CAreaSource::Constructor: endtime before starttime",BAD_DATA);}
  nspecies =NumSpecies;
	zone     =NULL;
	ConcArray=new double[nspecies];
	SorbArray=new double[nspecies];
	if (SorbArray==NULL){ExitGracefully("CAreaSource::Constructor",OUT_OF_MEMORY);}
	for (int s=0; s<nspecies; s++){
		ConcArray[s]=conc[s];
		SorbArray[s]=sorb[s];
		if ((ConcArray[s]<0.0) && (type==SPECIFIED_CONC         )){ConcArray[s]=NO_INFLUENCE;}
		if ((ConcArray[s]<0.0) && (type==SPECIFIED_CONC_RECHARGE)){ConcArray[s]=NO_INFLUENCE;}
		if ((SorbArray[s]<0.0)                                   ){SorbArray[s]=NO_INFLUENCE;}
	}
}
//----------------------------------------------------------------------------
CAreaSource::CAreaSource(sourcetype         ty,
						             CTimeSeries      **TimeSeries,
						             double            *sorb,
												 const int          NumSpecies){

	ExitGracefullyIf(TimeSeries==NULL,"CAreaSource::Constructor: NULL Time series",BAD_DATA);

  ConcArray =NULL;
	SorbArray =NULL;
	ConcSeries=NULL;
	zone      =NULL;
	type      =ty;	
	nspecies  =NumSpecies;
	
	SorbArray =new double       [nspecies];
	ConcSeries=new CTimeSeries *[nspecies];
	ExitGracefullyIf(ConcSeries==NULL,"CAreaSource::Constructor",OUT_OF_MEMORY);
	for (int s=0; s<nspecies; s++){
		SorbArray[s]=sorb[s];
		if (SorbArray[s]<0.0){SorbArray[s]=NO_INFLUENCE;}
		ConcSeries[s]=TimeSeries[s];
	}
}
//------------------------------------------------------------------------------
CAreaSource::~CAreaSource(){
	delete [] ConcArray;
	delete [] SorbArray;
	delete [] ConcSeries; //Just deletes pointers
}
/******************************************************************************
															MANIPULATORS
------------------------------------------------------------------------------
******************************************************************************/
void CAreaSource::SetTimeSeries(CTimeSeries **TimeSeries, const int numseries){
	ExitGracefullyIf(numseries!=nspecies,"CAreaSource::SetTimeSeries: inappropriate number of time series",BAD_DATA);
	ConcSeries=new CTimeSeries *[nspecies];
	for (int i=0;i<nspecies;i++){
		ConcSeries[i]=TimeSeries[i];
	}
}
/******************************************************************************
															ACCESSORS
------------------------------------------------------------------------------
******************************************************************************/
double        CAreaSource::GetArea() const {return zone->GetArea();}
//------------------------------------------------------------------------------
sourcetype    CAreaSource::GetType() const {return type;}
//------------------------------------------------------------------------------
void CAreaSource::GetRechargeConcentrations(const cmplex &z, const double &t, double *conc) const {
	//returns recharge concentration, mg/L
	if (ConcSeries==NULL){
		if ((type==SPECIFIED_CONC_RECHARGE) && (t<=endtime) && (t>=starttime))
		{
			//cout <<"Getting Recharge Concs!"<<endl;
			for (int s=0; s<nspecies; s++){conc[s]=ConcArray[s];}
		}
		else 
		{
			for (int s=0; s<nspecies; s++){conc[s]=NO_INFLUENCE;}
		}
	}
	else{
		if (type==SPECIFIED_CONC_RECHARGE){
			for (int s=0; s<nspecies; s++){
				conc[s]=ConcSeries[s]->GetValue(t);
				if (conc[s]<0){conc[s]=NO_INFLUENCE;}
			}
		}
		else                              {for (int s=0; s<nspecies; s++){conc[s]=NO_INFLUENCE;}}		
	}
}
//------------------------------------------------------------------------------
void CAreaSource::GetSpecifiedConcentrations(const cmplex &z, const double &t, double *conc) const {
	//returns concentrations: mg/L
	if (ConcSeries==NULL){
	if ((t<=endtime) && (t>=starttime)){
		if (type==SPECIFIED_CONC )        {for (int s=0; s<nspecies; s++){conc[s]=ConcArray[s];}}	
		else                              {for (int s=0; s<nspecies; s++){conc[s]=NO_INFLUENCE;}}
	}
	else                                {for (int s=0; s<nspecies; s++){conc[s]=NO_INFLUENCE;}}
	}
	else{
		if (type==SPECIFIED_CONC){
			for (int s=0; s<nspecies; s++){
				conc[s]=ConcSeries[s]->GetValue(t);
				if (conc[s]<0){conc[s]=NO_INFLUENCE;}
			}
		}
		else                              {for (int s=0; s<nspecies; s++){conc[s]=NO_INFLUENCE;}}		
	}
}
//------------------------------------------------------------------------------
void CAreaSource::GetInitialSorbedConcentrations(const cmplex &z, const double &t, double *sorb) const {
	//returns concentrations: mg/kg
	if ((t<=endtime) && (t>=starttime)){
		if (type==SPECIFIED_CONC )        {for (int s=0; s<nspecies; s++){sorb[s]=SorbArray[s];}}	
		else                              {for (int s=0; s<nspecies; s++){sorb[s]=NO_INFLUENCE;}}
	}
	else                                {for (int s=0; s<nspecies; s++){sorb[s]=NO_INFLUENCE;}}
}
//------------------------------------------------------------------------------
void CAreaSource::GetSpecifiedMassFlux(const cmplex &z, const double &t, double *flux) const{
	//Returns mass flux per unit area: mg/L *L/T
  if ((t<=endtime) && (t>=starttime)){
		if (type==SPECIFIED_FLUX)         {for (int s=0; s<nspecies; s++){flux[s]=ConcArray[s];}}
		else                              {for (int s=0; s<nspecies; s++){flux[s]=0.0;}}
	}
	else                                {for (int s=0; s<nspecies; s++){flux[s]=0.0;}}
}
/******************************************************************************
Member Functions
------------------------------------------------------------------------------
******************************************************************************/
bool CAreaSource::IsInside(const cmplex z) const{return zone->IsInside(z);}

/************************************************************************
           POLYGONAL SOURCE ZONE CLASS
************************************************************************/
CPolyAreaSource::CPolyAreaSource(cmplex						 *points, 
																 int								npoints, 
																 sourcetype         ty,
																 double						 *conc,
																 double            *sorb,
																 const double       startt,
																 const double       endt,
																 const int          NumSpecies):
                 CAreaSource(ty,conc,sorb,startt,endt,NumSpecies){
	//if type=constant, initial, or recharge, conc is the concentration 
	//if type=flux, then conc is the concentration flux C/t

	zone=new CPolyPropZone(source_zone,points,npoints,1.0);
}
//-----------------------------------------------------------------------
CPolyAreaSource::CPolyAreaSource(cmplex						*points, 
																 int               npoints, 
																 sourcetype        ty,
														     CTimeSeries		 **timeseries,
																 double						*sorb,
																 const int         NumSpecies):
                 CAreaSource(ty,timeseries,sorb,NumSpecies){
	//if type=constant, initial, or recharge, conc is the concentration 
	//if type=flux, then conc is the concentration flux C/t

	zone=new CPolyPropZone(source_zone,points,npoints,1.0);

}
/*****************************************************************************/
CPolyAreaSource *CPolyAreaSource::Parse(ifstream &input, sourcetype ty, int &l,const int NumSpecies){
	//string "SourceZone" or "InitConcZone" or some other variation
	//double starttime
	//double endtime
	//{double x double y}x(numlines+1)
	//&
	//{double conc double sorbconc}x(numspecies) ->if C<0, then NO_INFLUENCE
	//&

	CPolyAreaSource *pZone;

	int          nspecies(NumSpecies);
  bool         eof(false),done(false);
	int          Len,nlines(0),speccount(0);
	double			 startt,endt;
	char        *s[MAXINPUTITEMS];
  cmplex	     stringp[MAXPZLINES];
	double       initconc[MAX_SPECIES];
	double       initsorb[MAX_SPECIES];

  pZone=NULL;
	                                                      eof=TokenizeLine(input,s,Len); l++; 
	if      (Len==1){startt=s_to_d(s[0]);                 eof=TokenizeLine(input,s,Len); l++;}
	else if (Len==0){                                     eof=TokenizeLine(input,s,Len); l++;}
  else {cout <<"line"<< l << "is wrong length"<<endl;}
	if      (Len==1){endt=s_to_d(s[0]);                   eof=TokenizeLine(input,s,Len); l++;}
	else if (Len==0){                                     eof=TokenizeLine(input,s,Len); l++;}
  else {cout <<"line"<< l << "is wrong length"<<endl;}
  do {
	  if (Len==2){
      stringp[nlines]=s_to_c(s[0],s[1]);nlines++;       eof=TokenizeLine(input,s,Len); l++;}
    else if ((Len<=2) && (!strcmp(s[0],"&"))) {
			done=true;
		}		
		else if (Len==0){                                   eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while ((!done) && (!eof));
	done=false;                                           eof=TokenizeLine(input,s,Len); l++;
	do {
		if ((Len==2) && (speccount<nspecies)){
			initconc[speccount]=s_to_d(s[0]); 
			initsorb[speccount]=s_to_d(s[1]); 
			speccount++;                                      eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (speccount<nspecies) && (strcmp(s[0],"&"))){
			initconc[speccount]=s_to_d(s[0]); 
			initsorb[speccount]=-1; 
			speccount++;                                      eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (speccount==nspecies) && (!strcmp(s[0],"&"))) {//hit a &
			stringp[nlines]=stringp[0]; nlines++;
			pZone= new CPolyAreaSource(stringp,
																 nlines-1,
																 ty,
																 initconc,
																 initsorb,
																 startt,
																 endt,
																 NumSpecies);
			done=true;
		}
		else if (Len==1){ 
			cout <<"Too many species associated with source"<<endl;  return NULL;	
		}
		else if (Len==0){                                   eof=TokenizeLine(input,s,Len); l++;}
	}while ((!done) && (!eof));

  if (eof) {return NULL;}
	else     {return pZone;}
}
/************************************************************************
           ELLIPTICAL SOURCE ZONE CLASS
************************************************************************/
CEllAreaSource::CEllAreaSource  (const cmplex		  zcen,        //center of ellipse
																 const double     MajorAxis,   //major axis of ellipse
																 const double		  MinorAxis,   //minor axis of ellipse
															   const double     angle,			 //angle (in radians)  
																 sourcetype       ty,
																 double					 *conc,
																 double          *sorb,
																 const double     startt,
																 const double     endt,
																 const int        NumSpecies):
                 CAreaSource(ty,conc,sorb,startt,endt,NumSpecies){
	//if type=constant, initial, or recharge, conc is the concentration 
	//if type=flux, then conc is the concentration flux C/t

	zone=new CEllPropZone(source_zone,zcen,MajorAxis,MinorAxis,angle,1.0);	
}
/*****************************************************************************/
CEllAreaSource *CEllAreaSource::Parse(ifstream &input, sourcetype ty, int &l,const int NumSpecies){
	//string "EllSourceZone" or "EllInitConcZone" or some other variation
	//double starttime
	//double endtime
	//double x1 double y1 double A double B double angle (degrees)
	//&
	//{double conc}x(numspecies) ->if C<0, then NO_INFLUENCE
	//&

	CEllAreaSource *pZone;

	int          nspecies(NumSpecies);
  bool         eof(false),done(false);
	int          Len,nlines(0),speccount(0);
	double			 startt,endt;
	char        *s[MAXINPUTITEMS];	
	cmplex			 thisz;
	double       thisa,thisb,thisang;
	double       initconc[MAX_SPECIES];
	double       initsorb[MAX_SPECIES];
	
  pZone=NULL;
	                                                      eof=TokenizeLine(input,s,Len); l++; 
	if      (Len==1){startt=s_to_d(s[0]);                 eof=TokenizeLine(input,s,Len); l++;}
	else if (Len==0){                                     eof=TokenizeLine(input,s,Len); l++;}
  else {cout <<"line"<< l << "is wrong length"<<endl;}

	if      (Len==1){endt=s_to_d(s[0]);                   eof=TokenizeLine(input,s,Len); l++;}
	else if (Len==0){                                     eof=TokenizeLine(input,s,Len); l++;}
  else {cout <<"line"<< l << "is wrong length"<<endl;}

	done=false;  
	do {
		if (Len==5) {
      thisz   =s_to_c(s[0],s[1]);
			thisa   =s_to_d(s[2]);
			thisb   =s_to_d(s[3]);
			thisang =s_to_d(s[4]);                            
			thisang =PI/180.0*thisang;                         //eof=TokenizeLine(input,s,Len); l++; 
      done=true;
		}
		else if(Len==0) {                                    eof=TokenizeLine(input,s,Len); l++;}
    else {cout <<"line"<< l << "is wrong length"<<endl; break;}
	} while (!done);

	done=false;                                           eof=TokenizeLine(input,s,Len); l++;
	speccount=0;
	do {
		if (ProgramAborted()){return NULL;}
		if ((Len==2) && (speccount<nspecies)){
			initconc[speccount]=s_to_d(s[0]); 
			initsorb[speccount]=s_to_d(s[1]);
			//cout <<"species "<<speccount<<" added"<<endl;
			speccount++;                                      eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (speccount<nspecies) && (strcmp(s[0],"&"))){
			initconc[speccount]=s_to_d(s[0]); 
			initsorb[speccount]=-1;
			speccount++;                                      eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len>=1)  && (speccount==nspecies) && (strcmp(s[0],"&"))){ 
			cout <<"Too many species associated with source"<<endl;  return NULL;	
		}
		else if ((Len==1)  && (speccount<nspecies) && (!strcmp(s[0],"&"))){ 
			cout <<"Not enough species associated with source"<<endl;  return NULL;	
		}		
		else if ((Len==1) && (speccount==nspecies) && (!strcmp(s[0],"&")) ) {//hit a &
			pZone= new CEllAreaSource(thisz,
				                        thisa,
																thisb,
																thisang,
																ty,
																initconc,
																initsorb,
																startt,
																endt,NumSpecies);
			done=true;
		}
		else if (Len==0){                                   eof=TokenizeLine(input,s,Len); l++;}
	}while ((!done) && (!eof));

  if (eof) {cout <<"Reached end of file during parse of Elliptical source zone"<<endl;return NULL;}
	else     {return pZone;}
}
/***************************************************************************************
****************************************************************************************
			CLASS CPOINTSOURCE
***************************************************************************************/
CPointSource::CPointSource(sourcetype  ty,
													 cmplex	     z,  
													 double	    *initconc,
													 double      startt,
													 double      endt,
													 const int	 NumSpecies){
	type     =ty;
	starttime=startt;
	endtime  =endt;
  nspecies =NumSpecies;
	zs       =z;

	ConcArray   =new double    [nspecies];
	for (int s=0; s<nspecies; s++){
		ConcArray[s]=initconc[s];
		if (ConcArray[s]<0.0){ConcArray[s]=NO_INFLUENCE;}
	}
}
//------------------------------------------------------------------------------
CPointSource::~CPointSource(){
	delete [] ConcArray;
}
/******************************************************************************
            ACCESSORS
******************************************************************************/

cmplex CPointSource::GetLocation() const {return zs;}
//------------------------------------------------------------------------------
void   CPointSource::GetSpecifiedConcentrations(const double &t, double *conc) const {
  if ((t<=endtime) && (t>=starttime)){
		if (type==SPECIFIED_CONC){for (int s=0; s<nspecies; s++){conc[s]=ConcArray[s];}}	
	  else                     {for (int s=0; s<nspecies; s++){conc[s]=NO_INFLUENCE;}}
	}
	else                       {for (int s=0; s<nspecies; s++){conc[s]=NO_INFLUENCE;}}
}
//------------------------------------------------------------------------------
void   CPointSource::GetSpecifiedMassFlux(const double &t, double *flux) const {
	//returns specified mass flux [M/T]
  if ((t<=endtime) && (t>=starttime)){
		if (type==SPECIFIED_FLUX){for (int s=0; s<nspecies; s++){flux[s]=ConcArray[s];}}	
	  else                     {for (int s=0; s<nspecies; s++){flux[s]=0.0;}}
	}
	else                       {for (int s=0; s<nspecies; s++){flux[s]=0.0;}}
}
/******************************************************************************
             PARSE
******************************************************************************/
CPointSource *CPointSource::Parse(ifstream &input, int &l,const int NumSpecies){
	//string "PointSource" 
	//double x double y double starttime double endtime
	//{C}x[nspecies]
	//&
	CPointSource *pSource;
	int           nspecies(NumSpecies);
  bool          eof(false),done(false);
	int           Len,speccount(0);
	double        startt,endt;
	char         *s[MAXINPUTITEMS];
	double       *initconc=new double  [nspecies];
	cmplex z;

	pSource=NULL;
	if (parserdebug) {cout <<"Point Source"<< endl;}      eof=TokenizeLine(input,s,Len); l++;
  if      (Len==4){z=s_to_c(s[0],s[1]);startt=s_to_d(s[2]);endt=s_to_d(s[3]);}
	else if (Len==0){                                     eof=TokenizeLine(input,s,Len); l++;}
	else            {return NULL;}
	eof=TokenizeLine(input,s,Len); l++;
	do {
		if ((Len==1) && (speccount<nspecies) && (strcmp(s[0],"&"))){
			initconc[speccount]=s_to_d(s[0]); speccount++;    eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (speccount==nspecies) && (!strcmp(s[0],"&"))) {//hit a &
			pSource= new CPointSource(SPECIFIED_CONC, //TMP DEBUG
				                        z,
																initconc,
																startt,
																endt,
																NumSpecies);
			done=true;
		}
		else if (Len==1){ 
			cout <<"Too many species associated with source"<<endl;  return NULL;	
		}
		else if (Len==0){                                   eof=TokenizeLine(input,s,Len); l++;}
	}while ((!done) && (!eof));

	delete [] initconc;
	return pSource;
}
