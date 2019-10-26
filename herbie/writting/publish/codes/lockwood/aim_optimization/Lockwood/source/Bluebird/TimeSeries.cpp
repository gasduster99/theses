//TimeSeries.cpp

#include "MasterInclude.h"
#include "TimeSeries.h"

/**********************************************************************
     CONSTRUCTORS
**********************************************************************/
CTimeSeries::CTimeSeries(){
	starttime=0.0;
	endtime=0.0;
	novalue=0.0;
	pAllTimeSeries[TotalNumTimeSeries]=this;
	TotalNumTimeSeries++;
}
//---------------------------------------------------------------------
CTimeSeries::CTimeSeries(double startt, double endt, double no_val){
	starttime=startt;
	endtime  =endt;
	novalue =no_val;
	pAllTimeSeries[TotalNumTimeSeries]=this;
	TotalNumTimeSeries++;
}
//---------------------------------------------------------------------
CTimeSeries::~CTimeSeries(){}
//---------------------------------------------------------------------
double CTimeSeries::GetValue(const double &t) const{
	return novalue;
}
/**********************************************************************
							STATIC MEMBER FUNCTIONS
**********************************************************************/
int           CTimeSeries::TotalNumTimeSeries=0;
//---------------------------------------------------------------------
CTimeSeries  *CTimeSeries::pAllTimeSeries[];
//---------------------------------------------------------------------
void CTimeSeries::DestroyAllTimeSeries(){
	if (globaldebug){cout <<"DESTROYING ALL TIME SERIES"<<endl;}
	for (int i=0; i<TotalNumTimeSeries; i++){delete pAllTimeSeries[i];}
}
//========================================================================================================
/**********************************************************************
     HEAVISIDE CONSTRUCTORS
**********************************************************************/
CHeavisideTimeSeries::CHeavisideTimeSeries():
	                    CTimeSeries(){
   val=0;
}
//---------------------------------------------------------------------
CHeavisideTimeSeries::CHeavisideTimeSeries(const double    value,
											                     const double    startt,
											                     const double    endt,
																					 const double    novalue):
											CTimeSeries(startt,endt,novalue){
	val=value;						
}
//---------------------------------------------------------------------
CHeavisideTimeSeries::~CHeavisideTimeSeries(){}
/**********************************************************************
     ACCESSORS
**********************************************************************/
double CHeavisideTimeSeries::GetValue(const double &t) const{
	if ((t<=endtime) && (t>=starttime)){return val;}
	return novalue;
}
/**********************************************************************
     STATIC FUNCTIONS
**********************************************************************/
CHeavisideTimeSeries *CHeavisideTimeSeries::Parse(ifstream &input, int &l){
	return NULL;
}
//========================================================================================================
/**********************************************************************
     DISCRETE CONSTRUCTORS
**********************************************************************/
CDiscreteTimeSeries::CDiscreteTimeSeries():CTimeSeries(){
	ts=vs=NULL;N=0;
}
//---------------------------------------------------------------------
CDiscreteTimeSeries::CDiscreteTimeSeries(const double   *times,
										                     const double   *vals,
											                   const int       numvals,
																				 const double    novalue):
										 CTimeSeries(times[0],times[numvals-1],novalue){

	ExitGracefullyIf(times==NULL,"CDiscreteTimeSeries: NULL input(1)",BAD_DATA);
	ExitGracefullyIf(vals ==NULL,"CDiscreteTimeSeries: NULL input(2)",BAD_DATA);
	N=numvals;
	ts=new double [N];
	vs=new double [N];
	ExitGracefullyIf(vs==NULL,"CDiscreteTimeSeries: out of memory",OUT_OF_MEMORY);
	for (int i=0; i<N; i++){
		ts[i]=times[i];
		vs[i]=vals [i];
	}
}
//---------------------------------------------------------------------
CDiscreteTimeSeries::~CDiscreteTimeSeries(){
	delete [] ts;
	delete [] vs;
}
/**********************************************************************
     ACCESSORS
**********************************************************************/
double CDiscreteTimeSeries::GetValue(const double &t) const{
  if ((t<ts[0]) || (t>ts[N-1])){return novalue;}
	for (int i=0; i<N-1; i++){
		if ((t>ts[i]) && (t<=ts[i+1])){
			return LinInterp(vs[i],vs[i+1],ts[i],ts[i+1],t);
		}
	}
	return novalue;
}
/**********************************************************************
     STATIC FUNCTIONS
**********************************************************************/
CDiscreteTimeSeries *CDiscreteTimeSeries::Parse(ifstream &input, int &l){
	//string "DiscTimeSeries" (or variation)
	//{double t double v} x{num time series steps}
	//&
	char      *s[MAXINPUTITEMS];
	double times[MAX_DISCRETE_TIMES];
	double vals [MAX_DISCRETE_TIMES]; 
  int    numvals(0),Len;
  bool   eof(false), done;

	CDiscreteTimeSeries *pTS=NULL;

	done=false;                                          eof=TokenizeLine(input,s,Len); l++;
	do {
		if (Len==2){
			ExitGracefullyIf((numvals>=MAX_DISCRETE_TIMES),"Parse: Too many time series entries",BAD_DATA);
			times[numvals]=s_to_d(s[0]);
			vals [numvals]=s_to_d(s[1]);
			numvals++;                                       eof=TokenizeLine(input,s,Len); l++;
		}
		else if ((Len==1) && (!strcmp(s[0],"&"))) { 
      pTS=new CDiscreteTimeSeries(times,vals,numvals,-1);
			done=true;
		}
		else if (Len==0) {eof=TokenizeLine(input,s,Len); l++;}
		else {ImproperFormat(s,l); break;}
	} while ((!done) && (!eof));	

	return pTS;
}
//========================================================================================================
/**********************************************************************
     POLYNOMIAL CONSTRUCTORS
**********************************************************************/
CPolynomialTimeSeries::CPolynomialTimeSeries(){
	a=NULL; N=0;
}
//---------------------------------------------------------------------
CPolynomialTimeSeries::CPolynomialTimeSeries(const double   *coeff,
																						 const int       numterms,
																						 const double    startt,
																						 const double    endt,
																						 const double    novalue):
                       CTimeSeries(startt,endt,novalue){
	ExitGracefullyIf(a==NULL,"CPolynomialTimeSeries: NULL input",BAD_DATA);
												 
	N=numterms;
	a=new double[N];
	ExitGracefullyIf(a==NULL,"CPolynomialTimeSeries: out of memory",OUT_OF_MEMORY);
	for (int i=0; i<N; i++){
		a[i]=coeff[i];
	}
}
//---------------------------------------------------------------------
CPolynomialTimeSeries::~CPolynomialTimeSeries(){
	delete [] a;
}
/**********************************************************************
     ACCESSORS
**********************************************************************/
double CPolynomialTimeSeries::GetValue(const double &t) const{
	double sum(0.0);
	for (int i=0;i<N; i++){sum+=a[i]*t;} //slow way
	return sum;
}
/**********************************************************************
     STATIC FUNCTIONS
**********************************************************************/
CPolynomialTimeSeries *CPolynomialTimeSeries::Parse(ifstream &input, int &l){
	return NULL;
}
