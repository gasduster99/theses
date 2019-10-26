//LineFlowSink.cpp

#include "FlowSink.h"

/*********************************************************************
                        FLOW SOURCE/SINKS
**********************************************************************
                           CONSTRUCTORS
**********************************************************************/
CSourceSink::CSourceSink(){
	pLayer  =NULL;
	pDomain=NULL;
	concs  =NULL;
}
//--------------------------------------------------------------------
CSourceSink::CSourceSink(const CSingleLayerABC *pLay,
												 const C2DDomainABC    *pDom,
												 const double          *conc){
  pLayer =pLay;
	pDomain=pDom;
	if (pLay==NULL){ExitGracefully("CSourceSink::constructor: bad Layer ",RUNTIME_ERR);}
	if (pDom==NULL){ExitGracefully("CSourceSink::constructor: bad Domain",RUNTIME_ERR);}
	int nspecies=pDomain->GetNumSpecies();
	concs=new double [nspecies];
	for (int s=0;s<nspecies; s++){
		concs[s]=conc[s];
		if (concs[s]<0.0) {ExitGracefully("CSourceSink::constructor: bad concentration specified",BAD_DATA);}
	}

}
//--------------------------------------------------------------------
CSourceSink::~CSourceSink(){
	delete [] concs;
}
//--------------------------------------------------------------------
double  CSourceSink::GetConcentration(const double &t, const int s) const{
	return concs[s];
}
/*********************************************************************
                 POLYLINEAR FLOW SOURCE/SINKS
**********************************************************************
                           CONSTRUCTORS
**********************************************************************/
CLinearSourceSink::CLinearSourceSink():CSourceSink(NULL,NULL,NULL){
	pLeftFace =NULL;
	pRightFace=NULL;
	z1=z2=0.0;
}
//--------------------------------------------------------------------
CLinearSourceSink::CLinearSourceSink(const CSingleLayerABC *pLay,
																		 const C2DDomainABC		 *pDom,
																		 const cmplex						z1s, 
																		 const cmplex						z2s,
																		 const double					 *conc)
							    :CSourceSink(pLay,pDom,conc){
	cmplex offset;


	if (pLayer==NULL){ExitGracefully("CLinearSourceSink: Costructor- null layer",RUNTIME_ERR);}
	z1=z1s;
	z2=z2s;
  pLayer->GetDeltaPot(); //TMP DEBUG
	offset=-(z2-z1)*IM*MOVE_DIST;
	//cout <<"Creating Transect 1"<<endl;
	pLeftFace =new CTransect(-2,pLayer,z2+offset,z1+offset,100);

	//cout <<"Creating Transect 2"<<endl;
	pRightFace=new CTransect(-2,pLayer,z1-offset,z2-offset,100);

}
//--------------------------------------------------------------------
CLinearSourceSink::~CLinearSourceSink(){
	//delete pLeftFace;
	//delete pRightFace; /temporarily not working
}
//--------------------------------------------------------------------
void CLinearSourceSink::GetEndpoints(cmplex &z1s, cmplex &z2s) const{
	z1s=z1;
	z2s=z2;
}
/*********************************************************************
             IS INSIDE
**********************************************************************
//if within some tolerance distance from any segment, 
the point is close enough to the sink to be "inside"
--------------------------------------------------------------------*/
bool CLinearSourceSink::IsInside(const cmplex &z) const{
  cmplex Z=(z-0.5*(z2+z1))/(0.5*(z2-z1));
	if (((Z.imag())<sink_distance) || 
			 (pythag((Z-1.0).real(), Z.imag())<sink_distance) ||
			 (pythag((Z+1.0).real(), Z.imag())<sink_distance)){
		return true;
	}
	return false; 
}
/*********************************************************************
             CAPTURE
**********************************************************************
--------------------------------------------------------------------*/
bool CLinearSourceSink::Capture(const pt3D &pt) const{
	//checks if captured, depends upon height
	if (IsInside(c3Dto2D(pt))) {return true;} //tmp debug-only valid if flux is in both 
	return false; //TMP DEBUG
}
/*********************************************************************
             GET NET MASS FLUX
**********************************************************************
--------------------------------------------------------------------*/
double CLinearSourceSink::GetNetMassFlux(const double &t,const int s,double (*conc)(const cmplex &z, const double &t)) const{
	double net(0.0);
  
	net+=pLeftFace->CalculateMassFlux(t, conc);
	net+=pRightFace->CalculateMassFlux(t, conc);
	
	return net;
}
//--------------------------------------------------------------------
void CLinearSourceSink::CalculateFluxes(const double &t){
	pLeftFace ->CalculateStatistics(t);
	pRightFace->CalculateStatistics(t);
}
//--------------------------------------------------------------------
double CLinearSourceSink::GetInflux(const double &t) const{
  return pLeftFace->GetInflux() +pRightFace->GetInflux();
}
//--------------------------------------------------------------------
double CLinearSourceSink::GetOutflux(const double &t) const{
	return pLeftFace->GetOutflux()+pRightFace->GetOutflux();
}
//--------------------------------------------------------------------
double CLinearSourceSink::GetInfluxThruSegment(const cmplex &z1s, const cmplex &z2s,const double &t) const{
	double influx1,influx2, junk;

	pLeftFace ->GetFluxThroughSegment(z1s,z2s,influx1,junk);
	pRightFace->GetFluxThroughSegment(z1s,z2s,influx2,junk);	

	return influx1+influx2;
}
//--------------------------------------------------------------------
double CLinearSourceSink::GetOutfluxThruSegment(const cmplex &z1s, const cmplex &z2s,const double &t) const{
	double outflux1,outflux2, junk;

	pLeftFace ->GetFluxThroughSegment(z1s,z2s,junk,outflux1);
	pRightFace->GetFluxThroughSegment(z1s,z2s,junk,outflux2);	

	if ((outflux1<0) || (outflux2<0)){ExitGracefully("CLinearSourceSink::GetOutfluxThruSegment: negative outflux through transect", RUNTIME_ERR);}//TMP DEBUG
	return outflux1+outflux2;
}
/*********************************************************************
            PARSE
**********************************************************************
--------------------------------------------------------------------*/
CLinearSourceSink *CLinearSourceSink::Parse(ifstream &input, int &l, const CSingleLayerABC *pLay, C2DDomainABC *pDom){
  //string "LinearSourceSink"
  //double x1 double y1 double x2 double y2
  //&
  //[conc[s]]xnspecies
	//&
	CLinearSourceSink *pSink;
  pSink=NULL;
  bool     eof(false);
	int      Len;
	char   *s[MAXINPUTITEMS];
	cmplex   tmpz1,tmpz2;

	double   tmpc[MAX_SPECIES];
  int      currs(0);

	if (parserdebug) {cout << "Linear Source/Sink"<<endl;}  
																								eof=TokenizeLine(input,s,Len); l++;
	do{
		if (Len==4){
			tmpz1=s_to_c(s[0],s[1]);
			tmpz2=s_to_c(s[2],s[3]);
			                                          eof=TokenizeLine(input,s,Len); l++;
		}
	  else if (Len==0)                           {eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) {break;}
    else                                       {ImproperFormat(s,l); return NULL;}
	} while ((!strcmp(s[0],"&")) && (!eof));
																								eof=TokenizeLine(input,s,Len); l++;
	do{
		if ((Len==1) && (strcmp(s[0],"&")) ){
			if (currs>pDom->GetNumSpecies()){ExitGracefully("CPointSourceSink::Parse: too many concentrations associated with sink", BAD_DATA);}
			tmpc[currs]=s_to_d(s[0]);
			currs++;                                  eof=TokenizeLine(input,s,Len); l++;
		}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) { 		
			if (currs!=pDom->GetNumSpecies()){ExitGracefully("CPointSourceSink::Parse: not enough concentrations associated with sink", BAD_DATA);}
			pSink=new CLinearSourceSink(pLay,
																  pDom,
														      tmpz1,
		                              tmpz2,
																  tmpc);
			break; 
		}
    else                                       {ImproperFormat(s,l); return NULL;}
	} while ((!eof));

	if (eof){return NULL;}
	else    {return pSink;}
}

/*********************************************************************
                 POINT FLOW SOURCE/SINKS
**********************************************************************
                           CONSTRUCTORS
**********************************************************************/
CPointSourceSink::CPointSourceSink():CSourceSink(NULL,NULL,NULL){
	zsink=0;
	Q=0.0;//pTransect=NULL;
	radius=0.0;
}
//--------------------------------------------------------------------
CPointSourceSink::CPointSourceSink(const CSingleLayerABC *pLay,
																	 const C2DDomainABC		 *pDom,
																	 const cmplex						z,
																	 const double						r,
																	 const double					 *conc)
								 :CSourceSink(pLay,pDom,conc){
	zsink=z;
	radius=r;
	Q=0.0;//pTransect=new CCCirTransect(zsink,radius,10);
  
}
//--------------------------------------------------------------------
CPointSourceSink::~CPointSourceSink(){
}
/*********************************************************************
            CALCULATE FLUXES
**********************************************************************
--------------------------------------------------------------------*/
void CPointSourceSink::CalculateFluxes(const double &t){
	cmplex z1,z2,z3,z4;
	z1=zsink+radius;
	z2=zsink+IM*radius;
	z3=zsink-radius;
	z4=zsink-IM*radius;
	Q =pLayer->GetFluxThruFace(z1,z2,t).real();
	Q+=pLayer->GetFluxThruFace(z2,z3,t).real();
	Q+=pLayer->GetFluxThruFace(z3,z4,t).real();
	Q+=pLayer->GetFluxThruFace(z4,z1,t).real();

	//Q is negative if a sink, positive if a source
	//cout << "POINTSOURCEFLUX:" <<zsink<<" "<<radius<<" "<<Q<<endl;

}
/*********************************************************************
            ACCESSORS
**********************************************************************
--------------------------------------------------------------------*/
double  CPointSourceSink::GetOutfluxConcentration(const double &t,const int s ) const{
	return concs[s]; //time ignored for now, sink has constant concentration 
}
//--------------------------------------------------------------------
cmplex CPointSourceSink::GetZ() const{return zsink;}
//--------------------------------------------------------------------
bool CPointSourceSink::IsInside(const cmplex &z) const{return (abs(z-zsink)<=radius);}
//--------------------------------------------------------------------
bool CPointSourceSink::Capture(const pt3D &pt) const{
	if (IsInside(c3Dto2D(pt))) {return true;} //TMP DEBUG, depends upon flux
	return false; 
}
//--------------------------------------------------------------------
double CPointSourceSink::GetNetMassFlux(const double &t,const int s, double (*conc)(const cmplex &z, const double &t)) const{
	//Q is negative if a sink, positive if a source
	if (Q<0.0){return conc(zsink,t)*Q;}
	else      {return concs[s]*Q;     }
}
//--------------------------------------------------------------------
double CPointSourceSink::GetInflux(const double &t) const{
	//Q is negative if a sink, positive if a source
	if (Q<0.0){return -Q;}
	else      {return 0.0;}
}
//--------------------------------------------------------------------
double CPointSourceSink::GetOutflux(const double &t) const{
	//Q is negative if a sink, positive if a source
	if (Q>0.0){return Q;}
	else      {return 0.0;}
}
/*********************************************************************
            PARSE
**********************************************************************
--------------------------------------------------------------------*/
CPointSourceSink *CPointSourceSink::Parse(ifstream &input, int &l, const CSingleLayerABC *pLay, C2DDomainABC *pDom){
  //string "PointSourceSink"
  //double x double y double r 
  //&
  //[conc[s]]xnspecies
	//&
	CPointSourceSink *pSink;
  pSink=NULL;
  bool     eof(false);
	int      Len;
	char    *s[MAXINPUTITEMS];
	cmplex   tmpz;
	double   tmpr;

	double   tmpc[MAX_SPECIES];
  int      currs(0);

	if (parserdebug) {cout << "Point Source/Sink"<<endl;}  
																								eof=TokenizeLine(input,s,Len); l++;
	do{
		if (Len==3){
			tmpz=s_to_c(s[0],s[1]);
			tmpr=s_to_d(s[2]);
			                                          eof=TokenizeLine(input,s,Len); l++;
		}
	  else if (Len==0)                           {eof=TokenizeLine(input,s,Len); l++;}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) {break;}
    else                                       {ImproperFormat(s,l); return NULL;}
	} while ((!strcmp(s[0],"&")) && (!eof));
																								eof=TokenizeLine(input,s,Len); l++;
	do{
		if ((Len==1) && (strcmp(s[0],"&")) ){
			if (currs>pDom->GetNumSpecies()){ExitGracefully("CPointSourceSink::Parse: too many concentrations associated with sink", BAD_DATA);}
			tmpc[currs]=s_to_d(s[0]);
			//cout <<"currs"<<currs<<tmpc[currs]<<endl;
			currs++;                                  eof=TokenizeLine(input,s,Len); l++;
		}
		else if ((Len==1) && (!strcmp(s[0],"&")) ) { 		
			if (currs!=pDom->GetNumSpecies()){ExitGracefully("CPointSourceSink::Parse: not enough concentrations associated with sink", BAD_DATA);}
			pSink=new CPointSourceSink(pLay,
																 pDom,
														     tmpz,
		                             tmpr,
																 tmpc);
			break;
		}
    else                                       {ImproperFormat(s,l); return NULL;}
	} while ((!eof));

	if (eof){return NULL;}
	else    {return pSink;}
}