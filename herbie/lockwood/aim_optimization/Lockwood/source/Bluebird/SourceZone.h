//SourceZone.h
#ifndef SOURCE_H
#define SOURCE_H

#include "AbstractDomain.h" 
#include "AbstractLayer.h"
#include "PropertyZone.h"
#include "TimeSeries.h"

/****************************************************
 *  Class CAreaSource
 *  2D Chemical Source Zone Data Abstraction
 *    handles: dirichlet source zones
 *             constant concentration recharge zones
 *             specified areal flux zones
 *             initial concentrations 
 ***************************************************/
class CAreaSource{
 protected:/*----------------------------------------*/

	CPropZone         *zone;                //source zone geometry

	int								 nspecies;            //number of species in model
	double            *SorbArray;           //array of sorbed species concentrations (used for init conditions only)

	sourcetype         type;                //type of source (flux, constant conc,etc)
	
	CTimeSeries      **ConcSeries;          //nspecies time series for concentration profile

	//TMP DEBUG- should be able to handle any time series representation
	double            *ConcArray;           //array of species concentrations (or fluxes)
	double             starttime;           //starttime of source
	double             endtime;             //endtime of source 

 public:/*----------------------------------------*/
	//Constructors
	CAreaSource(sourcetype         ty,
							double						*conc,
							double            *sorb,
							double             startt,
							double             endt,
							const int          NumSpecies);
	CAreaSource(sourcetype         ty,
							CTimeSeries      **TimeSeries,
							double            *sorb,
							const int          NumSpecies);
  virtual ~CAreaSource();

	//manipulators
	void       SetTimeSeries(CTimeSeries **TimeSeries, const int numseries);

	//Accessors
  sourcetype         GetType() const; //TMP DEBUG - neccesary?
  double             GetArea() const; //TMP DEBUG - neccesary?



	void     GetSpecifiedConcentrations(const cmplex &z, const double &t, double *conc) const;
  void           GetSpecifiedMassFlux(const cmplex &z, const double &t, double *flux) const; 
  void      GetRechargeConcentrations(const cmplex &z, const double &t, double *conc) const;
	void GetInitialSorbedConcentrations(const cmplex &z, const double &t, double *sorb) const;

	//Member Functions
	bool              IsInside(const cmplex z) const;
};
/****************************************************
 *  Class CPolyAreaSource
 *  2D Polygonal Chemical Source Zone Data Abstraction
 ***************************************************/
class CPolyAreaSource:public CAreaSource{

public:
	CPolyAreaSource(cmplex						*points, 
									int                npoints, 
									sourcetype         ty,
									double						*conc,
									double						*sorb,
									double             startt,
									double             endt,
									const int          NumSpecies);
	CPolyAreaSource(cmplex						*points, 
									int                npoints, 
									sourcetype         ty,
									CTimeSeries			 **timeseries,
									double						*sorb,
									const int          NumSpecies);
	//Static Functions
	static CPolyAreaSource *Parse(ifstream &input,sourcetype ty, int &l,const int NumSpecies);
	
};
/****************************************************
 *  Class CEllAreaSource
 *  2D Elliptical Chemical Source Zone Data Abstraction
 ***************************************************/
class CEllAreaSource:public CAreaSource{

public:
	CEllAreaSource(const cmplex		  zcen,        //center of ellipse
								 const double     MajorAxis,   //major axis of ellipse
								 const double		  MinorAxis,   //minor axis of ellipse
								 const double     angle,			 //angle (in radians) 
								 sourcetype       ty,
								 double					 *conc, 
								 double          *sorb,
								 double           startt,
								 double           endt,
								 const int        NumSpecies);

	//Static Functions
	static CEllAreaSource *Parse(ifstream &input,sourcetype ty, int &l,const int NumSpecies);
	
};
/****************************************************
 *  Class CPointSource
 *  Chemical Point Source Data Abstraction
 ***************************************************/
class CPointSource{
 private:/*----------------------------------------*/

	sourcetype    type;
	cmplex        zs;                 //location

	double       *ConcArray;          //array of species concentrations
	int						nspecies;           //number of species in model

	double        starttime;          //starttime of source
	double        endtime;            //endtime of source 

 public:/*-----------------------------------------*/
	//Constructors
	CPointSource(sourcetype  ty,
		           cmplex			 z,  
							 double			*initconc,
							 double      startt,
							 double      endt,
							 const int NumSpecies);
 ~CPointSource();

	//Static Functions
	static CPointSource *Parse(ifstream &input, int &l,const int NumSpecies);

	//Accessors
  cmplex                GetLocation() const;
  void         GetSpecifiedMassFlux(const double &t, double *flux) const;
	void   GetSpecifiedConcentrations(const double &t, double *conc) const;

};
/****************************************************
 *  Class CLineSource
 *  Chemical Line Source Data Abstraction
 ***************************************************/
class CLineSource{
 private:/*----------------------------------------*/

	cmplex        z1;                 
  int           z2;              

	double       *ConcArray;          //array of species concentrations

	int						nspecies;           //number of species in model
	double        starttime;          //starttime of source
	double        endtime;            //endtime of source 

 public:/*-----------------------------------------*/
	//Constructors
	/*CLineSource();
	CLineSource(cmplex	z1,
		          cmplex  z2,
							double *initconc,
							double  startt,
							double  endt);
 ~CLineSource();

	//Static Functions
	static CLineSource *Parse(ifstream &input, int &l);

	//Accessors
  cmplex         GetClosestX(const cmplex &z) const;
  void               GetFlux(const double &X, const double &t, double *flux) const;
	bool     GetConcentrations(const double &X, const double &t, double *conc) const;
*/
};
#endif