//AnalysisLocation.h
#ifndef ANALYSISLOC_H
#define ANALYSISLOC_H

#include "AbstractAquifer.h"
#include "AbstractLayer.h"

/****************************************************
 *  Class CAnalysisLocation
 *  Locational Analysis Data Abstraction
 ***************************************************/
class CAnalysisLocation{
 private:/*----------------------------------------------------------*/
	 static int nLocations;
	 static CAnalysisLocation **pLocations;

 public:/*-----------------------------------------------------------*/
	CAnalysisLocation();
  virtual ~CAnalysisLocation();

  virtual void CalculateStatistics(const double &t)=0;
	virtual void WriteOutput        (const double &t) const=0;

	static void CalculateAndWriteAll(const double &t);
	static void DestroyAllAnalysisLocations();
};
/****************************************************
 *  Class CTransect
 *  Transect Data Abstraction
 *  a 2-D plane through which flow statistics are calculated
 ***************************************************/
class CTransect:public CAnalysisLocation{ 
 private:/*----------------------------------------------------------*/
	int										 ID;
	cmplex								 z1;
	cmplex								 z2;
	double								*x;
	double								*Qn;
	double								*Qt;
	int										 nx;
	const CSingleLayerABC *pLayer;
	double					       influx;
	double								 outflux;

  void TrapezoidIntegration(const double &Q1, const double &Q2, 
														const double &x1, const double &x2, 
														double &positive, double &negative);

 public:/*-----------------------------------------------------------*/
	CTransect();
	CTransect(const int my_ID,const CSingleLayerABC *pLay, cmplex zstart, cmplex zend, int numdiv);
 ~CTransect();

	double GetInflux   () const;
	double GetOutflux  () const;
	double GetTotalFlux() const;
  void   GetFluxThroughSegment(const cmplex &zint1, 
		                           const cmplex &zint2, 
														         double &influx, 
														         double &outflux);

	void   Subdivide    (cmplex *zc, Writeable1DArray fluxes, int &ns, const double fluxperdiv);
	double CalculateMassFlux(const double &t, double (*conc)(const cmplex &z, const double &t));

	void   CalculateStatistics(const double &t);
	void   WriteOutput        (const double &t) const;
};
/****************************************************
 *  Class CCirTransect
 *  Circular Transect Data Abstraction
 *  a circular boundary through which flow statistics are calculated
 *  may be used for zone budgets
 ***************************************************/
class CCirTransect:public CAnalysisLocation{ 
 private:/*----------------------------------------------------------*/
	int              ID;
	
	const CSingleLayerABC *pLayer;

	cmplex					 zcen;
	double           R;
	cmplex				  *zctrl;
	int						   ncontrol;
	
	double				  *Qn;
	double				  *Qt;
	double					 influx;
	double					 outflux;
	double           budget;

 public:/*-----------------------------------------------------------*/
	CCirTransect();
	CCirTransect(const int my_ID,const CSingleLayerABC *pLay, cmplex zc, double radius, int numdiv);
 ~CCirTransect();

	//void   Subdivide    (cmplex *zc, Writable1DArray fluxes, int &ns, const double fluxperdiv);
	
	void   CalculateStatistics(const double &t);
	void   WriteOutput        (const double &t) const;
};
/****************************************************
 *  Class CZoneBudget
 *  Zone Budget Polygon Data Abstraction
 *  a polygon through which mass balance statistics are calculated
 ***************************************************/
class CZoneBudget:public CAnalysisLocation{ 
 private:/*----------------------------------------------------------*/
	CTransect **pSides;
	int         nsides;
	int         ID;
	double			budget;
	double			influx;
	double			outflux;

  static int nZB; //total number of Zone budget polygons

 public:/*-----------------------------------------------------------*/
	CZoneBudget();
	CZoneBudget(const CSingleLayerABC *pLay,
		          const cmplex					*zp,
							const int							 pts, 
							const int							 numdiv);
 ~CZoneBudget();

	void   CalculateStatistics(const double &t);
	void   WriteOutput        (const double &t) const;

	double CalculateMassBudget(const double &t, double (*conc)(const cmplex &z, const double &t));
	
	static CZoneBudget *Parse(ifstream &input, int &l,CSingleLayerABC *pLay);
};


/***************************************************/

/****************************************************
 *  Class CHeadObsArray
 *  Observation point collection Data Abstraction
 ***************************************************/
class CHeadObsArray:public CAnalysisLocation{ 
 private:/*----------------------------------------------------------*/
	struct HeadObs{
		pt3D   pt;
		double head;
		double obshead;
		char  *name;
	};

	HeadObs					 **pHeadObservations;
	int							   nHeadObs;
	const CAquiferABC *pAq;

 public:/*-----------------------------------------------------------*/
	CHeadObsArray();
	CHeadObsArray(const CAquiferABC *pAq);
 ~CHeadObsArray();

	void AddObservation     (const cmplex &z, const double &h, char *name);
	
	void CalculateStatistics(const double &t);
	void WriteOutput        (const double &t) const;
};

/****************************************************
 *  Class CBFObsArray
 *  Base flow Observation point collection Data Abstraction
 ***************************************************/
class CBFObsArray:public CAnalysisLocation{ 
 private:/*----------------------------------------------------------*/
	struct BFObs{
		pt3D   pt;
		double flow;
		double obsflow;
		char  *name;
	};

	BFObs					   **pObsArray;
	int							   nObs;
	const CAquiferABC *pAq;

 public:/*-----------------------------------------------------------*/
	CBFObsArray();
	CBFObsArray(const CAquiferABC *pAq);
 ~CBFObsArray();

	void AddObservation     (const cmplex &z, const double &flow, char *name);
	
	void CalculateStatistics(const double &t);
	void WriteOutput        (const double &t) const;
};

/****************************************************
 *  Class CGradObsArray
 *  Gradient Observation point collection Data Abstraction
 ***************************************************/
class CGradObsArray:public CAnalysisLocation{ 
 private:/*----------------------------------------------------------*/
	struct GradObs{
		pt3D   pt;
		double gradient;
		double gradcomponent;
		double obsgradient;
		double angle;//radians
		double obsangle; //radians
		char  *name;
	};

	GradObs					 **pObsArray;
	int							   nObs;
	const CAquiferABC *pAq;

 public:/*-----------------------------------------------------------*/
	CGradObsArray();
	CGradObsArray(const CAquiferABC *pAq);
 ~CGradObsArray();

	void AddObservation     (const cmplex &z, const double &grad, const double &angle, char *name);//angle in radians
	
	void CalculateStatistics(const double &t);
	void WriteOutput        (const double &t) const;
};
	/****************************************************
 *  Class CLakeObsArray
 *  Lake Flux Observation Data Abstraction
 ***************************************************/
class CLakeObsArray:public CAnalysisLocation{ 
 private:/*----------------------------------------------------------*/
	struct LakeObs{
		pt3D   pt;
		double flux;
		double obsflux;
		char  *name;
	};

	LakeObs					 **pObsArray;
	int							   nObs;
	const CAquiferABC *pAq;

 public:/*-----------------------------------------------------------*/
	CLakeObsArray();
	CLakeObsArray(const CAquiferABC *pAq);
 ~CLakeObsArray();

	void AddObservation     (const cmplex &z, const double &flux,  char *name);
	
	void CalculateStatistics(const double &t);
	void WriteOutput        (const double &t) const;
};
#endif