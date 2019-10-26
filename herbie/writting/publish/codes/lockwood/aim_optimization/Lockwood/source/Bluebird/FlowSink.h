//FlowSink.h

#include "AnalysisLocation.h"
#include "AbstractDomain.h"

const double  sink_distance    =0.0001;   // local distance to linesink for sink capture

/****************************************************
 *  Class CSourceSink
 *  Source/Sink Data Abstraction
 *  behaves as both source and sink
 ***************************************************/
class CSourceSink{ 
 protected:
	const CSingleLayerABC  *pLayer;
	const C2DDomainABC     *pDomain;

	double *concs;         //concentration of source/sink

 public:
	CSourceSink();
  CSourceSink(const CSingleLayerABC    *pLay,
		          const C2DDomainABC *pDom,
		          const double *conc);
  virtual ~CSourceSink();

	double  GetConcentration(const double &t, const int s) const;

};

/****************************************************
 *  Class CLinearSourceSink
 *  Source/Sink Data Abstraction
 *  Used for capturing particles & contaminant
 *  behaves as both source and sink
 ***************************************************/
class CLinearSourceSink:public CSourceSink{ 
 private:/*----------------------------------------------------------*/

  cmplex            z1,z2;				//line endpoints (global coords) 

  CTransect        *pLeftFace;     //flow transects
	CTransect        *pRightFace;  

 public:/*-----------------------------------------------------------*/
	//Constructors
	CLinearSourceSink();
  CLinearSourceSink(const CSingleLayerABC *pLay,
		                const C2DDomainABC		*pDom,
		                const cmplex					 z1, 
								    const cmplex					 z2,
								    const double					*concs);
 ~CLinearSourceSink();

  //Initialization Functions
	void            CalculateFluxes(const double &t); 
	
	//Accessor functions
	void               GetEndpoints(cmplex &z1s, cmplex &z2s) const;

  double  GetOutfluxConcentration(const double &t, const int s) const;
	double           GetNetMassFlux(const double &t, const int s, double (*conc)(const cmplex &z, const double &t)) const;
  
	double                GetInflux(const double &t) const;
  double               GetOutflux(const double &t) const; 
	double     GetInfluxThruSegment(const cmplex &z1s, const cmplex &z2s,const double &t) const;
	double    GetOutfluxThruSegment(const cmplex &z1s, const cmplex &z2s,const double &t) const;

	//Neccessary?
  bool                   IsInside(const cmplex &z) const;
	bool                    Capture(const pt3D &pt) const;

	//Static Routines
  static CLinearSourceSink *CLinearSourceSink::Parse(ifstream &input, int &l, const CSingleLayerABC *pLay, C2DDomainABC *pDom);
};
/****************************************************
 *  Class CPointSourceSink
 *  Flow Sink Data Abstraction
 *  Used for capturing particles & contaminant 
 ***************************************************/
class CPointSourceSink:public CSourceSink{ 
 private:/*----------------------------------------------------------*/

  cmplex              zsink;				//location of source/sink
	double              radius; 
  double              Q;
	//CCirTransect     *pTransect;

 public:/*-----------------------------------------------------------*/
	//Constructors
	CPointSourceSink();
  CPointSourceSink(const CSingleLayerABC *pLay,
		               const C2DDomainABC		 *pDom,
		               const cmplex						z,
									 const double						r,
								   const double					 *concs);
 ~CPointSourceSink();

  //Initialization Functions
	void            CalculateFluxes(const double &t); 
	
	//Accessor Functions
	cmplex                     GetZ() const;

  double  GetOutfluxConcentration(const double &t, const int s) const;
  double           GetNetMassFlux(const double &t, const int s,double (*conc)(const cmplex &z, const double &t)) const;
	double                GetInflux(const double &t) const;
  double               GetOutflux(const double &t) const; 

	//Neccessary?
  bool                   IsInside(const cmplex &z) const;
	bool                    Capture(const pt3D &pt) const;

	//Static Routines
  static CPointSourceSink *CPointSourceSink::Parse(ifstream &input, int &l, const CSingleLayerABC *pLay, C2DDomainABC *pDom);
};