//AnalyticSolutions.h

#ifndef ANALYTICSOLS_H
#define ANALYTICSOLS_H

#include "AbstractDomain.h" 
//**************************************************************
class C2DAnalyticSol{

 protected:/*---------------------------------------------------*/
	C2DDomainABC          *pDomain;
	const CSingleLayerABC *pLayer;
	int nspecies;
	double starttime;

 public:/*----------------------------------------------------*/

	C2DAnalyticSol(C2DDomainABC      *pDom){pDomain=pDom;nspecies=0;starttime=0.0;}
  virtual ~C2DAnalyticSol(){}

	virtual void           Initialize (const double startt)=0;	
	virtual double    GetConcentration(const cmplex &z, const int s, const double &t) const=0;
  virtual double             GetMass(const int s, const double &t) const=0;
	virtual double        GetDecayLoss(const int s, const double &t) const=0;
};
//**************************************************************
class C2DPlaneSource:public C2DAnalyticSol{

 private:/*---------------------------------------------------*/


	cmplex z1,z2;  //endpoints of patch
	double alpha;  //diminishing factor of source
	double Co[MAX_SPECIES];     //initial concentration array
	double vx;     //normal velocity
	double vy;     //tangential velocity
	double Dx;     //Dispersion coeff
	double Dy;
  bool neg;

 public:/*----------------------------------------------------*/

	C2DPlaneSource(C2DDomainABC      *pDom,
		             const cmplex	     z1,
								 const cmplex      z2,
								 const double      alpha,
							   const double	    *ConcArray);
 ~C2DPlaneSource();

	void           Initialize (const double startt);
	double    GetConcentration(const cmplex &z, const int s, const double &t) const;
  double             GetMass(const int s, const double &t) const;
	double        GetDecayLoss(const int s, const double &t) const;

	static C2DPlaneSource *Parse(ifstream &INPUT, sourcetype ty, C2DDomainABC *pDom, int &l);
};
//**************************************************************
class C2DPtSource:public C2DAnalyticSol{

 private:/*---------------------------------------------------*/
	cmplex     zc;     //source location
	double     Mo[MAX_SPECIES];    //initial mass of species s or inital mass flow rate (depending)
	sourcetype type;   //constant or initial
	double     angle;  //orientation of flow
	double     v;      //velocity
	double     Dl;     //dispersion coeff
	double     Dt;

 public:/*----------------------------------------------------*/

	C2DPtSource(C2DDomainABC      *pDom,
		          const cmplex	     zcen,
							const double	    *MassArray,
							const sourcetype   ty);
 ~C2DPtSource();

	void           Initialize (const double startt);
	double    GetConcentration(const cmplex &z, const int s, const double &t) const;
  double             GetMass(const int s, const double &t) const;
	double        GetDecayLoss(const int s, const double &t) const;

	static C2DPtSource *Parse(ifstream &INPUT, sourcetype ty, C2DDomainABC *pDom, int &l);
};
#endif