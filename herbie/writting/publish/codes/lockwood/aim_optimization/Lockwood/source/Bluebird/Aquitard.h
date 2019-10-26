#ifndef AQUITARD_H
#define AQUITARD_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "PropertyZone.h"

const int     MAX_LEAK_ORDER=  100;        //maximum order of fourier leakage (one direction)
const int     MAX_LEAK_CTRL =  300;        //maximum fourier leakage ctrl pts (one direction)
const double  domain_buffer=     0.25;     //size of buffer domain around extents

/************************************************************************
  Class CAquiclude
  Analytic Aquiclude Data Abstraction
  parents: CAquicludeABC
***********************************************************************/

class CAquiclude : public CAquicludeABC { 
 private:/*----------------------------------------------------------*/
	
	static int				MinAquicludeIterations;  //minimum local solve iterations
	static int				MaxAquicludeIterations;  //maximum local solve iterations
  static double			AquicludeTolerance;      //maximum acceptable local tolerance

	//Member Variables:
	int               level;

	CLayerSetABC		 *pLayerBelow;  //pointer to layer below 
	CLayerSetABC     *pLayerAbove;  //pointer to layer above

	CAnalyticElem   **pElemArray;   //array of pointers to leakage elements 
	int               nElems;

	CPropZone       **pCondZoneArray;//array of pointers to conductance property zones
	int               nCondZones;

	double            conductance;   //background conductance of leaky layer (zero by default)
  double            elevation;		 //background elevation of leaky layer //noone should care

	double						DeltaLeakage;        //Range of leakage values in domain (Leakmax-Leakmin)
  window						Extents;						 //Extents of modeling domain 

	void							IdentifyDeltaLeakage(const double t);

 public:/*-----------------------------------------------------------*/
  //Constructors
  CAquiclude();
 ~CAquiclude();
  CAquiclude(CLayerSetABC *pLaybelow,CLayerSetABC *pLayabove, const double elev, const int lev);

  //Accessor Functions:
	double          GetConductance (const cmplex &z) const;
	double          GetElevation   (const cmplex &z) const;
	double          GetDeltaLeakage() const;

	CLayerSetABC   *GetLayerAbove	 () const;
	CLayerSetABC   *GetLayerBelow	 () const;

	//Manipulator Functions
	void						UpdateExtents        (const cmplex &z);
	void						UpdateExtents        (const cmplex &z, const double &r);

  //Assignment functions
  void            AddToLayer           (CAnalyticElem *Elemptr);
	void            AddToLayer           (CPropZone     *Propptr);

  //Member Functions 
  cmplex          GetDischargePotential(dir direct, const cmplex &z,const double &t) const;
  cmplex          GetW                 (dir direct, const cmplex &z,const double &t) const;
  double          GetLeakage           (const cmplex &z,const double &t) const;
	double          GetDesiredLeakage    (const cmplex &z,const double &t) const;
  double          GetNetFlux           (const double &t) const;

  void            SolveItself          (double &change, double &objective,const double t);
  void            WriteItself          (ofstream &SOL, const double &t) const;
  void            WriteOutput          (const double &t) const;
};

/************************************************************************
 *  Class CLeakyLayer
 *  Analytic Unbounded Area Sink Element Data Abstraction
 *  parents: CAnalyticElem
 ***********************************************************************/
class CLeakyLayer {
 private:/*----------------------------------------------------------*/
	//Static Variables:
  static int        order;        //order of fourier representation
	static double     fold;         //control points per order
  static int        nLeakCtrl;    //number of control points in one direction (N^0.5)

	static double   **SinTerm;       //2-D saved matrix of sin terms
  static double   **Buffer;        //2-D saved matrix of buffer function 

	//Member Variables:
	CAquicludeABC     *pAquiclude;    

	cmplex            zbl;          //bottom left corner of domain
  double            length;       //width/length of domain

  cmplex          **zctrl;        //pointer to Array of control points
  double          **leakctrl;     //pointer to Array of Leakage at given pts
  int               nleakcontrol; //number of control points (in 1 direction)

  double          **B;            //pointer to array of fourier coeff

 public:/*-----------------------------------------------------------*/
  //Constructors
  CLeakyLayer();
  CLeakyLayer(CAquicludeABC *pAq,int ord,  double OS); 
 ~CLeakyLayer();

	//Static member functions
	static void     Prepare();
	static void     Destroy();
	static void     SetPrecision (const int Precision,int &order, double &fold);

  //Accessor Functions:
	double          GetCond      () const;
	double          GetResistance() const;
	double          GetThickness () const;
	double          GetElevation () const;

  //Assignment functions
	void						SetBounds (window Extents);

  //Member Functions 
  cmplex          GetDischargePotential(const cmplex &z,const double &t) const;
  cmplex          GetW                 (const cmplex &z,const double &t) const;
  double          GetLeakage           (const cmplex &z,const double &t) const;
  double          GetNetDischarge      (const double &t) const;

  void            SolveItself          (double &change, double &objective,const double t);
  void            WriteItself          (ofstream &SOL, const double &t) const;
};
#endif