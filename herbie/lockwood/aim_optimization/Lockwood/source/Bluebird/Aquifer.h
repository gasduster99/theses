#ifndef AQUIFER_H
#define AQUIFER_H

#include "BluebirdLibrary.h"
#include "AbstractAquifer.h"
#include "AbstractLayer.h"
#include "AbstractMultiLayer.h"
#include "Layer.h"
#include "MultiLayer.h"
#include "Aquitard.h"

class CSingleLayer;
class CMultiLayer;
class CAquiclude;

/****************************************************
  Class CAquifer
  Aquifer Layer/LeakyLayer Container Data Abstraction
****************************************************/
class CAquifer:public CAquiferABC{  

 private:/*----------------------------------------------------------*/

	static int       MinAqIterations;  //minimum local solve iterations
	static int       MaxAqIterations;  //maximum local solve iterations
  static double    AqTolerance;      //maximum acceptable local tolerance

	char            *name;

	CSingleLayer    *pSLayerArray     [MAXAQLEVELS];
	CMultiLayer     *pMLayerArray     [MAXAQLEVELS];
	CAquiclude      *pAquicludeArray  [MAXAQLEVELS]; //Aquiclude between single or multi- layers
	
	bool             LevelIsMultiLayer[MAXAQLEVELS];

	//CLayerSet      **pLayerArray;//TMP DEBUG: ML MIGRATION to CLayerSet

  int              nLevels;         //Number of single or multilayers separated by aquicludes

 public:/*----------------------------------------------------------*/
  //Constructors:
	CAquifer(char * aqname);
	CAquifer(char * aqname, CSingleLayer **pLayers, CAquiclude **pLeakLays, int numlayers);//TMP DEBUG: ML MIGRATION to CLayerSet
	~CAquifer();
  
  
  static bool      Multilevel;      

  //accessor functions
  double           GetCond     (const cmplex &z, int L)const;
	double           GetCond     (const pt3D   &pt) const;
	double           GetPoro     (const cmplex &z, int L) const;
	double           GetPoro     (const pt3D   &pt) const;

	double           GetThick    ()               const;
  int              GetNumLayers()               const;
	int              GetNumLevels()               const;
	int              GetLevel    (const pt3D pt)  const;
	CSingleLayerABC *GetSLayer   (int L)          const; //TMP DEBUG: ML MIGRATION to CLayerSet
	CMultiLayerABC  *GetMLayer   (int L)          const;
	CAquiclude      *GetAquiclude(int L)          const; 
  CSingleLayer    *GetLayer    (const int l)    const; //TMP DEBUG: ML MIGRATION to CLayerSet - SHOULD REMOVE

  //assignment functions
  static void      SetSolveData(const int miniter, const int maxiter, const double tolerance);

//	void             AddLayerSetToAquifer(CLayerSet *pLayer); 
	void             AddSingleLayer(CSingleLayer *pLayer, CAquiclude *pAquicludeBeneath);
	void             AddMultiLayer (CMultiLayer  *pLayer, CAquiclude *pAquicludeBeneath);

  //member functions
  void             IterativeSolve       (double &t, ofstream &PROGRESS);
	void             IterExplicitSolve    (double &t, ofstream &PROGRESS);

	//inherited member functions
	double  				 GetHead              (const pt3D &pt, const double &t) const;
  double           GetLeakage           (const pt3D &pt,const double &t,const leak_type ltype) const;

  double  				 GetNetDischarge      (const double &t) const; 

	vector           GetVelocity2D        (const pt3D &pt,const double &t) const;
	vector           GetVelocity3D        (const pt3D &pt,const double &t) const;
 	vector           GetEffVelocity2D     (const pt3D &pt,const double &t, const disp_info &disp) const;
 	vector           GetEffVelocity3D     (const pt3D &pt,const double &t, const disp_info &disp) const;

	double           GetSaturatedThickness(const cmplex &z, int L,const double &t) const;
	double           GetSaturatedThickness(const pt3D &pt,const double &t) const;

	double           GetBaseflow          (const cmplex &z, int L, const double t) const;
  double           GetBaseflow          (const pt3D &pt, const double t) const;

	void             InterpolateQxQy      (const double &t);
	void             WriteItself          (ofstream &SOL, const double &t) const;
  void             WriteOutput          (const double &t) const;
};


#endif