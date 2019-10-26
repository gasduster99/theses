#ifndef AREASINK_H
#define AREASINK_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "PropertyZone.h"

class CAreaSink;
/************************************************************************
  Class CASinkDoublet
  Analytic Area Sink Doublet Boundary Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
***********************************************************************/
class CASinkDoublet: public CStringElem {
	friend class CAreasink;

 private:
	 CAreaSink *pAreaSink;

 public:
   CASinkDoublet(char									 *Name,
		             const CSingleLayerABC *pLay,                  //Layer of Element
								 const cmplex					 *points,                //Vertices of Polygon (first/last counted twice)
								 const int							NumOfLines,            //Number of lines
								 const int							prec,                  //precision of element
								      CAreaSink				 *pASink);               

  void            SolveItself          (double &change, double &objective,const double t);
	double          GetMaxError          (const double &t) const;
	void            WriteOutput          (const double &t) const;
};
/************************************************************************
  Class CASinkLinesink
  Analytic Area Sink Linesink Boundary Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
***********************************************************************/
class CASinkLinesink: public CStringElem {
 friend class CAreaSink;

 private:
	 CAreaSink *pAreaSink;

 public:
   CASinkLinesink(char									*Name,
		              const CSingleLayerABC *pLay,                //Layer of Element
									const cmplex					*points,              //Vertices of Polygon (first/last counted twice)
									const int							 NumOfLines,          //Number of lines
									const int							 prec,
									      CAreaSink				*pASink);             //precision of element

  void            SolveItself          (double &change, double &objective,const double t);
	double          GetMaxError          (const double &t) const;
	void            WriteOutput          (const double &t) const;
};
/************************************************************************
  Class CAreaSink
  Analytic Area Sink Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem
***********************************************************************/

class CAreaSink: public CAnalyticElem {
	friend class CASinkLinesink;
  friend class CASinkDoublet;

 protected:/*----------------------------------------------------------*/

	static int      RelaxStartIter;//Iteration at which easing in of recharge/leakage occurs
	static int      RelaxEndIter;  //Iteration at which easing in is completed- normal solve occurs
	static double   Relax;         //current relaxation of recharge/leakage

  double          * leakctrl;    //pointer to Array of given Leakage at control pts
  cmplex					*zleakctrl;    //pointer to Array of control Points
  int              nleakctrl;    //number of control points

  cmplex					*zMQbasis;     //Array of MQ Basis Points
  double          *MQcoeff;      //array of MQ interpolator coeff
  int              MQorder;      //order of MQ interpolator (number of basis points for given areasink) 

  double           ave_leak;     //average leakage in sink

	CASinkDoublet   *pDoubletBoundary;
	CASinkLinesink  *pLinesinkBoundary;

	virtual bool     CalcLeakage (const double t); //Calculates recharge/leakage at control points
																								 //returns true if leakage hasn't changed
	virtual void		 SolveMQCoeff(const double t); //solves for Multiquadric coefficients based upon leakage @ ctrl pts
	
	cmplex           GetInteriorPotential (const cmplex &z,const double &t) const;

	double           GetDivFluxThruFace   (const cmplex &z1, const cmplex &z2, const double &t) const;

 public:/*-------------------------------------------------------------*/
  //Constructors
  //specified leakage constuctor
  CAreaSink(char									*Name,									//Name of Element							 
						const CSingleLayerABC *pLay,									//Layer of Element
						const cmplex					*points,                //Vertices of Polygon (first/last counted twice)
						const int							 NumOfLines,            //Number of lines
						const double					*leakage,               //specified leakage at numlkpts points
						const cmplex					*givenlkpts,            //numlkpts locations of specified leakage
            const int							 numlkpts,              //precision of element
	          const int							 prec);                 //number of specified leakage points
  //constructor for subclasses
  CAreaSink(char									*Name,                  //Name of element           
					  const CSingleLayerABC	*pLay,                  //Layer of element
						const cmplex					*points,                //Vertices of polygon
						const int							 NumOfLines,            //Number of lines (first/last counted twice)
						const int							 prec);                 //precision of element
 ~CAreaSink();

  static CAreaSink *Parse              (ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);
	static void     SetPrecision         (const int Precision,int &order, double &fold);
	static void     SetRelaxation        (const int start, const int end);

  void            SetBlockOwner        (COwnerABC *BlockPtr, int seg, int IDinBlock);
	void            UpdateBlock          (const double &t) const;

  //Member Functions (Inherited from CAnalyticElem):

	bool            HasFlux              () const;

  cmplex					GetDischargePotential(const cmplex &z,const double &t) const;
  cmplex					GetW                 (const cmplex &z,const double &t) const;
  cmplex					GetGx                (const cmplex &z,const double &t) const;
  cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const double &t) const; 
	double          GetLeakage           (const cmplex &z,const double &t, const leak_type ltype) const; 
  double          GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, const leak_type ltype) const; 
	void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const;

  void            SolveItself          (double &change, double &objective,const double t);
	double          GetMaxError          (const double &t) const;
  void            WriteItself          (ofstream &SOL, const double &t) const;
  bool            ReadItself           (ifstream &SOL);
	void            WriteOutput          (const double &t) const;

  cmplex					Centroid						 () const;
  bool            IsInSquare					 (const cmplex &zc,const double w) const;
  bool            IsInCircle					 (const cmplex &zc,const double r) const;
  bool            PartInCircle				 (const cmplex &zc,const double r) const;
};

/************************************************************************
 *  Class CResLake
 *  Analytic Resistance-specified lake Element Data Abstraction
 *  parents: CAnalyticElem, CStringElem, CAreaSink
 ***********************************************************************/

class CResLake: public CAreaSink {
 private:/*----------------------------------------------------------*/
	bool CalcLeakage();//Calculates recharge/leakage at control points

 public:/*-----------------------------------------------------------*/
	//Constructor
  CResLake(char            *Name,
					 const CSingleLayerABC *pLay,
					 const cmplex    *points,
					 const int        NumOfLines,
					 const double     cond,
					 const double     thick, 
					 const double     depth,
					 const double     elev, 
					 const int        prec,
					 const int        Precision);
 ~CResLake();

	//Member Functions (Inherited from CAnalyticElem via StringElem via AreaSink):
	void            WriteOutput          (const double &t) const;

};

/************************************************************************
 *  Class CConductInhom
 *  Analytic Inhomogeneity ("hole") in conductance Element Data Abstraction
 *  parents: CAnalyticElem, CStringElem, CAreaSink
 ***********************************************************************/

class CConductInhom: public CAreaSink {
 private:/*----------------------------------------------------------*/

	const CAquicludeABC *pAquiclude;      //pointer to aquiclude in which it resides

	double               conduct;         //conductance of inhomogeneity

	CPolyPropZone       *pConductZone;    //conductance property zone

	double               relax;           //relaxation coefficient

	bool         CalcLeakage(const double t);   //Calculates recharge/leakage at control points
  void		    SolveMQCoeff(const double t);   //solves for multiquadric coefficients

 public:/*-----------------------------------------------------------*/
	//Constructor
  CConductInhom(char						*Name, 
								CAquicludeABC		*pAq,
								const CSingleLayerABC *pLay,						//pLayer only needed to calculate Black Hole (TMP DEBUG)
								const cmplex		*points,
								const int				 NumOfLines,
								const double		 cond, 
								const int				 prec,
								const int				 Precision);
	
 ~CConductInhom();

  static CConductInhom *Parse(ifstream  &input,int &l,CAquicludeABC *pAq, CSingleLayerABC *pLay, char *Name);

	//Accessor Functions
 	CPropZone      *GetConductZone() const;

	//Member Functions (Inherited from CAnalyticElem via StringElem via AreaSink):
  void            SolveItself   (double &change, double &objective,const double t);
	double          GetMaxError   (const double &t) const;
	void            WriteOutput   (const double &t) const;

};

#endif