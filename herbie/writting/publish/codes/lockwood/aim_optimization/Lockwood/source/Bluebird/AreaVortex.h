#ifndef AREAVORTEX_H
#define AREAVORTEX_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"
#include "StringElem.h"
#include "PropertyZone.h"

class CAreaVortex;
/************************************************************************
  Class CAVorDipole
  Analytic Area Vortex Dipole Boundary Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
***********************************************************************/
class CAVorDipole: public CStringElem {
 friend class CAreaVortex;

 private:
	 CAreaVortex *pAreaVortex;

 public:
   CAVorDipole(const CSingleLayerABC   *pLay,                  //Layer of Element
							 const cmplex						 *points,                //Vertices of Polygon (first/last counted twice)
							 const int								NumOfLines,            //Number of lines
							 const int								prec,                  //precision of element
								     CAreaVortex 			 *pAVor);               

  void            SolveItself          (double &change, double &objective,const double t);
	double          GetMaxError          (const double &t) const;
	void            WriteOutput          (const double &t) const;
};
/************************************************************************
  Class CAVorLineVortex
  Analytic Area Sink Linesink Boundary Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
***********************************************************************/
class CAVorLinevortex: public CStringElem {
 friend class CAreaVortex;

 private:
	 CAreaVortex *pAreaVortex;

 public:
   CAVorLinevortex(const CSingleLayerABC  *pLay,                //Layer of Element
									 const cmplex						*points,              //Vertices of Polygon (first/last counted twice)
									 const int							 NumOfLines,          //Number of lines
									 const int							 prec,                //precision of element
									       CAreaVortex			*pAVor);              

  void            SolveItself          (double &change, double &objective,const double t);
	double          GetMaxError          (const double &t) const;
	void            WriteOutput          (const double &t) const;
};

/************************************************************************
  Class CAreaVortex
  Analytic Area Vortex Element Data Abstraction
-------------------------------------------------------------------------
  parents: CAnalyticElem, CStringElem
***********************************************************************/

class CAreaVortex: public CAnalyticElem {
 friend class CAVorLinevortex;
 friend class CAVorDipole;
 friend class CSmoothInhom; //TMP DEBUG

 protected:/*----------------------------------------------------------*/

	//double         **JumpCoeff2;   //jump coefficients of element doublet
	//double         **FarFldCoeff2; //farfield coefficients of element doublet

  double          *curlctrl;     //Array of curl at control pts
  cmplex					*zcurlctrl;    //Array of control Points
  int              ncurlctrl;    //number of control points
  double           ave_curl;     //average curl in vortex-should be zero

  cmplex					*zMQbasis;     //Array of MQ Basis Points
  double          *MQcoeff;      //array of MQ interpolator coeff
  int              MQorder;      //order of MQ interpolator (number of basis points for given areasink) 

  CAVorLinevortex *pLinevortexBoundary;
  CAVorDipole     *pDipoleBoundary;   


	virtual bool     CalcCurl (const double t);    //Calculates curl at control points
																								 //returns true if curl hasn't changed
	virtual void		 SolveMQCoeff(const double t); //solves for Multiquadric coefficients based upon curl @ ctrl pts

	cmplex           GetInteriorPotential (const cmplex &z,const double &t) const;
	cmplex           GetInteriorW         (const cmplex &z,const double &t) const;

 public:/*-------------------------------------------------------------*/
  //Constructors
  //specified curl constuctor (not used)
  CAreaVortex(char									*Name,									
							const CSingleLayerABC *pLay,
							const cmplex					*points,							
							const int							 NumOfLines,
							const double					*curl,
							const cmplex					*givencurlpts,
							const int							 prec,
							const int							 numcurlpts);
  //constructor for subclasses
  CAreaVortex(char									*Name,            //Name of element                  
							const CSingleLayerABC	*pLay,            //Layer of element
							const cmplex					*points,          //vertices of polygon (first/last counted twice)
							const int							 NumOfLines,      //Number of polygon sides
							const int							 prec);           //precision of element
 ~CAreaVortex();

//  static CAreaVortex *Parse            (ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);
	static void     SetPrecision         (const int Precision,int &order, double &fold);
	static void     SetRelaxation        (const int start, const int end);

  void            SetBlockOwner        (COwnerABC *BlockPtr, int seg, int IDinBlock);
	void            UpdateBlock          (const double &t) const;

  //Member Functions (Inherited from CAnalyticElem via StringElem):
  cmplex					GetDischargePotential(const cmplex &z,const double &t) const;
  cmplex					GetW                 (const cmplex &z,const double &t) const;
	//cmplex          GetGx                (const cmplex &z,const double &t) const;
	//cmplex          GetFluxThruFace      (const cmplex &z1,const cmplex &z2, const double &t) const; 
  double          GetCurl              (const cmplex &z,const double &t) const;

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

#endif