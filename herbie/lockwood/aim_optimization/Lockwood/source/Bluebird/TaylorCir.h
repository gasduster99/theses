//TaylorCir.h
#ifndef TAYLORCIR_H
#define TAYLORCIR_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"

/************************************************************************
  Class CTaylorCir
  Taylor Series/Element Container Data Abstraction
  For grid, transport & track speedup
************************************************************************/

class CTaylorCir {
 private:/*----------------------------------------------------------*/
	//Member Variables:
	cmplex						zcen;					 //center of taylor expansion
	double						radius;        //radius of expansion
  
	cmplex					 *pTaylorCoeff;  //array of Taylor coeff for outer expansion 
	double            Q;             //net discharge from region

	int				        size;		       //number of contained elements
	CAnalyticElem   **pElemArray;    //an array of pointers to interior Analytic Elements 
	
	cmplex           *zctrl;         //array of control pts (global coords)  
  int               ncircontrol;   //number of control points 
  int               order;         //number of taylor coeff
  double            fold;          //overspecification ratio of ctrl pts

  CSingleLayerABC	 *pLayer;        //pointer to associated layer

	//Private Member Functions
  void              AddToCir        (CAnalyticElem *Elemptr);

 public:/*-----------------------------------------------------------*/
	//Constructors:
	CTaylorCir();
	CTaylorCir(CSingleLayerABC * pLay, cmplex zcenter, double rad, int ord, double OS);
 ~CTaylorCir();

	//Static Variables
	static bool        Active;

	//Static Member Functions
  static CTaylorCir *Parse                (ifstream &input, int &l,CSingleLayerABC *pLay, int defaultprecision);
	static void        SetPrecision         (const int Precision,int &order, double &fold);

	void							 FillCircle           (CAnalyticElem **ElemArray, int numElems);

	//Manipulator Functions
  void               SetTaylor            (const double &t);


	//Member Functions
	cmplex						 GetDischargePotential(const cmplex &z,const double &t) const;
  cmplex             GetW                 (const cmplex &z,const double &t) const;

	bool               IsInside             (const cmplex z) const;  
};
#endif