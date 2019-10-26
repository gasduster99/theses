#ifndef POINTELEMS_H
#define POINTELEMS_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"
#include "AnalyticElem.h"

/****************************************************
  Class CPointElem
  Analytic Point Element Data Abstraction
  parents: CAnalyticElem
****************************************************/

class CPointElem : public CAnalyticElem {
 protected:/*----------------------------------------------------------*/
	
	 //Enumerated Types----------------------------------------------------
	enum pointtype {WELL, PTDIPOLE, VORTEX};

	//Member Data----------------------------------------------------------
  cmplex					zp;           //complex coordinate of point
  double          Q;            //strength coefficient
	double          angle;        //orientation (in radians from positive x axis)
 	pointtype       Ptype;        //type of point element  

 public:/*-------------------------------------------------------------*/
  //Constructor:
  CPointElem(); 
  CPointElem(char						 *Name,  //name of element
						 CSingleLayerABC *pLay,  //layer of element
						 pointtype				PType, //point type (well/vortex/dipole)
						 double						ang,   //angle of force (if applicable)
						 double						str,   //strength of element
						 cmplex						z);    //location of point
 ~CPointElem();

  //Accessor Functions                                
  cmplex          GetZ                 () const;     
  double          GetArea              () const;   	
  bool            HasFlux              () const;

	//Manipulator Functions
	void            SetStrength          (double strength);

  //Member Functions (inherited from CAnalyticElem):
  cmplex					GetDischargePotential(const cmplex &z, const double &t) const;           
  cmplex					GetW                 (const cmplex &z, const double &t) const;
	cmplex          GetGx                (const cmplex &z, const double &t) const;
  cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const double &t) const; 
  double          GetLeakage           (const cmplex &z, const double &t,const leak_type ltype) const;
  double          GetNetDischarge      (const double &t) const; 
	void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, const cmplex &z3, const double &t, double &inQ, double &outQ) const;
  void            GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const double &t, double &Q1, double &Q2) const;
  cmplex					GetSingularStrength  (const cmplex &z, const double &t) const;

  void            SolveItself          (double &change, double &objective,const double t)=0;    
	void            WriteOutput          (const double &t) const=0;	
	double          GetMaxError          (const double &t) const=0;

  void            WriteItself          (ofstream &SOL, const double &t) const;
	bool            ReadItself           (ifstream &SOL);

	void            UpdateBlock          (const double &t) const;
	
  cmplex					Centroid             () const;          
  bool            IsInside						 (const cmplex &z) const; 
  bool            IsInSquare           (const cmplex &zc,const double w) const;
  bool            IsInCircle           (const cmplex &zc,const double r) const;
  bool            PartInCircle         (const cmplex &zc,const double r) const;
  bool            SharesNode           (const cmplex &zn) const; 
	void            WriteGeometry        (ofstream &BASEMAP) const;
};         

/****************************************************
 *  Class CDischargeWell
 *  Analytic Discharge Well Element Data Abstraction
 *  parents: CAnalyticElem,CPointElem
 ***************************************************/

class CDischargeWell:public CPointElem {
 private:/*----------------------------------------------------------*/
  
 public:/*-----------------------------------------------------------*/
  //Constructor
  CDischargeWell(char						 *Name,  
		             CSingleLayerABC *pLay, 
								 cmplex						z, 
								 double						prate);                                            
  
	//Static Functions
  static CDischargeWell *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Accessor Functions                                  
  double          GetDischarge() const;

  //Member Functions (inherited from CAnalyticElem via PointElem):
  void            SolveItself (double &change, double &objective,const double t);                                                       
	void            WriteOutput (const double &t) const;
	double          GetMaxError(const double &t) const;
};

/****************************************************
 *  Class CHeadSpecifiedWell
 *  Analytic Head-Specified Well Element Data Abstraction
 *  parents: CAnalyticElem, CPointElem
 ***************************************************/

class CHeadSpecifiedWell : public CPointElem {
 private:/*----------------------------------------------------------*/
  double          radius;            //radius of well
  double          H;                 //Desired head at reference point
  
 public:/*-----------------------------------------------------------*/
  //Constructor
  CHeadSpecifiedWell(char						 *Name, 
		                 CSingleLayerABC *pLay, 
										 cmplex						zw,
										 double						Head,  
										 double						r);     

	//Static Functions
  static CHeadSpecifiedWell *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);
	
  //Accessor Functions
  double          GetRadius() const;                                             
  double          GetHead() const;                                                                                                                      

  //Member Functions (inherited from CAnalyticElem via PointElem):   
  void            SolveItself(double &change, double & objective,const double t);                                            
	void            WriteOutput(const double &t) const;
	double          GetMaxError(const double &t) const;
};

/****************************************************
 *  Class CDryWell
 *  Analytic Discharge-Specified Dryable Extraction Well Element Data Abstraction
 *		-If Well pumps to base of aquifer, it behaves as head specified with head=base
 *  parents: CAnalyticElem, CPointElem
 ***************************************************/

class CDryWell : public CPointElem {
 private:/*----------------------------------------------------------*/
  double          radius;            //radius of well
	double          Qpump;             //desired pumping rate (not necc. actual strength, Q)
  
 public:/*-----------------------------------------------------------*/
  //Constructor
  CDryWell(char						 *Name, 
		       CSingleLayerABC *pLay, 
					 cmplex						zw, 
					 double						prate,  
					 double						r);     

	//Static Functions
  static CDryWell *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);
	
  //Accessor Functions
  double          GetRadius() const;                                                                                                                                                                 
  double          GetDischarge() const;

  //Member Functions (inherited from CAnalyticElem via PointElem):   
  void            SolveItself(double &change, double & objective,const double t);                                            
	void            WriteOutput(const double &t) const;
	double          GetMaxError(const double &t) const;
};
 
/****************************************************
 *  Class CVortex
 *  Analytic Vortex Element Data Abstraction
 *  parents: CAnalyticElem, CPointElem
 ***************************************************/

class CVortex : public CPointElem {
 private:/*----------------------------------------------------------*/
	
 public:/*-----------------------------------------------------------*/
  //Constructor
  CVortex(char						*Name, 
		      CSingleLayerABC *pLay, 
					cmplex					 zv, 
					double					 strength);      

	//Static Functions
  static CVortex *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Accessor Functions                                                                                
  double          GetStrength() const;                                   

  //Member Functions (inherited from CAnalyticElem via PointElem):   
  void            SolveItself (double &change, double & objective,const  double t);
	void            WriteOutput (const double &t) const;
	double          GetMaxError(const double &t) const;
};

/****************************************************
 *  Class CPtDipole
 *  Analytic Point Dipole Element Data Abstraction
 *  parents: CAnalyticElem, CPointElem
 ***************************************************/

class CPtDipole : public CPointElem {
 private:/*----------------------------------------------------------*/

 public:/*-----------------------------------------------------------*/
  //Constructor
  CPtDipole(char						*Name, 
		        CSingleLayerABC *pLay, 
						cmplex					 z, 
						double					 strength,
						double					 orientation);      
	
	//Static Functions
  static CPtDipole *Parse(ifstream &input, int &l,CSingleLayerABC *pLay, char * Name);

  //Accessor Functions                                                                                
  double          GetStrength() const;                                   
	double          GetOrientation() const;

  //Member Functions (inherited from CAnalyticElem via PointElem):   
  void            SolveItself (double &change, double & objective,const  double t);                                            
	void            WriteOutput (const double &t) const;
	double          GetMaxError(const double &t) const;
};

#endif

