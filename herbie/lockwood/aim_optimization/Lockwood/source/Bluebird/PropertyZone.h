//PropertyZone.h
#ifndef PROPZONE_H
#define PROPZONE_H

#include "MasterInclude.h"

//Property zone parameters---------------------------------------------
const double CLOSE_TO_LINE=  1.02;      // ellipse radius around line which determines closeness
const double NO_VALUE=-1.714893432;     // random unlikely number for property value
const int    MAXPZLINES=200;            //maximum number of sides for polypropzone
const int    MAXPZPTS=400;
const int    MAX_PROPZONES=20500;       //maximum number of property zones
 
enum proptype{notype,
							hydraulic_conductivity, 
							poro, 
							aq_hydraulic_conductivity, 
							aquitard_poro,
							base_elevation,
							hydraulic_conductance,
							layer_thickness,
							aquitard_thickness,
							soil_type,
							source_zone};

/****************************************************
 *  Class CPropZone
 *  Polygonal & Circular Property zone Data Abstraction
 *  used as area of constant (double) property value
 ***************************************************/

class CPropZone {
 protected:/*----------------------------------------------------------*/

	static CPropZone  *pAllPropZones[MAX_PROPZONES];
	static int         TotalPropZones;

	//Member Variables
	double           value;              //parameter value in zone	
	proptype         ptype;              //property type
	double           area;               //area of zone

 public:/*-------------------------------------------------------------*/
	 //Constructors
   CPropZone();
	 CPropZone(const proptype typ,
		         const double   val);
   virtual ~CPropZone();

   static void   DestroyAllPropZones();
	 static double NestedSift      (const cmplex &z, CPropZone **pZones, const int nzones, const double back_val);
	 static cmplex NestedSiftSlopes(const cmplex &z, CPropZone **pZones, const int nzones, const cmplex back_val);

	 //accessor functions
	 virtual double   GetValue        (const cmplex &z) const;
	 virtual cmplex   GetSlopes       (const cmplex &z) const;
	 double           GetValue        () const;
	 proptype         GetType         () const;
   double           GetArea         () const;

	 //Virtual Member Functions
	 virtual bool   IsInside  (const cmplex &z) const;
   virtual bool   IsInSquare(const cmplex &zc,const double w) const;
	
};

/****************************************************
 *  Class CPolyPropZone
 *  Polygonal Property zone Data Abstraction
 ***************************************************/
class CPolyPropZone: public CPropZone {
 private:/*----------------------------------------------------------*/
	//Member Variables
	cmplex          *ze;						//array of polygon endpoints (global coords)
  int              NLines;				//number of line segments in string
	window           boundingbox;		//bounding rectangle

 public:/*-----------------------------------------------------------*/
	//Constructors
	CPolyPropZone();
	CPolyPropZone(const proptype  typ,
		            const cmplex   *points,
								const int       NumOfLines,
								const double    val);
  ~CPolyPropZone();

	//Static Member Functions
	static CPropZone *Parse(ifstream &input, int &l,proptype typ);

	//Accessor Functions
	double   GetValue      (const cmplex &z) const;
	cmplex   GetSlopes     (const cmplex &z) const;
	double   Area          (bool &clockwise) const;
	void     GetPoly       (cmplex *points, int &npoints);
	window  *GetBoundingBox() const;

	//Member Functions (inherited from CPropzone)
	bool   IsInside        (const cmplex &z) const;
  bool   IsInSquare      (const cmplex &zc,const double w) const;
};

/****************************************************
 *  Class CMQPolyPropZone
 *  Polygonal Multiquadric Property zone Data Abstraction
 *  has smoothly varying values in polygon
 ***************************************************/
class CMQPolyPropZone: public CPolyPropZone {
 private:/*----------------------------------------------------------*/
	//Member Variables	
	double    *values;
	cmplex    *Zctrl;
	double    *MQcoeff;
	int        nbasispts;
	double     radius;

	double     anisotropy;
	double     anis_angle;
	cmplex     zref;


 public:/*-----------------------------------------------------------*/
	//Constructors
	CMQPolyPropZone();
	CMQPolyPropZone(const proptype  typ,
									const cmplex   *vertices,
									const int       NumOfLines,
									const cmplex   *basispts,
									const double   *values,
									const int       nbasis,
									const double    R, 
									const double    anis_ratio,
									const double    anis_ang);
  ~CMQPolyPropZone();

	//Static Member Functions
	static CMQPolyPropZone *Parse(ifstream &input, int &l,proptype typ);

	//Accessor Functions (inherited from CPropzone)
	double   GetValue      (const cmplex &z) const;
	cmplex   GetSlopes     (const cmplex &z) const;
};

/****************************************************
 *  Class CLinPolyPropZone
 *  Polygonal Linearly varying Property zone Data Abstraction
 *  has linearly varying values in polygon
 ***************************************************/
class CLinPolyPropZone: public CPolyPropZone {
 private:/*----------------------------------------------------------*/
	//Member Variables	
	cmplex    zcen;
	cmplex    slope;
	double    center_value;

 public:/*-----------------------------------------------------------*/
	//Constructors
	CLinPolyPropZone();
	CLinPolyPropZone(const proptype  typ,
									 const cmplex   *vertices,
									 const int       NumOfLines,
									 const cmplex    center,
									 const double    valcen,
									 const cmplex    slopes);
  ~CLinPolyPropZone();

	//Static Member Functions
	static CLinPolyPropZone *Parse(ifstream &input, int &l,proptype typ);

	//Accessor Functions (inherited from CPropzone)
	double   GetValue      (const cmplex &z) const;
	cmplex   GetSlopes     (const cmplex &z) const;
};

/****************************************************
 *  Class CCirPropZone
 *  Circular Property zone Data Abstraction
 ***************************************************/
class CCirPropZone: public CPropZone {
 private:/*----------------------------------------------------------*/
	//Member Variables
	cmplex           zcen;         //center of circle
	double           radius;       //radius of circle
 
 public:/*-----------------------------------------------------------*/
	//Constructors
	CCirPropZone(const proptype typ,
		           const cmplex   zc,
							 const double   rad,
							 const double   val);
  ~CCirPropZone();

	//Member Functions (inherited from CPropzone)
	bool   IsInside       (const cmplex &z) const;
  bool   IsInSquare     (const cmplex &zc,const double w) const;
};

/****************************************************
 *  Class CEllPropZone
 *  Elliptical Property zone Data Abstraction
 ***************************************************/
class CEllPropZone: public CPropZone {
 private:/*----------------------------------------------------------*/
	//Member Variables
	cmplex           zc;           //center of ellipse
	double           a;            //long axis
	double           b;            //short axis
 	double           angle;        //orientation (radians)

 public:/*-----------------------------------------------------------*/
	//Constructors
	CEllPropZone(const proptype typ,
		           const cmplex   zcen,
							 const double   a,
							 const double   b,
							 const double   orient,
							 const double   val);
  ~CEllPropZone();

	//Member Functions (inherited from CPropzone)
	bool   IsInside       (const cmplex &z) const;
  bool   IsInSquare     (const cmplex &zc,const double w) const;

	//Static Member Functions
	static CEllPropZone  *CEllPropZone::Parse(ifstream &input,int &l,proptype type);
};
#endif