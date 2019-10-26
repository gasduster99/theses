#ifndef MASTERINCLUDE_H
#define MASTERINCLUDE_H

/***********************************************************************************
This file includes declarations of common functions (mathematical/string functions)
Definitions are in:
   CommonFunctions.cpp
   MatrixSolvers.cpp
**********************************************************************************/

//neccesary libraries-----------------
#include <math.h>
#include <iostream>
#include <complex>
#include <string>
#include <stdio.h>
#include <typeinfo>
#include <strstream>
#include <time.h>

#include "GeneralInclude.h"
#include "MyVector.h"

//type definitions-------------------
typedef complex<double>              cmplex;

typedef const cmplex *               Unchangeable1DArray_z;
typedef       cmplex * const         Writeable1DArray_z;
typedef const cmplex * const         Ironclad1DArray_z;

typedef CVector                      pt3D;
typedef CVector                      vector;

const bool parserdebug=false;   //turn to true for debugging of parser
const bool really_noisy=false;  //turn to true for debugging of parser/tokenizeline
const int  MAX_MLAYERS=20;

const pt3D blankpoint;
struct window{double n,e,w,s;};
struct polygon{
	cmplex *vertices;
	int     nsides;
};
//Global Constants 
//************************************************************************
/*enum  badcode {//exit codes
	COMPLETED,
	RUNTIME_ERR,
	TOO_MANY,
	BAD_DATA,
	BAD_SOL_READ,
	SINGMAT,
	OUT_OF_MEMORY,
	NO_INPUT_FILE,
	VIRTUAL_ERROR,
	MESH_GEN_ERR,
	OTHER
};*/ 
enum  intercepttype {
	INTERSECTED,
	INTNODE,
	NO_INTERSECT,
	TWO_INTERSECTIONS
};

const int    MAXLINES       =    200;  // maximum number of lines in polylinear string (parsing limitation only)
   
//Mathematical
//************************************************************************
const double IOVER2PI    =0.15915494309189535;             //1/(2*pi)
const cmplex IOVER2PII   =cmplex(0,-0.15915494309189535);  //1/(2*pi*i) 
const cmplex IM          =cmplex(0,1);                     // i
const cmplex ZERO        =cmplex(0,0);                     // 0
const cmplex E1          =cmplex(1,0);                     //right end of line elem in Z coord
const cmplex E2          =cmplex(-1,0);                    //left end of line elem in Z coord
const double SURFERMAGIC =1.7014100091878E+038;            //no-grid in surfer
const double BIGGERIZE   =0.05;                            //size increase factor for bounding box/window
const double NEAR_FEATURE=0.01;                            //close to a feature in dimensionless coordinates

const int    MAXGAUSSPTS =20;                              //maximum number of gaussian integration points
const int    MAX_MQ_ORDER=45;
const int    DEFAULT_DIVS=200;                             //default number of divisions for numerical integration
//DXF format colors
//************************************************************************
const int		 DXF_RED=1;
const int		 DXF_YELLOW=2;
const int		 DXF_GREEN=3;
const int		 DXF_CYAN=4;
const int		 DXF_BLUE=5;
const int		 DXF_MAGENTA=6;
const int		 DXF_WHITE=7;
const int		 DXF_DKGRAY=8;
const int		 DXF_LTGRAY=9;
const int		 DXF_BROWN=26;
const int		 DXF_ORANGE=30;
const int		 DXF_VIOLET=202;
/************************************************************************

                  GLOBALLY ACCESSIBLE FUNCTIONS

************************************************************************/

/*--------------------------------------------------------------------------------
  Inline complex functions (simple)
--------------------------------------------------------------------------------*/
inline cmplex				 s_to_c (char *s1, char *s2)  {return cmplex (atof(s1),atof(s2));}
inline double        adjarg (const cmplex &z)     {double a=arg(z); if (a<0.0){a=2.0*PI+a;} return a;}
/********************************************************************************
	Zm1oZp1
	 optimized version of (Z-1.0)/(Z+1.0): 160% faster than direct call
-------------------------------------------------------------------------------*/
inline cmplex Zm1oZp1(const cmplex &Z){                            
	double den(Z.real()*Z.real()+Z.imag()*Z.imag()+Z.real()+Z.real()+1.0);
  return cmplex((den-Z.real()-Z.real()-2.0)/den,
								(Z.imag()+Z.imag())        /den);
}
/********************************************************************************
	c3Dto2D/c2Dto3D
	 conversion from 2D to 3D coords
-------------------------------------------------------------------------------*/
inline cmplex c3Dto2D  (const pt3D &pt) {return cmplex(pt.x,pt.y);}
inline pt3D   c2Dto3D  (const cmplex &z){return pt3D(z.real(),z.imag(),0.0);}

/********************************************************************************

  Generic Function declarations (functions located in CommonFunctions.cpp)

********************************************************************************/

//program exit strategy functions-------------------------------------------------

void   ExitGracefully(char *statement, badcode code);
void   ExitGracefullyIf(bool condition, char *statement, badcode code);
bool   ProgramAborted();
void   GlobalDelete  (bool pause);
//void   ImproperFormat(char **s, int l);

//array manipulation functions----------------------------------------------------

bool   DynArrayAppend(void **& pArr,void *xptr,int &size);
void   StraightSort  (Writeable1DArray a, const int size); 

//string manipulation functions---------------------------------------------------

//bool   TokenizeLine (ifstream &BBD, char **out, int &numwords);

//input/output functions----------------------------------------------------------

void   line_dxf			(ofstream &DXF,
										 double x1, double y1, double z1, 
										 double x2, double y2, double z2, 
										 char *layer,int color);

//complex power series functions--------------------------------------------------
cmplex Laurent     (const cmplex &Z, const int N, const cmplex * const &cf);
cmplex LaurentDer  (const cmplex &Z, const int N, const cmplex * const &cf, const double r);
cmplex LaurentDxx  (const cmplex &Z, const int N, const cmplex * const &cf, const double r);
cmplex LaurentRe   (const cmplex &Z, const int N, const double * const &cf);
cmplex LaurentDerRe(const cmplex &Z, const int N, const double * const &cf);
cmplex LaurentDxxRe(const cmplex &Z, const int N, const double * const &cf);
cmplex Taylor      (const cmplex &Z, const int N, const cmplex * const &cf);
cmplex TaylorDer   (const cmplex &Z, const int N, const cmplex * const &cf,const double &r);
cmplex TaylorDxx   (const cmplex &Z, const int N, const cmplex * const &cf,const double &r);

//random number functions---------------------------------------------------------

#define IA 16807               // Defined for random number generators-move to common functions ?? 
#define IMI 2147483647
#define AM (1.0/IMI)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IMI-1)/NTAB)
#define EPS 1.2e-7
#define RNMAX (1.0-EPS)
static long *idum=NULL;							

void   randinit   (long idum);
double ran1       ();
double gaussrand  ();
double quickrandom(const double low,const double high);

//geometric functions---------------------------------------------------------

bool					InSquare  (const cmplex &z,  const cmplex &squarecen, const double w);
bool					InEllipse (const cmplex &z,  const cmplex &z1, const cmplex &z2, const double dist);


intercepttype Intersect (const cmplex &z1a,const cmplex &z1b,const cmplex &z2a,const cmplex &z2b,cmplex &zint);
intercepttype LineIntersectCircle(const cmplex &z1,    const cmplex &z2, 
																	const cmplex &zcen,  const double  R, 
																	      cmplex &zint1,       cmplex &zint2);
double        IntegrateLine(const double &x1, const double &x2,
										        const double &y1, const double &y2, 
										        double &positive, double &negative);

//Basic Mathematical functions------------------------------------------------------

cmplex        arccosh  (const cmplex &z);
double        my_erf	 (const double &x);
double        my_erfc	 (const double &x);
double        hantush  (const double &u, const double b);
double        bessel_I0(const double &x);
double        bessel_K0(const double &x);
double        expint   (const double &x);
double    MQInterpolate(Ironclad1DArray   valctrl,
												Ironclad1DArray_z zctrl,
												int               numctrl,
												Writeable1DArray  MQcoeff,
												Ironclad1DArray_z zbasis,
												int               order, 
												double            R,
												double           &ave_val);
bool            PolyFit(Ironclad1DArray   x,
												Ironclad1DArray   y,
												Writeable1DArray  a,
												const int         N);
double         PolyEval(const double      x,
												Ironclad1DArray   a,
												const int         N);
//--------------------------------------------------------------------------
// Linear integration functions
enum lin_integrationtype{LG_1POINT, LG_2POINT,
                         LG_3POINT, LG_4POINT, LG_5POINT, LG_SIMPSONS, 
												 LG_TRAPEZOID1, LG_TRAPEZOID2, LG_TRAPEZOID3,LG_TRAPEZOID10,LG_TRAPEZOID100};
int GetNumLinIntegrationPts(const lin_integrationtype gausspts);

void GetLinIntegrationInfo(const double &x1,const double &x2,
											    	const lin_integrationtype gausspts, Writeable1DArray x, Writeable1DArray w);
//--------------------------------------------------------------------------
// Triangle functions, triangle integration functions
enum  trigausspoints{TG_1POINT,
                     TG_3POINT, TG_3POINTB, TG_3ENDPOINT, 
										 TG_4POINT, TG_6POINT,  TG_7POINT,
										 TG_9POINT, TG_13POINT, TG_16POINT};

inline cmplex TriLocalToGlobal(const cmplex &z1 , const cmplex &z2,  const cmplex &z3, 
												       const double &xii, const double &xij, const double &xik){
	return xii*z1 + xij*z2 + xik*z3;
}
void   TriGlobalToLocal(const cmplex &z1 , const cmplex &z2,  const cmplex &z3,
												const cmplex &z, 
															double &xii,       double &xij,       double &xik);
double				TriArea  (const cmplex &z1, const cmplex &z2, const cmplex &z3);
bool         InTriangle(const cmplex &z,  
												const cmplex &z1, const cmplex &z2, const cmplex &z3);

int   GetNumTriGaussPts(const trigausspoints gausspts);
void    GetTriGaussInfo(const cmplex &z1,const cmplex &z2,const cmplex &z3,
												const trigausspoints gausspts,cmplex *z, double *w);


#endif
