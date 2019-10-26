//GeneralInclude.h
//Basic types definitions, constants, and simple inline functions
//Used by Bluebird/Cardinal, Matrix Library, and Reaction Library
#ifndef GENERALINCLUDE_H
#define GENERALINCLUDE_H

#include <math.h>
#include <iostream>
#include <fstream>
#include<stdlib.h>
#include<string.h>

using namespace std;

//type definitions-------------------
typedef const double *               Unchangeable1DArray; //unchangeable but movable
typedef       double * const         Writeable1DArray;    //unmoveable but changeble
typedef const double * const         Ironclad1DArray;     //unmodifiable

typedef const double * const *       Unchangeable2DArray;
typedef       double *       * const Writeable2DArray;
typedef const double * const * const Ironclad2DArray;

//Global Constants 
//************************************************************************
enum  badcode {//error codes
	NO_ERROR,
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
}; 
const bool globaldebug=false;		//turn to true to track object destruction

//Strings/Parser
//************************************************************************
const int    MAXINPUTITEMS  =    256;  // maximum delimited input items per line
const int    MAXCHARINLINE  =   6000;  // maximum characters in line
const int    FILENAME_SIZE  =    256;  // characters in filename

//Mathematical
//************************************************************************
const double PI          =3.141592653589793238462643;      // pi 
const double EULER       =0.577215664901533;               //Euler's constant
const double LN_10       =2.30258509299;                   //ln(10)
const double REALSMALL   =1e-10;                           //almost zero
const double ALMOST_INF  =1e99;			                    	 //almost infinite or NaN number
const double MAXEXP      =300.0;                           //too large of an exponent for exp() function

/*--------------------------------------------------------------------------------
  Inline functions (simple)
--------------------------------------------------------------------------------*/
inline double        max    (double r, double l)  {return ((r>=l) ? r : l);}
inline double        min    (double r, double l)  {return ((r<=l) ? r : l);}
inline int           max    (int    r, int    l)  {return ((r>=l) ? r : l);}
inline int           min    (int    r, int    l)  {return ((r<=l) ? r : l);}
inline bool          oppsign(double r, double l)  {return ((r*l) < 0.0);   }
inline int           ipow   (int    l, int    r)  {return (int)(pow((double)(l),r)); }
inline int           iabs   (int i1)              {return (int)(fabs((double)(i1))); }
inline double        sqrt   (int i1)              {return sqrt((double)(i1));}
inline double        exp    (int i1)              {return exp((double)(i1));}
inline int					 s_to_i (char *s1)            {return (int)atof(s1);   }
inline double				 s_to_d (char *s1)            {return atof(s1);        }
inline bool				NotANumber(const double x)      {return ((x!=0.0) && (x!=1.0) && ((x*x)==x));}
inline void        lowerswap(double &u,const double v){if (v<u){u=v;}}
inline void        upperswap(double &u,const double v){if (v>u){u=v;}}
/********************************************************************************
	LinInterp
	 provides interpolated value between v1 and v2 
--------------------------------------------------------------------------------*/
inline double      LinInterp(const double v1,
														 const double v2, 
														 const double x1, 
														 const double x2, 
														 const double xint)   {
	if ((x2-x1)!=0.0){return (v1+(v2-v1)*(xint-x1)/(x2-x1));}
	else             {return v1;}
}
inline double    pythag(const double a,const double b){
	//not subject to overflow errs c=sqrt(a^2+b^2)
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa>absb){return absa*sqrt(1.0+(absb/absa)*(absb/absa));}
	else          {return (absb==0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));}
}
//string manipulation functions (found int ParseRoutines.cpp)-------------------

bool   TokenizeLine (ifstream &BBD, char **out, int &numwords);
void   ImproperFormat(char **s, int l);
void   WriteEllipsisToScreen(const int i,const int N, const int modN);

#endif
