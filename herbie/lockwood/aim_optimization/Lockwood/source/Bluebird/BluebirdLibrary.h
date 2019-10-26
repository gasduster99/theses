//BluebirdLibrary.h

#ifndef BLUEBIRDLIBRARY_H
#define BLUEBIRDLIBRARY_H

//Includes all files needed for AEM Groundwater/Surface Water Flow Solution
#include "MasterInclude.h"

/************************************************************************  
NAMING AND PROGRAMMING CONVENTIONS IN BLUEBIRD LIBRARY

all classes are named beginning with a capital C
pointers to OBJECTS (not basic data types) and Arrays of OBJECTS start with p

counters- i    loops through lines of a string or elements of an aggregate
							 also columns of a matrix
          j    loops through far field terms, rows of a matrix
					k    inner loop counter
					m    loops through control points along a boundary
					n    loops through terms of a power series
					s	   inner loops through terms of a power series 
					l    loops through lines of a file (for the parsers)
					L    loops through layers

************************************************************************/

//Global Solve and Identification
//************************************************************************

const int     MAXAQLEVELS     =3;       // maximum number of layers (parsing restriction only)
const int     MAXTIMES        =20;      // maximum number of transient times (parsing restriction only)
const int     MAXCONTROLPTS   =500;     // maximum number of ctrl points per element (explicit matrix restriction only)
const int     MAXDOF          =80;      // maximum # of degrees of freedom per element (explicit matrix restriction only)

const double  BHDIST          =500;     // relative dist from centroid to black hole (in net ext. far field)

const double  MOVE_DIST       =0.00001; // relative distance to move from singularities or boundaries 
const double  TEST_ERROR_RATIO=0.25;    // fraction of control points used in error evaluation

const double  MIN_LINEARIZED_POTENTIAL=0.0001; //Only used for stage specified linesink

const int     HEAD_ERROR_TAG  =0;
const int     FLUX_ERROR_TAG  =1;


enum dir       {ABOVE,BELOW};
enum leftright {DIR_LEFT, DIR_RIGHT};
enum leak_type {FROMTOP, FROMBOTTOM, FROMTOPANDBOTTOM};

struct anisotropy{
	double ratio; //kyy/kxx, less than one
	double dir;   //in radians, from 0 to 2*PI
};
//Strings
//************************************************************************
const int    NOT_A_STRING=    -1;       // indicator within an array of string segments 

//Surface Water Interaction
//************************************************************************


//Effective Velocity
//************************************************************************
struct disp_info{double al; double at; double D;};

/********************************************************************************
				Girinskii Potential Conversion Functions
********************************************************************************/

/********************************************************************************
	IsConfined
	 Returns true if local layer condition is confined, false otehrwise
-------------------------------------------------------------------------------*/
inline bool   IsConfined				(const double &potential, 
																 const double &K, 
																 const double &H){
	return (potential>=0.5*K*H*H);
}
/********************************************************************************
	ConvertToHead
	 Converts freshwater potential to LOCAL head (head-base)
   alpha is the ratio of saltwater density to freshwater ~1.035 (1.0 for non-intrusion problems)
   Hsea is the sea level in local head (<=0.0 for non-intrusion problems)
-------------------------------------------------------------------------------*/
inline double ConvertToHead     (const double &potential,
																 const double &K, 
																 const double &H,
																 const double &Hsea,
																 const double &alpha){

	if (K==0){return 0.0;}
	double Hs=max(Hsea,0.0);
	if				(potential>=0.5*K*H*H){									//confined--------------------------
		if			(potential>=K*H*(alpha*Hs-0.5*H)){			//freshwater head        
			return (potential/(K*H)+ 0.5*H);
		}
		else if (potential<=alpha*K*H*(Hs-0.5*H)){			//saltwater head
			return alpha*Hs-(alpha-1.0)*H;
		}
		else {																					//mixed head
			return sqrt(2.0*alpha*(potential-alpha*K*H*(Hs-0.5*H))/K*(1-alpha));
		}
	}
	else {																						//unconfined-------------------------
		if (potential>0.5*K*alpha*alpha*Hs*Hs){					//freshwater head
			return sqrt(2.0*max(potential,0.0)/K);
		}
		else if (potential<0.5*K*alpha*Hs*Hs){					//saltwater head
			return Hs;
		}
		else{																						//mixed head
			return sqrt(2.0*(alpha-1.0)*(max(potential,0.0)-0.5*K*alpha*Hs*Hs)/(K*alpha))+Hs;
		}
	}
}
/********************************************************************************
	ConvertToHead
	 Converts from LOCAL freshwater head to potential
   alpha is the ratio of saltwater density to freshwater ~1.035 (1.0 for non-intrusion problems)
   Hsea is the sea level in local head (<=0.0 for non-intrusion problems)
-------------------------------------------------------------------------------*/
inline double ConvertToPotential(const double &Head, 
																 const double &K, 
																 const double &H,
																 const double &Hsea,
																 const double &alpha){
	double Hs=max(Hsea,0.0);
	if (Head>H){																		//confined--------------------------
		if (Head>=alpha*Hs){													//freshwater Head
			return K*H*(Head-0.5*H);
		}
		else{																					//saltwater Head
			return 0.5*K*(1-alpha)/alpha*pow(Head+(alpha-1)*H-alpha*Hs,2)+alpha*K*H*(Hs-0.5*H); 
		}
	}
	else{																						//unconfined-------------------------
		if (Head>=alpha*Hs){													//freshwater Head
			return 0.5*K*Head*Head;
		}
		else{																					//saltwater Head
			return 0.5*K*(alpha/(alpha-1.0))*pow((Head-Hs),2)+0.5*K*alpha*Hs*Hs;
		}
	}
}
/********************************************************************************
 GetSaltwaterElev
	identifies Local elevation of saltwater/freshwater interface
	input: LOCAL head and LOCAL Sea Elevation
-------------------------------------------------------------------------------*/
inline double GetSaltwaterElev(const double &Head,
															 const double &Hsea,
															 const double &alpha){

	return max(min((1.0/(alpha-1.0))*(alpha*max(Hsea,0.0)-max(Head,0.0)),Head),0.0);
}
/********************************************************************************
 GetSaturatedThick
	identifies local saturated thickness of layer
-------------------------------------------------------------------------------*/
inline double GetSaturatedThick (const double &potential, 
																 const double &K, 
																 const double &H){
	/// need to change for saltwater problems
	if (K==0){return 0.0;}
	return ((potential<0.5*K*H*H) ? sqrt(2.0*max(potential,0.0)/K) : H);
}
/********************************************************************************
 GetTransmissivity
	identifies local transmissivity of layer
-------------------------------------------------------------------------------*/
inline double GetTransmissivity (const double &potential, 
																 const double &K, 
																 const double &H){
	return K*GetSaturatedThick(potential,K,H);
}

#endif