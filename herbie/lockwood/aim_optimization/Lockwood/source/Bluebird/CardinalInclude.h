#ifndef CARDINCLUDE_H
#define CARDINCLUDE_H

#include "MasterInclude.h"

/************************************************************************  
NAMING AND PROGRAMMING CONVENTIONS:

all classes are named beginning with a capital C

counters- i    loops through columns of the cell matrix
          j    loops through rows of the cell matrix
					k    loops through source zones
					m    loops through particles, pathlines
					s	   loops through species
					l    loops through lines of a file (for the parsers)

*///**********************************************************************

//Global Constants
//************************************************************************
const int			MAX_OBS_POINTS   =10;       // maximum number of observation points

//General Transport-------------------------------------
const int     MAX_SIA_ITER       =10;     //maximum sequential iterative approach iterations
const int     MAX_SOURCES_PER_CELL=3;     //maximum number of sources allowed in cell
const int     MAX_CONC_TIMES   =100;      //parsing limitation only

//Streamline approach-----------------------------------
const int     MAX_CELLS_JUMPED =200;
const int     MAX_CROSSES_OF_CELL=4;      //maximum times a cell boundary can be crossed by 1 streamline
const int     MAX_LINES_IN_CELL=50;       //maximum streamlines in a reaction cell (accounts for wells, hopefully)
const int     MAX_STREAMLINE_FACES=20;    //maximum number of streamline starting/injection faces



//MOC/MMOC----------------------------------------------
const double	MOC_ALLOC_TOLERANCE   =0.005;

//Sources-----------------------------------------------
const double  NO_INFLUENCE      =-0.0001;  //source concentration does not apply (either because of time or species)

const double  DEFAULT_DRY_DENSITY=1.0;


enum advtype    {  //advective stepping method
	CONSTANT_TIME_STEP, 
  ADAPTIVE_TIME_STEP,
	CONSTANT_SPACE_STEP
};
enum disptype   {  //dispersion representation
	EULERIAN,
  LAGRANGIAN,
	ANALYTIC
};
enum disp_dir   {  //dispersivity direction
	LONGITUDINAL,
  TRANSVERSE,
	VERTICAL
};
enum sourcetype {  //source zone boundary condition type
  SPECIFIED_CONC, 
	SPECIFIED_FLUX,
	SPECIFIED_CONC_RECHARGE,
	INITIAL_CONCENTRATION
};
enum trackdir   {  //tracking direction
	TRACK_FORWARD,
	TRACK_BACKWARD
};
enum orientation{  //3D tensor orientation
	OR_XX,OR_YY,OR_ZZ,
	OR_XY,OR_YX,OR_XZ,
	OR_ZX,OR_YZ,OR_ZY
};

#endif