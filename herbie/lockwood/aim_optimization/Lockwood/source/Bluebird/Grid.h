//Grid.h
#ifndef GRID_H
#define GRID_H

#include "MasterInclude.h"
#include "BluebirdLibrary.h"
#include "AnalyticElem.h"
#include "AbstractLayer.h"
#include "Aquifer.h"

class CAquifer;
/*--------------------------------------------------------------
      Header File for gridding routines
--------------------------------------------------------------*/
struct GridCommands{
	window W;                  //gridding window
	int    res;                //gridding resolution

	bool	 head;               //true if this variable is to be gridded for each layer
	bool   pot;
  bool   stream;
  bool   leak;
	bool   QxQy;
	bool   GxGy;
	bool   GxGynum;
	bool   cond;
	bool   satthick;
	bool   thick;
	bool   base;
	bool   top;
	bool   conduct;
	bool   saltH20;
	bool   vxvy;
	bool   effvel;
	bool   effvelnum;
	bool   curl;
	bool   DxxDx;
	bool   DxyDx;
	bool   VmagDer;
	bool   vder;

	bool   elem;							 //variables for gridding a single element
	CAnalyticElem *pGrdElement;
};
//--------------------------------------------------------------
void GridDriver(const CAquifer *Aq,         //pointer to aquifer
								GridCommands    G,          //Structure of gridding commands
								double         &t,          //time of gridding
								ofstream        &PROGRESS); //output file- progress.out
//--------------------------------------------------------------
void WriteProgress(const double pctdone);
void SetPlotLayer(const CSingleLayerABC *lay);
//--------------------------------------------------------------
// functions which create grid files of arbitrary functions (complex, real, or multivalued arrays)
bool cGrid(char         *filename1, 
					 char         *filename2, 
					 const window  W, 
					 const int     resolution, 
					 cmplex       (*funct)(const cmplex &z, const double &t),
					 const double &t);
//--------------------------------------------------------------
bool rGrid(char         *filename, 
					 const window  W, 
					 const int     resolution, 
					 double       (*funct)(const cmplex &z, const double &t),
					 const double &t);
//--------------------------------------------------------------
bool vGrid(       char		**filenames, 
					 const  window	 W, 
					 const  int			 resolution, 
					 const  int      numval,
					        void (*funct)(const cmplex &z, const double &t, double *v, const int nv),
					 const  double	&t);

//--------------------------------------------------------------
// dummy functions in the form cmplex f(const cmplex &z, const double &t) to be sent to cGrid
cmplex grdHeadStream		(const cmplex &z, const double &t);
cmplex grdQxQy					(const cmplex &z, const double &t);
cmplex grdGx  					(const cmplex &z, const double &t);
cmplex grdGy  					(const cmplex &z, const double &t);
cmplex grdPotStream			(const cmplex &z, const double &t);
cmplex grdElemPotStream (const cmplex &z, const double &t);
cmplex grdvxvy					(const cmplex &z, const double &t);
cmplex grdEffVel			  (const cmplex &z, const double &t);
//--------------------------------------------------------------
// dummy functions in the form double f(const cmplex &z, const double &t) to be sent to rGrid
double grdLeak					(const cmplex &z, const double &t);
double grdSaltwater     (const cmplex &z, const double &t);
double grdPot				  	(const cmplex &z, const double &t);
double grdCond					(const cmplex &z, const double &t);
double grdBase					(const cmplex &z, const double &t);
double grdThick			    (const cmplex &z, const double &t);
double grdLayerTop			(const cmplex &z, const double &t);
double grdCurl   			  (const cmplex &z, const double &t);
double grdConduct			  (const cmplex &z, const double &t);
double grdSatThick			(const cmplex &z, const double &t);
#endif