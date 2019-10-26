#ifndef BLUEBIRDMAIN_H
#define BLUEBIRDMAIN_H

//Header file for all-knowing drivers in Parse.cpp, CardinalParse.cpp, and MainDriver.cpp

//---GENERAL---------------------
#include "MasterInclude.h"
#include "BluebirdLibrary.h"
#include "AbstractMesh.h"
#include "Grid.h"
#include "TriagonalMesh.h"
#include "PropertyZone.h"
#include "AnalysisLocation.h"
#include "RectGrid.h"
#include "FluxGrid.h"
#include "TimeSeries.h"
//---CONTAINERS------------------
#include "AbstractLayer.h"
#include "AbstractMultiLayer.h"
#include "Aquifer.h"
#include "Layer.h"
#include "MultiLayer.h"
#include "Superblock.h"
#include "Aquitard.h"
//---ELEMENTS-GEOM---------------
#include "AnalyticElem.h"
#include "FarField.h"
#include "PointElems.h"
#include "StringElem.h"
#include "Circle.h"
#include "EllipseElem.h"
//---ELEMENTS-SPECIFIC-----------
#include "Inhom.h"
#include "River.h"
#include "Leaky.h"
#include "Drain.h"
#include "AreaSink.h"
#include "Stage.h"
#include "Qnorm.h"
#include "BaseInhom.h"
#include "HorWell.h"
#include "QSpecified.h"
#include "AreaVortex.h"
#include "SmoothInhom.h"
#include "SurfaceDrain.h"
//---SURFACE WATER---------------
#include "Flownode.h"
//---CARDINAL--------------------
#include "ChemDomain.h"
#include "Particle.h"
#include "Streamline.h"
#include "TransportScheme.h"
#include "FDEulerian.h"
#include "FEEulerian.h"
#include "CharacteristicMethods.h"
#include "RxNLib.h"
#include "RxNLib/ReactionScheme.h"
#include "RxNLib/SpeciesAndSoils.h"
#include "SourceZone.h"
#include "RxNLib/CationExchange.h"
//---MESHES---------------------
#include "AbstractMesh.h"
#include "RectGrid.h"
#include "TriagonalMesh.h"
/*----------------------------------------------------------*/
struct EngineCommands{
	bool   solve;              //true if solving is on
	bool   explicitsolve;      //true if explicit solve is on
	bool   grid;               //true if gridding is on
	bool   track;              //true if tracking is on
	bool   transport;          //true if transport is on

	bool   writesol;           //true if solution should be written
	bool   writeout;					 //true if output should be written
	bool   writesolvetime;     //true if solvetime should be written
	bool   warm;               //true if solve is false

	bool   transient;          //true if transience is on

	bool   debug;              //true if debug mode on 

	bool   interpolate;        //true if interpolation is on

	bool   meshgenerate;       //true if generating mesh
	bool   DXFMeshOn;          //true if mesh DXF files are to be created
	bool   BNAMeshOn;          //true if mesh BNA files are to be created
	
	bool   obsfileexists;      //true if observations.dat exists

	int    tag;                //additional debug info

	double warmtime;           //time of warm start
	double outtimes[MAXTIMES]; //Transient output times
	int    nTimes;						 //Number of transience output times
	double timestep;           //timestep of transient calculations

	bool   socket;             //true if a socket connection

  GridCommands Grid;         //contains important gridding commands
};
/*----------------------------------------------------------*/
bool	 Parse              (char            * filename,
											     CAquifer        *&aq,  
											     CPathline      ** particles, 
											     int              &numpart, 
											     EngineCommands   &engineinfo,
													 CTriMesh        *&pMesh);

bool   CardinalParse      (char *filename,
													 CChemDomain2D   *&TransDomain,
												   CAquifer        *&pAq,
													 CTriMesh        *pMesh);
bool ObservationsParse    (char             *filename,
									         CAquifer        *&pAq);
bool OldHeadFileParse     (CAquifer        *&pAq);
//void   SocketConnect      (CAquifer *pAquifer); //empty for now (sockets not used)

void   ClearOutputFiles   ();

void   ClearTransportFiles();

#endif 
