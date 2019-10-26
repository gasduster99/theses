//FluxGrid.h

#ifndef FLUXGRID_H
#define FLUXGRID_H

#include "MasterInclude.h"
#include "RectGrid.h"
#include "AbstractLayer.h"


class CFluxGrid{

 private:/*----------------------------------------------------------*/
	
	 CRectGrid *pGrid; 
	 double    **vx;
	 double    **vy;

 public:/*-----------------------------------------------------------*/
	//Constructors
	CFluxGrid();
  CFluxGrid(CRectGrid *grid);
	~CFluxGrid();

	void   Initialize   (const CSingleLayerABC *pLayer);

	pt3D   PollockTrack (const pt3D   &pt, 
		                   const double &tstep,
											 const bool    forward);

	void   WriteGeometry() const;

};

#endif 