
//FDEulerian.h
#ifndef FDEULERIAN_H
#define FDEULERIAN_H

#include "CardinalInclude.h"
#include "AbstractDomain.h"
#include "AbstractMesh.h"
#include "TransportScheme.h"
#include "SparseMatrix.h"

#include "RectGrid.h"
class CRectGrid;

/**************************************************************
	Class CFDEulerian
	Data Abstraction for implicit/explicit eulerian advection/dispersion on rectangular grid
	Has friend status with CRectGrid
**************************************************************/
class C2DFDEulerian:public C2DTransportScheme {
 private:/*------------------------------------------------*/
	const CSingleLayerABC *pLayer;

	const CRectGrid *pGrid;
	bool            *active;
	double          *sinkflux;
	double          *sourceflux;
	double          *rechargeflux;
  double          *borderoutflx;
	bool            *dirichlet;

	int              nCells;

	CSparseMatrix   *A; 
	CSparseMatrix   *AW;
	double          *B;      //Right hand side of system of eqns.

	double          *dC;
	double          *Cs;

	double           last_loc_tstep;

	double           w;      //temporal weighing term: 1.0 for implicit, 0.0 for explicit
	double           upwt;   //upstream weighing term: 1.0 for full upstream, 0.0 for central weighing

	bool             print_diagnostics;

	bool    IsActive(const int k) const;
	double  NinePointStencil(const double      D1,
		                       const double      D2,
													 const orientation dir,
													 const double      C     [3][3],
													 const bool        active[3][3],
													 const double      dx    [3],
													 const double      dy    [3]);
  void    GetLocalMatrix  (const int         k,       
		                             double      a[3][3], 
		                       const double      dx[3], 
													 const double      dy[3],
													 const cmplex      pts[5], 
													 const int         species, 
													 const CVector     x_orient,
													 const double      up_weighting, 
																 double     &vloc, 
																 double     &dloc);
  double     CellByCellMB (const int        s,
													 Ironclad2DArray  C,
													 Ironclad2DArray  Cend,
													 const double     up_weighting,
													 const double    &t,
													 const double    &tstep);
  void   InitializeBorderFlux(const double &t);

 public:/*---------------------------------------------------*/
	C2DFDEulerian(C2DDomainABC *pDom,CRectGrid *pConcGrid, const double time_weight, const double up_weight);
 ~C2DFDEulerian();

 	void          SetSubScheme(C2DTransportScheme *pSubScheme){}//Does nothing
	void           Initialize (const double         startt,
		                         const transport_type ty,
									           Ironclad1DArray      porosity, 
									           Ironclad1DArray      satthick);
  void             Transport(const double     t, 
											       const double     tstep,								   
											       Ironclad2DArray  C,  
											       Ironclad2DArray  Cend,   
											       Writeable2DArray Cnew);
  double CalculateSystemMass(const int        s,
		                         Ironclad2DArray  C,
														 Ironclad2DArray  Cs,
														 double  &SorbedMass) const; 
	double CalculateBorderFlux(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const;
	double CalculateSourceGain(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const;
	double   CalculateSinkLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const;
	double  CalculateDecayLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const;
	void           WriteOutput(const int outputstep);
  void  PrintGridDiagnostics();
	void       TranslateToMT3D() const;

};

#endif