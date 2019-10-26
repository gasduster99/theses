
#ifndef STREAMLINE_H
#define STREAMLINE_H

#include "AbstractLayer.h"
#include "AbstractDomain.h"
#include "AnalysisLocation.h"
#include "Particle.h"

struct strcell{
	double TOFa;     //starting time of flight
	double TOFb;     //ending time of flight
	double delTOF;   //duration of TOF in cell
	double v;        //average velocity in cell
	pt3D   pt;       //location of midpoint of cell     
	int    index;    //index in 1D grid
	int    ref_ind;  //reference index (W.R.T another grid)
	double q;        //streamline flux
	double *C;       //Concentration
};
enum str_update_type{UPDATE_FROM_UNEVEN,UPDATE_DIRECT,NO_UPDATE};

/****************************************************
 *  Class CStreamline
 *  Streamline Data Abstraction
 ***************************************************/
class CStreamline:public CPathline{
 private:/*----------------------------------------------------------*/

	static CTransect   **StartingFaces; 
	static int           nstartingfaces;

	const  C2DDomainABC *pDomain;       //TMP DEBUG- not necc. 2D
	const  CMesh        *pConcGrid;

	strcell             *UnevenCells;   //array of unevenly spaced 1D grid cells
	strcell							*EvenCells;     //array of evenly spaced 1D grid cells 
	int									 nEvenCells;    //number of 2D grid cells entered on traversal
	int									 nUnevenCells;  //number of 1D grid cells in evenly spaced grid

	double							 qstart;				//starting flux assoc w/ streamline
	int									 nspecies;      //number of species

	double							 maxtimestep;   //maximum time step to meet Courant constraint

	bool                 captured_out_of_grid;

	//Private Member functions-------------------------------------------
	bool HasBeenCaptured ();
  void Capture(const double &endtime);

	void TranslateEvenToUneven();
	void TranslateUnevenToEven();

	void UpdateImplicit   (const double &tstep,const double &t);
	void UpdateTVD        (const double tstep);

 public:/*-----------------------------------------------------------*/
	CStreamline();
  CStreamline(const int            ID, 
		          const CAquiferABC   *pAq, 
							const C2DDomainABC  *pDom,
							pt3D                 startpt, 
							double               flux);
 ~CStreamline();

	static void Destroy();
	static void CreateStreamlines(const CSingleLayerABC  *pLayer, 
																const C2DDomainABC *pDom, 
																CStreamline       **&lines, 
																int                 &nstreamlines, 
																double              t);
  static void AddStartingFace  (CTransect          *newface);

	double GetMaxTimeStep   () const ;
  void   Get2Dinfo				(const int i2D,int &ntimesincell,int* const pieces, Writeable1DArray TOFi) const;
	void   GetConcentrations(const int i, Writeable1DArray tmpconc) const;

	void	 SetFlux					(const double delQ, const int piece);
	void	 CleanFluxes			();
	void	 Create1DGrid     ();

	void   UpdateConcentrations(const double &tstep,const double &t,const str_update_type upty);
	void   WriteOutput();
};
#endif