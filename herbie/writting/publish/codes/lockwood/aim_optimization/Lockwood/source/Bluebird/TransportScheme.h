//TransportStyle.h
#ifndef TRANSPORTSCHEME_H
#define TRANSPORTSCHEME_H

#include "AbstractDomain.h" 
#include "AbstractMesh.h"
#include "MassParticle.h"      //needed solely for RW method
#include "Streamline.h"        //needed solely for streamline method
#include "AnalyticSolutions.h" //needed solely for the analytic method


enum transport_type{ADVECTION_ONLY,DISPERSION_ONLY,ADV_AND_DISP};
/**************************************************************
	Class C2DTransportScheme
	Data Abstraction for all vertically-averaged advective/dispersive transport schemes
**************************************************************/
class C2DTransportScheme{

 protected:/*------------------------------------------------*/
	 C2DDomainABC          *pDomain;          //transport domain
	 const CSingleLayerABC *pLayer;           //single aquifer layer

	 int										nspecies;         //number of species

	 const CMesh					 *pConcGrid;        //concentration mesh (generic)
	 int										nNodes;           //number of nodes

	 const double					 *poro;             //porosity on grid/mesh
	 const double					 *thick;            //thickness on grid/mesh

	 transport_type					process_type;     //advection, dispersion, or both

	 advtype								advectiontype;    //advection type (const space step, const time step, adaptive)
	 double									adv_timestep;     //advective time step
	 double									adv_spacestep;    //advective space step

	 double									sugg_tstep;       //dispersive time step

 public:/*---------------------------------------------------*/
	 C2DTransportScheme(C2DDomainABC *pDom){
		 pDomain      =pDom;
     pLayer       =pDom->GetLayer();
		 //pAquifer     =pDom->GetAquifer();
		 nspecies     =0;
		 pConcGrid    =NULL;
		 nNodes       =0;
		 poro         =NULL;
		 thick        =NULL;
		 process_type =ADV_AND_DISP;
		 advectiontype=ADAPTIVE_TIME_STEP;
		 adv_timestep =0.0;
		 adv_spacestep=0.0;
		 sugg_tstep   =ALMOST_INF;
	 }

	 virtual ~C2DTransportScheme(){if (globaldebug){cout<<"DESTROYING ABSTRACT 2D TRANSPORT STYLE"<<endl;}}

	 void     SetAdvectionParams(const advtype A, const double tstep, const double sstep){
		 advectiontype=A;
		 adv_timestep =tstep;
		 adv_spacestep=sstep;
	 }

	 virtual void          SetSubScheme(C2DTransportScheme *pSubScheme)=0;   
	 virtual void           Initialize (const double         startt,
		                                  const transport_type ty,
									                    Ironclad1DArray      porosity, 
														          Ironclad1DArray      satthick)=0;
	 virtual void             Transport(const double     t, 
													            const double     tstep,
													            Ironclad2DArray  C,
													            Ironclad2DArray  Cend,
													            Writeable2DArray Cnew)=0;
   virtual double CalculateSystemMass(const int        s, 
		                                  Ironclad2DArray  C,
														          Ironclad2DArray  Cs,
																			double  &SorbedMass) const=0; 
	 virtual double CalculateBorderFlux(const int        s,
														          Ironclad2DArray  C,
																			Ironclad2DArray  Cend,
														          const double    &t,
														          const double    &tstep) const=0; 
	 virtual double CalculateSourceGain(const int        s,
														          Ironclad2DArray  C,
																			Ironclad2DArray  Cend,
														          const double    &t,
														          const double    &tstep) const=0; 
	 virtual double   CalculateSinkLoss(const int        s,
														          Ironclad2DArray  C,
																			Ironclad2DArray  Cend,
														          const double    &t,
														          const double    &tstep) const=0;
	 virtual double  CalculateDecayLoss(const int        s,
														          Ironclad2DArray  C,
																			Ironclad2DArray  Cend,
														          const double    &t,
														          const double    &tstep) const=0; 

	 virtual void           WriteOutput(const int outputstep)=0;
};
//**************************************************************
class CRandomWalk:public C2DTransportScheme{
 //Doesnt neccesarily have to be 2D

 private:/*---------------------------------------------------*/
  const CAquiferABC    *pAquifer;
	 
	C2DTransportScheme   *pEulerian;      //pointer to associated eulerian scheme (for mass calculation)

  CMassParticle **pMP;                //Array of Mass particles
	int             nMassParticles;     //total number of particles
	int             nCells;             //number of Cells in Grid
	double         *CellVolume;         //array of cell volumes
	int            *CellIndices;        //array of cell indices of particles
	bool           *active;             //boolean [for each particle]
  bool            UseEPVA;            //boolean- true if EPVA formulation is used

	//int						 *numparts;           //number of particles per cell numparts[k]
	//int					  **pindex;             //indices of particles in cell pindex[k][j]

 public:/*----------------------------------------------------*/

	CRandomWalk(C2DDomainABC *pDom);
 ~CRandomWalk();

	void         SetParameters(const int numparts, const bool EPVA);

	void          SetSubScheme(C2DTransportScheme *pSubScheme){pEulerian=pSubScheme;}

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
};

//**************************************************************
class CStreamlineTrans:public C2DTransportScheme{
 //Doesnt neccesarily have to be 2D

 private:/*---------------------------------------------------*/

	C2DTransportScheme *pEulerian; //pointer to associated eulerian scheme (for dispersive transport)

	CStreamline **pStreamlines;  //array of pointers to streamlines
  int           nstreamlines;  //number of streamlines 

	int					**lineIDs;       //streamline IDs for 2D cells[nCells][MAX_LINES_IN_CELL];
	int					**linePart;      //streamline segment IDs for 2D cells [nCells][MAX_LINES_IN_CELL];
  double			**lineTOF;       //Time of flight of segment in cell [nCells][MAX_LINES_IN_CELL];
	int					 *nLinesInCell;  //number of line segments in cell [nCells]
  double			 *totalTOF;      //total time of flight of all segments in cell [nCells]
	double			 *Recharge;      //recharge to cell [nCells]

	int           nCells;        //number of Cells in 2D Grid

	str_update_type update_type; //type of updating of 1D Even grid concentrations

 public:/*----------------------------------------------------*/

	CStreamlineTrans(C2DDomainABC *pDom);
 ~CStreamlineTrans();

	void         SetParameters(const double ns, const str_update_type type);

	void          SetSubScheme(C2DTransportScheme *pSubScheme){pEulerian=pSubScheme;}

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
														 double  &SorbedMass)    const{return pEulerian->CalculateSystemMass(s,C,Cs,SorbedMass);}  
	double CalculateBorderFlux(const int        s,
														 Ironclad2DArray  C,
                             Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return pEulerian->CalculateBorderFlux(s,C,Cend,t,tstep);} 
	double CalculateSourceGain(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return pEulerian->CalculateSourceGain(s,C,Cend,t,tstep);} 
	double   CalculateSinkLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return pEulerian->CalculateSinkLoss(s,C,Cend,t,tstep);}
	double  CalculateDecayLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return pEulerian->CalculateDecayLoss(s,C,Cend,t,tstep);}
	void           WriteOutput(const int outputstep);
};
//**************************************************************
class C2DAnalytic:public C2DTransportScheme{
 
 //Superimposed analytic solutions (most appropriate in simply uniform flow)

 private:/*---------------------------------------------------*/
   C2DAnalyticSol **pSolutions;
	 int              nSolutions;

 public:/*----------------------------------------------------*/

	C2DAnalytic(C2DDomainABC *pDom);
 ~C2DAnalytic();

 	void           AddSolution(C2DAnalyticSol *pSol);

	void          SetSubScheme(C2DTransportScheme *pSubScheme){}
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
														 double  &SorbedMass)    const{SorbedMass=0.0;return 0.0;}//TMP DEBUG 
	double CalculateBorderFlux(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}
	double CalculateSourceGain(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}//TMP DEBUG 
	double   CalculateSinkLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}
	double  CalculateDecayLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}//TMP DEBUG
	void           WriteOutput(const int outputstep);
};
//**************************************************************
class C2DPatchSource:public C2DTransportScheme{

 private:/*---------------------------------------------------*/
	cmplex zc;     //center of patch
	double a;      //1/2 width of patch
	double alpha;  //diminishing factor of source
	double Co;     //initial concentration
	double vx;     //velocity
	double vy;
	double Dx;     //Dispersion coeff
	double Dy;
	double R;      //Retardation Factor

 public:/*----------------------------------------------------*/

	C2DPatchSource(C2DDomainABC *pDom);
 ~C2DPatchSource();
	void          SetSubScheme(C2DTransportScheme *pSubScheme){}
	void         SetParameters(const cmplex	&z1,
										         const cmplex	&z2,
										         const double &initconc,
										         const double &deteriorate);
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
														 double  &SorbedMass)    const{SorbedMass=0;return 0.0;}//TMP DEBUG 
	double CalculateBorderFlux(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}//TMP DEBUG 
	double CalculateSourceGain(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}//TMP DEBUG 
	double   CalculateSinkLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}
	double  CalculateDecayLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}//TMP DEBUG
	void           WriteOutput(const int outputstep);
};
//**************************************************************
class C2DPointSource:public C2DTransportScheme{

 private:/*---------------------------------------------------*/
	cmplex     zc;     //source location
	double     Mo;     //initial mass or inital mass flow rate (depending)
	sourcetype type;   //constant or initial
	double angle;      //orientation of flow
	double v;          //velocity
	double Dl;         //dispersion coeff
	double Dt;

 public:/*----------------------------------------------------*/

	C2DPointSource(C2DDomainABC *pDom);
 ~C2DPointSource();

	void          SetSubScheme(C2DTransportScheme *pSubScheme){}
	void         SetParameters(const cmplex	     zcen,
										         const double	     Mass,
										         const sourcetype  ty);
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
														 double  &SorbedMass) const{SorbedMass=0.0;return Mo;} //TMP DEBUG
	double CalculateBorderFlux(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}
	double CalculateSourceGain(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}//TMP DEBUG
	double   CalculateSinkLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}	
	double  CalculateDecayLoss(const int        s,
														 Ironclad2DArray  C,
														 Ironclad2DArray  Cend,
														 const double    &t,
														 const double    &tstep) const{return 0.0;}//TMP DEBUG
	void           WriteOutput(const int outputstep);
};
#endif