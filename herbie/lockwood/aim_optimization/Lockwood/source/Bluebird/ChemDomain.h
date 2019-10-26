//ChemDomain.h
#ifndef CHEMDOMAIN_H
#define CHEMDOMAIN_H

#include "CardinalInclude.h"
#include "AbstractDomain.h"  
#include "AbstractMesh.h"
#include "SourceZone.h"
#include "FlowSink.h"
#include "PropertyZone.h"
#include "TransportScheme.h"
#include "RxNLib.h"
//#include "SpeciesAndSoils.h"
//#include "ReactionScheme.h"
#include "Grid.h"

/****************************************************
  Class CChemDomain
  Arbitrary Dimensional Chemical Properties and Transport simulation Data Abstraction
***************************************************/
class CChemDomain: public C2DDomainABC{ //Temporary-these will be split, should inherit directly from CDomainABC
 protected://----------------------------------------------------------

	//Member Variables: 
	C2DTransportScheme	*pTransScheme;	      //pointer to associated advective(and perhaps dispersive) transport scheme
	C2DTransportScheme	*pDispScheme;		      //pointer to associated dispersive transport scheme
	CReactionScheme		 **pReactSchemes;	      //array of pointers to associated reaction schemes
	int                  nRxNs;               //Number of reactions

	double            *poro;                  //array of porosities at nodes/cell centers (i=0...nNodes)

	double           **C;											//current grid concentrations      [k][s] (k=0...nNodes)(s=0...nspecies)
	double           **Ctmp;									//intermediate grid concentrations [k][s] (k=0...nNodes)(s=0...nspecies)
	double           **Cnew;									//resultant grid concentrations    [k][s] (k=0...nNodes)(s=0...nspecies)
	double           **Csorb;									//sorbed concentrations            [k][s] (k=0...nNodes)(s=0...nspecies)
  double             Cback[MAX_SPECIES];    //Ambient concentrations  [s]
  double             Sback[MAX_SPECIES];    //Ambient sorbed concentrations  [s] //Not used much now

	double						 TimeStep;							//timestep for transport calculations
	double						 AdvTimeStep;						//timestep for advective transport 
  double						 AdvSpaceStep;					//spacestep for advective transport
	
	advtype            AdvectionType;					//type of advection for transport (constant time step, adaptive time step, etc)
	disptype           DispersionType;        //type of dispersion for transport (Eulerian, lagrangian, analytic, etc)

	double						 alpha_L;								//background longitudinal dispersivity;
	double						 alpha_TH;							//background transverse horizontal dispersivity;
	double						 alpha_TV;							//background transverse vertical dispersivity;

	CSpecies				 **pSpeciesArray;					//array[nspecies] of pointers to species in domain
	int								 nspecies;							//number of species

  CSoilType				 **pSoilArray;						//array [nSoils] of pointers to soils in domain
	int								 nSoils;								//number of soils

	CIsotherm				***IsothermTable;					//array[nSoils][nspecies] table of isotherm instances for each soil/species interaction
	bool							 sorbing[MAX_SPECIES];	//array[nspecies]:true if any isotherms exist for species s
	//bool              *noneq_sorpt;         //array[nspecies]:true if sorption is non-equilibrium
	bool               UsesIsotherm[MAX_SPECIES]; //-1 for equilibrium species with isotherm
  	
	double             StartTime;             //simulation start time 
	double            *OutputTimes;					  //times to output concentration plots
	int								 nOutputTimes;				  //number of times to plot
	int                WriteInterval;         //write output (e.g., MB) every this number of time steps

	pt3D              *ObsPoints;							//locations of observation points
	int								 nObsPoints;						//number of observation points	

	double           **NodalSourceFluxes;     //[k=0..nNodes][m=0..nNodalSources]
  CSourceSink     ***pNodalSources;         //[k=0..nNodes][m=0..nNodalSources] //holds index of source assoc with node i
	int               *nNodalSources;         //[k=0..nNodes]

	bool               UseInitialConditions;

	//Protected member functions
	double               GetConcentration(const int k,const double &t, const int s  ) const;
	void                GetConcentrations(const int k,const double &t, Writeable1DArray concs) const;
  double                  GetSorbedConc(const int k,const double &t, const int s  ) const;
//	void                   GetSorbedConcs(const int k,const double &t, Writeable1DArray concs) const;
 
	double               GetNewSorbedConc(const int k,const double &t, const int s  ) const;

  void           CopyConcentrationArray(Ironclad2DArray c1,Writeable2DArray cout);
	void         WriteConcentrationToFile(const int latesttimeplotted) const;
	void    ReadInitialConditionsFromFile();

	void                    AddSourceFlux(const int k, CSourceSink *pSource, const double &Q);

	//Protected virtual member functions (geometric in nature- filled out by child classes)
	virtual int               GetSoilType(const pt3D &pt) const=0;
 	virtual void			 PlotConcentrations(const int outputstep) const=0;
	virtual void			    PlotSorbedConcs(const int outputstep) const=0;	
	virtual void         CalculateMoments(const int s,const double &t) const=0;
	virtual void			     InitializeGrid()=0;
  virtual void        InitializeSources(const double &t)=0; 
	virtual void     AdjustForSourceTerms(const double t, const double tstep, Writeable1DArray massadded)=0;

 public://-----------------------------------------------------------
  //Constructors:
  CChemDomain();
 ~CChemDomain();
	
  //Static Manipulators
 	void    EnableEffectiveParams();

  //inherited accessor functions (inherited from DomainABC)
	int     GetNumSpecies				    () const;

  double	GetDispersivity			    (const pt3D &pt, disp_dir dir) const;
	void   	GetDispersivities		    (const pt3D &pt, const int s, disp_info &disp) const;
  double  GetDiffusionCoeff       (const int s) const;

  double  CalcDispersionCoeff     (vector &v, orientation dir, const CVector &alignment) const;	

	double  GetRetardation			    (const pt3D &pt, const double &t, const int s  ) const;	
  double  GetSpeciesDecay         (const int s) const;

	double  GetSorbedConc				    (const pt3D &pt, const double &t, const int s  ) const;	
	void    GetSorbedConcs  		    (const pt3D &pt, const double &t, Writeable1DArray concs) const;
	double  GetConcentration		    (const pt3D &pt, const double &t, const int s  ) const;
	void    GetConcentrations		    (const pt3D &pt, const double &t, Writeable1DArray concs) const;
	double  GetAmbientConc          (const int s) const;

	double  GetSourceConcentration  (const int k,    const double &t, const int s) const; 

  double  TranslateToSorbed       (const pt3D &pt, const double C, const int s) const;

	//purely virtual inherited accessor functions (inherited from DomainABC)
	double  GetDispersionCoeff	    (const pt3D &pt, const double &t, orientation dir, const CVector &alignment) const=0;
	
	double  GetDirichletConc        (const pt3D &pt, const double &t, const int s) const=0;
	double  GetRechargeConcentration(const pt3D &pt, const double &t, const int s) const=0;
	double  GetSpecifiedMassFlux    (const pt3D &pt, const double &t, const int s) const=0;

  double  GetBulkDryDensity       (const pt3D &pt) const=0;  
	
	//Manipulator functions 
	void    SetGrid						    	(CMesh *pGrid);

	void    SetInitialConditions    ();
  void    SetAmbientConcs         (Ironclad1DArray backC, Ironclad1DArray backS);

	void    SetVolumeRatio          (const double V);
	void    SetOutputTimes			    (double *conctimes, int numtimes, int interval);
	void    SetObservationPts		    (pt3D   *obspts, int nobs);

	void    SetTransportTimeStep    (const double tstep);

	void    SetTransportScheme	    (C2DTransportScheme *pTScheme);
	void    SetDispersionScheme	    (C2DTransportScheme *pDScheme);
	void    AddReactionScheme		    (CReactionScheme    *pRScheme);

	void    SetAdvectionParams	    (advtype   Atype, double tstep, double sstep);
  void    SetDispersionParams	    (const disptype  Dtype, const double aL, const double aTH, const double aTV);

	void    AddToDomain					    (CSpecies          *spec);
	void    AddToDomain					    (CSoilType         *soil);
	void    AddIsotherm					    (const int i, const int s, CIsotherm *iso);
	
	void    InitializeIsothermTable ();

	//member functions 
	//TEMPORARILY VIRTUAL 
	virtual void    Transport				(ofstream &PROGRESS){ExitGracefully("CChemDomain:Transport: temporarily virtual" ,VIRTUAL_ERROR);}
};
/****************************************************
 *  Class CChemDomain2D
 *  2-Dimensional Chemical Properties Storage Data Abstraction
 ***************************************************/
class CChemDomain2D: public CChemDomain{
//class CChemDomain2D: public C2DDomainABC: public CChemDomain{ //Eventually

 private:/*----------------------------------------------------------*/

	//Member Variables: 
	double            *satthick;              //array of saturated thickness at nodes/cell centers  (i=0...nNodes)

	CAreaSource			 **pAreaSources;					//array of pointers to 2D contaminant area sources
  int								 nsources;							//number of area sources
	//CLineSource      **pLineSources;        //array of pointers to contaminant polyline sources (not really necc)
	//int                nlinesources;        //number of contaminant polyline sources
	//CPointSource		 **pPointSources;				//array of pointers to contaminant point sources (not really necc)
	//int								 npointsources;				//number of contaminant point sources

	CLinearSourceSink**pLineSourceSinks;      //array of 2D contaminant polyline sinks/sources
	int                nlinesourcesinks;      //number of contaminant polyline sinks/sources
	CPointSourceSink **pPointSourceSinks;     //array of 2D contaminant point sinks (wells)
	int                npointsourcesinks;     //number of contaminant point sinks (wells)

	CPropZone				 **pSoilZones;						//pointer to 2D soiltype zones
  int								 nSoilZones;						//number of soiltype sources

	CPropZone         *pKdZone;               //pointer to smooth Kd zone (1)

	bool               MBEnabled;             //true if mass balance enabled
	bool               initialize_only;       //true if transport process is skipped-only initialization performed
  bool               PlotSorbed;            //true if sorbed mass should be plotted

	//Private Member Functions (inherited from CChemDomain)
  int               GetSoilType(const pt3D &pt) const;

	void			 PlotConcentrations(const int outputstep) const;
	void			    PlotSorbedConcs(const int outputstep) const;

	void         CalculateMoments(const int s,const double &t) const;

	void			       PlotKd() const;

	void			     InitializeGrid();
  void        InitializeSources(const double &t); 
	void     AdjustForSourceTerms(const double t, const double tstep, Writeable1DArray massadded);
	          //(should be used only if dispersion is lagrangian)

 public:/*-----------------------------------------------------------*/

  //Constructors:
  CChemDomain2D();
 ~CChemDomain2D();

  //inherited accessor functions (from CDomainABC via CChemDomain)
	double  GetDispersionCoeff	    (const pt3D &pt, const double &t, orientation dir, const CVector &alignment) const;
 
	double  GetDirichletConc        (const pt3D &pt, const double &t, const int s) const;
	double  GetRechargeConcentration(const pt3D &pt, const double &t, const int s) const;
	double  GetSpecifiedMassFlux    (const pt3D &pt, const double &t, const int s) const;

	double  GetBulkDryDensity       (const pt3D &pt) const;

	//overwritten for continuous Kd
	double  GetRetardation			    (const pt3D &pt, const double &t, const int s  ) const;	
  double  TranslateToSorbed       (const pt3D &pt, const double C, const int s) const;
	double  GetKd                   (const cmplex &z) const;

	//manipulator functions (fully external-invisible to all but drivers) 
	void    SetAquifer              (const CAquifer  *pAquifer); //TBR 
	void    SetLayer                (const CSingleLayerABC *pLay);

	void    AddToDomain					    (CAreaSource       *zone);
	void    AddToDomain					    (CPointSourceSink  *ptsource);
	void    AddToDomain					    (CLinearSourceSink *sink);
	void    AddToDomain					    (CPropZone         *zone);
	void    AddKdZone               (CPropZone         *kdzone);

	void    DisableMassBalance      ();
	void    InitializeOnly          (){initialize_only=true;}
	void    PlotSorbedMass          (){PlotSorbed=true;}

	//member functions 
	void    Transport					  (ofstream &PROGRESS);	
};


#endif