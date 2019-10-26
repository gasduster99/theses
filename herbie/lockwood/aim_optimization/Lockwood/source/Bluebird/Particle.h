//Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H

#include "MasterInclude.h"
#include "BluebirdLibrary.h"
#include "AbstractAquifer.h"
#include "AbstractDomain.h"  //for mass particles-to get dispersivity
//#include "FluxGrid.h"
class CFluxGrid;

//--Particle Parameters-----------------------------
const int    MAX_PARTICLES=   10000;         //maximum allowed particles	     
const int    MAX_STEPS=       5000;       //maximum increments allowed

const double DELTA0=         10.00;				//divisor of min_time_step if q=0
const double G_MIN_TIME_STEP= 0.0001;		//global minimum time step basis for 1st step
const double G_MIN_STEP=      0.00001;	//global minimum space step within aquifer
const double MAX_STEP=        5.00;				//maximum space step within aquifer
const double BAD_TURN=        0.707106781187; //90 degree turn -cos of 45

const double LIN_TOLERANCE=   0.01;        //tolerance for linearization distance
const double DIST_FRAC=       0.02;        //minimum fraction of path distance for space step in clean up

/****************************************************
 *  Struct pathnode
 *  Particle Tracking Location Data Abstraction
 ***************************************************/
struct pathnode{
	pt3D pt;														//location of pathline point
  double t;																//time of pathline traversal
	double s;																//distance along pathline (x-coord in 1-D transport)
  pathnode *next;													//next pathline point
	pathnode *last;													//previous pathline point
};
/****************************************************
 *  Class CParticle
 *  Particle/Pathline Tracking Data Abstraction
 *  only meant for single layer, SS particle tracking
 ***************************************************/
class CParticle{
 protected:/*----------------------------------------------------------*/
	//Static Variables
	static double			 sm_change;								//percent change considered too small
	static double			 lg_change;								//percent change considered too big
	static double			 bad_change;							//percent change considered way too big	
	static double			 sm_adjust;								//increase of time step if too small
	static double			 lg_adjust;								//reduction of time step if too large
	static double			 bad_adjust;							//reduction of time step if way too large

  //Member variables
  int                partID;									//particle ID
  int                dir;											//direction 1 for forward, -1 back
	const CAquiferABC *pAq;
  
  //Virtual Member Functions
	virtual void			 Advect          (const vector &movement, const double &tstep)=0;
	virtual bool			 HasBeenCaptured ()=0;
  virtual void			 BackTrack       ()=0;  
  virtual void			 Capture         (const double &endtime)=0;
  virtual void			 CleanPath       ()=0;

 public:/*-----------------------------------------------------------*/
  //Constructors
	CParticle();
  CParticle(int							 ID,
						const CAquiferABC  *pAquifer, 
						trackdir				 direct);
  virtual ~CParticle();
 
	//Static Functions
	static void				 SetPrecision     (int Precision);
	static void				 SetTrackPrecision(const double sm_ch,const double lg_ch,const double bd_ch,
														  	       const double sm_ad,const double lg_ad,const double bd_ad);

	//Manipulator functions
  void               SetAquifer       (const CAquiferABC  *pAquifer);  
	void               SetDirection     (trackdir	direct);

  //Accessor Functions
  int								 GetID            () const;

	//Purely virtual public functions
	virtual bool       IsCaptured       () const=0;  
  //virtual pt3D		   GetLocation      (const double &T) const=0;
  virtual pt3D		   GetLocation      () const=0;
	virtual double		 GetCurrentT      () const=0;
	virtual double		 GetLastT         () const=0;

	//Member Functions
  void							 Track            (const double    &timeperiod, 
																			 const advtype    atype, 
																			       double     tstep,
																			       double     sstepconst,
																			 const bool       eff_vel, 
																			 const disp_info &disp, 
																			 const bool       track3D); 
	void               TrackPollock     (const double    &timeperiod, 
		                                   const int intervals, 
																			 CFluxGrid       *pFluxGrid, 
																			 const bool       track3D,
																			 const bool       forward);

	//void

	virtual void			 WriteOutput      () const=0;
	virtual void			 WriteOutput      (const double &t) const=0;
};

const int BUFFERSIZE=200;
/****************************************************
 *  Class CPathline
 *  Contaminant Pathline Data Abstraction
 ***************************************************/
class CPathline:public CParticle{
 protected:/*----------------------------------------------------------*/
	//Static Variables 
	static pathnode *PATHEND;			 					 //blank pathnode- denotes end of path in pathline search
	static pathnode *PATHSTART;              //blank pathnode- denotes start of path in pathline search
	static int       Resolution;						 //resolution of pathline output- does not affect tracking accuracy
	static double    Duration;							 //duration of pathline travel times

	static CFluxGrid *pPollockGrid;          //Grid for pollocks method

  pathnode        *firstnode;              //pointer to first point on pathline
	pathnode        *lastnode;               //pointer to last point on pathline
	int              numsteps;               //number of pathline steps
	bool             captured;               //true if captured

  pathnode *buffarray[BUFFERSIZE];
	int       buffercount;


	//Member Functions (inherited from particle)				                     
	void   Advect          (const vector &movement, const double &tstep);
	void   Capture         (const double &endtime);
	bool	 HasBeenCaptured ();																											  
  void	 BackTrack       ();                               


 public:/*-----------------------------------------------------------*/
	CPathline();
  CPathline(int							    ID,
		        const CAquiferABC  *pAquifer, 
						pt3D					      startpt, 
						trackdir				    direct);
 ~CPathline();

	//Static Functions
  static void   SetResolution    (int res);
	static void   SetDuration      (double time);
	static double GetDuration      ();
  static void	  TrackAllPathlines(CPathline **Particles,const int numpart,
																	const double timeperiod, ofstream &PROGRESS);
  static void	  WriteAllPathlines(CPathline **Particles,const int numpart);
	static void   SetPollockGrid   (CFluxGrid *pFluxGrid);
	static CFluxGrid *GetPollockGrid(){return pPollockGrid;}
	static void   Destroy          ();

	//Accessor Functions (inherited from particle)
	bool          IsCaptured       () const;
  pt3D  				GetLocation      () const;                              
	double				GetCurrentT      () const;   
	double				GetLastT         () const;
  const CAquiferABC *GetAquifer() const {return pAq;}

  //Accessor Functions
  pt3D  				GetLocation			 (const double &T) const; 
  int						GetNumOfSteps		 () const;
	//void        GetGeometry      (pt3D *points, double *times, int numpoints) const;  

  //Member functions
  void          ReversePath      ();
  void	        CleanPath        ();
	void					WriteOutput	     () const;
	void					WriteOutput			 (const double &t) const;
};

#endif