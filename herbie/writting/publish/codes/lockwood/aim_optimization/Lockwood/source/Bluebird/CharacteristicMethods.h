//CharacteristicMethods.h

#include "AbstractDomain.h" 
#include "AbstractMesh.h"
#include "TransportScheme.h"
#include "MassParticle.h"
#include "Particle.h"    //needed for particle methods
#include "RectGrid.h" //needed for Pollocks method (not very elegant)

class CRectGrid;
/**************************************************************
	Class CMMOC
	Data Abstraction for all Modified MEthod of Characteristics scheme
**************************************************************/
//**************************************************************
class CMMOC:public C2DTransportScheme{
 //Doesnt neccesarily have to be 2D

 private:/*---------------------------------------------------*/
	const CAquiferABC    *pAquifer;

	C2DTransportScheme   *pEulerian;      //pointer to associated eulerian scheme (for dispersive transport)

	CFluxGrid            *pPollockGrid;   //pointer to grid for tracking w/ pollocks method

  CMassParticle **pMP;                //Array of Mass particles [one per species]
	int             nNodes;             //number of Nodes in Grid/Mesh
	pt3D					**BacktrackLocations; //locations after backtracking
	bool						tracked;            //true if backtracked locations have been identified

 public:/*----------------------------------------------------*/

	CMMOC(C2DDomainABC *pDom);
 ~CMMOC();

	void          SetSubScheme(C2DTransportScheme *pSubScheme){pEulerian=pSubScheme;}

	void         SetParameters(const int numpercell, CRectGrid *pGrid);

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
														 double  &SorbedMass)    const{return pEulerian->CalculateSystemMass(s,C,Cs,SorbedMass);}  //Not valid with lagrangian dispersion
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
class CMOC:public C2DTransportScheme{
 //Doesnt neccesarily have to be 2D

 private:/*---------------------------------------------------*/
	const CAquiferABC    *pAquifer;

	C2DTransportScheme   *pEulerian;        //pointer to associated eulerian scheme (for dispersive transport)

  int nparthigh;													//large number of particles per cell
	int npartlow;														//small number of particles per cell 

	CMassParticle			**pMP;	              //pointer to array of Mass Particles [one per species]
	int								 *nparts;             //number of particles per cell
	int                 nCells;             //number of Cells in Grid
	bool                tracked;
	pt3D             ***trackedpoints;      //pre-tracked point locations

	distribution_type		distribution;				//manner of distributing particles (random or structured)

  void AllocateParticles(Ironclad2DArray C);

 public:/*----------------------------------------------------*/

	CMOC(C2DDomainABC *pDom);
 ~CMOC();
	
	void         SetParameters(distribution_type ty, const int nplow, const int nphigh);

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