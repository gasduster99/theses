
//FEEulerian.h
#ifndef FEEULERIAN_H
#define FEEULERIAN_H

#include "CardinalInclude.h"
#include "AbstractDomain.h"
#include "AbstractMesh.h"
#include "TransportScheme.h"
#include "TriagonalMesh.h"
#include "SparseMatrix.h"

class CFEMesh;

/**************************************************************
	Class C2DFEEulerian
	Data Abstraction for eulerian advection/dispersion on triangular finite element mesh 
	Has friend status with CFEMesh
	Solves system of equations:

	([A]+w(tstep)[D]){Cn+1}=([A]-(1-w)*tstep*[D]){Cn}+tstep*((1-w){Fn} +w{Fn+1})
  or, equivalently
	[D]{Cn+1}={B}

	where w=0.0 (forward in time)
	      w=0.5 (crank-nicholson)
	      w=1.0 (backward in time)
**************************************************************/
class C2DFEEulerian:public C2DTransportScheme {

 private:/*------------------------------------------------*/

	bool            LinearParams;    //true if linear variation in parameters assumed
	bool            use_traditional; //true if velocities are obtained from heads at nodes and element saturated thickness is element avg.
	bool            SUPG;            //true if upstream weighting scheme uses streamline upwind petrov galerkin
	bool            lumped;          //true if sorption matrix is lumped


  const CSingleLayerABC *pLayer;

	const CFEMesh         *pMesh;

	double          last_loc_tstep;		//previous local time step
	double          last_Diff;				//previous species diffusion coeff

	double         *sinkflux;					//sink fluxes 
	double         *sourceflux;				//source fluxes
	double         *rechargeflux;     //recharge fluxes
	double         *borderflux;				//boundary fluxes
	double         *singularities;    //singularity strengths at nodes (e.g., Q/2pi for wells)

	CSparseMatrix  *A;								//Sorption Matrix  [A]
	CSparseMatrix  *Df;								//Diffusion Matrix [Df]
	CSparseMatrix  *M;								//Global matrix    [M]=[[A]+w(tstep)[[D]+[Df]+[S]] 
	double         *B;								//Master RHS      {B}= [[A]-(1-w)*tstep*[D]]{Cn}+tstep*((1-w){Fn} +w{Fn+1})	

	double         *wallflux;					//water flux through wall [L^3/T]
	double         *walllength;				//wall length at nodes
	double         *wallthick;				//wall thickness at nodes
	double         *dirflux;					//dirichlet flux estimate
	double         *tmpC;							//Temporary storage of concentrations (single species)

	bool           *dirichlet;

	double        **g_htheta;         //stores htheta at gauss points

	double          w;                //temporal weighing term
	double          upwt;             //upstream weighing term: 1.0 for full upstream, 0.0 for Galerkin weighing
  
	trigausspoints  gausstype;        //order of gauss integration

	bool            print_diagnostics;//true if diagnostics are printed to "diagnostics.csv"

  void        GetDispersionMatrix(const int e, const double &t, double De[3][3], double Df[3][3], double &vloc, double &dloc) const;
  void         GetSorptionMatrix (const int e, const double &t, double Ae[3][3], const int s) const;

 public:/*---------------------------------------------------*/

	C2DFEEulerian(C2DDomainABC *pDom,CFEMesh *pConcGrid, const double time_weight, const double up_weight);
 ~C2DFEEulerian();

 	void          SetSubScheme(C2DTransportScheme *pSubScheme){}//Does nothing
	void        SetParameters (trigausspoints  gausslevel, bool traditional, bool lump);
	void               UseSUPG();

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
  void  PrintMeshDiagnostics();
//	double   CalculateWallMass(const wallfluxstruct Wf) const;

/*	double CalculateWallOutflux(Ironclad1DArray Cin,
		                          Ironclad1DArray Cout,
															double Qn,
															const wallfluxstruct Wf);*/

};

#endif