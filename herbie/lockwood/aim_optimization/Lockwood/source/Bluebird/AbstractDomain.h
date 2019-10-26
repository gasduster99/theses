//Abstract Domain.h
#ifndef ABSTRACTDOMAIN_H
#define ABSTRACTDOMAIN_H

#include "RxNLib.h"
//#include "SpeciesAndSoils.h"
#include "CardinalInclude.h"
#include "Aquifer.h"
#include "AbstractLayer.h"
#include "AbstractMesh.h"

const int uses_isotherm=-1;
/****************************************************
 *  Class CDomainABC
 *  Chemical Domain Data Abstraction
 *  for minimal knowledge access for particles, streamlines, gridding
 ***************************************************/
class CDomainABC{  
 protected:/*--------------------------------------------------------*/

  static bool EffectiveParams;

	CMesh           *pConcGrid;
	double           VolumeRatio; //ratio of liters to flow units

 public:/*-----------------------------------------------------------*/
  //Constructors:
	CDomainABC(){pConcGrid=NULL;VolumeRatio=1000.0;} //assume cubic meters
	virtual ~CDomainABC(){}

	//Accessors
  virtual const CMesh  *GetGrid () const{if (pConcGrid==NULL){ExitGracefully("CDomainABC: no grid ",RUNTIME_ERR  );}return pConcGrid;}
	
	        double  Vratio                  () const {return VolumeRatio;}
	virtual int     GetNumSpecies           () const=0;
	
	virtual double	GetDispersivity		      (const pt3D &pt, disp_dir dir) const=0;
	virtual void   	GetDispersivities		    (const pt3D &pt, const int s, disp_info &disp) const=0;
	virtual double  GetDispersionCoeff      (const pt3D &pt, const double &t, orientation dir, const CVector &alignment) const=0;
  virtual double  GetDiffusionCoeff       (const int s) const=0;

  virtual double  CalcDispersionCoeff     (vector &v, orientation dir, const CVector &alignment) const=0;

	virtual double  GetRetardation		      (const pt3D &pt, const double &t, const int s) const=0;
	virtual double  GetBulkDryDensity       (const pt3D &pt) const=0;
  virtual double  GetSpeciesDecay         (const int s) const=0;

	virtual double  GetConcentration		    (const pt3D &pt, const double &t, const int s) const=0;
	virtual void    GetConcentrations		    (const pt3D &pt, const double &t, Writeable1DArray concs) const=0;	
	virtual double  GetSorbedConc				    (const pt3D &pt, const double &t, const int s  ) const=0;
	virtual void    GetSorbedConcs				  (const pt3D &pt, const double &t, Writeable1DArray concs) const=0;	

	virtual double  GetAmbientConc          (const int s) const=0;
//virtual CVector GetAdvectiveMassFlux    (const pt3D &pt, const double &t, const int s) const=0;
//virtual CVector GetDispersiveMassFlux   (const pt3D &pt, const double &t, const int s) const=0;
//virtual cmplex  GetVertIntFlux          (const pt3D &pt, const double &t, const int s) const=0;

	virtual double  GetDirichletConc        (const pt3D &pt, const double &t, const int s) const=0;
	virtual double  GetRechargeConcentration(const pt3D &pt, const double &t, const int s) const=0;
	virtual double  GetSpecifiedMassFlux    (const pt3D &pt, const double &t, const int s) const=0;

	virtual double  GetSourceConcentration  (const int i, const double &t, const int s) const=0;

	virtual double  TranslateToSorbed       (const pt3D &pt, const double C, const int s) const=0;

	bool            UseEffectiveParams()    {return EffectiveParams;}
};

class C2DDomainABC: public CDomainABC{  
 protected:/*--------------------------------------------------------*/

	const CAquiferABC        *pAq;       //so particles can access aquifer
	const CSingleLayerABC *pLayer;
	//C2DMesh           *pConcGrid;

 public:/*-----------------------------------------------------------*/
  //Constructors:
	C2DDomainABC():CDomainABC(){/*pLayer=NULL;*/}
	virtual ~C2DDomainABC(){}

	//Accessors
  virtual const CSingleLayerABC   *GetLayer()   const {if (pLayer==NULL){ExitGracefully("C2DDomainABC: no layer  ",RUNTIME_ERR);}return pLayer;}
  virtual const CAquiferABC       *GetAquifer() const {if (pAq==NULL)   {ExitGracefully("C2DDomainABC: no aquifer",RUNTIME_ERR);}return pAq;   }

	int     GetNumSpecies           () const=0;
	
	double	GetDispersivity		      (const pt3D &pt, disp_dir dir) const=0;
	void   	GetDispersivities		    (const pt3D &pt, const int s, disp_info &disp) const=0;	
	double  GetDispersionCoeff      (const pt3D &pt, const double &t, orientation dir, const CVector &alignment) const=0;
  double  GetDiffusionCoeff       (const int s) const=0;			

  double  CalcDispersionCoeff     (vector &v, orientation dir, const CVector &alignment) const=0;

	double  GetRetardation		      (const pt3D &pt, const double &t, const int s) const=0;
	double  GetBulkDryDensity       (const pt3D &pt) const=0;

	double  GetConcentration		    (const pt3D &pt, const double &t, const int s) const=0;
	void    GetConcentrations		    (const pt3D &pt, const double &t, Writeable1DArray concs) const=0;
	double  GetSorbedConc				    (const pt3D &pt, const double &t, const int s  ) const=0;	
	void    GetSorbedConcs				  (const pt3D &pt, const double &t, Writeable1DArray concs) const=0;	

  double  GetAmbientConc          (const int s) const=0;	
//virtual CVector GetAdvectiveMassFlux    (const pt3D &pt, const double &t, const int s) const=0;
//virtual CVector GetDispersiveMassFlux   (const pt3D &pt, const double &t, const int s) const=0;
//virtual cmplex  GetVertIntFlux          (const pt3D &pt, const double &t, const int s) const=0;

	double  GetDirichletConc        (const pt3D &pt, const double &t, const int s) const=0; 
	double  GetRechargeConcentration(const pt3D &pt, const double &t, const int s) const=0;	
	double  GetSpecifiedMassFlux    (const pt3D &pt, const double &t, const int s) const=0;	

	double  GetSourceConcentration  (const int i, const double &t, const int s) const=0;

	double  TranslateToSorbed       (const pt3D &pt, const double C, const int s) const=0;

};

#endif