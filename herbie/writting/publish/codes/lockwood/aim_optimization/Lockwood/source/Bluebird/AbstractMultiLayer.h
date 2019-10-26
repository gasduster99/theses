//AbstractMultiLayer.h

#ifndef MULTILAYER_ABC_H
#define MULTILAYER_ABC_H

#include "BluebirdLibrary.h"
#include "AbstractLayer.h"

/****************************************************
 *  Class CMultiLayerABC
 *  Multi-Layer Data Abstraction
 *  for minimal knowledge access from multilayer elements
 ***************************************************/
class CMultiLayerABC{  
 private:/*----------------------------------------------------------*/

 public:/*-----------------------------------------------------------*/
  //Constructors:
	CMultiLayerABC(){}
	virtual ~CMultiLayerABC(){if (globaldebug){cout <<" DESTROYING ABSTRACT SINGLELAYER "<<endl;}}

  //Accessor functions
	virtual int								GetNumLayers         () const=0;

	virtual double						GetBase							 (const cmplex &z) const=0;
  virtual double						GetTotalThickness 	 (const cmplex &z) const=0;

  virtual Ironclad1DArray		GetCondVect					 (const cmplex &z) const=0;
  virtual Ironclad1DArray		GetBottomVect				 (const cmplex &z) const=0;
  virtual Ironclad1DArray		GetThicknessVect  	 (const cmplex &z) const=0;
	virtual Ironclad1DArray		GetPoroVect					 (const cmplex &z) const=0;
	virtual Ironclad1DArray		GetConductVect       (const cmplex &z) const=0;
	virtual Ironclad1DArray		GetNormTransmissVect (const cmplex &z) const=0;
  virtual Ironclad1DArray		GetSystemEigenvector (const cmplex &z) const=0;

	virtual anisotropy			  GetAnisotropy        (const cmplex &z, const int lev) const=0;
	
	//Special Accessors (only really useful to elements)
	virtual double						GetDeltaPot       	 ()      const=0;														
  virtual cmplex						GetBlackHole      	 ()      const=0;
	virtual int								GetCurrentIter       ()      const=0;

  //Manipulator functions (extents are the only thing an element may change)
  virtual void							UpdateExtents        (const cmplex z)=0;												
  virtual void							UpdateExtents        (const cmplex z, const double r)=0;
  
  //member functions
  virtual Ironclad1DArray		GetHeadVect        (const cmplex &z, const double &t) const=0;
	virtual Ironclad1DArray   GetPotentialVect   (const cmplex &z, const double &t) const=0;
  virtual Ironclad1DArray_z GetQxQyVect        (const cmplex &z, const double &t) const=0;
	virtual Ironclad1DArray_z GetGxVect          (const cmplex &z, const double &t) const=0;
  virtual Ironclad1DArray   GetCurlVect        (const cmplex &z, const double &t) const=0;
	virtual Ironclad1DArray   GetLeakageVect     (const cmplex &z, const double &t,
		                                            const leak_type ltype) const=0;
	virtual double            GetCompPotential   (const cmplex &z, const double &t) const=0;

	virtual double						GetNetDischarge    (const int lev, const double &t) const=0;
  virtual double						GetNetDischarge    (const double &t)								const=0;

	//I can stub out these functions and later add transport capabilities
  virtual cmplex						GetFluxThruFace      (const cmplex &z1, const cmplex &z2, const int lev,const double &t) const=0;
	virtual double						GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, 
																									const cmplex &z3, const int lev, const double &t,
																									const leak_type ltype) const=0;
	virtual void							GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, 
																									const cmplex &z3, const int lev, const double &t, 
																									double &inQ, double &outQ) const=0;
  virtual void							GetFluxDistribution  (const cmplex &z1, const cmplex &z2, const int lev, 
																									const double &t, double &Q1, 
																									double &Q2) const=0;
/*	virtual cmplex						GetEffVelocity2D     (const cmplex &z, const int lev, const double &t, 
																									const disp_info &disp) const=0;                  
	virtual vector						GetEffVelocity3D     (const pt3D &pt,  const double &t,   
																									const disp_info &disp) const=0;
  */
	//must be able to calculate these for particle tracking
  virtual cmplex						GetVelocity2D        (const cmplex &z, const int lev, const double &t) const=0;
	virtual vector						GetVelocity3D        (const pt3D &pt,  const double &t) const=0;

//	virtual double          GetBaseflow          (const cmplex &z, const double t) const=0;

};

/****************************************************
 *  Class CMLSubLayer
 *  Single Layer Data Abstraction that uses information from Multi-layer
 *  to be used for vertically-averaged transport in a single sublayer of a multi-layer model
 ***************************************************/
class CMLSubLayer:public CSingleLayerABC{
 private:/*----------------------------------------------------------*/

	 CMultiLayerABC *pML;
	 int lev;

 public:/*-----------------------------------------------------------*/
  //Constructors:
	CMLSubLayer(CMultiLayerABC *pMultilayer, int level){
		ExitGracefullyIf(pMultilayer==NULL,"CMLSubLayer:Constructor: Null MultiLayer",RUNTIME_ERR);
		ExitGracefullyIf((level<0) || (level>=pML->GetNumLayers()),"CMLSubLayer:Constructor: bad sublayer index",RUNTIME_ERR);
		pML=pMultilayer;
		lev=level;
	}
	virtual ~CMLSubLayer(){if (globaldebug){cout <<" DESTROYING ABSTRACT ML SUBLAYER "<<endl;}}

  //Accessor functions
  double					GetCond							 (const cmplex &z) const{return pML->GetCondVect(z)[lev];}
  double					GetBase							 (const cmplex &z) const{return pML->GetBottomVect(z)[lev];}
  double					GetThick 						 (const cmplex &z) const{return pML->GetThicknessVect(z)[lev];}
	double          GetPoro							 (const cmplex &z) const{return pML->GetPoroVect(z)[lev];}
	anisotropy      GetAnisotropy        (const cmplex &z) const{static anisotropy anis; anis.ratio=1.0; anis.dir=0; return anis;}
	double          GetDeltaPot  				 ()      const{return 0.0;} //should never be called externally														
  cmplex          GetBlackHole      	 ()      const{return 0.0;} //should never be called externally
	double          GetSeaLevel          ()      const{return 0.0;}
	double          GetSaltwaterSG       ()      const{return 1.0;}
	int             GetCurrentIter       ()      const{return 0;}//should never be called externally

  //Manipulator functions (extents are the only thing an element may change)
  void            UpdateExtents        (const cmplex z){}												
  void            UpdateExtents        (const cmplex z, const double r){}
  
  //member functions
  double          GetHead              (const cmplex &z, const double &t) const{return pML->GetHeadVect(z,t)[lev];}
  double          GetSaturatedThickness(const cmplex &z, const double &t) const{return pML->GetThicknessVect(z)[lev];}
	void            GetHeadAndPotential  (const cmplex &z, cmplex &omega,
		                                    double &head,    const double &t) const{omega=pML->GetPotentialVect(z,t)[lev];head=pML->GetHeadVect(z,t)[lev];}
	cmplex          GetDischargePotential(const cmplex &z, const double &t) const{return pML->GetPotentialVect(z,t)[lev];}
  cmplex          GetW                 (const cmplex &z, const double &t) const{return pML->GetQxQyVect(z,t)[lev];}
	cmplex          GetGx                (const cmplex &z, const double &t) const{return pML->GetGxVect(z,t)[lev];}
  double          GetCurl              (const cmplex &z, const double &t) const{return pML->GetCurlVect(z,t)[lev];}
	double          GetLeakage           (const cmplex &z, const double &t,
																									 const leak_type ltype) const{return pML->GetLeakageVect(z,t,ltype)[lev];}
	double          GetNetDischarge      (const double &t)									const{return pML->GetNetDischarge(lev,t);}
  cmplex          GetFluxThruFace      (const cmplex &z1, const cmplex &z2,const double &t) const{return pML->GetFluxThruFace(z1,z2,lev,t);}
	double          GetIntegratedLeakage (const cmplex &z1, const cmplex &z2, 
		                                    const cmplex &z3, const double &t,
																				const leak_type ltype) const{return pML->GetIntegratedLeakage(z1,z2,z3,lev,t,ltype);}
	void            GetIntegratedBudget  (const cmplex &z1, const cmplex &z2, 
		                                    const cmplex &z3, const double &t, 
																				double &inQ, double &outQ) const{pML->GetIntegratedBudget(z1,z2,z3,lev,t,inQ,outQ);}
  void            GetFluxDistribution  (const cmplex &z1, const cmplex &z2, 
		                                    const double &t, double &Q1, 
																				double &Q2) const{pML->GetFluxDistribution(z1,z2,lev,t,Q1,Q2);}
  cmplex          GetVelocity2D        (const cmplex &z, const double &t) const{return pML->GetVelocity2D(z,lev,t);}
	vector          GetVelocity3D        (const pt3D &pt,  const double &t) const{return pML->GetVelocity3D(pt,t);}

	cmplex          GetEffVelocity2D     (const cmplex &z, const double &t, 
																				const disp_info &disp) const{return 0.0;}//TMP DEBUG                  
	vector          GetEffVelocity3D     (const pt3D &pt,  const double &t,   
		                                    const disp_info &disp) const{return 0.0;}//TMP DEBUG  

	double          GetBaseflow          (const cmplex &z, const double t) const{return 0.0;}
};

#endif