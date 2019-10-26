/******************************************************************************
File      : PumpAndTreat.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Defines a pump-and-treat optimization extension to the ObjectiveFunction class.

This class will support the following Pump-and-Treat objectives:
   Minimze the pumping rate
   Minimize the cost of pumping
   Minimize the cost of installation and pumping
   Minimize the cost of installation, pumping and treatment

Additionally, this class instantiates a set of constraint classes which can be 
added as a penalty to the objective function using a user-defined method 
(additive penalty, multiplicative penalty, etc.). The following constraints are
supported:
   Hydraulic gradient constraint which contain the plume
   Drawdown constraints which limit pumping rates
   Particle capture constraints which ensure plume cature
   Treatment capacity constraints which limit total pumping rate
   
Version History
05-07-04    lsm   created
01-10-05    lsm   Generalized the PatoConstraintABC class and modified to interface 
                  with abstract response variables (RespVarABC)
02-25-05    lsm   Added support for Mayer cost formulation
******************************************************************************/
#ifndef PUMP_AND_TREAT_H
#define PUMP_AND_TREAT_H

#include "GenConstrainedOpt.h"
#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include "ResponseVarGroup.h"
#include "MyTypes.h"
#include "GeometryUtility.h"

/******************************************************************************
class HydGradConstraint (hydraulic gradient constraint)

Hydraulic gradient constraints are composed of two head values which are stored
in the response variables group (ResponseVarGroup). The difference between these
two heads is the hydraulic gradient, which must be greater than or less than some
constraint value. The penalty is computed as the absolute value of the violation 
of the constraint multipled by a conversion factor which converts the units of 
the gradient violation (Length) to a cost unit (dollars). That is, the conversion
factor specifies the cost per unit length of gradient violation.
******************************************************************************/
class HydGradConstraint: public ConstraintABC
{
   public :
      ~HydGradConstraint(void){ Destroy();}
      HydGradConstraint(IroncladString name, RespVarABC * pHead1, 
		                RespVarABC * pHead2, double lwr, double upr, double conv);
      void Destroy(void);
      double CalcPenalty(void);
      ConstraintABC * GetNext(void){ return m_pNext;}
      void AddConstraint(ConstraintABC * pNxt);
	   void Write(FILE * pFile, int type); 
  	  double GetLowerLimit(void){ return m_Lwr;} 
	  double GetUpperLimit(void){ return m_Upr;} 
	  double GetResponseVar(void){ return m_pHead1->GetCurrentVal()-m_pHead2->GetCurrentVal();}
	  UnchangeableString GetName(void){ return m_Name;}

   private :
      ConstraintABC * m_pNext;
      StringType m_Name;
      StringType m_TypeStr;
	   //pointers into the response variable group
      RespVarABC * m_pHead1; 
	   RespVarABC * m_pHead2;
      double m_Lwr; //lower bound of the constraint
      double m_Upr; //upper bound of the constraint
	   double m_Conv; //conversion factor (cost per unit violation)
      double m_Viol; //constraint violation
}; /* end class HydGradConstraint */

/******************************************************************************
class DrawdownConstraint (drawdown constraint)

Drawdown constraints are composed of the initial and current head values and are
enforced at user-specified locations as specified in the response variables group 
(ResponseVarGroup). The difference between the initial and current heads is the 
drawdown, which must be greater than or less than some constraint value. The 
penalty is computed as the absolute value of the violation of the constraint 
multipled by a conversion factor which converts the units of the drawdown violation 
(Length) to a cost unit (dollars). That is, the conversion factor specifies the 
cost per unit length of drawdown violation.
******************************************************************************/
class DrawdownConstraint: public ConstraintABC
{
   public :
      ~DrawdownConstraint(void){ Destroy();}
      void Destroy(void);
      DrawdownConstraint(IroncladString name, RespVarABC * pLoc, double lwr, 
                         double upr, double conv);
      double CalcPenalty(void);      
      ConstraintABC * GetNext(void){ return m_pNext;}
      void AddConstraint(ConstraintABC * pNxt);
      void Write(FILE * pFile, int type);
  	  double GetLowerLimit(void){ return m_Lwr;} 
	  double GetUpperLimit(void){ return m_Upr;} 
	  double GetResponseVar(void){ return m_pLoc->GetCurrentVal();}
	  UnchangeableString GetName(void){ return m_Name;}

   private :
      ConstraintABC * m_pNext;
      StringType m_Name;
      StringType m_TypeStr;
      //pointer into the response variable group where the drawdown location is stored
      RespVarABC * m_pLoc; 
      double m_Lwr; //lower bound of the constraint
      double m_Upr; //upper bound of the constraint
	   double m_Conv; //conversion factor (cost per unit violation)
      double m_Viol; //constraint violation
}; /* end class DrawdownConstraint */

/******************************************************************************
class ParticleCaptureConstraint (particle capture constraint)

Particle capture constraints require that the location of a given particle be
within a well or within the original plume extents at the end of the planning 
horizon. These are specfied as (X,Y) pairs in the response variable group along
with a polygon that defines the plume geometry. At the end of the planning period,
a point-in-polygon test is performed to determine if the particle is in violation
of the capture/containment constraint. The penalty is computed as the square of 
the distance from the particle to the nearest plume boundary multiplied by a 
conversion factor which converts units from (Length^2) to cost (dollars).
Therefore, the conversion factor is the cost per unit violation of the particle
capture constraint.
******************************************************************************/
class ParticleCaptureConstraint: public ConstraintABC
{
   public :
      ~ParticleCaptureConstraint(void){ Destroy();}
      void Destroy(void);
      ParticleCaptureConstraint(IroncladString name, RespVarABC * pX, 
                                RespVarABC * pY, Point2D * pPlume, 
                                int nv, double conv);
      double CalcPenalty(void);
      ConstraintABC * GetNext(void){ return m_pNext;}
      void AddConstraint(ConstraintABC * pNxt);
      void Write(FILE * pFile, int type);
  	  double GetLowerLimit(void){ return NEARLY_ZERO;} 
	  double GetUpperLimit(void){ return NEARLY_HUGE;} 
	  double GetResponseVar(void){ return 0.00;}
	  UnchangeableString GetName(void){ return m_Name;}

   private :
      ConstraintABC * m_pNext;
      StringType m_Name;
      StringType m_TypeStr;
      //pointers to the x- and y-coordinates in the response variable group
      RespVarABC * m_pXcoord; 
      RespVarABC * m_pYcoord;
      Point2D * m_pPlume; //vertices in the plume geometry
      int m_NumVert;  //number of plume vertices
	   double m_Conv; //conversion factor (cost per unit violation)
      double m_Viol; //constraint violation
}; /* end class ParticleCaptureConstraint */

/* define a structure to contain plume vertices */
typedef struct PLUME_2D_STRUCT
{ 
   StringType name; // a name assigned to the plume
   Point2D * poly; //polygon vertices
   int nv; //number of vertices
}Plume2D;

/* define a structure to encapsulate well information */
typedef struct WELL_STRUCT
{
   StringType name;
   //design variables
   ParameterABC * pQ;
   ParameterABC * pXloc;
   ParameterABC * pYloc;
   //head at the well (or very near the well)
   RespVarABC * pHead;
   //surface topography at the well (can be a response variable or a constant)
   RespVarABC * pTopo;
   double Topo;
   //base of aquifer at the (can be a response variable or a constant)
   RespVarABC * pBase;
   double Base;
   //cost information
   double Cdrill, Cpump, Cnrg, Ctot;
}WellStruct;

/* define a structure that maps pump sizes */
typedef struct PUMP_LKUP_TABLE_STRUCT
{ 
   double Qmin, Qmax, Lmin, Lmax, cost;
}PumpLkupTableStruct;

/******************************************************************************
class PATO (pumo-and-treat optimization)
******************************************************************************/
class PATO : public ObjectiveFunction
{
   public :
      ~PATO(void);
      void Destroy(void);
      PATO(ParameterGroup * pParamGroup);      
      double CalcObjFunc(void);
      void WriteSetupToFile(FILE * pFile);
      void WriteWells(FILE * pFile, int type);
      void WriteCost(FILE * pFile, int type);
      void WriteConstraints(FILE * pFile, int type);
	  ConstraintABC * GetConstraintPtr(IroncladString pName);

   private :
      double CalcPumpingRate(void);
      double CalcOperationCost(void);
      double CalcCapitalCost(void);
      double CalcMayerCost(void);
      double CalcTreatmentCost(void);
      double LookupPumpCost(double rate, double lift);
      void InitFromFile(void);
      void InitResponseVars(void);
      void InitConstraints(void);
      void InitWells(void);
      void InitLookupTable(void);      
      void InitPlumes(void);
	  void ResizePlume(double x, double y, Plume2D * pPlume);

      PatoObjType m_ObjType;
      LmtPenType  m_PenType;

      double m_RateThresh;   //implicitly defines whether or not a well is active

      //cost factors for TOTQ formulation
      double m_ExtRateCF;    //extraction rate cost factor
      double m_InjRateCF;    //injection rate cost factor

      //captial cost factors for Mayer's formulation
      double m_MayerPumpCF;  //pump cost factor
      double m_MayerDrillCF; //drill cost factor

      //captial cost factors for RS Means formulation
      double m_FixWellCF;    //fixed well installation cost factor
      double m_VarWellCF;    //depth-dependent well installation cost factor

      //unit conversion factors for RS Means formulation
      double m_RateUCF;      //lookup table rate unit conversion factor
      double m_LiftUCF;      //lookup table lift unit conversion factor

      //cost factors for operational costs
      double m_TimeFrame;    //remediation time frame (years)
      double m_IntRate;      //interest rate
      double m_LaborRate;    //operational cost labor rate
      double m_ExtEnergyRate; //energy cost rate for extraction
      double m_InjEnergyRate; //energy cost rate for injection
      double m_AnalyticRate;  //analysis cost ($/sample)
      double m_SampleFreq;    //sample events per year
      double m_DisposalRate;  //disposal cost rate
      double m_MaintFactor;   //maintenance factor

      //cost factors for treatment costs (both operational and capital)
      double m_TreatCapCoeff;
      double m_TreatCapExpon;
      double m_TreatOpCoeff;
      double m_TreatOpExpon;

      double m_Costs[11];

      WellStruct * m_pWells;
      int m_MaxNumWells;

      PumpLkupTableStruct * m_pTbl;
      int m_TblSize;

      ParameterGroup * m_pParamGroup; // design variables
      ResponseVarGroup * m_pRespGroup; // response variables

      ConstraintABC * m_pConstraints; //linked-list of constraints

      Plume2D * m_pPlumes; //plume vertices
      int m_NumPlumes; //number of plumes
}; /* end class PATO */

#endif /* PUMP_AND_TREAT_H */
