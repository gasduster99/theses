/******************************************************************************
File      : GenConstrainedOpt.h
Author    : L. Shawn Matott
Copyright : 2005, L. Shawn Matott

Defines a general constrained optimization extension to the ObjectiveFunction class.

This class supports a veriety of cost and constraint formulations, allowing users to
define fairly generic objective functions without having to write a separate driver 
program.

This class instantiates a set of constraint classes which can be 
combined with the system cost using a user-selected penalty method (additive penalty, 
multiplicative penalty, etc.). Cost and constraints are made up of response variables
which are functions of model output and/or model parameters.
   
Version History
01-10-05    lsm   created
******************************************************************************/
#ifndef GCOP_H
#define GCOP_H

#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include "ResponseVarGroup.h"

/******************************************************************************
class GeneralConstraint

General constraints are imposed directly on the value a response variable 
specified in the response variables group (ResponseVarGroup). The penalty is 
computed as the absolute value of the violation of the constraint multipled by a 
conversion factor which converts the units of the constraint to a cost unit 
(dollars). That is, the conversion factor specifies the cost per unit of 
violation.
******************************************************************************/
class GeneralConstraint: public ConstraintABC
{
   public :
      ~GeneralConstraint(void){ Destroy();}
      void Destroy(void);
      GeneralConstraint(IroncladString name, RespVarABC * pVar, double lwr, 
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
      //pointer into the response variable group
      RespVarABC * m_pLoc; 
      double m_Lwr; //lower bound of the constraint
      double m_Upr; //upper bound of the constraint
      double m_Conv; //conversion factor (cost per unit violation)
      double m_Viol; //constraint violation
}; /* end class GeneralConstraint */

/******************************************************************************
class CapacityConstraint (capacity constraint)

Capacity constraints limit the summed value of a group of input parameters.
(for example, limits may be placed on the total pumping rate to ensure that an 
existing treatment plant is not overloaded). Constraint variables are 
stored in the ParameterGroup list and are identified by the name list. The 
penalty is computed as the absolute value of the violation of the constraint 
multiplied by a conversion factor which converts the units of the capacity violation 
(e.g. Length^3/Time for pumping rate) to a cost unit (dollars). That is, the 
conversion factor specifies the  cost per unit of capacity violation.
******************************************************************************/
class CapacityConstraint: public ConstraintABC
{
   public :
      ~CapacityConstraint(void){ Destroy();}
      void Destroy(void);
      CapacityConstraint(IroncladString name, IroncladString * nameList, 
                        int numNames, ParameterGroup * pGroup, double lwr, 
                        double upr, double conv);
      double CalcPenalty(void);
      ConstraintABC * GetNext(void){ return m_pNext;}
      void AddConstraint(ConstraintABC * pNxt);
      void Write(FILE * pFile, int type);
	  double GetLowerLimit(void){ return m_Lwr;} 
	  double GetUpperLimit(void){ return m_Upr;} 
	  double GetResponseVar(void){ return 0.00;}
	  UnchangeableString GetName(void){ return m_Name;}

   private :      
      ConstraintABC * m_pNext;
      StringType m_Name;
      StringType m_TypeStr;
      ParameterABC ** m_pParams; //array of parameters in the ParameterGroup
      int m_NumVars;  //number of parameters in capacity summation
      double m_Lwr; //lower bound of the constraint
      double m_Upr; //upper bound of the constraint
	   double m_Conv; //conversion factor (cost per unit violation)
      double m_Viol; //constraint violation
}; /* end class CapacityConstraint */

/******************************************************************************
class GCOP (General Constrained Optimization Problem)
******************************************************************************/
class GCOP : public ObjectiveFunction
{
   public :
      ~GCOP(void);
      void Destroy(void);
      GCOP(ParameterGroup * pParamGroup);      
      double CalcObjFunc(void);
      void WriteSetupToFile(FILE * pFile);
      void WriteConstraints(FILE * pFile, int type);
	  ConstraintABC * GetConstraintPtr(IroncladString pName);

   private :
      void InitFromFile(void);
      void InitResponseVars(void);
      void InitConstraints(void);

      LmtPenType  m_PenType;

      RespVarABC * m_pCostFunc;
      ParameterGroup * m_pParamGroup; // design variables
      ResponseVarGroup * m_pRespGroup; // response variables
      ConstraintABC * m_pConstraints; //linked-list of constraints
}; /* end class PATO */

//C-style functions
extern "C"
{
   IroncladString GetPenMethStr(LmtPenType i);
}

#endif /* GCOP */
