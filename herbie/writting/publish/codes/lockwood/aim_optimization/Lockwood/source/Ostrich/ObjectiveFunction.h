/******************************************************************************
File      : ObjectiveFunction.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Computes the objective function, which can either be weighted sum of squared 
errors (WSSE) or sum of the absolute weighted error (SAWE).

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added GetObjFuncStr()
07-08-04    lsm   added WriteSetup() to base class
08-17-04    lsm   added reporting of memory allocations
01-10-05    lsm   removed some unused member variables
******************************************************************************/
#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#include "ObservationGroup.h"
#include "MyTypes.h"

/******************************************************************************
class ConstraintABC (optimization constraints) base class
******************************************************************************/
class ConstraintABC
{
   public :
      virtual void Destroy(void)=0;
      virtual double CalcPenalty(void)=0;      
      virtual ConstraintABC * GetNext(void)=0;
      virtual void AddConstraint(ConstraintABC * pNxt)=0;
	  virtual void Write(FILE * pFile, int type)=0; 
	  virtual double GetLowerLimit(void)=0; 
	  virtual double GetUpperLimit(void)=0;
	  virtual double GetResponseVar(void)=0;
	  virtual UnchangeableString GetName(void)=0;
}; /* end class ConstraintABC */

/******************************************************************************
class ObjectiveFunction
******************************************************************************/
class ObjectiveFunction
{
   public:
      ObservationGroup * m_pObsGroup;
      virtual void Destroy(void)=0;
      virtual double CalcObjFunc(void)=0;      
      char m_ObjFuncStr[20];
      UnchangeableString GetObjFuncStr(void){ return m_ObjFuncStr;}
      virtual void WriteSetupToFile(FILE * pFile)=0;
	  virtual ConstraintABC * GetConstraintPtr(IroncladString pName)=0;
};/* end class ObjectiveFunction */

/******************************************************************************
class WSSE (weighted sum of squared errors)
******************************************************************************/
class WSSE : public ObjectiveFunction
{
   public :
      ~WSSE(void);
      void Destroy(void);
      WSSE(ObservationGroup * pObsGroup);
      double CalcObjFunc(void);
      void WriteSetupToFile(FILE * pFile) { return;}
	  ConstraintABC * GetConstraintPtr(IroncladString pName){ return NULL;}
}; /* end class WSSE */

/******************************************************************************
class SAWE (sum of average weighted errors)
******************************************************************************/
class SAWE : public ObjectiveFunction
{
   public :
      ~SAWE(void);
      void Destroy(void);
      SAWE(ObservationGroup * pObsGroup);
      double CalcObjFunc(void);
      void WriteSetupToFile(FILE * pFile) { return;}
	  ConstraintABC * GetConstraintPtr(IroncladString pName){ return NULL;}
}; /* end class SAWE */

/******************************************************************************
class UserObjFunc (user-defined objective function)

The user-defined objective function can be used for more general optimization
problems than the SAWE and WSSE objective functions. This class will be used
to encapsulate a user-supplied executable which is called upon to calculate
the objective function. 

In this scenario, the model executable input parameter is not the actual model, 
but is the user-supplied objective function program. Therefore, the user must 
write driver code that executes the model program and uses the model output 
to compute the objective function. Finally, output of the objective function 
program must be sent to stdout and should be in the following format:

   1) a line containing the syntax "ObjFuncVal <value>" must be present
      
   2) a line containing the syntax "ModelErrCode <err_code>" may be used to 
   signal the Ostrich code that some error occured in the model. If no errors 
   occured, this line may be omitted, or may be set to "ModelErrCode no_errors"

******************************************************************************/
class UserObjFunc : public ObjectiveFunction
{
   public :
      ~UserObjFunc(void);
      void Destroy(void);
      UserObjFunc(IroncladString pFileName);
      double CalcObjFunc(void);
      void WriteSetupToFile(FILE * pFile) { return;}
	  ConstraintABC * GetConstraintPtr(IroncladString pName){ return NULL;}
   private:
      void FileToString(void);

      StringType m_FileName;
      StringType m_FileStr;      
}; /* end class UserObjFunc */

#endif /* OBJECTIVE_FUNCTION_H */






