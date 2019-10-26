/******************************************************************************
File     : ModelABC.h
Author   : L. Shawn Matott
Copyright: 2006, L. Shawn Matott

Abstract Base Class for models

Version History
04-04-06    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef MODEL_ABC_H
#define MODEL_ABC_H

#include "ObservationGroup.h"
#include "ObjectiveFunction.h"
#include "ParameterGroup.h"
#include <stdio.h>

/******************************************************************************
class ModelABC

Abstract base class for Models (main [complex] model and surrogate models)
******************************************************************************/
class ModelABC
{
   public:
      virtual void Destroy(void)=0;
      virtual ObservationGroup * GetObsGroupPtr(void) = 0;
      virtual ParameterGroup *  GetParamGroupPtr(void) = 0;
      virtual ObjectiveFunction * GetObjFuncPtr(void) = 0;
      virtual double GetObjFuncVal(void) = 0;
      virtual void SetObjFuncVal(double curVal) = 0;
      virtual int GetCounter(void) = 0;
      virtual ObjFuncType GetObjFuncId(void) = 0;
      virtual UnchangeableString GetObjFuncStr(void) = 0;
      virtual UnchangeableString GetModelStr(void) = 0;
      virtual double Execute(void) = 0;
      virtual void Write(double objFuncVal) = 0;
      virtual void WriteMetrics(FILE * pFile) = 0;
      virtual void Bookkeep(bool bFinal) = 0;
      virtual int GetNumDigitsOfPrecision(void) = 0;
      virtual void PerformWarmStartCorrection(void) = 0;
}; /* end class ModelABC */

#endif /* MODEL_ABC_H */

