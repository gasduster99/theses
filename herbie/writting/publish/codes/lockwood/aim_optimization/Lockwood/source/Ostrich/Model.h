/******************************************************************************
File     : Model.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

The Model class encapsulates the interaction of the Ostrich optimization tools
with the externally executed groundwater modeling program. The class divides 
model components into three groups: the parameter group, the observation group 
and the objective function group. In addition to being able to execute the 
groundwater model, the Model class provides Ostrich algorithms with access to 
these groups.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added parallel model execution, GetObjFuncStr()
03-24-04    lsm   added internal model option
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added template file check, added support for PATO.
01-01-07    lsm   Added a ModelABC and created two Model sub-classes: a standard
                  Model class and the SurrogateModel class.
07-16-07    lsm   Added support for the EPA SuperMUSE cluster
******************************************************************************/
#ifndef MODEL_H
#define MODEL_H

#include "DecisionModule.h"
#include "SuperMUSE.h"
#include "ModelABC.h"
#include "FilePair.h"
#include "FileList.h"
#include "SurrogateParameterGroup.h"
#include <stdio.h>

/******************************************************************************
class Model

******************************************************************************/
class Model : public ModelABC
{
   public:
     Model(void);
	 ~Model(void);
     void Destroy(void);
          
     //retrieve member variables
     ObservationGroup *  GetObsGroupPtr(void);
     ParameterGroup   *  GetParamGroupPtr(void);     
     ObjectiveFunction * GetObjFuncPtr(void);
     double GetObjFuncVal(void) { return m_CurObjFuncVal;}
     void SetObjFuncVal(double curVal) { m_CurObjFuncVal = curVal;}
     int                 GetCounter(void);
     ObjFuncType GetObjFuncId(void) {return m_ObjFuncId;}
     UnchangeableString GetObjFuncStr(void){return m_pObjFunc->GetObjFuncStr();}
     UnchangeableString GetModelStr(void){return m_ExecCmd;}
     void PerformWarmStartCorrection(void);

     //misc. member functions     
     double Execute(void);
	 double Execute(double viol); //include parameter bounds violations in the objective function
     void   CheckGlobalSensitivity(void);
     void   Write(double objFuncVal);
     void   WriteMetrics(FILE * pFile);
     void   Bookkeep(bool bFinal);
     int GetNumDigitsOfPrecision(void) {return m_Precision;}
     bool WarmStart(double * val);
     bool CheckCache(double * val);
   private:
      ObjFuncType         m_ObjFuncId;
      ObservationGroup  * m_pObsGroup;
      ObjectiveFunction * m_pObjFunc;
      ParameterGroup    * m_pParamGroup;
      DecisionModule    * m_pDecision;

	   DatabaseABC * m_DbaseList;
      FilePair * m_FileList;
      int m_Counter;
      int m_NumCacheHits;
      int m_Precision;
      StringType  m_ExecCmd;
      FileList *  m_pFileCleanupList;
      bool m_SerialModelExec;
      bool m_InternalModel;
      bool m_bCheckGlobalSens;
      bool m_bUseSurrogates;
      bool m_bPreserveModelOutput;
      bool m_bWarmStart;
      bool m_bCaching;
      double m_CurObjFuncVal;

      void SetCmdToExecModel(IroncladString cmd);      
      void AddFilePair(FilePair * pFilePair);
      void AddDatabase(DatabaseABC * pDbase);

protected: //can be called by DecisionModule
      double StdExecute(double viol);
      friend class DecisionModule;

protected: //can be called by SuperMUSE class
      double GatherTask(char * pDir);
      friend class SuperMUSE;
}; /* end class Model */

/******************************************************************************
class SurrogateModel

******************************************************************************/
class SurrogateModel : public ModelABC
{
   public:
     SurrogateModel(UnmoveableString pFileName, ModelABC * pComplex, char * pType);
	 ~SurrogateModel(void);
     void Destroy(void);
          
     //retrieve member variables
     ObservationGroup *  GetObsGroupPtr(void);
     ParameterGroup   *  GetParamGroupPtr(void) { return NULL;}
     ObjectiveFunction * GetObjFuncPtr(void);
     double GetObjFuncVal(void) { return m_CurObjFuncVal;}
     void SetObjFuncVal(double curVal) { m_CurObjFuncVal = curVal;}
     int                 GetCounter(void);
     ObjFuncType GetObjFuncId(void) {return m_ObjFuncId;}
     UnchangeableString GetObjFuncStr(void){return m_pObjFunc->GetObjFuncStr();}
     UnchangeableString GetModelStr(void){return m_ExecCmd;}

     SurrogateParameterGroup * GetSurrogateParamGroupPtr(void);
     void Bookkeep(bool bFinal) { return;} 
     int GetNumDigitsOfPrecision(void) {return 6;}
 
     //misc. member functions     
     double Execute(void);
     void   Write(double objFuncVal);
     void   WriteMetrics(FILE * pFile);
      void PerformWarmStartCorrection(void) { return;}

   private:
      ObjFuncType         m_ObjFuncId;
      ObservationGroup  * m_pObsGroup;
      ObjectiveFunction * m_pObjFunc;
      SurrogateParameterGroup * m_pParamGroup;

      FilePair * m_FileList;
      int m_Counter;
      StringType  m_ExecCmd;
      StringType  m_pTypeStr;
      double m_CurObjFuncVal;

      void SetCmdToExecModel(IroncladString cmd);      
      void AddFilePair(FilePair * pFilePair);
}; /* end class SurrogateModel */
#endif /* MODEL_H */

