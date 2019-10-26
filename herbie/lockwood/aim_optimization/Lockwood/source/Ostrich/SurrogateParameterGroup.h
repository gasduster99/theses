/******************************************************************************
File      : SurrogateParameterGroup.h
Author    : L. Shawn Matott
Copyright : 2006, L. Shawn Matott

Encapsulates a group of surrogate parameters. These are parameters that are
tied to the parameters of the complex model.

Version History
04-06-06    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef SURROGATE_PARAMETER_GROUP_H
#define SURROGATE_PARAMETER_GROUP_H

#include "ParameterGroup.h"
#include "FilePair.h"
#include <stdio.h>

/******************************************************************************
class SurrogateParameterGroup

 Represents a collection of tied surrogate parameters.
******************************************************************************/
class SurrogateParameterGroup
{
   public:
     SurrogateParameterGroup(UnmoveableString pFileName, ParameterGroup * pComplex);
     ~SurrogateParameterGroup(void){ Destroy();}
      void Destroy(void);

     void SubIntoFile(FilePipe * pPipe);
     void Write(FILE * pFile, int type);     
     TiedParamABC * GetTiedParamPtr(IroncladString name);
     int GetNumTiedParams(void){ return m_NumTied;}
     void CheckTemplateFiles(FilePair * pList);

   private:      
      TiedParamABC ** m_pTied;
      int m_NumTied;
      void InitTiedParams(IroncladString pFileName, ParameterGroup * pComplex);
}; /* end class SurrogateParameterGroup */

#endif /* SURROGATE_PARAMETER_GROUP_H */

