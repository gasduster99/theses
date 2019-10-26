/******************************************************************************
File      : ObjectiveFunction.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Computes the objective function, which can either be weighted sum of squared 
errors (WSSE) or sum of the absolute weighted error (SAWE).

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added GetObjFuncStr()
08-17-04    lsm   added reporting of memory allocations
12-18-04    lsm   added more descriptive error messages to USER obj. func.
01-10-05    lsm   removed some unused member variables
******************************************************************************/
#include "ObjectiveFunction.h"
#include "Exception.h"
#include <string.h>
#include "Utility.h"

/******************************************************************************
WSSE::CTOR

Sets the observation group pointer.
******************************************************************************/
WSSE::WSSE(ObservationGroup * pObsGroup)
{
   m_pObsGroup = pObsGroup;
   strcpy(m_ObjFuncStr, "WSSE");
   IncCtorCount();
}/* end WSSE CTOR */

/******************************************************************************
WSSE::DTOR
******************************************************************************/
WSSE::~WSSE(void)
{
   Destroy();
}/* end WSSE DTOR */

/******************************************************************************
WSSE::Destroy()
******************************************************************************/
void WSSE::Destroy(void)
{
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WSSE::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double WSSE::CalcObjFunc(void)
{
   double measured;
   double computed;
   double weight;
   double error;
   double sum;
   int numObs;
   int i;
   Observation * pObs;

   sum = 0.00;
   numObs = m_pObsGroup->GetNumObs();

   for(i = 0; i < numObs; i++)
   {
      pObs = m_pObsGroup->GetObsPtr(i);
      measured  = pObs->GetMeasuredVal();
      computed  = pObs->GetComputedVal();
      weight    = pObs->GetWeight();
      error = measured - computed;
      /*----------------------------------------------------------
      The WSSE is calculated according to equation 2.8b of the 
      WinPest Manual. NOTE: weights are referenced to the standard 
      deviation and not the variance, so they must be squared 
      along with the residuals.
      ------------------------------------------------------------*/
      sum += ((error * error) * (weight * weight));
   } /* end for() */
   return sum;
} /* end WSSE::CalcObjFunc() */

/******************************************************************************
SAWE::DTOR
******************************************************************************/
SAWE::~SAWE(void)
{
   Destroy();
}/* end SAWE DTOR */

/******************************************************************************
SAWE::Destroy()
******************************************************************************/
void SAWE::Destroy(void)
{
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
WSSE::CTOR

Sets the observation group pointer.
******************************************************************************/
SAWE::SAWE(ObservationGroup * pObsGroup)
{
  m_pObsGroup = pObsGroup;
  strcpy(m_ObjFuncStr, "SAWE");
  IncCtorCount();
}/* end SAWE CTOR */

/******************************************************************************
SAWE::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double SAWE::CalcObjFunc(void)
{
   int i;
   int numObs;
   double measured;
   double computed;
   double weight;
   double error;   
   double sumOfErrors;
   Observation * pObs;

   numObs = m_pObsGroup->GetNumObs();
   sumOfErrors  = 0.00;

   for (i = 0; i < numObs; i++)
   {
      pObs = m_pObsGroup->GetObsPtr(i);
      measured  = pObs->GetMeasuredVal();
      computed  = pObs->GetComputedVal();
      weight    = pObs->GetWeight();      
      error = (measured - computed) * weight;
      if (error < 0) { error *= -1.0; }
      sumOfErrors += error;
   }/* end for() */

   return sumOfErrors;
} /* end SAWE::CalcObjFunc() */

/******************************************************************************
UserObjFunc::CTOR

Set the name of the output file of user-defined obj. function program, where
the objective function value is stored.
******************************************************************************/
UserObjFunc::UserObjFunc(IroncladString pFileName)
{
   int len;
   m_pObsGroup = NULL;
   strcpy(m_ObjFuncStr, "USER");

   len = (int)strlen(pFileName) + 1;
   NEW_PRINT("char", len);
   m_FileName = new char[len];
   MEM_CHECK(m_FileName);

   strcpy(m_FileName, pFileName);

   m_FileStr = NULL;

   IncCtorCount();
}/* end UserObjFunc CTOR */

/******************************************************************************
UserObjFunc::DTOR
******************************************************************************/
UserObjFunc::~UserObjFunc(void)
{
   Destroy();
}/* end UserObjFunc DTOR */

/******************************************************************************
UserObjFunc::Destroy()
******************************************************************************/
void UserObjFunc::Destroy(void)
{
   delete [] m_FileName;
   delete [] m_FileStr;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
UserObjFunc::CalcObjFunc()

Computes the objective function and returns the result.
******************************************************************************/
double UserObjFunc::CalcObjFunc(void)
{
   UnchangeableString curPos; //current position within the file string
   char tmp[DEF_STR_SZ], valStr[DEF_STR_SZ];
   double val;
   int i;

   //read the output of the user-defined executable
   FileToString();
   
   //verify required output
   curPos = strstr(m_FileStr, "OST_ObjFuncVal");
   if(curPos == NULL)
   {      
      LogError(ERR_FILE_IO, "Couldn't locate OST_ObjFuncVal tag-string in model output");
      ExitProgram(1);
   }/* end if() */

   //extract the objective function value, use the last occurence
   valStr[0] = (char)NULL;
   curPos = m_FileStr;
   while(curPos != NULL)
   {
      curPos = strstr(curPos, "OST_ObjFuncVal");
      if(curPos != NULL) 
      { 
         sscanf(curPos, "%s %s", tmp, valStr);
         if(valStr[0] == (char)NULL)
         {
            LogError(ERR_FILE_IO, "Couldn't locate objective function value for model output");
            ExitProgram(1);
         }
         val = atof(valStr);
         curPos++;
      }
   }/* end while() */

   //check for model errors
   curPos = strstr(m_FileStr, "OST_ModelErrCode");
   if(curPos != NULL) 
   {      
      if(strstr(curPos, "no_errors") == NULL)
      {   
         for(i = 0; curPos[i] != '\n'; i++){ tmp[i] = curPos[i];}
         tmp[i] = (char)NULL;
         LogError(ERR_MODL_EXE, tmp);
      }/* end if() */
   }/* end if() */

   return val;     
} /* end UserObjFunc::CalcObjFunc() */

/******************************************************************************
FileToString()

Reads a file into a string.
******************************************************************************/
void UserObjFunc::FileToString(void)
{
   int fileSize;
   int i;
   FILE * pFile;

   pFile = fopen(m_FileName, "r");
   if(pFile == NULL)
   {
      FileOpenFailure("UserObjFunc::FileToString", m_FileName);
   }/* end if() */

   /*
   count number of chars in file, 
   so that fileStr can be sized
   */
   fileSize = 0;
   while(feof(pFile) == 0) 
   {
      fileSize++;
      fgetc(pFile);
   }/* end while() */   
   fileSize--;

   //size fileStr
   if(m_FileStr != NULL)
   {
      delete [] m_FileStr;
   }/* end if() */
   NEW_PRINT("char", fileSize+1);
   m_FileStr = new char[fileSize+1];
   MEM_CHECK(m_FileStr);

   //fill fileStr
   rewind(pFile);
   for(i = 0; i < fileSize; i++)
   {
      m_FileStr[i] = (char)(fgetc(pFile));
   }/* end for() */
   m_FileStr[i] = 0;

   fclose(pFile);
} /* end FileToString() */

