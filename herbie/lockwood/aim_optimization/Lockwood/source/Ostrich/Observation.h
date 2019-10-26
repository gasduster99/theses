/******************************************************************************
File      : Observation.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

Encapsulates a single observation point.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
10-04-04    lsm   each observation can be assigned a different token
01-01-07    lsm   added copy CTOR and Reconfigure() routines to support Surrogate-
                  model approach
******************************************************************************/
#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <stdio.h>
#include "MyTypes.h"

#define OST_OBS_FILE "OstInterpolatedObs.txt"

/******************************************************************************
class Observation

 This class represents an observation data.
 Each observation has properties keyword, line, and column.
 The  program associates each observation with the value which is found on the 
 line and coloumn after the first occurence of the keyword in <fileName>.
******************************************************************************/
class Observation
{
   StringType m_Name;
   double m_MeasuredVal;
   double m_ComputedVal;
   double m_Weight;
   StringType m_FileName;
   StringType m_Keyword;
   int m_Line;
   int m_Column;
   char m_Tok;
   bool m_bAug; //if true, include observed values in augmented OstModel file

   public:
      Observation(IroncladString name, double value ,double weight, 
                  IroncladString fileName, IroncladString keyword, int line,
                  int column, char tok, bool bAug);

      Observation(Observation * pCopy);
      Observation(void);
      ~Observation(void){Destroy();}
      void Destroy(void);

      void Reconfigure(IroncladString fileName, IroncladString keyword, 
                       int line, int column, char tok, bool bAug);

      void Write(FILE * pFile, int type);
      void WriteSim(FILE * pFile, int type);
      double GetMeasuredVal(void);
      double GetComputedVal(void);
      double GetWeight(void);
      bool IsAugmented(void) { return m_bAug;}
      UnchangeableString GetKeyword(void);
      int GetLine(void);
      int GetColumn(void);     
      UnchangeableString GetFileName(void);
      UnchangeableString GetName(void);
      void SetComputedVal(double computedVal);
      double CalcResidual(void);
      char GetToken(void){return m_Tok;}
}; /* end class Observation */

#endif /* OBSERVATION_H */



