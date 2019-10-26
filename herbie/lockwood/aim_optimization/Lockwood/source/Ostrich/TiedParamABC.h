/******************************************************************************
File      : TiedParamABC.h
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Encapsulates a 'tied' parameter. Tied parameters are variables in the model 
which are computed from the values of one or more model parameters. The ABC
for the tied parameters encapsulates the interface used by other Ostrich 
modules, allowing various specific tied parameter relationships (linear, 
exponential, etc.) to be implemented as needed with minimal code change (just 
need to add the specific tied parameter class and some additional input file
parsing).

These specific tied-parameter classes are supported:
TiedParamLin1   : linear function of one parameter
TiedParamLin2   : linear function of two parameters
TiedParamExp    : base exponential function of one parameter
TiedParamLog    : logarithmic function of one parameter
TiedDistXY      : distance between two (x,y) parameters
TiedParamSimpleRatio  : simple ratio of two parameters (ax + b) / (cy + d)

TiedParamComplexRatio  : complex ratio of three parameters 
(Axyz + Bxy + Cxz + Dyz + Ex + Fy + Gz + H) / (Ixyz + Jxy + Kxz + Lyz + Mx + Ny + Oz + P)

Version History
07-07-04    lsm   Created
09-13-04    lsm   added distance (TiedDistXY)
01-01-07    lsm   added ratios (TiedParamSimpleRatio and TiedParamComplexRatio)
******************************************************************************/
#ifndef TIED_PARAM_ABC_H
#define TIED_PARAM_ABC_H

#include <stdio.h>
#include "MyTypes.h"
#include "ParameterABC.h"

/******************************************************************************
class TiedParamABC

Abstract base class of a tied-parameter.
******************************************************************************/
class TiedParamABC
{
   public:            
      virtual void   GetValAsStr(UnmoveableString valStr) = 0;
      virtual void   Write(FILE * pFile, int type) = 0;
      virtual double GetEstVal(void) = 0;      
      virtual UnchangeableString GetName(void) = 0;
      virtual void Destroy(void) = 0;
}; /* end class TiedParamABC */

/******************************************************************************
class TiedParamLin1

Represents a linear function of one parameter (F = aX+b)
******************************************************************************/
class TiedParamLin1 : public TiedParamABC
{
   public:      
      TiedParamLin1(void);
      TiedParamLin1(IroncladString name, ParameterABC * p1, 
                    UnmoveableString configStr);
      ~TiedParamLin1(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      ParameterABC * m_pTie;
      double m_C1, m_C0; //coefficients      
}; /* end class TiedParamLin */

/******************************************************************************
class TiedParamLin2

Represents a linear function of two parameters (F = aX + bY + cXY + d)
******************************************************************************/
class TiedParamLin2 : public TiedParamABC
{
   public:      
      TiedParamLin2(void);
      TiedParamLin2(IroncladString name, ParameterABC * p1, ParameterABC * p2,
                    UnmoveableString configStr);
      ~TiedParamLin2(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      ParameterABC * m_pTie1;
      ParameterABC * m_pTie2;
      double m_C3, m_C2, m_C1, m_C0; //coefficients
}; /* end class TiedParamLin2 */

/******************************************************************************
class TiedParamExp

Represents an exponential function of one parameter (F = a BASE^^(bX) + c)
******************************************************************************/
class TiedParamExp : public TiedParamABC
{
   public:      
      TiedParamExp(void);
      TiedParamExp(IroncladString name, ParameterABC * p1, 
                    UnmoveableString configStr);
      ~TiedParamExp(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      ParameterABC * m_pTie;
      double m_Base;
      double m_C2, m_C1, m_C0; //coefficients      
}; /* end class TiedParamExp */

/******************************************************************************
class TiedParamLog

Represents a logarithmic function of one parameter (F = a LOG(bX+c) + d)
******************************************************************************/
class TiedParamLog : public TiedParamABC
{
   public:      
      TiedParamLog(void);
      TiedParamLog(IroncladString name, ParameterABC * p1, 
                    UnmoveableString configStr);
      ~TiedParamLog(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      ParameterABC * m_pTie;
      double m_Base;
      double m_C3, m_C2, m_C1, m_C0; //coefficients
}; /* end class TiedParamLog */

/******************************************************************************
class TiedDistXY

Represents the distance between two (x,y) parameters.
******************************************************************************/
class TiedDistXY : public TiedParamABC
{
   public:      
      TiedDistXY(void);
      TiedDistXY(IroncladString name, ParameterABC * px1, ParameterABC * py1,
                                      ParameterABC * px2, ParameterABC * py2);
      ~TiedDistXY(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      ParameterABC * m_pX1;
      ParameterABC * m_pY1;
      ParameterABC * m_pX2;
      ParameterABC * m_pY2;
}; /* end class TiedDistXY */

/******************************************************************************
class TiedParamSimpleRatio

Represents a ratio of linear functions of two parameters :
   F = (aX + b) / (cY + d)
******************************************************************************/
class TiedParamSimpleRatio : public TiedParamABC
{
   public:      
      TiedParamSimpleRatio(void);
      TiedParamSimpleRatio(IroncladString name, ParameterABC * p1, ParameterABC * p2,
                     UnmoveableString configStr);
      ~TiedParamSimpleRatio(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      ParameterABC * m_pTie1;
      ParameterABC * m_pTie2;
      double m_C3, m_C2, m_C1, m_C0; //coefficients
}; /* end class TiedParamSimpleRatio */


/******************************************************************************
class TiedParamComplexRatio

Represents a complex ratio of linear functions of three parameters :
   (Axyz + Bxy + Cxz + Dyz + Ex + Fy + Gz + H) 
   -------------------------------------------
   (Ixyz + Jxy + Kxz + Lyz + Mx + Ny + Oz + P)
******************************************************************************/
class TiedParamComplexRatio : public TiedParamABC
{
   public:      
      TiedParamComplexRatio(void);
      TiedParamComplexRatio(IroncladString name, ParameterABC * p1, ParameterABC * p2,
                     ParameterABC * p3, UnmoveableString configStr);
      ~TiedParamComplexRatio(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      ParameterABC * m_pX;
      ParameterABC * m_pY;
      ParameterABC * m_pZ;
      double m_N[8]; //numerator coefficients
      double m_D[8]; //denominator coefficients
}; /* end class TiedParamComplexRatio */

/******************************************************************************
class TiedParamConstant

Represents a constant parameter.
******************************************************************************/
class TiedParamConstant : public TiedParamABC
{
   public:      
      TiedParamConstant(void);
      TiedParamConstant(IroncladString name, UnmoveableString pVal);
      ~TiedParamConstant(void){ Destroy();}
      void Destroy(void);

      void   GetValAsStr(UnmoveableString valStr){ GetPreciseValAsStr(valStr, GetEstVal());}
      void   Write(FILE * pFile, int type);
      double GetEstVal(void);
      UnchangeableString GetName(void){ return m_pName;}
            
   private:
      StringType m_pName;
      double m_val;
}; /* end class TiedParamConstant */

#endif /* TIED_PARAM_ABC_H */
