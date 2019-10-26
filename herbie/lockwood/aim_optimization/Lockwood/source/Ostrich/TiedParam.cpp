/******************************************************************************
File      : TiedParam.cpp
Author    : L. Shawn Matott
Copyright : 2004, L. Shawn Matott

Definition of various 'tied' parameters. Tied parameters are variables in the 
model which are computed from the values of one or more model parameters. The 
ABC for the tied parameters encapsulates the interface used by other Ostrich 
modules, allowing various specific tied parameter relationships (linear, 
exponential, etc.) to be implemented as needed with minimal code change (just 
need to add the specific tied parameter class and some additional input file
parsing).

These specific tied-parameter classes are supported:
TiedParamLin1         : linear function of one parameter
TiedParamLin2         : linear function of two parameters
TiedParamExp          : exponential function of one parameter
TiedParamLog          : logarithmic function of one parameter
TiedDistXY            : distance between two (x,y) parameters
TiedSimpleParamRatio  : ratio of two parameters (ax + b) / (cy + d)
TiedComplexParamRatio : complex ratio of three parameters
      (a*x*y*z + b*x*y + c*x*z + d*y*z + e*x + f*y + g*z + h) / 
      (i*x*y*z + j*x*y + k*x*z + l*y*z + m*x + n*y + o*z + p)
TiedParamConstant     : parameter is assigned a constant value

Version History
07-07-04    lsm   Created
09-13-04    lsm   added distance (TiedDistXY)
01-01-07    lsm   added ratios (TiedParamSimpleRatio and TiedParamComplexRatio)
08-06-09    lsm   added constants
******************************************************************************/
#include "TiedParamABC.h"
#include "Exception.h"
#include "Utility.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/******************************************************************************
TiedParamLin1::Destroy()
******************************************************************************/
void TiedParamLin1::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamLin1)
******************************************************************************/
TiedParamLin1::TiedParamLin1(void)
{
   m_pName = NULL;
   m_pTie  = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamLin1)
******************************************************************************/
TiedParamLin1::TiedParamLin1(IroncladString name, ParameterABC * p1, 
                             UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pTie = p1;  

   //parse the config string to determine value for C0 and C1
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin1()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamLin1::GetEstVal(void)
{
   double y, x;

   x = m_pTie->GetTransformedVal();
   y = m_C1*x + m_C0;
   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamLin1::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamLin1::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Param = %s\n", m_pTie->GetName());
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "Value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedParamLin1::WriteToFile() */

/******************************************************************************
TiedParamLin2::Destroy()
******************************************************************************/
void TiedParamLin2::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamLin2)
******************************************************************************/
TiedParamLin2::TiedParamLin2(void)
{
   m_pName = NULL;
   m_pTie1  = NULL;
   m_pTie2  = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   m_C3 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamLin2)
******************************************************************************/
TiedParamLin2::TiedParamLin2(IroncladString name, ParameterABC * p1, 
                             ParameterABC * p2, UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pTie1 = p1;
   m_pTie2 = p2;

   //parse the config string to determin valuese for C0, C1, C2 and C3
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin2()");
   m_C3 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin2()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLin2()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamLin2::GetEstVal(void)
{
   double y, x1, x2;

   x1 = m_pTie1->GetTransformedVal();
   x2 = m_pTie2->GetTransformedVal();
   y = m_C3*x1*x2 + m_C2*x2 + m_C1*x1 + m_C0;

   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamLin2::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamLin2::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Param #1 = %s\n", m_pTie1->GetName());
      fprintf(pFile, "Tied Param #2 = %s\n", m_pTie2->GetName());
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "C2 = %lf\n", m_C2);
      fprintf(pFile, "C3 = %lf\n", m_C3);
      fprintf(pFile, "value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedParamLin2::WriteToFile() */

/******************************************************************************
TiedParamExp::Destroy()
******************************************************************************/
void TiedParamExp::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamExp)
******************************************************************************/
TiedParamExp::TiedParamExp(void)
{
   m_pName = NULL;
   m_pTie  = NULL;
   m_Base = 0.00;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamExp)
******************************************************************************/
TiedParamExp::TiedParamExp(IroncladString name, ParameterABC * p1, 
                           UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];
   
   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pTie = p1;

   //parse the config string to determine values for base, C0, C1, and C2
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamExp()");
   if(strcmp(tmpStr, "exp") == 0){ m_Base = 2.718;}
   else{ m_Base = atof(tmpStr);}
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamExp()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamExp()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamExp::GetEstVal(void)
{
   double y, x;

   x = m_pTie->GetTransformedVal();

   if(m_Base == 2.718) //use exp
   {
      y = (m_C2 * exp(m_C1*x)) + m_C0;
   }
   else //use pow
   {
      y = (m_C2 * pow(m_Base, (m_C1*x))) + m_C0;
   }   

   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamExp::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamExp::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Param = %s\n", m_pTie->GetName());
      fprintf(pFile, "Exponent Base = %lf\n", m_Base);
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "C2 = %lf\n", m_C2);
      fprintf(pFile, "value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedParamExp::WriteToFile() */

/******************************************************************************
TiedParamLog::Destroy()
******************************************************************************/
void TiedParamLog::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamLog)
******************************************************************************/
TiedParamLog::TiedParamLog(void)
{
   m_pName = NULL;
   m_pTie  = NULL;
   m_Base = 0.00;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   m_C3 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamLog)
******************************************************************************/
TiedParamLog::TiedParamLog(IroncladString name, ParameterABC * p1, 
                           UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pTie = p1;

   //parse the config string to determine values for base, C0, C1, C2 and C3
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   if(strcmp(tmpStr, "ln") == 0){ m_Base = 2.718;}
   else{ m_Base = atof(tmpStr);}
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   m_C3 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamLog()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamLog::GetEstVal(void)
{
   double y, x, N;

   x = m_pTie->GetTransformedVal();
   N = (m_C2*x) + m_C1;

   if(m_Base == 2.718) //use exp
   {
      y = (m_C3 * log(N)) + m_C0;
   }
   else if(m_Base == 10.00) //use log10
   {
      y = (m_C3 * log10(N)) + m_C0;
   }   
   else //convert from log10, using loga(N) = log10(N)/log10(a)
   {
      y = (m_C3 * (log10(N)/log10(m_Base))) + m_C0;
   }

   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamLog::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamLog::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Param = %s\n", m_pTie->GetName());
      fprintf(pFile, "Log Base = %lf\n", m_Base);
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "C2 = %lf\n", m_C2);
      fprintf(pFile, "C3 = %lf\n", m_C3);
      fprintf(pFile, "value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedParamLog::WriteToFile() */

/******************************************************************************
TiedDistXY::Destroy()
******************************************************************************/
void TiedDistXY::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedDistXY)
******************************************************************************/
TiedDistXY::TiedDistXY(void)
{
   m_pName = NULL;
   m_pX1  = NULL;
   m_pY1  = NULL;
   m_pX2  = NULL;
   m_pY2  = NULL;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedDistXY)
******************************************************************************/
TiedDistXY::TiedDistXY(IroncladString name, ParameterABC * px1, 
                       ParameterABC * py1,  ParameterABC * px2, 
                       ParameterABC * py2)
{
   int len;

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pX1 = px1;
   m_pY1 = py1;
   m_pX2 = px2;
   m_pY2 = py2;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedDistXY::GetEstVal(void)
{
   double d, x1, x2, y1, y2;

   x1 = m_pX1->GetTransformedVal();
   x2 = m_pX2->GetTransformedVal();
   y1 = m_pY1->GetTransformedVal();
   y2 = m_pY2->GetTransformedVal();

   d = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));

   return d;
} /* end GetEstVal() */

/******************************************************************************
TiedDistXY::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedDistXY::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied X1 = %s\n", m_pX1->GetName());
      fprintf(pFile, "Tied Y1 = %s\n", m_pY1->GetName());
      fprintf(pFile, "Tied X2 = %s\n", m_pX2->GetName());
      fprintf(pFile, "Tied Y2 = %s\n", m_pY2->GetName());
      fprintf(pFile, "value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedDistXY::WriteToFile() */

/******************************************************************************
TiedParamSimpleRatio::Destroy()
******************************************************************************/
void TiedParamSimpleRatio::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamSimpleRatio)
******************************************************************************/
TiedParamSimpleRatio::TiedParamSimpleRatio(void)
{
   m_pName = NULL;
   m_pTie1  = NULL;
   m_pTie2  = NULL;
   m_C0 = 0.00;
   m_C1 = 0.00;
   m_C2 = 0.00;
   m_C3 = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamSimpleRatio)
******************************************************************************/
TiedParamSimpleRatio::TiedParamSimpleRatio(IroncladString name, ParameterABC * p1, 
                               ParameterABC * p2, UnmoveableString configStr)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pTie1 = p1;
   m_pTie2 = p2;

   //parse the config string to determine values for C0, C1, C2 and C3
   pTok = configStr;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamSimpleRatio()");
   m_C3 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamSimpleRatio()");
   m_C2 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   j = ValidateExtraction(j, 1, 1, "TiedParamSimpleRatio()");
   m_C1 = atof(tmpStr);
   pTok += j;
   j = ExtractString(pTok, tmpStr);
   m_C0 = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamSimpleRatio::GetEstVal(void)
{
   double y, x1, x2;

   x1 = m_pTie1->GetTransformedVal();
   x2 = m_pTie2->GetTransformedVal();
   y = (m_C3*x2 + m_C2) / (m_C1*x1 + m_C0);

   return y;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamSimpleRatio::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamSimpleRatio::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Param #1 = %s\n", m_pTie1->GetName());
      fprintf(pFile, "Tied Param #2 = %s\n", m_pTie2->GetName());
      fprintf(pFile, "C0 = %lf\n", m_C0);
      fprintf(pFile, "C1 = %lf\n", m_C1);
      fprintf(pFile, "C2 = %lf\n", m_C2);
      fprintf(pFile, "C3 = %lf\n", m_C3);
      fprintf(pFile, "value = %lf\n", val);
      fprintf(pFile, "function = (C3*P2 + C2)/(C1*P1 + C0)\n");
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedParamSimpleRatio::WriteToFile() */

/******************************************************************************
TiedParamComplexRatio::Destroy()
******************************************************************************/
void TiedParamComplexRatio::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamComplexRatio)
******************************************************************************/
TiedParamComplexRatio::TiedParamComplexRatio(void)
{
   m_pName = NULL;
   m_pX  = NULL;
   m_pY  = NULL;
   m_pZ  = NULL;
   for(int i = 0; i < 8; i++){ m_N[i] = m_D[i] = 0.00;}
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamComplexRatio)
******************************************************************************/
TiedParamComplexRatio::TiedParamComplexRatio
(
   IroncladString name, 
   ParameterABC * p1, 
   ParameterABC * p2, 
   ParameterABC * p3, 
   UnmoveableString configStr
)
{
   int i, j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   m_pX = p1;
   m_pY = p2;
   m_pZ = p3;

   //parse the config string to determine values for numer and denom coeffs
   pTok = configStr;
   for(i = 0; i < 8; i++)
   {
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, i, 8, "TiedParamSimpleRatio()");
      m_N[7-i] = atof(tmpStr);
      pTok += j;
   }
   for(i = 0; i < 8; i++)
   {
      j = ExtractString(pTok, tmpStr);
      j = ValidateExtraction(j, i, 8, "TiedParamSimpleRatio()");
      m_D[7-i] = atof(tmpStr);
      pTok += j;
   }

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamComplexRatio::GetEstVal(void)
{
   double x, y, z, num, den, val;

   x = m_pX->GetTransformedVal();
   y = m_pY->GetTransformedVal();
   z = m_pZ->GetTransformedVal();

   num = m_N[7]*x*y*z + m_N[6]*x*y + m_N[5]*x*z + m_N[4]*y*z + 
         m_N[3]*x     + m_N[2]*y   + m_N[1]*z   + m_N[0];

   den = m_D[7]*x*y*z + m_D[6]*x*y + m_D[5]*x*z + m_D[4]*y*z + 
         m_D[3]*x     + m_D[2]*y   + m_D[1]*z   + m_D[0];

   val = (num / den);
   return val;   
} /* end GetEstVal() */

/******************************************************************************
TiedParamComplexRatio::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamComplexRatio::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Tied Param #1 (X) = %s\n", m_pX->GetName());
      fprintf(pFile, "Tied Param #2 (Y) = %s\n", m_pY->GetName());
      fprintf(pFile, "Tied Param #3 (Z) = %s\n", m_pZ->GetName());
      fprintf(pFile, "A = %lf\n", m_N[7]);
      fprintf(pFile, "B = %lf\n", m_N[6]);
      fprintf(pFile, "C = %lf\n", m_N[5]);
      fprintf(pFile, "D = %lf\n", m_N[4]);
      fprintf(pFile, "E = %lf\n", m_N[3]);
      fprintf(pFile, "F = %lf\n", m_N[2]);
      fprintf(pFile, "G = %lf\n", m_N[1]);
      fprintf(pFile, "H = %lf\n", m_N[0]);
      fprintf(pFile, "I = %lf\n", m_D[7]);
      fprintf(pFile, "J = %lf\n", m_D[6]);
      fprintf(pFile, "K = %lf\n", m_D[5]);
      fprintf(pFile, "L = %lf\n", m_D[4]);
      fprintf(pFile, "M = %lf\n", m_D[3]);
      fprintf(pFile, "N = %lf\n", m_D[2]);
      fprintf(pFile, "O = %lf\n", m_D[1]);
      fprintf(pFile, "P = %lf\n", m_D[0]);
      fprintf(pFile, "value = %lf\n", val);
      fprintf(pFile, "function = (Axyz + Bxy + Cxz + Dyz + Ex + Fy + Gz + H) /");
      fprintf(pFile, "           (Ixyz + Jxy + Kxz + Lyz + Mx + Ny + Oz + P)\n");
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedParamComplexRatio::WriteToFile() */

/******************************************************************************
TiedParamConstant::Destroy()
******************************************************************************/
void TiedParamConstant::Destroy(void)
{
   delete [] m_pName;
   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
CTOR (TiedParamConstant)
******************************************************************************/
TiedParamConstant::TiedParamConstant(void)
{
   m_pName = NULL;
   m_val = 0.00;
   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
CTOR (TiedParamConstant)
******************************************************************************/
TiedParamConstant::TiedParamConstant(IroncladString name, UnmoveableString pVal)
{
   int j, len;
   char * pTok;
   char tmpStr[DEF_STR_SZ];

   len = (int)strlen(name) + 10;
   NEW_PRINT("char", len);
   m_pName = new char[len];
   MEM_CHECK(m_pName);

   strcpy(m_pName, name);

   //parse the config string to determine value for C0 and C1
   pTok = pVal;
   j = ExtractString(pTok, tmpStr);
   m_val = atof(tmpStr);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
GetEstVal()

Computes the value of the tied parameter.
******************************************************************************/
double TiedParamConstant::GetEstVal(void)
{
   return m_val;
} /* end GetEstVal() */

/******************************************************************************
TiedParamConstant::Write()

Writes formatted output to pFile argument.
******************************************************************************/
void TiedParamConstant::Write(FILE * pFile, int type)
{
   double val;

   val = GetEstVal();

   if(type == WRITE_SCI)
   {
      fprintf(pFile,"%E  ", val);
   }
   else if (type == WRITE_DEC)
   {
      fprintf(pFile,"%.6lf  ", val);
   }
   else if (type == WRITE_DBG)
   {
      fprintf(pFile, "Name = %s\n", m_pName);
      fprintf(pFile, "Value = %lf\n", val);
   }/* end else if() */
   else if (type == WRITE_TX_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
   else if (type == WRITE_OPT)
   {
      fprintf(pFile, "%-18s : %E\n", m_pName, val);
   }
   else // assuming (type == WRITE_BNR)
   {
      fprintf(pFile,"%-13s  ", m_pName);
   }/* end else() */
} /* end TiedParamConstant::WriteToFile() */
