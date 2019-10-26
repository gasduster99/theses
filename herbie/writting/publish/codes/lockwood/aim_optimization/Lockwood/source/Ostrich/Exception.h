/******************************************************************************
File     : Exception.h
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

This file defines a c-style interface for tracking and reporting any errors 
that cause the program to halt prematurely.

At the beginning of main(), InitErrors() must be called.

Prior to returning from main(), the contents of the error message and error 
code should be ouput to the user via a call to ReportErrors().

If a critical error occurs in the program, a call to LogError() should be 
made prior to exiting.

ReportErrors() will write errors to OstErrors.txt unless SetErrorFile() is 
called after InitErrors() is called.

FileOpenFailure() can be called whenever a file fails to open. This routine
will report the offending routine and file along with the current working 
directory, and then terminate the program.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-26-03    lsm   created version history field and updated comments.
11-19-03    lsm   Added WRITE_TX_BNR definition.
07-08-04    lsm   added Jacobian errror, added WRITE_OPT, doubled LINE_SIZE,
                  added GetParameterByName()
08-17-04    lsm   Now limiting the number of error messages stored in the log. 
                  User abort error was added, added GetTiedParameterByName()
10-19-05    lsm   Added support for BGA and GRID programs
01-01-07    lsm   Added five more error codes:
                     ERR_BAD_WGHT --> invalid observation weight
                     ERR_INS_PARM --> insensitive parameter
                     ERR_INS_OBS  --> insensitive observation
                     ERR_CONTINUE --> the error message is a continutation of a previous msg
                     ERR_NRINC    --> error inside of numerical recipes routine
07-13-07    lsm   Added SuperMUSE error code (ERR_SMUSE).
******************************************************************************/
#ifndef MY_EXCEPTION_H
#define MY_EXCEPTION_H

#include "ModelABC.h"
#include "AlgorithmABC.h"
#include "MyTypes.h"
#include "ParameterABC.h"
#include "TiedParamABC.h"
#include "GenConstrainedOpt.h"

/* VS7 deprecated access() and chdir() */
#ifdef WIN32
   #include <io.h>
   #define MY_ACCESS _access
   #define MY_CHDIR _chdir
#else
   #define MY_ACCESS access
   #define MY_CHDIR chdir
#endif

#define DEF_STR_SZ (320) //default string size (bytes)
#define ERR_VALUE (-999.999)

#define WRITE_SCI (1) //write values in scientific format (%E)
#define WRITE_DEC (2) //write values in decimal format (%.6lf)
#define WRITE_BNR (3) //write names/banner of values (%s)
//write names/banner of values (%s) also add information about transformations
#define WRITE_TX_BNR (4) 
#define WRITE_DBG (5) //write all variables
#define WRITE_OPT (6) //write using the 'optimal' format

/* 
An enum is defined for any error that can cause the program to halt
prematurely.
*/
typedef enum ERROR_CODE_TYPE
{
   ERR_NO_ERROR = 0,
   ERR_BAD_ARGS = 1,
   ERR_FILE_IO  = 2,
   ERR_MODL_EXE = 3,  /* execution of the model (split) */
   ERR_ARR_BNDS = 4,  /* array out of bounds */
   ERR_MISMATCH = 5,  /* parameter mismatch */
   ERR_SING_MAT = 6,  /* singular matrix */
   ERR_GRD_SIZE = 7,  /* grid size too large */
   ERR_SA_TEMP  = 8,  /* initial simulated annealing temperature */
   ERR_PRM_BNDS = 9,  /* parameter outside of bounds */
   ERR_BND_MIN  = 10, /* could not bound min */
   ERR_BND_UNK  = 11, /* unknown min bounding condition */  
   ERR_IN_PARSE = 12, /* failed to parse a line of input */
   ERR_MALLOC   = 13, /* couldn't allocate memory */
   ERR_JACOBIAN = 14, /* Jacobian insensitivity */
   ERR_ABORT    = 15, /* user abort */
   ERR_BGA      = 16, /* error in binary-coded GA */
   ERR_BAD_WGHT = 17, /* invalid observation weight */
   ERR_INS_PARM = 18, /* insensitive parameter */
   ERR_INS_OBS  = 19, /* insensitive observation */
   ERR_CONTINUE = 20, /* the error message is a continutation of a previous msg */
   ERR_NRINC    = 21, /* error inside of numerical recipes routine */
   ERR_SMUSE    = 22, /* error related to SuperMUSE parallel cluster */
   ERR_OVERFLOW = 23, /* overflow condition (typically caused by divide-by-zero) */
   ERR_NULL_PTR = 24, /* NULL pointer detected */
   ERR_STALL    = 25, /* algorithm stall detected */
   ERR_CLEANUP  = 26,
   ERR_PRM_NEST = 27  /* parameter is substring of another parameter */
}ErrorCodeType;

/* NUM_ERRORS is used to size ErrorMap defined in Exception.cpp */
#define NUM_ERRORS (28) 

/* 
These c-style routines allow access to the unique instance of the 
error structure defined in Exception.cpp.
*/
extern "C" 
{
   bool IsQuit(void);
   void InitErrors(void);
   void ReportErrors(void);   
   ErrorCodeType GetErrorCode(void);
   void SetErrorFile(IroncladString filename);
   void LogError(ErrorCodeType err, IroncladString msg);
   void FileOpenFailure(IroncladString routine, IroncladString file);
   void EndOfFileFailure(IroncladString routine, IroncladString file);
   void MissingTokenFailure(IroncladString token, IroncladString file);
   void ExitProgram(int code);

   void IncCtorCount(void);
   void IncDtorCount(void);

   void RegisterModelPtr(ModelABC * pModel);
   void RegisterAlgPtr(AlgorithmABC * pAlg);
   int GetNumDigitsOfPrecision(void);
   ParameterABC * GetParameterByName(IroncladString pName);
   TiedParamABC * GetTiedParameterByName(IroncladString pName);
   ConstraintABC * GetConstraintByName(IroncladString pName);
   IroncladString GetParameterName(int idx);
   IroncladString GetParameterValStr(int idx, double val);
   #define MEM_CHECK(a) MemCheck((void *)(a), __LINE__, __FILE__);

   //use this macro for memory debugging
   //#define NEW_PRINT(a, b) printf("new %s[%d]\n", a, b);
   //use this macro to turn off memory debugging
   #define NEW_PRINT(a, b) 

   void MemCheck(void * pMem, int line, char * file);
}

#endif /* MY_EXCEPTION_H */

