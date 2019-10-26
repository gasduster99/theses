/******************************************************************************
File     : MyTypes.h
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

This file contains various type defintions that are used by multiple header 
files.

Version History
02-09-04    lsm   added copyright information and initial comments.
03-05-04    lsm   added PSO
03-24-04    lsm   added PSO-LevMar hybrid
07-08-04    lsm   added HUGE log values
01-18-05    lsm   Addded support for Generalized Constrained Optimization (GCOP)
11-07-05    lsm   Added support for GRID, VSA, CSA, BGA, and status reporting
                  (used in grid computing). Also added PATO support for Mayer 
                  and RS Means cost functions.
01-01-07    lsm   Additional types for finte difference increments were defined.
                  The complete list of types are:
                     1. FD_RANGE_REL : step size is range-relative
                     2. FD_VALUE_REL : step size is relative to current value
                     3. FD_ABSOLUTE  : absolute step size
                     4. FD_OPTIMAL   : compute optimal at run-time
01-01-07    lsm   Added a population initialization type for Latin Hypercube Sampling
                  Added an n-dimensional point struct, used by RBF when in 
                  Surrogate-Model mode.
******************************************************************************/
#ifndef MY_TYPES_H
#define MY_TYPES_H
#include <float.h> //for DBL_MAX

#define MY_PI (3.1415926535897932384626433832795)                 
#define NEARLY_ZERO (1E-10)
#define NEARLY_HUGE (DBL_MAX)
#define NEARLY_HUGE_LN_EXP (log(NEARLY_HUGE))
#define NEARLY_HUGE_LOG10_EXP (log10(NEARLY_HUGE))

/* ------------------------------------------------
Typical precision of model output is 7 significant 
digits. This needs to be taken into account when 
computing numerical derivatives. Step sizes that 
are smaller than the model precision may trigger
a false appearance of insensitivity.
------------------------------------------------ */
#define MODEL_PRECISION (1E-7)

/* 
An enum is defined for each type of program that can be executed.
*/
typedef enum PROGRAM_TYPE
{
   SET_INFILE       = 0,
     GA_PROGRAM     = 1,
    BGA_PROGRAM     = 2,
     SA_PROGRAM     = 3,
    CSA_PROGRAM     = 4,
    VSA_PROGRAM     = 5,
    PSO_PROGRAM     = 6,
PSODESC_PROGRAM     = 7,
 PSO_LEV_PROGRAM    = 8,
PSODESC_LEV_PROGRAM = 9,
    LEV_PROGRAM     = 10,
   POWL_PROGRAM     = 11,
   BIS_PROGRAM      = 12,
   SMP_PROGRAM      = 13,
  STEEP_PROGRAM     = 14,
   FLRV_PROGRAM     = 15,
   STATS_PROGRAM    = 16,
   UTIL_PROGRAM     = 17,
   GRID_PROGRAM     = 18,
   EVAL_PROGRAM     = 19,
	 DDS_PROGRAM     = 20, /*JRC*/
  GMLMS_PROGRAM     = 21, 
   CCPSO_PROGRAM    = 22,
   CCPSO_LEV_PROGRAM = 23,
   CCRGA_PROGRAM    = 24,
   CCSA_PROGRAM     = 25,
   CCVSA_PROGRAM    = 26,
   SCEUA_PROGRAM    = 27,
   DDDS_PROGRAM     = 28,
   GLUE_PROGRAM     = 29,
   RJSMP_PROGRAM    = 30,
   METRO_PROGRAM    = 31,
   QUIT_PROGRAM     = 32
}ProgramType;

/* 
An enum for toggling debug print statements.
*/
typedef enum DEBUG_TYPE
{
   DBG_OFF    = 0,
   DBG_ON     = 1
}DebugType;

/* 
An enum for objective function type.
*/
typedef enum OBJ_FUNC_TYPE
{
   OBJ_FUNC_WSSE    = 0, /* minimize weighted sum of squared error */
   OBJ_FUNC_SAWE    = 1, /* minimize sum of absolute weighted error */ 
   OBJ_FUNC_USER    = 2, /* minimize user-defined objective function */
   OBJ_FUNC_PATO    = 3, /* pump-and-treat optimization */
   OBJ_FUNC_GCOP    = 4  /* general constrained optimization problem */
}ObjFuncType;

/* 
An enum for pump-and-treat objectives.
*/
typedef enum PATO_OBJ_TYPE
{
   PATO_OBJ_RATE       = 0, /* minimize pumping rate */
   PATO_OBJ_OP         = 1, /* minimize operational costs */ 
   PATO_OBJ_CAP_OP     = 2, /* minimize capital and operational costs (RS Means) */
   PATO_OBJ_MAYER      = 3, /* minimize capital and operational costs (Mayer) */
   PATO_OBJ_CAP_OP_TRE = 4  /* minimize capital, operational and treatment costs */
}PatoObjType;
#define NUM_COST_FUNCS (4)

/* 
An enum for constraint handling
APM : additive penalty method
MPM : multiplicative penalty method
EPM : exponential multiplier penalty method
*/
typedef enum LMT_PENALTY_TYPE
{
   PEN_TYPE_APM,
   PEN_TYPE_MPM,
   PEN_TYPE_EPM   
}LmtPenType;
#define NUM_PEN_METHS (3)

typedef enum FINITE_DIFF_TYPE
{
   FD_FORWARD = 0, /* forward differences */
   FD_OUT_CEN = 1, /* outside central differences */
   FD_PAR_CEN = 2, /* parabolic central differences */
   FD_FIT_CEN = 3  /* best-fit central differences */
}FiniteDiffType;

typedef enum FINITE_DIFF_INC_TYPE
{
   FD_RANGE_REL = 0, /* step size is range-relative */
   FD_VALUE_REL = 1, /* step size is relative to current value */
   FD_ABSOLUTE  = 2, /* absolute step size */
   FD_OPTIMAL   = 3  /* compute optimal at run-time */
}FiniteDiffIncType;

//Methods of population initialization
typedef enum POP_INIT_TYPE
{
  RANDOM_INIT = 0,
  QUAD_TREE_INIT = 1,
  LHS_INIT = 2
}PopInitType;

//an n-dimensional point
typedef struct POINT_ND_STRUCT{int ndim; double F; double * v;}MyPoint;

//a list of n-dimensional points
typedef struct PARAMETER_LIST_STRUCT{MyPoint p; struct PARAMETER_LIST_STRUCT * pNxt;} ParameterList;

//Typedefs for type-safe arrays, courtesy James Craig
typedef const double *       Unchangeable1DArray;
typedef       double * const Unmoveable1DArray;
typedef const double * const Ironclad1DArray;

typedef const int *       Unchangeable1DIntArray;
typedef       int * const Unmoveable1DIntArray;
typedef const int * const Ironclad1DIntArray;
 
typedef const double * const *       Unchangeable2DArray;
typedef       double *       * const Unmoveable2DArray;
typedef const double * const * const Ironclad2DArray;

//note that sscanf() does not trip up type violations
typedef char *             StringType;
typedef const char *       UnchangeableString;
typedef       char * const UnmoveableString;
typedef const char * const IroncladString;

typedef const void *       UnchangeableVoidPtr;
typedef       void * const UnmoveableVoidPtr;
typedef const void * const IroncladVoidPtr;

/* define a structure to encapsulate a grid */
typedef struct GRID_STRUCT
{
	double ** p;  //parameter sets
	double *  f;  //objective function values
   double *  dp; //parameter step sizes
   int nprm;     //number of parameters in each set
}GridStruct;

/* define a structure to encapsulate status information */
typedef struct STATUS_STRUCT
{
	float pct;   //percent done
	int maxIter;  //max iterations
   int curIter;  //current iteration
   int numRuns;
}StatusStruct;

#endif /* MY_TYPES_H */

