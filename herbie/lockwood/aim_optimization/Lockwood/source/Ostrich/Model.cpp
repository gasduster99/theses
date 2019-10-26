/******************************************************************************
File     : Model.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

The Model class encapsulates the interaction of the Ostrich optimization tools
with the externally executed modeling program. The class divides 
model components into three groups: the parameter group, the observation group 
and the objective function group. In addition to being able to execute the 
model, the Model class provides Ostrich algorithms with access to 
these groups. 

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
11-25-03    lsm   Modified to support long filenames in model executable.
03-05-04    lsm   added parallel model execution,
03-24-04    lsm   added support for internal models
07-08-04    lsm   added PATO hooks
08-17-04    lsm   RAM fragmentation fixes, metrics collection and reporting
                  Added template file check, added support for PATO.
12-02-04    lsm   Added support for Geometry parameters
01-18-05    lsm   Addded support for Generalized Constrained Optimization (GCOP)
01-18-05    lsm   Addded check for existence of model executable
10-22-05    lsm   Added support for SGI/IRIX environment, refined error handling
01-01-07    lsm   Added support for checking global sensitivity of parameters and
                  observations; items that are not sensitive over the entire range
                  of parameter values are removed from consideration. To select this
                  option, users enter the following line in the main configuration
                  section:
                     CheckSensitivities   yes
01-01-07    lsm   Added support for surrogate-model approach to calibration, ranking and
                  selection. To select this option, users enter the following line in the 
                  main configuration section:
                     SurrogateApproach   yes
01-01-07    lsm   Added model initialization and bookkepping calls.
07-13-07    lsm   Added support for the EPA SuperMUSE cluster. To select this option, 
                  users enter the following line in the main configuration section:
                     SuperMUSE yes
                  Users must also define various parameters in the 'SuperMUSE' 
                  section of the input file (See SuperMUSE.h for details):
                     BeginSuperMUSE
                       OstrichTaskerHostName <host_name>
                       TaskFile         <file_name>
                       TempFile         <file_name>
                       SuccessFile      <file_name>
                       ErrorFile        <file_name>
                       ScriptFile       <file_name>
                       MaxJobTime       <integer>
                     EndSuperMUSE
03-20-2010   lsm   Added support for warm starting Ostrich using previously
                   generated OstModel0.txt file.
******************************************************************************/
#include "mpi_stub.h"
#include "AccessConverter.h"
#include "NetCDFConverter.h"
#include "IsoParse.h"
#include "PumpAndTreat.h"
#include "Model.h"
#include "SuperMuseUtility.h"
#include "Utility.h"
#include "WriteUtility.h"
#include "Exception.h"
#include "MyPBS.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef WIN32
   #include <io.h>
   #include <direct.h>
#else
   #include <unistd.h>
#endif

#define JOB_SUCCEDDED (0)
#define JOB_FAILED    (1)
#define JOB_TIMED_OUT (2) 

/******************************************************************************
 default CTOR
******************************************************************************/
Model::Model(void)
{
   FilePair * pFilePair;
   char * line;
   char tmp1[DEF_STR_SZ];
   char tmp2[DEF_STR_SZ];
   char tmp3[DEF_STR_SZ];   
   FILE * pInFile;
   int id;
   int i;
   int j;
   bool quoteWrap; //if true, then executable must be wrapped in quotes

   RegisterModelPtr(this);

   IroncladString inFileName = GetInFileName();
   UnmoveableString pDirName = GetExeDirName();

   m_pDecision = NULL;
   m_ExecCmd = NULL;
   m_pParamGroup = NULL;
   m_pObsGroup = NULL;
   m_FileList = NULL;
   m_DbaseList = NULL;
   m_pFileCleanupList = NULL;
   m_Counter = 0;
   m_Precision = 6;
   m_pObjFunc = NULL;
   m_ObjFuncId = OBJ_FUNC_WSSE;
   m_SerialModelExec = true;
   m_bCheckGlobalSens = false;
   m_bPreserveModelOutput = false;
   m_bUseSurrogates = false;
   m_CurObjFuncVal = 0.00;
   m_bWarmStart = false;
   m_bCaching = false;
   m_NumCacheHits = 0;

   #ifdef WIN32
      m_pFileCleanupList = new FileList("Ostrich.exe");
   #else
      m_pFileCleanupList = new FileList("Ostrich");
   #endif
   
   pInFile = fopen(inFileName, "r");
   if(pInFile == NULL)
   {
      FileOpenFailure("Model::CTOR", inFileName);
   }/* end if() */

   //check for critical entries, entries which have no reasonable defaults
   FindToken(pInFile, "BeginFilePairs", inFileName);
   FindToken(pInFile, "EndFilePairs", inFileName);
   rewind(pInFile);
   FindToken(pInFile, "ModelExecutable", inFileName);   

   /*
   -------------------------------------------------------------
   Read in and create the directory from which the model will be 
   run.
   -------------------------------------------------------------
   */
   rewind(pInFile);
   if(CheckToken(pInFile, "ModelSubdir", inFileName) == true)
   {  
       line = GetCurDataLine(); 
       MyTrim(line);       
       if(strlen(line) < 12)
       { 
          LogError(ERR_IN_PARSE, "Bad ModelSubdir");
          ExitProgram(1);
       }       
       strcpy(tmp1, &line[11]);
       //strip whitespace
       MyTrim(tmp1);
       //strip quotes
       if(tmp1[0] == '"'){ tmp1[0] = ' ';}
       if(tmp1[strlen(tmp1)-1] == '"'){ tmp1[strlen(tmp1)-1] = ' ';}
       MyTrim(tmp1);
       strcpy(pDirName, tmp1);
   }

   if(pDirName[0] != '.')
   {
      MPI_Comm_rank(MPI_COMM_WORLD, &id);
      sprintf(tmp1, "%d", id);
      strcat(pDirName, tmp1);      
      sprintf(tmp1, "mkdir %s", pDirName);
      system(tmp1);      
   }/* end if() */

   /*
   -------------------------------------------
   Read in the model executable, modifying so 
   that output is redirected to a file. Also, 
   copy the executable file to the directory 
   of execution.
   -------------------------------------------
   */
   rewind(pInFile);
   FindToken(pInFile, "ModelExecutable", inFileName);
   line = GetCurDataLine();
   //line format = 'ModelExecutable   <var>'
   /*--------------------------------------------------------
   Read in executable, taking care to preserve full path, 
   even in the presence of long and space-separated filenames.
   --------------------------------------------------------*/   
   i = ExtractString(line, tmp2);
   i = ValidateExtraction(i, 1, 1, "Model()");
   i = ExtractFileName(&(line[i]), tmp1);

   //must wrap in quotes if there is whitespace in the execuable path
   quoteWrap = false;
   j = (int)strlen(tmp1);
   for(i = 0; i < j; i++){ if(tmp1[i] == ' ') {quoteWrap = true;}}     
   if(quoteWrap == true) { tmp1[j++] = '"';}
   tmp1[j] = (char)NULL;   
   if(quoteWrap == true) 
   { 
      MyStrRev(tmp1);
      tmp1[j++] = '"';
      tmp1[j] = (char)NULL;
      MyStrRev(tmp1);
   }

   /*-------------------------------------------------------------------------------------------
   Check to see if the model is internal. If it is, it will be a function call and end in a '()'
   --------------------------------------------------------------------------------------------*/
   if((strlen(tmp1) > 1) && (tmp1[strlen(tmp1)-2] == '(' ) && (tmp1[strlen(tmp1)-1] == ')' ))
   {
      m_InternalModel = true;
   }
   else
   {
      m_InternalModel = false;
   }

   /*------------------------------------------------------------------------------------------
   Extract executable file name from rest of path and save it. This is so that Ostrich can 
   cleanup after itself when copying files around.
   ------------------------------------------------------------------------------------------*/
   if(m_InternalModel == false)
   {
      j = 0;
      for(i = (int)(strlen(tmp1))-1; i > 0; i--)
      {
         if((tmp1[i] != '\\') && (tmp1[i] != '/'))
         {
            tmp2[j] = tmp1[i];
            j++;
         }
         else
         {
            if(quoteWrap == true)
            {
               tmp2[j] = '"';
               j++;
            }
            tmp2[j] = (char)NULL;
            MyStrRev(tmp2);
            m_pFileCleanupList->Insert(tmp2);
            break;
         }
      }
   }/* end if() */

   if((pDirName[0] != '.') && (m_InternalModel == false))
   {      
      #ifdef WIN32
         sprintf(tmp2, "copy %s %s", tmp1, pDirName);
      #else
         sprintf(tmp2, "cp %s %s", tmp1, pDirName);
      #endif
      system(tmp2);
   } /* end if() */

   if(m_InternalModel == false)
   {
      //make sure the executable exists
      strcpy(tmp2, tmp1);

      if(tmp2[0] == '"')
      {
         tmp2[0] = ' ';
         tmp2[strlen(tmp2)-1] = ' ';
         MyTrim(tmp2);
      }
      if(MY_ACCESS(tmp2, 0 ) == -1)
      {
         sprintf(tmp1, "Model executable (|%s|) not found", tmp2);
         LogError(ERR_FILE_IO, tmp1);
         ExitProgram(1);
      }

      #ifdef WIN32 //windows version
         strcat(tmp1, " > OstExeOut.txt");
      #else
         #ifdef IRIX_GMAKE  //SGI IRIX version
            strcpy(tmp2, tmp1);
            // '>&' redircts both output and error
            strcat(tmp2, " > OstExeOut.txt 2>&1"); 
            strcpy(tmp1, tmp2);
         #else // Linux version
            strcpy(tmp2, tmp1);
            // '>&' redircts both output and error
            strcat(tmp2, " > OstExeOut.txt"); 
            strcpy(tmp1, tmp2);
         #endif
      #endif
   }
   SetCmdToExecModel(tmp1);

   /*
   --------------------------------------------------------
   Read in the 'File Pairs': a set of template files and 
   their model equivalents.
   --------------------------------------------------------
   */
   rewind(pInFile);
   FindToken(pInFile, "BeginFilePairs", inFileName);
   line = GetNxtDataLine(pInFile, inFileName);
   while(strstr(line, "EndFilePairs") == NULL)
   {      
      if((strstr(line, ";") == NULL) && (strstr(line, "\t") == NULL))
      {
         LogError(ERR_FILE_IO, "Model::CTOR(): missing separator (;) in file pair.");
      }/* end if() */

      /*--------------------------------------------------------
      Read in file pairs, taking care to preserve full path, 
      even in the presence of long and space-separated filenames.
      --------------------------------------------------------*/
      //first file in pair.....
      i = ExtractFileName(line, tmp1);
      i = ExtractFileName(&(line[i]), tmp2);

      if(pDirName[0] != '.')
      {
         strcpy(tmp3, pDirName);
         #ifdef WIN32
            strcat(tmp3, "\\");
         #else
            strcat(tmp3, "/");
         #endif
         strcat(tmp3, tmp2);
         strcpy(tmp2, tmp3);
      }/* end if() */
      NEW_PRINT("FilePair", 1);
      pFilePair = new FilePair(tmp1, tmp2);
      MEM_CHECK(pFilePair);

      AddFilePair(pFilePair);

      line = GetNxtDataLine(pInFile, inFileName);
   }/* end while() */

   /*
   --------------------------------------------------------------------
   Read in any extra model files, these will need to be copied to the 
   model subdirectory.
   --------------------------------------------------------------------
   */
   rewind(pInFile);
   if(CheckToken(pInFile, "BeginExtraFiles", inFileName) == true)
   {
      //make sure end token exists
      FindToken(pInFile, "EndExtraFiles", inFileName);
      rewind(pInFile);
      FindToken(pInFile, "BeginExtraFiles", inFileName);

      line = GetNxtDataLine(pInFile, inFileName);
      while(strstr(line, "EndExtraFiles") == NULL)
      {
         //extra file
         ExtractFileName(line, tmp1);

         // add to cleanup list
         m_pFileCleanupList->Insert(tmp1);

         if(pDirName[0] != '.')
         {
            #ifdef WIN32
               sprintf(tmp2, "copy %s %s", tmp1, pDirName);
            #else
               sprintf(tmp2, "cp %s %s", tmp1, pDirName);
            #endif
            system(tmp2);

            strcpy(tmp2, pDirName);
            #ifdef WIN32
               strcat(tmp2, "\\");
            #else
               strcat(tmp2, "/");
            #endif
            strcat(tmp2, tmp1);
            strcpy(tmp1, tmp2);
         }/* end if() */

         line = GetNxtDataLine(pInFile, inFileName);
      }/* end while() */
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in parallel execution mode
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "ParallelModelExec", inFileName) == true)
   {
      line = GetCurDataLine();
      if(TestMyPBS() == true) //verify PBS is functional
      { 
         LogError(ERR_MODL_EXE, "Model will be run in parallel");
         m_SerialModelExec = false;
      }/* end if() */
      else //PBS not supported, don't allow parallel execution
      {
         LogError(ERR_MODL_EXE, "PBS not supported, can't run model in parallel");
         fclose(pInFile);
         ExitProgram(1);
      }/* end else() */
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Check for alternate objective function (default is WSSE)
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "ObjectiveFunction", inFileName) == true)
   {
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strstr(tmp2, "user") != NULL) {m_ObjFuncId = OBJ_FUNC_USER;}
      else if(strstr(tmp2, "sawe") != NULL) {m_ObjFuncId = OBJ_FUNC_SAWE;}
      else if(strstr(tmp2, "wsse") != NULL) {m_ObjFuncId = OBJ_FUNC_WSSE;}
      else if(strstr(tmp2, "pato") != NULL) {m_ObjFuncId = OBJ_FUNC_PATO;}
      else if(strstr(tmp2, "gcop") != NULL) {m_ObjFuncId = OBJ_FUNC_GCOP;}
   } /* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to check sensitivities
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "CheckSensitivities", inFileName) == true)
   {   
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) {m_bCheckGlobalSens = true;}
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to use surrogate models
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "SurrogateApproach", inFileName) == true)
   {   
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) {m_bUseSurrogates = true;}
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to use SuperMUSE
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   bool bSMUSE = false;
   if(CheckToken(pInFile, "SuperMUSE", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) 
      {
         EnableSuperMUSE();
         InitSuperMUSE(pInFile, (ModelABC *)this);
         bSMUSE = true;
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in flag to preserve model output files. This only applies if
   a model subdirectory is used.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "PreserveModelOutput", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) 
      {
         m_bPreserveModelOutput = true;
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in warm start flag. This only applies if we want to restart
   Ostrich from a previously aborted or interrupted run in which
   the OstModel0.txt file has been preserved.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "OstrichWarmStart", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) 
      {
         m_bWarmStart = true;
         RestoreRandomSeed();
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in caching flag. This only applies if we want to search the
   history of model evaluations (stored in OstModel0.txt) prior to
   running the model. If the candidate solution has already been
   evaluated, it will be found in the cache and so we won't have to
   run the model.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "OstrichCaching", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %s", tmp1, tmp2);
      MyStrLwr(tmp2);
      if(strncmp(tmp2, "yes", 3) == 0) 
      {
         m_bCaching = true;
      }
   }/* end if() */

   /*
   --------------------------------------------------------------------
   Read in number of digits of precision in I/O.
   --------------------------------------------------------------------
   */   
   rewind(pInFile);
   if(CheckToken(pInFile, "NumDigitsOfPrecision", inFileName) == true)
   {  
      line = GetCurDataLine(); 
      sscanf(line, "%s %d", tmp1, &m_Precision);
      if((m_Precision < 1) || (m_Precision > 32)) 
      {
         LogError(ERR_FILE_IO, "Invalid precision setting - defaulting to 6 digits.");
         m_Precision = 6;
      }
   }/* end if() */

   fclose(pInFile);

   if(bSMUSE) CleanSuperMUSE();

   /*
   --------------------------------------------------------
   Read in database conversion information.
   --------------------------------------------------------
   */
	DatabaseABC * access_DbaseList = new AccessConverter();
   DatabaseABC * netcdf_DbaseList = new NetCDFConverter();
	MEM_CHECK(access_DbaseList);
	MEM_CHECK(netcdf_DbaseList);
	if(access_DbaseList->ReadFromFile())
   {
		NEW_PRINT("AccessConverter", 1);
		m_DbaseList = access_DbaseList;
      delete netcdf_DbaseList;
		netcdf_DbaseList = NULL;
   }
	else if(netcdf_DbaseList->ReadFromFile())
	{
		NEW_PRINT("NetCDFConverter", 1);
		m_DbaseList = netcdf_DbaseList;
		delete access_DbaseList;
		access_DbaseList = NULL;
	}

   NEW_PRINT("ParameterGroup", 1);
   m_pParamGroup = new ParameterGroup();
   MEM_CHECK(m_pParamGroup);

   //check bounds on parameter group
   m_pParamGroup->CheckBounds();

   if((m_ObjFuncId == OBJ_FUNC_WSSE) || (m_ObjFuncId == OBJ_FUNC_SAWE))
   {
      NEW_PRINT("ObservationGroup", 1);
      m_pObsGroup = new ObservationGroup();
      MEM_CHECK(m_pObsGroup);
   }/* end if() */   

   //setup the objective function
   if(m_ObjFuncId == OBJ_FUNC_WSSE)
   { 
      NEW_PRINT("WSSE", 1);
      m_pObjFunc = new WSSE(m_pObsGroup);
   }
   else if(m_ObjFuncId == OBJ_FUNC_SAWE)
   { 
      NEW_PRINT("SAWE", 1);
      m_pObjFunc = new SAWE(m_pObsGroup);
   }
   else if(m_ObjFuncId == OBJ_FUNC_PATO)
   { 
      NEW_PRINT("PATO", 1);
      m_pObjFunc = new PATO(m_pParamGroup);
   }
   else if(m_ObjFuncId == OBJ_FUNC_GCOP)
   { 
      NEW_PRINT("GCOP", 1);
      m_pObjFunc = new GCOP(m_pParamGroup);
   }
   else /* OBJ_FUNC_USER) */ 
   { 
      NEW_PRINT("USER", 1);
      m_pObjFunc = new UserObjFunc("OstExeOut.txt");
   }
   MEM_CHECK(m_pObjFunc);

   /*-----------------------------------------------------------------------
   Read in special parameters.
   ------------------------------------------------------------------------*/
   m_pParamGroup->InitSpecialParams(inFileName);

   /*-----------------------------------------------------------------------
   Check template files against parameters, each parameter should appear in
   at least one template file or at least one database entry.
   ------------------------------------------------------------------------*/
   m_pParamGroup->CheckTemplateFiles(m_FileList);
   //m_pParamGroup->CheckDbaseFiles(m_DbaseList);

   /*-----------------------------------------------------------------------
   Check parameters for uniqueness, each parameter should be unique and 
   should not be a substring of another parameter.
   ------------------------------------------------------------------------*/
   m_pParamGroup->CheckMnemonics();

   /*-----------------------------------------------------------------------
   Initialize surrogate models, if the surrogate-based approach is enabled.
   ------------------------------------------------------------------------*/
   if(m_bUseSurrogates == true)
   {
      m_pDecision = new DecisionModule((ModelABC *)this);
   }

   CheckGlobalSensitivity();

   IncCtorCount();
} /* end default CTOR */

/******************************************************************************
destructor
******************************************************************************/
Model::~Model(void)
{
   Destroy();
} /* end DTOR */

/******************************************************************************
Free up memory.
******************************************************************************/
void Model::Destroy(void)
{
   delete m_pObsGroup;
   delete m_pParamGroup;
   if(m_pObjFunc != NULL){ m_pObjFunc->Destroy();}
   delete m_pObjFunc;
   delete m_FileList;
   delete m_DbaseList;
   delete [] m_ExecCmd;
   delete m_pDecision;

   if(m_pFileCleanupList != NULL)
   {
      IroncladString dirName = GetExeDirName(); 
      if(dirName[0] != '.')
      {
         m_pFileCleanupList->Cleanup(dirName);         
      }
      m_pFileCleanupList->Destroy();
   }
   IncDtorCount();
}/* end Destroy() */
   
/******************************************************************************
GetObjFuncPtr() 
   Returns a pointer to the objective function.
******************************************************************************/
ObjectiveFunction * Model::GetObjFuncPtr(void)
{
  return m_pObjFunc;
} /* end GetObjFuncPtr() */

/******************************************************************************
GetCounter()
   Returns the number of times the model has been executed
******************************************************************************/
int Model::GetCounter(void)
{
  return m_Counter;
} /* end GetCounter() */

/******************************************************************************
SetCmdToExecModel()
   Sets the syntax which is used to execute the model
*******************************************************************************/
void Model::SetCmdToExecModel(IroncladString cmd)
{
   int len;
   len = (int)strlen(cmd) + 1;
   NEW_PRINT("char", len);
   m_ExecCmd = new char[len];
   MEM_CHECK(m_ExecCmd);

   strcpy(m_ExecCmd, cmd);
} /* end SetCmdToExecModel() */

/******************************************************************************
AddFilePair()
   Adds a file pair to the model file pair list.
******************************************************************************/
void Model::AddFilePair(FilePair * pFilePair)
{
   if(m_FileList == NULL) { m_FileList = pFilePair;}
   else{m_FileList->InsertPair(pFilePair);}
} /* end AddFilePair() */

/******************************************************************************
AddDatabase()
   Adds a database conversion to the list.
******************************************************************************/
void Model::AddDatabase(DatabaseABC * pDbase)
{
   if(m_DbaseList == NULL) { m_DbaseList = pDbase;}
   else{m_DbaseList->InsertDbase(pDbase);}
} /* end AddDatabase() */

/******************************************************************************
GetObsGroupPtr()
   Returns the observation group pointer
******************************************************************************/
ObservationGroup * Model::GetObsGroupPtr(void)
{
  return m_pObsGroup;
} /* end GetObsGroupPtr() */

/*****************************************************************************
GetParamGroupPtr()
   Returns the parameter group pointer.
******************************************************************************/
ParameterGroup * Model::GetParamGroupPtr(void)
{
  return m_pParamGroup;
} /* end GetParamGroupPtr() */

/*****************************************************************************
Execute()
   Executes the model (or surrogate) and returns the objective function value.
******************************************************************************/
double Model::Execute(void)
{  
   if(m_bUseSurrogates == false)
   {
      return(StdExecute(0.00));
   } 
   else
   {
      return(m_pDecision->Execute());
   }
} /* end Execute() */

/*****************************************************************************
Execute()
   Executes the model (or surrogate) and returns the objective function value.
   Incorporates a penalty function for violation of parameter bounds.

   NOT compatible with surrogate-based appraoch.
******************************************************************************/
double Model::Execute(double viol)
{
	return(StdExecute(viol));
} /* end Execute() */

/******************************************************************************
PerformWarmStartCorrection()
   If WarmStarts are enabled and an algorithm has gotten off track somehow,
   then correct the current parameter set to match the parameter set in the
   OstModel0.txt file for the next evalaution.

   If WarmStarts are not enabled, then do nothing.
******************************************************************************/
void Model::PerformWarmStartCorrection(void)
{  
   if(m_bWarmStart == false)
   {
      return;
   } 
   else
   {
      int i, j, count;
      char line[DEF_STR_SZ], valStr[DEF_STR_SZ];
      char * pTmp;
      double objFunc, paramVal;
      bool bFound;

      ParameterABC * pParam;
      ParameterGroup * pGroup = GetParamGroupPtr();
      if(pGroup == NULL) return;
      int np = pGroup->GetNumParams();
      FILE * pIn = fopen("OstModel0.txt", "r");
      if(pIn == NULL) return;
      while(!feof(pIn))
      {
         fgets(line, DEF_STR_SZ, pIn);
         pTmp = line;
         MyTrim(pTmp);
         j = ExtractColString(pTmp, valStr, ' ');
         count = atoi(valStr);
         pTmp += j;
         MyTrim(pTmp);
         j = ExtractColString(pTmp, valStr, ' ');
         pTmp += j;
         MyTrim(pTmp);
         objFunc = atof(valStr);

         bFound = true;
         for(i = 0; i < np; i++)
         {
            if(count != (m_Counter+1))
            {
               bFound = false;
               break;
            }
            //found a match! Extract parameter values
            pParam = pGroup->GetParamPtr(i);
            j = ExtractColString(pTmp, valStr, ' ');
            paramVal = atof(valStr);
            pParam->SetEstVal(paramVal);

            pTmp += j;
            MyTrim(pTmp);
         }/* end for() */
         if(bFound == true)
         {
            fclose(pIn);
            return;
         }
      }/* end while() */
      fclose(pIn);
      return;
   }/* end else() */
} /* end PerformWarmStartCorrection() */

/*****************************************************************************
Bookkeep()
   Performs bookkeeping operations related to parallel executing.
   if bFinal == true, iteration is complete so collect metrics
   if bFinal == false, in the middle of iteration, share information between
    processors
******************************************************************************/
void Model::Bookkeep(bool bFinal)
{  
   int id, nprocs, temp, i;

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (nprocs == 1) return;

   if(m_bUseSurrogates == true)
   {
      m_pDecision->Bookkeep(bFinal);
   }   

   /* ---------------------------------------------------------------
   Collect total evals 
   --------------------------------------------------------------- */
   if(bFinal == true)
   {
      for(i = 1; i < nprocs; i++)
      {
         temp = m_Counter;
         MPI_Bcast(&temp, 1, MPI_INTEGER, i, MPI_COMM_WORLD);
         if(id == 0) m_Counter += temp;
      }
   }
} /* end Bookkeep() */

/*****************************************************************************
StdExecute()
   Executes the standard (complex) model and returns the objective function 
   value.
******************************************************************************/
double Model::StdExecute(double viol)
{   
   IroncladString dirName = GetExeDirName();
   FilePair * pCur;
   FilePipe * pPipe;
   double val;
   int num_procs, walltime, result;
   bool isGoodTopo;
   char tmp[DEF_STR_SZ];

   //exit early if the user has requested program termination
   if(IsQuit() == true){ return NEARLY_HUGE;}

   //inc. number of times model has been executed
   m_Counter++;

   /*
   adjust geometries to conform to topology rules, 
   if the topology cannot be 'fixed' then false will 
   be returned....
   */
   isGoodTopo = m_pParamGroup->FixGeometry();
   if(isGoodTopo == false)
   {
      LogError(ERR_MODL_EXE, "Could not correct model topology");
   }

   //make substitution of parameters into model input file
   pCur = m_FileList;
   while(pCur != NULL)
   {
      pPipe = pCur->GetPipe();
      m_pParamGroup->SubIntoFile(pPipe);
      pCur = pCur->GetNext();
   } /* end while() */

   //cd to model subdirectory, if needed
   if(dirName[0] != '.') { MY_CHDIR(dirName);}   

   //make substitution of parameters into model input databases
   if(m_DbaseList != NULL)
   {
      m_pParamGroup->SubIntoDbase(m_DbaseList);
   } /* end while() */

   /* -----------------------------------------------------------
   If warm start is enabled, attempt to read model evaluation
   from OstModel0.txt file.
   ----------------------------------------------------------- */
   if(m_bWarmStart == true)
   {      
      bool bCached = WarmStart(&val);
      if(bCached == true) //found previous model result
      {
         m_CurObjFuncVal = val;
         return (val);
      }
      else //couldn't find previous result, must be in new territory
      {
         Write(0.00);
         m_bWarmStart = false; //no more cached results
      }
   }

   /* -----------------------------------------------------------
   If caching is enabled, attempt to read model evaluation
   from OstModel0.txt file.
   ----------------------------------------------------------- */
   if(m_bCaching == true)
   {      
      bool bCached = CheckCache(&val);
      if(bCached == true) //found previous model result
      {
         m_NumCacheHits++;
         Write(val);
         m_CurObjFuncVal = val;
         return (val);
      }
   }/* end if() */

   if(m_InternalModel == true)
   {
      if(strcmp(m_ExecCmd, "Isotherm()") == 0){ Isotherm();}
      else if(strcmp(m_ExecCmd, "Orear()") == 0){ Orear();}
      else if(strcmp(m_ExecCmd, "McCammon()") == 0){ McCammon();}
      else if(strcmp(m_ExecCmd, "Kinniburgh()") == 0){ Kinniburgh();}
      else if(strcmp(m_ExecCmd, "AdvancedKinniburgh()") == 0){ AdvancedKinniburgh();}
      else{LogError(ERR_BAD_ARGS, "Unknown internal model"); ExitProgram(1);}
   }
   else
   {
      if(m_SerialModelExec == true) //run serial model
      {
         //invoke system command to execute the model   
         system(m_ExecCmd);
      }/* end if() */
      else //run model in parallel
      {         
         num_procs = 2;
         walltime = 10; //hours
         printf("Setting up script to request %d processors for %d hours\n", num_procs, walltime);
         //MyPBS_CreateParaScript(pSched);
         //MyPBS_qsub(pScript->ScriptFileName);
         //result = WaitForJobCompletion(pScript);
         result = JOB_FAILED;
         //error check
         if(result == JOB_TIMED_OUT)
         {
            LogError(ERR_MODL_EXE, "Parallel Job Timed Out");
            ExitProgram(1);
         }/* end if() */
         if(result == JOB_FAILED)
         {
            LogError(ERR_MODL_EXE, "Parallel Job Failed");
            ExitProgram(1);
         }/* end if() */
      }/* end else() */
   }/* end else (external model) */

   //extract computed reponses from model output database(s)
   if(m_DbaseList != NULL)
   {
      DatabaseABC * pCur;

      //first clean up ASCII files
      for(pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
      {
         pCur->DeleteASCIIFile();
      }

      //convert the responses to ASCII file
      for(pCur = m_DbaseList; pCur != NULL; pCur = pCur->GetNext())
      {
         pCur->ReadResponse();
      }
   }

   //extract computed observations from model output file(s)
   if(m_pObsGroup != NULL){ m_pObsGroup->ExtractVals();}

   //compute obj. func.
   val = m_pObjFunc->CalcObjFunc();

   //add in penalty for violation of parameter bounds
   val += viol*MyMax(1.00,val);

      //preserve model output, if desired
   if((dirName[0] != '.') && (m_bPreserveModelOutput == true))
   {
      #ifdef WIN32
         sprintf(tmp, "mkdir ..\\run%d", m_Counter); //use parent dir for temporary location, due to xcopy rules
         system(tmp);
         sprintf(tmp, "dir /B run* > Exclude.txt"); //need to exclude previous 'run' directories
         system(tmp); 
         sprintf(tmp, "xcopy * ..\\run%d /S /EXCLUDE:Exclude.txt >> OstExeOut.txt", m_Counter);  //perform copy
         system(tmp);
         sprintf(tmp, "move ..\\run%d .\\ >> OstExeOut.txt", m_Counter);  //relocate the directory
         system(tmp);
      #else
         sprintf(tmp, "mkdir run%d", m_Counter);
         system(tmp);
         sprintf(tmp, "cp * run%d 2>&1 | >> OstExeOut.txt", m_Counter);
         system(tmp);
      #endif

      sprintf(tmp, "run%d", m_Counter);
      m_pFileCleanupList->Cleanup(tmp);
   }

   //cd out of model subdirectory, if needed
   if(dirName[0] != '.') { MY_CHDIR("..");}
   
   //ouput results
   Write(val);

   m_CurObjFuncVal = val;

   return (val);
} /* end StdExecute() */

/******************************************************************************
WarmStart()
   Attempt to read previous result stored in OstModel0.txt file. If successful
   the corresponding obj. function value will be stored in val argument. Other-
   wise the function returns false and val is unchanged.

   Match must be exact for all parameters and for the model evaluation counter.
   Use CheckCache() if only a match on parameters (and not model eval. counter)
   is required.
******************************************************************************/
bool Model::WarmStart(double * val)
{
   int i, j, count;
   char line[DEF_STR_SZ], valStr[DEF_STR_SZ];
   char * pTmp;
   double objFunc, paramVal, pval, test, fmax;
   bool bFound;

   ParameterABC * pParam;
   ParameterGroup * pGroup = GetParamGroupPtr();
   if(pGroup == NULL) return false;
   int np = pGroup->GetNumParams();
   FILE * pIn = fopen("OstModel0.txt", "r");
   if(pIn == NULL) return false;

   //skip first m_Counter lines
   for(i = 0; i < m_Counter; i++)
   {
      fgets(line, DEF_STR_SZ, pIn); 
      if(feof(pIn)) break;
   }
  
   while(!feof(pIn))
   {
      fgets(line, DEF_STR_SZ, pIn);
      pTmp = line;
      MyTrim(pTmp);
      j = ExtractColString(pTmp, valStr, ' ');
      count = atoi(valStr);
      pTmp += j;
      MyTrim(pTmp);
      j = ExtractColString(pTmp, valStr, ' ');
      pTmp += j;
      MyTrim(pTmp);
      objFunc = atof(valStr);

      bFound = true;
      for(i = 0; i < np; i++)
      {
         if(count != m_Counter)
         {
            bFound = false;
            //break;
            fclose(pIn);
            return false;
         }
         pParam = pGroup->GetParamPtr(i);
         j = ExtractColString(pTmp, valStr, ' ');
         paramVal = atof(valStr);
         pval = pParam->GetEstVal();

         fmax = 1.00E-10;
         if(fabs(paramVal) > fmax) fmax = fabs(paramVal);
         if(fabs(pval) > fmax) fmax = fabs(pval);
         test=fabs(paramVal-pval)/fmax;
         if(test > 1E-6)
         {
            bFound = false;
            //break;
            fclose(pIn);
            return false;
         }/* end if() */
         pTmp += j;
         MyTrim(pTmp);
      }/* end for() */
      if(bFound == true)
      {
         fclose(pIn);
         *val = objFunc;
         return true;
      }
   }/* end while() */
   fclose(pIn);
   return false;
}/* end WarmStart() */

/******************************************************************************
CheckCache()
   Attempt to read previous result stored in OstModel0.txt file. If successful
   the corresponding obj. function value will be stored in val argument. Other-
   wise the function returns false and val is unchanged.
******************************************************************************/
bool Model::CheckCache(double * val)
{
   int i, j;
   char line[DEF_STR_SZ], valStr[DEF_STR_SZ];
   char * pTmp;
   double objFunc, paramVal;
   bool bFound;

   ParameterABC * pParam;
   ParameterGroup * pGroup = GetParamGroupPtr();
   if(pGroup == NULL) return false;
   int np = pGroup->GetNumParams();
   FILE * pIn = fopen("OstModel0.txt", "r");
   if(pIn == NULL) return false;
   while(!feof(pIn))
   {
      fgets(line, DEF_STR_SZ, pIn);
      pTmp = line;
      MyTrim(pTmp);
      j = ExtractColString(pTmp, valStr, ' ');
      pTmp += j;
      MyTrim(pTmp);
      j = ExtractColString(pTmp, valStr, ' ');
      pTmp += j;
      MyTrim(pTmp);
      objFunc = atof(valStr);

      bFound = true;
      for(i = 0; i < np; i++)
      {
         pParam = pGroup->GetParamPtr(i);
         j = ExtractColString(pTmp, valStr, ' ');
         paramVal = atof(valStr);
         if(fabs(paramVal - pParam->GetEstVal()) > 1E-10)
         {
            bFound = false;
            break;
         }/* end if() */
         pTmp += j;
         MyTrim(pTmp);
      }/* end for() */
      if(bFound == true)
      {
         fclose(pIn);
         *val = objFunc;
         return true;
      }
   }/* end while() */
   fclose(pIn);
   return false;
}/* end CheckCache() */

/******************************************************************************
GatherTask()
   Read output file of a SuperMUSE task (stored in pDir directory) and compute 
   the associated objective function.
******************************************************************************/
double Model::GatherTask(char * pDir)
{
  double val;

  //inc. number of times model has been executed
  m_Counter++;

  //cd to task subdirectory
  MY_CHDIR(pDir);

  //extract computed observations from model output file(s)
  if(m_pObsGroup != NULL){ m_pObsGroup->ExtractVals();}

  //compute obj. func.
  val = m_pObjFunc->CalcObjFunc();

  //cd out of task subdirectory
  MY_CHDIR("..");

  //ouput results
  Write(val);

  m_CurObjFuncVal = val;

  return (val);
}/* end GatherTask() */

/******************************************************************************
Write()
   Store parameter and objective function value to model output file.
******************************************************************************/
void Model::Write(double objFuncVal)
{
   FILE * pFile;
   char name[DEF_STR_SZ];
   int id;
   static bool firstCall = true;

   //Write() is only called once when WarmStart is active
   //and this is done to signal that subsequent Writes() 
   //should resume where WarmStart left off.
   if(m_bWarmStart == true)
   {
      firstCall = false;
      return;
   }

   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   sprintf(name, "OstModel%d.txt", id);

   if(firstCall == true) 
   {
      pFile = fopen(name, "r");
      if(pFile != NULL)
      {
         fclose(pFile);
         if(remove(name) != 0)
         {
            LogError(ERR_FILE_IO, "Write(): Couldn't delete OstModel.txt file");
            ExitProgram(1);
         }
      }
      firstCall = false;
      //write out banner.
      pFile = fopen(name, "a+");
      if(pFile == NULL)
      {
         LogError(ERR_FILE_IO, "Write(): Couldn't open OstModel.txt file");
         ExitProgram(1);
      }
      fprintf(pFile,"Run   obj.function   ");
      if(m_pObsGroup != NULL) m_pObsGroup->Write(pFile, WRITE_BNR);
      m_pParamGroup->Write(pFile, WRITE_BNR);
      fprintf(pFile,"\n");
      fclose(pFile);
   }

   pFile = fopen(name, "a+");
	fprintf(pFile, "%-4d  ", m_Counter);
   WritePreciseNumber(pFile, objFuncVal);
   fprintf(pFile, "  ");
   if(m_pObsGroup != NULL) m_pObsGroup->Write(pFile, WRITE_SCI);
   m_pParamGroup->Write(pFile, WRITE_SCI);
   fprintf(pFile, "\n");
   fclose(pFile);
} /* end Write() */

/******************************************************************************
WriteMetrics()
******************************************************************************/
void Model::WriteMetrics(FILE * pFile)
{  
   if(m_bUseSurrogates == true)
   {
      m_pDecision->WriteMetrics(pFile);
   }
   else
   {
      fprintf(pFile, "Total Evals             : %d\n", m_Counter);
      if(m_bCaching == true)
         fprintf(pFile, "Cache Hits              : %d\n", m_NumCacheHits);
   }
} /* end WriteMetrics() */

/******************************************************************************
CheckGlobalSensitivity()

Checks that each observation is sensitive to at least one parameter over the
range of possible parameter values. If an observation is not sensitive to any
parameters, a warning will be reported and the observation will be ignored in 
the calibration.

Also checks the sensitivity of each parameter. If a given parameter does not
affect the objective function over the entire parameter range, a warning will
be reported and the parameter will be ignored in the calibration.
******************************************************************************/
void Model::CheckGlobalSensitivity(void)
{
   if((m_pObsGroup == NULL) || (m_pParamGroup == NULL)) return;
   if(m_bCheckGlobalSens == false) return;

   int i, j, nobs, nprm; 
   double upr, lwr, Fupr, Flwr;
   double * pInit, * obs_sum, * prm_sum;
   double * ObsUpr, * ObsLwr;
   char tmp1[DEF_STR_SZ];
   UnchangeableString * obs_names, * prm_names;

   nobs = m_pObsGroup->GetNumObs();
   obs_sum = new double[nobs];
   ObsUpr = new double[nobs];
   ObsLwr = new double[nobs];
   obs_names = new UnchangeableString[nobs];

   nprm = m_pParamGroup->GetNumParams();
   pInit = new double[nprm];
   prm_sum = new double[nprm];
   prm_names = new UnchangeableString[nprm];

   //save intitial parameters
   m_pParamGroup->ReadParams(pInit);
   
   for(i = 0; i < nobs; i++)
   {
      obs_sum[i] = 0.00;
      obs_names[i] = m_pObsGroup->GetObsPtr(i)->GetName();
   }
   
   for(j = 0; j < nprm; j++)
   {
      prm_names[j] = m_pParamGroup->GetParamPtr(j)->GetName();
      upr = m_pParamGroup->GetParamPtr(j)->GetUprBnd();
      m_pParamGroup->GetParamPtr(j)->SetEstVal(upr);
      Fupr = Execute();

      //store computed observation values
      for(i = 0; i < nobs; i++)
      {
         ObsUpr[i] = m_pObsGroup->GetObsPtr(i)->GetComputedVal();
      }

      lwr = m_pParamGroup->GetParamPtr(j)->GetLwrBnd();
      m_pParamGroup->GetParamPtr(j)->SetEstVal(lwr);
      Flwr = Execute();

      //store computed observation values
      for(i = 0; i < nobs; i++)
      {
         ObsLwr[i] = m_pObsGroup->GetObsPtr(i)->GetComputedVal();
      }

      //store change in objective function
      prm_sum[j] = fabs(Fupr - Flwr);

      //accumulate changes in observation values
      for(i = 0; i < nobs; i++)
      {
         obs_sum[i] += fabs(ObsUpr[i] - ObsLwr[i]);
      }

      //restore parameter values
      m_pParamGroup->WriteParams(pInit);
   }/* end for() */

   //perform parameter sensitivity checks
   for(j = 0; j < nprm; j++)
   {
      if(prm_sum[j] <= NEARLY_ZERO)
      {
         sprintf(tmp1, "%s appears to be insensitive and has been set to a constant value", prm_names[j]);
         LogError(ERR_INS_PARM, tmp1);
         m_pParamGroup->ExcludeParam(prm_names[j]);
      }
   }

   //perform observation sensitivity checks: phase 1, log the warning messages
   for(i = 0; i < nobs; i++)
   {
      if(obs_sum[i] <= NEARLY_ZERO)
      {
         sprintf(tmp1, "%s appears to be insensitive and has been excluded from the calibration", obs_names[i]);
         LogError(ERR_INS_OBS, tmp1);
         m_pObsGroup->ExcludeObs(obs_names[i]);
      }
   }

   delete [] pInit;
   delete [] obs_sum;
   delete [] prm_sum;
   delete [] ObsUpr;
   delete [] ObsLwr;
   delete [] prm_names;
   delete [] obs_names;
} /* end CheckGlobalSensitivity() */

