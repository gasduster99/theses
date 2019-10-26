/******************************************************************************
File      : Utility.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

This file contains a bunch of useful c-style routines, ranging from matrix 
mathematics to string manipulation.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
10-27-03    lsm   Removed all input files except model input file.
11-25-03    lsm   Added string reverse function
03-05-04    lsm   added PSO, fixed filename parse error in obs. parse
                  added support for alternative tokens
03-24-04    lsm   added PSO-LevMar hybrid
07-08-04    lsm   GetNxtDataLine() will skip over blank lines
08-17-04    lsm   Added string replace and string occurence functions.
                  RAM fragmentation fixes.
11-30-04    lsm   Replaced calls to time() with MyTime()
11-07-05    lsm   Added MyRand() which takes into account processor type and
                  RAND_MAX value. Also added support for GRID, BGA, VSA and CSA
                  programs.
01-01-07    lsm   Replaced recursive calculation of determinant with LU decom-
                  position. Added support for temporary input files which store
                  copies of the surrogate sections of the input file.
******************************************************************************/
#include "mpi_stub.h"
#include "AccessConverter.h"
#include "ADOConnection.h"
#include "NumericalRecipes.h"
#include "StatUtility.h"
#include "Utility.h"
#include "Exception.h"
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/*
Global strings, these strings are set by command line arguments 
or by the Ostrich configuration file.
*/
char gOstFileName[DEF_STR_SZ]; //set by cmd line (or defaults to "ostIn.txt")
char gExeDirName[DEF_STR_SZ];  //set by cfg file (or defaults to ".")
char gInFileName[DEF_STR_SZ]; //set internally, main configuration file
char gSrgFileName[DEF_STR_SZ]; //set internally, surrogate module configuration file
char gDynFileName[DEF_STR_SZ]; //set internally, surrogate model configuration file

/* ----------------
Global random seed.
------------------ */
bool         gSeedInitialized = false;
unsigned int gRandomSeed = 258;

ProgramType gProgChoice;

/* 
Global storage for line of input.
*/
char * gLine = NULL;
int gLineSize = 0;

/******************************************************************************
SetProgramType()

Sets the type of the program.
******************************************************************************/
void SetProgramType(ProgramType progVal)
{
   gProgChoice = progVal;
}/* end SetProgramType() */

/******************************************************************************
GetProgramType()

Returns the program type.
******************************************************************************/
ProgramType GetProgramType(void)
{
   return gProgChoice;
}/* end GetProgramType() */

/******************************************************************************
GetRandomSeed()

Returns the seed for the random number generator.
******************************************************************************/
unsigned int GetRandomSeed(void)
{
   if(gSeedInitialized == false)
   {
      gRandomSeed = ReadRandomSeed();
      gSeedInitialized = true;
   }
   return gRandomSeed;
}/* end GetRandomSeed() */

/******************************************************************************
RestoreRandomSeed()

Restores the seed for the random number generator by reading an existing
OstOutput0.txt file.
******************************************************************************/
void RestoreRandomSeed(void)
{
   char line[DEF_STR_SZ];
   char * tok = "Seed for Random Nums.  :";
   char * pTmp;
   FILE * pIn = fopen("OstOutput0.txt", "r");
   if(pIn == NULL) return;
   while(!feof(pIn))
   {
      fgets(line, DEF_STR_SZ, pIn);
      
      if(strncmp(line, tok, strlen(tok)) == 0)
      {
         pTmp = &line[strlen(tok)];
         MyTrim(pTmp);
         gRandomSeed = atoi(pTmp);
         gSeedInitialized = true;
         fclose(pIn);
         return;
      }/* end if() */
   }/* end while() */
   fclose(pIn);   
}/* end RestoreRandomSeed() */

/******************************************************************************
ReadRandomSeed()

Reads the seed for the random number generator from the input file.
******************************************************************************/
unsigned int ReadRandomSeed(void)
{
   int id;
   FILE * pFile;
   char tmpStr[DEF_STR_SZ];
   char * line;

   pFile = fopen(GetInFileName(), "r");
   if(pFile == NULL) 
   {
      FileOpenFailure("ReadRandomSeed()", GetInFileName());
   }/* end if() */

   if(CheckToken(pFile, "RandomSeed", GetInFileName()) == true)
   {
      line = GetCurDataLine();
      sscanf(line, "%s %d", tmpStr, &gRandomSeed);
   }/* end if() */
   else
   {
      MPI_Comm_rank(MPI_COMM_WORLD, &id);
      gRandomSeed = MyTime() + id;
   }/* end else() */

   fclose(pFile);
   return gRandomSeed;
}/* end ReadRandomSeed() */

/******************************************************************************
ConvertToASCII()

Converts the designated files in ostIn.txt to their ASCII counterparts.
******************************************************************************/
void ConvertToASCII(void)
{
   int j;
   char tmpFileType[DEF_STR_SZ];
   char tmpFileName[DEF_STR_SZ];
   char * lineStr;
   FILE * pFile;

   IroncladString pFileName = GetOstFileName();
   pFile = fopen(pFileName, "r");

   if(pFile == NULL) 
   {
      FileOpenFailure("ConvertToASCII()", pFileName);
   }/* end if() */

   //make sure correct tokens are present
   FindToken(pFile, "BeginTypeConversion", pFileName);
   FindToken(pFile, "EndTypeConversion", pFileName);
   rewind(pFile);

   FindToken(pFile, "BeginTypeConversion", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);

   while(strstr(lineStr, "EndTypeConversion") == NULL)
   {      
      j = ExtractString(lineStr, tmpFileType);
	  lineStr += j;
      j = ExtractString(lineStr, tmpFileName);
      lineStr += j;
	   //deletes converted file, if it exists
      #ifdef WIN32
        string s_fileName(tmpFileName);
        s_fileName = s_fileName.substr(0, s_fileName.find_last_of('.'));
        s_fileName += ".txt";
        remove(s_fileName.c_str());
      #endif

      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */

   rewind(pFile);
   FindToken(pFile, "BeginTypeConversion", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);

   while(strstr(lineStr, "EndTypeConversion") == NULL)
   {      
      j = ExtractString(lineStr, tmpFileType);
	  lineStr += j;
	  //determine conversion type
	  if (strncmp(tmpFileType, "Access", 6) == 0)
	  {
	     AccessConverter * pConverter = new AccessConverter();
		  pConverter->Initialize(lineStr);
		  pConverter->Convert();
        delete pConverter;
	  }

      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */
   fclose(pFile);
}/* end ConvertToASCII */

/******************************************************************************
ReadProgramType()

Reads the program type from the input file.
******************************************************************************/
ProgramType ReadProgramType(void)
{
   char * line;
   FILE * pFile;
   char tmpType[DEF_STR_SZ];
   char tmpStr[DEF_STR_SZ];

   pFile = fopen(GetInFileName(), "r");
   if(pFile == NULL) 
   {
      FileOpenFailure("ReadProgramType()", GetInFileName());
   }/* end if() */

   if(CheckToken(pFile, "ProgramType", GetInFileName()) == true)
   {
      line = GetCurDataLine();
      sscanf(line, "%s %s", tmpStr, tmpType);
      MyStrLwr(tmpType);
      MyTrim(tmpType);
      if(strncmp(tmpType, "geneticalgorithm", 16) == 0)
      { gProgChoice = GA_PROGRAM;}
      else if(strncmp(tmpType, "binarygeneticalgorithm", 22) == 0)      
      { gProgChoice = BGA_PROGRAM;}
      else if(strncmp(tmpType, "shuffledcomplexevolution", 24) == 0)      
      { gProgChoice = SCEUA_PROGRAM;}
      else if(strncmp(tmpType, "bisectionalgorithm", 18) == 0)      
      { gProgChoice = BIS_PROGRAM;}
      else if(strncmp(tmpType, "samplingalgorithm", 17) == 0)      
      { gProgChoice = SMP_PROGRAM;}
      else if(strncmp(tmpType, "psodesc-gml", 11) == 0)      
      { gProgChoice = PSODESC_LEV_PROGRAM;}
      else if(strncmp(tmpType, "psodesc", 7) == 0)      
      { gProgChoice = PSODESC_PROGRAM;}
      else if(strstr(tmpType, "particleswarm") != NULL) 
      { gProgChoice = PSO_PROGRAM;}
      else if(strncmp(tmpType, "ccpso", 5) == 0) 
      { gProgChoice = CCPSO_PROGRAM;}
      else if(strncmp(tmpType, "ccrga", 5) == 0) 
      { gProgChoice = CCRGA_PROGRAM;}
      else if(strncmp(tmpType, "pso-gml",7) == 0) 
      { gProgChoice = PSO_LEV_PROGRAM;}
      else if(strncmp(tmpType, "ccpso-gml", 9) == 0) 
      { gProgChoice = CCPSO_LEV_PROGRAM;}
      else if(strncmp(tmpType, "simulatedannealing", 18) == 0) 
      { gProgChoice = SA_PROGRAM;}
      else if(strncmp(tmpType, "ccsa", 4) == 0) 
      { gProgChoice = CCSA_PROGRAM;}
      else if(strncmp(tmpType, "ccvsa", 5) == 0) 
      { gProgChoice = CCVSA_PROGRAM;}
      else if(strncmp(tmpType, "discretesimulatedannealing", 26) == 0) 
      { gProgChoice = CSA_PROGRAM;}
      else if(strncmp(tmpType, "vanderbiltsimulatedannealing", 28) == 0) 
      { gProgChoice = VSA_PROGRAM;}
      else if(strncmp(tmpType, "levenberg-marquardt", 19) == 0) 
      { gProgChoice = LEV_PROGRAM;}
      else if(strncmp(tmpType, "gml-ms", 6) == 0) 
      { gProgChoice = GMLMS_PROGRAM;}
      else if(strncmp(tmpType, "powell", 6) == 0) 
      { gProgChoice = POWL_PROGRAM;}
      else if(strncmp(tmpType, "steepest-descent", 16) == 0) 
      { gProgChoice = STEEP_PROGRAM;}
      else if(strncmp(tmpType, "fletcher-reeves", 15) == 0) 
      { gProgChoice = FLRV_PROGRAM;}
      else if(strncmp(tmpType, "regressionstatistics", 20) == 0) 
      { gProgChoice = STATS_PROGRAM;}
      else if(strncmp(tmpType, "modelevaluation", 15) == 0)      
      { gProgChoice = EVAL_PROGRAM;}
      else if(strncmp(tmpType, "gridalgorithm", 13) == 0) 
      { gProgChoice = GRID_PROGRAM;}
      else if(strncmp(tmpType, "dds", 3) == 0) 
      { gProgChoice = DDS_PROGRAM;}
      else if(strncmp(tmpType, "discretedds", 11) == 0) 
      { gProgChoice = DDDS_PROGRAM;}
      else if(strncmp(tmpType, "glue", 4) == 0) 
      { gProgChoice = GLUE_PROGRAM;}
      else if(strncmp(tmpType, "rejectionsampler", 16) == 0) 
      { gProgChoice = RJSMP_PROGRAM;}
      else if(strncmp(tmpType, "metropolissampler", 16) == 0) 
      { gProgChoice = METRO_PROGRAM;}
      else 
      { 
         LogError(ERR_FILE_IO, 
                  "Unknown program type, defaulting to Levenberg-Marquardt");
         gProgChoice = LEV_PROGRAM;
      }
   }/* end if() */
   else
   {
      LogError(ERR_FILE_IO, 
               "No program type, defaulting to Levenberg-Marquardt");
      gProgChoice = LEV_PROGRAM;
   }/* end else() */

   fclose(pFile);
   return gProgChoice;
}/* end ReadProgramType() */

/******************************************************************************
GetOstFileName()

Retrieves a pointer to a global string that is the name of the Ostrich input 
file (usually OstIn.txt).
******************************************************************************/
UnmoveableString GetOstFileName(void)
{
   return gOstFileName;
}/* end GetOstFileName() */

/******************************************************************************
GetExeDirName()

Retrieves a pointer to a global string that is the name of the directory from
which the model is to be executed.
******************************************************************************/
UnmoveableString GetExeDirName(void)
{
   return gExeDirName;
}/* end GetExeDirName() */

/******************************************************************************
CompDbl()

Compares 2 doubles, according to format specified by qsort()
******************************************************************************/
int CompDbl(UnchangeableVoidPtr arg1, UnchangeableVoidPtr arg2)
{
   double v1;
   double v2;

   v1 = *(double *)arg1;
   v2 = *(double *)arg2;

   if(v1 < v2) return -1;
   if(v1 > v2) return 1;
   return 0;
}/* end CompDbl() */

/******************************************************************************
MyStrLwr()

Converts a string tolower case
******************************************************************************/
void MyStrLwr(UnmoveableString pLine)
{
   int i, size;

   size = (int)strlen(pLine);

   for(i = 0; i < size; i++){ pLine[i] = (char)(tolower((int)(pLine[i])));}
}/* end MyStrLwr() */

/******************************************************************************
MyStrRev()

Reverses a string.
******************************************************************************/
void MyStrRev(UnmoveableString pLine)
{
   int i, j, len;
   char tmp[DEF_STR_SZ];

   len = (int)strlen(pLine);
   if(len > DEF_STR_SZ) len = DEF_STR_SZ;

   for(i = 0; i < len; i++)
   {
      j = len - 1 - i;
      tmp[i] = pLine[j];
   }/* end for() */  
   tmp[i] = (char)NULL;
   strcpy(pLine, tmp); 
}/* end MyStrRev() */

/******************************************************************************
MyTrim()

Removes leading and trailing whitespace from a string.
******************************************************************************/
void MyTrim(UnmoveableString pLine)
{
   int i, len;

   len = (int)strlen(pLine);

   if(len <= 0){ return;}

   for(i = (len-1); i >= 0; i--)
   {
      if(IsWhitespace(pLine[i]) == true){ pLine[i] = (char)NULL;}
      else{ break;}
   }/* end for() */

   if(IsWhitespace(pLine[0]) == true)
   {
      MyStrRev(pLine);
      MyTrim(pLine);
      MyStrRev(pLine);
   }/* end if() */    
}/* end MyTrim() */

/******************************************************************************
MyStrRep()

Replaces all occurences of find with replace in the provided string. If 
replacement string is larger than the find string, it is up to the user to have 
allocated enough space in pStr.

Returns the number of replacements made.
******************************************************************************/
int MyStrRep(UnmoveableString pStr, IroncladString pFind, IroncladString pRep)
{
   int findLen, repLen, diff, count, end, start, i;
   char * tmp;

   findLen = (int)strlen(pFind);
   repLen  = (int)strlen(pRep);
   diff = abs(findLen - repLen);

   /*---------------------------------------------------------------
   Actions taken depend on whether the replacement string is longer,
   shorter, or equal to the length of the find string.
   ---------------------------------------------------------------*/
   count = 0;
   //equal size, can make a one-for-one char replacement
   if(findLen == repLen)
   {
      while((tmp = strstr(pStr, pFind)) != NULL)
      {
         for(i = 0; i < findLen; i++)
         {
            *tmp = pRep[i];
            tmp++;
         }/* end for() */
         count++;
      }/* end while() */
   }/* end if() */
   //the replacement is smaller, must contract the string
   else if(findLen > repLen)
   {
      while((tmp = strstr(pStr, pFind)) != NULL)
      {
         for(i = 0; i < repLen; i++)
         {
            *tmp = pRep[i];
            tmp++;
         }/* end for() */
         while(*tmp != (char)NULL)
         {
            *tmp = *(tmp+diff);
            if(*tmp == (char)NULL){ break;}
            tmp++;
         }/* end for() */
         count++;
      }/* end while() */
   }
   //the replacement is bigger, must expand the string
   else// (findLen < repLen)
   {
      while((tmp = strstr(pStr, pFind)) != NULL)
      {
         end   = (int)strlen(pStr);
         start = (end - 1) - (int)strlen(tmp) + findLen;
         for(i = end; i > start; i--)
         {
            pStr[i+diff] = pStr[i];
         }
         for(i = 0; i < repLen; i++)
         {
            *tmp = pRep[i];
            tmp++;
         }/* end for() */
         count++;
      }/* end while() */
   }/* end else() */

   return count;
}/* end MyStrRep() */

/******************************************************************************
MyStrOccur()

Counts and returns the number of occurrences of pFind in the provided string 
(pStr). 
******************************************************************************/
int MyStrOccur(UnmoveableString pStr, IroncladString pFind)
{
   int count;
   char * pTmp;

   count = 0;
   pTmp = pStr;

   while((pTmp = strstr(pTmp, pFind)) != NULL)
   { 
      count++;
      pTmp += strlen(pFind);
   }

   return count;
}/* end MyStrOccur() */

/******************************************************************************
GetInFileName()

Retrieves a pointer to a global string that is the name of the main 
configuration file. This is a temporary file that does not include the 
surrogate model sections.
******************************************************************************/
UnmoveableString GetInFileName(void)
{
#ifdef ISOFIT_BUILD
   strcpy(gInFileName, GetOstFileName());
   return gInFileName;
#else
   FILE * pNew, * pOld;
   char * pStr;
   static bool bFirstTime = true;

   if(bFirstTime == false) return gInFileName;

   bFirstTime = false;

   pOld = fopen(GetOstFileName(), "r");
   if(pOld == NULL)
   {
      FileOpenFailure("GetInFileName()", GetOstFileName());
   }/* end if() */

   pStr = MyTempName(gInFileName);
   pNew = fopen(pStr, "w");
   if(pNew == NULL)
   {
      FileOpenFailure("GetInFileName()", pStr);
   }/* end if() */
   
   while(1) //search for beginning of surrogates section
   {  
      if(feof(pOld) != 0) //section doesn't exist
      {
         goto exit_function;
      }/* end if() */
      fgets(gLine, gLineSize, pOld);

      if(gLine[0] != '#') //not a comment line
      {
         if(strstr(gLine, "BeginSurrogateModels") == NULL) // not beginning
         {
            fputs(gLine, pNew);
         }
         else
         {
            break; //exit while loop
         }
      }/* end if() */
   }/* end while() */

   while(1) //search for end of surrogates section
   {  
      if(feof(pOld) != 0)
      {
         fclose(pNew);
         fclose(pOld);
         MissingTokenFailure("EndSurrogateModels", GetOstFileName());
         return NULL;
      }/* end if() */
      fgets(gLine, gLineSize, pOld);

      if(gLine[0] != '#') //not a comment line
      {
         if(strstr(gLine, "EndSurrogateModels") != NULL) // not beginning
         {
            break;
         }
      }/* end if() */
   }/* end while() */

   //copy remaining data lines
   while(1)
   {  
      fgets(gLine, gLineSize, pOld);

      if(feof(pOld) != 0) break;

      if(gLine[0] != '#') //not a comment line
      {
         fputs(gLine, pNew);
      }      
   }/* end while() */

exit_function :
   fclose(pNew);
   fclose(pOld);

   return gInFileName;
#endif
}/* end GetInFileName() */

/******************************************************************************
GetSrgFileName()

Retrieves a pointer to a global string that is the name of the surrogate models
configuration file. This is a temporary file.
******************************************************************************/
UnmoveableString GetSrgFileName(void)
{
   FILE * pNew, * pOld;
   char * pStr;
   static bool bFirstTime = true;

   if(bFirstTime == false) return gSrgFileName;

   bFirstTime = false;

   pOld = fopen(GetOstFileName(), "r");
   if(pOld == NULL)
   {
      FileOpenFailure("GetSrgFileName()", GetOstFileName());
   }/* end if() */

   pStr = MyTempName(gSrgFileName);
   pNew = fopen(pStr, "w");
   if(pNew == NULL)
   {
      FileOpenFailure("GetSrgFileName()", pStr);
   }/* end if() */
   
   while(1) //search for beginning of surrogates section
   {  
      if(feof(pOld) != 0) //section doesn't exist
      {
         fclose(pNew);
         fclose(pOld);
         MissingTokenFailure("BeginSurrogateModels", GetOstFileName());
         return NULL;
      }/* end if() */
      fgets(gLine, gLineSize, pOld);

      if(gLine[0] != '#') //not a comment line
      {
         if(strstr(gLine, "BeginSurrogateModels") != NULL) // not beginning
         {
            fputs(gLine, pNew);
            break;
         }
      }/* end if() */
   }/* end while() */

   while(1) //search for end of surrogates section
   {  
      if(feof(pOld) != 0)
      {
         fclose(pNew);
         fclose(pOld);
         MissingTokenFailure("EndSurrogateModels", GetInFileName());
         return NULL;
      }/* end if() */
      fgets(gLine, gLineSize, pOld);

      if(gLine[0] != '#') //not a comment line
      {
         fputs(gLine, pNew);
     
         if(strstr(gLine, "EndSurrogateModels") != NULL) //found end
         {
            break;
         }
      }/* end if() */
   }/* end while() */

   fclose(pNew);
   fclose(pOld);

   return gSrgFileName;
}/* end GetSrgFileName() */

/******************************************************************************
GetDynFileName()

Retrieves a pointer to a global string that is the name of one of the 
surrogate models configuration file. This is a temporary file.
******************************************************************************/
UnmoveableString GetDynFileName(IroncladString pTok)
{
   static bool first_time = true;
   char beg_tok[DEF_STR_SZ], end_tok[DEF_STR_SZ];
   FILE * pNew, * pOld;
   char * pStr;

   if(first_time == false){
      remove(gDynFileName);}
   first_time = false;
   if(pTok == NULL) return NULL;

   //generate dynamic tokens and filename
   sprintf(beg_tok, "Begin_%s_Model", pTok);
   sprintf(end_tok, "End_%s_Model", pTok);

   pOld = fopen(GetSrgFileName(), "r");
   if(pOld == NULL)
   {
      FileOpenFailure("GetDynFileName()", GetSrgFileName());
   }/* end if() */

   pStr = MyTempName(gDynFileName);
   pNew = fopen(pStr, "w");
   if(pNew == NULL)
   {
      FileOpenFailure("GetDynFileName()", pStr);
   }/* end if() */
   
   while(1) //search for beginning of surrogate section
   {  
      if(feof(pOld) != 0) //section doesn't exist
      {
         fclose(pNew);
         fclose(pOld);
         MissingTokenFailure(beg_tok, GetSrgFileName());
         return NULL;
      }/* end if() */
      fgets(gLine, gLineSize, pOld);

      if(gLine[0] != '#') //not a comment line
      {
         if(strstr(gLine, beg_tok) != NULL) // beginning
         {
            fputs(gLine, pNew);
            break;
         }
      }/* end if() */
   }/* end while() */

   while(1) //search for end of surrogate section
   {  
      if(feof(pOld) != 0)
      {
         fclose(pNew);
         fclose(pOld);
         MissingTokenFailure(end_tok, GetSrgFileName());
         return NULL;
      }/* end if() */
      fgets(gLine, gLineSize, pOld);

      if(gLine[0] != '#') //not a comment line
      {
         fputs(gLine, pNew);
     
         if(strstr(gLine, end_tok) != NULL) //found end
         {
            break;
         }
      }/* end if() */
   }/* end while() */

   fclose(pNew);
   fclose(pOld);

   return gDynFileName;
}/* end GetDynFileName() */

/******************************************************************************
FindToken()

Locates the token in the file, file pointer will be positioned at the 
following line.
******************************************************************************/
void FindToken(FILE * pFile, IroncladString token, IroncladString pName)
{
   gLine[0] = (char)NULL;

   do //search for token
   {  
      if(feof(pFile) != 0)
      {
         fclose(pFile);
         MissingTokenFailure(token, pName);
      }/* end if() */
      fgets(gLine, gLineSize, pFile);
   }while((gLine[0] == '#') || (strstr(gLine, token) == NULL));
} /* end FindTOken() */

/******************************************************************************
CheckToken()

Checks to see if the token exists in the file, file pointer will be 
positioned at the following line if token is found, else file will
be rewound.
******************************************************************************/
bool CheckToken(FILE * pFile, IroncladString token, IroncladString pName)
{
   gLine[0] = (char)NULL;
   do //search for token
   {  
      if(feof(pFile) != 0)
      {
         rewind(pFile);
         return(false);
      }/* end if() */
      fgets(gLine, gLineSize, pFile);
   }while((gLine[0] == '#') || (strstr(gLine, token) == NULL));

   return(true);
}/* end CheckToken() */

/******************************************************************************
InitDataLine()

Initialize the data line. Passing a NULL argument will free up the memory.
******************************************************************************/
void InitDataLine(IroncladString pName)
{
   int maxLineSize = 0;
   int lineSize = 0;
   char c;

   if(pName == NULL)
   {
      delete [] gLine;
      gLineSize = 0;
      return;
   }

   FILE * pFile = fopen(pName, "r");
   
   if(pFile == NULL)
   {
      printf("InitDataLine() : Couldn't open file |%s|\n", pName);
      return;
   }/* end if() */

   while(feof(pFile) == 0) 
   {
      c = (char)(fgetc(pFile));
      lineSize++;
      if(c == '\n')
      {
         if(lineSize > maxLineSize){
            maxLineSize = lineSize;
         }
         lineSize = 0;
      }
   }/* end while() */

   fclose(pFile);

   //last line might not end in a carriage return
   if(lineSize > maxLineSize){
      maxLineSize = lineSize;
   }
   
   maxLineSize *= 2;
   if(maxLineSize > gLineSize){
      delete [] gLine;
      gLineSize = maxLineSize;
      NEW_PRINT("char *", gLineSize);
      gLine = new char[gLineSize];
      MEM_CHECK(gLine);
   }
}/* end InitDataLine() */

/******************************************************************************
GetCurDataLine()

Retrieve current line of input (whatever is stored in gLine).
******************************************************************************/
char * GetCurDataLine(void)
{
   return gLine;
} /* end GetCurDataLine() */

/******************************************************************************
GetNxtDataLine()

Retrieve next line of input, skipping over comments and blank lines.
******************************************************************************/
char * GetNxtDataLine(FILE * pFile, IroncladString pName)
{
   gLine[0] = (char)NULL;

   do //skip over comments (lines whose first char is '#')
   {  
      if(feof(pFile) != 0)
      {
         fclose(pFile);
         EndOfFileFailure("GetNxtDataLine", pName);
      }/* end if() */
      fgets(gLine, gLineSize, pFile);
      MyTrim(gLine);
   }while((gLine[0] == '#') || (gLine[0] == (char)NULL));

  return gLine;
}/* end GetNxtDataLine() */

/******************************************************************************
IsWhitespace()

Returns true if char is whitespace, false otherwise.
******************************************************************************/
bool IsWhitespace(char c)
{
   if(c == ' ')  return true;
   if(c == '\t') return true;
   if(c == '\n') return true;
   if(c == '\r') return true;
   return false;
}/* end IsWhitespace() */

/******************************************************************************
ExtractFileName()

Parse a filename from pLine. Filenames can contain spaces in the middle and
so are terminated by TAB or ';' characters.

Returns the index of the beginning of the next string in pLine.
******************************************************************************/
int ExtractFileName(IroncladString pLine, UnmoveableString pOut)
{
   int i,j;

   i = 0;  
   j = 0;
   while((IsWhitespace(pLine[i]) == true) || (pLine[i] == ';')) { i++;} 
   while(((IsWhitespace(pLine[i]) == false) || (pLine[i] == ' ')) && (pLine[i] != ';'))
   { 
      if(pLine[i] == (char)NULL){break;}

      pOut[j] = pLine[i];
      i++;
      j++;
   }   
   pOut[j] = (char)NULL;
   //For Linux compatibility, must trim off any trailing whitespace
   j--;
   while(pOut[j] == ' '){pOut[j] = (char)NULL; j--;}

   if(pLine[i] == ';'){i++;}//skip past semi-colon to prevent future parse errors
   return i;
}/* end ExtractFileName() */

/******************************************************************************
ExtractString()

Parse a string from pLine. Strings are separated/terminated by whitespace
characters.

Returns the index of the beginning of the next string in pLine.

Returns -1 if no whitespace tokens could be find.
******************************************************************************/
int ExtractString(IroncladString pLine, UnmoveableString pOut)
{
   int i, j;

   i = 0;  
   j = 0;
   //skip any leading whitespace...
   while((pLine[i] == ' ') || (pLine[i] == '\t')) { i++;}
   //read until more whitespace is detected....
   while((pLine[i] != ' ') && (pLine[i] != '\t')) 
   { 
      if(pLine[i] == (char)NULL){pOut[j] = (char)NULL; return -1;} //end of string
      if(pLine[i] == '\n'){pOut[j] = (char)NULL; return -1;} //end of line
      if(pLine[i] == '\r'){pOut[j] = (char)NULL; return -1;} //end of line

      pOut[j] = pLine[i];
      i++;
      j++;
   }
   pOut[j] = (char)NULL;
   return i;
}/* end ExtractString() */

/******************************************************************************
ExtractColString()

Parse a column string from pLine. Column strings are separated/terminated by 
the 'tok' character.

Returns the index of the beginning of the next string in pLine.

Returns -1 if the specified token could not be find.
******************************************************************************/
int ExtractColString(IroncladString pLine, UnmoveableString pOut, char tok)
{
   int i, j;

   i = 0;  
   j = 0;
   //read until token is detected....
   while(pLine[i] != tok)
   { 
      if(pLine[i] == (char)NULL){pOut[j] = (char)NULL; return -1;} //end of string
      if(pLine[i] == '\n'){pOut[j] = (char)NULL; return -1;} //end of line
      if(pLine[i] == '\r'){pOut[j] = (char)NULL; return -1;} //end of line

      pOut[j] = pLine[i];
      i++;
      j++;
   }
   pOut[j] = (char)NULL;
   
   return (i+1); //inc. by 1 to skip over token
}/* end ExtractColString() */

/******************************************************************************
ValidateExtraction()

Tests the results of an ExtractColString() or ExtractString() call. Arguments 
are as follows:

j    - result of previous call to ExtractColString() or ExtractString()
cur  - the column entry that was expected to be extracted
last - the last possible column entry (the loop termination index)
pF   - the name of the calling routine (for inclusing in an error message)

Returns j or 0 if extraction is ok. Exits program gracefully otherwise.
******************************************************************************/
int ValidateExtraction(int j, int cur, int last, char * pF)
{
   if(j != -1) return j;
   if(cur == (last-1)) return 0; //last item in sequence, so ok
   char msg[DEF_STR_SZ];
   sprintf(msg, "%s : Unexpected end of input", pF);
   LogError(ERR_FILE_IO, msg);
   ExitProgram(1);
   return 0;
}/* end ValidateExtraction() */


/******************************************************************************
CheckExtraction()

Tests the results of an ExtractColString() or ExtractString() call. Arguments 
are as follows:

j    - result of previous call to ExtractColString() or ExtractString()
cur  - the column entry that was expected to be extracted
last - the last possible column entry (the loop termination index)
pF   - the name of the calling routine (for inclusing in an error message)

Returns j or 0 if extraction is ok. Logs and error and returns -1 otherwise.
******************************************************************************/
int CheckExtraction(int j, int cur, int last, char * pF)
{
   if(j != -1) return j;
   if(cur == (last-1)) return 0; //last item in sequence, so ok
   char msg[DEF_STR_SZ];
   sprintf(msg, "%s : Unexpected end of input", pF);
   LogError(ERR_FILE_IO, msg);
   return -1;
}/* end CheckExtraction() */

/******************************************************************************
SortInc()

Sorts a list of numbers in increasing order using quick sort algorithm.
******************************************************************************/
void SortInc(Unmoveable1DArray v, int size)
{
   qsort((void *)v, (size_t)size, sizeof(double), CompDbl);
} /* end SortInc() */

/******************************************************************************
MatMult()

Multiplies two matrices and stores result in mOut.
******************************************************************************/
void MatMult(Ironclad2DArray m1, Ironclad2DArray m2, Unmoveable2DArray mOut, 
             int row1, int row2, int col2)
{
   int i;
   int j;
   int k;
   double sum;

   for(i = 0; i < row1; i++)
   {
      for(j = 0; j < col2; j++)
      {
         sum = 0.00;
         for(k = 0; k < row2; k++) 
         {
            sum += (m1[i][k] * m2[k][j]);
         }/* end for() */
         mOut[i][j] = sum;
      }/* end for() */
   }/* end for() */
}/* end MatMult() */

/******************************************************************************
VectMult()

Multiples a vector and a matrix and stores resulting vector in vOut.
******************************************************************************/
void VectMult(Ironclad2DArray m, Ironclad1DArray v, Unmoveable1DArray vOut, 
              int rows, int cols)
{
   double sum;
   int i;
   int j;
   
   for(i = 0; i < rows; i++)
   {
      sum = 0.00;
      for(j = 0; j < cols; j++)
      {
         sum += (m[i][j] * v[j]);
      }/* end for() */
      vOut[i] = sum;
   }/* end for() */
}/* end VectMult() */

/******************************************************************************
MatInv()

Inverts a matrix and stores result in inv without altering m.

Returns true if successful, false otehrwise.
******************************************************************************/
bool MatInv(Ironclad2DArray m, Unmoveable2DArray inv, int size)
{
   int i;
   int col;
   int row;
   int pivRow;
   int pivCol;
   double max;
   double val;
   double sf;
   //declared static to avoid memory fragmentation
   static double ** copy = NULL; 
   static int lastSize = 0;

   //empty input args signifies deletion of copy matrix
   if((m == NULL) && (inv == NULL) && (size == 0))
   {
      for(i = 0; i < lastSize; i++){ delete [] copy[i];}
      delete [] copy;
      copy = NULL;
      lastSize = 0;
      return true;
   }

   /*-------------------------------------------------
   Copy input matrix and work off of the copy, because 
   the inversion algorithm will destroy the input matrix.
   ---------------------------------------------------*/
   if(size > lastSize) //resize copy, if necessary
   {
      for(i = 0; i < lastSize; i++){ delete [] copy[i];}
      delete [] copy;

      NEW_PRINT("double *", size);
      copy = new double *[size];
      MEM_CHECK(copy);

      for(row = 0; row < size; row++)
      {
         NEW_PRINT("double", size);
         copy[row] = new double [size];
         MEM_CHECK(copy[row]);
      }   
      lastSize = size;
   }/* end if() */

   //size only input arg signifies creation of copy matrix only
   if((m == NULL) && (inv == NULL) && (size != 0))
   {
      return true;
   }

   for(row = 0; row < size; row++)
   {
      for(col = 0; col < size; col++)
      {
         copy[row][col] = m[row][col];
      }
   }   

   //init. inv. matrix to be identity 
   for(row = 0; row < size; row++){
    for(col = 0; col < size; col++){
      inv[row][col] = 0.00;}}
   for(row = 0; row < size; row++){
      inv[row][row] = 1.00;}

   //use gauss-jordan with partial pivoting   
   for(i = 0; i < size; i++)
   {
      //determine pivot row
      max = fabs(copy[i][i]);
      pivRow = i;
      for(row = i; row < size; row++)
      {
         val = fabs(copy[row][i]);
         if(val > max)
         {
            max = val;
            pivRow = row;
         }/* end if() */
      }/* end for() */

      //check for singular matrix
      if(max <= NEARLY_ZERO) 
      { 
         LogError(ERR_SING_MAT, "MatInv(): pivot too small");
         return false;
      }/* end if() */  

      //perform pivot
      if(pivRow != i)
      {
         //swap rows of both matrices
         for(pivCol = i; pivCol < size; pivCol++)
         {
            val = copy[i][pivCol];
            copy[i][pivCol] = copy[pivRow][pivCol];
            copy[pivRow][pivCol] = val;
         }/* end for() */
         for(pivCol = 0; pivCol < size; pivCol++)
         {
            val = inv[i][pivCol];
            inv[i][pivCol] = inv[pivRow][pivCol];
            inv[pivRow][pivCol] = val;
         }/* end for() */
      }/* end if() */    

      //perform gaussian elimination
      for(row = 0; row < size; row++)
      {
         if(row != i)
         {
            sf = (copy[row][i] / copy[i][i]);
            for(col = i; col < size; col++){copy[row][col] -= (sf*copy[i][col]);}
            for(col = 0; col < size; col++){inv[row][col]-=(sf*inv[i][col]);}
         }/* end if() */
      }/* end for() */
   }/* end for() */

   //divide through by diagonals
   for(row = 0; row < size; row++)
   {
      sf = copy[row][row];
      for(col = 0; col < size; col++)
      {  
         copy[row][col]   /= sf;
         inv[row][col] /= sf;
      }/* end for() */
   }/* end for() */

   return true;
}/* end MatInv() */

/******************************************************************************
CholeskyDecomp()

Decomposes matrix (A) into the product of a matrix (L) and it's transpose (LT),
using Cholesky decomposition method. If matrix A is not symmetric, the 
decomposition proceeds, but symmetry is faked by using the average of Aij and Aji 
in place of non-symmetrical values. This allows the routine to be used in 
situations where only approximate symmetry is required. Likewise, if matrix A is
not positive definite, this too will be faked by using absolute values.

Returns 
    0  if A is symmetric and positive definite
    1  if A is positive definite, but non-symmetric
    2  if A is not positive definite, but is symmetric
    3  if A is not positive definite and non-symmetric
******************************************************************************/
int CholeskyDecomp(Ironclad2DArray A, Unmoveable2DArray L, Unmoveable2DArray LT, int size)
{
   bool isSym, isPos;
   int errCode;
   int i, j, k;
   double Aij, Aji, sum;

   //check symmetry
   isSym = true;
   for(i = 0; i < size; i++)
   {
      for(j = 0; j < size; j++)
      {
         Aij = A[i][j];
         Aji = A[j][i];
         if(fabs(Aij - Aji) > NEARLY_ZERO)
         {
            isSym = false;
         }
      }/* end for() */
   }/* end for() */

   //Setup initial L matrix
   for(i = 0; i < size; i++)
   {
      for(j = 0; j < size; j++)
      {         
         Aij = A[i][j];
         Aji = A[j][i];
         L[i][j] = 0.5*(Aij + Aji);
      }/* end for() */
   }/* end for() */
   
   //perform decomposition
   isPos = true;
   for(i = 0; i < size; i++)
   {
      for(j = i; j < size; j++)
      {
         for(sum = L[i][j], k = i - 1; k >= 0; k--)
         {
            sum -= L[i][k]*L[j][k];
         }/* end for() */
         if(i == j)
         {
            if(sum <= 0.00)
            { 
               isPos = false; 
               sum = fabs(sum)+NEARLY_ZERO;
            }/* end if() */
            LT[j][i] = sqrt(sum);               
         }/* end if() */
         else
         {
            L[j][i] = sum / LT[i][i];
         }/* end else() */
      }/* end for() */
   }/* end for() */

   //copy diagonal from LT to L
   for(i = 0; i < size; i++)
   {
      L[i][i] = LT[i][i];
   }

   //zero the upper triangle
   for(i = 0; i < size; i++)
   {
      for(j = i; j < size; j++)
      {
         if(i != j){ L[i][j] = 0.00;}
      }
   }

   //fill transpose
   for(i = 0; i < size; i++)
   {
      for(j = 0; j < size; j++)
      {
         LT[i][j] = L[j][i];
      }
   }
   
   //set error code and return
   errCode = 0;
   if(isSym == false){ errCode += 1;}
   if(isPos == false){ errCode += 2;}
   return errCode;
}/* end CholeskyDecomp() */

/******************************************************************************
CalcDeterminant()

Computes the determinant of a matrix. Matrix m must be a square matrix.
******************************************************************************/
double CalcDeterminant(Ironclad2DArray m, int size)
{
   bool isOK;
   double ** pScratch; //scratch  matrix
   double det;
   int i, j, *indx;

   if(size == 1)
   {
      det = m[0][0];
   }
   else if(size == 2)
   {
      det = (m[0][0]*m[1][1]) - (m[0][1]*m[1][0]);
   }/* end if() */
   else//utilize LU decomposition, then multiply diagonals
   {
      //allocate space for scratch matrix
      NEW_PRINT("double *", size);
      pScratch = new double*[size];
      MEM_CHECK(pScratch);

      for(i = 0; i < size; i++)
      { 
         NEW_PRINT("double", size);
         pScratch[i] = new double[size];
      }
      MEM_CHECK(pScratch[i-1]);

      NEW_PRINT("int", size);
      indx = new int[size];
      MEM_CHECK(indx);

      for(i = 0; i < size; i++){
         for(j = 0; j < size; j++){
            pScratch[i][j] = m[i][j];}}

      isOK = ludcmp(pScratch, size, indx, &det);
      if(isOK == true)
      {
         for(i = 0; i < size; i++) { det *= pScratch[i][i];}
      }
      else
      {
         det = 0.00;
      }

      for(i = 0; i < size; i++){ delete [] pScratch[i];}
      delete [] pScratch;
      delete [] indx;
   }/* end else() */

   return det;
}/* end CalcDeterminant() */

/******************************************************************************
DotProduct()

Compute the dot product of two vectors.
******************************************************************************/
double DotProduct(Ironclad1DArray v1, Ironclad1DArray v2, int size)
{
   int i;
   double sum;

   sum = 0.00;
   for(i = 0; i < size; i++){ sum += v1[i]*v2[i];}
   return sum;
}/* end DotProduct() */

/******************************************************************************
CheckOverflow()

Checks to see if a number is too big. If overflow occurs, returns true.
******************************************************************************/
bool CheckOverflow(double num)
{
   if((num == 0.00) || (num == 1.00)){return false;}
   if((num*num) == num) {return true;}
   if(num != num) {return true;}
   return false;
}/* end CheckOverflow() */

/******************************************************************************
MyTime()

Get the time elapsed (in seconds).
******************************************************************************/
unsigned int MyTime(void)
{
#ifdef USE_MY_TIME
   return((unsigned)time(NULL));
#else
   return((unsigned)time(NULL));
#endif
}/* end MyTime() */

/******************************************************************************
GetElapsedTime()

Get the time elapsed (in seconds) from the start of the program.
******************************************************************************/
unsigned int GetElapsedTime(void)
{
   static unsigned int start = MyTime();
   return (MyTime() - start);
}/* end GetElapsedTime() */

/******************************************************************************
MyRand()

Generates a 32-bit random number by successive calls to rand(), which may be 
limited to a maximum value less than 2^32-1
******************************************************************************/
unsigned int MyRand(void)
{
   static bool set_seed = true;
   unsigned int r;
   unsigned int t1, t2, t3, t4;

   if(set_seed == true)
   {
      set_seed = false;
	   srand(GetRandomSeed());
   }
   if(RAND_MAX >= MY_RAND_MAX) return rand();

   t1 = (rand() & 0x0000007F) << 24;
   t2 = (rand() & 0x000000FF) << 16;
   t3 = (rand() & 0x000000FF) << 8;
   t4 = (rand() & 0x000000FF);
   r = (t1 | t2 | t3 | t4);
   return r;   
}/* end MyRand() */

/**********************************************************************
MyGaussRand()

Returns a random number sampled from a normal distribution with 
mean of m and a standard deviation of s.
**********************************************************************/
double MyGaussRand(double m, double s)
{
   double p = UniformRandom(); //probability
   double r = StdNormInvCDF(p); //random number of std. dev. from mean
   return (m + r*s);
}/* end MyGaussRand() */

/**********************************************************************
		UniformRandom
-----------------------------------------------------------------------
	return uniformly distributed random number between 0 and 1
	coded by James Craig
**********************************************************************/
double UniformRandom()
{
	return (double)(MyRand()) / (double)(MY_RAND_MAX);
}
/**********************************************************************
		GaussRandom
-----------------------------------------------------------------------
	returns a standard Gaussian random number 
	based upon Numerical recipes (gasdev) and Marsagalia-Bray Algorithm
	coded by James Craig
**********************************************************************/
double GaussRandom()
{
	double  ranval,zvalue;
	double  Work3,Work2,Work1;

	Work3=2.0; 
	while ((Work3>=1.0) || (Work3==0.0))
	{
		ranval=UniformRandom(); 
		Work1 = 2.0 * ranval - 1.0;
		ranval=UniformRandom(); 
		Work2 = 2.0 * ranval - 1.0;
		Work3 = Work1 * Work1 + Work2 * Work2;
	}
	Work3 = pow((-2.0*log(Work3))/Work3,0.5);  // natural log
	// pick one of two deviates at random: (don't worry about trying to use both)
	ranval=UniformRandom(); 

	if (ranval<0.5) {zvalue = Work1 * Work3;}
	else            {zvalue = Work2 * Work3;}

	return zvalue;
}
/******************************************************************************
MyTempName()

Generates a temporary file name
******************************************************************************/
char * MyTempName(UnmoveableString pStr)
{
   static int count = 0;
   int id;

   MPI_Comm_rank(MPI_COMM_WORLD, &id); 

   sprintf(pStr, "OstTemp_%02d_%02d.txt", id, count++);

   return pStr;
}/* end MyTempName() */

double MyMin(double a, double b)
{
   return (a < b ? a : b);
}
double MyMax(double a, double b)
{
   return (a > b ? a : b);
}
int		 iMax(const int a,const int b)      
{
	return (a > b ? a : b);
}

/******************************************************************************
GetPreciseValAsStr()

Write a number (x) to the specified string (valStr) using the requested number
of digits of precision.
******************************************************************************/
void GetPreciseValAsStr(UnmoveableString pOut, double x)
{
   int precision = GetNumDigitsOfPrecision();

   switch(precision)
   {
      case (1) : sprintf(pOut, "%.1E", x); break;
      case (2) : sprintf(pOut, "%.2E", x); break;
      case (3) : sprintf(pOut, "%.3E", x); break;
      case (4) : sprintf(pOut, "%.4E", x); break;
      case (5) : sprintf(pOut, "%.5E", x); break;
      case (6) : sprintf(pOut, "%.6E", x); break;
      case (7) : sprintf(pOut, "%.7E", x); break;
      case (8) : sprintf(pOut, "%.8E", x); break;
      case (9) : sprintf(pOut, "%.9E", x); break;
      case (10) : sprintf(pOut, "%.10E", x); break;
      case (11) : sprintf(pOut, "%.11E", x); break;
      case (12) : sprintf(pOut, "%.12E", x); break;
      case (13) : sprintf(pOut, "%.13E", x); break;
      case (14) : sprintf(pOut, "%.14E", x); break;
      case (15) : sprintf(pOut, "%.15E", x); break;
      case (16) : sprintf(pOut, "%.16E", x); break;
      case (17) : sprintf(pOut, "%.17E", x); break;
      case (18) : sprintf(pOut, "%.18E", x); break;
      case (19) : sprintf(pOut, "%.19E", x); break;
      case (20) : sprintf(pOut, "%.20E", x); break;
      case (21) : sprintf(pOut, "%.21E", x); break;
      case (22) : sprintf(pOut, "%.22E", x); break;
      case (23) : sprintf(pOut, "%.23E", x); break;
      case (24) : sprintf(pOut, "%.24E", x); break;
      case (25) : sprintf(pOut, "%.25E", x); break;
      case (26) : sprintf(pOut, "%.26E", x); break;
      case (27) : sprintf(pOut, "%.27E", x); break;
      case (28) : sprintf(pOut, "%.28E", x); break;
      case (29) : sprintf(pOut, "%.29E", x); break;
      case (30) : sprintf(pOut, "%.30E", x); break;
      case (31) : sprintf(pOut, "%.31E", x); break;
      case (32) : sprintf(pOut, "%.32E", x); break;
      default  : sprintf(pOut, "%.6E", x); break;        
   }
}/* end GetPreciseValAsStr() */

/******************************************************************************
SampleWithReplacement()

Sample from a range of integer values with replcament.
******************************************************************************/
int SampleWithReplacement(int opFlag, int range)
{
   int val, r;
   static int * pDeck = NULL; //the 'deck' that is to be sampled
   static int size = 0; //size of the deck

   if(opFlag == -1) //initialize
   {
      delete [] pDeck;
      size = range;
      pDeck = new int[range];
      for(int i = 0; i < range; i++) pDeck[i] = i;
      return 0;
   }
   else if(opFlag == -2) //re-initialize
   {
      size = range;
      for(int i = 0; i < range; i++) pDeck[i] = i;
      return 0;
   }
   else if(opFlag == -3) //free up memory
   {
      delete [] pDeck;
      pDeck = NULL;
      size = 0;
      return 0;
   }
   else //sample with replacement
   {
      r = MyRand()%size;
      val = pDeck[r];
      size--;
      pDeck[r] = pDeck[size];
      pDeck[size] = val;

      //FILE * pWgt = fopen("SampleWithReplacement.txt", "a");
      //fprintf(pWgt, "Sampled Value = %d\n", val);
      //fclose(pWgt);
      
      return val;
   }
}/* end SampleWithReplacement() */

/******************************************************************************
ExecuteCommandLine()

Execute the the given command line.
******************************************************************************/
void ExecuteCommandLine(IroncladString cmd, bool isRead, IroncladString fileName, IroncladString paramName)
{
   char psBuffer[DEF_STR_SZ];
	char output[DEF_STR_SZ];
   FILE * pPipe;

   #ifdef WIN32
      pPipe = _popen(cmd, "rt");
	   if (isRead)
	   {
		   while(!feof(pPipe))
		   {
			   if(fgets(psBuffer, DEF_STR_SZ, pPipe) != NULL)
			   {
				   strcpy(output, paramName);
				   strcat(output, " ");
				   strcat(output, psBuffer);
				   FILE * pFile;
				   pFile = fopen(fileName, "a");
				   fprintf(pFile, output);
				   fclose(pFile);
			   }
		   }
	   }
	   _pclose(pPipe);
   #else
      system(cmd);
   #endif

}/* end ExecuteCommandLine() */
