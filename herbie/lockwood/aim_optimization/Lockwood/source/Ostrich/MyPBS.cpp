/******************************************************************************
File     : MyPBS.cpp
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

This module contains PBS commands which can be executed from within the Ostrich
program. This will allow Ostrich to run parallel versions of the model executable.

Version History
02-09-04    lsm   added copyright information and initial comments.
08-17-04    lsm   added reporting of memory allocations
******************************************************************************/
#include "MyPBS.h"
#include "Exception.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char gPbsCmdLine[DEF_STR_SZ];
const char * gPbsOutFile = "pbs_result.txt";
const char * gPbsErrFile = "pbs_error.txt";
char * gPbsCmdResult = NULL;
bool gPbsCmdError = false;

void Get_PBS_Result(FILE * pFile);
bool Check_PBS_Error_File(void);
char * Run_PBS_Command(char * cmd, char * args);

//Verifies that all PBS functions are supported by the shell.
bool TestMyPBS(void)
{
   char * myStr;
   myStr = MyPBS_Qstat("");
   printf("qstat result = %s", myStr);
   return (!gPbsCmdError);
}/* end TestMyPBS() */

char * Run_PBS_Command(char * cmd, char * args)
{
   FILE * pFile;

   strcpy(gPbsCmdLine, cmd);   
   strcat(gPbsCmdLine, " ");
   strcat(gPbsCmdLine, args);
   strcat(gPbsCmdLine, " 2> ");
   strcat(gPbsCmdLine, gPbsErrFile);
   strcat(gPbsCmdLine, " 1> ");
   strcat(gPbsCmdLine, gPbsOutFile);

   printf("%s\n", gPbsCmdLine);
   system(gPbsCmdLine);

   //check to see if error file contains anything
   if(Check_PBS_Error_File() == true)
   {
      pFile = fopen(gPbsErrFile, "r");
      gPbsCmdError = true;
   }
   else
   {
      pFile = fopen(gPbsOutFile, "r");
      gPbsCmdError = false;
   }
   if(pFile == NULL)
   {
      printf("Run_PBS_Command(): couldn't retrieve result\n");
      gPbsCmdError = true;
      return(NULL);
   }
   else
   {
      Get_PBS_Result(pFile);
      fclose(pFile);      
   }
   return gPbsCmdResult;
}/* end Run_PBS_Command() */

char * MyPBS_Qsub(char * args)
{
   return(Run_PBS_Command("qsub", args));
}/* end MyQsub() */

char * MyPBS_Qstat(char * args)
{
   return(Run_PBS_Command("qstat", args));
}
char * MyPBS_Qdel(char * args)
{
   return(Run_PBS_Command("qdel", args));
}
char * MyPBS_Qselect(char * args)
{
   return(Run_PBS_Command("qselect", args));
}
char * MyPBS_Qrerun(char * args)
{
   return(Run_PBS_Command("qrerun", args));
}
char * MyPBS_Qorder(char * args)
{
   return(Run_PBS_Command("qorder", args));
}
char * MyPBS_Qmove(char * args)
{
   return(Run_PBS_Command("qmove", args));
}
char * MyPBS_Qhold(char * args)
{
   return(Run_PBS_Command("qhold", args));
}
char * MyPBS_Qalter(char * args)
{
   return(Run_PBS_Command("qalter", args));
}
char * MyPBS_Qmsg(char * args)
{
   return(Run_PBS_Command("qmsg", args));
}
char * MyPBS_Qrls(char * args)
{
   return(Run_PBS_Command("qrls", args));
}

char * MyPBS_Showbf(char * args)
{
   return(Run_PBS_Command("showbf", args));
}

char * MyPBS_Sleep(char * args)
{
   return(Run_PBS_Command("sleep", args));
}

/******************************************************************************
Get_PBS_Result()

Reads the file containing the result of the PBS command is into a string.
******************************************************************************/
void Get_PBS_Result(FILE * pFile)
{
   int fileSize;
   int i;

   /*
   count number of chars in file, so that the string can be sized
   */
   fileSize = 0;
   while(feof(pFile) == 0) 
   {
      fileSize++;
      fgetc(pFile);
   }/* end while() */
   fileSize--;

   //size fileStr
   if(gPbsCmdResult != NULL)
   {
      delete [] gPbsCmdResult;
   }/* end if() */
   NEW_PRINT("char", fileSize+1);
   gPbsCmdResult = new char[fileSize + 1];
   MEM_CHECK(gPbsCmdResult);

   //fill fileStr
   rewind(pFile);   
   for(i = 0; i < fileSize; i++)
   {
      gPbsCmdResult[i] = (char)(fgetc(pFile));
   }/* end for() */
   gPbsCmdResult[i] = 0;
} /* end Get_PBS_Result() */

/******************************************************************************
Check_PBS_Error_File()

Checks to see if there is anyhting in the error file. If so, return true, else
return false.
******************************************************************************/
bool Check_PBS_Error_File(void)
{
   FILE * pFile;
   int fileSize;

   pFile = fopen(gPbsErrFile, "r");
   if(pFile == NULL){ return false;}
   
   //count number of chars in file
   fileSize = 0;
   while(feof(pFile) == 0) 
   {
      fileSize++;
      fgetc(pFile);
   }/* end while() */
   fileSize--;
   fclose(pFile);

   if(fileSize > 0) { return true;}
   return false;
}/* end Check_PBS_Error_File() */

