/******************************************************************************
File     : MyPBS.h
Author   : L. Shawn Matott
Copyright: 2004, L. Shawn Matott

This module contains PBS commands which can be executed from within the Ostrich
program. This will allow Ostrich to run parallel versions of the model executable.

Version History
02-09-04    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef MY_PBS_H
#define MY_PBS_H

extern "C" {
   bool TestMyPBS(void);
   char * MyPBS_Qsub(char * args);
   char * MyPBS_Qstat(char * args);
   char * MyPBS_Qdel(char * args);
   char * MyPBS_Qselect(char * args);
   char * MyPBS_Qrerun(char * args);
   char * MyPBS_Qorder(char * args);
   char * MyPBS_Qmove(char * args);
   char * MyPBS_Qhold(char * args);
   char * MyPBS_Qalter(char * args);
   char * MyPBS_Qmsg(char * args);
   char * MyPBS_Qrls(char * args);
   char * MyPBS_Showbf(char * args);
   char * MyPBS_Sleep(char * args);
}

#endif /* MY_PBS_H */


