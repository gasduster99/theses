/******************************************************************************
File      : SCEUA.cpp
Author    : L. Shawn Matott - Converted from Origincal SCEUA Fortran code
Copyright : 2009, L. Shawn Matott

The SCE-UA method is a general purpose global optimization
program - Shuffled Complex Evolution.  It was originally developed 
by Dr. Qingyun Duan as part of his doctoral dissertation work at
University of Arizona, Tucson, AZ 85721, USA. 

The dissertation is entitled "A Global Optimization Strategy for
Efficient and Effective Calibration of Hydrologic Models".  The
program has since been modified to make it easier for use on
problems of users' interests.  

The algorithm has been described in detail in an article entitled :
"Effective and Efficient Global Optimization for Conceptual Rainfall-Runoff Models", 
Water Resources Research, Vol 28(4), pp.1015-1031, 1992

And in an article entitled:
"A Shuffled Complex Evolution Approach for Effective and Efficient Global Minimization",
 by Q. Duan, V.K. Gupta and S. Sorooshian, Journal of Optimization Theory and its 
Applications, Vol 76(3), pp 501-521, 1993.  

A paper entitled "Optimal Use of the SCE-UA Global Optimization Method for Calibrating Watershed Models", 
by Q. Duan, S. Sorooshian, & V.K. Gupta, Journal of Hydrology, Vol.158, 265-284, 1994, 
discussed how to use the SCE-UA Method in an efficient and  effective manner.

Input Summary for the SCEUA algorithm (adapted from original Fortran-based
description):
==========================================================================
variable   type     description
MAXN       integer  Maximum number of trials allowed before
                    optimization is terminated.  The purpose of
                    MAXN is to stop an optimization search before
                    too much computer time is expended.  MAXN
                    should be set large enough so that optimization
                    is generally completed before MAXN trials are
                    performed. Recommended value is 10,000 (increase or
                    decrease as necessary).
---> this parameter is called m_Budget within Ostrich

KSTOP      integer  Number of shuffling loops in which the 
                    criterion must improve by the specified
                    percentage or else optimization will be
                    terminated. Recommended value: 5.
---> this parameter is called m_Kstop within Ostrich

PCENTO     double   Percentage by which the criterion value must
                    change in the specified number of shuffling 
                    loops or else optimization is terminated
                    (Use decimal equivalent: Percentage/100).
                    Recommended value: 0.01.
---> this parameter is called m_Pcento within Ostrich

NGS        integer  Number of complexes used for optimization
                    search.  Minimum value is 1.
                    Recommended value is between 2 and 20 depending
                    on the number of parameters to be optimized and
                    on the degree of difficulty of the problem.
---> this parameter is called m_NumComplexes within Ostrich

ISEED      integer  Random seed used in optimization search.  Enter
                    any integer number.  Default value (=1969) is
                    assumed if this field is left blank.
                    Recommended value: any large integer.
---> this parameter is called m_Seed within Ostrich

IDEFLT     integer  Flag for setting the control variables of the
                    SCE-UA algorithm.  Enter false or leave the field
                    blank for default setting.  Enter true for user
                    specified setting.
                    If option true is chosen, user must specify alg. 
                    parameters.
---> this parameter is called m_bUserConfig within Ostrich

NPG        integer  Number of points in each complex.  NPG should
                    be greater than or equal to 2.  The default
                    value is equal to (2 * number of optimized
                    parameters + 1).
---> this parameter is called m_PtsPerComplex within Ostrich

NPS        integer  Number of points in each sub-complex.  NPS
                    should be greater than or equal to 2 and less
                    than NPG.  The default value is equal to 
                    (number of optimized parameters + 1).
---> this parameter is called m_PtsPerSubComplex within Ostrich

NSPL       integer  Number of evolution steps taken by each complex
                    before next shuffling.  Default value is equal
                    to NPG.
---> this parameter is called m_NumEvoSteps within Ostrich

MINGS      integer  Minimum number of complexes required for
                    optimization search, if the number of complexes
                    is allowed to reduce as the optimization search
                    proceeds.  The default value is equal to NGS.
---> this parameter is called m_MinComplexes within Ostrich

INIFLG     integer  Flag on whether to include an initial point in
                    the starting population.  Enter true if the initial 
                    point is to be included.  The default value is equal to false.
---> this parameter is called m_bUseInitPt within Ostrich

IPRINT    integer   Print-out control flag.  Enter '0' to print out
                    the best estimate of the global optimum at the
                    end of each shuffling loop.  Enter '1' to print
                    out every point in the entire sample population
                    at the end of each shuffling loop.  The default
                    value is equal to 0. Enter 2 to ignore this variable
                    and use conventional Ostrich output.
---> this parameter is called m_OutputMode within Ostrich

PARAMS     double * Initial estimates of the parameters to be optimized.
---> this parameter is called m_pParams within Ostrich

LOWER      double * Lower bounds of the parameters to be optimized.
---> this parameter is called m_pLower within Ostrich

UPPER      double * Upper bounds of the parameters to be optimized.
---> this parameter is called m_pUpper within Ostrich

Version History
10-31-09    lsm   Created
******************************************************************************/
#include <math.h>
#include "SCEUA.h"
#include "WriteUtility.h"
#include "Model.h"

/******************************************************************************
CTOR

Registers the algorithm pointer and creates instances of member variables.
******************************************************************************/
SCEUA::SCEUA(ModelABC * pModel)
{
   RegisterAlgPtr(this);
   m_pModel = pModel;
   m_pParams = NULL;
   m_pUpper = NULL;
   m_pLower = NULL;
   m_bUseInitPt = false;

   IncCtorCount();
}/* end CTOR() */

/******************************************************************************
Destroy()
******************************************************************************/
void SCEUA::Destroy(void)
{ 
   delete [] m_pParams;
   delete [] m_pLower;
   delete [] m_pUpper;
}/* end Destroy() */

/******************************************************************************
Calibrate()

Solve the Least-Squares minimization problem using SCEUA.
******************************************************************************/
void SCEUA::Calibrate(void)
{ 
   int id;
   char fileName[DEF_STR_SZ];
   FILE * pFile;

   NEW_PRINT("StatsClass", 1);
   m_pStats = new StatsClass(m_pModel);
   MEM_CHECK(m_pStats);
   
   Optimize();

   //compute statistics (variance and covariance)
   m_pStats->CalcStats();

   id = 0;
   sprintf(fileName, "OstOutput%d.txt", id);

   //write statistics of best parameter set to output file
   pFile = fopen(fileName, "a");   
   m_pStats->WriteStats(pFile);
   fclose(pFile);

   //write statistics of best parameter set to output file
   m_pStats->WriteStats(stdout);
} /* end Calibrate() */

/******************************************************************************
Optimize()

Minimize the objective function using SCEUA.
******************************************************************************/
void SCEUA::Optimize(void)
{
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();

   InitFromFile(GetInFileName());

   WriteSetup(m_pModel, "Shuffled Complex Evolution - University of Arizona");
   //write banner
   WriteBanner(m_pModel, "gen   best value     ", "Pct. Complete");
  
   scemain(); //main SCE implemenation, converted from FORTRAN

   //place model at optimal prameter set
   pGroup->WriteParams(m_pParams);
   m_pModel->Execute();

   WriteOptimal(m_pModel, m_Best);
   m_pStatus.pct = 100.00;
   m_pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&m_pStatus);
   //write algorithm metrics
   WriteAlgMetrics(this);
} /* end Optimize() */

/******************************************************************************
scemain()

THIS IS THE MAIN PROGRAM CALLING SUBROUTINES SCEIN AND SCEUA
******************************************************************************/
void SCEUA::scemain(void)
{
   //implicit real*8 (a-h,o-z)
   double * a, * bl, * bu;//dimension a(16), bl(16), bu(16)
   int jseed[10] = {2,3,5,7,11,13,17,19,23,29}; //jseed(10) data jseed/2,3,5,7,11,13,17,19,23,29/

   a = m_pParams;
   bl = m_pLower;
   bu = m_pUpper;

   if(m_OutputMode != 2)
      printf(" ENTER THE MAIN PROGRAM --- \n"); //write (*,*) ' ENTER THE MAIN PROGRAM --- '

   int nopt, kstop, iseed, ngs, npg, nps, nspl, mings, iprint, maxn;
   int iniflg;
   double pcento;
   nopt = m_np = m_pModel->GetParamGroupPtr()->GetNumParams();
   scein(a,bl,bu,nopt,&maxn,&kstop,&pcento,&iseed,
         &ngs,&npg,&nps,&nspl,&mings,&iniflg,&iprint);

   /*
   if (iseed .gt. 0) then
      nrun = min(iseed,10)
   else
      nrun = 1
   end if
   */
   int nrun = 1;

   /*
   do i=1,nrun
      if (nrun .ne. 1) iseed = jseed(i)
      write (*,*) '@ SCE-UA Run Number',i,' Random Seed Value',iseed
      call sceua(a,bl,bu,nopt,maxn,kstop,pcento,iseed,
      &           ngs,npg,nps,nspl,mings,iniflg,iprint)
   end do
   */
   int i;
   for(i = 0; i < nrun; i++)
   {
      if (nrun != 1) iseed = jseed[i];
      if(m_OutputMode != 2)
         printf("@ SCE-UA Run Number %d Random Seed Value %d\n",i,iseed);
      sceua(a,bl,bu,nopt,maxn,kstop,pcento,iseed,
            ngs,npg,nps,nspl,mings,iniflg,iprint);
   }/* end for() */
}/* end scemain() */

/******************************************************************************
scein()

THIS SUBROUTINE READS AND PRINTS THE INPUT VARIABLES FOR
SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
 -- Version 2.1

ADAPTED FROM FORTRAN CODE WRITTEN BY 
    QINGYUN DUAN - UNIVERSITY OF ARIZONA, APRIL 1992
******************************************************************************/
void SCEUA::scein
(
   double * a,
   double * bl,
   double * bu,
   int nopt,
   int *maxn,
   int *kstop,
   double *pcento,
   int *iseed,
   int *ngs,
   int *npg,
   int *nps,
   int *nspl,
   int *mings,
   int *iniflg,
   int *iprint
)
{
   //implicit real*8 (a-h,o-z)
   //dimension a(16),bl(16),bu(16)
   //common /iopar/ in,ipr
   char pcntrl[100],deflt[100],usrsp[100]; //character*10 pcntrl,deflt,usrsp
   char reduc[40],initl[40],ysflg[40],noflg[40],**xname; //character*4 reduc,initl,ysflg,noflg,xname(16)

   strcpy(deflt, "DEFAULT"); //data deflt/' DEFAULT  '/
   strcpy(usrsp, "USER SPEC."); //data usrsp/'USER SPEC.'/
   strcpy(ysflg, "YES"); //data ysflg/'YES '/
   strcpy(noflg, "NO");  //data noflg/'NO  '/
   //data xname /'  X1','  X2','  X3','  X4','  X5','  X6','  X7',
   //&'  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X16'/
   xname = new char*[nopt];
   for(int i = 0; i < nopt; i++)
   {
     xname[i] = new char[50];   
     //sprintf(xname[i], "X%d", i+1);
     strncpy(xname[i], m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetName(), 9);
   }

   //write (*,*) ' ENTER THE SCEIN SUBROUTINE --- '
   if(m_OutputMode != 2) printf("ENTER THE SCEIN SUBROUTINE --- \n");

   //INITIALIZE I/O VARIABLES
   /*
   in = 2
   ipr = 5
   open(unit=in,file='scein.dat',status='old')
   open(unit=ipr,file='sceout.dat',status='unknown')
   */
   FILE * pIn  = fopen("sce.in", "r");
   FILE * pOut = fopen("sce.out", "w");

   /*
   ierror = 0
   iwarn = 0
   write(ipr,700)
   700 format(10x,'SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION',
     &       /,10x,46(1h=))
   */
   int ierror = 0;
   int iwarn = 0;
   if(m_OutputMode != 2)
   {
      fprintf(pOut, "\
          SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION\n\
          ==============================================\n\n\n");
   }

   // READ THE SCE CONTROL PARAMETERS
   /*
       ideflt = 0
       read(in,800) maxn,kstop,pcento,ngs,iseed,ideflt
   800 format(2i5,f5.2,3i5)
       if (iseed .eq. 0) iseed = 1969
   */
   int ideflt = 0;
   char line[160];
   fgets(line, 1600, pIn);
   sscanf(line, "%d %d %lf %d %d %d", maxn, kstop, pcento, ngs, iseed, &ideflt);
   if (*iseed == 0) *iseed = 1969;

  //IF ideflt IS EQUAL TO 1, READ THE SCE CONTROL PARAMETERS
  /*
      if (ideflt .eq. 1) Then
        read(in,810) npg,nps,nspl,mings,iniflg,iprint
  810   format(6i5)
        pcntrl = usrsp
      else
        read(in,*)
        pcntrl = deflt
      end if
  */
   if (ideflt == 1)
   {
      fgets(line, 1600, pIn);
      sscanf(line, "%d %d %d %d %d %d", npg, nps, nspl, mings, iniflg, iprint);
      strcpy(pcntrl, usrsp);
   }
   else
   {
      strcpy(pcntrl, deflt);
   }

   //READ THE INITIAL PARAMETER VALUES AND THE PARAMETER BOUNDS
   /*
      iopt = 0
   820 iopt = iopt + 1
      read(in,830,end=840) a(iopt), bl(iopt), bu(iopt)
   830 format(3f10.3)
      go to 820
   840 nopt = iopt - 1
   */
   int iopt;
   for(iopt = 0; iopt < nopt; iopt++)
   {
      fgets(line, 1600, pIn);
      sscanf(line,"%lf %lf %lf", &(a[iopt]), &(bl[iopt]), &(bu[iopt]));
   }

   //IF ideflt IS EQUAL TO 0, SET THE SCE CONTROL PARAMETERS TO THE DEFAULT VALUES
   /*
   if (ideflt .eq. 0) then
      npg = 2*nopt + 1
      nps = nopt + 1
      nspl = npg
      mings = ngs
      iniflg = 0
      iprint = 0
   end if
   */
   if (ideflt == 0)
   {
      *npg = 2*nopt + 1;
      *nps = nopt + 1;
      *nspl = *npg;
      *mings = *ngs;
      *iniflg = 0;
      *iprint = 0;
   }/* end if() */

   //CHECK IF THE SCE CONTROL PARAMETERS ARE VALID
   /*
   if (ngs .lt. 1 .or. ngs .ge. 1320) then
      write(ipr,900) ngs
900   format(//,1x,'**ERROR** NUMBER OF COMPLEXES IN INITIAL ',
   *         ' POPULATION ',i5,' IS NOT A VALID CHOICE')
      ierror = ierror + 1
   end if
   */
   if ((*ngs < 1) || (*ngs >= 1320))
   {
      fprintf(pOut, 
              "**ERROR** NUMBER OF COMPLEXES IN INITIAL POPULATION (%d) IS NOT A VALID CHOICE\n",
              *ngs);
      ierror = ierror + 1;
   }

/*
   if (kstop .lt. 0 .or. kstop .ge. 20) then
      write(ipr,901) kstop
   901   format(//,1x,'**WARNING** THE NUMBER OF SHUFFLING LOOPS IN',
   *  ' WHICH THE CRITERION VALUE MUST CHANGE ',/,13x,'SHOULD BE',
   *  ' GREATER THAN 0 AND LESS THAN 10.  ','kstop = ',i2,
   *  ' WAS SPECIFIED.'/,13x,'BUT kstop = 5 WILL BE USED INSTEAD.')
      iwarn = iwarn + 1
      kstop=5
   end if
*/

   if ((*kstop < 0) || (*kstop >= 20))
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF SHUFFLING LOOPS IN \
WHICH THE CRITERION VALUE MUST CHANGE SHOULD BE \
GREATER THAN 0 AND LESS THAN 10. kstop = %d WAS SPECIFIED. \
BUT kstop = 5 WILL BE USED INSTEAD.\n", 
      *kstop);
      iwarn = iwarn + 1;
      *kstop = 5;
   }/* end if() */

   /*
   if (mings .lt. 1 .or. mings .gt. ngs) then
      write(ipr,902) mings
   902   format(//,1x,'**WARNING** THE MINIMUM NUMBER OF COMPLEXES ',
   *         i2,' IS NOT A VALID CHOICE. SET IT TO DEFAULT')
      iwarn = iwarn + 1
      mings = ngs
   end if
   */
   if((*mings < 1) || (*mings > *ngs))
   {
      fprintf(pOut, 
"**WARNING** THE MINIMUM NUMBER OF COMPLEXES (%d) \
IS NOT A VALID CHOICE. SET IT TO DEFAULT \n", 
      *mings);
      iwarn = iwarn + 1;
      *mings = *ngs;
   }/* end if() */

   /*
   if (npg .lt. 2 .or. npg .gt. 1320/max(ngs,1)) then
      write(ipr,903) npg
   903   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A COMPLEX ',
   *         I4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
      iwarn = iwarn + 1
      npg = 2*nopt+1
   end if
   */
   if ((*npg < 2) || (*npg > 1320/MyMax(*ngs,1))) 
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF POINTS IN A COMPLEX (%d) \
IS NOT A VALID CHOICE, SET IT TO DEFAULT\n", 
      *npg);
      iwarn = iwarn + 1;
      *npg = 2*nopt+1;
   }/* end if() */

   /*
   if (nps.lt.2 .or. nps.gt.npg .or. nps.gt.50) then
      write(ipr,904) nps
   9 04   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A SUB-',
   *  'COMPLEX ',i4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
      iwarn = iwarn + 1
      nps = nopt + 1
   end if
   */

   if ((*nps < 2) || (*nps > *npg) || (*nps > 50))
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF POINTS IN A SUB-COMPLEX (%d) \
IS NOT A VALID CHOICE, SET IT TO DEFAULT\n",
      *nps);
      iwarn = iwarn + 1;
      *nps = nopt + 1;
   }/* end if() */

   /*
   if (nspl .lt. 1) then
      write(ipr,905) nspl
   905   format(//,1x,'**WARNING** THE NUMBER OF EVOLUTION STEPS ',
   *         'TAKEN IN EACH COMPLEX BEFORE SHUFFLING ',I4,/,13x,
   *         'IS NOT A VALID CHOICE, SET IT TO DEFAULT')
      iwarn = iwarn + 1
      nspl = npg
   end if
   */
   if (*nspl < 1)
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF EVOLUTION STEPS \
TAKEN IN EACH COMPLEX BEFORE SHUFFLING (%d) \
IS NOT A VALID CHOICE, SET IT TO DEFAULT\n", 
      *nspl);
      iwarn = iwarn + 1;
      *nspl = *npg;
   }/* end if() */
   
   // COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPULATION
   int npt;
   npt = (*ngs) * (*npg); //npt = ngs * npg

   /*
   if (npt .gt. 1320) then
      write(ipr,906) npt
   906   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN INITIAL ',
   *         'POPULATION ',i5,' EXCEED THE POPULATION LIMIT,',/,13x,
   *         'SET NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS')
      iwarn = iwarn + 1
      ngs = 2
      npg = 2*nopt + 1
      nps = nopt + 1
      nspl = npg
   end if
   */

   if (npt > 1320)
   {
      fprintf(pOut, 
"**WARNING** THE NUMBER OF POINTS IN INITIAL \
POPULATION (%d) EXCEED THE POPULATION LIMIT \
SET NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS\n", 
      npt);
      iwarn = iwarn + 1;
      *ngs = 2;
      *npg = 2*nopt + 1;
      *nps = nopt + 1;
      *nspl = *npg;
   } /* end if() */

   // PRINT OUT THE TOTAL NUMBER OF ERROR AND WARNING MESSAGES
   /* if (ierror .ge. 1) write(ipr,907) ierror
      907 format(//,1x,'*** TOTAL NUMBER OF ERROR MESSAGES IS ',i2) */
   if (ierror >= 1)
      fprintf(pOut, "*** TOTAL NUMBER OF ERROR MESSAGES IS %d\n",ierror);

   /* if (iwarn .ge. 1) write(ipr,908) iwarn
      908 format(//,1x,'*** TOTAL NUMBER OF WARNING MESSAGES IS ',i2) */
   if (iwarn >= 1)
      fprintf(pOut, "*** TOTAL NUMBER OF WARNING MESSAGES IS %d\n",iwarn);

   /*
   if (mings .lt. ngs) then
      reduc = ysflg
   else
      reduc = noflg
   end if
   */
   if (*mings < *ngs)
      strcpy(reduc, ysflg);
   else
      strcpy(reduc, noflg);

   /*   
   if (iniflg .ne. 0) then
      initl = ysflg
   else
      initl = noflg
   end if
   */

   if (iniflg != 0)
      strcpy(initl, ysflg);
   else
      strcpy(initl, noflg);

   //PRINT SHUFFLED COMPLEX EVOLUTION OPTIMIZATION OPTIONS
  /* 104 write(ipr,910)
     910 format(//,2x,'SCE CONTROL',5x,'MAX TRIALS',5x,
     &'REQUIRED IMPROVEMENT',5x,'RANDOM',/,3x,'PARAMETER',8x,
     &'ALLOWED',6x,'PERCENT',4x,'NO. LOOPS',6x,'SEED',/,
     &2x,11(1h-),5x,10(1H-),5x,7(1h-),4x,9(1h-),5x,6(1h-)) */
   fprintf(pOut,"\
  SCE CONTROL     MAX TRIALS     REQUIRED IMPROVEMENT     RANDOM\n\
   PARAMETER        ALLOWED      PERCENT    NO. LOOPS      SEED\n\
  -----------     ----------     -------    ---------     ------\n");

   /*
   pcenta=pcento*100.
   write(ipr,912) pcntrl,maxn,pcenta,kstop,iseed
   912 format(3x,a10,7x,i5,10x,f3.1,9x,i2,9x,i5)
   */
   double pcenta;
   pcenta=(*pcento)*100.;
   fprintf(pOut,"  %-11s     %-10d     %7.2lf    %-9d     %-6d\n\n\n",
   pcntrl, *maxn, pcenta, *kstop, *iseed);

   /*
      write(ipr,914) ngs,npg,npt,nps,nspl
      914 format(//,18x,'SCE ALGORITHM CONTROL PARAMETERS',/,18x,32(1H=),
     &//,2x,'NUMBER OF',5x,'POINTS PER',5x,'POINTS IN',6x,'POINTS PER',
     &4x,'EVOL. STEPS',/,2x,'COMPLEXES',6X,'COMPLEX',6x,'INI. POPUL.',
     &5x,'SUB-COMPLX',4x,'PER COMPLEX',/,2x,9(1h-),5x,10(1h-),4x,
     &11(1h-),5x,10(1h-),4x,11(1h-),5x,/,2x,5(i5,10x))
   */
   fprintf(pOut,"\
                  SCE ALGORITHM CONTROL PARAMETERS\n\
                  ================================\n\n\
  NUMBER OF     POINTS PER     POINTS IN      POINTS PER    EVOL. STEPS\n\
  COMPLEXES      COMPLEX      INI. POPUL.     SUB-COMPLX    PER COMPLEX\n\
  ---------     ----------    -----------     ----------    -----------\n");
   fprintf(pOut,"  %-9d     %-10d    %-11d     %-10d    %-11d\n\n\n", 
           *ngs, *npg, npt, *nps, *nspl);

   /*
   write(ipr,915) reduc,mings,initl
   915 format(//,15x,'COMPLX NO.',5x,'MIN COMPLEX',5x,'INI. POINT',/,
   &15x,'REDUCTION',6x,'NO. ALLOWED',6x,'INCLUDED',/,
   &15x,10(1h-),5x,11(1h-),5x,10(1h-),/,18x,a4,6x,i8,13x,a4)
   */
   fprintf(pOut,"\
               COMPLX NO.     MIN COMPLEX     INI. POINT\n\
               REDUCTION      NO. ALLOWED      INCLUDED\n\
               ----------     -----------     ----------\n");
   fprintf(pOut,"               %-10s     %-11d     %-10s\n\n\n", 
           reduc, *mings, initl);

   /*
   write(ipr,916)
   916 format(//,8x,'INITIAL PARAMETER VALUES AND PARAMETER BOUNDS',/,
     &       8x,45(1h=),//,2x,'PARAMETER',5x,'INITIAL VALUE',5x,
     &       'LOWER BOUND',5x,'UPPER BOUND',/,2x,9(1h-),5x,13(1h-),5x,
     &       11(1h-),5x,11(1h-))
      do 920 i = 1, nopt
        write(ipr,918) xname(i),a(i),bl(i),bu(i)
   918   format(4x,a4,4x,3(6x,f10.3))
   920 continue
   */
   fprintf(pOut,"\
        INITIAL PARAMETER VALUES AND PARAMETER BOUNDS\n\
        =============================================\n\n\
  PARAMETER     INITIAL VALUE     LOWER BOUND     UPPER BOUND\n\
  ---------     -------------     -----------     -----------\n");
   for(int i = 0; i < nopt; i++)
   {
      fprintf(pOut, "  %-9s     %13.3lf     %11.3lf     %11.3lf\n",
              xname[i], a[i], bl[i], bu[i]);
   }/* end for() */
   fprintf(pOut,"\n\n");

   for(int i = 0; i < nopt; i++)
   {
     delete [] xname[i];   
   }
   delete [] xname;

   /*
   if (ierror .ge. 1) then
   write(ipr,922)
   922 format(//,'*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE',
   &       ' OF INPUT DATA ERROR ***')
   stop
   end if
   */
   if (ierror >= 1)
   {
      fprintf(pOut, 
"*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE OF INPUT DATA ERROR ***\n");
      fclose(pIn);
      fclose(pOut);
      ExitProgram(1);
   }/* end if() */

   fclose(pIn);
   fclose(pOut);
}/* end scein() */

/******************************************************************************
sceua()

  SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
     -- Version 2.1

  by QINGYUN DUAN (adpated to C++ by L. Shawn Matott)
  DEPARTMENT OF HYDROLOGY & WATER RESOURCES
  UNIVERSITY OF ARIZONA, TUCSON, AZ 85721
  (602) 621-9360, email: duan@hwr.arizona.edu

  WRITTEN IN OCTOBER 1990.
  REVISED IN AUGUST 1991
  REVISED IN APRIL 1992

  STATEMENT BY AUTHOR:
  --------------------

     This general purpose global optimization program is developed at
     the Department of Hydrology & Water Resources of the University
     of Arizona.  Further information regarding the SCE-UA method can
     be obtained from Dr. Q. Duan, Dr. S. Sorooshian or Dr. V.K. Gupta
     at the address and phone number listed above.  We request all
     users of this program make proper reference to the paper entitled
     'Effective and Efficient Global Optimization for Conceptual
     Rainfall-runoff Models' by Duan, Q., S. Sorooshian, and V.K. Gupta,
     Water Resources Research, Vol 28(4), pp.1015-1031, 1992.

  LIST OF INPUT ARGUEMENT VARIABLES
     a(.) = initial parameter set
     bl(.) = lower bound on parameters
     bu(.) = upper bound on parameters
     nopt = number of parameters to be optimized

  LIST OF SCE ALGORITHMIC CONTROL PARAMETERS:
     ngs = number of complexes in the initial population
     npg = number of points in each complex
     npt = total number of points in initial population (npt=ngs*npg)
     nps = number of points in a sub-complex
     nspl = number of evolution steps allowed for each complex before
         complex shuffling
     mings = minimum number of complexes required, if the number of
         complexes is allowed to reduce as the optimization proceeds
     iseed = initial random seed
     iniflg = flag on whether to include the initial point in population
         = 0, not included
         = 1, included
     iprint = flag for controlling print-out after each shuffling loop
         = 0, print information on the best point of the population
         = 1, print information on every point of the population

  CONVERGENCE CHECK PARAMETERS
     maxn = max no. of trials allowed before optimization is terminated
     kstop = number of shuffling loops in which the criterion value must
         chang by the given percentage before optimization is terminated
     pcento = percentage by which the criterion value must change in
         given number of shuffling loops
     ipcnvg = flag indicating whether parameter convergence is reached
         (i.e., check if gnrng is less than 0.001)
         = 0, parameter convergence not satisfied
         = 1, parameter convergence satisfied

  LIST OF LOCAL VARIABLES
     x(.,.) = coordinates of points in the population
     xf(.) = function values of x(.,.)
     xx(.) = coordinates of a single point in x
     cx(.,.) = coordinates of points in a complex
     cf(.) = function values of cx(.,.)
     s(.,.) = coordinates of points in the current simplex
     sf(.) = function values of s(.,.)
     bestx(.) = best point at current shuffling loop
     bestf = function value of bestx(.)
     worstx(.) = worst point at current shuffling loop
     worstf = function value of worstx(.)
     xnstd(.) = standard deviation of parameters in the population
     gnrng = normalized geometric mean of parameter ranges
     lcs(.) = indices locating position of s(.,.) in x(.,.)
     bound(.) = bound on ith variable being optimized
     ngs1 = number of complexes in current population
     ngs2 = number of complexes in last population
     iseed1 = current random seed
     criter(.) = vector containing the best criterion values of the last
         10 shuffling loops
******************************************************************************/
void SCEUA::sceua
(
   double * a,
   double * bl,
   double * bu,
   int nopt,
   int maxn,
   int kstop,
   double pcento,
   int iseed,
   int ngs,
   int npg,
   int nps,
   int nspl,
   int mings,
   int iniflg,
   int iprint
)
{
   FILE * pOut = fopen("sce.out", "a");
   //implicit real*8 (a-h,o-z)

   //LOCAL ARRAYS
   /* dimension x(2000,16),xx(16),bestx(16),worstx(16),xf(2000)
   dimension s(50,16),sf(50),lcs(50),cx(2000,16),cf(2000)
   dimension xnstd(16),bound(16),criter(20),unit(16) */
   double ** x, * xx, * bestx, * worstx, xf[2000];
   double ** s, sf[50], ** cx, cf[2000];
   double * xnstd, * bound, criter[20], * unit;
   int lcs[50];

   //allocate memory
   x = new double * [2000];
   cx = new double * [2000];
   int i;
   for(i = 0; i < 2000; i++)
   {
      x[i] = new double[nopt];
      cx[i] = new double[nopt];
   }
   xx = new double[nopt];
   bestx = new double[nopt];
   worstx = new double[nopt];
   xnstd = new double[nopt];
   bound = new double[nopt];
   unit = new double[nopt];
   s = new double *[50];
   for(i = 0; i < 50; i++)
   {
      s[i] = new double[nopt];
   }

   //common /iopar/ in,ipr

   char **xname; //character*4 xname(16)      
   /* data xname /'  X1','  X2','  X3','  X4','  X5','  X6','  X7',
   &'  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X16'/ &/ */
   xname = new char*[nopt];
   for(i = 0; i < nopt; i++)
   {
     xname[i] = new char[50];   
     //sprintf(xname[i], "X%d", i+1);
     strncpy(xname[i], m_pModel->GetParamGroupPtr()->GetParamPtr(i)->GetName(), 9);
   }

   //write (*,*) ' ENTER THE SCEUA SUBROUTINE --- '     
   if(m_OutputMode != 2) printf("ENTER THE SCEUA SUBROUTINE --- \n");

   // INITIALIZE VARIABLES
   /*
   nloop = 0
   loop = 0
   igs = 0
   nopt1 = 8
   if (nopt.lt.8) nopt1 = nopt
   nopt2 = 12
   if (nopt.lt.12) nopt2 = nopt
   */
   int nloop, loop, igs, nopt1, nopt2;
   nloop = 0;
   loop = 0;
   igs = 0;
   nopt1 = 8;
   if (nopt < 8) nopt1 = nopt;
   nopt2 = 12;
   if (nopt < 12) nopt2 = nopt;

   //INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER
   int iseed1;
   iseed1 = -abs(iseed); //iseed1 = -abs(iseed)

   //COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
   int npt, ngs1, npt1;
   npt = ngs * npg; //npt = ngs * npg
   ngs1 = ngs; //ngs1 = ngs
   npt1 = npt; //npt1 = npt

   //write(ipr,400)
   fprintf(pOut, "\
  ==================================================\n\
  ENTER THE SHUFFLED COMPLEX EVOLUTION GLOBAL SEARCH\n\
  ==================================================\n\n\n");

   //write (*,*) ' ***  Evolution Loop Number ',nloop
   if(m_OutputMode != 2) printf(" ***  Evolution Loop Number %d\n", nloop);

   //COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
   /*
   do j = 1, nopt
      bound(j) = bu(j) - bl(j)
      unit(j) = 1.0
   end do
   */
   int j;
   for(j = 1; j <= nopt; j++)
   {
      bound[j-1] = bu[j-1] - bl[j-1];
      unit[j-1] = 1.0;
   }

   //COMPUTE THE FUNCTION VALUE OF THE INITIAL POINT
   //fa = functn(nopt,a)
   double fa;
   m_pModel->GetParamGroupPtr()->WriteParams(a);
   fa = m_pModel->Execute();

   //write initial config.
   int nleft = m_Budget - m_pModel->GetCounter();
   
   m_pStatus.curIter = 0;
   m_pStatus.maxIter = m_Budget;
   m_pStatus.pct = (float)100.00*((float)1.00-(float)nleft/(float)m_Budget);
   m_pStatus.numRuns = m_pModel->GetCounter();
   m_pModel->GetParamGroupPtr()->WriteParams(m_pParams);
   WriteRecord(m_pModel, 0, fa, m_pStatus.pct);
   WriteStatus(&m_pStatus);

   //PRINT THE INITIAL POINT AND ITS CRITERION VALUE
   /*
   write(ipr,500)
   write(ipr,510) (xname(j),j=1,nopt2)
   write(ipr,520) fa,(a(j),j=1,nopt2)
   if (nopt.lt.13) go to 101
   write(ipr,530) (xname(j),j=13,nopt)
   write(ipr,540) (a(j),j=13,nopt)
   101 continue
   */
   fprintf(pOut, "\
*** PRINT THE INITIAL POINT AND ITS CRITERION VALUE ***\n\n\
 CRITERION    ");
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut, "%-9s    ", xname[i]);
   }
   fprintf(pOut, "\n  %8.0lf     ", fa);
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut, "%5.3lf     ", a[i]);
   }
   fprintf(pOut, "\n\n\n");

   // GENERATE AN INITIAL SET OF npt1 POINTS IN THE PARAMETER SPACE
   // IF iniflg IS EQUAL TO 1, SET x(1,.) TO INITIAL POINT a(.)
   /* if (iniflg .eq. 1) then
      do j = 1, nopt
         x(1,j) = a(j)
      end do
      xf(1) = fa */
   if (iniflg == 1)
   {
      for(j = 0; j < nopt; j++)
         x[0][j] = a[j];
      xf[0] = fa;
      WriteInnerEval(WRITE_SCE, npt, '.');
      WriteInnerEval(1, npt, '.');
   }/* end if() */
   // ELSE, GENERATE A POINT RANDOMLY AND SET IT EQUAL TO x(1,.)
   /* 
   else
      call getpnt(nopt,1,iseed1,xx,bl,bu,unit,bl)
      do j=1, nopt
         x(1,j) = xx(j)
      end do
      xf(1) = functn(nopt,xx)
   end if
   icall = 1
   if (icall .ge. maxn) go to 9000
   */
   else
   {
      getpnt(nopt,1,&iseed1,xx,bl,bu,unit,bl);
      for(j = 0; j < nopt; j++)
         x[0][j] = xx[j];
      WriteInnerEval(WRITE_SCE, npt, '.');
      WriteInnerEval(1, npt, '.');
      m_pModel->GetParamGroupPtr()->WriteParams(xx);
      xf[0] = m_pModel->Execute();
   }
   int icall;
   icall = 1;
   if (icall >= maxn) goto label_9000;

   //GENERATE npt1-1 RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
   //SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES
   /*
   do i = 2, npt1
      call getpnt(nopt,1,iseed1,xx,bl,bu,unit,bl)
      do j = 1, nopt
         x(i,j) = xx(j)
      end do
      xf(i) = functn(nopt,xx)
      icall = icall + 1
      if (icall .ge. maxn) then
         npt1 = i
         go to 45
      end if
   end do
   */
   for(i = 1; i < npt1; i++)
   {
      getpnt(nopt,1,&iseed1,xx,bl,bu,unit,bl);
      for(j = 0; j < nopt; j++)
      {
         x[i][j] = xx[j];
      }
      WriteInnerEval(i+1, npt, '.');
      m_pModel->GetParamGroupPtr()->WriteParams(xx);
      xf[i] = m_pModel->Execute();
      icall = icall + 1;
      if (icall >= maxn)
      {
         npt1 = i;
         break;
      }
   }/* end for() */
   WriteInnerEval(WRITE_ENDED, npt, '.');

   // ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE
   //45 call sort(npt1,nopt,x,xf)
   sort(npt1,nopt,x,xf);

   // RECORD THE BEST AND WORST POINTS
   /*
   do j = 1, nopt
      bestx(j) = x(1,j)
      worstx(j) = x(npt1,j)
   end do
   bestf = xf(1)
   worstf = xf(npt1)
   */
   for(j = 0; j < nopt; j++)
   {
      bestx[j] = x[0][j];
      worstx[j] = x[npt1-1][j];
   }
   double bestf, worstf;
   bestf = xf[0];
   worstf = xf[npt1-1];

   // COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
   // call parstt(npt1,nopt,x,xnstd,bound,gnrng,ipcnvg)
   double gnrng;
   int ipcnvg;
   parstt(npt1,nopt,x,xnstd,bound,&gnrng,&ipcnvg);

   // PRINT THE RESULTS FOR THE INITIAL POPULATION
   /*
   write(ipr,600)
   write(ipr,610) (xname(j),j=1,nopt1)
   if (nopt .lt. 9) go to 201
   write(ipr,620) (xname(j),j=9,nopt)
   201 continue
   write(ipr,630) nloop,icall,ngs1,bestf,worstf,gnrng,
   &               (bestx(j),j=1,nopt1)
   if (nopt .lt. 9) go to 301
   write(ipr,640) (bestx(j),j=9,nopt)
   301 continue
   */
   fprintf(pOut,"\
**** PRINT THE RESULTS OF THE SCE SEARCH ***\n\n\
 LOOP TRIALS COMPLXS  BEST F   WORST F   PAR RNG         ");
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut,"%-9s ", xname[i]);
   }
   fprintf(pOut,"\n");
   fprintf(pOut," %4d %6d %7d  %6.2lf  %9.3E  %8.3lf      ",
           nloop,icall,ngs1,bestf,worstf,gnrng);
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut,"%6.3lf    ", bestx[i]);
   }
   fprintf(pOut,"\n");

   /*
   if (iprint .eq. 1) then
   write(ipr,650) nloop
      do i = 1, npt1
         write(ipr,660) xf(i),(x(i,j),j=1,nopt1)
         if (nopt .lt. 9) go to 401
         write(ipr,640) (x(i,j),j=9,nopt)
      401   end do
   end if
   */
   if (iprint == 1)
   {
      fprintf(pOut, "POPULATION AT LOOP (%d)\n", nloop);
      for(i = 0; i < npt1; i++)
      {
         fprintf(pOut, "%8.3lf    ", xf[i]);
         for(j = 0; j < nopt; j++)
         {
            fprintf(pOut, "%8.3lf    ", x[i][j]);
         }
         fprintf(pOut, "\n");
      }/* end for() */
   }/* end if() */

   if (icall >= maxn) goto label_9000; //if (icall .ge. maxn) go to 9000
   if (ipcnvg == 1) goto label_9200; //if (ipcnvg .eq. 1) go to 9200

   // BEGIN THE MAIN LOOP ----------------
 label_1000:

   nleft = m_Budget - m_pModel->GetCounter();
   m_pStatus.curIter = nloop+1;
   if(IsQuit() == true){ goto label_9999;}
   if(nleft <= 0){ m_pStatus.pct = 100.00;  goto label_9000;}

   nloop = nloop + 1; //nloop = nloop + 1

   //write (*,*) ' ***  Evolution Loop Number ',nloop
   if(m_OutputMode != 2) printf(" ***  Evolution Loop Number %d\n",nloop); 
   
   //BEGIN LOOP ON COMPLEXES
   for(igs = 1; igs <= ngs1; igs++) //do 3000 igs = 1, ngs1
   {
      // ASSIGN POINTS INTO COMPLEXES
      int k1, k2;
      for(k1 = 1; k1 <= npg; k1++) //do k1 = 1, npg
      {
        k2 = (k1-1) * ngs1 + igs; //k2 = (k1-1) * ngs1 + igs
        for(j = 1; j <= nopt; j++) //do j = 1, nopt
        {
            cx[k1-1][j-1] = x[k2-1][j-1]; //cx(k1,j) = x(k2,j)
        } //end do
        cf[k1-1] = xf[k2-1]; // cf(k1) = xf(k2)
      } //end do
      // BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
      int tmp = 0;
      WriteInnerEval(WRITE_SCE, m_NumEvoSteps, '.');
      //icall = 0;
      for(loop = 0; loop < nspl; loop++) // do 2000 loop = 1, nspl
      {
         // CHOOSE A SUB-COMPLEX (nps points) ACCORDING TO A LINEAR
         // PROBABILITY DISTRIBUTION
         if (nps == npg) // if (nps .eq. npg) then
         {
            int k;
            for(k = 0; k < nps; k++) //do k = 1, nps
            {
               lcs[k] = k; //lcs(k) = k
            } // end do
            goto label_85; //go to 85
         } // end if

         double myrand;
         myrand = ran1(&iseed1); //rand = ran1(iseed1)
         //lcs(1) = 1 + dint(npg + 0.5 - dsqrt( (npg+.5)**2 -
         //&         npg * (npg+1) * rand ))
         lcs[0] = (int)(npg + 0.5 - sqrt(pow((npg+0.5),2.00) - npg*(npg+1.00)*myrand));

         int k, lpos;
         for(k = 2; k <= nps; k++) // do k = 2, nps
         {
label_60:
            myrand = ran1(&iseed1); //rand = ran1(iseed1)
            //lpos = 1 + dint(npg + 0.5 - dsqrt((npg+.5)**2 -
            //&  npg * (npg+1) * rand ))
            lpos = (int)(npg + 0.5 - sqrt(pow((npg+0.5),2.00) - npg*(npg+1.00)*myrand));

            for(k1 = 1; k1 <= k-1; k1++) //do k1 = 1, k-1
            {
               if (lpos == lcs[k1-1]) goto label_60; //if (lpos .eq. lcs(k1)) go to 60
            } // end do
            lcs[k-1] = lpos; // lcs(k) = lpos
         } // end do

         // ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
         sort1(nps,lcs); //call sort1(nps,lcs)

         // CREATE THE SUB-COMPLEX ARRAYS
label_85: 
         for(k = 1; k <= nps; k++) // do k = 1, nps
         {
            for(j = 1; j <= nopt; j++) //do j = 1, nopt
            {
               s[k-1][j-1] = cx[lcs[k-1]][j-1]; //s(k,j) = cx(lcs(k),j)
            } // end do
            sf[k-1] = cf[lcs[k-1]]; //sf(k) = cf(lcs(k))
         } // end do

         // USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
         //call cce(nopt,nps,s,sf,bl,bu,xnstd,icall,maxn,iseed1)
         cce(nopt,nps,s,sf,bl,bu,xnstd,&tmp,maxn,&iseed1);

         // IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
         // INTO THE COMPLEX
         for(k = 1; k <= nps; k++) //do k = 1, nps
         {
            for(j = 1; j <= nopt; j++) //do j = 1, nopt
            {
               cx[lcs[k-1]][j-1] = s[k-1][j-1]; //cx(lcs(k),j) = s(k,j)
            } // end do
            cf[lcs[k-1]] = sf[k-1]; // cf(lcs(k)) = sf(k)
         } //end do

         // SORT THE POINTS
         sort(npg,nopt,cx,cf); //call sort(npg,nopt,cx,cf)

         //IF MAXIMUM NUMBER OF RUNS EXCEEDED, BREAK OUT OF THE LOOP
         if (icall >= maxn) break; //if (icall .ge. maxn) go to 2222
         // END OF INNER LOOP ------------
      } /* end for() */ // 2000 continue
      //2222 continue
      WriteInnerEval(WRITE_ENDED, m_NumEvoSteps, '.');
      icall += tmp;

      // REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)
      for(k1 = 1; k1 <= npg; k1++) //do k1 = 1, npg
      {
         k2 = (k1-1) * ngs1 + igs; //k2 = (k1-1) * ngs1 + igs
         for(j = 1; j <= nopt; j++) //do j = 1, nopt
         {
            x[k2-1][j-1] = cx[k1-1][j-1]; //x(k2,j) = cx(k1,j)
         } // end do
         xf[k2-1] = cf[k1-1]; //xf(k2) = cf(k1)
      } // end do
      if (icall >= maxn) break; //if (icall .ge. maxn) go to 3333
      //END LOOP ON COMPLEXES
   } /* end for() */ //3000 continue

   // RE-SORT THE POINTS
   sort(npt1,nopt,x,xf); //3333 call sort(npt1,nopt,x,xf)

   // RECORD THE BEST AND WORST POINTS
   for(j = 0; j < nopt; j++) //do j = 1, nopt
   {
      m_pParams[j] = bestx[j] = x[0][j]; //bestx(j) = x(1,j)
      worstx[j] = x[npt1-1][j]; //worstx(j) = x(npt1,j)
   } //end do
   m_Best = bestf = xf[0]; //bestf = xf(1)
   worstf = xf[npt1-1]; //worstf = xf(npt1)

   // TEST THE POPULATION FOR PARAMETER CONVERGENCE
   parstt(npt1,nopt,x,xnstd,bound,&gnrng,&ipcnvg); // call parstt(npt1,nopt,x,xnstd,bound,gnrng,ipcnvg)

   // PRINT THE RESULTS FOR CURRENT POPULATION
   m_pModel->GetParamGroupPtr()->WriteParams(m_pParams);
   nleft = m_Budget - m_pModel->GetCounter();
   m_pStatus.pct = (float)100.00*((float)1.00-(float)nleft/(float)m_Budget);
   m_pStatus.numRuns = m_pModel->GetCounter();
   WriteStatus(&m_pStatus);
   WriteRecord(m_pModel, nloop, m_Best, m_pStatus.pct);

   //if (mod(nloop,5) .ne. 0) go to 501
   //write(ipr,610) (xname(j),j=1,nopt1)
   //if (nopt .lt. 9) go to 501
   //write(ipr,620) (xname(j),j=9,nopt)
   //501 continue
   //write(ipr,630) nloop,icall,ngs1,bestf,worstf,gnrng,
   //&               (bestx(j),j=1,nopt1)
   //if (nopt.lt.9) go to 601
   //write(ipr,640) (bestx(j),j=9,nopt)
   //601 continue
   if((nloop%5) == 0)
   {
      fprintf(pOut,"\
 LOOP TRIALS COMPLXS  BEST F   WORST F   PAR RNG         ");
      for(i = 0; i < nopt; i++)
      {
         fprintf(pOut,"%-9s ", xname[i]);
      }
      fprintf(pOut,"\n");
   }
   fprintf(pOut," %4d %6d %7d  %6.2lf  %9.3E  %8.3lf      ",
           nloop,icall,ngs1,bestf,worstf,gnrng);
   for(i = 0; i < nopt; i++)
   {
      fprintf(pOut,"%6.3lf    ", bestx[i]);
   }
   fprintf(pOut,"\n");

   //if (iprint .eq. 1) then
   //write(ipr,650) nloop
   //do i = 1, npt1
   //write(ipr,660) xf(i),(x(i,j),j=1,nopt1)
   //if (nopt .lt. 9) go to 701
   //write(ipr,640) (x(i,j),j=9,nopt)
   //701   end do
   //end if
   if (iprint == 1)
   {
      fprintf(pOut, "POPULATION AT LOOP (%d)\n", nloop);
      for(i = 0; i < npt1; i++)
      {
         fprintf(pOut, "%8.3lf    ", xf[i]);
         for(j = 0; j < nopt; j++)
         {
            fprintf(pOut, "%8.3lf    ", x[i][j]);
         }
         fprintf(pOut, "\n");
      }/* end for() */
   }/* end if() */

   // TEST IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED
   if (icall >= maxn) goto label_9000; //if (icall .ge. maxn) go to 9000

   // COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
   criter[19] = bestf; //criter(20) = bestf
   if (nloop < (kstop+1)) goto label_132; //if (nloop .lt. (kstop+1)) go to 132
   double denomi, timeou;
   denomi = fabs(criter[19-kstop] + criter[19]) / 2.0; //denomi = dabs(criter(20-kstop) + criter(20)) / 2.
   timeou = fabs(criter[19-kstop] - criter[19]) / denomi; //timeou = dabs(criter(20-kstop) - criter(20)) / denomi
   if (timeou < pcento) goto label_9100; //if (timeou .lt. pcento) go to 9100
label_132: 
   int l;
   for(l=0; l < 19; l++) //do l = 1, 19
   {
        criter[l] = criter[l+1]; //criter(l) = criter(l+1)
   } //end do


   //IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
   if (ipcnvg == 1) goto label_9200; //if (ipcnvg .eq. 1) go to 9200

   //NONE OF THE STOPPING CRITERIA IS SATISFIED, CONTINUE SEARCH

   //CHECK FOR COMPLEX NUMBER REDUCTION
   int ngs2;
   if (ngs1 > mings) //if (ngs1 .gt .mings) then
   {
        ngs2 = ngs1; //ngs2 = ngs1
        ngs1 = ngs1 - 1; //ngs1 = ngs1 - 1
        npt1 = ngs1 * npg; //npt1 = ngs1 * npg
        comp(nopt,npt1,ngs1,ngs2,npg,x,xf,cx,cf); //call comp(nopt,npt1,ngs1,ngs2,npg,x,xf,cx,cf)
   } //end if

   // END OF MAIN LOOP -----------
   goto label_1000;

   // SEARCH TERMINATED
label_9000: 
      //write(ipr,800) maxn,loop,igs,nloop
      fprintf(pOut, "\
*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE\n\
    LIMIT ON THE MAXIMUM NUMBER OF TRIALS (%d)\n\
    WAS EXCEEDED.  SEARCH WAS STOPPED AT %d SUB-COMPLEX\n\
    OF COMPLEX %d IN SHUFFLING LOOP %d ***\n\n",
    maxn, loop, igs, nloop);

      goto label_9999; //go to 9999

label_9100:
      //write(ipr,810) pcento*100.,kstop
      fprintf(pOut, "\
*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION\n\
    VALUE HAS NOT CHANGED %5.2lf PERCENT IN %d\n\
    SHUFFLING LOOPS ***\n\n", pcento*100.0,kstop);

      goto label_9999; //go to 9999

 label_9200:
      //write(ipr,820) gnrng*100.
      fprintf(pOut, "\
 *** OPTIMIZATION TERMINATED BECAUSE THE POPULATION HAS\n\
     CONVERGED INTO %5.2lf PERCENT OF THE FEASIBLE SPACE ***\n\n",
      gnrng*100.0);

 label_9999:
      // PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
      //write(ipr,830)
      fprintf(pOut, "\
*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS CRITERION VALUE ***\n\n\
 CRITERION        ");

      //write(ipr,510) (xname(j),j=1,nopt2)
      //write(ipr,520) bestf,(bestx(j),j=1,nopt2)
      //if (nopt.lt.13) go to 801
      //write(ipr,530) (xname(j),j=13,nopt)
      //write(ipr,540) (bestx(j),j=13,nopt)
      for(i = 0; i < nopt; i++)
      {
         fprintf(pOut,"%-9s ", xname[i]);
      }
      fprintf(pOut,"\n%6.3lf    ", bestf);
      for(i = 0; i < nopt; i++)
      {
         fprintf(pOut,"%6.3lf    ", bestx[i]);
      }
      fprintf(pOut,"\n");

//label_801:
   // END OF SUBROUTINE SCEUA
   // return
   //  400 format(//,2x,50(1h=),/,2x,'ENTER THE SHUFFLED COMPLEX EVOLUTION',
   //     &       ' GLOBAL SEARCH',/,2x,50(1h=))
   //  500 format(//,'*** PRINT THE INITIAL POINT AND ITS CRITERION ',
   //     &       'VALUE ***')
   //  510 format(/,' CRITERION',12(6x,a4),/1x,60(1h-))
   //  520 format(g10.3,12f10.3)
   //  530 format(10x,12(6x,a4))
   //  540 format(10x,12f10.3)
   //  600 format(//,1x,'*** PRINT THE RESULTS OF THE SCE SEARCH ***')
   //  610 format(/,1x,'LOOP',1x,'TRIALS',1x,'COMPLXS',2x,'BEST F',3x,
   //     &       'WORST F',3x,'PAR RNG',1x,8(6x,a4))
   //  620 format(49x,8(6x,a4))
   //  630 format(i5,1x,i5,3x,i5,3g10.3,8(f10.3))
   //  640 format(49x,8(f10.3))
   //  650 format(/,1x,'POPULATION AT LOOP ',i3,/,1x,22(1h-))
   //  660 format(15x,g10.3,20x,8(f10.3))
   //  800 format(//,1x,'*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE',
   //     &       ' LIMIT ON THE MAXIMUM',/,5x,'NUMBER OF TRIALS ',i5,
   //     &       ' EXCEEDED.  SEARCH WAS STOPPED AT',/,5x,'SUB-COMPLEX ',
   //     &       i3,' OF COMPLEX ',i3,' IN SHUFFLING LOOP ',i3,' ***')
   //810 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION',
   //   &       ' VALUE HAS NOT CHANGED ',/,5x,f5.2,' PERCENT IN',i3,
   //   &       ' SHUFFLING LOOPS ***')
   //820 format(//,1x,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION',
   //   &       ' HAS CONVERGED INTO ',/,4x,f5.2,' PERCENT OF THE',
   //   &       ' FEASIBLE SPACE ***')
   //830 format(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS',
   //   &       ' CRITERION VALUE ***')

   fclose(pOut);
   /* clean up memory */
   for(i = 0; i < nopt; i++)
   {
     delete [] xname[i];   
   }
   delete [] xname;
   for(i = 0; i < 2000; i++)
   {
      delete [] x[i];
      delete [] cx[i];
   }
   delete [] x;
   delete [] cx;
   delete [] xx;
   delete [] bestx;
   delete [] worstx;
   delete [] xnstd;
   delete [] bound;
   delete [] unit;
   for(i = 0; i < 50; i++)
   {
      delete [] s[i];
   }
   delete [] s;
} /* end sceua() */

/******************************************************************************
cce()

ALGORITHM TO GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
******************************************************************************/
void SCEUA::cce
(
   int nopt,
   int nps,
   double ** s,
   double * sf,
   double * bl,
   double * bu,
   double * xnstd,
   int * icall,
   double maxn,
   int * iseed
)
{
   //SUB-COMPLEX VARIABLES
   //implicit real*8 (a-h,o-z)
   const double c1=0.8;
   const double c2=0.4;
   //dimension s(50,16),sf(50),bu(16),bl(16),xnstd(16)

   /* ----------------------------------------------------
   LIST OF LOCAL VARIABLES
      sb(.) = the best point of the simplex
      sw(.) = the worst point of the simplex
      w2(.) = the second worst point of the simplex
      fw = function value of the worst point
      ce(.) = the centroid of the simplex excluding wo
      snew(.) = new point generated from the simplex
      iviol = flag indicating if constraints are violated
            = 1 , yes
            = 0 , no
   ----------------------------------------------------- */
   double * sw, * sb, *ce, *snew; //dimension sw(16),sb(16),ce(16),snew(16)
   sw = new double[nopt];
   sb = new double[nopt];
   ce = new double[nopt];
   snew = new double[nopt];

   //EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
   int n = nps;
   int m = nopt;
   double alpha = 1.0;
   double beta = 0.5;

   /* ---------------------------------------------------
   IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
   COMPUTE THE CENTROID ce OF THE REMAINING POINTS
   COMPUTE step, THE VECTOR BETWEEN wo AND ce
   IDENTIFY THE WORST FUNCTION VALUE fw
   --------------------------------------------------- */
   int i, j;
   for(j = 0; j < m; j++) //do j = 1, m
   {
      sb[j] = s[0][j]; //sb(j) = s(1,j)
      sw[j] = s[n-1][j]; //sw(j) = s(n,j)
      ce[j] = 0.0; //ce(j) = 0.0
      for(i = 0; i < n-1; i++) //do i = 1, n-1
      {
         ce[j] = ce[j] + s[i][j]; //ce(j) = ce(j) + s(i,j)
      } //end do
     ce[j] = ce[j]/(double)(n-1); //ce(j) = ce(j)/dble(n-1)
   } //end do
   double fw;
   fw = sf[n-1]; //fw = sf(n)

   //COMPUTE THE NEW POINT snew
   //FIRST TRY A REFLECTION STEP
   for(j = 0; j < m; j++) //do j = 1, m
   {
      //snew(j) = ce(j) + alpha * (ce(j) - sw(j))
      snew[j] = ce[j] + alpha * (ce[j] - sw[j]);
   } //end do

   //printf("(1) %E  %E\n", snew[0], snew[1]);
   //CHECK IF snew SATISFIES ALL CONSTRAINTS
   int ibound;
   chkcst(nopt,snew,bl,bu,&ibound); //call chkcst(nopt,snew,bl,bu,ibound)

   /* ------------------------------------------------------------------
   snew IS OUTSIDE THE BOUND,
   CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
   A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
   AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
   ------------------------------------------------------------------ */
   //if (ibound .ge. 1) call getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb)
   if (ibound >= 1) getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb);
   //printf("(2) %E  %E\n", snew[0], snew[1]);

   //COMPUTE THE FUNCTION VALUE AT snew
   double fnew;
   WriteInnerEval(*icall+1, m_NumEvoSteps, '.');
   m_pModel->GetParamGroupPtr()->WriteParams(snew);
   fnew = m_pModel->Execute(); //fnew = functn(nopt,snew)
   *icall = *icall + 1;

   //COMPARE fnew WITH THE WORST FUNCTION VALUE fw
   //fnew IS LESS THAN fw, ACCEPT THE NEW POINT snew AND RETURN
   if (fnew <= fw) goto label_2000;   // if (fnew .le. fw) go to 2000
   if (*icall >= maxn) goto label_3000; //if (icall .ge. maxn) go to 3000

   //fnew IS GREATER THAN fw, SO TRY A CONTRACTION STEP
   for(j = 0; j < m; j++) //do j = 1, m
   {
      //snew(j) = ce(j) - beta * (ce(j) - sw(j))
      snew[j] = ce[j] - beta * (ce[j] - sw[j]);
   } //end do

   //COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
   WriteInnerEval(*icall+1, m_NumEvoSteps, '.');
   m_pModel->GetParamGroupPtr()->WriteParams(snew);
   fnew = m_pModel->Execute(); //fnew = functn(nopt,snew)
   *icall = *icall + 1;

   //COMPARE fnew TO THE WORST VALUE fw
   //IF fnew IS LESS THAN OR EQUAL TO fw, THEN ACCEPT THE POINT AND RETURN
   if (fnew <= fw) goto label_2000; //if (fnew .le. fw) go to 2000
   if (*icall >= maxn) goto label_3000; //if (icall .ge. maxn) go to 3000

   /* ---------------------------------------------------------------------
   IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
   ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
   AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
   --------------------------------------------------------------------- */
//label_1000:
   //call getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb)
   getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb);

   //COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
   WriteInnerEval(*icall+1, m_NumEvoSteps, '.');
   m_pModel->GetParamGroupPtr()->WriteParams(snew);
   fnew = m_pModel->Execute(); //fnew = functn(nopt,snew)
   *icall = *icall + 1;

   //REPLACE THE WORST POINT BY THE NEW POINT
label_2000:
   for(j = 0; j < m; j++) //do j = 1, m
   {
      s[n-1][j] = snew[j]; //s(n,j) = snew(j)
   } //end do
   sf[n-1] = fnew; //sf(n) = fnew

label_3000:
   //free up memory  
   delete [] sw;
   delete [] sb;
   delete [] ce;
   delete [] snew;
   //c  END OF SUBROUTINE CCE
   // return
   // end
} /* end cce() */

/******************************************************************************
getpnt()

This subroutine generates a new point within feasible region

x(.) = new point
xi(.) = focal point
bl(.) = lower bound
bu(.) = upper bound
std(.) = standard deviation of probability distribution
idist = probability flag
      = 1 - uniform distribution
      = 2 - Gaussian distribution
******************************************************************************/
void SCEUA::getpnt
(
   int nopt,
   int idist,
   int * iseed,
   double * x,
   double * bl,
   double * bu,
   double * std,
   double * xi
)
{
   //implicit real*8 (a-h,o-z)
   //dimension x(16),bl(16),bu(16),std(16),xi(16)
   int ibound;
   int j;
   double myrand;

   int icount = m_pModel->GetCounter();
label_1:
   for(j = 0; j < nopt; j++) //1   do j=1, nopt
   {
label_2:
      //2   if (idist .eq. 1) rand = ran1(iseed)
      if (idist == 1)
      {
         myrand = ran1(iseed);
      }
      //if (idist .eq. 2) rand = gasdev(iseed)
      else if (idist == 2)
      {
         myrand = gasdev(iseed);
      }
      else
      {
         printf("unknown distribution!\n");
      }

      //x(j) = xi(j) + std(j) * myrand * (bu(j) - bl(j))
      x[j] = xi[j] + std[j] * myrand * (bu[j] - bl[j]);

      FILE * pFile = fopen("getpnt.txt", "a+");
      fprintf(pFile, "%d\tx[%d]:%E\txi[%d]:%e\tstd[%d]:%E\tmyrand : %E\tbu[%d]:%E\tbl[%d]:%E\n",
                     icount, j+1, x[j], j+1, xi[j], j+1, std[j], myrand, j+1, bu[j], j+1, bl[j]);
      fclose(pFile);

      //Check explicit constraints
      //call chkcst(1,x(j),bl(j),bu(j),ibound)
      chkcst(1,&x[j],&bl[j],&bu[j],&ibound);
      //if (ibound .ge. 1) go to 2
      if (ibound >= 1)
      {
         goto label_2;
      }
   } //end do

   //Check implicit constraints    
   chkcst(nopt,x,bl,bu,&ibound); //call chkcst(nopt,x,bl,bu,ibound)
   if (ibound >= 1)
   {
      goto label_1; //if (ibound .ge. 1) go to 1
   }
   //return
   //end
}/* end getpnt() */

/******************************************************************************
parstt()

SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
******************************************************************************/
void SCEUA::parstt
(
   int npt,
   int nopt,
   double ** x,
   double * xnstd,
   double * bound,
   double * gnrng,
   int * ipcnvg
)
{
   //implicit real*8 (a-h,o-z)
   //dimension x(2000,16),xmax(16),xmin(16)
   //dimension xmean(16),xnstd(16),bound(16)
   double * xmax, * xmin, * xmean;
   xmax = new double[nopt];
   xmin = new double[nopt];
   xmean = new double[nopt];

   //parameter (delta = 1.0d-20,peps=1.0d-3)
   const double delta = 1.0e-20;
   const double peps=1.0E-3;

   //COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
   double gsum, xsum1, xsum2;
   int i,k;
   gsum = 0.0; //gsum = 0.d0
   for(k = 0; k < nopt; k++) //do k = 1, nopt
   {
      xmax[k] = -1.0E+20; //xmax(k) = -1.0d+20
      xmin[k] = 1.0E+20; //xmin(k) = 1.0d+20
      xsum1 = 0.0; //xsum1 = 0.d0
      xsum2 = 0.0; //xsum2 = 0.d0

      for(i = 0; i < npt; i++) //do i = 1, npt
      {
         xmax[k] = MyMax(x[i][k], xmax[k]); //xmax(k) = dmax1(x(i,k), xmax(k))
         xmin[k] = MyMin(x[i][k], xmin[k]); //xmin(k) = dmin1(x(i,k), xmin(k))
         xsum1 = xsum1 + x[i][k]; //xsum1 = xsum1 + x(i,k)
         xsum2 = xsum2 + x[i][k]*x[i][k]; //xsum2 = xsum2 + x(i,k)*x(i,k)
      } //end do

      xmean[k] = xsum1/(double)npt; //xmean(k) = xsum1 / dble(npt)
      //xnstd(k) = (xsum2 / dble(npt) - xmean(k)*xmean(k))
      xnstd[k] = (xsum2 / (double)npt - xmean[k]*xmean[k]);

      //if (xnstd(k) .le. delta) xnstd(k) = delta
      if (xnstd[k] <= delta) xnstd[k] = delta; 
      xnstd[k] = sqrt(xnstd[k]); //xnstd(k) = dsqrt(xnstd(k))
      xnstd[k] = xnstd[k] / bound[k]; //xnstd(k) = xnstd(k) / bound(k)
      //gsum = gsum + dlog( delta + (xmax(k)-xmin(k))/bound(k) )
      gsum = gsum + log(delta + (xmax[k]-xmin[k])/bound[k]);
   } //end do
   *gnrng = exp(gsum/(double)nopt); //gnrng = dexp(gsum/dble(nopt))

   //CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
   *ipcnvg = 0;
   if (*gnrng <= peps) //if (gnrng .le. peps) then
   {
      *ipcnvg = 1;
   } //end if

   delete [] xmax;
   delete [] xmin;
   delete [] xmean;

   //END OF SUBROUTINE PARSTT
   //return
   //end
}/* end parstt() */

/******************************************************************************
comp()

THIS SUBROUTINE REDUCE INPUT MATRIX a(n,ngs2*npg) TO MATRIX
b(n,ngs1*npg) AND VECTOR af(ngs2*npg) TO VECTOR bf(ngs1*npg)
******************************************************************************/
void SCEUA::comp
(
   int n,
   int npt,
   int ngs1,
   int ngs2,
   int npg,
   double ** a,
   double * af,
   double ** b,
   double * bf
)
{
   //implicit real*8 (a-h,o-z)
   //dimension a(2000,16),af(2000),b(2000,16),bf(2000)
   int i, igs, ipg, k1,  k2;
   for(igs = 1; igs <= ngs1; igs++) //do igs=1, ngs1
   {
      for(ipg = 1; ipg <= npg; ipg++) //do ipg=1, npg
      {
         k1=(ipg-1)*ngs2 + igs; //k1=(ipg-1)*ngs2 + igs
         k2=(ipg-1)*ngs1 + igs; //k2=(ipg-1)*ngs1 + igs
         for(i = 1; i <= n; i++) //do i=1, n
         {
            b[k2-1][i-1] = a[k1-1][i-1]; //b(k2,i) = a(k1,i)
         } //end do
         bf[k2-1] = af[k1-1]; //bf(k2) = af(k1)
      } //end do
   } //end do

   int j;
   for(j = 0; j < npt; j++) //do j=1, npt
   {
      for(i = 0; i < n; i++) //do i=1, n
      {
         a[j][i] = b[j][i]; //a(j,i) = b(j,i)
      } //end do
      af[j] = bf[j];
   } //end do

   //END OF SUBROUTINE COMP
   //return
   //end
} /* end comp() */

/******************************************************************************
sort()

SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
BY W.H. PRESS ET AL., pp. 233-234

LIST OF VARIABLES
   ra(.) = array to be sorted
   rb(.,.) = arrays ordered corresponding to rearrangement of ra(.)
   wk(.,.), iwk(.) = local varibles
******************************************************************************/
void SCEUA::sort
(
   int n,
   int m,
   double ** rb,
   double * ra
)
{
   //implicit real*8 (a-h,o-z)
   //dimension ra(2000),rb(2000,16),wk(2000,16),iwk(2000)
   double ** wk;
   int iwk[2000];

   wk = new double *[2000];
   int i;
   for(i = 0; i < 2000; i++) wk[i] = new double[n];

   indexx(n, ra, iwk); //call indexx(n, ra, iwk)

   for(i = 0; i < n; i++) //do 11 i = 1, n
   {
      wk[i][0] = ra[i]; //wk(i,1) = ra(i)
   }
//label_11:
   for(i = 0; i < n; i++) //do 12 i = 1, n
   {
      ra[i] = wk[iwk[i]][0]; //ra(i) = wk(iwk(i),1)
   }
//label_12:
   int j;
   for(j = 0; j < m; j++) //do 14 j = 1, m
   {
      for(i = 0; i < n; i++) //do 13 i = 1, n
      {
         wk[i][j] = rb[i][j]; //wk(i,j) = rb(i,j)
      }
   }
//label_13:
//label_14:
   for(j = 0; j < m; j++) //do 16 j = 1, m
   {
      for(i = 0; i < n; i++) //do 15 i = 1, n
      {
         rb[i][j] = wk[iwk[i]][j]; //rb(i,j) = wk(iwk(i),j)
      }
   }
//label_15:
//label_16:
   for(i = 0; i < 2000; i++) delete [] wk[i];
   delete [] wk;
   //END OF SUBROUTINE SORT
   //return
   //end
} /* end sort() */

/******************************************************************************
sort1()

SORTING SUBROUTINE ADAPTED FROM "NUMERICAL RECIPES"
BY W.H. PRESS ET AL., pp. 231

LIST OF VARIABLES
   ra(.) = integer array to be sorted
******************************************************************************/
void SCEUA::sort1
(
   int n,
   int * ra
)
{
   int i, j;
   int nn = n + 1;
   int * nra = new int[nn];
   for(i = 1; i <= n; i++)
   {
      nra[i] = ra[i-1];
   }

   //implicit real*8 (a-h,o-z)
   //dimension ra(n)

   int ir, rra; //integer ra, rra
   int l;
   
   l = (n / 2) + 1;
   ir = n;

label_10:
   if (l > 1) //if (l .gt. 1) then
   {
      l = l - 1;
      rra = nra[l]; //rra = ra(l)
   }
   else
   {
      rra = nra[ir]; //rra = ra(ir)
      nra[ir] = nra[1]; //ra(ir) = ra(1)
      ir = ir - 1;
      if (ir == 1) //if (ir .eq. 1) then
      {
         nra[1] = rra; //ra(1) = rra

         for(i = 1; i <= n; i++)
         {
            ra[i-1] = nra[i];
         }
         delete [] nra;

         return;
      } //end if
   } //end if
   i = l;
   j = l + l;

label_20:
   if (j <= ir) //if (j .le. ir) then
   {
      if (j < ir) //if (j .lt. ir) then
      {
         if (nra[j] < nra[j + 1]) //if (ra(j) .lt. ra(j + 1)) j = j + 1
         {
            j = j + 1;
         }
      } //end if
      if (rra < nra[j]) //if (rra .lt. ra(j)) then
      {
         nra[i] = nra[j]; //ra(i) = ra(j)
         i = j;
         j = j + j;
      }
      else
      {
         j = ir + 1;
      } //end if
      goto label_20;
   } //end if
   nra[i] = rra; //ra(i) = rra
   goto label_10;

   //END OF SUBROUTINE SORT1
   //end
} /* end sort1() */

/******************************************************************************
indexx()

THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
******************************************************************************/
void SCEUA::indexx
(
   int n, 
   double * arrin, 
   int * indx
)
{
   int i;
   double q;
   int nn = n + 1;
   double * narrin = new double[nn];
   int * nindx = new int[nn];
   for(i = 1; i <= n; i++)
   {
      narrin[i] = arrin[i-1];
      nindx[i] = i;
   }

   //implicit real*8 (a-h,o-z)
   //dimension arrin(n), indx(n)

   int j, l, ir, indxt;
   for(j = 1; j <= n; j++) //do 11 j = 1, n
   {
      nindx[j] = j; //indx(j) = j
   }
//label_11:
   l = (n / 2) + 1; //l = (n / 2) + 1
   ir = n; //ir = n
label_10:
   if (l > 1) //if (l .gt. 1) then
   {
      l = l - 1;
      indxt = nindx[l]; //indxt = indx(l)
      q = narrin[indxt]; //q = arrin(indxt)
   }
   else
   {
      indxt = nindx[ir];//indxt = indx(ir)
      q = narrin[indxt]; //q = arrin(indxt)
      nindx[ir] = nindx[1]; //indx(ir) = indx(1)
      ir = ir - 1;
      if (ir == 1) //if (ir .eq. 1) then
      {
         nindx[1] = indxt; //indx(1) = indxt

         for(i = 1; i <= n; i++)
         {
            indx[i-1] = nindx[i]-1;
         }
         delete [] nindx;
         delete [] narrin;
         return;
      } //end if
   }// end if
   i = l;
   j = l + l;
label_20:
   if (j <= ir)//if (j .le. ir) then
   {
      if (j < ir) //if (j .lt. ir) then
      {
         //if (arrin(indx(j)) .lt. arrin(indx(j + 1))) j = j + 1
         if (narrin[nindx[j]] < narrin[nindx[j + 1]])
         {
            j = j + 1;
         }
      } //end if
      if (q < narrin[nindx[j]]) //if (q .lt. arrin(indx(j))) then
      {
         nindx[i] = nindx[j]; //indx(i) = indx(j)
         i = j;
         j = j + j;
      }
      else
      {
         j = ir + 1;
      } //end if
      goto label_20;
   } //end if
   nindx[i] = indxt; //indx(i) = indxt
   goto label_10;

   //  END OF SUBROUTINE INDEXX
   //end
} /* end indexx() */

/******************************************************************************
ran1()

THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
******************************************************************************/
double SCEUA::ran1(int * idum)
{
   static int ix1, ix2, ix3, j;
   static double myran1;
   //implicit real*8 (a-h,o-z)
   static double r[97]; //dimension r(97)
   //parameter (m1 = 259200, ia1 = 7141, ic1 = 54773, rm1 = 3.8580247e-6)
   int m1 = 259200;
   int ia1 = 7141;
   int ic1 = 54773;
   double rm1 = 3.8580247e-6;

   //parameter (m2 = 134456, ia2 = 8121, ic2 = 28411, rm2 = 7.4373773e-6)
   int m2 = 134456; 
   int ia2 = 8121;
   int ic2 = 28411;
   double rm2 = 7.4373773e-6;

   //parameter (m3 = 243000, ia3 = 4561, ic3 = 51349)
   int m3 = 243000;
   int ia3 = 4561;
   int ic3 = 51349;

   //save
   //data iff / 0 /
   static int iff = 0;

   //if ((idum .lt. 0) .or. (iff .eq. 0)) then
   if ((*idum < 0) || (iff == 0))
   {
      iff = 1;
      ix1 = (ic1 - (*idum))%m1; //ix1 = mod(ic1 - idum,m1)
      ix1 = ((ia1 * ix1) + ic1)%m1; //ix1 = mod((ia1 * ix1) + ic1,m1)
      ix2 = ix1%m2; //ix2 = mod(ix1,m2)
      ix1 = ((ia1 * ix1) + ic1)%m1; //ix1 = mod((ia1 * ix1) + ic1,m1)
      ix3 = ix1%m3; //ix3 = mod(ix1,m3)
      for(j = 0; j < 97; j++) //do 11 j = 1, 97
      {
         ix1 = ((ia1 * ix1) + ic1)%m1; //ix1 = mod((ia1 * ix1) + ic1,m1)
         ix2 = ((ia2 * ix2) + ic2)%m2; //ix2 = mod((ia2 * ix2) + ic2,m2)
         //r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
         r[j] = ((double)ix1 + ((double)ix2 * rm2)) * rm1;
      }
//label_11:
      *idum = 1;
   } //end if
   ix1 = ((ia1 * ix1) + ic1)%m1; //ix1 = mod((ia1 * ix1) + ic1,m1)
   ix2 = ((ia2 * ix2) + ic2)%m2; //ix2 = mod((ia2 * ix2) + ic2,m2)
   ix3 = ((ia3 * ix3) + ic3)%m3; //ix3 = mod((ia3 * ix3) + ic3,m3)
   j = ((97 * ix3) / m3);
   //if ((j .gt. 97) .or. (j .lt. 1)) pause
   if ((j >= 97) || (j < 0)) 
   {
     j = 1+MyRand()%97; /* sleep(1000); */
   }
   myran1 = r[j]; //ran1 = r(j)
   //r(j) = (dble(ix1) + (dble(ix2) * rm2)) * rm1
   r[j] = ((double)ix1 + ((double)ix2 * rm2)) * rm1; 

   // END OF SUBROUTINE RAN1
   FILE * pOut = fopen("randoms.txt", "a");
   fprintf(pOut, "%14.7E\n", myran1);
   fclose(pOut);
   return myran1;
   //end
} /* end ran1() */ 

/******************************************************************************
gasdev()

THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
******************************************************************************/
double SCEUA::gasdev(int * idum)
{
   int icount = m_pModel->GetCounter();
   double mygasdev, v1, v2, r, fac;
   //implicit real*8 (a-h,o-z)
   //common /gasblk/ iset, gset
   static double gset = 0.00;
   static int iset = 0;
   //data iset / 0 /
   if (iset == 0) //if (iset .eq. 0) then
   {
label_1:
      v1 = (2.0 * ran1(idum)) - 1.0;
      v2 = (2.0 * ran1(idum)) - 1.0;
      //r = (v1 ** 2) + (v2 ** 2)
      r = (v1 * v1) + (v2 * v2);
      if (r >= 1.0) goto label_1; //if (r .ge. 1.) goto 1

      fac = sqrt(- ((2.0 * log(r)) / r));
      gset = v1 * fac;
      mygasdev = v2 * fac;
      iset = 1;
   }
   else
   {
      mygasdev = gset;
      iset = 0;
   } //end if

   //END OF SUBROUTINE GASDEV
   return mygasdev;
   //end
} /* end gasdev() */

/******************************************************************************
chkcst()


     This subroutine check if the trial point satisfies all
     constraints.

     ibound - violation indicator
            = -1 initial value
            = 0  no violation
            = 1  violation
     nopt = number of optimizing variables
     ii = the ii'th variable of the arrays x, bl, and bu
******************************************************************************/
void SCEUA::chkcst
(
   int nopt,
   double * x,
   double * bl,
   double * bu,
   int * ibound
)
{
   /*
   int i;
   *ibound = 0;
   for(i = 0; i < nopt; i++)
   {
      if(x[i] > bu[i]) (*ibound)++;
      if(x[i] < bl[i]) (*ibound)++;
   }
   */

   //implicit real*8 (a-h,o-z)
   //dimension x(nopt),bl(nopt),bu(nopt)

   *ibound = -1;

   //Check if explicit constraints are violated
   //do ii=1, nopt
   //  if (x(ii) .lt. bl(ii) .or. x(ii) .gt. bu(ii)) go to 10
   //end do
   //if (nopt .eq. 1) go to 9

   for(int ii=1; ii<=nopt; ii++)
   {
      if ((x[ii-1] < bl[ii-1]) || (x[ii-1] > bu[ii-1])) goto label10;
   }
   if (nopt == 1) goto label9;


//     Check if implicit constraints are violated
//     (no implicit constraints for this function)
//
//     No constraints are violated
//      
label9:    *ibound = 0;
      return;

//    At least one of the constraints are violated
label10:   *ibound = 1;
      return;
}/* end chkcst() */

/******************************************************************************
InitFromFile()

Read configuration information from the given filename, then write the 
configuration info. to the file "sce.in" (maintains compatibility with SCE
fortran implementation).
******************************************************************************/
void SCEUA::InitFromFile(IroncladString pFileName)
{
   FILE * pFile;
   int i;
   char * line;
   char tmp[DEF_STR_SZ], tmp2[DEF_STR_SZ];

   //assign defaults
   m_np = m_pModel->GetParamGroupPtr()->GetNumParams();
   m_Budget = 10000; //MAXN
   m_Kstop = 5; //KSTOP
   m_Pcento = 0.01; //PCENTO
   m_NumComplexes = (int)(sqrt((double)m_np)); //NGS
   m_Seed = 1969; //ISEED
   m_UserConfig = 1; //IDEFLT
   m_PtsPerComplex = 2*m_np + 1; //NPG
   m_PtsPerSubComplex = m_np+1; //NPS
   m_NumEvoSteps = m_PtsPerComplex; //NSPL
   m_MinComplexes = m_NumComplexes; //MINGS
   m_UseInitPt = 1; //INIFLG
   m_OutputMode = 2; //IPRINT
   m_bUseInitPt = false; //INIFLG

   //allocate initial parameter configuration
   ParameterGroup * pGroup = m_pModel->GetParamGroupPtr();
   m_np = pGroup->GetNumParams();
   NEW_PRINT("double", m_np);
   m_pParams = new double[m_np];
   MEM_CHECK(m_pParams);

   NEW_PRINT("double", m_np);
   m_pLower = new double[m_np];
   MEM_CHECK(m_pParams);

   NEW_PRINT("double", m_np);
   m_pUpper = new double[m_np];
   MEM_CHECK(m_pParams);

   for(i = 0; i < m_np; i++)
   {
      m_pParams[i] = pGroup->GetParamPtr(i)->GetEstVal();
      m_pLower[i] = pGroup->GetParamPtr(i)->GetLwrBnd();
      m_pUpper[i] = pGroup->GetParamPtr(i)->GetUprBnd();
   }/* end for() */

   //read in SCEUA configuration
   pFile = fopen(pFileName, "r");
   if(pFile == NULL) 
   {
      //couldn't open file, use defaults and log the error.
      LogError(ERR_FILE_IO, "Couldn't open SCEUA config. file. Using Defaults");      
      return;
   }/* end if() */   

   if(CheckToken(pFile, "RandomSeed", GetInFileName()) == true)
   {
      fclose(pFile);
      m_Seed = GetRandomSeed();
      pFile = fopen(pFileName, "r");
   }
   rewind(pFile);

   //make sure correct tokens are present
   if(CheckToken(pFile, "BeginSCEUA", pFileName) == true)
   {
      FindToken(pFile, "EndSCEUA", pFileName);
      rewind(pFile);

      FindToken(pFile, "BeginSCEUA", pFileName);
      line = GetNxtDataLine(pFile, pFileName);
      while(strstr(line, "EndSCEUA") == NULL)
      {         
         if(strstr(line, "Budget") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_Budget); 
            if(m_Budget < 100)
            {
               LogError(ERR_FILE_IO, "Invalid SCEUA budget. Defaulting to 100.");
               m_Budget = 100;
            }
         }/*end if() */
         else if(strstr(line, "LoopStagnationCriteria") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_Kstop); 
         }
         else if(strstr(line, "PctChangeCriteria") != NULL)
         {
            sscanf(line, "%s %lf", tmp, &m_Pcento); 
         }
         else if(strstr(line, "NumComplexes") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumComplexes); 
         }
         else if(strstr(line, "NumPointsPerComplex") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_PtsPerComplex); 
         }
         else if(strstr(line, "NumPointsPerSubComplex") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_PtsPerSubComplex); 
         }
         else if(strstr(line, "NumEvolutionSteps") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_NumEvoSteps); 
         }
         else if(strstr(line, "MinNumOfComplexes") != NULL)
         {
            sscanf(line, "%s %d", tmp, &m_MinComplexes); 
         }
         else if(strstr(line, "UseInitialPoint") != NULL)
         {
            sscanf(line, "%s %s", tmp, tmp2); 
            MyStrLwr(tmp2);
            if(strcmp(tmp2, "yes") == 0) m_bUseInitPt = true;
         }

         else
         {
            sprintf(tmp, "Unknown token: %s", line);
            LogError(ERR_FILE_IO, tmp);
         }/* end else() */
         line = GetNxtDataLine(pFile, pFileName);
      } /* end while() */
   }/* end if() */   
   fclose(pFile);

   //create sce.in file
   int iniflg = 0;
   if(m_bUseInitPt) iniflg = 1;
   FILE * pOut = fopen("sce.in", "w");
   fprintf(pOut, "%d  %d  %lf  %d  %d  1\n", 
           m_Budget, m_Kstop, m_Pcento, m_NumComplexes, m_Seed);
   fprintf(pOut, "%d  %d  %d  %d  %d  2\n", 
           m_PtsPerComplex, m_PtsPerSubComplex, m_NumEvoSteps, m_MinComplexes, iniflg);
   for(i = 0; i < m_np; i++)
   {
      fprintf(pOut, "%E %E %E\n", m_pParams[i], m_pLower[i], m_pUpper[i]);
   }
   fclose(pOut);
} /* end InitFromFile() */

/******************************************************************************
WriteMetrics()

Write out algorithm metrics and setup.
******************************************************************************/
void SCEUA::WriteMetrics(FILE * pFile) 
{
   fprintf(pFile, "\nAlgorithm Metrics\n");
   fprintf(pFile, "Algorithm                : Shuffled Complex Evolution (SCE)\n");
   fprintf(pFile, "Budget                   : %d\n", m_Budget);
   fprintf(pFile, "Loop Stagnation Criteria : %d\n", m_Kstop); 
   fprintf(pFile, "Pct Change Criteria      : %lf\n", m_Pcento); 
   fprintf(pFile, "Number of Complexes      : %d\n", m_NumComplexes); 
   fprintf(pFile, "Points Per Complex       : %d\n", m_PtsPerComplex); 
   fprintf(pFile, "Points Per Sub-Complex   : %d\n", m_PtsPerSubComplex); 
   fprintf(pFile, "Num. of Evolution Steps  : %d\n", m_NumEvoSteps); 
   fprintf(pFile, "Min. Num. of Complexes   : %d\n", m_MinComplexes); 
  
   m_pModel->WriteMetrics(pFile);
}/* end WriteMetrics() */

/******************************************************************************
SCEUA_Program()

Calibrate or optimize the model using SCE.
******************************************************************************/
void SCEUA_Program(int argC, StringType argV[])
{
   NEW_PRINT("Model", 1);
   ModelABC * model = new Model();

   NEW_PRINT("SCEUA", 1);
   SCEUA * SCE = new SCEUA(model);
   MEM_CHECK(SCE);

   if(model->GetObjFuncId() == OBJ_FUNC_WSSE) { SCE->Calibrate(); }
   else { SCE->Optimize(); }

   delete SCE;
   model->Destroy();
} /* end SCEUA_Program() */
