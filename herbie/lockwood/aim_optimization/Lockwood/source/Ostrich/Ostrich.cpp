/******************************************************************************
File     : Ostrich.cpp
Author   : L. Shawn Matott
Copyright: 2003, L. Shawn Matott

Main program execution. Provides a text interface for the set of optimization
and gridding algorithms that make up the Ostrich program.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added PSO
03-24-04    lsm   added PSO-LevMar hybrid, added ISOFIT_BUILD option
11-07-05    lsm   added support for BGA, GRID, VSA and CSA programs
03-03-07    jrc   added DDS program
******************************************************************************/
#include "mpi_stub.h"
#include <stdio.h>
#include <string.h>
#include "Exception.h"
#include "Utility.h"
#include "OptMathClass.h"
#include "OptSearchClass.h"
#include "SAAlgorithm.h"
#include "VandSA.h"
#include "CCSA.h"
#include "CCVSA.h"
#include "SCEUA.h"
#include "ComboSA.h"
#include "BisectionAlgorithm.h"
#include "PowellAlgorithm.h"
#include "SamplingAlgorithm.h"
#include "SteepDescAlgorithm.h"
#include "FletchReevesAlgorithm.h"
#include "LevenbergAlgorithm.h"
#include "GridAlgorithm.h"
#include "StatsClass.h"
#include "GeneticAlgorithm.h"
#include "BinaryGA.h"
#include "ParticleSwarm.h"
#include "ParticleSwarmDESC.h"
#include "DDSAlgorithm.h"/*JRC*/
#include "DiscreteDDSAlgorithm.h"
#include "MyPBS.h"
#include "StatUtility.h"
#include "QuadTree.h"
#include "IsoParse.h"
#include "CCPSO.h"
#include "CCRGA.h"
#include "GLUE.h"
#include "RejectionSampler.h"

#ifdef ISOFIT_BUILD
int Ostrich(int argc, StringType argv[])
#else
int main(int argc, StringType argv[])
#endif
{
   ProgramType program;
   UnmoveableString pInFile = GetOstFileName();

   //initialize time tracker
   GetElapsedTime();

   MPI_Init(&argc,&argv);

#ifndef ISOFIT_BUILD
   InitErrors();   
#endif

   //initialize input files (assume only one input file)
   strcpy(GetOstFileName(), "ostIn.txt");
   strcpy(GetExeDirName(), ".");
   InitDataLine("ostIn.txt");
   
   program = ReadProgramType();
   SetProgramType(program);

   //execute desired operation
   switch(program)
   {
      case(GA_PROGRAM) : 
      {
         GA_Program(argc, argv);
         break;
      }/* end case(GA_PROGRAM) */
      case(BGA_PROGRAM) : 
      {
         BGA_Program(argc, argv);
         break;
      }/* end case(BGA_PROGRAM) */
      case(GRID_PROGRAM) : 
      {
         GRID_Program(argc, argv);
         break;
      }/* end case(GRID_PROGRAM) */
      case(SA_PROGRAM) : 
      {
         SA_Program(argc, argv);
         break;
      }/* end case(SA_PROGRAM) */
      case(CSA_PROGRAM) : //combinatorial simulated annealing
      {
         CSA_Program(argc, argv);
         break;
      }/* end case(CSA_PROGRAM) */
      case(CCSA_PROGRAM) : //computation-constrained simulated annealing
      {
         CCSA_Program(argc, argv);
         break;
      }/* end case(CCSA_PROGRAM) */
      case(CCVSA_PROGRAM) : //computation-constrained Vanderbilt-Louie simulated annealing
      {
         CCVSA_Program(argc, argv);
         break;
      }/* end case(CCVSA_PROGRAM) */
      case(VSA_PROGRAM) : //vanderbilt-louie simulated annealing
      {
         VSA_Program(argc, argv);
         break;
      }/* end case(VCSA_PROGRAM) */
      case(PSO_PROGRAM) : 
      {
         PSO_Program(argc, argv);
         break;
      }/* end case(PSO_PROGRAM) */
      case(PSODESC_PROGRAM) : 
      {
         PSODESC_Program(argc, argv);
         break;
      }/* end case(PSODESC_PROGRAM) */
      case(PSO_LEV_PROGRAM) : 
      {
         PSO_LEVMAR_Program(argc, argv);
         break;
      }/* end case(PSO_LEV_PROGRAM) */
      case(PSODESC_LEV_PROGRAM) : 
      {
         PSODESC_LEVMAR_Program(argc, argv);
         break;
      }/* end case(PSODESC_LEV_PROGRAM) */
      case(CCPSO_PROGRAM) : 
      {
         CCPSO_Program(argc, argv);
         break;
      }/* end case(CCPSO_PROGRAM) */
      case(SCEUA_PROGRAM) : 
      {
         SCEUA_Program(argc, argv);
         break;
      }/* end case(SCEUA_PROGRAM) */
      case(CCRGA_PROGRAM) : 
      {
         CCRGA_Program(argc, argv);
         break;
      }/* end case(CCRGA_PROGRAM) */
      case(CCPSO_LEV_PROGRAM) : 
      {
         CCPSO_LEVMAR_Program(argc, argv);
         break;
      }/* end case(CCPSO_LEV_PROGRAM) */
      case(LEV_PROGRAM) : 
      {
         LEV_Program(argc, argv);
         break;
      }/* end case(LEV_PROGRAM) */
      case(GMLMS_PROGRAM) : 
      {
         GMLMS_Program(argc, argv);
         break;
      }/* end case(GMLMS_PROGRAM) */
      case(POWL_PROGRAM) : 
      {
         PWL_Program(argc, argv);
         break;
      }/* end case(POWL_PROGRAM) */
      case(STEEP_PROGRAM) : 
      {
         STPDSC_Program(argc, argv);
         break;
      }/* end case(STEEP_PROGRAM) */
      case(FLRV_PROGRAM) : 
      {
         FLRV_Program(argc, argv);
         break;
      }/* end case(FLRV_PROGRAM) */
      case(BIS_PROGRAM) : 
      {
         BIS_Program(argc, argv);
         break;
      }/* end case(BIS_PROGRAM) */
      case(SMP_PROGRAM) : 
      {
         SMP_Program(argc, argv);
         break;
      }/* end case(SMP_PROGRAM) */
      case(STATS_PROGRAM) : 
      {
         STATS_Program(argc, argv);
         break;
      }/* end case(STATS_PROGRAM) */
      case(EVAL_PROGRAM) : 
      {
         EVAL_Program(argc, argv);
         break;
      }/* end case(EVAL_PROGRAM) */
      case(UTIL_PROGRAM) : 
      {
		 ConvertToASCII();
         //STATS_TestFdist();
         //STATS_TestStudentDist();
         //STATS_TestStdNormDist();
         break;
      }/* end case(UTIL_PROGRAM) */
      case(DDS_PROGRAM) :/*JRC*/
      {
         DDS_Program(argc,argv);
			break;
      }
      case(DDDS_PROGRAM) :
      {
         DiscreteDDS_Program(argc,argv);
			break;
      }
      case(GLUE_PROGRAM) :
      {
         GLUE_Program(argc,argv);
			break;
      }
      case(RJSMP_PROGRAM) :
      {
         RJSMP_Program(argc,argv);
			break;
      }
      case(METRO_PROGRAM) :
      {
         METRO_Program(argc,argv); //Metropolis MCMC
			break;
      }
      case(QUIT_PROGRAM) : 
      default:
      {            
         break;
      }/* end case(QUIT_PROGRAM) */      
   }/* end switch() */   

   MPI_Finalize();

   #ifndef ISOFIT_BUILD
      ExitProgram(0);
   #endif

   return (0);
} /* end main() */
