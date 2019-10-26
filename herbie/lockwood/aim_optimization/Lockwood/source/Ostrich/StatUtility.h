/******************************************************************************
File     : StatUtility.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

StatUtility contains C routines that assist in statistical calculations.

Version History
06-19-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   replaced magic # with NEARLY_ZERO
01-01-07    lsm   Added Durbin-Watson and runs tests.
******************************************************************************/
#ifndef STAT_UTILITY_H
#define STAT_UTILITY_H

#include "MyTypes.h"

double GammaLn(double val);
double CalcMean(Ironclad1DArray  v, int size);
double CalcMedian(Unmoveable1DArray v, int size);
double CalcStdDev(Ironclad1DArray v, int size);

double FdistCDF(int u, int v, double xLwr, double xUpr);
double FdistPDF(int u, int v, double x);
double FdistInvCDF(int u, int v, double p);

double StudentCDF(int dof, double x);
double StudentPDF(int dof, double x);
double StudentInvCDF(int dof, double p);

double StdNormCDF(double x);
double StdNormPDF(double x);
double StdNormInvCDF(double p);

bool DurbinWatsonTest(double * res, double ** jac, int n, int k, double * D, 
                      double * plwr, double * pupr);
double ImhofSum(double d, double * wr, int n);
double Imhof(double u);
bool RunsTest(double * residuals, int num, int * nPos, int * nNeg, int * nRuns, 
              double * plwr, double * pupr);
double nCr(int n, int r);

void STATS_TestFdist(void);
void STATS_TestStudentDist(void);
void STATS_TestStdNormDist(void);
void STATS_TestAutocorrelation(void);

#endif /* STAT_UTILITY_H */

