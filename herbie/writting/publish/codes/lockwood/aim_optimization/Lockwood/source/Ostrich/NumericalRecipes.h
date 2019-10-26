/******************************************************************************
File     : NumericalRecipes.cpp
Authors  : L. Shawn Matott, William H. Press, Saul A. Teukolsky,
           William T. Vetterling, and Brian P. Flannery
Copyright: 2004, L. Shawn Matott and 1986-92, Numerical Recipes Software o!7.]41y.

Macros and routines that are adapted from Numerical Recipes in C.

Version History
05-07-04    lsm   created file
01-01-07    lsm   Added helper functions for new statistics (Runs Test and 
                  Durbin-Watson)
******************************************************************************/
#ifndef NUMERICAL_RECIPES_H
#define NUMERICAL_RECIPES_H

extern "C" {
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SWAP(g,h)     {y=(g);(g)=(h);(h)=y;}

int    IMAX(int a, int b);
double FMAX(double a, double b);
double FMIN(double a, double b);
double SIGN(double a, double b);
double SQR(double a);
void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n, double * wr, double * wi);

double qromo(double (*func)(double), double a, double b,
             double (*choose)(double(*)(double), double, double, int));

double midpnt(double (*func)(double), double a, double b, int n);
double midinf(double (*funk)(double), double aa, double bb, int n);

void polint(double * xa, double * ya, int n, double x, double *y, double *dy);
void nrerror(char * error_text);
bool ludcmp(double **a, int n, int *indx, double *d);
}

#endif /* NUMERICAL_RECIPES_H */

