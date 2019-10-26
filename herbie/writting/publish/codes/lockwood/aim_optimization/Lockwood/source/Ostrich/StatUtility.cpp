/******************************************************************************
File     : StatUtility.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

StatUtility contains C routines that assist in statistical calculations.

Version History
06-19-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   replaced magic # with NEARLY_ZERO
01-01-07    lsm   Added Durbin-Watson and runs tests.
******************************************************************************/
#include "StatUtility.h"
#include "Utility.h"
#include <math.h>
#include "MyErf.h"
#include "NumericalRecipes.h"

/*------------------------------------------------------------------
Imhof globals, these are passed between ImhofSum() and Imhof() via 
globals so that Imhof() is compatible with the Numerical Recipes 
calling convention.
------------------------------------------------------------------*/
int gImhofN; 
double * gImhofTau;
double gImhofR0;

/******************************************************************************
GammaLn()

Compute natural log of the gamma function. Code adapted from Numerical 
Recipes in C++, 2nd Edition, page 219.
******************************************************************************/
double GammaLn(double val)
{
   int j;
   double x, y, tmp, ser;
   
   static double coeff[6] = 
   { 76.18009172947146, -86.50532032941677, 24.01409824083091, 
    -1.231739572450155, 0.001208650973866179, -0.5395239384953e-5};

   y = x = val;
   tmp = x + 5.5;
   tmp -= (x + 0.5) * log(tmp);
   ser = 1.000000000190015;

   for(j = 0; j < 6; j++)
   {
      y++;
      ser += coeff[j] / y;
   }/* end for() */       

   return (-tmp + log(2.5066282746310005 * ser / x));
}/* end GammaLn() */

/******************************************************************************
CalcStdDev()

Returns the standard deviation of a list of numbers.
******************************************************************************/
double CalcStdDev(Ironclad1DArray v, int size)
{
   int i;
   double mean;
   double sum;
   double v1;

   mean = CalcMean(v, size);

   sum = 0.00;

   for(i = 0; i < size; i++)
   {
      v1 = v[i];
      sum += ((v1 - mean) * (v1 - mean));
   } /* end for() */
  
   return sqrt(sum / (double)size);
} /* end CalcStdDev() */

/******************************************************************************
CalcMean()

Returns the mean of a list of numbers
******************************************************************************/
double CalcMean(Ironclad1DArray v, int size)
{
  int i;
  double sum;

  sum = 0.00;
  for(i = 0; i < size; i++){ sum += v[i];}
  return (sum / (double)size);
} /* end CalcMean() */

/******************************************************************************
CalcMedian()

Returns the median of a list of numbers
******************************************************************************/
double CalcMedian(Unmoveable1DArray v, int size)
{
  int i = size/2;

  SortInc(v,size);

  //even
  if((size % 2) == 0){ return 0.5*(v[i]+v[i-1]); }
  //odd
  return v[i];
} /* end CalcMedian() */

/******************************************************************************
FdistPDF()

Calculates the probability density function of the F-distribution. 

Code is based on formula presented in "Applied Statistics and Probability 
for Engineers", Montgomery and Runger, 1994. Page 316, Theorem 6-5.
******************************************************************************/
double FdistPDF(int u, int v, double x)
{
   double u1 = (double)u;
   double v1 = (double)v;
   double val;   

   double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

   if(x < 0) { return 0;}   
   
   tmp1 = GammaLn((u1 + v1) / 2.00);
   tmp2 = GammaLn(u1 / 2.00);
   tmp3 = GammaLn(v1 / 2.00);
   tmp4 = exp(tmp1 - (tmp2 + tmp3));
   tmp5 = pow((u1 / v1), (u1 / 2.00));
   tmp6 = tmp4 * tmp5 * pow(x, ((u1 / 2.00) - 1.00));
   tmp7 = ((u1 / v1) * x) + 1.00;
   tmp7 = pow(tmp7, ((u1 + v1) / 2.00));

   val = (tmp6 / tmp7);

   return val;
}/* end FdistPDF() */

/******************************************************************************
FdistCDF()

Calculates the cumulative density function of the F-distribution. The args are
u = d.o.f. in numerator
v = d.o.f. in denominator
xLwr = lower limit
xUpr = upper limit

Returns the probability that a F-distributed variable with u and v degress of 
freedom is less than or equal to xUpr and greter than or equal to xLwr 
(i.e. P[xLwr <= X <= xUpr]).
******************************************************************************/
double FdistCDF(int u, int v, double xLwr, double xUpr)
{
   int j;
   double lwr, upr, mid; //integration limits
   double dt;
   double sum, old;
   double Flwr, Fupr, Fmid;
   double stop_val = 1E-6;
      
   /*--------------------------------------
   determine integration limits.
   ---------------------------------------*/
   if(xLwr <= 0.00){xLwr = 0.00;}
   if(xUpr <= 0.00){return 0.00;}
   if(xLwr == xUpr){return 0.00;}
   if((xLwr == 0.00) && (xUpr == 1.00) && (u == v)){return 0.50;}

   if(u == 1)
   { 
      Fupr = (2.00*StudentCDF(v, sqrt(xUpr)) - 1.00);
      Flwr = (2.00*StudentCDF(v, sqrt(xLwr)) - 1.00);

      return (Fupr - Flwr);
   }/* end if() */
   
   lwr = xLwr; 
   upr = xUpr;

   /*-----------------------------------------------------
   integrate over the limits, the integration technique
   is an iterative form of the trapezoidal rule. The 
   number of intervals is doubled each iteration until 
   the area calculation has converged.
   ------------------------------------------------------*/
   Flwr = FdistPDF(u, v, lwr);
   Fupr = FdistPDF(u, v, upr);
   dt = (upr - lwr);
   //initially apply trapezoidal rule at the end points
   sum = 0.5*dt*(Flwr+Fupr); 
   old = sum;
   j = 0;   
   //iterate until converged or at least 5 times
   while((fabs(sum - old) > stop_val) || (j < 5))
   {
      j++;      
      dt /= 2; //halve step size
      mid = lwr + dt; //initial value of points
   
      old = sum;
      sum = 0.00;
      //calculate F() of points
      while(mid <= upr)
      {
         Fmid = FdistPDF(u, v, mid);
         sum +=  Fmid;
         //a step of 2*dt skips previously computed points    
         mid += (2.00*dt); 
      }/* end while() */
      sum *= dt;      
      sum += (old * 0.5);
   }/* end while() */   
      
   return (sum);
}/* end FdistCDF() */

/******************************************************************************
FdistInvCDF()

Calculates the F-distribution upper-tail percentage point, (i.e. given  a 
probability p, compute the value of x such that the P[X < x] = p).
******************************************************************************/
double FdistInvCDF(int u, int v, double p)
{
   int j;
   double x, upr, lwr;
   double F;
   double stopVal = 1E-6;

   if(p <= 0.00){ return 0.00;}
   if(p >= 1.00){ p = 1.00;} //max. probability is 1

   //initial guess is at center of distribution
   x = 1.00;
   F = FdistCDF(u, v, 0.00, x);

   upr = 1.00;
   lwr = 0.00;
   
   //find upper limit
   while(F < p) 
   { 
      lwr = x;
      x *= 2;
      //add in the area between lwr and x
      F += FdistCDF(u, v, lwr, x); 
   }/* end while() */   

   //converge on answer
   j = 0;
   while((fabs(F - p) > stopVal) || (j < 5))
   {
      if(F >= p) 
      { 
         upr = x;
         x = 0.5*(upr + lwr);
         //subtract the area between x and upr
         F -= FdistCDF(u, v, x, upr);
      }/* end if() */
      else 
      { 
         lwr = x;
         x = 0.5*(upr + lwr);
         //add the area between lwr and x
         F += FdistCDF(u, v, lwr, x);
      }/* end else() */

      j++;
   }/* end while() */
   
   return x;
}/* end FdistInvCDF() */

/******************************************************************************
StudentPDF()

Calculates the probability density function of the Student's t-distribution, 
given dof (degrees of freedom) and x.
******************************************************************************/
double StudentPDF(int dof, double x)
{
   double coeff; //coefficient, only depends on dof
   double e; //exponent, only depends on dof
   double tmp1, tmp2, g1, g2;
   double val;

   /*-----------------------------------
   Calculate terms that rely only on dof
   -----------------------------------*/
   tmp1 = 0.5*((double)dof + 1.00);
   e = -tmp1; 
   tmp2 = 0.5*(double)dof;
   g1 = GammaLn(tmp1);
   g2 = GammaLn(tmp2);
   tmp1 = exp(g1 - g2);
   coeff = tmp1 / sqrt((double)dof*MY_PI);

   /*-----------------------------------
   calculate the x-dependent term
   ------------------------------------*/
   tmp1  = (1.00 + ((x*x)/(double)dof));
   tmp2 = pow(tmp1, e);
   val = (coeff*tmp2);
      
   return (val);
}/* end StudentPDF() */

/******************************************************************************
StudentCDF()

Calculates the Student's cumulative t-distribution, given dof (degrees of 
freedom) and x.

Returns the probability that a t-distributed variable with dof degress of 
freedom is less than or equal to x (i.e. P[X <= x]).
******************************************************************************/
double StudentCDF(int dof, double x)
{
   double P; //desired probability
   double lwr, upr, mid; //integration limits
   double dt;
   double sum, old;
   double Flwr, Fupr, Fmid;   
   double stop_val = 1E-6;
   int j;

   //area under PDF from -infinity to 0 is 0.5
   P = 0.5; 

   /*--------------------------------------
   determine integration limits.
   ---------------------------------------*/
   if(x == 0.00) { return P;}
   else if(x > 0.00){lwr = 0.00; upr = x;}
   else {lwr = x; upr = 0.00;}

   /*-----------------------------------------------------
   integrate over the limits, the integration technique
   is an iterative form of the trapezoidal rule. The 
   number of intervals is doubled each iteration until 
   the area calculation has converged.
   ------------------------------------------------------*/
   Flwr = StudentPDF(dof, lwr);
   Fupr = StudentPDF(dof, upr);
   dt = (upr - lwr);
   //initially apply trapezoidal rule at the end points
   sum = 0.5*dt*(Flwr+Fupr); 
   old = sum;
   j = 0;   
   //iterate until converged or at least 5 times
   while((fabs(sum - old) > stop_val) || (j < 5))
   {
      j++;
      dt /= 2; //halve step size
      mid = lwr + dt; //initial value of points
   
      old = sum;   
      sum = 0.00;
      //calculate F() of points
      while(mid <= upr)
      {
         Fmid = StudentPDF(dof, mid);
         sum +=  Fmid;
         //a step of 2*dt skips previously computed points    
         mid += (2.00*dt); 
      }/* end while() */
      sum *= dt;
      sum += (old * 0.5);
   }/* end while() */

   if(x > 0.00) { P += sum;}
   else { P -= sum;}
      
   return (P);
}/* end StudentCDF() */

/******************************************************************************
StudentInvCDF()

Calculates the student-distribution upper-tail percentage point, (i.e. given  a 
probability p, compute the value of x such that P[X < x] = p).
******************************************************************************/
double StudentInvCDF(int dof, double p)
{
   double x, upr, lwr;
   double F;
   double stopVal = 1E-6;
   bool flip = false;

   if(p <= 0.00){ p = 0.00;} //min. probability is 0
   if(p >= 1.00){ p = 1.00;} //max. probability is 1
   if(p == 0.50){ return 0.00;}   
   if(p < 0.50){ flip = true; p = 1 - p;}

   //initial guess is at center of distribution
   x = 0.00; F = 0.50;

   upr = 1.00; lwr = 0.00;
   
   //find upper limit
   while(F < p) 
   { 
      lwr = x;
      x = upr;
      upr *= 2;
      F = StudentCDF(dof, x);
   }/* end while() */   

   //converge on answer
   while(fabs(F - p) > stopVal)
   {      
      if(F >= p) { upr = x;}   
      else { lwr = x;}

      x = 0.5*(upr + lwr);
      F = StudentCDF(dof, x);
   }/* end while() */
   
   if(flip == true){ return -x;}
   return x;
}/* end StudentInvCDF() */

/******************************************************************************
StdNormPDF()

Calculates the probability density function of the standard normal distribution.
This is a normal distribution with mean = 0 and std. dev = 1.
******************************************************************************/
double StdNormPDF(double x)
{
   double val;
   
   val = (exp(-(x*x)*0.5)/(sqrt(2.00 * MY_PI)));
   
   return (val);
}/* end StdNormPDF() */

/******************************************************************************
StdNormCDF()

Calculates the cumulative density function of the standard normal distribution.
This is a normal distribution with mean = 0 and std. dev = 1.

Returns the probability that a normally variable is less than or equal to x 
(i.e. P[X <= x]).

NOTE: Utilizes the relationship of the std. normal distribution to the error
function.
******************************************************************************/
double StdNormCDF(double x)
{
   double val;
   
   val = 0.5*(MyErf(x/sqrt(2.00)) + 1.00);
   
   return (val);
}/* end StdNormCDF() */

/******************************************************************************
StdNormInvCDF()

Calculates the standard normal distribution upper-tail percentage point, 
(i.e. given  a probability p, compute the value of x such that P[X < x] = p).
******************************************************************************/
double StdNormInvCDF(double p)
{
   double x, upr, lwr;
   double F;
   double stopVal = NEARLY_ZERO;
   bool flip = false;

   if(p <= 0.00){ p = 0.00;} //min. probability is 0
   if(p >= 1.00){ p = 1.00;} //max. probability is 1
   if(p == 0.50){ return 0.00;}
   if(p < 0.50){ flip = true; p = 1 - p;}

   //initial guess is at center of distribution
   x = 0.00; F = 0.50;
   lwr = 0.00; upr = 1.00;
   
   //find upper limit
   while(F < p) 
   { 
      lwr = x;
      x = upr;
      upr *= 2;

      F = StdNormCDF(x);
      if(fabs(F - p) < stopVal){ break;}
   }/* end while() */   

   //converge on answer
   while(fabs(F - p) > stopVal)
   {      
      if(F > p) { upr = x;}   
      else { lwr = x;}

      x = 0.5*(upr + lwr);
      F = StdNormCDF(x);
   }/* end while() */

   if(flip == true){ return -x;}   
   return x;
}/* end StdNormInvCDF() */

/******************************************************************************
DurbinWatsonTest()

Perform a Durbin-Watson test for autocorrelation of residuals. The residuals 
should be placed in an appropriate order by the calling routine (e.g. ordered by
increasing value of observation).

Computation is based on material presented in:

   White KJ. 1992. "The Durbin-Watson Test for Autocorrelation in 
   Nonlinear Models", The Review of Economics and Statistics, vol. 74, no. 
   2, pg. 370-373.

   and

   Durbin J. Watson GS. 1971. "Testing for Serial COrrelation in Least Squares
   Regression III", Biometrika, vol. 58, no. 1, pg. 1-19.

Input:
   res : Residuals vector
   jac : Jacobian matrix
   n   : number of observations (residuals)
   k   : number of parameters
   
Output:
   Returns true if successful (n > 2) and sets the following output values:
      D    : Durbin-Watson statistic (between 0 and 4)
      plwr : p-value of the lower tail Prob(D < d)
      pupr : p-value of the upper tail Prob(D > d)
   Return false if n is less than or equal to 2 (too few samples)
******************************************************************************/
bool DurbinWatsonTest(double * res, double ** jac, int n, int k, double * D, 
                      double * plwr, double * pupr)
{
   int i, j;
   double numer, denom;
   double ** A, ** M, ** MA, ** JT, ** TMP, ** INV;
   double * wr, * wi;

   *plwr = *pupr = 0.50;
   
   if(n <= 2){ *D = 2.00; return false;}

   numer = denom = 0.00;
   for(i = 0; i < n; i++){ denom += (res[i] * res[i]);}
   for(i = 1; i < n; i++){ numer += (res[i] - res[i-1])*(res[i] - res[i-1]);}
   *D = (numer / denom);

   if((jac == NULL) || (n < 4)){ return true;}   

   //allocate and initialize the A matrix (see White. 1992)
   A   = new double * [n];
   M   = new double * [n];
   MA  = new double * [n];
   JT  = new double * [k];
   INV  = new double * [k];
   TMP = new double * [n];
   wr = new double[n];
   wi = new double[n];
   for(i = 0; i < n; i++)
   { 
      A[i] = new double[n];
      M[i] = new double[n];
      MA[i] = new double[n];
      if(i < k){ 
         JT[i] = new double[n];
         INV[i] = new double[k];}
      TMP[i] = new double[n];
   }
   for(i = 0; i < n; i++){ for(j = 0; j < n; j++){ 
      A[i][j] = 0.00e00; 
      if(i < k) JT[i][j] = jac[j][i]; }}
   for(i = 0; i < n; i++)
   { 
      A[i][i] = 2.00;
      if(i < n - 1) A[i][i+1] = A[i+1][i] = -1.00e00;
   }
   A[0][0] = A[n-1][n-1] = 1.00;

   //Compute M = I - J[J^T J]^-1 J^T  and then MA = M * A (see White. 1992)
   MatMult(JT, jac, INV, k, n, k);
   if(MatInv(INV, INV, k) == false)
   {
      return false;
   }
   MatMult(jac, INV, TMP, n, k, k);
   MatMult(TMP, JT, M, n, k, n);
   for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
         if(i == j) M[i][j] = 1.00 - M[i][j];
         else       M[i][j] = - M[i][j]; }}
   MatMult(M, A, MA, n, n, n);

    /*-------------------------------------------
   compute the eigenvalues of MA (see White. 1992) 
   using canned numerical recipes in C routines
   --------------------------------------------*/
   balanc(MA, n);
   elmhes(MA, n);
   hqr(MA, n, wr, wi);

   //integrate the CDF [see equation (1-11) and (1-12) of Durbin-Watson. 1971]
   *pupr = ImhofSum(*D, wr, n-k);
   if(*pupr < 0.00) *pupr = 0.00;
   if(*pupr > 1.00) *pupr = 1.00;
   *plwr = 1.00 - *pupr;
   
   for(i = 0; i < n; i++)
   { 
      delete A[i];
      delete M[i];
      delete MA[i];
      delete TMP[i];
      if(i < k){
         delete JT[i];
         delete INV[i]; }
   }
   delete [] A;
   delete [] M;
   delete [] MA;
   delete [] TMP;
   delete [] JT;
   delete [] INV;
   delete [] wr;
   delete [] wi;

   return true;
}/* end DurbinWatsonTest() */

/******************************************************************************
ImhofSum()

Compute the probability that D > d by integration of equation 1-11 of 
Durbin-Watson (1971).

Input:
   d   : the Durbin-Watson test statistic
   wr  : the vector of eigenvalues
   n   : number of eigenvalues
******************************************************************************/
double ImhofSum(double d, double * wr, int n)
{
   double p, a, b;

   //'pass' globals to Imhof() prior to invoking NR integration routines
   gImhofN = n;
   gImhofTau = wr;
   gImhofR0 = d;

   /*---------------------------------------------
   First integrate from a = 0 to b = 2 * PI
   ----------------------------------------------*/
   a = qromo(Imhof, 0.00, 2.00 * MY_PI, midpnt);

   /*----------------------------------------------
   Now integrate from a = 2 * Pi to b = infinity 
   using change of variables suggested by NR in C : 
   a' = 0 (1/infinity) and b' = 1/a
   ----------------------------------------------*/
   b = qromo(Imhof, 2.00 * MY_PI, NEARLY_HUGE, midinf);

   p = 0.5 + (a + b)/MY_PI;

   return p;
}/* end ImhofSum() */

/******************************************************************************
Imhof()

Compute the integration function given in (1-11) of Durbin and Watson (1971).
******************************************************************************/
double Imhof(double u)
{
   int i;
   double theta = 0.00, rho = 1.00, f;
   int n = gImhofN;
   double * tau = gImhofTau;
   double r0 = gImhofR0;

   for(i = 0; i < n; i++)
   {
      theta += atan(u * (tau[i] - r0));
      rho   *= pow((1.00 + ((tau[i]-r0)*(tau[i]-r0))*(u*u)), 0.25);
   }
   theta *= 0.5;

   f = sin(theta)/(u * rho);

   return f;
}/* end Imhof() */

/******************************************************************************
RunsTest()

Perform a runs test for autocorrelation of residuals. The residuals should be 
placed in an appropriate order by the calling routine (e.g. ordered by
increasing value of observation).

Code is based on the formulation presented by Sweed and Eisenhart in "Tables for 
testing randomness of grouping in a sequence of alternatives", The Annals of 
Mathematical Statistics, vol. 14, no. 1, March 1943, pg. 66-87.

NOTE: a value of zero is counted as positive.

Input:
   residuals : Array of residuals, runs test examines the sign of the residuals.
   n : size of residuals array
Output:
   Returns true if successful (npos and nneg > 1) and sets the following output 
   values:
      nPos   : number of positive values
      nNeg   : number of negative values
      nRuns  : number of runs
      plwr   : p-value for lower tail Prob(u < nruns)
      pupr   : p-value for upper tail Prob(u > nruns)
   Return false if npos or nneg is less than or equal to 1 (too few samples)
******************************************************************************/
bool RunsTest(double * residuals, int num, int * nPos, int * nNeg, int * nRuns, 
              double * plwr, double * pupr)
{
   int u, i, s1, s2, k, m, n;
   double sum, fu;
   
   *nPos = *nNeg = 0; 
   *nRuns = 1; 

   for(i = 0; i < num; i++)
   {
      if(residuals[i] < 0.00) (*nNeg)++;
      else                    (*nPos)++;
   }
   if((*nPos < 2) || (*nNeg < 2)){*plwr = *pupr = 0.50; return false;}

   for(i = 1; i < num; i++)
   {
      if(residuals[i] < 0.00) s1 = -1;
      else s1 = +1;
      if(residuals[i-1] < 0.00) s2 = -1;
      else s2 = +1;
      if(s1 != s2) (*nRuns)++;
   }

   sum = 0.00;
   if(*nPos <= *nNeg) {m = *nPos; n = *nNeg;}
   else               {m = *nNeg; n = *nPos;}

   for(u = 2; u <= *nRuns; u++)
   {
      //u is even
      if((u % 2) == 0)
      {
         k = u / 2;
         fu = 2.00* nCr(m-1, k-1) * nCr(n-1, k-1);         
      }
      //u is odd
      else
      {
         k = (u + 1) / 2;
         fu = (nCr(m-1, k-1) * nCr(n-1, k-2) + nCr(m-1, k-2) * nCr(n-1, k-1));
      }
      sum += fu;
   }/* end for() */
   *plwr = sum / nCr(m+n, n);
   *pupr = 1.00 - ((sum - fu)/nCr(m+n, n));

   return true;
}/* end RunsTest() */

/******************************************************************************
nCr()

Computes the number of combinations of n things taken r at a time.
   nCr = n! / [ k! (n-k)! ]
******************************************************************************/
double nCr(int n, int r)
{
   int i;
   double numer, denom;

   numer = denom = 1.00e00;

   if((n - r) < r) r = (n - r);

   for(i = 0; i < r; i++) numer *= ((double)(n - i));
   for(i = 0; i < r; i++) denom *= ((double)(r - i));

   return (numer / denom);
}/* end nCr() */

/******************************************************************************
STATS_TestAutocorrelation()

Test out the Durban-Watson and Runs autocorrelation tests.
******************************************************************************/
void STATS_TestAutocorrelation(void)
{
   int i, j;

   //double r[16] = {1,-1,1,1,-1,1,1,1,-1,1,1,1,-1,1,-1,1};
   //double r[20] = {1,-1,1,1,-1,1,-1,-1,-1,1,-1,1,-1,1,1,-1,-1,1,1,1};
    double r[30] = /*{2.86548,1.59472,-0.18592,3.27336,3.98284,1.78204,1.66188,
                  -3.35868,-3.799,-3.52944,-6.62036,-7.18124,-5.352,-2.98244,
                  1.08712,1.18672,2.56656,3.0962,3.45588,6.35548};*/

                  /* Example 2, p-val = 0.025, D = 1.08 
                  {-32.325,-26.598,2.22,-16.962,-1.144,-2.508,-1.963,11.673,
                   -0.509,27.036,-4.419,40.035,23.58,33.943,-2.785,-8.604,
                   0.577,6.849,-18.97,-29.062}; */

                  /* Example 2, p-val = 0.01, D = 0.95 
                  {-27.85108549,-21.94514121,6.992154982,-12.07054883,
                   3.866747361,2.74133974,3.584580214,17.45917259,
                   5.396468783,33.23970926,2.082949731,46.8948383,
                   30.73807877,41.51861544,5.267800196,-0.014366951,
                   9.703465902,16.45265066,-8.829516486,-18.20573935}; */

                   /* example 4, p-val = 0.05, D = 1.08 
                   {1.2911,0.0408,-2.1098,-1.1591,-0.0484,
                    0.2477,0.4824,0.5296,0.9292,0.5075,
                    0.5628,0.4638,-0.2064,-1.5155,-0.0941,
                    0.0000,0.0000,0.0000,0.0000,0.0000}; */

                   /* example 4, p-val = 0.025, D = 0.95 
                  {1.5852,0.3367,-1.7999,-0.8425,0.2645,
                   0.5602,0.8019,0.8535,1.2550,0.8425,
                   0.8908,0.7759,0.1282,-1.1640,0.2173,
                   0.0000,0.0000,0.0000,0.0000,0.0000}; */

                  /* example 4, p-val = 0.01, D = 0.81 
                  {1.7522,0.5047,-1.6238,-0.6626,0.4422,
                  0.7377,0.9834,1.0376,1.4401,1.0329,
                  1.0772,0.9532,0.3184,-0.9641,0.3942,
                  0.0000,0.0000,0.0000,0.0000,0.0000}; */

                  /* SHAZAM example 2, D = 0.57, lwr = 0.00, upr = 1.00 */
                  {-0.166860,-0.141225,-0.056331,-0.097397,-0.218790,
                   -0.180189,0.013461,-0.140583,-0.052185,0.223997,
                   0.246285,0.266461,0.242079,0.073455,0.022708,
                   0.059364,0.097253,0.155865,0.065346,0.073757,
                   0.094268,-0.045548,-0.081874,-0.155034,0.008652,
                  -0.093809,-0.129256,-0.139787,-0.057169,0.113163}; 

                  /* SHAZAM example 1, D = 2.01, lwr = 0.301270, upr = 0.698730 
                  {-5.507930127,-2.576793164,-1.421168928,5.181466085,0.251704188,
                    5.310204304, 1.945793402,-0.574298765,-4.395641673,-1.542614199,
                   -4.59467751,  4.956713133, 8.897180132,-6.415885462, 2.56132361,
                    7.288520277,-9.365001685, 0.00,       0.00,         0.00,
                    0.00,       0.00,         0.00,       0.00,         0.00,
                     0.00,       0.00,         0.00,       0.00,         0.00}; */
                   
   //jacobian/sensitivity matrix
   double J[30][6] = /* Example 1
                      {{1.00,159.3},{1.00,161.2},{1.00,162.8},{1.00,164.6},
                      {1.00,165.9},{1.00,167.9},{1.00,168.3},{1.00,169.7},
                      {1.00,170.5},{1.00,171.6},{1.00,173.9},{1.00,176.1},
                      {1.00,178.0},{1.00,179.1},{1.00,180.2},{1.00,181.2},
                      {1.00,181.6},{1.00,182.5},{1.00,183.3},{1.00,184.3}}; */
                      /* Example 2 
                      {{1.00,75},{1.00,78},{1.00,80},{1.00,82},{1.00,84},
                      {1.00,88},{1.00,93},{1.00,97},{1.00,99},{1.00,104},
                      {1.00,109},{1.00,115},{1.00,120},{1.00,127},{1.00,135},
                     {1.00,144},{1.00,153},{1.00,161},{1.00,170},{1.00,182}}; */
                      /* example 4 
                     {{1.00,794},{1.00,799},{1.00,837},{1.00,855},{1.00,845},
                      {1.00,844},{1.00,863},{1.00,875},{1.00,880},{1.00,905},
                      {1.00,886},{1.00,843},{1.00,904},{1.00,950},{1.00,841},
                      {1.00,0.00},{1.00,0.00},{1.00,0.00},{1.00,0.00},{1.00,0.00}}; */
                     /* SHAZAM example 2 */
                     {{1.00, 25.4, 9.90, 17, 1, 0},
                     {1.00, 26.70, 4.70, 18, 1, 0},
                     {1.00, 29.10, 1.90, 23, 1, 0},
                     {1.00, 29.20, 3.20, 28, 1, 0},
                     {1.00, 29.20, 1.90, 30, 1, 0},
                     {1.00, 27.80, 3.9, 27, 0, 0},
                     {1.00, 27.40, 3.9, 24, 0, 1},
                     {1.00, 28.00, 3.80, 23, 0, 1},
                     {1.00, 28.30, 5.9, 25, 0, 1},
                     {1.00, 28.80, 5.3, 23, 1, 0},
                     {1.00, 29.30, 3.3, 24, 1, 0},
                     {1.00, 29.40, 3.00, 25, 1, 0},
                     {1.00, 29.20, 2.90, 25, 1, 0},
                     {1.00, 29.40, 5.50, 25, 0, 0},
                     {1.00, 30.20, 4.40, 25, 0, 1},
                     {1.00, 31.00, 4.10, 24, 0, 1},
                     {1.00, 31.20, 4.30, 25, 0, 1},
                     {1.00, 31.5, 6.80, 25, 0, 0},
                     {1.00, 31.70, 5.50, 26, 0, 0},
                     {1.00, 32.30, 5.50, 27, 0, 0},
                     {1.00, 32.60, 6.70, 26, 0, 0},
                     {1.00, 32.70, 5.5, 26, 0, 0},
                     {1.00, 33.20, 5.70, 26, 0, 0},
                     {1.00, 33.60, 5.20, 26, 0, 0},
                     {1.00, 34.00, 4.5, 27, 1, 0},
                     {1.00, 34.60, 3.80, 27, 1, 0},
                     {1.00, 35.10, 3.80, 27, 1, 0},
                     {1.00, 35.50, 3.60, 28, 1, 0},
                     {1.00, 36.30, 3.50, 30, 1, 0},
                     {1.00, 36.70, 4.90, 33, 1, 0}}; 
                     /* SHAZAM example 1 
                     {{1,96.7,101},
                     {1,98.1,100.1},
                     {1,100,100},
                     {1,104.9,90.6},
                     {1,104.9,86.5},
                     {1,109.5,89.7},
                     {1,110.8,90.6},
                     {1,112.3,82.8},
                     {1,109.3,70.1},
                     {1,105.3,65.4},
                     {1,101.7,61.3},
                     {1,95.4,62.5},
                     {1,96.4,63.6},
                     {1,97.6,52.6},
                     {1,102.4,59.7},
                     {1,101.6,59.5},
                     {1,103.8,61.3}}; */
                   
   double ** JAC;
   JAC = new double *[30];
   for(i = 0; i < 30; i++) JAC[i] = new double[6];
   for(i = 0; i < 30; i++) for(j = 0; j < 6; j++) JAC[i][j] = J[i][j];

   double D, plwr, pupr;
   int pos, neg, runs;

   DurbinWatsonTest(r, JAC, 30, 6, &D, &plwr, &pupr);
   RunsTest(r, 30, &pos, &neg, &runs, &plwr, &pupr);

   for(i = 0; i < 6; i++) delete [] JAC[i];
   delete [] JAC;

}/* end STATS_TestAutocorrelation() */

/******************************************************************************
STATS_TestFdist()

Computes the F-distribution inverse CDF for numerous values of u,v and p.
The results can be compared to those generated in another program (e.g. Excel)
to verify proper operation of the F-distribution code.
******************************************************************************/
void STATS_TestFdist(void)
{
   int i, j, k;
   double val;
   int u[5] = {1,2,4,8,16};
   int v[7] = {1,2,4,8,16,32,64};
   double p[5] ={ 0.250, 0.500, 0.750, 0.975, 0.99};

   printf("***** F-distribution Inverse CDF *****\n");
   printf("u  v  0.25      0.50      0.75      0.975     0.99    \n");
   for(i = 0; i < 5; i++)
   {
      for(j = 0; j < 7; j++)
      {
         printf("%d  %d  ",u[i],v[j]);
         for(k = 0; k < 5; k++)
         {
            val = FdistInvCDF(u[i], v[j], p[k]);
            printf("%lf  ", val);
         }/* end for() */
         printf("\n");
      }/* end for() */
   }/* end for() */
   
}/* end STATS_TestFdist() */

/******************************************************************************
STATS_TestStudentDist()

Computes the t-distribution inverse CDF for numerous values of n,and p.
The results can be compared to those generated in another program (e.g. Excel)
to verify proper operation of the t-distribution code.
******************************************************************************/
void STATS_TestStudentDist(void)
{
   int i, j;
   double val;
   int u[11] = {1,2,4,8,16,32,64,128,256,512,1024};
   double p[5] ={ 0.550, 0.750, 0.800, 0.975, 0.995};

   printf("***** Student's t-distribution Inverse CDF *****\n");
   printf("u  0.550      0.750      0.800      0.975     0.995\n");
   for(i = 0; i < 11; i++)
   {
      printf("%d  ",u[i]);
      for(j = 0; j < 5; j++)
      {         
         val = StudentInvCDF(u[i], p[j]);
         printf("%lf  ", val);
      }/* end for() */
      printf("\n");
   }/* end for() */
}/* end STATS_TestStudentDist() */

/******************************************************************************
STATS_TestStdNormDist()

Computes the std. normal distribution inverse CDF for numerous values of p.
The results can be compared to those generated in another program (e.g. Excel)
to verify proper operation of the code.
******************************************************************************/
void STATS_TestStdNormDist(void)
{
   int i;
   double val;   

   printf("***** Standard Normal Inverse CDF *****\n");
   printf("Probability  X\n");
   for(i = 0; i < 19; i++)
   {
      printf("%lf  ",0.05*(i+1));
      val = StdNormInvCDF(0.05*(i+1));
      printf("%lf\n", val);
   }/* end for() */

   printf("%lf  ",0.99999);
   val = StdNormInvCDF(0.99999);
   printf("%lf\n", val);
}/* end STATS_TestStdNormDist() */
