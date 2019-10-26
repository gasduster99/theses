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
#include <math.h>
#include "NumericalRecipes.h"
#include "Exception.h"

/*--------------- Numerical Recipes definitions---------------*/
#define FPNT(x) ((*func)(x))
#define FINF(x) ((*funk)(1.0/(x))/((x)*(x)))

//a Numerical Recipes routine, returns square of 'a'
double SQR(double a)
{
   if(a == 0.00){ return 0.00;}
   return (a*a); 
}/* end SQR() */

//return max of two doubles
double FMAX(double a, double b)
{
   if(a > b) return a;
   return b;
}/* end FMAX() */

//return max of two integers
int IMAX(int a, int b)
{
   if(a > b) return a;
   return b;
}/* end IMAX() */

//return min of two doubles
double FMIN(double a, double b)
{
   if(a < b) return a;
   return b;
}/* end FMAX() */

//match the sign of a to the sign of b
double SIGN(double a, double b) 
{
  if(b >= 0.0) return (fabs(a));
  return (-fabs(a));
}/* end SIGN() */

/*--------------- End Numerical Recipes definitions------------*/

void balanc(double **a, int n)
{
   double RADIX = 2.0;
	int last,j,i;
	double s,r,g,f,c,sqrdx;

	sqrdx=RADIX*RADIX;
	last=0;
	while (last == 0) {
		last=1;
		for (i=1;i<=n;i++) {
			r=c=0.0;
			for (j=1;j<=n;j++)
				if (j != i) {
					c += fabs(a[j-1][i-1]);
					r += fabs(a[i-1][j-1]);
				}
			if (c && r) {
				g=r/RADIX;
				f=1.0;
				s=c+r;
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g=r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last=0;
					g=1.0/f;
					for (j=1;j<=n;j++) a[i-1][j-1] *= g;
					for (j=1;j<=n;j++) a[j-1][i-1] *= f;
				}
			}
		}
	}
}/* end balanc() */

void elmhes(double **a, int n)
{
	int m,j,i;
	double y,x;

	for (m=2;m<n;m++) {
		x=0.0;
		i=m;
		for (j=m;j<=n;j++) {
			if (fabs(a[j-1][m-1-1]) > fabs(x)) {
				x=a[j-1][m-1-1];
				i=j;
			}
		}
		if (i != m) {
			for (j=m-1;j<=n;j++) SWAP(a[i-1][j-1],a[m-1][j-1])
			for (j=1;j<=n;j++) SWAP(a[j-1][i-1],a[j-1][m-1])
		}
		if (x) {
			for (i=m+1;i<=n;i++) {
				if ((y=a[i-1][m-1-1]) != 0.0) {
					y /= x;
					a[i-1][m-1-1]=y;
					for (j=m;j<=n;j++)
						a[i-1][j-1] -= y*a[m-1][j-1];
					for (j=1;j<=n;j++)
						a[j-1][m-1] += y*a[j-1][i-1];
				}
			}
		}
	}
}/* end elmhes() */

void hqr(double **a, int n, double * wr, double * wi)
{
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm;
    double EPS = 1.0e-6;
    double diff;

	anorm=0.0;
	for (i=1;i<=n;i++)
		for (j=IMAX(i-1,1);j<=n;j++)
			anorm += fabs(a[i-1][j-1]);
	nn=n;
	t=0.0;
	while (nn >= 1) {
		its=0;
		do {
			for (l=nn;l>=2;l--) {
				s=fabs(a[l-1-1][l-1-1])+fabs(a[l-1][l-1]);
				if (s == 0.0) s=anorm;
				diff = fabs((float)(fabs(a[l-1][l-1-1]) + s) - s);
				if (diff <= EPS){ break;}
			}
			x=a[nn-1][nn-1];
			if (l == nn) {
				wr[nn-1]=x+t;
				wi[nn-1]=0.0; nn--;
			} else {
				y=a[nn-1-1][nn-1-1];
				w=a[nn-1][nn-1-1]*a[nn-1-1][nn-1];
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1-1]=wr[nn-1]=x+z;
						if (z) wr[nn-1]=x-w/z;
						wi[nn-1-1]=wi[nn-1]=0.0;
					} else {
						wr[nn-1-1]=wr[nn-1]=x+p;
						wi[nn-1-1]= -(wi[nn-1]=z);
					}
					nn -= 2;
				} else {
					if (its == 30){ LogError(ERR_SING_MAT, "Too many iterations in hqr()"); return;}
					if (its == 10 || its == 20) {
						t += x;
						for (i=1;i<=nn;i++) a[i-1][i-1] -= x;
						s=fabs(a[nn-1][nn-1-1])+fabs(a[nn-1-1][nn-2-1]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2);m>=l;m--) {
						z=a[m-1][m-1];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1-1][m-1]+a[m-1][m+1-1];
						q=a[m+1-1][m+1-1]-z-r-s;
						r=a[m+2-1][m+1-1];
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a[m-1][m-1-1])*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a[m-1-1][m-1-1])+fabs(z)+fabs(a[m+1-1][m+1-1]));
						diff = fabs((float)(u+v) - v);
						if (diff <= EPS) break;
					}
					for (i=m+2;i<=nn;i++) {
						a[i-1][i-2-1]=0.0;
						if (i != (m+2)) a[i-1][i-3-1]=0.0;
					}
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a[k-1][k-1-1];
							q=a[k+1-1][k-1-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2-1][k-1-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a[k-1][k-1-1] = -a[k-1][k-1-1];
							} else
								a[k-1][k-1-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<=nn;j++) {
								p=a[k-1][j-1]+q*a[k+1-1][j-1];
								if (k != (nn-1)) {
									p += r*a[k+2-1][j-1];
									a[k+2-1][j-1] -= p*z;
								}
								a[k+1-1][j-1] -= p*y;
								a[k-1][j-1] -= p*x;
							}
							mmin = nn<k+3 ? nn : k+3;
							for (i=l;i<=mmin;i++) {
								p=x*a[i-1][k-1]+y*a[i-1][k+1-1];
								if (k != (nn-1)) {
									p += z*a[i-1][k+2-1];
									a[i-1][k+2-1] -= p*r;
								}
								a[i-1][k+1-1] -= p*q;
								a[i-1][k-1] -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}
}/* end hqr() */

#define JMAX   14
#define JMAXP (JMAX+1)
double qromo(double (*func)(double), double a, double b,
	double (*choose)(double(*)(double), double, double, int))
{
   double EPS = 1.0e-6;
   int K=5;
	int j;
	double ss,dss,h[JMAXP+1],s[JMAXP];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo()");
	return 0.0;
}/* end qromo() */

double midpnt(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FPNT(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FPNT(x);
			x += ddel;
			sum += FPNT(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}/* end midpnt() */

double midinf(double (*funk)(double), double aa, double bb, int n)
{
	double x,tnm,sum,del,ddel,b,a;
	static double s;
	int it,j;

	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		return (s=(b-a)*FINF(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FINF(x);
			x += ddel;
			sum += FINF(x);
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}/* end midinf() */

void polint(double * xa, double * ya, int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c = new double[1+n];
	d = new double[1+n];
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint()");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	delete [] d;
   delete [] c;
}/* end polint() */

void nrerror(char * error_text)
{
   LogError(ERR_NRINC, error_text);
}/* end() */

bool ludcmp(double **a, int n, int *indx, double *d)
{
	bool isOK = true;
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=new double[n];
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) 
		{ 
			nrerror("ludcmp() : singular matrix"); 
			isOK = false;
			goto exit;
		}
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
      imax = j;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0)
		{ 
			nrerror("ludcmp() : singular matrix"); 
			isOK = false;
			goto exit;
		}			
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
exit:
	delete [] vv;
	return isOK;
}/* end ludcmp() */
