#include <mex.h>
#include <limits>
#include <cmath>
#include "private/Gamma.h"

/*
Confidence level estimators for a Binomial distribution.

For a detailed discussion see:
    "Interval Estimation for a Binomial Proportion",
    L. D. Brown, T. T. Cai and A. DasGupta,
    2001, Statistical Sciences 16 (2) 101-133,
    http://www-stat.wharton.upenn.edu/~tcai/paper/Binomial-StatSci.pdf

The D0 implementation is detailed in the source code and documentation for the class TGraphAsymmErrors.

Implementation by Andrea Bocci <andrea.bocci@cern.ch>
Last modified on Wed Jun 25 19:23:52 CEST 2008

This is free software, licenced under the GNU LGPL version 2.1, or (at your option) any later version.
*/



class BinomialConfidence
{
protected:
  double GLOBAL_k;   // used to pass k[i] into equations
  double GLOBAL_N;   // used to pass N[i] into equations
  double CONFLEVEL;  // confidence level for the interval


  double Sign(double a, double b)
  {
    return (b >= 0) ? std::fabs(a) : -std::fabs(a);
  }

  //______________________________________________________________________________
  double BetaCf(double x, double a, double b)
  {
     // Continued fraction evaluation by modified Lentz's method
     // used in calculation of incomplete Beta function.

     int itmax = 500;
     double eps = 3.e-14;
     double fpmin = 1.e-30;

     int m, m2;
     double aa, c, d, del, qab, qam, qap;
     double h;
     qab = a+b;
     qap = a+1.0;
     qam = a-1.0;
     c = 1.0;
     d = 1.0 - qab*x/qap;
     if (std::fabs(d)<fpmin) d=fpmin;
     d=1.0/d;
     h=d;
     for (m=1; m<=itmax; m++) {
        m2=m*2;
        aa = m*(b-m)*x/((qam+ m2)*(a+m2));
        d = 1.0 +aa*d;
        if(std::fabs(d)<fpmin) d = fpmin;
        c = 1 +aa/c;
        if (std::fabs(c)<fpmin) c = fpmin;
        d=1.0/d;
        h*=d*c;
        aa = -(a+m)*(qab +m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if(std::fabs(d)<fpmin) d = fpmin;
        c = 1.0 +aa/c;
        if (std::fabs(c)<fpmin) c = fpmin;
        d=1.0/d;
        del = d*c;
        h*=del;
        if (std::fabs(del-1)<=eps) break;
     }
     if (m>itmax) {
       mexErrMsgIdAndTxt("BinomialConfidence:BetaCf", "a or b too big, or itmax too small, a=%g, b=%g, x=%g, h=%g, itmax=%d",
             a,b,x,h,itmax);
     }
     return h;
  }

  double Ibetai(double a, double b, double x)
  {
     // Calculates the incomplete beta function  I_x(a,b); this is
     // the incomplete beta function divided by the complete beta function

     double bt;
     if (x < 0.0 || x > 1.0) {
        mexErrMsgIdAndTxt("BinomialConfidence:Ibetai", "[Ibetai] Illegal x in routine Ibetai: x = %.5g", x);
        return 0;
     }
     if (x == 0.0 || x == 1.0)
        bt=0.0;
     else
        bt=std::exp(LogGamma(a+b)-LogGamma(a)-LogGamma(b)+a*log(x)+b*log(1.0-x));

     if (x < (a+1.0)/(a+b+2.0))
        return bt*BetaCf(x,a,b)/a;
     else
        return 1.0-bt*BetaCf(1-x,b,a)/b;
  }

  double Beta_ab(double a, double b, double k, double N)
  {
     // Calculates the fraction of the area under the
     // curve x^k*(1-x)^(N-k) between x=a and x=b

     if (a == b) return 0;    // don't bother integrating over zero range
     double c1 = k+1;
     double c2 = N-k+1;
     return Ibetai(c1,c2,b)-Ibetai(c1,c2,a);
  }




  double SearchLower(double high, double k, double N, double c)
  {
     // Integrates the binomial distribution with
     // parameters k,N, and determines what is the lower edge of the
     // integration region which ends at high, and which contains
     // probability content c. If a lower limit is found, the value is
     // returned. If no solution is found, the -1 is returned.
     // check to see if there is any solution by verifying that the integral down
     // to the minimum lower limit (0) is greater than c

     double integral = Beta_ab(0.0, high, k, N);
     if (integral == c) return 0.0;      // lucky -- this is the solution
     if (integral < c) return -1.0;      // no solution exists
     double too_low = 0.0;               // lower edge estimate
     double too_high = high;
     double test;

     // use a bracket-and-bisect search
     // LM: looping 20 times might be not enough to get an accurate precision.
     // see for example bug https://savannah.cern.ch/bugs/?30246
     // now break loop when difference is less than 1E-15
     // t.b.d: use directly the beta distribution quantile

     for (int loop=0; loop<50; loop++) {
        test = 0.5*(too_high + too_low);
        integral = Beta_ab(test, high, k, N);
        if (integral > c)  too_low = test;
        else too_high = test;
        if ( std::fabs(integral - c) <= 1.E-15) break;
     }
     return test;
  }

  double SearchUpper(double low, double k, double N, double c)
  {
     // Integrates the binomial distribution with
     // parameters k,N, and determines what is the upper edge of the
     // integration region which starts at low which contains probability
     // content c. If an upper limit is found, the value is returned. If no
     // solution is found, -1 is returned.
     // check to see if there is any solution by verifying that the integral up
     // to the maximum upper limit (1) is greater than c

     double integral = Beta_ab(low, 1.0, k, N);
     if (integral == c) return 1.0;    // lucky -- this is the solution
     if (integral < c) return -1.0;    // no solution exists
     double too_high = 1.0;            // upper edge estimate
     double too_low = low;
     double test;

     // use a bracket-and-bisect search
     // LM: looping 20 times might be not enough to get an accurate precision.
     // see for example bug https://savannah.cern.ch/bugs/?30246
     // now break loop when difference is less than 1E-15
     // t.b.d: use directly the beta distribution quantile

     for (int loop=0; loop<50; loop++) {
        test = 0.5*(too_low + too_high);
        integral = Beta_ab(low, test, k, N);
        if (integral > c)  too_high = test;
        else too_low = test;
        if ( std::fabs(integral - c) <= 1.E-15) break;
     }
     return test;
  }



  double Brent(double ax, double bx, double cx, double tol, double *xmin)
  {
     // Implementation file for the numerical equation solver library.
     // This includes root finding and minimum finding algorithms.
     // Adapted from Numerical Recipes in C, 2nd edition.
     // Translated to C++ by Marc Paterno

     const int    kITMAX = 100;
     const double kCGOLD = 0.3819660;
     const double kZEPS  = 1.0e-10;

     int iter;
     double a,b,d=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
     double e=0.0;

     a=(ax < cx ? ax : cx);
     b=(ax > cx ? ax : cx);
     x=w=v=bx;
     fw=fv=fx=Interval(x);
     for (iter=1;iter<=kITMAX;iter++) {
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*std::fabs(x)+kZEPS);
        if (std::fabs(x-xm) <= (tol2-0.5*(b-a))) {
           *xmin=x;
           return fx;
        }
        if (std::fabs(e) > tol1) {
           r=(x-w)*(fx-fv);
           q=(x-v)*(fx-fw);
           p=(x-v)*q-(x-w)*r;
           q=2.0*(q-r);
           if (q > 0.0) p = -p;
           q=std::fabs(q);
           etemp=e;
           e=d;
           if (std::fabs(p) >= std::fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) d=kCGOLD*(e=(x >= xm ? a-x : b-x));
           else {
              d=p/q;
              u=x+d;
              if (u-a < tol2 || b-u < tol2) d=Sign(tol1,xm-x);
           }
        } else {
           d=kCGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(std::fabs(d) >= tol1 ? x+d : x+Sign(tol1,d));
        fu=Interval(u);
        if (fu <= fx) {
           if (u >= x) a=x; else b=x;
           v  = w;
           w  = x;
           x  = u;
           fv = fw;
           fw = fx;
           fx = fu;
        } else {
           if (u < x) a=u; else b=u;
           if (fu <= fw || w == x) {
              v=w;
              w=u;
              fv=fw;
              fw=fu;
           } else if (fu <= fv || v == x || v == w) {
              v=u;
              fv=fu;
           }
        }
     }
     mexErrMsgIdAndTxt("BinomialConfidence:Brent", "[Brent] Too many interations");
     *xmin=x;
     return fx;
  }


  double Interval(double low)
  {
     // Return the length of the interval starting at low
     // that contains CONFLEVEL of the x^GLOBAL_k*(1-x)^(GLOBAL_N-GLOBAL_k)
     // distribution.
     // If there is no sufficient interval starting at low, we return 2.0

     double high = SearchUpper(low, GLOBAL_k, GLOBAL_N, CONFLEVEL);
     if (high == -1.0) return 2.0; //  so that this won't be the shortest interval
     return (high - low);
  }


public:
  void Efficiency(double k, double N, double conflevel,
        double& mode, double& low, double& high)
  {
     // Calculate the shortest central confidence interval containing the required
     // probability content.
     // Interval(low) returns the length of the interval starting at low
     // that contains CONFLEVEL probability. We use Brent's method,
     // except in two special cases: when k=0, or when k=N
     // Main driver routine
     // Author: Marc Paterno

     //If there are no entries, then we know nothing, thus return the prior...
     if (0==N) {
        mode = .5; low = 0.0; high = 1.0;
        return;
     }

     // Calculate the most probable value for the posterior cross section.
     // This is easy, 'cause it is just k/N
     double efficiency = (double)k/N;

     double low_edge;
     double high_edge;

     if (k == 0) {
        low_edge = 0.0;
        high_edge = SearchUpper(low_edge, k, N, conflevel);
     } else if (k == N) {
        high_edge = 1.0;
        low_edge = SearchLower(high_edge, k, N, conflevel);
     } else {
        GLOBAL_k = k;
        GLOBAL_N = N;
        CONFLEVEL = conflevel;
        Brent(0.0, 0.5, 1.0, 1.0e-9, &low_edge);
        high_edge = low_edge + Interval(low_edge);
     }

     // return output
     mode = efficiency;
     low = low_edge;
     high = high_edge;
  }
};



///////////////////////////////////////////////////////////////////////////
// Main entry point to a MEX function
///////////////////////////////////////////////////////////////////////////


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  
  // Check inputs to mex function
  if (nrhs < 2 || nrhs > 3 || nlhs < 1 || nlhs > 2) {
    mexErrMsgIdAndTxt ( "binointerval:usage", "Usage:  [phat, pci] = binointerval(x, n, [alpha = 0.05])" );
  }

  const int             nPoints     = mxGetNumberOfElements(prhs[0]);
  const double*         k           = mxGetPr(prhs[0]);
  const double*         n           = mxGetPr(prhs[1]);
  const double          level       = ( nrhs > 2 ? 1 - mxGetScalar(prhs[2]) : 0.954499736103642 );

  if (mxGetNumberOfElements(prhs[1]) != nPoints)
    mexErrMsgIdAndTxt("binointerval:arguments", "Inputs x and n must have the same number of elements, encountered %d != %d.", nPoints, mxGetNumberOfElements(prhs[1]));


  // Create output structures
  if (nlhs > 0)         plhs[0]     = mxCreateDoubleMatrix(nPoints, 1, mxREAL);
  if (nlhs > 1)         plhs[1]     = mxCreateDoubleMatrix(nPoints, 2, mxREAL);

  double*               phat        = ( nlhs > 0 ? mxGetPr(plhs[0])           : 0 );
  double*               pcilo       = ( nlhs > 1 ? mxGetPr(plhs[1])           : 0 );
  double*               pciup       = ( nlhs > 1 ? mxGetPr(plhs[1]) + nPoints : 0 );


  // Loop over data points
  BinomialConfidence    computer;
  for (int iPoint = 0; iPoint < nPoints; ++iPoint)
  {
    if (k[iPoint] <= n[iPoint])
      computer.Efficiency(k[iPoint], n[iPoint], level, phat[iPoint], pcilo[iPoint], pciup[iPoint]);
    else
      phat[iPoint] = pcilo[iPoint] = pciup[iPoint] = mxGetNaN();
  }
}
