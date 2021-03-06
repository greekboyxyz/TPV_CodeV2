#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

int main(int argc, char** argv)
{
    int i;
    double lam, n, k;
    FILE *epsfp, *lamfp, *wfp;
    epsfp = fopen("Fe3O4.txt","r");
    lamfp = fopen("W_Palik.txt","r");
    //wfp   = fopen("W_ald_interpolated","w");
    int dim = 609;
    std::vector<double> X(dim), Y(dim), Z(dim);
    std::vector<double> LAM(4000);

    for (i=0; i<dim; i++) {

      fscanf(epsfp,"%lf",&lam);
      fscanf(epsfp,"%lf",&n);
      fscanf(epsfp,"%lf",&k);
      X[i] = lam;
      Y[i] = n;
      Z[i] = k;
      //if (i>0) printf(" %i %12.10f\n",i,X[i]-X[i-1]); 
    }    
    for (i=0; i<4000; i++) {
    

      fscanf(lamfp,"%lf",&lam);
      fscanf(lamfp,"%lf",&n);
      fscanf(lamfp,"%lf",&k);
      LAM[i] = lam;

    }
    tk::spline s;
    tk::spline t;
    s.set_points(X,Y);    // currently it is required that X is already sorted
    t.set_points(X,Z);

    double x;

    for (i=0; i<4000; i++) {

      printf("%12.10f %12.10f  %12.10f\n", LAM[i], s(LAM[i]), t(LAM[i]));
   
    }

    return EXIT_SUCCESS;
}
