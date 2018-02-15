#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include</usr/include/complex.h>

// Function Prototypes
int *VEC_INT(int dim);
double *VEC_DOUBLE(int dim);
char *VEC_CHAR(int dim);
double complex *VEC_CDOUBLE(int dim);
void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C);
double PrintSpectrum(char *filename, int NBR, int NML, double *d_array, double complex *nlow, double complex *nhi, double complex *abs_n, double complex *sub_eps, int NumLam, double *LamList, double Temp);
void TransferMatrix(int Nlayer,double thetaI, double k0, double complex *rind, double *d,
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21);
double SpectralEfficiency(double *emissivity, int N, double *lambda, double lambdabg, double Temperature, double *P);
void Bruggenman(double f, double complex epsD, double complex epsM, double *eta, double *kappa);
void MaxwellGarnett(double f, double epsD, double complex epsM, double *eta, double *kappa);
void Lorentz(double we, double de, double w, double *epsr, double *epsi);
int ReadDielectric(char *file, double *lambda, double complex *epsM);
int IsDominated(int idx, int LENGTH, double *O1,double *O2);
void ReadBRRind(int numBR, double lambda, double *BRlambda, double complex *BRind, double *n, double *k); 
void BuildML(int N_BR, int N_ML, double *d_layer, double complex *n, double complex *rind, double *d);
// Global variables
int polflag;
int alloyflag;
double c = 299792458;
double pi=3.141592653589793;

int main(int argc, char* argv[]) {
  // complex double precision variables
  double complex m11, m21, r, t, st, cosL;
  double complex *rind, nlow, nhi, *n_array;
  double complex *absorber_n, *BR_HighN, *BR_LowN, *substrate_eps;
  double SE, PU;

  double h=6.626e-34;
  double kb = 1.38064852e-23;
  double rho;

  // real double precision variables
  double R, T, A, *d, *d_array, thetaI, lambda, k0, alpha, beta;
  double sti, n1, n2, thetaT, rp, Rp, Tangle;
  double eta, kappa;
  double we, de, w;
  int Nlayer, N_BR, N_ML;
  // Lists for spectral efficiency
  double *LamList, *Emiss, *clam;
  // Variables for Spectral Efficiency
  // This is intentionally larger than number of wavelength values in data files
  int NumLam=10000;


  //  Allocate arrays for spectral efficiency
  LamList = VEC_DOUBLE(NumLam);
  Emiss   = VEC_DOUBLE(NumLam);
  clam    = VEC_DOUBLE(NumLam);

  //  Arrays for material constants
  absorber_n      = VEC_CDOUBLE(NumLam);
  BR_HighN        = VEC_CDOUBLE(NumLam);
  BR_LowN         = VEC_CDOUBLE(NumLam);

  substrate_eps   = VEC_CDOUBLE(NumLam);

  FILE *fp;

  // Character string(s)
  char *write, *line, *subfile, *absorberfile, *BR_HighFile, *BR_LowFile, *prefix, *pfile;

  write   = VEC_CHAR(1000);
  line    = VEC_CHAR(1000);
  pfile   = VEC_CHAR(1000);
  prefix = VEC_CHAR(1000);

  subfile        = VEC_CHAR(1000);
  absorberfile   = VEC_CHAR(1000);
  BR_HighFile    = VEC_CHAR(1000);
  BR_LowFile     = VEC_CHAR(1000);

  //  Did we pass a filename to the program?
  if (argc==1) {
    exit(0);
  }

  strcpy(write,argv[1]);

  
  // initialize variables to be varied over
  int N_min=0;
  int N_max=0;
  int NumVars=0;
  double d1 = 0.;
  double d1min=0.0;
  double d1max=0.0;
  double d1_delta = 0.;
  double d2 = 0.;
  double d2min=0.0;
  double d2max=0.0;
  double d2_delta = 0.;
  double d3 = 0.0;
  double d3min=0.0;
  double d3max=0.0;
  double d3_delta = 0.;
  double d4 = 0.0;
  double d4min = 0.0;
  double d4max = 0.0;
  double d4_delta = 0;
  double Temp = 0.;
  double Tempmin=0.;
  double Tempmax=0.;
  double T_delta = 0.;

  // Fixed PV lambda_bg
  double lbg=2254e-9;;
  // Open the file for writing!
  fp = fopen(write,"r");
  printf("  going to read from file %s\n",write);
  fflush(stdout);
 
  // read info about min and max number of layers
  fscanf(fp,"%s",line);
  printf("%s\n",line);  // Nlayer
  fscanf(fp,"%i",&N_min);
  fscanf(fp,"%i",&N_max);
  // read info about min and max thickness of low-RI materials
  fscanf(fp,"%s",line);
  fscanf(fp,"%lf",&d1min);
  fscanf(fp,"%lf",&d1max);
  // read info about min and max thickness of high-RI materials
  fscanf(fp,"%s",line);
  fscanf(fp,"%lf",&d2min);
  fscanf(fp,"%lf",&d2max);
  // read info about min and max thickness of dielectric in MultiLayer
  fscanf(fp,"%s",line);
  fscanf(fp,"%lf",&d3min);
  fscanf(fp,"%lf",&d3max);
  // read info about min and max thickness of metal in MultiLayer
  fscanf(fp,"%s",line);
  fscanf(fp,"%lf",&d4min);
  fscanf(fp,"%lf",&d4max);
  // read info about min and max temperature
  fscanf(fp,"%s",line);
  fscanf(fp,"%lf",&Tempmin);
  fscanf(fp,"%lf",&Tempmax);
  // File prefix for output files (Pareto front and spectra)
  fscanf(fp,"%s",line);
  fscanf(fp,"%s",prefix);
  strcpy(pfile,prefix);
  strcat(pfile,"_Pareto.txt");
  // File name to read absorber data from 
  fscanf(fp,"%s",line);
  fscanf(fp,"%s",absorberfile);
  // Indicate if we will use MaxwellGarnett or Bruggenman
  fclose(fp);

  if (alloyflag!=1 && alloyflag!=2) {
	  printf("  INVALID CHOICE FOR ALLOY!\n");
	  printf("  DEFAULTING TO MAXWELL-GARNETT\n");
          alloyflag=1;
  }
  
  int matflag=0;
  if (strcmp(absorberfile,"DIEL/W_Palik.txt")==0) {

          printf("  READING FROM W\n");
 	  matflag = 1;
  }
  else {
	  matflag = 0;
	  printf("  NOT W\n");

  }

  // Now define the specific file names
  // Substrate is Tungsten
  strcpy(subfile,"DIEL/W_Palik.txt");
  // Dielectric in Alloy
  strcpy(BR_LowFile,"DIEL/Al2O3_Spline.txt");
  strcpy(BR_HighFile,"DIEL/HfO2_Spline.txt");
  int CheckNum;
  // How many data points are in the file W_Palik.txt?  This function  will tell us

  // Substrate W data - dielectric function
  NumLam =   ReadDielectric(subfile, LamList, substrate_eps);
  // Alloy Materials
  CheckNum = ReadDielectric(absorberfile, clam, absorber_n);
  CheckNum = ReadDielectric(BR_LowFile, clam, BR_LowN);
  CheckNum = ReadDielectric(BR_HighFile, clam, BR_HighN);

  // Refractive index of alumina
  nlow = 1.76+0.*I;
  // Refractive index of ZrO2
  nhi  = 2.15+0.*I;

  printf("#  Read from files!\n");
  printf("#  Read %i entries from file %s\n",NumLam,subfile);
  printf("#  Read %i entries from file %s\n",CheckNum,absorberfile);
  printf("#  Read %i entries from file %s\n",CheckNum,BR_LowFile);
  // How many variations will we try?
  NumVars = N_max-N_min;
  printf("  %i\n",N_min);
  printf("  %i\n",N_max);
  printf("  %i\n",NumVars);
  // What will be the delta on d1?
  d1_delta = (d1max-d1min)/(NumVars-1);
  printf("  d1 between %f and %f in increments of %f\n",d1min,d1max,d1_delta);
  // What will be the delta on d2?
  d2_delta = (d2max-d2min)/(NumVars-1);
  printf("  d2 between %f and %f in increments of %f\n",d2min,d2max,d2_delta);
  // What will be the delta on vf?
  d3_delta = (d3max-d3min)/(NumVars-1);
  printf("  vf between %f and %f in increments of %f\n",d3min,d3max,d3_delta);
  // What will be the delta on T?
  T_delta = (Tempmax-Tempmin)/(NumVars-1);
  printf("  T between %f and %f in increments of %f\n",Tempmin,Tempmax,T_delta);

  Temp = Tempmin; 
  polflag=1;

  int *NLA, *PF;
  double *D1A, *D2A, *D3A, *D4A, *TA, *SEA, *SDA;

  int totalVars = pow(NumVars,5.0);
  // Arrays that are only NumVars long
  NLA = (int*)malloc(NumVars*sizeof(int));
  D1A = (double*)malloc(NumVars*sizeof(double));
  D2A = (double*)malloc(NumVars*sizeof(double));
  D3A = (double*)malloc(NumVars*sizeof(double));
  D4A = (double*)malloc(NumVars*sizeof(double));
  TA  = (double*)malloc(NumVars*sizeof(double));

  // Arrays that are totalVars long
  SEA = (double*)malloc(totalVars*sizeof(double));
  SDA = (double*)malloc(totalVars*sizeof(double));
  PF  = (int*)malloc(totalVars*sizeof(int));

  // Arrays for TMM - 1000 is excessively long, but safe
  d = VEC_DOUBLE(1000);
  rind = VEC_CDOUBLE(1000);
  d_array = VEC_DOUBLE(4);
  n_array = VEC_CDOUBLE(5);

  int varcount=-1;
  printf("  NL d1        d2       d3        d4         Temp         SE                SD\n");
  // Loop over Nlayer - keeping fixed for now!
  for (int i=0; i<NumVars; i++) {

    Nlayer = N_min+i;
    N_ML = 3;
    N_BR = Nlayer - N_ML - 3;
    NLA[i] = Nlayer;

    // Loop over d1
    for (int j=0; j<NumVars; j++) {

      d1 = d1min + j*d1_delta;
      D1A[j] = d1; 
      d_array[0] = d1;
      // Loop over d2
      for (int k=0; k<NumVars; k++) {
      
        d2 = d2min + k*d2_delta;
	D2A[k] = d2;
        d_array[1] = d2;
        // Loop over d3
	for (int l=0; l<NumVars; l++) {
       
          d3 = d3min + l*d3_delta;
	  D3A[l] = d3;
          d_array[2] = d3;
	  // Loop over d4
	  for (int m=0; m<NumVars; m++) {

            d4 = d4min + m*d4_delta;
            D4A[m] = d4;
            d_array[3] = d4;
	    varcount++;
   
            // Normal incidence
            thetaI = 0;
            // Structure is established, now analayze it for its spectrum 
            for (int ii=0; ii<NumLam; ii++) {

	     	     
             lambda = LamList[ii];    // Lambda in meters
             k0 = 2*pi*1e-6/lambda;  // k0 in inverse microns - verified
             w=2*pi*c/lambda;        // angular frequency 

             // Fill n_array
             n_array[0] = BR_LowN[ii];
	     n_array[1] = BR_HighN[ii];
	     n_array[2] = BR_LowN[ii];
	     n_array[3] = absorber_n[ii]; 
	     n_array[4] = csqrt(substrate_eps[ii]);

	     // Now build ML rind and d arrays
	     BuildML(N_BR, N_ML, d_array, n_array, rind, d);

	     // Solve the Transfer Matrix Equations
	     TransferMatrix(Nlayer, thetaI, k0, rind, d, &cosL, &beta, &alpha, &m11, &m21);

	     rho = (2*h*c*c/pow(lambda,5))*(1/(exp(h*c/(lambda*kb*Temp))-1));
 
	     // Fresnel reflection coefficient (which is complex if there are absorbing layers)
	     r = m21/m11; 
 
             // Fresnel transmission coefficient (also complex if there are absorbing layers)
             t = 1./m11;
 
             // Reflectance, which is a real quantity between 0 and 1
             R = creal(r*conj(r));
             Tangle =  1.0*creal(cosL)/(1.0*cos(thetaI));
             T = creal(t*conj(t))*Tangle;
             A = 1 - R - T;
 
             // Store absorbance/emissivity in array Emiss
             Emiss[ii] = A;
           }
           SE = SpectralEfficiency(Emiss, NumLam, LamList, lbg, Temp, &PU);
           printf(" %i %8.6f %8.6f %8.6f %8.6f %8.6f %12.10e %12.10e\n",Nlayer, d1, d2, d3, d4, Temp, SE, PU);
           fflush(stdout);
	   SEA[varcount] = SE;
           SDA[varcount] = PU;

         }
       }
     }
   }
 }

  
 FILE *pf;
 pf = fopen(pfile,"w");
 int id;
 varcount=-1;
 int po=0;
 char *specfile;
 specfile = VEC_CHAR(1000);
 for (int i=0; i<NumVars; i++) {

   for (int j=0; j<NumVars; j++) {

     for (int k=0; k<NumVars; k++) {

       for (int l=0; l<NumVars; l++) {
  
         for (int m=0; m<NumVars; m++) {

           varcount++;
           //printf("  varcount is %i\n",varcount);
	   //fflush(stdout);
	   // id is 1 if member i is dominated by at least one other member j!=i
	   // a member is pareto optimal only if it is NOT dominated 

	   id = IsDominated(varcount, totalVars, SEA, SDA);

	   if (id) PF[varcount] = 0;

	   else {

	      	   PF[varcount] = 1;
                   po++;
		   char lab[10];
		   sprintf(lab, "%d", po);
		   strcpy(specfile,prefix);
		   strcat(specfile,lab);
		   strcat(specfile,"_spectra.txt");
                   // Note VFA currently holds d3... might want to change the name!
		   fprintf(pf,"  %i  %f  %f     %f       %f   %12.10f  %12.10e\n",
                                   
		      		   NLA[i],D1A[j],D2A[k],D3A[l],D4A[m],SEA[varcount],SDA[varcount]);
                   d_array[0] = D1A[j];
		   d_array[1] = D2A[k];
		   d_array[2] = D3A[l];
		   d_array[3] = D4A[m];
		   //PrintSpectrum(char *filename, int NBR, int NML, double d_array, double complex *nlow, double complex *nhi, double complex *abs_n, double complex *sub_eps, int NumLam, double *LamList, double Temp)
		   double RT = PrintSpectrum(specfile,NLA[i]-6, 3, d_array, BR_LowN, BR_HighN, absorber_n, substrate_eps, NumLam, LamList, Temp);
	   }
       	 }
       }
     }
   }
 }
 fclose(fp);
 fclose(pf);

 return 0;

}

// Functions
//
double PrintSpectrum(char *filename, int NBR, int NML, double *d_array, double complex *nlow, double complex *nhi, double complex *abs_n, double complex *sub_eps, int NumLam, double *LamList, double Temp) {
           int Nlayer = NBR + NML + 3;
           FILE *fp;
           double Tangle;
      	   double eta, kappa;
           double h=6.626e-34;
           double kb = 1.38064852e-23;
	   double *Emiss;
	   double lbg = 2254e-9;
	   double PU, SE;
	   Emiss = (double *)malloc(NumLam*sizeof(double));
           fp = fopen(filename,"w");

           //  Top/Bottom layer RI for Transmission calculation
           double n1 = 1.;
           double n2 = 1.;
           double complex *rind, *n_array;
	   double *d;
	   rind    = VEC_CDOUBLE(1000);
	   n_array = VEC_CDOUBLE(5); 
	   d       = VEC_DOUBLE(1000);
           // Normal incidence
           double thetaI = 0;
           // Structure is established, now analayze it for its spectrum
           double complex m21, m11, cosL, r, t;
	   double beta, alpha, R, T, A;
	   for (int ii=0; ii<NumLam; ii++) {

             double lambda = LamList[ii];    // Lambda in meters
             double k0 = 2*pi*1e-6/lambda;  // k0 in inverse microns - verified
             double w=2*pi*c/lambda;        // angular frequency

	                  // Weak absorber and dielectric data stored as RI, want them in eps
             // for effective medium theory
             n_array[0] = nlow[ii];
	     n_array[1] = nhi[ii];
	     n_array[2] = nlow[ii];
	     n_array[3] = abs_n[ii];
	     n_array[4] = csqrt(sub_eps[ii]);


	     BuildML(NBR, NML, d_array, n_array, rind, d);

             // Solve the Transfer Matrix Equations
             TransferMatrix(Nlayer, thetaI, k0, rind, d, &cosL, &beta, &alpha, &m11, &m21);

	     double rho = (2*h*c*c/pow(lambda,5))*(1/(exp(h*c/(lambda*kb*Temp))-1));


             // Fresnel reflection coefficient (which is complex if there are absorbing layers)
             r = m21/m11;

             // Fresnel transmission coefficient (also complex if there are absorbing layers)
             t = 1./m11;

             // Reflectance, which is a real quantity between 0 and 1
             R = creal(r*conj(r));
             Tangle =  n2*creal(cosL)/(n1*cos(thetaI));
             T = creal(t*conj(t))*Tangle;
             A = 1 - R - T;

             // Store absorbance/emissivity in array Emiss
             Emiss[ii] = A;
             fprintf(fp,"%8.6e  %8.6f  %8.6f  %8.6f  %8.6f\n",LamList[ii],R,A,rho*A,rho);
           }
           SE = SpectralEfficiency(Emiss, NumLam, LamList, lbg, Temp, &PU);
           fprintf(fp,"# %i  %8.6f %8.6f  %8.6f  %8.6f  %8.6f  %12.10e  %12.10e\n",Nlayer, d_array[0], d_array[1], d_array[2], d_array[3],Temp, SE, PU);
           free(Emiss);
	   fclose(fp);
	   free(rind);
	   free(d);
           return R;
}


int *VEC_INT(int dim){
  int *v,i;
  v = (int *)malloc(dim*sizeof(int));
  if (v==NULL) {
     printf("\n\nVEC_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0;
  return v;
}
double *VEC_DOUBLE(int dim){
  int i;
  double *v;
  v = (double *)malloc(dim*sizeof(double));
  if (v==NULL) {
     printf("\n\nVEC_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0;
  return v;
}

double complex *VEC_CDOUBLE(int dim) {
  int i;
  double complex *v;
  v = (double complex *)malloc(dim*sizeof(double complex));
  if (v==NULL) {
     printf("\n\nVEC_CDOUBLE:  Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0 + I*0.;
  return v;
}

char *VEC_CHAR(int dim){
  char *v;
  v = (char *)malloc(dim*sizeof(char));
  if (v==NULL) {
     printf("\n\nVEC_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  return v;
}

void TransferMatrix(int Nlayer, double thetaI, double k0, double complex *rind, double *d, 
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21) { 
  
  int i, j, k, indx;
  double complex *kz, *phiL, *D, *Dinv, *Pl, *phil;
  double complex *EM, ctheta, tmp, *tmp2, *tmp3, c0, c1, ci, kx;

  kz   = VEC_CDOUBLE(Nlayer);
  phil = VEC_CDOUBLE(Nlayer);
  D    = VEC_CDOUBLE(4*Nlayer);
  Dinv = VEC_CDOUBLE(4*Nlayer);
  Pl   = VEC_CDOUBLE(4*Nlayer);
  EM   = VEC_CDOUBLE(4);
  tmp2 = VEC_CDOUBLE(4); 
  tmp3 = VEC_CDOUBLE(4);

  c0 = 0. + I*0.;
  c1 = 1. + I*0.;
  ci = 0. + I*1.;

  //  x-component of incident wavevector...
  //  should be in dielectric material, so the imaginary 
  //  component should be 0.
  kx = k0*rind[0]*sin(thetaI);

  //  Now get the z-components of the wavevector in each layer
  for (i=0; i<Nlayer; i++) {
     kz[i] = (rind[i]*k0)*(rind[i]*k0) - kx*kx;
     kz[i] = csqrt(kz[i]);
     // Want to make sure the square root returns the positive imaginary branch
     if (cimag(kz[i])<0.)  {
        kz[i] = creal(kz[i]) - cimag(kz[i]);
     }
   }



   //  Calculate the P matrix
   for (i=1; i<Nlayer-1; i++) {
     phil[i]=kz[i]*d[i];

     //  Upper left (diagonal 1)
     Pl[i*4] = cexp(-ci*phil[i]);  
     //  upper right (off diagonal 1)
     Pl[i*4+1] = c0;
     //  lower left (off diagonal 2)
     Pl[i*4+2] = c0;
     //  lower right (diagonal 2)
     Pl[i*4+3] = cexp(ci*phil[i]);

   }

 
   //  Calculate the D and Dinv matrices
   for (i=0; i<Nlayer; i++) {
     ctheta = kz[i]/(rind[i]*k0);
     //  p-polarized incident waves
     if (polflag==1) {  

       //  Upper left (diagonal 1)
       D[i*4] = ctheta;
       // upper right
       D[i*4+1] = ctheta;
       // lower left
       D[i*4+2] = rind[i];
       // lower right
       D[i*4+3] = -rind[i];

     } 
     //  s-polarized incident waves
     if (polflag==2) {

       // upper left
       D[i*4] = 1;
       // upper right
       D[i*4+1] = 1;
       // lower left
       D[i*4+2] = rind[i]*ctheta;
       // lower right
       D[i*4+3] = -1*rind[i]*ctheta;

     }
     //  Now compute inverse
     //  Compute determinant of each D matrix
     tmp = D[i*4]*D[i*4+3]-D[i*4+1]*D[i*4+2];
     tmp = 1./tmp;

     //printf("  tmp is %12.10f  %12.10f\n",creal(tmp),cimag(tmp));    
     Dinv[i*4]=tmp*D[i*4+3];
     Dinv[i*4+1]=-1*tmp*D[i*4+1];
     Dinv[i*4+2]=-1*tmp*D[i*4+2];
     Dinv[i*4+3]=tmp*D[i*4];
 
   }


   // Initial EM matrix
   EM[0] = c1;
   EM[1] = c0;
   EM[2] = c0;
   EM[3] = c1;
   for (i=Nlayer-2; i>0; i--) {
     CMatMult2x2(i, Pl  , i,  Dinv, 0, tmp2);
     CMatMult2x2(i, D   , 0, tmp2,  0, tmp3);
     CMatMult2x2(0, tmp3, 0, EM  ,  0, tmp2); 

     for (j=0; j<2; j++) {
       for (k=0; k<2; k++) {
          EM[2*j+k] = tmp2[2*j+k];
       }
     }
   }
   CMatMult2x2(0, EM  , Nlayer-1, D   , 0, tmp2);
   CMatMult2x2(0, Dinv, 0, tmp2, 0, EM); 


   //  Finally, collect all the quantities we wish 
   //  to have available after this function is called
   *m11 = EM[0*2+0];  //  
   *m21 = EM[1*2+0];
   *beta = creal(kx);
   *alpha = cimag(kx);
   *cosL = ctheta;

   free(kz);  
   free(phil);
   free(D);
   free(Dinv);
   free(Pl);
   free(EM);
   free(tmp2);
   free(tmp3);

}

 

void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C) {
     int i, j, k, m, n;

     double complex sum;

     for (k=0; k<2; k++) {
       for (i=0; i<2; i++) {

          sum = 0. + 0*I;
          for (j=0; j<2; j++) {

            m = 2*i + j;
            n = 2*j + k;           
            sum += A[Aidx*4+m]*B[Bidx*4+n];

          }
          
          C[Cidx*4 + (2*i+k)] = sum;
        }
      }

}

double SpectralEfficiency(double *emissivity, int N, double *lambda, double lbg, double T, double *P){
    int i;
    double dlambda, sumD, sumN;
    double l, em;
    double h = 6.626e-34;
    double kb = 1.38064852e-23;
    double rho;

    sumD = 0;
    sumN = 0;

    for (i=1; i<N-1; i++) {

      em = emissivity[i];
      l = lambda[i];
      rho = (2.*pi*h*c*c/pow(l,5))*(1/(exp(h*c/(l*kb*T))-1));

      dlambda = fabs((lambda[i+1]- lambda[i-1])/(2));
      sumD += em*rho*dlambda;

      if (l<=lbg) {

        sumN += (l/lbg)*em*rho*dlambda;

      }
    }

    *P = sumN;
 
    return sumN/sumD;

}


void Bruggenman(double f, double complex epsD, double complex epsM, double *eta, double *kappa) {
  // medium 1 is surrounding medium (dielectric)
  // medium 2 is inclusion (W) - f passed to function is volume fraction of inclusion
  double f1, f2;
  double complex b, eps1, eps2, epsBG;
  eps1 = epsD;
  eps2 = epsM;


  f1 = (1 - f);
  f2 = f;
  b = (2*f1 - f2)*eps1 + (2*f2 - f1)*eps2;

  // The real part of epsBG can be positive or negative
  // but the imaginary part needs to be positive
  epsBG = (b + csqrt(8.*eps1*eps2 + b*b))/4.;

  // test to see that epsBG satisfy Bruggenman condition
  double complex test;
  // both real and imaginary part of n should be positive
   *eta   = fabs(creal(csqrt(epsBG)));
   *kappa = fabs(cimag(csqrt(epsBG)));

}


void MaxwellGarnett(double f, double epsD, double complex epsM, double *eta, double *kappa) {
   double complex num, denom;

   num   = epsD*(2*f*(epsM - epsD) + epsM + 2*epsD);
   denom = 2*epsD + epsM + f*(epsD-epsM); 

   *eta   = creal(csqrt(num/denom));
   *kappa = cimag(csqrt(num/denom));

}

//  Evaluates real and imaginary part of refractive index from 
//  the Lorent oscillator model given omega_0, gamma_0, and omega
void Lorentz(double we, double de, double w, double *nreal, double *nimag) {

  double complex epsilon;
  double complex n;

  
  epsilon = 1 + pow(we,2)/(pow(we,2) - 2*I*de*w - pow(w,2));

  //printf("  w:  %12.10e  we:  %12.10f  de:  %12.10f  epsr:  %12.10f  epsi:  %12.10f\n",w,we,de,creal(epsilon),cimag(epsilon));
  n = csqrt(epsilon);

  *nreal = creal(n);
  *nimag = cimag(n);

}

int ReadDielectric(char *file, double *lambda, double complex *epsM) {
   int i;
   FILE *fp;
   double lam, epsr, epsi;

   fp = fopen(file,"r");

   i=0;
   while(!feof(fp)) {

     fscanf(fp, "%lf",&lam);
     fscanf(fp, "%lf",&epsr);
     fscanf(fp, "%lf",&epsi);

     lambda[i] = lam;
     epsM[i]   = epsr + I*epsi;

     i++;
   }

   printf("#  There are %i elements in file %s\n",i,file);
   fflush(stdout);
   return i;
   fclose(fp);
}

void ReadBRRind(int numBR, double lambda, double *BRlambda, double complex *BRind, double *n, double *k) {
  int i, fdx, bdx, die;
  double temp, eta, kappa;

  // The wavelength we are interested in is smaller than any in the range of data
  if (lambda<BRlambda[0]) {

    *n = creal(BRind[0]) + (lambda - BRlambda[0])*((creal(BRind[1]) - creal(BRind[0]))/(BRlambda[1] - BRlambda[0]));
    *k = cimag(BRind[0]) + (lambda - BRlambda[0])*((cimag(BRind[1]) - cimag(BRind[0]))/(BRlambda[1] - BRlambda[0]));


  }
  // The wavelength we are interested in is larger than any in the range of data
  else if (lambda>BRlambda[numBR-2]) {

    *n = creal(BRind[numBR-2]) +(lambda - BRlambda[numBR-2])*((creal(BRind[numBR-2]) - creal(BRind[numBR-3]))/(BRlambda[numBR-2] - BRlambda[numBR-3]));
    *k = cimag(BRind[numBR-2]) +(lambda - BRlambda[numBR-2])*((cimag(BRind[numBR-2]) - cimag(BRind[numBR-3]))/(BRlambda[numBR-2] - BRlambda[numBR-3]));


  }
  // We need to scan the data to find the BRlambda for two lambdas that straddle the lambda of interest
  else {

    i=0; 
    die=1;
    do {

      temp = BRlambda[i];
      if (temp>lambda) {
      
        die=0;
        fdx = i;
        bdx = i-1; 

      }
      else i++; 

    }while(die);

    *n = creal(BRind[bdx]) + (lambda - BRlambda[fdx])*((creal(BRind[fdx]) - creal(BRind[bdx]))/(BRlambda[fdx] - BRlambda[bdx]));
    *k = cimag(BRind[bdx]) + (lambda - BRlambda[fdx])*((cimag(BRind[fdx]) - cimag(BRind[bdx]))/(BRlambda[fdx] - BRlambda[bdx]));
  
  }

}


int IsDominated(int idx, int LENGTH, double *O1,double *O2) {
  int i, is, rval;
  double Val1, Val2;

  Val1 = O1[idx];
  Val2 = O2[idx];

  // start by assuming solution is NOT dominated
  rval = 0;
  for (i=0; i<LENGTH; i++)

      if (i!=idx) {

        // Trying to maximize the function, xi dominates xidx if 
        // fj(xi) >= fj(xidx) for all j and fj(xi) < fj(xidx) for at least one j
        if ((O1[i]>=Val1 && O2[i]>=Val2) && (O1[i]>Val1 || O2[i]>Val2)) {

          //printf("  x%i is dominated by x%i\n",idx,i);
          //printf("  f1(%i):  %12.10f  f1(%i): %12.10f  f2(%i): %12.10f  f2(%i):  %12.10f\n",
          //idx,Val1,i,O1[i],idx,Val2,i,O2[i]);
          //printf("  terminating early!  i is %i out of %i\n",i,LENGTH);
          i=LENGTH;
          rval = 1;

        }
      }
  return rval;
}

// Will call this function for each wavelength in the spectrum.
// It will take as arguments: 
// The number of layers in the Bragg Reflector (N_BR)
// The number of layers in the MultiLayer absorber (N_ML) (total number of layer N_Layer = N_BR + N_ML + 3)
// The array of thicknesses of d1 (low-RI layer in BR), d2 (high-RI in BR), d3 (dielectric in ML), d4 (metal in MOL)
// The array of complex refractive indices of each layer (n1 -> low RI in BR, n2 -> high RI in BR,
// n3 -> dielectric in ML, n4 -> metal in ML, n5 -> metal in substrate)
// The function will return an array *rind filled with the complex RI of each layer in the target
// stack and  and array *d filled with the thickness of each layer in the target stack
void BuildML(int N_BR, int N_ML, double *d_layer, double complex *n, double complex *rind, double *d) {

	int Numlayer = N_BR + N_ML + 3;
        int ML_idx = N_BR + N_ML;
	double d1, d2, d3, d4;
	double complex n1, n2, n3, n4, n5;

	// thickness of low-RI layer in BR
	d1 = d_layer[0];
	// RI of low_RI layer in BR
	n1 = n[0];
	// thickness of high-RI layer in BR
	d2 = d_layer[1];
	// RI of high_RI layer in BR
	n2 = n[1];
	// thickness of dielectric in ML
	d3 = d_layer[2];
	// RI of dielectric in ML
	n3 = n[2];
	// thickness of metal in ML
	d4 = d_layer[3];
	// RI of metal in ML
	n4 = n[3];
	// RI of W substrate
	n5 = n[4];

	// First and last layers of structure are always air with thickness d=0
	rind[0] = 1.0+0*I;
	d[0] = 0.;

	rind[Numlayer-1] = 1.0+0.*I;
	d[Numlayer-1] = 0.;

	// Second to last layer of structure is always optically-thick Tungsten
	rind[Numlayer-2] = n5; 
	d[Numlayer-2] = 0.9;

	// Take care of Bragg Reflector
	for (int i=1; i<=N_BR; i++) {
                // Even index -> high RI layers
		if (i%2==0) {
			rind[i] = n2;
			d[i] = d2;
		}
		// Odd index -> low RI layers
		else {
			rind[i] = n1;
			d[i] = d2;
		}
	}

	// Take care of MultiLayer
	for (int i=N_BR+1; i<=ML_idx; i++) {

		// Even index -> dielectric
		if (i%2==0) {
			rind[i] = n3;
			d[i] = d3;
		}
		// Odd index -> metal
		else {
			rind[i] = n4;
			d[i] = d4;
		}
	}

/*	for (int i=0; i<Numlayer; i++) {
		printf("  %12.10e (%12.10e + i*%12.10e)\n",d[i],creal(rind[i]),cimag(rind[i]));
	}
*/

}

