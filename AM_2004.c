/* standard header files */
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>

/* calculate expected SAD using equations 6 and 8 in Alonso and McKane(2004) Sampling Hubbells neutral theory of biodiversity */

void start_time(void);
void prn_time(void);
void prn_total_time(void);
int binomialCoeff(int, int);

int main(int argc, char *argv[])
{


  int JM; /* argv[1], metacommunity size,*/
  double theta; /* argv[2], speciation rate*/
  int J; /* argv[3], sample size*/
  double *SM;
  double *x;
  double *S;
  
  FILE *out;
    
  JM = atoi(argv[1]);
  theta = atof(argv[2]);
  J = atoi (argv[3]);

  printf("Will calculate expected species abundance distribution in a sample of size %d from a metacommunity of size %d based on Eq. 6 and 8 in Alonso and McKane 2004. Note the probability of speciation per birth nu=%f\n", J, JM, theta/((double) JM-1) );

  out = fopen ("out.dat", "w");
  
  /* double SM[JM+1]; /\*metacommunity SAD*\/ */
  /* double x[JM+1]; /\*relative abundance in metacommunity*\/ */
  /* double S[J+1]; /\*sample SAD*\/ */

  SM = calloc(JM+1, sizeof(double));
  x = calloc(JM+1, sizeof(double));
  S = calloc(J+1, sizeof(double));


  int i, j;
  
  for(i=1; i<=JM; ++i) {
    SM[i] = (theta)/((double) i)*exp(lgamma(JM+1)+lgamma(JM+theta-i)-lgamma(JM+1-i)-lgamma(JM+theta));
    x[i] = ((double) i)/((double) JM);
    fprintf(out, "SM[%d]=\t%f\n", i, SM[i]);
  }
  for(j=1; j<=J; ++j) {
    S[j] = 0;
    for(i=1; i<=JM; ++i) {
      S[j] = S[j] +(SM[i]*exp(lgamma(J+1)-lgamma(J-j+1)-lgamma(j+1))*pow(x[i],((double)j))*pow(1-x[i],((double) J-j)));
     }
    fprintf(out, "S[%d]=\t%f\n", j, S[j]);
  }

}

