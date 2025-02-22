/* standard header files */
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>

/* Based on the algorithm described in: Rosindell, J., Wong, Y. & Etienne, R.S. (2008). A coalescence approach to spatial neutral ecology. */


int get_rand_integ_intvl(int, int);
double get_rand_unit(void);
void start_time(void);
void prn_time(void);
void prn_total_time(void);

/*  defining function to sort spabund  */
int comparator(const void *p1, const void *p2)
{
  const int *dp1 = (const int *) p1;
  const int *dp2 = (const int *) p2;
  return (*dp1 > *dp2) - (*dp1 < *dp2);
}


int main(int argc, char *argv[])
{

  int* spabund; 
  int* metacomm;
  time_t t; /*declares a time variable*/

  int N; /* argv[1], metacommunity size*/
  double nu; /* argv[2], speciation rate*/
  int num_runs; /*argv [3], number of simulations */
  char* outfile = argv[4]; /*name of output files as an argument -- will be the slurm ID*/
  int AR_ID; /*argv[5], array task ID to set the seed */

  FILE *MMout;
    
  N = atoi(argv[1]);
  nu = atof(argv[2]);
  num_runs = atoi(argv[3]);
  AR_ID = atoi(argv[5]);

  /* set seed by array ID number*/
  srand(AR_ID);

  /* Create arrays */
  spabund=(int*)calloc(N+1,sizeof(int));
  metacomm=(int*)calloc(N+1,sizeof(int));

  /* Set up spabund and metacomm arrays */
  int i, j;
  for (i = 1; i <= N; i++)
    spabund [i] = 0;
  for (j = 1; j <= N; j++)
    metacomm [j] = 1;
        
  start_time();
  
  MMout = fopen(outfile, "w");
  int numzero; /* will hold the number of metacomm slots that are currently 0*/
  int whichslot; /* which nonzero patch is sampled */
  int parent;
  double toss; /* sampling to see if speciation or not */
  int numsp;
  int curr_num_ind;
  int r;
  for (r=1; r<=num_runs; ++r) {

    int i, j;
  for (i = 1; i <= N; i++)
    spabund [i] = 0;
  for (j = 1; j <= N; j++)
    metacomm [j] = 1;
    numzero=0;
    numsp=0;
    curr_num_ind=0;
    while(numzero < N-1) {
      whichslot = get_rand_integ_intvl(1, N-numzero); /* will keep metacomm organized so has all non-zero slots first, then zero slots*/
      toss = get_rand_unit();
      if (toss < nu) {
	numsp = numsp+1;
	spabund[numsp] = metacomm[whichslot];
	curr_num_ind = curr_num_ind + metacomm[whichslot];
	for (i=whichslot; i < N - numzero; ++i) {
	  metacomm[i] = metacomm[i+1];
	}
	metacomm[N-numzero] = 0;
	numzero = numzero + 1;
      }
      else { /*toss > nu*/
	parent = get_rand_integ_intvl(1, N);
	if(metacomm[parent]!=0){
	  metacomm[parent] = metacomm[parent] + metacomm[whichslot];
	  for (i=whichslot; i < N - numzero; ++i) {
	    metacomm[i] = metacomm[i+1];
	  }
	  metacomm[N-numzero]=0;
	  numzero = numzero +1;
	}
      }
    }
    numsp = numsp+1;
    spabund[0] = 0;
    spabund[numsp] = N-curr_num_ind;
    qsort(spabund, numsp+1, sizeof(int), comparator);  
    fprintf(MMout, "%d", numsp);
    for (i=numsp;i>=1; --i)
      fprintf(MMout, ",%d", spabund[i]);
    fprintf(MMout,"\n");
  }
    prn_total_time();
    fclose(MMout);
    return 0;
}

int get_rand_integ_intvl(int x, int y)
{

  double toss;

  toss = get_rand_unit();
  
  while (toss==0) toss=get_rand_unit();

  return (int) ceil(toss*(y-x+1))-1 + x;

}

double get_rand_unit()
{

  double toss;

  toss = (double) rand();
  toss = toss/((double) RAND_MAX);
  
  return toss;

}

# define MAXSTRING 100

typedef struct {
  clock_t begin_clock, save_clock;
  time_t begin_time, save_time;
} time_keeper;

static time_keeper tk;

void start_time(void)
{
  tk.begin_clock = tk.save_clock = clock();
  tk.begin_time = tk.save_time = time(NULL);
}

void prn_time(void)
{
  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, user_time, real_time;
  
  user_time = (clock() - tk.save_clock) / clocks_per_second;
  real_time = difftime(time(NULL), tk.save_time);
  tk.save_clock = clock();
  tk.save_time = time(NULL);

  /* print the values found, and do it neatly */

  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf("SINCE LAST PRN_TIME:\n");
  printf("%s%*.1f%s%f%s%f%s\n%s%*.1f%s%f%s%f%s\n", 
	 "User time: ", field_width, user_time, " seconds, ", user_time/60, " minutes, ", user_time/3600, " hours ", 
	 "Real time: ", field_width, real_time, " seconds, ", real_time/60, " minutes, ", real_time/3600, " hours ");

}

void prn_total_time(void)
{
  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, user_time, real_time;
  
  user_time = (clock() - tk.begin_clock) / clocks_per_second;
  real_time = difftime(time(NULL), tk.begin_time);

  /* print the values found, and do it neatly */

  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf("SINCE START:\n");
  printf("%s%*.1f%s%f%s%f%s\n%s%*.1f%s%f%s%f%s\n\n", 
	 "User time: ", field_width, user_time, " seconds, ", user_time/60, " minutes, ", user_time/3600, " hours ", 
	 "Real time: ", field_width, real_time, " seconds, ", real_time/60, " minutes, ", real_time/3600, " hours ");

}



