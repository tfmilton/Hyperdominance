/* standard header files */
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>
# include <string.h>

/* time keeper stuff*/
void start_time(void);
void prn_time(void);
void prn_total_time(void);

/* general functions*/
int comparator(const void *, const void *);
int myround(double);
void modular(int, int, int*);
double get_rand_unit();
int get_rand_integ_intvl(int, int);
double get_rand_real_intvl(double, double);

/* lattice model functions*/
void convert_to_coords(int, int, int*);
void choose_random_location(int, int, int*);
void choose_randloc_atdist(double, int*);
void choose_random_nearest_neighbor(int*, int, int, int*);
double gaussian(double, double*);
double gaussianxy(double, double, double*);
double fat_tailed(double, double*);
double fat_tailedxy(double, double, double*);
double *discrete_normalized_cummulative(double f(double, double*), double*, double, int, int*);  
void select_randloc_frm_cum_disp(int*, int, int, double*, int, int*);
double distance(double, double, double, double);

/* macrodynamic functions */
int *speciate(int**, int*, int*, int*, int*, int*);
void replace(int**, int, int, int*, int*, int*, double*, int, int*);


/* main function#*/

int main(int argc, char *argv[])
{

  /* input parameters */
  int size_x, size_y, maxtime, start_out_time;
  int intev_btwn_specrich_time, intev_btwn_specabund_time, intev_btwn_landscape_time;
  int num_runs; /*argv [15], number of simulations */
  double nu;
  double sigma, p, u, parameters[2], kern_cutoff, max_dist;

  /* variables for initial conditions of landscape*/
  int numspec, *ABUND, splab_array_size = 100, larg_splab_inuse, **LANDSCAPE;

  /* other variables */
  int t, ind1_coords[2], ind2_coords[2], loc_kern_cutoff, r, testvar, AR_ID;
  int i,j, xmody[2], num_mutants_per_timestep, num_specevent=0;
  double toss, *disc_norm_cum, doublemask;
  char landscape_outfile[200], specrich_outfile[200], specabund_outfile[200],all_abund_outfile[200] ; 
  char  curr_time[30], countr[30], inID[30];
  FILE *specrich_ofp, *specabund_ofp, *landscape_ofp, *MMout;
  

  /* argv[0] = program name */
  /* argv[1] = lattice size in x direction */
  /* argv[2] = lattice size in y direction */
  /* argv[3] = number of time steps */
  /* argv[4] = time step to start ouputting lanscape and spec abund*/
  /* argv[5] = interval of time to wait between outputting spec richness */
  /* argv[6] = interval of time to wait between outputting spec abund*/
  /* argv[7] = interval of time to wait between outputting landscape*/
  /* argv[8] = speciation rate nu */
  /* argv[9] = type of dispersal (gaussian or fat_tailed or global) */
  /* argv[10] = guassian sigma or fat_tailed u*/
  /* argv[11] = fat_tailed p*/
  /* argv[12] = kernel cutoff */
  /* argv[13] = beginning of output file names*/
  /* argv[14] = max dispersal distance*/
  /* argv[15] = number of runs*/
  /* argv[16] = array task ID */

  /* start program timer*/
  start_time();

  /* seed random number generator */
    srand(atoi(argv[16]));
    AR_ID = atoi(argv[16]);

    printf("The array ID is %d.\n", AR_ID);

  /* get input variables */
  size_x = atoi(argv[1]);
  size_y = atoi(argv[2]);
  maxtime = atoi(argv[3]);
  start_out_time = atoi(argv[4]);
  intev_btwn_specrich_time = atoi(argv[5]);
  intev_btwn_specabund_time = atoi(argv[6]);
  intev_btwn_landscape_time = atoi(argv[7]);
  nu = (double) atof(argv[8]);
  num_runs = atoi (argv[15]);
  

  /* prepare dispersal kernel */
  printf("\n\nThe dispersal kernel is %s.\n", argv[9]);
  if (strcmp(argv[9], "gaussian") == 0) { /* if argv[9] == "gaussian"*/
    kern_cutoff = (double) atof(argv[12]);
    sigma = (double) atof(argv[10]);
    parameters[0] = sigma;
    max_dist = (double) atof(argv[14]);
    disc_norm_cum = discrete_normalized_cummulative(gaussian, parameters, kern_cutoff, max_dist, &loc_kern_cutoff);
  }
  else if (strcmp(argv[9], "fat_tailed") == 0) {
    kern_cutoff = (double) atof(argv[12]);
    u = (double) atof(argv[10]);
    p = (double) atof(argv[11]);
    parameters[0] = u;
    parameters[1] = p;
    max_dist = (double) atof(argv[14]);
    disc_norm_cum = discrete_normalized_cummulative(fat_tailed, parameters, kern_cutoff, max_dist, &loc_kern_cutoff);
  }
  else if (strcmp(argv[9], "global") ==0) {
    disc_norm_cum = calloc(1, sizeof(double));
    disc_norm_cum[0] = 2;
  }
  else if (strcmp(argv[9], "nearest_neighbor") ==0) {
    disc_norm_cum = calloc(1, sizeof(double));
    disc_norm_cum[0] = 3;
  }
  else {
    printf("No dispersal kernel chosen!!\n");
    exit(1);
  }
  /* if you change this remember to change at end up dynamic loop */
  /* Move these down to within the r runs loop*/
  /* strcpy(specrich_outfile, argv[13]); */
  /* strcpy(specabund_outfile, argv[13]); */
  /* strcpy(landscape_outfile, argv[13]); */

  /* for array, all file names need to be based on array task ID*/
  strcpy(all_abund_outfile, argv[13]);
  sprintf(inID, "%d", AR_ID);
  strcpy(all_abund_outfile, "SADs_");
  strcat (all_abund_outfile, inID);
  strcat (all_abund_outfile, ".dat");
   MMout = fopen(all_abund_outfile, "w");

  printf("The lattice size is %d x %d.\nThe number of time steps the model will be run is %d.\nRecording of the landscape and species abundances will begin at %d.\nThe species richness will be recorded every %d.\nThe number of speciation events per individual per time step is %f.\nThe cutoff for calculating the dispersal kernel is %f.\nThe beginning of the landscape filename is %s.\n", size_x, size_y, maxtime, start_out_time, intev_btwn_specrich_time, nu, kern_cutoff, landscape_outfile);

  /* allocate necessary arrays */  
  ABUND = calloc(splab_array_size+1, sizeof(int));

  LANDSCAPE = calloc(size_x, sizeof(int *));
  for (i = 0; i<size_x; ++i)
    LANDSCAPE[i] = calloc(size_y, sizeof(int));

  /*Start loop here for many runs */
  for (r=1; r<=num_runs; r++){
  
  /* if you change this remember to change at end up dynamic loop */
  strcpy(specrich_outfile, argv[13]);
  strcpy(specabund_outfile, argv[13]);
  strcpy(landscape_outfile, argv[13]);

  
  /* start landscape with all one species */
  numspec = 1;
  larg_splab_inuse = 1;
  ABUND[1] = size_x*size_y;
  for (i=0; i<size_x; ++i)
    for (j=0; j<size_y; ++j)
      LANDSCAPE[i][j] = 1;

      testvar = r;
      sprintf(countr, "%d", testvar);
      sprintf(inID, "%d", AR_ID);
      strcat(specrich_outfile, "sr_");
      strcat(specrich_outfile, countr);
      strcat(specrich_outfile, "_");
      strcat(specrich_outfile, inID);
      strcat(specrich_outfile, ".dat");
      specrich_ofp = fopen(specrich_outfile, "w");

  num_mutants_per_timestep = myround(nu*size_x*size_y);
  printf("The number of mutants that will be created per time step is %d.\n", num_mutants_per_timestep); 


  for (t=1; t<=maxtime; ++t) {
    
  
    for (i=1; i<= size_x*size_y; ++i) {
     
     
      toss=get_rand_unit();
      if (toss < nu) {
      	ABUND = speciate(LANDSCAPE, ind1_coords, ABUND, &numspec,
			 &splab_array_size, &larg_splab_inuse);
    	num_specevent = num_specevent+1;
      }
      else{
      choose_random_location(size_x, size_y, ind1_coords);
      replace(LANDSCAPE, size_x, size_y, ind1_coords, ABUND, &numspec, disc_norm_cum, loc_kern_cutoff, ind2_coords);
      }
    }

    /* the output */
    modular(t, intev_btwn_specrich_time, xmody);
    if (xmody[1] == 0) {
      fprintf(specrich_ofp, "%d\t%d\n", t, numspec);
    }

    modular(t, intev_btwn_specabund_time, xmody);
    if (t >= start_out_time && xmody[1] == 0) {
      doublemask = t;
      testvar = r;
      sprintf(curr_time, "%.1e", doublemask);
      sprintf(countr, "%d", testvar);
      sprintf(inID, "%d", AR_ID);
      strcat(specabund_outfile, "sa_");
      strcat(specabund_outfile, curr_time);
      strcat(specabund_outfile, "_");
      strcat(specabund_outfile, countr);
      strcat(specabund_outfile, "_");
      strcat(specabund_outfile, inID);
      strcat(specabund_outfile, ".dat");
      specabund_ofp = fopen(specabund_outfile, "w");
      for (i=1; i<=larg_splab_inuse; ++i)
	fprintf(specabund_ofp,"  %d      %d\n", i, ABUND[i]);
      fclose(specabund_ofp);
      strcpy(specabund_outfile, argv[13]);
    }

    modular(t, intev_btwn_landscape_time, xmody);
    if (t >= start_out_time && xmody[1] == 0) {
      doublemask = t;
      testvar = r;
      sprintf(curr_time, "%.1e", doublemask);
      sprintf(countr, "%d", testvar);
      sprintf(inID, "%d", AR_ID);
      strcat(landscape_outfile, "lscp_");
      strcat(landscape_outfile, curr_time);
      strcat(landscape_outfile, "_");
      strcat(landscape_outfile, countr);
      strcat(landscape_outfile, "_");
      strcat(landscape_outfile, inID);
      strcat(landscape_outfile, ".dat");
      landscape_ofp = fopen(landscape_outfile, "w");
      for (i=0; i<size_x; ++i) {
	for (j=0; j<size_y; ++j)
	  fprintf(landscape_ofp," %d  ", LANDSCAPE[i][j]);
	fprintf(landscape_ofp,"\n");
      }
      fclose(landscape_ofp);
      strcpy(landscape_outfile, argv[13]);

    }
    /* printf ("time is %d \n", t); */
  }
  /*Here ABUND it is the final Abundance array after all timesteps*/
  /*So sort ABUND, then print it to MMout, for the length of ABUND, new line, then loop back for next run*/ 
  qsort (ABUND, larg_splab_inuse+1, sizeof(int), comparator);
  for (i=larg_splab_inuse; i>=1; --i)
    fprintf(MMout,"%d,", ABUND[i]);
  fprintf(MMout, "\n");
  }
  
  printf("The number of speciation events that occurred is %d.\n\n\n", num_specevent);

  prn_time();
  printf("\n\n");
}


/* time_keeper stuff */
   
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

/* various general functions*/

/* defining function to sort spabund (ascending order) */
int comparator(const void *p1, const void *p2)
{
  const int *dp1 = (const int *) p1;
  const int *dp2 = (const int *) p2;

  return (*dp1 > *dp2) - (*dp1 < *dp2);
}

/* get_rand_integ_intvl returns a random integer in the interval [x,y], that
   is, including x and y as possibilities. */
/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this
   function! --srand(time(NULL)) */

int get_rand_integ_intvl(int x, int y)
{

  double toss;

  toss = get_rand_unit();
  /*  printf("toss=%f, y-x+1=%d, toss*(y-x+1)=%f, ceil(toss*(y-x+1))=%f, returning %d.\n", toss, y-x+1, toss*(y-x+1), ceil(toss*(y-x+1)), (int) ceil(toss*(y-x+1))-1 + x);*/

  /* technically should throw out zeros, because should choose a particular integer based on a region non-inclusive of the one before that, even when choosing the lowest value in the integer range */
  while (toss==0) toss=get_rand_unit();

  return (int) ceil(toss*(y-x+1))-1 + x;

  

}

/* get_rand_real_intvl returns a random real number in the interval
   [x,y], that is, including x and y as possibilities. */
/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this
   function! --srand(time(NULL)) */
double get_rand_real_intvl(double x, double y)
{

  double toss;

  toss = rand();
  return ((toss/RAND_MAX)*(y-x) + x);

}

/* get_rand_unit returns a real random number in the interval [0,1],*/
/*   that is, including 0 and 1 as possiblities.  rand() returns a */
/*   random integer between 0 and RAND_MAX.  Hence there is a discrete */
/*   set of real numbers that this function can return--this set */
/*   contains about 2 billion (2X10^9) numbers.  */

/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this function! srand(time(NULL)) */
double get_rand_unit()
{

  double toss;

  toss = (double) rand();
  toss = toss/((double) RAND_MAX);
  
  return toss;

}

/*note: this returns first the answer to x/y, then the remainder*/
void modular(int x, int y, int *xmody)
{

  xmody[0] = floor(((double) x)/((double) y));
  xmody[1] = x - xmody[0]*y;
  
}

/* take care not to use this on doubles which are out of the INT_MAX range */
/* I've read that it is more efficient to use a #def command instead of a function here, so if wanting to find a way to speed things up, may want to do that instead.  I'm not sure how much it really matters.  Most of my lattice model functions don't use quite this anyway, which is were speed would matter. */
int myround(double x)
{
  int intcast;

  if (fabs(((double) ceil(x)) - x) > 0.5) 
    return (int) floor(x);
  else
    return (int) ceil(x);
}


/* lattice model stuff*/

# define PI 3.14159265358979323846


/* Picture a lattice with y columns.  Label each site by one number.
   Site (0,0) is labeled by 1, (0,1) by 2, (0,2) by 3, ..., (1,0) by
   y+1, etc.  So the idea is to label the lattice sites by a number
   that increases as you go in the y direction on the lattice until
   the end of the row, and then keeps increasing on the next row.
   This function will convert from the one number label to the
   coords (x,y)*/
   
void convert_to_coords(int label, int size_y, int *coords) {

  modular(label-1,size_y, coords);

}

int convert_to_label(int size_y, int *coords) {
  
  int x, y;

  x=coords[0];
  y=coords[1];

  return x*size_y + y+1;

}


/* choose_random_location chooses at random one of the sites on the
   lattice of size size_x by size_y.  The coordinates of the randomly chosen
   site are returned in the ind1_coords array.  Only one random number
   generation is used by using the convert_to_coords function.  The
   coordinates go from 0 to size_x-1 and 0 to size_y-1. */ 
/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this function! --srand(time(NULL)) */
void choose_random_location(int size_x, int size_y, int *coords)
{

  int chosen_location;
  
  chosen_location = get_rand_integ_intvl(1, size_x*size_y);
  convert_to_coords(chosen_location, size_y, coords);

}


/* choose_randloc_atdist chooses a random location on the lattice
   approximately a distance r from the origin.  It does this by
   picking an angle at random and rounding the coordinats xcos(angle)
   and xsin(angle).  This rounding picks the closest point on the
   lattice to these real coordinates, and returns the integer
   coordinates of that closest point in coords.*/

void choose_randloc_atdist(double r, int *coords)  
{
  double angle;

  angle=get_rand_real_intvl(0,2*PI);
  coords[0] = myround(r*cos(angle));
  coords[1] = myround(r*sin(angle));

  /*  printf("angle chosen is %f radians\n", angle); */
  /*  printf("r*cos(angle)=%f, coords[0]=%d, r*sin(angle)=%f, coords[1]=%d\n", r*cos(angle), coords[0], r*sin(angle), coords[1]);*/

}


void choose_random_nearest_neighbor(int *ind1_coords, int size_x, int size_y, int *ind2_coords)
{

  int chosen_location;
  
  
  chosen_location = get_rand_integ_intvl(1, 4);
  if (chosen_location == 1) {
    ind2_coords[0] = ind1_coords[0] + 0;
    ind2_coords[1] = ind1_coords[1] + 1;
  }
  if (chosen_location == 2) {
    ind2_coords[0] = ind1_coords[0] + 1;
    ind2_coords[1] = ind1_coords[1] + 0;
  }
  if (chosen_location == 3) {
    ind2_coords[0] = ind1_coords[0] + 0;
    ind2_coords[1] = ind1_coords[1] + -1;
  }
  if (chosen_location == 4) {
    ind2_coords[0] = ind1_coords[0] + -1;
    ind2_coords[1] = ind1_coords[1] + 0;
  }
  while (ind2_coords[0] >= size_x || ind2_coords[0] < 0 ||
	 ind2_coords[1] >= size_y || ind2_coords[1] < 0) {
      chosen_location = get_rand_integ_intvl(1, 4);
      chosen_location = get_rand_integ_intvl(1, 4);
      if (chosen_location == 1) {
	ind2_coords[0] = ind1_coords[0] + 0;
	ind2_coords[1] = ind1_coords[1] + 1;
      }
      if (chosen_location == 2) {
	ind2_coords[0] = ind1_coords[0] + 1;
	ind2_coords[1] = ind1_coords[1] + 0;
      }
      if (chosen_location == 3) {
	ind2_coords[0] = ind1_coords[0] + 0;
	ind2_coords[1] = ind1_coords[1] + -1;
      }
      if (chosen_location == 4) {
	ind2_coords[0] = ind1_coords[0] + -1;
	ind2_coords[1] = ind1_coords[1] + 0;
      }
  }
}


/*guassian(r, parameters) returns the gaussian probability for a seed
  to land in a ring between r and r+dr.  It is 2*pi*r times the one
  dimensional guassian probability distribution. */
double gaussian(double r, double *parameters)
{
  double sigma;

  sigma = parameters[0];

  return ((r/pow(sigma,2)) * (exp(-pow(r,2)/(2*pow(sigma,2)))));

}

double gaussianxy(double x, double y, double *parameters)
{
  double sigma;

  sigma = parameters[0];

  /*  printf("sigma=%f.\n", sigma);
  printf("pow(sigma,2)=%f.\n", pow(sigma,2));
  printf("(1/(2*PI*pow(sigma,2)))=%f.\n", (1/(2*PI*pow(sigma,2))));
  printf("distance(x,y,0,0)=%f.\n", distance(x,y,0,0));
  printf("(-pow(distance(x,y,0,0),2)/(2*pow(sigma,2)))=%f.\n", (-pow(distance(x,y,0,0),2)/(2*pow(sigma,2))) );
  printf("exp(-pow(distance(x,y,0,0),2)/(2*pow(sigma,2)))=%f.\n", exp(-pow(distance(x,y,0,0),2)/(2*pow(sigma,2))) );*/

  return (1/(2*PI*pow(sigma,2)))*(exp(-pow(distance(x,y,0,0),2)/(2*pow(sigma,2))));
}


/* fat_tailed(r, parameters) returns the probability for a seed to
   land in a ring between r and r+dr according to a fat-tailed
   dispersal kernel developped in Clark et al. 1999 Ecology
   80:1475. It is 2*pi*r times the one dimensional probability
   distribution. */
double fat_tailed(double r, double *parameters)
{
  double u, p;

  u = parameters[0];
  p = parameters[1];

  return ((2*p*r)/(u*pow(1 + (pow(r,2)/u),p+1)));

}

double fat_tailedxy(double x, double y, double *parameters)
{

  double u, p;

  u = parameters[0];
  p = parameters[1];

  /*  printf("x=%f, y=%f.\n", x, y);
  printf("distance(x,y,0,0)=%f.\n", distance(x,y,0,0));
  printf("pow(1 + (pow(distance(x,y,0,0),2)/u),p+1)=%f.\n", pow(1 + (pow(distance(x,y,0,0),2)/u),p+1) );
  printf("1/(PI*u*pow(1 + (pow(distance(x,y,0,0),2)/u),p+1))=%f.\n", 1/(PI*u*pow(1 + (pow(distance(x,y,0,0),2)/u),p+1)) );*/

  return (p/(PI*u*pow(1 + (pow(distance(x,y,0,0),2)/u),p+1)));

}

/* discrete_normalized_cummulative will return the discrete, normalized
   cummulative probability distribution for the dispersal kernel f(r,
   parameters) out to the point where f(r, parameters) is smaller than
   the cutoff.  It will also return the location of this cutoff point
   in loc_cutoff. */
double *discrete_normalized_cummulative(double f(double, double*), 
				    double *parameters, 
				    double cutoff, int max_disp_dist, 
					int *point_loc_cutoff)  
{ 

  int curr_arraysize=100, pastmax=0;
  double *disc_norm_cum;
  
  int i, j;

  disc_norm_cum = calloc(curr_arraysize, sizeof(double));
  printf("\nCalculating dispersal kernel...\n");
  fflush(NULL);

  i=1;/*note starting at 1 because not allowing self-replacement*/
  /* while we are either not past the maximum or the value to add is
     not less than the cutoff, keep filling array*/
  while (!pastmax || (1-disc_norm_cum[i-1] >= cutoff && i<= max_disp_dist)) {
    /* if haven't already passed the maximum, check if you just did. */
    if (!pastmax)
      if (f(i,parameters) <= f(i-1,parameters)) pastmax=1;
    /* fill if at point less than array size, otherwise reallocate */
    if (i==curr_arraysize) {
      curr_arraysize = curr_arraysize + 100;
      disc_norm_cum = realloc(disc_norm_cum,
			      curr_arraysize*sizeof(double)); 
    }
    disc_norm_cum[i] = disc_norm_cum[i-1];
    if (i==0) /*this not being used due to no self-replacement, but keep it here in case change that later*/
      for(j=0; j<=4; ++j) 
	disc_norm_cum[i] = disc_norm_cum[i] + 0.1*f(i+j*0.1, parameters);
    else
      for(j=0; j<=9; ++j) 
	disc_norm_cum[i] = disc_norm_cum[i] + 0.1*f(i-0.5+j*0.1, parameters);
    /*    printf("i=%d, disc_norm_cum=%f\n", i, disc_norm_cum[i]); */
    ++i;
  }
  
  /* return the location where the cutoff was met */
  *point_loc_cutoff = i-1;

  /* normalize the cummulative probabilities, making the largest
     probability be 1 */
  printf("Normalizing w/ loc_cutoff=%d, remaining prob=%f...\n", *point_loc_cutoff, 1-disc_norm_cum[*point_loc_cutoff]);
  for (i=1; i<=*point_loc_cutoff; ++i) {
    disc_norm_cum[i] = disc_norm_cum[i]/disc_norm_cum[*point_loc_cutoff];
    /*    printf("i=%d, disc_norm_cum=%f\n", i, disc_norm_cum[i]);*/
  }

  /* return the pointer to the discrete normalized cummulative
     probability array */
  return disc_norm_cum;

}

/* select_randloc_frm_cum_disp selects a random location according to
   the cummulative disperal probability for a distance r in
   disc_norm_cum.  It first selects a distance according to this
   cummulative distribution, and then calls a function that selects a
   lattice site approximately that distance away at random. */
/* In order to test if this way of picking sites according to the
   dispersal kernel works well, I */
void select_randloc_frm_cum_disp(int *ind1_coords, int size_x, int size_y, double *disc_norm_cum, int loc_cutoff, int *ind2_coords)
{

  int i;
  double toss;

  toss = get_rand_unit();
  i=1;/*note we don't look in the 0th because we are not allowing self-replacement and the array was created and normalized assuming no self-replacement*/
  while((toss > disc_norm_cum[i]) && (i<loc_cutoff)) {
    ++i;
  }
  choose_randloc_atdist(i, ind2_coords);
  ind2_coords[0] = ind2_coords[0] + ind1_coords[0];
  ind2_coords[1] = ind2_coords[1] + ind1_coords[1];
  while (ind2_coords[0] >= size_x || ind2_coords[0] < 0 ||
	 ind2_coords[1] >= size_y || ind2_coords[1] < 0) {
    toss = get_rand_unit();
    i=0;
    while((toss > disc_norm_cum[i]) && (i<loc_cutoff)) {
      ++i;
    }
    choose_randloc_atdist(i, ind2_coords);
    ind2_coords[0] = ind2_coords[0] + ind1_coords[0];
    ind2_coords[1] = ind2_coords[1] + ind1_coords[1];
  }
}

double distance(double x_1, double y_1, double x_2, double y_2)
{

  double x_dist, y_dist;

  x_dist = x_2-x_1;
  y_dist = y_2-y_1;

  /*  printf("x_dist=%f, pow(x_dist,2)=%f.\n", x_dist, pow(x_dist,2));
      printf("y_dist=%f, pow(y_dist,2)=%f.\n", y_dist, pow(y_dist,2));*/

  return sqrt(pow(x_dist,2) + pow(y_dist, 2));

}
  
void disperse(int *coords1, int size_x, int size_y, double *disc_norm_cum, int loc_cutoff, int *coords2)
{

  int i;

  if (disc_norm_cum[0] == 2) {
    choose_random_location(size_x, size_y, coords2);
    while ((coords2[0] == coords1[0]) && (coords2[1] ==
					  coords1[1])) /* no self-replacement*/
      choose_random_location(size_x, size_y, coords2);
  }
  else if (disc_norm_cum[0] == 3) 
    choose_random_nearest_neighbor(coords1, size_x, size_y, coords2);
  else 
    select_randloc_frm_cum_disp(coords1, size_x, size_y, disc_norm_cum, loc_cutoff, coords2);

  
}


/* macrodynamic functions*/

int *speciate(int **LANDSCAPE, int *ind1_coords, int *ABUND, int
	      *point_numspec, int *point_splab_array_size, int *point_larg_splab_inuse)
{

  int i;

  int species, *ABUND_TEMP, newspecies=1;

  /* decrease abundance of the species of ind1 by 1 */
  species = LANDSCAPE[ind1_coords[0]][ind1_coords[1]];
  ABUND[species] = ABUND[species] - 1;

  /* use the lowest species label not in use */
  i=1;
  while (ABUND[i] != 0 && i<=*point_larg_splab_inuse) ++i;
  newspecies = i;
  if (i>*point_larg_splab_inuse) *point_larg_splab_inuse = i;

  /* increase numspec if the first species didn't go extinct */
  if (ABUND[species] != 0) 
    *point_numspec = *point_numspec +1;

  /* put an individual of the new species where ind1 was*/
  LANDSCAPE[ind1_coords[0]][ind1_coords[1]] = newspecies;

  /* printf ("new species %d  at %d, %d\n", newspecies, ind1_coords[0], ind1_coords[1]); */


  /* allocate more space in the species abundance array if necessary */
  if (newspecies > *point_splab_array_size) {
    ABUND=realloc(ABUND, (*point_splab_array_size+101)*sizeof(int));
    *point_splab_array_size = *point_splab_array_size+100;
    /* initialize new part of abundance array because realloc doesn't do it. */
    for(i=newspecies; i<=*point_splab_array_size; ++i)
	  ABUND[i] = 0;
  }

  /* set the abundance of the newspecies to be 1 */
  ABUND[newspecies] = 1;

  /* return the pointer to the abundance array, because it may be a new pointer and the program outside doesn't seem to pick this up unless the value is actually returned. */
  return ABUND;

}



void replace(int **LANDSCAPE, int size_x, int size_y, int
	     *ind1_coords, int *ABUND, int *point_numspec, double
	     *disc_norm_cum, int loc_cutoff, int *ind2_coords)
{
  int species1, species2;

  int i;


  if (disc_norm_cum[0] == 2) {
    choose_random_location(size_x, size_y, ind2_coords);
    while ((ind2_coords[0] == ind1_coords[0]) && (ind2_coords[1] ==
	   ind1_coords[1])) 
      choose_random_location(size_x, size_y, ind2_coords);
  }
  else if (disc_norm_cum[0] == 3) 
    choose_random_nearest_neighbor(ind1_coords, size_x, size_y, ind2_coords);
  else 
    select_randloc_frm_cum_disp(ind1_coords, size_x, size_y, disc_norm_cum, loc_cutoff, ind2_coords);

 /* printf("ind1_coords = (%d, %d), ind2_coords = (%d,%d)\n", ind1_coords[0], ind1_coords[1], ind2_coords[0], ind2_coords[1]); */
  
  species1 = LANDSCAPE[ind1_coords[0]][ind1_coords[1]];
  species2 = LANDSCAPE[ind2_coords[0]][ind2_coords[1]];

   /* printf("species1 is %d, species2 is %d\n", species1, species2); */

  /* if the species identity of the first individual is different than
     the species identity of the individual chosen to give offspring
     to the first individual's site, decrease the abundance of species
     1 by one, increase the abundance of species 2 by one, and put an
     individual of species 1 at the first individual's site */
  if (species1 != species2) {
    ABUND[species1] = ABUND[species1] - 1;
    ABUND[species2] = ABUND[species2] + 1;
    LANDSCAPE[ind1_coords[0]][ind1_coords[1]] = species2;

    /* printf ("replacement: species 2 is %d at coords %d, %d\n", species2, ind1_coords[0], ind1_coords[1]); */
  }
  
  /* if species1 goes extinct, update the number of species.  note
     that another variable in the main program will keep track of the
     largest possible species label with a non-zero abundance, which is
     different than numspec because there may be intermediate labels
     not currently in use.  do things this way to avoid having to go
     through the entire lattice updating all of the species labels.  */
  if (ABUND[species1] == 0)
    *point_numspec = *point_numspec -1;

}






