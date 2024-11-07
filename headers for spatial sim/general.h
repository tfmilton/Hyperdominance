# define MAX_FLNAME_SIZE 200
# define MAX_STR_SIZE 200
# define TRUE 1
# define FALSE 0
# define PI 3.14159265358979323846

int is_x_lower_than_y(void*, void*);

int myround(double);

void modular(int, int, int*);

void modular_dbl(double, double, double*);

double abs_dbl(double);

double get_rand_unit();

int get_rand_integ_intvl(int, int);

double get_rand_real_intvl(double, double);

double area_ring(double, double);

void descending_order(int*, int*, int);

double deg_to_radians(double);

double frac_circle_in_region(double, double, double, double, double);

double frac_circle_in_region_faster(double, double, double, double, double);

double max(double, double);

int intmax3(int, int, int);

int intmax4(int, int, int, int);

int intmax_array(int*, int);

void printout_2d_dist(double**, int, double, double, double, double, char*, char*, char*);

double sample_var(int, double*, double*, double*, double*);

void get_90perc_CL(int num_samples, double *samples, double *point_mean, double *point_lower_CL, double *point_upper_CL, double *p_sum_of_samples, double *p_sum_of_sqrd_samples);

double get_median(int num_samples, double *samples);

void get_quantile_boundaries(int, double*, double, double*, double*, double*);

double get_mean(int, double*);

void delete_lowest(int, double*);

double get_lowest(int, double*);

void delete_highest(int, double*);

double get_highest(int, double*);

double correlation_coefficient(int num_samples, double *X, double *Y);

int conv_rad_to_deg(double);

double conv_deg_to_rad(int);

int mymax(int*, int);
