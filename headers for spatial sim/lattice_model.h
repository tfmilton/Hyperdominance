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
