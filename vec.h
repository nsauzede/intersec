
#define SMALL 0.001

typedef double v3[3];

double* cross3( v3 dest, const v3 src1, const v3 src2);
double* sum3( v3 dest, const v3 src1, const v3 src2);
double* diff3( v3 dest, const v3 src1, const v3 src2);
double* div3( v3 n, double nn);
double* mult3( v3 n, double nn);
double* add3( v3 n, double nn);
double dot3( v3 p, v3 n);
double norm3( v3 n);
int intersec_plane( v3 p0, v3 p1, v3 p2, v3 l0, v3 l, double *pt);
int solvetri( double a, double b, double c, double *t1, double *t2);
int intersec_sphere( v3 cs, double sr, v3 e, v3 v, double *tmin, double *tmax);
