#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>

#ifdef WIN32
#include <windows.h>
#define MAP_FAILED ((HANDLE)-1)
#else
#include <sys/mman.h>
#endif

#define W 60
#define H 60
#define SMALL	0.001
#define BIG		10.0

#ifdef DEBUG
#define dprintf(...) do{printf(__VA_ARGS__);}while(0)
#else
#define dprintf(...) do{}while(0)
#endif

#if 1
#include "vec.h"
#else
int solvetri( double a, double b, double c, double *t1, double *t2)
{
	int result = 0;
	double d, sd;
	d = b * b - 4 * a * c;
	sd = sqrt( d);
	
	if (d > 0)
	{
		*t1 = (-b - sd) / 2 / a;
		*t2 = (-b + sd) / 2 / a;
		result = 2;
	}
	else if (d == 0)
	{
		*t1 = -b / 2 / a;
		result = 1;
	}
	else
	{
		result = 0;
	}
	return result;
}

int intersec_sphere( double cx, double cy, double cz, double sr, double ex, double ey, double ez, double vx, double vy, double vz, double *t)
{
	int result = 0;
	double a, b, c;
	double t1, t2;
	
	double tx, ty, tz;
	double n;
	n = sr * sr;
	tx = ex - cx;
	ty = ey - cy;
	tz = ez - cz;
	a = vx * vx + vy * vy + vz * vz;
	b = 2 * ( vx * tx + vy * ty + vz * tz);
	c = (tx * tx + ty * ty + tz * tz) - n;
	int sol;
	sol = solvetri( a, b, c, &t1, &t2);
	if (sol == 2)
	{
//		printf( "two solutions : %f and %f\n", t1, t2);
		if (fabs( t1) > fabs( t2))
			t1 = t2;
	}
	else if (sol == 1)
	{
//		printf( "one solution : %f\n", t1);
	}
	else
	{
//		printf( "no solutions\n");
	}
	if (sol && (t1 > SMALL))
	{
		result = 1;
//		printf( "returning %f\n", t1);
		if (t)
		{
			*t = t1;
		}
	}

	return result;
}
#endif

typedef struct {
	double cx, cy, cz, sr;
	double r, g, b;
} sphere_t;

sphere_t spheres[] = {
	{ .cx = 1.0, .cy = 0.0, .cz = 1.0, .sr = 2.0, .r = 0.0, .g = 0.0, .b = 1.0 },
};
int nsph = sizeof (spheres) / sizeof (spheres[0]);

int main( int argc, char *argv[])
{
	v3 e = {7,0,5};
	v3 v = {-1.0,0.0,-1.0};
	v3 p = {9,0,0};
	double t;
	int ni;
	double sr;
	v3 cs;

	cs[0] = spheres[0].cx; cs[1] = spheres[0].cy; cs[2] = spheres[0].cz; sr = spheres[0].sr;
	double n;
	n = norm3( v);
	div3( v, n);
	ni = intersec_sphere( cs, sr, e, v, &t, 0);
	sum3( p, e, mult3( copy3( p, v), t));
	v3 r;
	v3 nn;
	double dot;
	diff3( nn, p, cs);
	div3( nn, norm3( nn));
	dot = dot3( nn, v);
	diff3( r, v, mult3( nn, 2 * dot));
	div3( r, norm3( r));
	disp3( "e", e);
	disp3( "v", v);
	disp3( "cs", cs);
	printf( "ni=%d t=%f\n", ni, t);
	disp3( "p", p);
	disp3( "nn", nn);
	printf( "dot=%f\n", dot);
	disp3( "r", r);
//	printf( "vx=%f vy=%f vz=%f ni=%d t=%f x=%f y=%f z=%f nx=%f ny=%f nz=%f dot=%f rx=%f ry=%f rz=%f\n", vx, vy, vz, ni, t, x, y, z, nx, ny, nz, dot, rx, ry, rz);

	return 0;
}
