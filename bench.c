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

#include "vec.h"

v3 spheres[] = {
	{1, 0, 1},
	{2},
};
int nsph = sizeof (spheres) / sizeof (spheres[0]);

int main( int argc, char *argv[])
{
	v3 p;
	v3 e = {7,0,5};
	v3 v = {-1.0,0.0,-1.0};
	double *cs = spheres[0];
	double sr = spheres[1][0];
	double t;
	int ni;

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
	diff3( r, v, mult3( copy3( r, nn), 2 * dot));
	div3( r, norm3( r));
	disp3( "e", e);
	disp3( "v", v);
	disp3( "cs", cs);
	printf( "ni=%d t=%f\n", ni, t);
	disp3( "p", p);
	disp3( "nn", nn);
	printf( "dot=%f\n", dot);
	disp3( "r", r);

	return 0;
}
