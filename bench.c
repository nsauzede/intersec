#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

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
int nsph = sizeof (spheres) / sizeof (spheres[0]) / 2;

v3 boxes[] = {
	{ -1, -1, -1},
	{ 2, 2, 2},
};

int main( int argc, char *argv[])
{
	v3 p;
	v3 e = {7,0,5};
	v3 v = {-1.0,0.0,-1.0};
	double t;
	int ni;

#if 0
	double *cs = spheres[0];
	double sr = spheres[1][0];
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
{
	int i, num = 2000000;
	printf( "generating %d rand spheres\n", num);
	v3 *spheres2 = malloc( num * 2 * sizeof (v3));
	memset( spheres2, 0, num * 2 * sizeof (v3));
	int *res = malloc( num * sizeof( int));
	double *t = malloc( num * sizeof( double));
	srand( time( 0));
	double old_tm, tm;
	struct timeval old_tv, tv;
	for (i = 0; i < num; i++)
	{
		v3 v;
		v[0] = 0.01 * (rand() % 100);
		v[1] = 0.01 * (rand() % 100);
		v[2] = 0.01 * (rand() % 100);
		sum3( spheres2[i * 2], spheres[0], v);
		spheres2[i * 2 + 1][0] = spheres[1][0] + 0.01 * (rand() % 100);
	}
	printf( "intersecting %d rand spheres\n", num);
	gettimeofday( &old_tv, 0);
	for (i = 0; i < num; i++)
	{
		res[i] = intersec_sphere( spheres2[i * 2], spheres2[i * 2 + 1][0], e, v, &t[i], 0);
	}
	gettimeofday( &tv, 0);
	old_tm = (double)old_tv.tv_sec * 1000000 + old_tv.tv_usec;
	tm = (double)tv.tv_sec * 1000000 + tv.tv_usec;
	double del = tm - old_tm;
	double interss = 1000000 * del / num;
	printf( "delay=%f inter/s=%f\n", del, interss);
	printf( "t=%f\n", t[num - 1]);
}
#else
	double *lower = boxes[0];
	double *upper = boxes[1];
	double n;
	n = norm3( v);
	div3( v, n);
	double t2;
	v3 p2;
	ni = intersec_box( lower, upper, e, v, &t, &t2);
	sum3( p, e, mult3( copy3( p, v), t));
	sum3( p2, e, mult3( copy3( p2, v), t2));
	disp3( "e", e);
	disp3( "v", v);
	disp3( "lower", lower);
	disp3( "upper", lower);
	printf( "ni=%d t=%f t2=%f\n", ni, t, t2);
	disp3( "p", p);
	disp3( "p2", p2);
{
	int i, num = 2000000;
	printf( "generating %d rand boxes\n", num);
	v3 *boxes2 = malloc( num * 2 * sizeof (v3));
	memset( boxes2, 0, num * 2 * sizeof (v3));
	int *res = malloc( num * sizeof( int));
	double *t = malloc( num * sizeof( double));
	srand( time( 0));
	double old_tm, tm;
	struct timeval old_tv, tv;
	for (i = 0; i < num; i++)
	{
		v3 v;
		v[0] = 0.01 * (rand() % 100);
		v[1] = 0.01 * (rand() % 100);
		v[2] = 0.01 * (rand() % 100);
		sum3( boxes2[i * 2], boxes[0], v);
		boxes2[i * 2 + 1][0] = boxes[1][0] + 0.01 * (rand() % 100);
	}
	printf( "intersecting %d rand boxes\n", num);
	gettimeofday( &old_tv, 0);
	for (i = 0; i < num; i++)
	{
		res[i] = intersec_box( boxes2[i * 2], boxes2[i * 2 + 1], e, v, &t[i], 0);
	}
	gettimeofday( &tv, 0);
	old_tm = (double)old_tv.tv_sec * 1000000 + old_tv.tv_usec;
	tm = (double)tv.tv_sec * 1000000 + tv.tv_usec;
	double del = tm - old_tm;
	double interss = 1000000 * del / num;
	printf( "delay=%f inter/s=%f\n", del, interss);
	printf( "t=%f\n", t[num - 1]);
}
#endif

	return 0;
}
