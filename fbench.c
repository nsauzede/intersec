#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <string.h>
#include <sys/time.h>

#ifdef DEBUG
#define dprintf(...) do{printf(__VA_ARGS__);}while(0)
#else
#define dprintf(...) do{}while(0)
#endif

#include "fvec.h"

v3 boxes[] = {
	{ S(-1), S(-1), S(-1)},
	{ S(2), S(2), S(2)},
};

int main( int argc, char *argv[])
{
	v3 e = {S(7),S(0),S(5)};
	v3 v = {S(-1),S(0),S(-1)};

	scalar_t n;
	n = norm3( v);
	div3( v, n);

#if 1
#define NUM 2000000 
#else
#define NUM 1
#endif
	int i, num = NUM;
	printf( "generating %d rand boxes\n", num);
	v3 *boxes2 = malloc( num * 2 * sizeof (v3));
	memset( boxes2, 0, num * 2 * sizeof (v3));
	int *res = malloc( num * sizeof (int));
	scalar_t *t = malloc( num * sizeof (scalar_t));
	memset( t, 0, num * sizeof (scalar_t));
	srand( time( 0));
	double old_tm, tm;
	struct timeval old_tv, tv;
	for (i = 0; i < num; i++)
	{
		v3 tv = { S(0), S(0), S(0)};
#if 0
		tv[0] = S(0.001 * (rand() % 10));
		tv[1] = S(0.001 * (rand() % 10));
		tv[2] = S(0.001 * (rand() % 10));
#endif
		sum3( boxes2[i * 2], boxes[0], tv);
		sum3( boxes2[i * 2 + 1], boxes[1], tv);
	}
	printf( "intersecting %d rand boxes\n", num);
	gettimeofday( &old_tv, 0);
	for (i = 0; i < num; i++)
	{
		res[i] = intersec_box( boxes2[i * 2], boxes2[i * 2 + 1], e, v, &t[i], 0);
#if NUM == 1
		printf( "res=%d\n", res[i]);
#endif
	}
	gettimeofday( &tv, 0);
	old_tm = (double)old_tv.tv_sec * 1000000 + old_tv.tv_usec;
	tm = (double)tv.tv_sec * 1000000 + tv.tv_usec;
	double del = tm - old_tm;
	del += 1; // in case we get 0 us
	double interss = (double)num * 1000000 / del;
	char *unit = "";
	double mult = 1;
	if (interss > 1e9)
	{
		mult = 1e9;
		unit = "G";
	}
	else if (interss > 1e6)
	{
		mult = 1e6;
		unit = "M";
	}
	else if (interss > 1e3)
	{
		mult = 1e3;
		unit = "K";
	}
	interss /= mult;
	printf( "delay=%f rate=%lu %sinters/s\n", del, (unsigned long)interss, unit);
	printf( "t=%f\n", D(t[0]));

	return 0;
}
