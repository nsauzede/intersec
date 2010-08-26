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

typedef struct {
	double cx, cy, cz, sr;
	double r, g, b;
} sphere_t;

sphere_t spheres[] = {
	{ .cx = 1.0, .cy = 0.0, .cz = 1.0, .sr = 2.0, .r = 0.0, .g = 0.0, .b = 1.0 },
};
int nsph = sizeof (spheres) / sizeof (spheres[0]);

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

int main( int argc, char *argv[])
{
	double ex, ey, ez, vx, vy, vz;

	ex = 7; ey = 0; ez = 5; vx = -1.0; vy = 0.0; vz = -1.0;

	double x = 0.0, y = 0.0, z = 0.0;
	double t;

	int ni;
	double cx, cy, cz, sr;
	cx = spheres[0].cx; cy = spheres[0].cy; cz = spheres[0].cz; sr = spheres[0].sr;
	double n;
	n = sqrt( vx * vx + vy * vy + vz * vz);
	vx /= n;
	vy /= n;
	vz /= n;
	ni = intersec_sphere( cx, cy, cz, sr, ex, ey, ez, vx, vy, vz, &t);
	x = ex + t * vx;
	y = ey + t * vy;
	z = ez + t * vz;
	double rx, ry, rz;
	double nx, ny, nz, dot;
	nx = x - cx;
	ny = y - cy;
	nz = z - cz;
	n = sqrt( nx * nx + ny * ny + nz * nz);
	nx /= n;
	ny /= n;
	nz /= n;
	dot = nx * vx + ny * vy + nz * vz;
#if 1
	rx = vx - 2 * dot * nx;
	ry = vy - 2 * dot * ny;
	rz = vz - 2 * dot * nz;
#else
	rx = 2 * dot * nx - vx;
	ry = 2 * dot * ny - vy;
	rz = 2 * dot * nz - vz;
#endif
	n = sqrt( rx * rx + ry * ry + rz * rz);
	rx /= n;
	ry /= n;
	rz /= n;
	printf( "vx=%f vy=%f vz=%f ni=%d t=%f x=%f y=%f z=%f nx=%f ny=%f nz=%f dot=%f rx=%f ry=%f rz=%f\n", vx, vy, vz, ni, t, x, y, z, nx, ny, nz, dot, rx, ry, rz);

	return 0;
}
