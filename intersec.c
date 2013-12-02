/*
 * intersec : small naive ray tracer
 * copyright 2012 Nicolas Sauzede (nsauzede@laposte.net)
 * This code is GPLv3
 *
 */
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#ifdef USETIMEOFDAY
#include <sys/time.h>
#else
#include <time.h>
#endif

#ifdef WIN32
#include <windows.h>
#define MAP_FAILED ((HANDLE)-1)
#else
#include <sys/mman.h>
#endif

#ifdef USE_SDL
#include <SDL.h>
#endif

#if 0
#define W 60
#define H 60
#else
#define W 640
#define H 480
#endif

#define SMALL	0.001
#define BIG		10.0

#define WINSCALE 1.0

#ifdef DEBUG
#define dprintf(...) do{printf(__VA_ARGS__);}while(0)
#else
#define dprintf(...) do{}while(0)
#endif

typedef struct {
	double cx, cy, cz, sr;
	double r, g, b;
	double rdatt, gdatt, bdatt;		// diffuse
	double rratt, gratt, bratt;		// reflection
	double rfatt, gfatt, bfatt;		// refraction
	double rindex;					// refraction index
} sphere_t;

sphere_t spheres[] = {
#if 0
#if 1
	{ .cx = 0.05, .cy = 0.3, .cz = -0.8, .sr = 0.2, .r = 1.0, .g = 0.0, .b = 0.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#define D1 1.0
	{ .cx = -0.2, .cy = -0.1, .cz = -0.6, .sr = 0.1, .r = 0.0, .g = 1.0, .b = 0.0, .rdatt = D1, .gdatt = D1, .bdatt = D1, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#define D 0.8
#define R 0.0
	{ .cx = 0.5, .cy = -0.05, .cz = -0.7, .sr = 0.2, .r = 0.1, .g = 0.1, .b = 1.0, .rdatt = D, .gdatt = D, .bdatt = D, .rratt = R, .gratt = R, .bratt = R },
	{ .cx = 0.0, .cy = 0.0, .cz = 0.0, .sr = 0.03, .r = 1.0, .g = 1.0, .b = 1.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#endif
#define RFL 0.0
#define RFR 1.0
#define RINDEX 1.0
	{ .cx = 0.0, .cy = 0.0, .cz = 0.3, .sr = 0.1, .r = 1.0, .g = 1.0, .b = 1.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = RFL, .gratt = RFL, .bratt = RFL, .rfatt = RFR, .gfatt = RFR, .bfatt = RFR, .rindex = RINDEX },
#elif 1
#define SR 0.05
#define R 1.0
#define G 1.0
#define B 1.0
#define Z 0
#define CX 0.05*WINSCALE
#define CY 0.08*WINSCALE
#define CZ 0.1*WINSCALE
//	{ .cx = 0*CX, .cy = 0*CY, .cz = -50*CZ, .sr = 5.0*WINSCALE, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#if 1 /* ioccc ray */
	{ .cx = -7*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -3*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 5*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 8*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -7*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -4*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 1*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 4*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -7*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -4*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 1*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 4*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -7*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -3*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 5*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 8*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

#undef CX
#undef CY
#undef CZ
#undef R
#undef G
#undef B
#define R 0.0
#define G 0.0
#define B 0.0
#define CX 0.025*WINSCALE
#define CY 0.013*WINSCALE
#define CZ 0.05*WINSCALE
#undef Z
#define Z 5
	{ .cx = -11*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -9*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 11*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -11*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -7*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 11*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -11*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -9*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 9*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -11*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -7*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 9*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#endif
#else
	{ .cx = 0.0, .cy = -0.1, .cz = 0.0, .sr = 0.2, .r = 1.0, .g = 0.0, .b = 0.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0.0, .cy = 0.1, .cz = 0.0, .sr = 0.2, .r = 0.0, .g = 1.0, .b = 0.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0.0, .cy = 0.4, .cz = 0.0, .sr = 0.2, .r = 0.0, .g = 0.0, .b = 1.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#endif
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

int intersec_sphere( double cx, double cy, double cz, double sr, double ex, double ey, double ez, double vx, double vy, double vz, double *tmin, double *tmax)
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
		{
			double temp = t1;
			t1 = t2;
			t2 = temp;
		}
	}
	else if (sol == 1)
	{
//		printf( "one solution : %f\n", t1);
		t2 = t1;
	}
	else
	{
//		printf( "no solutions\n");
	}
	if (sol && (t1 > SMALL))
	{
		result = 1;
//		printf( "returning %f\n", t1);
		if (tmin)
			*tmin = t1;
		if (tmax)
			*tmax = t2;
	}

	return result;
}

#define ATT_MIN 0.001
#define LEV_MAX -1

int sky_color( int level, double ex, double ey, double ez, double _vx, double _vy, double _vz, double *_r, double *_g, double *_b)
{
	double rmin = 0, gmin = 0, bmin = 0;
#ifdef USE_SKY
	double coef2;
	double th;
//	double ph;
	double n = 1;
	n = sqrt( _vx * _vx + _vy * _vy + _vz * _vz);
	th = asin( _vy / n) * 180 / M_PI;
//	th = acos( _vx / n) * 180 / M_PI;
//	th = atan( _vy / _vx) * 180 / M_PI;
#define THR 180
	double thr = THR;
#if 1
	if (level == 0)
	{
		thr /= 3;
		coef2 = (fmod( fabs( th - thr / 2) , thr)) / thr;
	}
	else
#endif
	coef2 = (fmod( fabs( th - thr / 2) , thr)) / thr;
	double a, b;
	double x1, x2, y1, y2;
	if (coef2 <= 0.333)
	{
		x1 = 0; x2 = 0.333;
		coef2 -= x1;

		y1 = 6; y2 = 255;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		rmin = a * coef2 + b;

		y1 = 105; y2 = 255;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		gmin = a * coef2 + b;

		y1 = 155; y2 = 255;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		bmin = a * coef2 + b;
	}
	else if (coef2 <= 0.666)
	{
		x1 = 0.333; x2 = 0.666;
		coef2 -= x1;

		y1 = 255; y2 = 218;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		rmin = a * coef2 + b;

		y1 = 255; y2 = 178;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		gmin = a * coef2 + b;

		y1 = 255; y2 = 127;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		bmin = a * coef2 + b;
	}
	else if (coef2 <= 1.0)
	{
		x1 = 0.666; x2 = 1.0;
		coef2 -= x1;

		y1 = 218; y2 = 103;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		rmin = a * coef2 + b;

		y1 = 178; y2 = 55;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		gmin = a * coef2 + b;

		y1 = 127; y2 = 26;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		bmin = a * coef2 + b;
	}
	else
	{
		printf( "boom sky color coef2 out of bound coef2=%f\n", coef2);fflush( stdout);
		getchar();
		exit( 3);
	}
	rmin /= 255;
	gmin /= 255;
	bmin /= 255;
				if ((rmin > 1.0) || (gmin > 1.0) || (bmin > 1.0) || (rmin < 0.0) || (gmin < 0.0) || (bmin < 0.0))
				{
					printf( "boom sky color overflow r=%f g=%f b=%f\n", rmin, gmin, bmin);fflush( stdout);
					getchar();
					exit( 3);
				}
#endif
	*_b = bmin;
	*_g = gmin;
	*_r = rmin;

	return 0;
}

unsigned long traced = 0;
unsigned long intersected = 0;
unsigned long reflected = 0;
unsigned long refracted = 0;
int level_max_reached = -1;
int level_max = LEV_MAX;
int traceray( int level, double ex, double ey, double ez, double _vx, double _vy, double _vz, double *r, double *g, double *b, char *pix, double ratt, double gatt, double batt)
{
	if (level == 0)
		traced++;
	if (level > level_max_reached)
		level_max_reached = level;
	unsigned long smin = -1, s;
	double tmin = BIG, tmax = 0, srmin = 0;
	double rmin = 0.0, gmin = 0.0, bmin = 0.0;		// material
	double rdatt = 0.0, gdatt = 0.0, bdatt = 0.0;	// diffuse
	double rratt = 0.0, gratt = 0.0, bratt = 0.0;	// reflected
	double rfatt = 0.0, gfatt = 0.0, bfatt = 0.0;	// refracted
	double rindex = 0.0;
	
	for (s = 0; s < nsph; s++)
	{
		double cx, cy, cz, sr;
		cx = spheres[s].cx; cy = spheres[s].cy; cz = spheres[s].cz; sr = spheres[s].sr;
//		printf( "sphere: c(%f;%f;%f) r(%f)\n", cx, cy, cz, sr);
		int res;
		double t = 0, t2 = 0;
		res = intersec_sphere( cx, cy, cz, sr, ex, ey, ez, _vx, _vy, _vz, &t, &t2);
//		printf( "%c", res ? 's' : '.');
		if (res)
		{
			if (t < tmin)
			{
				smin = s;
				tmin = t;
				tmax = t2;
				srmin = spheres[s].sr;
				rmin = spheres[s].r;
				gmin = spheres[s].g;
				bmin = spheres[s].b;
				rdatt = spheres[s].rdatt;
				gdatt = spheres[s].gdatt;
				bdatt = spheres[s].bdatt;
				rratt = spheres[s].rratt;
				gratt = spheres[s].gratt;
				bratt = spheres[s].bratt;
				rfatt = spheres[s].rfatt;
				gfatt = spheres[s].gfatt;
				bfatt = spheres[s].bfatt;
				rindex = spheres[s].rindex;
//				pix = '0' + s;
				if (spheres[s].r > 0)
					*pix = 'r';
				else if (spheres[s].g > 0)
					*pix = 'g';
				else if (spheres[s].b > 0)
					*pix = 'b';
				else
					*pix = 'k';
			}
		}
	}
	if (tmin >= BIG)	// sky
	{
		sky_color( level, ex, ey, ez, _vx, _vy, _vz, &rmin, &gmin, &bmin);
#if 1
		if (level > 0)
			dprintf( "sky color %f %f %f\n", rmin, gmin, bmin);
#endif
	}
	else				// object
	{
		double rdiff = 0.0, gdiff = 0.0, bdiff = 0.0;
		double rrefl = 0.0, grefl = 0.0, brefl = 0.0;
		double rrefr = 0.0, grefr = 0.0, brefr = 0.0;
	
		if (level > 0)
			intersected++;
//		printf( "object\n");
// compute intersections here
		double rx, ry, rz, rx2, ry2, rz2;		// coord of intersec
		rx = ex + _vx * tmin;
		ry = ey + _vy * tmin;
		rz = ez + _vz * tmin;
		rx2 = ex + _vx * tmax;
		ry2 = ey + _vy * tmax;
		rz2 = ez + _vz * tmax;
		double n;
		dprintf( "ray is (%f,%f,%f)\n", _vx, _vy, _vz);
		dprintf( "inters with sphere %lu at (%f,%f,%f)\n", smin, rx, ry, rz);
		double nvx, nvy, nvz;	// normal vect

// compute distance between two intersec points
		nvx = rx - rx2;
		nvy = ry - ry2;
		nvz = rz - rz2;
		n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
		double grad;
		grad = n / srmin / 2;
		grad = sqrt( grad);
		rdiff = rmin * rdatt * grad;
		gdiff = gmin * gdatt * grad;
		bdiff = bmin * bdatt * grad;

		double rvx, rvy, rvz;	// vect of reflected/refracted ray

		if (1
			&& (((ratt * rratt) >= ATT_MIN) || ((gatt * gratt) >= ATT_MIN) || ((batt * bratt) >= ATT_MIN))
			&& ((level < level_max) || (level_max == -1))
		)
		{
// loose energy
			double loss = 0.9;
			ratt *= loss;
			gatt *= loss;
			batt *= loss;
		
			if (level > 0)
				reflected++;
// normal at intersec
		nvx = rx - spheres[smin].cx;
		nvy = ry - spheres[smin].cy;
		nvz = rz - spheres[smin].cz;
		n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
		nvx /= n;
		nvy /= n;
		nvz /= n;

		double dot = 1.0;
		dot = (_vx * nvx + _vy * nvy + _vz * nvz);
		rvx = _vx - 2 * dot * nvx;
		rvy = _vy - 2 * dot * nvy;
		rvz = _vz - 2 * dot * nvz;
		n = sqrt( rvx * rvx + rvy * rvy + rvz * rvz);
		rvx /= n;
		rvy /= n;
		rvz /= n;
		dprintf( "refl vec (%f,%f,%f) (s=%lu)\n", rvx, rvy, rvz, smin);
			traceray( level + 1, rx, ry, rz, rvx, rvy, rvz, &rrefl, &grefl, &brefl, pix, ratt * rratt, gatt * gratt, batt * bratt);
		}

		if (1
			&& (((ratt * rfatt) >= ATT_MIN) || ((gatt * gfatt) >= ATT_MIN) || ((batt * bfatt) >= ATT_MIN))
			&& ((level < level_max) || (level_max == -1))
		)
		{
// loose energy
			double loss = 0.9;
			ratt *= loss;
			gatt *= loss;
			batt *= loss;
		
			if (level > 0)
				refracted++;
#if 1
// normal at intersec
		nvx = rx - spheres[smin].cx;
		nvy = ry - spheres[smin].cy;
		nvz = rz - spheres[smin].cz;
		n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
		nvx /= n;
		nvy /= n;
		nvz /= n;

		double dot = 1.0;
		dot = (-_vx * nvx + -_vy * nvy + -_vz * nvz);
		rvx = -dot * nvx - _vx;
		rvy = -dot * nvy - _vy;
		rvz = -dot * nvz - _vz;
		rvx *= 1.0 - rindex;
		rvy *= 1.0 - rindex;
		rvz *= 1.0 - rindex;
		rvx += _vx;
		rvy += _vy;
		rvz += _vz;
		n = sqrt( rvx * rvx + rvy * rvy + rvz * rvz);
		rvx /= n;
		rvy /= n;
		rvz /= n;
#else
		rvx = _vx;
		rvy = _vy;
		rvz = _vz;
#endif
		dprintf( "refr vec (%f,%f,%f) (s=%lu)\n", rvx, rvy, rvz, smin);
			traceray( level + 1, rx, ry, rz, rvx, rvy, rvz, &rrefr, &grefr, &brefr, pix, ratt * rfatt, gatt * gfatt, batt * bfatt);
		}

		if ((rdatt + rratt + rfatt + gdatt + gratt + gfatt + bdatt + bratt + bfatt) > SMALL)
		{
			rmin = (rdiff * rdatt + rrefl * rratt + rrefr * rfatt) / (rdatt + rratt + rfatt);
			gmin = (gdiff * gdatt + grefl * gratt + grefr * gfatt) / (gdatt + gratt + gfatt);
			bmin = (bdiff * bdatt + brefl * bratt + brefr * bfatt) / (bdatt + bratt + bfatt);
		}

		if ((rmin > 1.0) || (gmin > 1.0) || (bmin > 1.0) || (rmin < 0.0) || (gmin < 0.0) || (bmin < 0.0))
		{
			printf( "boom object color overflow r=%f g=%f b=%f\n", rmin, gmin, bmin);fflush( stdout);
			getchar();
			exit( 3);
		}
	}
	*r = rmin;
	*g = gmin;
	*b = bmin;
	return 0;
}

int stats()
{
	double cur_t;
#ifndef USETIMEOFDAY
	cur_t = time( 0);
#else
	struct timeval tv;
	gettimeofday( &tv, 0);
	cur_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
	static double old_t = SMALL / 2;
	if (old_t <= SMALL)
		old_t = cur_t;
	printf( " trcd=%lu intr=%lu ", traced, intersected);
	printf( "refl=%lu refr=%lu MAX=%d max=%d ", reflected, refracted, level_max, level_max_reached);
	double duration = cur_t - old_t;
	unsigned long last_traced = 0;
	unsigned long perf = (traced - last_traced) / duration;
	if (duration > 0)
	{
		printf( "- %lu ray/s dur=%.2fs", perf, duration);
	}
	printf( "\n");
//	old_t = cur_t;
	last_traced = traced;
	return 0;
}

#define USE_CAMERA
#ifdef USE_CAMERA
typedef double v3[3];

void cross3( v3 dest, v3 src1, v3 src2)
{
	dest[0] = src1[1] * src2[2] - src1[2] * src2[1];
	dest[1] = src1[2] * src2[0] - src1[0] * src2[2];
	dest[2] = src1[0] * src2[1] - src1[1] * src2[0];
}

void sum3( v3 dest, v3 src1, v3 src2)
{
	dest[0] = src1[0] + src2[0];
	dest[1] = src1[1] + src2[1];
	dest[2] = src1[2] + src2[2];
}

void diff3( v3 dest, v3 src1, v3 src2)
{
	dest[0] = src1[0] - src2[0];
	dest[1] = src1[1] - src2[1];
	dest[2] = src1[2] - src2[2];
}
#endif

int main( int argc, char *argv[])
{
#ifdef USE_SDL
	SDL_Surface *screen = 0;
	int do_sdl = 1;
#endif
	double ex, ey, ez, vx, vy, vz;
	
	ex = 0*WINSCALE; ey = 0*WINSCALE; ez = 1*WINSCALE; vx = 0.0; vy = 0.0; vz = -1*WINSCALE;
	unsigned long w, h;
	w = W; h = H;
	double winw, winh;
	double winscale = WINSCALE;
	int do_tga = 1;
	int do_txt = 0;
	int bpp = 24;
	int tga_fd = -1;
	unsigned char *tga_map = MAP_FAILED;
	size_t tga_size = 0;
	size_t tga_index = 0;

	int arg = 1;
	if (argc > arg)
	{
		sscanf( argv[arg++], "%lu", &w);
		if (argc > arg)
		{
			sscanf( argv[arg++], "%lu", &h);
			if (argc > arg)
			{
				sscanf( argv[arg++], "%d", &level_max);
			}
		}
	}
#ifdef USE_SDL
	if (do_sdl)
	{
		SDL_Init( SDL_INIT_VIDEO);
		screen = SDL_SetVideoMode( w, h, bpp, 0);
		if (screen)
			atexit( SDL_Quit);
	}
#endif
	winw = 1.0 * winscale; winh = winw * h / w;
	if (do_tga)
	{
		tga_fd = open( "out.tga", O_CREAT | O_RDWR 
#ifdef WIN32
			| O_BINARY
#endif
			, S_IRUSR | S_IWUSR);
		if (tga_fd == -1)
		{
			perror( "open");
			exit( 2);
		}
		tga_size = w * h * bpp / 8 + 18;
		dprintf( "tga_size=%lu\n", (unsigned long)tga_size);
		if (ftruncate( tga_fd, tga_size))
			perror( "ftruncate");
#ifdef WIN32
		tga_map = malloc( tga_size);
		if (!tga_map)
			tga_map = MAP_FAILED;
		lseek( tga_fd, 0, SEEK_SET);
		read( tga_fd, tga_map, tga_size);
#else
		tga_map = mmap( 0, tga_size, PROT_READ | PROT_WRITE, MAP_SHARED, tga_fd, 0);
#endif
#define TGA_BYTE(b) do {tga_map[tga_index++] = b;}while(0)
#define TGA_SHORT(s) do{*(unsigned short *)(&tga_map[0] + tga_index) = s;tga_index+=sizeof( unsigned short);}while(0)
		TGA_BYTE( 0);
		TGA_BYTE( 0);
		TGA_BYTE( 2);
		TGA_BYTE( 0);
		TGA_BYTE( 0);
		TGA_BYTE( 0);
		TGA_BYTE( 0);
		TGA_BYTE( 0);
		TGA_SHORT( 0);
		TGA_SHORT( 0);
		TGA_SHORT( w);
		TGA_SHORT( h);
		TGA_BYTE( bpp);
		TGA_BYTE( 0x20);
	}
	if (do_txt)
	{
		printf( "w=%lu h=%lu\n", w, h);
		printf( "eye: e(%f;%f;%f) v(%f;%f;%f)\n", ex, ey, ez, vx, vy, vz);
	}
	
	printf( "w=%lu h=%lu winw=%.2f winh=%.2f vx=%f vy=%f vz=%f level_max=%d\n", w, h, winw, winh, vx, vy, vz, level_max);
	unsigned long i, j;
	double old_t;
	double last_t;
	double cur_t;
#ifndef USETIMEOFDAY
	old_t = time( 0);
#else
	struct timeval tv;
	gettimeofday( &tv, 0);
	old_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
	last_t = old_t;
	for (j = 0; j < h; j++)
	{
		double _vx, _vy, _vz;
//		_vy = vy + (((double)h - 1) - j) / (h - 1) - winh / 2;
		_vy = vy - winh / 2 + winh * (double)(h - j - 1) / (h - 1);
		_vz = vz;

		for (i = 0; ; i++)
		{
#ifndef USETIMEOFDAY
			cur_t = time( 0);
#else
			gettimeofday( &tv, 0);
			cur_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
			int percent = (double)100.0 * (j * w + i) / (w * h);
			static int old_percent = -1;
			if ((percent != old_percent) || (cur_t >= (last_t + 1)))
			{
				if (cur_t >= (last_t + 1))
					last_t = cur_t;
				old_percent = percent;
				if (!(old_percent % 10))
				{
					printf( " %d%%", old_percent);
					stats();
					fflush( stdout);
				}
#ifdef USE_SDL
				if (screen)
				{
					SDL_UpdateRect( screen, 0, 0, 0, 0);
				}
#endif
			}
			if (i >= w)
				break;
			char pix = '.';
#ifndef USE_CAMERA
			_vx = vx - winw / 2 + winw * (double)i / (w - 1);
#else
			v3 v;
			v[0] = vx;
			v[1] = vy;
			v[2] = vz;
			v3 e;
			e[0] = ex;
			e[1] = ey;
			e[2] = ez;
			// p=p0+(p1-p0)u+(p2-p0)v
			v3 up = { 0, 1, 0};
			v3 right;
			cross3( right, v, up);
			cross3( up, right, v);
			v3 s0, s1, s2, s;
			sum3( s0, e, v);
			sum3( s1, s0, up);
			sum3( s2, s0, right);
			double u, vv;
			vv = -(double)0.5 + (double)1.0 * (double)i / ((double)w - 1);
			u = -(double)0.5 + (double)1.0 * (double)(h - j - 1) / ((double)h - 1);
			s[0] = s0[0] + (s1[0] - s0[0]) * u + (s2[0] - s0[0]) * vv;
			s[1] = s0[1] + (s1[1] - s0[1]) * u + (s2[1] - s0[1]) * vv;
			s[2] = s0[2] + (s1[2] - s0[2]) * u + (s2[2] - s0[2]) * vv;
			v3 _v;
			diff3( _v, s, e);
			_vx = _v[0];
			_vy = _v[1];
			_vz = _v[2];
#endif

//			printf( "j=%lu i=%lu vy=%f vx=%f\n", j, i, _vy, _vx);
			double r = 0.0, g = 0.0, b = 0.0;

			traceray( 0, ex, ey, ez, _vx, _vy, _vz, &r, &g, &b, &pix, 1.0, 1.0, 1.0);
			unsigned char cr, cg, cb;
			cr = (double)255 * r;
			cg = (double)255 * g;
			cb = (double)255 * b;
#ifdef USE_SDL
			if (screen)
			{
				Uint32 col;
				SDL_Rect rect;
				rect.x = i;
				rect.y = j;
				rect.w = rect.h = 1;
				col = SDL_MapRGB( screen->format, cr, cg, cb);
				SDL_FillRect( screen, &rect, col);
			}
#endif
			//dprintf( "  (r=%f g=%f b=%f)", r, g, b);
			if (do_tga)
			{
#if 1
				static unsigned long count = 0;
				if ((r > 1.0) || (g > 1.0) || (b > 1.0) || (r < 0.0) || (g < 0.0) || (b < 0.0))
				{
					printf( "boom color overflow %lu at (%lu,%lu) r=%f g=%f b=%f\n", count, i, j, r, g, b);fflush( stdout);
					getchar();
					exit( 3);
				}
				count++;
#endif
			//	dprintf( "{%02x:%02x:%02x}", cr, cg, cb);
				TGA_BYTE( cb);
				TGA_BYTE( cg);
				TGA_BYTE( cr);
//				TGA_BYTE( 0);
			}
			if (do_txt)
				printf( "%c", pix);
		}
		if (do_txt)
			printf( "\n");
		//dprintf( "\n");
	}
#ifndef USETIMEOFDAY
	cur_t = time( 0);
#else
	gettimeofday( &tv, 0);
	cur_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
	printf( "\n");
	if (do_tga)
	{
#ifdef WIN32
		lseek( tga_fd, 0, SEEK_SET);
		write( tga_fd, tga_map, tga_size);
		free( tga_map);
#else
		munmap( tga_map, tga_size);
#endif
		close( tga_fd);
	}
//	dprintf( "tga_size=%lu\n", (unsigned long)tga_size);
//	dprintf( "tga_index=%lu\n", (unsigned long)tga_index);
	printf( "traced=%lu intersected=%lu\n", traced, intersected);
	printf( "reflected=%lu refracted=%lu level_max=%d level_max_reached=%d\n", reflected, refracted, level_max, level_max_reached);
	double duration = cur_t - old_t;
	if (!duration)
		duration = 1;
	unsigned long perf = (w * h) / duration;
	printf( "perf : %lu ray/s duration=%.2fs (old_t=%.2f cur_t=%.2f)\n", perf, duration, old_t, cur_t);
#ifdef USE_SDL
	if (screen)
	{
		SDL_Event event;
		int end = 0;
		while (!end)
		{
			SDL_WaitEvent( &event);
			switch (event.type)
			{
				case SDL_QUIT:
					end = 1;
					break;
				case SDL_KEYDOWN:
					end = 1;
					break;
				case SDL_MOUSEBUTTONDOWN:
					end = 1;
					break;
				default:
					break;
			}
		}
	}
#endif
	
	return 0;
}
