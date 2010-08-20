#include <stdio.h>
#include <math.h>

#define W 10
#define H 10
#define SMALL	0.001
#define BIG		10.0

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

int intersec_sphere( double cx, double cy, double cz, double sr, 
double ex, double ey, double ez, double vx, double vy, double vz,
double *t)
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
		if (abs( t1) > abs( t2))
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

int main()
{
	printf( "hello intersec\n");
	
	double cx, cy, cz, sr;
	double ex, ey, ez, vx, vy, vz;
	
	cx = 0.05; cy = 0.1; cz = 0.5; sr = 0.2;
	ex = 0; ey = 0; ez = 1; vx = 0.05; vy = 0.1; vz = -1;
	int w, h;
	w = 8; h = 6;
	double winw, winh;
	winw = 1.0; winh = 1.0;
	
	printf( "sphere: c(%f;%f;%f) r(%f)\n", cx, cy, cz, sr);
	printf( "eye: e(%f;%f;%f) v(%f;%f;%f)\n", ex, ey, ez, vx, vy, vz);
	
	int i, j;
	for (j = 0; j < h; j++)
	{
	for (i = 0; i < w; i++)
	{
		int res;
		double _vx, _vy, _vz;
		_vx = vx - winw / 2 + winw * (double)i / (w - 1);
		_vy = vy + winh / 2 - winh * (double)j / (h - 1);
		_vz = vz;
		double t = 0;
		res = intersec_sphere( cx, cy, cz, sr, ex, ey, ez, _vx, _vy, _vz, &t);
//		printf( "%c", res ? 's' : '.');
		if (res)
		{
#define COLS 9
			int pix = (int)((BIG - (t - (double)SMALL)) * COLS / BIG);
			if (pix < 0)
				pix = -pix;
			pix %= COLS;
			printf( "%d", pix);
		}
		else
			printf( ".");
	}
	printf( "\n");
	}
	
	return 0;
}
