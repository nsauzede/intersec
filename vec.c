#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "vec.h"

#define dprintf(...) do{}while(0)

double* cross3( v3 dest, const v3 src1, v3 const src2)
{
	dest[0] = src1[1] * src2[2] - src1[2] * src2[1];
	dest[1] = src1[2] * src2[0] - src1[0] * src2[2];
	dest[2] = src1[0] * src2[1] - src1[1] * src2[0];
	return dest;
}

double* sum3( v3 dest, v3 const src1, v3 const src2)
{
	dest[0] = src1[0] + src2[0];
	dest[1] = src1[1] + src2[1];
	dest[2] = src1[2] + src2[2];
	return dest;
}

double* diff3( v3 dest, v3 const src1, v3 const src2)
{
	dest[0] = src1[0] - src2[0];
	dest[1] = src1[1] - src2[1];
	dest[2] = src1[2] - src2[2];
	return dest;
}

double* div3( v3 n, double nn)
{
	n[0] /= nn;
	n[1] /= nn;
	n[2] /= nn;
	return n;
}

double* mult3( v3 n, double nn)
{
	n[0] *= nn;
	n[1] *= nn;
	n[2] *= nn;
	return n;
}

double* add3( v3 n, double nn)
{
	n[0] += nn;
	n[1] += nn;
	n[2] += nn;
	return n;
}

double dot3( v3 p, v3 n)
{
	return p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
}

double norm3( v3 n)
{
	return sqrt( dot3( n, n));	
}

double *copy3( v3 dest, const v3 src)
{
	memcpy( dest, src, sizeof( v3));
	return dest;
}

void disp3( char *s, const v3 v)
{
	printf( "%s: %f,%f,%f\n", s, v[0], v[1], v[2]);
}

// plane : (p - po) . n = 0
// line : p = dl + l0
int intersec_plane( v3 p0, v3 p1, v3 p2, v3 l0, v3 l, double *pt)
{
	int result = 0;
	dprintf( "l0 is %f,%f,%f l is %f,%f,%f\n", l0[0], l0[1], l0[2], l[0], l[1], l[2]);
	dprintf( "p0 is %f,%f,%f p1 is %f,%f,%f p2 is %f,%f,%f\n", p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
	// plane normal :
	v3 v1, v2;
	diff3( v1, p1, p0);
	diff3( v2, p2, p0);
	dprintf( "v1 is %f,%f,%f v2 is %f,%f,%f\n", v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
	v3 n;
	cross3( n, v1, v2);
	double nn;
	nn = norm3( n);
	div3( n, nn);
	dprintf( "normal is %f,%f,%f\n", n[0], n[1], n[2]);
	// t = ((p0 - l0) . n) / (l . n)
	double num;
	v3 p;
	diff3( p, p0, l0);
	num = dot3( p, n);
	dprintf( "num is %f\n", num);
	double den;
	den = dot3( l, n);
	dprintf( "den is %f\n", den);
	if (den == 0)
	{
		if (num == 0)
		{
			dprintf( "all intersections\n");
		}
		else
		{
			dprintf( "no intersections\n");
		}
	}
	else
	{
		double det;
		double a, b, c, d, e, f, g, h, i;
		v3 lb;
		sum3( lb, l0, l);
		a = l0[0] - lb[0];
		b = p1[0] - p0[0];
		c = p2[0] - p0[0];
		d = l0[1] - lb[1];
		e = p1[1] - p0[1];
		f = p2[1] - p0[1];
		g = l0[2] - lb[2];
		h = p1[2] - p0[2];
		i = p2[2] - p0[2];
		det = a * (e * i - f * h) - b * (i * d - f * g) + c * (d * h - e * g);
		dprintf( "det=%f\n", det);
		if (det > 0.001) // matrix is invertible
		{
			double A, B, C, D, E, F, G, H, I;
			A = (e * i - f * h);
			B = -(d * i - f * g);
			C = (d * h - e * g);
			D = -(b * i - c * h);
			E = (a * i - c * g);
			F = -(a * h - b * g);
			G = (b * f - c * e);
			H = -(a * f - c * d);
			I = (a * e - b * d);
			double t, u = 0, v = 0;
			// p=p0+u*v1+v*v2
			t = (A * (l0[0] - p0[0]) + D * (l0[1] - p0[1]) + G * (l0[2] - p0[2])) / det;
			u = (B * (l0[0] - p0[0]) + E * (l0[1] - p0[1]) + H * (l0[2] - p0[2])) /det;
			v = (C * (l0[0] - p0[0]) + F * (l0[1] - p0[1]) + I * (l0[2] - p0[2])) / det;
			dprintf( "t=%f u=%f v=%f\n", t, u, v);
			if ((u >= 0 && u <= 1) && (v >= 0 && v <= 1) && ((u + v) <= 1))
			{
				dprintf( "inters inside triangle, yay !\n");
				if (t > 1) // _t must not be too small
				{
					if (pt)
						*pt = t;
					result = 1;
				}
			}
		}
	}
	return result;
}

#define SMALL 0.001
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

int intersec_sphere( v3 cs, double sr, v3 e, v3 v, double *tmin, double *tmax)
{
	double ex = e[0];
	double ey = e[1];
	double ez = e[2];
	double cx = cs[0];
	double cy = cs[1];
	double cz = cs[2];
	double vx = v[0];
	double vy = v[1];
	double vz = v[2];
	
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
		if (fabs( t1) > fabs( t2))
		{
			double temp = t1;
			t1 = t2;
			t2 = temp;
		}
	}
	else if (sol == 1)
	{
		t2 = t1;
	}
	else
	{
//		printf( "no solutions\n");
	}
	if (sol && (t1 > SMALL))
	{
		result = 1;
		if (tmin)
			*tmin = t1;
		if (tmax)
			*tmax = t2;
	}

	return result;
}
