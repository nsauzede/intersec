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

void set3( v3 dest, double x, double y, double z)
{
	dest[0] = x;
	dest[1] = y;
	dest[2] = z;
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

int intersec_parallelogram( v3 p0, v3 p1, v3 p2, v3 l0, v3 l, double *pt)
{
	v3 delta, location;
	double dot,pos1,pos2,gu[4],gv[4];
	int i,j,crossing_no; 
	
	// plane normal :
	v3 vect1, vect2;
	diff3( vect1, p1, p0);
	diff3( vect2, p2, p0);
	v3 norm;
	cross3( norm, vect1, vect2);
	double nn;
	nn = norm3( norm);
	div3( norm, nn);				// should be done at constr

	dot = dot3( norm, l);
	
	double n1 = dot3( norm, p0);	// should be done at obj constr ?

	if (fabs(dot)<SMALL)
		return 0;
	pos1 = n1;
	pos2 = dot3( norm, l0);
	*pt=(pos1-pos2)/dot;
	sum3( location, l0, mult3( copy3( location, l), *pt));
	diff3( delta, location, p0);
	if ((fabs(norm[0]) > fabs(norm[1])) && (fabs(norm[0]) > fabs(norm[2])))
	{ 
		gu[0] = - delta[1];
		gv[0] = - delta[2];
		gu[1] = vect1[1] - delta[1];
		gv[1] = vect1[2] - delta[2];
		gu[2] = vect2[1] + vect1[1] - delta[1];
		gv[2] = vect2[2] + vect1[2] - delta[2];
		gu[3] = vect2[1] - delta[1];
		gv[3] = vect2[2] - delta[2];
	}
	else
	{
		if (fabs(norm[1]) >= fabs(norm[2]))
		{ 
			gu[0] = - delta[0];
			gv[0] = - delta[2];
			gu[1] = vect1[0] - delta[0];
			gv[1] = vect1[2] - delta[2];
			gu[2] = vect2[0] + vect1[0] - delta[0];
			gv[2] = vect2[2] + vect1[2] - delta[2]; 
			gu[3] = vect2[0] - delta[0];
			gv[3] = vect2[2] - delta[2]; 
		}
		else
		{
			gu[0] = - delta[0];
			gv[0] = - delta[1];
			gu[1] = vect1[0] - delta[0];
			gv[1] = vect1[1] - delta[1];
			gu[2] = vect2[0] + vect1[0] - delta[0];
			gv[2] = vect2[1] + vect1[1] - delta[1];
			gu[3] = vect2[0] - delta[0];
			gv[3] = vect2[1] - delta[1]; 
		}
	}
	crossing_no = 0;
	for (i=0; i<4; i++)
	{ 
		j = (i + 1) % 4;
		if (((gv[i] < 0) && (gv[j] >= 0)) || ((gv[j] < 0) && (gv[i]>= 0))) 
		{
			if ((gu[i]>=0) && (gu[j] >= 0))
				crossing_no++; 
			else
			{
				if ((gu[i]>=0) || (gu[j] >= 0))
				{
					if (gu[i] - gv[i] * (gu[j] + gu[i]) / (gv[j] - gv[i]) >0)
						crossing_no++;
				}
			}
		}
	}
	if ((crossing_no % 2) == 0)
		return 0;
	return 1;
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

#if 0
// taken from "fractal prog & ray tracing with c++" book
int intersec_sphere( v3 cs, double sr, v3 e, v3 v, double *tmin, double *tmax)
{
	// already computed ?
	double n1 = sr * sr;

	double a, b, c, d;
	v3 temp;

	double t1, t2;

	diff3( temp, cs, e);
	c = dot3( temp, temp) - n1;
	b = -2 * dot3( v, temp);
	a = dot3( v, v);
	d = b * b - 4 * a * c;
	if (d <= 0)
		return 0;
	d = sqrt( d);
	t1 = (-b + d) / (a + a);
	t2 = (-b - d) / (a + a);
	if (t1 <= SMALL && t2 <= SMALL)
		return 0;
	if (t1 < t2)
	{
		*tmin = t1;
	}
	else
	{
		*tmin = t2;
	}

	return 1;
}
#else
int intersec_sphere( v3 cs, double sr, v3 e, v3 v, double *tmin, double *tmax)
{
	int result = 0;
	double a, b, c;
	double t1, t2;
	
	v3 t;
	double n;
	n = sr * sr;
	diff3( t, e, cs);
	a = dot3( v, v);
	b = 2 * dot3( v, t);
	c = dot3( t, t) - n;
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
#endif

#define BIG 3e30
#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
#define MIN3(a,b,c) (a < b ? (a < c ? a : c) : b < c ? b : c)
#define MAX3(a,b,c) (a > b ? (a > c ? a : c) : b > c ? b : c)
int intersec_box( v3 lower, v3 upper, v3 e, v3 v, double *_tmin, double *_tmax, int *num)
{
	double tmin, tmax, tminx, tmaxx, tminy, tmaxy, tminz, tmaxz, t1, t2;

	int numx = 0, numy = 0, numz = 0;

	if (fabs(v[0]) < SMALL)
	{
		if ((lower[0] < e[0]) && (upper[0] > e[0]))
		{
			tminx = -BIG;
			tmaxx = BIG;
		}
		else
			return 0;
	}
	else
	{
		t1 = (lower[0] - e[0]) / v[0];
		t2 = (upper[0] - e[0]) / v[0];
		tminx = MIN( t1, t2);
		tmaxx = MAX( t1, t2);
		if (tmaxx < 0)
			return 0;

		if (tminx == t1)
			numx = 1;
		else
			numx = 2;
	}
	if (fabs(v[1]) < SMALL)
	{
		if ((lower[1] < e[1]) && (upper[1] > e[1]))
		{
			tminy = -BIG;
			tmaxy = BIG;
		}
		else
			return 0;
	}
	else
	{
		t1 = (lower[1] - e[1]) / v[1];
		t2 = (upper[1] - e[1]) / v[1];
		tminy = MIN( t1, t2);
		tmaxy = MAX( t1, t2);
		if (tmaxy < 0)
			return 0;

		if (tminy == t1)
			numy = 1;
		else
			numy = 2;
	}
	if (fabs(v[2]) < SMALL)
	{
		if ((lower[2] < e[2]) && (upper[2] > e[2]))
		{
			tminz = -BIG;
			tmaxz = BIG;
		}
		else
			return 0;
	}
	else
	{
		t1 = (lower[2] - e[2]) / v[2];
		t2 = (upper[2] - e[2]) / v[2];
		tminz = MIN( t1, t2);
		tmaxz = MAX( t1, t2);
		if (tmaxz < 0)
			return 0;

		if (tminz == t1)
			numz = 1;
		else
			numz = 2;
	}
	tmin = MAX3( tminx, tminy, tminz);
	tmax = MIN3( tmaxx, tmaxy, tmaxz);
	if (tmax < tmin)
		return 0;

	if (num)
	{
	if (tmin == tminx)
		*num = numx;
	else if (tmin == tminy)
		*num = numy << 2;
	else if (tmin == tminz)
		*num = numz << 4;
	}

	if (_tmin)
		*_tmin = tmin;
	if (_tmax)
		*_tmax = tmax;
	return 2;
}

int use_solid = 0;
int intersec_quad( v3_t obj, v3_t e, v3_t v, v1_t *tmin, v1_t *tmax)
{
	int result = 0;

	v1_t a, b, c;
	v1_t t1 = 0, t2 = 0;
	int sol;

	v1_t z0 = 0;
	v1_t x0 = 0;
	v1_t y0 = 0;
	x0 = obj[0];
	y0 = obj[1];
	z0 = obj[2];
	v1_t rx = obj[3 * 1 + 0] * (double)M_PI / 180.0;
	v1_t ry = obj[3 * 1 + 1] * (double)M_PI / 180.0;
	v1_t rz = obj[3 * 1 + 2] * (double)M_PI / 180.0;

	v1_t A = obj[3 * 2 + 0];
	v1_t B = obj[3 * 2 + 1];
	v1_t C = obj[3 * 2 + 2];
	v1_t D = obj[3 * 2 + 3];
	v1_t E = obj[3 * 2 + 4];
	v1_t lim = obj[3 * 2 + 5];
	v1_t Vx = v[0];
	v1_t Vy = v[1];
	v1_t Vz = v[2];
	v1_t lx0 = e[0];
	v1_t ly0 = e[1];
	v1_t lz0 = e[2];
//translate
	lx0 -= x0;
	ly0 -= y0;
	lz0 -= z0;

	v1_t xx = 0, yy = 0, zz = 0;
	rx = rx;
	rz = rz;
	ry = ry;
	xx = xx;
	zz = zz;
	yy = yy;

	//rotate
	yy = ly0 * cos(rx) - lz0 * sin(rx);
	zz = ly0 * sin(rx) + lz0 * cos(rx);
	ly0 = yy;
	lz0 = zz;
	yy = Vy * cos(rx) - Vz * sin(rx);
	zz = Vy * sin(rx) + Vz * cos(rx);
	Vy = yy;
	Vz = zz;

	xx = lx0 * cos(ry) + lz0 * sin(ry);
	zz = -lx0 * sin(ry) + lz0 * cos(ry);
	lx0 = xx;
	lz0 = zz;
	xx = Vx * cos(ry) + Vz * sin(ry);
	zz = -Vx * sin(ry) + Vz * cos(ry);
	Vx = xx;
	Vz = zz;

	xx = lx0 * cos(rz) - ly0 * sin(rz);
	yy = lx0 * sin(rz) + ly0 * cos(rz);
	lx0 = xx;
	lz0 = yy;
	xx = Vx * cos(rz) - Vy * sin(rz);
	yy = Vx * sin(rz) + Vy * cos(rz);
	Vx = xx;
	Vy = yy;

	a = A * Vx * Vx + B * Vy * Vy + C * Vz * Vz;
	b = 2 * A * lx0 * Vx + 2 * B * ly0 * Vy + 2 * C * lz0 * Vz + E * Vy;
	c = A * lx0 * lx0 + B * ly0 * ly0 + C * lz0 * lz0 + E * ly0 - D;

	sol = solvetri( a, b, c, &t1, &t2);
	if (sol == 2)
	{
		if (t1 > t2)
		{
			v1_t temp = t1;
			t1 = t2;
			t2 = temp;
		}
	}
	else if (sol == 1)
	{
		t2 = t1;
	}
		v1_t z1 = 0;
		v1_t z2 = 0;

	if (sol 
		//&& (fabs(t1) > SMALL)
	)
	{
		int hit = 0;
		z1 = lz0 + Vz * t1;
		z2 = lz0 + Vz * t2;

		if (fabs( lim) < SMALL)
		{
			hit = 1;
		}
		else
		{
			v1_t zmin = z0;
			v1_t zmax = z0 + lim;
			if (use_solid)
			{
//				v1_t z1a = z1 - zmax;
//				v1_t z2a = z2 - zmax;

				if ((z1 < zmin) && (z2 > zmin))
				{
					t1 = (zmin - lz0) / Vz;
					hit = 1;
					result = 1;
				}
#if 0
				else
					if ((z1a * z2a) < 0)
					{
						t1 = (zmax - lz0) / Vz;
						hit = 1;
						result = 1;
					}
#endif
			}
			if (!hit)
			{
				hit = (z1 >= zmin) && (z1 <= zmax);
				if (!hit && ((sol > 1) && (fabs(t2) > SMALL)
							))
				{
					hit = (z2 >= zmin) && (z2 <= zmax);
					if (hit)
					{
						v1_t temp = t2;
						t1 = t2;
						t2 = temp;
					}
				}
			}
		}
		if (hit)
		{
			result = sol;
			if (tmin)
				*tmin = t1;
			if (tmax)
				*tmax = t2;
		}
	}

	return result;
}



