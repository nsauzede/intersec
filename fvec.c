#include <math.h>

#include "fvec.h"

scalar_t dot3( v3 p, v3 n)
{
	scalar_t res;
	scalar_t t1 = p[0] * n[0] / COEF;
	scalar_t t2 = p[1] * n[1] / COEF;
	scalar_t t3 = p[2] * n[2] / COEF;
	res = (t1 + t2 + t3);
	return res;
}

scalar_t norm3( v3 n)
{
	scalar_t res = dot3( n, n);
	double temp = D(res);
	temp = sqrt( temp);
	res = S( temp);
	return res;
}

scalar_t* sum3( v3 dest, v3 const src1, v3 const src2)
{
	dest[0] = src1[0] + src2[0];
	dest[1] = src1[1] + src2[1];
	dest[2] = src1[2] + src2[2];
	return dest;
}

scalar_t* div3( v3 n, scalar_t nn)
{
	n[0] *= COEF;
	n[1] *= COEF;
	n[2] *= COEF;

	n[0] /= nn;
	n[1] /= nn;
	n[2] /= nn;

	return n;
}

#define abs(s) (s < 0 ? -s : s)
int intersec_box( v3 lower, v3 upper, v3 e, v3 v, scalar_t *_tmin, scalar_t *_tmax)
{
	scalar_t tmin, tmax, tminx, tmaxx, tminy, tmaxy, tminz, tmaxz, t1, t2;
	if (abs(v[0]) < SMALL)
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
		t1 = (lower[0] - e[0]) * COEF / v[0];
		t2 = (upper[0] - e[0]) * COEF / v[0];
		tminx = MIN( t1, t2);
		tmaxx = MAX( t1, t2);
		if (tmaxx < 0)
			return 0;
	}
	if (abs(v[1]) < SMALL)
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
		t1 = (lower[1] - e[1]) * COEF / v[1];
		t2 = (upper[1] - e[1]) * COEF / v[1];
		tminy = MIN( t1, t2);
		tmaxy = MAX( t1, t2);
		if (tmaxy < 0)
			return 0;
	}
	if (abs(v[2]) < SMALL)
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
		t1 = (lower[2] - e[2]) * COEF / v[2];
		t2 = (upper[2] - e[2]) * COEF / v[2];
		tminz = MIN( t1, t2);
		tmaxz = MAX( t1, t2);
		if (tmaxz < 0)
			return 0;
	}
	tmin = MAX3( tminx, tminy, tminz);
	tmax = MIN3( tmaxx, tmaxy, tmaxz);
	if (tmax < tmin)
		return 0;

	if (_tmin)
		*_tmin = tmin;
	if (_tmax)
		*_tmax = tmax;
	return 2;
}

