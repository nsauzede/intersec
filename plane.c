#include <stdio.h>
#include <math.h>

typedef double v3[3];

// plane : (p - po) . n = 0
// line : p = dl + l0
int intersec_plane( v3 p0, v3 p1, v3 p2, v3 l0, v3 l, double *t)
{
	int result = 0;
	printf( "l0 is %f,%f,%f l is %f,%f,%f\n", l0[0], l0[1], l0[2], l[0], l[1], l[2]);
	printf( "p0 is %f,%f,%f p1 is %f,%f,%f p2 is %f,%f,%f\n", p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
	// plane normal :
	v3 v1, v2;
	v1[0] = p1[0] - p0[0];
	v1[1] = p1[1] - p0[1];
	v1[2] = p1[2] - p0[2];
	v2[0] = p2[0] - p0[0];
	v2[1] = p2[1] - p0[1];
	v2[2] = p2[2] - p0[2];
	printf( "v1 is %f,%f,%f v2 is %f,%f,%f\n", v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
	v3 n;
	n[0] = v1[1] * v2[2] - v1[2] * v2[1];
	n[1] = v1[2] * v2[0] - v1[0] * v2[2];
	n[2] = v1[0] * v2[1] - v1[1] * v2[0];
	double nn;
	nn = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	n[0] /= nn;
	n[1] /= nn;
	n[2] /= nn;
	printf( "normal is %f,%f,%f\n", n[0], n[1], n[2]);
	// t = ((p0 - l0) . n) / (l . n)
	double num;
	v3 p;
	p[0] = p0[0] - l0[0];
	p[1] = p0[1] - l0[1];
	p[2] = p0[2] - l0[2];
	num = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
	printf( "num is %f\n", num);
	double den;
	den = l[0] * n[0] + l[1] * n[1] + l[2] * n[2];
	printf( "den is %f\n", den);
	if (den == 0)
	{
		if (num == 0)
		{
			printf( "all intersections\n");
		}
		else
		{
			printf( "no intersections\n");
		}
	}
	else
	{
		double _t = num / den;
		printf( "one intersection at %f\n", _t);
		double x = l0[0] + _t * l[0];
		double y = l0[1] + _t * l[1];
		double z = l0[2] + _t * l[2];
		printf( "inters is %f,%f,%f\n", x, y, z);
		double det;
		double a, b, c, d, e, f, g, h, i;
		v3 lb;
		lb[0] = l0[0] + l[0];
		lb[1] = l0[1] + l[1];
		lb[2] = l0[2] + l[2];
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
		printf( "det=%f\n", det);
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
			double tt, u = 0, v = 0;
			// p=p0+u*v1+v*v2
			tt = (A * (l0[0] - p0[0]) + D * (l0[1] - p0[1]) + G * (l0[2] - p0[2])) / det;
			u = (B * (l0[0] - p0[0]) + E * (l0[1] - p0[1]) + H * (l0[2] - p0[2])) /det;
			v = (C * (l0[0] - p0[0]) + F * (l0[1] - p0[1]) + I * (l0[2] - p0[2])) / det;
			printf( "tt=%f u=%f v=%f\n", tt, u, v);
			if ((u >= 0 && u <= 1) && (v >= 0 && v <= 1) && ((u + v) <= 1))
			{
				printf( "inters inside triangle, yay !\n");
				if (_t > 1) // _t must not be too small
				{
					if (t)
						*t = _t;
					result = 1;
				}
			}
		}
	}
	return result;
}

int main()
{
	printf( "hello plane\n");
	v3 p0 = { 1, 0, 0 };
	v3 p1 = { 0, 1, 0 };
	v3 p2 = { 0, 0, 1 };
	v3 e = { 4, 5, 5 };
	v3 v = { -1, -1, -1 };
	double t = 0;
	int res = intersec_plane( p0, p1, p2, e, v, &t);
	printf( "result=%d t=%f\n", res, t);
	if (res)
	{
		double x = e[0] + t * v[0];
		double y = e[1] + t * v[1];
		double z = e[2] + t * v[2];
		printf( "inters is %f,%f,%f\n", x, y, z);
	}
	return 0;
}
