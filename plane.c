#include <stdio.h>
#include <math.h>

#define dprintf(...) do{}while(0)

typedef double v3[3];

// plane : (p - po) . n = 0
// line : p = dl + l0
int intersec_plane( v3 p0, v3 p1, v3 p2, v3 l0, v3 l, double *pt)
{
	int result = 0;
	dprintf( "l0 is %f,%f,%f l is %f,%f,%f\n", l0[0], l0[1], l0[2], l[0], l[1], l[2]);
	dprintf( "p0 is %f,%f,%f p1 is %f,%f,%f p2 is %f,%f,%f\n", p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
	// plane normal :
	v3 v1, v2;
	v1[0] = p1[0] - p0[0];
	v1[1] = p1[1] - p0[1];
	v1[2] = p1[2] - p0[2];
	v2[0] = p2[0] - p0[0];
	v2[1] = p2[1] - p0[1];
	v2[2] = p2[2] - p0[2];
	dprintf( "v1 is %f,%f,%f v2 is %f,%f,%f\n", v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
	v3 n;
	n[0] = v1[1] * v2[2] - v1[2] * v2[1];
	n[1] = v1[2] * v2[0] - v1[0] * v2[2];
	n[2] = v1[0] * v2[1] - v1[1] * v2[0];
	double nn;
	nn = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	n[0] /= nn;
	n[1] /= nn;
	n[2] /= nn;
	dprintf( "normal is %f,%f,%f\n", n[0], n[1], n[2]);
	// t = ((p0 - l0) . n) / (l . n)
	double num;
	v3 p;
	p[0] = p0[0] - l0[0];
	p[1] = p0[1] - l0[1];
	p[2] = p0[2] - l0[2];
	num = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
	dprintf( "num is %f\n", num);
	double den;
	den = l[0] * n[0] + l[1] * n[1] + l[2] * n[2];
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
#if 0
		double _t = num / den;
		dprintf( "one intersection at %f\n", _t);
		double x = l0[0] + _t * l[0];
		double y = l0[1] + _t * l[1];
		double z = l0[2] + _t * l[2];
		dprintf( "inters is %f,%f,%f\n", x, y, z);
#endif
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

int main()
{
	printf( "hello plane\n");
	v3 p0 = { 1, 0, 0 };
	v3 p1 = { 0, 1, 0 };
	v3 p2 = { 0, 0, 1 };
	v3 e = { 5, 5, 5 };
	v3 v = { -1, -1, -1 };
	double t = 0;
	int i, j, w, h;
	w = 20;
	h = 20;
	for (j = 0; j < h; j++)
	{
		for (i = 0; i < w; i++)
		{
			v3 _v;
			//_v[0] = v[0] + ((double)i - w / 2) / w;
			_v[0] = v[0] - 0.5 + (double)1.0 * (double)i / (w - 1);
			//_v[1] = v[1] + ((double)j - h / 2) / h;
			_v[1] = v[1] - 0.5 + 1.0 * (double)(h - j - 1) / (h - 1);
			_v[2] = v[2];
			int res = intersec_plane( p0, p1, p2, e, _v, &t);
			dprintf( "result=%d t=%f\n", res, t);
			if (res)
			{
#if 0
				double x = e[0] + t * v[0];
				double y = e[1] + t * v[1];
				double z = e[2] + t * v[2];
				dprintf( "inters is %f,%f,%f\n", x, y, z);
#endif
			}
			printf( "%c", res ? '1' : '0');
		}
		printf( "\n");
	}
	return 0;
}
