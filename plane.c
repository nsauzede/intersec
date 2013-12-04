#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#define dprintf(...) do{}while(0)

typedef double v3[3];

double* cross3( v3 dest, v3 src1, v3 src2)
{
	dest[0] = src1[1] * src2[2] - src1[2] * src2[1];
	dest[1] = src1[2] * src2[0] - src1[0] * src2[2];
	dest[2] = src1[0] * src2[1] - src1[1] * src2[0];
	return dest;
}

double* sum3( v3 dest, v3 src1, v3 src2)
{
	dest[0] = src1[0] + src2[0];
	dest[1] = src1[1] + src2[1];
	dest[2] = src1[2] + src2[2];
	return dest;
}

double* diff3( v3 dest, v3 src1, v3 src2)
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

double dot3( v3 p, v3 n)
{
	return p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
}

double norm3( v3 n)
{
	return sqrt( dot3( n, n));	
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

// camera is : eye coordinate (vector e)..
#define E 80
	v3 e = { 0, 0, E };
// ..a direction (vector v)..
#define V -1
	v3 v = { 0, 0, V };
// ..and an "up" (vector up) (camera "head" rotation, default pointing to the "sky")
v3 up = { 0, 1, 0};
// camera screen size
int w = 80, h = 40;
// 3d scene :
#if 0
// a simple triangle, 3 points3
v3 p0 = { 1, 0, 0 };
v3 p1 = { 0, 1, 0 };
v3 p2 = { 0, 0, 1 };
// a sphere
v3 cs = { 0, 0, 0 };
double sr = 0.3;
#else
#define LOAD_SCENE
// multiple facets (3 vertex3 + 1 color3 each)
#ifndef LOAD_SCENE
v3 facets[] = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 },
		{ 1, 0, 0 },	// color

		{ 0, 0, 1 },
		{ 0, 1, 0 },
		{ -1, 0, 0 },
		{ 0, 1, 0 },	// color
};
int nfacets = sizeof(facets) / sizeof(facets[0]) / 4;
#else
v3 *facets = 0;
int nfacets = 0;
#endif
// multiple spheres
v3 spheres[][3] = {
#if 0
	{
		{ 0, -0.5, 0 },
		{ 0.5, 0, 0 },	// radius
		{ 0, 0, 1 },	// color
	},
#endif
};
int nspheres = sizeof(spheres) / sizeof(spheres[0]);
#endif

#define EPS 0.1
#define BIG 1000.0
#define COMP_EPS(x,val) (x >= ((val) - EPS) && x <= ((val) + EPS))
void traceray( v3 e, v3 _v, v3 col)
{
	double *pcol = 0;
	char c = '.';
	double tmin = BIG;
	double t;
	double *p0, *p1, *p2;
	int res;
	int i;
	for (i = 0; i < nfacets; i++)
	{
		p0 = facets[i * 4 + 0];
		p1 = facets[i * 4 + 1];
		p2 = facets[i * 4 + 2];
		t = 0;
		res = intersec_plane( p0, p1, p2, e, _v, &t);
		dprintf( "result=%d t=%f\n", res, t);
		if (res && (t < tmin))
		{
			tmin = t;
			c = '0' + i;
			pcol = &facets[i * 4][3];
#if 0
			double x = e[0] + t * _v[0];
			double y = e[1] + t * _v[1];
			double z = e[2] + t * _v[2];
			dprintf( "inters is %f,%f,%f\n", x, y, z);
			if (COMP_EPS( x, 1.0) && COMP_EPS( y, 0.0) && COMP_EPS( z, 0.0))
			{
				c = 'R';
				col[0] = 1;
			}
			else if (COMP_EPS( x, 0.0) && COMP_EPS( y, 1.0) && COMP_EPS( z, 0.0))
			{
				c = 'G';
				col[1] = 1;
			}
			else if (COMP_EPS( x, 0.0) && COMP_EPS( y, 0.0) && COMP_EPS( z, 1.0))
			{
				c = 'B';
				col[2] = 1;
			}
#endif
		}
	}
	for (i = 0; i < nspheres; i++)
	{
		p0 = spheres[i][0];		// sphere center
		p1 = spheres[i][1];		// sphere radius
		t = 0;
		res = intersec_sphere( p0, *p1, e, _v, &t, 0);
		dprintf( "result=%d t=%f\n", res, t);
		if (res && (t < tmin))
		{
			tmin = t;
			c = '0' + i;
			pcol = &spheres[i][2][0];
		}
	}
	if (pcol)
		memcpy( col, pcol, sizeof( col));
	printf( "%c", c);
}

int main( int argc, char *argv[])
{
	char *scene = "scene.stl";
	
	printf( "hello plane\n");
	int arg = 1;
	if (arg < argc)
	{
		scene = argv[arg++];
	}
	
	int i, j;
#ifdef LOAD_SCENE
	// load scene ?
	if (scene)
	{
	FILE *in = fopen( scene, "rt");
	if (in)
	{
		char buf[1024];
		char *ptr;
		while (!feof( in))
		{
			fgets( buf, sizeof( buf), in);
			ptr = strstr( buf, "endfacet");
			if (ptr)
				nfacets++;
		}
		printf( "read nfacets=%d\n", nfacets);
		rewind( in);
		if (nfacets)
		{
		facets = malloc( nfacets * sizeof(v3[4]));
		i = 0;
		j = 0;
		while (!feof( in))
		{
			fgets( buf, sizeof( buf), in);
			ptr = strstr( buf, "vertex");
			if (ptr)
			{
				float p[3];
				sscanf( ptr, "%*s %f %f %f", &p[0], &p[1], &p[2]);
				printf( "read p: %f,%f,%f\n", p[0], p[1], p[2]);
				facets[i * 4 + j][0] = p[0];
				facets[i * 4 + j][1] = p[1];
				facets[i * 4 + j][2] = p[2];
				if (++j >= 3)
				{
					i++;
					j = 0;
				}
			}
		}
		}
		fclose( in);
	}
	}
#endif
	printf( "nfacets=%d\n", nfacets);
	// compute camera screen as three vectors of a plane : s0, s1 and s2
	v3 s0, s1, s2;
	v3 right;
	cross3( right, v, up);
	cross3( up, right, v);
	sum3( s0, e, v);
	sum3( s1, s0, up);
	sum3( s2, s0, right);

	for (j = 0; j < h; j++)
	{
		for (i = 0; i < w; i++)
		{

			double u, v;
			v = -(double)0.5 + (double)1.0 * (double)i / ((double)w - 1);
			u = -(double)0.5 + (double)1.0 * (double)(h - j - 1) / ((double)h - 1);
			// compute a point s in camera screen plane based on {u,v}
			v3 s;
			// p=p0+(p1-p0)u+(p2-p0)v
			s[0] = s0[0] + (s1[0] - s0[0]) * u + (s2[0] - s0[0]) * v;
			s[1] = s0[1] + (s1[1] - s0[1]) * u + (s2[1] - s0[1]) * v;
			s[2] = s0[2] + (s1[2] - s0[2]) * u + (s2[2] - s0[2]) * v;
			// compute final camera=> pixel vector _v
			v3 _v;
			diff3( _v, s, e);

			v3 col;
			traceray( e, _v, col);

			char c = ' ';
			if (col[0] > col[1]) // R>G
			{
				if (col[1] > col[2]) // G>B
				{
					c = 'R';	// R>G>B
				}
				else if (col[0] > col[2]) // R>B
				{
					c = 'r';	// R>B>G
				}
				else // B>G, B>R
				{
					c = 'b';	// B>R>G
				}
			}
			else if (col[1] > col[2]) // G>R, G>B
			{
				if (col[0] > col[2]) // R>B
				{
					c = 'G';	// G>R>B
				}
				else	// B>R
				{
					c = 'g';	// G>B>R
				}
			}
			else // G>R, B>G
			{
				c = 'B'; // B>G>R
			}
			c = c;
//			printf( "%c", c);
		}
		printf( "\n");
	}
	return 0;
}
