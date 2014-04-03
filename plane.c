#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#define dprintf(...) do{}while(0)

typedef double v1_t;
typedef v1_t v3_t[3];

v1_t* cross3( v3_t dest, v3_t src1, v3_t src2)
{
	dest[0] = src1[1] * src2[2] - src1[2] * src2[1];
	dest[1] = src1[2] * src2[0] - src1[0] * src2[2];
	dest[2] = src1[0] * src2[1] - src1[1] * src2[0];
	return dest;
}

v1_t* sum3( v3_t dest, v3_t src1, v3_t src2)
{
	dest[0] = src1[0] + src2[0];
	dest[1] = src1[1] + src2[1];
	dest[2] = src1[2] + src2[2];
	return dest;
}

v1_t* diff3( v3_t dest, v3_t src1, v3_t src2)
{
	dest[0] = src1[0] - src2[0];
	dest[1] = src1[1] - src2[1];
	dest[2] = src1[2] - src2[2];
	return dest;
}

v1_t* div3_t( v3_t n, v1_t nn)
{
	n[0] /= nn;
	n[1] /= nn;
	n[2] /= nn;
	return n;
}

v1_t* mult3( v3_t n, v1_t nn)
{
	n[0] *= nn;
	n[1] *= nn;
	n[2] *= nn;
	return n;
}

v1_t dot3( v3_t p, v3_t n)
{
	return p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
}

v1_t norm3( v3_t n)
{
	return sqrt( dot3( n, n));	
}

// plane : (p - po) . n = 0
// line : p = dl + l0
int intersec_plane( v3_t p0, v3_t p1, v3_t p2, v3_t l0, v3_t l, v1_t *pt)
{
	int result = 0;
	dprintf( "l0 is %f,%f,%f l is %f,%f,%f\n", l0[0], l0[1], l0[2], l[0], l[1], l[2]);
	dprintf( "p0 is %f,%f,%f p1 is %f,%f,%f p2 is %f,%f,%f\n", p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
	// plane normal :
	v3_t v1, v2;
	diff3( v1, p1, p0);
	diff3( v2, p2, p0);
	dprintf( "v1 is %f,%f,%f v2 is %f,%f,%f\n", v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
	v3_t n;
	cross3( n, v1, v2);
	v1_t nn;
	nn = norm3( n);
	div3_t( n, nn);
	dprintf( "normal is %f,%f,%f\n", n[0], n[1], n[2]);
	// t = ((p0 - l0) . n) / (l . n)
	v1_t num;
	v3_t p;
	diff3( p, p0, l0);
	num = dot3( p, n);
	dprintf( "num is %f\n", num);
	v1_t den;
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
		v1_t det;
		v1_t a, b, c, d, e, f, g, h, i;
		v3_t lb;
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
			v1_t A, B, C, D, E, F, G, H, I;
			A = (e * i - f * h);
			B = -(d * i - f * g);
			C = (d * h - e * g);
			D = -(b * i - c * h);
			E = (a * i - c * g);
			F = -(a * h - b * g);
			G = (b * f - c * e);
			H = -(a * f - c * d);
			I = (a * e - b * d);
			v1_t t, u = 0, v = 0;
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
int _solvetri( v1_t a, v1_t b, v1_t c, v1_t *t1, v1_t *t2)
{
	int result = 0;
	v1_t d, sd;
	d = b * b - 4 * a * c;
	sd = sqrt( d);
	
	if (d > 0)
	{
		*t1 = (-b - sd) / 2 / a;
		*t2 = (-b + sd) / 2 / a;
#if 0
		if ((*t1 < 0) || (*t2 < 0))
		asm volatile( "int $3");
#endif
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
#define solvetri _solvetri
#else
#define solvetri __solvetri
#endif

int __solvetri( v1_t a, v1_t b, v1_t c, v1_t *_t1, v1_t *_t2)
{
	v1_t t1 = 0, t2 = 0;
	int sol = _solvetri( a, b, c, &t1, &t2);
	if (sol > 0)
	{
		if (t1 < 0)
		{
			sol--;
			if (sol)
			{
				if (t2 >= 0)
				{
					t1 = t2;
				}
				else
					sol--;
			}
		}
		else if (sol > 1)
		{
			if (t2 < 0)
				sol--;
		}
	}
	if (_t1)
		*_t1 = t1;
	if (_t2)
		*_t2 = t2;
	return sol;
}

int intersec_sphere( v3_t cs, v1_t sr, v3_t e, v3_t v, v1_t *tmin, v1_t *tmax)
{
	v1_t ex = e[0];
	v1_t ey = e[1];
	v1_t ez = e[2];
	v1_t cx = cs[0];
	v1_t cy = cs[1];
	v1_t cz = cs[2];
	v1_t vx = v[0];
	v1_t vy = v[1];
	v1_t vz = v[2];
	
	int result = 0;
	v1_t a, b, c;
	v1_t t1, t2;
	
	v1_t tx, ty, tz;
	v1_t n;
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
			v1_t temp = t1;
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

#define USE_SOLID
#ifdef USE_SOLID
int use_solid = 1;
#endif

int intersec_cyl( v3_t cy, v3_t e, v3_t v, v1_t *tmin, v1_t *tmax)
{
	int result = 0;

	v1_t a, b, c;
	v1_t t1, t2;
	int sol;

	v1_t z0 = 0;
#define DO_TRANSLATE
#ifdef DO_TRANSLATE
	v1_t x0 = 0;
	v1_t y0 = 0;
	x0 = cy[0];
	y0 = cy[1];
	z0 = cy[2];
#endif

#define DO_ROTATE
#ifdef DO_ROTATE
	//v1_t rx = cy[3 * 1 + 0] * (double)M_PI / 180.0;
	v1_t rx = +0*45.0 * (double)M_PI / 180.0;
	//v1_t ry = cy[3 * 1 + 1] * (double)M_PI / 180.0;
	v1_t ry = -0*45.0 * (double)M_PI / 180.0;
	//v1_t rz = cy[3 * 1 + 2] * (double)M_PI / 180.0;
	v1_t rz = +0*45.0 * (double)M_PI / 180.0;
#endif

	v1_t Hc = cy[3 * 2 + 0];
	v1_t Rc = cy[3 * 2 + 1];
	v1_t Vx = v[0];
	v1_t Vy = v[1];
	v1_t Vz = v[2];
	v1_t lx0 = e[0];
	v1_t ly0 = e[1];
	v1_t lz0 = e[2];

#ifdef DO_TRANSLATE
	lx0 -= x0;
	ly0 -= y0;
	lz0 -= z0;
#endif

#ifdef DO_ROTATE
	v1_t xx = 0, yy = 0, zz = 0;
	rx = rx;
	rz = rz;
	ry = ry;
	xx = xx;
	zz = zz;
	yy = yy;

#if 1
	yy = ly0 * cos(rx) - lz0 * sin(rx);
	zz = ly0 * sin(rx) + lz0 * cos(rx);
	ly0 = yy;
	lz0 = zz;
	yy = Vy * cos(rx) - Vz * sin(rx);
	zz = Vy * sin(rx) + Vz * cos(rx);
	Vy = yy;
	Vz = zz;
#endif

#if 1
	xx = lx0 * cos(ry) + lz0 * sin(ry);
	zz = -lx0 * sin(ry) + lz0 * cos(ry);
	lx0 = xx;
	lz0 = zz;
	xx = Vx * cos(ry) + Vz * sin(ry);
	zz = -Vx * sin(ry) + Vz * cos(ry);
	Vx = xx;
	Vz = zz;
#endif

#if 0
	xx = lx0 * cos(rz) - ly0 * sin(rz);
	yy = lx0 * sin(rz) + ly0 * cos(rz);
	lx0 = xx;
	lz0 = yy;
	xx = Vx * cos(rz) - Vy * sin(rz);
	yy = Vx * sin(rz) + Vy * cos(rz);
	Vx = xx;
	Vy = yy;
#endif

#endif

	a = Vx * Vx + Vy * Vy;
	b = 2 * lx0 * Vx + 2 * ly0 * Vy;
	c = lx0 * lx0 + ly0 * ly0 - Rc * Rc;

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
	v1_t foo1 = 0.0;
	v1_t foo2 = 0.0;
	v1_t foo = 0.0;
	if (sol 
			//&& (fabs(t1) > SMALL)
			)
	{
		int hit = 0;

		v1_t zmin = z0;
		v1_t zmax = z0 + Hc;
		v1_t z1 = lz0 + Vz * t1;
		v1_t z2 = lz0 + Vz * t2;
#ifdef USE_SOLID
		if (use_solid)
		{
		v1_t z1a = z1 - zmax;
		v1_t z2a = z2 - zmax;

		if ((z1 < zmin) && (z2 > zmin))
		{
			t1 = (zmin - lz0) / Vz;
			hit = 1;
			foo = lz0 + Vz * t1;
			result = 1;
		}
		else
		if ((z1a * z2a) < 0)
		{
			t1 = (zmax - lz0) / Vz;
			hit = 1;
			foo = lz0 + Vz * t1;
			result = 1;
		}
		}
		if (!hit)
#endif
		{
		hit = (z1 >= zmin) && (z1 <= zmax);
		if (hit)
			foo = z1;
		foo1 = z1;
		foo2 = z2;
		if (!hit && ((sol > 1) && (fabs(t2) > SMALL)))
		{
			hit = (z2 >= zmin) && (z2 <= zmax);
			if (hit)
			{
				foo = z2;
				v1_t temp = t2;
				t1 = t2;
				t2 = temp;
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
#if 0
	printf( "%3.0f/", foo1);
	printf( "%3.0f", foo2);
#else
	foo1 = foo1;
	foo2 = foo2;
#endif
	foo = foo;
	printf( "%3.0f", foo);
	//printf( "%d", result);
	//printf( "%d", sol);

	return result;
}

int intersec_quad( v3_t obj, v3_t e, v3_t v, v1_t *tmin, v1_t *tmax)
{
	int result = 0;

	v1_t a, b, c;
	v1_t t1 = 0, t2 = 0;
	int sol;

	v1_t z0 = 0;
#define DO_TRANSLATE
#ifdef DO_TRANSLATE
	v1_t x0 = 0;
	v1_t y0 = 0;
	x0 = obj[0];
	y0 = obj[1];
	z0 = obj[2];
#endif

#define DO_ROTATE
#ifdef DO_ROTATE
	//v1_t rx = cy[3 * 1 + 0] * (double)M_PI / 180.0;
	v1_t rx = +0*45.0 * (double)M_PI / 180.0;
	//v1_t ry = cy[3 * 1 + 1] * (double)M_PI / 180.0;
	v1_t ry = +0*25.0 * (double)M_PI / 180.0;
	//v1_t rz = cy[3 * 1 + 2] * (double)M_PI / 180.0;
	v1_t rz = +0*45.0 * (double)M_PI / 180.0;
#endif

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

#ifdef DO_TRANSLATE
	lx0 -= x0;
	ly0 -= y0;
	lz0 -= z0;
#endif

#ifdef DO_ROTATE
	v1_t xx = 0, yy = 0, zz = 0;
	rx = rx;
	rz = rz;
	ry = ry;
	xx = xx;
	zz = zz;
	yy = yy;

#if 1
	yy = ly0 * cos(rx) - lz0 * sin(rx);
	zz = ly0 * sin(rx) + lz0 * cos(rx);
	ly0 = yy;
	lz0 = zz;
	yy = Vy * cos(rx) - Vz * sin(rx);
	zz = Vy * sin(rx) + Vz * cos(rx);
	Vy = yy;
	Vz = zz;
#endif

#if 1
	xx = lx0 * cos(ry) + lz0 * sin(ry);
	zz = -lx0 * sin(ry) + lz0 * cos(ry);
	lx0 = xx;
	lz0 = zz;
	xx = Vx * cos(ry) + Vz * sin(ry);
	zz = -Vx * sin(ry) + Vz * cos(ry);
	Vx = xx;
	Vz = zz;
#endif

#if 0
	xx = lx0 * cos(rz) - ly0 * sin(rz);
	yy = lx0 * sin(rz) + ly0 * cos(rz);
	lx0 = xx;
	lz0 = yy;
	xx = Vx * cos(rz) - Vy * sin(rz);
	yy = Vx * sin(rz) + Vy * cos(rz);
	Vx = xx;
	Vy = yy;
#endif

#endif

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
	v1_t foo1 = 0.0;
	v1_t foo2 = 0.0;
	v1_t foo = 0.0;
		v1_t z1 = 0;
		v1_t z2 = 0;

#if 0
		if (t1 < 0)
			asm volatile("int $3");
#endif
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
			foo = z1;
		}
		else
		{
			v1_t zmin = z0;
			v1_t zmax = z0 + lim;
#ifdef USE_SOLID
			if (use_solid)
			{
				v1_t z1a = z1 - zmax;
				v1_t z2a = z2 - zmax;

				if ((z1 < zmin) && (z2 > zmin))
				{
					t1 = (zmin - lz0) / Vz;
					hit = 1;
					foo = lz0 + Vz * t1;
					result = 1;
				}
				else
					if ((z1a * z2a) < 0)
					{
						t1 = (zmax - lz0) / Vz;
						hit = 1;
						foo = lz0 + Vz * t1;
						result = 1;
					}
			}
			if (!hit)
#endif
			{
				hit = (z1 >= zmin) && (z1 <= zmax);
				if (hit)
					foo = z1;
				foo1 = z1;
				foo2 = z2;
				if (!hit && ((sol > 1) && (fabs(t2) > SMALL)
							))
				{
					hit = (z2 >= zmin) && (z2 <= zmax);
					if (hit)
					{
						foo = z2;
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
#if 0
	printf( "%3.0f/", foo1);
	printf( "%3.0f", foo2);
#else
	foo1 = foo1;
	foo2 = foo2;
#endif
	foo = foo;
//	foo = sol;
//	foo = z1;
	printf( "%3.0f", foo);
//	printf( "[%4.0f/%4.0f]", t1, t2);
	//printf( "%d", result);
	//printf( "%d", sol);

	return result;
}

#define EPS 0.1
#define BIG 1000.0
#define COMP_EPS(x,val) (x >= ((val) - EPS) && x <= ((val) + EPS))
// camera is : eye coordinate (vector e)..
#if 0
#define EE 5
	v3_t e = { -2*EE, 0*EE, 16*EE };
// ..a direction (vector v)..
#define V EE/2
	v3_t v = { -0*V, -0*V, -2*V };
#elif 0
#define EE 20
	v3_t e = { 0*EE, 0*EE, -2*EE };
// ..a direction (vector v)..
#define V EE/2
	v3_t v = { -0*V, -0*V, 2*V };
#elif 0
#define EE 20
	v3_t e = { 0*EE, 0*EE, 2*EE };
// ..a direction (vector v)..
#define V EE/2
	v3_t v = { -0*V, -0*V, -2*V };
#elif 0
#define EE 20
	v3_t e = { 2*EE, -2*EE, 2*EE };
// ..a direction (vector v)..
#define V EE/2
	v3_t v = { -2*V, 2*V, -2*V };
#elif 1
#define EE 20
	v3_t e = { 2*EE, 2*EE, 2*EE };
// ..a direction (vector v)..
#define V EE/2
	v3_t v = { -2*V, -2*V, -2*V };
#else
#define EE 50
	v3_t e = { 0*EE, 0*EE, 1*EE };
// ..a direction (vector v)..
#define V EE/2
	v3_t v = { -0*V, -0*V, -1*V };
#endif
// ..and an "up" (vector up) (camera "head" rotation, default pointing to the "sky")
v3_t up = { -0, +1, 0};
// camera screen size
#if 1
int w = 44, h = 80;
#elif 1
int w = 22, h = 40;
#elif 1
int w = 180, h = 80;
#else
//int w = 640, h = 480;
//int w = 1024, h = 768;
int w = 1920, h = 1080;
#endif

// 3d scene :
//#define LOAD_SCENE
// multiple facets (3 vertex3 + 1 color3 each)
#define F 40
v3_t _facets[] = {
#if 0
		{ 0, 0, 0 },
		{ F, 0, 0 },
		{ 0, F, 0 },
		{ 1, 0, 0 },	// color

		{ 0, 0, 0 },
		{ 0, F, 0 },
		{ 0, 0, F },
		{ 0, 1, 0 },	// color

		{ 0, 0, 0 },
		{ 0, 0, F },	// WARNING : visible faces must have normal pointing to eye !
		{ F, 0, 0 },
		{ 0, 0, 1 },	// color
#endif
};
int nfacets = sizeof( _facets) / sizeof( _facets[0]) / 4;
v3_t *facets = _facets;

// multiple spheres
#define S 20
v3_t _spheres[] = {
#if 0
	{ S, S, S },
	{ S/4, 0, 0 },	// radius
	{ 1, 0, 1 },	// color
#endif
};
int nspheres = sizeof( _spheres) / sizeof( _spheres[0]) / 3;
v3_t *spheres = _spheres;

// multiple cyls
#define CH 10
#define CR 5
v3_t _cyls[] = {
#if 1
	{ 0, 0, 0 },	// center
	{ 0, 0, 0 },	// axis
	{ CH, CR, 0 },	// CH, CR
	{ 0, 1, 0 },	// color
#endif
};
int ncyls = sizeof( _cyls) / sizeof( _cyls[0]) / 4;
v3_t *cyls = _cyls;

// multiple quads
#if 0 // Z-cone
#define A +1		// +
#define B +1		// +
#define C -1		// -
#define D +0		// 0
#define E 0			// N/A
#define LIM 10
#elif 0 // Y-cone
#define A +1		// +
#define B -1		// -
#define C +1		// +
#define D +0		// 0
#define E 0			// N/A
#define LIM 0
#elif 0 // Y-cylinder
#define A +1		// +
#define B +0		// 0
#define C +1		// +
#define D +25		// +
#define E 0			// N/A
#define LIM 0
#elif 1 // Z-cylinder
#define A +1		// +
#define B +1		// +
#define C +0		// 0
#define D +25		// +
#define E 0			// N/A
#define LIM 10
#endif
v3_t _quads[] = {
#if 0
	{ 0, 0, 0 },	// translation
	{ 0, 0, 0 },	// rotation
	{ A, B, C },	// 
	{ D, E, LIM },	// 
	{ 0, 1, 0 },	// color
#endif
};
int nquads = sizeof( _quads) / sizeof( _quads[0]) / 5;
v3_t *quads = _quads;

void traceray( v3_t e, v3_t _v, v3_t col)
{
	v1_t *pcol = 0;
	v1_t tmin = BIG;
	v1_t t;
	v1_t *p0, *p1, *p2;
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
			pcol = &facets[i * 4 + 3][0];
		}
	}
	for (i = 0; i < nspheres; i++)
	{
		p0 = spheres[i * 3 + 0];		// sphere center
		p1 = spheres[i * 3 + 1];		// sphere radius
		t = 0;
		res = intersec_sphere( p0, *p1, e, _v, &t, 0);
		dprintf( "result=%d t=%f\n", res, t);
		if (res && (t < tmin))
		{
			tmin = t;
			pcol = &spheres[i * 3 + 2][0];
		}
	}
	for (i = 0; i < ncyls; i++)
	{
		p0 = cyls[i * 4 + 0];		// cyl center
//		p1 = cyls[i * 4 + 1];		// cyl axis
//		p2 = cyls[i * 4 + 2];		// cyl H and R
		t = 0;
		res = intersec_cyl( p0, e, _v, &t, 0);
		dprintf( "result=%d t=%f\n", res, t);
		dprintf( " %d", res);
		if (res && (t < tmin))
		{
			tmin = t;
			pcol = &cyls[i * 4 + 3][0];
		}
	}
	for (i = 0; i < nquads; i++)
	{
		p0 = quads[i * 5 + 0];		// quad
		t = 0;
		res = intersec_quad( p0, e, _v, &t, 0);
		dprintf( "result=%d t=%f\n", res, t);
		dprintf( " %d", res);
		if (res && (t < tmin))
		{
			tmin = t;
			pcol = &quads[i * 5 + 4][0];
		}
	}
	if (pcol)
		memcpy( col, pcol, 3 * sizeof( col[0]));
}

int main( int argc, char *argv[])
{
	char *scene = 0;
	
	printf( "hello plane\n");
	int arg = 1;
	if (arg < argc)
	{
#ifdef USE_SOLID
		if (!strcmp( "--hollow", argv[arg]))
		{
			use_solid = 0;
			arg++;
		}
#endif
		scene = argv[arg++];
	}
#ifdef USE_SOLID
	printf( "use_solid=%d\n", use_solid);
#endif
	
	int i, j;
	if (scene)	// load scene ?
	{
	FILE *in = fopen( scene, "rt");
	if (in)
	{
		nfacets = 0;
		nspheres = 0;
		char buf[1024];
		char *ptr;
		while (!feof( in))
		{
			if (!fgets( buf, sizeof( buf), in))
				break;
			ptr = strstr( buf, "endfacet");
			if (ptr)
				nfacets++;
		}
		printf( "read nfacets=%d\n", nfacets);
		rewind( in);
		if (nfacets)
		{
		facets = malloc( nfacets * sizeof(v3_t[4]));
		i = 0;
		j = 0;
#ifdef USE_NORMAL
		float normal[3] = { 0, 0, 0};
#endif
		while (!feof( in))
		{
			if (!fgets( buf, sizeof( buf), in))
				break;
#ifdef USE_NORMAL
			ptr = strstr( buf, "facet");
			if (ptr)
			{
				float p[3];
				sscanf( ptr, "facet normal %f %f %f", &p[0], &p[1], &p[2]);
				printf( "read normal: %f,%f,%f\n", p[0], p[1], p[2]);
				normal[0] = p[0];
				normal[1] = p[1];
				normal[2] = p[2];
			}
#endif
			ptr = strstr( buf, "vertex");
			if (ptr)
			{
				float p[3];
				sscanf( ptr, "%*s %f %f %f", &p[0], &p[1], &p[2]);
//				printf( "read p: %f,%f,%f\n", p[0], p[1], p[2]);
				facets[i * 4 + j][0] = p[0];
				facets[i * 4 + j][1] = p[1];
				facets[i * 4 + j][2] = p[2];
				if (++j >= 3)
				{
#ifndef USE_NORMAL
					facets[i * 4 + j][0] = !(i % 3);	// fake facet color with facet index r
					facets[i * 4 + j][1] = !((i + 1) % 3);	// g
					facets[i * 4 + j][2] = !((i + 2) % 3);	// b
#else
					facets[i * 4 + j][0] = normal[0];	// fake facet color with facet normal r
					facets[i * 4 + j][1] = normal[1];	// g
					facets[i * 4 + j][2] = normal[2];	// b
#endif
					i++;
					j = 0;
				}
			}
		}
		}
		fclose( in);
	}
	}
	printf( "nfacets=%d\n", nfacets);
	printf( "nspheres=%d\n", nspheres);
	printf( "ncyls=%d\n", ncyls);
	printf( "nquads=%d\n", nquads);
	// compute camera screen as three vectors of a plane : s0, s1 and s2
	v3_t s0, s1, s2;
	v3_t right;
	// normalize v
	div3_t( v, norm3( v));
	cross3( right, v, up);
	cross3( up, right, v);
	sum3( s0, e, v);
	sum3( s1, s0, up);
	sum3( s2, s0, right);

	for (j = 0; j < h; j++)
	{
		for (i = 0; i < w; i++)
		{

			v1_t u, v;
			v = -(v1_t)0.5 + (v1_t)1.0 * (v1_t)i / ((v1_t)w - 1);
			u = -(v1_t)0.5 + (v1_t)1.0 * (v1_t)(h - j - 1) / ((v1_t)h - 1);
			// compute a point s in camera screen plane based on {u,v}
			v3_t s;
			// p=p0+(p1-p0)u+(p2-p0)v
			s[0] = s0[0] + (s1[0] - s0[0]) * u + (s2[0] - s0[0]) * v;
			s[1] = s0[1] + (s1[1] - s0[1]) * u + (s2[1] - s0[1]) * v;
			s[2] = s0[2] + (s1[2] - s0[2]) * u + (s2[2] - s0[2]) * v;
			// compute final camera=> pixel vector _v
			v3_t _v;
			diff3( _v, s, e);

			v3_t col;
			memset( col, 0, sizeof( col));
			traceray( e, _v, col);

			char c = '?';
			if (COMP_EPS( col[0], 0.0) && COMP_EPS( col[1], 0.0) && COMP_EPS( col[2], 0.0))
				c = '.';
			else if (COMP_EPS( col[0], 1.0) && COMP_EPS( col[1], 0.0) && COMP_EPS( col[2], 0.0))
				c = 'R';
			else if (COMP_EPS( col[0], 0.0) && COMP_EPS( col[1], 1.0) && COMP_EPS( col[2], 0.0))
				c = 'G';
			else if (COMP_EPS( col[0], 0.0) && COMP_EPS( col[1], 0.0) && COMP_EPS( col[2], 1.0))
				c = 'B';
			else
			{
				if (col[0] > col[1])
				{
					if (col[0] > col[2])
						c = 'r';
					else
						c = 'b';
				}
				else if (col[1] > col[2])
					c = 'g';
				else
					c = 'b';
			}
			c = c;
			printf( "%c", c);
			fflush( stdout);
		}
		printf( "\n");
	}
	return 0;
}
