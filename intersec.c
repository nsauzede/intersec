#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <time.h>

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
#if 0
	{ .cx = 0.05, .cy = 0.1, .cz = 0.5, .sr = 0.1, .r = 1.0, .g = 0.0, .b = 0.0 },
	{ .cx = -0.05, .cy = 0.0, .cz = 0.6, .sr = 0.05, .r = 0.0, .g = 1.0, .b = 0.0 },
	{ .cx = 0.1, .cy = 0.0, .cz = 0.7, .sr = 0.03, .r = 0.0, .g = 0.0, .b = 1.0 },
#elif 1
#define SR 0.05
#define R 0.1
#define G 0.1
#define B 0.1
#define CX 0.05
#define CY 0.08
#define CZ 0.1
#if 1 /* ioccc ray */
	{ .cx = -7*CX, .cy = 5*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = G, .b = B },
	{ .cx = -3*CX, .cy = 5*CY, .cz = 0*CZ, .sr = SR, .r = 1*R, .g = 10*G, .b = B },
	{ .cx = 2*CX, .cy = 5*CY, .cz = 0*CZ, .sr = SR, .r = 1*R, .g = G, .b = 10*B },
	{ .cx = 5*CX, .cy = 5*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = G, .b = B },
	{ .cx = 8*CX, .cy = 5*CY, .cz = 0*CZ, .sr = SR, .r = 1*R, .g = 10*G, .b = B },

	{ .cx = -7*CX, .cy = 4*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = 1*G, .b = 10*B },
	{ .cx = -4*CX, .cy = 4*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = 1*G, .b = B },
	{ .cx = -2*CX, .cy = 4*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = 10*G, .b = B },
	{ .cx = 1*CX, .cy = 4*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = 1*G, .b = 10*B },
	{ .cx = 4*CX, .cy = 4*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = 1*G, .b = B },
	{ .cx = 7*CX, .cy = 4*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = 10*G, .b = B },

	{ .cx = -7*CX, .cy = 3*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = G, .b = 10*B },
	{ .cx = -4*CX, .cy = 3*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = G, .b = 1*B },
	{ .cx = -2*CX, .cy = 3*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = 10*G, .b = 1*B },
	{ .cx = 1*CX, .cy = 3*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = G, .b = 10*B },
	{ .cx = 4*CX, .cy = 3*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = G, .b = 1*B },
	{ .cx = 7*CX, .cy = 3*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = 10*G, .b = 1*B },

	{ .cx = -7*CX, .cy = 2*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = G, .b = 10*B },
	{ .cx = -3*CX, .cy = 2*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = G, .b = B },
	{ .cx = 2*CX, .cy = 2*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = 10*G, .b = B },
	{ .cx = 5*CX, .cy = 2*CY, .cz = 0*CZ, .sr = SR, .r = R, .g = G, .b = 10*B },
	{ .cx = 8*CX, .cy = 2*CY, .cz = 0*CZ, .sr = SR, .r = 10*R, .g = G, .b = B },

	{ .cx = 0*CX, .cy = 0*CY, .cz = 0*CZ, .sr = 1*SR, .r = R, .g = G, .b = B },
#endif
#else
	{ .cx = 0.0, .cy = -0.1, .cz = 0.0, .sr = 0.2, .r = 1.0, .g = 0.0, .b = 0.0 },
	{ .cx = 0.0, .cy = 0.1, .cz = 0.0, .sr = 0.2, .r = 0.0, .g = 1.0, .b = 0.0 },
	{ .cx = 0.0, .cy = 0.4, .cz = 0.0, .sr = 0.2, .r = 0.0, .g = 0.0, .b = 1.0 },
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

#define ATT_MIN 0.001
#define LEV_MAX -1

int sky_color( double ex, double ey, double ez, double _vx, double _vy, double _vz, double *r, double *g, double *b)
{
	double th, ph, n, coef2;
	n = sqrt( _vx * _vx + _vy * _vy + _vz * _vz);
	th = acos( _vz / n);
	ph = acos( _vx / n) * 180 / M_PI;
	coef2 = (ph + 90) / 180;
	if (coef2 > 1.0)
		coef2 = 1.0;
	if (coef2 < 0.0)
		coef2 = 0.0;
	*r = 1.0 - coef2;
	*g = 0.5;
	*b = coef2;

	return 0;
}

unsigned long traced = 0;
unsigned long intersected = 0;
unsigned long reflected = 0;
int level_max_reached = -1;
int level_max = LEV_MAX;
int traceray( int level, double ex, double ey, double ez, double _vx, double _vy, double _vz, double *r, double *g, double *b, char *pix, double att)
{
	traced++;
	if (level > level_max_reached)
		level_max_reached = level;
	unsigned long smin = -1, s;
	double tmin = BIG;
	double rmin = 0.0, gmin = 0.0, bmin = 0.0;
	for (s = 0; s < nsph; s++)
	{
		double cx, cy, cz, sr;
		cx = spheres[s].cx; cy = spheres[s].cy; cz = spheres[s].cz; sr = spheres[s].sr;
//		printf( "sphere: c(%f;%f;%f) r(%f)\n", cx, cy, cz, sr);
		int res;
		double t = 0;
		res = intersec_sphere( cx, cy, cz, sr, ex, ey, ez, _vx, _vy, _vz, &t);
//		printf( "%c", res ? 's' : '.');
		if (res)
		{
			if (t < tmin)
			{
				smin = s;
				tmin = t;
				rmin = spheres[s].r;
				gmin = spheres[s].g;
				bmin = spheres[s].b;
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
		sky_color( ex, ey, ez, _vx, _vy, _vz, &rmin, &gmin, &bmin);
//		printf( "sky color %f %f %f\n", rmin, gmin, bmin);
				if ((rmin > 1.0) || (gmin > 1.0) || (bmin > 1.0) || (rmin < 0.0) || (gmin < 0.0) || (bmin < 0.0))
				{
					printf( "boom sky color overflow r=%f g=%f b=%f\n", rmin, gmin, bmin);fflush( stdout);
					getchar();
					exit( 3);
				}
	}
	else				// object -> ambient
	{
		intersected++;
//		printf( "object\n");
		if (att > 1.0)
			att = 1.0;
		rmin *= att;
		gmin *= att;
		bmin *= att;
				if ((rmin > 1.0) || (gmin > 1.0) || (bmin > 1.0) || (rmin < 0.0) || (gmin < 0.0) || (bmin < 0.0))
				{
					printf( "boom object color overflow r=%f g=%f b=%f\n", rmin, gmin, bmin);fflush( stdout);
					getchar();
					exit( 3);
				}
// compute reflexions here
		double rx, ry, rz;		// coord of intersec
		rx = ex + _vx * tmin;
		ry = ey + _vy * tmin;
		rz = ez + _vz * tmin;
		double n;
//		dprintf( "inters with sphere %d at (%f,%f,%f)", smin, rx, ry, rz);
		double nvx, nvy, nvz;	// normal vect
		nvx = rx - spheres[smin].cx;
		nvy = ry - spheres[smin].cy;
		nvz = rz - spheres[smin].cz;
		n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
		nvx /= n;
		nvy /= n;
		nvz /= n;
		double rvx, rvy, rvz;	// vect of reflected ray
		double dot = 1.0;
		dot = (_vx *nvx + _vy * nvy + _vz * nvz);
		rvx = _vx - 2 * dot * nvx;
		rvy = _vy - 2 * dot * nvy;
		rvz = _vz - 2 * dot * nvz;
		n = sqrt( rvx * rvx + rvy * rvy + rvz * rvz);
		rvx /= n;
		rvy /= n;
		rvz /= n;
//		dprintf( "refl vec (%f,%f,%f)\n", rvx, rvy, rvz);
		double rr = 0.0, rg = 0.0, rb = 0.0;		// reflected color to add
		if (1
			&& (att >= ATT_MIN)
			&& ((level < level_max) || (level_max == -1))
		)
		{
			reflected++;
			traceray( level + 1, rx, ry, rz, rvx, rvy, rvz, &rr, &rg, &rb, pix, att * 0.9);
#if 0
			double ratt = 0.5;
			rmin = (rmin + rr * ratt) / (1.0 + ratt);
			gmin = (gmin + rg * ratt) / (1.0 + ratt);
			bmin = (bmin + rb * ratt) / (1.0 + ratt);
#elif 1
			double iatt = 1.0, ratt = 0.9;
			rmin = (rmin * iatt + rr * ratt) / (iatt + ratt);
			gmin = (gmin * iatt + rg * ratt) / (iatt + ratt);
			bmin = (bmin * iatt + rb * ratt) / (iatt + ratt);
#else
			rmin = rr;
			gmin = rg;
			bmin = rb;
#endif
		}
	}
	*r = rmin;
	*g = gmin;
	*b = bmin;
	return 0;
}

int main( int argc, char *argv[])
{
	double ex, ey, ez, vx, vy, vz;
	
	ex = 0; ey = 0; ez = 1; vx = 0.0; vy = 0.0; vz = -1;
	unsigned long w, h;
	w = W; h = H;
	double winw, winh;
	winw = 1.0; winh = 1.0;
	int do_tga = 1;
	int do_txt = 0;
	int do_att = 1;
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
	unsigned long old_t = time( 0);
	for (j = 0; j < h; j++)
	{
		double _vx, _vy, _vz;
		_vy = vy + winh / 2 - winh * (double)j / (h - 1);
		_vz = vz;

		for (i = 0; i < w; i++)
		{
			int percent = (double)100.0 * (j * w + i) / (w * h);
			static int old_percent = -1;
			if (percent != old_percent)
			{
				old_percent = percent;
				if (!(old_percent % 10))
				{
					printf( " %d%%", old_percent);fflush( stdout);
				}
			}
			char pix = '.';
			_vx = vx - winw / 2 + winw * (double)i / (w - 1);

			double r = 0.0, g = 0.0, b = 0.0;

			traceray( 0, ex, ey, ez, _vx, _vy, _vz, &r, &g, &b, &pix, do_att);
			dprintf( "  (r=%f g=%f b=%f)", r, g, b);
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
				unsigned char cr, cg, cb;
				cr = (double)255 * r;
				cg = (double)255 * g;
				cb = (double)255 * b;
				dprintf( "{%02x:%02x:%02x}", cr, cg, cb);
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
		dprintf( "\n");
	}
	unsigned long cur_t = time( 0);
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
	dprintf( "tga_size=%lu\n", (unsigned long)tga_size);
	dprintf( "tga_index=%lu\n", (unsigned long)tga_index);
	printf( "traced=%lu intersected=%lu\n", traced, intersected);
	printf( "reflected=%lu level_max=%d level_max_reached=%d\n", reflected, level_max, level_max_reached);
	unsigned long duration = cur_t - old_t;
	if (!duration)
		duration = 1;
	unsigned long perf = (w * h) / duration;
	printf( "perf : %lu ray/s duration=%lus (old_t=%lu cur_t=%lu)\n", perf, duration, old_t, cur_t);
	
	return 0;
}
