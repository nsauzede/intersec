#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>

#ifdef WIN32
#include <windows.h>
#define MAP_FAILED ((HANDLE)-1)
#else
#include <sys/mman.h>
#endif

#define W 600
#define H 600
#define SMALL	0.001
#define BIG		10.0

#ifdef DEBUG
#define dprintf(...) do{printf(__VA_ARGS__);}while(0)
#else
#define dprintf(...) do{}while(0)
#endif

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

typedef struct {
	double cx, cy, cz, sr;
	double r, g, b;
} sphere_t;

sphere_t spheres[] = {
	{ .cx = 0.05, .cy = 0.1, .cz = 0.5, .sr = 0.2, .r = 1.0, .g = 0.0, .b = 0.0 },
	{ .cx = 0.0, .cy = 0.0, .cz = 0.6, .sr = 0.1, .r = 0.0, .g = 1.0, .b = 0.0 },
	{ .cx = 0.1, .cy = 0.0, .cz = 0.7, .sr = 0.03, .r = 0.0, .g = 0.0, .b = 1.0 },
};
int nsph = sizeof (spheres) / sizeof (spheres[0]);

int traceray( double ex, double ey, double ez, double _vx, double _vy, double _vz, double *r, double *g, double *b, char *pix, int do_att)
{
	int s;
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
		double coef1, coef2;
		coef1 = 1.0 * _vx;
		coef2 = 1.0 * _vy;
		bmin = 1.0 - coef2;
		gmin = 0.5;
		rmin = coef2;
	}
	else				// object -> ambient
	{
		double att = 1.0;	// should be a vector(rgb) ?
		if (do_att)
			att = 1.0 / sqrt(tmin);
		rmin *= att;
		gmin *= att;
		bmin *= att;
	}
	*r = rmin;
	*g = gmin;
	*b = bmin;
	return 0;
}

int main( int argc, char *argv[])
{
	double ex, ey, ez, vx, vy, vz;
	
	ex = 0; ey = 0; ez = 1; vx = 0.05; vy = 0.1; vz = -1;
	int w, h;
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
		sscanf( argv[arg++], "%d", &w);
		if (argc > arg)
		{
			sscanf( argv[arg++], "%d", &h);
			if (argc > arg)
			{
				sscanf( argv[arg++], "%d", &do_tga);
			}
		}
	}
	if (do_tga)
	{
		tga_fd = open( "out.tga", O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
		if (tga_fd == -1)
		{
			perror( "open");
			exit( 2);
		}
		tga_size = w * h * bpp / 8 + 18;
		dprintf( "tga_size=%lu\n", (unsigned long)tga_size);
		ftruncate( tga_fd, tga_size);
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
		printf( "w=%d h=%d\n", w, h);
		printf( "eye: e(%f;%f;%f) v(%f;%f;%f)\n", ex, ey, ez, vx, vy, vz);
	}
	
	int i, j;
	for (j = 0; j < h; j++)
	{
		double _vx, _vy, _vz;
		_vy = vy + winh / 2 - winh * (double)j / (h - 1);
		_vz = vz;

		for (i = 0; i < w; i++)
		{
			char pix = '.';
			_vx = vx - winw / 2 + winw * (double)i / (w - 1);

			double r = 0.0, g = 0.0, b = 0.0;

			traceray( ex, ey, ez, _vx, _vy, _vz, &r, &g, &b, &pix, do_att);
			dprintf( "  (r=%f g=%f b=%f)", r, g, b);
			if (do_tga)
			{
				unsigned char cr, cg, cb;
				cr = ((double)255 * r);
				cg = (unsigned char)(255 * g);
				cb = (unsigned char)(255 * b);
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
	dprintf( "tga_index=%lu", (unsigned long)tga_index);
	
	return 0;
}
