#include <stdio.h>
#include <math.h>
#include <unistd.h>

#define W 80
#define H 60
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
};
int nsph = sizeof (spheres) / sizeof (spheres[0]);

int main( int argc, char *argv[])
{
	double ex, ey, ez, vx, vy, vz;
	
	ex = 0; ey = 0; ez = 1; vx = 0.05; vy = 0.1; vz = -1;
	int w, h;
	w = W; h = H;
	double winw, winh;
	winw = 1.0; winh = 1.0;
	int do_tga = 0;
	int do_att = 1;
	int bpp = 24;
	unsigned char tga_bdata;
	unsigned short tga_sdata;

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
#define TGA_BYTE(b) tga_bdata = b; write( 1, &tga_bdata, sizeof( tga_bdata))
#define TGA_SHORT(s) tga_sdata = s; write( 1, &tga_sdata, sizeof( tga_sdata))
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
	else
	{
		printf( "w=%d h=%d\n", w, h);
		printf( "eye: e(%f;%f;%f) v(%f;%f;%f)\n", ex, ey, ez, vx, vy, vz);
	}
	
	int i, j;
	for (j = 0; j < h; j++)
	{
		for (i = 0; i < w; i++)
		{
			int s;
			double tmin = BIG;
			char pix = '.';
			double rmin = 0.0, gmin = 0.0, bmin = 0.0;
			for (s = 0; s < nsph; s++)
			{
				double cx, cy, cz, sr;
				cx = spheres[s].cx; cy = spheres[s].cy; cz = spheres[s].cz; sr = spheres[s].sr;
//				printf( "sphere: c(%f;%f;%f) r(%f)\n", cx, cy, cz, sr);
				int res;
				double _vx, _vy, _vz;
				_vx = vx - winw / 2 + winw * (double)i / (w - 1);
				_vy = vy + winh / 2 - winh * (double)j / (h - 1);
				_vz = vz;
				double t = 0;
				res = intersec_sphere( cx, cy, cz, sr, ex, ey, ez, _vx, _vy, _vz, &t);
//				printf( "%c", res ? 's' : '.');
				if (res)
				{
					if (t < tmin)
					{
						tmin = t;
						rmin = spheres[s].r; gmin = spheres[s].g; bmin = spheres[s].b;
//						pix = '0' + s;
						if (spheres[s].r > 0)
							pix = 'r';
						else if (spheres[s].g > 0)
							pix = 'g';
						else if (spheres[s].b > 0)
							pix = 'b';
						else
							pix = 'k';
					}
				}
			}
			if (do_tga)
			{
				double att = 1.0;
				if (do_att)
					att = 1.0 / sqrt(tmin);
				TGA_BYTE( (unsigned char)(255 * rmin * att));
				TGA_BYTE( (unsigned char)(255 * gmin * att));
				TGA_BYTE( (unsigned char)(255 * bmin * att));
//				TGA_BYTE( 0);
			}
			else
				printf( "%c", pix);
		}
		if (!do_tga)
			printf( "\n");
	}
	
	return 0;
}
