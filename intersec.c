/*
 * intersec : small naive ray tracer
 * copyright 2012 Nicolas Sauzede (nsauzede@laposte.net)
 * This code is GPLv3
 *
 */
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#ifdef USETIMEOFDAY
#include <sys/time.h>
#else
#include <time.h>
#endif

#ifdef WIN32
#include <windows.h>
#define MAP_FAILED ((HANDLE)-1)
#else
#include <sys/mman.h>
#endif

#ifdef USE_SDL
#include <SDL.h>
#endif

#include "vec.h"

#if 0
#define W 60
#define H 60
#else
#define W 640
#define H 480
#endif

#define SMALL	0.001
#define BIG		10.0

#define WINSCALE 1.0

#ifdef DEBUG
#define dprintf(...) do{printf(__VA_ARGS__);}while(0)
#else
#define dprintf(...) do{}while(0)
#endif

typedef struct {
	double cx, cy, cz, sr;
	double r, g, b;
	double rdatt, gdatt, bdatt;		// diffuse
	double rratt, gratt, bratt;		// reflection
	double rfatt, gfatt, bfatt;		// refraction
	double rindex;					// refraction index
} sphere_t;

typedef double v3[3];

typedef v3 lamp_t;

typedef struct {
	sphere_t *spheres;
	int nspheres;
	v3 *facets;
	int nfacets;
	lamp_t *lamps;
	int nlamps;
} scene_t;

int load_scene( scene_t *scene, char *file)
{
	int i, j;
	FILE *in = fopen( file, "rt");
	if (!in)
	{
		printf( "couln't open scene file '%s' !\n", file);
		exit( 1);
	}
	else
	{
		enum { UNKNOWN = 0, STL = 1, POV = 2 } type= UNKNOWN;
		scene->nfacets = 0;
		scene->nspheres = 0;
		char buf[1024];
		char *ptr;
		while (!feof( in))
		{
			if (!fgets( buf, sizeof( buf), in))
				break;
			ptr = strstr( buf, "endfacet");
			if (ptr)
			{
				if (type &= ~STL)
				{
					printf( "bad scene file : type=%x and STL token found\n", type);
					exit( 1);
				}
				type = STL;
				scene->nfacets++;
			}
			ptr = strstr( buf, "camera");
			if (ptr)
			{
				if (type &= ~POV)
				{
					printf( "bad scene file : type=%x and POV token found\n", type);
					exit( 1);
				}
				type = POV;
			}
			ptr = strstr( buf, "sphere");
			if (ptr)
			{
				if (type &= ~POV)
				{
					printf( "bad scene file : type=%x and POV token found\n", type);
					exit( 1);
				}
				type = POV;
				scene->nspheres++;
			}
			ptr = strstr( buf, "light_source");
			if (ptr)
			{
				if (type &= ~POV)
				{
					printf( "bad scene file : type=%x and POV token found\n", type);
					exit( 1);
				}
				type = POV;
				scene->nlamps++;
			}
		}
		printf( "type=%x (%s)\n", type, type == STL ? "STL" : type == POV ? "POV" : "unknown");
		printf( "read nfacets=%d\n", scene->nfacets);
		printf( "read nspheres=%d\n", scene->nspheres);
		printf( "read nlamps=%d\n", scene->nlamps);
		rewind( in);
		if ((type == POV) && scene->nspheres)
		{
			if (scene->nlamps)
			{
				scene->lamps = malloc( scene->nlamps * sizeof(*scene->lamps));
				memset( scene->lamps, 0, scene->nlamps * sizeof(*scene->lamps));
			}
			scene->spheres = malloc( scene->nspheres * sizeof(*scene->spheres));
			memset( scene->spheres, 0, scene->nspheres * sizeof(*scene->spheres));
			i = -1;
			int l = 0;
			while (!feof( in))
			{
				float p[4];
				if (!fgets( buf, sizeof( buf), in))
					break;
				ptr = strstr( buf, "light_source");
				if (ptr)
				{
					ptr = strchr( buf, '<');
					printf( "light_source=[%s]\n", ptr);
					sscanf( ptr, "<%f, %f, %f>", &p[0], &p[1], &p[2]);
					printf( "read l: %f,%f,%f\n", p[0], p[1], p[2]);
					scene->lamps[l][0] = p[0];
					scene->lamps[l][1] = p[1];
					scene->lamps[l][2] = p[2];
					l++;
					continue;
				}
				ptr = strstr( buf, "sphere");
				if (ptr)
				{
					if (!fgets( buf, sizeof( buf), in))
						break;
					ptr = strchr( buf, '<');
					printf( "sphere=[%s]\n", ptr);
					sscanf( ptr, "<%f, %f, %f>, %f", &p[0], &p[1], &p[2], &p[3]);
					printf( "read s: %f,%f,%f %f\n", p[0], p[1], p[2], p[3]);
					i++;
					scene->spheres[i].cx = p[0];
					scene->spheres[i].cy = p[1];
					scene->spheres[i].cz = p[2];
					scene->spheres[i].sr = p[3];
					continue;
				}
				ptr = strstr( buf, "rgb");
				if (ptr)
				{
					ptr = strchr( buf, '<');
					printf( "rgb=[%s]\n", ptr);
					sscanf( ptr, "<%f, %f, %f>", &p[0], &p[1], &p[2]);
					printf( "read c: %f,%f,%f\n", p[0], p[1], p[2]);
					scene->spheres[i].r = p[0];
					scene->spheres[i].g = p[1];
					scene->spheres[i].b = p[2];
					continue;
				}
				ptr = strstr( buf, "diffuse");
				if (ptr)
				{
					printf( "diffuse=[%s]\n", ptr);
					sscanf( ptr, "diffuse %f", &p[0]);
					printf( "read d: %f\n", p[0]);
					scene->spheres[i].rdatt = p[0];
					scene->spheres[i].gdatt = p[0];
					scene->spheres[i].bdatt = p[0];
					continue;
				}
				ptr = strstr( buf, "reflection");
				if (ptr)
				{
					printf( "reflection=[%s]\n", ptr);
					sscanf( ptr, "reflection %f", &p[0]);
					printf( "read r: %f\n", p[0]);
					scene->spheres[i].rratt = p[0];
					scene->spheres[i].gratt = p[0];
					scene->spheres[i].bratt = p[0];
					continue;
				}
				ptr = strstr( buf, "specular");
				if (ptr)
				{
					printf( "specular=[%s]\n", ptr);
					sscanf( ptr, "specular %f", &p[0]);
					printf( "read r: %f\n", p[0]);
					scene->spheres[i].rfatt = p[0];
					scene->spheres[i].gfatt = p[0];
					scene->spheres[i].bfatt = p[0];
					continue;
				}
			}
		}
		if ((type == STL) && scene->nfacets)
		{
		scene->facets = malloc( scene->nfacets * sizeof(v3[4]));
		memset( scene->spheres, 0, scene->nspheres * sizeof(*scene->spheres));
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
				scene->facets[i * 4 + j][0] = p[0];
				scene->facets[i * 4 + j][1] = p[1];
				scene->facets[i * 4 + j][2] = p[2];
				if (++j >= 3)
				{
#ifndef USE_NORMAL
					scene->facets[i * 4 + j][0] = !(i % 3);	// fake facet color with facet index r
					scene->facets[i * 4 + j][1] = !((i + 1) % 3);	// g
					scene->facets[i * 4 + j][2] = !((i + 2) % 3);	// b
#else
					scene->facets[i * 4 + j][0] = normal[0];	// fake facet color with facet normal r
					scene->facets[i * 4 + j][1] = normal[1];	// g
					scene->facets[i * 4 + j][2] = normal[2];	// b
#endif
					i++;
					j = 0;
				}
			}
		}
		}
		fclose( in);
	}
	printf( "nfacets=%d\n", scene->nfacets);
	printf( "nspheres=%d\n", scene->nspheres);
////
	return 0;
}
scene_t scene;

sphere_t _spheres[] = {
#if 1
#if 1
	{ .cx = 0.05, .cy = 0.3, .cz = -0.8, .sr = 0.2, .r = 1.0, .g = 0.0, .b = 0.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#define D1 1.0
	{ .cx = -0.2, .cy = -0.1, .cz = -0.6, .sr = 0.1, .r = 0.0, .g = 1.0, .b = 0.0, .rdatt = D1, .gdatt = D1, .bdatt = D1, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#define D 0.8
#define R 0.0
	{ .cx = 0.5, .cy = -0.05, .cz = -0.7, .sr = 0.2, .r = 0.1, .g = 0.1, .b = 1.0, .rdatt = D, .gdatt = D, .bdatt = D, .rratt = R, .gratt = R, .bratt = R },
	{ .cx = 0.0, .cy = 0.0, .cz = 0.0, .sr = 0.03, .r = 1.0, .g = 1.0, .b = 1.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#endif
#define RFL 0.0
#define RFR 1.0
#define RINDEX 1.0
	{ .cx = 0.0, .cy = 0.0, .cz = 0.3, .sr = 0.1, .r = 1.0, .g = 1.0, .b = 1.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = RFL, .gratt = RFL, .bratt = RFL, .rfatt = RFR, .gfatt = RFR, .bfatt = RFR, .rindex = RINDEX },
#elif 1
#define SR 0.05
#define R 1.0
#define G 1.0
#define B 1.0
#define Z 0
#define CX 0.05*WINSCALE
#define CY 0.08*WINSCALE
#define CZ 0.1*WINSCALE
//	{ .cx = 0*CX, .cy = 0*CY, .cz = -50*CZ, .sr = 5.0*WINSCALE, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#if 1 /* ioccc ray */
	{ .cx = -7*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -3*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 5*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 8*CX, .cy = 4*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -7*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -4*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 1*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 4*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = 3*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -7*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -4*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 1*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 4*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = 2*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -7*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -3*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = 1*G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 5*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = R, .g = G, .b = 1*B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 8*CX, .cy = 1*CY, .cz = Z*CZ, .sr = SR, .r = 1*R, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

#undef CX
#undef CY
#undef CZ
#undef R
#undef G
#undef B
#define R 0.0
#define G 0.0
#define B 0.0
#define CX 0.025*WINSCALE
#define CY 0.013*WINSCALE
#define CZ 0.05*WINSCALE
#undef Z
#define Z 5
	{ .cx = -11*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -9*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 11*CX, .cy = -5*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -11*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -7*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 7*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 11*CX, .cy = -10*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -11*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -9*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 9*CX, .cy = -15*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },

	{ .cx = -11*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -7*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = 1, .g = G, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = -2*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 2*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = 1, .b = B, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 9*CX, .cy = -20*CY, .cz = Z*CZ, .sr = 1*SR, .r = R, .g = G, .b = 1, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#endif
#else
	{ .cx = 0.0, .cy = -0.1, .cz = 0.0, .sr = 0.2, .r = 1.0, .g = 0.0, .b = 0.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0.0, .cy = 0.1, .cz = 0.0, .sr = 0.2, .r = 0.0, .g = 1.0, .b = 0.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
	{ .cx = 0.0, .cy = 0.4, .cz = 0.0, .sr = 0.2, .r = 0.0, .g = 0.0, .b = 1.0, .rdatt = 1.0, .gdatt = 1.0, .bdatt = 1.0, .rratt = 1.0, .gratt = 1.0, .bratt = 1.0 },
#endif
};
int __nsph = sizeof (_spheres) / sizeof (_spheres[0]);
sphere_t *__spheres = _spheres;

#define ATT_MIN 0.001
#define LEV_MAX -1

int sky_color( int level, double ex, double ey, double ez, double _vx, double _vy, double _vz, double *_r, double *_g, double *_b)
{
	double rmin = 0, gmin = 0, bmin = 0;
#ifdef USE_SKY
	double coef2;
	double th;
//	double ph;
	double n = 1;
	n = sqrt( _vx * _vx + _vy * _vy + _vz * _vz);
	th = asin( _vy / n) * 180 / M_PI;
//	th = acos( _vx / n) * 180 / M_PI;
//	th = atan( _vy / _vx) * 180 / M_PI;
#define THR 180
	double thr = THR;
#if 1
	if (level == 0)
	{
		thr /= 3;
		coef2 = (fmod( fabs( th - thr / 2) , thr)) / thr;
	}
	else
#endif
	coef2 = (fmod( fabs( th - thr / 2) , thr)) / thr;
	double a, b;
	double x1, x2, y1, y2;
	if (coef2 <= 0.333)
	{
		x1 = 0; x2 = 0.333;
		coef2 -= x1;

		y1 = 6; y2 = 255;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		rmin = a * coef2 + b;

		y1 = 105; y2 = 255;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		gmin = a * coef2 + b;

		y1 = 155; y2 = 255;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		bmin = a * coef2 + b;
	}
	else if (coef2 <= 0.666)
	{
		x1 = 0.333; x2 = 0.666;
		coef2 -= x1;

		y1 = 255; y2 = 218;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		rmin = a * coef2 + b;

		y1 = 255; y2 = 178;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		gmin = a * coef2 + b;

		y1 = 255; y2 = 127;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		bmin = a * coef2 + b;
	}
	else if (coef2 <= 1.0)
	{
		x1 = 0.666; x2 = 1.0;
		coef2 -= x1;

		y1 = 218; y2 = 103;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		rmin = a * coef2 + b;

		y1 = 178; y2 = 55;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		gmin = a * coef2 + b;

		y1 = 127; y2 = 26;
		a = (y2 - y1) / (x2 - x1);
		b = y1;
		bmin = a * coef2 + b;
	}
	else
	{
		printf( "boom sky color coef2 out of bound coef2=%f\n", coef2);fflush( stdout);
		getchar();
		exit( 3);
	}
	rmin /= 255;
	gmin /= 255;
	bmin /= 255;
				if ((rmin > 1.0) || (gmin > 1.0) || (bmin > 1.0) || (rmin < 0.0) || (gmin < 0.0) || (bmin < 0.0))
				{
					printf( "boom sky color overflow r=%f g=%f b=%f\n", rmin, gmin, bmin);fflush( stdout);
					getchar();
					exit( 3);
				}
#endif
	*_b = bmin;
	*_g = gmin;
	*_r = rmin;

	return 0;
}

unsigned long traced = 0;
unsigned long intersected = 0;
unsigned long reflected = 0;
unsigned long refracted = 0;
int level_max_reached = -1;
int level_max = LEV_MAX;
int traceray( int level, double ex, double ey, double ez, double _vx, double _vy, double _vz, double *r, double *g, double *b, char *pix, double ratt, double gatt, double batt)
{
	v3 e, _v;
	set3( e, ex, ey, ez);
	set3( _v, _vx, _vy, _vz);
	int nsph = scene.nspheres;
	sphere_t *spheres = scene.spheres;
	int nlamps = scene.nlamps;
	lamp_t *lamps = scene.lamps;
	if (level == 0)
		traced++;
	if (level > level_max_reached)
		level_max_reached = level;
	unsigned long smin = -1, s;
	double tmin = BIG, tmax = 0;
//#define CHEAPO_LIGHTING
#ifdef CHEAPO_LIGHTING
	double srmin = 0;
#endif

	double rmin = 0.0, gmin = 0.0, bmin = 0.0;		// material
	double rdatt = 0.0, gdatt = 0.0, bdatt = 0.0;	// diffuse
	double rratt = 0.0, gratt = 0.0, bratt = 0.0;	// reflected
	double rfatt = 0.0, gfatt = 0.0, bfatt = 0.0;	// refracted
	double rindex = 0.0;
	
	for (s = 0; s < nsph; s++)
	{
//		double cx, cy, cz, sr;
//		cx = spheres[s].cx; cy = spheres[s].cy; cz = spheres[s].cz; sr = spheres[s].sr;
		v3 c;
		double sr;
		c[0] = spheres[s].cx; c[1] = spheres[s].cy; c[2] = spheres[s].cz; sr = spheres[s].sr;
//		printf( "sphere: c(%f;%f;%f) r(%f)\n", cx, cy, cz, sr);
		int res;
		double t = 0, t2 = 0;
//		res = intersec_sphere( cx, cy, cz, sr, ex, ey, ez, _vx, _vy, _vz, &t, &t2);
		res = intersec_sphere( c, sr, e, _v, &t, &t2);
//		printf( "%c", res ? 's' : '.');
		if (res)
		{
			if (t < tmin)
			{
				smin = s;
				tmin = t;
				tmax = t2;
#ifdef CHEAPO_LIGHTING
				srmin = spheres[s].sr;
#endif

				rmin = spheres[s].r;
				gmin = spheres[s].g;
				bmin = spheres[s].b;
				rdatt = spheres[s].rdatt;
				gdatt = spheres[s].gdatt;
				bdatt = spheres[s].bdatt;
				rratt = spheres[s].rratt;
				gratt = spheres[s].gratt;
				bratt = spheres[s].bratt;
				rfatt = spheres[s].rfatt;
				gfatt = spheres[s].gfatt;
				bfatt = spheres[s].bfatt;
				rindex = spheres[s].rindex;
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
		sky_color( level, ex, ey, ez, _vx, _vy, _vz, &rmin, &gmin, &bmin);
#if 1
		if (level > 0)
			dprintf( "sky color %f %f %f\n", rmin, gmin, bmin);
#endif
	}
	else				// object
	{
		double rdiff = 0.0, gdiff = 0.0, bdiff = 0.0;
		double rrefl = 0.0, grefl = 0.0, brefl = 0.0;
		double rrefr = 0.0, grefr = 0.0, brefr = 0.0;
	
		if (level > 0)
			intersected++;
//		printf( "object\n");
// compute intersections here
		double rx, ry, rz, rx2, ry2, rz2;		// coord of intersec
		rx = ex + _vx * tmin;
		ry = ey + _vy * tmin;
		rz = ez + _vz * tmin;
		rx2 = ex + _vx * tmax;
		ry2 = ey + _vy * tmax;
		rz2 = ez + _vz * tmax;
		double n;
		dprintf( "ray is (%f,%f,%f)\n", _vx, _vy, _vz);
		dprintf( "inters with sphere %lu at (%f,%f,%f)\n", smin, rx, ry, rz);
		double nvx, nvy, nvz;	// normal vect

		double grad = 1;
#ifdef CHEAPO_LIGHTING
		// this is a cheapo way to make ambient light/shadows
// compute distance between two intersec points
		nvx = rx - rx2;
		nvy = ry - ry2;
		nvz = rz - rz2;
		n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
		grad = n / srmin / 2;
		grad = sqrt( grad);
		rdiff = rmin * rdatt * grad;
		gdiff = gmin * gdatt * grad;
		bdiff = bmin * bdatt * grad;
#else
// compute distance between lamp and intersec point
		int i;
		for (i = 0; i < nlamps; i++)
		{
			rx2 = lamps[i][0];
			ry2 = lamps[i][1];
			rz2 = lamps[i][2];
			nvx = rx2 - rx;
			nvy = ry2 - ry;
			nvz = rz2 - rz;
			n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
#define LIGHT_CONST 0.001
			n += LIGHT_CONST;
			grad = 1 / n / n;
			break;
		}
		if (grad > 1)
		{
			rdiff = 1;
			gdiff = 1;
			bdiff = 1;
		}
		else
		{
		rdiff = rmin * rdatt * grad;
		gdiff = gmin * gdatt * grad;
		bdiff = bmin * bdatt * grad;
		}
#endif

		double rvx, rvy, rvz;	// vect of reflected/refracted ray

		if (1
			&& (((ratt * rratt) >= ATT_MIN) || ((gatt * gratt) >= ATT_MIN) || ((batt * bratt) >= ATT_MIN))
			&& ((level < level_max) || (level_max == -1))
		)
		{
// loose energy
			double loss = 0.9;
			ratt *= loss;
			gatt *= loss;
			batt *= loss;
		
			if (level > 0)
				reflected++;
// normal at intersec
		nvx = rx - spheres[smin].cx;
		nvy = ry - spheres[smin].cy;
		nvz = rz - spheres[smin].cz;
		n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
		nvx /= n;
		nvy /= n;
		nvz /= n;

		double dot = 1.0;
		dot = (_vx * nvx + _vy * nvy + _vz * nvz);
		rvx = _vx - 2 * dot * nvx;
		rvy = _vy - 2 * dot * nvy;
		rvz = _vz - 2 * dot * nvz;
		n = sqrt( rvx * rvx + rvy * rvy + rvz * rvz);
		rvx /= n;
		rvy /= n;
		rvz /= n;
		dprintf( "refl vec (%f,%f,%f) (s=%lu)\n", rvx, rvy, rvz, smin);
			traceray( level + 1, rx, ry, rz, rvx, rvy, rvz, &rrefl, &grefl, &brefl, pix, ratt * rratt, gatt * gratt, batt * bratt);
		}

		if (1
			&& (((ratt * rfatt) >= ATT_MIN) || ((gatt * gfatt) >= ATT_MIN) || ((batt * bfatt) >= ATT_MIN))
			&& ((level < level_max) || (level_max == -1))
		)
		{
// loose energy
			double loss = 0.9;
			ratt *= loss;
			gatt *= loss;
			batt *= loss;
		
			if (level > 0)
				refracted++;
#if 1
// normal at intersec
		nvx = rx - spheres[smin].cx;
		nvy = ry - spheres[smin].cy;
		nvz = rz - spheres[smin].cz;
		n = sqrt( nvx * nvx + nvy * nvy + nvz * nvz);
		nvx /= n;
		nvy /= n;
		nvz /= n;

		double dot = 1.0;
		dot = (-_vx * nvx + -_vy * nvy + -_vz * nvz);
		rvx = -dot * nvx - _vx;
		rvy = -dot * nvy - _vy;
		rvz = -dot * nvz - _vz;
		rvx *= 1.0 - rindex;
		rvy *= 1.0 - rindex;
		rvz *= 1.0 - rindex;
		rvx += _vx;
		rvy += _vy;
		rvz += _vz;
		n = sqrt( rvx * rvx + rvy * rvy + rvz * rvz);
		rvx /= n;
		rvy /= n;
		rvz /= n;
#else
		rvx = _vx;
		rvy = _vy;
		rvz = _vz;
#endif
		dprintf( "refr vec (%f,%f,%f) (s=%lu)\n", rvx, rvy, rvz, smin);
			traceray( level + 1, rx, ry, rz, rvx, rvy, rvz, &rrefr, &grefr, &brefr, pix, ratt * rfatt, gatt * gfatt, batt * bfatt);
		}

		if ((rdatt + rratt + rfatt + gdatt + gratt + gfatt + bdatt + bratt + bfatt) > SMALL)
		{
			rmin = (rdiff * rdatt + rrefl * rratt + rrefr * rfatt) / (rdatt + rratt + rfatt);
			gmin = (gdiff * gdatt + grefl * gratt + grefr * gfatt) / (gdatt + gratt + gfatt);
			bmin = (bdiff * bdatt + brefl * bratt + brefr * bfatt) / (bdatt + bratt + bfatt);
		}

		if ((rmin > 1.0) || (gmin > 1.0) || (bmin > 1.0) || (rmin < 0.0) || (gmin < 0.0) || (bmin < 0.0))
		{
			printf( "boom object color overflow r=%f g=%f b=%f\n", rmin, gmin, bmin);fflush( stdout);
			getchar();
			exit( 3);
		}
	}
	*r = rmin;
	*g = gmin;
	*b = bmin;
	return 0;
}

int stats()
{
	double cur_t;
#ifndef USETIMEOFDAY
	cur_t = time( 0);
#else
	struct timeval tv;
	gettimeofday( &tv, 0);
	cur_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
	static double old_t = SMALL / 2;
	if (old_t <= SMALL)
		old_t = cur_t;
	printf( " trcd=%lu intr=%lu ", traced, intersected);
	printf( "refl=%lu refr=%lu MAX=%d max=%d ", reflected, refracted, level_max, level_max_reached);
	double duration = cur_t - old_t;
	unsigned long last_traced = 0;
	unsigned long perf = (traced - last_traced) / duration;
	if (duration > 0)
	{
		printf( "- %lu ray/s dur=%.2fs", perf, duration);
	}
	printf( "\n");
//	old_t = cur_t;
	last_traced = traced;
	return 0;
}

#define USE_CAMERA

int main( int argc, char *argv[])
{
#ifdef USE_SDL
	SDL_Surface *screen = 0;
	int do_sdl = 1;
#endif
	double ex, ey, ez, vx, vy, vz;
	char *scene_file = 0;	
	
	ex = 0*WINSCALE; ey = 0*WINSCALE; ez = 1*WINSCALE; vx = 0.0; vy = 0.0; vz = -1*WINSCALE;
	unsigned long w, h;
	w = W; h = H;
	double winw, winh;
	double winscale = WINSCALE;
	int do_tga = 1;
	int do_txt = 0;
	int bpp = 24;
	int tga_fd = -1;
	unsigned char *tga_map = MAP_FAILED;
	size_t tga_size = 0;
	size_t tga_index = 0;

	int arg = 1;
	while (arg < argc)
	{
		if (!strcmp( argv[arg], "-w"))
		{
			if (++arg >= argc)
			{
				printf( "missing width argument\n");
				exit( 1);
			}
			sscanf( argv[arg++], "%lu", &w);
		}
		else if (!strcmp( argv[arg], "-h"))
		{
			if (++arg >= argc)
			{
				printf( "missing height argument\n");
				exit( 1);
			}
			sscanf( argv[arg++], "%lu", &h);
		}
		else if (!strcmp( argv[arg], "-d"))
		{
			if (++arg >= argc)
			{
				printf( "missing depth argument\n");
				exit( 1);
			}
			sscanf( argv[arg++], "%d", &level_max);
		}
		else
		{
			if (scene_file)
			{
				printf( "scene file can't be set twice (already set as '%s')\n", scene_file);
				exit( 1);
			}
			scene_file = argv[arg++];
		}
	}
	scene.spheres = __spheres;
	scene.nspheres = __nsph;
	if (scene_file)
	{
		load_scene( &scene, scene_file);
	}
	printf( "nspheres=%d\n", scene.nspheres);
	printf( "nlamps=%d\n", scene.nlamps);
#ifdef USE_SDL
	if (do_sdl)
	{
		SDL_Init( SDL_INIT_VIDEO);
		screen = SDL_SetVideoMode( w, h, bpp, 0);
		if (screen)
			atexit( SDL_Quit);
	}
#endif
	winw = 1.0 * winscale; winh = winw * h / w;
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
	double old_t;
	double last_t;
	double cur_t;
#ifndef USETIMEOFDAY
	old_t = time( 0);
#else
	struct timeval tv;
	gettimeofday( &tv, 0);
	old_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
	last_t = old_t;
	for (j = 0; j < h; j++)
	{
		double _vx, _vy, _vz;
//		_vy = vy + (((double)h - 1) - j) / (h - 1) - winh / 2;
		_vy = vy - winh / 2 + winh * (double)(h - j - 1) / (h - 1);
		_vz = vz;

		for (i = 0; ; i++)
		{
#ifndef USETIMEOFDAY
			cur_t = time( 0);
#else
			gettimeofday( &tv, 0);
			cur_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
			int percent = (double)100.0 * (j * w + i) / (w * h);
			static int old_percent = -1;
			if ((percent != old_percent) || (cur_t >= (last_t + 1)))
			{
				if (cur_t >= (last_t + 1))
					last_t = cur_t;
				old_percent = percent;
				if (!(old_percent % 10))
				{
					printf( " %d%%", old_percent);
					stats();
					fflush( stdout);
				}
#ifdef USE_SDL
				if (screen)
				{
					SDL_UpdateRect( screen, 0, 0, 0, 0);
				}
#endif
			}
			if (i >= w)
				break;
			char pix = '.';
#ifndef USE_CAMERA
			_vx = vx - winw / 2 + winw * (double)i / (w - 1);
#else
			v3 v;
			v[0] = vx;
			v[1] = vy;
			v[2] = vz;
			v3 e;
			e[0] = ex;
			e[1] = ey;
			e[2] = ez;
			// p=p0+(p1-p0)u+(p2-p0)v
			v3 up = { 0, 1, 0};
			v3 right;
			cross3( right, v, up);
			cross3( up, right, v);
			v3 s0, s1, s2, s;
			sum3( s0, e, v);
			sum3( s1, s0, up);
			sum3( s2, s0, right);
			double u, vv;
			vv = -(double)0.5 + (double)1.0 * (double)i / ((double)w - 1);
			u = -(double)0.5 + (double)1.0 * (double)(h - j - 1) / ((double)h - 1);
			s[0] = s0[0] + (s1[0] - s0[0]) * u + (s2[0] - s0[0]) * vv;
			s[1] = s0[1] + (s1[1] - s0[1]) * u + (s2[1] - s0[1]) * vv;
			s[2] = s0[2] + (s1[2] - s0[2]) * u + (s2[2] - s0[2]) * vv;
			v3 _v;
			diff3( _v, s, e);
			_vx = _v[0];
			_vy = _v[1];
			_vz = _v[2];
#endif

//			printf( "j=%lu i=%lu vy=%f vx=%f\n", j, i, _vy, _vx);
			double r = 0.0, g = 0.0, b = 0.0;

			traceray( 0, ex, ey, ez, _vx, _vy, _vz, &r, &g, &b, &pix, 1.0, 1.0, 1.0);
			unsigned char cr, cg, cb;
			cr = (double)255 * r;
			cg = (double)255 * g;
			cb = (double)255 * b;
#define LO 3
#define HI 253
			if ((cr <= LO) && (cg <= LO) && (cb <= LO))
			{
				static int count = 0;
				if (!count)
				{
					printf( "got black ! x=%lu y=%lu\n", i, j);
					count++;
				}
			}
			if ((cr >= HI) && (cg >= HI) && (cb >= HI))
			{
				static int count = 0;
				if (!count)
				{
					printf( "got white ! x=%lu y=%lu r=%u g=%u b=%u\n", i, j, cr, cg, cb);
					count++;
				}
			}
			if ((cr >= HI) && (cg <= LO) && (cb <= LO))
			{
				static int count = 0;
				if (!count)
				{
					printf( "got red ! x=%lu y=%lu r=%u g=%u b=%u\n", i, j, cr, cg, cb);
					count++;
				}
			}
			if ((cr <= LO) && (cg >= HI) && (cb <= LO))
			{
				static int count = 0;
				if (!count)
				{
					printf( "got green ! x=%lu y=%lu r=%u g=%u b=%u\n", i, j, cr, cg, cb);
					count++;
				}
			}
			if ((cr <= LO) && (cg <= LO) && (cb >= HI))
			{
				static int count = 0;
				if (!count)
				{
					printf( "got blue ! x=%lu y=%lu r=%u g=%u b=%u\n", i, j, cr, cg, cb);
					count++;
				}
			}
#ifdef USE_SDL
			if (screen)
			{
				Uint32 col;
				SDL_Rect rect;
				rect.x = i;
				rect.y = j;
				rect.w = rect.h = 1;
				col = SDL_MapRGB( screen->format, cr, cg, cb);
				SDL_FillRect( screen, &rect, col);
			}
#endif
			//dprintf( "  (r=%f g=%f b=%f)", r, g, b);
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
			//	dprintf( "{%02x:%02x:%02x}", cr, cg, cb);
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
		//dprintf( "\n");
	}
#ifndef USETIMEOFDAY
	cur_t = time( 0);
#else
	gettimeofday( &tv, 0);
	cur_t = tv.tv_sec + tv.tv_usec / 1000000.0;
#endif
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
//	dprintf( "tga_size=%lu\n", (unsigned long)tga_size);
//	dprintf( "tga_index=%lu\n", (unsigned long)tga_index);
	printf( "traced=%lu intersected=%lu\n", traced, intersected);
	printf( "reflected=%lu refracted=%lu level_max=%d level_max_reached=%d\n", reflected, refracted, level_max, level_max_reached);
	double duration = cur_t - old_t;
	if (!duration)
		duration = 1;
	unsigned long perf = (w * h) / duration;
	printf( "perf : %lu ray/s duration=%.2fs (old_t=%.2f cur_t=%.2f)\n", perf, duration, old_t, cur_t);
#ifdef USE_SDL
	if (screen)
	{
		SDL_Event event;
		int end = 0;
		while (!end)
		{
			SDL_WaitEvent( &event);
			switch (event.type)
			{
				case SDL_QUIT:
					end = 1;
					break;
				case SDL_KEYDOWN:
					end = 1;
					break;
				case SDL_MOUSEBUTTONDOWN:
					end = 1;
					break;
				default:
					break;
			}
		}
	}
#endif
	
	return 0;
}
