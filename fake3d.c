#include <stdio.h>

#include <SDL.h>
#include <SDL_thread.h>
#include <SDL_mutex.h>

//#define DEBUG

#include "vec.h"

#ifdef DEBUG
#define dprintf printf
#else
#define dprintf(...) do{}while(0)
#endif

typedef struct {
// mutable, protected by sem_init
	int num;
	SDL_sem *sem_init;

// non-mutable, common, or protected data (workers_done, start..)
	int tot;
	void *payload;		// opaque user ptr to store useful common payload data
	SDL_mutex *mutex;
	SDL_cond *cond;
	int *workers_done;

	SDL_mutex *mutex_go;
	SDL_cond *cond_go;
	int *start;			// incrementing
	int *quit;
} worker_t;

typedef struct {
	//	objs
	//	{
	//		{
	//			s type		sphere=0; facet=1; box=2
	//			s s0		sphere: radius
	//			s N/A
	//		}
	//		v3 loc0		sphere: center;	facet: vertex0; box: lower
	//		v3 loc1		facet: vertex1; box: upper
	//		v3 loc2		facet: vertex2
	//		v3 color
	//		}
	//	}
#define LEN_OBJ 5
#define OBJ_LOC0 1
#define OBJ_LOC1 2
#define OBJ_LOC2 3
#define OBJ_COLOR 4
	v3 *objs;
	int nobjs;

	v3 *facets;
	int nfacets;
	v3 *spheres;
	int nspheres;
	v3 *boxes;
	int nboxes;

// camera is : eye coordinate (vector e)..
// ..a direction (vector v)..
// ..and an "up" (vector up) (camera "head" rotation, default pointing to the "sky")
	double *e, *v, *up;
// camera screen size
	int w;
	int h;
	// camera screen as three vectors of a plane : s0, s1 and s2
	v3 s0, s1, s2;
} scene_t;

typedef struct {
	SDL_Surface *screen;
	int w;
	int h;
	scene_t *scene;
} payload_t;

#define EPS 0.1
#define BIG 1000.0
#define COMP_EPS(x,val) (x >= ((val) - EPS) && x <= ((val) + EPS))

int intersec_obj( v3 *obj, v3 e, v3 v, double *tmin)
{
	int res = 0;
	if (obj[0][0] == 0) // sphere
		res = intersec_sphere( obj[OBJ_LOC0], obj[0][1], e, v, tmin, 0);
	else if (obj[0][0] == 1) // facet
		res = intersec_plane( obj[OBJ_LOC0], obj[OBJ_LOC1], obj[OBJ_LOC2], e, v, tmin);
	else if (obj[0][0] == 2) // box
		res = intersec_box( obj[OBJ_LOC0], obj[OBJ_LOC1], e, v, tmin, 0);
	return res;
}

void traceray( scene_t *scene, v3 _e, v3 _v, v3 col)
{
	double *pcol = 0;
	double tmin = BIG;
	double t;
	int res;
	int i;
	for (i = 0; i < scene->nobjs; i++)
	{
		t = 0;
		res = intersec_obj( &scene->objs[i * LEN_OBJ], _e, _v, &t);
//		dprintf( "result=%d t=%f\n", res, t);
		if (res && (t < tmin))
		{
			tmin = t;
			pcol = scene->objs[i * LEN_OBJ + OBJ_COLOR];
		}
	}
	if (pcol)
		memcpy( col, pcol, 3 * sizeof( col[0]));
}

// input : s0,s1,s2,e
inline int compute_scene( scene_t *scene, int i, int j, int *r, int *g, int *b)
{
	int w = scene->w;
	int h = scene->h;
	double u, v;
	v = -(double)0.5 + (double)1.0 * (double)i / ((double)w - 1);
	u = -(double)0.5 + (double)1.0 * (double)(h - j - 1) / ((double)h - 1);
	// compute a point s in camera screen plane based on {u,v}
	v3 s;
	// p=p0+(p1-p0)u+(p2-p0)v
	s[0] = scene->s0[0] + (scene->s1[0] - scene->s0[0]) * u + (scene->s2[0] - scene->s0[0]) * v;
	s[1] = scene->s0[1] + (scene->s1[1] - scene->s0[1]) * u + (scene->s2[1] - scene->s0[1]) * v;
	s[2] = scene->s0[2] + (scene->s1[2] - scene->s0[2]) * u + (scene->s2[2] - scene->s0[2]) * v;
	// compute final camera=> pixel vector _v
	v3 _v;
	diff3( _v, s, scene->e);
	v3 col;
	memset( col, 0, sizeof( col));
	traceray( scene, scene->e, _v, col);

	if ((col[0] < 0) || (col[0] > 1) || (col[1] < 0) || (col[1] > 1) || (col[2] < 0) || (col[2] > 1))
	{
		printf( "boom\n");
		exit( 0);
	}
	*r = col[0] * 255;
	*g = col[1] * 255;
	*b = col[2] * 255;

	return 0;
}

// input : v,up,e
// output : s0,s1,s2,up,v
inline int init_scene( scene_t *scene)
{
	// compute camera screen as three vectors of a plane : s0, s1 and s2
	v3 right;
	// normalize v
	div3( scene->v, norm3( scene->v));
	cross3( right, scene->v, scene->up);
	cross3( scene->up, right, scene->v);
	sum3( scene->s0, scene->e, scene->v);
	sum3( scene->s1, scene->s0, scene->up);
	sum3( scene->s2, scene->s0, right);
	return 0;
}

inline int load_scene( scene_t *scene, char *file)
{
	int i, j;
	FILE *in = fopen( file, "rt");
	if (in)
	{
		scene->nobjs = 0;
		scene->nfacets = 0;
		scene->nspheres = 0;
		scene->nboxes = 0;
		char buf[1024];
		char *ptr;
		while (!feof( in))
		{
			if (!fgets( buf, sizeof( buf), in))
				break;
			ptr = strstr( buf, "endfacet");
			if (ptr)
				scene->nfacets++;
		}
		printf( "read nfacets=%d\n", scene->nfacets);
		rewind( in);
		if (scene->nfacets)
		{
		scene->facets = malloc( scene->nfacets * sizeof(v3[4]));
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
////
	return 0;
}

int thr( void *opaque)
{
	worker_t *w = opaque;
	payload_t *p = 0;
	int num = -1;

	if (!w)
		return -1;
	num = w->num;
	SDL_SemPost( w->sem_init);		// signal master we have done with reading worker data

	p = w->payload;
	int x1, y1, x2 = 0, y2 = 0;
	int x = 0;
	int y = 0;
	int running = 0;
	int start = *w->start;
	dprintf( "thr %d ready to work\n", num);	
	while (1)
	{
		// now wait for master to signal go
		int end = 0;
		while (!end && !running && !*w->quit)
		{
			SDL_LockMutex( w->mutex_go);
			while (((*w->start == start) && !running && !*w->quit) && SDL_CondWait( w->cond_go, w->mutex_go) == 0)				// this blocks without hogging cpu
				continue;
			if ((*w->start != start) || running || *w->quit)
				end = 1;
			if (*w->start != start)
			{
				running = 0;
				start = *w->start;
			}
			SDL_UnlockMutex( w->mutex_go);
		}
		if (*w->quit)
		{
			dprintf( "thr %d will quit\n", num);	
			break;
		}
		if (!running)
		{
			dprintf( "thr %d starting %d\n", num, start);
			x1 = 0;
			x2 = p->w - 1;
			y1 = num * p->h / w->tot;
			y2 = (num + 1) * p->h / w->tot - 1;
			x = x1;
			y = y1;
			// start working
			running = 1;
		}

		// working
		Uint32 col;
		SDL_Rect rect;
		rect.x = x;
		rect.y = y;
		rect.w = 1;
		rect.h = 1;
		int r, g, b;
		compute_scene( p->scene, x, y, &r, &g, &b);
		col = SDL_MapRGB( p->screen->format, r, g, b);
		SDL_FillRect( p->screen, &rect, col);

//		SDL_Delay( 1);

		if (x < x2)
			x++;
		else
		{
			x = x1;
			if (y < y2)
				y++;
			else
			{
				running = 0;
				// now signal master we're done
				SDL_LockMutex( w->mutex);
				dprintf( "thr %d done %d\n", num, start);
				(*w->workers_done)++;
				SDL_UnlockMutex( w->mutex);
				SDL_CondSignal( w->cond);
			}
		}
	}
	dprintf( "thr %d quitting\n", num);
	return 0;
}

// 3d scene :
#define F 40
#define S 20
// multiple objects
v3 __objs[] = {

	{ 1, 0, 0 },	// 1=facet
	{ 0, 0, 0 },
	{ F, 0, 0 },
	{ 0, F, 0 },
	{ 1, 0, 0 },	// color

	{ 1, 0, 0 },	// 1=facet
	{ 0, 0, 0 },
	{ 0, F, 0 },
	{ 0, 0, F },
	{ 0, 1, 0 },	// color

	{ 1, 0, 0 },	// 1=facet
	{ 0, 0, 0 },
	{ 0, 0, F },	// WARNING : visible faces must have normal pointing to eye !
	{ F, 0, 0 },
	{ 0, 0, 1 },	// color

	{ 2, 0, 0 },	// 2=box
	{ -S, -S, -S },	// lower
	{ S, S, S },	// upper
	{},
	{ 1, 1, 0 },	// color
};
int _nobjs = sizeof( __objs) / sizeof( __objs[0]) / LEN_OBJ;
// multiple facets (3 vertex3 + 1 color3 each)
v3 __facets[] = {
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
int _nfacets = sizeof( __facets) / sizeof( __facets[0]) / 4;
// multiple spheres
v3 __spheres[] = {
#if 0
	{ S, S, S },
	{ S/4, 0, 0 },	// radius
	{ 1, 0, 1 },	// color
#endif
};
int _nspheres = sizeof( __spheres) / sizeof( __spheres[0]) / 3;

v3 __boxes[] = {
#if 0
	{ -S, -S, -S },	// lower
	{ S, S, S },	// upper
	{ 1, 1, 0 },	// color
#endif
};
int _nboxes = sizeof( __boxes) / sizeof( __boxes[0]) / 3;

v3 *_objs = __objs;
v3 *_facets = __facets;
v3 *_spheres = __spheres;
v3 *_boxes = __boxes;

// camera is : eye coordinate (vector e)..
#if 1
#define E 20
	v3 _e = { 2*E, 2*E, 2*E };
// ..a direction (vector v)..
#define V E/2
	v3 _v = { -2*V, -2*V, -2*V };
#else
#define E 10
	v3 _e = { 0*E, 0*E, 1*E };
// ..a direction (vector v)..
#define V E/2
	v3 _v = { -0*V, -0*V, -1*V };
#endif
// ..and an "up" (vector up) (camera "head" rotation, default pointing to the "sky")
v3 _up = { 0, 1, 0};

#define W 800
#define H 600
int main( int argc, char *argv[])
{
	char *scene_file = 0;
	SDL_Surface *screen = 0;
	int w = W;
	int h = H;
	int bpp = 32;
	int n = 1;
	int arg = 1;
	if (arg < argc)
	{
		scene_file = argv[arg++];
		if (arg < argc)
			sscanf( argv[arg++], "%d", &n);
	}
	int i;
	SDL_mutex *mutex;
	SDL_cond *cond;
	SDL_mutex *mutex_go;
	SDL_cond *cond_go;
	SDL_sem *sem_init;
	int workers_done = 0;
	int start = 0;
	int quit = 0;

	printf( "hello fake3d\n");

	SDL_Init( SDL_INIT_VIDEO);
	atexit( SDL_Quit);
	screen = SDL_SetVideoMode( w, h, bpp, 0);
	SDL_EnableKeyRepeat( SDL_DEFAULT_REPEAT_DELAY / 2, SDL_DEFAULT_REPEAT_INTERVAL);
	mutex = SDL_CreateMutex();
	cond = SDL_CreateCond();
	mutex_go = SDL_CreateMutex();
	cond_go = SDL_CreateCond();
	sem_init = SDL_CreateSemaphore( 0);
	worker_t wo;
	payload_t p;
	scene_t scene;
	scene.objs = _objs;
	scene.nobjs = _nobjs;
	scene.facets = _facets;
	scene.nfacets = _nfacets;
	scene.spheres = _spheres;
	scene.nspheres = _nspheres;
	scene.boxes = _boxes;
	scene.nboxes = _nboxes;
	scene.w = w;
	scene.h = h;
	scene.e = _e;
	scene.v = _v;
	scene.up = _up;
	p.scene = &scene;
	if (scene_file)	// load scene ?
		load_scene( &scene, scene_file);
	p.screen = screen;
	p.w = w;
	p.h = h;
	wo.tot = n;
	wo.mutex = mutex;
	wo.cond = cond;
	wo.workers_done = &workers_done;
	wo.start = &start;
	wo.quit = &quit;
	wo.mutex_go = mutex_go;
	wo.cond_go = cond_go;
	wo.payload = &p;
	wo.sem_init = sem_init;
	SDL_Thread **thrs = malloc( n * sizeof( SDL_Thread*));
	for (i = 0; i < n; i++)
	{
		wo.num = i;
		thrs[i] = SDL_CreateThread( thr, &wo);
		SDL_SemWait( sem_init);
	}
	SDL_DestroySemaphore( sem_init);
	int end = 0;
	while (!end)
	{
		static int go = 1;
		static int done = 1;
		int just_done = 0;
		SDL_Event event;
		while (SDL_PollEvent( &event))
		{
			switch (event.type)
			{
				case SDL_QUIT:
					end = 1;
					break;
				case SDL_KEYUP:
					switch (event.key.keysym.sym)
					{
						case SDLK_ESCAPE:
							end = 1;
							break;
						case SDLK_UP:
							go = 1;
							add3( scene.e, -0.5);
							break;
						case SDLK_DOWN:
							go = 1;
							add3( scene.e, +0.5);
							break;
						default:
							break;
					}
					break;
				case SDL_MOUSEBUTTONUP:
	//				if (done)
						go = 1;
					break;
				default:
					break;
			}
		}
		if (end)
			break;

		static Uint32 old_t = 0;
		Uint32 t;
		t = SDL_GetTicks();
#define DELAY_REFRESH	500
		if (t >= (old_t + DELAY_REFRESH))
		{
			SDL_UpdateRect( screen, 0, 0, 0, 0);
			old_t = t;
		}
		
		if (!done)
		{
			// see if workers are done
			SDL_LockMutex( mutex);
#define DELAY_WAIT_WORKERS 100
			while ((workers_done != n) && SDL_CondWaitTimeout( cond, mutex, DELAY_WAIT_WORKERS) == 0)
				continue;
			if (workers_done == n)
			{
				done = 1;
				just_done = 1;
			}
			SDL_UnlockMutex( mutex);
		}
		else
			SDL_Delay( DELAY_WAIT_WORKERS);

		int r, g, b;
		Uint32 col;
		SDL_Rect rect;
		int _w = 10;
		int _h = 10;

		if (just_done)
		{
			rect.x = w - _w;
			rect.y = h - _h;
			rect.w = _w;
			rect.h = _h;
			g = b = 0;
			r = 255;
			col = SDL_MapRGB( screen->format, r, g, b);
			SDL_FillRect( screen, &rect, col);
			SDL_UpdateRect( screen, 0, 0, 0, 0);
		}
		
//		if (done)
		{
			if (go)
			{
#if 0
				rect.x = 0;
				rect.y = 0;
				rect.w = w;
				rect.h = h;
				r = g = b = 128;
				col = SDL_MapRGB( screen->format, r, g, b);
				SDL_FillRect( screen, &rect, col);
#endif

				rect.x = w - _w;
				rect.y = h - _h;
				rect.w = _w;
				rect.h = _h;
				r = b = 0;
				g = 255;
				col = SDL_MapRGB( screen->format, r, g, b);
				SDL_FillRect( screen, &rect, col);

				SDL_UpdateRect( screen, 0, 0, 0, 0);

				init_scene( &scene);
				// start workers
				go = 0;
				done = 0;
				SDL_LockMutex( mutex_go);
				start++;
				dprintf( "master say start %d\n", start);
				workers_done = 0;
				SDL_CondBroadcast( cond_go);
				SDL_UnlockMutex( mutex_go);
			}
		}
	}
	SDL_LockMutex( mutex_go);
	dprintf( "master say quit\n");
	quit = 1;
	SDL_CondBroadcast( cond_go);
	SDL_UnlockMutex( mutex_go);
	for (i = 0; i < n; i++)
	{
		SDL_WaitThread( thrs[i], 0);
	}
	free( thrs);
	SDL_DestroyCond( cond);
	SDL_DestroyMutex( mutex);
	SDL_DestroyCond( cond_go);
	SDL_DestroyMutex( mutex_go);
	return 0;
}
