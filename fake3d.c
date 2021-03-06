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
//		{
//			s type		sphere=0; facet=1; box=2; paral=3 ; cyl=4
//			s s0		sphere: radius; cyl: radius
//			s c0		cyl: height
//		}
//		v3 loc0		sphere: center;	facet: vertex0; box: lower; cyl: axis
//		v3 loc1		facet: vertex1; box: upper
//		v3 loc2		facet: vertex2
//		v3 color
#define LEN_OBJ 5
#define OBJ_LOC0 1
#define OBJ_LOC1 2
#define OBJ_LOC2 3
#define OBJ_COLOR 4
	v3 *objs;
	int nobjs;

	v3 *lamps;
	int nlamps;

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

int intersec_obj( v3 *obj, v3 e, v3 v, double *tmin, v3 normal, v3 col)
{
	int res = 0;

	mult3( normal, 0);	// default no normal
	if (obj[0][0] == 0) // sphere
	{
		res = intersec_sphere( obj[OBJ_LOC0], obj[0][1], e, v, tmin, 0);
		if (res)
		{
			v3 pinter; // intersection pos
			copy3( pinter, v);
			mult3( pinter, *tmin);
			sum3( pinter, e, pinter);
			diff3( normal, pinter, obj[OBJ_LOC0]);
			div3( normal, norm3( normal));
		}
	}
	else if (obj[0][0] == 1) // facet
	{
		res = intersec_plane( obj[OBJ_LOC0], obj[OBJ_LOC1], obj[OBJ_LOC2], e, v, tmin);
		if (res)
		{
			v3 n1, n2;
			diff3( n1, obj[OBJ_LOC0], obj[OBJ_LOC1]);
			diff3( n2, obj[OBJ_LOC0], obj[OBJ_LOC2]);
			cross3( normal, n1, n2);
		}
	}
	else if (obj[0][0] == 3) // parallelogram
	{
		res = intersec_parallelogram( obj[OBJ_LOC0], obj[OBJ_LOC1], obj[OBJ_LOC2], e, v, tmin);
		if (res)
		{
			v3 n1, n2;
			diff3( n1, obj[OBJ_LOC0], obj[OBJ_LOC1]);
			diff3( n2, obj[OBJ_LOC0], obj[OBJ_LOC2]);
			cross3( normal, n1, n2);
		}
	}
	else if (obj[0][0] == 2) // box
	{
		int num = 0;
		res = intersec_box( obj[OBJ_LOC0], obj[OBJ_LOC1], e, v, tmin, 0, &num);
		if (res)
		{
#if 0
//		asm volatile( "int $3");
		if (num == 1)
		{
			col[0] = 1;
			col[1] = 0;
			col[2] = 0;
		}
		else if (num == 2)
		{
			col[0] = 1;
			col[1] = 1;
			col[2] = 0;
		}
		else if (num == 4)
		{
			col[0] = 1;
			col[1] = 1;
			col[2] = 1;
		}
		else if (num == 8)
		{
			col[0] = 1;
			col[1] = 1;
			col[2] = 1;
		}
		else if (num == 16)
		{
			col[0] = 0;
			col[1] = 1;
			col[2] = 1;
		}
		else if (num == 32)
		{
			col[0] = 0;
			col[1] = 0;
			col[2] = 1;
		}
#endif
			if (num == 1)
			{
				normal[0] = 0.0;
				normal[1] = 0.0;
				normal[2] = 0.0;
			}
			else if (num == 2)
			{
				normal[0] = 0.0;
				normal[1] = 0.9;
				normal[2] = 0.0;
			}
			else if (num == 4)
			{
				normal[0] = 0.0;
				normal[1] = 0.0;
				normal[2] = 0.0;
			}
			else if (num == 8)
			{
				normal[0] = 0.0;
				normal[1] = 0.9;
				normal[2] = 0.0;
			}
			else if (num == 16)
			{
				normal[0] = 0.0;
				normal[1] = 0.0;
				normal[2] = 0.0;
			}
			else if (num == 32)
			{
				normal[0] = 0.0;
				normal[1] = 0.0;
				normal[2] = 0.0;
			}
		}
	}
	else if (obj[0][0] == 4) // quad
	{
		res = intersec_quad( obj[OBJ_LOC0], e, v, tmin, 0);
	}
	return res;
}

void traceray( scene_t *scene, v3 _e, v3 _v, v3 col)
{
	double tmin = BIG;
	v3 objnormal;
	int res;
	int i;
	v3 *obj = 0;
	v3 col0;
	mult3( col0, 0);
	for (i = 0; i < scene->nobjs; i++)
	{
		v3 *cur_obj = &scene->objs[i * LEN_OBJ];
		v3 cur_objnormal;
		double t = 0;
		res = intersec_obj( cur_obj, _e, _v, &t, cur_objnormal, col0);
		if (res && (t < tmin))
		{
			tmin = t;
			obj = cur_obj;
			copy3( objnormal, cur_objnormal);
		}
	}
	if (!obj)
		return;

	v3 pinter; // intersection pos
	copy3( pinter, _v);
	mult3( pinter, tmin);
	sum3( pinter, _e, pinter);

	v3 vcol;
	copy3( vcol, obj[OBJ_COLOR]);
#if 1
	//	copy3( vcol, col0);		// start with black
	for (i = 0; i < scene->nlamps; i++)
	{
		v3 vdist;	// this is shadow (or enlightment) vec for this lamp
		diff3( vdist, scene->lamps[i], pinter);

		double colin = dot3( vdist, objnormal);
		if (colin <= 0)
		{
			mult3( vcol, 0);
			break;
		}

		int j;
//		double tshad = BIG;
		res = 0;
		v3 nvdist;
		copy3( nvdist, vdist);
		div3( nvdist, norm3( nvdist)); // normalize
		for (j = 0; j < scene->nobjs; j++)
		{
			v3 *cur_obj = &scene->objs[i * LEN_OBJ];
			double t = 0;
			v3 objnormal2;
			v3 col1;
			res = intersec_obj( cur_obj, pinter, vdist, &t, objnormal2, col1);
			if (res)
				break;
		}
		if (res)	// obj found that shadows lamp
			mult3( vcol, 0);
		else
		{
			double att = 900;
			double dist = norm3( vdist) + 0.01;
			att /= dist * dist;
			if (att > 1.0)
				att = 1.0;
			mult3( vcol, att);
		}
	}
#endif
	copy3( col, vcol);
}

// input : s0,s1,s2,e
int compute_scene( scene_t *scene, int i, int j, int *r, int *g, int *b)
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
	// normalize up
	div3( scene->up, norm3( scene->up));
	cross3( right, scene->v, scene->up);
	cross3( scene->up, right, scene->v);
	sum3( scene->s0, scene->e, scene->v);
	sum3( scene->s1, scene->s0, scene->up);
	sum3( scene->s2, scene->s0, right);
	return 0;
}

int load_scene( scene_t *scene, char *file)
{
	int i, j;
	FILE *in = fopen( file, "rt");
	if (in)
	{
		scene->nobjs = 0;
		int nspheres = 0;
		int nfacets = 0;
		int nboxes = 0;
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
		scene->nobjs = nspheres + nfacets + nboxes;
		if (scene->nobjs)
		{
		scene->objs = malloc( scene->nobjs * sizeof(v3) * LEN_OBJ);
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
				if (j == 0)
					scene->objs[i * LEN_OBJ][0] = 1; // type
				scene->objs[i * LEN_OBJ + OBJ_LOC0 + j][0] = p[0];
				scene->objs[i * LEN_OBJ + OBJ_LOC0 + j][1] = p[1];
				scene->objs[i * LEN_OBJ + OBJ_LOC0 + j][2] = p[2];
				if (++j >= 3)
				{
#ifndef USE_NORMAL
					scene->objs[i * LEN_OBJ + OBJ_COLOR][0] = !(i % 3);	// fake facet color with facet index r
					scene->objs[i * LEN_OBJ + OBJ_COLOR][1] = !((i + 1) % 3);	// g
					scene->objs[i * LEN_OBJ + OBJ_COLOR][2] = !((i + 2) % 3);	// b
#else
					scene->objs[i * LEN_OBJ + OBJ_COLOR][0] = normal[0];	// fake facet color with facet normal r
					scene->objs[i * LEN_OBJ + OBJ_COLOR][1] = normal[1];	// g
					scene->objs[i * LEN_OBJ + OBJ_COLOR][2] = normal[2];	// b
#endif
					i++;
					j = 0;
				}
			}
		}
		}
		fclose( in);
	}
	printf( "nobjs=%d\n", scene->nobjs);
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
			while (((*w->start == start) && !running && !*w->quit) && SDL_CondWait( w->cond_go, w->mutex_go) == 0)
				// this blocks without hogging cpu
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

#if 1
// 3d scene :
#define F 40
#define S 20
// objects
v3 __objs[] = {
#if 1
	// a sphere
	{ 0, S/4, 0 },	// 0=sphere
	{ S, S, -S/2 },
	{  },
	{  },
	{ 1, 1, 1 },	// color
// a facet
	{ 1, 0, 0 },	// 1=facet
	{ 0, 0, 0 },
	{ F, 0, 0 },
	{ 0, F, 0 },
	{ 1, 0, 0 },	// color
// a facet
	{ 1, 0, 0 },	// 1=facet
	{ 0, 0, 0 },
	{ 0, F, 0 },
	{ 0, 0, F },
	{ 0, 1, 0 },	// color
// a facet
	{ 1, 0, 0 },	// 1=facet
	{ 0, 0, 0 },
	{ 0, 0, F },	// WARNING : visible faces must have normal pointing to eye !
	{ F, 0, 0 },
	{ 0, 0, 1 },	// color
#if 1
	// a box
	{ 2, 0, 0 },	// 2=box
	{ -S, -S, -S },	// lower
	{ S, S, S },	// upper
	{},
	{ 1, 1, 0 },	// color
#endif
#else
	// a quad
	{ 2, 0, 0 },	// 2=box
	{ -S, -S, -S },	// lower
	{ S, S, S },	// upper
	{},
	{ 1, 1, 0 },	// color
#endif
};
int _nobjs = sizeof( __objs) / sizeof( __objs[0]) / LEN_OBJ;
v3 *_objs = __objs;

// lamps
v3 __lamps[] = {
// a lamp
	{ F, F, -F },
};
int _nlamps = sizeof( __lamps) / sizeof( __lamps[0]);
v3 *_lamps = __lamps;

// camera is : eye coordinate (vector e)..
#if 1
#define E 30
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

#else
// plate4
// objects
v3 __objs[] = {
	// a sphere
	{ 0, 40, 0 },	// 0=sphere
	{ 210, 55, -80 },
	{  },
	{  },
	{ .5, .3, .5 },	// color
	// a sphere
	{ 0, 35, 0 },	// 0=sphere
	{ 200, 90, 70 },
	{  },
	{  },
	{ .3, .5, .9 },	// color
#if 0
	// a parallelogram
	{ 3, 0, 0 },	// 3=parallelogram
	{ -10000, 0, -10000 },
	{ 10000, 0, -10000 },
	{ -10000, 0, 10000 },
	{ 1, 0, 0 },	// color
#endif
#if 1
#define T 3
#else
#define T 1
#endif	// a parallelogram
	{ T, 0, 0 },	// 3=parallelogram
	{ 150, 0, -15 },
	{ 150, 0, 15 },
	{ 150, 60, -15 },
	{ 1, .5, .5 },	// color
	// a parallelogram
	{ T, 0, 0 },	// 3=parallelogram
	{ 150, 0, -15 },
	{ 300, 0, -15 },
	{ 150, 60, -15 },
	{ .5, 1, .5 },	// color
	// a parallelogram
	{ T, 0, 0 },	// 3=parallelogram
	{ 150, 0, 15 },
	{ 300, 0, 15 },
	{ 150, 60, 15 },
	{ .5, .5, 1 },	// color
	// a parallelogram
	{ T, 0, 0 },	// 3=parallelogram
	{ 150, 60, -15 },
	{ 150, 60, 15 },
	{ 300, 60, -15 },
	{ 1, .5, 1 },	// color

	// a parallelogram
	{ T, 0, 0 },	// 3=parallelogram
	{ 200, 60, -15 },
	{ 200, 60, 15 },
	{ 200, 120, -15 },
	{ 1, 1, .5 },	// color
	// a parallelogram
	{ T, 0, 0 },	// 3=parallelogram
	{ 200, 60, -15 },
	{ 250, 60, -15 },
	{ 200, 120, -15 },
	{ .5, .5, .5 },	// color
	// a parallelogram
	{ T, 0, 0 },	// 3=parallelogram
	{ 200, 60, 15 },
	{ 250, 60, 15 },
	{ 200, 120, 15 },
	{ 1, 1, 1 },	// color
};
int _nobjs = sizeof( __objs) / sizeof( __objs[0]) / LEN_OBJ;
v3 *_objs = __objs;

// lamps
v3 __lamps[] = {
#if 1
	// a lamp
	{ 120, 120, -50 },
#endif
#if 0
	// a lamp
	{ 120, 150, 80 },
#endif
};
int _nlamps = sizeof( __lamps) / sizeof( __lamps[0]);
v3 *_lamps = __lamps;

// camera is : eye coordinate (vector e)..
v3 _e = { -20, 70, -40 };
// ..a direction (vector v)..
v3 _v = { 200, 50, 0 };
// ..and an "up" (vector up) (camera "head" rotation, default pointing to the "sky")
v3 _up = { -190, -70, 40};
//v3 _up = { 0, 1, 0};

#endif

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
	
	SDL_Rect rect;
	Uint32 col;
	int r, g, b;
	rect.x = 0;
	rect.y = 0;
	rect.w = w;
	rect.h = h;
	r = 64;
	g = 64;
	b = 64;
	col = SDL_MapRGB( screen->format, r, g, b);
	SDL_FillRect( screen, &rect, col);
//	SDL_UpdateRect( screen, 0, 0, 0, 0);
	
	SDL_EnableKeyRepeat( SDL_DEFAULT_REPEAT_DELAY, SDL_DEFAULT_REPEAT_INTERVAL);
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
	scene.lamps = _lamps;
	scene.nlamps = _nlamps;
	scene.w = w;
	scene.h = h;
	scene.e = _e;
	scene.v = _v;
	scene.up = _up;
	p.scene = &scene;
	if (scene_file)	// load scene ?
		load_scene( &scene, scene_file);
	printf( "using %d objs\n", scene.nobjs);
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
				case SDL_KEYDOWN:
					switch (event.key.keysym.sym)
					{
						case SDLK_ESCAPE:
							end = 1;
							break;
						case SDLK_PAGEUP:
							if (done)
							{
								v3 tmp;
								sum3( scene.e, scene.e, mult3( copy3( tmp, scene.v), 10));
								go = 1;
							}
							break;
						case SDLK_PAGEDOWN:
							if (done)
							{
								v3 tmp;
								diff3( scene.e, scene.e, mult3( copy3( tmp, scene.v), 10));
								go = 1;
							}
							break;
#define DELT 10
						case SDLK_LEFT:
							if (done)
							{
								scene.e[2] += DELT;
								go = 1;
							}
							break;
						case SDLK_RIGHT:
							if (done)
							{
								scene.e[2] -= DELT;
								go = 1;
							}
							break;
						case SDLK_UP:
							if (done)
							{
								scene.e[1] += DELT;
								go = 1;
							}
							break;
						case SDLK_DOWN:
							if (done)
							{
								scene.e[1] -= DELT;
								go = 1;
							}
							break;
						case SDLK_HOME:
							if (done)
							{
								scene.e[0] += DELT;
								go = 1;
							}
							break;
						case SDLK_END:
							if (done)
							{
								scene.e[0] -= DELT;
								go = 1;
							}
							break;
						default:
							break;
					}
					break;
				case SDL_MOUSEBUTTONUP:
					if (done)
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
