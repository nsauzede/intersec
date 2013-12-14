#include <stdio.h>

#include <SDL.h>
#include <SDL_thread.h>
#include <SDL_mutex.h>

#if 0
#define dprintf printf
#else
#define dprintf(...) do{}while(0)
#endif

// non-mutable, or protected (workers_done, start..)
typedef struct {
	SDL_Surface *screen;
	int w;
	int h;
	int tot;

	SDL_mutex *mutex;
	SDL_cond *cond;
	int *workers_done;

	SDL_mutex *mutex_go;
	SDL_cond *cond_go;
	int *start;
	int *quit;
} payload_t;

typedef struct {
	int num;
	void *payload;
	SDL_sem *sem_init;
} thr_init_t;

int thr( void *opaque)
{
	payload_t *p = 0;
	int num = -1;

	if (!opaque)
		return -1;
	else
	{
		thr_init_t *ti = opaque;
		p = ti->payload;
		num = ti->num;
		SDL_SemPost( ti->sem_init);		// signal master we have done with reading thr init
	}

	int x1, y1, x2, y2;
	int x;
	int y;
	int c = 0;
	int running = 0;
	int start = *p->start;
	dprintf( "thr %d ready to work\n", num);	
	while (1)
	{
		// now wait for master to signal go
		int end = 0;
		while (!end)
		{
			SDL_LockMutex( p->mutex_go);
			while (((*p->start == start) && !running && !*p->quit) && SDL_CondWait( p->cond_go, p->mutex_go) == 0)				// this blocks without hogging cpu
				continue;
			if ((*p->start != start) || running || *p->quit)
				end = 1;
			if (*p->start != start)
			{
				running = 0;
				start = *p->start;
			}
			SDL_UnlockMutex( p->mutex_go);
		}
		if (*p->quit)
		{
			dprintf( "thr %d will quit\n", num);	
			break;
		}
		if (!running)
		{
			dprintf( "thr %d starting %d\n", num, start);
			c++;
			x1 = 0;
			x2 = p->w - 1;
			y1 = num * p->h / p->tot;
			y2 = (num + 1) * p->h / p->tot - 1;
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
		col = SDL_MapRGB( p->screen->format, (c+num)&1?2555:0, (c+num)&2?255:0, (c+num)&4?255:0);
		SDL_FillRect( p->screen, &rect, col);
		SDL_Delay( 1);
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
				SDL_LockMutex( p->mutex);
				dprintf( "thr %d done %d\n", num, start);
				(*p->workers_done)++;
				SDL_UnlockMutex( p->mutex);
				SDL_CondSignal( p->cond);
			}
		}
	}
	dprintf( "thr %d quitting\n", num);
	return 0;
}

int main( int argc, char *argv[])
{
	SDL_Surface *screen = 0;
	SDL_Init( SDL_INIT_VIDEO);
	atexit( SDL_Quit);
	int w = 100;
	int h = 100;
	int bpp = 32;
	screen = SDL_SetVideoMode( w, h, bpp, 0);
	int n = 2;
	int arg = 1;
	if (arg < argc)
		sscanf( argv[arg++], "%d", &n);
	int i;
	SDL_mutex * mutex = SDL_CreateMutex();
	SDL_cond *cond = SDL_CreateCond();
	SDL_mutex * mutex_go = SDL_CreateMutex();
	SDL_cond *cond_go = SDL_CreateCond();
	int workers_done = 0;
	int start = 0;
	int quit = 0;
	SDL_sem *sem_init = SDL_CreateSemaphore( 0);
	payload_t p;
	p.screen = screen;
	p.w = w;
	p.h = h;
	p.tot = n;
	p.mutex = mutex;
	p.cond = cond;
	p.workers_done = &workers_done;
	p.start = &start;
	p.quit = &quit;
	p.mutex_go = mutex_go;
	p.cond_go = cond_go;
	thr_init_t ti;
	ti.payload = &p;
	ti.sem_init = sem_init;
	SDL_Thread **thrs = malloc( n * sizeof( SDL_Thread*));
	for (i = 0; i < n; i++)
	{
		ti.num = i;
		thrs[i] = SDL_CreateThread( thr, &ti);
		SDL_SemWait( sem_init);
	}
	SDL_DestroySemaphore( sem_init);
	int end = 0;
	while (!end)
	{
		static int go = 1;
		static int done = 1;
		SDL_Event event;
		while (SDL_PollEvent( &event))
		{
			switch (event.type)
			{
				case SDL_QUIT:
				case SDL_KEYUP:
					end = 1;
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
		SDL_UpdateRect( screen, 0, 0, 0, 0);
		
#ifdef USE_TIME
		static Uint32 old_t = 0;
#endif
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
#ifdef USE_TIME
				old_t = SDL_GetTicks();
#endif
			}
			SDL_UnlockMutex( mutex);
		}
		else
			SDL_Delay( DELAY_WAIT_WORKERS);

//		if (done)
		{
#ifdef USE_TIME
			Uint32 t = SDL_GetTicks();
#define DELAY_RESTART_WORKERS 1000
			if (t > (old_t + DELAY_RESTART_WORKERS))
#else
			if (go)
#endif
			{
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
