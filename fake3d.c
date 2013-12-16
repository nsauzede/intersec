#include <stdio.h>

#include <SDL.h>
#include <SDL_thread.h>
#include <SDL_mutex.h>

#if 0
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
	SDL_Surface *screen;
	int w;
	int h;
} payload_t;

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
	int x1, y1, x2, y2;
	int x;
	int y;
	int c = 0;
	int running = 0;
	int start = *w->start;
	dprintf( "thr %d ready to work\n", num);	
	while (1)
	{
		// now wait for master to signal go
		int end = 0;
		while (!end)
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
			c++;
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
	worker_t wo;
	payload_t p;
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
