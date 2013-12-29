
#ifdef USE_DOUBLE
//floating point
typedef double scalar_t;
#define COEF 1
#define S(d) d
#define D(s) s
#define SMALL S(0.0001)
//#define BIG S(3e30)
#define BIG S(1e10)
#else
// fixed point
#include <inttypes.h>
#if 0
typedef int64_t scalar_t;
#define S(d) ((scalar_t)((double)d * 0x8000000000000000))
#define PRIscalar PRIx64
#define SMALL S(0.0001)
#define BIG S(3e30)
#else
typedef int32_t scalar_t;
#define COEF 1000
#define S(d) ((scalar_t)((double)d * COEF))
#define D(s) ((double)((scalar_t)s / COEF))
#define PRIscalar PRIx32
#define SMALL S((double)1.0/COEF)
#define BIG S(1e6)
#endif
#endif


typedef scalar_t v3[3];

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
#define MIN3(a,b,c) (a < b ? (a < c ? a : c) : b < c ? b : c)
#define MAX3(a,b,c) (a > b ? (a > c ? a : c) : b > c ? b : c)

scalar_t norm3( v3 n);
scalar_t* div3( v3 n, scalar_t nn);
scalar_t* sum3( v3 dest, v3 const src1, v3 const src2);
int intersec_box( v3 lower, v3 upper, v3 e, v3 v, scalar_t *_tmin, scalar_t *_tmax);

