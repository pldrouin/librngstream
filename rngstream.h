/***********************************************************************\
 *
 * Original implementation:      Pierre L'Ecuyer, University of Montreal
 *
\***********************************************************************/

//Modified by Pierre-Luc Drouin (pldrouin@pldrouin.net) in February 2012:
//-Created rng_rand24, that generates true 24 bit pseudo random unsigned integers.
//-Now using more inline functions for increased performance.
//-Now using 64 bit integers instead of double precision variables for
// increased performance on 64 bit CPUs.

#ifndef RNGSTREAM_H
#define RNGSTREAM_H
 
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>

static const int64_t m1   =       4294967087;
static const int64_t m2   =       4294944443;
static const int64_t a12  =       1403580;
static const int64_t a13n =       810728;
static const int64_t a21  =       527612;
static const int64_t a23n =       1370589;
static const int64_t corr1 =	   (m1 * a13n);
static const int64_t corr2 =	   (m2 * a23n);
static const int64_t two17 =      131072;
static const double two53 =      9007199254740992.0;

// The following are the transition matrices of the two MRG components
// (in matrix form), raised to the powers -1, 1, 2^76, and 2^127, resp.

static const int64_t InvA1[3][3] = {          // Inverse of A1p0
       { 184888585,   0,  1945170933 },
       {         1,   0,           0 },
       {         0,   1,           0 }
       };

static const int64_t InvA2[3][3] = {          // Inverse of A2p0
       {      0,  360363334,  4225571728 },
       {      1,          0,           0 },
       {      0,          1,           0 }
       };

static const int64_t A1p0[3][3] = {
       {       0,        1,       0 },
       {       0,        0,       1 },
       { -810728,  1403580,       0 }
       };

static const int64_t A2p0[3][3] = {
       {        0,        1,       0 },
       {        0,        0,       1 },
       { -1370589,        0,  527612 }
       };

static const int64_t A1p76[3][3] = {
       {      82758667, 1871391091, 4127413238 },
       {    3672831523,   69195019, 1871391091 },
       {    3672091415, 3528743235,   69195019 }
       };

static const int64_t A2p76[3][3] = {
       {    1511326704, 3759209742, 1610795712 },
       {    4292754251, 1511326704, 3889917532 },
       {    3859662829, 4292754251, 3708466080 }
       };

static const int64_t A1p127[3][3] = {
       {    2427906178, 3580155704,  949770784 },
       {     226153695, 1230515664, 3580155704 },
       {    1988835001,  986791581, 1230515664 }
       };

static const int64_t A2p127[3][3] = {
       {    1464411153,  277697599, 1610723613 },
       {      32183930, 1464411153, 1022607788 },
       {    2824425944,   32183930, 2093834863 }
       };

//-------------------------------------------------------------------------
// The default seed of the package; will be the seed of the first
// declared RNGStream, unless SetPackageSeed is called.
//
static uint32_t rng_nextseed[6] =
{
   12345, 12345, 12345, 12345, 12345, 12345
};

typedef struct
{
	uint32_t Cg[6], Bg[6], Ig[6];
	uint32_t fill;
	uint8_t favail;
} rng_stream;

uint32_t rng_multmodm (int64_t a, const int64_t s, const uint32_t c, const uint32_t m);
void rng_matvecmodm (const int64_t A[3][3], uint32_t const* const s, uint32_t* const v, const uint32_t m);
void rng_matvecmodmll (const int64_t A[3][3], int64_t const* const s, int64_t* const v, const uint32_t m);
void rng_matmatmodm (const int64_t A[3][3], const int64_t B[3][3], int64_t C[3][3], const uint32_t m);
void rng_mattwopowmodm (const int64_t A[3][3], int64_t B[3][3], const uint32_t m, const long e);
void rng_matpowmodm (const int64_t A[3][3], int64_t B[3][3], const uint32_t m, long n);
int rng_checkseed (const uint32_t seed[6]);

void rng_init(rng_stream* s);
bool rng_setseed(rng_stream* s, uint32_t const* const seed);
inline static void rng_advanceseed(uint32_t const* seedin, uint32_t* seedout){rng_matvecmodm (A1p127,seedin,seedout,m1); rng_matvecmodm (A2p127,&seedin[3],&seedout[3],m2);}
inline static bool rng_setpackageseed(uint32_t const* const seed){if(rng_checkseed(seed)) return false; memcpy(rng_nextseed,seed,6*sizeof(uint32_t)); return true;}
inline static void rng_resetstartstream(rng_stream* s){int i; for(i = 5; i >=0; --i) s->Cg[i] = s->Bg[i] = s->Ig[i];}
inline static void rng_resetstartsubstream(rng_stream* s){memcpy(s->Cg,s->Bg,6*sizeof(uint32_t));}
inline static void rng_resetnextsubstream(rng_stream* s){rng_matvecmodm(A1p76, s->Bg, s->Bg, m1); rng_matvecmodm(A2p76, s->Bg+3, s->Bg+3, m2); memcpy(s->Cg,s->Bg,6*sizeof(uint32_t));}
void rng_advancestate(rng_stream* s, const long e, const long c);
inline static void rng_getstate(rng_stream* s, uint32_t* const seed){memcpy(seed,s->Cg,6*sizeof(uint32_t));}
void rng_writestate(rng_stream* s);
void rng_writestatefull(rng_stream* s);

//This function returns a uniform deviate in the interval [0,m1-1]. The original
//U01 function was returning a uniform deviate in the interval [1,m1], before
//multiplying by norm equal to 1/(m1+1), so U01 was returning a value in the
//interval ]0,1[. 
inline static uint32_t rng_rand32weak(rng_stream* s)
{
  int64_t r=(int64_t)s->Cg[2]-s->Cg[5];
  r-=m1*(r>>63);

  //printf("%16f %16f %16f %16f %16f %16f\n",(double)Cg[0],(double)Cg[1],(double)Cg[2],(double)Cg[3],(double)Cg[4],(double)Cg[5]);

  /* Component 1 */
  uint32_t p=(uint32_t)((a12 * s->Cg[1] - a13n * s->Cg[0] + corr1)%m1);

  s->Cg[0]=s->Cg[1];
  s->Cg[1]=s->Cg[2];
  s->Cg[2]=p;

  /* Component 2 */
  p=(uint32_t)((a21 * s->Cg[5] - a23n * s->Cg[3] + corr2)%m2);

  s->Cg[3]=s->Cg[4];
  s->Cg[4]=s->Cg[5];
  s->Cg[5]=p;

  //printf("%16f %16f\n",(double)p1,(double)p2);

  /* Combination */
  return r;
}

inline static uint32_t rng_rand24(rng_stream *s){return rng_rand32weak(s)>>8;}

inline static uint32_t rng_rand32(rng_stream *s){
  if(!s->favail) {
    s->fill=rng_rand24(s);
    s->favail=2;
    uint32_t ret=rng_rand24(s)|(s->fill<<24);
    s->fill>>=8;
    return ret;
  }
  uint32_t ret=rng_rand24(s)|(s->fill<<24);
  s->fill>>=8;
  s->favail-=1;
  return ret;
}

inline static uint64_t rng_rand64(rng_stream *s){
  if(!s->favail) {
    s->fill=rng_rand24(s);
    s->favail=1;
    uint64_t ret=rng_rand24(s)|((rng_rand24(s)|(((uint64_t)s->fill)<<24))<<24);
    s->fill>>=16;
    return ret;

  } else if(s->favail==1) {
    s->fill|=(rng_rand24(s)<<8);
    s->favail=2;
    uint64_t ret=rng_rand24(s)|((rng_rand24(s)|(((uint64_t)s->fill)<<24))<<24);
    s->fill>>=16;
    return ret;
  }
  uint64_t ret=rng_rand24(s)|((rng_rand24(s)|(((uint64_t)s->fill)<<24))<<24);
  s->fill>>=16;
  s->favail-=2;
  return ret;
}

#endif


