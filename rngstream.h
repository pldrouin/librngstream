/**
 * @file rngstream.h
 * @brief Optimised C implementation of the RNG Stream algorithm.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 * Original public domain algorithm and implementation from Pierre L'Ecuyer, University of
 * Montreal. Also based on public domain optimised implementation of MRG32k3a by Sebastiano Vigna (vigna@acm.org).
 * Indicated generation times were obtained on computer with a i7-6600U CPU @
 * 2.60GHz with a static library compiled using GCC 9.30 with compiler flags
 * "-O3 -march=native -mieee-fp -pipe -Wall -g -Wcast-align" on Ubuntu 20.04.1
 * LTS x86_64.
 */

#ifndef RNGSTREAM_H
#define RNGSTREAM_H
 
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>

#define __rngstream_m1    INT64_C(4294967087)
#define __rngstream_m2    INT64_C(4294944443)
#define __rngstream_a12   INT64_C(1403580)
#define __rngstream_a13n  INT64_C(810728)
#define __rngstream_a21   INT64_C(527612)
#define __rngstream_a23n  INT64_C(1370589)
#define __rngstream_corr1 INT64_C(3482050076509336) //m1 * a13n
#define __rngstream_corr2 INT64_C(5886603609186927) //m2 * a23n
#define __rngstream_m1m2  INT64_C(4294967085) //m1 - 2
#define __rngstream_two17 INT64_C(131072)
#define __rngstream_two53 INT64_C(9007199254740992)

// The following are the transition matrices of the two MRG components
// (in matrix form), raised to the powers -1, 1, 2^76, and 2^127, resp.

static const int64_t __rngstream_InvA1[3][3] = {          // Inverse of A1p0
       { 184888585,   0,  1945170933 },
       {         1,   0,           0 },
       {         0,   1,           0 }
       };

static const int64_t __rngstream_InvA2[3][3] = {          // Inverse of A2p0
       {      0,  360363334,  4225571728 },
       {      1,          0,           0 },
       {      0,          1,           0 }
       };

static const int64_t __rngstream_A1p0[3][3] = {
       {       0,        1,       0 },
       {       0,        0,       1 },
       { -810728,  1403580,       0 }
       };

static const int64_t __rngstream_A2p0[3][3] = {
       {        0,        1,       0 },
       {        0,        0,       1 },
       { -1370589,        0,  527612 }
       };

static const int64_t __rngstream_A1p76[3][3] = {
       {      82758667, 1871391091, 4127413238 },
       {    3672831523,   69195019, 1871391091 },
       {    3672091415, 3528743235,   69195019 }
       };

static const int64_t __rngstream_A2p76[3][3] = {
       {    1511326704, 3759209742, 1610795712 },
       {    4292754251, 1511326704, 3889917532 },
       {    3859662829, 4292754251, 3708466080 }
       };

static const int64_t __rngstream_A1p127[3][3] = {
       {    2427906178, 3580155704,  949770784 },
       {     226153695, 1230515664, 3580155704 },
       {    1988835001,  986791581, 1230515664 }
       };

static const int64_t __rngstream_A2p127[3][3] = {
       {    1464411153,  277697599, 1610723613 },
       {      32183930, 1464411153, 1022607788 },
       {    2824425944,   32183930, 2093834863 }
       };

//-------------------------------------------------------------------------
// The default seed of the package; will be the seed of the first
// declared RNGStream, unless SetPackageSeed is called.
//
static uint64_t rng_nextseed[6] =
{
   12345, 12345, 12345, 12345, 12345, 12345
};

typedef struct
{
	uint64_t Cg[6], Bg[6], Ig[6];  //!< Values can fit in uint32_t, but using uint64_t for faster processing speed on 64bit CPUs
	uint64_t fill64;
	uint32_t fill32;
	uint8_t favail32;
	uint8_t favail64;
} rng_stream;

uint64_t rng_multmodm (int64_t a, const int64_t s, const uint64_t c, const uint64_t m);
void rng_matvecmodm (const int64_t A[3][3], uint64_t const* const s, uint64_t* const v, const uint64_t m);
void rng_matvecmodmll (const int64_t A[3][3], int64_t const* const s, int64_t* const v, const uint64_t m);
void rng_matmatmodm (const int64_t A[3][3], const int64_t B[3][3], int64_t C[3][3], const uint64_t m);
void rng_mattwopowmodm (const int64_t A[3][3], int64_t B[3][3], const uint64_t m, const long e);
void rng_matpowmodm (const int64_t A[3][3], int64_t B[3][3], const uint64_t m, long n);
int rng_checkseed (const uint64_t seed[6]);

void rng_init(rng_stream* s);
bool rng_setseed(rng_stream* s, uint64_t const* const seed); //!< Values can fit in uint32_t, but using uint64_t for faster processing speed on 64bit CPUs
inline static void rng_advanceseed(uint64_t const* seedin, uint64_t* seedout){rng_matvecmodm (__rngstream_A1p127,seedin,seedout,__rngstream_m1); rng_matvecmodm (__rngstream_A2p127,&seedin[3],&seedout[3],__rngstream_m2);}
inline static bool rng_setpackageseed(uint64_t const* const seed){if(rng_checkseed(seed)) return false; memcpy(rng_nextseed,seed,6*sizeof(uint64_t)); return true;} //!< Values can fit in uint32_t, but using uint64_t for faster processing speed on 64bit CPUs
inline static void rng_resetstartstream(rng_stream* s){int i; for(i = 5; i >=0; --i) s->Cg[i] = s->Bg[i] = s->Ig[i];}
inline static void rng_resetstartsubstream(rng_stream* s){memcpy(s->Cg,s->Bg,6*sizeof(uint64_t));}
inline static void rng_resetnextsubstream(rng_stream* s){rng_matvecmodm(__rngstream_A1p76, s->Bg, s->Bg, __rngstream_m1); rng_matvecmodm(__rngstream_A2p76, s->Bg+3, s->Bg+3, __rngstream_m2); memcpy(s->Cg,s->Bg,6*sizeof(uint64_t));}
void rng_advancestate(rng_stream* s, const long e, const long c);
inline static void rng_getstate(rng_stream* s, uint64_t* const seed){memcpy(seed,s->Cg,6*sizeof(uint64_t));}
void rng_writestate(rng_stream* s);
void rng_writestatefull(rng_stream* s);

/**
 * @brief Uniform deviate in the interval [0,m1-1], with m1=4294967087.
 *
 * This function returns a uniform deviate in the interval [0,m1-1]. The original
 * U01 function was returning a uniform deviate in the interval [1,m1], before
 * multiplying by norm equal to 1/(m1+1), so U01 was returning a value in the
 * interval ]0,1). Generation time is less than 4.294 ns.
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,m1-1].
 */
inline static uint64_t rng_rand_m1(rng_stream* s)
{
  int64_t r=(int64_t)s->Cg[2]-s->Cg[5];
  r-=__rngstream_m1*(r>>63);

  //printf("%16f %16f %16f %16f %16f %16f\n",(double)Cg[0],(double)Cg[1],(double)Cg[2],(double)Cg[3],(double)Cg[4],(double)Cg[5]);

  /* Component 1 */
  uint64_t p=(__rngstream_a12 * s->Cg[1] + __rngstream_corr1 - __rngstream_a13n * s->Cg[0])%__rngstream_m1;

  s->Cg[0]=s->Cg[1];
  s->Cg[1]=s->Cg[2];
  s->Cg[2]=p;

  /* Component 2 */
  p=(__rngstream_a21 * s->Cg[5] + __rngstream_corr2 - __rngstream_a23n * s->Cg[3])%__rngstream_m2;

  s->Cg[3]=s->Cg[4];
  s->Cg[4]=s->Cg[5];
  s->Cg[5]=p;

  //printf("%16f %16f\n",(double)p1,(double)p2);

  /* Combination */
  return r;
}

/**
 * @brief Uniform deviate in the interval [1,m1], with m1=4294967087.
 *
 * This function returns a uniform deviate in the interval [1,m1]. Generation
 * tims is less than 4.294 ns.
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [1,m1].
 */
inline static uint64_t rng_rand_pm1(rng_stream* s)
{
  int64_t r=(int64_t)s->Cg[2]-s->Cg[5];
  r-=__rngstream_m1*((r-1)>>63);

  //printf("%16f %16f %16f %16f %16f %16f\n",(double)Cg[0],(double)Cg[1],(double)Cg[2],(double)Cg[3],(double)Cg[4],(double)Cg[5]);

  /* Component 1 */
  uint64_t p=(__rngstream_a12 * s->Cg[1] + __rngstream_corr1 - __rngstream_a13n * s->Cg[0])%__rngstream_m1;

  s->Cg[0]=s->Cg[1];
  s->Cg[1]=s->Cg[2];
  s->Cg[2]=p;

  /* Component 2 */
  p=(__rngstream_a21 * s->Cg[5] + __rngstream_corr2 - __rngstream_a23n * s->Cg[3])%__rngstream_m2;

  s->Cg[3]=s->Cg[4];
  s->Cg[4]=s->Cg[5];
  s->Cg[5]=p;

  //printf("%16f %16f\n",(double)p1,(double)p2);

  /* Combination */
  return r;
}

/**
 * @brief Uniform deviate in the interval [0,2^24-1].
 *
 * This function returns a uniform deviate in the interval [0,2^24-1].
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,2^24-1].
 */
#define RNG_RAND24(s) (rng_rand_m1(s)>>8)

/**
 * @brief Uniform deviate in the interval [0,72057590531489791].
 *
 * This function returns a uniform deviate in the interval [0,72057590531489791]. Generation time is less than 7.975 ns (same as rng_rand32 without bit recycling).
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,72057590531489791].
 */
inline static uint64_t rng_rand_m1_24(rng_stream *s){return (rng_rand_m1(s)<<24)|(rng_rand_m1(s)>>8);}

/**
 * @brief Uniform deviate in the interval [0,2^32-1].
 *
 * This function returns a uniform deviate in the interval [0,2^32-1].
 * Generation time is less than 5.680 ns. Generation time without bit recycling
 * is less than 7.933 ns.
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,2^32-1].
 */
inline static uint32_t rng_rand32(rng_stream *s){
  if(!s->favail32) {
    s->fill32=RNG_RAND24(s);
    uint32_t ret=RNG_RAND24(s)|(s->fill32<<24);
    s->favail32=2;
    s->fill32>>=8;
    return ret;
  }
  uint32_t ret=RNG_RAND24(s)|(s->fill32<<24);
  s->favail32-=1;
  s->fill32>>=8;
  return ret;
}
/*
inline static uint32_t rng_rand32(rng_stream *s){
  return (RNG_RAND24(s)<<8)|(rng_rand_m1(s)>>24);
}
*/

/**
 * @brief Uniform deviate in the interval [0,2^64-1].
 *
 * This function returns a uniform deviate in the interval [0,2^64-1].
 * Generation time is less than 14.08 ns. Generation time with previous
 * algorithm without bit recycling was less than 15.56 ns.
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,2^64-1].
 */
inline static uint64_t rng_rand64(rng_stream *s){
  if(!s->favail64) {
    s->fill64=RNG_RAND24(s);
    uint64_t ret=RNG_RAND24(s)|(RNG_RAND24(s)<<24)|(s->fill64<<48);
    s->favail64=1;
    s->fill64>>=16;
    return ret;

  } else if(s->favail64==1) {
    s->fill64|=(RNG_RAND24(s)<<8);
    uint64_t ret=RNG_RAND24(s)|(RNG_RAND24(s)<<24)|(s->fill64<<48);
    s->favail64=2;
    s->fill64>>=16;
    return ret;
  }
  uint64_t ret=RNG_RAND24(s)|(RNG_RAND24(s)<<24)|(s->fill64<<48);
  s->favail64=0;
  s->fill64>>=16;
  return ret;
}
/*
inline static uint64_t rng_rand64(rng_stream *s){
  return (RNG_RAND24(s)<<48)|(RNG_RAND24(s)<<24)|RNG_RAND24(s);
}
*/

/**
 * @brief Uniform deviate in the interval [0,1), with a non-truncated minimum spacing of 1/4294967087.
 *
 * This function returns a uniform deviate in the interval [0,1). The minimum
 * non-truncated spacing between two distinct values is 1/4294967087. Generation
 * time is less than 4.292 ns (very similar to rng_rand_m1).
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,1).
 */
inline static double rng_rand_u01(rng_stream *s){return rng_rand_m1(s)*0x1.000000d10000bp-32;}

/**
 * @brief Uniform deviate in the interval (0,1], with a non-truncated minimum spacing of 1/4294967087.
 *
 * This function returns a uniform deviate in the interval (0,1]. The minimum
 * non-truncated spacing between two distinct values is 1/4294967087. Generation
 * time is less than 4.292 ns (very similar to rng_rand_pm1).
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,1).
 */
inline static double rng_rand_pu01(rng_stream *s){return rng_rand_pm1(s)*0x1.000000d10000bp-32;}

/**
 * @brief Uniform deviate in the interval [0,1], with a non-truncated minimum spacing of 1/72057590531489792.
 *
 * This function returns a uniform deviate in the interval [0,1]. The minimum
 * non-truncated spacing between two distinct values is 1/72057590531489792.
 * Upper bound is included only due to double precision truncation error.
 * Generation time is less than 7.959 ns.
 * No advantage of using this function over rng_rand_u01d as it is not any
 * faster (maybe it is even slightly slower) and it has a larger minimum
 * spacing.
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,1).
 */
inline static double rng_rand_u01e(rng_stream *s){return rng_rand_m1_24(s)*0x1.000000d10000bp-56;}

/**
 * @brief Uniform deviate in the interval [0,1], with a non-truncated minimum spacing of 1/18446742278413265569.
 *
 * This function returns a uniform deviate in the interval [0,1]. The minimum
 * non-truncated spacing between two distinct values is 1/18446742278413265569.
 * Upper bound is included only due to double precision truncation error.
 * Generation time is less than 7.938 ns.
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,1].
 */
inline static double rng_rand_u01d(rng_stream *s){return rng_rand_m1(s)*0x1.000001a200020p-64+rng_rand_m1(s)*0x1.000000d10000bp-32;}
//inline static double rng_rand_u01d(rng_stream *s){return (rng_rand_m1(s)+rng_rand_m1(s)*0x1.000000d10000bp-32)*0x1.000000d10000bp-32;}

/**
 * @brief Uniform deviate in the interval [0,1), with a minimum spacing of
 * 2^85/2097150.
 *
 * This function returns a uniform deviate in the interval [0,1). The minimum
 * non-truncated spacing between two distinct values is 2^85/2097150, used only
 * in the interval [4294967086/4294967087,0), with a non-truncated spacing of
 * 1/18446742278413265569 elsewhere.
 * It is a slightly modified version of rng_rand_u01dm, where the multiplicating
 * factor for the 2nd order random number is slightly biased down in the
 * mentioned interval to ensure that the upper bound remains excluded.
 * Generation time of less than 7.977 ns was measured (very similar to
 * rng_rand_u01d).
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval [0,1).
 */
inline static double rng_rand_u01dm(rng_stream *s){int64_t first=rng_rand_m1(s); return rng_rand_m1(s)*(0x1.000001a200020p-64+((__rngstream_m1m2-first)>>63)*0x1.1a200020p-84)+first*0x1.000000d10000bp-32;}

/**
 * @brief Uniform deviate in the interval (0,1], with a non-truncated minimum spacing of 1/18446742278413265569.
 *
 * This function returns a uniform deviate in the interval (0,1]. The minimum
 * non-truncated spacing between two distinct values is 1/18446742278413265569.
 * Generation time is less than 8.032 ns (barely slower than rng_rand_u01d).
 *
 * @param s: Handle to rng_stream.
 * @return uniform deviate in the interval (0,1].
 */
inline static double rng_rand_pu01d(rng_stream *s){return rng_rand_m1(s)*0x1.000001a200020p-64+rng_rand_pm1(s)*0x1.000000d10000bp-32;}
//inline static double rng_rand_pu01d(rng_stream *s){double ret=rng_rand_u01d(s); while(ret==0) ret=rng_rand_u01d(s); return ret;}
//inline static double rng_rand_pu01d(rng_stream *s){return (rng_rand_m1(s)+rng_rand_pm1(s)*0x1.000000d10000bp-32)*0x1.000000d10000bp-32;}

#endif


