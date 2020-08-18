/***********************************************************************\
 *
 * Original implementation:      Pierre L'Ecuyer, University of Montreal
 *
 \***********************************************************************/

//Modified by Pierre-Luc Drouin (pldrouin@pldrouin.net) in February 2012:
//-Created Rand24, that generates true 24 bit pseudo random unsigned integers.
//-Now using more inline functions for increased performance.
//-Now using 64 bit integers instead of double precision variables for
// increased performance on 64 bit CPUs.

//Different instances can be used safely in a multi-threaded application (initialization of static tables is thread-safe).

#include <stdlib.h>
#include <stdio.h>
#include "rngstream.h"

//*************************************************************************
// Public members of the class start here

void rng_init(rng_stream *s)
{
  /* Information on a stream. The arrays {Cg, Bg, Ig} contain the current
     state of the stream, the starting state of the current SubStream, and the
     starting state of the stream. rng_nextseed
     will be the seed of the next declared RNGStream. */

  int i;

  for (i = 5; i >=0; --i) {
    s->Bg[i] = s->Cg[i] = s->Ig[i] = rng_nextseed[i];
  }
  s->favail=0;

  rng_matvecmodm (A1p127, rng_nextseed, rng_nextseed, m1);
  rng_matvecmodm (A2p127, &rng_nextseed[3], &rng_nextseed[3], m2);
  rng_rand32weak(s);
}

//-------------------------------------------------------------------------
bool rng_setseed(rng_stream* s, uint32_t const* const seed)
{
  if (rng_checkseed (seed)) return false; // FAILURE     
  int i;

  for (i = 5; i >=0; --i) s->Cg[i] = s->Bg[i] = s->Ig[i] = seed[i];
  return true; // SUCCESS
}

//-------------------------------------------------------------------------
// if e > 0, let n = 2^e + c;
// if e < 0, let n = -2^(-e) + c;
// if e = 0, let n = c.
// Jump n steps forward if n > 0, backwards if n < 0.
//
void rng_advancestate(rng_stream* s, const long e, const long c)
{
  int64_t B1[3][3], C1[3][3], B2[3][3], C2[3][3];

  if (e > 0) {
    rng_mattwopowmodm (A1p0, B1, m1, e);
    rng_mattwopowmodm (A2p0, B2, m2, e);
  } else if (e < 0) {
    rng_mattwopowmodm (InvA1, B1, m1, -e);
    rng_mattwopowmodm (InvA2, B2, m2, -e);
  }

  if (c >= 0) {
    rng_matpowmodm (A1p0, C1, m1, c);
    rng_matpowmodm (A2p0, C2, m2, c);
  } else {
    rng_matpowmodm (InvA1, C1, m1, -c);
    rng_matpowmodm (InvA2, C2, m2, -c);
  }

  if (e) {
    rng_matmatmodm (B1, C1, C1, m1);
    rng_matmatmodm (B2, C2, C2, m2);
  }

  rng_matvecmodm (C1, s->Cg, s->Cg, m1);
  rng_matvecmodm (C2, s->Cg+3, s->Cg+3, m2);
}

//-------------------------------------------------------------------------
void rng_writestate(rng_stream* s)
{
  printf("The current state of the Rngstream:\n   Cg = { ");
  int i;

  for (i = 0; i < 5; i++) printf("%" PRIu32 ", ",s->Cg[i]);
  printf("%" PRIu32 " }\n\n",s->Cg[5]);
}


//-------------------------------------------------------------------------
void rng_writestatefull(rng_stream* s)
{
  int i;

  printf("The RNGStream   Ig = { ");

  for (i = 0; i < 5; i++) printf("%" PRIu32 ", ",s->Ig[i]);
  printf("%" PRIu32 " }\n",s->Ig[5]);

  printf("   Bg = { ");
  for (i = 0; i < 5; i++) printf("%" PRIu32 ", ",s->Bg[i]);
  printf("%" PRIu32 " }\n",s->Bg[5]);

  printf("   Cg = { ");
  for (i = 0; i < 5; i++) printf("%" PRIu32 ", ",s->Cg[i]);
  printf("%" PRIu32 " }\n\n",s->Cg[5]);
}

//-------------------------------------------------------------------------
// Compute the vector v = A*s MOD m. Assume that -m < s[i] < m.
// Works also when v = s.
//
void rng_matvecmodm (const int64_t A[3][3], uint32_t const* const s, uint32_t* const v, const uint32_t m)
{
  int i;
  uint32_t x[3];               // Necessary if v = s

  for (i = 2; i >=0; --i) {
    x[i] = rng_multmodm (A[i][0], s[0], 0, m);
    x[i] = rng_multmodm (A[i][1], s[1], x[i], m);
    x[i] = rng_multmodm (A[i][2], s[2], x[i], m);
  }
  memcpy(v,x,3*sizeof(uint32_t));
}

//-------------------------------------------------------------------------
// Compute the vector v = A*s MOD m. Assume that -m < s[i] < m.
// Works also when v = s.
//
void rng_matvecmodmll (const int64_t A[3][3], int64_t const* const s, int64_t* const v, const uint32_t m)
{
    int i;
    uint32_t x[3];               // Necessary if v = s

    for (i = 2; i >=0; --i) {
        x[i] = rng_multmodm(A[i][0], s[0], 0, m);
        x[i] = rng_multmodm(A[i][1], s[1], x[i], m);
        x[i] = rng_multmodm(A[i][2], s[2], x[i], m);
    }
    for (i = 2; i >=0; --i) v[i] = x[i];
}

//-------------------------------------------------------------------------
// Return (a*s + c) MOD m; a, s, c and m must be < 2^35
//
uint32_t rng_multmodm (int64_t a, const int64_t s, const uint32_t c, const uint32_t m)
{
  double v;
  long a1;

  v = (double)(a) * s + c;

  if (v >= two53 || v <= -two53) {
    a1 = (long)(a / two17);    a -= a1 * two17;
    v  = a1 * s;
    a1 = (long)(v / m);     v -= a1 * m;
    v = v * two17 + a * s + c;
  }

  a1 = (long)(v / m);
  /* in case v < 0)*/
  if ((v -= a1 * m) < 0.0) return v += m;   else return (uint32_t)v;
}

//-------------------------------------------------------------------------
// Compute the matrix C = A*B MOD m. Assume that -m < s[i] < m.
// Note: works also if A = C or B = C or A = B = C.
//
void rng_matmatmodm (const int64_t A[3][3], const int64_t B[3][3], int64_t C[3][3], const uint32_t m)
{
  int i, j;
  int64_t V[3], W[3][3];

  for (i = 2; i >=0; --i) {

    for (j = 2; j >=0; --j) V[j] = B[j][i];
    rng_matvecmodmll(A, V, V, m);

    for (j = 2; j >=0; --j) W[j][i] = V[j];
  }
  memcpy(C,W,9*sizeof(int64_t));
}

//-------------------------------------------------------------------------
// Compute the matrix B = (A^(2^e) Mod m);  works also if A = B. 
//
void rng_mattwopowmodm (const int64_t A[3][3], int64_t B[3][3], const uint32_t m, const long e)
{
  long i;

  /* initialize: B = A */
  if (A != B) memcpy(B,A,9*sizeof(int64_t));

  /* Compute B = A^(2^e) mod m */
  for (i = e-1; i >=0; --i) rng_matmatmodm (B, B, B, m);
}


//-------------------------------------------------------------------------
// Compute the matrix B = (A^n Mod m);  works even if A = B.
//
void rng_matpowmodm (const int64_t A[3][3], int64_t B[3][3], const uint32_t m, long n)
{
  int j;
  int64_t W[3][3];

  /* initialize: W = A; B = I */
  memcpy(W,A,9*sizeof(int64_t));
  memset(B,0,9*sizeof(int64_t));

  for (j = 2; j >=0; --j) B[j][j] = 1.0;

  /* Compute B = A^n mod m using the binary decomposition of n */
  while (n > 0) {
    if (n % 2) rng_matmatmodm (W, B, B, m);
    rng_matmatmodm (W, W, W, m);
    n /= 2;
  }
}

//-------------------------------------------------------------------------
// Check that the seeds are legitimate values. Returns 0 if legal seeds,
// -1 otherwise.
//
int rng_checkseed (const uint32_t seed[6])
{
  int i;

  for (i = 0; i < 3; ++i) {
    if (seed[i] >= m1) {
      fprintf(stderr,"****************************************\n"
	"ERROR: Seed[%i] >= 4294967087, Seed is not set."
	"\n****************************************\n\n",i);
      return (-1);
    }
  }
  for (i = 3; i < 6; ++i) {
    if (seed[i] >= m2) {
      fprintf(stderr,"*****************************************\n"
	"ERROR: Seed[%i] >= 4294944443, Seed is not set."
	"\n*****************************************\n\n",i);
      return (-1);
    }
  }
  if (seed[0] == 0 && seed[1] == 0 && seed[2] == 0) {
    fprintf(stderr,"****************************\n"
      "ERROR: First 3 seeds = 0.\n"
      "****************************\n\n");
    return (-1);
  }
  if (seed[3] == 0 && seed[4] == 0 && seed[5] == 0) {
    fprintf(stderr,"****************************\n"
      "ERROR: Last 3 seeds = 0.\n"
      "****************************\n\n");
    return (-1);
  }

  return 0;
}
