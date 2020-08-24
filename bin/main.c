#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "rngstream.h"

int main(int nargs, char* args[])
{
  //int64_t sum=0;
  //int64_t i;
  uint32_t r;
  rng_stream s;
  rng_init(&s);

  //double sum=0;
  /*
  int64_t sump=0;
  int64_t sumn=0;
  int64_t sumpd=0;
  int64_t sumnd=0;
  int64_t sumpdp=0;
  int64_t sumndp=0;
  int64_t ret;
  */
  //const int64_t niters=10000000000;
  //const int64_t niters=1000000;
  //const int64_t niters=1000000000;
  //const uint32_t pdomain=m1-m2;

  //rng_writestate(&s);

  //for(i=0; i<niters; ++i) {
    //rng_writestate(&s);
  for(;;) {
    //sum+=rng_rand_m1(&s);
    //sum+=rng_rand_pm1(&s);
    //sum+=rng_rand_m1_24(&s);
    //sum+=rng_rand32(&s);
    //sum+=rng_rand64(&s);
    //sum+=rng_rand_m1_24(&s);
    //sum+=rng_rand_u01(&s);
    //sum+=rng_rand_pu01(&s);
    //sum+=rng_rand_u01e(&s);
    //sum+=rng_rand_u01d(&s);
    //sum+=rng_rand_u01dm(&s);
    //sum+=rng_rand_pu01d(&s);
    //sum+=rng_rand_pu01d(&s);
    //r=(rng_rand_m1(&s)<<16)|(rng_rand_m1(&s)&0xFFFF);;
    //r=rng_rand_m1(&s);
    //r=rng_rand32(&s);
    r=(uint32_t)(rng_rand_u01d(&s)*4294967296);
    write(1, &r,4);
    /*
    ret=rng_rand_m1(&s);

    if(ret>=0) {
      ++sump;
      
      if(ret<=pdomain) ++sumpd;
      else if(ret>=m2) ++sumpdp;
      
    } else {
     ++sumn;
    
     if(ret<=-m2) ++sumnd;
     else if(ret>=m2-m1) ++sumndp;
    //sum+=(rng_rand_m1(&s)<=pdomain);
    }
    */
  }
  //sum/=niters;
  //printf("%22.15e +/- %22.15e vs %22.15e\n",sum,sqrt(((double)pdomain+1)/m1*((double)(m1-pdomain+1))/((double)m1*niters)),((double)pdomain+1)/m1);
  //printf("pos: %20" PRIi64 "/%20" PRIi64 " (%22.15e)\tneg: %20" PRIi64 "/%20" PRIi64 " (%22.15e) vs %22.15e\n",sumpd,sump,((double)sumpd)/sump,sumnd,sumn,((double)sumnd)/sumn,((double)pdomain)/m1);
  //printf("pos: %20" PRIi64 "/%20" PRIi64 " (%22.15e)\tneg: %20" PRIi64 "/%20" PRIi64 " (%22.15e) vs %22.15e\n",sumpdp,sump,((double)sumpdp)/sump,sumndp,sumn,((double)sumndp)/sumn,((double)pdomain)/m1);
  //printf("Sum is %22.15e\n",sum);
  fflush(stdout);
  return 0;
}
