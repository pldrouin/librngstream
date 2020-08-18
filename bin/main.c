#include <stdio.h>
#include <unistd.h>

#include "rngstream.h"

int main(int nargs, char* args[])
{
  //int64_t sum=0;
  int64_t i;
  uint32_t r;
  rng_stream s;
  rng_init(&s);

  //for(i=0; i<100000000; ++i) {
  for(;;) {
    //sum+=rng_rand32weak(&s);
    r=rng_rand32weak(&s);
    //r=rng_rand32(&s);
    write(1, &r,4);
  }
  fflush(stdout);
  //printf("Sum is %" PRIi64 "\n",sum);
  return 0;
}
