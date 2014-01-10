#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#include <mex.h>
#include "utils.h"
extern int interrupted;

static const  int     n_timers = 8;
static struct timeval tstart[8];
static int error(const char *, ...);

/*---------------------------------------------------------------------------*/
void tic (int timernum)
/*---------------------------------------------------------------------------*/
{
  if (timernum > n_timers)
  {
    error ("Timer number %d outside range\n", timernum);
    exit(1);
  }
  gettimeofday(&tstart[timernum], NULL);
}

/*---------------------------------------------------------------------------*/
double toc (int timernum)
/*---------------------------------------------------------------------------*/
{
  if (timernum > n_timers)
  {
    error ("Timer number %d outside range\n", timernum);
    exit(1);
  }
  struct timeval tend;
  gettimeofday(&tend, NULL);

  double seconds =
    difftime(tend.tv_sec, tstart[timernum].tv_sec)
    +1.0e-6*(tend.tv_usec-tstart[timernum].tv_usec);

  return seconds;
}

/*---------------------------------------------------------------------------*/
static int error (const char *format, ...)
/*---------------------------------------------------------------------------*/
{
  int retval = 0;
  if (1) /* Always show error messages */
  {
    const int size = 1024;
    char *message = mxMalloc(size*sizeof(char));
    int n;

    va_list ap;
    va_start(ap, format);
    n  = snprintf(message, size, "Error:\n");
    n += vsnprintf (message+n, size, format,  ap);
    va_end(ap);

    mexErrMsgTxt(message);
    mxFree(message);
  }
  return retval;
}
