#ifndef _SystemHeaders_h_
#define _SystemHeaders_h_


#if(1)

//////////////////
// UNIX section //
//////////////////

#include <unistd.h>


#else

//////////////////
// VC++ section //
//////////////////

#include <math.h>
static double rint(double x)
{
    return floor(x);
}


#define M_PI ((double)3.1415926535897932384626433832795)

#include <windows.h>
#include <winbase.h>
static void sleep(unsigned int seconds)
{
    Sleep(1000 * seconds);
}


#include <float.h>
static int isnan(double x)
{
    return _isnan(x);
}


static void srand48(long seedval)
{
    //srand(seedval);
}
#include "rank.h"
#include "IntRandom.h" 
static double drand48(void)
{
    return ((double)rand_r(&rank2))/((double)RAND_MAX+1);;
}

#endif


#endif
