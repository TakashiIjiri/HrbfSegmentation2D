#ifndef __CLAPACK_HEAD__
#define __CLAPACK_HEAD__


extern "C"
{
#include "./clapack/f2c.h"
#include "./clapack/clapack.h"
}

#ifdef _WIN64

#ifdef _DEBUG

#pragma comment(lib, "BLASd.lib")
#pragma comment(lib, "libf2cd.lib")
#pragma comment(lib, "clapackd.lib")

#else
#pragma comment(lib, "BLAS.lib")
#pragma comment(lib, "libf2c.lib")
#pragma comment(lib, "clapack.lib")

#endif

#else

#ifdef _DEBUG

#pragma comment(lib, "BLASd.lib")
#pragma comment(lib, "libf2cd.lib")
#pragma comment(lib, "clapackd.lib")

#else

#pragma comment(lib, "BLAS.lib")
#pragma comment(lib, "libf2c.lib")
#pragma comment(lib, "clapack.lib")

#endif

#endif

//extern "C" long _ftol( double ); //defined by VC6 C libs
//extern "C" long _ftol2( double dblSource ) { return _ftol( dblSource ); }

#endif //__UMFPACK_HEAD__ 

