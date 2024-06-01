#ifndef __UMFPACK_HEAD__
#define __UMFPACK_HEAD__

// This file should be used in the context of extern "C"{}.
//extern "C"
//{
//#include "cblas.h"
//}

extern "C"
{
#include "./umfpack/umfpack.h"
}

#ifdef _WIN64

#pragma comment(lib, "libamd_woBlas64.lib")
#pragma comment(lib, "libumfpack_woBlas64.lib")

#else

#pragma comment(lib, "libamd_woBlas32.lib")
#pragma comment(lib, "libumfpack_woBlas32.lib")

#endif

//extern "C" long _ftol( double ); //defined by VC6 C libs
//extern "C" long _ftol2( double dblSource ) { return _ftol( dblSource ); }

#endif //__UMFPACK_HEAD__ 

