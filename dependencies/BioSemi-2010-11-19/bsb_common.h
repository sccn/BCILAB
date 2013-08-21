#include "mex.h"

/* library handling */
#ifdef _WIN32
    #include <windows.h>
    #define LOAD_LIBRARY(lpath) LoadLibrary(lpath)
    #define LOAD_FUNCTION(lhandle,fname) GetProcAddress(lhandle,fname)
    #define FREE_LIBRARY(lhandle) FreeLibrary(lhandle)
#else
    #include <dlfcn.h>
    #define LOAD_LIBRARY(lpath) dlopen(lpath,RTLD_NOW | RTLD_LOCAL)
    #define LOAD_FUNCTION(lhandle,fname) dlsym(lhandle,fname)
    #define FREE_LIBRARY(lhandle) dlclose(lhandle)
#endif

/* 32/64 bit handling */
#ifdef __LP64__
    #define uintptr_t uint64_T
    #define PTR_CLASS mxUINT64_CLASS
#else
    #ifdef _WIN64
        #define uintptr_t uint64_T
        #define PTR_CLASS mxUINT64_CLASS
    #else
        #define uintptr_t uint32_T
        #define PTR_CLASS mxUINT32_CLASS
    #endif
#endif
        
/* length of the buffer in samples */
#define BUFFER_SAMPLES (60*512)

#ifdef _WIN32
    #define LINKAGE __cdecl
#else
    #define LINKAGE
#endif

/* function handle types */
typedef void* (LINKAGE *OPEN_DRIVER_ASYNC_t)(void);
typedef int (LINKAGE *USB_WRITE_t)(void*, char*);
typedef int (LINKAGE *READ_MULTIPLE_SWEEPS_t)(void*,char*,int);
typedef int (LINKAGE *READ_POINTER_t)(void*,unsigned*);
typedef int (LINKAGE *CLOSE_DRIVER_ASYNC_t)(void*);
