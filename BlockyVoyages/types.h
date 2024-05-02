/*
\file          basictypes.h
\author        Dwayne Robinson
\since         2005-04-08
\date          2005-09-18

Very basic data types that ensure field size consistency between compilers
and languages. Consistent size is necessary if using share libraries or
mixing with other languages like assembly.

Also declares several little constants that should have been added to the
core language long ago. Putting these in here reduces dependencies on
system specific header files.
*/

#pragma once

#pragma warning (disable : 4996)

// essential primitives every version of C should have had from day 1
#if defined(_MSC_VER)
    typedef          unsigned int                uint;
    typedef          unsigned __int8            byte;
    typedef             signed __int8            sbyte;
    typedef             signed __int8            int8;
    typedef             signed __int16          int16;
    typedef             signed __int32          int32;
    typedef             signed __int64          int64;
    typedef          unsigned __int8            ubyte;
    typedef          unsigned __int8            uint8;
    typedef          unsigned __int16          uint16;
    typedef          unsigned __int32          uint32;
    typedef          unsigned __int64          uint64;
#elif defined(__MWERKS__)
    typedef          unsigned int                uint;
    typedef          unsigned int8_t            byte;
    typedef             signed int8_t            sbyte;
    typedef             signed int8_t            int8;
    typedef             signed int16_t          int16;
    typedef             signed int32_t          int32;
    typedef             signed int64_t          int64;
    typedef          unsigned int8_t            ubyte;
    typedef          unsigned int8_t            uint8;
    typedef          unsigned int16_t          uint16;
    typedef          unsigned int32_t          uint32;
    typedef          unsigned int64_t          uint64;
#elif defined(__GNUC__)
    typedef          unsigned int                uint;
    typedef          unsigned char              byte;
    typedef             signed char              sbyte;
    typedef             signed char              int8;
    typedef             signed short             int16;
    typedef             signed int                int32;
    typedef             signed long              int64;
    typedef          unsigned char              ubyte;
    typedef          unsigned char              uint8;
    typedef          unsigned short             uint16;
    typedef          unsigned int                uint32;
    typedef          unsigned long              uint64;
#endif

// added by Gregory Esch 2007-3-2
typedef                      float                  float32;
typedef                      double                 float64;

// These two should also have been defined from the very beginning of C.
#if !defined(__cplusplus_cli) && !defined(__cplusplus)
    #ifndef true
        #define true  1
        #define false 0
    #endif
#endif

// Surprisingly, this essential primitive is sometimes undefined on Unix systems.
// Retardly, Visual Studio 2005 highlights the reserved name but only compiles
// correctly if the CLR is included. Bjarne Stroustrup advocates it in C++.
#if !defined(__cplusplus_cli)
    #if defined(__GNUC__)
        #define nullptr __null
    #elif !defined(__cplusplus)
        #define nullptr ((void*)0)
    #endif
#endif

#define ElemCount(elementlist, T)  sizeof(elementlist) / sizeof(T)