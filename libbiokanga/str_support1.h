// Str_Support1.h
//
// Version 2.4.0
//
// This file is a part of the Str Library, and is subject to the
// clauses of the License Agreement accompanying the software.

#ifdef STR_WIN32
#include <windows.h>
#if !defined(STR_NO_UNICODE) && !defined(STR_NO_WINSTUFF)
#include <wtypes.h>
#include <oleauto.h>
#endif
#endif

#ifdef STR_NO_RTLIBRARY
#ifndef STR_WIN32
#error Cannot define STR_NO_RTLIBRARY under Linux/gcc!
#endif
#endif

#ifndef STR_NO_RTLIBRARY
#include <stdlib.h>		// malloc/free
#include <stdio.h>		// For sprintf
#include <memory.h>		// For memcpy, memmove, etc
#include <stdarg.h>		// for FormatVA
#endif

class Str;

// CPOS is a character index or size (expressed in characters).  This one we
//   want outside of the StrImplementation namespace.  It MUST be signed,
//   otherwise multiple compiler conflicts with LPCWSTR will occur.
typedef int CPOS;


// Check platform declaration
#if !defined(STR_WIN32)  &&  !defined(STR_LINUX)
#error Please define either STR_WIN32 or STR_LINUX
#endif
#if defined(STR_WIN32)  &&  defined(STR_LINUX)
#error Please define only one of STR_WIN32 or STR_LINUX
#endif
#if defined(UNDER_CE)  &&  (UNDER_CE>0)
#define STR_WINDOWS_CE
#if defined(STR_LINUX)
#error STR_LINUX cannot be defined when compiling for Windows CE
#endif
#ifdef STR_NO_RTLIBRARY
#error STR_NO_RTLIBRARY cannot be defined when compiling for Windows CE
#endif
// STR_NO_WINSTUFF always turned on for EVC
#undef  STR_NO_WINSTUFF
#define STR_NO_WINSTUFF
#endif

// Borland C++ builder support
#if defined(__BORLANDC__)  &&  (__BORLANDC__ >= 0x550)
#ifdef STR_LINUX
#error STR_LINUX is incompatible with Borland C++; use STR_WIN32
#endif
#ifdef STR_NO_RTLIBRARY
#error STR_NO_RTLIBRARY is incompatible with Borland C++
#endif
#ifdef _UNICODE
#error _UNICODE is incompatible with Borland C++
#endif
#ifdef STR_UNICODE
#error STR_UNICODE is incompatible with Borland C++
#endif
#define STR_BORCPPBUILDER
#endif


#if defined(STR_DEBUG)  ||  defined(STR_NODEBUG)
//
// Host application has manually declared the library debug mode
#if defined(STR_DEBUG)  &&  defined(STR_NODEBUG)
#error Cannot define STR_DEBUG and STR_NODEBUG at the same time!
#endif
#else
//
// Determine debug mode.  As a first step, assume DEBUG means _DEBUG
#if defined(DEBUG) && !defined(_DEBUG)
#define _DEBUG
#endif
#if !defined(STR_WINDOWS_CE)  &&  !defined(STR_BORCPPBUILDER)
// Under Embedded VC, NDEBUG is always defined and we can't use it as a marker.
// With Borland C++ Builder, NDEBUG is not used - the presence or absence of
// _DEBUG tells us what to do
// Do some checks on other compilers, though (require one of _DEBUG / NDEBUG)
#if !defined(_DEBUG)  &&  !defined(NDEBUG)
#error Please define either _DEBUG or NDEBUG
#endif
#if defined(_DEBUG)  &&  defined(NDEBUG)
#error Please define either _DEBUG or NDEBUG, but not both
#endif
#endif
// Use STR_DEBUG or STR_NODEBUG as our own marker
#ifdef _DEBUG
#define STR_DEBUG
#else
#define STR_NODEBUG
#endif
#endif


// Unicode support?
#ifdef STR_NO_UNICODE
#ifdef STR_UNICODE
#error Cannot define STR_UNICODE when STR_NO_UNICODE is defined
#endif
#endif


// Get a few Linux typedef-s outside of the StrImplementation namespace
#ifdef STR_LINUX

#include <stdarg.h>	// Do NOT use varargs.h!!!
#include <ctype.h>	// toupper(), tolower()

typedef int BOOL;
typedef unsigned char BYTE;
typedef void* LPVOID;
typedef const void* LPCVOID;
typedef unsigned int UINT;
typedef unsigned int DWORD;

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

#endif

// No C RTL automatically implies No Exceptions
#ifdef STR_NO_RTLIBRARY
#ifndef STR_NO_EXCEPTIONS
#define STR_NO_EXCEPTIONS
#endif
#endif

// No Exceptions not compatible with any of the advanced functionality
#ifdef STR_NO_EXCEPTIONS
#ifdef STR_USE_REGEX
#error Regular expressions functionality cannot be used when STR_NO_EXCEPTIONS is turned on
#endif
#ifdef STR_USE_EXTRAS
#error StrArray and other extra features cannot be used when STR_NO_EXCEPTIONS is turned on
#endif
#endif

// Unicode/ANSI helpers
#if (!defined(STR_UNICODE))  &&  (!defined(STR_ANSI))
#ifdef _UNICODE
#define STR_UNICODE
#else
#define STR_ANSI
#endif
#else
#if defined(STR_UNICODE)  &&  defined(STR_ANSI)
#error Cannot define STR_UNICODE and STR_ANSI simultaneously
#endif
#endif
#if defined(STR_UNICODE)  &&  defined(STR_NO_UNICODE)
#error Cannot define STR_UNICODE and STR_NO_UNICODE at the same time
#endif

#if defined(STR_LINUX)  &&  !defined(_T)
#ifdef STR_UNICODE
#define _T(x) L##x
#else
#define _T(x) x
#endif
#endif

// Managed C++ extension features
#ifdef STR_MCPP_FEATURES

#ifdef STR_NO_RTLIBRARY
#error Managed C++ extensions in Str Library are not compatible with STR_NO_RTLIBRARY
#endif
#ifdef STR_NO_UNICODE
#error Managed C++ extensions in Str Library are not compatible with STR_NO_UNICODE
#endif

#using <mscorlib.dll>

#endif

// Unicode support in Linux
#if defined(STR_LINUX)

#include <wchar.h>		// wchar_t type, if not built in, and some wcsXXX functions
#include <stdlib.h>		// wcsXXX functions

// Note: Str Library requires that sizeof(short)==2 for Unicode support on this platform
// Note: Str Library requires that sizeof(wchar_t)==sizeof(short) for Unicode support on this platform

#endif


// DLL exported version?
#ifdef STR_EXPORT
#ifndef STR_WIN32
#error STR_EXPORT can be used only in Win32 projects!
#endif
#else
#define STR_EXPORT		/* No-op */
#endif


// Limit the number of free blocks kept in blocks pool.  Thanks to
//   Mr. Diener from OR Soft Janicke for contributing this mechanism.
#ifdef STR_MAX_FREE_CACHE
#ifdef STR_NO_BLOCK_CACHE
#error Cannot define STR_MAX_FREE_CACHE together with STR_NO_BLOCK_CACHE
#endif
#else
#define STR_MAX_FREE_CACHE 5
#endif


// STRCHAR and related types
#ifdef STR_LINUX
#ifndef STR_NO_UNICODE
typedef wchar_t STRWCHAR;
#endif
typedef char    STRACHAR;
#else
#ifndef STR_NO_UNICODE
typedef WCHAR	STRWCHAR;
#endif
typedef CHAR	STRACHAR;
#endif

#ifdef STR_UNICODE
typedef STRWCHAR STRCHAR;
#else
typedef STRACHAR STRCHAR;
#endif


// If calling LoadString and not using MFC, caller must implement STR_get_stringres.
#ifdef STR_WIN32
HANDLE STR_get_stringres();
#endif

// If using non-MFC or Linux app, caller must implement STR_get_debugname
const STRCHAR* STR_get_debugname();


// Additional sanity checks for regular expressions support
#if defined(STR_USE_REGEX)  &&  defined(STR_NO_RTLIBRARY)
#error Regular expressions cannot be used when STR_NO_RTLIBRARY is turned on
#endif


// Very basic platform-dependent typedefs and includes
#ifdef STR_WIN32

#include <windows.h>	// For HANDLE type
#include <tchar.h>

#else

#define STR_NO_WINSTUFF		/* Always required under Linux */

#ifdef STR_THREAD_SAFE
#error No support for thread safety exists under Linux.  Please do not define STR_THREAD_SAFE
#endif

#endif

#define STR_max(x,y)		\
	((x) > (y) ? (x) : (y))
#define STR_min(x,y)		\
	((x) < (y) ? (x) : (y))


namespace StrImplementation {

// MSVC compiler
#if defined(STR_WIN32)  &&  !defined(STR_BORCPPBUILDER)
#pragma intrinsic (memset)
#pragma warning (disable:4514 4127 4511 4512 4710)
#endif

// Borland C++ compiler
#if defined(STR_BORCPPBUILDER)
#pragma warn -8027
#endif

// GCC compiler
#ifdef STR_LINUX
#ifdef STR_GCC_295
#define memset(x,y,z) memset((x),(y),(z))
#else
#define memset(x,y,z) __builtin_memset((x),(y),(z))
#endif
#endif


struct verifymt;		// Internal -- implemented in .cpp
struct syncstr;			// Internal -- implemented in .cpp


// Internal malloc/free replacements
#if defined(STR_DEBUG)  &&  defined(STR_WIN32)  &&  !defined(STR_BORCPPBUILDER)

#define STR_malloc(x)	\
	StrImplementation::STR_malloc_dbg((x), __LINE__, __FILE__)
void* STR_malloc_dbg(size_t size, int line, const char* file);

#else

void* STR_malloc(size_t size);

#endif

void STR_free(void *memblock);


// Check for structure endian-ness
#if defined(STR_BIG_ENDIAN)  &&  defined(STR_LITTLE_ENDIAN)
#error Cannot define STR_BIG_ENDIAN and STR_LITTLE_ENDIAN at the same time
#endif

#ifdef STR_WIN32
#ifdef STR_BIG_ENDIAN
#error On WIN32 the architecture is always little endian!  Please undefine STR_BIG_ENDIAN
#endif
#ifndef STR_LITTLE_ENDIAN
#define STR_LITTLE_ENDIAN
#endif
#else	// ifdef STR_WIN32
#if !defined(STR_BIG_ENDIAN)  &&  !defined(STR_LITTLE_ENDIAN)
#error You must define STR_BIG_ENDIAN or STR_LITTLE_ENDIAN, depending on your machine CPU!
#endif
#endif


// Check for name duplication of two important macros
#ifdef null_b
#error The macro null_b is used by Str Library and should not be defined here
#endif
#ifdef null_str
#error The macro null_str is used by Str Library and should not be defined here
#endif
#ifdef null_str_m
#error The macro null_str_m is used by Str Library and should not be defined here
#endif

#ifndef STR_STR_FACTOR
#define STR_STR_FACTOR   5
#else
#if (STR_STR_FACTOR < 3)
#error The absolute minimum value for STR_STR_FACTOR is 3
#endif
#endif
#define ANSI_BAE_MASK	((1 << STR_STR_FACTOR) - 1)
#define WIDE_BAE_MASK	((1 << (STR_STR_FACTOR - 1)) - 1)


// Check if a bytesize value is a multiple of the BAE factor
inline BOOL bae_aligned(int size)
{ return !(size & ((1 << STR_STR_FACTOR)-1)); }


// The structure below is implementation specific, do not
//   attempt to use it directly!
#pragma pack(4)
typedef struct SBlock {
	// Data fields
#ifdef STR_THREAD_SAFE
	long	m_Locked;	// Instance locked by another thread
	DWORD   m_OwnerT;	// Owner thread
#endif
	long	m_Ref;		// Usage counter, 0 in empty slots, always 2 in null str
	DWORD	m_Flags;	// SFLAG_ combination
	CPOS	m_Alloc;	// Allocated length, excl. 0
	CPOS	m_Length;	// Actual string length, excl. 0
	// 0-terminated string starts here

// Flag management
	BOOL    GetFlag(DWORD flagmask)     { return (m_Flags & flagmask) == flagmask; }
	void    SetFlag(DWORD flagmask);		// Inlined later
	void    ClearFlag(DWORD flagmask);		// Inlined later

} SBlock;
#pragma pack()

#define SFLAG_MT		1		/* Marked as multithreaded */

// Wrapper around SBlock, used for the single shared empty string instance
#pragma pack(1)
typedef struct SBlock_Empty {
	SBlock   data;
	STRACHAR emptys;
} SBlock_Empty;
#pragma pack()


/*********************************************************************
* verifymt is a helper class instantiated in the beginning of all
*	methods that _might_ call syncstr.  In release and non-mt builds,
*   it does nothing.  Otherwise, it verifies that the operation is
*   not trying to access a string object that is not MT marked and
*   belongs to a different thread.
*********************************************************************/

struct verifymt
{
	Str* the_object;
	verifymt(Str* obj);
	verifymt(Str* obj, int dummy);
};


// STR_DITEMS: biggest item to be kept in blocks pool
#ifndef STR_DITEMS
#if STR_STR_FACTOR > 5
#define STR_DITEMS ((1 << STR_STR_FACTOR) * 10)
#else
#define STR_DITEMS 320
#endif
#endif


/*********************************************************************
* Proc:		STR_Inc and STR_Dec
* Purpose:	Add or subtract one from the given variable in a thread-
*			safe way
*********************************************************************/

#ifdef STR_THREAD_SAFE
#define STR_Inc(val) InterlockedIncrement(&(val))
#define STR_Dec(val) InterlockedDecrement(&(val))
#else
#define STR_Inc(val) ++(val)		/* Thread unsafe version */
#define STR_Dec(val) --(val)		/* Thread unsafe version */
#endif


/*********************************************************************
* Debugging facilities
*********************************************************************/

#ifdef STR_DEBUG

void doAssert(int line, const char* file);

#if defined(STR_WIN32)  && !defined(STR_BORCPPBUILDER)  &&  !defined(STR_WINDOWS_CE)
#define BKPT DebugBreak()
#else
#define BKPT ((void) 0)
#endif

#ifndef STR_ASSERT
#define STR_ASSERT(x)		\
	if ((x) == 0) { StrImplementation::doAssert(__LINE__, __FILE__); 	BKPT; } \
	else { (void) 0; }
#endif
#define STR_ASSERT_FALSE()

#else

#ifndef STR_ASSERT
#define STR_ASSERT(x)		/* no-op */
#define STR_ASSERT_FALSE()	/* no-op */
#endif

#endif


STR_EXPORT void STR_rls_block_realrelease(SBlock* block);
STR_EXPORT int  STR_ReplaceCore(Str& obj, const STRCHAR *oldstr, const STRCHAR *newstr, BOOL nocase);


// Works only for positive numbers; factor must be power of 2
inline int up_factor (int value, int factor)
{
	int cutdown = (int) ((value+(factor-1)) / factor);
	return cutdown * factor;
}

}

