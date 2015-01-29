// Str_Implement.cpp
//
// Version 2.4.0
//
// This file is a part of the Str Library, and is subject to the
// clauses of the License Agreement accompanying the software.

//!!stdafx

#include "str.h"

#ifdef STR_BORCPPBUILDER
#define Char StrChar
#pragma warn -8004
#pragma warn -8008
#pragma warn -8066
#endif

#ifdef STR_MCPP_FEATURES
#pragma warning (disable:4189)
#include <vcclr.h>
#pragma warning (default:4189)
#endif
#ifdef STR_WIN32
#pragma warning (disable:4702 4710)
#endif

using namespace StrImplementation;

struct StrMTType {
	int dummy;
};
StrMTType gStrMT = {1};
const StrMTType& StrMT = gStrMT;


#if defined(__AFX_H__)  &&  !defined(STR_NO_WINSTUFF)
HANDLE STR_get_stringres()		{ return AfxGetResourceHandle(); }
#endif


/*********************************************************************
* Global variables
*********************************************************************/

// Shutdown flag, effectively disables keep-free-mem mechanism
namespace StrImplementation {
BOOL STR_Shutting = FALSE;
}


/*********************************************************************
* Static instance of null string (single-threaded & multi-threaded)
*   and other on-startup initializations
*********************************************************************/

namespace StrImplementation {

STR_Initializer STR_EXPORT oInitializer;

#ifndef STR_NO_RTLIBRARY
STR_Initializer::STR_Initializer()
#else
void STR_Initializer::InitializeStr()
#endif
{
	memset (&m_Empty, 0, sizeof(m_Empty));
	m_Empty.data.m_Ref = 2;
	m_Empty.data.m_Flags = 0;
	m_Empty.data.m_Alloc = ANSI_BAE_MASK;
}

// Instantiate empty string helper instance
#ifndef STR_NO_RTLIBRARY
Str STR_EXPORT gEmptyStr;
#endif

}

#ifndef STR_NO_RTLIBRARY
const Str& Str::EmptyStr = StrImplementation::gEmptyStr;
#endif

#ifdef STR_NO_RTLIBRARY
#define REFER_EmptyStr _T("")
#else
#define REFER_EmptyStr Str::EmptyStr
#endif


/*********************************************************************
* Proc:		Str::CompactFreeMemory, Str::SetShutdownFlag
* Purpose:	Releases any unused blocks we might be caching at the moment.
*********************************************************************/

/*static*/ void Str::CompactFreeMemory()
{
#ifndef STR_NO_BLOCK_CACHE
	for (UINT i=0; i <= (STR_DITEMS / (1 << STR_STR_FACTOR)); i++)
		STR_rls_chain(i);
#endif
}

/*static*/ void Str::SetShutdownFlag(BOOL shutdown /*= TRUE*/)
{
	if (shutdown) {
		CompactFreeMemory();
		STR_Shutting = TRUE;
	}
	else
		STR_Shutting = FALSE;
}


/*********************************************************************
* Proc:		doAssert
* Purpose:	Handler for STR_ASSERT() macro in debug builds
*********************************************************************/

#ifdef __AFX_H__
const STRCHAR* STR_get_debugname()
{
	return AfxGetAppName();
}
#endif

#ifdef STR_DEBUG

void StrImplementation::doAssert(int line, const char* file)
{
	if (file == NULL)
		file = "?";

#ifdef STR_UNICODE
	const STRACHAR* f = "[%S] assertion failed!\n\nFile:\t%s\nLine:\t%d\n";
#else
	const STRACHAR* f = "[%s] assertion failed!\n\nFile:\t%s\nLine:\t%d\n";
#endif

	Str temp;
	temp.Format(f, STR_get_debugname(), file, line);
#if defined(STR_WIN32)  &&  !defined(STR_NO_WINSTUFF)
	::MessageBox (NULL, temp, STR_get_debugname(), MB_SYSTEMMODAL | MB_ICONSTOP);
#else
	fprintf (stderr, temp.GetStringA());	// Will always leak
#endif
	STR_abort();
}

#endif


/*********************************************************************
* Proc:		Str::Error, ErrorNoObject and supporting functions
*********************************************************************/

#ifndef STR_USER_ERRORS

#if defined(STR_NO_EXCEPTIONS)
static void STR_abort_application(int fatal_type)
{
	// Display message and abort
	Str msg;
	msg.Format(_T("Str fatal %d"), fatal_type);
#ifdef STR_UNICODE
	STRACHAR* msg1 = msg.GetStringA();		// Will produce memory leak
#else
	const STRACHAR* msg1 = msg;
#endif
#if defined(STR_WIN32)  &&  !defined(STR_NO_WINSTUFF)
	::MessageBoxA (NULL, msg1, NULL, MB_SYSTEMMODAL | MB_ICONSTOP);
#else
	fprintf (stderr, msg1);
#endif
	STR_abort();
}
#endif

#if !defined(STR_NO_EXCEPTIONS)  &&  defined(STR_DEBUG)  &&  defined(STR_WIN32)  &&  !defined(STR_NO_WINSTUFF)
inline void LogThrowingException(StrException* eobject)
{
	STRACHAR buf[128];
	wsprintfA(buf, "Throwing StrException (obj=%p code=%d)\n", &(eobject->m_Object), eobject->m_Error);
	OutputDebugStringA(buf);
	eobject;
}
#endif


void Str::Error(StrEC ecode) const
{
#ifdef STR_NO_EXCEPTIONS
	STR_abort_application((int) ecode);
#else
	StrException eobject;
	eobject.m_Error = ecode;
	eobject.m_Object = *this;
#if defined(STR_DEBUG)  &&  defined(STR_WIN32)  &&  !defined(STR_NO_WINSTUFF)
	LogThrowingException(&eobject);
#endif
	throw eobject;
#endif
}

/*static*/ void Str::ErrorNoObject(StrEC ecode)
{
#ifdef STR_NO_EXCEPTIONS
	STR_abort_application((int) ecode);
#else
	StrException eobject;
	eobject.m_Error = ecode;
#if defined(STR_DEBUG)  &&  defined(STR_WIN32)  &&  !defined(STR_NO_WINSTUFF)
	LogThrowingException(&eobject);
#endif
	throw eobject;
#endif
}

#endif


/*********************************************************************
* Proc:		Str_ValidatePtr
* Purpose:	Make sure the value given looks like a ptr to R/W memory
*********************************************************************/

#ifdef STR_DEBUG

void Str_ValidatePtr(const void* x, BOOL ornull);	// CodeWarrior fix

void Str_ValidatePtr(const void* x, BOOL ornull)
{
	if (ornull  &&  x == NULL)
		return;
#ifdef STR_WIN32
	if (IsBadReadPtr(x, 1) || IsBadWritePtr((void*) x, 1))
		DebugBreak();
#else
	/* AWAITING IMPLEMENTATION! */
#endif
}
#endif


/*********************************************************************
* Proc:		STR_blkalloc / STR_blkverify / STR_blkrealfree
* Purpose:	Helpers
*********************************************************************/

#define PRESERVE_malloc malloc
#define PRESERVE_free   free
#undef  malloc
#undef  free

#ifdef STR_DEBUG

#ifdef _CRTDBG_MAP_ALLOC
#error _CRTDBG_MAP_ALLOC cannot be defined at this point
#endif

#if defined(STR_WIN32)  &&  !defined(STR_NO_RTLIBRARY)

#if !defined(STR_BORCPPBUILDER)  &&  !defined(STR_WINDOWS_CE)

// Win32, regular C RTL
#include <crtdbg.h>
#define DECLARE_STR_blkalloc	\
	BYTE* STR_blkalloc (UINT bytes, int line, const char* file);  /* CodeWarrior fix */	\
	BYTE* STR_blkalloc (UINT bytes, int line, const char* file)

#else

// Win32, Borland C++ Builder or Embedded Visual C++
#define _malloc_dbg(a,b,c,d) malloc(a)
#define _free_dbg(a,b)       free(a)
#define DECLARE_STR_blkalloc	\
	BYTE* STR_blkalloc (UINT bytes, int /*line*/, const char* /*file*/)

#endif

#else

#define _malloc_dbg(a,b,c,d) Str::OS_malloc((a))
#define _free_dbg(a,b)       Str::OS_free(a)

#define DECLARE_STR_blkalloc	\
	BYTE* STR_blkalloc (UINT bytes, int /*line*/, const char* /*file*/);  /* CodeWarrior fix */ \
	BYTE* STR_blkalloc (UINT bytes, int /*line*/, const char* /*file*/)

#endif

// Debug alloc: prepend byte length, append protective mask
DECLARE_STR_blkalloc
{
	BYTE* x = (BYTE*) _malloc_dbg(bytes+8, _NORMAL_BLOCK, file, line);
	((UINT*)x)[0] = bytes;
	x += 4;
	memset (x+bytes, 0x78, 4);
	return x;
}

// Debug verify: make sure byte length matches and protective mask is OK

void STR_blkverify (BYTE* block, UINT bytes);		// CodeWarrior fix

void STR_blkverify (BYTE* block, UINT bytes)
{
	BYTE* realblock = block-4;
	STR_ASSERT (((UINT*)realblock)[0] == bytes);
	STR_ASSERT (block[bytes+0] == 0x78  &&  block[bytes+1] == 0x78);
	STR_ASSERT (block[bytes+2] == 0x78  &&  block[bytes+3] == 0x78);
}

// Debug real free: call free() 4 bytes earlier

void STR_blkrealfree(BYTE* block, UINT bytes);

void STR_blkrealfree(BYTE* block, UINT bytes)
{
	if (block == NULL)
		return;
	if (bytes != (UINT) -1)
		STR_blkverify(block, bytes);
	BYTE* realblock = block-4;
	_free_dbg(realblock, _NORMAL_BLOCK);
}

#else

// Release alloc: just do malloc
inline BYTE* STR_blkalloc (UINT bytes)
{
	return (BYTE*) Str::OS_malloc (bytes);
}

// Release verify: do nothing
inline void STR_blkverify (BYTE*, UINT)
{ }

// Release real free: just call free()
inline void STR_blkrealfree(BYTE* block, UINT)
{
	Str::OS_free(block);
}

#endif


/*********************************************************************
* Proc:		STR_malloc
*********************************************************************/

#if defined(STR_DEBUG)  &&  defined(STR_WIN32)  &&  !defined(STR_BORCPPBUILDER)
void* StrImplementation::STR_malloc_dbg(size_t siz, int line, const char* file)
{
#if defined(STR_NO_RTLIBRARY)  ||  defined(STR_WINDOWS_CE)  ||  defined(STR_BORCPPBUILDER)
	line; file;
#endif
#ifdef STR_NO_RTLIBRARY
	return Str::OS_malloc(siz);
#else
	return _malloc_dbg(siz, _NORMAL_BLOCK, file, line);
#endif
}
#else
void* StrImplementation::STR_malloc(size_t siz)
{
	return Str::OS_malloc(siz);
}
#endif


/*********************************************************************
* Proc:		STR_free
*********************************************************************/

void StrImplementation::STR_free(void *memblock)
{
	Str::OS_free(memblock);
}


/////////////////////////////////////
// Start of Str block repository code

namespace StrImplementation {
inline SBlock* blocks_default(SBlock* obtained, CPOS numchars)
{
	obtained->m_Ref    = 1;
	obtained->m_Flags  = 0;
#ifdef STR_THREAD_SAFE
	obtained->m_Locked = 0;
	obtained->m_OwnerT = STR_GetCurrentThreadID();
#endif
	obtained->m_Alloc  = numchars;
	return obtained;
}
}


#ifndef STR_NO_BLOCK_CACHE

/*********************************************************************
* Free block repository support (single-threaded config).  The array
*    contains a linked list of available memory blocks of exactly
*    the required size (8, 16, 24, etc.)
*********************************************************************/

#define STR_REPOS_ASIZE ((STR_DITEMS / (1 << STR_STR_FACTOR))+1)

// The following two arrays are automatically initialized to zero bytes
//   as per C/C++ language specification.  Thus we avoid an on-startup-
//   static-constructor to fill them in with zeros.
static BYTE* Str_Repos[STR_REPOS_ASIZE]; // = { 0 }
static int   Str_ReposCnt[STR_REPOS_ASIZE]; // = { 0 }

#ifndef STR_THREAD_SAFE

namespace StrImplementation {
struct cache_arb
{
	cache_arb(int)			{ }
	~cache_arb()			{ }
};
}

#else

/*********************************************************************
* Helper class structure -- faster (and more efficient) CRIT_SECT,
*   but WITHOUT reentrancy support!
*********************************************************************/

namespace StrImplementation {
typedef struct {
	long lock;
	void Enter();
	void Leave();
protected:
	void GetStructAccess();
	void RlsStructAccess();
} fastlock;
}

inline void fastlock::Enter()
{
	for (;;) {
		if (STR_Inc(lock) == 1)
			return;
		STR_Dec(lock);
		Sleep(0);
	}
}

inline void fastlock::Leave()
{
	STR_Dec(lock);
}


/*********************************************************************
* Free block repository support (multithreaded config).  The
*   additional array contains a number of Win32 critical section
*   objects used to protect individual items within the linked list.
*********************************************************************/

namespace StrImplementation {
struct cache_arb
{
	int m_cindex;
	static fastlock ReoProt[];
	cache_arb(int cindex)
		{
			m_cindex = cindex;
			ReoProt[cindex].Enter();
		}
	~cache_arb()	
		{
			ReoProt[m_cindex].Leave(); 
		}
};
}

fastlock cache_arb::ReoProt[STR_REPOS_ASIZE];

#endif


/*********************************************************************
* Proc:		get_block
* Purpose:	Helper.  Allocates (from the cache or memory) a block
*				with a specified amount of chars for a string.
*				Also fills the Ref, Locked and Alloc fields
*********************************************************************/

#ifdef STR_DEBUG
SBlock* StrImplementation::STR_get_block_dbg(CPOS numchars, int line, const char* file)
#else
SBlock* StrImplementation::STR_get_block_rls(CPOS numchars)
#endif
{
	// Make sure we don't get invalid charcounts
	STR_ASSERT_UPPED(numchars);
	// Compute number of bytes; check for large blocks
	int bytes = sizeof(SBlock) + ((numchars+1) * sizeof(STRCHAR));
	SBlock* obtained;
	if (bytes >= STR_DITEMS)
		goto DoObtain;
	// Check for block from cache
	{
		UINT cindex = (bytes-sizeof(SBlock)) / (1 << STR_STR_FACTOR);
		{
			cache_arb keeper (cindex);
			BYTE* x = Str_Repos[cindex];
			if (x == NULL)
				goto DoObtain;
			// Get the block at the head
			BYTE** x_fp = (BYTE**) x;
			Str_Repos[cindex] = *x_fp;
			STR_ASSERT(Str_ReposCnt[cindex] > 0);
			--Str_ReposCnt[cindex];
			obtained = (SBlock*) x;
			goto DoReturn;
		}
	}
DoObtain:
	{ }
#ifdef STR_DEBUG
	obtained = (SBlock*) STR_blkalloc (bytes, line, file);
#else
	obtained = (SBlock*) STR_blkalloc (bytes);
#endif
	// Put some simple defaults and return to caller
DoReturn:
	return blocks_default(obtained, numchars);
}


/*********************************************************************
* Proc:		STR_rls_block / STR_rls_block_realrelease
* Purpose:	Decrement the usage count of the passed block, and when
*			the count reaches zero, frees the block.
*********************************************************************/

void StrImplementation::STR_rls_block_realrelease(SBlock* block)
{
	// Release it
	UINT bytes = sizeof(SBlock) + ((block->m_Alloc+1) * sizeof(STRCHAR));
	STR_ASSERT((bytes & 3) == 0);
	STR_blkverify ((BYTE*) block, bytes);
	// Large blocks are returned directly to the memory manager.
	//   This is also done automatically for all blocks if
	//   the "shutting down" flag is TRUE
	if (bytes >= STR_DITEMS  ||  STR_Shutting) {
		STR_blkrealfree ((BYTE*) block, bytes);
		return;
	}
	// Debug only: set block data to garbage to detect past-lifetime
	//   usage
#ifdef STR_DEBUG
	memset (block, 0xCC, bytes);
#endif
	// Put reference to next free block in first 4 bytes
	UINT cindex = (bytes-sizeof(SBlock)) / (1 << STR_STR_FACTOR);
	BYTE** block_fp = (BYTE**) block;
	BOOL do_real_free = TRUE;
	{
		cache_arb keeper (cindex);
		// And add block as head of cache, if appropriate (assume linear distribution
		//   of maximum pools size, i.e. bigger blocks less often needed)
		const UINT cindex_MAX = STR_DITEMS / (1 << STR_STR_FACTOR);
		int per_cindex_maxcnt = (int)(cindex_MAX - cindex) * STR_MAX_FREE_CACHE;
		if (Str_ReposCnt[cindex] < per_cindex_maxcnt) {
			*block_fp = Str_Repos[cindex];
			Str_Repos[cindex] = (BYTE*) block;		// Add at head of cache
			++Str_ReposCnt[cindex];
			do_real_free = FALSE;
		}
	}
	if (do_real_free)
		STR_blkrealfree ((BYTE*) block, bytes);
}


#ifndef STR_NO_BLOCK_CACHE

void StrImplementation::STR_rls_block(verifymt& mtcheck)
{
	STRCHAR* blk_s = Str::block_s(mtcheck.the_object->block());
	if (blk_s == null_str())
		return;
	SBlock* block = (SBlock*) (((BYTE*) blk_s) - sizeof(SBlock));
#if defined(STR_THREAD_SAFE)  &&  defined(STR_DEBUG)
	if (!block->GetFlag(SFLAG_MT)) {
		STR_ASSERT(block->m_Locked == 0);
	}
#endif
	// Decrement usage counter
	long x = --(block->m_Ref);
	if (x)
		return;
	STR_rls_block_realrelease(block);
}

#endif


/*********************************************************************
* Proc:		STR_rls_chain
* Purpose:	Helper for CompactFree(); releases the entire chain of
*			free buffers with a specified size
*********************************************************************/

void StrImplementation::STR_rls_chain(UINT cindex)
{
	// Make a fast check.  We're allowed to do this OUTSIDE the
	//   arbitration code -- it may fail very, very rarely, and
	//   even in those cases it would be completely harmless.
	//   It's not advisable to combine the check with the ASSERT
	//   at the end of the method, though :)
	if (Str_Repos[cindex] == NULL)
		return;
	// Lock this chain
	cache_arb keeper (cindex);
	while (Str_Repos[cindex] != NULL) {
		BYTE** x = (BYTE**) Str_Repos[cindex];
		Str_Repos[cindex] = *x;
		--Str_ReposCnt[cindex];
		int to_be_sized = cindex * (1 << STR_STR_FACTOR) + sizeof(SBlock);
		STR_blkrealfree ((BYTE*) x, to_be_sized);
	}
	STR_ASSERT(Str_ReposCnt[cindex] == 0);
}

#else	// !defined(STR_NO_BLOCK_CACHE)

/*********************************************************************
* Proc:		get_block
* Purpose:	Version used when the block-repository mechanism is not
*			active.
*********************************************************************/


#ifdef STR_DEBUG
SBlock* StrImplementation::STR_get_block_dbg(CPOS numchars, int line, const char* file)
#else
SBlock* StrImplementation::STR_get_block_rls(CPOS numchars)
#endif
{
	// Make sure we don't get invalid charcounts
	STR_ASSERT_UPPED(numchars);

	// Compute number of bytes
	int bytes = sizeof(SBlock) + ((numchars+1) * sizeof(STRCHAR));
#ifdef STR_DEBUG
	SBlock* obtained = (SBlock*) STR_blkalloc (bytes, line, file);
#else
	SBlock* obtained = (SBlock*) STR_blkalloc (bytes);
#endif
	// Put some simple defaults and return to caller
	return blocks_default(obtained, numchars);
}


/*********************************************************************
* Proc:		rls_block
* Purpose:	Version used when the block repository is turned off
*********************************************************************/

void StrImplementation::STR_rls_block(verifymt& mtcheck)
{
	STRCHAR* blk_s = Str::block_s(mtcheck.the_object->block());
	if (blk_s == null_str())
		return;

	SBlock* block = (SBlock*) (((BYTE*) blk_s) - sizeof(SBlock));
#if defined(STR_THREAD_SAFE)  &&  defined(STR_DEBUG)
	if (!block->GetFlag(SFLAG_MT)) {
		STR_ASSERT(block->m_Locked == 0);
	}
#endif

	// Decrement usage counter
	long x = --(block->m_Ref);
	if (x)
		return;

	// Release it
	UINT bytes = sizeof(SBlock) + ((block->m_Alloc+1) * sizeof(STRCHAR));
	STR_ASSERT((bytes & 3) == 0);
	STR_blkverify ((BYTE*) block, bytes);

	// Debug only: set block data to garbage to detect past-lifetime
	//   usage

#ifdef STR_DEBUG
	memset (block, 0xCC, bytes);
#endif

	// Return to memory manager
	STR_blkrealfree ((BYTE*) block, bytes);
}


#endif	// !defined(STR_NO_BLOCK_CACHE)


/*********************************************************************
* Proc:	STR_strtol, STR_strtod
*********************************************************************/

#if defined(STR_NO_RTLIBRARY)  ||  defined(STR_LINUX)

long StrImplementation::STR_strtol(const STRCHAR *nptr, STRCHAR **endptr, int ibase)
{
	//#define FL_UNSIGNED   1
	const int FL_NEG        = 2;      // negative sign found
	//const int FL_OVERFLOW   = 4;      // overflow occured
	const int FL_READDIGIT  = 8;      // we've read at least one correct digit

	const STRCHAR *p = nptr;        // p is our scanning pointer
	unsigned long number = 0;       // start with zero
	int flags = 0;

	STRCHAR c = *p++;               // read char
	while (STR_CharIsSpace(c))
		c = *p++;               // skip whitespace
    if (c == '-') {
		flags |= FL_NEG;        // remember minus sign
		c = *p++;
	}
	else if (c == '+')
		c = *p++;               // skip sign
	if (ibase < 0 || ibase == 1 || ibase > 36) {
		if (endptr)
			*endptr = (STRCHAR*) nptr;
		return 0;
	}
    else if (ibase == 0) {
		// determine base free-lance, based on first two chars of
		// string
	if (c != '0')
		ibase = 10;
	else if (*p == 'x' || *p == 'X')
		ibase = 16;
	else
		ibase = 8;
    }
	if (ibase == 16) {
		// we might have 0x in front of number; remove if there
		if (c == '0' && (*p == 'x' || *p == 'X')) {
			++p;
			c = *p++;
		}
	}
	for (;;) {
		// convert c to value
		unsigned int digval;
		if (STR_CharIsDigit(c))
			digval = c - '0';
		else if (STR_CharIsAlpha(c))
			digval = STR_CharToUpper(c) - _T('A') + 10;
		else
			break;
		if (digval >= (unsigned)ibase)
			break;          // exit loop if bad digit found
		// record the fact we have read one digit
		flags |= FL_READDIGIT;
		// we now need to compute number = number * base + digval
		number = number * ibase + digval;
		c = *p++;               // read next digit
	}
	--p;                            // point to place that stopped scan
	if (!(flags & FL_READDIGIT)) {
		// no number there; return 0 and point to beginning of string
		if (endptr)
			p = nptr;
		number = 0;
	}
	if (endptr != NULL)
		// store pointer to char that stopped the scan
		*endptr = (STRCHAR*) p;
	if (flags & FL_NEG)
		// negate result if there was a neg sign
		number = (unsigned long)(-(long)number);

        return number;
}

#endif

#if defined(STR_LINUX)
double StrImplementation::STR_strtod(const STRCHAR *s, STRCHAR **endptr)
{
#ifdef STR_UNICODE
	STRACHAR target[80];
	STR_ToAnsi(s, target, sizeof(target));
	STRACHAR* endptr_ansi;
	double v = strtod(target, &endptr_ansi);
	if (endptr_ansi)
		*endptr = (STRCHAR*) (s + (endptr_ansi - target));
	else
		*endptr = NULL;
	return v;
#else
	return strtod(s, endptr);
#endif
}
#endif


#ifdef STR_NO_RTLIBRARY

void StrImplementation::STR_copy (LPVOID dest, LPCVOID src, int bytes)
{
	BYTE* d = (BYTE*) dest;
	BYTE* s = (BYTE*) src;
	for (int i=0; i<bytes; i++)
		d[i] = s[i];
}

void StrImplementation::STR_move (LPVOID dest, LPCVOID src, int bytes)
{
	BYTE* d = (BYTE*) dest;
	BYTE* s = (BYTE*) src;
	if (dest < s) {
		for (int i=bytes-1; i>=0; i--)
			d[i] = s[i];
	}
	else
		STR_copy(dest, src, bytes);
}

#endif

int STR_EXPORT StrImplementation::STR_strcmp_ex(const STRCHAR *s, const STRCHAR *s2)
{
	int i = STR_strcmp (s, s2);
	if (i == 0)
		return 0;
	return (i < 0) ? -1 : 1;
}

int STR_EXPORT StrImplementation::STR_stricmp_ex(const STRCHAR *s, const STRCHAR *s2)
{
	int i = STR_stricmp (s, s2);
	if (i == 0)
		return 0;
	return (i < 0) ? -1 : 1;
}


/*********************************************************************
* Proc:	STR_sprintf_X
*********************************************************************/

// Important: this function cannot process %s or %S properly because
//   it has no features to reverse the polarity of string parameter
//   charsets!
static void STR_dosprintf_X(STRCHAR* dest, size_t destSZC, const STRCHAR* fmt, ...)
{
	va_list marker;
	va_start(marker, fmt);

#ifdef STR_DEBUG
	const STRCHAR* d_fmtp = fmt;
	while(d_fmtp) {
		d_fmtp = STR_scan(d_fmtp, (STRCHAR) '%');
		if (d_fmtp) {
			if (d_fmtp[1] == 0)
				break;				// Last char was %, ignore
			STR_ASSERT(d_fmtp[1] != (STRCHAR) 's');
			STR_ASSERT(d_fmtp[1] != (STRCHAR) 'S');
			STR_ASSERT(d_fmtp[1] != (STRCHAR) '!');		// Our special formatting symbol
			d_fmtp++;
		}
	}
#endif

#ifdef STR_UNICODE

#ifdef STR_WIN32

#ifdef STR_NO_RTLIBRARY
	// The wvsprintfX API does not support floating point numbers
	wvsprintfW(dest, fmt, marker);
#else
#if _MSC_VER >= 1400
	vswprintf_s(dest, destSZC, fmt, marker);		// VS2005 - safe buffer handling
#else
	(void) destSZC;
	vswprintf(dest, fmt, marker);			// Pre-VS2005
#endif
#endif

#else
	// Wide-character version of sprintf on Linux platforms; since we
	//   do not need to support %s / %S, a double-charset-conversion
	//   method will work just fine and we can also afford to use stack space
	int chars = STR_wstrlen(fmt);
	STRACHAR* temp_ansi_fmt = (STRACHAR*) alloca(chars+1);
	STR_ToAnsi(fmt, temp_ansi_fmt, chars+1);
	const int tr_size = 96;		// Should be more than enough even for a very long double
	STRACHAR* temp_result = (STRACHAR*) alloca(tr_size);
	vsprintf(temp_result, temp_ansi_fmt, marker);
	int result_chars = STR_astrlen(temp_result);	// Pseudo-max buffer size
	STR_ToWide(temp_result, dest, result_chars+1);
#endif

#else

#ifdef STR_NO_RTLIBRARY
	// The wvsprintfX API does not support floating point numbers
	wvsprintfA(dest, fmt, marker);
#else
#if _MSC_VER >= 1400
	vsprintf_s(dest, destSZC, fmt, marker);			// VS2005 safe buffer handling
#else
	(void) destSZC;
	vsprintf(dest, fmt, marker);			// Pre-VS2005
#endif
#endif

#endif
}

void StrImplementation::STR_sprintf_i(STRCHAR* dest, size_t destSZC, const STRCHAR* fmt, int value)
{
	STR_dosprintf_X(dest, destSZC, fmt, value);
}

void StrImplementation::STR_sprintf_p(STRCHAR* dest, size_t destSZC, const STRCHAR* fmt, void* value)
{
	STR_dosprintf_X(dest, destSZC, fmt, value);
}


/*********************************************************************
* Proc:	STR_ToUpper/STR_ToLower
*********************************************************************/

#ifndef STR_WIN32

void StrImplementation::STR_ToUpper (STRCHAR* buffer, DWORD blen)
{
	STRCHAR* x = buffer;
	while (blen > 0  &&  x[0] != 0) {
		STRCHAR ch = x[0];
		blen--;
		if (STR_CharIsAlpha(ch))
			x[0] = STR_CharToUpper(ch);
		x++;
	}
}

void StrImplementation::STR_ToLower (STRCHAR* buffer, DWORD blen)
{
	STRCHAR* x = buffer;
	while (blen > 0  &&  x[0] != 0) {
		STRCHAR ch = x[0];
		blen--;
		if (STR_CharIsAlpha(ch))
			x[0] = STR_CharToLower(ch);
		x++;
	}
}

#endif


/*********************************************************************
* Proc:		STR_stricmp
*********************************************************************/

#ifdef STR_LINUX

int StrImplementation::STR_stricmp(const STRCHAR* s, const STRCHAR* s2)
{
	for (int i=0;; i++) {
		// Get two chars, handle EOS
		STRCHAR c1 = s[i];
		STRCHAR c2 = s2[i];
		if (c1 == 0) {
			if (c2 == 0)
				return 0;
			// s ends before s2, so it must be 'smaller'
			return -1;
		}
		// s ends after s2, so it must be 'greater'
		if (c2 == 0)
			return 1;
		// Compare them insensitively
		c1 = STR_CharToUpper(c1);
		c2 = STR_CharToUpper(c2);
		if (c1 == c2)
			continue;
		return (c1 > c2) ? 1 : -1;
	}
}

#ifdef STR_UNICODE

int StrImplementation::STR_strcmp(const STRCHAR* s, const STRCHAR* s2)
{
	for (int i=0;; i++) {
		// Get two chars, handle EOS
		STRWCHAR c1 = s[i];
		STRWCHAR c2 = s2[i];
		if (c1 == 0) {
			if (c2 == 0)
				return 0;
			// s ends before s2, so it must be 'smaller'
			return -1;
		}
		// s ends after s2, so it must be 'greater'
		if (c2 == 0)
			return 1;
		// Compare them
		if (c1 == c2)
			continue;
		return (c1 > c2) ? 1 : -1;
	}
}

int StrImplementation::STR_strncmp(const STRCHAR *s, const STRCHAR *s2, int nchars)
{
	// Unicode AND Linux version of strncmp
	for (int i=0; i<nchars; i++) {
		// Get two chars, handle EOS
		STRWCHAR c1 = s[i];
		STRWCHAR c2 = s2[i];
		if (c1 == 0  &&  c2 == 0)
			return 0;
		// No need to handle case 'just one string ends' because it will
		//   automatically result in the correct 'mismatch' comparision
		if (c1 != c2)
			return (c1 > c2) ? 1 : -1;
	}
	// We used up all characters for the comparision
	return 0;
}

#endif

int StrImplementation::STR_strnicmp(const STRCHAR *s, const STRCHAR *s2, int nchars)
{
	// Linux version of strnicmp (unicode / ansi)
	for (int i=0; i<nchars; i++) {
		// Get two chars, handle EOS
		STRCHAR c1 = s[i];
		STRCHAR c2 = s2[i];
		if (c1 == 0  &&  c2 == 0)
			return 0;
		// No need to handle case 'just one string ends' because it will
		//   automatically result in the correct 'mismatch' comparision
		STRCHAR c1_u = STR_CharToUpper(c1);
		STRCHAR c2_u = STR_CharToUpper(c2);
		if (c1_u != c2_u)
			return (c1_u > c2_u) ? 1 : -1;
	}
	// We used up all characters for the comparision
	return 0;
}

#endif


/*********************************************************************
* Proc:		STR_scan / STR_rscan (equivalent to strchr / strrchr)
*********************************************************************/

#if !defined(STR_UNICODE)  &&  !defined(STR_NO_RTLIBRARY)

inline STRCHAR* StrImplementation::STR_scan(const STRCHAR* s, STRCHAR c)
{
	STR_ASSERT(c != 0);
	return (STRCHAR*)((s[0] == 0) ? NULL : strchr(s, c));
}

inline STRCHAR* StrImplementation::STR_rscan(const STRCHAR* s, STRCHAR c)
{
	STR_ASSERT(c != 0);
	return (STRCHAR*)((s[0] == 0) ? NULL : strrchr(s, c));
}

inline STRCHAR* StrImplementation::STR_scan2(const STRCHAR* s, STRCHAR c, int nchars)
{
	STR_ASSERT(nchars > 0);
	STR_ASSERT(c != 0);
	return (STRCHAR*)((s[0] == 0) ? NULL : (STRCHAR*) memchr (s, c, nchars));
}

#else

inline STRCHAR* StrImplementation::STR_scan(const STRCHAR* s, STRCHAR c)
{
	STR_ASSERT(c != 0);
	while (s[0] != c  &&  s[0] != 0)
		s++;
	if (s[0] == 0)
		return NULL;
	return (STRCHAR*) s;
}

STRCHAR* StrImplementation::STR_rscan(const STRCHAR* s, STRCHAR c)
{
	STR_ASSERT(c != 0);
	// Find end of string
	const STRCHAR* endptr = s;
	while (endptr[0] != 0)
		endptr++;
	if (endptr == s)
		return NULL;
	// Scan backwards
	while (--endptr >= s) {
		if (endptr[0] == c)
			return (STRCHAR*) endptr;
	}
	return NULL;
}

#if defined(STR_LINUX)  ||  defined(STR_NO_RTLIBRARY)  ||  defined(STR_WINDOWS_CE)

STRCHAR* StrImplementation::STR_scan2(const STRCHAR* s, STRCHAR c, int nchars); // CodeWarrior fix

STRCHAR* StrImplementation::STR_scan2(const STRCHAR* s, STRCHAR c, int nchars)
{
	STR_ASSERT(nchars > 0);
	STR_ASSERT(c != 0);
	while (nchars > 0  &&  s[0] != c) {
		nchars--;
		s++;
	}
	return (nchars == 0) ? NULL : (STRCHAR*) s;
}

#else

inline STRCHAR* StrImplementation::STR_scan2(const STRCHAR* s, STRCHAR c, int nchars)
{
	STR_ASSERT(nchars > 0);
	STR_ASSERT(c != 0);
	return (s[0] == 0) ? NULL : (STRCHAR*) wmemchr (s, c, nchars);
}

#endif

#endif


#ifdef STR_UNICODE

STRACHAR* STR_ansiscan(const STRACHAR *s, STRCHAR c)
{
	STR_ASSERT(c != 0);
	// Handle bad unicode characters
	if (c > 255)		// No need to check for c<0
		return NULL;
	// The rest is simple
	while (s[0] != (STRACHAR) c  &&  s[0] != 0)
		s++;
	if (c != 0  &&  s[0] == 0)
		return NULL;
	return (STRACHAR*) s;
}

#else

#define STR_ansiscan STR_scan

#endif


/*********************************************************************
* Proc:		Str::copy constructor
*********************************************************************/

Str::Str(const Str& source)
{
	STR_ASSERT(this != &source);
	if (source.dptr == null_str()) {
		dptr = null_str();		// source is always single-threaded here
		return;
	}
	// Spawn a non-MT copy if source is MT-marked
#ifdef STR_THREAD_SAFE
	if (source.IsMT()) {
		dptr = null_str();
		*this = (const STRCHAR*) source;
		return;
	}
#endif
	// Regular copy: just increase refcount
	STR_Inc(source.block()->m_Ref);
	dptr = source.dptr;
}


/*********************************************************************
* Proc:		Str::NewEmptyMT
*********************************************************************/

#ifdef STR_THREAD_SAFE
STRCHAR* StrImplementation::NewEmptyMT(CPOS prealloc)
#else
STRCHAR* StrImplementation::NewEmptyMTForceNI(CPOS prealloc)
#endif
{
	// Round up to blkalloc() granularity
	up_alloc (prealloc);
	// Get memory for string
	SBlock* target = STR_get_block (prealloc);
	target->SetFlag(SFLAG_MT);
	STRCHAR* dptr = Str::block_s(target);
	// Fill in length, place a null terminator
	target->m_Length = 0;
	*dptr = 0;
	
	return dptr;
}


/*********************************************************************
* Proc:		Str::NewFromString
* Purpose:	Core code for 'copy char* constructor' and '= operator'.
*			Assumes our 'data' field is garbage on entry.
*********************************************************************/

void Str::NewFromString(const STRCHAR* s, CPOS slen, CPOS prealloc, BOOL want_mt)
{
	// Determine actual size of buffer that needs to be allocated
	if (slen > prealloc)
		prealloc = slen;
	// Empty string?
	if (slen == 0  &&  prealloc == 0) {
		dptr = want_mt ? NewEmptyMT(0) : null_str();
		return;
	}
	// Round up to blkalloc() granularity
	up_alloc (prealloc);
	// Get memory for string
	SBlock* target = STR_get_block (prealloc);
	if (want_mt)
		target->SetFlag(SFLAG_MT);
	dptr = block_s(target);
	// Fill new descriptor
	target->m_Length = slen;
	// Copy the string bytes, including the null.  Use copy, not
	//   move, because we're guaranteed we won't overlap with source
	STR_copy (dptr, s, (slen+1) * sizeof(STRCHAR));
}

void Str::NewFromString(const STRCHAR* s, CPOS prealloc, BOOL want_mt)
{
	NewFromString(s, STR_strlen(s), prealloc, want_mt);
}


/*********************************************************************
* Proc:		Str::Str (CPOS prealloc)
* Purpose:	Instantiates an empty string, but with the specified
*			number of bytes in the preallocated buffer.
* In:		prealloc - number of characters to reserve
*********************************************************************/

#if 0
// Obsolete - removed because it was often a source of confusion
Str::Str(CPOS prealloc)
{
	up_alloc (prealloc);
	// Get memory for string
	StrImplementation::SBlock* target = STR_get_block (prealloc);
	dptr = block_s(target);
	// Fill in descriptor values left out by get_block
	target->m_Length = 0;
	dptr[0] = 0;
}
#endif


/*********************************************************************
* Proc:		Str::Str(const CString&)
* Purpose:	Implemented only when using MFC.  Instantiates the CString
*			from an MFC-provided one.
*********************************************************************/

#if defined(__AFX_H__)  &&  !defined(STR_NO_WINSTUFF)
Str::Str(const CString& source, CPOS prealloc /*= 0*/)
{
	NewFromString((const STRCHAR*) source, (CPOS) source.GetLength(), prealloc, FALSE);
}
#endif


/*********************************************************************
* Proc:		Str::operator = [const CString&]
* Purpose:	Implemented only when MFC is used.  
*********************************************************************/

#if defined(__AFX_H__)  &&  !defined(STR_NO_WINSTUFF)
const Str& Str::operator=(const CString& source)
{
	*this = (const STRCHAR*) source;
	return *this;
}
#endif


#if defined(STR_THREAD_SAFE)  &&  defined(STR_DEBUG)
StrImplementation::verifymt::verifymt(Str* obj)
{
	the_object = obj;
	if (Str::block_s(obj->block()) == null_str())
		return;
	DWORD cur = STR_GetCurrentThreadID();
	if (obj->block()->m_OwnerT == cur)
		return;
	if (obj->block()->GetFlag(SFLAG_MT))
		return;
	obj->Error(Str::SE_BadMTUse);
}
#endif


/*static*/ void Str::ConstructAt(Str* location, int count)
{
#ifndef STR_BORCPPBUILDER
	STR_ASSERT(sizeof(Str) == sizeof(char*));
#endif
	for (int i=0; i<count; i++)
		location[i].dptr = null_str();
}

/*static*/ void Str::ConstructAt(Str* location, int count, const StrMTType&)
{
#ifndef STR_BORCPPBUILDER
	STR_ASSERT(sizeof(Str) == sizeof(char*));
#endif
	for (int i=0; i<count; i++)
#ifdef STR_THREAD_SAFE
		location[i].dptr = StrImplementation::NewEmptyMT(0);
#else
		location[i].dptr = null_str();
#endif
}


/*static*/ void Str::DestructAt(Str* location, int count /*= 1*/)
{
#ifndef STR_BORCPPBUILDER
	STR_ASSERT(sizeof(Str) == sizeof(STRCHAR*));
#endif
	for (int i=0; i<count; i++) {
		StrImplementation::verifymt mtcheck(&(location[i]));
		StrImplementation::STR_rls_block(mtcheck);
	}
}


/*********************************************************************
* syncstr is a helper class for serializing write access to the
*	string data across all application threads.
*********************************************************************/

struct StrImplementation::syncstr
{
#ifdef STR_THREAD_SAFE
	SBlock*		m_ToUnlock;		// NULL means "nothing to unlock"
#endif
	SBlock*		keep_around;	// Non-null means "deref + maybe free in dtor"

	syncstr(verifymt& mtobj, CPOS newlength, BOOL copyspawn);
	~syncstr();
};


// li_spawn is a helper routine that actually spawns a copy
//   of the string object's data.  it can be called only by
//   syncstr; the copyspawn parameter determines whether
//   the old string contents will be copied to the new place.
void Str::li_spawn(CPOS newlength, BOOL copyspawn, syncstr* caller)
{
	// Get a fixed length block belonging to the current thread
	SBlock* prevdata = block();
	BOOL was_mt = prevdata->GetFlag(SFLAG_MT);
	if (newlength < prevdata->m_Length)
		newlength = prevdata->m_Length;
	up_alloc(newlength);
	// Get memory for new block
	SBlock* target = STR_get_block (newlength);
	dptr = block_s(target);
	if (was_mt)
		target->SetFlag(SFLAG_MT);
	// Copy data bytes, or ignore previous content?
	if (copyspawn) {
		CPOS slen = prevdata->m_Length;
		target->m_Length = slen;
		target->m_Alloc  = newlength;
		STR_copy (dptr, block_s(prevdata), (slen+1) * sizeof(STRCHAR));
	}
	else {
		target->m_Length = 0;
		dptr[0] = 0;
	}
	// Keep previous buffer until the end of the lifetime of syncstr
	//   (unless it's the empty string instance)
	STR_ASSERT(caller->keep_around == NULL);
	STRCHAR* prevdata_s = Str::block_s(prevdata);
	if (prevdata_s != null_str())
		caller->keep_around = prevdata;
}


#ifdef STR_THREAD_SAFE

static void GetMtObjectLock(SBlock* objblock)
{
	for (;;) {
		if (STR_Inc(objblock->m_Locked) == 1)
			break;
		STR_Dec(objblock->m_Locked);
		Sleep(0);
	}
}

inline static void FreeMtObjectLock(SBlock* objblock)
{
	STR_Dec(objblock->m_Locked);
}

#endif



#if !defined(STR_THREAD_SAFE)  &&  defined(STR_NODEBUG)
inline
#endif
syncstr::syncstr(verifymt& mtobj, CPOS newlength, BOOL copyspawn)
{
	Str* obj = mtobj.the_object;
	//
	// Set basic structure elements
	keep_around = NULL;
	SBlock* objblock = obj->block();
	long ob_Ref = objblock->m_Ref;
	CPOS ob_Alloc = objblock->m_Alloc;
#ifdef STR_THREAD_SAFE
	m_ToUnlock = NULL;
#endif
	if (newlength == 0)
		newlength = ob_Alloc;
	STR_ASSERT_UPPED(newlength);

	//
	// Trying to work with null strings always results in expansion
	if (Str::block_s(objblock) == null_str()) {
		obj->li_spawn(newlength, copyspawn, this);
		return;
	}

	//
	// Get a lock for the string data block if necessary
#ifdef STR_THREAD_SAFE
	BOOL is_mt = objblock->GetFlag(SFLAG_MT);
	if (is_mt)
		GetMtObjectLock(objblock);
#endif

	//
	// If it's a single instance, and MT-enabled, remain in locked state
	STR_ASSERT(objblock->m_Ref > 0);
	if (ob_Ref == 1  &&  newlength <= ob_Alloc) {
#ifdef STR_THREAD_SAFE
		if (is_mt)
			m_ToUnlock = objblock;
#endif
		return;
	}

	//
	// Either the instance is not only one, or we need to expand
	obj->li_spawn(newlength, copyspawn, this);
#ifdef STR_THREAD_SAFE
	STR_ASSERT(m_ToUnlock == NULL);		// No lock to decrement
#endif
}

inline syncstr::~syncstr()
{
	SBlock* cache_keep_around = keep_around;	// Avoid refetches by the compiler
#ifdef STR_THREAD_SAFE
	if (m_ToUnlock)
		FreeMtObjectLock(m_ToUnlock);
	BOOL is_mt = cache_keep_around  &&  cache_keep_around->GetFlag(SFLAG_MT);
	BOOL must_unlock_keep = is_mt ? (cache_keep_around->m_Ref > 1) : FALSE;
#endif
	if (cache_keep_around) {
		STRCHAR* dummy_content = Str::block_s(cache_keep_around);
		Str* dummy_content2 = (Str*) (&dummy_content);
		verifymt mtcheck_dummy (dummy_content2, -1);
		STR_rls_block(mtcheck_dummy);
#ifdef STR_THREAD_SAFE
		if (must_unlock_keep)
			FreeMtObjectLock(cache_keep_around);
#endif
	}
}


/*********************************************************************
* Proc:		SetMT (thread safe implementation only)
*********************************************************************/

#ifdef STR_THREAD_SAFE
void Str::SetMT()
{
	if (!IsMT()) 
		SetMT_Internal();
	else {
		verifymt mtcheck(this);
		syncstr locker (mtcheck, 0, TRUE);
	}
}
#endif


/*********************************************************************
* Proc:		ChangeOwnerThread (thread safe implementation only)
*********************************************************************/

#ifdef STR_THREAD_SAFE
void Str::ChangeOwnerThread()
{
	STR_ASSERT(IsMT());
	SBlock* blk = block();
	GetMtObjectLock(blk);
	block()->m_OwnerT = STR_GetCurrentThreadID();
    FreeMtObjectLock(blk);
}
#endif


/*********************************************************************
* Proc:		Str::Compact
* Purpose:	If m_Alloc is bigger than m_Length, shrinks the
*			buffer to hold only as much memory as necessary.
* In:		only_above - will touch the object only if the difference
*				is greater than this value.  By default, 0 is passed,
*				but other good values might be 7 and up.
* Rems:		Will compact even the buffer for a string shared between
*			several Str instances.
*********************************************************************/

void Str::Compact(CPOS only_above /*= 0*/)
{
	verifymt mtcheck(this);
	// Nothing to do?
	if (dptr == null_str())
		return;
	// The algorithm is as follows:
	//   1) Do nothing if this is a multithreaded object
	//   2) Do nothing if the reference is >1 (because we don't
	//      have the list of instances referring the block and
	//      therefore cannot reallocate and patch all of them;
	//		also, in a multithreaded env we might be interrupted
	//		in a bad place)
	//   3) Do nothing if the difference in bytes would be
	//      too small after the compaction
	//   4) Do the actual compaction
#ifdef STR_THREAD_SAFE
	if (IsMT())
		return;
#endif
	int slen = block()->m_Length;
	if (block()->m_Ref > 1)
		goto CantCompact;
	if ((block()->m_Alloc - slen) <= only_above)
		goto CantCompact;
	{
		CPOS to_shrink = slen;
		up_alloc (to_shrink);
		if ((to_shrink - slen) <= only_above)
			goto CantCompact;
		// Reallocate safely; first, check for empty string
		SBlock* prevb = block();
		STRCHAR*  prevs = dptr;
		if (slen == 0) {
			Reset();
			if (block() != prevb)
				STR_rls_block(mtcheck);
			return;
		}
		// Get new buffer for string
		SBlock* target = STR_get_block (to_shrink);
		dptr = block_s(target);
		// Fill new descriptor, copy old data
		target->m_Length = slen;
		STR_copy (dptr, prevs, (slen+1) * sizeof(STRCHAR));
		// Release old block
		STR_rls_block(mtcheck);
		return;
	}
	// Couldn't do anything
CantCompact:
	(void) 0;
}


/*********************************************************************
* Proc:		Str::Copy_XXX methods
*********************************************************************/

void Str::Copy_Native(const STRCHAR* s)
{
	verifymt mtcheck(this);
	// Check for zero length string.
	CPOS slen = STR_strlen(s);
	BOOL ismt = IsMT();
	if (slen == 0) {
		if (dptr != null_str()) {
			STR_rls_block(mtcheck);
			dptr = ismt ? NewEmptyMT(0) : null_str();
		}
		return;
	}
	// Do the copying
	CPOS alen = slen;
	up_alloc(alen);
	syncstr keeper (mtcheck, alen, FALSE);
	block()->m_Length = slen;
	STR_copy (dptr, s, (slen+1) * sizeof(STRCHAR));
}

#ifdef STR_UNICODE
void Str::Copy_FromAnsi(const STRACHAR *value)
{
	verifymt mtcheck(this);
	// Check how many characters will the outstring have
	int need_ch = STR_ToWideSlots(value);
	STR_ASSERT(need_ch > 0);
	BOOL ismt = IsMT();
	if (need_ch == 1) {
		// Empty string
		if (dptr != null_str()) {
			STR_rls_block(mtcheck);
			dptr = ismt ? NewEmptyMT(0) : null_str();
		}
		return;
	}
	// Perform actual conversion.
	CPOS need_slots = need_ch-1;	// Because of terminating null
	up_alloc(need_slots);
	syncstr keeper (mtcheck, need_slots, FALSE);
	block()->m_Length = need_ch-1;
	STR_ToWide(value, dptr, need_ch);
}
#endif

#if !defined(STR_NO_UNICODE)  &&  defined(STR_ANSI)
void Str::Copy_FromUni(const STRWCHAR* value)
{
	verifymt mtcheck(this);
	// Check how many characters will the outstring have
	int need_ch = STR_ToAnsiBytes(value);
	STR_ASSERT(need_ch > 0);
	BOOL ismt = IsMT();
	if (need_ch == 1) {
		// Empty string
		if (dptr != null_str()) {
			STR_rls_block(mtcheck);
			dptr = ismt ? NewEmptyMT(0) : null_str();
		}
		return;
	}
	// Perform actual conversion.
	CPOS need_bytes = need_ch-1;	// Because of terminating null
	up_alloc(need_bytes);
	syncstr keeper (mtcheck, need_bytes, FALSE);
	block()->m_Length = need_ch-1;
	STR_ToAnsi(value, dptr, need_ch);
}
#endif


/*********************************************************************
* Proc:		CsCopyCore
* Purpose:	Fast helper for Str::operator= Str&
*********************************************************************/

STRCHAR* StrImplementation::STR_CsCopyCore(Str* this_obj, const Str& src)
{
	verifymt mtcheck (this_obj);
	STRCHAR* this_dptr = Str::block_s(this_obj->block());
	const STRCHAR* src_dptr = (const STRCHAR*) src.GetString();
	// Self-assignments must be filtered out
	if (this_dptr == src_dptr)
		return this_dptr;
	// If source is MT and we're not, force a spawn
#ifdef STR_THREAD_SAFE
	if (src_dptr != null_str()  &&  src.IsMT()  &&  !this_obj->IsMT()) {
		*this_obj = (const STRCHAR*) src_dptr;		// Forcibly decouple
		return (STRCHAR*) this_obj->GetString();
	}
#endif
	// Increase source string's ref count
	if (src_dptr != null_str()) {
		SBlock* foreign = (SBlock*) (((BYTE*) src_dptr) - sizeof(SBlock));
		STR_Inc(foreign->m_Ref);
	}
	// Dereference our own buffer, maybe release it
	if (this_dptr != null_str())
		STR_rls_block(mtcheck);
	return (STRCHAR*) src_dptr;
}


/*********************************************************************
* Proc:		Str::Empty
* Purpose:	Sets the string to null.  However, the allocated buffer
*			is not released.
*********************************************************************/

Str& Str::Empty()
{
	verifymt mtcheck(this);
	// Already empty, and with buffer zero?
	if (dptr == null_str())
		return *this;
	// More than one instance served?  Alleviate situation.
	syncstr keeper (mtcheck, 0, FALSE);
	dptr[0] = 0;
	block()->m_Length = 0;
	return *this;
}


/*********************************************************************
* Proc:		Str::Reset
* Purpose:	Sets the string to null, deallocates buffer
*********************************************************************/

Str& Str::Reset()
{
	verifymt mtcheck(this);
	if (dptr != null_str()) {
		BOOL ismt = IsMT();
		STR_rls_block(mtcheck);
		dptr = ismt ? NewEmptyMT(0) : null_str();
	}
	return *this;
}


/*********************************************************************
* Proc:		Str::AppendChar
* Purpose:	Appends a single character to the end of the string
*********************************************************************/

#define IMPLEMENT_APPENDCHAR_NATIVECHARSET				\
	CPOS new_len = block()->m_Length + 1;				\
	CPOS new_alen = new_len;							\
	up_alloc(new_alen);									\
	syncstr keeper (mtcheck, new_alen, TRUE);			\
	dptr[new_len-1] = ch;								\
	dptr[new_len  ] = 0;								\
	block()->m_Length = new_len;

Str& Str::AppendChar(STRACHAR ch)
{
#if !defined(STR_LINUX)  ||  !defined(STR_UNICODE)
	verifymt mtcheck(this);
#endif
#ifdef STR_ANSI
	IMPLEMENT_APPENDCHAR_NATIVECHARSET
#else
	STRACHAR buf[2];
	buf[0] = ch;
	buf[1] = 0;
	AppendString(buf);
#endif
	return *this;
}

#ifndef STR_NO_UNICODE
Str& Str::AppendChar(STRWCHAR ch)
{
#if !defined(STR_LINUX)  ||  defined(STR_UNICODE)
	verifymt mtcheck(this);
#endif
#ifdef STR_ANSI
	STRWCHAR buf[2];
	buf[0] = ch;
	buf[1] = 0;
	AppendString(buf);
#else
	IMPLEMENT_APPENDCHAR_NATIVECHARSET
#endif
	return *this;
}
#endif


/*********************************************************************
* Proc:		Str::AppendInt
* Purpose:	Appends a decimal signed integer, possibly with - sign
*********************************************************************/

Str& Str::AppendInt(int value)
{
	STRCHAR buf[32];
	STR_sprintf_i(buf, sizeof(buf)/sizeof(STRCHAR), _T("%d"), value);
	AppendString (buf);
	return *this;
}


/*********************************************************************
* Proc:		Str::AppendFloat, AppendDouble
* Purpose:	Append a signed value, use specified # of digits
*********************************************************************/

#ifndef STR_NO_RTLIBRARY

Str& Str::AppendFloat(float value, UINT after_dot)
{
	if (after_dot > 32)
		after_dot = 32;
	STRCHAR fmt[16];
	STRCHAR buf[64];
	STR_sprintf_i (fmt, sizeof(fmt)/sizeof(STRCHAR), _T("%%.%df"), (int) after_dot);
	STR_dosprintf_X(buf, sizeof(buf)/sizeof(STRCHAR), fmt, (double) value);
	AppendString (buf);
	return *this;
}

Str& Str::AppendDouble(double value, UINT after_dot)
{
	if (after_dot > 48)
		after_dot = 48;
	STRCHAR fmt[16];
	STRCHAR buf[128];
	STR_sprintf_i (fmt, sizeof(fmt)/sizeof(STRCHAR), _T("%%.%df"), (int) after_dot);
	STR_dosprintf_X(buf, sizeof(buf)/sizeof(STRCHAR), fmt, value);
	AppendString (buf);
	return *this;
}

#endif


/*********************************************************************
* Proc:		Str::ToInt, ToFloat, ToDouble
* Purpose:	Convert the string to number
*********************************************************************/

static BOOL ToNumCore(const Str* obj, BOOL* error, STRCHAR* endptr)
{
	// No error?
	if (endptr == NULL  ||  endptr[0] == 0) {
		if (error)
			*error = FALSE;
		return TRUE;
	}
	// Oops.  Either set *error, or throw an exception
	if (error == NULL)
		obj->Error(Str::SE_BadNumber);
	*error = TRUE;
	return FALSE;
}

int Str::ToInt(BOOL* error /*= NULL*/) const
{
	STRCHAR* endptr = NULL;
	long l = STR_strtol(dptr, &endptr, 10);
	return ToNumCore(this, error, endptr) ? l : 0;
}

#ifndef STR_NO_RTLIBRARY

float Str::ToFloat(BOOL* error /*= NULL*/) const
{
	return (float) ToDouble(error);
}

double Str::ToDouble(BOOL* error /*= NULL*/) const
{
	STRCHAR* endptr = NULL;
	double v = STR_strtod(dptr, &endptr);
	return ToNumCore(this, error, endptr) ? v : 0.0;
}

#endif


/*********************************************************************
* Proc:		Str::CoreAppendChars
* Purpose:	Core code for AppendChars() and operators +=
*********************************************************************/

void Str::CoreAppendChars(const STRCHAR* s, CPOS howmany, verifymt& mtcheck)
{
	STR_ASSERT(mtcheck.the_object == this);
	if (howmany == 0)
		return;
	// Prepare big enough buffer.
	CPOS to_alloc = GetLength() + howmany;
	up_alloc(to_alloc);
	syncstr keeper (mtcheck, to_alloc, TRUE);
	// And copy the bytes
	STRCHAR* dest = dptr + block()->m_Length;
	STR_copy (dest, s, howmany*sizeof(STRCHAR));
	dest[howmany] = 0;
	block()->m_Length += howmany;
}


/*********************************************************************
* Proc:		Str::SetAt
*********************************************************************/

void Str::SetAt(CPOS idx, STRCHAR newch)
{
	verifymt mtcheck(this);
#ifdef STR_DEBUG
	if (idx >= GetLength())
		Error(SE_BadCharPos);
	if (newch == 0)
		Error(SE_BadCharPos);
#endif
	syncstr keeper (mtcheck, 0, TRUE);	// copy-on-write?
	dptr[idx] = newch;
}


/*********************************************************************
* Proc:		Str::operator += (both from const char* and from Str)
* Purpose:	Append a string to what we already contain.
*********************************************************************/

void Str::operator += (const STRCHAR* s)
{
	verifymt mtcheck (this);
	CPOS slen = STR_strlen(s);
	CoreAppendChars (s, slen, mtcheck);
}

#ifdef STR_UNICODE
void Str::operator += (const STRACHAR* s)
{
	Str temp(s);
	*this += temp;
}
#endif
#if !defined(STR_NO_UNICODE)  &&  defined(STR_ANSI)
void Str::operator += (const STRWCHAR* s)
{
	Str temp(s);
	*this += temp;
}
#endif


/*********************************************************************
* Proc:		Str::AppendChars
* Purpose:	Catenate a number of characters to our internal data.
*********************************************************************/

void Str::AppendChars(const STRCHAR* s, CPOS howmany /*= -1*/)
{
	verifymt mtcheck(this);
	if (howmany == (CPOS) -1)
		howmany = STR_strlen(s);
	CoreAppendChars (s, howmany, mtcheck);
}


/*********************************************************************
* Proc:		Str::Insert
* Purpose:	Insert a substring at a specified character position
*********************************************************************/

Str& Str::Insert(const STRCHAR* s, CPOS startat, CPOS howmany /*= -1*/)
{
	verifymt mtcheck(this);
	int s_len = STR_strlen(s);
	if (howmany == (CPOS) -1)
		howmany = s_len;
	else
		howmany = STR_min (howmany, s_len);
	if (howmany == 0)
		return *this;
	if (startat > GetLength())		// startat=Len means basically "append"
		Error(SE_BadScanPos);
	// Enlarge the buffer
	CPOS to_alloc = GetLength() + howmany;
	up_alloc(to_alloc);
	syncstr keeper (mtcheck, to_alloc, TRUE);
	// Open a slot for insertion
	STRCHAR* str_ptr = dptr + startat;
	STR_move(str_ptr+howmany, str_ptr, (GetLength()-startat+1) * sizeof(STRCHAR));
	// Copy the bytes for the to-be-inserted string
	STR_copy (str_ptr, s, howmany * sizeof(STRCHAR));
	// Adjust length
	block()->m_Length += howmany;
	return *this;
}


/*********************************************************************
* Proc:		operator +LPCSTR for Str
*********************************************************************/

STR_EXPORT Str operator+(const STRCHAR* lpsz, const Str& s)
{
	Str temp;
	if (s.IsMT())
		temp.SetMT();
	temp.Preallocate(s.GetLength() + STR_strlen(lpsz));
	temp  = lpsz;
	temp += s;
	return temp;
}


/*********************************************************************
* Proc:		operator +char for Str
*********************************************************************/

STR_EXPORT Str operator+(STRCHAR ch, const Str& s)
{
	Str temp;
	if (s.IsMT())
		temp.SetMT();
	temp.Preallocate(s.GetLength() + 1);
	temp.AppendChar (ch);
	temp += s;
	return temp;
}


/*********************************************************************
* Proc:		Str::GetLastChar
*********************************************************************/

STRCHAR Str::GetLastChar() const
{
	CPOS l = GetLength();
#ifdef STR_DEBUG
	if (l < 1)
		Error(SE_StringEmpty);
#endif
	return dptr[l-1];
}


/*********************************************************************
* Proc:		Str::TruncateAt
* Purpose:	Cuts off the string at the character with the specified
*			index.  The allocated buffer remains the same.
*********************************************************************/

void Str::TruncateAt (CPOS idx)
{
	verifymt mtcheck(this);
	if (idx >= GetLength())
		return;
	// Spawn a copy if necessary (preserve buffer size, though)
	syncstr keeper (mtcheck, 0, TRUE);
	// And do the truncation
	dptr[idx] = 0;
	block()->m_Length = idx;
}


/*********************************************************************
* Proc:		Str::Find and ReverseFind
* Purpose:	Scan the string for a particular character (must not
*			be 0); return the index where the character is found
*			first, or -1 if cannot be met
*********************************************************************/

int Str::Find (STRCHAR ch, CPOS startat /*= 0*/) const
{
	STR_ASSERT(ch != 0);
	// Start from middle of string?
	if (startat > 0) {
		if (startat > GetLength())
			Error(SE_BadScanPos);
	}
	const STRCHAR* scan = STR_scan (dptr+startat, ch);
	if (scan == NULL)
		return -1;
	else
		return (CPOS) (scan - dptr);
}

int Str::ReverseFind (STRCHAR ch, CPOS startat /*= -1*/) const
{
	if (startat == (CPOS) -1) {
		// Scan entire string
		const STRCHAR* scan = STR_rscan(dptr, ch);
		if (scan)
			return (CPOS) (scan - dptr);
	}
	else {
		// Make sure the index is OK
		if (startat >= GetLength())
			Error(SE_BadScanPos);
		for (int findex = (int) startat; findex >= 0; findex--) {
			if (dptr[findex] == ch)
				return findex;
		}
	}
	return -1;
}


/*********************************************************************
* Proc:		Str::FindNoCase and ReverseFindNoCase
*********************************************************************/

int Str::FindNoCase(STRCHAR ch, CPOS startat /*= 0*/) const
{
	Str self_copy (*this);
	self_copy.MakeLower();
	STRCHAR ch_copy = ch;
	STR_ToLower(&ch_copy, 1);
	return self_copy.Find(ch_copy, startat);
}

int Str::ReverseFindNoCase(STRCHAR ch, CPOS startat /*= -1*/) const
{
	Str self_copy (*this);
	self_copy.MakeLower();
	STRCHAR ch_copy = ch;
	STR_ToLower(&ch_copy, 1);
	return self_copy.ReverseFind(ch_copy, startat);
}


/*********************************************************************
* Proc:		Str::FindOneOf
*********************************************************************/

int Str::FindOneOf (const STRCHAR* charset, CPOS startat /*= 0*/) const
{
	// Bad input params?
	if (charset[0] == 0)
		Error(SE_StringEmpty);
	if (startat > 0) {
		if (startat >= GetLength())
			Error(SE_BadScanPos);
	}
	CPOS foundat = StringSpanExcluding(dptr+startat, charset);
	if (foundat >= 0)
		foundat += startat;
	return foundat;
}


/*********************************************************************
* Proc:		Str::Preallocate
* Purpose:	If the buffer is smaller than the amount of characters
*			specified, reallocates the buffer.  This function cannot
*			reallocate to a buffer smaller than the existing one.
*********************************************************************/

void Str::Preallocate(CPOS size)
{
	verifymt mtcheck(this);
	up_alloc(size);
	syncstr keeper(mtcheck, size, TRUE);
}


/*********************************************************************
* Proc:		Str::operator == (basic forms, the rest are inline)
*********************************************************************/

BOOL operator ==(const Str& s1, const Str& s2)
{
	CPOS slen = s2.GetLength();
	if (s1.GetLength() != slen)
		return FALSE;
	return STR_memcmp(s1.operator const STRCHAR*(), s2, slen*sizeof(STRCHAR)) == 0;
}

BOOL operator ==(const Str& s1, const STRCHAR* s2)
{
	CPOS slen = STR_strlen(s2);
	if (s1.GetLength() != slen)
		return FALSE;
	return STR_memcmp(s1.operator const STRCHAR*(), s2, slen*sizeof(STRCHAR)) == 0;
}


/*********************************************************************
* Proc:		Str::FmtPriorBuffer
* Purpose:	Helper for Str::Format -- estimates the size of the
*			buffer which Format() will preallocate before it will
*			attempt to do the formatting.
*********************************************************************/

void Str::FmtPriorBuffer (const STRCHAR* st, verifymt& mtcheck)
{
	// The whole thing below is done only if we're an
	//    empty string with an empty buffer.  Otherwise,
	//    we presume the caller has allocated enough bytes.
	if (dptr != null_str()) {
		Empty();			// bug fix Mar 15, 99
		return;
	}
	// Compute how much space we need to allocate
	CPOS bl = 0;
	const STRCHAR* x = st;
	while (TRUE) {
		// Locate next % sign
		const STRCHAR* x0 = x;
		while (x[0] != 0  &&  x[0] != '%')
			x++;
		bl += (CPOS) (x-x0);
		if (x[0] == 0)
			break;
		// Add some more space
		if (x[1] == 's'  ||  x[1] == 'S')
			bl += 4 * (WIDE_BAE_MASK+1);
		else
			bl += (WIDE_BAE_MASK+1);
		// Go on
		x++;
	}
	up_alloc(bl);
	// Provide suitable minimum buffer size, erase string
	syncstr keeper(mtcheck, bl, FALSE);
	dptr[0] = 0;
	block()->m_Length = 0;
}


/*********************************************************************
* Proc:		Str_MakeUL
* Purpose:	Local helper for MakeUpper() and MakeLower()
* Rems:		Note the trick of using a temporary string and
*			assigning only if no changes were found.  In this way,
*			we can assure that multiple Str instances do not get
*			'severed' from each other (syncstr() would do just that).
*********************************************************************/

void Str::Str_MakeUL(BOOL do_upr)
{
#ifdef STR_THREAD_SAFE
	verifymt mtcheck(this);
#endif
	// Get a temporary copy of the string and work on it
	int slen   = GetLength();
	int sbytes = (slen+1) * sizeof(STRCHAR);
	const int max_onstack = 128;
	BYTE tempcopy_local[max_onstack];
	STRCHAR* tempcopy = (STRCHAR*) tempcopy_local;
	if (sbytes > max_onstack)
		tempcopy = (STRCHAR*) STR_malloc (sbytes);
	STR_copy (tempcopy, GetString(), sbytes);
	do_upr ? STR_ToUpper(tempcopy, slen) : STR_ToLower(tempcopy, slen);
	// Compare it to the original; if a mismatch is found, reassign
	if (STR_memcmp(GetString(), tempcopy, sbytes))
		*this = tempcopy;
	if (tempcopy != (STRCHAR*) tempcopy_local)
		STR_free(tempcopy);
}


/*********************************************************************
* Proc:		Str::GetStringA
* Purpose:	Returns the string content as an ANSI string.  Result
*			should be freed with OS_free().
*********************************************************************/

#ifndef STR_NO_UNICODE
STRACHAR* Str::GetStringA() const
{
#ifdef STR_UNICODE
	// Compute how many bytes are needed, allocate buffer.
	int r = STR_ToAnsiBytes(dptr);
	if (r == 0)
		Error(SE_UnicodeError);
	STRACHAR* result = (STRACHAR*) OS_malloc (r);
	if (r == 1) {
		result[0] = 0;
		return result;
	}
	// Perform the actual conversion (will also append
	//   terminating zero as appropriate)
	STR_ToAnsi(dptr, result, r);
#else
	// No conversion: just do malloc/copy
	int r = GetLength() + 1;
	STRACHAR* result = (STRACHAR*) OS_malloc (r);
	STR_copy(result, GetString(), r);
#endif
	return result;
}
#endif


/*********************************************************************
* Proc:		Str::GetStringW
* Purpose:	Available only in non-UNICODE apps -- returns the string
*			content in UNICODE
*********************************************************************/

#ifndef STR_NO_UNICODE
STRWCHAR* Str::GetStringW() const
{
#ifdef STR_UNICODE
	// No conversion: just do malloc/copy
	int r = (GetLength() + 1) * sizeof(STRCHAR);
	STRWCHAR* result = (STRWCHAR*) OS_malloc (r);
	STR_copy(result, GetString(), r);
#else
	// Compute how many bytes are needed, allocate buffer.
	//   The function returns number of wide chars needed,
	//   _including_ the space for the NULL
	int r = STR_ToWideSlots(dptr);
	if (r == 0)
		Error(SE_UnicodeError);
	STRWCHAR* result = (STRWCHAR*) OS_malloc (r * sizeof(STRWCHAR));
	if (r == 1) {
		result[0] = 0;
		return result;
	}
	// Perform the actual conversion (will also append
	//   terminating zero as appropriate)
	STR_ToWide(dptr, result, r);
#endif
	return result;
}
#endif


/*********************************************************************
* Proc:		operator + (Str and Str, Str and LPCSTR)
*********************************************************************/

STR_EXPORT Str operator+(const Str& s1, const Str& s2)
{
	Str out;
	if (s1.IsMT()  ||  s2.IsMT())
		out.SetMT();
	out.Preallocate(CPOS(s1.GetLength() + s2.GetLength()));
	verifymt mtcheck(&out, -1);
	out.CoreAppendChars(s1.dptr, s1.GetLength(), mtcheck);
	out.CoreAppendChars (s2.dptr, s2.block()->m_Length, mtcheck);	// Bugfix
	return out;
}

STR_EXPORT Str operator+(const Str& s, const STRCHAR* lpsz)
{
	CPOS slen = STR_strlen(lpsz);
	Str out;
	if (s.IsMT())
		out.SetMT();
	out.Preallocate(s.GetLength() + slen);
	verifymt mtcheck(&out, -1);
	out.CoreAppendChars(s.dptr, s.GetLength(), mtcheck);
	out += lpsz;
	return out;
}

STR_EXPORT Str operator+(const Str& s, STRCHAR ch)
{
	Str out;
	if (s.IsMT())
		out.SetMT();
	out.Preallocate(s.GetLength() + 1);
	verifymt mtcheck(&out, -1);
	out.CoreAppendChars(s.dptr, s.GetLength(), mtcheck);
	out += ch;
	return out;
}


/*********************************************************************
* Proc:		Str::LoadString
* Purpose:	Loads a string from the stringtable.
* In:		resid - resource ID
* Out:		True if OK, FALSE if could not find such a string
*********************************************************************/

#if defined(STR_WIN32) &&  !defined(STR_NO_WINSTUFF)
BOOL Str::LoadString(UINT resid)
{
	verifymt mtcheck(this);
	HANDLE hLoad = STR_get_stringres();
	// Try smaller resources first
	STRCHAR buffer[96];
	if (::LoadString((HINSTANCE) hLoad, resid, 
		buffer, sizeof(buffer) / sizeof(STRCHAR)) == 0)
	{
		return FALSE;
	}
	CPOS slen = STR_strlen(buffer);
	if (slen < ((sizeof(buffer)/sizeof(STRCHAR))-1)) {
		Empty();
		CoreAppendChars(buffer, slen, mtcheck);
		return TRUE;
	}
	// We have a large string, use a big buffer
	const UINT s_big = 16384;
	STRCHAR* bigbuf = (STRCHAR*) STR_malloc (s_big * sizeof(STRCHAR));
	if (::LoadString((HINSTANCE) hLoad, resid, bigbuf, s_big) == 0) {
		STR_free(bigbuf);
		return FALSE;
	}
	*this = (const STRCHAR*) bigbuf;
	STR_free (bigbuf);
	return TRUE;
}
#endif


/*********************************************************************
* Proc:		AllocSysString, SetSysString
*********************************************************************/

#if defined(STR_WIN32)  &&  !defined(STR_NO_UNICODE) && !defined(STR_NO_WINSTUFF)

BSTR Str::AllocSysString() const
{
	if (IsEmpty())
		return NULL;		// Proper way to represent 0-character string
#ifdef STR_UNICODE
	const STRCHAR* widestr = dptr;
	int charlen = GetLength();
#else
	STRWCHAR* widestr = GetStringW();
	int charlen = (int) lstrlenW(widestr);
#endif
	STR_ASSERT(charlen > 0);
	BSTR result = ::SysAllocStringLen(NULL, charlen);
	if (result == NULL) {
		this->Error(SE_OutOfMemory);
		return NULL;		// Keep compiler happy
	}
	STR_copy(result, widestr, charlen*sizeof(WCHAR));
#ifndef STR_UNICODE
	Str::OS_free(widestr);
#endif
	return result;
}

BSTR Str::SetSysString(BSTR* str) const
{
#ifdef STR_UNICODE
	const STRCHAR* widestr = dptr;
	int charlen = GetLength();
#else
	STRWCHAR* widestr = GetStringW();
	int charlen = (int) lstrlenW(widestr);
#endif
	::SysReAllocStringLen(str, NULL, charlen);
	if (charlen)
		STR_copy(*str, widestr, charlen*sizeof(WCHAR));
#ifndef STR_UNICODE
	Str::OS_free(widestr);
#endif
	return *str;
}

const Str& Str::operator=(BSTR s)
{
	if (s)
		*this = (const STRWCHAR*) s;
	else
		Empty();
	return *this;
}

#endif


/*********************************************************************
* Proc:		Str::DeleteLeft, Delete
*********************************************************************/

Str& Str::DeleteLeft(CPOS count)
{
	verifymt mtcheck(this);
	if (GetLength() <= count) {
		Empty();
		return *this;
	}
	if (count == 0)
		return *this;
	syncstr keeper (mtcheck, 0, TRUE);		// Preserve buffer size
	STR_move (dptr, 
		((BYTE*) dptr)+(count*sizeof(STRCHAR)), 
		(GetLength()-count+1) * sizeof(STRCHAR));
	block()->m_Length -= count;
	return *this;
}

Str& Str::Delete(CPOS start, CPOS count /* = 1*/)
{
	verifymt mtcheck(this);
	if (GetLength() <= start) {
		Empty();
		return *this;
	}
	syncstr keeper (mtcheck, 0, TRUE);		// Preserve buffer size
	STRCHAR* pstart = dptr + start;
	if (GetLength() <= (start+count)) {
		pstart[0] = 0;
		block()->m_Length = start;
		return *this;
	}
	STR_move (pstart, 
		((BYTE*) pstart)+(count*sizeof(STRCHAR)), 
		(GetLength()-(start+count)+1) * sizeof(STRCHAR));
	block()->m_Length -= count;
	return *this;
}


/*********************************************************************
* Proc:		Str::FmtOneValue
* Purpose:	Helper for Str::Format, formats one %??? item
* In:		x - ptr to the '%' sign in the specification; on exit,
*				will point to the first char after the spec.
* Out:		True if OK, False if should end formatting (but also copy
*			the remaining characters at x)
*********************************************************************/

#ifndef STR_NO_RTLIBRARY
const char Str_dblletters[] = "eEfgG";
#endif

namespace StrImplementation {
static int CountEmbeddedStars(STRCHAR* fsbuf) {
	int cnt = 0;
	for (CPOS at=1; TRUE; at++) {
		if (fsbuf[at] == 0)
			return cnt;
		if (fsbuf[at] == '*')
			cnt++;
	}
}
}

BOOL Str::FmtOneValue (const STRCHAR*& x, verifymt& mtcheck, va_list& marker, BOOL rev_polarity)
{
	const STRACHAR fprefix[]    = "-+0 #*.123456789hlL";
	const STRACHAR intletters[] = "cCdiouxX";
	const STRACHAR nulltext[]   = "(null)";
	// Start copying format specifier to a local buffer
	const int chars_fsbuf = 64;
	STRCHAR fsbuf[chars_fsbuf];
	fsbuf[0] = '%';
	int fsp = 1;
GetMoreSpecifiers:
	// Get one character
#ifdef STR_DEBUG
	if (fsp >= chars_fsbuf) {
		Error(SE_BadFmtSpecifier);
		return FALSE;
	}
#endif
	STRCHAR ch = x[0];
	if (ch == 0)
		return FALSE;		// unexpected end of format string
	x++;
	// Chars that may exist in the format prefix
	if (STR_ansiscan (fprefix, ch) != NULL) {
		fsbuf[fsp] = ch;
		fsp++;
		goto GetMoreSpecifiers;
	}
	// 's' and S are the most important parameter specifier types
	if (ch == 's'  ||  ch == 'S') {
		// Determine whether we should fetch an ANSI or UNICODE parameter
#ifdef STR_UNICODE
		int takeansi = 0;
#else
		int takeansi = 1;
#endif
		if (ch == 'S')
			takeansi ^= 1;
		if (rev_polarity)
			takeansi ^= 1;
#ifdef STR_NO_UNICODE
		if (takeansi == 0)
			Error(SE_BadUnicodeUse);
#endif
		// Get a properly casted pointer to the actual string bits.
		//   Also handle the case with a null pointer.
		void* val = va_arg (marker, void*);
		const STRACHAR* val_ansi = NULL;
#ifndef STR_NO_UNICODE
		const STRWCHAR* val_uni  = NULL;
#endif
		if (val == NULL) {
			val_ansi = nulltext;
			takeansi = 1;
		}
		else {
			if (takeansi)
				val_ansi = (const STRACHAR*) val;
#ifndef STR_NO_UNICODE
			else
				val_uni = (const STRWCHAR*) val;
#endif
		}
		// Find out how many characters should we actually print.
		//   To do this, get the string length, but also try to
		//   detect a .precision field in the format specifier prefix.
		CPOS slen = 0;
		if (takeansi)
			slen = STR_astrlen(val_ansi);
#ifndef STR_NO_UNICODE
		else
			slen = STR_wstrlen(val_uni);
#endif
		fsbuf[fsp] = 0;
		const STRCHAR* precis = STR_scan (fsbuf, '.');
		if (precis != NULL  &&  precis[1] != 0) {
			// Convert value after dot, put within 0 and slen
			STRCHAR* endptr;
			int result = (int) STR_strtol (precis+1, &endptr, 10);
			if (result >= 0  &&  result < int(slen))
				slen = (UINT) result;
		}
		// Copy the appropriate number of characters
		if (slen > 0) {
#ifdef STR_UNICODE
#ifndef STR_NO_UNICODE
			if (takeansi) {
				Str temp;
				STRCHAR* temp_buf = temp.GetBuffer(slen+1);
				STR_ToWide_Counted(val_ansi, temp_buf, slen+1, slen);
				temp_buf[slen] = 0;
				temp.ReleaseBuffer(slen);
				CoreAppendChars (temp, slen, mtcheck);
			}
			else
#endif
				CoreAppendChars (val_uni, slen, mtcheck);
#else
#ifndef STR_NO_UNICODE
			if (takeansi)
#endif
				CoreAppendChars (val_ansi, slen, mtcheck);
#ifndef STR_NO_UNICODE
			else {
				Str temp;
				STRCHAR* temp_buf = temp.GetBuffer(slen+1);
				STR_ToAnsi_Counted(val_uni, temp_buf, slen+1, slen);
				temp_buf[slen] = 0;
				temp.ReleaseBuffer(slen);
				CoreAppendChars (temp, slen, mtcheck);
			}
#endif
#endif
		}
		return TRUE;
	}
	// '!' is our private extension, allows direct passing of Str*
	if (ch == '!') {
		// No precision characters taken into account here.
		const Str* value = va_arg (marker, const Str*);
		*this += *value;
		return TRUE;
	}
	// Chars that format an integer value
	if (STR_ansiscan (intletters, ch) != NULL) {
		fsbuf[fsp] = ch;
		fsbuf[fsp+1] = 0;
		STRCHAR valbuf[64];
		int value = va_arg (marker, int);
		STR_sprintf_i (valbuf, sizeof(valbuf)/sizeof(STRCHAR), fsbuf, value);
		*this += valbuf;
		return TRUE;
	};
	// Chars that format a double value
#ifndef STR_NO_RTLIBRARY
	if (STR_ansiscan (Str_dblletters, ch) != NULL) {
		fsbuf[fsp] = ch;
		fsbuf[fsp+1] = 0;
		int cntS = CountEmbeddedStars(fsbuf);
		STRCHAR valbuf[128];
		if (cntS == 0) {
			double value = va_arg (marker, double);
			STR_dosprintf_X(valbuf, sizeof(valbuf)/sizeof(STRCHAR), fsbuf, value);
		}
		else if (cntS == 1) {
			int digs1 = va_arg (marker, int);
			double value = va_arg (marker, double);
			STR_dosprintf_X(valbuf, sizeof(valbuf)/sizeof(STRCHAR), fsbuf, digs1, value);
		}
		else if (cntS == 2) {
			int digs1 = va_arg (marker, int);
			int digs2 = va_arg (marker, int);
			double value = va_arg (marker, double);
			STR_dosprintf_X(valbuf, sizeof(valbuf)/sizeof(STRCHAR), fsbuf, digs1, digs2, value);
		}
		else {
			Error(SE_BadFmtSpecifier);
			return FALSE;
		}
		*this += valbuf;
		return TRUE;
	}
#endif
	// 'Print pointer' is supported
	if (ch == 'p') {
		fsbuf[fsp] = ch;
		fsbuf[fsp+1] = 0;
		STRCHAR valbuf[64];
		void* value = va_arg (marker, void*);
		STR_sprintf_p (valbuf, sizeof(valbuf)/sizeof(STRCHAR), fsbuf, value);
		*this += valbuf;
		return TRUE;
	};
	// 'store # written so far' is obscure and unsupported
	if (ch == 'n') {
		Error(SE_UnsupFmtSpecifier);
#if defined(STR_NO_EXCEPTIONS)
		return FALSE;
#endif
	}
	// If we fall here, the character does not represent an item
	AppendChar (ch);
	return TRUE;
}

void Str::FormatCore (const STRCHAR* x, verifymt& mtcheck, va_list& marker, BOOL rev_polarity)
{
	for (;;) {
		// Locate next % sign, copy chunk, and exit if no more
		const STRCHAR* next_p = STR_scan (x, '%');
		if (next_p == NULL)
			break;
		if (next_p > x)
			CoreAppendChars (x, (CPOS) (next_p-x), mtcheck);
		x = next_p+1;
		// We're at a parameter
		if (!FmtOneValue (x, mtcheck, marker, rev_polarity))
			break;		// Copy rest of string and return
	}
	if (x[0] != 0)
		*this += x;
}


/*********************************************************************
* Proc:		Str::Mid, Left, Right
*********************************************************************/

Str Str::Mid (CPOS start, CPOS chars) const
{
	Str result;
	verifymt mtcheck(&result, -1);
	// Nothing to return?
	CPOS l = GetLength();
	if (l == 0  ||  (start+chars) == 0)
		return result;
	// Do not return data beyond the end of the string
	if (start >= l)
		return result;
	if ((start+chars) >= l)
		chars = l - start;
	// Copy bytes
	result.CoreAppendChars((this->operator const STRCHAR*())+start, chars, mtcheck);
	return result;
}

Str Str::Right (CPOS chars) const
{
	Str result;
	if (chars >= GetLength()) {
		result = *this;
		return result;
	}
	result = Mid(GetLength()-chars, chars);
	return result;
}


/*********************************************************************
* Proc:		Str::DeleteRight
*********************************************************************/

Str& Str::DeleteRight(CPOS count)
{
#ifdef STR_THREAD_SAFE
	verifymt mtcheck (this);
#endif
	if (GetLength() <= count)
		Empty();
	else
		TruncateAt (GetLength() - count);
	return *this;
}


/*********************************************************************
* Proc:		Str::Format
* Purpose:	sprintf-like method
*********************************************************************/

Str& Str::Format(const STRACHAR* fmt, ...)
{
	verifymt mtcheck(this);
	va_list marker;
	va_start(marker, fmt);
#ifdef STR_UNICODE
	Str temp_fmt(fmt);
	const STRCHAR* fmt1 = temp_fmt;
	BOOL rev_polarity = TRUE;
#else
	const STRCHAR* fmt1 = fmt;
	BOOL rev_polarity = FALSE;
#endif
	// Provide suitable minimum buffer size, erase string
	FmtPriorBuffer(fmt1, mtcheck);
	// Walk the string
	FormatCore (fmt1, mtcheck, marker, rev_polarity);
	va_end(marker);
	return *this;
}

#ifndef STR_NO_UNICODE
Str& Str::Format(const STRWCHAR* fmt, ...)
{
	verifymt mtcheck(this);
	va_list marker;
	va_start(marker, fmt);
#ifdef STR_UNICODE
	const STRCHAR* fmt1 = fmt;
	BOOL rev_polarity = FALSE;
#else
	Str temp_fmt(fmt);
	const STRCHAR* fmt1 = temp_fmt;
	BOOL rev_polarity = TRUE;
#endif
	// Provide suitable minimum buffer size, erase string
	FmtPriorBuffer(fmt1, mtcheck);
	// Walk the string
	FormatCore (fmt1, mtcheck, marker, rev_polarity);
	va_end(marker);
	return *this;
}
#endif

void Str::FormatVaArg(const STRCHAR *fmt, va_list& marker, BOOL rev_polarity)
{
	verifymt mtcheck(this);
	FmtPriorBuffer(fmt, mtcheck);
	FormatCore (fmt, mtcheck, marker, rev_polarity);
}

Str& Str::AppendFormat(const STRCHAR* fmt, ...)
{
	verifymt mtcheck(this);
	// Do not attempt to expand string buffer before
	//   doing the formatting; do not erase existing content
	va_list marker;
	va_start(marker, fmt);
	FormatCore (fmt, mtcheck, marker, FALSE);
	va_end(marker);
	return *this;
}

#if defined(STR_WIN32)  &&  !defined(STR_NO_WINSTUFF)

Str& Str::FormatRes(UINT resid, ...)
{
	verifymt mtcheck(this);
	// Get a copy of the format specifier
	Str formatter;
	if (!formatter.LoadString(resid)) {
		Empty();
		return *this;
	}
	// Provide suitable minimum buffer size, erase string
	FmtPriorBuffer(formatter, mtcheck);
	// Walk the string
	va_list marker;
	va_start(marker, resid);
	FormatCore (formatter, mtcheck, marker, FALSE);
	va_end(marker);
	return *this;
}
#endif


/*********************************************************************
* Proc:		Str::UrlEncode / UrlDecode
* Purpose:	Encode/decode a URL string in a format suitable for
*			sending over the Internet
*********************************************************************/

#ifndef STR_NO_RTLIBRARY

inline bool UrlCharNeedsEscaping(STRACHAR ch)
{
	if ((ch >= 'A'  &&  ch <= 'Z')  ||  (ch >= 'a'  &&  ch <= 'z')  ||  (ch >= '0'  &&  ch <= '9'))
		return false;		// No special encoding
	if (ch == '-'  ||  ch == '.'  ||  ch == '_'  ||  ch == '~')
		return false;		// "Unreserved" characters are not encoded 
	// All other characters need encoding.  These may be either characters that
	// are "reserved" - such as + & { } etc. - or characters that simply do not
	// have a good ASCII representation -- all characters <32 and >=127 fall in
	// this category.  Also note that we return true for the space character
	// (codepoint 32).  Our caller needs to decide whether to encode the space
	// as +, or as %20
	return true;
}

inline STRCHAR UrlHexDig(unsigned int value) 
{
	unsigned int vv = value & 0x0F;
	if (vv < 10)
		return ((STRCHAR) '0') + (STRCHAR) vv;
	else
		return ((STRCHAR) 'A') + (STRCHAR) vv - 10;
}

void Str::UrlEncode(Str& dest, BOOL spaceAsPlus /*= FALSE*/) 
{
	dest.Empty();
	dest.Preallocate(GetLength());		// Minimum possible buffer size, may grow
	verifymt mtdest(&dest);
	CPOS at;

#ifdef STR_UNICODE
	// Do a quick scan for characters above 127.  If at least one is found, we need
	// to jump through hoops a little bit -- convert to UTF8 first
	STRACHAR* utfs = NULL;
	for (at=0; TRUE; at++) {
		STRCHAR ch = dptr[at];
		if (ch == 0)
			break;
		if ((unsigned int) ch > 127) {
			utfs = ToUTF8(NULL, -1);
			if (utfs == NULL)
				Error(SE_OutOfMemory);
			break;
		}
	}
#endif

	for (at=0; TRUE; at++) {
		STRACHAR ch;
#ifdef STR_UNICODE
		if (utfs != NULL)
			ch = utfs[at];
		else
			ch = (STRACHAR) dptr[at];		// Guaranteed to fit in STRACHAR range
#else
		ch = dptr[at];
#endif
		if (ch == 0)
			break;
		if (spaceAsPlus  &&  (ch == ' '))
			ch = '+';
		else if (UrlCharNeedsEscaping(ch)) {
			STRCHAR tocat[3];
			tocat[0] = '%';
			tocat[1] = UrlHexDig((unsigned int) ch >> 4);
			tocat[2] = UrlHexDig((unsigned int) ch);
			dest.CoreAppendChars(tocat, 3, mtdest);
			continue;
		}
		STRCHAR onechar = (STRCHAR) ch;
		dest.CoreAppendChars(&onechar, 1, mtdest);
	}

#ifdef STR_UNICODE
	if (utfs != NULL)
		OS_free(utfs);
#endif
}

#endif	//ndef STR_NO_RTLIBRARY


/*********************************************************************
* Proc:		Str::GetWindowText
* Purpose:	Helper for Windows apps
*********************************************************************/

#ifndef STR_NO_WINSTUFF
void Str::GetWindowText (HWND wnd)
{
	verifymt mtcheck(this);
	Empty();
	// How big a buffer should we have?
	CPOS tl = (CPOS) ::SendMessage (wnd, WM_GETTEXTLENGTH, 0, 0L);
	if (tl == 0)
		return;
	CPOS tl_alloc = tl+2;
	up_alloc (tl_alloc);
	syncstr keeper (mtcheck, tl_alloc, TRUE);
	// And get the text
	::SendMessage (wnd, WM_GETTEXT, (WPARAM) tl+1, (LPARAM) dptr);
	block()->m_Length = tl;
}
#endif


/*********************************************************************
* Proc:		Str::TrimLeft
* Purpose:	Remove leading characters from a Str.  All characters
*			to be excluded are passed as a parameter; NULL means
*			'truncate spaces, tabs and CR/LF'
*********************************************************************/

namespace StrImplementation {
static STRCHAR* trim_defaultcharset = _T(" \t\r\n");
}

void Str::ImpTrimLeft(const STRCHAR* charset)
{
	CPOS good = 0;
	if (charset == NULL)
		charset = trim_defaultcharset;
	// Bugfix Apr 2 / 05: would ASSERT when the loop below reached the
	//   terminating null character
	while (dptr[good] != 0  &&  STR_scan (charset, dptr[good]) != NULL)
		good++;
	if (good > 0)
		DeleteLeft (good);
}

Str& Str::TrimLeft(const STRCHAR* charset /*= NULL*/)
{
#ifdef STR_THREAD_SAFE
	verifymt mtcheck(this);
#endif
	ImpTrimLeft(charset);
	return *this;
}


/*********************************************************************
* Proc:		Str::TrimRight
* Purpose:	Remove trailing characters; see TrimLeft
*********************************************************************/

void Str::ImpTrimRight(const STRCHAR* charset /*= NULL*/)
{
	CPOS good = block()->m_Length;
	if (good == 0)
		return;
	if (charset == NULL)
		charset = trim_defaultcharset;
	while (good > 0  &&  STR_scan (charset, dptr[good-1]) != NULL)
		--good;
	TruncateAt (good);		// Also works well with good == 0
}

Str& Str::TrimRight(const STRCHAR* charset /*= NULL*/)
{
#ifdef STR_THREAD_SAFE
	verifymt mtcheck(this);
#endif
	ImpTrimRight(charset);
	return *this;
}

Str& Str::Trim(const STRCHAR* charset /*= NULL*/)
{
#ifdef STR_THREAD_SAFE
	verifymt mtcheck(this);
#endif
	ImpTrimRight(charset);
	ImpTrimLeft(charset);
	return *this;
}


/*********************************************************************
* Proc:		Str::Remove
* Purpose:	Remove all occurences of a character from a string
*			(case-sensitive compare). Return # of chars removed.
*********************************************************************/

int Str::Remove(STRCHAR ch)
{
	verifymt mtcheck(this);
	if (ch == 0)
		Error(SE_NoNullChar);
	STRCHAR* next_at = STR_scan(dptr, ch);
	if (!next_at)
		return 0;
	CPOS start = (CPOS) (next_at - dptr);
	int count = 1;
	syncstr keeper (mtcheck, 0, TRUE);
	next_at = dptr + start;		// Need to reconstruct because of keeper{}
	for (;;) {
		// Cut it out here
		start = (CPOS) (next_at - dptr);
		STR_move (next_at,
			((BYTE*) next_at)+sizeof(STRCHAR), 
			(GetLength()-(start+1)+1) * sizeof(STRCHAR));
		block()->m_Length -= 1;
		// Find next occurence
		next_at = STR_scan(next_at, ch);
		if (!next_at)
			return count;
		count++;
	}
}


/*********************************************************************
* Proc:		Str::Find
* Purpose:	Scan the string for a particular substring; return the 
*			index where found first, or -1 if cannot be found
*********************************************************************/

CPOS Str::Find (const STRCHAR* substr, CPOS startat /*= 0*/) const
{
	int s_len = STR_strlen(substr);
	if (s_len == 0)
		return (CPOS) -1;
	int s_at = (int) startat;
	int e_at = (int) GetLength() - s_len;
	while (s_at <= e_at) {
		// Find next occurrence of substr first character
		const STRCHAR* charscan = STR_scan (dptr+s_at, substr[0]);
		if (charscan == NULL)
			return (CPOS) -1;
		s_at = (int) (charscan - dptr);
		// Single-character string?
		if (s_len > 1) {
			// Try to match at this location
			const BYTE* b1 = ((const BYTE*) GetString()) + ((s_at+1) * sizeof(STRCHAR));
			const BYTE* b2 = ((const BYTE*) substr) + sizeof(STRCHAR);
			if (STR_memcmp(b1, b2, (s_len-1)*sizeof(STRCHAR)) == 0)
				return (CPOS) s_at;
		}
		else
			return (CPOS) s_at;
		// Nope, go on
		s_at++;
	}
	return (CPOS) -1;
}

CPOS Str::FindNoCase (const STRCHAR* substr, CPOS startat /*= 0*/) const
{
	Str temp_obj (*this);
	temp_obj.MakeLower();
	Str temp_str (substr);
	temp_str.MakeLower();
	return temp_obj.Find(temp_str, startat);
}


/*********************************************************************
* Proc:		Str::ReplaceAt
* Purpose:	Given a starting position and a number of characters,
*			replaces that portion of the string with another string
*********************************************************************/

void Str::ReplaceAt(CPOS pos, CPOS numchars, const STRCHAR* newstr)
{
	verifymt mtcheck(this);
	if (numchars == 0)
		Error(SE_BadScanPos);
	if ((pos + numchars) > GetLength())
		Error(SE_BadScanPos);
	int repl_len = STR_strlen(newstr);
	if (repl_len == numchars) {
		// Replace "in situ"
		syncstr keeper (mtcheck, 0, TRUE);
		BYTE* dest = ((BYTE*) dptr) + (pos * sizeof(STRCHAR));
		STR_copy(dest, newstr, repl_len*sizeof(STRCHAR));
	}
	else {
		Delete(pos, numchars);
		Insert(newstr, pos, repl_len);
	}
}


/*********************************************************************
* Proc:		Str::Replace
* Purpose:	Replaces all occurences of a given substring with another
*			substring
*********************************************************************/

STR_EXPORT int StrImplementation::STR_ReplaceCore
	(Str& obj, const STRCHAR* oldstr, const STRCHAR* newstr, BOOL nocase)
{
#ifdef STR_THREAD_SAFE
	verifymt mtcheck(&obj);
#endif
	int count = 0;
	int olds_len = STR_strlen(oldstr);
	if (olds_len == 0)
		obj.Error(Str::SE_StringEmpty);
	int news_len = STR_strlen(newstr);
	// Bugfix July 5, 2002: was look_at < maxl_at
	// Bugfix Aug 3, 2003: if the replacement string was longer than the source one,
	//   some of the last occurrences would not be found due to a wrong look_at <= maxl_at
	//   check (which would be correct only if <oldstr> string length was <= <newstr>
	//   length
	for (CPOS look_at = 0;; ) {
		// End of scan?
		if (look_at > (obj.GetLength() - olds_len))
			break;
		// Find next occurrence
		look_at = nocase ? obj.FindNoCase(oldstr, look_at) : obj.Find(oldstr, look_at);
		if (look_at == (CPOS) -1)
			return count;
		// Replace here
		obj.ReplaceAt(look_at, olds_len, newstr);
		count++;
		// Advance 'scan' position to skip just replaced string :)
		look_at += news_len;
	}
	return count;
}


/*********************************************************************
* Proc:		Str::EndInSlash
* Purpose:	Makes sure the string has a \ or / (depending on the OS)
*			at its end, unless it is completely empty.  Useful when
*			forming path specifications.
*********************************************************************/

Str& Str::EndInSlash()
{
#ifdef STR_THREAD_SAFE
	verifymt mtcheck(this);
#endif
#ifdef STR_LINUX
	const STRCHAR mustbe = '/';
#else
	const STRCHAR mustbe = '\\';
#endif
	if (!IsEmpty()) {
		if (GetLastChar() != mustbe)
			AppendChar(mustbe);
	}
	return *this;
}


/*********************************************************************
* Proc:		Str::StringSpanIncluding (static)
* Purpose:	Find the first character in the string that is NOT one
*			of the tokens passed. May return the position of the 
*			terminating null (complete match, or empty string if 0)
*********************************************************************/

/*static*/ CPOS Str::StringSpanIncluding(const STRCHAR* string, const STRCHAR* tokens)
{
	int tklen = STR_strlen(tokens);
	int pos;
	if (tklen == 1) {
		// Optimize for one token character
		pos = 0;
		while (string[pos] == tokens[0])
			pos++;
	}
	else {
		// Several token characters: loop char by char
		pos = 0;
		while (string[pos] != 0) {
			const STRCHAR* found = STR_scan2(tokens, string[pos], tklen);
			if (!found)
				return (CPOS) pos;
			pos++;
		}
	}
	return (CPOS) pos;
}


/*********************************************************************
* Proc:		Str::StringSpanExcluding (static)
* Purpose:	Find the first character in the string that is one
*			of the tokens passed. Similar to FindOneOf.  May return
*			-1 if a match was never found.
*********************************************************************/

/*static*/ CPOS Str::StringSpanExcluding(const STRCHAR* string, const STRCHAR* tokens)
{
	CPOS minpos = (CPOS) -1;
	for (int idx=0; tokens[idx] != 0; idx++) {
		STRCHAR ch = tokens[idx];
		const STRCHAR* scan = STR_scan (string, ch);
		if (scan) {
			CPOS atpos = (CPOS) (scan - string);
			if (atpos < minpos  ||  minpos == (CPOS) -1)
				minpos = atpos;
		}
	}
	return minpos;
}


/*********************************************************************
* Proc:		Str::SpanIncluding
* Purpose:	Returns the beginning of the string, up to (but excluding)
*			the first character NOT in the passed tokens list
*********************************************************************/

Str Str::SpanIncluding(const STRCHAR* tokens) const
{
	STR_ASSERT(STR_strlen(tokens) > 0);
	CPOS endat = StringSpanIncluding(dptr, tokens);
	return (endat > 0) ? Left(endat) : REFER_EmptyStr;
}


/*********************************************************************
* Proc:		Str::SpanExcluding
* Purpose:	Returns the beginning of the string, up to (but excluding)
*			the first character matching one of the passed tokens
*********************************************************************/

Str Str::SpanExcluding(const STRCHAR* tokens) const
{
	STR_ASSERT(STR_strlen(tokens) > 0);
	CPOS endat = StringSpanExcluding(dptr, tokens);
	return (endat > 0) ? Left(endat) : REFER_EmptyStr;
}


/*********************************************************************
* Proc:		Str::Tokenize
* Purpose:	Similar to strtok()
*********************************************************************/

Str Str::Tokenize(const STRCHAR* tokens, CPOS& startat) const
{
	STR_ASSERT(STR_strlen(tokens) > 0);
	// Nothing more to do?
	if (startat < 0)
		return REFER_EmptyStr;
	const STRCHAR* scanStart = dptr + startat;
	int slength = GetLength();
	if (scanStart < (dptr+slength)) {
		// Skip all characters in the tokens list
		CPOS nSkip = StringSpanIncluding(scanStart, tokens);
		if (nSkip >= 0  &&  scanStart[nSkip] != 0) {
			scanStart += nSkip;
			CPOS nEnd = StringSpanExcluding(scanStart, tokens);
			if (nEnd < 0)
				nEnd = (CPOS) STR_strlen(scanStart);
			int i_from = startat + nSkip;
			int i_to   = (i_from + nEnd) - 1;
			startat = i_to + 1;
			if (startat >= slength)
				startat = (CPOS) -1;
			return Mid(i_from, (i_to-i_from)+1);
		}
	}
	// Not found, return an empty string
	startat = (CPOS) -1;
	return REFER_EmptyStr;
}


/*********************************************************************
* Proc:		GetBuffer / ReleaseBuffer family
*********************************************************************/

STRCHAR* Str::GetBuffer()
{
	verifymt mtcheck(this);
	syncstr keeper (mtcheck, 0, TRUE);
	return dptr;
}

STRCHAR* Str::GetBufferSetLength(CPOS newlen)
{
	STRCHAR* buf = GetBuffer(newlen);
	STR_ASSERT(buf == dptr);
	block()->m_Length = newlen;
	buf[newlen] = 0;
	return buf;
}

void Str::ReleaseBuffer(CPOS newlen /*= -1*/)
{
	// Fix (or maybe "improvement" is better) 2.2.1 (5-Jul-04) ::
	//     always null-terminate, even when newlen != -1
	if (newlen == (CPOS) -1)
		newlen = STR_strlen(dptr);
	STR_ASSERT(newlen <= block()->m_Alloc);
	block()->m_Length = newlen;
	dptr[newlen] = 0;
}


#ifndef STR_NO_RTLIBRARY

Str& Str::FormatVA(const STRACHAR *fmt, va_list marker) 
{
	verifymt mtcheck(this);
#ifdef STR_UNICODE
	Str temp_fmt(fmt);
	const STRCHAR* fmt1 = temp_fmt;
	BOOL rev_polarity = TRUE;
#else
	const STRCHAR* fmt1 = fmt;
	BOOL rev_polarity = FALSE;
#endif
	// Provide suitable minimum buffer size, erase string
	FmtPriorBuffer(fmt1, mtcheck);
	// Walk the string
	FormatCore (fmt1, mtcheck, marker, rev_polarity);
	return *this;
}

#ifndef STR_NO_UNICODE
Str& Str::FormatVA(const STRWCHAR *fmt, va_list marker)
{
	verifymt mtcheck(this);
#ifdef STR_UNICODE
	const STRCHAR* fmt1 = fmt;
	BOOL rev_polarity = FALSE;
#else
	Str temp_fmt(fmt);
	const STRCHAR* fmt1 = temp_fmt;
	BOOL rev_polarity = TRUE;
#endif
	// Provide suitable minimum buffer size, erase string
	FmtPriorBuffer(fmt1, mtcheck);
	// Walk the string
	FormatCore (fmt1, mtcheck, marker, rev_polarity);
	return *this;
}

#endif

#endif


/*********************************************************************
* Proc:		Managed C++ Features
*********************************************************************/

#ifdef STR_MCPP_FEATURES

Str::Str(System::String __gc* source)
{
	dptr = null_str();
	const __wchar_t __pin* s_str = PtrToStringChars(source);
#ifdef STR_UNICODE
	AppendChars(s_str);
#else
	AppendString(s_str);
#endif
}

const Str& Str::operator=(System::String __gc* source)
{
	Empty();
	const __wchar_t __pin* s_str = PtrToStringChars(source);
#ifdef STR_UNICODE
	AppendChars(s_str);
#else
	AppendString(s_str);
#endif
	return *this;
}

#endif


/*********************************************************************
* Proc:		CArchive friend functions
*********************************************************************/

#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)

void Str::operator +=(const CString& obj)
{
	const STRCHAR* objText = (const STRCHAR*) obj;
	*this += objText;
}

void StrImplementation::StmWriteStringLength(CArchive& ar, int len, BOOL is_uni)
{
	if (is_uni)
	{
		// Tag Unicode strings
		ar << (BYTE) 0xff;
		ar << (WORD) 0xfffe;
	}

	if (len < 255)
	{
		ar << (BYTE) len;
	}
	else if (len < 0xfffe)
	{
		ar << (BYTE) 0xff;
		ar << (WORD) len;
	}
	else {
		ar << (BYTE) 0xff;
		ar << (WORD) 0xffff;
		ar << (DWORD) len;
	}
}

int StrImplementation::StmReadStringLength(CArchive& ar, int& charsize)
{
	DWORD dwLength;
	WORD wLength;
	BYTE bLength;

	charsize = sizeof(STRACHAR);

	ar >> bLength;
	if (bLength < 0xff)
		return bLength;

	ar >> wLength;
	if (wLength == 0xfffe)
	{
#ifdef STR_NO_UNICODE
		Str::ErrorNoObject(Str::SE_BadUnicodeUse);
#else
		charsize = sizeof(STRWCHAR);

		ar >> bLength;
		if (bLength < 0xff)
			return bLength;

		ar >> wLength;
#endif
	}
	if (wLength < 0xffff)
		return wLength;

	ar >> dwLength;
	return dwLength;
}

/*static*/ void Str::WriteStringLP(CArchive& ar, const STRCHAR* lpsz)
{
	int len = STR_strlen(lpsz);
	BOOL use_uni;
#ifdef STR_NO_UNICODE
	use_uni = TRUE;
#else
	use_uni = sizeof(STRCHAR) == sizeof(STRWCHAR);
#endif
	StmWriteStringLength(ar, len, use_uni);
	ar.Write(lpsz, len * sizeof(STRCHAR));
}

/*static*/ BOOL Str::ReadStringLP(CArchive& ar, Str& dest)
{
	int charsize;  // 1 = STRACHAR, 2 = STRWCHAR
	int len = (int) StmReadStringLength(ar, charsize);
	int lenbytes = len * charsize;
	int oursize = sizeof(STRCHAR);
	if (charsize == oursize) {
		// Matching character set; use internal buffer
		STRCHAR* buf = dest.GetBuffer(len);
		int read = ar.Read(buf, lenbytes);
		if (read != lenbytes)
			AfxThrowArchiveException(CArchiveException::endOfFile);
		dest.ReleaseBuffer(len);
	}
	else {
#ifdef STR_NO_UNICODE
		Str::ErrorNoObject(Str::SE_BadUnicodeUse);
		return FALSE;
#else
		// First, load source data in a temporary buffer
		BYTE* buf = (BYTE*) OS_malloc(lenbytes + charsize);
		int read = ar.Read(buf, lenbytes);
		if (read != lenbytes)
			AfxThrowArchiveException(CArchiveException::endOfFile);
		buf[lenbytes] = 0;
		if (charsize == sizeof(STRWCHAR))
			buf[lenbytes + 1] = 0;
		// Then use assignment operator to perform proper conversion
		if (charsize == sizeof(STRACHAR))
			dest = (const STRACHAR*) buf;
		else
			dest = (const STRWCHAR*) buf;
		OS_free(buf);
#endif
	}

	return TRUE;
}

/*static*/ void Str::WriteString(CArchive& ar, const STRCHAR* lpsz)
{
	ar.Write(lpsz, StrImplementation::STR_strlen(lpsz) * sizeof(STRCHAR));
}

/*static*/ BOOL Str::ReadString(CArchive& ar, Str& dest)
{
	STRCHAR buf[64];
	int gotbuf = 0;
	BOOL chunk1 = TRUE;
	dest.Empty();
	for (;;) {
		STRCHAR ch;
		if (ar.Read(&ch, sizeof(ch)) != sizeof(ch))
			break;
		// Stop at end of line (\r is discarded, if present)
		if (ch == '\n'  ||  ch == '\r') {
			if (ch == '\r')
				if (ar.Read(&ch, sizeof(ch)) != sizeof(ch))	// Go past \n also
					break;
			goto GotIt;
		}
		buf[gotbuf++] = ch;
		if (gotbuf == (sizeof(buf) / sizeof(STRCHAR))) {
			dest.AppendChars(buf, gotbuf);
			gotbuf = 0;
			chunk1 = FALSE;
		}
	}
	// We come here if EOF condition is found
	if (gotbuf == 0  &&  chunk1)
		return FALSE;
GotIt:
	if (gotbuf > 0)
		dest.AppendChars(buf, gotbuf);

	return TRUE;
}

#endif


/*********************************************************************
* Proc:		Str::ToUTF8 / Str::FromUTF8
*********************************************************************/

#if defined(STR_UNICODE)  &&  !defined(STR_NO_RTLIBRARY)

#ifdef STR_WIN32
#pragma warning (disable:4244)
#endif

#define BIN_00000001 0x01
#define BIN_00000011 0x03
#define BIN_00000111 0x07
#define BIN_00001111 0x0F
#define BIN_00011111 0x1F
#define BIN_00111111 0x3F
#define BIN_10000000 0x80
#define BIN_11000000 0xC0
#define BIN_11100000 0xE0
#define BIN_11110000 0xF0
#define BIN_11111000 0xF8
#define BIN_11111100 0xFC
#define BIN_11111110 0xFE

inline int UtfMakeBytes(UtfChType ch, BYTE* toadd)
{
	if (ch <= 0x7F) {
		toadd[0] = (BYTE) ch;
		return 1;
	}
	if (ch <= 0x7FF) {
		toadd[0] = (BYTE) (ch >> 6) | BIN_11000000;
		toadd[1] = (BYTE) (ch & BIN_00111111) | BIN_10000000;
		return 2;
	}
	if (ch <= 0xFFFF) {
		toadd[0] = (BYTE) (ch >> 12) | BIN_11100000;
		toadd[1] = ((BYTE) (ch >> 6) & BIN_00111111) | BIN_10000000;
		toadd[2] = (BYTE) (ch & BIN_00111111) | BIN_10000000;
		return 3;
	}
#ifdef STR_WCIS16BITS
	// Should never be able to come here...
	return -1;
#else
	if (ch <= 0x1FFFFF) {
		toadd[0] = (BYTE) (ch >> 18) | BIN_11110000;
		toadd[1] = ((BYTE) (ch >> 12) & BIN_00111111) | BIN_10000000;
		toadd[2] = ((BYTE) (ch >>  6) & BIN_00111111) | BIN_10000000;
		toadd[3] = (BYTE) (ch & BIN_00111111) | BIN_10000000;
		return 4;
	}
	if (ch <= 0x3FFFFFF) {
		toadd[0] = (BYTE) (ch >> 24) | BIN_11111000;
		toadd[1] = ((BYTE) (ch >> 18) & BIN_00111111) | BIN_10000000;
		toadd[2] = ((BYTE) (ch >> 12) & BIN_00111111) | BIN_10000000;
		toadd[3] = ((BYTE) (ch >>  6) & BIN_00111111) | BIN_10000000;
		toadd[4] = (BYTE) (ch & BIN_00111111) | BIN_10000000;
		return 5;
	}
	toadd[0] = (BYTE) (ch >> 30) | BIN_11111100;
	toadd[1] = ((BYTE) (ch >> 24) & BIN_00111111) | BIN_10000000;
	toadd[2] = ((BYTE) (ch >> 18) & BIN_00111111) | BIN_10000000;
	toadd[3] = ((BYTE) (ch >> 12) & BIN_00111111) | BIN_10000000;
	toadd[4] = ((BYTE) (ch >>  6) & BIN_00111111) | BIN_10000000;
	toadd[5] = (BYTE) (ch & BIN_00111111) | BIN_10000000;
	return 6;
#endif
}

inline bool UtfNonEncodable(UtfChType ch)
{
    return (ch >= 0xD800  &&  ch <= 0xDFFF)  ||  ch == 0xFFFE  ||  ch == 0xFFFF;
}

STRACHAR* Str::ToUTF8(STRACHAR* buffer, int buflen /*= -1*/) const
{
	// Need to maintain our own buffer?
	BOOL ownbuf = FALSE;
	CPOS ourlen = GetLength();
	if (buflen < 0) {
		STR_ASSERT(buffer == NULL);
		buflen = (int) (ourlen * 1.25);			// Empirical
		buflen = STR_max(buflen, 16);			// No point in allocating less
		buffer = (STRACHAR*) OS_malloc(buflen);
		ownbuf = TRUE;
	}
	else {
		STR_ASSERT(buffer != NULL);
	}
	// Start doing the conversion
	int inbuffer = 0;
	for (int i=0; i<ourlen; i++) {
		//
		// Prepare the representation of the current character
		STRACHAR toadd[6];
		UtfChType ch = (UtfChType) dptr[i];
		if (UtfNonEncodable(ch))
			Error(SE_MalformedUtf8Char);	// Surrogates not allowed in UTF-8 sequences 
		int bcount = UtfMakeBytes(ch, (BYTE*) toadd);
		if (bcount < 0)
			Error(SE_MalformedUtf8Char);	// Can't really happen

		//
		// Manage buffer size
		if ((inbuffer + bcount) >= buflen) {
			// Caller's buffer cannot be expanded
			if (!ownbuf)
				return (STRACHAR*) -1;
			// Do the expansion
            int newbuflen;
			if (buflen <= 32)
				newbuflen = 128;
			else
				newbuflen = (int) (buflen * 1.5);
			STRACHAR* newbuffer = (STRACHAR*) OS_malloc(newbuflen);
			if (inbuffer)
				memcpy(newbuffer, buffer, inbuffer);
			free(buffer);
			buffer = newbuffer;
			buflen = newbuflen;
		}

		// Append and continue to next character
		memcpy(buffer + inbuffer, toadd, bcount);
		inbuffer += bcount;
	}
	// Append terminating null and return
	buffer[inbuffer] = 0;
	return buffer;
}

inline bool GoodUtfFollowByte(BYTE b)
{
	return (b & BIN_11000000) == BIN_10000000;
}

inline int UtfAssembleUcs(const STRACHAR* buffer, UtfChType& result)
{
	// First character
	BYTE b0 = buffer[0];
	if (b0 == 0)
		return 0;			// Null character??
	if (b0 <= 0x7F) {
		result = (unsigned short) b0;
		return 1;
	}
	// Is it a two-byte sequence?
	BYTE b1 = buffer[1];
	if (!GoodUtfFollowByte(b1))
		return -1;
	if ((b0 & BIN_11100000) == BIN_11000000) {
		// We have two good bytes, assemble result
		result = (b1 & BIN_00111111);
		result |= (b0 & BIN_00011111) << 6;
		return 2;
	}
	// Is it a three-byte sequence? (BUGFIX: Jul 2003 - was checking the wrong UtfFollowByte)
	BYTE b2 = buffer[2];
	if (!GoodUtfFollowByte(b2))
		return -1;
	if ((b0 & BIN_11110000) == BIN_11100000) {
		// We have three good bytes, assemble result
		result = (b2 & BIN_00111111);
		result |= (b1 & BIN_00111111) << 6;
		result |= (b0 & BIN_00001111) << 12;
		return 3;
	}
#ifdef STR_WCIS16BITS
	// Longer sequences are not supported
	return -1;
#else
	// Is it a 4-byte sequence?
	BYTE b3 = buffer[3];
	if (!GoodUtfFollowByte(b3))
		return -1;
	if ((b0 & BIN_11111000) == BIN_11110000) {
		// We have 4 good bytes, assemble result
		result = (b3 & BIN_00111111);
		result |= (b2 & BIN_00111111) << 6;
		result |= (b1 & BIN_00111111) << 12;
		result |= (b0 & BIN_00000111) << 18;
		return 4;
	}
	// Is it a 5-byte sequence?
	BYTE b4 = buffer[4];
	if (!GoodUtfFollowByte(b4))
		return -1;
	if ((b0 & BIN_11111100) == BIN_11111000) {
		// We have 5 good bytes, assemble result
		result = (b4 & BIN_00111111);
		result |= (b3 & BIN_00111111) << 6;
		result |= (b2 & BIN_00111111) << 12;
		result |= (b1 & BIN_00111111) << 18;
		result |= (b0 & BIN_00000011) << 24;
		return 4;
	}
	// 6-byte sequence is what's left to check
	BYTE b5 = buffer[5];
	if (!GoodUtfFollowByte(b5))
		return -1;
	if ((b0 & BIN_11111110) == BIN_11111100) {
		// We have 6 good bytes, assemble result
		result = (b5 & BIN_00111111);
		result |= (b4 & BIN_00111111) << 6;
		result |= (b3 & BIN_00111111) << 12;
		result |= (b2 & BIN_00111111) << 18;
		result |= (b1 & BIN_00111111) << 24;
		result |= (b0 & BIN_00000001) << 30;
		return 4;
	}
	return -1;
#endif
}

void Str::FromUTF8(const STRACHAR* buffer, int length /*= -1*/)
{
	// Preparations
	if (length < 0)
		length = STR_astrlen(buffer);
	Empty();
	int preallocated = length;						// Empirical
	STRWCHAR* direct = GetBuffer(preallocated+1);	// Will also handle mtcheck()
	// Convert characters one by one
	int newlen = 0;
#ifndef STR_NO_EXCEPTIONS
	try {
#endif
		int total_in = 0;
		while (length > 0) {
			// Get one UTF-8 character and assemble a UCS-2 (-4) value
			UtfChType result = (UtfChType) 0;		// Keep VC6 happy
			int go_past = UtfAssembleUcs(buffer+total_in, result);
			if (go_past == 0)
				Error(SE_MalformedUtf8Char);	// Cannot have end of string in the middle!
			if (go_past < 0)
				Error(SE_MalformedUtf8Char);	// Bad sequence in the decoder
			// Decoded OK, but is this sequence allowed?
			if (UtfNonEncodable(result))
				Error(SE_MalformedUtf8Char);	// Surrogates not allowed in UTF-8 sequences 
			// Store the result and go on scanning
			direct[newlen] = (STRWCHAR) result;
			newlen++;
			STR_ASSERT(newlen <= preallocated);	// Should never happen
			total_in += go_past;
			length -= go_past;
			STR_ASSERT(length >= 0);
		}
#ifndef STR_NO_EXCEPTIONS
	}
	catch (...) {
		ReleaseBuffer(0);		// Make sure we clean up after ourselves
		throw;
	}
#endif		// We'll leak here but that's life
	// Done
	ReleaseBuffer(newlen);
}

#ifdef STR_WIN32
#pragma warning (default:4244)
#endif

#endif


/*********************************************************************
* Proc:	Soundex algorithm
*********************************************************************/

#ifndef STR_NO_RTLIBRARY

inline BOOL SoundexIsAlpha(STRCHAR ch)
{ return (ch >= 'a'  &&  ch <= 'z')  ||  (ch >= 'A'  &&  ch <= 'Z'); }

inline STRCHAR SoundexToUpper(STRCHAR ch)
{ return (ch >= 'a'  &&  ch <= 'z') ? (STRCHAR) (ch - ('a' - 'A')) : ch; }

Str Str::Soundex() const
{
    static int s_tab[] = {
        // A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
           0,1,2,3,0,1,2,0,0,2,2,4,5,5,0,1,2,6,2,3,0,1,0,2,0,2
    };
	// Skip leading non-alpha characters
	CPOS pos = 0;
	CPOS ourlen = GetLength();
	while (pos < ourlen  &&  !SoundexIsAlpha(dptr[pos]))
		pos++;
	if (pos == ourlen)
		return Str::EmptyStr;	// No meaningful characters
	// Get initial letter as first result item
	STRCHAR result[4+1];
	result[0] = SoundexToUpper(dptr[pos++]);
	int count = 1;
	while ((pos < ourlen)  &&  (count <= 3)) {
		STRCHAR ch = dptr[pos++];
		if (SoundexIsAlpha(ch)  &&  ch != dptr[pos-2]) {
			int s_val = s_tab[SoundexToUpper(ch) - 'A'];
			if (s_val)
				// Another result character
				result[count++] = (STRCHAR) ('0' + s_val);
		}
	}
	while (count <= 3)
		result[count++] = '0';
	// We're done
	result[4] = 0;
	return Str(result);
}

#endif


/*********************************************************************
* Proc:		Char::ToANSI, ToUnicode
*********************************************************************/

#ifndef STR_NO_UNICODE
/*static*/ STRACHAR Char::ToANSI(STRWCHAR ch)
{
	STR_ASSERT(ch == 0  ||  STR_ToAnsiBytes_Counted(&ch, 1) == 1);
	STRACHAR buffer;
    STR_ToAnsi_Counted(&ch, &buffer, 1, 1);
	return buffer;
}

/*static*/ STRWCHAR Char::ToUnicode(STRACHAR ch)
{
	STRWCHAR buffer;
	STR_ToWide_Counted(&ch, &buffer, 1, 1);
	return buffer;
}
#endif


#ifdef STR_BORCPPBUILDER
#undef Char
#endif
