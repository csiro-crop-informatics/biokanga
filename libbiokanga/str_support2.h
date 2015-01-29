// STR_Support2.h
//
// Version 2.4.0
//
// This file is a part of the Str Library, and is subject to the
// clauses of the License Agreement accompanying the software.

// Static methods to call OS-provided memory management

#ifndef STR_USER_ALLOCATOR

#ifdef STR_NO_RTLIBRARY

inline /*static*/ void* Str::OS_malloc (size_t size)
{
	return ::HeapAlloc(::GetProcessHeap(), HEAP_GENERATE_EXCEPTIONS, size); 
}

inline /*static*/ void Str::OS_free (void* mem)
{
	if (mem)
		::HeapFree(::GetProcessHeap(), 0, mem); 
}

#else

inline /*static*/ void* Str::OS_malloc (size_t size)
{ return malloc(size); }

inline /*static*/ void Str::OS_free (void* mem)
{ free(mem); }

#endif

#else

void* Str_OS_malloc_handler(size_t size);	// Must be supplied by user code
void Str_OS_free_handler(void* mem);		// Must be supplied by user code

inline /*static*/ void* Str::OS_malloc (size_t size)
{ return Str_OS_malloc_handler(size); }

inline /*static*/ void Str::OS_free (void* mem)
{ Str_OS_free_handler(mem); }

#endif


//
// Everything below is in a separate namespace

namespace StrImplementation {

#ifdef STR_DEBUG
StrImplementation::SBlock* STR_get_block_dbg(CPOS, int line, const char* file);
#define STR_get_block(x)		\
	StrImplementation::STR_get_block_dbg((x), __LINE__, __FILE__)
#else
#define STR_get_block(x)		\
	StrImplementation::STR_get_block_rls(x)
StrImplementation::SBlock* STR_get_block_rls(CPOS);
#endif

STR_EXPORT void     STR_rls_block(verifymt& mtcheck);
STR_EXPORT STRCHAR* STR_scan(const STRCHAR *s, STRCHAR c);
STR_EXPORT STRCHAR* STR_rscan(const STRCHAR *s, STRCHAR c);
STR_EXPORT STRCHAR* STR_scan2(const STRCHAR *s, STRCHAR c, int nchars);
STR_EXPORT STRCHAR* STR_CsCopyCore(Str*, const Str&);
STR_EXPORT void     STR_rls_chain(UINT cindex);

class STR_Initializer {
public:
	SBlock_Empty m_Empty;
public:
#ifndef STR_NO_RTLIBRARY
	STR_Initializer();
#else
	void InitializeStr();
#endif
};

extern STR_Initializer STR_EXPORT oInitializer;

#ifdef STR_THREAD_SAFE
STR_EXPORT STRCHAR* NewEmptyMT(CPOS prealloc);
#else
STR_EXPORT STRCHAR* NewEmptyMTForceNI(CPOS prealloc);
#endif

#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)
void STR_EXPORT StmWriteStringLength(CArchive& ar, int len, BOOL is_uni);
int  STR_EXPORT StmReadStringLength(CArchive& ar, int& charsize);
#endif

}


/*********************************************************************
* Very simple utility macros and Str inlines
*********************************************************************/

// Round up to the byte alignment desired (minimum 4, always power of 2)
// When using unicode, this returns divides the alignment by 2, because
//   unicode characters are two bytes long
inline void up_alloc (CPOS& value)
{
#ifdef STR_UNICODE
	value = CPOS(value | WIDE_BAE_MASK);
#else
	value = CPOS(value | ANSI_BAE_MASK);
#endif
}

// ASSERT that something is rounded up with up_alloc
#ifdef STR_DEBUG
inline void STR_ASSERT_UPPED(CPOS bytes)
{
#ifdef STR_UNICODE
	STR_ASSERT ((bytes & WIDE_BAE_MASK) == WIDE_BAE_MASK);
#else
	STR_ASSERT ((bytes & ANSI_BAE_MASK) == ANSI_BAE_MASK);
#endif
}
#else
#define STR_ASSERT_UPPED(x)  { }
#endif


// Subtract sizeof(StrImplementation::SBlock) from dptr and return address properly casted
inline StrImplementation::SBlock* Str::block() const
{
	BYTE* dptr_bytes = (BYTE*) dptr;
	return (StrImplementation::SBlock*) (dptr_bytes - sizeof(StrImplementation::SBlock));
}

// Add sizeof(StrImplementation::SBlock) to given pointer and cast as LPTSTR (static method)
inline STRCHAR* Str::block_s (StrImplementation::SBlock* b)
{
	BYTE* b_bytes = (BYTE*) b;
	return (STRCHAR*) (b_bytes + sizeof(StrImplementation::SBlock));
}

// Return StrImplementation::SBlock() and string addresses for null block

#define null_b()		\
	(&StrImplementation::oInitializer.m_Empty.data)
#define null_str()		\
	((STRCHAR*) (&StrImplementation::oInitializer.m_Empty.emptys))

#if defined(STR_UNICODE)  &&  !defined(STR_NO_RTLIBRARY)
namespace StrImplementation {
#ifdef STR_LINUX
#define STR_WCIS32BITS
typedef unsigned int UtfChType;
#else
#define STR_WCIS16BITS
typedef unsigned short UtfChType;
#endif
}
#endif


/*********************************************************************
* Additional simple debugging tools
*********************************************************************/

// Additional debugging tools
namespace StrImplementation {

inline void STR_abort()
{
#if defined(STR_NO_RTLIBRARY)  ||  defined(STR_WINDOWS_CE)
	TerminateProcess(GetCurrentProcess(), 255);
#else
	abort();
#endif
}

void Str_ValidatePtr(const void* x, BOOL ornull);
}
#if !defined(STR_ASSERT_PTR) && !defined(STR_ASSERT_PTR_OR_NULL)
#ifdef STR_DEBUG
#define STR_ASSERT_PTR(x)			Str_ValidatePtr(x, FALSE)
#define STR_ASSERT_PTR_OR_NULL(x)	Str_ValidatePtr(x, TRUE)
#else
#define STR_ASSERT_PTR(x)
#define STR_ASSERT_PTR_OR_NULL(x)
#endif
#endif


/*********************************************************************
* UNICODE and ANSI compatible versions of various CRT functions
*********************************************************************/

namespace StrImplementation {

inline void SBlock::SetFlag(DWORD flagmask)
{
	STR_ASSERT(Str::block_s(this) != null_str());
	m_Flags |= flagmask;
}

inline void SBlock::ClearFlag(DWORD flagmask)
{
	STR_ASSERT(Str::block_s(this) != null_str());
	m_Flags &= ~flagmask;
}


void STR_EXPORT STR_sprintf_i(STRCHAR* dest, size_t destSZC, const STRCHAR* fmt, int value);
void STR_EXPORT STR_sprintf_p(STRCHAR* dest, size_t destSZC, const STRCHAR* fmt, void* value);
#ifndef STR_NO_RTLIBRARY
void STR_EXPORT STR_sprintf_d(STRCHAR* dest, size_t destSZC, const STRCHAR* fmt, double value);
#endif

#ifdef STR_WIN32
inline void STR_ToUpper (STRCHAR* buffer, DWORD blen)
{ CharUpperBuff(buffer, blen); }
inline void STR_ToLower (STRCHAR* buffer, DWORD blen)
{ CharLowerBuff(buffer, blen); }
#else
void STR_EXPORT STR_ToUpper (STRCHAR* buffer, DWORD blen);
void STR_EXPORT STR_ToLower (STRCHAR* buffer, DWORD blen);
#endif


int STR_EXPORT STR_strcmp_ex(const STRCHAR *s, const STRCHAR *s2);
int STR_EXPORT STR_stricmp_ex(const STRCHAR *s, const STRCHAR *s2);

//
// Avoid using C RTL as much as possible
#ifdef STR_NO_RTLIBRARY

inline int  STR_strcmp(const STRCHAR *s, const STRCHAR *s2)		{ return lstrcmp (s, s2); }
inline int  STR_stricmp(const STRCHAR *s, const STRCHAR *s2)	{ return lstrcmpi (s, s2); }
int STR_EXPORT STR_strncmp(const STRCHAR *s, const STRCHAR *s2, int nchars);
int STR_EXPORT STR_strnicmp(const STRCHAR *s, const STRCHAR *s2, int nchars);

long STR_EXPORT STR_strtol(const STRCHAR *s, STRCHAR **endptr, int ibase);
// strtod not available

// strlen
inline CPOS STR_astrlen(const STRACHAR *s)		{ return (CPOS) lstrlenA(s); }
inline CPOS STR_wstrlen(const STRWCHAR *s)		{ return (CPOS) lstrlenW(s); }
inline CPOS STR_strlen(const STRCHAR *s)		{ return (CPOS) lstrlen(s); }
// memcpy and memmove
void STR_EXPORT STR_copy (LPVOID dest, LPCVOID src, int bytes);
void STR_EXPORT STR_move (LPVOID dest, LPCVOID src, int bytes);

// Uppercase / lowercase char helpers
inline int STR_CharIsUpper (STRCHAR value)
{ return (value >= 'A' && value <= 'Z'); }
inline STRCHAR STR_CharToUpper (STRCHAR value)
{ return (value >= 'a' && value <= 'z') ? (STRCHAR) (value - ('a' - 'A')) : value; }
inline int STR_CharIsLower (STRCHAR value)
{ return (value >= 'a' && value <= 'z'); }
inline STRCHAR STR_CharToLower (STRCHAR value)
{ return (value >= 'A' && value <= 'Z') ? (STRCHAR) (value + ('a' - 'A')) : value; }
inline BOOL STR_CharIsAlpha (STRCHAR value)
{ return STR_CharIsUpper(value)  ||  STR_CharIsLower(value); }
inline BOOL STR_CharIsSpace (STRCHAR value)
{ return (value >= 0x09 && value <= 0x0D)  ||  (value == 0x20); }

#else  // ifdef STR_NO_RTLIBRARY

#ifdef STR_UNICODE

#ifdef STR_WIN32
inline int STR_strcmp(const STRCHAR* s, const STRCHAR* s2)				{ return wcscmp (s, s2); }
inline int STR_strncmp(const STRCHAR *s, const STRCHAR *s2, int nchars)	{ return wcsncmp (s, s2, nchars); }
#else
int STR_strcmp(const STRCHAR* s, const STRCHAR* s2);				// Non-inline
int STR_strncmp(const STRCHAR *s, const STRCHAR *s2, int nchars);	// Non-inline
#endif

#else

inline int  STR_strcmp(const STRCHAR *s, const STRCHAR *s2)		{ return strcmp (s, s2); }
inline int  STR_strncmp(const STRCHAR *s, const STRCHAR *s2, int nchars)	{ return strncmp (s, s2, nchars); }

#endif

#ifdef STR_WIN32
#ifdef STR_UNICODE
inline int  STR_stricmp(const STRCHAR *s, const STRCHAR *s2)	{ return _wcsicmp (s, s2); }
inline int  STR_strnicmp(const STRCHAR *s, const STRCHAR *s2, int nchars)	{ return _wcsnicmp(s, s2, nchars); }
#else
#ifndef STR_BORCPPBUILDER
inline int  STR_stricmp(const STRCHAR *s, const STRCHAR *s2)	{ return _stricmp (s, s2); }
inline int  STR_strnicmp(const STRCHAR *s, const STRCHAR *s2, int nchars)	{ return _strnicmp(s, s2, nchars); }
#else
inline int  STR_stricmp(const STRCHAR *s, const STRCHAR *s2)	{ return stricmp (s, s2); }
inline int  STR_strnicmp(const STRCHAR *s, const STRCHAR *s2, int nchars)	{ return strnicmp(s, s2, nchars); }
#endif
#endif
#else
int STR_stricmp(const STRCHAR *s, const STRCHAR *s2);
int STR_strnicmp(const STRCHAR *s, const STRCHAR *s2, int nchars);
#endif

#if defined(STR_WIN32)

inline long STR_strtol(const STRCHAR *s, STRCHAR **endptr, int ibase)
#ifdef STR_UNICODE
{ return wcstol(s, endptr, ibase); }
#else
{ return strtol(s, endptr, ibase); }
#endif

inline double STR_strtod(const STRCHAR *s, STRCHAR **endptr)
#ifdef STR_UNICODE
{ return wcstod(s, endptr); }
#else
{ return strtod(s, endptr); }
#endif

#else

long   STR_strtol(const STRCHAR *s, STRCHAR **endptr, int ibase);	// Linux: in cpp
double STR_strtod(const STRCHAR *s, STRCHAR **endptr);				// Linux: in cpp

#endif


// strlen
inline CPOS STR_astrlen(const STRACHAR *s)		{ return (CPOS) strlen(s); }
#ifndef STR_NO_UNICODE
inline CPOS STR_wstrlen(const STRWCHAR *s)		{ return (CPOS) wcslen(s); }
#endif
#ifdef STR_UNICODE
#define STR_strlen STR_wstrlen
#else
#define STR_strlen STR_astrlen
#endif


// memcpy, memmove
inline void STR_copy (LPVOID dest, LPCVOID src, int bytes)
{
#if defined(STR_WIN32)  ||  defined(STR_GCC_295)
	memmove (dest, src, bytes);
#else
	__builtin_memmove (dest, src, bytes);
#endif
}

inline void STR_move (LPVOID dest, LPCVOID src, int bytes)
{ memmove (dest, src, bytes); }


// Uppercase / lowercase low level helpers for Linux / Unicode
#if defined(STR_LINUX)  &&  defined(STR_UNICODE)

#define iswalpha(x)		\
	( ((x) >= 'a'  &&  (x) <= 'z')  ||  ((x) >= 'A'  &&  (x) <= 'Z') )
#define iswalnum(x)		\
	( iswalpha(x)  ||  (x >= '0'  &&  x <= '9') )
#define iswupper(x)		\
	( (x) >= 'A'  &&  (x) <= 'Z' )
#define iswlower(x)		\
	( (x) >= 'a'  &&  (x) <= 'z' )
#define towupper(x)		\
	( iswlower(x) ? ((int) ((x)-32)) : (x) )
#define towlower(x)		\
	( iswupper(x) ? ((int) ((x)+32)) : (x) )
#define iswspace(x)		\
	( (x) == 32  ||  ( (x) >= 0x09  &&  (x) <= 0x0D ) )

#endif


// Uppercase / lowercase char helpers

inline int STR_CharIsUpper (STRCHAR value)
{
#ifdef STR_ANSI
	return isupper((int) (unsigned char) value);
#else
	return iswalpha(value) && iswupper(value);
#endif
}

inline int STR_CharIsLower (STRCHAR value)
{
#ifdef STR_ANSI
	return islower((int) (unsigned char) value);
#else
	return iswalpha(value) && iswlower(value);
#endif
}

inline BOOL STR_CharIsAlpha (STRCHAR value)
{
#ifdef STR_ANSI
	return isalpha((int) (unsigned char) value) ? TRUE : FALSE;
#else
	return iswalpha(value) ? TRUE : FALSE;
#endif
}

inline STRCHAR STR_CharToUpper (STRCHAR value)
{
#ifdef STR_ANSI
	return STR_CharIsLower(value) ? (STRCHAR) toupper((int) (unsigned char) value) : value;
#else
	return towupper(value);
#endif
}

inline STRCHAR STR_CharToLower (STRCHAR value)
{
#ifdef STR_ANSI
	return STR_CharIsUpper(value) ? (STRCHAR) tolower((int) (unsigned char) value) : value;
#else
	return towlower(value);
#endif
}

inline BOOL STR_CharIsSpace (STRCHAR value)
{
#ifdef STR_ANSI
	return isspace((int) (unsigned char) value) ? TRUE : FALSE;
#else
	return iswspace(value) ? TRUE : FALSE;
#endif
}

#endif

inline BOOL STR_CharIsDigit(STRCHAR value)
{ return (value >= '0'  &&  value <= '9') ? TRUE : FALSE; }

// memcmp: no alternative Windows function, so we'll keep it this
//   way even when STR_NO_RTLIBRARY is defined
inline int STR_memcmp(LPCVOID b1, LPCVOID b2, int bytes)
	{ return memcmp(b1, b2, bytes); }

}


/*********************************************************************
* Proc:		STR_ToAnsi, STR_ToAnsiBytes
* Purpose:	Helpers
*********************************************************************/

namespace StrImplementation {

#ifndef STR_NO_UNICODE

#ifdef STR_LINUX
// Return # of characters in the passed wide string
inline int STR_wcslen(const STRWCHAR* source)
{
	int count=0;
	while (source[0] != 0) {
		count++;
		source++;
	}
	return count;
}
#endif

inline int STR_ToAnsiBytes(const STRWCHAR *source)
{
#ifdef STR_WIN32
	return WideCharToMultiByte(CP_ACP, 0, source, -1, NULL, 0, NULL, NULL);
#else
	// wcstombs doesn't work well with deststring=NULL in some environments,
	//   so we'll do our calculations assuming SBCS
	return STR_wcslen(source) + sizeof(STRACHAR);	// Account for terminating null
#endif
}

inline int STR_ToAnsiBytes_Counted(const STRWCHAR *source, int source_len)
{
	STR_ASSERT(source_len > 0);
#ifdef STR_WIN32
	return WideCharToMultiByte(CP_ACP, 0, source, source_len, NULL, 0, NULL, NULL);
#else
	// Always exactly equal to the # of chars in the wide string
	return source_len;		
#endif
}

inline void STR_ToAnsi(const STRWCHAR *source, STRACHAR *target, int targbytes)
{
#ifdef STR_WIN32
	WideCharToMultiByte(CP_ACP, 0, source, -1, target, targbytes-1, NULL, NULL);
	target[targbytes-1] = 0;
#else
	wcstombs(target, source, targbytes-1);
	target[targbytes-1] = 0;
#endif
}

inline void STR_ToAnsi_Counted(const STRWCHAR *source, STRACHAR *target, int targbytes, int source_len)
{
	STR_ASSERT(source_len > 0);
#ifdef STR_WIN32
	WideCharToMultiByte(CP_ACP, 0, source, source_len, target, targbytes, NULL, NULL);
#else
	wcstombs(target, source, source_len);
#endif
}

inline int STR_ToWideSlots(const STRACHAR *source)
{
#ifdef STR_WIN32
	return MultiByteToWideChar(CP_ACP, 0, source, -1, NULL, 0);
#else
	// mbstowcs doesn't work well with deststring=NULL in some environments,
	//   so we'll do our calculations assuming SBCS
	return strlen(source) + sizeof(STRACHAR);	// Account for terminating null
#endif
}

inline int STR_ToWideSlots_Counted(const STRACHAR *source, int source_len)
{
	STR_ASSERT(source_len > 0);
#ifdef STR_WIN32
	return MultiByteToWideChar(CP_ACP, 0, source, source_len, NULL, 0);
#else
	// Always exactly equal to the # of chars in the sbcs string
	return source_len;
#endif
}

inline void STR_ToWide(const STRACHAR *source, STRWCHAR *target, int targchars)
{
#ifdef STR_WIN32
	MultiByteToWideChar(CP_ACP, 0, source, -1, target, targchars-1);
	target[targchars-1] = 0;
#else
	size_t fresult = mbstowcs(target, source, targchars-1);
	if (fresult == (size_t) (targchars-1))
		target[targchars-1] = 0;
#endif
}

inline void STR_ToWide_Counted(const STRACHAR *source, STRWCHAR *target, int targchars, int source_len)
{
	STR_ASSERT(source_len > 0);
#ifdef STR_WIN32
	MultiByteToWideChar(CP_ACP, 0, source, source_len, target, targchars);
#else
	mbstowcs(target, source, source_len);
#endif
}

#endif

}


/*********************************************************************
* Proc:		STR_GetCurrentThreadID
*********************************************************************/

namespace StrImplementation {

#ifdef STR_THREAD_SAFE
inline DWORD STR_GetCurrentThreadID()
{
	return GetCurrentThreadId();
}
#endif

#ifndef STR_THREAD_SAFE
inline STRCHAR* NewEmptyMT(CPOS prealloc)
{
	if (prealloc == 0)
		return null_str();
	else
		return NewEmptyMTForceNI(prealloc);
}
#endif

}


#if !defined(STR_THREAD_SAFE)  ||  defined(STR_NODEBUG)
inline StrImplementation::verifymt::verifymt(Str* obj)
{
	the_object = obj;
}
#endif

// This implementation _never_ checks for correctness
inline StrImplementation::verifymt::verifymt(Str* obj, int /*dummy*/)
{
	the_object = obj;
}


/*********************************************************************
* Proc:		Str::Str (STRCHARCHAR char)
* Purpose:	Instantiates a string consisting of a single character
*********************************************************************/

inline Str::Str(STRCHAR source)
{
	// Handle null character case
	if (source == 0) {
		dptr = null_str();
		return;
	}
	// Construct 1-character string
	STRCHAR source_str[2];
	source_str[0] = source;
	source_str[1] = 0;
	NewFromString(source_str, 1, 0, FALSE);
}


/*********************************************************************
* Proc:		Str::IsMT(), SetMT()
*********************************************************************/

inline BOOL Str::IsMT() const
{
#ifdef STR_THREAD_SAFE
	return block()->GetFlag(SFLAG_MT);
#else
	return FALSE;
#endif
}

inline void Str::SetMT_Internal()
{
#ifdef STR_THREAD_SAFE
	//STR_ASSERT(block()->m_OwnerT == StrImplementation::STR_GetCurrentThreadID()); Bugfix Aug 09, 2002
	STR_ASSERT(block()->m_OwnerT == 0  ||  block()->m_OwnerT == StrImplementation::STR_GetCurrentThreadID());
	if (dptr == null_str())
		// Must move to separate instance with proper thread marker
		dptr = StrImplementation::NewEmptyMT(block()->m_Alloc);
	else
		block()->SetFlag(SFLAG_MT);
#endif
}

#ifndef STR_THREAD_SAFE

inline void Str::SetMT()
{ }

inline void Str::ChangeOwnerThread()
{ }

#endif

inline Str::Str()
{
	dptr = null_str();
}

inline Str::Str(const StrMTType&)
{
#ifdef STR_THREAD_SAFE
	dptr = StrImplementation::NewEmptyMT(0);
#else
	dptr = null_str();
#endif
}


/*********************************************************************
* Proc:		Str::Str(const char*)
* Purpose:	Construct an instance that exactly copies the specified
*			string.  An optional second parameter specifies the buffer
*			size (will be ignored if it is less than what's needed)
*********************************************************************/

#ifdef STR_UNICODE
#define IMP_COPY_MT_FOREIGNCHARSET(wantmt)			\
	NewFromString(L"", 0, prealloc, wantmt);		\
	*this = s
#else
#define IMP_COPY_MT_FOREIGNCHARSET(wantmt)			\
	NewFromString("", 0, prealloc, wantmt);			\
	*this = s
#endif

#define IMP_COPY_MT_NATIVECHARSET(wantmt)				\
	NewFromString(s, prealloc, wantmt)

inline Str::Str(const STRACHAR *s)
{
	CPOS prealloc = 0;
#ifdef STR_UNICODE
	IMP_COPY_MT_FOREIGNCHARSET(FALSE);
#else
	IMP_COPY_MT_NATIVECHARSET(FALSE);
#endif
}

inline Str::Str(CPOS prealloc, const STRACHAR *s)
{
#ifdef STR_UNICODE
	IMP_COPY_MT_FOREIGNCHARSET(FALSE);
#else
	IMP_COPY_MT_NATIVECHARSET(FALSE);
#endif
}

inline Str::Str(const STRACHAR* s, const StrMTType&)
{
	CPOS prealloc = 0;
#ifdef STR_UNICODE
	IMP_COPY_MT_FOREIGNCHARSET(TRUE);
#else
	IMP_COPY_MT_NATIVECHARSET(TRUE);
#endif
}

inline Str::Str(CPOS prealloc, const STRACHAR* s, const StrMTType&)
{
#ifdef STR_UNICODE
	IMP_COPY_MT_FOREIGNCHARSET(TRUE);
#else
	IMP_COPY_MT_NATIVECHARSET(TRUE);
#endif
}

#ifndef STR_NO_UNICODE

inline Str::Str(const STRWCHAR *s)
{
	CPOS prealloc = 0;
#ifdef STR_UNICODE
	IMP_COPY_MT_NATIVECHARSET(FALSE);
#else
	IMP_COPY_MT_FOREIGNCHARSET(FALSE);
#endif
}

inline Str::Str(CPOS prealloc, const STRWCHAR *s)
{
#ifdef STR_UNICODE
	IMP_COPY_MT_NATIVECHARSET(FALSE);
#else
	IMP_COPY_MT_FOREIGNCHARSET(FALSE);
#endif
}

inline Str::Str(const STRWCHAR *s, const StrMTType&)
{
	CPOS prealloc = 0;
#ifdef STR_UNICODE
	IMP_COPY_MT_NATIVECHARSET(TRUE);
#else
	IMP_COPY_MT_FOREIGNCHARSET(TRUE);
#endif
}

inline Str::Str(CPOS prealloc, const STRWCHAR *s, const StrMTType&)
{
#ifdef STR_UNICODE
	IMP_COPY_MT_NATIVECHARSET(TRUE);
#else
	IMP_COPY_MT_FOREIGNCHARSET(TRUE);
#endif
}

#endif

#undef IMP_COPY_MT_FOREIGNCHARSET
#undef IMP_COPY_MT_NATIVECHARSET


/*********************************************************************
* Proc:		Str::operator = Str
* Purpose:	Copies a Str into another Str, destroying the
*			previous content.
*********************************************************************/

inline const Str& Str::operator=(const Str& source)
{
	BOOL wasmt = IsMT();
	dptr = StrImplementation::STR_CsCopyCore(this, source);
	if (wasmt)
		SetMT();
	return *this;
}


/*********************************************************************
* Proc:		Str::operator = LPCSTR and LPCWSTR
* Purpose:	Wrappers around internal methods
*********************************************************************/

inline const Str& Str::operator=(const STRACHAR *s)
{
#ifdef STR_UNICODE
	Copy_FromAnsi(s);
#else
	Copy_Native(s);
#endif
	return *this;
}

#ifndef STR_NO_UNICODE

inline const Str& Str::operator=(const STRWCHAR* s)
{
#ifdef STR_UNICODE
	Copy_Native(s);
#else
	Copy_FromUni(s);
#endif
	return *this;
}

#endif


/*********************************************************************
* Proc:		Str::destructor
*********************************************************************/

inline Str::~Str()
{
	StrImplementation::verifymt mtcheck(this);
	StrImplementation::STR_rls_block(mtcheck);
}


/*********************************************************************
* Proc:		Str::GetFirstChar, IsEmpty, GetLength
*********************************************************************/

inline STRCHAR Str::GetFirstChar() const
{
	return dptr[0];
}

inline BOOL Str::IsEmpty() const
{
	return dptr[0] == 0;
}

inline CPOS Str::GetLength() const
{
	return block()->m_Length;
}

inline CPOS Str::GetBufferLength() const
{
	return block()->m_Alloc;
}

inline Str Str::Mid(CPOS start) const
{
	return Mid(start, GetLength());
}

/*********************************************************************
* Proc:		Str::operators LPCTSTR, [], methods GetString and GetAt
*********************************************************************/

inline STRCHAR Str::operator[](const CPOS idx) const
{
#ifdef STR_DEBUG
	if (idx < 0  ||  idx >= GetLength())
		return (STRCHAR) 0;
#endif
	return dptr[idx];
}

inline STRCHAR Str::GetAt(CPOS idx) const
{
	return operator[] (idx);
}

inline Str::operator const STRCHAR*() const
{ return dptr; }

inline const STRCHAR* Str::GetString() const
{ return dptr; }

#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)

inline CString Str::GetCString() const
{ return CString(GetString()); }

inline Str::operator CString() const
{ return CString(GetString()); }

#endif

/*********************************************************************
* Proc:		Str::Compare and CompareNoCase
*********************************************************************/

inline int Str::Compare (const STRCHAR* match) const
{ return StrImplementation::STR_strcmp_ex(dptr, match); }

inline int Str::CompareNoCase (const STRCHAR* match) const
{ return StrImplementation::STR_stricmp_ex(dptr, match); }


/*********************************************************************
* Proc:		Str::[operator ==] and [operator !=] inlined forms
*********************************************************************/

inline BOOL operator ==(const STRCHAR* s1, const Str& s2)
{ return (s2 == s1); }

inline BOOL operator !=(const Str& s1, const Str& s2)
{ return !(s1 == s2); }
inline BOOL operator !=(const Str& s1, const STRCHAR* s2)
{ return !(s1 == s2); }
inline BOOL operator !=(const STRCHAR* s1, const Str& s2)
{ return !(s2 == s1); }

inline BOOL operator< (const Str& s1, const Str& s2)
{ return (s1.Compare(s2) < 0); }
inline BOOL operator< (const Str& s1, const STRCHAR* s2)
{ return (s1.Compare(s2) < 0); }
inline BOOL operator< (const STRCHAR* s1, const Str& s2)
{ return (s2.Compare(s1) >= 0); }

inline BOOL operator<= (const Str& s1, const Str& s2)
{ return (s1.Compare(s2) <= 0); }
inline BOOL operator<= (const Str& s1, const STRCHAR* s2)
{ return (s1.Compare(s2) <- 0); }
inline BOOL operator<= (const STRCHAR* s1, const Str& s2)
{ return (s2.Compare(s1) > 0); }

inline BOOL operator> (const Str& s1, const Str& s2)
{ return (s1.Compare(s2) > 0); }
inline BOOL operator> (const Str& s1, const STRCHAR* s2)
{ return (s1.Compare(s2) > 0); }
inline BOOL operator> (const STRCHAR* s1, const Str& s2)
{ return (s2.Compare(s1) <= 0); }

inline BOOL operator>= (const Str& s1, const Str& s2)
{ return (s1.Compare(s2) >= 0); }
inline BOOL operator>= (const Str& s1, const STRCHAR* s2)
{ return (s1.Compare(s2) >= 0); }
inline BOOL operator>= (const STRCHAR* s1, const Str& s2)
{ return (s2.Compare(s1) < 0); }


/*********************************************************************
* Proc:		Str::operator += (Str)
* Purpose:	Append a string to what we already contain.
*********************************************************************/

inline void Str::operator += (const Str& obj)
{
	StrImplementation::verifymt mtcheck(this);
	CoreAppendChars (obj.dptr, obj.block()->m_Length, mtcheck);
}


/*********************************************************************
* Proc:		Str::AppendString - synonyms for operators +=
*********************************************************************/

inline Str& Str::AppendString(const Str& obj)
{ *this += obj;  return *this; }

inline Str& Str::AppendString(const STRACHAR *s)
{ *this += s;  return *this; }

#ifndef STR_NO_UNICODE
inline Str& Str::AppendString(const STRWCHAR* s)
{ *this += s;  return *this; }
#endif


inline void Str::operator += (STRACHAR ch)
{ AppendChar (ch); }

#ifndef STR_NO_UNICODE
inline void Str::operator += (STRWCHAR ch)
{ AppendChar (ch); }
#endif


/*********************************************************************
* Inline methods that compile to a single helper function call
*********************************************************************/

inline int Str::Replace (const STRCHAR* oldstr, const STRCHAR* newstr)
{ return StrImplementation::STR_ReplaceCore(*this, oldstr, newstr, FALSE); }

inline int Str::ReplaceNoCase (const STRCHAR* oldstr, const STRCHAR* newstr)
{ return StrImplementation::STR_ReplaceCore(*this, oldstr, newstr, TRUE); }

inline Str& Str::MakeUpper()
{ Str_MakeUL(TRUE); return *this; }

inline Str& Str::MakeLower()
{ Str_MakeUL(FALSE); return *this; }

inline Str Str::Left (CPOS chars) const
{ return Mid(0, chars); }

inline Str& Str::InsertAtLeft (const STRCHAR* s, CPOS howmany /*= -1*/)
{ Insert(s, 0, howmany);  return *this; }


/*********************************************************************
* Char constructors and assignment operators
*********************************************************************/

inline Char::Char()
{ data = 0; }

#ifdef STR_ANSI

inline Char::Char(STRACHAR source)
{ data = source; }

#ifndef STR_NO_UNICODE
inline Char::Char(STRWCHAR source)
{ data = ToANSI(source); }
#endif

inline const Char& Char::operator=(STRACHAR source)
{ data = source; return *this; }

inline void Char::Set(STRACHAR source)
{ data = source; }

#ifndef STR_NO_UNICODE

inline const Char& Char::operator=(STRWCHAR source)
{ data = ToANSI(source); return *this; }

inline void Char::Set(STRWCHAR source)
{ data = ToANSI(source); }

#endif

#else // ifdef STR_ANSI

inline Char::Char(STRWCHAR source)
{ data = source; }

inline Char::Char(STRACHAR source)
{ data = ToUnicode(source); }

inline const Char& Char::operator=(STRACHAR source)
{ data = ToUnicode(source); return *this; }

inline void Char::Set(STRACHAR source)
{ data = ToUnicode(source); }

inline const Char& Char::operator=(STRWCHAR source)
{ data = source; return *this; }

inline void Char::Set(STRWCHAR source)
{ data = source; }

#endif


/*********************************************************************
* Char get / compare operations
*********************************************************************/

inline BOOL Char::IsNull() const
{ return data == 0; }

inline Char::operator const STRCHAR() const
{ return data; }

inline const STRCHAR Char::GetChar() const
{ return data; }

#ifdef STR_ANSI

inline const STRACHAR Char::GetCharA() const
{ return data; }

#ifndef STR_NO_UNICODE
inline const STRWCHAR Char::GetCharW() const
{ return ToUnicode(data); }
#endif

#else

inline const STRACHAR Char::GetCharA() const
{ return ToANSI(data); }

inline const STRWCHAR Char::GetCharW() const
{ return data; }

#endif

inline BOOL operator ==(const Char& c1, const Char& c2)
{ return c1.GetChar() == c2.GetChar(); }

inline BOOL operator ==(STRCHAR c1, const Char& c2)
{ return c1 == c2.GetChar(); }

inline BOOL operator ==(const Char& c1, STRCHAR c2)
{ return c1.GetChar() == c2; }

inline BOOL operator !=(const Char& c1, const Char& c2)
{ return c1.GetChar() != c2.GetChar(); }

inline BOOL operator !=(STRCHAR c1, const Char& c2)
{ return c1 != c2.GetChar(); }

inline BOOL operator !=(const Char& c1, STRCHAR c2)
{ return c1.GetChar() != c2; }


/*********************************************************************
* Uppercase / lowercase / IsXXX family
*********************************************************************/

inline Char& Char::MakeUpper()
{ data = StrImplementation::STR_CharToUpper(data); return *this; }

inline Char& Char::MakeLower()
{ data = StrImplementation::STR_CharToLower(data); return *this; }

inline BOOL Char::IsUpper() const
{ return StrImplementation::STR_CharIsUpper(data) ? TRUE : FALSE; }

inline BOOL Char::IsLower() const
{ return StrImplementation::STR_CharIsLower(data) ? TRUE : FALSE; }

inline BOOL Char::IsAlpha() const
{ return StrImplementation::STR_CharIsAlpha(data); }

inline BOOL Char::IsDigit() const
{ return StrImplementation::STR_CharIsDigit(data); }

inline BOOL Char::IsAlphaNumeric() const
{ return IsAlpha()  ||  IsDigit(); }


/*********************************************************************
* Buffer management and miscellaneous
*********************************************************************/

inline STRCHAR* Str::GetBuffer(CPOS min_chars)
{
	Preallocate(min_chars);
	return GetBuffer();
}

#ifdef STR_MCPP_FEATURES
inline System::String __gc* Str::ToManaged() const
{
	System::String __gc* retS;
#ifdef STR_UNICODE
	retS = new System::String(*this);
#else
	STRWCHAR* wide = GetStringW();
	retS = new System::String(wide);
	OS_free(wide);
#endif
	return retS;
}
#endif


/*********************************************************************
* Serialization
*********************************************************************/

#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)

inline CArchive& operator<<(CArchive& ar, const Str& source)
{
	Str::WriteStringLP(ar, (const STRCHAR*) source);
	return ar;
}

inline CArchive& operator>>(CArchive& ar, Str& dest)
{
	Str::ReadStringLP(ar, dest);
	return ar;
}

#endif


/*********************************************************************
* Extra COM/Automation-related features
*********************************************************************/

#if defined(STR_WIN32) && !defined(STR_NO_WINSTUFF) && !defined(STR_NO_UNICODE)

inline VARIANT Str::ToVariant() const
{
	VARIANT ret;
	ret.vt = VT_BSTR;
	ret.bstrVal = AllocSysString();
	return ret;
}

#endif


/*********************************************************************
* Regular expressions
*********************************************************************/

#ifdef STR_USE_REGEX

inline void RegexpMatch::GetCompleteMatch(const STRCHAR **szStart, const STRCHAR **szEnd)
{
	STR_ASSERT(szStart != NULL);
	STR_ASSERT(szEnd != NULL);
	*szStart = m_CompleteMatch.szStart;
	*szEnd = m_CompleteMatch.szEnd;
}

inline Str RegexpMatch::GetCompleteMatch()
{
	Str dest;
	dest.AppendChars(m_CompleteMatch.szStart, (CPOS) (m_CompleteMatch.szEnd-m_CompleteMatch.szStart));
	return dest;
}

#endif


/*********************************************************************
* StrArray
*********************************************************************/

#ifdef STR_USE_EXTRAS

inline int StrArray::GetUpperBound() const
{ return GetCount() - 1; }

inline const Str& StrArray::operator[](int index) const
{ return GetAt(index); }

inline Str& StrArray::operator[](int index)
{ return ElementAt(index); }

inline int StrArray::GrowRel (int additionalCount)
{ int count = GetUpperBound(); MaybeGrow(count + additionalCount); return count+1; } 

inline void StrArray::SetAtGrow(int index, const STRCHAR* elem)
{ MaybeGrow(index); SetAt(index, elem); }

inline void StrArray::SetAtGrow(int index, const Str& elem)
{ MaybeGrow(index); SetAt(index, elem); }

inline int StrArray::Add(const STRCHAR* elem)
{ int newpos = GrowRel(1); SetAt(newpos, elem); return newpos; }

inline int StrArray::Add(const Str& elem)
{ int newpos = GrowRel(1); SetAt(newpos, elem); return newpos; }

inline void StrArray::InsertAt(int index, const STRCHAR* elem, int nCount /*= 1*/)
{ Str toIns(elem); InsertAt(index, toIns, nCount); }

#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)

inline void StrArray::SetAt(int index, const CString& elem)
{ SetAt (index, (const STRCHAR*) elem); }

inline void StrArray::SetAtGrow(int index, const CString& elem)
{ SetAtGrow (index, (const STRCHAR*) elem); }

inline int StrArray::Add(const CString& elem)
{ return Add ((const STRCHAR*) elem); }

inline void InsertAt(int index, const CString& newElement, int nCount /*= 1*/)
{ InsertAt (index, (const STRCHAR*) newElement, nCount); }

#endif // STR_USE_EXTRAS and MFC incorporated

#endif // STR_USE_EXTRAS


#if defined(STR_USE_REGEX)  ||  defined(STR_USE_EXTRAS)

#include <new>			// Placement new

namespace StrImplementation {

#define VECLITE_MINCAPACITY 4

template<class T> class veclite
{
private:
	T*		_v;
	size_t	_size;
	size_t	_capacity;
public:
	typedef veclite<T> _MyT;
	veclite();
	veclite(const _MyT& copyFrom);
	~veclite();
	_MyT& operator=(const _MyT& copyFrom);
public:
	void clear();
	T& at(size_t idx);
	T& operator[] (size_t idx);
	const T& at(size_t idx) const;
	const T& operator[] (size_t idx) const;
	T& back();
	const T& back() const;
	size_t size() const							{ return _size; }
	size_t capacity() const						{ return _capacity; }
	size_t max_size() const;
	bool empty() const							{ return (_size == 0); }
	void pop_back();
	void push_back(const T& val)				{ doResize(_size+1, val); }
	void resize(size_t newSize, T initTo);
	void resize(size_t newSize)					{ resize (newSize, T()); }
	void reserve(size_t newCapacity);
	size_t insert(size_t idx, size_t count, const T& val);
	size_t insert(size_t idx, const T& val)		{ return insert(idx, 1, val); }
	void   erase_range(size_t idx, size_t count);
	void   erase (size_t idx)					{ erase_range(idx, 1); }
	T* direct_access()							{ return _v; }
	const T* direct_access() const				{ return _v; }
private:
	void increaseCapacity(size_t newCapacity, bool beFlexible);
	void doResize(size_t newSize, const T& initTo);
	void copyFromCore(const _MyT& copyFrom);
	// Replace the implementations of the following with "throw
	// MyTypeOfException"; these default implementations are very, very crude
	void throwMaxLengthExceeded()				{ STR_abort(); }
	void throwNoMemory()						{ STR_abort(); }
};

template<class T>
StrImplementation::veclite<T>::veclite()
{
	_v = NULL;
	_size = 0;
	_capacity = 0;
}

template<class T>
StrImplementation::veclite<T>::~veclite()
{
	clear();
}

template<class T>
void veclite<T>::copyFromCore(const _MyT& copyFrom)
{
	size_t copyFromSize = copyFrom._size;
	if (copyFromSize == 0)
		return;
	increaseCapacity(copyFromSize, true);
	for (size_t i=0; i<copyFromSize; i++) {
		T* tempPtr = &(_v[i]);
		new (tempPtr) T(copyFrom._v[i]);
	}
	_size = copyFromSize;
}

template<class T>
veclite<T>& veclite<T>::operator=(const _MyT& copyFrom)
{
	clear();
	copyFromCore(copyFrom);
	return *this;
}

template<class T>
StrImplementation::veclite<T>::veclite(const _MyT& copyFrom)
{
    _v = NULL;
	_size = 0;
	_capacity = 0;
	copyFromCore(copyFrom);
}

template<class T>
void veclite<T>::clear()
{
	if (_v == NULL)
		return;
	for (size_t i=0; i<_size; i++)
		_v[i].~T();
	Str::OS_free(_v);
	_v = NULL;
	_size = 0;
	_capacity = 0;
}

template<class T>
T& veclite<T>::operator[] (size_t idx)
{
	STR_ASSERT(idx < _size);
	return _v[idx];
}

template<class T>
T& veclite<T>::at (size_t idx)
{ return operator[] (idx); }

template<class T>
const T& veclite<T>::operator[] (size_t idx) const
{
	STR_ASSERT(idx < _size);
	return _v[idx];
}

template<class T>
const T& veclite<T>::at (size_t idx) const
{ return operator[] (idx); }

template<class T>
T& veclite<T>::back()
{
	STR_ASSERT(_size > 0);
	return at(_size-1);
}

template<class T>
const T& veclite<T>::back() const
{
	STR_ASSERT(_size > 0);
	return at(_size-1);
}

template<class T>
size_t veclite<T>::max_size() const
{
	// (2^31 - 16) / sizeof(T)
	return (size_t) (2147483632 / sizeof(T));
}

template<class T>
void veclite<T>::pop_back()
{
	STR_ASSERT(_size > 0);
	_v[_size-1].~T();
	_size--;
}

template<class T>
void veclite<T>::increaseCapacity(size_t newCapacity, bool beFlexible)
{
	STR_ASSERT(newCapacity > _capacity);
	if (newCapacity > max_size())
		throwMaxLengthExceeded();
	//
	// In fixed mode we should always allocate exactly as much as requested
	if (!beFlexible)
		goto DoRealloc;
	//
	// Flexible mode code follows.  First, establish a minimum
	if (newCapacity < VECLITE_MINCAPACITY)
		newCapacity = VECLITE_MINCAPACITY;
	// If we've gone past 1/2 of our maximum size, ignore heuristics
	// and increase capacity to the exact value requested.  This may
	// cause a great deal of thrashing but is the only decent thing to
	// do (and callers shouldn't allocate a 1G array without explicitly
	// doing a capacity change first!)
	if (_capacity >= (size_t) (max_size() / 2))
		goto DoRealloc;		// Fall back to exact requested value
	// Compute current capacity + 50%, compare it to requested capacity,
	// and take whichever is greater
	{
		size_t newHeuristical = _capacity + (size_t) (_capacity / 2);
		if (newHeuristical > newCapacity)
			newCapacity = newHeuristical;
	}
DoRealloc:
	T* newV = (T*) Str::OS_malloc(newCapacity * sizeof(T));
	if (newV == NULL)
		throwNoMemory();
	if (_size > 0) {
		// Current size >0, so be sure to copy out existing elements
		for (size_t i=0; i<_size; i++) {
			T* tempPtrDest = &(newV[i]);
			new (tempPtrDest) T(_v[i]);
			_v[i].~T();
		}
		// And release old block
		Str::OS_free(_v);
	}
	_v = newV;
	_capacity = newCapacity;
}

template<class T>
void veclite<T>::resize(size_t newSize, T initTo)
{
	if (newSize != _size) {
		if (newSize == 0)
			clear();
		else
			doResize(newSize, initTo);
	}
}

//
// Must NOT be called with newSize==0
template<class T>
void veclite<T>::doResize(size_t newSize, const T& initTo)
{
	if (_size == newSize)
		return;
	if (_size < newSize) {
		if (newSize > _capacity)
			increaseCapacity(newSize, true);
		for (size_t i=_size; i<newSize; i++) {
			T* tempPtr = &(_v[i]);
			new (tempPtr) T(initTo);
		}
	}
	else {
		for (size_t i=newSize; i<_size; i++)
			_v[i].~T();
	}
	_size = newSize;
}

template<class T>
void veclite<T>::reserve(size_t newCapacity)
{
	STR_ASSERT(newCapacity > 0);
	if (_v == NULL)
		_capacity = 0;		// So our capacity increment logic doesn't go haywire
	if (newCapacity <= _capacity)
		return;
	increaseCapacity(newCapacity, true);
}


template<class T>
size_t veclite<T>::insert(size_t idx, size_t count, const T& val)
{
	if (count == 0)
		return idx;
	if (idx == size()) {
		// Special case: append
		doResize(size()+count, val);
		return idx;
	}
	STR_ASSERT(idx < size());
	size_t newCapacity = size()+count;
	if (newCapacity > _capacity)
		increaseCapacity(newCapacity, true);
	size_t i;
	for (i=size()-1; i>=idx  &&  i != (size_t) -1; i--) {
		T* tempPtr = &(_v[i+count]);
		new (tempPtr) T(_v[i]);
		_v[i].~T();
	}
	for (i=0; i<count; i++) {
		T* tempPtr = &(_v[idx+i]);
		new (tempPtr) T(val);
	}
	_size += count;
	return idx;
}

template<class T>
void veclite<T>::erase_range(size_t idx, size_t count)
{
	if (count == 0)
		return;
	STR_ASSERT(idx <= (size()-count));
	if (size() == count  &&  idx == 0) {
		clear();				// Ensure block is released completely
		return;
	}
	size_t i;
	for (i=size()-1; i>=(idx+count)  &&  i!=(size_t)-1; i--)
		_v[i-count] = _v[i];
	for (i=0; i<count; i++)
		_v[size()-(i+1)].~T();
	_size -= count;
}


} // namespace StrImplementation

#endif	// defined(STR_USE_REGEX)  ||  defined(STR_USE_EXTRAS)
