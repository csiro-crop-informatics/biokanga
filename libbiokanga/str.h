// Str.h
//
// Version 2.4.0
//
// This file is a part of the Str Library, and is subject to the
// clauses of the License Agreement accompanying the software.

#ifndef _STR_H
#define _STR_H

// User-defined header file: contains configuration declarations
#ifndef STR_DEFINES_DONE
#include "str_config.h"
#endif

// Low-level definitions and miscellaneous internal helpers
#include "str_support1.h"

// Additional parameter for constructors to automatically call SetMT
struct STR_EXPORT StrMTType;
extern STR_EXPORT const StrMTType& StrMT;


/*********************************************************************
* Class:	Char
* Purpose:	All-purpose single character handling
* Remarks:  Called StrChar under Borland C++ Builder
*********************************************************************/

#ifdef STR_BORCPPBUILDER
#define Char StrChar
#endif

class STR_EXPORT Char
{
private:
	STRCHAR	data;

// Construction, copying, assignment
public:
	Char();
	Char(STRACHAR source);
	void Set(STRACHAR source);
	const Char& operator = (STRACHAR source);
#ifndef STR_NO_UNICODE
	Char(STRWCHAR source);
	const Char& operator = (STRWCHAR source);
	void Set(STRWCHAR source);
#endif

// ANSI <-> Unicode
#ifndef STR_NO_UNICODE
	static STRACHAR ToANSI(STRWCHAR ch);
	static STRWCHAR ToUnicode(STRACHAR ch);
#endif

// Get attributes, get data, compare
	BOOL	 IsNull() const;
	operator const STRCHAR() const;
	const STRCHAR  GetChar() const;
	const STRACHAR GetCharA() const;
#ifndef STR_NO_UNICODE
	const STRWCHAR GetCharW() const;
#endif

// Miscellaneous operations
	Char& MakeUpper();
	Char& MakeLower();
	BOOL  IsUpper() const;
	BOOL  IsLower() const;
	BOOL  IsAlpha() const;
	BOOL  IsDigit() const;
	BOOL  IsAlphaNumeric() const;

};

BOOL operator ==(const Char& c1, const Char& c2);
BOOL operator ==(STRCHAR c1, const Char& c2);
BOOL operator ==(const Char& c1, STRCHAR c2);
BOOL operator !=(const Char& c1, const Char& c2);
BOOL operator !=(STRCHAR c1, const Char& c2);
BOOL operator !=(const Char& c1, STRCHAR c2);


/*********************************************************************
* Class:	Str
* Purpose:	All-purpose string handling
*********************************************************************/

class STR_EXPORT Str
{
// Data block and internal access methods
private:
	STRCHAR *dptr;

// Construction, copying, assignment
public:
	Str();
	Str(const StrMTType&);
	Str(const Str& source);
	Str(CPOS prealloc, const STRACHAR *s);
	Str(CPOS prealloc, const STRACHAR *s, const StrMTType&);
	Str(const STRACHAR *s);
	Str(const STRACHAR *s, const StrMTType&);
	Str(STRCHAR source);
	static void ConstructAt(Str* location, int count);
	static void ConstructAt(Str* location, int count, const StrMTType&);
	const Str& operator = (const Str& source);
	const Str& operator = (const STRACHAR *s);
#ifndef STR_NO_UNICODE
	Str(CPOS prealloc, const STRWCHAR *s);
	Str(CPOS prealloc, const STRWCHAR *s, const StrMTType&);
	Str(const STRWCHAR *s);
	Str(const STRWCHAR *s, const StrMTType&);
	const Str& operator = (const STRWCHAR *s);
#endif
#if defined(STR_WIN32) && !defined(STR_NO_WINSTUFF) && !defined(STR_NO_UNICODE)
	const Str& operator = (BSTR s);
	VARIANT ToVariant() const;
#endif
	~Str();
	static void DestructAt(Str* location, int count);

// MFC-specific features
#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)
	Str(const CString& source, CPOS prealloc = 0);
	const Str& operator = (const CString& source);
	CString GetCString() const;
	operator CString() const;
	void operator += (const CString& obj);
	static BOOL ReadString(CArchive& ar, Str& dest);
	static void WriteString(CArchive& ar, const STRCHAR* lpsz);
	static BOOL ReadStringLP(CArchive& ar, Str& dest);
	static void WriteStringLP(CArchive& ar, const STRCHAR* lpsz);
	friend CArchive& operator<<(CArchive& ar, const Str& source);
	friend CArchive& operator>>(CArchive& ar, Str& dest);
#endif

// Managed C++ features
#ifdef STR_MCPP_FEATURES
	Str(System::String __gc* source);
	const Str& operator = (System::String __gc* source);
	System::String __gc* ToManaged() const;
#endif

// Get attributes, get data, compare
	BOOL	 IsEmpty() const;
	CPOS	 GetLength() const;
	operator const STRCHAR*() const;
	STRCHAR GetFirstChar() const;
	STRCHAR GetLastChar() const;
	STRCHAR operator[] (const CPOS idx) const;
	STRCHAR GetAt(CPOS idx) const;		// Same as operator []
	Str   Left (CPOS chars) const;
	Str   Right (CPOS chars) const;
	Str   Mid (CPOS start, CPOS chars) const;
	Str   Mid (CPOS start) const;
	STRCHAR* GetBuffer();
	STRCHAR* GetBuffer(CPOS min_chars);
	STRCHAR* GetBufferSetLength (CPOS newlen);
	void     ReleaseBuffer(CPOS newlen = (CPOS) -1);

// Compare, find, replace
	int   Compare (const STRCHAR* match) const;		// -1, 0 or 1
	int   CompareNoCase (const STRCHAR* match) const;	// -1, 0 or 1
	CPOS  Find (STRCHAR ch, CPOS startat = 0) const;
	CPOS  FindNoCase (STRCHAR ch, CPOS startat = 0) const;
	CPOS  ReverseFind (STRCHAR ch, CPOS startat = (CPOS) -1) const;
	CPOS  ReverseFindNoCase (STRCHAR ch, CPOS startat = (CPOS) -1) const;
	CPOS  Find (const STRCHAR *substr, CPOS startat = 0) const;
	CPOS  FindNoCase (const STRCHAR *substr, CPOS startat = 0) const;
	CPOS  FindOneOf (const STRCHAR *charset, CPOS startat = 0) const;
	int   Replace(const STRCHAR *, const STRCHAR *newstr);
	int   ReplaceNoCase(const STRCHAR *oldstr, const STRCHAR *newstr);
	void  ReplaceAt(CPOS pos, CPOS numchars, const STRCHAR *newstr);

// Memory management and reset calls
	Str& Empty();			// Sets length to 0, but keeps buffer around
	Str& Reset();			// This also releases the buffer
	void Preallocate(CPOS size);
	CPOS GetBufferLength() const;
	void Compact(CPOS only_above_b = 0);
	static void CompactFreeMemory();
	static void SetShutdownFlag(BOOL shutdown = TRUE);

// Various
	void SetMT();
	BOOL IsMT() const;
	void ChangeOwnerThread();
	Str& EndInSlash();
	Str& Format(const STRACHAR *fmt, ...);
#ifndef STR_NO_UNICODE
	Str& Format(const STRWCHAR *fmt, ...);
#endif
#ifndef STR_NO_RTLIBRARY
	Str& FormatVA(const STRACHAR *fmt, va_list marker);
#ifndef STR_NO_UNICODE
	Str& FormatVA(const STRWCHAR *fmt, va_list marker);
#endif
#endif
	Str& AppendFormat(const STRCHAR *fmt, ...);
#if defined(STR_WIN32)  && !defined(STR_NO_WINSTUFF)
	Str& FormatRes(UINT resid, ...);
	BOOL LoadString(UINT resid);
#endif
#if defined(STR_WIN32)  &&  !defined(STR_NO_UNICODE) && !defined(STR_NO_WINSTUFF)
	BSTR AllocSysString() const;
	BSTR SetSysString(BSTR* str) const;
#endif
	Str& MakeUpper();
	Str& MakeLower();

// Catenation, truncation
	void operator += (const Str& obj);
	void operator += (const STRACHAR* s);
	void operator += (STRACHAR ch);		// Same as AppendChar
#ifndef STR_NO_UNICODE
	void operator += (const STRWCHAR *s);
	void operator += (STRWCHAR ch);		// Same as AppendChar
	Str& AppendString(const STRWCHAR *s);			// Same as +=
#endif
	void SetAt(CPOS idx, STRCHAR newch);		// Can't be beyond current end of string
	Str& AppendString(const Str& obj);				// Same as +=
	Str& AppendString(const STRACHAR* s);			// Same as +=
	Str& AppendChar(STRACHAR ch);
#ifndef STR_NO_UNICODE
	Str& AppendChar(STRWCHAR ch);
#endif
	void AppendChars(const STRCHAR *s, CPOS howmany = (CPOS) -1);
	Str& Insert(const STRCHAR *s, CPOS startat, CPOS howmany = (CPOS) -1);
	Str& InsertAtLeft(const STRCHAR *s, CPOS howmany = (CPOS) -1);
	Str& AppendInt(int value);
	int  ToInt(BOOL* error = NULL) const;
#ifndef STR_NO_RTLIBRARY
	Str& AppendFloat(float value, UINT after_dot);
	Str& AppendDouble(double value, UINT after_dot);
	float  ToFloat(BOOL* error = NULL) const;
	double ToDouble(BOOL* error = NULL) const;
#endif
	Str& DeleteLeft(CPOS count);
	Str& DeleteRight(CPOS count);
	Str& Delete(CPOS start, CPOS count = 1);
	void TruncateAt(CPOS idx);
	friend STR_EXPORT Str operator+(const Str& s1, const Str& s2);
	friend STR_EXPORT Str operator+(const Str& s, const STRCHAR *lpsz);
	friend STR_EXPORT Str operator+(const STRCHAR *lpsz, const Str& s);
	friend STR_EXPORT Str operator+(const Str& s, STRCHAR ch);
	friend STR_EXPORT Str operator+(STRCHAR ch, const Str& s);
	Str& TrimLeft(const STRCHAR *charset = NULL);		// Remove all trailing
	Str& TrimRight(const STRCHAR *charset = NULL);		// Remove all leading
	Str& Trim(const STRCHAR *charset = NULL);			// Remove both
	int  Remove(STRCHAR ch);

// String span and tokenization
	static CPOS StringSpanIncluding(const STRCHAR* string, const STRCHAR* tokens);
	static CPOS StringSpanExcluding(const STRCHAR* string, const STRCHAR* tokens);
	Str  SpanIncluding(const STRCHAR* tokens) const;
	Str  SpanExcluding(const STRCHAR* tokens) const;
	Str  Tokenize(const STRCHAR* tokens, CPOS& startat) const;	// Similar to strtok(), initialize startat to 0 in the beginning
#ifndef STR_NO_RTLIBRARY
	Str  Soundex() const;
#endif

// Window operations and other utilities
#ifndef STR_NO_WINSTUFF
	void GetWindowText (HWND wnd);
#endif
#ifndef STR_NO_RTLIBRARY
	void UrlEncode(Str& dest, BOOL spaceAsPlus = FALSE);
#endif

// Convert from internal repres. to ANSI or UNICODE
#if defined(STR_UNICODE)  &&  !defined(STR_NO_RTLIBRARY)
	STRACHAR* ToUTF8(STRACHAR* buffer, int buflen) const;
	void FromUTF8(const STRACHAR* buffer, int length = -1);
#endif
#ifndef STR_NO_UNICODE
	STRACHAR* GetStringA() const;		// Free result with OS_free()
	STRWCHAR* GetStringW() const;		// Free result with OS_free()
#endif
// No conversion -- exactly like typecast
	const STRCHAR* GetString() const;		// Do NOT free result!

// A helper instance of a single-threaded empty string
#ifndef STR_NO_RTLIBRARY
	static const Str& EmptyStr;
#endif

// Miscellaneous implementation methods
public:
	// These two might be used on occasion by callers
	static void* OS_malloc(size_t size);
	static void  OS_free(void* mem);
	// Access private block
	StrImplementation::SBlock* block() const;
	// Block manipulation helpers
	static STRCHAR* block_s (StrImplementation::SBlock* b);
	void li_spawn(CPOS newlength, BOOL copyspawn, StrImplementation::syncstr* caller);
	void FormatVaArg(const STRCHAR *fmt, va_list& marker, BOOL rev_polarity);
protected:
	// Higher level -- miscellaneous
	void NewFromString(const STRCHAR* s, CPOS slen, CPOS prealloc, BOOL want_mt);
	void NewFromString(const STRCHAR* s, CPOS prealloc, BOOL want_mt);
	void CoreAppendChars(const STRCHAR* s, CPOS howmany, StrImplementation::verifymt& mtcheck);
	void Copy_Native(const STRCHAR* s);
	void Str_MakeUL(BOOL do_upr);
	void ImpTrimLeft(const STRCHAR *charset);
	void ImpTrimRight(const STRCHAR *charset);
	// Higher level -- formatting
	void FormatCore (const STRCHAR* x, StrImplementation::verifymt& mtcheck, va_list& marker, BOOL rev_polarity);
	BOOL FmtOneValue (const STRCHAR*& x, StrImplementation::verifymt& mtcheck, va_list& marker, BOOL rev_polarity);
	void FmtPriorBuffer (const STRCHAR* x, StrImplementation::verifymt& mtcheck);
	// UNICODE <-> ANSI stuff
#ifndef STR_NO_UNICODE
#ifdef STR_UNICODE
	void Copy_FromAnsi(const STRACHAR* s);
#else
	void Copy_FromUni (const STRWCHAR* s);
#endif
#endif
	// Low level
	void SetMT_Internal();
public:
	typedef enum {
		SE_OK = 0,
		//
		// Str and Char class basic errors
		SE_BadCharPos,				// SetAt() with bad char index
		SE_BadScanPos,				// Bad char index / char count in Find, Replace, Insert...
		SE_StringEmpty,				// GetFirst/LastChar when string is empty
		SE_UnicodeError,			// Unicode conv functions failed for some reason
		SE_BadFmtSpecifier,			// Bad format/precision specifier in Format/FormatRes
		SE_UnsupFmtSpecifier,		// Unsupported specifier in Format() family
		SE_NoNullChar,				// Null character passed where it is not allowed
		SE_BadMTUse,				// Attempt to modify non-MT object on different thread 
		SE_BadUnicodeUse,			// Need Unicode at runtime, but STR_NO_UNICODE turned on
		SE_BadNumber,				// ToInt / ToFloat / ToDouble could not convert properly
		SE_MalformedUtf8Char,		// ToUTF8 / FromUTF8 cannot perform a conversion
		SE_Reserved3,
		SE_OutOfMemory,				// An API or method reported 'no memory'
		//
		// Regexp and RegexpMatch errors
		SE_RX_Brace_Expected,		// A closing brace was expected
		SE_RX_Paren_Expected,		// A closing parenthesis was expected
		SE_RX_Bracket_Expected,		// A closing bracket was expected
		SE_RX_Unexpected,			// An unspecified fatal error occurred
		SE_RX_Empty_Range,			// A range expression was empty
		SE_RX_Invalid_Group,		// A backreference was made to a non-existing group
		SE_RX_Invalid_Range,		// An invalid range was specified
		SE_RX_Empty_RepeatOp,		// A possibly empty * or + was detected
		//
		// StrArray errors
		SE_ARR_BadFormat,			// Bad VT_ format specified in conversion function
		SE_ARR_BadElementIndex,		// Bad array element index or count
		SE_ARR_COMError				// Miscellaneous error reported by COM API function
		//
		} StrEC;
	// This may be reimplemented in user code -- see documentation
	void Error(StrEC ecode) const;
	static void ErrorNoObject(StrEC ecode);

// Own user methods support
#ifdef STR_OWN_METHODS
public:
#include "str_own_methods.h"
#endif
};

BOOL STR_EXPORT operator ==(const Str& s1, const Str& s2);
BOOL STR_EXPORT operator ==(const Str& s1, const STRCHAR* s2);
BOOL operator ==(const STRCHAR* s1, const Str& s2);
BOOL operator !=(const Str& s1, const Str& s2);
BOOL operator !=(const Str& s1, const STRCHAR* s2);
BOOL operator !=(const STRCHAR* s1, const Str& s2);
BOOL operator <(const Str& s1, const Str& s2);
BOOL operator <(const Str& s1, const STRCHAR* s2);
BOOL operator <(const STRCHAR* s1, const Str& s2);
BOOL operator <=(const Str& s1, const Str& s2);
BOOL operator <=(const Str& s1, const STRCHAR* s2);
BOOL operator <=(const STRCHAR* s1, const Str& s2);
BOOL operator >(const Str& s1, const Str& s2);
BOOL operator >(const Str& s1, const STRCHAR* s2);
BOOL operator >(const STRCHAR* s1, const Str& s2);
BOOL operator >=(const Str& s1, const Str& s2);
BOOL operator >=(const Str& s1, const STRCHAR* s2);
BOOL operator >=(const STRCHAR* s1, const Str& s2);

#ifndef STR_NO_EXCEPTIONS
class STR_EXPORT StrException {
public:
	Str			m_Object;	// A copy of the string that caused the error, or empty string if N/A
	Str::StrEC	m_Error;	// Error code
};
#endif


// System shutdown flag (effectively disables keep-free-mem mechanism)
namespace StrImplementation {
extern BOOL STR_EXPORT STR_Shutting;
}


/*********************************************************************
* Classes:	Regexp, RegexpMatch
* Purpose:	Regular expressions support
*********************************************************************/

#ifdef STR_USE_REGEX

class RegexpInner;
class RegexpMatchInner;

class STR_EXPORT Regexp;

class STR_EXPORT RegexpMatch
{
public:
	friend class Regexp;
	friend class RegexpInner;
	friend class RegexpMatchInner;

	RegexpMatch(int stack = -1);
	~RegexpMatch();

	void GetCompleteMatch(const STRCHAR **szStart, const STRCHAR **szEnd);
	Str  GetCompleteMatch();
	int GetCount() const					{ return m_uNumGroups; }
	void GetMatch(int index, const STRCHAR **szStart, const STRCHAR **szEnd);
	Str  GetMatch(int index);

protected:
	struct MatchGroup
	{
		const STRCHAR *szStart;
		const STRCHAR *szEnd;
	};

	UINT m_uNumGroups;
	MatchGroup m_CompleteMatch;

	RegexpMatchInner* m_Inner;
	int stackptr;

	void Initialize(UINT uRequiredMem, UINT uNumGroups);
	BOOL Push(void *p);
	BOOL Push(int n);
	void *Pop();
};

class STR_EXPORT Regexp
{
public:
	Regexp();
	~Regexp();
	void Parse(const STRCHAR *re, BOOL bCaseSensitive=TRUE);
	BOOL Match(const STRCHAR *str, RegexpMatch* dest, const STRCHAR **end_at=NULL);

protected:
	Str::StrEC m_Error;
	void ThrowError(const STRCHAR** rerror = NULL);
	void Reset();

	enum REInstructionType { 
		RX_ADVANCE,
		RX_ANY,
		RX_CALL,
		RX_FAIL,
		RX_GET_CHARPOS,
		RX_GET_STACKPOS,
		RX_GROUP_END, 
		RX_GROUP_START,
		RX_JMP,
		RX_MATCH,
		RX_NG_PLUS,
		RX_NG_QUESTION,
		RX_NG_STAR_BEGIN, 
		RX_NOP,
		RX_NOTRANGE,
		RX_NOTRANGE_EX,
		RX_PLUS,
		RX_POP_CHARPOS,
		RX_POP_GROUP,
		RX_POP_MEMORY,
		RX_PREVIOUS,
		RX_PUSH_CHARPOS,
		RX_PUSH_GROUP,
		RX_PUSH_MEMORY,
		RX_QUESTION,
		RX_RANGE,
		RX_RANGE_EX,
		RX_RET_NOMATCH,
		RX_RETURN,
		RX_STAR_BEGIN,
		RX_STORE_CHARPOS,
		RX_STORE_STACKPOS,
		RX_SYMBOL
	};

	struct INSTRUCTION_SYMBOL
	{
		int nSymbol;
	};

	struct INSTRUCTION_JMP
	{
		int nTarget;	
	};

	struct INSTRUCTION_GROUP
	{
		int nGroup;
	};

	struct INSTRUCTION_CALL
	{
		int nTarget;
	};

	struct INSTRUCTION_MEMORY
	{
		int nIndex;
	};

	struct INSTRUCTION_PREVIOUS
	{
		int nGroup;
	};

	struct INSTRUCTION_RANGE_EX
	{
		int nTarget;
	};

	int InstructionsPerRangeBitField();

public:
	struct INSTRUCTION
	{
		REInstructionType type;
		union
		{
			INSTRUCTION_SYMBOL symbol;
			INSTRUCTION_JMP jmp;
			INSTRUCTION_GROUP group;
			INSTRUCTION_CALL call;
			INSTRUCTION_MEMORY memory;
			INSTRUCTION_PREVIOUS prev;
			INSTRUCTION_RANGE_EX range;
		};
	};

protected:
	RegexpInner* m_Inner;
	UINT m_uNumGroups;
	UINT m_uRequiredMem;
	BOOL m_bCaseSensitive;

	// class used internally to restore parsing state when unwinding
	class ParseBacktrack
	{
	public:
		int m_nNumInstructions;
		UINT m_uNumGroups;
		UINT m_uRequiredMem;

		ParseBacktrack(Regexp *pRegExp);
		void Restore(Regexp *pRegExp);
	};
	friend class ParseBacktrack;

	INSTRUCTION& GetInstruction(int nIndex);
	int AddInstruction(REInstructionType type);
	int AddInstructions(int nNumInstructions);
	int AddMemInstruction(REInstructionType type);
	BOOL PeekToken(const STRCHAR **ppszRE, int ch);
	BOOL MatchToken(const STRCHAR **ppszRE, int ch);
	STRCHAR GetEscapedChar(STRCHAR ch);
	int ParseArg(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseGroup(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseCharItem(const STRCHAR **ppszRE, STRCHAR *pchStartChar, STRCHAR *pchEndChar);
	int ParseCharSet(const STRCHAR **ppszRE, BOOL bNot);
	int ParseCharClass(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseNot(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseAbbrev(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseSE(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseE(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseAltE(const STRCHAR **ppszRE, BOOL &bEmpty);
	int ParseRE(const STRCHAR **re, BOOL &bEmpty);
	void FixupMatchPos(RegexpMatch* pContext, const STRCHAR *szOrig, const STRCHAR *szNew);
	BOOL MatchCore(const STRCHAR *szInput, RegexpMatch* dest);
};

#endif		// ifdef STR_USE_REGEX


/*********************************************************************
* Classes:	StrArray
* Purpose:	Various utility classes for working with Str
*********************************************************************/

#ifdef STR_USE_EXTRAS

class STR_EXPORT StrArray
{
public:
	StrArray();
	~StrArray();
	int GetSize() const;
	int GetCount() const;
	BOOL IsEmpty() const;
	int GetUpperBound() const;
	void SetSize(int newsize, int growby = -1);
	void FreeExtra();
	const Str& GetAt(int index) const;
	void SetAt(int index, const STRCHAR* elem);
	void SetAt(int index, const Str& elem);
	Str& ElementAt(int index);
	void SetAtGrow(int index, const STRCHAR* elem);
	void SetAtGrow(int index, const Str& elem);
#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)
	void SetAt(int index, const CString& elem);
	void SetAtGrow(int index, const CString& elem);
	int Add(const CString& elem);
	int Append(const CStringArray& src);
	void InsertAt(int index, const CString& newElement, int nCount = 1);
#endif
	int Add(const STRCHAR* elem);
	int Add(const Str& elem);
	int Append(const StrArray& src);
	void Copy(const StrArray& src);
	const Str& operator[](int index) const;
	Str& operator[](int index);
	void InsertAt(int index, const STRCHAR* elem, int nCount = 1);
	void InsertAt(int index, const Str& newElement, int nCount = 1);
	void InsertAt(int startIndex, const StrArray& newArray);
	void RemoveAt(int index, int nCount = 1);
	void RemoveAll();
	void CopyElements(int destStartIndex, int srcStartIndex, int count, const StrArray& srcArray);
	/* Unsupported methods:
	Str* GetData();				- no equivalent facility in <veclite>
	const Str* GetData() const; - no equivalent facility in <veclite>
	*/
#if defined(STR_WIN32) && !defined(STR_NO_WINSTUFF) && !defined(STR_NO_UNICODE)
	SAFEARRAY* ToSafeArray(VARTYPE destVtType, int startIndex = 0, int nCount = -1);
#endif

// Implementation
protected:
	void MaybeGrow (int toMaxIndex);
	int  GrowRel (int additionalCount);
	int  m_GrowBy;		// Not directly supported by <veclite> so we need to emulate it
public:
	// This member is public for implementation reasons.  Never use from client code!
	void* m_Data;		// Actually veclite<Str>
};

#endif // STR_USE_EXTRAS


/////////////////////////////////////////////
// Inline implementations of numerous methods
#include "str_support2.h"

#ifdef STR_BORCPPBUILDER
#undef Char
#endif

#endif  // _STR_H
