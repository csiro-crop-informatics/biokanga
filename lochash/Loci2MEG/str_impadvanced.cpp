// Str_ImpAdvanced.cpp
//
// Version 2.4.0
//
// This file is a part of the Str Library, and is subject to the
// clauses of the License Agreement accompanying the software.

//!!stdafx

#include "str.h"

#if !defined(STR_USE_REGEX)  &&  !defined(STR_USE_EXTRAS)
#error You must define STR_USE_REGEX or STR_USE_EXTRAS to use these features
#endif

#ifdef STR_BORCPPBUILDER
#define Char StrChar
#pragma warn -8004
#pragma warn -8008
#pragma warn -8066
#endif


/*********************************************************************
* Classes:	TAutoVectorPtr
* Purpose:	Portable equivalent to ATL CAutoVectorPtr
*********************************************************************/

namespace StrImplementation {

template< typename T >
class TAutoVectorPtr
{
public:
	TAutoVectorPtr() : m_p(NULL)				{ }
	TAutoVectorPtr( TAutoVectorPtr< T >& p )	{ m_p = p.Detach(); }
	explicit TAutoVectorPtr( T* p ) : m_p(p)	{ }
	~TAutoVectorPtr()							{ Free(); }
	operator T*() const							{ return m_p; }

	TAutoVectorPtr< T >& operator=( TAutoVectorPtr< T >& p )
	{
		Free();
		Attach( p.Detach() );
		return *this;
	}

	// Allocate the vector
	void Allocate( size_t nElements )
	{
		STR_ASSERT (m_p == NULL);
		m_p = new T[nElements];
	}

	// Attach to an existing pointer (takes ownership)
	void Attach( T* p )
	{
		STR_ASSERT (m_p == NULL);
		m_p = p;
	}

	// Detach the pointer (releases ownership)
	T* Detach()
	{
		T* p;

		p = m_p;
		m_p = NULL;

		return p;
	}

	// Delete the vector pointed to, and set the pointer to NULL
	void Free()
	{
		delete[] m_p;
		m_p = NULL;
	}

public:
	T* m_p;
};

}		// namespace


using namespace StrImplementation;


/*********************************************************************
* Regex global variables
*********************************************************************/

#ifdef STR_USE_REGEX

#define STR_REGEXP_STACKDEPTH 256

static const STRCHAR* RegexpAbbrevs[] = 
{
	_T("a([a-zA-Z0-9])"),	// alpha numeric
	_T("b([ \\t])"),		// white space (blank)
	_T("c([a-zA-Z])"),		// alpha
	_T("d([0-9])"),			// digit
	_T("h([0-9a-fA-F])"),	// hex digit
	_T("n(\r|(\r?\n))"),	// newline
	_T("q(\"[^\"]*\")|(\'[^\']*\')"),	// quoted string
	_T("w([a-zA-Z]+)"),		// simple word
	_T("z([0-9]+)"),		// unsigned integer
	NULL
};

class RegexpMatchInner
{
public:
	veclite<void *> m_stack;
	TAutoVectorPtr<void *> m_Mem;
	TAutoVectorPtr<RegexpMatch::MatchGroup> m_Matches;
};

RegexpMatch::RegexpMatch(int nInitStackSize /* = -1*/)
{
	if (nInitStackSize == -1)
		nInitStackSize = STR_REGEXP_STACKDEPTH;
	m_uNumGroups = 0;
	stackptr = 0;
	m_Inner = new RegexpMatchInner;
	m_Inner->m_stack.resize(nInitStackSize);
	m_CompleteMatch.szStart = NULL;
	m_CompleteMatch.szEnd = NULL;
}

RegexpMatch::~RegexpMatch()
{
	delete m_Inner;
}

inline int Regexp::InstructionsPerRangeBitField()
{
	return (256/8) / sizeof(INSTRUCTION) + ((256/8) % sizeof(INSTRUCTION) ? 1 : 0);
}

inline STRCHAR Regexp::GetEscapedChar(STRCHAR ch)
{
	if (ch == 't')
		return '\t';
	return ch;
}

inline BOOL Regexp::PeekToken(const STRCHAR **ppszRE, int ch)
{
	if (**ppszRE != (STRCHAR) ch)
		return FALSE;
	return TRUE;
}

inline BOOL Regexp::MatchToken(const STRCHAR **ppszRE, int ch)
{
	if (!PeekToken(ppszRE, ch))
		return FALSE;
	*ppszRE += 1;
	return TRUE;
}

void RegexpMatch::GetMatch(int nIndex, const STRCHAR **szStart, const STRCHAR **szEnd)
{
	STR_ASSERT(szStart != NULL);
	STR_ASSERT(szEnd != NULL);
	*szStart = m_Inner->m_Matches[nIndex].szStart;
	*szEnd = m_Inner->m_Matches[nIndex].szEnd;
}

Str RegexpMatch::GetMatch(int nIndex)
{
	const STRCHAR* start;
	const STRCHAR* end;
	GetMatch(nIndex, &start, &end);
	Str dest;
	dest.AppendChars(start, (CPOS) (end-start));
	return dest;
}

void RegexpMatch::Initialize(UINT uRequiredMem, UINT uNumGroups)
{
	stackptr = 0;
	m_uNumGroups = 0;
	m_Inner->m_Matches.Free();
	m_Inner->m_Matches.Allocate(uNumGroups);
	m_uNumGroups = uNumGroups;
	m_Inner->m_Mem.Free();
	m_Inner->m_Mem.Allocate(uRequiredMem);
	memset(m_Inner->m_Mem.m_p, 0x00, uRequiredMem*sizeof(void *));
	memset(m_Inner->m_Matches, 0x00, m_uNumGroups * sizeof(MatchGroup));
}

BOOL RegexpMatch::Push(int n)
{
	stackptr++;
	if (m_Inner->m_stack.size() <= (UINT) stackptr)
		m_Inner->m_stack.resize((stackptr+1)*2);
	m_Inner->m_stack[stackptr] = (void*) (size_t) n;
	return TRUE;
}

// Can't inline in H file because of warn pragmas 4311 / 4312
inline BOOL RegexpMatch::Push(void* p)
{
	return Push((int) (size_t) p);
}

void* RegexpMatch::Pop()
{
	STR_ASSERT(stackptr > 0);	// Stack underflow - should never happen
	void* p = m_Inner->m_stack[stackptr];
	stackptr--;
	return p;
}

class RegexpInner
{
public:
	veclite<Regexp::INSTRUCTION> m_Instructions;
};

Regexp::ParseBacktrack::ParseBacktrack(Regexp *re)
{
	m_nNumInstructions = (int) re->m_Inner->m_Instructions.size();
	m_uNumGroups = re->m_uNumGroups;
	m_uRequiredMem = re->m_uRequiredMem;
}

void Regexp::ParseBacktrack::Restore(Regexp *re)
{
	re->m_Inner->m_Instructions.resize(m_nNumInstructions);
	re->m_uNumGroups = m_uNumGroups;
	re->m_uRequiredMem = m_uRequiredMem;
}

Regexp::Regexp()
{
	m_uNumGroups = 0;
	m_uRequiredMem = 0;
	m_bCaseSensitive = TRUE;
	m_Error = Str::SE_OK;
	m_Inner = new RegexpInner;
}

Regexp::~Regexp()
{
	delete m_Inner;
}

// Can't be inline because of STR_EXPORT
Regexp::INSTRUCTION& Regexp::GetInstruction(int nIndex)
{
	return m_Inner->m_Instructions[nIndex];
}

void Regexp::Parse(const STRCHAR *re, BOOL bCaseSensitive/*=TRUE*/)
{
	Reset();

	m_bCaseSensitive = bCaseSensitive;

	BOOL release_szin = FALSE;
	const STRCHAR *szInput = re;

	try {
		if (!bCaseSensitive) {
			// copy the string
			int len = STR_strlen(re);
			int bytesize = (len + 1) * sizeof(STRCHAR);
			szInput = (const STRCHAR *) Str::OS_malloc(bytesize);
			release_szin = TRUE;
			STR_copy((char *) szInput, re, bytesize);
			// lowercase the copy
			STR_ToLower ((STRCHAR*) szInput, len);
		}
		const STRCHAR *sz = szInput;

		int nCall = AddInstruction(RX_CALL);
		if (*sz == '^') {
			AddInstruction(RX_FAIL);
			sz++;
		}
		else
			AddInstruction(RX_ADVANCE);

		BOOL bEmpty = true;
		ParseRE(&sz, bEmpty);
		STR_ASSERT(m_Error == Str::SE_OK);
		GetInstruction(nCall).call.nTarget = 2;
		AddInstruction(RX_MATCH);
	}
	catch (...) {
		if (release_szin)
			Str::OS_free((void *) szInput);
		throw;
	}

	if (release_szin)
		Str::OS_free((void *) szInput);

	if (m_Error)
		ThrowError(&re);
}

void Regexp::ThrowError(const STRCHAR** rerror)
{
	STR_ASSERT(m_Error);
	Str error_in;
	if (rerror)
		error_in = *rerror;
	error_in.Error(m_Error);
}

BOOL Regexp::Match(const STRCHAR *str, RegexpMatch* con, const STRCHAR **end_at /*=NULL*/)
{
	STR_ASSERT(str);
	STR_ASSERT(con);

	if (end_at)
		*end_at = NULL;

	BOOL release_szin = FALSE;
	const STRCHAR *szInput = str;
	BOOL core_result = FALSE;		// Keep certain compilers happy - assign a value!

	try {
		if (!m_bCaseSensitive) {
			// copy the string
			int len = STR_strlen(str);
			int bytesize = (len + 1) * sizeof(STRCHAR);
			szInput = (const STRCHAR *) Str::OS_malloc(bytesize);
			release_szin = TRUE;
			STR_copy((char *) szInput, str, bytesize);
			// lowercase the copy
			STR_ToLower ((STRCHAR*) szInput, len);
		}
		// Call the real matcher
		core_result = MatchCore(szInput, con);
		// Fixup output results
		const STRCHAR* sz = con->m_CompleteMatch.szEnd;		// MatchCore left it here
		if (!m_bCaseSensitive)
			FixupMatchPos(con, str, szInput);
		if (end_at)
			*end_at = str + (sz - szInput);
		goto RetResult;
	}
	catch (...) {
		if (release_szin)
			Str::OS_free((void*) szInput);
		throw;
	}
RetResult:
	if (release_szin)
		Str::OS_free((void*) szInput);
	return core_result;
}

BOOL Regexp::MatchCore(const STRCHAR *szInput, RegexpMatch* con)
{
	con->Initialize(m_uRequiredMem, m_uNumGroups);
	int ip = 0;

	const STRCHAR *sz = szInput;
	const STRCHAR *szCurrInput = szInput;

	for (;;) {
		if (ip == 0)
			con->m_CompleteMatch.szStart = sz;

		switch (GetInstruction(ip).type)
 		{
		case RX_NOP:
			ip++;
			break;

		case RX_SYMBOL:
			if (GetInstruction(ip).symbol.nSymbol == (unsigned char)*sz)
			{
				sz++;
				ip++;
			}
			else
			{
				ip = (int) (size_t) con->Pop();
			}
			break;

		case RX_ANY:
			if (*sz)
			{
				sz++;
				ip++;
			}
			else
			{
				ip = (int) (size_t) con->Pop();
			}
			break;

		case RX_GROUP_START:
			con->m_Inner->m_Matches[GetInstruction(ip).group.nGroup].szStart = sz;
			ip++;
			break;

		case RX_GROUP_END:
			con->m_Inner->m_Matches[GetInstruction(ip).group.nGroup].szEnd = sz;
			ip++;
			break;

		case RX_PUSH_CHARPOS:
			con->Push((void *) sz);
			ip++;
			break;

		case RX_POP_CHARPOS:
			sz = (STRCHAR *) con->Pop();
			ip++;
			break;

		case RX_CALL:
			con->Push(ip+1);
			ip = GetInstruction(ip).call.nTarget;
			break;

		case RX_JMP:
			ip = GetInstruction(ip).jmp.nTarget;
			break;

		case RX_RETURN:
			ip = (int) (size_t) con->Pop();
			break;

		case RX_PUSH_MEMORY:
			con->Push((void *) (con->m_Inner->m_Mem[GetInstruction(ip).memory.nIndex]));
			ip++;
			break;

		case RX_POP_MEMORY:
			con->m_Inner->m_Mem[GetInstruction(ip).memory.nIndex] = con->Pop();
			ip++;
			break;

		case RX_STORE_CHARPOS:
			con->m_Inner->m_Mem[GetInstruction(ip).memory.nIndex] = (void *) sz;
			ip++;
			break;

		case RX_GET_CHARPOS:
			sz = (STRCHAR *) con->m_Inner->m_Mem[GetInstruction(ip).memory.nIndex];
			ip++;
			break;

		case RX_STORE_STACKPOS:
			con->m_Inner->m_Mem[GetInstruction(ip).memory.nIndex] = (void *) (size_t) con->stackptr;
			ip++;
			break;

		case RX_GET_STACKPOS:
			con->stackptr = (int) (size_t) con->m_Inner->m_Mem[GetInstruction(ip).memory.nIndex];
			ip++;
			break;

		case RX_RET_NOMATCH:
			if (sz == (STRCHAR *) con->m_Inner->m_Mem[GetInstruction(ip).memory.nIndex])
			{
				// do a return
				ip = (int) (size_t) con->Pop();
			}
			else
				ip++;
			break;

		case RX_ADVANCE:
			sz = szCurrInput + 1;
			szCurrInput = sz;
			if (*sz == '\0')
				goto Error;
			ip = 0;
			con->stackptr = 0;
			break;

		case RX_FAIL:
			goto Error;

		case RX_RANGE:
			{
				if (*sz == '\0')
				{
					ip = (int) (size_t) con->Pop();
					break;
				}

				unsigned char *pBits = (unsigned char *) (&m_Inner->m_Instructions[ip]+1);
				size_t u = (size_t) *sz;
				if (pBits[u >> 3] & 1 << (u & 0x7))
				{
					ip += InstructionsPerRangeBitField();
					ip++;
					sz++;
				}
				else
				{
					ip = (int) (size_t) con->Pop();
				}
			}
			break;

		case RX_NOTRANGE:
			{
				if (*sz == '\0')
				{
					ip = (int) (size_t) con->Pop();
					break;
				}

				unsigned char *pBits = (unsigned char *) (&m_Inner->m_Instructions[ip]+1);
				size_t u = (size_t) * ((unsigned char *) sz);
				if (pBits[u >> 3] & 1 << (u & 0x7))
				{
					ip = (int) (size_t) con->Pop();
				}
				else
				{
					ip += InstructionsPerRangeBitField();
					ip++;
					sz++;
				}
			}
			break;

		case RX_RANGE_EX:
			{
				if (*sz == '\0')
				{
					ip = (int) (size_t) con->Pop();
					break;
				}

				BOOL bMatch = FALSE;
				int inEnd = GetInstruction(ip).range.nTarget;
				ip++;

				while (ip < inEnd)
				{						
					if ((unsigned char)*sz >= GetInstruction(ip).memory.nIndex && 
						(unsigned char)*sz <= GetInstruction(ip+1).memory.nIndex)
					{
						// if we match, we jump to the end
						sz++;
						ip = inEnd;
						bMatch = TRUE;
					}
					else
					{
						ip += 2;
					}
				}
				if (!bMatch)
				{
					ip = (int) (size_t) con->Pop();
				}
			}
			break;

		case RX_NOTRANGE_EX:
			{
				if (*sz == '\0')
				{
					ip = (int) (size_t) con->Pop();
					break;
				}

				BOOL bMatch = TRUE;
				int inEnd = GetInstruction(ip).range.nTarget;
				ip++;

				while (ip < inEnd)
				{
					if ((unsigned char)*sz >= GetInstruction(ip).memory.nIndex && 
						(unsigned char)*sz <= GetInstruction(ip+1).memory.nIndex)
					{
						ip = (int) (size_t) con->Pop();
						bMatch = FALSE;
						break;
					}
					else
					{
						// if we match, we jump to the end
						ip += 2;
					}
				}
				if (bMatch)
					sz++;
			}
			break;

		case RX_PREVIOUS:
			{
				BOOL bMatch = FALSE;
				const STRCHAR* ssta = con->m_Inner->m_Matches[GetInstruction(ip).prev.nGroup].szStart;
				const STRCHAR* send = con->m_Inner->m_Matches[GetInstruction(ip).prev.nGroup].szEnd;
				if (m_bCaseSensitive)
				{
					bMatch = !STR_strncmp(sz, ssta, (int) (send - ssta));
				}
				else
				{
					bMatch = !STR_strnicmp(sz, ssta, (int) (send - ssta));
				}
				if (bMatch)
				{
					sz += send - ssta;
					ip++;
					break;
				}
				ip = (int) (size_t) con->Pop();
			}
			break;

		case RX_MATCH:
			con->m_CompleteMatch.szEnd = sz;
			return TRUE;

		case RX_PUSH_GROUP:
			con->Push((void *) con->m_Inner->m_Matches[GetInstruction(ip).group.nGroup].szStart);
			con->Push((void *) con->m_Inner->m_Matches[GetInstruction(ip).group.nGroup].szEnd);
			ip++;
			break;

		case RX_POP_GROUP:
			con->m_Inner->m_Matches[GetInstruction(ip).group.nGroup].szEnd = (const STRCHAR *) con->Pop();
			con->m_Inner->m_Matches[GetInstruction(ip).group.nGroup].szStart = (const STRCHAR *) con->Pop();
			ip++;
			break;

		default:
			STR_ASSERT_FALSE();
			break;
		}

	}

Error:
	con->m_CompleteMatch.szEnd = sz;
	return FALSE;
}

void Regexp::Reset()
{
	m_Inner->m_Instructions.resize(0);
	m_uRequiredMem = 0;
	m_bCaseSensitive = TRUE;
	m_uNumGroups = 0;
	m_Error = Str::SE_OK;
}

int Regexp::AddInstruction(REInstructionType type)
{
	m_Inner->m_Instructions.resize(m_Inner->m_Instructions.size()+1);
	m_Inner->m_Instructions[m_Inner->m_Instructions.size()-1].type = type;
	return (int) m_Inner->m_Instructions.size()-1;
}

int Regexp::ParseArg(const STRCHAR **ppszRE, BOOL &bEmpty)
{
	int nPushGroup = AddInstruction(RX_PUSH_GROUP);
	if (nPushGroup < 0)
		return -1;

	GetInstruction(nPushGroup).group.nGroup = m_uNumGroups;

	int p = AddInstruction(RX_GROUP_START);
	if (p < 0)
		return -1;
	GetInstruction(p).group.nGroup = m_uNumGroups++;

	int nCall = AddInstruction(RX_CALL);
	if (nCall < 0)
		return -1;

	int nPopGroup = AddInstruction(RX_POP_GROUP);
	if (nPopGroup < 0)
		return -1;
	GetInstruction(nPopGroup).group.nGroup = GetInstruction(nPushGroup).group.nGroup;

	if (AddInstruction(RX_RETURN) < 0)
		return -1;

	int nAlt = ParseRE(ppszRE, bEmpty);
	if (nAlt < 0)
	{
		if (!PeekToken(ppszRE, '}')) {
			m_Error = Str::SE_RX_Brace_Expected;
			return -1;
		}

		// in the case of an empty group, we add a nop
		nAlt = AddInstruction(RX_NOP);
		if (nAlt < 0)
			return -1;
	}

	GetInstruction(nCall).call.nTarget = nAlt;

	if (!MatchToken(ppszRE, '}'))
	{
		m_Error = Str::SE_RX_Brace_Expected;
		return -1;
	}

	int nEnd = AddInstruction(RX_GROUP_END);
	if (nEnd < 0)
		return -1;
	GetInstruction(nEnd).group.nGroup = GetInstruction(p).group.nGroup;
	return nPushGroup;
}

// ParseGroup: parse grammar rule Group
int Regexp::ParseGroup(const STRCHAR **ppszRE, BOOL &bEmpty)
{
	int nCall = AddInstruction(RX_CALL);
	if (nCall < 0)
		return -1;

	if (AddInstruction(RX_RETURN) < 0)
		return -1;

	int nAlt = ParseRE(ppszRE, bEmpty);
	if (nAlt < 0) {
		if (!PeekToken(ppszRE, ')'))
		{
			m_Error = Str::SE_RX_Paren_Expected;
			return -1;
		}

		// in the case of an empty group, we add a nop
		nAlt = AddInstruction(RX_NOP);
		if (nAlt < 0)
			return -1;
	}

	GetInstruction(nCall).call.nTarget = nAlt;

	if (!MatchToken(ppszRE, ')'))
	{
		m_Error = Str::SE_RX_Paren_Expected;
		return -1;
	}

	return nCall;
}

// ParseCharItem: parse grammar rule CharItem
int Regexp::ParseCharItem(const STRCHAR **ppszRE, STRCHAR *pchStartChar, STRCHAR *pchEndChar)
{
	if (**ppszRE == '\\')
	{
		*ppszRE += 1;
		*pchStartChar = GetEscapedChar(**ppszRE);
	}
	else
		*pchStartChar = **ppszRE;
	*ppszRE += 1;

	if (!MatchToken(ppszRE, '-'))
	{
		*pchEndChar = *pchStartChar;
		return 0;
	}

	// check for unterminated range
	if (!**ppszRE || PeekToken(ppszRE, ']'))
	{
		m_Error = Str::SE_RX_Bracket_Expected;
		return -1;
	}

	*pchEndChar = **ppszRE;
	*ppszRE += 1;

	if (*pchEndChar < *pchStartChar)
	{
		m_Error = Str::SE_RX_Invalid_Range;
		return -1;
	}
	return 0;
}

int Regexp::AddInstructions(int nNumInstructions)
{
	size_t nCurr = m_Inner->m_Instructions.size();
	m_Inner->m_Instructions.resize(nCurr+nNumInstructions);
	return (int) nCurr;
}

// ParseCharSet: parse grammar rule CharSet
int Regexp::ParseCharSet(const STRCHAR **ppszRE, BOOL bNot)
{
	int p = AddInstruction(bNot ? RX_NOTRANGE_EX : RX_RANGE_EX);
	if (p < 0)
		return -1;

	STRCHAR chStart;
	STRCHAR chEnd;

	while (**ppszRE && **ppszRE != ']')
	{
		if (ParseCharItem(ppszRE, &chStart, &chEnd))
			return -1;

		int nStart = AddInstruction(RX_NOP);
		if (nStart < 0)
			return -1;

		int nEnd = AddInstruction(RX_NOP);
		if (nEnd < 0)
			return -1;

		GetInstruction(nStart).memory.nIndex = (int) chStart;
		GetInstruction(nEnd).memory.nIndex = (int) chEnd;
	}

	GetInstruction(p).range.nTarget = (int) m_Inner->m_Instructions.size();

	return p;
}

// ParseCharClass: parse grammar rule CharClass
int Regexp::ParseCharClass(const STRCHAR **ppszRE, BOOL &bEmpty)
{
	bEmpty = false;
	if (MatchToken(ppszRE, ']'))
	{
		m_Error = Str::SE_RX_Empty_Range;
		return -1;
	}

	BOOL bNot = FALSE;
	if (MatchToken(ppszRE, '^'))
		bNot = TRUE;

	if (MatchToken(ppszRE, ']'))
	{
		m_Error = Str::SE_RX_Empty_Range;
		return -1;
	}

	int p = ParseCharSet(ppszRE, bNot);
	if (p < 0)
		return p;
	if (!MatchToken(ppszRE, ']'))
	{
		m_Error = Str::SE_RX_Bracket_Expected;
		return -1;
	}

	return p;
}

int Regexp::AddMemInstruction(REInstructionType type)
{
	int p = AddInstruction(type);
	if (p < 0)
		return p;
	GetInstruction(p).memory.nIndex = m_uRequiredMem++;
	return p;
}

// helper for parsing !SE
int Regexp::ParseNot(const STRCHAR **ppszRE, BOOL &bEmpty)
{
	int nStoreCP = AddMemInstruction(RX_STORE_CHARPOS);
	int nStoreSP = AddMemInstruction(RX_STORE_STACKPOS);

	int nCall = AddInstruction(RX_CALL);
	if (nCall < 0)
		return -1;

	int nGetCP = AddInstruction(RX_GET_CHARPOS);
	if (nGetCP < 0)
		return -1;
	GetInstruction(nGetCP).memory.nIndex = GetInstruction(nStoreCP).memory.nIndex;

	int nGetSP = AddInstruction(RX_GET_STACKPOS);
	if (nGetSP < 0)
		return -1;
	GetInstruction(nGetSP).memory.nIndex = GetInstruction(nStoreSP).memory.nIndex;

	int nJmp = AddInstruction(RX_JMP);
	if (nJmp < 0)
		return -1;

	int nSE = ParseSE(ppszRE, bEmpty);
	if (nSE < 0)
		return nSE;

	// patch the call
	GetInstruction(nCall).call.nTarget = nSE;

	int nGetCP1 = AddInstruction(RX_GET_CHARPOS);
	if (nGetCP1 < 0)
		return -1;
	GetInstruction(nGetCP1).memory.nIndex = GetInstruction(nStoreCP).memory.nIndex;

	int nGetSP1 = AddInstruction(RX_GET_STACKPOS);
	if (nGetSP1 < 0)
		return -1;
	GetInstruction(nGetSP1).memory.nIndex = GetInstruction(nStoreSP).memory.nIndex;

	int nRet = AddInstruction(RX_RETURN);
	if (nRet < 0)
		return -1;

	GetInstruction(nJmp).jmp.nTarget = nRet+1;

	return nStoreCP;
}

// ParseAbbrev: parse grammar rule Abbrev
int Regexp::ParseAbbrev(const STRCHAR **ppszRE, BOOL &bEmpty)
{
	const STRCHAR **szAbbrevs = RegexpAbbrevs;

	while (*szAbbrevs)
	{
		if (**ppszRE == **szAbbrevs)
		{
			const STRCHAR *szAbbrev = (*szAbbrevs)+1;
			int p = ParseE(&szAbbrev, bEmpty);
			if (p < 0)
			{
				m_Error = Str::SE_RX_Unexpected;
				return p;
			}
			*ppszRE += 1;
			return p;
		}
		szAbbrevs++;
	}
	return -1;
}

// ParseSE: parse grammar rule SE (simple expression)
int Regexp::ParseSE(const STRCHAR **ppszRE, BOOL &bEmpty)
{

	if (MatchToken(ppszRE, '{'))
		return ParseArg(ppszRE, bEmpty);
	if (MatchToken(ppszRE, '('))
		return ParseGroup(ppszRE, bEmpty);
	if (MatchToken(ppszRE, '['))
		return ParseCharClass(ppszRE, bEmpty);

	if (MatchToken(ppszRE, '\\'))
	{
		Char temp(**ppszRE);
		if (!temp.IsDigit())
		{
			// check for abbreviations
			int p;
			p = ParseAbbrev(ppszRE, bEmpty);
			if (p >= 0)
				return p;

			if (m_Error)
				return -1;

			// escaped char
			p = AddInstruction(RX_SYMBOL);
			if (p < 0)
				return -1;
			GetInstruction(p).symbol.nSymbol = (int) **ppszRE;
			*ppszRE += 1;
			return p;
		}
		// previous match
		bEmpty = false;
		int nPrev = AddInstruction(RX_PREVIOUS);
		if (nPrev < 0)
			return -1;

		UINT uValue = (UINT) STR_strtol(*ppszRE, (STRCHAR **) ppszRE, 10);
		if (uValue >= m_uNumGroups)
		{
			m_Error = Str::SE_RX_Invalid_Group;
			return -1;
		}
		GetInstruction(nPrev).prev.nGroup = (int) uValue;
		return nPrev;
	}

	if (MatchToken(ppszRE, '!'))
		return ParseNot(ppszRE, bEmpty);

	if (**ppszRE == '}' || **ppszRE == ']' || **ppszRE == ')')
	{
		return -1;
	}

	if (**ppszRE == '\0')
	{
		return -1;
	}

	int p;
	if (**ppszRE == '.')
	{
		p = AddInstruction(RX_ANY);
		if (p < 0)
			return -1;
		bEmpty = false;
	}
	else if (**ppszRE == '$' && (*ppszRE)[1] == '\0')
	{
		p = AddInstruction(RX_SYMBOL);
		if (p < 0)
			return -1;
		GetInstruction(p).symbol.nSymbol = 0;
		bEmpty = false;
	}
	else
	{
		p = AddInstruction(RX_SYMBOL);
		if (p < 0)
			return -1;
		GetInstruction(p).symbol.nSymbol = (int) **ppszRE;
		bEmpty = false;
	}
	*ppszRE += 1;
	return p;
}

// ParseE: parse grammar rule E (expression)
int Regexp::ParseE(const STRCHAR **ppszRE, BOOL &bEmpty)
{
	ParseBacktrack ParseBacktrack(this);
	const STRCHAR *sz = *ppszRE;

	int nSE;

	int nFirst = ParseSE(ppszRE, bEmpty);
	if (nFirst < 0)
		return nFirst;

	REInstructionType type = RX_MATCH;

	if (MatchToken(ppszRE, '*'))
		if(MatchToken(ppszRE, '?'))
			type = RX_NG_STAR_BEGIN;
		else
			type = RX_STAR_BEGIN;


	else if (MatchToken(ppszRE, '+'))
		if(MatchToken(ppszRE, '?'))
			type = RX_NG_PLUS;
		else
			type = RX_PLUS;

	else if (MatchToken(ppszRE, '?'))
		if(MatchToken(ppszRE, '?'))
			type = RX_NG_QUESTION;
		else
			type = RX_QUESTION;


	if (type == RX_MATCH)
		return nFirst;

	if (type == RX_STAR_BEGIN || type == RX_QUESTION|| type == RX_NG_STAR_BEGIN || type == RX_NG_QUESTION)
	{
		ParseBacktrack.Restore(this);
	}
	else
	{
		m_uNumGroups = ParseBacktrack.m_uNumGroups;
	}
	*ppszRE = sz;

	int nE;

	if (type == RX_NG_STAR_BEGIN || type == RX_NG_PLUS || type == RX_NG_QUESTION) // Non-Greedy
	{			
		int nCall = AddInstruction(RX_CALL);
		if (nCall < 0)
			return -1;

		bEmpty = false;

		nSE = ParseSE(ppszRE, bEmpty);
		if (nSE < 0)
			return nSE;

		if (bEmpty && (type == RX_NG_STAR_BEGIN || type == RX_NG_PLUS))
		{
			m_Error = Str::SE_RX_Empty_RepeatOp;
			return -1;
		}
		bEmpty = true;

		*ppszRE += 1;
		*ppszRE += 1;

		if (type == RX_NG_STAR_BEGIN || type == RX_NG_PLUS)
		{
			int nJmp = AddInstruction(RX_JMP);
			if (nJmp < 0)
				return -1;
			GetInstruction(nCall).call.nTarget = nJmp+1;
			GetInstruction(nJmp).jmp.nTarget = nCall;
		}
		else
			GetInstruction(nCall).call.nTarget = nSE+1;

		if (type == RX_NG_PLUS)
			nE = nFirst;
		else
			nE = nCall;
	}
	else // Greedy
	{

		int nPushMem = AddInstruction(RX_PUSH_MEMORY);
		if (nPushMem < 0)
			return -1;

		int nStore = AddInstruction(RX_STORE_CHARPOS);
		if (nStore < 0)
			return -1;

		if (AddInstruction(RX_PUSH_CHARPOS) < 0)
			return -1;

		int nCall = AddInstruction(RX_CALL);
		if (nCall < 0)
			return -1;

		if (AddInstruction(RX_POP_CHARPOS) < 0)
			return -1;

		int nPopMem = AddInstruction(RX_POP_MEMORY);
		if (nPopMem < 0)
			return -1;

		int nJmp = AddInstruction(RX_JMP);
		if (nJmp < 0)
			return -1;

		GetInstruction(nPushMem).memory.nIndex = m_uRequiredMem++;
		GetInstruction(nStore).memory.nIndex = GetInstruction(nPushMem).memory.nIndex;
		GetInstruction(nCall).call.nTarget = nJmp+1;
		GetInstruction(nPopMem).memory.nIndex = GetInstruction(nPushMem).memory.nIndex;

		bEmpty = false;

		nSE = ParseSE(ppszRE, bEmpty);
		if (nSE < 0)
			return nSE;

		if (bEmpty && (type == RX_STAR_BEGIN || type == RX_PLUS))
		{
			m_Error = Str::SE_RX_Empty_RepeatOp;
			return -1;
		}

		if (type != RX_PLUS && type != RX_NG_PLUS)
			bEmpty = true;

		*ppszRE += 1;


		int nRetNoMatch = AddInstruction(RX_RET_NOMATCH);
		if (nRetNoMatch < 0)
			return -1;

		int nStore1 = AddInstruction(RX_STORE_CHARPOS);
		if (nStore1 < 0)
			return -1;

		GetInstruction(nRetNoMatch).memory.nIndex = GetInstruction(nPushMem).memory.nIndex;
		GetInstruction(nStore1).memory.nIndex = GetInstruction(nPushMem).memory.nIndex;

		if (type != RX_QUESTION)
		{
			int nJmp1 = AddInstruction(RX_JMP);
			if (nJmp1 < 0)
				return -1;
			GetInstruction(nJmp1).jmp.nTarget = nPushMem;
		}

		GetInstruction(nJmp).jmp.nTarget = (int) m_Inner->m_Instructions.size();
		if (type == RX_PLUS)
			nE = nFirst;
		else
			nE = nPushMem;
	}

	return nE;
}

// ParseAltE: parse grammar rule AltE
int Regexp::ParseAltE(const STRCHAR **ppszRE, BOOL &bEmpty)
{
	const STRCHAR *sz = *ppszRE;
	ParseBacktrack ParseBacktrack(this);

	int nPush = AddInstruction(RX_PUSH_CHARPOS);
	if (nPush < 0)
		return -1;

	int nCall = AddInstruction(RX_CALL);
	if (nCall < 0)
		return -1;

	GetInstruction(nCall).call.nTarget = nPush+4;
	if (AddInstruction(RX_POP_CHARPOS) < 0)
		return -1;

	int nJmpNext = AddInstruction(RX_JMP);
	if (nJmpNext < 0)
		return -1;

	int nE = ParseE(ppszRE, bEmpty);
	if (nE < 0)
	{
		if (m_Error)
			return -1;
		ParseBacktrack.Restore(this);
		return nE;
	}

	int nJmpEnd = AddInstruction(RX_JMP);
	if (nJmpEnd < 0)
		return -1;

	GetInstruction(nJmpNext).jmp.nTarget = nJmpEnd+1;

	if (!MatchToken(ppszRE, '|'))
	{
		ParseBacktrack.Restore(this);
		*ppszRE = sz;

		return ParseE(ppszRE, bEmpty);
	}

	BOOL bEmptyAltE;
	int nAltE = ParseAltE(ppszRE, bEmptyAltE);
	GetInstruction(nJmpEnd).jmp.nTarget = (int) m_Inner->m_Instructions.size();
	GetInstruction(nJmpNext).jmp.nTarget = nAltE;
	if (nAltE < 0)
	{
		if (m_Error)
			return -1;
		ParseBacktrack.Restore(this);
		return nAltE;
	}
	bEmpty = bEmpty | bEmptyAltE;
	return nPush;
}

// ParseRE: parse grammar rule RE (regular expression)
int Regexp::ParseRE(const STRCHAR **re, BOOL &bEmpty)
{
	if (**re == '\0')
		return -1;

	int p = ParseAltE(re, bEmpty);
	if (p < 0) {
		if (m_Error)
			ThrowError();
		return p;
	}

	BOOL bEmptyRE = true;
	ParseRE(re, bEmptyRE);
	if (m_Error)
		ThrowError();
	bEmpty = bEmpty && bEmptyRE;
	return p;
}

void Regexp::FixupMatchPos(RegexpMatch* con, const STRCHAR *original_s, const STRCHAR *new_s)
{
	con->m_CompleteMatch.szStart = original_s + (con->m_CompleteMatch.szStart - new_s);
	con->m_CompleteMatch.szEnd = original_s + (con->m_CompleteMatch.szEnd - new_s);
	for (UINT i=0; i<con->m_uNumGroups; i++) {
		con->m_Inner->m_Matches[i].szStart = original_s + (con->m_Inner->m_Matches[i].szStart - new_s);
		con->m_Inner->m_Matches[i].szEnd = original_s + (con->m_Inner->m_Matches[i].szEnd - new_s);
	}
}

#endif // STR_USE_REGEX


/*********************************************************************
* Classes:	StrArray
* Purpose:	Portable equivalent of MFC CStringArray
*********************************************************************/

#ifdef STR_USE_EXTRAS

typedef veclite<Str> StrArray_Vector;

#ifdef STR_WIN32
#define VECTOR_REF_AT(x)	at((x))
#else
#define VECTOR_REF_AT(x)	operator[]((x))
#endif

inline StrArray_Vector* SaData(const StrArray* obj)
{
	return (StrArray_Vector*) obj->m_Data;
}

StrArray::StrArray()
{
	m_Data = NULL;
	m_GrowBy = -1;
	StrArray_Vector* vec = new StrArray_Vector();
	m_Data = vec;
}

StrArray::~StrArray()
{
	if (m_Data) {
		StrArray_Vector* vec = SaData(this);
		vec->clear();
		delete SaData(this);
	}
}

int StrArray::GetSize() const
{
	StrArray_Vector* vec = SaData(this);
	return static_cast<int> (vec->size());
}

int StrArray::GetCount() const
{
	StrArray_Vector* vec = SaData(this);
	return static_cast<int> (vec->size());
}

BOOL StrArray::IsEmpty() const
{
	StrArray_Vector* vec = SaData(this);
	return vec->empty() ? TRUE : FALSE;
}

void StrArray::RemoveAll()
{
	StrArray_Vector* vec = SaData(this);
	// Bugfix 15 Jan 03: return statement with an expression
	//   evaluating to void not supported in vc6
	vec->clear();
}

void StrArray::SetSize(int newsize, int growby /*= -1*/)
{
	StrArray_Vector* vec = SaData(this);
	vec->resize(newsize);
	m_GrowBy = growby;
}

void StrArray::FreeExtra()
{
	StrArray_Vector* vec = SaData(this);
	// The following may actually be a no-op with certain STL implementations
	vec->reserve(vec->size());
}

const Str& StrArray::GetAt(int index) const
{
	StrArray_Vector* vec = SaData(this);
	return vec->operator [](index);
}

Str& StrArray::ElementAt(int index)
{
	StrArray_Vector* vec = SaData(this);
	return vec->VECTOR_REF_AT(index);
}

void StrArray::SetAt(int index, const STRCHAR* elem)
{
	StrArray_Vector* vec = SaData(this);
	vec->VECTOR_REF_AT(index) = elem;
}

void StrArray::SetAt(int index, const Str& elem)
{
	StrArray_Vector* vec = SaData(this);
	vec->VECTOR_REF_AT(index) = elem;
}

void StrArray::MaybeGrow (int toMaxIndex)
{
	StrArray_Vector* vec = SaData(this);
	size_t vec_size = vec->size();
	size_t toMaxCount = static_cast<size_t> (toMaxIndex) + 1;
	if (vec_size < toMaxCount) {
		size_t vec_capacity = vec->capacity();
		if (vec_capacity < toMaxCount) {
			if (m_GrowBy < 0)
				vec->reserve(toMaxCount);
			else
				vec->reserve(STR_max(vec_capacity + m_GrowBy, toMaxCount));
		}
		vec->resize(toMaxCount);
	}
}

int StrArray::Append(const StrArray& src)
{
	StrArray_Vector* vec = SaData(this);
	int src_count = src.GetCount();
	int append_at = GrowRel(src_count);
	StrArray_Vector* src_vec = SaData(&src);
	for (int i=0; i<src_count; i++)
		vec->VECTOR_REF_AT(append_at + i) = src_vec->operator [](i);
	return append_at;
}

#if defined(__AFX_H__) && !defined(STR_NO_WINSTUFF)
int StrArray::Append(const CStringArray& src)
{
	StrArray_Vector* vec = SaData(this);
	int src_count = static_cast<int> (src.GetSize());
	int append_at = GrowRel(src_count);
	for (int i=0; i<src_count; i++)
		vec->VECTOR_REF_AT(append_at + i) = (const STRCHAR*) src[i];
	return append_at;
}
#endif

void StrArray::CopyElements(int destStartIndex, int srcStartIndex, int count, const StrArray& srcArray)
{
	StrArray_Vector* vec = SaData(this);
	StrArray_Vector* src_vec = SaData(&srcArray);
	for (int i=0; i<count; i++)
		vec->VECTOR_REF_AT(destStartIndex+i) = src_vec->operator [](srcStartIndex+i);
}

void StrArray::InsertAt(int startIndex, const StrArray& newArray)
{
	StrArray_Vector* vec = SaData(this);
	int newArrayCount = newArray.GetCount();
	{
		Str emptyEl;
		vec->insert(startIndex, newArrayCount, emptyEl);
	}
	this->CopyElements(startIndex, 0, newArrayCount, newArray);
}

void StrArray::Copy(const StrArray& src)
{
	StrArray_Vector* vec = SaData(this);
	StrArray_Vector* src_vec = SaData(&src);
	size_t newsize = src_vec->size();
	if (vec->size() != newsize)
		vec->resize(newsize);
	for (size_t i=0; i<newsize; i++)
		vec->VECTOR_REF_AT(i) = src_vec->operator [](i);
}

void StrArray::InsertAt(int index, const Str& newElement, int nCount /*= 1*/)
{
	StrArray_Vector* vec = SaData(this);
	vec->insert(index, nCount, newElement);
}

void StrArray::RemoveAt(int index, int nCount /*= 1*/)
{
	StrArray_Vector* vec = SaData(this);
	vec->erase_range(index, nCount);
}

#if defined(STR_WIN32) && !defined(STR_NO_WINSTUFF) && !defined(STR_NO_UNICODE)

SAFEARRAY* StrArray::ToSafeArray(VARTYPE destVtType, int startIndex /*= 0*/, int nCount /*= -1*/)
{
	if (destVtType != VT_VARIANT  &&  destVtType != VT_BSTR) {
		Str::ErrorNoObject(Str::SE_ARR_BadFormat);
		return NULL;	// Keep compiler happy
	}
	if (GetSize() == 0) {
		// Special case: return empty array
		if (startIndex != 0  ||  (nCount != -1  &&  nCount != 0)) {
			Str::ErrorNoObject(Str::SE_ARR_BadElementIndex);
			return NULL;
		}
		SAFEARRAYBOUND bounds0[1];
		bounds0[0].lLbound = 0;
		bounds0[0].cElements = 0;
		SAFEARRAY* ret0 = ::SafeArrayCreate(destVtType, 1, bounds0);
		if (ret0 == NULL) {
			Str::ErrorNoObject(Str::SE_OutOfMemory);
			return NULL;
		}
		return ret0;
	}
	if (startIndex < 0  ||  startIndex >= GetSize()) {
		Str::ErrorNoObject(Str::SE_ARR_BadElementIndex);
		return NULL;
	}
	if (nCount == -1)
		nCount = GetSize() - startIndex;
	else {
		if (nCount > (GetSize() - startIndex)) {
			Str::ErrorNoObject(Str::SE_ARR_BadElementIndex);
			return NULL;
		}
	}
	SAFEARRAYBOUND bounds[1];
	bounds[0].lLbound = 0;
	bounds[0].cElements = nCount;
	SAFEARRAY* ret = ::SafeArrayCreate(destVtType, 1, bounds);
	if (ret == NULL) {
		Str::ErrorNoObject(Str::SE_OutOfMemory);
		return NULL;
	}
	void* pvData = NULL;
	try {
		if (::SafeArrayAccessData(ret, &pvData) != 0) {
			Str::ErrorNoObject(Str::SE_ARR_COMError);
			return NULL;
		}
		if (destVtType == VT_BSTR) {
			BSTR* vtData = (BSTR*) pvData;
			for (int i=0; i<nCount; i++)
				vtData[i] = this->GetAt(startIndex+i).AllocSysString();
		}
		else {
			VARIANT* vtData = (VARIANT*) pvData;
			for (int i=0; i<nCount; i++) {
				vtData[i].vt = VT_BSTR;
				vtData[i].bstrVal = this->GetAt(startIndex+i).AllocSysString();
			}
		}
		::SafeArrayUnaccessData(ret);
		return ret;
	}
	catch (...) {
		if (pvData)
			::SafeArrayUnaccessData(ret);
		::SafeArrayDestroy(ret);
		throw;
	}
}

#endif

#endif	// STR_USE_EXTRAS

#ifdef STR_BORCPPBUILDER
#undef Char
#endif
