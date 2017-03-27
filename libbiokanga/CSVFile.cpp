/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CCSVFile::CCSVFile(void)
{
m_hFile = -1;
m_pBuffer = NULL;
m_pFields = NULL;
Reset(true);
}

CCSVFile::~CCSVFile(void)
{
Reset();
}


// CSVEstSizes
// Assumes CSV rows are reasonably inter-row field populated
const int cEstChrsBuff = 500000;		// 1st 500K chars from file should contain sufficient rows to give reasonable estimate for total rows in file

UINT32									// returns estimated number of rows, 0 if unable to estimate
CCSVFile::CSVEstSizes(char *pszFile,	// CSV file path+name to estimate sizes
			  INT64 *pFileSize,			// file is this size on disk
			  INT32 *pMaxNumFields,		// with this max number of fields in any row
			  INT32 *pMeanNumFields,	// mean number of fields (should be same as MaxFields)
			  INT32 *pMaxNumChrsRow,	// and this maximal number of chars in any row
			  INT32 *pMeanNumChrsRow)	// mean number of chars per row
{
int hFile;
INT64 FileSize;
int NumInBuff;
int BuffSize;
UINT8 *pBuff;
char *pChr;
int NumChrsParsed;
char Chr;
UINT32 TotNumRows;
UINT32 EstTotNumRows;
int MaxNumFields;
int NumFieldsThisRow;
int TotNumFieldsAllRows;
int NumChrsThisRow;
int TotNumChrsAllRows;
int MaxChrsRow;
bool bInDoubleQuotes;
bool bInSingleQuotes;
int NumFieldsQuoted;
int StatRslt;

if(pMaxNumFields != NULL)
	*pMaxNumFields = 0;
if(pMeanNumFields != NULL)
	*pMeanNumFields = 0;
if(pMaxNumChrsRow != NULL)
	*pMaxNumChrsRow = 0;
if(pMeanNumChrsRow != NULL)
	*pMeanNumChrsRow = 0;
if(pFileSize != NULL)
	*pFileSize = 0;

#ifdef _WIN32
struct _stat64 st;
if((StatRslt=_stat64(pszFile,&st))==0)
#else
struct stat64 st;
if((StatRslt = stat64(pszFile,&st)) == 0)
#endif
	FileSize = (INT64)st.st_size;
else
	FileSize = 0;

if(pFileSize != NULL)
	*pFileSize = FileSize;

if(FileSize == 0)		// 0 if file not readable or if 0 length
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"CSVEstSizes: Unable to estimate sizes for file '%s', does file exist and not 0 length, or is it readable (error == %d)",pszFile, StatRslt);
	return(0);
	}

if(FileSize < 10)		// arbitrary minimum CSV file size...
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"CSVEstSizes: Unable to estimate sizes for file '%s', file exists but is only %d long",pszFile,FileSize);
	return(0);
	}

BuffSize = (int)min(FileSize+100,(INT64)cEstChrsBuff);
if((pBuff = new UINT8 [BuffSize])==NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"CSVEstSizes: Unable to allocate %d memory for file '%s'",BuffSize,pszFile);
	return(0);
	}
// now can try to actually open file
hFile = open(pszFile,O_READSEQ);
if(hFile == -1)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"CSVEstSizes: Unable to open file '%s'",pszFile);
	delete pBuff;
	return(0);
	}

NumInBuff = (int)read(hFile,pBuff,BuffSize - 1);
close(hFile);

if(NumInBuff < 10)			// what kind of CSV file would be this small :-)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"CSVEstSizes: Unable to estimate sizes for file '%s', file exists but is only %d long",pszFile,FileSize);
	return(0);
	}
pBuff[NumInBuff] = 0;

TotNumRows = 0;
MaxNumFields = 0;
NumFieldsThisRow = 0;
TotNumFieldsAllRows = 0;
TotNumChrsAllRows = 0;
NumChrsThisRow = 0;
MaxChrsRow = 0;
NumFieldsQuoted = 0;
bInDoubleQuotes = false;
bInSingleQuotes = false;
pChr = (char *)pBuff;
NumChrsParsed = 0;
while((Chr = *pChr++) != '\0')
	{
	NumChrsParsed += 1;
	NumChrsThisRow += 1;
	if(Chr == '\r' || Chr == '\n')				// end of current row?
		{
		while(*pChr == '\r' || *pChr == '\n')	// slough repeated CR/LF's
			{
			pChr++;
			NumChrsParsed += 1;
			}
		bInDoubleQuotes = false;				// if CR/LFs are within quotes then most parsers will fail!
		bInSingleQuotes = false;
		if(NumFieldsThisRow > 0)					// accept as row if contains at least one field
			{
			if(NumFieldsThisRow > MaxNumFields)
				MaxNumFields = NumFieldsThisRow;
			TotNumFieldsAllRows += NumFieldsThisRow;
			TotNumRows += 1;
			}
		if(NumChrsThisRow > MaxChrsRow)
			MaxChrsRow = NumChrsThisRow;
		TotNumChrsAllRows += NumChrsThisRow;
		NumFieldsThisRow = 0;
		NumChrsThisRow = 0;
		NumFieldsQuoted = 0;
		continue;
		}
	if(Chr <= ' ')							// not interested in any whitespace or control chrs
		continue;
	if(NumFieldsThisRow == 0)						// if any non-whitespace then there is at least one col
		NumFieldsThisRow = 1;
	
	if(Chr == ',' && !(bInDoubleQuotes || bInSingleQuotes))					// could be a col separator, but not if contained within quotes
		{
		NumFieldsThisRow += 1;
		continue;
		}

	// double quotes within single, and single within double processing
	if(Chr == '\'')								
		{
		if(bInDoubleQuotes)
			continue;
		if(!bInSingleQuotes)
			NumFieldsQuoted += 1;
		bInSingleQuotes = !bInSingleQuotes;
		continue;
		}

	if(Chr == '"')
		{
		if(bInSingleQuotes)
			continue;
		if(!bInDoubleQuotes)
			NumFieldsQuoted += 1;
		bInDoubleQuotes = !bInDoubleQuotes;
		continue;
		}
	}

delete pBuff;

if(TotNumRows)
	{
	if(pMeanNumChrsRow)
		*pMeanNumChrsRow = TotNumChrsAllRows/TotNumRows;
	if(pMeanNumFields)
		*pMeanNumFields = TotNumFieldsAllRows/TotNumRows;
	}

if(NumFieldsThisRow > 0)
	{
	if(NumFieldsThisRow > MaxNumFields)
		MaxNumFields = NumFieldsThisRow;
	TotNumRows += 1;
	}
if(NumChrsThisRow > MaxChrsRow)
	MaxChrsRow = NumChrsThisRow;
if(pMaxNumFields != NULL)
	*pMaxNumFields = MaxNumFields;
if(pMaxNumChrsRow != NULL)
	*pMaxNumChrsRow = MaxChrsRow;

if(NumChrsParsed == FileSize)
	EstTotNumRows = TotNumRows;
else
	EstTotNumRows = (UINT32)(((TotNumRows + 1) * FileSize)/(UINT64)NumChrsParsed);	// better to slightly over estimate than underestimate..
return(EstTotNumRows);
}

UINT32						// returns an estimate of the number of rows in currently opened CSV file
CCSVFile::EstNumRows(void)
{
return(m_EstNumRows);
}

int 
CCSVFile::Open(char *pszFileName)		// file to open
{
int Idx;
int Rslt;
Reset();

// estimate number of rows, max fields and max line length
m_EstNumRows = CSVEstSizes(pszFileName,&m_StatFileSize,&m_EstMaxFields,NULL,&m_EstMaxChrsRow);
if(m_EstNumRows == 0)
	{
	if(m_StatFileSize < 5)
		AddErrMsg("CCSVFile::Open","Unable to estimate sizes for '%s', does it exist and have appropriate read permissions or is empty?",pszFileName);
	else
		AddErrMsg("CCSVFile::Open","Unable to estimate sizes for '%s'",pszFileName);
	Rslt = eBSFerrOpnFile;
	return(Rslt);
	}

// now can try to actually open file
m_hFile = open(pszFileName,O_READSEQ);
if(m_hFile == -1)
	{
	AddErrMsg("CCSVFile::Open","Unable to open %s - %s",pszFileName,strerror(errno));
	Rslt = eBSFerrOpnFile;
	return(Rslt);
	}

m_MaxBuffered = (int)min(m_StatFileSize+1,(INT64)cMaxAllocInBuff);

strcpy(m_szFileName,pszFileName);
m_pBuffer = new char[m_MaxBuffered + 16];		// additional 16 as a safety margin
if(m_pBuffer == NULL)
	{
	AddErrMsg("CCSVFile::Open","Memory allocation of %d bytes for %s- %s",m_MaxBuffered,pszFileName,strerror(errno));
	Reset();			// closes opened file..
	return(eBSFerrMem);
	}

m_pFields = new tsCSVField [m_MaxFields];
if(m_pFields == NULL)
	{
	AddErrMsg("CCSVFile::Open","Memory allocation (m_MaxFields) of %d bytes for %s- %s",m_MaxFields * sizeof(tsCSVField),pszFileName,strerror(errno));
	Reset();			// closes opened file..
	return(eBSFerrMem);
	}
memset(m_pFields,0,m_MaxFields * sizeof(tsCSVField));
for(Idx=0;Idx < m_MaxFields; Idx++)
	{
	m_pFields[Idx].pValue = new char [m_MaxFieldLen];
	if(m_pFields[Idx].pValue == NULL)
		{
		AddErrMsg("CCSVFile::Open","Memory allocation (m_MaxFieldLen) of %d bytes for %s- %s",m_MaxFieldLen,pszFileName,strerror(errno));
		Reset();			// closes opened file..
		return(eBSFerrMem);
		}
	m_pFields[Idx].FieldID = Idx+1;
	}

return(eBSFSuccess);
}

int 
CCSVFile::Close(void)
{
Reset();
return(eBSFSuccess);
}

void
CCSVFile::Reset(bool bInit)
{
int Idx;
if(m_hFile != -1)
	{
	close(m_hFile);
	m_hFile = -1;
	}

m_szFileName[0]='\0';	// name of file
m_bEOF = false;

if(m_pBuffer!=NULL)	// line buffer
	{
	delete m_pBuffer;
	m_pBuffer = NULL;
	}
m_NumBuffered = 0;		// total number of chars currently buffered in m_pBuffer
m_MaxBuffered = 0;		// maximum number that can be buffered
m_BuffParseOfs = 0;		// offset within m_pBuffer at which parsing will next start
m_LastProcLineLen = 0;

m_EstNumRows = 0;	// file is estimated to hold this many rows
m_EstMaxFields = 0;	// estimated to hold this max number of fields for any row
m_EstMaxChrsRow = 0;	// longest row estimated to be of this length

m_CurLineNum = 0;		// current line number
m_CurNumFields = 0;		// number of fields parsed from current line
m_CurMaxFieldLen = 0;	// max len of any field in current line
if(m_pFields != NULL)	// pts to allocated array sized to hold m_MaxFields fields
	{
	for(Idx = 0; Idx < m_MaxFields; Idx++)
		delete m_pFields[Idx].pValue;
	delete m_pFields;
	m_pFields = NULL;
	}

if(bInit)
	{
	m_FieldSepChar = cDfltFieldSepChar;	// character used to as field separator
	m_TextQuote = cDfltTextChar;		// character used to quote text strings
	m_CommentChar = cDfltCommentChar;	// character to treat as starting a comment line
	m_MaxFields = cCSVDfltFields;		// maximum expected fields per line
	m_MaxFieldLen = cCSVDfltFieldLen;	// maximum expected length of any parsed field
	m_MaxLineLen = cCSVMaxLineLen;		// maximum expected number of characters in any line
	}
}

int 
CCSVFile::SetMaxFields(int NumFields)	// set maximum number of fields - Note: must be set prior to Open() or default used
{
if(NumFields > cCSVMaxFields)
   return(eBSFerrParams);
if(m_hFile != -1)
	return(eBSFerrFileOpened);
if(NumFields < cCSVMinFields)
	m_MaxFields = cCSVMinFields;
else
	m_MaxFields = NumFields;
return(eBSFSuccess);
}
int 
CCSVFile::SetMaxFieldLen(int FieldLen)	// set maximum length of any field  - Note: must be set prior to Open() or default used
{
if(FieldLen > cCSVMaxFieldLen)
   return(eBSFerrParams);
if(m_hFile != -1)
	return(eBSFerrFileOpened);
if(FieldLen < cCSVMinFieldLen)
	m_MaxFieldLen = cCSVMinFieldLen;
else
	m_MaxFieldLen = FieldLen;
m_MaxFieldLen = FieldLen;
return(eBSFSuccess);
}
int 
CCSVFile::SetMaxLineLen(int MaxLineLen)		// set maximum length of any line - Note: must be set prior to Open() or default used
{
if(MaxLineLen > cCSVMaxLineLen)
   return(eBSFerrParams);
if(m_hFile != -1)
	return(eBSFerrFileOpened);
if(MaxLineLen < cCSVMaxLineLen)
	m_MaxLineLen = cCSVMaxLineLen;
else
	m_MaxLineLen = MaxLineLen;
return(eBSFSuccess);
}

int 
CCSVFile::GetMaxFields(void)			// get maximum number of fields
{
return(m_MaxFields);
}

int 
CCSVFile::GetMaxFieldLen(void)			// get maximum length of any field
{
return(m_MaxFieldLen);
}

int 
CCSVFile::GetMaxLineLen(void)			// get maximum length of any line
{
return(m_MaxLineLen);
}

int 
CCSVFile::GetCurFields(void)			// get current number of fields
{
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
return(m_CurNumFields);
}

int 
CCSVFile::GetCurFieldLen(void)			// get current maximum length of any field
{
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
return(m_CurMaxFieldLen);
return(0);
}

int 
CCSVFile::GetCurLineLen(void)			// get current line length
{
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
return(m_CurLineLen);
}

// 0: field was parsed and terminated by field separator
// 1: field was parsed and terminated by 0x0a (NL) or 0x0d (CR)
// 2: field was parsed and terminated by 0x00 ('\00')
int 
CCSVFile::ParseField(char *pszField,		// field to parse 
					 int *pChrsParsed,		// returns number of chars actually parsed from pszLine
					 tsCSVField *pField)	// field to parse into
{
char szField[cCSVMaxFieldLen];				// field is initially parsed into this buffer, trimmed for whitespace, and then copied into field value pField->pValue
char Chr;
char *pValue = szField;
int ParseState = 0;
int FieldState = 2;

pField->CurFieldLen = 0;                                                                                                                                                                                              
pField->bQuoted = false;
pField->bQuotedErrs = false;
pField->bQuotedIncomplete = false;
pField->bTruncated = false;

pField->pValue[0] = '\0';
*pChrsParsed = 0;

while(FieldState == 2 && (Chr = *pszField++) != '\0')
	{
	*pChrsParsed += 1;
	if(pField->CurFieldLen >= cCSVMaxFieldLen)
		break;
	switch(ParseState) {
		case 0:					// haven't started to process field
			if(Chr == ' ' || Chr == '\t')	// slough leading whitespace
				continue;
			if(Chr == m_TextQuote)
				{
				pField->bQuoted = true;
				ParseState = 2;		// quoted text field starting
				}
			else
				{
				if(Chr == 0x0d || Chr == 0x0a)
					{
					if(Chr == 0x0a)
						m_CurLineNum++;
					FieldState = 1;
					continue;
					}

				if(Chr == m_FieldSepChar)	// empty field
					{
					FieldState = 0;
					continue;
					}
				
				*pValue++ = Chr;			// unquoted or numeric field
				pField->CurFieldLen=1;			
				ParseState = 1;		
				}
			continue;

		case 1:							// parsing non-quoted field
			if(Chr == m_FieldSepChar ||	// end of current field?
				Chr == 0x0d || Chr == 0x0a)
				{
				// trim any trailing whitespace
				while(pValue[-1] == ' ' || pValue[-1] == '\t')
					{
					pValue -= 1;
					pField->CurFieldLen--;
					}
				if(Chr == 0x0a)
					m_CurLineNum++;

				FieldState = (Chr == m_FieldSepChar) ? 0 : 1;
				continue;
				}
			*pValue++ = Chr;
			pField->CurFieldLen++;
			continue;

		case 2:						// parsing textual field
			if(Chr == m_TextQuote)	// end of current field?
				{
				if(*pszField == m_TextQuote) // double quotes allowed to escape a single quote
					{
					pszField++;				// slough second double quote
					*pChrsParsed += 1;
					}
				else				// single quote, must be end of quoted text
					{
					ParseState = 3;	// slough whitespace until field separator
					continue;
					}
				}
			*pValue++ = Chr;
			pField->CurFieldLen++;
			continue;

		case 3:						// finished quoted string, looking for field separator or line terminator
			if(Chr == m_FieldSepChar)
				{
				FieldState = 0;
				continue;
				}
			if(Chr == 0x0d || Chr == 0x0a)
				{
				if(Chr == 0x0a)
					m_CurLineNum++;
				FieldState = 1;
				continue;
				}
			if(Chr == ' ' || Chr == '\t') // simply slough whitespace
				continue;
			pField->bQuotedErrs = true;	// additional chars before field separator
			continue;
		}
	}

if(FieldState == 2)
	{
	switch(ParseState) {
		case 0:				// haven't started to process field
			break;

		case 1:			// processing non-quoted field
			while(pValue[-1] == ' ' || pValue[-1] == '\t')	// trim any trailing whitespace
				{
				pValue--;
				pField->CurFieldLen--;
				}
			break;

		case 2:			// parsing quoted textual field
			pField->bQuotedIncomplete = true;
			break;
			
		case 3:			// finished quoted string, looking for field separator or line terminator
			break;
		}
	}

if(pField->CurFieldLen >= m_MaxFieldLen)
	{
	pField->bTruncated = true;
	pField->CurFieldLen = m_MaxFieldLen-1;
	}
if(pField->CurFieldLen)
	memmove(pField->pValue,szField,pField->CurFieldLen);
pField->pValue[pField->CurFieldLen] = '\0';

return(FieldState);
}

int 
CCSVFile::NextLine(void)				// move to next line in CSV file, skips blank lines and comment lines starting with comment char
{
char *pChr;
int LenRead;
bool bInComment;
int Ofs;
int ChrsParsed;
int Idx;
int ParseRslt;
char Chr;
bool bMore = true;
int NumFields = 0;


if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);

m_CurNumFields = 0;
m_CurLineLen = 0;
m_CurMaxFieldLen = 0;
bInComment = false;

while(bMore) {
	// ensure that at all times there is sufficient buffered for at least one line of m_MaxLineLen chars
	if(!m_bEOF && (!m_NumBuffered || (m_NumBuffered - m_BuffParseOfs) < m_MaxLineLen))
		{
		if(m_BuffParseOfs > 0)
			{
			if(m_NumBuffered && ((m_NumBuffered - m_BuffParseOfs) > 0))
				{
				memmove(m_pBuffer,&m_pBuffer[m_BuffParseOfs],m_NumBuffered - m_BuffParseOfs);	// replaced memcpy with memmove because of gcc overlap problem
				}
			m_NumBuffered -= m_BuffParseOfs;
			m_pBuffer[m_NumBuffered] = '\0';
			m_BuffParseOfs = 0;
			}
		LenRead = read(m_hFile,&m_pBuffer[m_NumBuffered],m_MaxBuffered - m_NumBuffered - 1); // -1 to allow for a terminating '\0'
		if(LenRead <= 0)
			m_bEOF = true;
		else
			m_NumBuffered += LenRead;
		m_pBuffer[m_NumBuffered] = '\0';
		}
	
	if(m_bEOF && (m_NumBuffered == 0 || m_BuffParseOfs == m_NumBuffered))
		{
		bMore = false;
		continue;
		}

	// slough blank lines and those which are comments
	pChr = &m_pBuffer[m_BuffParseOfs];
	if(*pChr < ' ' &&  !(*pChr == '\t' || *pChr == '\n' || *pChr == '\r'))
		{
		bMore = false;
		continue;
		}

	if(*pChr == 0x0d || *pChr == 0x0a)
		{
		m_BuffParseOfs++;
		if(*pChr == 0x0d && pChr[1] == 0x0a)
			{
			m_CurLineNum++;
			m_BuffParseOfs++;
			}
		else
			if(*pChr == 0x0a)
				m_CurLineNum++;
		bInComment = false;
		continue;
		}
	if(bInComment || *pChr == ' ' || *pChr == '\t')
		{
		m_BuffParseOfs++;
		continue;
		}
	if(*pChr == m_CommentChar)
		{
		bInComment = true;
		m_BuffParseOfs++;
		continue;
		}

	// make a copy of the current line
	m_LastProcLineLen = 0;
	pChr = &m_pBuffer[m_BuffParseOfs];
	while(Chr = *pChr++)
		{
		m_szLastProcLine[m_LastProcLineLen++] = Chr;
		if((m_LastProcLineLen+1) == cCSVMaxLineLen || (*pChr < ' ' && *pChr != '\t'))
			break;
		}
	m_szLastProcLine[m_LastProcLineLen] = '\0';

	pChr = &m_pBuffer[m_BuffParseOfs];
	Ofs = 0;
	ChrsParsed = 0;
	for(Idx = 0; Idx < m_MaxFields; Idx++)
		{
		ParseRslt = ParseField(&pChr[Ofs],  // field to parse 
					 &ChrsParsed,			// returns number of chars actually parsed from pszLine
					 &m_pFields[Idx]);  // field to parse into
		Ofs += ChrsParsed;
		NumFields++;
		if(m_pFields[Idx].CurFieldLen > m_CurMaxFieldLen)
			m_CurMaxFieldLen = m_pFields[Idx].CurFieldLen;
		if(ParseRslt != 0)
			break;
		}

	if(ParseRslt == 0)	//0: field was parsed and terminated by field separator, means not EOL or EOT
		{
		// slough everything upto EOL or EOT
		pChr = &m_pBuffer[m_BuffParseOfs + Ofs];
		while(Chr = *pChr++)
			{
			Ofs += 1;
			if(*pChr == 0x0d || *pChr == 0x0a)
				{
				if(Chr == 0x0a)
					m_CurLineNum++;
				break;
				}
			}
		}

	m_BuffParseOfs += Ofs;
	m_CurNumFields = NumFields;
	m_CurLineLen = Ofs;
	return(NumFields);
	}
return(eBSFSuccess);
}

// returns a copy of raw line as was processed by NextLine 
int 
CCSVFile::GetLine(int Max2Ret,			// max to return, includes trailing '\0'
				  char *pRetLine)
{
if(Max2Ret < 1 || pRetLine == NULL)
	return(0);
Max2Ret -= 1;	
if(Max2Ret > m_LastProcLineLen)
	Max2Ret = m_LastProcLineLen;
if(Max2Ret > 0)
	memmove(pRetLine,m_szLastProcLine,Max2Ret);
pRetLine[Max2Ret] = '\0';
return(Max2Ret);
}

// IsLikelyHeaderLine
// Heuristically determines if the currently parsed CSV record is a likely header line
// If line contains at most 2 empty fields with the remainder all quoted values or non-numerics then assume likely to be a header line
// Caller makes the final decision as to if the current line really is a header line!
bool
CCSVFile::IsLikelyHeaderLine(void)
{
int Idx;
int NumEmpty;
tsCSVField *pField;
char *pTerm;
double Val;
pField = m_pFields;
if(m_hFile == -1 || !m_CurNumFields)
	return(false);
NumEmpty = 0;
for(Idx = 0; Idx < m_CurNumFields; Idx++,pField++)
	{
	if(pField->bQuoted)
		continue;
	if(!pField->CurFieldLen)
		{
		if(++NumEmpty > 2)
			return(false);
		continue;
		}
	pTerm = NULL;
	Val = strtod(pField->pValue,&pTerm);
	if(pTerm != NULL && *pTerm == '\0')
		return(false);
	}
return(true);
}

int 
CCSVFile::GetInt(int FieldID,int *pRetInt) 
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
*pRetInt = atoi(m_pFields[FieldID-1].pValue);
return(errno == ERANGE ? eBSFerrNumRange : eBSFSuccess);
}

int 
CCSVFile::GetLong(int FieldID,long *pRetLong) 
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
*pRetLong = atol(m_pFields[FieldID-1].pValue);
return(errno == ERANGE ? eBSFerrNumRange : eBSFSuccess);
}

int 
CCSVFile::GetInt64(int FieldID,INT64 *pRetInt64) 
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
#ifdef _WIN32
*pRetInt64 = _atoi64(m_pFields[FieldID-1].pValue);
#else
*pRetInt64 = atoll(m_pFields[FieldID-1].pValue);
#endif
return(errno == ERANGE ? eBSFerrNumRange : eBSFSuccess);
}

int 
CCSVFile::GetDouble(int FieldID,double *pRetDouble) 
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
*pRetDouble = atof(m_pFields[FieldID-1].pValue);
return(errno == ERANGE ? eBSFerrNumRange : eBSFSuccess);
}

int
CCSVFile::GetText(int FieldID,char **ppszRetText)			
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
*ppszRetText = m_pFields[FieldID-1].pValue;
return(eBSFSuccess);
}


int 
CCSVFile::GetLineNumber(void)			// get current line number
{
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
return(m_CurLineNum);
}

int 
CCSVFile::GetQuoted(int FieldID)	
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
return(m_pFields[FieldID-1].bQuoted ? 1 : eBSFSuccess);
}

int 
CCSVFile::GetQuotedErrs(int FieldID)		// quoted field had additional non-whitespace characters before next field separator
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
return(m_pFields[FieldID-1].bQuotedErrs ? eBSFerrQuoteErrs : eBSFSuccess);
}

int 
CCSVFile::GetQuotedIncomplete(int FieldID) // quoted field was terminated by '\0' before closing quote parsed
{
if(FieldID < 1 || FieldID > m_CurNumFields)
	return(eBSFerrFieldID);
if(m_hFile == -1)						// file has to be opened!
	return(eBSFerrFileClosed);
return(m_pFields[FieldID-1].bQuotedIncomplete ? eBSFerrQuoteIncomplete : eBSFSuccess);
}

