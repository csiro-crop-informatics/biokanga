#pragma once

const char cDfltFieldSepChar = ',';	// default character used to separate fields in CSV file
const char cDfltCommentChar = '#';	// default character used as a comment line indicator in CSV file
const char cDfltTextChar = '"';		// default character used to quote text in CSV file


const int cCSVMaxFields   = 2048;	// max number of fields
const int cCSVDfltFields  = 100;	// default number of fields
const int cCSVMinFields   = 10;		// can always handle at least 10 fields


const int cCSVMaxFieldLen = 0x03fff;	// maximum length (chars) of any single field supported
const int cCSVDfltFieldLen = 200;		// default length of any single field supported
const int cCSVMinFieldLen  = 10;		// can always handle fields of at least this size

// Note: the '10' in the following two consts is to allow for comma separators plus some whitespace
const int cCSVMaxLineLen  = (cCSVMaxFieldLen * 4);	// max number of characters which can be buffered in any line
const int cCSVDfltLineLenx = cCSVDfltFields * (cCSVDfltFieldLen + 10);// default number of characters which can be buffered
const int cCSVMinLineLenx  = cCSVMinFields * (cCSVMinFieldLen + 10);// minimum number of characters which can be buffered

const int cMaxAllocInBuff = 0x0ffffff; // allocate upto this sized input buffer

#pragma pack(1)
typedef struct TAG_sCSVField {
	int FieldID;			// which field (1..n)
	int CurFieldLen;		// number of chars currently ptd at by pValue
	bool bQuoted;			// true if field was a quoted value
	bool bQuotedErrs;		// quoted field had additional non-whitespace characters before next field separator
	bool bQuotedIncomplete; // quoted field was terminated by '\0' before closing quote parsed
	bool bTruncated;		// true if field value was truncated to m_MaxFieldLen
	char *pValue;			// field value
	} tsCSVField;



#pragma pack()


class CCSVFile  : public CErrorCodes
{
	int m_hFile;			// opened file
	INT64 m_StatFileSize;		// file size as returned by stat() when file initially opened
	char m_szFileName[_MAX_PATH+1];	// name of file
	bool m_bEOF;			// true when CSV file read upto EOF

	UINT32 m_EstNumRows;	// file is estimated to hold this many rows
	INT32 m_EstMaxFields;	// estimated to hold this max number of fields for any row
	INT32 m_EstMaxChrsRow;	// longest file row estimated to be of this length

	char m_FieldSepChar;	// character used to as field separator
	char m_TextQuote;		// character used to quote text strings
	char m_CommentChar;		// character to treat as starting a comment line

	int m_MaxFields;		// maximum expected fields per line
	int m_MaxFieldLen;		// maximum expected length of any parsed field
	int m_MaxLineLen;		// maximum expected line length

	int m_NumBuffered;		// total number of chars currently buffered in m_pBuffer
	int m_MaxBuffered;		// maximum number that can be buffered
	int m_BuffParseOfs;		// offset within m_pBuffer at which parsing will next start
	char *m_pBuffer;		// pts to buffered chars read from m_hFile


	int m_CurLineNum;		// current line number
	int m_CurLineLen;		// current line length
	int m_CurNumFields;		// number of fields parsed from current line
	int m_CurMaxFieldLen;	// maximum field length of any field in current line
	tsCSVField *m_pFields;	// pts to allocated array sized to hold m_MaxFields fields

	char m_szLastProcLine[cCSVMaxLineLen+1];	// copy of raw current line as just processed by NextLine
	int m_LastProcLineLen;		// strlen(m_szLastProcLine)

	void Reset(bool bInit = false);
	int ParseField(char *pszField,int *pChrsParsed,tsCSVField *pField); // parse from pszLine into pField, sets *pChrsParsed to number of chars parsed

public:
	CCSVFile(void);
	~CCSVFile(void);

	UINT32									// returns estimated number of rows, 0 if unable to estimate
		CSVEstSizes(char *pszFile,	// CSV file path+name to estimate sizes
			  INT64 *pFileSize=NULL,			// file is this size on disk
			  INT32 *pMaxNumFields=NULL,		// with this max number of fields in any row
			  INT32 *pMeanNumFields=NULL,	// mean number of fields (should be same as MaxFields)
			  INT32 *pMaxNumChrsRow=NULL,	// and this maximal number of chars in any row
			  INT32 *pMeanNumChrsRow=NULL);	// mean number of chars per row

	UINT32	EstNumRows(void);				// returns an estimate of the number of rows in currently opened CSV file
	
	int Open(char *pszFileName);	// file to open
	int Close(void);

	int SetMaxFields(int NumFields);	// set maximum number of fields - Note: must be set prior to Open() or default used
	int SetMaxFieldLen(int FieldLen);	// set maximum length of any field  - Note: must be set prior to Open() or default used
	int SetMaxLineLen(int MaxLineLen);	// set maximum length of any line - Note: must be set prior to Open() or default used
	int GetMaxFields(void);				// get maximum number of fields
	int GetMaxFieldLen(void);			// get maximum length of any field
	int GetMaxLineLen(void);			// get maximum length of any line

	int GetCurFields(void);				// get current number of fields
	int GetCurFieldLen(void);			// get current maximum length of any field
	int GetCurLineLen(void);			// get current line length

	int NextLine(void);					// move to next line in CSV file, skips blank lines and comment lines starting with comment char
	int GetLineNumber(void);			// get current line number
	bool IsLikelyHeaderLine(void);		// returns true if current line is heuristically likely to be a header line
	int GetQuoted(int FieldID);			  // 0 if unquoted, 1 if quoted, < 0 if errors
	int GetQuotedErrs(int FieldID);		  // 0 if noerrs, < 0 if errors
	int GetQuotedIncomplete(int FieldID); // 0 if noerrs, < 0 if errors

	int GetInt(int FieldID,int *pRetInt);	// parses (atoi) and returns specified field as int
	int GetLong(int FieldID,long *pRetLong);		// parses (atol) and returns specified field as long
	int GetInt64(int FieldID,INT64 *pRetInt64);	// parses (atoi64) and returns specified field as INT64
	int GetDouble(int FieldID,double *pRetDouble); // parses (atof) and returns specified field as double
	int GetText(int FieldID,char **ppRetText);			// returns text from specified field in current line 
	int GetLine(int Max2Ret,char *pRetLine);	// returns a copy of raw line as processed by NextLine 
};
