#pragma once
#include "./commdefs.h"

// diagnostics output levels
typedef enum {
	eDLNone = 0,			// do not write any diagnostics to file/screen
	eDLFatal,				// fatal diagnostic - always output to file/screen
	eDLWarn,				// warning + fatal write to file/screen
	eDLInfo,				// info + warning + fatal write to file/screen
	eDLDiag,				// diag + info + warning + fatal write to file/screen
	eDLDebug				// all diagnostics write to file/screen
} etDiagLevel;

const etDiagLevel cDfltScreenDiagLevel = eDLInfo; // default to screen diagnostics level
const etDiagLevel cDfltFileDiagLevel = eDLInfo;	  // default to file diagnostics level

const int cMaxDiagLen = cMaxReadLen+200;			// max line length (includes '\0' terminator) of any diagnostics message line
class CDiagnostics
{
	bool m_Opened;						// true if diagnostics output to be processed - Open() sets true, Close() sets false
	etDiagLevel m_FileDiagLevel;		// cutoff level for output to file diagnostics, levels above this are sloughed
	etDiagLevel m_ScreenDiagLevel;		// cutoff level for output to screen diagnostics, levels above this are sloughed
	int m_hFile;


public:
	CDiagnostics(void);
	CDiagnostics(const char *pszFile,etDiagLevel ScreenDiagLevel,etDiagLevel FileDiagLevel);
	CDiagnostics(const char *pszFile,etDiagLevel ScreenDiagLevel,etDiagLevel FileDiagLevel,bool bAppend);
	~CDiagnostics(void);

	bool Open(const char *pszFile,			// file to write diagnostics to - can be NULL if diagnostics to be written to file
		   etDiagLevel ScreenDiagLevel=cDfltScreenDiagLevel, // cutoff level for diagnostics to screen - if not eDLNone..eDLDebug then existing level is used 
		   etDiagLevel FileDiagLevel=cDfltFileDiagLevel,   // cutoff level for diagnostics to file - if not eDLNone..eDLDebug then existing level is used 
		   bool bAppend=false);			// true if file to be opened in append mode, default is to truncate existing log file 

	void Close(void);					// closes diagnostics output file and stops screen writes


	bool			// true if diagnostics was actually output
	DiagOut(etDiagLevel DiagLevel,					// diagnostics level
				    const char *pszSource,							// identifies diagnostics message source (can be NULL if source not to be output)
					const char *pszFormat,...);						// message format

	bool			// true if diagnostics was actually output
	DiagOutMsgOnly(etDiagLevel DiagLevel,				// diagnostics level
					const char *pszFormat,...);					// message format

	etDiagLevel SetScreenDiagLevel(etDiagLevel DiagLevel); // cutoff level for diagnostics to screen - if not eDLNone..eDLDebug then existing level is used 
	etDiagLevel SetFileDiagLevel(etDiagLevel DiagLevel); // cutoff level for diagnostics to file - if not eDLNone..eDLDebug then existing level is used 
	etDiagLevel GetScreenDiagLevel(void);
	etDiagLevel GetFileDiagLevel(void);

};
