/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CDiagnostics::CDiagnostics(void)
{
m_Opened = false;
m_ScreenDiagLevel = eDLNone;
m_FileDiagLevel = eDLNone;
m_hFile = -1;
}

CDiagnostics::CDiagnostics(const char *pszFile,etDiagLevel ScreenDiagLevel,etDiagLevel FileDiagLevel)
{
m_hFile = -1;
m_Opened = false;
if(pszFile != NULL && pszFile[0] != '\0')
	Open(pszFile,ScreenDiagLevel,FileDiagLevel,false);
else
	{
	m_ScreenDiagLevel=ScreenDiagLevel;
	m_FileDiagLevel=FileDiagLevel;
	}
}

CDiagnostics::CDiagnostics(const char *pszFile,etDiagLevel ScreenDiagLevel,etDiagLevel FileDiagLevel,bool bAppend)
{
m_Opened = false;
m_hFile = -1;
if(pszFile != NULL && pszFile[0] != '\0')
	Open(pszFile,ScreenDiagLevel,FileDiagLevel,bAppend);
else
	{
	m_ScreenDiagLevel=ScreenDiagLevel;
	m_FileDiagLevel=FileDiagLevel;
	}
}

CDiagnostics::~CDiagnostics(void)
{
if(m_hFile != -1)
	close(m_hFile);
}

bool
CDiagnostics::Open(const char *pszFile,			// file to write diagnostics to - can be NULL if diagnostics not to be written to file
				   etDiagLevel ScreenDiagLevel, // cutoff level for diagnostics to screen - if not eDLNone..eDLDebug then existing level is used 
				   etDiagLevel FileDiagLevel,   // cutoff level for diagnostics to file - if not eDLNone..eDLDebug then existing level is used 
				   bool bAppend)			// true if file to be opened in append mode, default is to truncate existing log file 
{
if(m_hFile != -1)
	{
	close(m_hFile);
	m_hFile = -1;
	}
m_Opened = false;
SetScreenDiagLevel(ScreenDiagLevel);
SetFileDiagLevel(FileDiagLevel);
if(pszFile != NULL && pszFile[0] != '\0')
	{
#ifdef _WIN32
        m_hFile = open(pszFile,(bAppend ? _O_APPEND : _O_TRUNC) | _O_CREAT | _O_TEXT | _O_RDWR,_S_IREAD | _S_IWRITE);
#else
    if(bAppend)
        {
        if(( m_hFile = open64(pszFile,O_RDWR | O_APPEND))==-1 && errno == 2)
          m_hFile = open64(pszFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
        }
    else
          if((m_hFile = open64(pszFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			{
			if(ftruncate(m_hFile,0)!=0)
				{
		        printf("\nUnable to truncate log file '%s' - error: '%s'\n",pszFile,strerror(errno));
		        fflush(stdout);
				return(false);
				}
			}
#endif
    if(m_hFile == -1)
		{
        printf("\nUnable to open log file '%s' - error: '%s'\n",pszFile,strerror(errno));
        fflush(stdout);
		return(false);
        }
	}
m_Opened = true;
return(true);
}

void
CDiagnostics::Close(void)
{
if(m_hFile != -1)
	{
	close(m_hFile);
	m_hFile = -1;
	}
m_Opened = false;
}

etDiagLevel		// returns level as set 
CDiagnostics::SetScreenDiagLevel(etDiagLevel DiagLevel)	// cutoff level for diagnostics to screen - if not eDLNone..eDLDebug then existing level is used 
{
if(DiagLevel >= eDLNone && DiagLevel <= eDLDebug)
	m_ScreenDiagLevel=DiagLevel;
return(m_ScreenDiagLevel);
}

etDiagLevel		// returns level as set 
CDiagnostics::SetFileDiagLevel(etDiagLevel DiagLevel)	// cutoff level for diagnostics to file - if not eDLNone..eDLDebug then existing level is used 
{
if(DiagLevel >= eDLNone && DiagLevel <= eDLDebug)
	m_FileDiagLevel = DiagLevel;
return(m_FileDiagLevel);
}

// returns current diagnostics level
etDiagLevel 
CDiagnostics::GetScreenDiagLevel(void)
{
return(m_ScreenDiagLevel);
}

// returns current diagnostics level
etDiagLevel 
CDiagnostics::GetFileDiagLevel(void)
{
return(m_FileDiagLevel);
}

// DiagOut
// Writes diagnostic message to screen/file
bool					// returns false if not Open()'d
CDiagnostics::DiagOut(etDiagLevel DiagLevel,					// diagnostics level
				    const char *pszSource,							// identifies diagnostics message source (can be NULL)
					const char *pszFormat,...)						// message format
{
va_list Args;
int LineLen;
char szDiag[cMaxDiagLen];
char szLine[cMaxDiagLen+128];			// need to allow for timestamp and message source
char szTimestamp[200];
#ifdef _WIN32
struct _timeb timebuffer;
#else
struct timeb timebuffer;
#endif
char *timeline;

if(!m_Opened)
	return(false);
if((DiagLevel <= eDLNone && DiagLevel > eDLDebug) ||
   (m_ScreenDiagLevel == eDLNone && (m_hFile == -1 || m_FileDiagLevel == eDLNone)))
	return(true);
#ifdef _WIN32
_ftime(&timebuffer);
#else
ftime(&timebuffer);
#endif
timeline = ctime(&timebuffer.time);

va_start(Args, pszFormat );
#ifdef _WIN32
_vsnprintf(szDiag,cMaxDiagLen,pszFormat,Args);
#else
vsnprintf(szDiag,cMaxDiagLen,pszFormat,Args);
#endif
szDiag[cMaxDiagLen-1] = '\0';

LineLen = sprintf(szTimestamp,"\n[%.15s.%03d %.4s]",&timeline[4],(int)timebuffer.millitm, &timeline[20]);
if(pszSource != NULL)
	sprintf(&szTimestamp[LineLen],"(%s) ",pszSource);
else
	sprintf(&szTimestamp[LineLen],"   ");

LineLen=sprintf(szLine,"%s%s",szTimestamp,szDiag);
if(m_hFile != -1 && DiagLevel <= m_FileDiagLevel)
	{
#ifdef _WIN32
	write(m_hFile,szLine,LineLen);
	_commit(m_hFile);
#else
	if(write(m_hFile,szLine,LineLen));
	fsync(m_hFile);
#endif
	}
if(DiagLevel <= m_ScreenDiagLevel)
	{
	fputs(szLine,stdout);
	fflush(stdout);
	}
return(true);
}

bool			// true if diagnostics was actually output
CDiagnostics::DiagOutMsgOnly(etDiagLevel DiagLevel,				// diagnostics level
					const char *pszFormat,...)					// message format
{
va_list Args;
int LineLen;
char szDiag[cMaxDiagLen];
char szLine[cMaxDiagLen+2];			// need to allow for '\n'

if(!m_Opened)
	return(false);
if((DiagLevel <= eDLNone && DiagLevel > eDLDebug) ||
   (m_ScreenDiagLevel == eDLNone && (m_hFile == -1 || m_FileDiagLevel == eDLNone)))
	return(true);

va_start(Args, pszFormat );
#ifdef _WIN32
_vsnprintf(szDiag,cMaxDiagLen,pszFormat,Args);
#else
vsnprintf(szDiag,cMaxDiagLen,pszFormat,Args);
#endif

szDiag[cMaxDiagLen-1] = '\0';

LineLen=sprintf(szLine,"\n        %s",szDiag);
if(m_hFile != -1 && DiagLevel <= m_FileDiagLevel)
	{
#ifdef _WIN32
	write(m_hFile,szLine,LineLen);
	_commit(m_hFile);
#else
	if(write(m_hFile,szLine,LineLen));
	fsync(m_hFile);
#endif
	}
if(DiagLevel <= m_ScreenDiagLevel)
	{
	printf("%s",szLine);
	fflush(stdout);
	}
return(true);
}

