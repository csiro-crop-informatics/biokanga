// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License
//   Please contact stuart.stephen@csiro.au for support or 
//   to submit modifications to this source

#include "stdafx.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

#include "SQLiteSummaries.h"

// Following database schema is utilised
// Tables

//  TblExprs
//	    ExprID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this experiment
//      ExprName VARCHAR(50) UNIQUE,	-- Experiment name
//      ExprTitle VARCHAR(50),			-- Experiment title
//      ExprDescr VARCHAR(1000)			-- describes Experiment

//	TblProcess
//		ProcessID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this process
//      ProcessName VARCHAR(50),			-- process name
//      ProcessTitle VARCHAR(50),			-- process title
//      ProcessDescr VARCHAR(1000)			-- describes process

//	TblProcessing
//		ProcessingID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this processing instance
//	    ExprID INTEGER,							-- processing in this experiment
//		ProcessID INTEGER,					    -- using this process
//      ProcessVer VARCHAR(20)					-- process version
//		Start TIMESTAMP,					    -- processing started
//      Finish TIMESTAMP,					    -- processing completed
//      ResultCode INTEGER					    -- processing end result code

//	TblProcessingLog
//		ProcessingLogID INTEGER PRIMARY KEY ASC, -- Uniquely identifies this processing log instance
//		ProcessingID INTEGER,					-- identifies this processing instance
//      LogText VARCHAR(2000)					-- log text

// TblParams
//		ParamID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this parameter instance
//		ProcessingID INTEGER,				-- identifies this processing instance
//      ParamType INTEGER,					-- parameter value type - see teSQLliteSummParamTypes
//		ParmName VARCHAR(50),				-- parameter name
//      ParmValue VARCHAR(2000),			-- parameter value

// TblResults
//		ResultID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this result instance
//		ProcessingID INTEGER,				-- identifies this processing instance
//      GroupAs VARCHAR(50),				-- result is part of this grouping
//      ResultType INTEGER,					-- result value type - see teSQLliteSummParamTypes
//		ResultName VARCHAR(50),				-- result name
//      ResultValue VARCHAR(500),			-- result value

// TblXYResults
//		ResultID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this result instance
//		ProcessingID INTEGER,				-- identifies this processing instance
//      GroupAs VARCHAR(50),				-- result is part of this grouping
//      ResultXType INTEGER,				-- result X value type - see teSQLliteSummParamTypes
//		ResultXName VARCHAR(50),			-- result X name
//      ResultXValue VARCHAR(500),			-- result X value
//      ResultYType INTEGER,				-- result Y value type - see teSQLliteSummParamTypes
//		ResultYName VARCHAR(50),			-- result Y name
//      ResultYValue VARCHAR(500),			-- result Y value

tsSummStmSQL CSQLiteSummaries::m_StmSQL[7] = {
	{(char *)"TblExprs",
		(char *)"CREATE TABLE IF NOT EXISTS TblExprs (ExprID INTEGER PRIMARY KEY ASC,ExprName VARCHAR(50) UNIQUE, ExprTitle VARCHAR(50) UNIQUE,ExprDescr VARCHAR(1000) DEFAULT '')",
		(char *)"INSERT INTO TblExprs (ExprName,ExprTitle,ExprName,ExprDescr) VALUES(?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC)",
		NULL },
	{ (char *)"TblProcess",
		(char *)"CREATE TABLE IF NOT EXISTS TblProcess ( ProcessID INTEGER PRIMARY KEY ASC,ProcessName VARCHAR(50) UNIQUE,ProcessTitle VARCHAR(50),ProcessDescr VARCHAR(1000))",
		(char *)"INSERT INTO TblProcess (ProcessName,ProcessTitle,ProcessDescr) VALUES(?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblProcess_ProcessName' ON 'TblProcess' ('ProcessName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblProcess_ProcessName' ON 'TblProcess' ('ProcessName' ASC)",
		NULL },
	{ (char *)"TblProcessing",
		(char *)"CREATE TABLE IF NOT EXISTS TblProcessing (ProcessingID INTEGER PRIMARY KEY ASC,ExprID INTEGER,ProcessID INTEGER, ProcessVer VARCHAR(20), Start VARCHAR(24), Finish VARCHAR(24), ResultCode INTEGER)", 
		(char *)"INSERT INTO TblProcessing (ExprID,ProcessID,ProcessVer,Start,Finish,ResultCode) VALUES(?,?,?,?,?,?)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblProcessing_ExprIDProcessID'",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblProcessing_ExprIDProcessID' ON 'TblProcessing' ('ExprID' ASC,'ProcessID' ASC)",
		NULL },
	{ (char *)"TblProcessingLog",
		(char *)"CREATE TABLE IF NOT EXISTS TblProcessingLog (ProcessingLogID INTEGER PRIMARY KEY ASC,ProcessingID INTEGER,Timestamp VARCHAR(24),LogText VARCHAR(2000))",
		(char *)"INSERT INTO TblProcessingLog (ProcessingID,Timestamp,LogText) VALUES(?,?,?)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblProcessingLog_ProcessingID'",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblProcessingLog_ProcessingID' ON 'TblProcessingLog' ('ProcessingID' ASC)",
		NULL },

	{ (char *)"TblParams",
		(char *)"CREATE TABLE IF NOT EXISTS TblParams (ParamID INTEGER PRIMARY KEY ASC,ProcessingID INTEGER,ParamType INTEGER,ParmName VARCHAR(50),ParmValue VARCHAR(2000))",
		(char *)"INSERT INTO TblParams (ProcessingID,ParamType,ParmName,ParmValue) VALUES(?,?,?,?)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblParams_ParmName'",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblParams_ParmName' ON 'TblParams' ('ParmName' ASC)",
		NULL},

	{ (char *)"TblResults",
		(char *)"CREATE TABLE IF NOT EXISTS TblResults (ResultID INTEGER PRIMARY KEY ASC,ProcessingID INTEGER,GroupAs VARCHAR(50), ResultType INTEGER, ResultName VARCHAR(50),ResultValue VARCHAR(500))",
		(char *)"INSERT INTO TblResults (ProcessingID,GroupAs,ResultType,ResultName,ResultValue) VALUES(?,?,?,?,?)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblResults_ProcessingID'",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblResults_ProcessingID' ON 'TblResults' ('ProcessingID' ASC)",
		NULL },
	{ (char *)"TblXYResults",
		(char *)"CREATE TABLE IF NOT EXISTS TblXYResults (ResultID INTEGER PRIMARY KEY ASC,ProcessingID INTEGER,GroupAs VARCHAR(50), ResultXType INTEGER, ResultXName VARCHAR(50),ResultXValue VARCHAR(500),ResultYType INTEGER, ResultYName VARCHAR(50),ResultYValue VARCHAR(500))",
		(char *)"INSERT INTO TblXYResults (ProcessingID,GroupAs,ResultXType,ResultXName,ResultXValue,ResultYType,ResultYName,ResultYValue) VALUES(?,?,?,?,?,?,?,?)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblXYResults_ProcessingID'",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblXYResults_ProcessingID' ON 'TblXYResults' ('ProcessingID' ASC)",
		NULL },
};


CSQLiteSummaries::CSQLiteSummaries(void)
{
m_pDB = NULL;
m_bInitialised = Init() == eBSFSuccess ? true : false;
}


CSQLiteSummaries::~CSQLiteSummaries(void)
{
if(m_pDB != NULL)
	{
	sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	}
#ifndef _WIN32
if(m_bInitialised)
	pthread_spin_destroy(&m_hSpinLock);
#endif
}

int
CSQLiteSummaries::Init(void)
{
#ifdef _WIN32
if(!InitializeCriticalSectionAndSpinCount(&m_hSCritSect,1000))
	{
#else
if(pthread_spin_init(&m_hSpinLock,PTHREAD_PROCESS_PRIVATE)!=0)
	{
#endif
	return(eBSFerrInternal);
	}
sqlite3_initialize();
return(eBSFSuccess);
}

// serialise access to base flags in the upper nibble with base in lower 4 bits
inline void
CSQLiteSummaries::SQLiteSerialise(void)
{
int SpinCnt = 5000;
#ifdef _WIN32
while(!TryEnterCriticalSection(&m_hSCritSect))
	{
	if(SpinCnt -= 1)
		continue;
	SwitchToThread();
	SpinCnt = 500;
	}
#else
while(pthread_spin_trylock(&m_hSpinLock)==EBUSY)
	{
	if(SpinCnt -= 1)
		continue;
	pthread_yield();
	SpinCnt = 500;
	}
#endif
}

inline void
CSQLiteSummaries::SQLiteRelease(void)
{
#ifdef _WIN32
LeaveCriticalSection(&m_hSCritSect);
#else
pthread_spin_unlock(&m_hSpinLock);
#endif
}


sqlite3 *
CSQLiteSummaries::OpenDatabase(char *pszDatabase,		// database to open, if not already existing then will be created
							bool bReplace)			// if database already exists then replace
{
tsSummStmSQL *pStms;
int TblIdx;
int sqlite_error;

// note if database already exists in case bReplace is requested
struct stat TargStat;
int StatRslt = stat(pszDatabase,&TargStat);
if(StatRslt >= 0 && bReplace)
	remove(pszDatabase);

// try opening/creating the database 
if((sqlite_error = sqlite3_open_v2(pszDatabase, &m_pDB,SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't open database: %s", sqlite3_errmsg(m_pDB));
	sqlite3_shutdown();
	m_pDB = NULL;
	return(NULL);
	}

// create all tables which may not already exist
pStms = m_StmSQL;
for(TblIdx = 0; TblIdx < 7; TblIdx++,pStms++)
	{
	pStms->pPrepInsert = NULL;
	if(pStms->pszCreateTbl == NULL)
		continue;
	if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszCreateTbl,0,0,0))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszCreateTbl);   
		sqlite3_close_v2(m_pDB);
		sqlite3_shutdown();
		m_pDB = NULL;
		return(NULL);
		}
	}


pStms = m_StmSQL;
for(TblIdx = 0; TblIdx < 7; TblIdx++,pStms++)
	{
	if(pStms->pszOpenCreateIndexes == NULL)
		continue;
	if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszOpenCreateIndexes,0,0,0))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create indexes on table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszOpenCreateIndexes);   
		sqlite3_close_v2(m_pDB);
		sqlite3_shutdown();
		m_pDB = NULL;
		return(NULL);
		}
	}

// prepare all insert statements
pStms = m_StmSQL;
for(TblIdx = 0; TblIdx < 7; TblIdx++,pStms++)
	{
	if(pStms->pszInsert == NULL)
		{
		pStms->pPrepInsert = NULL;
		continue;
		}
	if((sqlite_error = sqlite3_prepare_v2(m_pDB,pStms->pszInsert,-1,&pStms->pPrepInsert,NULL))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't prepare insert statement on table %s: %s", pStms->pTblName, sqlite3_errmsg(m_pDB));
		while(TblIdx > 0)
			{
			pStms -= 1;
			if(pStms->pPrepInsert != NULL)
				{
				sqlite3_finalize(pStms->pPrepInsert);
				pStms->pPrepInsert = NULL;
				}
			}
		sqlite3_close_v2(m_pDB);
		sqlite3_shutdown();
		m_pDB = NULL;
		return(NULL);
		}
	}
return(m_pDB);
}

int
CSQLiteSummaries::CloseDatabase(void)
{
int TblIdx;
int Rslt = 0;
tsSummStmSQL *pStms;
pStms = m_StmSQL;
if(m_pDB != NULL)
	{
	for(TblIdx = 0; TblIdx < 7; TblIdx++,pStms++)
		{
		if(pStms->pPrepInsert == NULL)
			continue;
		sqlite3_finalize(pStms->pPrepInsert);
		pStms->pPrepInsert = NULL;
		}
	Rslt = sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	}
return(Rslt);
}



char *
CSQLiteSummaries::RemoveQuotes(char *pszRawText)
{
char *pSrcChr;
char *pDstChr;
bool bInSpace;
char Chr;
CUtility::TrimQuotedWhitespcExtd(pszRawText);
pSrcChr = pszRawText;
pDstChr = pSrcChr;
bInSpace = false;
while((Chr = *pSrcChr++)!= '\0')
	{
	if(Chr == '\'' || Chr == '"')
		continue;
	if(Chr == ' ' || Chr == '\t')
		{
		if(bInSpace)
			continue;
		bInSpace = true;
		}
	else
		bInSpace = false;
	*pDstChr++ = Chr;
	}
*pDstChr = '\0';
return(pszRawText);
}


// callbacks from sqlite3_exec returning an identifier
int CSQLiteSummaries::ExecCallbackID(void *pCallP1, // callback function processing identifier (4th arg to sqlite3_exec())
					int NumCols,			// number of result columns 
					char **ppColValues,		// array of ptrs to column values 
					char **ppColName)		// array of ptrs to column names
{
int ValChars;
int *pID;
char *ppEnd;

// some basic validation of call back parameter values
if(pCallP1 == NULL || NumCols != 1 || ppColValues == NULL || ppColValues[0] == NULL || *ppColValues[0] == '\0')
	return(1);

pID = (int *)pCallP1;
ValChars = (int)strlen(ppColValues[0]);
*pID = strtol(ppColValues[0],&ppEnd,10);
return(0);
}

int														// returned experiment identifier
CSQLiteSummaries::StartExperiment(char *pszDatabase,	// summary results to this SQLite database
							bool bReplace,				// if false then append to existing, if true then replace any existing database
							bool bContinue,				// if true then try to reuse any existing experiment identifier with same experiment name
							char *pszExprimentName,		// experiment name
							char *pszExperimentTitle,	// experiment title
							char *pszExperimentDescr)   // describes experiment
{
tsSummStmSQL *pStm;
int sqlite_error;
int ExprID;
char szSQLStatement[cMaxSQLStatement];
char szExprName[cMaxNameLen+1];
char szExprTitle[cMaxNameLen+1];
char szExprDescr[cMaxDescrText+1];

if(pszExprimentName == NULL || pszExperimentTitle == NULL || pszExperimentDescr == NULL ||
   pszExprimentName[0] == '\0' || pszExperimentTitle[0] == '\0' || pszExperimentDescr[0] == '\0')
	return(eBSFerrParams);

strncpy(szExprName,pszExprimentName,sizeof(szExprName));
szExprName[sizeof(szExprName)-1] = '\0';
strncpy(szExprTitle,pszExperimentTitle,sizeof(szExprTitle));
szExprTitle[sizeof(szExprTitle)-1] = '\0';
strncpy(szExprDescr,pszExperimentDescr,sizeof(szExprDescr));
szExprDescr[sizeof(szExprDescr)-1] = '\0';

SQLiteSerialise();
if(OpenDatabase(pszDatabase,bReplace)== NULL)
	{
	SQLiteRelease();
	return(0);
	}
if(m_pDB == NULL)
	{
	SQLiteRelease();
	return(eBSFerrInternal);
	}

char *pszBeginTransaction = (char *)"BEGIN TRANSACTION";
char *pszEndTransaction = (char *)"END TRANSACTION";

char *pszPragmaSyncOff = (char *)"PRAGMA synchronous = OFF";
char *pszPragmaSyncOn = (char *)"PRAGMA synchronous = ON";
char *pszPragmaJournMem = (char *)"PRAGMA journal_mode = MEMORY";

// synchronous writes off
if((sqlite_error = sqlite3_exec(m_pDB,pszPragmaSyncOff,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't turn synchronous writes off: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

// bracket inserts as a single transaction
if((sqlite_error = sqlite3_exec(m_pDB,pszBeginTransaction,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't begin transactions: %s",sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

pStm = &m_StmSQL[0];

// experiment already known?
ExprID = -1;
sprintf(szSQLStatement,"select ExprID from TblExprs where ExprName LIKE '%s'",szExprName);
if((sqlite_error = sqlite3_exec(m_pDB,szSQLStatement,ExecCallbackID,&ExprID,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite3_exec - getting ExprID: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if(ExprID >= 1)					// >= 1 if already known 
	{
	SQLiteRelease();
	return(ExprID);
	}

// experiment previously unknown, create it...
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 1, szExprName,(int)strlen(szExprName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, szExprTitle,(int)strlen(szExprTitle)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, szExprDescr,(int)strlen(szExprDescr)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

ExprID = (int)sqlite3_last_insert_rowid(m_pDB);
SQLiteRelease();
return(ExprID);
}

int					// returned process identifier 
CSQLiteSummaries::AddProcess(char *pszProcessName,		// process name
							char *pszProcessTitle,		// process title
							char *pszProcessDescr)		// describes process
{
tsSummStmSQL *pStm;
int sqlite_error;
int ProcessID;

char szSQLStatement[cMaxSQLStatement];
char szProcessName[cMaxNameLen+1];
char szProcessTitle[cMaxNameLen+1];
char szProcessDescr[cMaxDescrText+1];

if(pszProcessName == NULL || pszProcessTitle == NULL || pszProcessDescr == NULL ||
   pszProcessName[0] == '\0' || pszProcessTitle[0] == '\0' || pszProcessDescr[0] == '\0')
	return(eBSFerrParams);

strncpy(szProcessName,pszProcessName,sizeof(szProcessName));
szProcessName[sizeof(szProcessName)-1] = '\0';
strncpy(szProcessTitle,pszProcessTitle,sizeof(szProcessTitle));
szProcessTitle[sizeof(szProcessTitle)-1] = '\0';
strncpy(szProcessDescr,pszProcessDescr,sizeof(szProcessDescr));
szProcessDescr[sizeof(szProcessDescr)-1] = '\0';


if(m_pDB == NULL)
	return(eBSFerrInternal);

pStm = &m_StmSQL[1];
SQLiteSerialise();
// process already known?
ProcessID = -1;
sprintf(szSQLStatement,"select ProcessID from TblProcess where ProcessName LIKE '%s'",szProcessName);
if((sqlite_error = sqlite3_exec(m_pDB,szSQLStatement,ExecCallbackID,&ProcessID,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite3_exec - getting ProcessID: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if(ProcessID >= 1)					// >= 1 if already known
	{
	SQLiteRelease();
	return(ProcessID);
	}
// process previously unknown, create it...
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 1, szProcessName,(int)strlen(szProcessName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, szProcessTitle,(int)strlen(szProcessTitle)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, szProcessDescr,(int)strlen(szProcessDescr)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

ProcessID = (int)sqlite3_last_insert_rowid(m_pDB);
SQLiteRelease();
return(ProcessID);
}

int					// length - strlen() - of timestamp string 
CSQLiteSummaries::GetTimeStamp(char *pszTimeStamp)	// copy timestamp into this
{
int LineLen;
#ifdef _WIN32
struct _timeb timebuffer;
#else
struct timeb timebuffer;
#endif
char *timeline;
#ifdef _WIN32
_ftime(&timebuffer);
#else
ftime(&timebuffer);
#endif
timeline = ctime(&timebuffer.time);
LineLen = sprintf(pszTimeStamp,"%.15s.%03d %.4s",&timeline[4],(int)timebuffer.millitm, &timeline[20]);
return(LineLen);
}

int													// uiniqely identifies this starting experiment process instance
CSQLiteSummaries::StartProcessing(int ExperimentID,	// identifier returned by StartExperiment()
						 int ProcessID,				// identifier as returned by AddProcess()
						 char *pszProcessVersion)	// process version
{
tsSummStmSQL *pStm;
int sqlite_error;
int ProcessingID;

if(m_pDB == NULL)
	return(eBSFerrInternal);

pStm = &m_StmSQL[2];
SQLiteSerialise();
char szTimestamp[cMaxNameLen];
GetTimeStamp(szTimestamp);
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExperimentID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, ProcessID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, pszProcessVersion,(int)strlen(pszProcessVersion)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, szTimestamp,(int)strlen(szTimestamp)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 5, szTimestamp,(int)strlen(szTimestamp)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}


if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 6, -1))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

ProcessingID = (int)sqlite3_last_insert_rowid(m_pDB);
SQLiteRelease();
return(ProcessingID);
}
				
int													// uiniqely identifies this starting experiment process log instance
CSQLiteSummaries::AddLog(int ProcessingID,			// identifier returned by StartProcessing()
						 const char *pszFormat,...) // printf style format
{
va_list Args;
char szLogText[cMaxDiagLen+1];
tsSummStmSQL *pStm;
int sqlite_error;
int ProcessingLogID;
if(m_pDB == NULL)
	return(eBSFerrInternal);

pStm = &m_StmSQL[3];
va_start(Args, pszFormat );
#ifdef _WIN32
_vsnprintf(szLogText,cMaxDiagLen,pszFormat,Args);
#else
vsnprintf(szLogText,cMaxDiagLen,pszFormat,Args);
#endif
szLogText[cMaxDiagLen-1] = '\0';

SQLiteSerialise();
char szTimestamp[cMaxNameLen];
GetTimeStamp(szTimestamp);

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ProcessingID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, szTimestamp,(int)strlen(szTimestamp)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, szLogText,(int)strlen(szLogText)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

ProcessingLogID = (int)sqlite3_last_insert_rowid(m_pDB);
SQLiteRelease();
return(ProcessingLogID);
}

int
CSQLiteSummaries::AddParameter(int ProcessingID,		// identifier returned by StartProcessing()
						 teSQLliteSummParamTypes ParamType,	// parameter type
						 int ValueSize,					// parameter value is of this byte length (if text then excludes terminating '\0')
						 const char *pszParamName,			// parameter name
						 void *pParmValue)				// parameter value
{
int Rslt;
int ParameterID;
tsSummStmSQL *pStm;
int sqlite_error;
char szParamValue[cMaxTextBuff];

if(m_pDB == NULL || 
   ProcessingID == 0 ||
   (ParamType != ePTText && ValueSize < 1) ||
   pszParamName == NULL || pszParamName[0] == '\0')
	return(eBSFerrInternal);

if(ParamType == ePTText && (pParmValue == NULL || ValueSize == 0 || *(char *)pParmValue == '\0'))
	{
	pParmValue = (void *)"N/A";
	ValueSize = 3;
	}

pStm = &m_StmSQL[4];
if((Rslt = ValueToText(ParamType,ValueSize,pParmValue,sizeof(szParamValue)-1,szParamValue)))
   return(Rslt);
SQLiteSerialise();

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ProcessingID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, ParamType))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, pszParamName,(int)strlen(pszParamName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, szParamValue,(int)strlen(szParamValue)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

ParameterID = (int)sqlite3_last_insert_rowid(m_pDB);

SQLiteRelease();
return(ParameterID);
}



int
CSQLiteSummaries::AddResult(int ProcessingID,			// identifier returned by StartProcessing()
						 const char *GroupAs,			// result is part of this grouping
						 teSQLliteSummParamTypes ParamType,	// result value type
						 int ValueSize,					// result value is of this byte length
						 const char *pszResultName,		// result name
						 void *pResultValue)			// result value
{
int Rslt;
int ResultID;
tsSummStmSQL *pStm;
int sqlite_error;
char szResultValue[cMaxTextBuff];

if(m_pDB == NULL || 
   ProcessingID == 0 ||
   (ParamType != ePTText && ValueSize < 1) ||
   GroupAs == NULL || GroupAs[0] == '\0' ||
   pszResultName == NULL || pszResultName[0] == '\0')
	return(eBSFerrInternal);

if(ParamType == ePTText && (pResultValue == NULL || ValueSize == 0 || *(char *)pResultValue == '\0'))
	{
	pResultValue = (void *)"N/A";
	ValueSize = 3;
	}

pStm = &m_StmSQL[5];

if((Rslt = ValueToText(ParamType,ValueSize,pResultValue,sizeof(szResultValue)-1,szResultValue)))
   return(Rslt);
SQLiteSerialise();
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ProcessingID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, GroupAs,(int)strlen(GroupAs)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, ParamType))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, pszResultName,(int)strlen(pszResultName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 5, szResultValue,(int)strlen(szResultValue)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);
ResultID = (int)sqlite3_last_insert_rowid(m_pDB);
SQLiteRelease();
return(ResultID);
}


int
CSQLiteSummaries::AddResultXY(int ProcessingID,			// identifier returned by StartProcessing()
						 const char *GroupAs,			// result is part of this grouping
						 teSQLliteSummParamTypes ParamXType,	// result value X type
						 int ValueXSize,				// result value is of this byte length
						 const char *pszResultXName,	// result X name
						 void *pResultXValue,			// result X value
 						 teSQLliteSummParamTypes ParamYType,	// result value Y type
						 int ValueYSize,				// result value is of this byte length
						 const char *pszResultYName,	// result Y name
						 void *pResultYValue)			// result Y value
{
int Rslt;
int ResultID;
tsSummStmSQL *pStm;
int sqlite_error;
char szResultXValue[cMaxTextBuff];
char szResultYValue[cMaxTextBuff];

if(m_pDB == NULL || 
   ProcessingID == 0 ||
   (ParamXType != ePTText && ValueXSize < 1) ||
   (ParamYType != ePTText && ValueYSize < 1) ||
   GroupAs == NULL || GroupAs[0] == '\0' ||
   pszResultXName == NULL || pszResultXName[0] == '\0' ||
   pszResultYName == NULL || pszResultYName[0] == '\0')
	return(eBSFerrInternal);

if(ParamXType == ePTText && (pResultXValue == NULL || ValueXSize == 0 || *(char *)pResultXValue == '\0'))
	{
	pResultXValue = (void *)"N/A";
	ValueXSize = 3;
	}

if(ParamYType == ePTText && (pResultYValue == NULL || ValueYSize == 0 || *(char *)pResultYValue == '\0'))
	{
	pResultYValue = (void *)"N/A";
	ValueYSize = 3;
	}

pStm = &m_StmSQL[6];

if((Rslt = ValueToText(ParamXType,ValueXSize,pResultXValue,sizeof(szResultXValue)-1,szResultXValue)))
   return(Rslt);
if((Rslt = ValueToText(ParamYType,ValueYSize,pResultYValue,sizeof(szResultYValue)-1,szResultYValue)))
   return(Rslt);

SQLiteSerialise();
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ProcessingID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, GroupAs,(int)strlen(GroupAs)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, ParamXType))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, pszResultXName,(int)strlen(pszResultXName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 5, szResultXValue,(int)strlen(szResultXValue)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 6, ParamYType))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 7, pszResultYName,(int)strlen(pszResultXName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 8, szResultYValue,(int)strlen(szResultYValue)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);
ResultID = (int)sqlite3_last_insert_rowid(m_pDB);
SQLiteRelease();
return(ResultID);
}



int					// returned process identifier
CSQLiteSummaries::EndProcessing(int ProcessingID,	// identifier returned by StartProcessing()
						  int ResultCode)			// processing completed result code
{
char szUpdate[cMaxSQLStatement];
int sqlite_error;
if(m_pDB == NULL)
	return(eBSFerrInternal);
SQLiteSerialise();
char szTimestamp[cMaxNameLen];
GetTimeStamp(szTimestamp);
sprintf(szUpdate,"UPDATE TblProcessing SET ResultCode = %d, Finish = \"%s\" WHERE ProcessingID = %d",ResultCode,szTimestamp,ProcessingID);
if((sqlite_error = sqlite3_exec(m_pDB,szUpdate,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't update processing result code: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}
SQLiteRelease();
return(ProcessingID);
}

int
CSQLiteSummaries::EndExperiment(int ExperimentID)			// identifier returned by StartExperiment()
{
int Rslt;
int sqlite_error;

char *pszEndTransaction = (char *)"END TRANSACTION";
char *pszPragmaSyncOn = (char *)"PRAGMA synchronous = ON";

SQLiteSerialise();
	// end transaction
if((sqlite_error = sqlite3_exec(m_pDB,pszEndTransaction,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't end transactions: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

tsSummStmSQL *pStms;
pStms = m_StmSQL;
int TblIdx;
for(TblIdx = 0; TblIdx < 6; TblIdx++,pStms++)
	{
	if(pStms->pszCreateIndexes == NULL)
		continue;
	if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszCreateIndexes,0,0,0))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create indexes on table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszCreateIndexes);   
		CloseDatabase();
		SQLiteRelease();
		return(eBSFerrInternal);
		}
	}
// synchronous writes off
if((sqlite_error = sqlite3_exec(m_pDB,pszPragmaSyncOn,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't turn synchronous writes on: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase();
	SQLiteRelease();
	return(eBSFerrInternal);
	}

Rslt = CloseDatabase();
SQLiteRelease();
return(Rslt);
}


int							// eBSFSuccess if pszValueText contains textified value
CSQLiteSummaries::ValueToText(teSQLliteSummParamTypes ParamType,	// result value type
			int ValueSize,			// result value is of this byte length
			void *pValue,			// result value
			int MaxTextLen,			// truncate returned value text at this many chars - remember to allow for terminating '\0'
			char *pszValueText)		// user allocated to hold returned value as text
{
char szBuff[cMaxTextBuff];

if(ValueSize < 1 || pszValueText == NULL || MaxTextLen < 1)
	return(eBSFerrParams);
pszValueText[0] = '\0';
switch(ParamType) {
	case ePTBool:					// text as 'False' or 'True'
		if(ValueSize != sizeof(bool))
			return(eBSFerrCvrtType);
		if(*(bool *)pValue)
			strcpy(szBuff,"True");
		else
			strcpy(szBuff,"False");
		break;
	case ePTInt32:					// signed int32
		if(ValueSize != sizeof(INT32))
			return(eBSFerrCvrtType);
		sprintf(szBuff,"%d",*(INT32 *)pValue);
		break;
	case ePTUint32:					// unsigned int32
		if(ValueSize != sizeof(UINT32))
			return(eBSFerrCvrtType);
		sprintf(szBuff,"%u",*(UINT32 *)pValue);
		break;
	case ePTInt64:					// signed int64
		if(ValueSize != sizeof(INT64))
			return(eBSFerrCvrtType);
		sprintf(szBuff,"%lld",*(INT64 *)pValue);
		break;
	case ePTUint64:					// unsigned int64
		if(ValueSize != sizeof(UINT64))
			return(eBSFerrCvrtType);
		sprintf(szBuff,"%llu",*(UINT64 *)pValue);
		break;
	case ePTDouble:					// double
		if(ValueSize != sizeof(double))
			return(eBSFerrCvrtType);
		sprintf(szBuff,"%0.f",*(double *)pValue);
		break;
	case ePTText:					// text
		strncpy(pszValueText,(char *)pValue,MaxTextLen);
		pszValueText[MaxTextLen] = '\0';
		return(eBSFSuccess);
	default:
		return(eBSFerrParams);
	}
strncpy(pszValueText,szBuff,MaxTextLen);
pszValueText[MaxTextLen] = '\0';
return(eBSFSuccess);
}


