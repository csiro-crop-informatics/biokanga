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

#include "SQLiteDE.h"

// Following database schema is utilised
// Tables
//	TblExprs		One row for each experiment
//  TblTrans        One row for each transcript
//  TblExpr         One row for each transcript differentially expressed
//  TblBins			One row for each bin count instance

// In each table the following columns are defined
//	TblExprs		One row for each experiment
//     ExprID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this experiment instance
//     ExprType INTEGER                     -- type, DE 2
//     ExprInFile VARCHAR(200),             -- Input CSV filename
//	   ExprName VARCHAR(50) UINQUE,			-- Short name of this experiment
//     ExprDescr VARCHAR(200),				-- Describes the experiment
//	   CtrlConditions VARCHAR(1000),		-- Describes the control conditions for the control transcriptome sequencing
//	   ExprConditions VARCHAR(1000),		-- Describes the experimental conditions for the experimental transcriptome sequencing
//     NumBins INTEGER                      -- Number of bins for Pearson analysis

//  TblTrans         One row for each transcript 
//     TransID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this transcript instance
//     ExprID INTEGER,						-- in this experiment
//	   TransName VARCHAR(80) UNIQUE,		-- Short name of this transcript
//     Exons INTEGER,					    -- number of exons
//	   TransLen  INTEGER,				    -- transcript length
//	   TransAnnotation VARCHAR(1000)		-- any known annotation for transcript

//  TblExpres         One row for each transcript differentially expressed
//     ExpresID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this expression instance
//     ExprID INTEGER,						-- in this experiment
//     TransID INTEGER,						-- expression is on this transcript
//     Class INTEGER                        -- classification
//	   Score INTEGER						-- Overall DE score
//     DECntsScore INTEGER,                 -- Counts only DE score
//     PearsonScore INTEGER,                -- Pearson only DE score
//     CtrlUniqueLoci INTEGER,				-- number of unique loci within control transcript at which read starts observed
//     ExprUniqueLoci INTEGER,				-- number of unique loci within experiment transcript at which read starts observed
//     CtrlExprLociRatio REAL				-- control/experiment unique loci ratio
//     PValueMedian REAL
//     PValueLow95 REAL
//     PValueHi95 REAL
//     TotCtrlCnts INTEGER
//     TotExprCnts INTEGER
//     TotCtrlExprCnts INTEGER
//     ObsFoldChange REAL
//     FoldMedian REAL
//     FoldLow95 REAL
//     FoldHi95 REAL
//     ObsPearson REAL
//     PearsonMedian REAL
//     PearsonLow95 REAL
//     PearsonHi95 REAL
//     CtrlAndExprBins INTEGER              -- intersect of control and experiment bins having at least one attributed count
//     CtrlOnlyBins INTEGER  				-- control bins having at least one attributed count with same bin in experiment having no counts
//     ExprOnlyBins INTEGER					-- experiment bins having at least one attributed count with same bin in control having no counts


//  TblBins     One row for each bin count set
//     BinID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this bin instance
//     ExprID INTEGER,						-- in this experiment
//		TransID INTEGER,					-- bins for this transcript
//		ExpresID INTEGER,					-- and this expression
//      NthBin  INTEGER,					-- Nth bin along the transcript
//      CtrlCounts INTEGER,				    -- control counts in this bin
//      ExprCounts INTEGER				    -- experiment counts in this bin

tsDEStmSQL CSQLiteDE::m_StmSQL[4] = {
	{(char *)"TblExprs",
		(char *)"CREATE TABLE TblExprs (ExprID INTEGER PRIMARY KEY ASC,ExprType Integer, ExprInFile VARCHAR(200), ExprName VARCHAR(50) UNIQUE,ExprDescr VARCHAR(200) DEFAULT '', CtrlConditions VARCHAR(1000),ExprConditions VARCHAR(1000),NumBins INTEGER)",
		(char *)"INSERT INTO TblExprs (ExprType,ExprInFile,ExprName,ExprDescr,CtrlConditions,ExprConditions,NumBins) VALUES(?,?,?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC)",
		NULL,
		NULL },
	{ (char *)"TblTrans",
		(char *)"CREATE TABLE TblTrans ( TransID INTEGER PRIMARY KEY ASC,ExprID INTEGER,TransName VARCHAR(80) UNIQUE,Exons Integer,TransLen INTEGER,TransAnnotation VARCHAR(1000))",
		(char *)"INSERT INTO TblTrans (ExprID,TransName,Exons,TransLen,TransAnnotation) VALUES(?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblTrans_TransName' ON 'TblTrans' ('TransName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblTrans_TransName' ON 'TblTrans' ('TransName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblTrans_TransName' ON 'TblTrans' ('TransName' ASC);CREATE INDEX IF NOT EXISTS 'TblTrans_TransLen' ON 'TblTrans' ('TransLen' ASC)",
		NULL },
	{ (char *)"TblExpres",
		(char *)"CREATE TABLE TblExpres (ExpresID INTEGER PRIMARY KEY ASC,ExprID INTEGER,TransID INTEGER," \
				"Class INTEGER,Score INTEGER,DECntsScore INTEGER, PearsonScore INTEGER,CtrlUniqueLoci INTEGER,ExprUniqueLoci INTEGER,CtrlExprLociRatio REAL," \
				"PValueMedian REAL,PValueLow95 REAL,PValueHi95 REAL,TotCtrlCnts INTEGER,TotExprCnts INTEGER,TotCtrlExprCnts INTEGER,ObsFoldChange REAL," \
				"FoldMedian REAL,FoldLow95 REAL,FoldHi95 REAL,ObsPearson REAL,PearsonMedian REAL,PearsonLow95 REAL, PearsonHi95 REAL," \
				"CtrlAndExprBins INTEGER ,CtrlOnlyBins INTEGER,ExprOnlyBins INTEGER)",
		(char *)"INSERT INTO TblExpres (ExprID,TransID,Class,Score,DECntsScore, PearsonScore,CtrlUniqueLoci,ExprUniqueLoci,CtrlExprLociRatio," \
				"PValueMedian,PValueLow95,PValueHi95,TotCtrlCnts,TotExprCnts,TotCtrlExprCnts,ObsFoldChange," \
				"FoldMedian,FoldLow95,FoldHi95,ObsPearson,PearsonMedian,PearsonLow95, PearsonHi95," \
				"CtrlAndExprBins,CtrlOnlyBins,ExprOnlyBins) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExpres_ExprIDTransID' ON 'TblExpres' ('ExprID' ASC,'TransID' ASC)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblExpres_ExprIDTransID';CREATE INDEX IF NOT EXISTS 'TblExpres_ExprID' ON 'TblExpres' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblExpres_TransID' ON 'TblExpres' ('TransID' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblExpres_ExprIDTransID';DROP INDEX IF EXISTS 'TblExpres_ExprID';DROP INDEX IF EXISTS 'TblExpres_TransID'"},
	{ (char *)"TblBins",
		(char *)"CREATE TABLE TblBins (BinID INTEGER PRIMARY KEY ASC,ExprID INTEGER,TransID INTEGER, ExpresID INTEGER, NthBin INTEGER,CtrlCounts INTEGER, ExprCounts INTEGER)",
		(char *)"INSERT INTO TblBins (ExprID,TransID,ExpresID,NthBin,CtrlCounts,ExprCounts) VALUES(?,?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblBins_ExprIDTransIDExpresIDNthBin' ON 'TblBins' ('ExprID' ASC,'TransID' ASC,'ExpresID' ASC, 'NthBin' ASC)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblBins_ExprIDTransIDExpresIDNthBin' ON 'TblBins' ('ExprID' ASC,'TransID' ASC,'ExpresID' ASC,'NthBin' ASC);CREATE INDEX IF NOT EXISTS 'TblBins_TransID' ON 'TblBins' ('TransID' ASC);CREATE INDEX IF NOT EXISTS 'TblBins_ExpresID' ON 'TblBins' ('ExpresID' ASC);CREATE INDEX IF NOT EXISTS 'TblBins_NthBin' ON 'TblBins' ('NthBin' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblBins_ExprIDTransIDExpresIDNthBin';DROP INDEX IF EXISTS 'TblBins_TransID';DROP INDEX IF EXISTS 'TblBins_ExpresID',DROP INDEX IF EXISTS 'TblBins_NthBin'"},
	};


char *
CSQLiteDE::RemoveQuotes(char *pszRawText)
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

CSQLiteDE::CSQLiteDE(void)
{
m_pDB = NULL;
m_NumTransMRA = 0;
m_NumTrans = 0;
m_NumExpres = 0;
m_bSafe = true;
}


CSQLiteDE::~CSQLiteDE(void)
{
if(m_pDB != NULL)
	{
	sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	}
}


sqlite3 *
CSQLiteDE::CreateDatabase(bool bSafe,				// true if sourcing from input CSV of unknown origin which may contain duplicates etc..
				char *pszDatabase)		// database to create (any existing database is deleted then clean created)
{
tsDEStmSQL *pStms;
int TblIdx;
int sqlite_error;
// note if database already exists in case bReplace is requested
struct stat TargStat;
int StatRslt = stat(pszDatabase,&TargStat);
if(StatRslt >= 0)
	remove(pszDatabase);

m_bSafe = bSafe;

// try creating the database
if((sqlite_error = sqlite3_open_v2(pszDatabase, &m_pDB,SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't open database: %s", sqlite3_errmsg(m_pDB));
	sqlite3_shutdown();
	m_pDB = NULL;
	return(NULL);
	}

// create all tables
pStms = m_StmSQL;
for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
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
if(bSafe)
	{
	for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
		{
		if(pStms->pszOpenCreateSafeIndexes == NULL)
			continue;
		if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszOpenCreateSafeIndexes,0,0,0))!=SQLITE_OK)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create safe indexes on table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszOpenCreateIndexes);   
			sqlite3_close_v2(m_pDB);
			sqlite3_shutdown();
			m_pDB = NULL;
			return(NULL);
			}
		}
	}
else
	{
	for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
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
	}

// prepare all insert statements
pStms = m_StmSQL;
for(TblIdx = 0; TblIdx <4; TblIdx++,pStms++)
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
CSQLiteDE::CloseDatabase(bool bNoIndexes)
{
int TblIdx;
int Rslt = 0;
tsDEStmSQL *pStms;
pStms = m_StmSQL;
if(m_pDB != NULL)
	{
	if(!bNoIndexes)
		{
		for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
			{
			if(pStms->pPrepInsert == NULL)
				continue;
			sqlite3_finalize(pStms->pPrepInsert);
			pStms->pPrepInsert = NULL;
			}
		}
	Rslt = sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	}
return(Rslt);
}

// callbacks from sqlite3_exec returning an identifier
int CSQLiteDE::ExecCallbackID(void *pCallP1, // callback function processing identifier (4th arg to sqlite3_exec())
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

int												// errors if < eBSFSuccess, if positive then the ExprID
CSQLiteDE::CreateExperiment(int CSVtype,		// 0 if short form, 1 if including individual bin counts
				char *pszInFile,				// CSV file containing expression analysis results
				char *pszExprName,				// experiment identifier
				char *pszExprDescr,				// describes experiment
				char *pszCtrlConditions,		// control conditions
				char *pszExprConditions,		// experiment conditions
				int NumBins)
{
int sqlite_error;
int ExprID;
char szExprName[128];
tsDEStmSQL *pStm;

if(m_pDB == NULL)
	return(eBSFerrInternal);

m_NumTransMRA = 0;
m_NumTrans = 0;
m_NumExpres = 0;

pStm = &m_StmSQL[0];
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, CSVtype))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, pszInFile,(int)strlen(pszInFile)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, pszExprName,(int)strlen(pszExprName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, pszExprDescr,(int)strlen(pszExprDescr)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 5, pszCtrlConditions,(int)strlen(pszCtrlConditions)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 6, pszExprConditions,(int)strlen(pszExprConditions)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 7, NumBins))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

	// find out the ExprID assigned to this experiment
if(m_bSafe)
	{
	ExprID = -1;
	sprintf(szExprName,"select ExprID from TblExprs where ExprName LIKE '%s'",pszExprName);
	if((sqlite_error = sqlite3_exec(m_pDB,szExprName,ExecCallbackID,&ExprID,NULL))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite3_exec - getting ExprID: %s", sqlite3_errmsg(m_pDB));   
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	}
else
	ExprID = (int)sqlite3_last_insert_rowid(m_pDB);

return(ExprID);
}

int										// returned identifier for transcript
CSQLiteDE::AddTrans(int ExprID,			// experiment identifier
		char *pszTransName,				// transcript name
		int NumExons,				    // number of exons in this transcript
		int TransLen,					// transcript length
		char *pszAnnotation)			// any associated annotations
{
int sqlite_error;
tsDEStmSQL *pStm;
int Idx;
int TransID;
tsMRATrans *pMRATrans;
tsMRATrans *pLRATrans;
char szQueryTransName[200];

if(m_pDB == NULL)
	return(eBSFerrInternal);

if(!m_NumTransMRA)
	memset(m_MRATrans,0,sizeof(m_MRATrans));

// quickly check if sequence is a recently accessed sequence and if so then return the identifier
pMRATrans = m_MRATrans;
for(Idx = 0; Idx < m_NumTransMRA; Idx++, pMRATrans++)
	{
	if(!stricmp(pszTransName,pMRATrans->szTransName))
		{
		if(pMRATrans->AccessCnt < 0x7ffffff)
			pMRATrans->AccessCnt += 10;
		return(pMRATrans->TransID);
		}
	if(pMRATrans->AccessCnt > 0)
		pMRATrans->AccessCnt -= 1;
	}

TransID = -1;
if(m_bSafe)
	{
	// not a recently accessed sequence so need to check if already known to SQLite
	sprintf(szQueryTransName,"select TransID from TblTrans where ExprID = %d AND TransName LIKE '%s'",ExprID,pszTransName);
	sqlite3_exec(m_pDB,szQueryTransName,ExecCallbackID,&TransID,NULL);
	}

if(TransID == -1)	// will be -1 if not already in database so need to add
	{
	pStm = &m_StmSQL[1];								// access sequence statements
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}

	if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, pszTransName,(int)strlen(pszTransName)+1,SQLITE_STATIC))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}

	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, NumExons))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 4, TransLen))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}

	if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 5, pszAnnotation,(int)strlen(pszAnnotation)+1,SQLITE_STATIC))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}

	if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	sqlite3_reset(pStm->pPrepInsert);

	if(m_bSafe)
		{
		sprintf(szQueryTransName,"select TransID from TblTrans where ExprID = %d AND TransName LIKE '%s'",ExprID,pszTransName);
		sqlite3_exec(m_pDB,szQueryTransName,ExecCallbackID,&TransID,NULL);
		}
	else
		TransID = (int)sqlite3_last_insert_rowid(m_pDB);
	m_NumTrans += 1;						// number of seqs added to TblSeqs
	}

// replace lRA sequence
if(m_NumTransMRA < cMaxMRATrans)
	pMRATrans = &m_MRATrans[m_NumTransMRA++];
else
	{
	pMRATrans = m_MRATrans;
	pLRATrans = pMRATrans++;
	for(Idx = 1; Idx < m_NumTransMRA; Idx++, pMRATrans++)
		{
		if(pMRATrans->AccessCnt < pLRATrans->AccessCnt)
			pLRATrans = pMRATrans;
		}
	pMRATrans = pLRATrans;
	}
pMRATrans->AccessCnt = 1000;
pMRATrans->TransID = TransID;
strcpy(pMRATrans->szTransName,pszTransName);
return(TransID);
}


int										// returned identifier for transcript
CSQLiteDE::AddExpres(int ExprID,		// experiment identifier
		int TransID,					// transcript identifier
		int Class,
		int Score,
		int DECntsScore,
		int PearsonScore,
		int CtrlUniqueLoci,
		int ExprUniqueLoci,
		double CtrlExprLociRatio,
		double PValueMedian,
		double PValueLow95,
		double PValueHi95,
		int TotCtrlCnts,
		int TotExprCnts,
		int TotCtrlExprCnts,
		double ObsFoldChange,
		double FoldMedian,
		double FoldLow95,
		double FoldHi95,
		double ObsPearson,
		double PearsonMedian,
		double PearsonLow95,
		double PearsonHi95,
		int CtrlAndExprBins,
		int CtrlOnlyBins,
		int ExprOnlyBins)
{
tsDEStmSQL *pStm;
int sqlite_error;
int ExpresID;
char szQueryExpresID[200];

pStm = &m_StmSQL[2];								// access sequence statements
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, TransID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, Class))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 4, Score))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 5, DECntsScore))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 6, PearsonScore))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 7, CtrlUniqueLoci))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 8, ExprUniqueLoci))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 9, CtrlExprLociRatio))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 10, PValueMedian))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 11, PValueLow95))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 12, PValueHi95))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 13, TotCtrlCnts))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 14, TotExprCnts))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 15, TotCtrlExprCnts))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 16, ObsFoldChange))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 17, FoldMedian))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 18, FoldLow95))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 19, FoldHi95))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 20, ObsPearson))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 21, PearsonMedian))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 22, PearsonLow95))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_double(pStm->pPrepInsert, 23, PearsonHi95))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 24, CtrlAndExprBins))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 25, CtrlOnlyBins))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 26, ExprOnlyBins))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

if(m_bSafe)
	{
	sprintf(szQueryExpresID,"select ExpresID from TblExpres where ExprID = %d AND TransID = %d",ExprID,TransID);
	sqlite3_exec(m_pDB,szQueryExpresID,ExecCallbackID,&ExpresID,NULL);
	}
else
	ExpresID = (int)sqlite3_last_insert_rowid(m_pDB);
m_NumExpres += 1;						// number of expressions added to TblExpres
return(ExpresID);
}

int
CSQLiteDE::AddBin(int ExprID,
				int TransID,
				int ExpresID,
				int NthBin,
				int CtrlCounts,
				int ExprCounts)
{
int sqlite_error;
tsDEStmSQL *pStm;
int BinID;
char szQueryBinID[200];
pStm = &m_StmSQL[3];								// access sequence statements
if(m_pDB == NULL)
	return(eBSFerrInternal);
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, TransID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, ExpresID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 4, NthBin))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 5, CtrlCounts))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 6, CtrlCounts))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}



if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

if(m_bSafe)
	{
	BinID = -1;
	sprintf(szQueryBinID,"select BinID from TblBins where ExprID = %d AND TransID = %d AND ExpresID = %d AND NthBin = %d",ExprID,TransID,ExpresID,NthBin);
	sqlite3_exec(m_pDB,szQueryBinID,ExecCallbackID,&BinID,NULL);
	}
else
	BinID = (int)sqlite3_last_insert_rowid(m_pDB);	
return(BinID);
}

int
CSQLiteDE::ProcessCSV2SQLite(int PMode,			// currently just the one mode...default is to parse from CSV and create/populate SQLite database
				  bool bSafe,					// if true then use indexing on all tables whilst inserting... much slower but perhaps safer if multiple threads ...
				  char *pszExprName,			// name by which this experiment is identified
				  char *pszExprDescr,			// describes experiment
				  char *pszCtrlConditions,		// control conditions
				  char *pszExprConditions,		// experiment conditions
				  char *pszInFile,				// parse from this input CSV file
				  char *pszDatabase)			// SQLite database file
{
int Rslt;
bool bExtdBins;			// false if CSV file does not contain bins, true if determined that CSV contains bin counts
int sqlite_error;
sqlite3_stmt *prepstatement = NULL;
int ExprID;		// experiment identifier
int TransID;
int ExpresID;
char *pszTransName;
char *pszContOrExpr;
int TransLen;
int NumExons;
int Class;
int Score;
int DECntsScore;
int PearsonScore;
int CtrlUniqueLoci;
int ExprUniqueLoci;
double CtrlExprLociRatio;
double PValueMedian;
double PValueLow95;
double PValueHi95;
int TotCtrlCnts;
int TotExprCnts;
int TotCtrlExprCnts;
double ObsFoldChange;
double FoldMedian;
double FoldLow95;
double FoldHi95;
double ObsPearson;
double PearsonMedian;
double PearsonLow95;
double PearsonHi95;
int CtrlAndExprBins;
int CtrlOnlyBins;
int ExprOnlyBins;
int TotBins;
int BinIdx;
int BinValues[100 * 2];		// sized to contain 100 bin counts for both control and experiment
int *pBinValue;
int NumFields;
int NumElsRead;
int NumBins;
int BinID;
int ExpNumBins;
int ExpNumFields;


// load CSV file and determine number of fields, from these it can be determined if the file also contains the bin counts
CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	CloseDatabase(true);
	return(eBSFerrObj);
	}
pCSV->SetMaxFields(cDESummMaxBinFields + 10);	// could be upto 100 bins, add 10 in case more fields than expected!
if((Rslt=pCSV->Open(pszInFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInFile);
	delete pCSV;
	return(Rslt);
	}
if((Rslt=pCSV->NextLine()) < 1)	// have to be able to read at least one!
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read any lines from %s",pszInFile);
	delete pCSV;
	return(Rslt);
	}
NumFields = pCSV->GetCurFields();		// expecting at least 24, if 24 then summary, if 27+ then summary plus bin counts
if(NumFields < cDESummFields)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CSV transcipt expression file '%s' are expected to contain at least %d fields, parsed %d",pszInFile,cDESummFields,NumFields);
	delete pCSV;
	return(eBSFerrFieldCnt);
	}
if(NumFields > cDESummFields && NumFields < cDESummMinBinFields)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected CSV file '%s' to contain either %d (no bins) or at least %d (with bins) fields, file has %d fields",pszInFile,cDESummFields,cDESummMinBinFields,NumFields);
	delete pCSV;
	return(eBSFerrFieldCnt);
	}
if(NumFields > cDESummMaxBinFields)					// if summary plus bins then expecting at most 100 bins
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected CSV file '%s' to contain no more than 100 bins",pszInFile);
	delete pCSV;
	return(eBSFerrFieldCnt);
	}
if(NumFields == cDESummFields)
	{
	bExtdBins = false;
	NumBins = 0;
	}
else
	{
	bExtdBins = true;
	NumBins = NumBins = 5 + NumFields - cDESummMinBinFields; 
	}
pCSV->Close();

sqlite3_initialize();

if((CreateDatabase(bSafe,pszDatabase))==NULL)
	{
	sqlite3_shutdown();
	return(eBSFerrInternal);
	}

if((Rslt = CreateExperiment(bExtdBins ? 1 : 0,pszInFile,pszExprName,pszExprDescr,pszCtrlConditions,pszExprConditions,NumBins)) < 1)
	{
	CloseDatabase(true);
	return(Rslt);
	}
ExprID = Rslt;

char *pszBeginTransaction = (char *)"BEGIN TRANSACTION";
char *pszEndTransaction = (char *)"END TRANSACTION";
char *pszDropIndexes = (char *)"DROP INDEX IF EXISTS 'Markers_LociID'";

char *pszPragmaSyncOff = (char *)"PRAGMA synchronous = OFF";
char *pszPragmaSyncOn = (char *)"PRAGMA synchronous = ON";
char *pszPragmaJournMem = (char *)"PRAGMA journal_mode = MEMORY";

gDiagnostics.DiagOut(eDLInfo,gszProcName,"sqlite - populating tables");


// synchronous writes off
if((sqlite_error = sqlite3_exec(m_pDB,pszPragmaSyncOff,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't turn synchronous writes off: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

// bracket inserts as a single transaction
if((sqlite_error = sqlite3_exec(m_pDB,pszBeginTransaction,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't begin transactions: %s",sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}


// load CSV file and start populating the SQLite database
if((Rslt=pCSV->Open(pszInFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInFile);
	delete pCSV;
	CloseDatabase(true);
	return(Rslt);
	}


bExtdBins = false;
NumElsRead = 0;
ExpNumBins = 0;
ExpNumFields = 0;
while((Rslt=pCSV->NextLine()) > 0)			// onto next line containing fields
	{
	if(!(NumElsRead % (bSafe ? 5000 : 100000)) && NumElsRead > 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d CSV lines - Transcripts: %d",NumElsRead, m_NumTrans);

	NumFields = pCSV->GetCurFields();		// expecting at least 24, if 24 then summary, if 27+ then summary plus bin counts
	if(ExpNumFields > 0 && NumFields != ExpNumFields)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CSV transcipt expression file '%s' has varying number of fields : expected %d, parsed %d at line %d",pszInFile,ExpNumFields,NumFields,NumElsRead);
		delete pCSV;
		CloseDatabase(true);
		return(eBSFerrFieldCnt);
		}
	if(ExpNumFields == 0)
		{
		ExpNumFields = NumFields;
		if(NumFields < cDESummFields)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"CSV transcipt expression file '%s' are expected to contain at least %d fields, parsed %d at line %d",pszInFile,cDESummFields,NumFields,NumElsRead);
			delete pCSV;
			CloseDatabase(true);
			return(eBSFerrFieldCnt);
			}
		if(NumFields > cDESummFields)							// summary plus bins?
			{
			if(NumFields < cDESummMinBinFields)					// if summary plus bins then expecting at least 5 bins
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected CSV file '%s' to contain at least 5 bins at line %d",pszInFile,NumElsRead);
				delete pCSV;
				CloseDatabase(true);
				return(eBSFerrFieldCnt);
				}
			if(NumFields > cDESummMaxBinFields)					// if summary plus bins then expecting at most 100 bins
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected CSV file '%s' to contain no more than 100 bins at line %d",pszInFile,NumElsRead);
				delete pCSV;
				CloseDatabase(true);
				return(eBSFerrFieldCnt);
				}

			NumBins = 5 + NumFields - cDESummMinBinFields;
			ExpNumBins = NumBins;
			bExtdBins = true;
			}
		}

	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;
	NumElsRead += 1;

	// first 20 fields are common to both CSV file types
	pCSV->GetInt(1,&Class);
	pCSV->GetText(2,&pszTransName);
	pCSV->GetInt(3,&TransLen);
	pCSV->GetInt(4,&NumExons);

	pCSV->GetInt(5,&Score);
	pCSV->GetInt(6,&DECntsScore);
	pCSV->GetInt(7,&PearsonScore);
	pCSV->GetInt(8,&CtrlUniqueLoci);
	pCSV->GetInt(9,&ExprUniqueLoci);
	
	pCSV->GetDouble(10,&CtrlExprLociRatio);
	pCSV->GetDouble(11,&PValueMedian);
	pCSV->GetDouble(12,&PValueLow95);
	pCSV->GetDouble(13,&PValueHi95);

	pCSV->GetInt(14,&TotCtrlCnts);
	pCSV->GetInt(15,&TotExprCnts);
	pCSV->GetInt(16,&TotCtrlExprCnts);

	pCSV->GetDouble(17,&ObsFoldChange);
	pCSV->GetDouble(18,&FoldMedian);
	pCSV->GetDouble(19,&FoldLow95);
	pCSV->GetDouble(20,&FoldHi95);
	// file formats diverge dependent on if containing bin counts
	if(!bExtdBins)
		{
		pszContOrExpr = NULL;
		pCSV->GetDouble(21,&ObsPearson);
		pCSV->GetDouble(22,&PearsonMedian);
		pCSV->GetDouble(23,&PearsonLow95);
		pCSV->GetDouble(24,&PearsonHi95);
		TotBins = 0;
		CtrlAndExprBins = 0;
		CtrlOnlyBins = 0;
		ExprOnlyBins = 0;
		}
	else
		{
		pCSV->GetText(21,&pszContOrExpr);
		pCSV->GetDouble(22,&ObsPearson);
		pCSV->GetDouble(23,&PearsonMedian);
		pCSV->GetDouble(24,&PearsonLow95);
		pCSV->GetDouble(25,&PearsonHi95);
		pCSV->GetInt(26,&TotBins);
		pCSV->GetInt(27,&CtrlAndExprBins);
		pCSV->GetInt(28,&CtrlOnlyBins);
		pCSV->GetInt(29,&ExprOnlyBins);
		if(!stricmp(pszContOrExpr,"Control"))
			pBinValue = &BinValues[0];
		else
			pBinValue = &BinValues[1];
		for(BinIdx = 0; BinIdx < NumBins; BinIdx++,pBinValue += 2)
			*pBinValue = pCSV->GetInt(BinIdx + 30,pBinValue);
		}
	if(!bExtdBins || (bExtdBins && stricmp(pszContOrExpr,"Control")))
		{
		TransID = AddTrans(ExprID,pszTransName,NumExons,TransLen,(char *)"N/A");
		ExpresID = AddExpres(ExprID,TransID,Class,Score,DECntsScore,PearsonScore,CtrlUniqueLoci,ExprUniqueLoci,CtrlExprLociRatio,PValueMedian,PValueLow95,PValueHi95,
					TotCtrlCnts,TotExprCnts,TotCtrlExprCnts,ObsFoldChange,FoldMedian,FoldLow95,FoldHi95,ObsPearson,PearsonMedian,PearsonLow95,PearsonHi95,
					CtrlAndExprBins,CtrlOnlyBins,ExprOnlyBins);
		}
	if(bExtdBins && stricmp(pszContOrExpr,"Control"))
		{
		pBinValue = &BinValues[0];
		for(BinIdx = 1; BinIdx <= NumBins; BinIdx++,pBinValue += 2)
			BinID = AddBin(ExprID,TransID,ExpresID,BinIdx,*pBinValue,pBinValue[1]);
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d CSV lines - transcripts: %d",NumElsRead, m_NumTrans);

	// end transaction
if((sqlite_error = sqlite3_exec(m_pDB,pszEndTransaction,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't end transactions on '%s': %s", "Markers",sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed populating the sqlite database");

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating indexes ...");

tsDEStmSQL *pStms;
pStms = m_StmSQL;
int TblIdx;
for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
	{
	if(pStms->pszCreateIndexes == NULL)
		continue;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating indexes on table %s ...", pStms->pTblName);
	if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszCreateIndexes,0,0,0))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create indexes on table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszCreateIndexes);   
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Indexes generated");
// synchronous writes off
if((sqlite_error = sqlite3_exec(m_pDB,pszPragmaSyncOn,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't turn synchronous writes on: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

CloseDatabase();
sqlite3_shutdown();
gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database ready for use");
return(eBSFSuccess);
}

