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

#include "SQLiteMarkers.h"

// Following database schema is utilised
// Tables
//	TblExprs		One row for each experiment
//  TblCults        One row for each cultivar or species
//  TblSeqs         One row for each sequence (DNA chromosome or assembly contig if genomic, transcript if RNA) 
//  TblLoci			One row for each loci identified on a sequence plus the cannonical base at that loci 
//  TblSnps			One row for each SNP at an identified loci for each cultivar and relavant experiment
//  TblMarkers      One row for each marker at an identified loci plus cultivar specific base and score
//  TblMarkerSnps   One row for each SNP processed and contributing to a specific identified marker

// In each table the following columns are defined
//	TblExprs		One row for each experiment
//     ExprID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this experiment instance
//     ExprType INTEGER                     -- type, markers 0 or SNPs 1
//     ExprInFile VARCHAR(200),             -- Input CSV filename
//	   ExprName VARCHAR(50) UINQUE,			-- Short name of this experiment
//     ExprDescr VARCHAR(200),				-- Describes the experiment
//	   CultName VARCHAR(50),				-- Short name of target cultivar against which alignments were made

//  TblCults        One row for each cultivar or species
//     CultID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this cultivar instance
//	   CultName VARCHAR(50) UNIQUE,			-- Short name of this cultivar

//  TblSeqs         One row for each sequence (DNA chromosome or assembly contig if genomic, transcript if RNA) 
//     SeqID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this sequence instance
//     ExprID INTEGER,						-- Sequence was target in this experiment
//	   SeqName VARCHAR(50),					-- Short name of this sequence

//  TblLoci			One row for each loci identified on a sequence plus the cannonical base at that loci
//     LociID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this loci instance
//     ExprID INTEGER,						-- Loci in this experiment
//     SeqID INTEGER,						-- Loci is on this sequence instance
//     Offset INTEGER,                      -- And at this offset on the sequence
//     Base VARCHAR(1),						-- With the sequence offset having this cannonical base

//  TblSnps			One row for each SNP at an identified loci for each cultivar and relavant experiment
//     SnpID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this SNP instance
//     ExprID INTEGER,						-- SNP identified within this experiment
//     CultID INTEGER,						-- SNP is in this cultivar relative to the experiment target cultivar
//     LociID INTEGER,						-- Identifies the loci instance at which the SNP is being called
//     Acnt INTEGER DEFAULT 0,			    -- Number of bases A in relative cultivar reads covering the SNP loci
//     Ccnt INTEGER DEFAULT 0,			    -- Number of bases C in relative cultivar reads covering the SNP loci
//     Gcnt INTEGER DEFAULT 0,			    -- Number of bases G in relative cultivar reads covering the SNP loci
//     Tcnt INTEGER DEFAULT 0,			    -- Number of bases T in relative cultivar reads covering the SNP loci
//     Ncnt INTEGER DEFAULT 0,			    -- Number of bases N in in relative cultivar reads covering the SNP loci
//     TotCovCnt INTEGER DEFAULT 0,			-- Total number of bases in relative cultivar reads covering the SNP loci
//     TotMMCnt INTEGER DEFAULT 0,			-- Total number of mismatches bases in relative cultivar reads covering the SNP loci

//  TblMarkers      One row for each marker at an identified loci plus cultivar specific base and score
//     MarkerID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this marker instance
//     ExprID INTEGER,		                -- marker identified within this experiment
//     CultID INTEGER,			            -- marker is in this cultivar relative to the experiment target cultivar
//     LociID INTEGER,			            -- Identifies the loci instance at which the marker is being called
//     Base VARCHAR(1),						-- Called marker base
//     Score INTEGER,						-- Called marker score

//  TblMarkerSnps   One row for each SNP processed and contributing to a specific identified marker
//     MarkerSnpsID INTEGER PRIMARY KEY ASC,-- Uniquely identifies this marker SNP instance
//     SnpID INTEGER,						-- Identifies SNP instance
//     MarkerID INTEGER,					-- Used to generate this marker instance

tsStmSQL CSQLiteMarkers::m_StmSQL[7] = {
	{(char *)"TblExprs",
		(char *)"CREATE TABLE TblExprs (ExprID INTEGER PRIMARY KEY ASC,ExprType Integer, ExprInFile VARCHAR(200), ExprName VARCHAR(50) UNIQUE,ExprDescr VARCHAR(200) DEFAULT '', CultName VARCHAR(50),CultDescr VARCHAR(1000) DEFAULT '')",
		(char *)"INSERT INTO TblExprs (ExprType,ExprInFile,ExprName,ExprDescr,CultName,CultDescr) VALUES(?,?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC);CREATE INDEX IF NOT EXISTS 'TblExprs_CultName' ON 'TblExprs' ('CultName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC);CREATE INDEX IF NOT EXISTS 'TblExprs_CultName' ON 'TblExprs' ('CultName' ASC)",
		NULL,
		NULL },
	{ (char *)"TblCults",
		(char *)"CREATE TABLE TblCults ( CultID INTEGER PRIMARY KEY ASC,CultName VARCHAR(50) UNIQUE)",
		(char *)"INSERT INTO TblCults (CultName) VALUES(?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblCults_CultName' ON 'TblCults' ('CultName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblCults_CultName' ON 'TblCults' ('CultName' ASC)",
		NULL,
		NULL },
	{ (char *)"TblSeqs",
		(char *)"CREATE TABLE TblSeqs (SeqID INTEGER PRIMARY KEY ASC,ExprID INTEGER,SeqName VARCHAR(50))",
		(char *)"INSERT INTO TblSeqs (ExprID,SeqName) VALUES(?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblSeqs_ExprIDSeqName' ON 'TblSeqs' ('ExprID' ASC,'SeqName' ASC)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblSeqs_ExprIDSeqName';CREATE INDEX IF NOT EXISTS 'TblSeqs_ExprID' ON 'TblSeqs' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblSeqs_SeqName' ON 'TblSeqs' ('SeqName' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblSeqs_ExprIDSeqName';DROP INDEX IF EXISTS 'TblSeqs_ExprID';DROP INDEX IF EXISTS 'TblSeqs_SeqName'"},
	{ (char *)"TblLoci",
		(char *)"CREATE TABLE TblLoci (LociID INTEGER PRIMARY KEY ASC,ExprID INTEGER, SeqID INTEGER,Offset INTEGER, Base VARCHAR(1))",
		(char *)"INSERT INTO TblLoci (ExprID,SeqID,Offset,Base) VALUES(?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblLoci_ExprIDSeqIDOffset' ON 'TblLoci' ('ExprID' ASC,'SeqID' ASC,'Offset' ASC)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblLoci_ExprIDSeqIDOffset';CREATE INDEX IF NOT EXISTS 'TblLoci_ExprIDSeqIDOffset' ON 'TblLoci' ('ExprID' ASC,'SeqID' ASC,'Offset' ASC);CREATE INDEX IF NOT EXISTS 'TblLoci_SeqID' ON 'TblLoci' ('SeqID' ASC);CREATE INDEX IF NOT EXISTS 'TblLoci_SeqIDOffset' ON 'TblLoci' ('SeqID' ASC,'Offset' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblLoci_SeqID';DROP INDEX IF EXISTS 'TblLoci_SeqIDOffset'" },
	{ (char *)"TblSnps",	
		(char *)"CREATE TABLE TblSnps (SnpID INTEGER PRIMARY KEY ASC,ExprID INTEGER,CultID INTEGER,LociID INTEGER,Acnt INTEGER DEFAULT 0,Ccnt INTEGER DEFAULT 0,Gcnt INTEGER DEFAULT 0,Tcnt INTEGER DEFAULT 0,Ncnt INTEGER DEFAULT 0,TotCovCnt INTEGER DEFAULT 0,TotMMCnt INTEGER DEFAULT 0)",
		(char *)"INSERT INTO TblSnps (ExprID,CultID,LociID,Acnt,Ccnt,Gcnt,Tcnt,Ncnt,TotCovCnt,TotMMCnt) VALUES(?,?,?,?,?,?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblSnps_ExprIDCultIDLociID' ON 'TblSnps' ('ExprID' ASC,'CultID' ASC, 'LociID' ASC)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblSnps_ExprIDCultIDLociID';CREATE INDEX IF NOT EXISTS 'TblSnps_ExprID' ON 'TblSnps' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblSnps_CultID' ON 'TblSnps' ('CultID' ASC);CREATE INDEX IF NOT EXISTS 'TblSnps_LociID' ON 'TblSnps' ('LociID' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblSnps_ExprID'; DROP INDEX IF EXISTS 'TblSnps_CultID';DROP INDEX IF EXISTS 'TblSnps_LociID'"},
	{ (char *)"TblMarkers",	
		(char *)"CREATE TABLE TblMarkers (MarkerID INTEGER PRIMARY KEY ASC,ExprID INTEGER,CultID INTEGER,LociID INTEGER,Base VARCHAR(1),Score INTEGER)",
		(char *)"INSERT INTO TblMarkers (ExprID, CultID,LociID,Base,Score) VALUES(?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblMarkers_ExprIDCultIDLociID' ON 'TblMarkers' ('ExprID' ASC,'CultID' ASC,'LociID' ASC)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblMarkers_ExprIDCultIDLociID';CREATE INDEX IF NOT EXISTS 'TblMarkers_ExprID' ON 'TblMarkers' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblMarkers_CultID' ON 'TblMarkers' ('CultID' ASC);CREATE INDEX IF NOT EXISTS 'TblMarkers_LociID' ON 'TblMarkers' ('LociID' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblMarkers_ExprIDCultIDLociID';DROP INDEX IF EXISTS 'TblMarkers_ExprID'; DROP INDEX IF EXISTS 'TblMarkers_CultID';DROP INDEX IF EXISTS 'TblMarkers_LociID'"},
	{ (char *)"TblMarkerSnps",	
		(char *)"CREATE TABLE TblMarkerSnps (MarkerSnpsID INTEGER PRIMARY KEY ASC,ExprID INTEGER, MarkerID INTEGER, SnpID INTEGER)",
		(char *)"INSERT INTO TblMarkerSnps (ExprID,MarkerID,SnpID) VALUES(?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblMarkerSnps_ExprIDMarkerIDSnpID' ON 'TblMarkerSnps' ('ExprID' ASC,'MarkerID' ASC, 'SnpID' ASC)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblMarkerSnps_ExprIDMarkerIDSnpID';CREATE INDEX IF NOT EXISTS 'TblMarkerSnps_SnpID' ON 'TblMarkerSnps' ('SnpID' ASC);CREATE INDEX IF NOT EXISTS 'TblMarkerSnps_MarkerID' ON 'TblMarkerSnps' ('MarkerID' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblMarkerSnps_ExprIDMarkerIDSnpID';DROP INDEX IF EXISTS 'TblMarkerSnps_SnpID';DROP INDEX IF EXISTS 'TblMarkerSnps_MarkerID'"}
	};


char *
CSQLiteMarkers::RemoveQuotes(char *pszRawText)
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

CSQLiteMarkers::CSQLiteMarkers(void)
{
m_pDB = NULL;
m_NumSeqMRA = 0;
m_NumSeqs = 0;
m_NumSNPLoci = 0;
m_NumSNPs = 0;
m_NumMarkers = 0;
m_bSafe = true;
}


CSQLiteMarkers::~CSQLiteMarkers(void)
{
if(m_pDB != NULL)
	{
	sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	}
}



sqlite3 *
CSQLiteMarkers::CreateDatabase(bool bSafe,				// true if sourcing from input CSV of unknown origin which may contain duplicates etc..
				char *pszDatabase)		// database to create (any existing database is deleted then clean created)
{
tsStmSQL *pStms;
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
if(bSafe)
	{
	for(TblIdx = 0; TblIdx < 7; TblIdx++,pStms++)
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
CSQLiteMarkers::CloseDatabase(bool bNoIndexes)
{
int TblIdx;
int Rslt = 0;
tsStmSQL *pStms;
pStms = m_StmSQL;
if(m_pDB != NULL)
	{
	if(!bNoIndexes)
		{
		for(TblIdx = 0; TblIdx < 7; TblIdx++,pStms++)
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

// callbacks from sqlite3_exec
int CSQLiteMarkers::ExecCallbackCultID(void *pCallP1,	// callback function processing identifier (4th arg to sqlite3_exec())
					int NumCols,			// number of result columns 
					char **ppColValues,		// array of ptrs to column values 
					char **ppColName)		// array of ptrs to column names
{
int ValChars;
tsCultivar *pCult;
char *ppEnd;

// some basic validation of call back parameter values
if(pCallP1 == NULL || NumCols != 1 || ppColValues == NULL || ppColValues[0] == NULL || *ppColValues[0] == '\0')
	return(1);

pCult = (tsCultivar *)pCallP1;
ValChars = (int)strlen(ppColValues[0]);

pCult->CultID = strtol(ppColValues[0],&ppEnd,10);
return(0);
}


// callbacks from sqlite3_exec returning an identifier
int CSQLiteMarkers::ExecCallbackID(void *pCallP1, // callback function processing identifier (4th arg to sqlite3_exec())
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
CSQLiteMarkers::CreateExperiment(int CSVtype,	// 0 if markers, 1 if SNPs
				char *pszInFile,				// file containing SNPs or markers
				char *pszExprName,				// experiment identifier
				char *pszExprDescr,				// describes experiment
				char *pszAssembName)			// targeted assembly
{
int sqlite_error;
int ExprID;
char szExprName[128];
tsStmSQL *pStm;

if(m_pDB == NULL)
	return(eBSFerrInternal);

m_NumSeqs = 0;
m_NumSNPLoci = 0;
m_NumSNPs = 0;
m_NumMarkers = 0;

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
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 5, pszAssembName,(int)strlen(pszAssembName)+1,SQLITE_STATIC))!=SQLITE_OK)
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

int															// errors if < eBSFSuccess
CSQLiteMarkers::CreateCultivars(int NumCultivars,				// number of cultivars to add
				tsCultivar *pCultivars)			// pts to array of cultivars
{
int sqlite_error;
int CultIdx;
tsStmSQL *pStm;
tsCultivar *pCult;
char szCultTarg[200]; 
if(m_pDB == NULL)
	return(eBSFerrInternal);

pCult = pCultivars;
pStm = &m_StmSQL[1];								// access cultivar statements
for(CultIdx = 0; CultIdx < NumCultivars; CultIdx++,pCult++)
	{
	if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 1, pCult->szCultivarName,(int)strlen(pCult->szCultivarName)+1,SQLITE_STATIC))!=SQLITE_OK)
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


	// find out the CultID assigned to this cultivar
	if(m_bSafe)
		{
		sprintf(szCultTarg,"select CultID from TblCults where CultName LIKE '%s'",pCult->szCultivarName);
		sqlite3_exec(m_pDB,szCultTarg,ExecCallbackCultID,pCult,NULL);
		}
	else
		pCult->CultID = (int)sqlite3_last_insert_rowid(m_pDB);
	}

return(NumCultivars);
}





int										// returned sequence identifier for sequence
CSQLiteMarkers::AddSeq(	int ExprID,						// experiment identifier
		char *pszSeqName)				// target assembly sequence name
{
int sqlite_error;
tsStmSQL *pStm;
int Idx;
int SeqID;
tsMRASeq *pMRASeq;
tsMRASeq *pLRASeq;
char szSeqTarg[200];

if(m_pDB == NULL)
	return(eBSFerrInternal);

if(!m_NumSeqMRA)
	memset(m_MRASeqs,0,sizeof(m_MRASeqs));

// quickly check if sequence is a recently accessed sequence and if so then return the identifier
pMRASeq = m_MRASeqs;
for(Idx = 0; Idx < m_NumSeqMRA; Idx++, pMRASeq++)
	{
	if(!stricmp(pszSeqName,pMRASeq->szSeqName))
		{
		if(pMRASeq->AccessCnt < 0x7ffffff)
			pMRASeq->AccessCnt += 10;
		return(pMRASeq->SeqID);
		}
	if(pMRASeq->AccessCnt > 0)
		pMRASeq->AccessCnt -= 1;
	}

SeqID = -1;
if(m_bSafe)
	{
	// not a recently accessed sequence so need to check if already known to SQLite
	sprintf(szSeqTarg,"select SeqID from TblSeqs where ExprID = %d AND SeqName LIKE '%s'",ExprID,pszSeqName);
	sqlite3_exec(m_pDB,szSeqTarg,ExecCallbackID,&SeqID,NULL);
	}

if(SeqID == -1)	// will be -1 if not already in database so need to add
	{
	pStm = &m_StmSQL[2];								// access sequence statements
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}

	if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, pszSeqName,(int)strlen(pszSeqName)+1,SQLITE_STATIC))!=SQLITE_OK)
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
		sprintf(szSeqTarg,"select SeqID from TblSeqs where ExprID = %d AND SeqName LIKE '%s'",ExprID,pszSeqName);
		sqlite3_exec(m_pDB,szSeqTarg,ExecCallbackID,&SeqID,NULL);
		}
	else
		SeqID = (int)sqlite3_last_insert_rowid(m_pDB);
	m_NumSeqs += 1;						// number of seqs added to TblSeqs
	}

// replace lRA sequence
if(m_NumSeqMRA < cMaxMRASeqs)
	pMRASeq = &m_MRASeqs[m_NumSeqMRA++];
else
	{
	pMRASeq = m_MRASeqs;
	pLRASeq = pMRASeq++;
	for(Idx = 1; Idx < m_NumSeqMRA; Idx++, pMRASeq++)
		{
		if(pMRASeq->AccessCnt < pLRASeq->AccessCnt)
			pLRASeq = pMRASeq;
		}
	pMRASeq = pLRASeq;
	}
pMRASeq->AccessCnt = 1000;
pMRASeq->SeqID = SeqID;
strcpy(pMRASeq->szSeqName,pszSeqName);
return(SeqID);
}

int				// returned loci identifier 
CSQLiteMarkers::AddLoci(int ExprID,			// experiment identifier
		int SeqID,			// target assembly sequence identifier
		int Offset,			// offset on sequence			
		char Base)			// cannonical base at loci				
{
int sqlite_error;
tsStmSQL *pStm;
int LociID;
char szLoci[200];
char szBase[2];

szBase[0] = Base;			// SQLite seems to treat chars as 1byte integers and the command line SQLite shell displays as a numeric
szBase[1] = '\0';			// so easiest to use a VARCHAR(1) text string

if(m_pDB == NULL)
	return(eBSFerrInternal);

pStm = &m_StmSQL[3];								// access sequence statements
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, SeqID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, Offset))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, szBase,2,SQLITE_STATIC))!=SQLITE_OK)
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
	sprintf(szLoci,"select LociID from TblLoci where ExprID = %d AND SeqID = %d AND Offset = %d and Base = %d",ExprID,SeqID,Offset,(int)Base);
	sqlite3_exec(m_pDB,szLoci,ExecCallbackID,&LociID,NULL);
	}
else
	LociID = (int)sqlite3_last_insert_rowid(m_pDB);
m_NumSNPLoci += 1;					// number of SNP loci added to TblLoci
return(LociID);
}


int							// returned SNP indentifier
CSQLiteMarkers::AddSNP(int ExprID,			// experiment identifier
		int CultID,			// SNP in this cultivar relative to target sequence
		int LociID,			// identifies target sequence loci at which SNP has been identified
		int Acnt,			// count of A's in reads covering loci
		int Ccnt,			// count of C's in reads covering loci
		int Gcnt,			// count of G's in reads covering loci
		int Tcnt,			// count of T's in reads covering loci
		int Ncnt,			// count of N's in reads covering loci
		int TotCovCnt,		// total count of reads covering loci
		int TotMMCnt)		// of which there were this many reads with mismatches
{
int sqlite_error;
tsStmSQL *pStm;
int SnpID;
char szSNP[200];
pStm = &m_StmSQL[4];								// access sequence statements

if(m_pDB == NULL)
	return(eBSFerrInternal);

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, CultID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, LociID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 4, Acnt))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 5, Ccnt))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 6, Gcnt))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 7, Tcnt))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 8, Ncnt))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 9, TotCovCnt))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 10, TotMMCnt))!=SQLITE_OK)
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
	SnpID = -1;
	sprintf(szSNP,"select SnpID from TblSnps where ExprID = %d AND LociID = %d AND CultID = %d",ExprID,LociID,CultID);
	sqlite3_exec(m_pDB,szSNP,ExecCallbackID,&SnpID,NULL);
	}
else
	SnpID = (int)sqlite3_last_insert_rowid(m_pDB);
m_NumSNPs += 1;						// number of SNPs added to TblSnps
return(SnpID);
}

int						// returned marker identifier
CSQLiteMarkers::AddMarker(int ExprID,		// experiment identifier
		int CultID,			// marker in this cultivar relative to other cultivars
		int LociID,			// identifies target sequence loci at which marker has been identified
		char MarkerBase,		// marker base
		int MarkerScore)	// score
{
int sqlite_error;
tsStmSQL *pStm;
int MarkerID;
char szMarker[200];
char szBase[2];
szBase[0] = MarkerBase;		// SQLite seems to treat chars as 1byte integers and the command line SQLite shell displays as a numeric
szBase[1] = '\0';			// so easiest to use a VARCHAR(1) text string
pStm = &m_StmSQL[5];								// access sequence statements

if(m_pDB == NULL)
	return(eBSFerrInternal);

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, CultID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, LociID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, szBase,2,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 5, MarkerScore))!=SQLITE_OK)
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
	MarkerID = -1;
	sprintf(szMarker,"select MarkerID from TblMarkers where ExprID = %d AND LociID = %d AND CultID = %d AND Base = %d",ExprID,LociID,CultID,MarkerBase);
	sqlite3_exec(m_pDB,szMarker,ExecCallbackID,&MarkerID,NULL);
	}
else
	MarkerID = (int)sqlite3_last_insert_rowid(m_pDB);	
m_NumMarkers += 1;					// number of markers add to TblMarkers
return(MarkerID);
}

int
CSQLiteMarkers::AddMarkerSnp(int ExprID,
				int MarkerID,
				int SnpID)
{
int sqlite_error;
tsStmSQL *pStm;
int MarkerSnpID;
char szMarkerSnp[200];
pStm = &m_StmSQL[6];								// access sequence statements
if(m_pDB == NULL)
	return(eBSFerrInternal);
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, MarkerID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, SnpID))!=SQLITE_OK)
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
	MarkerSnpID = -1;
	sprintf(szMarkerSnp,"select MarkerSnpsID from TblMarkerSnps where ExprID = %d AND MarkerID = %d AND SnpID = %d",ExprID,MarkerID,SnpID);
	sqlite3_exec(m_pDB,szMarkerSnp,ExecCallbackID,&MarkerSnpID,NULL);
	}
else
	MarkerSnpID = (int)sqlite3_last_insert_rowid(m_pDB);	
return(MarkerSnpID);
}

int
CSQLiteMarkers::ProcessCSV2SQLite(int PMode,	// currently just the one mode...default is to parse from CSV and create/populate SQLite database
				  bool bSafe,					// if true then use indexing on all tables whilst inserting... much slower but perhaps safer if multiple threads ...
			      int CSVtype,					// input CSV file has this format (0: markers, 1: SNPs)
				  char *pszExprName,			// name by which this experiment is identified
				  char *pszExprDescr,			// describes experiment
				  char *pTargAssemb,			// assembly against which aligments for SNP discovery
				  int NumSpecies,				// number of species used in alignments
				  char *pszSpeciesNames[],		// names of species
				  char *pszInFile,				// parse from this input CSV file
				  char *pszDatabase)			// SQLite database file
{
int Rslt;
int ExprID;

int sqlite_error;
sqlite3_stmt *prepstatement = NULL;
tsCultivar *pCultivar;
int CultIdx;

sqlite3_initialize();

if((CreateDatabase(bSafe,pszDatabase))==NULL)
	{
	sqlite3_shutdown();
	return(eBSFerrInternal);
	}

if((Rslt = CreateExperiment(CSVtype,pszInFile,pszExprName,pszExprDescr,pTargAssemb)) < 1)
	{
	CloseDatabase(true);
	return(Rslt);
	}
ExprID = Rslt;

pCultivar = Cultivars;
for(CultIdx = 0; CultIdx < NumSpecies; CultIdx++, pCultivar++)
	{
	pCultivar->CultID = 0;
	if(CSVtype == 0)
		pCultivar->CultIdx = (4 + CultIdx * 8);
	else
		pCultivar->CultIdx = 0;
	strcpy(pCultivar->szCultivarName,pszSpeciesNames[CultIdx]);
	}

if((Rslt = CreateCultivars(NumSpecies,Cultivars)) < 1)
	{
	CloseDatabase(true);
	return(Rslt);
	}

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
CCSVFile *pCSV = new CCSVFile;
if(pCSV == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to instantiate CCSVfile");
	CloseDatabase(true);
	return(eBSFerrObj);
	}

if((Rslt=pCSV->Open(pszInFile))!=eBSFSuccess)
	{
	while(pCSV->NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,pCSV->GetErrMsg());
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open file: %s",pszInFile);
	delete pCSV;
	CloseDatabase(true);
	return(Rslt);
	}

int NumFields;
int NumElsRead;
int NumCultivars;
char *pszSeqName;
int SeqID;
NumElsRead = 0;
while((Rslt=pCSV->NextLine()) > 0)			// onto next line containing fields
	{
	if(!(NumElsRead % (bSafe ? 5000 : 100000)) && NumElsRead > 0)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d CSV lines - unique sequences: %d, SNP Loci: %d, SNPs: %d, Markers: %d",NumElsRead, m_NumSeqs,m_NumSNPLoci, m_NumSNPs, m_NumMarkers);

	NumFields = pCSV->GetCurFields();		// SNP files have 21, Marker CSV have 4 + (8 * NumCultivars) fields
	switch(CSVtype) {
		case 0:					// markers
			NumCultivars = (NumFields - 4)/8;
			if(NumFields < 12 || NumFields != (NumCultivars * 8) + 4)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected Marker CSV file number of fields to be ((NumCultivars * 8) + 4) in '%s', GetCurFields() returned '%d'",pszInFile,NumFields);
				delete pCSV;
				CloseDatabase(true);
				return(eBSFerrFieldCnt);
				}

			if(NumSpecies != NumCultivars)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected Marker CSV file to contain %d SNP species, NumCultivars in '%s' is %d",NumSpecies,pszInFile,NumCultivars);
				delete pCSV;
				CloseDatabase(true);
				return(eBSFerrFieldCnt);
				}

			if(!NumElsRead && pCSV->IsLikelyHeaderLine())
				continue;
			break;

		case 1:					// SNPs
			if(NumFields != 21)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Expected Marker CSV file number of fields to be 21 at line %d in '%s', GetCurFields() returned '%d'",NumElsRead,pszInFile,NumFields);
				delete pCSV;
				CloseDatabase(true);
				return(eBSFerrFieldCnt);
				}
			NumCultivars = 1;
			break;
		}

	if(!NumElsRead && pCSV->IsLikelyHeaderLine())
		continue;
	NumElsRead += 1;

	int Loci;
	int LociID;
	int SnpID;
	int CultID;
	char *pszLociBase;
	char *pszMarkerBase;
	int MarkerScore;
	int MarkerID;
	char LociBase;
	char MarkerBase;
	int TotCovCnt;
	int TotMMCnt;
	int SpeciesIdx;
	int Acnt;
	int Ccnt;
	int Gcnt;
	int Tcnt;
	int Ncnt;
	int FieldIdx;
	int NumCovd;
	
	switch(CSVtype) {
		case 0:					// markers
			pCSV->GetText(1,&pszSeqName);		// sequence containing marker
			RemoveQuotes(pszSeqName);
			SeqID = AddSeq(ExprID,pszSeqName);

			pCSV->GetInt(2,&Loci);		// loci on sequence at which marker has been determined
			pCSV->GetText(3,&pszLociBase); // cannonical target sequence base at the marker loci
			switch(*pszLociBase) {
				case 'a': case 'A':
					LociBase = 'A';
					break;
				case 'c': case 'C':
					LociBase = 'C';
					break;
				case 'g': case 'G':
					LociBase = 'G';
					break;
				case 't': case 'T':
					LociBase = 'T';
					break;
				default:
					LociBase = 'N';
					break;
				}
			LociID = AddLoci(ExprID,SeqID,Loci,LociBase);
			FieldIdx = 5;
			NumCovd = 0;
			for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++,FieldIdx += 8)
				{
				pCSV->GetText(FieldIdx,&pszMarkerBase);		// marker base
				switch(*pszMarkerBase) {
					case 'a': case 'A':
						MarkerBase = 'A';
						break;
					case 'c': case 'C':
						MarkerBase = 'C';
						break;
					case 'g': case 'G':
						MarkerBase = 'G';
						break;
					case 't': case 'T':
						MarkerBase = 'T';
						break;
					default:
						MarkerBase = 'N';
						break;
					}
				pCSV->GetInt(FieldIdx+1,&MarkerScore);		// cannonical target sequence base at the marker loci
				pCSV->GetInt(FieldIdx+2,&TotCovCnt);	
				pCSV->GetInt(FieldIdx+3,&Acnt);
				pCSV->GetInt(FieldIdx+4,&Ccnt);
				pCSV->GetInt(FieldIdx+5,&Gcnt);
				pCSV->GetInt(FieldIdx+6,&Tcnt);
				pCSV->GetInt(FieldIdx+7,&Ncnt);
				switch(LociBase) {
					case 'A':
						TotMMCnt = TotCovCnt - Acnt;
						break;
					case 'C':
						TotMMCnt = TotCovCnt - Ccnt;
						break;
					case 'G':
						TotMMCnt = TotCovCnt - Gcnt;
						break;
					case 'T':
						TotMMCnt = TotCovCnt - Tcnt;
						break;
					default:
						TotMMCnt = TotCovCnt - Ncnt;
						break;
					}
				CultID = (int)Cultivars[SpeciesIdx].CultID;
				if(TotCovCnt > 0 || MarkerScore > 0)
					{
					SnpID = AddSNP(ExprID,CultID,LociID,Acnt,Ccnt,Gcnt,Tcnt,Ncnt,TotCovCnt,TotMMCnt);
					Cultivars[SpeciesIdx].SnpID = SnpID;
					}
				else
					Cultivars[SpeciesIdx].SnpID = 0;

				if(MarkerScore > 0)
					{
					MarkerID = AddMarker(ExprID,CultID,LociID,MarkerBase,MarkerScore);
					Cultivars[SpeciesIdx].MarkerID = MarkerID;
					}
				else
					Cultivars[SpeciesIdx].MarkerID = 0;
				}
			for(SpeciesIdx = 0; SpeciesIdx < NumSpecies; SpeciesIdx++)
				{
				if(Cultivars[SpeciesIdx].MarkerID != 0)
					{
					MarkerID = Cultivars[SpeciesIdx].MarkerID;
					for(int Idx = 0; Idx < NumSpecies; Idx++)
						if(Cultivars[Idx].SnpID > 0)
							AddMarkerSnp(ExprID,MarkerID,Cultivars[Idx].SnpID);
					}
				}
			break;

		case 1:					// SNPs
			pCSV->GetText(4,&pszSeqName);		// sequence containing SNP
			RemoveQuotes(pszSeqName);
			SeqID = AddSeq(ExprID,pszSeqName);

			pCSV->GetInt(5,&Loci);		// loci on sequence at which marker has been determined
			pCSV->GetInt(11,&TotCovCnt);
			pCSV->GetInt(12,&TotMMCnt);
			pCSV->GetText(13,&pszLociBase); // cannonical target sequence base at the marker loci
			switch(*pszLociBase) {
				case 'a': case 'A':
					LociBase = 'A';
					break;
				case 'c': case 'C':
					LociBase = 'C';
					break;
				case 'g': case 'G':
					LociBase = 'G';
					break;
				case 't': case 'T':
					LociBase = 'T';
					break;
				default:
					LociBase = 'N';
					break;
				}
			LociID = AddLoci(ExprID,SeqID,Loci,LociBase);

			pCSV->GetInt(14,&Acnt);
			pCSV->GetInt(15,&Ccnt);
			pCSV->GetInt(16,&Gcnt);
			pCSV->GetInt(17,&Tcnt);
			pCSV->GetInt(18,&Ncnt);

			if(TotCovCnt) {
				switch(LociBase) {
					case 'A':
						if(Acnt == 0)
							Acnt = TotCovCnt - TotMMCnt;
						break;
					case 'C':
						if(Ccnt == 0)
							Ccnt = TotCovCnt - TotMMCnt;
						break;
					case 'G':
						if(Gcnt == 0)
							Gcnt = TotCovCnt - TotMMCnt;
						break;
					case 'T':
						if(Tcnt == 0)
							Tcnt = TotCovCnt - TotMMCnt;
						break;
					case 'N':
						if(Ncnt == 0)
							Ncnt = TotCovCnt - TotMMCnt;
						break;
					}
				CultID = (int)Cultivars[0].CultID;
				SnpID = AddSNP(ExprID,CultID,LociID,Acnt,Ccnt,Gcnt,Tcnt,Ncnt,TotCovCnt,TotMMCnt);
				}
			break;
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Parsed %d CSV lines - unique sequences: %d, SNP Loci: %d, SNPs: %d, Markers: %d",NumElsRead, m_NumSeqs,m_NumSNPLoci, m_NumSNPs, m_NumMarkers);

	// end transaction
if((sqlite_error = sqlite3_exec(m_pDB,pszEndTransaction,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't end transactions on '%s': %s", "Markers",sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed populating the sqlite database");

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating indexes ...");

tsStmSQL *pStms;
pStms = m_StmSQL;
int TblIdx;
for(TblIdx = 0; TblIdx < 7; TblIdx++,pStms++)
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
