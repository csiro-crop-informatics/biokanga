#pragma once
#include "./commdefs.h"

const unsigned int cBSGOAssocVersion = 10;		// increment each time the file strucure is changed
const unsigned int cBSGOAssocVersionBack = 10;	// can handle all versions starting from this minimum supported version

const int cAllocNumGenes = 50000;				// allocate tsGOGenes in this size chunks (number of tsGOGenes's) 
const int cAllocNumAssoc = cAllocNumGenes * 5;	// allocate tsGOAssoc in this size chunks (number of tsGOAssoc's) 
const int cAllocIdents = cAllocNumAssoc * 20;   // allocate for identifers in this size chunks (number of chars)
const int cMaxGOIDsPerGene = 2000;				// max number of GO IDs per gene

const int cMaxLenGO_DB = 20;					// max len GO Annotation File DB field
const int cMaxGOFieldLen = 256;					// max length of any other GO Annotation File DB field

const int cMaxMapFroms = 1000;					// max number of map from genes when parsing mapped gene file

const int cMinIsonameLen = 6;					// min element name length to be considered as being an isoform name

typedef enum eGOAssocParseType {
	eGOAPMUCSU = 0,								// parse file as UCSC 
	eGOAPGO,									// parse file as GO annotations
	eGOAPGOTAIR,								// parse file as TAIR GO annotations
    eGOAPGOFB									// parse file as GO Flybase annotations
} etGOAssocParseType;

typedef enum eGOaspect {
	eGOaspectUndefined = 0,
	eGOaspectP,									// biological process P
	eGOaspectF,									// molecular function F
	eGOaspectC									// cellular process C
} etGOaspect;

// GO association qualifiers (bit fields)
typedef enum eGOQual {
		eGOQnone = 0,				// no qualifiers specified
		eGOQNOT = 1,				//'NOT' is used to make an explicit note that the gene product is not associated with the GO term
		eGOQcontributes_to = 2,		//'Contributes_to' may be used only with molecular function terms
		eGOQcolocalizes_with = 4	//'Colocalizes_with' may be used only with cellular component terms
} etGOQual; 

// GO association object type
typedef enum eGOobjType {
	eGOTundefined = 0,
	eGOTgene, 
	eGOTtranscript, 
	eGOTprotein, 
	eGOTprotein_structure, 
	eGOTcomplex
}etGOobjType;

// GO association evidence (bit fields)
typedef enum eGOEvidence {
	eGOEnone = 0,
	eGOEIMP = 1,				//inferred from mutant phenotype
	eGOEIGI = 2,				//inferred from genetic interaction [with <database:gene_symbol[allele_symbol]>]
	eGOEIPI = 4,				//inferred from physical interaction [with <database:protein_name>]
	eGOEISS = 8,				//inferred from sequence similarity [with <database:sequence_id>]
	eGOEIDA = 16,				//inferred from direct assay
	eGOEIEP = 32,				//inferred from expression pattern
	eGOEIEA = 64,				//inferred from electronic annotation [with <database:id>]
	eGOETAS = 128,				//traceable author statement
	eGOENAS = 256,				//non-traceable author statement
	eGOEND = 512,				//no biological data available
	eGOEIC = 1024,				//inferred from reviewed computational analysis
	eGOERC = 2048				//inferred by curator [from <GO:id>]
}etGOEvidence;

#pragma pack(4)

typedef struct TAG_sGeneMap {
	char szFromGene[cMaxGeneNameLen];	// gene name to map from (sorted on szFromGene then szToGene)
	char szToGene[cMaxGeneNameLen];		// gene name to map to
	int nThInstance;					// nTh instance of szFromGene
} tsGeneMap;

typedef struct TAG_sGeneFilter {
	char szGene[cMaxGeneNameLen];		// retain this gene
} tsGeneFilter;

typedef struct TAG_sGeneGOitem {
	etGOobjType ObjType;				// object type
	char *pszGeneName;					// gene name
    char *pszQuals;						// list of GO qualifiers
	char *pszAspects;					// list of GO aspects
	char *pszEvidences;					// list of GO evidence
	char *pszGOIDs;						// list of GO:id
} tsGeneGOitem;

typedef union uIntPtr {
	INT32 Idx;
	INT64 Pad64;				// simply to ensure that ptr can be either 32 or 64bit so sizeof(tuGOIntPtr) is constant
	void *ptr;
	} tuIntPtr;

typedef struct TAG_sGOGene {
	tuIntPtr Name;					// gene name
	tGeneID GeneID;					// unique identifier for this gene or transcriptional loci
	INT32 NumGOAssoc;				// number of associated GO terms
	INT32 GOAssocIdx;				// index at which first GOAssoc (tsGOAssoc) for this gene starts
	UINT8 Type;						// etGOobjType
	} tsGOGene;

typedef struct TAG_sGOAssoc {
	tGeneID GeneID;					// associated with which gene
	INT32 IdentIdx;					// index at which associated GO:term identifier starts
	UINT16 Evidence;				// eGOEvidence
	UINT8 Qual;						// etGOQual
	UINT8 Aspect;					// etGOaspect
	} tsGOAssoc;
#pragma pack()

#pragma pack(8)
typedef struct TAG_sGOAssocFileHdr {
	unsigned char Magic[4];			// magic chars to identify this file as a biosequence file
	UINT64 FileLen;					// current file length
	INT64 GOGenesOfs;				// file offset to genes
	INT64 GOAssocsOfs;				// file offset to associations
	INT64 GeneIdentsOfs;			// file offset to gene identifiers
	INT64 GOIdentsOfs;				// file offset to GO identifiers

	UINT32 Type;					// GOAssoc file type 
	UINT32 Version;					// header version, incremented if structure changes with later releases
	UINT32 SizeOfHdr;				// total size of this header

	INT32 GOGenesCnt;				// number of genes
	INT32 GOGenesSize;				// size (bytes) genes

	INT32 GOAssocsCnt;				// number of associations
	INT32 GOAssocsSize;				// size (bytes) associations

	INT32 GeneIdentsCnt;			// number of gene identifiers
	INT32 GeneIdentsSize;			// size (bytes) gene identifiers
	
	INT32 GOIdentsCnt;				// number of GO identifiers
	INT32 GOIdentsSize;				// size (bytes) GO identifiers

	char szDescription[cMBSFFileDescrLen];// describes contents of file
	char szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc
}tsGOAssocFileHdr;

#pragma pack()

class CGOAssocs : protected CEndian, public CErrorCodes
{
	int m_hFile;				// opened/created file handle
	char m_szFile[_MAX_PATH+1];	// file name as opened/created
	bool m_bCreate;				// TRUE if file opened for create 
	bool m_bHdrDirty;			// TRUE if header needs to be written to file
	bool m_bGOAssocAvail;		// TRUE if GO associations are loaded and can be accessed

	tsGOAssocFileHdr m_FileHdr;

	bool m_bGeneNamesAsPtrs;	// true if the name union in tsGOGene instances ptd at by m_pGenes is ptr

	int m_GOIdentsCnt;			// number of GO identifiers
	int m_NxtGOIdentOfs;		// offset into m_pszGOIdents where to next write GO:term identifiers
	int m_MaxAllocGOIdents;		// size of memory currently allocated to m_pszGOIdents
	char *m_pszGOIdents;		// memory allocated to hold conactenated '\0' separated GO:term identifiers

	int m_GeneIdentsCnt;		// current number of gene name identifiers
	int m_NxtGeneIdentOfs;		// offset into m_pszGOIdents where to next write gene names
	int m_MaxAllocGeneIdents;	// size of memory currently allocated to m_pszGeneIdents
	char *m_pszGeneIdents;		// memory allocated to hold conactenated '\0' separated gene names

	int m_GOBIOGenesCnt;		// current number of genes associated with at least one BIO GO term
	int m_GOCELGenesCnt;		// current number of genes associated with at least one CEL GO term
	int m_GOMOLGenesCnt;		// current number of genes associated with at least one MOL GO term
	int m_GOGenesCnt;			// current number of genes associated with at least one GO term
	int m_MaxAllocGOGene;		// max number of tsGOGene that m_pGenes can hold
	tsGOGene *m_pGenes;			// memory allocated to hold tsGOGene structures

	int m_GOAssocCnt;			// current number of associations
	int m_MaxAllocGOAssoc;		// max number of tsGOAssoc that m_pAssocs can hold
	tsGOAssoc *m_pAssocs;		// memory allocated to hold association structures

	etGOEvidence m_DfltEvidence; // evidence to assume if none supplied when adding associations

	bool m_bDoGeneMaps;			// true if gene mapping required
	int m_NumGeneMaps;			// number of gene name mappings
	int m_AllocdGeneMaps;		// number of allocated gene name mappings
	tsGeneMap *m_pGeneMappings;	// pts to optional gene name mappings

	bool m_bDoGeneFilters;		// true if gene filtering required
	int m_NumGeneFilters;		// number of gene name filterings
	int m_AllocdGeneFilters;	// number of allocated gene name filterings
	tsGeneFilter *m_pGeneFilters; // pts to optional gene names to filter (retain)

	bool m_bGOFlybase;			// true if processing GO flybase, means need to map any BED filter genes (these include transcript suffix) onto GO cannonical gene names
	etGOAssocParseType m_FileParseType;  // expected association file type being parsed
	int m_NumGeneGOitems;		// number of gene GO items
	int m_AllocdGeneGOitems;	// number of allocated gene GO items
	tsGeneGOitem *m_pGeneGOitems; // pts to gene GO items
	int m_NumGeneGOitemsTxt;	// chrs used in m_pszGeneGoItemsTxt
	int m_AllocdGeneGoItemsTxt;	// current allocation
	char *m_pszGeneGoItemsTxt;	// allocated to hold all gene GO items text		

	teBSFrsltCodes Disk2Hdr(char *pszGoFile,int FileType);
	teBSFrsltCodes Hdr2Disk(void);

	void Reset(bool bFlush = false);
	void ClearAssociations(void);
	void InitHdr(void);
	teBSFrsltCodes Flush2Disk(void);
	
	
	int ParseGOAnnotation(etGOAssocParseType Type,  // expected file type
			  char *pszGOAnnotation,	// GO annotation file to parse
			  char *pszGeneMap = NULL,	// optional gene mapping file
			  char *pszGeneFilters = NULL); // optional gene filters file
	int ParseUCSC(char *pszUCSCGoAssoc,char *pszGeneFilters = NULL);			// UCSC file to parse
	teBSFrsltCodes LoadAssociations(void);
	teBSFrsltCodes ReadDisk(INT64 DiskOfs,int Len,void *pTo);
	int	ParseNxtField(char *pszLine,				// parse start
						 bool bReqValue,			// value is mandatory
						 int MaxLen,				// max allowed ret value
						 char *pszRetValue,			// where to return parsed value	void SwitchGeneNameIDX2PTR(void);
						 bool bAcceptComma=false);	// accept comma delimited as well as usual tab

	void SwitchGeneNameIDX2PTR(void);
	void SwitchGeneNamePTR2IDX(void);

	static int SortGeneGoItems( const void *arg1, const void *arg2);
	static int SortMapGenes( const void *arg1, const void *arg2);
	static int SortGeneFilters( const void *arg1, const void *arg2);
	static int SortGenes( const void *arg1, const void *arg2);
	static int SortDedupeGOIDs( const void *arg1, const void *arg2);
	tsGOGene *LocateGene(char *pszGene);
	tsGeneFilter  *LocateGeneFilter(const char *pszGene);
	tsGeneMap *LocateGeneMap(char *pszGene,int nThInstance=0);
	bool MapGene(char *pszGene,int nThInstance = 0);	// gene to be mapped

	int Add(etGOobjType ObjType,
   			   const char *pszGene,			// gene name
			   const char *pszQuals,			// list of GO qualifiers
			   const char *pszAspects,		// list of GO aspects
			   const char *pszEvidences,		// list of GO evidence
			   const char *pszGOIDs);			// list of GO:id

	int AddGeneGOlist(etGOobjType ObjType,
   			   const char *pszGene,			// gene name
			   const char *pszQuals,			// list of GO qualifiers
			   const char *pszAspects,		// list of GO aspects
			   const char *pszEvidences,		// list of GO evidence
			   const char *pszGOIDs);			// list of GO:id

	int ProvAddGeneGOlist(bool bDoGeneMaps,
							 char *pszGOID,
							 char *pszEvidence,
							 char *pszQualifier,
							 char *pszAspect,
							 char *pszObjectType,
							 char *pszObjectID,
							 char *pszObjectSymbol,
							 char *pszObjectName,
							 char *pszSynonyms);

	int ParseNxtSynonym(char *pszToParse,char *pszRetField);

	int ProcessGeneGOlist(void);		// process gene GO list, dedupe and combine
	int DeDupeGOIDs(char *pszQuals,char *pszAspects,char *pszEvidences,char *pszGOIDs); // remove duplicate GOIDs
	int GenGONumGenes(void);			// generate gene to GO counts for each root GO term

public:
	CGOAssocs(void);
	~CGOAssocs(void);
	teBSFrsltCodes Open(char *pszFileName,bool bCreate=false);

	int Parse(etGOAssocParseType Type,   // GO annotation file type
			  char *pszGOAnnotation,	// GO annotation file to parse
			  char *pszGeneMap = NULL,	// optional gene mapping file
			  char *pszGeneFilters = NULL); // optional gene filters file
	int ParseGeneMapping(char *pszGeneMapping);			// gene mapping file to parse
	int ParseGeneFilters(char *pszGeneFilters);			// gene filter file to parse
	teBSFrsltCodes Close(bool bFlush2Disk = false);

	int GetNumGeneGOitems(void);						// returns total number of gene GO items
	int GetNumGeneMaps(void);							// returns total number of gene mappings

	int GetNumGenes(void);					// returns total number of genes
	char *GetGene(int NthGene);				// returns Nth gene name

	int GetNumGOIDs(char *pszGeneName);				// returns number of GO identifiers associated with gene
	char *GetGOID(char *pszGeneName,int NthTerm);	// returns Nth GO identifier associated with gene or NULL if none
	int GetGOAttribs(char *pszGeneName,int NthTerm, etGOEvidence *pEvidence,etGOobjType *pType,etGOQual *pQual,etGOaspect *pAspect);

	int TrimWhitespace(char *pTxt);		// inplace trims any leading/trailing whitespace
	int StripQuotesWS(char *pszTxt);	// inplace strips any bracketing double or single quotes plus trims leading/trailing whitespace
	int ParseNthField(const char *pszToParse,char *pszRetField,int nthField);
	etGOQual MapGOQual(char *pszTxt);
	const char *MapGOaspect2Txt(etGOaspect GOaspect);
	etGOaspect MapGOaspect(char *pszTxt);
	const char *MapGOQual2Txt(etGOQual GOqual);
	etGOobjType MapGOobjType(char *pszTxt);
	const char *MapGOobjType2Txt(etGOobjType GOobjType);
	etGOEvidence MapGOEvidence(char *pszTxt);
	char *MapGOEvidence2Txt(etGOEvidence GOEvidence);

	int GetNumGOGenes(etOntologies OntologyClass);	// returns number of genes with at least one GO term (could be redundant) for specified Ontology
	bool TrimNameIso(char *pszName);				// Inplace remove any name suffix of the form '.[0-99]'
};

