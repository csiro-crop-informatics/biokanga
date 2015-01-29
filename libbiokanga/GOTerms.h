#pragma once

const unsigned int cBSGOTermsVersion = 11;		// increment each time the file strucure is changed
const unsigned int cBSGOTermsVersionBack = 11;	// can handle all versions starting from this minimum supported version
const int cMaxOBOLineLen = 2048;				// maximum expected length of any individual .obo file line before joining
const int cMaxOBOJoinedLineLen = cMaxOBOLineLen * 20;	// maximum expected length of any joined line when lines are joined via '\' indicator
const int cMaxOBOTagLen = 64;					// maximum expected length of any tag or stanza type name
const int cMaxOBOGOID = 32;						// maximum expected length of any GO:term identifier
const int cMaxOBOGOname = 128;					// maximum expected length of any GO:term name
const int cMaxOBOGOnamespace = 64;				// maximum expected length of any GO:term namespace
const int cMaxOBOGOdef = 2048;					// maximum expected length of any GO:term definition
const int cMaxOBOGOcomment = 2048;				// maximum expected length of any GO:term name
const int cMaxOBOGOdDbxref = 128;				// maximum expected length of any GO:term dbxref
const int cMaxOBOParentGOIDs = 100;				// allow any single GO:term on average to have upto this many parent terms
const int cMaxOBOAltGOIDs = 100;				// allow any single GO:term on average to have upto this many alternative IDs
const int cMaxOBOPartOfGOIDs = 100;				// allow any single GO:term on average to have upto this many part_of IDs
const int cAllocTerms = 15000;					// allocation size chunks (tsGOTerms) when reallocing m_pGOTerms
const int cAllocGOIDlists = 20000;				// allocation size chunks (ints) when reallocing alternative and parent GOIDofs
const int cAllocGOIDs = 0x0fffff;				// allocation size chunks (bytes) when reallocing m_pGOIDs
const int cAllocGOTagVals = 0x03fffff;			// allocation size chunks (bytes) when reallocing m_pGOTagVals

const unsigned int cClampCountVal = 0x7fffffff;	// counts are clamped to be <= this value

// tag definitions extracted from the following documents
// http://www.geneontology.org/GO.format.shtml (The OBO Flat File Format)
// http://www.godatabase.org/dev/doc/obo_format_spec.html (The OBO Flat File Format Specification, version 1.2)
// Note that there are slight differences in the tags specified by each document - I've gone for the union of all tags!



// stanza type enumerations
typedef enum eStanzaType {
	eTSTDocHdr,		// virtual stanza type for document header
	eTSTTerm,			// term stanza
	eTSTTypeDef,		// typedef stanza
	eTSTInstance		// instance stanza
	} teStanzaType;
typedef UINT8 tStanzaType;

// stanza tag enumerations
typedef enum eStanzaTags {
// tags in a [term] stanza
	eTSTid,				//The unique id of the current term.
	eTSTname,			//The term name
	eTSTis_anonymous,	//Whether or not the current object has an anonymous id
	eTSTalt_id,			//Defines an alternate id for this term	
	eTSTnamespace,		//The namespace in which the term belongs
	eTSTdef,			//The definition of the current term
	eTSTcomment,		//A comment for this term.
	eTSTsubset,			//This tag indicates a term subset to which this term belongs
	eTSTsynonym,		//This tag gives a synonym for this term, some xrefs to describe the origins of the synonym, and may indicate a synonym category or scope information.
	eTSTrelated_synonym,// Depreciated. Related alias for the synonym tag with the scope modifier set to RELATED.
	eTSTexact_synonym,	//Deprecated. An alias for the synonym tag with the scope modifier set to EXACT.
	eTSTnarrow_synonym,	//Deprecated. An alias for the synonym tag with the scope modifier set to NARROW.
	eTSTbroad_synonym,	//Deprecated. An alias for the synonym tag with the scope	modifier set to BROAD.
	eTSTxref,			//A dbxref that describes an analagous term in another vocabulary 	
	eTSTxref_analog,	//Deprecated. An alias for the xref tag.
	eTSTxref_unknown,	//Deprecated. An alias for the xref tag.
	eTSTis_a,			//This tag describes a subclassing relationship between one term and another
	eTSTintersection_of,//This tag indicates that this term represents the intersection of several other terms
	eTSTunion_of,		//This tag indicates that this term represents the union of several other terms
	eTSTdisjoint_from,	//This tag indicates that a term is disjoint from another, meaning that the two terms have no instances or subclasses in common.
	eTSTrelationship,	//This tag describes a typed relationship between this term and another term.
	eTSTis_obsolete,	//Whether or not this term is obsolete. Allowable values are	"true" and "false"
	eTSTreplaced_by,	//Gives a term which replaces an obsolete term.	
	eTSTconsider,		//Gives a term which may be an appropriate substitute for an obsolete term
	eTSTuse_term,		//Deprecated. Equivalent to consider.
	eTSTbuiltin,		//Whether or not this term or relation is builtin to the obo format
// tags in [typedef] stanza 
// also include all [term] stanza tags with exception of union_of, intersection_of and disjoint_from
	eTSTdomain,  			//The id of a term, or a special reserved identifier, which indicates the domain for this relationship type.
	eTSTrange, 				//The id of a term, or a special reserved identifier, which indicates acceptable range for this relationship type.
	eTSTinverse_of, 		//The id of another relationship type that is the inverse of this relationship type.
	eTSTtransitive_over,	//The id of another relationship type that this	relationship type is transitive over.
	eTSTis_cyclic,			//Whether or not a cycle can be made from this relationship type.
	eTSTis_reflexive,		//Whether this relationship is reflexive. 
	eTSTis_symmetric,		//Whether this relationship is symmetric
	eTSTis_anti_symmetric,	//Whether this relationship is anti-symmetric.
	eTSTis_transitive,		//Whether this relationship is transitive
	eTSTis_metadata_tag,	//Whether this relationship is a metadata tag. Properties that are marked as metadata tags are used to record object metadata. Object metadata is additional information about an object that is useful to track, but does not impact the definition of the object or how it should be treated by a reasoner. Metadata tags might be used to record special term synonyms or structured notes about a term, for example.
// tags in [instance] stanza
// also include [term] stanza tags id,name,is_anonymous,namespace,alt_id,comment,xref,synonym,is_obsolete,replaced_by and consider
	eTSTinstance_of,			//The term id that gives the class of which this is an instance.
	eTSTproperty_value,  		//This tag binds a property to a value in this instance
// tags in document header
	eTSTformat_version,		//Gives the obo specification version that this file uses
	eTSTtyperef,			//A url pointing to a type description document
	eTSTdefault_namespace,	//default namespace
	eTSTdata_version,		//Gives the version of the current ontology.
	eTSTversion,			//Deprecated. Use data-version instead.
	eTSTdate,				//The current date in dd:MM:yyyy HH:mm format.
	eTSTsaved_by,			//The username of the person to last save this file.
	eTSTauto_generated_by,	//The program that generated the file.
	eTSTsubsetdef,			//A description of a term subset	
	eTSTimport,				//A url pointing to another obo document to be appended to current obo document
	eTSTsynonymtypedef,		//A description of a user-defined synonym type
	eTSTidspace,			// 	A mapping between a "local" ID space and a "global" ID space.
	eTSTdefault_relationship_id_prefix,	// Any relationship lacking an ID space will be	prefixed with the value of this tag
	eTSTid_mapping,			// 	maps a Term or Typedef ID to another Term or Typedef ID
	eTSTremark,				// 	General comments for this file
// end of accepted tags in .obo files
	eTSTEndMarker			//used to flag the number of enums
} teStanzaTags;
typedef UINT8 tStanzaTag;

typedef enum eTSTValueTypes {
	eTSTValNone,				// no value to parse
	eTSTTermID,					// term identifier, e.g. GO:123457
	eTSTQuoteString,			// double quoted string e.g. "text string"
	eTSTWord,					// unquoted single word e.g. text_string
	eTSTDbxref,					// dbxref
	eTSTBool,					// boolean, e.g true, false, 1, 0
	eTSTText					// catchall general text which may be quoted/unquoted
} teTSTValueType;
typedef UINT8 tTSTValueType;

typedef enum eTSTError {
	eTSTErrNone,				// no error
	eTSTErrKeyTerm,				// A line does not contain a colon to indicate the end of the tag name	
	eTSTErrTagValue,			// There is no value to the right of the key tag colon
	eTSTErrEsc,					// An unrecognized escape sequence has been used
	eTSTErrEOL,					// More data was expected, but the input line ended
	eTSTErrQuote,				// A quoted string was expected, but something else was found
	eTSTErrQuoteEnd,			// A quoted string was found in an appropriate place, but no closing quote was encountered.
	eTSTErrDbxref,				// A dbxref list was expected, but something else was found.
	eTSTErrDbxrefMal,			// A dbxref list was found in an appropriate place, but it was not well formed
	eTSTErrDbxrefEnd,			// A dbxref list was found in an appropriate place, but no closing bracket was found
	eTSTErrGenErr				// catchall error
} etTSTError;
typedef UINT8 tTSTError;

#pragma pack(1)

// following is a fixed size union (32bit int, 32bit ptr or 64bit ptr) for either 32 or 64bit compiles
// as this union could be created on a 32bit machine, written to disk and processed on a 64bit machine
// thus size needs to be independent
typedef union uGOIntPtr {
	int Idx;
	INT64 Pad64;				// simply to ensure that ptr can be either 32 or 64bit so sizeof(tuGOIntPtr) is constant
	void *ptr;
	} tuGOIntPtr;

typedef struct TAG_sTSTerrMap {
	tTSTError	Error;
	const char *pszErr;
} tsTSTerrMap;

typedef struct TAG_sStanzaTag {
	tStanzaTag TagID;		// uniquely identifies this stanza tag
	int Min;			    // minimum number of instances of this tag allowed in stanza
	int Max;			    // maximum number (-1 == unlimited) of instances of this tag allowed in stanza
	const char *pszTagName;		// stanza tag name
	tTSTValueType Type1;	// expected value 1 type
	tTSTValueType Type2;	// expected value 2 type
} tsStanzaTag;

typedef struct TAG_sStanzaType {
	tStanzaType StanzaType;			// type of stanza
	const char *pszStanza;					// stanza type without bracketing '[' and ']'
	int NumTags;						// number of tags in this stanza
	tsStanzaTag *pTags;					// tags for this stanza
} tsStanzaType;

// GO element identifiers are maintained as a linear list
typedef struct TAG_sGOID {
	UINT16 Len;	// strlen(Txt)
	UINT16 Hash;		// hash over element text to quickly determine if element could match a probe element
	char Txt[1];			// element text placeholder, text to be '\0' terminated
} tsGOID;

// GO tag values are maintained as a linear list
typedef struct TAG_sGOTagVal {
	UINT16 Len;			// strlen(Txt)
	UINT16 Hash;		// hash over tag value text to quickly determine if text could match a probe text
	tStanzaTag Tag;		// tag
	char Val[1];		// tag value text placeholder, text to be '\0' terminated
} tsGOTagVal;

typedef struct TAG_sGOTerm {
	tuGOIntPtr TermGOID;		// GO:Ident
	UINT32 GOTermID;			// uniquely identifies this term instance (1..n)
	INT32 TermNameIdx;			// GO term name
	INT32 TermDefIdx;				// GO term definition
	INT32 ReplacedByIdx;			// GO term replacement GO:Ident
	INT32 GOSupTermIdx;			// GO:term superceding identifier
	INT32 GOAltTermsIdx;			// list of GO:term alternative identifiers
	INT32 GOParentIDsIdx;			// list of IsA GO:Idents
	INT32 GOPartOfIDsIdx;			// list of part_of GO:Idents (treated as if parents) 
	INT32 GOChildIDsIdx;			// list of child GO:Idents 
	INT32 GOTermCntsIdx;			// associated counts for this term or 0 if no counts

	UINT16 NumParents;  // number of parents to this term as a child (GOParentIDs) 
	UINT16 NumPartOfs; // number of part_of identifiers for this term (treated as if parents)
	UINT16 NumChildren; // number of children to this term as a parent (GOChildIDs)
	UINT16 NumAltTermIDs; // number of alternative identifiers for this term

	UINT8 RootOntology;	// to which root class (etOntologies) this term belongs 
	UINT8 bIsObsolete:1;     // if set then this term is obsolete and the ReplacedByIdx term should be used
	UINT8 bNoClass:1;		  // if set then term is classless - currently just bIsObsolete
	} tsGOTerm;

typedef struct TAG_sGOTermCnts {
	double HGProbK;					// probabilites from hypergeometric cdf
	double HGWeightedProbK;			// probabilites from hypergeometric cdf with gene weightings
	double HGMTCProbK;				// MTC probabilites from hypergeometric cdf
	double HGMTCWeightedProbK;		// MTC probabilites from hypergeometric cdf with gene weightings
	int TermID;						// to which term these counts apply
	int UpdateSeq;					// update sequence, used to track updates to counts
	unsigned int NumBkgdGenes;		// number of background genes contributing to this term
	unsigned int BkgdCnt;			// background counts 
	unsigned int NumSampleGenes;	// number of sample genes
	unsigned int SampleCnt;			// sample counts
	unsigned char RootOntology;		// to which root class (etOntologies) this term count belongs 
} tsGOTermCnts;
#pragma pack()

#pragma pack(8)
typedef struct TAG_sGOTermsFileHdr {
	unsigned char Magic[4];				// magic chars to identify this file as a biosequence file
	UINT64 FileLen;						// current file length
	INT64 GOTagValOfs;					// file offset to GO tag values
	INT64 ChildGOIDsOfs;				// file offset to GO child references
	INT64 ParentGOIDsOfs;				// file offset to GO parent references
	INT64 AltGOIDsOfs;					// file offset to GO alternative references
	INT64 GOTermOfs;					// file offset to GO terms
	INT64 GOIDOfs;						// file offset to GO identifiers

	INT32 Type;							// GOTerm file type 
	INT32 Version;						// header version, incremented if structure changes with later releases
	INT32 SizeOfHdr;					// total size of this header

	INT32 CellTermID;					// identifies root cellular component term
	INT32 BioTermID;					// identifies root biological process term 
	INT32 MolTermID;					// identifies root molecular function term

	INT32 POAnatomicalID;				// identifies plant plant anatomical entity term  (currently PO:0025131 as at June 2012)
	INT32 POdevID;						// identified plant structure growth term (currently PO:0009012  as at June 2012)

	INT32 GOIDcnt;						// number of GO identifiers
	INT32 GOIDSize;						// size (bytes) on disk of GO identifiers

	INT32 GOTagValcnt;					// number of GO tag values
	INT32 GOTagValSize;					// size (bytes) on disk of GO tag values

	INT32 ChildGOIDscnt;				// number of GO child references
	INT32 ChildGOIDsSize;				// size (bytes) on disk of GO child references

	INT32 ParentGOIDscnt;				// number of GO parent references
	INT32 ParentGOIDsSize;				// size (bytes) on disk of GO parent references

	INT32 AltGOIDscnt;					// number of GO alternative references
	INT32 AltGOIDsSize;					// size (bytes) on disk of GO alternative references

	INT32 GOCellTermcnt;				// number of terms classed as cellular
	INT32 GOBioTermcnt;					// number of terms classed as biological
	INT32 GOMolTermcnt;					// number of terms classed as molecular

	INT32 POAnatomyTermCnt;				// number of terms classed as plant anatomical
	INT32 PODevTermCnt;					// number of terms classed as plant structure developmental growth

	INT32 GOTermcnt;					// number of GO terms
	INT32 GOTermSize;					// size (bytes) on disk of GO terms

	char szDescription[cMBSFFileDescrLen];// describes contents of file
	char szTitle[cMBSFShortFileDescrLen];	// short title by which this file can be distingished from other files in dropdown lists etc
}tsGOTermsFileHdr;
#pragma pack()

class CGOTerms : protected CEndian, public CErrorCodes
{
	int m_hFile;						// opened/created file handle
	char m_szFile[_MAX_PATH+1];			// file name as opened/created
	bool m_bCreate;					    // TRUE if file opened for create 
	bool m_bHdrDirty;					// TRUE if header needs to be written to file
	bool m_bGOTermsAvail;				// TRUE if GO terms are loaded and can be accessed

	tsGOTermsFileHdr m_FileHdr;

	bool m_bTermGOIDsAsPtrs;			// TRUE if TermGOID currently set as ptrs instead of idx into m_pGOIDs (required when sorting)
	
	int m_GOCellTermcnt;				// current number of terms classed as cellular
	int m_GOBioTermcnt;					// current number of terms classed as biological
	int m_GOMolTermcnt;					// current number of terms classed as molecular

	int m_POAnatomyTermCnt;				// number of terms classed as plant anatomical
	int m_PODevTermCnt;					// number of terms classed as plant structure growth

	int m_GOTermcnt;					// current number of terms loaded into m_pGOTerms
	int m_AllocdGOTerms;				// total number elements (tsGOTerm) allocated into m_pGOTerms
	tsGOTerm *m_pGOTerms;				// memory allocated to hold array of GO terms

	int m_GOIDcnt;					// number of GOIDs currently in m_pGOIDs
	int m_NxtGOIDofs;				// offset into m_pGOIDs at which to write next GOID
	int m_AllocdGOIDs;				// total memory (bytes)  currently allocated for m_pGOIDs
	tsGOID *m_pGOIDs;				// memory allocated to hold list of GO identifiers

	int m_GOTagValcnt;				// number of tag + values currently in m_pGOTagVals
	int m_NxtGOTagValOfs;			// offset into m_pGOTagVals at which to write next tag + value
	int m_AllocdGOTagVals;			// total memory (bytes) currently allocated for m_pGOTagVals
	tsGOTagVal *m_pGOTagVals;		// memory (m_AllocdGOTagVals) allocated to hold tag + values

	int m_ChildGOIDscnt;			// currrent number of GOIDs in m_pGOChildIDs;
	int m_NxtGOChildIDIdx;			// index into m_pGOChildIDs of where to write next child GOIDofs
	int m_AllocdGOChildIDs;			// total memory (ints) currently allocated for m_pGOChildIDs
	INT32 *m_pGOChildIDs;			// array of term child GOIDofs

	int m_ParentGOIDscnt;			// currrent number of GOIDs in m_pParentGOIDs;
	int m_NxtGOParentIDIdx;			// index into m_pGOParentIDs of where to write next parent GOIDofs
	int m_AllocdGOParentIDs;		// total memory (ints) currently allocated for m_pGOParentIDs
	INT32 *m_pGOParentIDs;			// array of term parent GOIDofs

	int m_AltGOIDscnt;				// currrent number of GOIDs in m_pAltGOIDs;
	int m_NxtAltGOIDIdx;			// index into m_pAltGOIDs of where to write next alternative GOIDofs
	int m_AllocdAltGOIDs;			// total memory (ints) currently allocated for m_pAltGOIDs
	INT32 *m_pAltGOIDs;				// array of alternative GOIDofs

	int m_CurOBOLineNum;			// line number in current OBO file being processed

	int m_GOTermCntsCnt;			// currrent number of GOTermCnts in m_pGOTermCnts;
	int m_NxtGOTermCntsIdx;			// index into m_pGOTermCnts of where to write next tsGOTermCnts
	int m_AllocdGOTermCnts;			// total elements (tsGOTermCnts) currently allocated for m_pGOTermCnts
	tsGOTermCnts *m_pGOTermCnts;	// array of GOTermCnts

	int m_CurUpdateSeq;				// current update sequence counter - incremented each time a new AddCount() call is made to track propagation of counts into parent terms
	INT64 m_TotbBkgndLen;			// total of all backgound counts


	int	GetJoinedLine(char *pszRetLine,int MaxRetLineLen,FILE *pOBOstream);
	int TrimWhitespace(char *pTxt);		// inplace trims any leading/trailing whitespace
	int StripQuotesWS(char *pszTxt);	// inplace strips any bracketing double or single quotes plus trims leading/trailing whitespace
	tsStanzaTag *LocateStanzaTag(tsStanzaType *pCurStanza,char *pszTagName);

	bool Unescape(char *pszTxt);
	int	ParseTagName(char *pszTxt,int MaxValLen,char *pszTagName,bool bUnescape);
	int ParseTagValue(teTSTValueType Type,char *pszTxt,int MaxValLen,char *pszTagValue,bool bUnescape);
	int ParseTagDbxref(char *pszTxt,int MaxValLen,char *pszTagDbxref,bool bUnescape);

	int GetGOID(char *pszGOID);			// get byte offset into m_pGOIDs for tsGOElItem with matching GO:Term identifier
	int AddGOID(char *pszGOID);				// adds new tsGOElItem to m_pGOIDs or reuses existing if GO:Term identifier already in m_pGOIDs
	char *GetGOID(int GOIDOfs);				// returns pszGOID for GOIDOfs as returned from GetGOID() or AddGOID()

	int LocateGOTagVal(teStanzaTags Tag,char *pszVal);	// get byte offset into m_pGOTagVals for tsGOTagVal with matching tag + value
	int AddGOTagVal(teStanzaTags Tag,char *pszVal);		// adds new tsGOTagVal to m_pGOTagVals or reuses existing if tag + value already in m_pGOTagVals
	char *GetGOVal(int GOTagValOfs);// returns pszVal for offset as returned from AddGOTagVal() or LocateGOTagVal()
	teStanzaTags GetGOTag(int GOTagValOfs);// returns Tag for offset as returned from AddGOTagVal() or LocateGOTagVal()

	int AddTerm(char *pszGOID,char *pszTermName,char *pszTermNamespace,char *pszTermDef,char *pszTermDefDbxref,
			char *pszTermComment,bool bIsObsolete,char *pszGOIDreplacedby,
			int NumPartOfGOIDs,char *pszPartOfGOIDs,int NumGOIDAltIDs,char *pszGOIDAltIDs,int NumParentGOIDs,char *pszParentGOIDs);



	teBSFrsltCodes Disk2Hdr(char *pszGoFile,int FileType);
	teBSFrsltCodes Hdr2Disk(void);
	void Reset(bool bFlush = false);
	void ClearTerms(void);
	void InitHdr(void);
	teBSFrsltCodes Flush2Disk(void);
	teBSFrsltCodes LoadTerms(void);
	teBSFrsltCodes ReadDisk(INT64 DiskOfs,int Len,void *pTo);
	int GenBackRefs(void);
	
	void SwitchTermGOIDIDX2PTR(void);
	void SwitchTermGOIDPTR2IDX(void);

	static int SortTermsByGOID(const void *arg1, const void *arg2);
	static int SortTermCntsByTermID(const void *arg1, const void *arg2);
	static int SortTermCntsByHGProbK(const void *arg1, const void *arg2);
	static int SortTermCntsByHGWeightedProbK(const void *arg1, const void *arg2);
	static int SortTermCntsByHGMTCProbK(const void *arg1, const void *arg2);
	static int SortTermCntsByHGMTCWeightedProbK(const void *arg1, const void *arg2);

	int SetOntologyClass4Terms(void);
	void SetChildClass(tsGOTerm *pCurTerm,etOntologies OntologyClass); // NOTE: recursive function

	void ResetUpdateSeqs(void);			// reset all term.UpdateSeqs back to 0
	int AddCountRecurse(tsGOTerm *pTerm, // recurse into parents of this term
				   bool bSample,	// if true then update SampleCnts. if false then update background counts
				   int Count);		// count to increment by (1..n)

	unsigned int ClampCount(unsigned int Count,unsigned int Increment); // Returns Count + Increment, clamping to ensure <= cClampCountVal 
public:
	CGOTerms(void);
	~CGOTerms(void);
	teBSFrsltCodes Open(char *pszFileName,bool bCreate=false);
	teBSFrsltCodes Close(bool bFlush2Disk = false);
	int Parse(FILE *pOBOstream);
	int NumGOTerms(etOntologies RootOntology); // returns total number of GO:Terms for each root ontology
	int ClearStats(void);				 // clears all term stats
	int ClearSampleCounts(void);		 // clears all sample counts/stats
	int AddCount(etOntologies OntologyClass, // which class of ontologies to add count to
				bool bProp,			 // if true then propagate counts into parents
				   bool bSample,		 // if true then update SampleCnts. if false then update background counts 
				   int Count,		     // count to increment by (1..n)
				   int NumGOIDs,	     // number of GO:Term identifiers in pszGOID[]
				   char *pszGOID[]);     // array of ptrs to GO:Term identifiers

	int SetBkgndCnts(etOntologies OntologyClass, // which class of ontologies to set counts for
		bool bProp,						// propagate counts from GO:Term into parent terms
		bool bBkgndLen,					// background counts proportional to gene lengths
		int Strand,						// background counts are for which strand 0==both,1=='+' and 2=='-'
		CBEDfile *pBED,					// gene BED file
		CGOAssocs *pAssocs);			// gene to ID association file
		
	tsGOTerm *LocateGOID(char *pszGOID); // LocateGOID searches for specified GOID and returns a ptr to term identified by specified GOID or NULL if unable to locate term
	tsGOTerm *LocateNthParent(tsGOTerm *pCurTerm, int NthParent); // Returns a ptr to GO:Term which is the Nth (1..n) parent relative to the specified term
	tsGOTerm *LocateNthChild(tsGOTerm *pCurTerm, int NthChild);   // Returns a ptr to GO:Term which is the Nth (1..n) child relative to the specified term
	tsGOTerm *LocateRoot(tsGOTerm *pCurTerm);					  // Returns a ptr to root GO:Term of specified term
	tsGOTermCnts *GetTermCnts(tsGOTerm *pTerm);					  // returns ptr to counts associated with specified term, allocates new counts if none previously associated
	tsGOTermCnts *GetExistingTermCnts(tsGOTerm *pTerm);	  // returns ptr to counts associated with specified term or NULL if no counts associated
	tsGOTermCnts *GetTermCnts(int Ith); // returns the Ith, 1..GetNumTermCnts(), term counts
	
	int GetNumTermCnts(etOntologies RootOntology);			// returns number of term counts which have been associated with terms
	int GetNumSampledTermCnts(etOntologies RootOntology,unsigned int MinGenes=1); // returns number of term counts which are associated back to sample MinGenes	
	int GetNumBkgndTermCnts(etOntologies RootOntology,unsigned int MinGenes=1); // returns number of term counts which are associated back to population MinGenes	
	int SortCntsTermID(void);			// sorts counts by their TermID
	int SortCntsHGProbK(void);			// sorts counts by HGProbK
	int SortCntsHGWeightedProbK(void);	// sorts counts by HGWeightedProbK
	int SortCntsHGMTCProbK(void);		// sorts counts by HGMTCProbK
	int SortCntsHGMTCWeightedProbK(void);// sorts counts by HGMTCWeightedProbK

	tsGOTerm *GetGOTermByID(int TermID); // returns ptr to term identified by TermID
	char *GetGOIDforTerm(tsGOTerm *pTerm); // returns GO:Term identifier for specified term
	char *GetGONameforTerm(tsGOTerm *pTerm); // returns GO:Term name for specified term

	char *RootOntology2Txt(etOntologies Ontology); // returns root class ontology term
	int GetRootGOTermID(etOntologies Ontology);	   // returns TermID for root class ontology term

};
