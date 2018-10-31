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

// tags allowed in [Term] stanzas
tsStanzaTag TermStanzaTags[] = {
	{eTSTid,1,1,			"id",eTSTTermID,eTSTValNone},	//The unique id of the current term.
	{eTSTname,1,1,			"name",eTSTText,eTSTValNone},	//The term name
	{eTSTis_anonymous,0,1,	"is_anonymous",eTSTBool,eTSTValNone},	//Whether or not the current object has an anonymous id
	{eTSTalt_id,0,-1,		"alt_id",eTSTTermID,eTSTValNone},		//Defines an alternate id for this term	
	{eTSTnamespace,0,1,		"namespace",eTSTWord,eTSTValNone},	//The namespace in which the term belongs
	{eTSTdef,0,1,			"def",eTSTQuoteString,eTSTDbxref},	//The definition of the current term
	{eTSTcomment,0,1,		"comment",eTSTText,eTSTValNone},		//A comment for this term.
	{eTSTsubset,0,-1,		"subset",eTSTWord,eTSTValNone},		//This tag indicates a term subset to which this term belongs
	{eTSTsynonym,0,-1,		"synonym",eTSTQuoteString,eTSTDbxref},	//This tag gives a synonym for this term, some xrefs to describe the origins of the synonym, and may indicate a synonym category or scope information.
	{eTSTrelated_synonym,0,-1,"related_synonym",eTSTQuoteString,eTSTDbxref},// Depreciated. Related alias for the synonym tag with the scope modifier set to RELATED.
	{eTSTexact_synonym,0,-1,"exact_synonym",eTSTQuoteString,eTSTDbxref},	//Deprecated. An alias for the synonym tag with the scope	modifier set to EXACT.
	{eTSTnarrow_synonym,0,-1,"narrow_synonym",eTSTQuoteString,eTSTDbxref},	//Deprecated. An alias for the synonym tag with the scope modifier set to NARROW.
	{eTSTbroad_synonym,0,-1,"broad_synonym",eTSTQuoteString,eTSTDbxref},	//Deprecated. An alias for the synonym tag with the scope	modifier set to BROAD.
	{eTSTxref,0,-1,			"xref",eTSTText,eTSTValNone},			//A dbxref that describes an analagous term in another vocabulary 	
	{eTSTxref_analog,0,-1,	"xref_analog",eTSTText,eTSTValNone},		//Deprecated. An alias for the xref tag.
	{eTSTxref_unknown,0,-1,	"xref_unknown",eTSTText,eTSTValNone},	//Deprecated. An alias for the xref tag.
	{eTSTis_a,0,-1,			"is_a",eTSTTermID,eTSTValNone},			//This tag describes a subclassing relationship between one term and another
	{eTSTintersection_of,0,-1,"intersection_of",eTSTText,eTSTValNone},//This tag indicates that this term represents the intersection of several other terms
	{eTSTunion_of,0,-1,		"union_of",eTSTText,eTSTValNone},		//This tag indicates that this term represents the union of several other terms
	{eTSTdisjoint_from,0,-1,"disjoint_from",eTSTText,eTSTValNone},	//This tag indicates that a term is disjoint from another, meaning that the two terms have no instances or subclasses in common.
	{eTSTrelationship,0,-1,	"relationship",eTSTText,eTSTValNone},	//This tag describes a typed relationship between this term and another term.
	{eTSTis_obsolete,0,1,	"is_obsolete",eTSTBool,eTSTValNone},	//Whether or not this term is obsolete. Allowable values are	"true" and "false"
	{eTSTreplaced_by,0,1,	"replaced_by",eTSTTermID,eTSTValNone},	//Gives a term which replaces an obsolete term.	
	{eTSTconsider,0,-1,		"consider",eTSTTermID,eTSTValNone},		//Gives a term which may be an appropriate substitute for an obsolete term
	{eTSTuse_term,0,-1,		"use_term",eTSTTermID,eTSTValNone},		//Deprecated. Equivalent to consider.
	{eTSTbuiltin,0,1,		"builtin",eTSTBool,eTSTValNone}			//Whether or not this term or relation is builtin to the obo format
};
const int cNumTermStanzaTags = sizeof(TermStanzaTags)/sizeof(TermStanzaTags[1]);

// tags allowed in [typedef] stanzas
tsStanzaTag TypedefStanzaTags[] = {
	{eTSTid,1,1,			"id",eTSTTermID,eTSTValNone},			//The unique id of the current term.
	{eTSTname,1,1,			"name",eTSTText,eTSTValNone},			//The term name
	{eTSTis_anonymous,0,1,	"is_anonymous",eTSTBool,eTSTValNone},	//Whether or not the current object has an anonymous id
	{eTSTalt_id,0,-1,		"alt_id",eTSTTermID,eTSTValNone},			//Defines an alternate id for this term	
	{eTSTnamespace,0,1,		"namespace",eTSTWord,eTSTValNone},		//The namespace in which the term belongs
	{eTSTdef,0,1,			"def",eTSTQuoteString,eTSTDbxref},			//The definition of the current term
	{eTSTcomment,0,1,		"comment",eTSTText,eTSTDbxref},		//A comment for this term.
	{eTSTsubset,0,-1,		"subset",eTSTWord,eTSTValNone},			//This tag indicates a term subset to which this term belongs
	{eTSTsynonym,0,-1,		"synonym",eTSTQuoteString,eTSTDbxref},		//This tag gives a synonym for this term, some xrefs to describe the origins of the synonym, and may indicate a synonym category or scope information.
	{eTSTrelated_synonym,0,-1,"related_synonym",eTSTQuoteString,eTSTDbxref},// Depreciated. Related alias for the synonym tag with the scope modifier set to RELATED.
	{eTSTexact_synonym,0,-1,"exact_synonym",eTSTQuoteString,eTSTDbxref},	//Deprecated. An alias for the synonym tag with the scope	modifier set to EXACT.
	{eTSTnarrow_synonym,0,-1,"narrow_synonym",eTSTQuoteString,eTSTDbxref},	//Deprecated. An alias for the synonym tag with the scope modifier set to NARROW.
	{eTSTbroad_synonym,0,-1,"broad_synonym",eTSTQuoteString,eTSTDbxref},	//Deprecated. An alias for the synonym tag with the scope	modifier set to BROAD.
	{eTSTxref,0,-1,			"xref",eTSTText,eTSTValNone},			//A dbxref that describes an analagous term in another vocabulary 	
	{eTSTxref_analog,0,-1,	"xref_analog",eTSTText,eTSTValNone},	//Deprecated. An alias for the xref tag.
	{eTSTxref_unknown,0,-1,	"xref_unknown",eTSTText,eTSTValNone},	//Deprecated. An alias for the xref tag.
	{eTSTis_a,0,-1,			"is_a",eTSTTermID,eTSTValNone},			//This tag describes a subclassing relationship between one term and another
	{eTSTrelationship,0,-1,	"relationship",eTSTText,eTSTValNone},	//This tag describes a typed relationship between this term and another term.
	{eTSTis_obsolete,0,1,	"is_obsolete",eTSTBool,eTSTValNone},	//Whether or not this term is obsolete. Allowable values are	"true" and "false"
	{eTSTreplaced_by,0,1,	"replaced_by",eTSTTermID,eTSTValNone},	//Gives a term which replaces an obsolete term.	
	{eTSTconsider,0,-1,		"consider",eTSTTermID,eTSTValNone},		//Gives a term which may be an appropriate substitute for an obsolete term
	{eTSTuse_term,0,-1,		"use_term",eTSTTermID,eTSTValNone},		//Deprecated. Equivalent to consider.
	{eTSTbuiltin,0,1,		"builtin",eTSTBool,eTSTValNone},		//Whether or not this term or relation is builtin to the obo format
	{eTSTdomain,0,-1,  		"domain",eTSTTermID,eTSTValNone},		//The id of a term, or a special reserved identifier, which indicates the domain for this relationship type.
	{eTSTrange, 0,-1,		"range",eTSTTermID,eTSTValNone},		//The id of a term, or a special reserved identifier, which indicates acceptable range for this relationship type.
	{eTSTinverse_of, 0,-1,	"inverse_of",eTSTTermID,eTSTValNone},	//The id of another relationship type that is the inverse of this relationship type.
	{eTSTtransitive_over,0,-1,"transitive_over",eTSTTermID,eTSTValNone},	//The id of another relationship type that this	relationship type is transitive over.
	{eTSTis_cyclic,0,1,	"is_cyclic",eTSTBool,eTSTValNone},			//Whether or not a cycle can be made from this relationship type.
	{eTSTis_reflexive,0,1,	"is_reflexive",eTSTBool,eTSTValNone},	//Whether this relationship is reflexive. 
	{eTSTis_symmetric,0,1,	"is_symmetric",eTSTBool,eTSTValNone},	//Whether this relationship is symmetric
	{eTSTis_anti_symmetric,0,1,"is_anti_symmetric",eTSTBool,eTSTValNone},	// Whether this relationship is anti-symmetric.
	{eTSTis_transitive,0,1,	"is_transitive",eTSTBool,eTSTValNone},	// Whether this relationship is transitive
	{eTSTis_metadata_tag,0,1,"is_metadata_tag",eTSTBool,eTSTValNone}	// Whether this relationship is a metadata tag. Properties that are marked as metadata tags are used to record object metadata. Object metadata is additional information about an object that is useful to track, but does not impact the definition of the object or how it should be treated by a reasoner. Metadata tags might be used to record special term synonyms or structured notes about a term, for example.
};
const int cNumTypedefStanzaTags = sizeof(TypedefStanzaTags)/sizeof(TypedefStanzaTags[1]);
 
// tags allowed in [instance] stanzas
tsStanzaTag InstanceStanzaTags[] = {
	{eTSTid,1,1,			"id",eTSTTermID,eTSTValNone},			//The unique id of the current term.
	{eTSTname,1,1,			"name",eTSTText,eTSTValNone},			//The term name
	{eTSTinstance_of,1,1,	"instance_of",eTSTTermID,eTSTValNone},	// The term id that gives the class of which this is an instance.
	{eTSTis_anonymous,0,1,	"is_anonymous",eTSTBool,eTSTValNone},	// Whether or not the current object has an anonymous id
	{eTSTalt_id,0,-1,		"alt_id",eTSTTermID,eTSTValNone},		// Defines an alternate id for this term	
	{eTSTcomment,0,1,		"comment",eTSTText,eTSTValNone},		// A comment for this term.
	{eTSTxref,0,-1,			"xref",eTSTText,eTSTValNone},			// A dbxref that describes an analagous term in another vocabulary 	
	{eTSTsynonym,0,-1,		"synonym",eTSTQuoteString,eTSTDbxref},	// This tag gives a synonym for this term, some xrefs to describe the origins of the synonym, and may indicate a synonym category or scope information.
	{eTSTis_obsolete,0,1,	"is_obsolete",eTSTBool,eTSTValNone},	// Whether or not this term is obsolete. Allowable values are	"true" and "false"
	{eTSTreplaced_by,0,1,	"replaced_by",eTSTTermID,eTSTValNone},	// Gives a term which replaces an obsolete term.	
	{eTSTconsider,0,-1,		"consider",eTSTTermID,eTSTValNone},		// Gives a term which may be an appropriate substitute for an obsolete term
	{eTSTproperty_value,0,-1,"property_value",eTSTText,eTSTValNone} // This tag binds a property to a value in this instance
};
const int cNumInstanceStanzaTags = sizeof(InstanceStanzaTags)/sizeof(InstanceStanzaTags[1]);


// tags allowed in document headers
tsStanzaTag DocHdrTags[] = {
	{eTSTformat_version,1,1,"format-version",eTSTText,eTSTValNone},	//Gives the obo specification version that this file uses
	{eTSTtyperef,1,1,"typeref",eTSTText,eTSTValNone},				//A url pointing to a type description document
	{eTSTdefault_namespace,0,1,"default-namespace",eTSTWord,eTSTValNone},	//default namespace
	{eTSTdata_version,0,1,"data-version",eTSTText,eTSTValNone},		//Gives the version of the current ontology.
	{eTSTversion,0,1,"version",eTSTText,eTSTValNone},				//Deprecated. Use data-version instead.
	{eTSTdate,0,1,"date",eTSTText,eTSTValNone},						//The current date in dd:MM:yyyy HH:mm format.
	{eTSTsaved_by,0,1,"saved-by",eTSTText,eTSTValNone},				//The username of the person to last save this file.
	{eTSTauto_generated_by,0,1,"auto-generated-by",eTSTText,eTSTValNone},	//The program that generated the file.
	{eTSTsubsetdef,0,-1,"subsetdef",eTSTWord,eTSTQuoteString},			//A description of a term subset	
	{eTSTimport,0,-1,"import",eTSTText,eTSTValNone},					//A url pointing to another obo document to be appended to current obo document
	{eTSTsynonymtypedef,0,-1,"synonymtypedef",eTSTText,eTSTValNone},	//A description of a user-defined synonym type
	{eTSTidspace,0,-1,"idspace",eTSTText,eTSTValNone},				// 	A mapping between a "local" ID space and a "global" ID space.
	{eTSTdefault_relationship_id_prefix,0,1,"default-relationship-id-prefix"},	// Any relationship lacking an ID space will be	prefixed with the value of this tag
	{eTSTid_mapping,0,-1,"id-mapping",eTSTText,eTSTValNone},			// 	maps a Term or Typedef ID to another Term or Typedef ID
	{eTSTremark,0,1,"remark",eTSTText,eTSTValNone}					// 	General comments for this file
 };
const int cNumDocHdrStanzaTags = sizeof(DocHdrTags)/sizeof(DocHdrTags[1]);


tsStanzaType StanzaTypes[] = {
	{eTSTDocHdr,"DocHeader",cNumDocHdrStanzaTags,DocHdrTags},		// virtual stanza type for document header
	{eTSTTerm,"Term",cNumTermStanzaTags,TermStanzaTags},			// term stanza
	{eTSTTypeDef,"TypeDef",cNumTypedefStanzaTags,TypedefStanzaTags},	// typedef stanza
	{eTSTInstance,"Instance",cNumInstanceStanzaTags,InstanceStanzaTags}		// instance stanza
	};
const int cNumDocStanzaTypes = sizeof(StanzaTypes)/sizeof(StanzaTypes[1]);

tsTSTerrMap TSTerrMaps[] = {
	{eTSTErrNone,		"No Errors"},	
	{eTSTErrKeyTerm,	"Cannot find key-terminating colon"},  
	{eTSTErrTagValue,	"Tag has no value"}, 
	{eTSTErrEsc,		"Unrecognized escape character"}, 
	{eTSTErrEOL,		"Unexpected end of line"}, 
	{eTSTErrQuote,		"Expected quoted string"}, 
	{eTSTErrQuoteEnd,	"Unclosed quoted string"}, 	
	{eTSTErrDbxref,		"Expected dbxref list"}, 	
	{eTSTErrDbxrefMal,	"Malformed dbxref list"}, 
	{eTSTErrDbxrefEnd,	"Unclosed dbxref list"},
	{eTSTErrGenErr,		"General Error"}
};
const int cNumParseErrs = sizeof(TSTerrMaps)/sizeof(TSTerrMaps[1]);

// root terms
const char *pszRootTerms[] = { "cellular_component","biological_process","molecular_function", "plant anatomical entity","plant structure development stage"};

void
ReportParseError(etTSTError Error,int LineNum,const char *pszSupInfo)
{
const char *pszErrTxt;
int Idx;
for(Idx = 0; Idx < cNumParseErrs; Idx++)
	if(TSTerrMaps[Idx].Error == Error)
		break;
if(Idx < cNumParseErrs)
	pszErrTxt = TSTerrMaps[Idx].pszErr;
else
	pszErrTxt = "General Error";
	
printf("\nError: %s at line: %d %s",pszErrTxt,LineNum,pszSupInfo);
}

 CGOTerms::CGOTerms(void)
{
m_hFile = -1;
m_pGOTerms = NULL;
m_pGOIDs = NULL;
m_pGOTagVals = NULL;
m_pGOParentIDs = NULL;
m_pGOChildIDs = NULL;
m_pAltGOIDs =NULL;
m_pGOTermCnts = NULL;	
Reset();
}

CGOTerms::~CGOTerms(void)
{
Reset();
}


void 
CGOTerms::Reset(bool bFlush)
{
if(bFlush)
	Flush2Disk();

if(m_hFile != -1)
	{
	close(m_hFile);
	m_hFile = -1;
	}

ClearTerms();

m_szFile[0] = '\0';			// file name as opened/created
m_bCreate=false;		    // TRUE if file opened for create 
m_bHdrDirty=false;			// TRUE if header needs to be written to file
m_bGOTermsAvail=false;		// TRUE if GO terms are loaded and can be accessed
m_bTermGOIDsAsPtrs=false;	// TRUE if TermGOID currently set as ptrs instead of idx into m_pGOIDs (required when sorting) 
m_CurOBOLineNum=0;			// line number in current OBO file being processed
InitHdr();
}

teBSFrsltCodes
CGOTerms::Disk2Hdr(char *pszGoFile,int FileType)
{
if(_lseeki64(m_hFile,0,SEEK_SET)!=0)			// read in header..
	{
	AddErrMsg("CGOTerms::Open","Seek failed to offset 0 on %s - %s",pszGoFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsGOTermsFileHdr) != read(m_hFile,&m_FileHdr,sizeof(tsGOTermsFileHdr)))
	{
	AddErrMsg("CGOTerms::Open","Read of file header failed on %s - %s",pszGoFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

	// header read, validate it as being a GOTerms file header
if(tolower(m_FileHdr.Magic[0]) != 'b' || 
	tolower(m_FileHdr.Magic[1]) != 'i' || 
	tolower(m_FileHdr.Magic[2]) != 'o' || 
	tolower(m_FileHdr.Magic[3]) != 's')
	{
	AddErrMsg("CGOTerms::Open","%s opened but no magic signature - not a GOTerms file",pszGoFile);
	Reset(false);			// closes opened file..
	return(eBSFerrNotBioseq);
	}



if(m_bIsBigEndian)	// file was written with little-endian ordering
	{
	m_FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);						// current file length
	m_FileHdr.GOTagValOfs = SwapUI64Endians(m_FileHdr.GOTagValOfs);					// file offset to GO tag values
	m_FileHdr.ChildGOIDsOfs = SwapUI64Endians(m_FileHdr.ChildGOIDsOfs);				// file offset to GO child references
	m_FileHdr.ParentGOIDsOfs = SwapUI64Endians(m_FileHdr.ParentGOIDsOfs);				// file offset to GO parent references
	m_FileHdr.AltGOIDsOfs = SwapUI64Endians(m_FileHdr.AltGOIDsOfs);					// file offset to GO alternative references
	m_FileHdr.GOTermOfs = SwapUI64Endians(m_FileHdr.GOTermOfs);					// file offset to GO terms

	m_FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);							// GOTerm file type 
	m_FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);						// header version, incremented if structure changes with later releases
	m_FileHdr.SizeOfHdr = SwapUI32Endians(m_FileHdr.SizeOfHdr);					// total size of this header

	m_FileHdr.CellTermID = SwapUI32Endians(m_FileHdr.CellTermID);				// identifies root cellular component term
	m_FileHdr.BioTermID = SwapUI32Endians(m_FileHdr.BioTermID);					// identifies root biological process term 
	m_FileHdr.MolTermID = SwapUI32Endians(m_FileHdr.MolTermID);					// identifies root molecular function term

	m_FileHdr.POAnatomicalID = SwapUI32Endians(m_FileHdr.POAnatomicalID);				// identifies plant anotomical term  (currently PO:0025131 as at June 2012) 
	m_FileHdr.POdevID = SwapUI32Endians(m_FileHdr.POdevID);						// identifies plant structure growth term (currently PO:0009012  as at June 2012)

	m_FileHdr.GOIDcnt = SwapUI32Endians(m_FileHdr.GOIDcnt);						// number of GO identifiers
	m_FileHdr.GOIDSize = SwapUI32Endians(m_FileHdr.GOIDSize);						// size (bytes) on disk of GO identifiers
	m_FileHdr.GOIDOfs = SwapUI64Endians(m_FileHdr.GOIDOfs);						// file offset to GO identifiers

	m_FileHdr.GOTagValcnt = SwapUI32Endians(m_FileHdr.GOTagValcnt);					// number of GO tag values
	m_FileHdr.GOTagValSize = SwapUI32Endians(m_FileHdr.GOTagValSize);					// size (bytes) on disk of GO tag values

	m_FileHdr.ChildGOIDscnt = SwapUI32Endians(m_FileHdr.ChildGOIDscnt);				// number of GO child references
	m_FileHdr.ChildGOIDsSize = SwapUI32Endians(m_FileHdr.ChildGOIDsSize);				// size (bytes) on disk of GO child references

	m_FileHdr.ParentGOIDscnt = SwapUI32Endians(m_FileHdr.ParentGOIDscnt);				// number of GO parent references
	m_FileHdr.ParentGOIDsSize = SwapUI32Endians(m_FileHdr.ParentGOIDsSize);				// size (bytes) on disk of GO parent references

	m_FileHdr.AltGOIDscnt = SwapUI32Endians(m_FileHdr.AltGOIDscnt);					// number of GO alternative references
	m_FileHdr.AltGOIDsSize = SwapUI32Endians(m_FileHdr.AltGOIDsSize);					// size (bytes) on disk of GO alternative references

	m_FileHdr.GOCellTermcnt = SwapUI32Endians(m_FileHdr.GOCellTermcnt);				// number of terms classed as cellular
	m_FileHdr.GOBioTermcnt = SwapUI32Endians(m_FileHdr.GOBioTermcnt);				// number of terms classed as biological
	m_FileHdr.GOMolTermcnt = SwapUI32Endians(m_FileHdr.GOMolTermcnt);				// number of terms classed as molecular

	m_FileHdr.POAnatomyTermCnt = SwapUI32Endians(m_FileHdr.POAnatomyTermCnt);			// number of terms classed as plant structural
	m_FileHdr.PODevTermCnt = SwapUI32Endians(m_FileHdr.PODevTermCnt);				// number of terms classed as plant structural developmental

	m_FileHdr.GOTermcnt = SwapUI32Endians(m_FileHdr.GOTermcnt);					// number of GO terms
	m_FileHdr.GOTermSize = SwapUI32Endians(m_FileHdr.GOTermSize);					// size (bytes) on disk of GO terms
	}

	// check GOTerms file is the type we are expecting
if(m_FileHdr.Type != FileType)
	{
	AddErrMsg("CGOTerms::Open","%s opened as a GOTerms file - expected type %d, file type is %d",pszGoFile,FileType,m_FileHdr.Type);
	Reset(false);			// closes opened file..
	return(eBSFerrFileType);
	}

	// can we handle this version?
if(m_FileHdr.Version > cBSGOTermsVersion || m_FileHdr.Version < cBSGOTermsVersionBack)
	{
	AddErrMsg("CBEDfile::Open","%s opened as a GOTerms file - can only handle versions %d to %d, file version is %d",pszGoFile,
			cBSGOTermsVersionBack,cBSGOTermsVersion,m_FileHdr.Version);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}
return(eBSFSuccess);
}

teBSFrsltCodes
CGOTerms::Hdr2Disk(void)
{
tsGOTermsFileHdr FileHdr;
tsGOTermsFileHdr *pHdr;
int WrtLen;

WrtLen = sizeof(tsGOTermsFileHdr);

if(m_bIsBigEndian)	// if on a big-endian machine then need to make little endian as that is our native file format
	{
	memmove(&FileHdr,&m_FileHdr,WrtLen);
	FileHdr.FileLen = SwapUI64Endians(m_FileHdr.FileLen);						// current file length
	FileHdr.GOTagValOfs = SwapUI64Endians(m_FileHdr.GOTagValOfs);					// file offset to GO tag values
	FileHdr.ChildGOIDsOfs = SwapUI64Endians(m_FileHdr.ChildGOIDsOfs);				// file offset to GO child references
	FileHdr.ParentGOIDsOfs = SwapUI64Endians(m_FileHdr.ParentGOIDsOfs);				// file offset to GO parent references
	FileHdr.AltGOIDsOfs = SwapUI64Endians(m_FileHdr.AltGOIDsOfs);					// file offset to GO alternative references
	FileHdr.GOTermOfs = SwapUI64Endians(m_FileHdr.GOTermOfs);					// file offset to GO terms

	FileHdr.Type = SwapUI32Endians(m_FileHdr.Type);							// GOTerm file type 
	FileHdr.Version = SwapUI32Endians(m_FileHdr.Version);						// header version, incremented if structure changes with later releases
	FileHdr.SizeOfHdr = SwapUI32Endians(m_FileHdr.SizeOfHdr);					// total size of this header

	FileHdr.CellTermID = SwapUI32Endians(m_FileHdr.CellTermID);					// identifies root cellular component term
	FileHdr.BioTermID = SwapUI32Endians(m_FileHdr.BioTermID);				// identifies root biological process term 
	FileHdr.MolTermID = SwapUI32Endians(m_FileHdr.MolTermID);				// identifies root molecular function term

	FileHdr.POAnatomicalID = SwapUI32Endians(m_FileHdr.POAnatomicalID);				// identifies plant structure term  (currently PO:0025131 as at June 2012) 
	FileHdr.POdevID = SwapUI32Endians(m_FileHdr.POdevID);					// identified plant structure growth term (currently PO:0009012  as at June 2012)


	FileHdr.GOIDcnt = SwapUI32Endians(m_FileHdr.GOIDcnt);						// number of GO identifiers
	FileHdr.GOIDSize = SwapUI32Endians(m_FileHdr.GOIDSize);						// size (bytes) on disk of GO identifiers
	FileHdr.GOIDOfs = SwapUI64Endians(m_FileHdr.GOIDOfs);						// file offset to GO identifiers

	FileHdr.GOTagValcnt = SwapUI32Endians(m_FileHdr.GOTagValcnt);					// number of GO tag values
	FileHdr.GOTagValSize = SwapUI32Endians(m_FileHdr.GOTagValSize);					// size (bytes) on disk of GO tag values

	FileHdr.ChildGOIDscnt = SwapUI32Endians(m_FileHdr.ChildGOIDscnt);				// number of GO child references
	FileHdr.ChildGOIDsSize = SwapUI32Endians(m_FileHdr.ChildGOIDsSize);				// size (bytes) on disk of GO child references

	FileHdr.ParentGOIDscnt = SwapUI32Endians(m_FileHdr.ParentGOIDscnt);				// number of GO parent references
	FileHdr.ParentGOIDsSize = SwapUI32Endians(m_FileHdr.ParentGOIDsSize);				// size (bytes) on disk of GO parent references

	FileHdr.AltGOIDscnt = SwapUI32Endians(m_FileHdr.AltGOIDscnt);					// number of GO alternative references
	FileHdr.AltGOIDsSize = SwapUI32Endians(m_FileHdr.AltGOIDsSize);					// size (bytes) on disk of GO alternative references

	FileHdr.GOCellTermcnt = SwapUI32Endians(m_FileHdr.GOCellTermcnt);				// number of terms classed as cellular
	FileHdr.GOBioTermcnt = SwapUI32Endians(m_FileHdr.GOBioTermcnt);					// number of terms classed as biological
	FileHdr.GOMolTermcnt = SwapUI32Endians(m_FileHdr.GOMolTermcnt);					// number of terms classed as molecular

	m_FileHdr.POAnatomyTermCnt = SwapUI32Endians(m_FileHdr.POAnatomyTermCnt);			// number of terms classed as plant anatomical entity
	m_FileHdr.PODevTermCnt = SwapUI32Endians(m_FileHdr.PODevTermCnt);				// number of terms classed as plant structural developmental


	FileHdr.GOTermcnt = SwapUI32Endians(m_FileHdr.GOTermcnt);					// number of GO terms
	FileHdr.GOTermSize = SwapUI32Endians(m_FileHdr.GOTermSize);					// size (bytes) on disk of GO terms
	pHdr = &FileHdr;
	}
else
	pHdr = &m_FileHdr;

if(_lseeki64(m_hFile,0,SEEK_SET) ||
			write(m_hFile,pHdr,WrtLen)!=WrtLen)
	{
	AddErrMsg("CBEDfile::Flush2Disk","Unable to write file header to disk on file %s - error %s",m_szFile,strerror(errno));
	Reset(false);
	return(eBSFerrFileAccess);
	}
m_bHdrDirty = false;
return(eBSFSuccess);
}



void
CGOTerms::ClearTerms(void)
{
if(m_pGOTerms != NULL)
	{
	delete m_pGOTerms;
	m_pGOTerms = NULL;
	}
if(m_pGOTagVals != NULL)
	{
	delete m_pGOTagVals;
	m_pGOTagVals = NULL;
	}
if(m_pGOIDs != NULL)
	{
	delete m_pGOIDs;
	m_pGOIDs = NULL;
	}
if(m_pGOChildIDs!=NULL)
	{
	delete m_pGOChildIDs;
	m_pGOChildIDs = NULL;
	}
if(m_pGOParentIDs != NULL)
	{
	delete m_pGOParentIDs;
	m_pGOParentIDs = NULL;
	}
if(m_pAltGOIDs != NULL)
	{
	delete m_pAltGOIDs;
	m_pAltGOIDs = NULL;
	}
if(m_pGOTermCnts != NULL)
	{	
	delete m_pGOTermCnts;
	m_pGOTermCnts = NULL;
	}

m_CurUpdateSeq=0;
m_GOCellTermcnt=0;				// current number of terms classed as cellular
m_GOBioTermcnt=0;					// current number of terms classed as biological
m_GOMolTermcnt=0;					// current number of terms classed as molecular
m_GOTermcnt=0;				// current number of terms loaded into m_pGOTerms
m_AllocdGOTerms=0;			// number of terms allocated into m_pGOTerms
m_GOIDcnt=0;				// number of GOIDs currently in m_pGOIDs
m_NxtGOIDofs=0;				// offset into m_pGOIDs at which to write next GOID
m_AllocdGOIDs=0;			// total memory (bytes)  currently allocated for m_pGOIDs
m_GOTagValcnt=0;			// number of tag + values currently in m_pGOTagVals
m_NxtGOTagValOfs=0;			// offset into m_pGOTagVals at which to write next tag + value
m_AllocdGOTagVals=0;		// total memory (bytes) currently allocated for m_pGOTagVals
m_ChildGOIDscnt=0;			// currrent number of GOIDs in m_pGOChildIDs;
m_NxtGOChildIDIdx=0;		// index into m_pGOChildIDs of where to write next child GOIDofs
m_AllocdGOChildIDs=0;		// total memory (ints) currently allocated for m_pGOChildIDs
m_ParentGOIDscnt=0;			// currrent number of GOIDs in m_pParentGOIDs;
m_AllocdGOParentIDs=0;		// total memory (ints) currently allocated for m_pGOParentIDs
m_AltGOIDscnt=0;			// currrent number of GOIDs in m_pAltGOIDs;
m_AllocdAltGOIDs=0;			// total memory (ints) currently allocated for m_pAltGOIDs
m_GOCellTermcnt=0;			// current number of terms classed as cellular
m_GOBioTermcnt=0;			// current number of terms classed as biological
m_GOMolTermcnt=0;			// current number of terms classed as molecular

m_POAnatomyTermCnt = 0;		// number of terms classed as plant structure
m_PODevTermCnt = 0;			// number of terms classed as plant structure growth

m_GOTermCntsCnt=0;			// currrent number of GOTermCnts in m_pGOTermCnts
m_NxtGOTermCntsIdx=0;		// index into m_pGOTermCnts of where to write next tsGOTermCnts
m_AllocdGOTermCnts=0;		// total elements (tsGOTermCnts) currently allocated for m_pGOTermCnts
}

// InitHdr
// Initiates the file header as containing no GO terms or assocated relationships
void 
CGOTerms::InitHdr(void)
{
memset(&m_FileHdr,0,sizeof(m_FileHdr));
m_FileHdr.Magic[0] = 'b';
m_FileHdr.Magic[1] = 'i';
m_FileHdr.Magic[2] = 'o';
m_FileHdr.Magic[3] = 's';
m_FileHdr.Type = cBSFTypeGOTerms;		// biosequence file type 
m_FileHdr.Version = cBSGOTermsVersion;	// header version, incremented if structure changes with later releases
m_FileHdr.FileLen = sizeof(m_FileHdr);	// current file length
m_FileHdr.SizeOfHdr = sizeof(m_FileHdr);// total size of this header
m_FileHdr.szDescription[0] = '\0';
m_FileHdr.szTitle[0] = '\0';
m_bHdrDirty = true;						// TRUE if header needs to be written to file
}
	

teBSFrsltCodes 
CGOTerms::Open(char *pszBioGO,bool bCreate)
{
teBSFrsltCodes Rslt;
if(pszBioGO == NULL || *pszBioGO == '\0') // validate parameters
	{
	AddErrMsg("CGOTerms::Open","Parameter errors - %s",pszBioGO);
	return(eBSFerrParams);
	}


Reset(false);						// reset context in case file was previously opened

#ifdef _WIN32
if(!bCreate)
	m_hFile = open(pszBioGO, O_READSEQ ); // file access is normally sequential..
else
	m_hFile = open(pszBioGO, O_CREATETRUNC);
#else
if(!bCreate)
	m_hFile = open64(pszBioGO, O_READSEQ); // file access is normally sequential..
else
	{
     if((m_hFile = open64(pszBioGO,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
        if(ftruncate(m_hFile,0)!=0)
			{
			AddErrMsg("CGOTerms::Open","Unable to truncate %s - %s",pszBioGO,strerror(errno));
			Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
			return(Rslt);
			}
	}
#endif

if(m_hFile == -1)					// check if file open succeeded
	{
	AddErrMsg("CGOTerms::Open","Unable to open %s - %s",pszBioGO,strerror(errno));
	Rslt = bCreate ? eBSFerrCreateFile : eBSFerrOpnFile;
	return(Rslt);
	}

strncpy(m_szFile,pszBioGO,_MAX_PATH);
m_szFile[_MAX_PATH-1] = '\0';
if(bCreate)
	{
	m_bCreate = true;
	m_bGOTermsAvail = false;
	InitHdr();
	if((Rslt = Flush2Disk()) != eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(Rslt);
		}
	}
else // else opening existing file
	{
	if((Rslt=Disk2Hdr(pszBioGO,cBSFTypeGOTerms))!=eBSFSuccess)
		{
		Reset(false);			// closes opened file..
		return(Rslt);
		}

	// if not empty then load terms + parent/child relationships...
	if(m_FileHdr.GOTermcnt)
		{
		if((Rslt=LoadTerms())!=eBSFSuccess)	// close file after terms loaded
			{
			AddErrMsg("CGOTerms::Open","Error loading GO:Terms from %s",pszBioGO);
			Reset(false);
			return(Rslt);
			}
		}
	m_bGOTermsAvail = true;		// terms are in memory and can now be accessed
	}
return(eBSFSuccess);
}

// LoadTerms
// Load GO:terms from previously opened file
// Any already loaded GO:terms are deleted
teBSFrsltCodes
CGOTerms::LoadTerms(void)
{
teBSFrsltCodes Rslt;
int Idx;
if(m_hFile == -1)
	return(eBSFerrClosed);


ClearTerms();

m_bTermGOIDsAsPtrs=false;	// TRUE if TermGOID currently set as ptrs instead of idx into m_pGOIDs (required when sorting) 

if(!m_FileHdr.GOTermcnt)	// any terms in file to load?
	return(eBSFerrNoFeatures);

// allocate all memory up front before loading from file
if(m_FileHdr.GOIDcnt)
	{
	if((m_pGOIDs = (tsGOID *)new unsigned char [m_FileHdr.GOIDSize])==NULL)
		return(eBSFerrMem);
	}
m_GOIDcnt =	m_FileHdr.GOIDcnt;
m_NxtGOIDofs = m_FileHdr.GOIDSize;
m_AllocdGOIDs = m_FileHdr.GOIDSize;

if(m_FileHdr.GOTagValcnt)
	{
	if((m_pGOTagVals = (tsGOTagVal *)new unsigned char [m_FileHdr.GOTagValSize])==NULL)
		{
		ClearTerms();
		return(eBSFerrMem);
		}
	}
m_GOTagValcnt = m_FileHdr.GOTagValcnt;
m_NxtGOTagValOfs = m_FileHdr.GOTagValSize;
m_AllocdGOTagVals = m_FileHdr.GOTagValSize;
if(m_FileHdr.ChildGOIDscnt)
	{
	if((m_pGOChildIDs = (int *)new unsigned char [m_FileHdr.ChildGOIDsSize])==NULL)
		{
		ClearTerms();
		return(eBSFerrMem);
		}
	}
m_ChildGOIDscnt = m_FileHdr.ChildGOIDscnt;
m_NxtGOChildIDIdx = m_FileHdr.ChildGOIDscnt;
m_AllocdGOChildIDs = m_FileHdr.ChildGOIDscnt;


if(m_FileHdr.ParentGOIDscnt)
	{
	if((m_pGOParentIDs = (int *)new unsigned char [m_FileHdr.ParentGOIDsSize])==NULL)
		{
		ClearTerms();
		return(eBSFerrMem);
		}
	}
m_ParentGOIDscnt = m_FileHdr.ParentGOIDscnt;
m_NxtGOParentIDIdx = m_FileHdr.ParentGOIDscnt;
m_AllocdGOParentIDs = m_FileHdr.ParentGOIDscnt;


if(m_FileHdr.AltGOIDscnt)
	{
	if((m_pAltGOIDs = (int *)new unsigned char [m_FileHdr.AltGOIDsSize])==NULL)
		{
		ClearTerms();
		return(eBSFerrMem);
		}
	}
m_AltGOIDscnt = m_FileHdr.AltGOIDscnt;
m_NxtAltGOIDIdx = m_FileHdr.AltGOIDscnt;
m_AllocdAltGOIDs = m_FileHdr.AltGOIDscnt;

if(m_FileHdr.GOTermcnt)
	{
	if((m_pGOTerms = (tsGOTerm *)new unsigned char [m_FileHdr.GOTermSize])==NULL)
		{
		ClearTerms();
		return(eBSFerrMem);
		}
	}

m_GOCellTermcnt= m_FileHdr.GOCellTermcnt;
m_GOBioTermcnt= m_FileHdr.GOBioTermcnt;	
m_GOMolTermcnt= m_FileHdr.GOMolTermcnt;	
m_POAnatomyTermCnt = m_FileHdr.POAnatomyTermCnt;
m_PODevTermCnt = m_FileHdr.PODevTermCnt;


m_GOTermcnt = m_FileHdr.GOTermcnt;
m_AllocdGOTerms = m_FileHdr.GOTermcnt;

// all required memory has been allocated, now read from disk
if(m_pGOIDs != NULL)
	{
	Rslt = ReadDisk(m_FileHdr.GOIDOfs,m_FileHdr.GOIDSize,m_pGOIDs);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		tsGOID *pGOID = m_pGOIDs;
		UINT8 *pByte = (UINT8 *)pGOID;
		for(Idx = 0; Idx < m_GOIDcnt; Idx++)
			{
			pByte += pGOID->Len + sizeof(tsGOID);
			pGOID->Len = SwapUI32Endians(pGOID->Len);
			pGOID = (tsGOID *)pByte;
			}
		}
	}

if(Rslt == eBSFSuccess && m_pGOTagVals != NULL)
	{
	Rslt = ReadDisk(m_FileHdr.GOTagValOfs,m_FileHdr.GOTagValSize,m_pGOTagVals);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		tsGOTagVal *pGOTagVal = m_pGOTagVals;
		UINT8 *pByte = (UINT8 *)pGOTagVal;
		for(Idx = 0; Idx < m_GOTagValcnt; Idx++)
			{
			pByte += pGOTagVal->Len + sizeof(tsGOTagVal);
			pGOTagVal->Len = SwapUI32Endians(pGOTagVal->Len);
			pGOTagVal = (tsGOTagVal *)pByte;
			}
		}
	}

if(Rslt == eBSFSuccess && m_pGOChildIDs != NULL)
	{
	Rslt = ReadDisk(m_FileHdr.ChildGOIDsOfs,m_FileHdr.ChildGOIDsSize,m_pGOChildIDs);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		INT32 *pGOChildID = m_pGOChildIDs;
		for(Idx = 0; Idx < m_ChildGOIDscnt; Idx++,pGOChildID++)
			*pGOChildID = SwapUI32Endians(*pGOChildID);
		}
	}

if(Rslt == eBSFSuccess && m_pGOParentIDs != NULL)
	{
	Rslt = ReadDisk(m_FileHdr.ParentGOIDsOfs,m_FileHdr.ParentGOIDsSize,m_pGOParentIDs);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		INT32 *pGOParentID = m_pGOParentIDs;
		for(Idx = 0; Idx < m_ParentGOIDscnt; Idx++,pGOParentID++)
			*pGOParentID = SwapUI32Endians(*pGOParentID);
		}
	}

if(Rslt == eBSFSuccess && m_pAltGOIDs != NULL)
	{
	Rslt = ReadDisk(m_FileHdr.AltGOIDsOfs,m_FileHdr.AltGOIDsSize,m_pAltGOIDs);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		INT32 *pAltGOIDs = m_pAltGOIDs;
		for(Idx = 0; Idx < m_AltGOIDscnt; Idx++,pAltGOIDs++)
			*pAltGOIDs = SwapUI32Endians(*pAltGOIDs);
		}
	}

if(Rslt == eBSFSuccess && m_pGOTerms != NULL)
	{
	Rslt = ReadDisk(m_FileHdr.GOTermOfs,m_FileHdr.GOTermSize,m_pGOTerms);
	if(Rslt == eBSFSuccess && m_bIsBigEndian)
		{
		tsGOTerm *pGOTerm = m_pGOTerms;
		for(Idx = 0; Idx < m_GOTermcnt; Idx++,pGOTerm++)
			{
			pGOTerm->TermGOID.Pad64 = SwapUI64Endians(pGOTerm->TermGOID.Pad64);		// GO:Ident
			pGOTerm->GOTermID = SwapUI32Endians(pGOTerm->GOTermID);			// uniquely identifies this term instance (1..n)
			pGOTerm->TermNameIdx = SwapUI32Endians(pGOTerm->TermNameIdx);			// GO term name
			pGOTerm->TermDefIdx = SwapUI32Endians(pGOTerm->TermDefIdx);				// GO term definition
			pGOTerm->ReplacedByIdx = SwapUI32Endians(pGOTerm->ReplacedByIdx);			// GO term replacement GO:Ident
			pGOTerm->GOSupTermIdx = SwapUI32Endians(pGOTerm->GOSupTermIdx);			// GO:term superceding identifier
			pGOTerm->GOAltTermsIdx = SwapUI32Endians(pGOTerm->GOAltTermsIdx);			// list of GO:term alternative identifiers
			pGOTerm->GOParentIDsIdx = SwapUI32Endians(pGOTerm->GOParentIDsIdx);			// list of IsA GO:Idents
			pGOTerm->GOPartOfIDsIdx = SwapUI32Endians(pGOTerm->GOPartOfIDsIdx);			// list of part_of GO:Idents (treated as if parents) 
			pGOTerm->GOChildIDsIdx = SwapUI32Endians(pGOTerm->GOChildIDsIdx);			// list of child GO:Idents 
			pGOTerm->GOTermCntsIdx = SwapUI32Endians(pGOTerm->GOTermCntsIdx);			// associated counts for this term or 0 if no counts

			pGOTerm->NumParents = SwapUI16Endians(pGOTerm->NumParents);  // number of parents to this term as a child (GOParentIDs) 
			pGOTerm->NumPartOfs = SwapUI16Endians(pGOTerm->NumPartOfs); // number of part_of identifiers for this term (treated as if parents)
			pGOTerm->NumChildren = SwapUI16Endians(pGOTerm->NumChildren); // number of children to this term as a parent (GOChildIDs)
			pGOTerm->NumAltTermIDs = SwapUI16Endians(pGOTerm->NumAltTermIDs); // number of alternative identifiers for this term
			}
		}
	}

if(Rslt != eBSFSuccess)
	ClearTerms();
return(Rslt);
}

// ReadDisk
// Reads block of size 'Len' from disk starting at 'DiskOfs' into preallocated memory at 'pTo'
teBSFrsltCodes
CGOTerms::ReadDisk(INT64 DiskOfs,int Len,void *pTo)
{
if(_lseeki64(m_hFile,DiskOfs,SEEK_SET)!=DiskOfs)
	{
	AddErrMsg("CGOTerms::ReadDisk","Seek failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
if(read(m_hFile,pTo,Len)!=Len)
	{
	AddErrMsg("CGOTerms::ReadDisk","Read failed on %s - %s",m_szFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}
return(eBSFSuccess);
}

teBSFrsltCodes 
CGOTerms::Close(bool bFlush2Disk)
{
if(bFlush2Disk) {
	if(m_GOTermcnt > 1 && m_pGOTerms != NULL)
		{
		SwitchTermGOIDIDX2PTR();
		qsort(m_pGOTerms,m_GOTermcnt,sizeof(tsGOTerm),SortTermsByGOID);
		SwitchTermGOIDPTR2IDX();
		}
	Flush2Disk();
	}
Reset();
return(eBSFSuccess);
}

teBSFrsltCodes 
CGOTerms::Flush2Disk(void)
{
int Rslt;
int WrtLen;
int Idx;

if(m_hFile != -1 && m_bCreate)		// if file opened for write because creating..
	{
	if(m_GOTermcnt)					// are there terms to write?
		{
		// ensure terms are sorted by identifier
		SwitchTermGOIDIDX2PTR();
		qsort(m_pGOTerms,m_GOTermcnt,sizeof(tsGOTerm),SortTermsByGOID);
		SwitchTermGOIDPTR2IDX();

		// write out m_GOIDcnt GO identifiers
		m_FileHdr.GOIDcnt = m_GOIDcnt;
		m_FileHdr.GOIDSize = m_NxtGOIDofs;
		m_FileHdr.GOIDOfs = m_FileHdr.FileLen;

		if(m_bIsBigEndian)
			{
			tsGOID *pGOID = m_pGOIDs;
			UINT8 *pByte = (UINT8 *)pGOID;
			for(Idx = 0; Idx < m_GOIDcnt; Idx++)
				{
				pByte += pGOID->Len + sizeof(tsGOID);
				pGOID->Len = SwapUI32Endians(pGOID->Len);
				pGOID = (tsGOID *)pByte;
				}
			}

		WrtLen = m_FileHdr.GOIDSize;	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,(char *)m_pGOIDs,WrtLen)!=WrtLen)
			{
			AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO identifiers to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FileLen += WrtLen;

		// now write out the GO term tag values
		m_FileHdr.GOTagValcnt = m_GOTagValcnt;
		m_FileHdr.GOTagValSize = m_NxtGOTagValOfs;
		m_FileHdr.GOTagValOfs = m_FileHdr.FileLen;


		if(m_bIsBigEndian)
			{
			tsGOTagVal *pGOTagVal = m_pGOTagVals;
			UINT8 *pByte = (UINT8 *)pGOTagVal;
			for(Idx = 0; Idx < m_GOTagValcnt; Idx++)
				{
				pByte += pGOTagVal->Len + sizeof(tsGOTagVal);
				pGOTagVal->Len = SwapUI32Endians(pGOTagVal->Len);
				pGOTagVal = (tsGOTagVal *)pByte;
				}
			}

		WrtLen = m_FileHdr.GOTagValSize ;	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,(char *)m_pGOTagVals,WrtLen)!=WrtLen)
			{
			AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO Tag values to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FileLen += WrtLen;


		// now write out the GO term child references
		m_FileHdr.ChildGOIDscnt = m_ChildGOIDscnt;
		m_FileHdr.ChildGOIDsSize = m_ChildGOIDscnt * sizeof(int);
		m_FileHdr.ChildGOIDsOfs = m_FileHdr.FileLen;
		WrtLen = m_FileHdr.ChildGOIDsSize;

		if(m_bIsBigEndian)
			{
			INT32 *pGOChildID = m_pGOChildIDs;
			for(Idx = 0; Idx < m_ChildGOIDscnt; Idx++,pGOChildID++)
				*pGOChildID = SwapUI32Endians(*pGOChildID);
			}

		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,(char *)m_pGOChildIDs,WrtLen)!=WrtLen)
			{
			AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO child references to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FileLen += WrtLen;

		// now write out the GO term parent references
		m_FileHdr.ParentGOIDscnt = m_ParentGOIDscnt;
		m_FileHdr.ParentGOIDsSize = m_ParentGOIDscnt * sizeof(int);
		m_FileHdr.ParentGOIDsOfs = m_FileHdr.FileLen;
		WrtLen = m_FileHdr.ParentGOIDsSize;	

		if(m_bIsBigEndian)
			{
			INT32 *pGOParentID = m_pGOParentIDs;
			for(Idx = 0; Idx < m_ParentGOIDscnt; Idx++,pGOParentID++)
				*pGOParentID = SwapUI32Endians(*pGOParentID);
			}

		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,(char *)m_pGOParentIDs,WrtLen)!=WrtLen)
			{
			AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO parent references to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FileLen += WrtLen;


		// now write out the GO term alternative references
		m_FileHdr.AltGOIDscnt = m_AltGOIDscnt;
		m_FileHdr.AltGOIDsSize = m_AltGOIDscnt * sizeof(int);
		m_FileHdr.AltGOIDsOfs = m_FileHdr.FileLen;
		WrtLen = m_FileHdr.AltGOIDsSize;	

		if(m_bIsBigEndian)
			{
			INT32 *pAltGOIDs = m_pAltGOIDs;
			for(Idx = 0; Idx < m_AltGOIDscnt; Idx++,pAltGOIDs++)
				*pAltGOIDs = SwapUI32Endians(*pAltGOIDs);
			}

		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,(char *)m_pAltGOIDs,WrtLen)!=WrtLen)
			{
			AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO parent references to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FileLen += WrtLen;

		// now write out the GO terms
		m_FileHdr.GOCellTermcnt = m_GOCellTermcnt;
		m_FileHdr.GOBioTermcnt = m_GOBioTermcnt;	
		m_FileHdr.GOMolTermcnt = m_GOMolTermcnt;

		m_FileHdr.POAnatomyTermCnt = m_POAnatomyTermCnt;	
		m_FileHdr.PODevTermCnt = m_PODevTermCnt;

		m_FileHdr.GOTermcnt = m_GOTermcnt;
		m_FileHdr.GOTermSize = m_GOTermcnt * sizeof(tsGOTerm);
		m_FileHdr.GOTermOfs = m_FileHdr.FileLen;

		if(m_bIsBigEndian)
			{
			tsGOTerm *pGOTerm = m_pGOTerms;
			for(Idx = 0; Idx < m_GOTermcnt; Idx++,pGOTerm++)
				{
				pGOTerm->TermGOID.Pad64 = SwapUI64Endians(pGOTerm->TermGOID.Pad64);		// GO:Ident
				pGOTerm->GOTermID = SwapUI32Endians(pGOTerm->GOTermID);			// uniquely identifies this term instance (1..n)
				pGOTerm->TermNameIdx = SwapUI32Endians(pGOTerm->TermNameIdx);			// GO term name
				pGOTerm->TermDefIdx = SwapUI32Endians(pGOTerm->TermDefIdx);				// GO term definition
				pGOTerm->ReplacedByIdx = SwapUI32Endians(pGOTerm->ReplacedByIdx);			// GO term replacement GO:Ident
				pGOTerm->GOSupTermIdx = SwapUI32Endians(pGOTerm->GOSupTermIdx);			// GO:term superceding identifier
				pGOTerm->GOAltTermsIdx = SwapUI32Endians(pGOTerm->GOAltTermsIdx);			// list of GO:term alternative identifiers
				pGOTerm->GOParentIDsIdx = SwapUI32Endians(pGOTerm->GOParentIDsIdx);			// list of IsA GO:Idents
				pGOTerm->GOPartOfIDsIdx = SwapUI32Endians(pGOTerm->GOPartOfIDsIdx);			// list of part_of GO:Idents (treated as if parents) 
				pGOTerm->GOChildIDsIdx = SwapUI32Endians(pGOTerm->GOChildIDsIdx);			// list of child GO:Idents 
				pGOTerm->GOTermCntsIdx = SwapUI32Endians(pGOTerm->GOTermCntsIdx);			// associated counts for this term or 0 if no counts

				pGOTerm->NumParents = SwapUI16Endians(pGOTerm->NumParents);  // number of parents to this term as a child (GOParentIDs) 
				pGOTerm->NumPartOfs = SwapUI16Endians(pGOTerm->NumPartOfs); // number of part_of identifiers for this term (treated as if parents)
				pGOTerm->NumChildren = SwapUI16Endians(pGOTerm->NumChildren); // number of children to this term as a parent (GOChildIDs)
				pGOTerm->NumAltTermIDs = SwapUI16Endians(pGOTerm->NumAltTermIDs); // number of alternative identifiers for this term
				}
			}

		WrtLen = m_FileHdr.GOTermSize;	
		if(_lseeki64(m_hFile,m_FileHdr.FileLen,SEEK_SET) != m_FileHdr.FileLen ||
			write(m_hFile,(char *)m_pGOTerms,WrtLen)!=WrtLen)
			{
			AddErrMsg("CGOTerms::Flush2Disk","Unable to write GO terms to disk on file %s - error %s",m_szFile,strerror(errno));
			Reset(false);
			return(eBSFerrFileAccess);
			}
		m_FileHdr.FileLen += WrtLen;
		}

		// now write the header to disk
	if((Rslt = Hdr2Disk())!=eBSFSuccess)
		{
		Reset(false);
		return((teBSFrsltCodes)Rslt);
		}

	m_bHdrDirty = false;
	}
return(eBSFSuccess);
}


// TrimWhitespace
// Inplace trim any leading/trailing whitespace
// Returns number of chrs in trimmed string excluding terminating '\0'
int
CGOTerms::TrimWhitespace(char *pTxt)
{
int Len = 0;
char *pStart = pTxt;
char Chr;

if(pTxt == NULL || *pTxt == '\0')
	return(eBSFerrParams);

	// strip any leading whitespace
while((Chr = *pTxt) && isspace(Chr))
	pTxt++;
if(Chr == '\0')					// empty line?
	{
	*pStart = '\0';
	return(0);
	}

while(*pTxt)			// fast forward to string terminator '\0' 
	{
	*pStart++ = *pTxt++;
	Len++;
	}
pStart--;				// backup to last chr 
while(isspace(*pStart))
	{
	pStart--;
	Len--;
	}
pStart[1] = '\0';
return(Len);
}

//LocateStanzaTag
//Does a linear search for tsStanzaTag with tag name same as that requested
tsStanzaTag *
CGOTerms::LocateStanzaTag(tsStanzaType *pCurStanza,char *pszTagName)
{
tsStanzaTag *pCurTag;
int Idx;
if(pCurStanza == NULL || pszTagName == NULL || *pszTagName == '\0')
	return(NULL);
pCurTag = pCurStanza->pTags;
for(Idx = 0; Idx < pCurStanza->NumTags; Idx++,pCurTag++)
	if(!stricmp(pCurTag->pszTagName,pszTagName))
		return(pCurTag);
return(NULL);
}

int				// count of chrs returned in pszRetLine
CGOTerms::GetJoinedLine(char *pszRetLine,
						int MaxRetLineLen,
						FILE *pOBOstream)
{
char *pTxt;
char Chr;
char *pStart;
int CurLen = 0;
if(pszRetLine == NULL || MaxRetLineLen <= 0 || pOBOstream == NULL)
	return(eBSFerrParams);

while(fgets(&pszRetLine[CurLen],MaxRetLineLen - CurLen,pOBOstream)!= NULL)
	{
	m_CurOBOLineNum++;
	
	// simply slough lines which are just whitespace or which start with '!' as comment introducer
	pTxt = &pszRetLine[CurLen];
	while((Chr = *pTxt) && Chr != '!' && isspace(Chr))
		pTxt++;
	if(Chr == '\0' || Chr == '!')		// if an empty or completely commented out line
		continue;						// simply slough line
	
	pStart = &pszRetLine[CurLen];
	while((Chr = *pTxt++) && Chr != '\n' && Chr != '\r' && (Chr != '!' || pTxt[-2] == '\\'))			// copy down rest of line
		{
		*pStart++ = Chr;
		CurLen++;
		}
	*pStart = '\0';

	if((Chr == '\n' || Chr == '\r') && pStart[-1] == '\\')		// '\' at end of line means that line continues over onto next
		{
		CurLen--;
		continue;
		}
		
	break;
	}		
return(CurLen);	
}

//Unescape
//Inplace unescape string which could contain escaped chars
/* 
	\n  	 newline
	\W 	single space
	\t 	tab
	\: 	colon
	\, 	comma
	\" 	double quote
	\\ 	backslash
	\( 	open parenthesis
	\) 	close parenthesis
	\[ 	open bracket
	\] 	close bracket
	\{ 	open brace
	\} 	close brace
*/

bool
CGOTerms::Unescape(char *pszTxt)
{
char *pDst;
char Chr;
pDst = pszTxt;
if(pszTxt == NULL || *pszTxt == '\0')
	return(false);
while(Chr = *pszTxt++) {
	if(Chr == '\\')
		{
		switch(Chr = *pszTxt++) {
			case 'n':	//		\n  	 newline
				*pDst++ = '\n';
				break;
			case 'W':	//		\W 	single space
				*pDst++ = ' ';
				break;
			case 't':	//		\t 	tab
				*pDst++ = '\t';
				break;
			case ':':	//		\: 	colon
				*pDst++ = ':';
				break;
			case ',':	//		\, 	comma
				*pDst++ = ',';
				break;
			case '"':	//		\" 	double quote
				*pDst++ = '"';
				break;
			case '\\':	//		\\ 	backslash
				*pDst++ = '\\';
				break;
			case '(':	//		\( 	open parenthesis
				*pDst++ = '(';
				break;
			case ')':	//	 	\)  close parenthesis
				*pDst++ = ')';
				break;
			case '[':	//		\[ 	open bracket
				*pDst++ = '[';
				break;
			case ']':	//		\] 	close bracket
				*pDst++ = ']';
				break;
			case '{':	//		\{ 	open brace
				*pDst++ = '{';
				break;
			case '}':	//		\} 	close brace
				*pDst++ = '}';
				break;
			case '\0':			// special handling required for string terminator
				*pDst++ = '\\';
				*pDst = '\0';
				return(true);

			default:			// any other escaped chrs are treated as literal '\' followed by Chr
				*pDst++ = '\\';
				*pDst++ = Chr;
				break;
			}
		continue;
		}
	*pDst++ = Chr;
	}
return(true);
}

// ParseTagName
// Parses tag from pszTxt into pszTag
// If bUnescape then any escape chars will be translated to their single char equivalent
int						// number of chars parsed from pszTxt including ':' tag name terminator or one of 
CGOTerms::ParseTagName(char *pszTxt,int MaxTagNameLen,char *pszTagName,bool bUnescape)
{
int NumChrs = 0;
int TagNameLen = 0;
char *pDst;

if(pszTxt == NULL || *pszTxt == '\0' ||  pszTagName == NULL || MaxTagNameLen <= 0)
	return(eBSFerrParams);

*pszTagName = '\0';

// slough all preceding whitespace
while(*pszTxt && isspace(*pszTxt))
	{
	pszTxt++;
	NumChrs++;
	}
if(*pszTxt == '\0')
	return(0);

// copy into szTag until whitespace or ':'
pDst = pszTagName;
while(*pszTxt && !isspace(*pszTxt) && *pszTxt != ':')
	{
	if(TagNameLen < (MaxTagNameLen-1))
		{
		*pDst++ = *pszTxt++;
		TagNameLen++;
		}
	else
		pszTxt++;
	NumChrs++;
	}
*pDst = '\0';

// slough trailing whitespace upto terminating ':'
while(*pszTxt && isspace(*pszTxt))
	{
	pszTxt++;
	NumChrs++;
	}
if(*pszTxt != ':')
	{
	ReportParseError(eTSTErrKeyTerm,m_CurOBOLineNum,"");
	return(-1);
	}
if(bUnescape)
	Unescape(pszTagName);
return(NumChrs+1);
}

// StripQuotesWS
// Inplace strips any bracketing double quotes plus leading/trailing whitespace
int							// returns length of quote and whitespace stripped *pszTxt
CGOTerms::StripQuotesWS(char *pszTxt)
{
char QuoteChr = '\0';
int Len = 0;
char *pDst = pszTxt;

if(pszTxt == NULL || *pszTxt == '\0')
	return(0);

// first slough any leading whitespace
while(*pszTxt && isspace(*pszTxt))
	pszTxt++;
// check if left with empty string
if(*pszTxt == '\0')
	{
	*pDst = '\0';
	return(0);
	}

// check if single or double quote
if(*pszTxt == '"' || *pszTxt == '\'')
	{
	QuoteChr = *pszTxt;		// note which quote was used
	pszTxt++;
	// could be more whitespace, slough 
	while(*pszTxt && isspace(*pszTxt))
		pszTxt++;
	// check if left with empty string
	if(*pszTxt == '\0')
		{
		*pDst = '\0';
		return(0);
		}
	}

// any leading whitespace plus quote char have been trimed
while(*pszTxt)
	{
	Len++;
	*pDst++ = *pszTxt++;
	}

pDst--;	// backup to pt at last chr, not at terminating '\0'
while(isspace(*pDst))	// slough any trailing whitespace
	{
	pDst--;
	Len--;
	}

// check for quote chr of same type which started string
if(Len && *pDst == QuoteChr)	
	{
	pDst--;				// slough trailing quote chr
	Len--;
	// could be more whitespace, slough 
	while(Len && isspace(*pDst))
		{
		Len--;
		pDst--;
		}
	}
if(!Len)
	*pDst = '\0';
else
	pDst[1] = '\0';
return(Len);
}

// ParseTagValue
// Parses value from pszTxt into pszTagValue
// If bUnescape then any escape chars will be translated to their single char equivalent
int						// number of chars parsed from pszTxt 
CGOTerms::ParseTagValue(teTSTValueType Type,char *pszTxt,int MaxValLen,char *pszTagValue,bool bUnescape)
{
int ValLen;
int NumChrs = 0;
char *pChr;
*pszTagValue = '\0';

if(Type == eTSTValNone)
	return(0);
if(pszTxt == NULL || *pszTxt == '\0' || pszTagValue == NULL)
	return(eBSFerrParams);
while(*pszTxt && isspace(*pszTxt))
	{
	pszTxt++;
	NumChrs++;
	}
if(pszTxt[0] == '\0')
	return(NumChrs);

pChr = pszTagValue;
ValLen = 0;
switch(Type) {
	case eTSTTermID:			// term identifier, e.g. GO:123457
	case eTSTWord:				// unquoted single word e.g. text_string
	case eTSTBool:				// boolean, e.g true, false
		// copy into szTag until whitespace
		while(*pszTxt && !isspace(*pszTxt))
			{
			if(ValLen < (MaxValLen-1))
				{
				*pChr++ = *pszTxt++;
				ValLen++;
				}
			else
				pszTxt++;
			NumChrs++;
			}
		*pChr = '\0';
		break;
	
	case eTSTQuoteString:		// double quoted string e.g. "text string"
		// expecting double quote
		if(*pszTxt++ != '\"')
			return(-1);
		NumChrs++;

		// copy into szTag until string terminator or unescaped double quote
		while(*pszTxt)
			{
			NumChrs++;
			if(*pszTxt == '"' && pszTxt[-1] != '\\')
				break;
			if(ValLen < (MaxValLen-1))
				{
				*pChr++ = *pszTxt++;
				ValLen++;
				}
			else
				pszTxt++;
			}
		*pChr = '\0';
		TrimWhitespace(pszTagValue);
		break;

	case eTSTDbxref:			// dbxref
	case eTSTText:				// catchall general text which may be quoted/unquoted
		while(*pszTxt)
			{
			if(ValLen < (MaxValLen-1))
				{
				*pChr++ = *pszTxt++;
				ValLen++;
				}
			else
				pszTxt++;

			NumChrs++;
			}
		*pChr = '\0';
		// remove any bracketing quotes
		StripQuotesWS(pszTagValue);
		break;

	default:
		break;
	}

if(bUnescape)
	Unescape(pszTagValue);
return(NumChrs);
}

// ParseTagDbxref
// Parses value from pszTxt into pszTagDbxref with '[' and ']' stripped off
// If bUnescape then any escape chars will be translated to their single char equivalent
int						// number of chars parsed from pszTxt 
ParseTagDbxref(char *pszTxt,int MaxDbxrefLen,char *pszTagDbxref,bool bUnescape)
{
return(0);
}

// GetGOID(pszGOID)
// Does a linear search for matching GOID
// Returns eBSFerrGOID if no matching GOID or an offset 0..n if GOID located
int
CGOTerms::GetGOID(char *pszGOID)
{
unsigned short Len;
UINT16 Hash;
tsGOID *pID;
if(pszGOID == NULL || *pszGOID == '\0')
	return(eBSFerrParams);
if(m_pGOIDs == NULL || !m_GOIDcnt)
	return(eBSFerrGOID);
Hash = CUtility::GenHash16(pszGOID);
Len = (unsigned short)strlen(pszGOID);
pID = m_pGOIDs;
for(int Idx=0; Idx < m_GOIDcnt; Idx++)
	{
	if(pID->Len == Len && pID->Hash == Hash && !stricmp(pszGOID,pID->Txt))
		return((int)((char *)pID - (char *)m_pGOIDs));
	pID = (tsGOID *)((char *)pID + pID->Len + sizeof(tsGOID));
	}
return(eBSFerrGOID);	// couldn't locate
}

// AddGOID
// Adds specified GOID to list of known IDs
// If GOID already known then returns its unique identifier
// If GOID not known then adds to list and returns new identifier
// Identifiers returned are in range 0..n, note that they are offsets into m_pGOIDs and thus not contiguous
// < 0 returned if errors
int
CGOTerms::AddGOID(char *pszGOID)
{
UINT16 Hash;
unsigned char *pTmp;
tsGOID *pID;
unsigned short Len;
int GOIDofs;
int Mem2Alloc;

if(pszGOID == NULL || *pszGOID == '\0')
	return(eBSFerrParams);

// if identifier already known then return existing index
Hash = CUtility::GenHash16(pszGOID);
Len = (unsigned short)strlen(pszGOID);
if(m_pGOIDs != NULL && m_GOIDcnt)
	{
	pID = m_pGOIDs;
	for(int Idx=0; Idx < m_GOIDcnt; Idx++)
		{
		if(pID->Len == Len && pID->Hash == Hash && !stricmp(pszGOID,pID->Txt))
			return((int)((char *)pID - (char *)m_pGOIDs));
		pID = (tsGOID *)((char *)pID + pID->Len + sizeof(tsGOID));
		}
	}
// couldn't locate
if(m_pGOIDs == NULL || (m_NxtGOIDofs + (int)Len + (int)sizeof(tsGOID)) >= m_AllocdGOIDs)
	{
	if(m_pGOIDs == NULL)
		{
		m_NxtGOIDofs = 0;
		m_AllocdGOIDs = 0;
		m_GOIDcnt = 0;
		Mem2Alloc = 2 * cAllocGOIDs;
		}
	else
		Mem2Alloc = m_AllocdGOIDs + cAllocGOIDs;
	pTmp = new unsigned char[Mem2Alloc];
	if(pTmp == NULL)
		return(eBSFerrMem);

	if(m_pGOIDs != NULL)
		{
		if(m_NxtGOIDofs)
			memmove(pTmp,m_pGOIDs,m_NxtGOIDofs);
		delete m_pGOIDs;
		}
	else
		m_AllocdGOIDs = 0;
	m_pGOIDs = (tsGOID *)pTmp;
	m_AllocdGOIDs = Mem2Alloc;
	}

GOIDofs = m_NxtGOIDofs;
pID = (tsGOID *)((char *)m_pGOIDs + m_NxtGOIDofs);
pID->Hash = Hash;
pID->Len = Len;
strcpy(pID->Txt,pszGOID);
m_GOIDcnt++;
m_NxtGOIDofs += (Len + sizeof(tsGOID));
return(GOIDofs);
}

// GetGOID
// Returns ptr to GOID pTxt at specified GOIDOfs as returned by AddGOID(pTxt) or GetGOID(pTxt)
// Returns NULL if GOIDOfs not in range 0..m_NxtGOIDofs-1
char *
CGOTerms::GetGOID(int GOIDOfs)
{
tsGOID *pID;
if(m_pGOIDs == NULL || !m_GOIDcnt || GOIDOfs < 0 || GOIDOfs >= m_NxtGOIDofs)
	return(NULL);
pID = (tsGOID *)((char *)m_pGOIDs + GOIDOfs);
return(pID->Txt);
}


// LocateGOTagVal
// Does a linear search for matching Tag + pszVal
// Identifiers returned are in range 0..n, note that they are actually offsets into m_pGOTagVals and thus not contiguous
// Returns eBSFerrGOTagVal if no matching Tag + pszVal
int 
CGOTerms::LocateGOTagVal(teStanzaTags Tag,char *pszVal)
{
unsigned short Len;
UINT16 Hash;
tsGOTagVal *pTagVal;

if(pszVal == NULL || *pszVal == '\0')
	return(eBSFerrParams);

if(m_pGOTagVals == NULL || !m_GOTagValcnt)
	return(eBSFerrGOID);
Hash = CUtility::GenHash16(pszVal);
Len = (unsigned short)strlen(pszVal);
pTagVal = m_pGOTagVals;
for(int Idx=0; Idx < m_GOTagValcnt; Idx++)
	{
	if(pTagVal->Tag == Tag && pTagVal->Len == Len && pTagVal->Hash == Hash && !stricmp(pszVal,pTagVal->Val))
		return((int)((char *)pTagVal - (char *)m_pGOTagVals));
	pTagVal = (tsGOTagVal *)((char *)pTagVal + pTagVal->Len + sizeof(tsGOTagVal));
	}
return(eBSFerrGOTagVal);	// couldn't locate
}

// AddGOTagVal
// adds new tsGOTagVal to m_pGOTagVals or reuses existing if tag + value already in m_pGOTagVals
// If Tag + pszVal then returns its unique identifier
// If Tag + pszVal not known then adds to list and returns new identifier
// Identifiers returned are in range 0..n, note that they are actually offsets into m_pGOTagVals and thus not contiguous
// < 0 (eBSFerrMem) returned if errors
int 
CGOTerms::AddGOTagVal(teStanzaTags Tag,char *pszVal)
{
UINT16 Hash;
unsigned char *pTmp;
tsGOTagVal *pTagVal;
unsigned short Len;
int GOTagValOfs;
int Mem2Alloc;
if(pszVal == NULL || *pszVal == '\0')
	return(eBSFerrParams);

// if identifier already known then return existing index
Hash = CUtility::GenHash16(pszVal);
Len = (unsigned short)strlen(pszVal);
if(m_pGOTagVals != NULL && m_GOTagValcnt)
	{
	pTagVal = m_pGOTagVals;
	for(int Idx=0; Idx < m_GOTagValcnt; Idx++)
		{
		if(pTagVal->Tag == Tag && pTagVal->Len == Len && pTagVal->Hash == Hash && !stricmp(pszVal,pTagVal->Val))
			return((int)((char *)pTagVal - (char *)m_pGOTagVals));
		pTagVal = (tsGOTagVal *)((char *)pTagVal + pTagVal->Len + sizeof(tsGOTagVal));
		}
	}
// couldn't locate
if(m_pGOTagVals == NULL || (m_NxtGOTagValOfs + (int)Len + (int)sizeof(tsGOTagVal)) >= m_AllocdGOTagVals)
	{
	if(m_pGOTagVals == NULL)
		{
		m_NxtGOTagValOfs = 0;
		m_AllocdGOTagVals = 0;
		m_GOTagValcnt = 0;
		Mem2Alloc = 2 * cAllocGOTagVals;
		}
	else
		Mem2Alloc = m_AllocdGOTagVals + cAllocGOTagVals;
	pTmp = new unsigned char[Mem2Alloc];
	if(pTmp == NULL)
		return(eBSFerrMem);

	if(m_pGOTagVals != NULL)
		{
		if(m_NxtGOTagValOfs)
			memmove(pTmp,(char *)m_pGOTagVals,m_NxtGOTagValOfs);
		delete m_pGOTagVals;
		}
	m_pGOTagVals = (tsGOTagVal *)pTmp;
	m_AllocdGOTagVals = Mem2Alloc;
	}

GOTagValOfs = m_NxtGOTagValOfs;
pTagVal = (tsGOTagVal *)((char *)m_pGOTagVals + m_NxtGOTagValOfs);
pTagVal->Hash = Hash;
pTagVal->Len = Len;
pTagVal->Tag = Tag;
strcpy(pTagVal->Val,pszVal);
m_GOTagValcnt++;
m_NxtGOTagValOfs += (Len + sizeof(tsGOTagVal));
return(GOTagValOfs);
}

// GetGOVal
// returns pszVal for offset as returned from AddGOTagVal() or LocateGOTagVal()
// Returns NULL if GOTagValOfs not in range 0..m_NxtGOTagValOfs-1
char *
CGOTerms::GetGOVal(int GOTagValOfs)
{
tsGOTagVal *pTagVal;
if(m_pGOTagVals == NULL || !m_GOTagValcnt || GOTagValOfs < 0 || GOTagValOfs >= m_NxtGOTagValOfs)
	return(NULL);
pTagVal = (tsGOTagVal *)((char *)m_pGOTagVals + GOTagValOfs);
return(pTagVal->Val);
}


// GetGOTag
// returns Tag for offset as returned from AddGOTagVal() or LocateGOTagVal()
// Returns -1 if GOTagValOfs not in range 0..m_NxtGOTagValOfs-1
teStanzaTags
CGOTerms::GetGOTag(int GOTagValOfs)
{
tsGOTagVal *pTagVal;
if(m_pGOTagVals == NULL || !m_GOTagValcnt || GOTagValOfs < 0 || GOTagValOfs >= m_NxtGOTagValOfs)
	return((teStanzaTags)-1);
pTagVal = (tsGOTagVal *)((char *)m_pGOTagVals + GOTagValOfs);
return((teStanzaTags)pTagVal->Tag);
}

// AddTerm
// Adds specified term
int
CGOTerms::AddTerm(char *pszGOID,
		char *pszTermName,
		char *pszTermNamespace,
		char *pszTermDef,
		char *pszTermDefDbxref,
		char *pszTermComment,
		bool bIsObsolete,
		char *pszGOIDreplacedby,
		int NumPartOfGOIDs,    // 
		char *pszPartOfGOIDs,
		int NumGOIDAltIDs,
		char *pszGOIDAltIDs,
		int NumParentGOIDs,
		char *pszParentGOIDs)
{
tsGOTerm *pTmp;
int *piTmp;
int Cnt;
tsGOTerm *pTerm;
int *pID;
int Terms2Alloc;

if(pszGOID == NULL || pszGOID[0] == '\0' || 
   NumGOIDAltIDs < 0 || (NumGOIDAltIDs > 0 && (pszGOIDAltIDs == NULL || *pszGOIDAltIDs == '\0')) ||
   NumParentGOIDs < 0 || (NumParentGOIDs > 0 && (pszParentGOIDs == NULL || *pszParentGOIDs == '\0')))
	return(eBSFerrParams);

// time to allocate or realloc?
if(m_pGOTerms == NULL || (m_GOTermcnt + 1) >= m_AllocdGOTerms)
	{
	if(m_pGOTerms == NULL)
		{
		m_GOCellTermcnt=0;				// current number of terms classed as cellular
		m_GOBioTermcnt=0;					// current number of terms classed as biological
		m_GOMolTermcnt=0;					// current number of terms classed as molecular

		m_POAnatomyTermCnt=0;					// current number of terms classed as plant structure
		m_PODevTermCnt=0;					// current number of terms classed as plant structure developmental growth


		m_GOTermcnt = 0;
		m_AllocdGOTerms = 0;
		Terms2Alloc = 2 * cAllocTerms;
		}
	else
		Terms2Alloc = m_AllocdGOTerms + cAllocTerms;
	pTmp = new tsGOTerm[Terms2Alloc];
	if(pTmp == NULL)
		return(eBSFerrMem);
	if(m_pGOTerms != NULL)
		{
		if(m_GOTermcnt)
			memmove((char *)pTmp,(char *)m_pGOTerms,m_GOTermcnt * sizeof(tsGOTerm));
		delete m_pGOTerms;
		}
	m_pGOTerms = pTmp;
	m_AllocdGOTerms = Terms2Alloc;
	}

if(NumGOIDAltIDs && (m_pAltGOIDs == NULL || (m_AltGOIDscnt + NumGOIDAltIDs + 1) >= m_AllocdAltGOIDs))
	{
	if(m_pAltGOIDs == NULL)
		{
		m_AltGOIDscnt = 0;
		m_AllocdAltGOIDs = 0;
		Terms2Alloc = 2 * cAllocGOIDlists;
		}
	else
		Terms2Alloc = m_AllocdAltGOIDs + NumGOIDAltIDs + cAllocGOIDlists;
	piTmp = new int [Terms2Alloc];
	if(piTmp == NULL)
		return(eBSFerrMem);
	if(m_pAltGOIDs != NULL)
		{
		if(m_AltGOIDscnt)
			memmove((char *)piTmp,(char *)m_pAltGOIDs,m_AltGOIDscnt * sizeof(int));
		delete m_pAltGOIDs;
		}
	m_pAltGOIDs = piTmp;
	m_AllocdAltGOIDs = Terms2Alloc;
	}

if((NumParentGOIDs + NumPartOfGOIDs) && (m_pGOParentIDs == NULL || (m_ParentGOIDscnt + NumParentGOIDs + NumPartOfGOIDs + 1) >= m_AllocdGOParentIDs))
	{
	if(m_pGOParentIDs == NULL)
		{
		m_ParentGOIDscnt = 0;
		m_AllocdGOParentIDs = 0;
		Terms2Alloc = 2 * cAllocGOIDlists;
		}
	else
		Terms2Alloc = m_AllocdGOParentIDs + NumParentGOIDs + NumPartOfGOIDs + cAllocGOIDlists;
	piTmp = new int [Terms2Alloc];
	if(piTmp == NULL)
		return(eBSFerrMem);
	if(m_pGOParentIDs != NULL)
		{
		if(m_ParentGOIDscnt)
			memmove((char *)piTmp,(char *)m_pGOParentIDs,m_ParentGOIDscnt * sizeof(int));
		delete m_pGOParentIDs;
		}
	m_pGOParentIDs = piTmp;
	m_AllocdGOParentIDs = Terms2Alloc;
	}

pTerm = &m_pGOTerms[m_GOTermcnt++];
pTerm->GOTermID = m_GOTermcnt;
pTerm->TermNameIdx = eBSFerrGOTagVal;
pTerm->TermDefIdx = eBSFerrGOTagVal;
pTerm->ReplacedByIdx = eBSFerrGOID;
pTerm->GOAltTermsIdx = eBSFerrGOID;
pTerm->NumPartOfs = 0;
pTerm->GOPartOfIDsIdx = eBSFerrGOID;
pTerm->NumAltTermIDs = 0;
pTerm->GOParentIDsIdx = eBSFerrGOID;
pTerm->NumChildren = 0;
pTerm->GOChildIDsIdx = eBSFerrGOID;
pTerm->RootOntology = 0;
pTerm->NumParents = 0;
pTerm->NumPartOfs = 0;
pTerm->GOTermCntsIdx = 0;
pTerm->bIsObsolete = false;
pTerm->bNoClass = false;

pTerm->TermGOID.Idx = AddGOID(pszGOID);

pTerm->bIsObsolete = bIsObsolete;
if(bIsObsolete)
	pTerm->bNoClass = true;

if(pszTermName != NULL && pszTermName[0] != '\0')
	pTerm->TermNameIdx = AddGOTagVal(eTSTname,pszTermName);
if(pszTermDef != NULL && pszTermDef[0] != '\0')
	pTerm->TermDefIdx = AddGOTagVal(eTSTdef,pszTermDef);
if(bIsObsolete && (pszGOIDreplacedby != NULL && pszGOIDreplacedby[0] != '\0'))
	pTerm->ReplacedByIdx = AddGOID(pszGOIDreplacedby);

if(NumGOIDAltIDs)
	{
	pID = &m_pAltGOIDs[m_AltGOIDscnt];
	for(Cnt = 0; Cnt < NumGOIDAltIDs; Cnt++)
		{
		*pID++ = AddGOID(pszGOIDAltIDs);
		pszGOIDAltIDs += strlen(pszGOIDAltIDs)+1;
		}
	pTerm->GOAltTermsIdx = m_AltGOIDscnt;
	pTerm->NumAltTermIDs = NumGOIDAltIDs;
	m_AltGOIDscnt += NumGOIDAltIDs;
	}

if(pTerm->bNoClass)
	return(pTerm->GOTermID);

if(NumParentGOIDs)
	{
	pID = &m_pGOParentIDs[m_ParentGOIDscnt];
	for(Cnt = 0; Cnt < NumParentGOIDs; Cnt++)
		{
		*pID++ = AddGOID(pszParentGOIDs);
		pszParentGOIDs += strlen(pszParentGOIDs)+1;
		}
	pTerm->GOParentIDsIdx = m_ParentGOIDscnt;
	pTerm->NumParents = NumParentGOIDs;
	m_ParentGOIDscnt += NumParentGOIDs;
	}

if(NumPartOfGOIDs)
	{
	pID = &m_pGOParentIDs[m_ParentGOIDscnt];
	for(Cnt = 0; Cnt < NumPartOfGOIDs; Cnt++)
		{
		*pID++ = AddGOID(pszPartOfGOIDs);
		pszPartOfGOIDs += strlen(pszPartOfGOIDs)+1;
		}
	pTerm->GOPartOfIDsIdx = m_ParentGOIDscnt;
	pTerm->NumPartOfs = NumPartOfGOIDs;
	m_ParentGOIDscnt += NumPartOfGOIDs;
	}


if(!NumParentGOIDs && !NumPartOfGOIDs)	// if no parents then should be a root term
	{
	// check if plant ontology
	if(!stricmp(pszTermName,pszRootTerms[3]))	// plant anatomical?
		{
		if(m_FileHdr.POAnatomicalID)
			{
			AddErrMsg("CGOTerms::AddTerm","Root term - %s - already parsed",pszTermName);
			return(eBSFerrDupGOTerm);
			}
		m_FileHdr.POAnatomicalID = pTerm->GOTermID;	// identifies plant structure term
		return(pTerm->GOTermID);
		}

	if(!stricmp(pszTermName,pszRootTerms[4]))	// plant structure development stage?
		{
		if(m_FileHdr.POdevID)
			{
			AddErrMsg("CGOTerms::AddTerm","Root term - %s - already parsed",pszTermName);
			return(eBSFerrDupGOTerm);
			}
		m_FileHdr.POdevID = pTerm->GOTermID;	// identifies  plant structure development stage term
		return(pTerm->GOTermID);
		}

	if(!stricmp(pszTermName,pszRootTerms[0]))	// cellular component?
		{
		if(m_FileHdr.CellTermID)
			{
			AddErrMsg("CGOTerms::AddTerm","Root term - %s - already parsed",pszTermName);
			return(eBSFerrDupGOTerm);
			}
		m_FileHdr.CellTermID = pTerm->GOTermID;	// identifies root cellular component term
		return(pTerm->GOTermID);
		}

	if(!stricmp(pszTermName,pszRootTerms[1]))	// biological process?
		{
		if(m_FileHdr.BioTermID)
			{
			AddErrMsg("CGOTerms::AddTerm","Root term - %s - already parsed",pszTermName);
			return(eBSFerrDupGOTerm);
			}
		m_FileHdr.BioTermID = pTerm->GOTermID;	// identifies root biological process term 
		return(pTerm->GOTermID);
		}
	
	if(!stricmp(pszTermName,pszRootTerms[2]))	// molecular function?
		{
		if(m_FileHdr.MolTermID)
			{
			AddErrMsg("CGOTerms::AddTerm","Root term - %s - already parsed",pszTermName);
			return(eBSFerrDupGOTerm);
			}
		m_FileHdr.MolTermID = pTerm->GOTermID;	// identifies root molecular function term 
		return(pTerm->GOTermID);
		}
	
	AddErrMsg("CGOTerms::AddTerm","Root term - %s - not recognised",pszTermName);
	return(eBSFerrGOID);
	}
	
return(pTerm->GOTermID);
}



// Parse
// Parses opened OBO format file
int
CGOTerms::Parse(FILE *pOBOstream)
{
char szTagName[cMaxOBOTagLen+1];		// to hold each tag as parsed from szLineBuff
char szStanzaType[cMaxOBOTagLen+1];		// to hold each stanza type as parsed from "[stanzatype]"
char szGOID[cMaxOBOGOID+1];
char szTermName[cMaxOBOGOname+1];
char szTermNamespace[cMaxOBOGOnamespace+1];
char szTermDef[cMaxOBOGOdef+1];
char szTermDefDbxref[cMaxOBOGOdDbxref+1];
char szTermComment[cMaxOBOGOcomment+1];
char szGOIDreplacedby[cMaxOBOGOID+1];


char szTmp[cMaxOBOGOcomment+1];		// shorter intermediate parsing
char szTmpLong[cMaxOBOGOcomment+1];	// to hold longer intermediate parse values
bool bIsObsolete;

char szParentGOIDs[cMaxOBOParentGOIDs*(cMaxOBOGOID+1)];
int NumParentGOIDs;
char *pszNxtParentGOID;

char szGOIDAltIDs[cMaxOBOAltGOIDs*(cMaxOBOGOID+1)];
int NumGOIDAltIDs;
char *pszNxtGOIDAltID;

char szPartOfGOIDs[cMaxOBOPartOfGOIDs*(cMaxOBOGOID+1)];
char *pszNxtPartOfGOID;
int NumPartOfGOIDs;

tsStanzaType *pCurStanza;
tsStanzaType *pStanza;
tsStanzaTag *pTag;

bool bSloughStanza;
int Rslt = -1;
char szJoinedLine[cMaxOBOJoinedLineLen+1];	// to hold each joined line, this is the line which is parsed

char *pszTxt;							// ptr into szLineBuff
char *pszTag;							// ptr into szTag
int LineNum;
int LineLen;
int Idx;
int ParseLen;

LineNum = 0;
pCurStanza = &StanzaTypes[0];		// start by defaulting stanza type to document header
bSloughStanza = false;
szTagName[0]='\0';		// to hold each tag as parsed from szLineBuff
szGOID[0]= '\0';
szTermName[0] = '\0';
szTermNamespace[0]= '\0';
szTermDef[0]= '\0';
szTermDefDbxref[0]= '\0';
szTermComment[0]= '\0';
szGOIDreplacedby[0]='\0';
szGOIDAltIDs[0] = '\0';
pszNxtGOIDAltID = szGOIDAltIDs;
NumGOIDAltIDs = 0;
szParentGOIDs[0] = '\0';
pszNxtParentGOID = szParentGOIDs;
NumParentGOIDs = 0;
szPartOfGOIDs[0] = '\0';
pszNxtPartOfGOID = szPartOfGOIDs;
NumPartOfGOIDs = 0;
bIsObsolete = false;


while((LineLen = GetJoinedLine(szJoinedLine,cMaxOBOJoinedLineLen,pOBOstream))>0)
	{
	pszTxt = szJoinedLine;
	// check if new stanza type starting or can continue with existing
	if(*pszTxt == '[')	// '[' introduces new stanza terminating the current stanza
		{
		if(szGOID[0] != '\0')
			{
			Rslt = AddTerm(szGOID,szTermName,szTermNamespace,szTermDef,szTermDefDbxref,szTermComment,bIsObsolete,szGOIDreplacedby,NumPartOfGOIDs,szPartOfGOIDs,NumGOIDAltIDs,szGOIDAltIDs,NumParentGOIDs,szParentGOIDs);
			if(Rslt < eBSFSuccess)
				{
				ReportParseError(eTSTErrGenErr,m_CurOBOLineNum,"AddTerm error");
				return(Rslt);
				}
			szTagName[0]='\0';		// to hold each tag as parsed from szLineBuff
			szGOID[0]= '\0';
			szTermName[0] = '\0';
			szTermNamespace[0]= '\0';
			szTermDef[0]= '\0';
			szTermDefDbxref[0]= '\0';
			szTermComment[0]= '\0';
			szGOIDreplacedby[0]='\0';
			szGOIDAltIDs[0] = '\0';
			pszNxtGOIDAltID = szGOIDAltIDs;
			NumGOIDAltIDs = 0;
			szParentGOIDs[0] = '\0';
			pszNxtParentGOID = szParentGOIDs;
			NumParentGOIDs = 0;
			szPartOfGOIDs[0] = '\0';
			pszNxtPartOfGOID = szPartOfGOIDs;
			NumPartOfGOIDs = 0;
			bIsObsolete = false;
			}	
		pszTxt++;
		// slough any leading whitespace
		while(*pszTxt && isspace(*pszTxt))
			pszTxt++;
		if(*pszTxt == '\0' || *pszTxt == ']')	// error if no or empty stanza specified
			{
			ReportParseError(eTSTErrGenErr,m_CurOBOLineNum,"Empty stanza type '[]'");
			return(-1);
			}
		pszTag = szStanzaType;
		while(*pszTxt && !isspace(*pszTxt) && *pszTxt != ']')
			*pszTag++ = *pszTxt++;
		if(*pszTxt == '\0')
			{
			ReportParseError(eTSTErrEOL,m_CurOBOLineNum,"Expected stanza type closing ']'");
			return(-1);
			}
		*pszTag = '\0';
		while(*pszTxt && isspace(*pszTxt))
			pszTxt++;
		if(*pszTxt != ']')
			{
			ReportParseError(eTSTErrEOL,m_CurOBOLineNum,"Expected stanza type closing ']'");
			return(-1);
			}
		pStanza = StanzaTypes;
		for(Idx = 0; Idx < cNumDocStanzaTypes; Idx++,pStanza++)
			{
			if(!stricmp(pStanza->pszStanza,szStanzaType))
				break;
			}
		if(Idx >= cNumDocStanzaTypes)
			{
			bSloughStanza = true;
			continue;
			}
		pCurStanza = pStanza;	// know how to parse this stanza
		bSloughStanza = false;
		szTagName[0]='\0';		// to hold each tag as parsed from szLineBuff
		szGOID[0]= '\0';
		szTermName[0] = '\0';
		szTermNamespace[0]= '\0';
		szTermDef[0]= '\0';
		szTermDefDbxref[0]= '\0';
		szTermComment[0]= '\0';
		szGOIDreplacedby[0]='\0';
		szGOIDAltIDs[0] = '\0';
		pszNxtGOIDAltID = szGOIDAltIDs;
		NumGOIDAltIDs = 0;
		szParentGOIDs[0] = '\0';
		pszNxtParentGOID = szParentGOIDs;
		NumParentGOIDs = 0;
		szPartOfGOIDs[0] = '\0';
		pszNxtPartOfGOID = szPartOfGOIDs;
		NumPartOfGOIDs = 0;		
		bIsObsolete = false;
		continue;
		}
	if(bSloughStanza)		// if can't handle current stanza then slough lines until [stanza] we can handle
		continue;

	// assuming that tag name is now being parsed as it is not a stanza starting
	ParseLen = ParseTagName(pszTxt,sizeof(szTagName),szTagName,true);
	if(ParseLen < 0)
		{
		ReportParseError(eTSTErrKeyTerm,m_CurOBOLineNum,"");
		return(-1);
		}
	else
		if(!ParseLen)
			continue;
	pszTxt+=ParseLen;

	pTag = LocateStanzaTag(pCurStanza,szTagName);
	if(pTag == NULL)		// simply slough unrecognised tags
		continue;

	switch(pCurStanza->StanzaType) {
		case eTSTDocHdr:					// parsing doc header
			switch(pTag->TagID) {
				case eTSTformat_version:	//Gives the obo specification version that this file uses
				case eTSTdata_version:		//Gives the version of the current ontology.
				case eTSTversion:			//Deprecated. Use data-version instead.
				case eTSTdefault_namespace:	//default namespace
					break;
				default:					// not currently interested in other header tags
					break;
				}
			break;

		case eTSTTerm:					// parsing [Term] stanza
			switch(pTag->TagID) {
				case eTSTid:			//The unique id of the current term
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szGOID),szGOID,false);
					break;

				case eTSTname:			//The term name
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szTermName),szTermName,false);
					break;

				case eTSTnamespace:		//The namespace in which the term belongs
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szTermNamespace),szTermNamespace,false);
					break;

				case eTSTdef:			//The definition of the current term
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szTermDef),szTermDef,true);
					if(ParseLen > 0)
						ParseLen=ParseTagValue((teTSTValueType)pTag->Type2,pszTxt+ParseLen,sizeof(szTermDefDbxref),szTermDefDbxref,true);
					break;

				case eTSTcomment:		//A comment for this term
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szTermComment),szTermComment,true);
					break;

				case eTSTis_a:			//This tag describes a subclassing relationship between one term and another
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,(int)sizeof(szParentGOIDs) - (int)(pszNxtParentGOID-szParentGOIDs),pszNxtParentGOID,false);
					if(ParseLen > 0)
						{
						NumParentGOIDs++;
						pszNxtParentGOID += strlen(pszNxtParentGOID) + 1;
						}
					break;

				case eTSTalt_id:	// an alternative id for this term
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,(int)sizeof(szGOIDAltIDs) - (int)(pszNxtGOIDAltID-szGOIDAltIDs),pszNxtGOIDAltID,false);
					if(ParseLen > 0)
						{
						NumGOIDAltIDs++;
						pszNxtGOIDAltID += strlen(pszNxtGOIDAltID) + 1;
						}
					break;

				case eTSTis_obsolete:
					bIsObsolete = false;
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szTmp),szTmp,false);
					if(ParseLen > 0 && !stricmp(szTmp,"true"))
						bIsObsolete = true;
					break;

				case eTSTreplaced_by:	// Term is replaced by this GO:ID
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szGOIDreplacedby),szGOIDreplacedby,false);
					break;

				case eTSTrelationship:
					*pszNxtPartOfGOID = '\0';
					ParseLen=ParseTagValue((teTSTValueType)pTag->Type1,pszTxt,sizeof(szTmpLong),szTmpLong,false);
					if(ParseLen > 20)
						{
						if(2 == sscanf(szTmpLong,"%s %s ",szTmp,pszNxtPartOfGOID))
							{
							if(!stricmp(szTmp,"part_of") && !strnicmp(pszNxtPartOfGOID,"GO:",3))
								{
								NumPartOfGOIDs++;
								pszNxtPartOfGOID += strlen(pszNxtPartOfGOID) + 1;
								}
							}
						}
					break;

				default:
					break;
				}
			break;
		
		case eTSTTypeDef:		// parsing [TypeRef] stanza
		case eTSTInstance:		// parsing [Instance] stanza
		default:
			break;
		}
	if(ParseLen <= 0)
		printf("\nProblems!!");
	}
if(szGOID[0] != '\0')
	AddTerm(szGOID,szTermName,szTermNamespace,szTermDef,szTermDefDbxref,szTermComment,bIsObsolete,szGOIDreplacedby,NumPartOfGOIDs,szPartOfGOIDs,NumGOIDAltIDs,szGOIDAltIDs,NumParentGOIDs,szParentGOIDs);

// check if all referenced root terms identified
if(!(m_FileHdr.POdevID || m_FileHdr.POAnatomicalID || m_FileHdr.CellTermID || m_FileHdr.BioTermID || m_FileHdr.MolTermID))
	{
	AddErrMsg("CGOTerms::Parse","Root term - Expected at least one root term, no root terms parsed");
	return(eBSFerrGOID);
	}

if(m_FileHdr.CellTermID || m_FileHdr.BioTermID || m_FileHdr.MolTermID)
	{
	if(!m_FileHdr.CellTermID)			// identifies root cellular component term
		{
		AddErrMsg("CGOTerms::Parse","Root term - %s - not parsed",pszRootTerms[0]);
		return(eBSFerrGOID);
		}

	if(!m_FileHdr.BioTermID)				// identifies root biological process term 
		{
		AddErrMsg("CGOTerms::Parse","Root term - %s - not parsed",pszRootTerms[1]);
		return(eBSFerrGOID);
		}

	if(!m_FileHdr.MolTermID)				// identifies root molecular function term 
		{
		AddErrMsg("CGOTerms::Parse","Root term - %s - not parsed",pszRootTerms[2]);
		return(eBSFerrGOID);
		}
	}

GenBackRefs();

if(m_GOTermcnt > 1  && m_pGOTerms != NULL) 
	{
	SwitchTermGOIDIDX2PTR();
	qsort(m_pGOTerms,m_GOTermcnt,sizeof(tsGOTerm),SortTermsByGOID);
	SwitchTermGOIDPTR2IDX();
	
	}

Rslt = SetOntologyClass4Terms();

return(Rslt);
}


// GenBackRefs
// Generate back references to all children of each GO:Term node
// Processes by iterating over all terms, for each term iterating over the parents and
// adding the term as being a child of the parent
// Back references are not generated for obsoleted terms
// 
int
CGOTerms::GenBackRefs(void)
{
int Parent;
tsGOTerm *pParent;
int Child;
tsGOTerm *pChild;
int NumParents;
int *pIdx;
int Terms2Alloc;
int *piTmp;

pParent = m_pGOTerms;
for(Parent = 0; Parent < m_GOTermcnt; Parent++,pParent++)
	{
	if(pParent->bNoClass)	// classless terms have no children
		continue;
	
	pParent->NumChildren = 0;
	pParent->GOChildIDsIdx = eBSFerrGOID;
	pChild = m_pGOTerms; 
	for(Child = 0; Child < m_GOTermcnt; Child++,pChild++)
		{
		NumParents = pChild->NumParents + pChild->NumPartOfs;

		if(Parent == Child ||			// parent can't be child of same!
			pChild->bNoClass ||			// if obsolete then no parents
			!NumParents)				// or child is a root
			continue;
		if(pChild->NumParents)
			pIdx = &m_pGOParentIDs[pChild->GOParentIDsIdx];
		else
			pIdx = &m_pGOParentIDs[pChild->GOPartOfIDsIdx];

		while(NumParents--)
			{
			if(*pIdx++ == pParent->TermGOID.Idx)
				{
				// bingo, add the child to the parent
				if(m_pGOChildIDs == NULL || (m_ChildGOIDscnt + 1) >= m_AllocdGOChildIDs)
					{
					if(m_pGOChildIDs == NULL)
						{
						m_ChildGOIDscnt = 0;
						m_AllocdGOChildIDs = 0;
						Terms2Alloc = 2 * cAllocGOIDlists;
						}
					else
						Terms2Alloc = m_AllocdGOChildIDs + 1 + cAllocGOIDlists;
					piTmp = new int [Terms2Alloc];
					if(piTmp == NULL)
						return(eBSFerrMem);
					if(m_pGOChildIDs != NULL)
						{
						if(m_ChildGOIDscnt)
							memmove((char *)piTmp,(char *)m_pGOChildIDs,m_ChildGOIDscnt * sizeof(int));
						delete m_pGOChildIDs;
						}
					m_pGOChildIDs = piTmp;
					m_AllocdGOChildIDs = Terms2Alloc;
					}

				if(!pParent->NumChildren)
					pParent->GOChildIDsIdx = m_ChildGOIDscnt;
				m_pGOChildIDs[m_ChildGOIDscnt++] = pChild->TermGOID.Idx;
				pParent->NumChildren++;
				break;
				}
			}
		}
	}
return(eBSFSuccess);
}

// NumGOTerms
// Returns total number of GO:terms in specified RootOntology
// Note that this includes obsoleted terms if RootOntology == eONTAny
//
int 
CGOTerms::NumGOTerms(etOntologies RootOntology)		
{
if(!m_bGOTermsAvail || m_pGOTerms == NULL)
	return(0);
switch(RootOntology) {
	case eONTCellular:		// Cellular component
		return(m_GOCellTermcnt);

	case eONTBiological:	// Biological process
		return(m_GOBioTermcnt);

	case eONTMolecular:		// Molecular function
		return(m_GOMolTermcnt);

	case eONTPlantAnatomical:	// plant structure
		return(m_POAnatomyTermCnt);

	case eONTPlantDev:		// plant structure developmental growth
		return(m_PODevTermCnt);

	default:
		break;
	}
return(m_GOTermcnt);
}

// GetNumTermCnts
// Returns number of term counts which have been associated with terms
// in specified RootOntology
// Note that as only non-obsoleted terms can have counts then returned numbers exclude obsoleted terms
int 
CGOTerms::GetNumTermCnts(etOntologies RootOntology)			
{
int NumTerms;
int Idx;
tsGOTermCnts *pCnt;

if((pCnt=m_pGOTermCnts) == NULL || !m_GOTermCntsCnt)
	return(0);

for(NumTerms = 0,Idx = 0; Idx < m_GOTermCntsCnt; Idx++,pCnt++)
	if(RootOntology == pCnt->RootOntology)
		NumTerms++;

return(NumTerms);
}

// GetNumSampledTermCnts
// Returns number of termcounts which have sampled gene hits of MinGenes or greater
int 
CGOTerms::GetNumSampledTermCnts(etOntologies RootOntology,unsigned int MinGenes)	
{
int NumTerms;
int Idx;
tsGOTermCnts *pCnt;
if((pCnt=m_pGOTermCnts) == NULL || !m_GOTermCntsCnt)
	return(0);
for(NumTerms = 0,Idx = 0; Idx < m_GOTermCntsCnt; Idx++,pCnt++)
	if(RootOntology == pCnt->RootOntology && pCnt->NumSampleGenes >= MinGenes)
		NumTerms++;
return(NumTerms);
}

// GetNumBkgndTermCnts
// Returns number of termcounts which have background (population) gene hits of MinGenes or greater
int 
CGOTerms::GetNumBkgndTermCnts(etOntologies RootOntology,unsigned int MinGenes)	
{
int NumTerms;
int Idx;
tsGOTermCnts *pCnt;
if((pCnt=m_pGOTermCnts) == NULL || !m_GOTermCntsCnt)
	return(0);
for(NumTerms = 0,Idx = 0; Idx < m_GOTermCntsCnt; Idx++,pCnt++)
	if(RootOntology == pCnt->RootOntology && pCnt->NumBkgdGenes >= MinGenes)
		NumTerms++;
return(NumTerms);
}

// GetTermCnts
// returns the Ith, 1..GetNumTermCnts(), term counts
tsGOTermCnts *
CGOTerms::GetTermCnts(int Ith) 
{
if(m_pGOTermCnts == NULL || Ith < 1 || Ith > m_GOTermCntsCnt)
	return(NULL);
return(&m_pGOTermCnts[Ith-1]);
}

// GetTermCnts
// Returns ptr to counts associated with specified term
// If term is obsolete (or allocation error) then returns NULL
// Allocates new counts if none exist
tsGOTermCnts *
CGOTerms::GetTermCnts(tsGOTerm *pTerm)
{
tsGOTermCnts *pCnts;
if(pTerm->bNoClass)	// can't associate counts with an obsolete or classless term
	return(NULL);

if(pTerm->GOTermCntsIdx < 1)
	{
	if(m_pGOTermCnts == NULL)
		{
		m_pGOTermCnts = new tsGOTermCnts [m_GOTermcnt];
		if(m_pGOTermCnts == NULL)
			return(NULL);
		m_GOTermCntsCnt = 0;
		m_NxtGOTermCntsIdx = 0;
		m_AllocdGOTermCnts = m_GOTermcnt;
		}
	pCnts = &m_pGOTermCnts[m_NxtGOTermCntsIdx++];
	memset(pCnts,0,sizeof(tsGOTermCnts));
	pCnts->UpdateSeq = -1;
	pCnts->TermID = pTerm->GOTermID;
	pCnts->RootOntology = pTerm->RootOntology;
	pTerm->GOTermCntsIdx = m_NxtGOTermCntsIdx;
	m_GOTermCntsCnt = m_NxtGOTermCntsIdx;
	}
pCnts = &m_pGOTermCnts[pTerm->GOTermCntsIdx-1];
return(pCnts);
}


// GetExistingTermCnts
// Returns ptr to counts associated with specified term
// If term is obsolete or no counts associated then returns NULL
tsGOTermCnts *
CGOTerms::GetExistingTermCnts(tsGOTerm *pTerm)
{
tsGOTermCnts *pCnts;
if(pTerm->bNoClass)	// can't associate counts with an obsolete or classless term
	return(NULL);
if(pTerm->GOTermCntsIdx < 1 || pTerm->GOTermCntsIdx > m_NxtGOTermCntsIdx)
	return(NULL);
pCnts = &m_pGOTermCnts[pTerm->GOTermCntsIdx-1];
return(pCnts);
}

// ClearStats
// Clears all GOTerm count references
int 
CGOTerms::ClearStats(void)		 
{
int Idx;
tsGOTermCnts *pCnts;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt == 0)
	return(eBSFSuccess);

pCnts = m_pGOTermCnts;
for(Idx=0; Idx < m_GOTermCntsCnt; Idx++,pCnts++)
	m_pGOTerms[pCnts->TermID-1].GOTermCntsIdx = 0;
m_GOTermCntsCnt = 0;
m_NxtGOTermCntsIdx = 0;
return(eBSFSuccess);
}


int 
CGOTerms::ClearSampleCounts(void)		 // clears all sample counts/stats
{
int Idx;
tsGOTermCnts *pCnts;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt == 0)
	return(eBSFSuccess);
pCnts = m_pGOTermCnts;
for(Idx=0; Idx < m_GOTermCntsCnt; Idx++,pCnts++)
	{
	pCnts->UpdateSeq = -1;
	pCnts->SampleCnt = 0;
	pCnts->NumSampleGenes = 0;
	}
m_CurUpdateSeq = 0;
return(eBSFSuccess);

}


// ResetUpdateSeqs
// resets all UpdateSeqs back to -1 without resetting the counts or references from GOTerms to counts
void 
CGOTerms::ResetUpdateSeqs(void)			
{
int Idx;
tsGOTermCnts *pCnts;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt == 0)
	return;
pCnts = m_pGOTermCnts;
for(Idx=0; Idx < m_GOTermCntsCnt; Idx++,pCnts++)
	pCnts->UpdateSeq = -1;
m_CurUpdateSeq = 0;
}


// ClampCount
// Adds Increment to Count and returns
// If Increment + Count would overflow cMaxCountVal then returns cMaxCountVal
unsigned int
CGOTerms::ClampCount(unsigned int Count,unsigned int Increment)
{
if(((INT64)Count + Increment) > cClampCountVal)
	return(cClampCountVal);
return(Count + Increment);
}

int
CGOTerms::AddCountRecurse(tsGOTerm *pTerm, // recurse into parents of this term
				   bool bSample,	// if true then update SampleCnts. if false then update background counts
				   int Count)		// count to increment by (1..n)
{
int Rslt;
int ParentIdx;
tsGOID *pCurID;
tsGOTerm *pParent;
tsGOTermCnts *pCnts;
int NumParents;

if(pTerm == NULL || Count == 0)
	return(eBSFerrParams);

NumParents = pTerm->NumParents + pTerm->NumPartOfs;
if(NumParents < 1)		// check if any parents, if none then already at root node so simple return
	return(eBSFSuccess);

for(ParentIdx = 0; ParentIdx < NumParents; ParentIdx++)
	{
	if(pTerm->NumParents)
		pCurID = (tsGOID *)((char *)m_pGOIDs + m_pGOParentIDs[pTerm->GOParentIDsIdx + ParentIdx]);
	else
		pCurID = (tsGOID *)((char *)m_pGOIDs + m_pGOParentIDs[pTerm->GOPartOfIDsIdx + ParentIdx]);

	if((pParent = LocateGOID(pCurID->Txt))==NULL)
		return(eBSFerrGOID);
	
	pCnts = GetTermCnts(pParent);
	if(pCnts == NULL)
		return(eBSFerrMem);

	if(pCnts->UpdateSeq != m_CurUpdateSeq)	// already updated this term?
		{
		if(bSample)
			{
			pCnts->SampleCnt = ClampCount(pCnts->SampleCnt,Count);
			pCnts->NumSampleGenes = ClampCount(pCnts->NumSampleGenes,1);
			}
		else
			{
			pCnts->BkgdCnt = ClampCount(pCnts->BkgdCnt,Count);
			pCnts->NumBkgdGenes = ClampCount(pCnts->NumBkgdGenes,1);
			}

		pCnts->UpdateSeq = m_CurUpdateSeq;		// term updated
		if((pParent->NumParents + pParent->NumPartOfs) > 0)
			if((Rslt=AddCountRecurse(pParent,bSample,Count))!=eBSFSuccess)
				return(Rslt);
		}
	}
return(eBSFSuccess);
}


int				// number of terms in pszGOID to which count was actually incremented
CGOTerms::AddCount(etOntologies OntologyClass, // which class of ontologies to add count to
				   bool bProp,		// if true then propagate counts into parents
				   bool bSample,	// if true then update SampleCnts. if false then update background counts
				   int Count,		// count to increment by (1..n)
				   int NumGOIDs,	// number of GO:Term identifiers in pszGOID[]
				   char *pszGOID[]) // array of ptrs to GO:Term identifiers 
{
int Rslt;
tsGOTerm *pTerm;
char *pszCurChildGOID;
int ChildIdx;
tsGOTermCnts *pCnts;
int NumTerms = 0;

if(++m_CurUpdateSeq > 0x3fffffff)
	{
	ResetUpdateSeqs();	
	m_CurUpdateSeq = 1;
	}

for(ChildIdx = 0; ChildIdx < NumGOIDs; ChildIdx++)
	{
	pszCurChildGOID = *pszGOID++;
	pTerm = LocateGOID(pszCurChildGOID);
	if(pTerm == NULL)		// if GO:term not present then should be an error, but if treated as 
		continue;			// a hard error then too many GOs would fail

	if(pTerm->bNoClass)	// not part of any class 
		continue;

	if(OntologyClass != pTerm->RootOntology)
		continue;

	pCnts = GetTermCnts(pTerm);
	if(pCnts == NULL)
		return(eBSFerrMem);

	if(pCnts->UpdateSeq == m_CurUpdateSeq)	// already updated this term?
		continue;

	if(bSample)
		{
		pCnts->SampleCnt = ClampCount(pCnts->SampleCnt,Count);
		pCnts->NumSampleGenes = ClampCount(pCnts->NumSampleGenes,1);
		}
	else
		{
		pCnts->BkgdCnt = ClampCount(pCnts->BkgdCnt,Count);
		pCnts->NumBkgdGenes = ClampCount(pCnts->NumBkgdGenes,1);
		}

	pCnts->UpdateSeq = m_CurUpdateSeq;		// term updated

	if(bProp)
		if((Rslt=AddCountRecurse(pTerm,bSample,Count))!=eBSFSuccess)
			return(Rslt);
	NumTerms++;
	}
return(NumTerms);
}

int
CGOTerms::SetBkgndCnts(etOntologies OntologyClass, // which class of ontologies to set counts for
		bool bProp,	// propagate counts from GO:Term into parent terms
		bool bBkgndLen,			// background counts proportional to gene lengths
		int Strand,				// background counts are for which strand 0==both,1=='+' and 2=='-'
		CBEDfile *pBED,			// gene BED file
		CGOAssocs *pAssocs)		// gene to ID association file
{
int Rslt;
int TotbBkgndLen;
int CurFeatureID;
char szGene[100];
char szChrom[100];
int GeneStart;
int GeneEnd;
int GeneLen;
char GeneStrand;
int NumGOIDs;
int Cnt;
char *pszGOID;
int NumFeatures;
int LocatedTermCnt;
int NotLocatedTermCnt;
int LocatedGOIDs;
int NotLocatedGOIDs;

char *pszGOIDs[1000];
char szGOIDs[10000];
char *pszConcatGOIDs;
int CntGOIDs;

LocatedTermCnt = 0;
NotLocatedTermCnt = 0;
LocatedGOIDs = 0;
NotLocatedGOIDs = 0;
CurFeatureID = 0;
TotbBkgndLen = 0;
NumFeatures = pBED->GetNumFeatures(); 
while((CurFeatureID = pBED->GetNextFeatureID(CurFeatureID)) > 0)
	{
	Rslt = pBED->GetFeature(CurFeatureID,szGene,szChrom,&GeneStart,&GeneEnd,NULL,&GeneStrand);
	if(Strand != 0 && ((Strand == 1 && GeneStrand != '+') || (Strand == 2 && GeneStrand != '-')))
		continue;
	GeneLen = GeneEnd - GeneStart;
	NumGOIDs = pAssocs->GetNumGOIDs(szGene);

	if(NumGOIDs > 0)
		{
		LocatedGOIDs++;
		pszConcatGOIDs = szGOIDs;
		*pszConcatGOIDs = '\0';
		CntGOIDs = 0;
		for(Cnt = 0; Cnt < NumGOIDs; Cnt++)
			{
			pszGOID = pAssocs->GetGOID(szGene,Cnt+1);
			if(pszGOID != NULL)
				{
				strcpy(pszConcatGOIDs,pszGOID);
				pszGOIDs[CntGOIDs++] = pszConcatGOIDs;
				pszConcatGOIDs += strlen(pszGOID) + 1;
				}
			}
		if(CntGOIDs)
			{
			TotbBkgndLen += (bBkgndLen ? GeneLen : 1);
			if((Rslt=AddCount(OntologyClass,bProp,false,bBkgndLen ? GeneLen : 1,CntGOIDs,pszGOIDs))!=eBSFSuccess)
				break;
			}
		}
	else
		NotLocatedGOIDs++;
	}
m_TotbBkgndLen = TotbBkgndLen;
return(Rslt);
}

// LocateGOID performs a binary search for specified GOID
// Returns a ptr to term identified by specified GOID or NULL if unable to locate term
// Uses binary search
tsGOTerm *
CGOTerms::LocateGOID(char *pszGOID)
{
tsGOTerm *pProbe;
char *pszID;
int *pID;
int Rslt;
if(pszGOID == NULL || *pszGOID == '\0')
	return(NULL);

int Lo,Mid,Hi;	// search limits
Lo = 0; Hi = m_GOTermcnt-1;
while(Hi >= Lo) {
	Mid = (Hi + Lo)/2;
	pProbe = &m_pGOTerms[Mid];
	pszID = ((tsGOID *)((char *)m_pGOIDs + pProbe->TermGOID.Idx))->Txt;
	Rslt = stricmp(pszGOID,pszID);
	if(Rslt < 0)	
		{
		Hi = Mid - 1;
		continue;
		}
	if(Rslt > 0)	
		{
		Lo = Mid + 1;
		continue;
		}
	return(pProbe);
	}

// could be that the sought GOID is an alternative identifier
// do a linear search, won't often be required so cost shouldn't be an issue
pProbe = m_pGOTerms;
for(Lo=0;Lo < m_GOTermcnt; Lo++, pProbe++)
	{
	if(!pProbe->NumAltTermIDs)
		continue;
	pID = &m_pAltGOIDs[pProbe->GOAltTermsIdx];
	for(Hi = 0; Hi < pProbe->NumAltTermIDs; Hi++,pID++)
		{
		pszID = ((tsGOID *)((char *)m_pGOIDs + *pID))->Txt;
		if(!stricmp(pszGOID,pszID))
			return(pProbe);
		}
	}
return(NULL);
}

// LocateRoot
// Returns a ptr to root GO:Term of specified term or 
// a) NULL if already root
// b) NULL if term marked as obsolete
tsGOTerm *
CGOTerms::LocateRoot(tsGOTerm *pCurTerm)
{
tsGOID *pCurID;
// already a root node with no parents?
if(pCurTerm == NULL || pCurTerm->bNoClass || !(pCurTerm->NumParents + pCurTerm->NumPartOfs))
	return(NULL);

// follow links back through 1st parent to root term
// all parents share common root so it doesn't matter which parent is used
do {
	if(pCurTerm->NumParents)
		pCurID = (tsGOID *)((char *)m_pGOIDs + m_pGOParentIDs[pCurTerm->GOParentIDsIdx]);
	else
		pCurID = (tsGOID *)((char *)m_pGOIDs + m_pGOParentIDs[pCurTerm->GOPartOfIDsIdx]);

	pCurTerm=LocateGOID(pCurID->Txt);
	}
while(pCurTerm != NULL && (pCurTerm->NumParents + pCurTerm->NumPartOfs));
return(pCurTerm);
}


// LocateNthParent
// Returns a ptr to GO:Term which is the Nth (1..n) parent relative to the specified term
// Returns NULL if the Nth parent not located
tsGOTerm *
CGOTerms::LocateNthParent(tsGOTerm *pCurTerm, int NthParent)
{
tsGOID *pCurID;
int NumParents;

if(pCurTerm == NULL || pCurTerm->bNoClass || NthParent < 1)
	return(NULL);
NumParents = pCurTerm->NumParents + pCurTerm->NumPartOfs;
if(!NumParents || NthParent > NumParents)
	return(NULL);
if(pCurTerm->NumParents)
	pCurID = (tsGOID *)((char *)m_pGOIDs + m_pGOParentIDs[pCurTerm->GOParentIDsIdx + NthParent - 1]);
else
	pCurID = (tsGOID *)((char *)m_pGOIDs + m_pGOParentIDs[pCurTerm->GOPartOfIDsIdx + NthParent - 1]);
return(LocateGOID(pCurID->Txt));
}

// LocateNthChild
// Returns a ptr to GO:Term which is the Nth (1..n) child relative to the specified term
// Returns NULL if the Nth child not located
tsGOTerm *
CGOTerms::LocateNthChild(tsGOTerm *pCurTerm, int NthChild)
{
tsGOID *pCurID;
if(pCurTerm == NULL || pCurTerm->bNoClass || NthChild < 1)
	return(NULL);
if(!pCurTerm->NumChildren || NthChild > pCurTerm->NumChildren)
	return(NULL);
pCurID = (tsGOID *)((char *)m_pGOIDs + m_pGOChildIDs[pCurTerm->GOChildIDsIdx + NthChild - 1]);
return(LocateGOID(pCurID->Txt));
}

// SetTermsOntologyClass
int
CGOTerms::SetOntologyClass4Terms(void)
{
tsGOTerm *pCurTerm;

if(m_FileHdr.POAnatomicalID || m_FileHdr.POdevID)
	{
	if(m_FileHdr.POAnatomicalID)
		{
		if((pCurTerm = GetGOTermByID(m_FileHdr.POAnatomicalID))==NULL)
			{
			AddErrMsg("CGOTerms::SetOntologyClass4Terms","Root term - %s - not located",pszRootTerms[1]);
			return(eBSFerrGOID);
			}
		SetChildClass(pCurTerm,eONTPlantAnatomical);
		}
	if(m_FileHdr.POdevID)
		{
		if((pCurTerm = GetGOTermByID(m_FileHdr.POdevID))==NULL)
			{
			AddErrMsg("CGOTerms::SetOntologyClass4Terms","Root term - %s - not located",pszRootTerms[1]);
			return(eBSFerrGOID);
			}
		SetChildClass(pCurTerm,eONTPlantDev);
		}
	}
if(m_FileHdr.MolTermID || m_FileHdr.BioTermID || m_FileHdr.CellTermID)
	{
	if((pCurTerm = GetGOTermByID(m_FileHdr.MolTermID))==NULL)
		{
		AddErrMsg("CGOTerms::SetOntologyClass4Terms","Root term - %s - not located",pszRootTerms[1]);
		return(eBSFerrGOID);
		}
	SetChildClass(pCurTerm,eONTMolecular);

	if((pCurTerm = GetGOTermByID(m_FileHdr.BioTermID))==NULL)
		{
		AddErrMsg("CGOTerms::SetOntologyClass4Terms","Root term - %s - not located",pszRootTerms[2]);
		return(eBSFerrGOID);
		}
	SetChildClass(pCurTerm,eONTBiological);

	if((pCurTerm = GetGOTermByID(m_FileHdr.CellTermID))==NULL)
		{
		AddErrMsg("CGOTerms::SetOntologyClass4Terms","Root term - %s - not located",pszRootTerms[0]);
		return(eBSFerrGOID);
		}
	SetChildClass(pCurTerm,eONTCellular);
	}
m_bGOTermsAvail = true;
return(eBSFSuccess);
}

void
CGOTerms::SetChildClass(tsGOTerm *pCurTerm,etOntologies OntologyClass)
{
tsGOTerm *pChildTerm;
int NthChild;

if(pCurTerm->bNoClass)
	return;

if(!pCurTerm->RootOntology)
	{
	pCurTerm->RootOntology = OntologyClass;
	switch(OntologyClass) {
		case eONTCellular:	// Cellular component
			m_GOCellTermcnt++;
			break;
		case eONTBiological:// Biological process
			m_GOBioTermcnt++;
			break;
		case eONTMolecular:	// Molecular function
			m_GOMolTermcnt++;
			break;

		case eONTPlantAnatomical:// number of terms classed as plant structure
			m_POAnatomyTermCnt++;
			break;
		case eONTPlantDev:	// number of terms classed as plant structure developmental growth
			m_PODevTermCnt++;
			break;

		default:
			break;
		}
	}

if(!pCurTerm->NumChildren)
	return;
NthChild = 1;
while((pChildTerm = LocateNthChild(pCurTerm,NthChild++))!=NULL)
	SetChildClass(pChildTerm,OntologyClass);
}


// sorts counts by by their TermID
int 
CGOTerms::SortCntsTermID(void)			
{
int Idx;
tsGOTermCnts *pCnt;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt < 1)
	return(eBSFerrInternal);
qsort(m_pGOTermCnts,m_GOTermCntsCnt,sizeof(tsGOTermCnts),SortTermCntsByTermID);
// need to readjust indexes in m_pGOTerms
pCnt = m_pGOTermCnts;
for(Idx = 0; Idx < m_GOTermCntsCnt; Idx++, pCnt++)
	m_pGOTerms[pCnt->TermID-1].GOTermCntsIdx = Idx+1;
return(eBSFSuccess);
}


// sorts counts by HGProbK
int 
CGOTerms::SortCntsHGProbK(void)			
{
int Idx;
tsGOTermCnts *pCnt;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt < 1)
	return(eBSFerrInternal);
qsort(m_pGOTermCnts,m_GOTermCntsCnt,sizeof(tsGOTermCnts),SortTermCntsByHGProbK);
// need to readjust indexes in m_pGOTerms
pCnt = m_pGOTermCnts;
for(Idx = 0; Idx < m_GOTermCntsCnt; Idx++, pCnt++)
	m_pGOTerms[pCnt->TermID-1].GOTermCntsIdx = Idx+1;
return(eBSFSuccess);
}

// sorts counts by HGWeightedProbK
int 
CGOTerms::SortCntsHGWeightedProbK(void)	
{
int Idx;
tsGOTermCnts *pCnt;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt < 1)
	return(eBSFerrInternal);
qsort(m_pGOTermCnts,m_GOTermCntsCnt,sizeof(tsGOTermCnts),SortTermCntsByHGWeightedProbK);
// need to readjust indexes in m_pGOTerms
pCnt = m_pGOTermCnts;
for(Idx = 0; Idx < m_GOTermCntsCnt; Idx++, pCnt++)
	m_pGOTerms[pCnt->TermID-1].GOTermCntsIdx = Idx+1;
return(eBSFSuccess);
}

// sorts counts by HGMTCProbK
int 
CGOTerms::SortCntsHGMTCProbK(void)		
{
int Idx;
tsGOTermCnts *pCnt;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt < 1)
	return(eBSFerrInternal);
qsort(m_pGOTermCnts,m_GOTermCntsCnt,sizeof(tsGOTermCnts),SortTermCntsByHGMTCProbK);
// need to readjust indexes in m_pGOTerms
pCnt = m_pGOTermCnts;
for(Idx = 0; Idx < m_GOTermCntsCnt; Idx++, pCnt++)
	m_pGOTerms[pCnt->TermID-1].GOTermCntsIdx = Idx+1;
return(eBSFSuccess);
}

// sorts counts by HGMTCWeightedProbK
int 
CGOTerms::SortCntsHGMTCWeightedProbK(void)
{
int Idx;
tsGOTermCnts *pCnt;
if(m_pGOTermCnts == NULL || m_GOTermCntsCnt < 1)
	return(eBSFerrInternal);
qsort(m_pGOTermCnts,m_GOTermCntsCnt,sizeof(tsGOTermCnts),SortTermCntsByHGMTCWeightedProbK);
// need to readjust indexes in m_pGOTerms
pCnt = m_pGOTermCnts;
for(Idx = 0; Idx < m_GOTermCntsCnt; Idx++, pCnt++)
	m_pGOTerms[pCnt->TermID-1].GOTermCntsIdx = Idx+1;
return(eBSFSuccess);
}


tsGOTerm *
CGOTerms::GetGOTermByID(int TermID) // returns ptr to term identified by TermID
{
if(TermID < 1 || TermID > m_GOTermcnt)
	return(NULL);
return(&m_pGOTerms[TermID-1]);
}

char *
CGOTerms::GetGOIDforTerm(tsGOTerm *pTerm)
{
char *pszID;
if(pTerm == NULL)
	return(NULL);
pszID = ((tsGOID *)((char *)m_pGOIDs + pTerm->TermGOID.Idx))->Txt;
return(pszID);
}

char *
CGOTerms::GetGONameforTerm(tsGOTerm *pTerm)
{
char *pszName;
if(pTerm == NULL)
	return(NULL);
pszName = GetGOVal(pTerm->TermNameIdx);
return(pszName);
}

int
CGOTerms::GetRootGOTermID(etOntologies Ontology)
{
switch(Ontology) {
	case eONTCellular:	// Cellular component
		return(m_FileHdr.CellTermID);

	case eONTBiological:	// Biological process
		return(m_FileHdr.BioTermID);

	case eONTMolecular:	// Molecular function
		return(m_FileHdr.MolTermID);

	case eONTPlantAnatomical:	// plant anatomical entity
		return(m_FileHdr.POAnatomicalID);

	case eONTPlantDev:	// plant developmental stage
		return(m_FileHdr.POdevID);

	default:
		break;
	}
return(0);
}

// returns root class ontology
char *
CGOTerms::RootOntology2Txt(etOntologies Ontology) 
{
const char *pszOntology;
	switch(Ontology) {
		case eONTCellular:	// Cellular component
			pszOntology = "Cellular Component";
			break;

		case eONTBiological:	// Biological process
			pszOntology = "Biological Process";
			break;

		case eONTMolecular:	// Molecular function
			pszOntology = "Molecular Function";
			break;

		case eONTPlantAnatomical:	// plant structure
			pszOntology = "Plant Anatomical Entity";
			break;

		case eONTPlantDev:	// plant developmental stage
			pszOntology = "Plant Structure Development Stage";
			break;

		default:
			pszOntology = "Undefined";
			break;
	}
return((char *)pszOntology);
}

// SwitchTermGOIDIDX2PTR
// Switch - map - m_pGOIDs relative offsets to pointers because SortTermsByGOID() depends on IDs being pointers
// Call SwitchTermGOIDPTR2IDX after sorting has completed
void 
CGOTerms::SwitchTermGOIDIDX2PTR(void)
{
tsGOTerm *pTerm;

/// if already ptrs, or no terms then easy return
if(m_bTermGOIDsAsPtrs)
	return;
if((pTerm = m_pGOTerms)==NULL || !m_GOTermcnt)
	{
	m_bTermGOIDsAsPtrs = true;
	return;
	}
for(int Cnt=0;Cnt < m_GOTermcnt; Cnt++,pTerm++)
	pTerm->TermGOID.ptr = (char *)m_pGOIDs + pTerm->TermGOID.Idx;
m_bTermGOIDsAsPtrs = true;
}


// SwitchTermGOIDPTR2IDX
// Switch - map - term IDs from pointers to offsets relative to m_pGOIDs
// pointers were set by SwitchTermGOIDIDX2PTR() because SortTermsByGOID() depends on IDs being pointers
void 
CGOTerms::SwitchTermGOIDPTR2IDX(void)
{
tsGOTerm *pTerm;

// if already Idx's, or no terms then easy return
if(!m_bTermGOIDsAsPtrs)
	return;
if((pTerm = m_pGOTerms)==NULL || !m_GOTermcnt)
	{
	m_bTermGOIDsAsPtrs = false;
	return;
	}
for(int Cnt=0;Cnt < m_GOTermcnt; Cnt++,pTerm++)
	pTerm->TermGOID.Idx = (int)((char *)pTerm->TermGOID.ptr - (char *)m_pGOIDs);
m_bTermGOIDsAsPtrs = false;
}


// SortTermsByGOID
// Sort GO:terms by their ID
int 
CGOTerms::SortTermsByGOID(const void *arg1, const void *arg2)
{
tsGOTerm *pEl1 = (tsGOTerm *)arg1;
tsGOTerm *pEl2 = (tsGOTerm *)arg2;
char *pszName1 = (char *)((tsGOID *)pEl1->TermGOID.ptr)->Txt;
char *pszName2 = (char *)((tsGOID *)pEl2->TermGOID.ptr)->Txt;
char c1 = tolower(*pszName1);
char c2 = tolower(*pszName2);
if(c1 < c2)
	return(-1);
if(c1 > c2)
	return(1);
return(stricmp(pszName1,pszName2));
}

// SortTermCntsByGOID
// Sort term cnts by their TermID
int 
CGOTerms::SortTermCntsByTermID(const void *arg1, const void *arg2)
{
tsGOTermCnts *pEl1 = (tsGOTermCnts *)arg1;
tsGOTermCnts *pEl2 = (tsGOTermCnts *)arg2;
if(pEl1->TermID < pEl2->TermID)
	return(-1);
if(pEl1->TermID > pEl2->TermID)
	return(1);
return(0);
}

// SortTermCntsByHGProbK
// Sort term cnts by their HGProbK (probabilites from hypergeometric cdf)
int 
CGOTerms::SortTermCntsByHGProbK(const void *arg1, const void *arg2)
{
tsGOTermCnts *pEl1 = (tsGOTermCnts *)arg1;
tsGOTermCnts *pEl2 = (tsGOTermCnts *)arg2;
if(pEl1->HGProbK < pEl2->HGProbK)
	return(-1);
if(pEl1->HGProbK > pEl2->HGProbK)
	return(1);
return(0);
}

// SortTermCntsByHGWeightedProbK
// Sort term cnts by their HGWeightedProbK (probabilites from weighted hypergeometric cdf)
int 
CGOTerms::SortTermCntsByHGWeightedProbK(const void *arg1, const void *arg2)
{
tsGOTermCnts *pEl1 = (tsGOTermCnts *)arg1;
tsGOTermCnts *pEl2 = (tsGOTermCnts *)arg2;
if(pEl1->HGWeightedProbK < pEl2->HGWeightedProbK)
	return(-1);
if(pEl1->HGWeightedProbK > pEl2->HGWeightedProbK)
	return(1);
return(0);
}

// SortTermCntsByHGMTCProbK
// Sort term cnts by their HGMTCProbK (MTC probabilites from hypergeometric cdf)
int 
CGOTerms::SortTermCntsByHGMTCProbK(const void *arg1, const void *arg2)
{
tsGOTermCnts *pEl1 = (tsGOTermCnts *)arg1;
tsGOTermCnts *pEl2 = (tsGOTermCnts *)arg2;
if(pEl1->HGMTCProbK < pEl2->HGMTCProbK)
	return(-1);
if(pEl1->HGMTCProbK > pEl2->HGMTCProbK)
	return(1);
return(0);
}

// SortTermCntsByHGMTCWeightedProbK
// Sort term cnts by their HGMTCWeightedProbK (MTC probabilites from hypergeometric cdf with gene weightings)
int 
CGOTerms::SortTermCntsByHGMTCWeightedProbK(const void *arg1, const void *arg2)
{
tsGOTermCnts *pEl1 = (tsGOTermCnts *)arg1;
tsGOTermCnts *pEl2 = (tsGOTermCnts *)arg2;
if(pEl1->HGMTCWeightedProbK < pEl2->HGMTCWeightedProbK)
	return(-1);
if(pEl1->HGMTCProbK > pEl2->HGMTCWeightedProbK)
	return(1);
return(0);
}




