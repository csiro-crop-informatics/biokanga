#pragma once

// General Feature Format V1
// file record format:
// comments are preceded by a single '#', the rest of line following the '#' is a comment
// Each record consists of 
//<seqname><tab><source><tab><feature><tab><start><tab><end><tab><score><tab><strand><tab><frame>[<tab><group>][<comments>]<EOL>
//seqname	name of sequence or chromosome
//source    feature source
//feature   feature type name
//start		feature start locus (1..end)
//end		feature end locus (inclusive and >= start)
//score		floating point value, if no value then use '0'
//strand	'+','-', or '.' if no associated strand
//frame		'0','1','2' or '.' if 1st base not in a codon frame
//group		optional string-valued field that can be used as a name to group together a set of records

// General Feature Format V2 (http://www.sanger.ac.uk/resources/software/gff/spec.html)
// Also known as Gene Transfer Format - GTF
// file record format:
// optional metadata at start of file, lines starts with '##' followed by metadata
// common metadata:
// ##gff-version 2
// ##source-version <source> <version text>
// ##date <date>
//
// Each record consists of 
//<seqname><tab><source><tab><feature><tab><start><tab><end><tab><score><tab><strand><tab><frame>[<tab><attributelist>][<comments>]<EOL>
// Except for score (use '.' instead of '0' if no associated score), 1st 8 fields are same as fir GFF V1
//
// The attribute list consists of type/value pairs, each pair (except for last) is terminated by a ';' and is
// separated from the next attribute by exactly one space char ' '. Textual values should be bracketed by double quotes.
// There are two mandatory attributes:
// gene_id value			a globally unique identifier for genomic source of the sequence
// transcript_id value		a globally unique identifier for the predicted transcript
// Example GTF attributes: gene_id "Em:U62317.C22.6.mRNA"; transcript_id "Em:U62317.C22.6.mRNA"; exon_number 1

// General Feature Format V3 (http://www.sequenceontology.org/gff3.shtml)
// first line of file must be '##gff-version 3'
// strand may be '?' if strand unknown
// attribute/value pairs are in the format 'attribute=value', multiple pairs are separated ';'s.
// there can be multiple values associated to a single attribute, these values are comma separated.
// defined attributes:
// ID, Name, Alias, Parent, Target, Gap, Derives_from, Note, Dbxref, Ontology_term



const int cGFFMaxAttribs	= 50;	// max number of attributes
const int cGFFMaxValLen		= 1024;		// maximum length (chars) of any single attribute value supported
const int cGFFMaxLineLen  = (1024 + (cGFFMaxAttribs * cGFFMaxValLen));	// max number of characters which can be buffered

typedef enum TAG_etGFFGeneType {
	eGGFFany = 0,		// any GFF record
	eGGFFgene,			// only gene related records
	eGGFFtransposon,	// only transposon related records
	eGGFFmiRNA,			// only miRNA related records
	eGGFFsnoRNA,		// only snoRNA related records
	eGGFFtRNA,			// only tRNA related records
	eGGFFpseudogene,	// only pseudogene related records
	eGGFFncRNA,			// only ncRNA (also referenced as other_RNA)
	eGGFFplaceholder
} etGGFFGeneType;

#pragma pack(1)
typedef struct TAG_sGFFFields {
	int GFFVersion;									// GFF format is this version
	bool bIsMetadata;								// true if record line contains metadata
	char szRawLine[cGFFMaxLineLen];					// raw record line
	char szSeqName[cMaxDatasetSpeciesChrom];		// sequence name
	char szSeqSource[cMaxDatasetSpeciesChrom];		// sequence source
	char szFeature[cMaxDatasetSpeciesChrom];		// feature type
	UINT32 Start;									// feature start locus (1..End)
	UINT32 End;										// feature end locus (inclusive and >= start)
	bool bDfltScore;								// true if no score specified and default (0.0) is used
	float Score;									// if none then 0.0
	bool bDfltStrand;								// true if no strand specified and default '?' is used
	char Strand;									// '+', '-' or '?'
	bool bDfltFrame;								// true if no frame specified and default (0) is used
	int Frame;										// 0..2
	int AttribOfs;									// offset into szRawLine at which attributes start
	} tsGFFFields;

#pragma pack()


class CGFFFile  : public CErrorCodes
{
	FILE *m_pGFFStream;		// GFF file opened as a stream
	char m_szInFile[_MAX_PATH+1];	// input gff file

	bool m_bEOF;			// set TRUE when last record has been parsed

	int m_CurLineNum;		// current line number
	int m_CurLineLen;		// current line length
	int m_CurNumFields;		// number of fields parsed from current line
	int m_CurMaxFieldLen;	// maximum field length of any field in current line
	tsGFFFields m_Fields;	// parsed fields

	void Reset(void);
	int ParseLine(void);	// read next line and parse
	char *TrimWhitespace(char *pTxt); // inplace trim leading/trailing whitespace

public:
	CGFFFile(void);
	~CGFFFile(void);
	int Open(char *pszFileName);	// file to open
	int Close(void);

	int GetCurFields(void);				// get current number of attributes
	int GetCurFieldLen(void);			// get current maximum length of any attribute value
	int GetCurLineLen(void);			// get current line length

	int NextRecordOfType(etGGFFGeneType GGFFGeneType); // moves to next record of specified type

	int NextLine(void);					// move to next line in GFF file, skips blank lines and comment lines starting with comment char
	int GetLineNumber(void);			// get current line number
	bool IsMetadataLine(void);			// returns true if current line contains metadata
	int GetGFFversion(void);
	char *GetMetadata(void);
	char *GetRecord(void);				// returns complete current record
	char *GetSeqName(void);
	char *GetSource(void);
	char *GetFeature(void);
	int GetStart(void);
	int GetEnd(void);
	float GetScore(void);
	char GetStrand(void);
	char GetFrame(void);
	char *GetNoteValue(char *pszAttribute);

};
