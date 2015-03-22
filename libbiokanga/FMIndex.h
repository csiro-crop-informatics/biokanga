#pragma once

const int cMinL2BlockSize  = 1;						// minimum accepted level 2 block size in 1K increments
const int cDfltL2BlockSize = 1;						// default level 2 block size in 1K increments
const int cMaxL2BlockSize  = 4;						// maximum accepted level 2 block size in 1K increments
const int cMinL1BlockSize  = cMinL2BlockSize * 4;	// minimum accepted level 1 block size (in 1K increments)
const int cDfltL1BlockSize = cMinL1BlockSize * 4;	// default level 1 block size (in 1K increments)
const int cMaxL1BlockSize  = 0x7fff;				// maximum accepted level 1 block size (in 1K increments)
const double cDfltMarkerFreq = 0.02;				// marker frequency (0.0 <= Freq <= 1.0)

const int BZ_N_GROUPS = 6;
const int BZ_MAX_CODE_LEN = 23;
const int BZ_MAX_ALPHA_SIZE = 258;
const int ALPHASIZE = 256;

		
#define EXT ".fmi"



class CFMIndex
{
	typedef struct TAG_sBucket_lev1 {
	  UINT32 *occ;             /* occ chars of compact alph in prev. superbuc */
	  UINT16 alpha_size;
	  unsigned char *bool_char_map;   /* boolean map of chars occurring in this superbucket */
	} bucket_lev1;
	 
	typedef struct TAG_sFm_index {
	  unsigned char *text;					/* input text */
	  unsigned char *oldtext;				/* A copy of the input text */
	  unsigned char *compress;				/* compress text */
	  unsigned char *bwt;					/* BWT of input text */
	  UINT32 *lf;					/* lf-mapping or suffix-array*/
	  UINT32 compress_size;			/* size of compressed file */
	  UINT32 text_size;				/* size of text */
	  bucket_lev1 *buclist_lev1; 	/* array of num_bucs buckets */
		
	  /* Info readed/writed on index prologue */
	  UINT32 bucket_size_lev1;		/* size of level 1 buckets */	
	  UINT32 bucket_size_lev2;		/* size of level 2 buckets */
	  UINT32 num_bucs_lev1;			/* number buckets lev 1 */
	  UINT32 num_bucs_lev2;			/* number buckets lev 2 */
	  UINT32 bwt_eof_pos;			/* position of EOF within BWT */
	  UINT16 alpha_size;				/* actual size of alphabet in input text */
	  UINT16 type_compression;		/* buckets lev 1 type of compression */
	  double freq;					/* frequency of marked chars */
	  UINT16 owner;					/* == 0 frees the text and alloc a new with overshoot */	
	  UINT16 compress_owner;         /* == 1 frees the compress */
	  UINT16 smalltext;				/* If == 1 stores plain text without compression */

	  /* Starting position (in byte) of each group of info */
	  UINT32 start_prologue_info_sb;	/* byte inizio info sui superbuckets */
	  UINT32 start_prologue_info_b;	/* byte inizio posizioni buckets */
	  unsigned char *start_prologue_occ;	/* byte inizio posizioni marcate */
	  UINT32 start_positions;        /* byte inizio posizioni marcate per build */
	  UINT32 *start_lev2;			/* starting position of each buckets in compr file */

	  /* Chars remap info */
	  UINT16 bool_char_map[ALPHASIZE];/* is 1 if char i appears on text */
	  unsigned char char_map[ALPHASIZE];	 /* cm[i]=j say that chars i is remapped on j */
	  unsigned char inv_char_map[ALPHASIZE]; /* icm[i] = j say that j was remapped on i */

	  /* Multiple locate. Chars substitution */
	  unsigned char specialchar;			/* carattere speciale che indica marcamento */
	  unsigned char subchar;				/* carattere sostituito dal carattere speciale */

	  // Da spostare in una struct working_space
	  /* Running temp info of actual superbucket and bucket */
	  UINT32 pfx_char_occ[ALPHASIZE];	/* i stores # of occ of chars 0.. i-1 in the text */
	  UINT16	bool_map_sb[ALPHASIZE]; 	/* info alphabet for superbucket to be read */
	  unsigned char inv_map_sb[ALPHASIZE]; 	/* inverse map for the current superbucket */
	  UINT16 alpha_size_sb;	  					/* current superbucket alphasize */
	  
	  unsigned char	bool_map_b[ALPHASIZE];		/* info alphabet for the bucket to be read */
	  unsigned char inv_map_b[ALPHASIZE];		/* inverse map for the current superbucket */
	  UINT16 alpha_size_b;     					/* actual size of alphabet in bucket */
	    
	  unsigned char mtf[ALPHASIZE];  		/* stores MTF-picture of bucket to be decompressed */
	  unsigned char *mtf_seq;				/* store bucket decompressed */
	  UINT32 occ_bucket[ALPHASIZE];  /* number chars occurences in the actual bucket needed by Mtf2 */
	  UINT16 int_dec_bits; 			/* log2(log2(text_size)) */
	  UINT16 log2textsize; 			/* int_log2(s.text_size-1) */
	  UINT16 var_byte_rappr;   		/* variable byte-length repr. (log2textsize+7)/8;*/
	  UINT32 num_marked_rows;		/* number of marked rows */
	  UINT32 bwt_occ[ALPHASIZE];     /* entry i stores # of occ of chars 0..i-1 */
	  UINT32 char_occ[ALPHASIZE];	/* (useful for mtf2) entry i stores # of occ of char i == bwt_occ[i]-bwt_occ[i-1] */
	  UINT16 skip;			/* 0 no marked pos, 1 all pos are marked, 2 only a char is marked*/
	  UINT16 sb_bitmap_size;			/* size in bytes of the bitmap of superbuckets */	
	  UINT32 occcharinf; 			/* Number occs of the chars  < index->specialchar */
	  UINT32 occcharsup; 			/* Number occs of the chars  <= index->specialchar */
	  
	  /* Needed by fm_build */
	  UINT32 *loc_occ;				/* Positions of marked rows */
	  
	} fm_index;





	/* Report rows from mulri_count */
	typedef struct TYPE_sMulti_count {
		UINT32 first_row; 			/* riga inizio occorrenza */
		UINT32 elements;   			/* numero occorrenze */	
		} multi_count;


	int __Fm_Verbose;

	fm_index *m_pIndex;							// contains complete context for compressed index

	UINT32 * __Num_Bytes;				/* numero byte letti/scritti */
	unsigned char * __MemAddress;				/* indirizzo della memoria dove scrivere */
	int __Bit_buffer_size;						/* number of unread/unwritten bits in Bit_buffer */
	UINT32 __pos_read;
	UINT32 __Bit_buffer;

	unsigned char *pmtf_start;
	int gmtflen;

	int allocated;
	int used;	/* var usate dalla count multipla */
	multi_count *lista;


	int count_row_mu (unsigned char * pattern, UINT32 len, UINT32 sp,UINT32 ep);
	inline void get_pos (UINT32 first_row, UINT32 element, UINT16 step, UINT32 * pos);
	int multi_locate (UINT32 sp, UINT32 element, UINT32 * positions);
	void preBmBc(unsigned char *x, int m, int bmBc[]);
	void suffixes(unsigned char *x, int m, int *suff);
	void preBmGs(unsigned char *x, int m, int bmGs[]);
	int fm_boyermoore(unsigned char * pattern, UINT32 length, UINT32 ** occ, UINT32 * numocc);
	int open_file(char * filename, unsigned char ** file, UINT32 * size);
	int fm_read_basic_prologue (void);
	int fm_multi_count (unsigned char * pattern, UINT32 len,multi_count ** list);

	int occ_all (UINT32 sp, UINT32 ep, UINT32 * occsp, UINT32 * occep,unsigned char * char_in);
	int get_info_sb (UINT32 pos, UINT32 * occ);
	int get_info_b (unsigned char ch, UINT32 pos, UINT32 *occ, int flag);
	unsigned char get_b_multihuf(UINT32 k, UINT32 * occ,int is_odd);
	inline void unmtf_unmap (unsigned char * mtf_seq, int len_mtf);
	int compress_bucket(unsigned char *in, UINT32 len, UINT16 alphasize);
	void mtf_string(unsigned char *in, unsigned char *out, UINT32 len, UINT16 mtflen);
	int fm_bwt_compress(void);
	void fm_unmtf(unsigned char *in, unsigned char *out, int length);
	int fm_bwt_uncompress(void);

	int fm_multihuf_compr (unsigned char * in, int len, int alpha_size);
	inline int decode_unary(void);
	int fm_multihuf_decompr (unsigned char * dest, int alpha_size, int limit);
	int fm_uncompress_bucket_multihuf (unsigned char * dest, int len, int alpha_size);

	unsigned char huf_len[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// coding and decoding
	int huf_code[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// coding
	int rfreq[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// coding
	int mtf_freq[BZ_MAX_ALPHA_SIZE];	// coding

	unsigned char huf_minLens[BZ_N_GROUPS];	// decoding
	int huf_limit[BZ_N_GROUPS][BZ_MAX_CODE_LEN];	// decoding
	int huf_base[BZ_N_GROUPS][BZ_MAX_CODE_LEN];	// decoding
	int huf_perm[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];	// decoding

	void fm_init_bit_writer (unsigned char * mem, UINT32 * pos_mem);
	void fm_bit_write (int n, UINT32 vv);
	int fm_bit_read (int n);
	void fm_bit_write24(int bits, UINT32 num);
	void fm_uint_write (UINT32 uu);
	UINT32 fm_uint_read (void);
	void fm_init_bit_reader (unsigned char * mem);
	void fm_bit_flush (void);
	UINT32 fm_integer_decode (unsigned short int headbits);

	int parse_options(char *optionz);
	int build_index(unsigned char *text, UINT32 length, char *build_options);
	int build_sa(void);
	void count_occ(void);
	int build_bwt(void);
	int compute_locations(void);
	int compute_info_superbuckets(void);
	int compute_info_buckets(void);
	void write_prologue(void);
	int compress_superbucket(UINT32 num);
	void write_susp_infos(void);
	void write_locations(void);
	int select_subchar(void);
	void dealloc(void);
	void dealloc_bucketinfo(void);
	int errore(int error);
	int save_index(char *filename);
	int fm_read_file(char *filename,		// file to read from 
				 unsigned char **textt,			// returned buffer allocated using malloc() containing contents of file + space for suffix sorting 
				 UINT32 *length);

	const char *error_index(int e);
	int int_log2(int u);
	int int_pow2(int u);

	UINT32 go_back(UINT32 row, UINT32 len, unsigned char *dest);
	UINT32 go_forw(UINT32 row, UINT32 len, unsigned char *dest);
	UINT32 fl_map(UINT32 row, unsigned char ch);
	unsigned char get_firstcolumn_char(UINT32 row);

	int fm_snippet(UINT32 row, UINT32 plen, UINT32 clen, unsigned char *dest, 
			   UINT32 *snippet_length);
	int read_prologue(void);
	int uncompress_data(void);
	int uncompress_superbucket(UINT32 numsb, unsigned char *out);
	int fm_compute_lf(void);
	int fm_invert_bwt(void);
	void free_unbuild_mem(void);
	int fm_unbuild(unsigned char ** text, UINT32 *length);

	void hbMakeCodeLengths(unsigned char *len,UINT32 *freq,UINT32 alphaSize,UINT32 maxLen);
	void hbAssignCodes(UINT32 *code,unsigned char *length,UINT32 minLen,UINT32 maxLen,UINT32 alphaSize);
	void hbCreateDecodeTables(UINT32 *limit,UINT32 *base,UINT32 *perm,unsigned char *length,UINT32 minLen,UINT32 maxLen,UINT32 alphaSize);
	

public:
	CFMIndex(void);
	~CFMIndex(void);
	int CreateIndex(char *pszInFile,		// create from contents of this file
					char *pszOutFile,		// write index into this file
					unsigned int *pText_len = NULL, // returned file content length before index created	
					unsigned int *pIndex_len = NULL, // created index length
					int bsl1 = cDfltL1BlockSize,	// level 1 block size (in 1K increments) will be forced to be multiple of bsl2
					int bsl2 = cDfltL2BlockSize,	// level 2 block size (in byte increments) will be forced to be multiple of 256
					double Freq = cDfltMarkerFreq);	// marker frequency (0.0 <= Freq <= 1.0)

	int load_index (char * filename);
	int extract(UINT32 from, UINT32 to, unsigned char **dest,UINT32 *snippet_length);
	int free_index (void);
	int display(unsigned char *pattern, UINT32 length, UINT32 nums, UINT32 *numocc, 
			unsigned char **snippet_text, UINT32 **snippet_len);
	int locate (unsigned char * pattern, UINT32 length, UINT32 ** occ,UINT32 * numocc);
	int count (unsigned char * pattern, UINT32 length, UINT32 * numocc);
	int get_length (UINT32 * length);

	int fm_build_config(double freq, UINT32 bsl1, UINT32 bsl2, UINT16 owner);
	int fm_build(unsigned char *text, UINT32 length);
	int index_size(UINT32 *size);


	int load_index_mem(unsigned char *compress, UINT32 size);		// loads compressed content from user supplied memory
	int save_index_mem(unsigned char *compress);					// saves compressed content into user supplied memory	

	char *GetErrText(int error);

};

