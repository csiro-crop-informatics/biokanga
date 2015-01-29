#pragma once

const int cMaxComparedLCP = 300000;					// maximum length LCP compared, if longer then treat LCP+1 as mismatch

const int cDSSFREESIZE=5000;
const int cDSSMax_thresh=30;



class CDSsort
{
/* ------- node of blind trie -------- */ 
typedef struct nodex {
  INT32 skip;
  unsigned char key;  
  struct nodex  *down;      // first child
  struct nodex *right;      // next brother
	} node;


		void *m_pFreeArray[cDSSFREESIZE];
		node *m_pBufn;
		INT32 m_bufn_num;
		INT32 m_free_num;
		INT32 m_Aux_written;
		INT32 *m_pAux;
		node **m_ppStack;
		int m_Stack_size;
		INT32 m_Cmp_done;
		INT32  m_Text_size;           // size of input string 
		unsigned char  *m_pText;        // input string+ overshoot
		INT32  *m_pSuffixArray;       // suffix array
		unsigned char  *m_pEndOfText;   // m_pText+m_Text_size
		INT32  *m_pAnchorRank;        // rank (in the sorted suffixes of the  
									// anchor points (-1 if rank is unknown))
		UINT16  *m_pAnchorOffset;     // offset (wrt to the anchor) of the suffix
									// whose rank is in m_pAnchorRank. 
		INT32 m_NumAnchorPts;           // number of anchor points
		INT32 m_FtabArray[65537];   
		INT32 m_RunningOrderArray[256];

		unsigned char m_BucketRankedArray[65536];

		INT32 m_AnchorDist;                // distance between anchors
		INT32 m_DsVerbose;                // how verbose it the algorithm?
		INT32 m_DsWordSize;              // # of bytes in word in mkqs
		INT32 m_MkQsThresh;               // recursion limit for mk quicksort:
		INT32 m_MaxPseudoAnchorOffset; // maximum offset considered when 
										// searching a pseudo anchor
		INT32 m_B2gRratio;                // maximum ratio bucket_size/group_size
										// accepted for pseudo anchor_sorting
		INT32 m_UpdateAnchorRanks;      // if!=0 update anchor ranks when determining
										// rank for pseudo-sorting
		INT32 m_BlindSortRatio;         // blind sort is used for groups of size 
										// <= m_Text_size/m_BlindSortRatio

		INT32 m_Calls2HelpedSort;     
		INT32 m_Calls2AnchorSortForw;     
		INT32 m_Calls2AnchorSortBbackw;    
		INT32 m_Calls2PseudoAnchorSortForw;      
		INT32 m_Calls2DeepSort;     
		INT32 m_ShallowLimit;						// Max depth for shallow sorting
		unsigned char *m_pShallowTextLimit;         // m_pText+m_ShallowLimit
		INT32 m_CmpLeft;
		INT32 m_LCPAuxArray[1+cDSSMax_thresh];
		INT32 *m_pLCP; 



	void blind_ssort(INT32 *a, INT32 n, INT32 depth);
	node *find_companion(node *head, unsigned char *s);
	node *get_leaf(node *head);
	inline node *new_node__blind_ssort(void);
	void insert_suffix(node *h, INT32 suf, int n, unsigned char mmchar);
	void traverse_trie(node *h);
	inline INT32 get_lcp_unrolled(unsigned char *b1, unsigned char *b2, INT32 cmp_limit);
	INT32 compare_suffixes(INT32 suf1, INT32 suf2, INT32 depth);
	static int neg_integer_cmp(const void *a, const void *b);
	static int integer_cmp(const void *a, const void *b);
	inline INT32 cmp_unrolled_lcp(unsigned char *b1, unsigned char *b2);
	void qs_unrolled_lcp(INT32 *a, int n, int depth, int blind_limit);
	void deep_sort(INT32 *a, INT32 n, INT32 depth);
	void calc_running_order(void);
	void set_global_variables(void);
	int compute_overshoot(void);
	void pretty_putchar(int c);
	int scmp3(unsigned char *p, unsigned char *q, int *l, int maxl);
	void helped_sort(INT32 *a, int n, int depth);
	void pseudo_or_deep_sort(INT32 *a, INT32 n, INT32 depth);
	void pseudo_anchor_sort(INT32 *a,INT32 n,INT32 pseudo_anchor_pos, INT32 offset);
	void general_anchor_sort(INT32 *a, INT32 n,INT32 anchor_pos, INT32 anchor_rank, INT32 offset);
	INT32 get_rank(INT32 pos);
	INT32 get_rank_update_anchors(INT32 pos);
	void update_anchors(INT32 *a, INT32 n);
	INT32 split_group(INT32 *a, int n, int depth,int offset,INT32 pivot,int *first);
	void shallow_sort(INT32 *a, int n);
	inline void vecswap2(INT32 *a, INT32 *b, int n);
	inline INT32 *med3func(INT32 *a, INT32 *b, INT32 *c, unsigned char *text_depth);
	void shallow_mkq(INT32 *a, int n, unsigned char *text_depth);
	void shallow_mkq16(INT32 *a, int n, unsigned char *text_depth);
	void shallow_mkq32(INT32 *a, int n, unsigned char *text_depth);
	int check_global_variables(void);
	inline INT32 cmp_unrolled_shallow_lcp(unsigned char *b1, unsigned char *b2);
	void shallow_inssort_lcp(INT32 *a, INT32 n, unsigned char *text_depth);
	void free_node_mem(void);
public:
	CDSsort(void);
	~CDSsort(void);
	void ds_ssort(unsigned char *pToSort,	// string to sort suffixes on (NOTE: must have additional allocated memory past end for holding overshoot)
				  INT32 *pSuffixArray,	// sort into this prealloc'd suffix array
				  INT32 ToSortLen);		// string length (NOTE: excludes additional memory at end of pToSort allocated for overshoot)
	int init_ds_ssort(int adist = 500, int bs_ratio = 2000); // returns required overshoot memory size
};
