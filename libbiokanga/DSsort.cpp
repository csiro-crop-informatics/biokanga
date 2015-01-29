#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "../hdrs/commhdrs.h"
#endif

#define True   ((bool)1)
#define False  ((bool)0)
#define Cmp_overshoot 16


#ifndef min
#define min(a, b) ((a)<=(b) ? (a) : (b))
#endif

#ifndef max
#define max(a, b) ((a)>=(b) ? (a) : (b))
#endif


#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a, b) ((a)>=(b) ? (a) : (b))


// constant and macro for marking groups
#define SETMASK (1 << 30)
#define CLEARMASK (~(SETMASK))
#define IS_SORTED_BUCKET(sb) (m_FtabArray[sb] & SETMASK)
#define BUCKET_FIRST(sb) (m_FtabArray[sb]&CLEARMASK)
#define BUCKET_LAST(sb) ((m_FtabArray[sb+1]&CLEARMASK)-1)
#define BUCKET_SIZE(sb) ((m_FtabArray[sb+1]&CLEARMASK)-(m_FtabArray[sb]&CLEARMASK))

#define BUFSIZE 1000
#define BIGFREQ(b) (m_FtabArray[((b)+1) << 8] - m_FtabArray[(b) << 8])


CDSsort::CDSsort(void)
{
set_global_variables();
}

CDSsort::~CDSsort(void)
{
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   string sorting routine based on a blind tree
   26-jun-01 ver 1.0
   03-jul-01 ver 1.1 (get rid of node header)
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */ 

/* ==================================================================
   comment: it is a little tricky how we handle the case in which 
   there are two strings s1 s2 such that s1 is a prefix of s2.
   (the correct ordering is that s1 preceeds lexicographically s2)
   We proceed as follows. We insert the strings in order of increasing 
   length so that s1 is inserted before s2. When s2 is inserted we put
   it in a leaf which is to the left of s1's leaf (obviously they have
   the same parent node). This is wrong acording to the alphabetic
   ordering but is done so that is there is a third string s3 which
   has s1 as a prefix we are certain that s3 meets s2's leaf and not s1's.
   When we traverse the trie to get the sorted string we check if there 
   are two sibling with the same key and if so we invert them 
   to get the correct ordering
   =================================================================== */


/* ****************************************************************
   routine for deep-sorting the suffixes a[0] ... a[n-1]
   knowing that they have a common prefix of length "depth"
  **************************************************************** */   
void 
CDSsort::blind_ssort(INT32 *a, INT32 n, INT32 depth)
{
  INT32 i,j,m_pLCP;
  INT32 aj;
  node nh, *root, *h;

  // ---- sort suffixes in order of increasing length
  qsort(a,n, sizeof(INT32), neg_integer_cmp);

  // --- skip suffixes which have already reached the end-of-text
  for(j=0;j<n;j++)
    if(a[j]+depth < m_Text_size)
      break;
  if(j>=n-1) return;  // everything is already sorted!

  // ------ init stack -------
  m_ppStack = (node **) malloc(n*sizeof(node *));
  if(m_ppStack==NULL) {
    fprintf(stderr,"Out of memory! (blind_ssort)\n");
    exit(1);
  }

  // ------- init root with the first unsorted suffix
  nh.skip = -1;   
  nh.right = NULL; 
  nh.down = (nodex *) (intptr_t)a[j]; 
  root = &nh;

  // ------- insert suffixes a[j+1] ... a[n-1]
  for(i=j+1;i<n;i++) {
    h=find_companion(root, m_pText+a[i]);
    assert(h->skip==-1);
    assert(m_Stack_size<=i-j);
    aj= (INT32)(intptr_t)h->down;
    assert(aj > a[i]);
    m_pLCP = compare_suffixes(aj,a[i],depth);
    insert_suffix(root, a[i], m_pLCP, m_pText[aj+m_pLCP]);
  }

  // ---- traverse the trie and get suffixes in lexicographic order  
  m_pAux=a;  m_Aux_written = j;
  traverse_trie(root);
  assert(m_Aux_written==n);
 
  free_node_mem();
  free(m_ppStack);
}

/* ***********************************************************************
   this function traverses the trie rooted at head following the string s. 
   Returns the leaf "corresponding" to the string s
   *********************************************************************** */
CDSsort::node *
CDSsort::find_companion(node *head, unsigned char *s)
{
  unsigned char c;
  node *p;
  int t;

  m_Stack_size = 0;                // init stack
  while(head->skip >= 0) {
    m_ppStack[m_Stack_size++] = head;
    t = head->skip;
    if(s+t>=m_pEndOfText)    // s[t] does not exist: mismatch 
      return get_leaf(head);
    c = s[t]; p = head->down;
  repeat:
    if(c==p->key) {              // found branch corresponding to c
      head = p;
      continue;
    }
    else if(c<p->key)            // no branch corresponding to c: mismatch
      return get_leaf(head);
    if((p=(p->right))==NULL)     // no other branches: mismatch
      return get_leaf(head);
    goto repeat;                 // look at next branch
  }
  m_ppStack[m_Stack_size++] = head;
  return head;
}


// this function returns a leaf below "head". 
// any leaf will do for the algorithm: we take the easiest to reach
CDSsort::node *
CDSsort::get_leaf(node *head)
{
  assert(head->skip>=0);

  do {
    head = head->down;
  } while(head->skip>=0);
  return head;
}



CDSsort::node *
CDSsort::new_node__blind_ssort(void)
{
  if(m_bufn_num-- == 0) {
    m_pBufn = (node *) malloc(BUFSIZE * sizeof(node));
    if(m_pBufn==NULL) {
      fprintf(stderr,"Out of mem (new_node1)\n"); exit(1);}
    m_pFreeArray[m_free_num++] = (void *) m_pBufn; 
    if(m_free_num>=cDSSFREESIZE) {
      fprintf(stderr,"Out of mem (new_node2)\n"); exit(1);}
   m_bufn_num = BUFSIZE-1;
  }
  return m_pBufn++;
}


/* *****************************************************
   insert a suffix in the trie rooted at *p.
   we know that the trie already contains a string
   which share the first n chars with suf
   ***************************************************** */
void 
CDSsort::insert_suffix(node *h, INT32 suf, int n, unsigned char mmchar)
{
  INT32 t;
  unsigned char c, *s;
  node *p, **pp;

  s = m_pText + suf;

  for(t=0;t<m_Stack_size;t++) {
    h=m_ppStack[t];
    if(h->skip<0 || h->skip>=n) break;
  }  
  
  assert(s[n]!=mmchar || h->skip==-1 || h->skip==n);

  // --------- insert a new node before node *h if necessary
  if(h->skip!=n) {
    p = new_node__blind_ssort();     // create and init new node
    p->key = mmchar;
    p->skip = h->skip;  // p inherits skip and children of *h
    p->down = h->down;   
    p->right = NULL;
    h->skip = n;
    h->down = p;        // now *h has p as the only child 
  }
  assert(h->skip==n);

  // -------- search the position of s[n] among *h offsprings
  c=s[n]; pp = &(h->down);
  while((*pp)!=NULL) {
    if((*pp)->key>=c) 
      break;
    pp = &((*pp)->right);
  }
  // ------- insert new node containing suf
  p = new_node__blind_ssort();
  p->skip = -1;
  p->key = c; 
  p->right = *pp; *pp = p;
  p->down = (nodex *)(intptr_t)suf;
  return;
}

/* ************************************************************
   this procedures traverse the trie in depth first order
   so that the suffixes (stored in the leaf) are recovered
   in lexicographic order
   ************************************************************ */
void 
CDSsort::traverse_trie(node *h)
{
  node *p, *nextp;

  if(h->skip < 0)
    m_pAux[m_Aux_written++] = (INT32)(intptr_t)h->down;
  else {
    p = h->down;
    assert(p!=NULL);
    do {
      nextp = p->right;
      if(nextp!=NULL) {
	assert(nextp->key>=p->key);
	// if there are 2 nodes with equal keys 
	// they must be considered in inverted order
	if(nextp->key==p->key) {
	  traverse_trie(nextp);
	  traverse_trie(p);
	  p = nextp->right;
	  continue;
	}
      }
      traverse_trie(p);
      p=nextp;
    } while(p!=NULL);
  }
}



/* ***********************************************************************
   Function to compute the m_pLCP of two strings originating from the *b1 and *b2
   the parameter is the length of s1 (which is shortest than s2)
   if s1 is a prefix of s2 we return the length of s1 -1
   The size of the unrolled loop must be at most equal to the costant 
   Cmp_overshoot defined in common.h
   the function return the result of the comparison (+ or -) and writes 
   in m_Cmp_done the number of comparisons done
   *********************************************************************** */ 
INT32 
CDSsort::get_lcp_unrolled(unsigned char *b1, unsigned char *b2, INT32 cmp_limit)
{
  INT32 cmp2do; 
  assert(b1 != b2);
  assert(cmp_limit > 0);



#ifdef USETHISOLDCODE
    unsigned char c1, c2;

  // execute blocks of 16 comparisons untill a difference
  // is found or we reach cmp_limit comparisons
  cmp2do = cmp_limit;
  do {
    // 1
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      break;}
    b1++; b2++; 
    // 2
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  1; break; }
    b1++; b2++; 
    // 3
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  2; break; }
    b1++; b2++; 
    // 4
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  3; break; }
    b1++; b2++; 
    // 5
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  4; break; }
    b1++; b2++; 
    // 6
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  5; break; }
    b1++; b2++; 
    // 7
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  6; break; }
    b1++; b2++; 
    // 8
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  7; break; }
    b1++; b2++; 
    // 9
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  8; break; }
    b1++; b2++; 
    // 10
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -=  9; break; }
    b1++; b2++; 
    // 11
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 10; break; }
    b1++; b2++; 
    // 12
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 11; break; }
    b1++; b2++; 
    // 13
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 12; break; }
    b1++; b2++; 
    // 14
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 13; break; }
    b1++; b2++; 
    // 15
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 14; break; }
    b1++; b2++; 
    // 16
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      cmp2do -= 15; break; }
    b1++; b2++; 

    cmp2do -= 16;
  } while(cmp2do>0);

assert(cmp2do >= 0);

  if(cmp_limit - cmp2do < cmp_limit)
    return cmp_limit-cmp2do;

  return cmp_limit-1;

#else
unsigned char *pT1 = b1;
unsigned char *pT2 = b2;
cmp2do = 0;
while(cmp2do < cmp_limit && (*pT1++ == *pT2++))
	cmp2do++;
if(cmp2do == cmp_limit)
	cmp2do--;
return(cmp2do);
#endif
} 



/* ************************************************************************
   this function returns the m_pLCP between suf1 and suf2 (that is returns n 
   such that suf1[n]!=suf2[n] but suf1[i]==suf2[i] for i=0..n-1
   However, it is possible that suf1 is a prefix of suf2 (not vice-versa
   because of the initial sorting of suffixes in order of descreasing length)
   in this case the function returns n=length(suf1)-1. So in this case 
   suf1[n]==suf2[n] (and suf1[n+1] does not exists). 
   ************************************************************************ */
 INT32 
CDSsort::compare_suffixes(INT32 suf1, INT32 suf2, INT32 depth)
{
  int limit;
  unsigned char *s1, *s2;

  assert(suf1 > suf2);
  s1  = m_pText + depth + suf1;
  s2  = m_pText + depth + suf2;
  limit = (int)(m_Text_size - suf1 - depth);
  return depth + get_lcp_unrolled(s1 ,s2, limit);
}


  
/* ****************************************************************** 
   comparison function used to sort suffixes in order of 
   increasing length. Since suffixes are represented by their offset
   in the array, we sort these offsets in order of decreasing length.
   ****************************************************************** */
int 
CDSsort::neg_integer_cmp(const void *a, const void *b)
{
  return *((INT32 *) b) -  *((INT32 *) a); 
}


// free memory used for trie nodes
void 
CDSsort::free_node_mem(void)
{
  int i;

  for(i=m_free_num-1;i>=0;i--) {
    assert(m_pFreeArray[i]!=NULL);
    free(m_pFreeArray[i]);
  }
  // clear counters
  m_bufn_num=m_free_num=0;
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   deep.c

   "deep" sorting routines
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */



/* ***********************************************************************
   Function to compare two strings originating from the *b1 and *b2
   The size of the unrolled loop must be at most equal to the costant 
   Cmp_overshoot defined in common.h
   the function return the result of the comparison (+ or -) and writes 
   in m_Cmp_done the number of successfull comparisons done
   *********************************************************************** */ 

INT32 
CDSsort::cmp_unrolled_lcp(unsigned char *b1, unsigned char *b2)
{
  unsigned char c1, c2;
  assert(b1 != b2);
  m_Cmp_done=0;

#ifdef USETHISOLDCODE
  // execute blocks of 16 comparisons untill a difference
  // is found or we run out of the string 
  do {
    // 1
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 2
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  1; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 3
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  2; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 4
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  3; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 5
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  4; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 6
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  5; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 7
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  6; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 8
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  7; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 9
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  8; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 10
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done +=  9; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 11
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done += 10; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 12
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done += 11; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 13
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done += 12; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 14
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done += 13; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 15
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done += 14; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 16
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_Cmp_done += 15; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 

    m_Cmp_done += 16;

  } while( b1<m_pEndOfText && b2<m_pEndOfText);
assert(b1 <= m_pEndOfText && b2 <= m_pEndOfText);
  //return (b2-m_pText) - (b1-m_pText);   // we have  b2>b1 <=> *b2<*b1
  return (int)(ptrdiff_t)(b2 - b1);
#else
unsigned char *pT1 = b1;
unsigned char *pT2 = b2;
unsigned int Max2Cmp = (unsigned int)min((ptrdiff_t)(m_pEndOfText - b1),(ptrdiff_t)(m_pEndOfText - b2));
while(Max2Cmp-- && ((c1=*pT1++) == (c2=*pT2++)))
	m_Cmp_done++;
if(c1 != c2)
	return(c1 > c2 ? 1 : -1);
return (int)(ptrdiff_t)(b2 - b1);
#endif
} 

/* **************************************************************
   ternary quicksort (seward-like) with m_pLCP information
   ************************************************************** */
#define STACK_SIZE 100
#define Swap(i,j) {tmp=a[i]; a[i]=a[j]; a[j]=tmp;}
#define Pushd(x,y,z) {stack_lo[sp]=x; stack_hi[sp]=y; stack_d[sp]=z; sp++;}
#define Popd(x,y,z)  {sp--; x=stack_lo[sp]; y=stack_hi[sp]; z=stack_d[sp];} 

void 
CDSsort::qs_unrolled_lcp(INT32 *a, int n, int depth, int blind_limit)
{ 
  unsigned char *text_depth, *text_pos_pivot;
  INT32 stack_lo[STACK_SIZE];
  INT32 stack_hi[STACK_SIZE];
  INT32 stack_d[STACK_SIZE];
  INT32 sp,r,r3,med,tmp;
  INT32 i, j, lo, hi,ris,lcp_lo,lcp_hi;

  // ----- init quicksort --------------
  r=sp=0;
  Pushd(0,n-1,depth);

  // ----- repeat untill stack is empty ------
  while (sp > 0) {
    assert ( sp < STACK_SIZE );
    Popd(lo,hi,depth);
    text_depth = m_pText+depth;

    // --- use shellsort for small groups
    if(hi-lo<blind_limit) { 
       blind_ssort(a+lo,hi-lo+1,depth);
       continue;
    }

    /* Random partitioning. Guidance for the magic constants 
       7621 and 32768 is taken from Sedgewick's algorithms
       book, chapter 35.
    */
    r = ((r * 7621) + 1) % 32768;
    r3 = r % 3;
    if (r3 == 0) med = lo; else
    if (r3 == 1) med = (lo+hi)>>1; else
                 med = hi;

    // --- partition ----
    Swap(med,hi);  // put the pivot at the right-end
    text_pos_pivot=text_depth+a[hi];
    i=lo-1; j=hi; 
    lcp_lo=lcp_hi=INT_MAX;
    while(1) {
      while(++i<hi) {
		ris=cmp_unrolled_lcp(text_depth+a[i], text_pos_pivot);
        if(ris>0) 
			{
			if(m_Cmp_done < lcp_hi) 
				lcp_hi=m_Cmp_done; 
			break;
			} 
		else 
			if(m_Cmp_done < lcp_lo) 
				lcp_lo=m_Cmp_done;
		}
      while(--j>lo) {
	ris=cmp_unrolled_lcp(text_depth+a[j], text_pos_pivot);
    if(ris<0) 
		{ 
		if(m_Cmp_done < lcp_lo) 
			lcp_lo=m_Cmp_done; 
		break; 
		}
	else 
		if(m_Cmp_done < lcp_hi) 
			lcp_hi=m_Cmp_done;
      }
   if (i >= j) break; 
   Swap(i,j);
   }
    Swap(i,hi);  // put pivot at the middle

    // ---- testing ---------
    assert(lcp_lo<INT_MAX || i==lo);
    assert(lcp_hi<INT_MAX || i==hi);

    // --------- insert subproblems in stack; smallest last
    if(i-lo < hi-i) {
      Pushd(i+1,hi,depth+lcp_hi);
      if(i-lo>1) Pushd(lo,i-1,depth+lcp_lo);
    }
    else {
      Pushd(lo,i-1,depth+lcp_lo);
      if(hi-i>1) Pushd(i+1,hi,depth+lcp_hi);
    }
  }
}



/* ****************************************************************
   routine for deep-sorting the suffixes a[0] ... a[n-1]
   knowing that they have a common prefix of length "depth"
  **************************************************************** */   
void 
CDSsort::deep_sort(INT32 *a, INT32 n, INT32 depth)
{
  int blind_limit;

  m_Calls2DeepSort++;    
  assert(n>1);    // test to discover useless calls

  blind_limit=m_Text_size/m_BlindSortRatio;
  if(n<=blind_limit)
    blind_ssort(a,n,depth);  // small_group
  else 
    qs_unrolled_lcp(a,n,depth,blind_limit);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   ds.c
   deep-shallow sorting algorithm. these routines are taken mainly by 
   Seward's d_copyEQ_u12.c 
   Buckets are sorted one at a time; when  bucket Y is completely sorted
   (except for YY) "pointer copying" is used to sort small buckets of 
   the form XY, X=A..Z
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */ 


/* ------------------------------------------------------------------------
   The use of m_pAnchorRank[] and m_pAnchorOffset is the following: 
   m_pAnchorRank[i] is either -1 or contains the rank (in the list of
   the sorted suffixes of the suffix starting at position 
     i*m_AnchorDist + m_pAnchorOffset[i].
   Initially m_pAnchorRank[i] = -1 and m_pAnchorOffset[i]=m_AnchorDist,
   then, if a suffix in position t (i*m_AnchorDist <= t < (i+1)*m_AnchorDist)
   appears to be in a large group which is sorted, the rank of
   t is stored in m_pAnchorRank[i], and the value t-(i*m_AnchorDist) 
   is written to m_pAnchorOffset[i]. Both vaulues can be later updated,
   but the value in m_pAnchorOffset[i] can only decrease, so no further
   changes are done when m_pAnchorOffset[i] is = 0. The invariant is:
   if m_pAnchorRank[i]>=0 then 
       m_pSuffixArray[m_pAnchorRank[i]]=i*m_AnchorDist+m_pAnchorOffset[i]
   -------------------------------------------------------------------------*/
 

  


/* ************************************************************
   This is the main deep/shallow suffix sorting routines
   It divides the suffixes in buckets according to the 
   first two characters. Some of the buckets are then sorted 
   calling shallow_sort(). When all buckets of kind ax (pToSort!=a) are 
   sorted, we use this ordering to sort all suffixes in the 
   buckets ya (for any y including y=a).
    ************************************************************* */
// NOTE: Expects additional writeable memory at end of pToSort to have been allocated for holding 'overshoot area'
// The size of this 'overshoot area' is expected to be that returned by a previous call to init_ds_ssort()

void 
CDSsort::ds_ssort(unsigned char *pToSort,	// string to sort suffixes on. NOTE: Expected to have space allocated at end for holding overshoot
				  INT32 *pSuffixArray,	// sort into this prealloc'd suffix array
				  INT32 ToSortLen)		// string length
{
#ifdef USEOLDSIZES
  int overshoot;
  INT32  i, j, ss, sb;
  INT32 k;
  unsigned char  c1, c2;
  bool   bigDone[256];
  INT32  copyStart[256];
  INT32  copyEnd  [256];
  INT32  numQSorted = 0;
#else
  int overshoot;
  size_t  i, j, ss, sb;
  size_t k;
  unsigned char  c1, c2;
  bool   bigDone[256];
  size_t  copyStart[256];
  size_t  copyEnd  [256];
  size_t  numQSorted = 0;
#endif


  // ------ set some global variables ------
  m_pText=pToSort;
  m_Text_size=ToSortLen;
  m_pSuffixArray = pSuffixArray;
  m_pEndOfText = m_pText + m_Text_size;
  // ------ fill overshoot area
  overshoot = compute_overshoot();
  for(i=ToSortLen;i<((size_t)ToSortLen+overshoot);i++) m_pText[i]=0; 

  // ------ init array containing positions of anchors
  if(m_AnchorDist==0) {
    m_NumAnchorPts=0; m_pAnchorRank=NULL; m_pAnchorOffset=NULL;
  }
  else {
    m_NumAnchorPts = 2 + (ToSortLen-1)/m_AnchorDist;  // see comment for helped_sort() 
    m_pAnchorRank = (INT32 *) malloc(m_NumAnchorPts*sizeof(INT32));
    m_pAnchorOffset = (UINT16 *) malloc(m_NumAnchorPts*sizeof(UINT16));
    if(!m_pAnchorRank || !m_pAnchorOffset) {
      fprintf(stderr, "malloc failed (ds_sort)\n");
      exit(1);
    }
    for(i=0;i<(size_t)m_NumAnchorPts;i++) {
      m_pAnchorRank[i]= -1;               // pos of anchors is initially unknown
      m_pAnchorOffset[i] = m_AnchorDist;   // maximum possible value
    }
  }

  // ---------- init m_FtabArray ------------------
  for (i = 0; i <= 65536; i++) m_FtabArray[i] = 0;
  c1 = m_pText[0];
  for (i = 1; i <= (size_t)m_Text_size; i++) {
    c2 = m_pText[i];
    m_FtabArray[(c1 << 8) + c2]++;
    c1 = c2;
  }
  for (i = 1; i <= 65536; i++) 
	  m_FtabArray[i] += m_FtabArray[i-1];

  // -------- sort suffixes considering only the first two chars 
  c1 = m_pText[0];
  for (int i2 = 0; i2 < m_Text_size; i2++) {
    c2 = m_pText[i2+1];
    j = (c1 << 8) + c2;
    c1 = c2;
    m_FtabArray[j]--;
    m_pSuffixArray[m_FtabArray[j]] = i2;
  }

  /* decide on the running order */
  calc_running_order();
  for (i = 0; i < 256; i++) bigDone[i] = False;

   /* Really do the suffix sorting */
  for (i = 0; i <= 255; i++) {

    /*--
      Process big buckets, starting with the least full.
      --*/
    ss = m_RunningOrderArray[i];
    if(m_DsVerbose>2)
      fprintf(stderr,"group %3d;  size %d\n",(int)ss,(int)(BIGFREQ(ss)&CLEARMASK)); 

    /*--
      Complete the big bucket [ss] by sorting
      any unsorted small buckets [ss, j].  Hopefully
      previous pointer-scanning phases have already
      completed many of the small buckets [ss, j], so
      we don't have to sort them at all.
      --*/
    for (j = 0; j <= 255; j++) {
      if (j != ss) {
	sb = (ss << 8) + j;
	if ( ! (m_FtabArray[sb] & SETMASK) ) {
	  INT32 lo = m_FtabArray[sb]   & CLEARMASK;
	  INT32 hi = (m_FtabArray[sb+1] & CLEARMASK) - 1;
	  if (hi > lo) {
	    if (m_DsVerbose>2)
	      fprintf(stderr,"sorting [%02x, %02x], done %d "
			"this %d\n", (unsigned int)ss, (unsigned int)j, (unsigned int)numQSorted, hi - lo + 1 );
	    shallow_sort(m_pSuffixArray+lo, hi-lo+1);
            #if 0
	    check_ordering(lo, hi);
            #endif
	    numQSorted += ( hi - lo + 1 );
	  }
	}
	m_FtabArray[sb] |= SETMASK;
      }
    }
    assert (!bigDone[ss]);
    // ------ now order small buckets of type [xx,ss]  --------
    {
      for (j = 0; j <= 255; j++) {
	copyStart[j] =  m_FtabArray[(j << 8) + ss]     & CLEARMASK;
	copyEnd  [j] = (m_FtabArray[(j << 8) + ss + 1] & CLEARMASK) - 1;
      }
      // take care of the virtual -1 char in position m_Text_size+1
      if(ss==0) {
	k=m_Text_size-1;
	c1 = m_pText[k];
	if (!bigDone[c1])
	  m_pSuffixArray[ copyStart[c1]++ ] = (INT32)k;
      }
      for (j = m_FtabArray[ss << 8] & CLEARMASK; j < copyStart[ss]; j++) {
	k = m_pSuffixArray[j]-1; 
	if (k < 0) 
		continue;  
	c1 = m_pText[k];
	if (!bigDone[c1])
	  m_pSuffixArray[ copyStart[c1]++ ] = (INT32)k;
      }
      for (j = (m_FtabArray[(ss+1) << 8] & CLEARMASK) - 1; j > copyEnd[ss]; j--) {
	k = m_pSuffixArray[j]-1; 
	if (k < 0) 
		continue;
	c1 = m_pText[k];
	if (!bigDone[c1]) 
	  m_pSuffixArray[ copyEnd[c1]-- ] = (INT32)k;
      }
    }
    assert (copyStart[ss] - 1 == copyEnd[ss]);
    for (j = 0; j <= 255; j++) m_FtabArray[(j << 8) + ss] |= SETMASK;
    bigDone[ss] = True;
  }
  if (m_DsVerbose) {
    fprintf(stderr, "\t %d pointers, %d sorted, %d scanned\n",
	      m_Text_size, (unsigned int)numQSorted, (unsigned int)(m_Text_size - numQSorted ));
    fprintf(stderr, "\t %d calls to helped_sort\n",m_Calls2HelpedSort);      
    fprintf(stderr, "\t %d calls to anchor_sort (forward)\n",
	    m_Calls2AnchorSortForw);      
    fprintf(stderr, "\t %d calls to anchor_sort (backward)\n",
	    m_Calls2AnchorSortBbackw);      
    fprintf(stderr, "\t %d calls to pseudo_anchor_sort (forward)\n",
    	    m_Calls2PseudoAnchorSortForw);      
    fprintf(stderr, "\t %d calls to deep_sort\n",m_Calls2DeepSort);      
  }
  // ---- done! ---------------------------------------- 
  free(m_pAnchorOffset);
  free(m_pAnchorRank);
}



/* ****************************************************************
   compute running =(sorting) order for big buckets: start with 
   the least full and proceed to the largest one.
   The sorting is done using shellsort
   **************************************************************** */ 

void 
CDSsort::calc_running_order ( void )
{
   INT32 i, j;
   for (i = 0; i <= 255; i++) m_RunningOrderArray[i] = i;

   {
      INT32 vv;
      INT32 h = 1;
      do h = 3 * h + 1; while (h <= 256);
      do {
         h = h / 3;
         for (i = h; i <= 255; i++) {
            vv = m_RunningOrderArray[i];
            j = i;
            while ( BIGFREQ(m_RunningOrderArray[j-h]) > BIGFREQ(vv) ) {
               m_RunningOrderArray[j] = m_RunningOrderArray[j-h];
               j = j - h;
               if (j <= (h - 1)) goto zero;
            }
            zero:
            m_RunningOrderArray[j] = vv;
         }
      } while (h != 1);
   }
}



/* *******************************************************************
   globals.c
   Ver 1.0   14-oct-02
   This file contains the definition of the global variables 
   which can be defined by the user + some relate procedures
   ******************************************************************* */

/* *******************************************************************
   procedure to be called by external program before calling ds_ssort()
   using this procedure external programs can choose
   the parameters m_AnchorDist and m_BlindSortRatio.
   The procedure returns 0 if something goes wrong, otherwise 
   it returns the overshhot, that is the amount of extra space
   required at the end of the array contanining the text
   ******************************************************************** */
int  
CDSsort::init_ds_ssort(int adist, 
					   int bs_ratio)
{
  set_global_variables();
  m_AnchorDist = adist;
  m_BlindSortRatio=bs_ratio;
  m_ShallowLimit =  m_AnchorDist + 50;
  if(check_global_variables())
    return 0;
  return compute_overshoot();
}


// set default values for the global variables
void 
CDSsort::set_global_variables(void)
{
m_BlindSortRatio=2000;
m_AnchorDist = 500;
m_ShallowLimit = 550;
m_DsVerbose = 0;
m_DsWordSize = 4;
m_MkQsThresh=20; 
m_MaxPseudoAnchorOffset=0;
m_B2gRratio=1000;
m_UpdateAnchorRanks=0;
m_bufn_num=0;
m_free_num=0;
m_Calls2HelpedSort=0;     
m_Calls2AnchorSortForw=0;     
m_Calls2AnchorSortBbackw=0;    
m_Calls2PseudoAnchorSortForw=0;      
m_Calls2DeepSort=0; 
m_pLCP=&m_LCPAuxArray[1]; 
}

// check if the global variables passed as parameters
// are in the valid range
int 
CDSsort::check_global_variables(void)
{
  if((m_AnchorDist<100) && (m_AnchorDist!=0)) {
    fprintf(stderr,"Anchor distance must be 0 or greater than 99\n");
    return 1;
  }
  if(m_AnchorDist>65535) {
    fprintf(stderr,"Anchor distance must be less than 65536\n");
    return 1;
  }
  if(m_ShallowLimit<2) {
    fprintf(stderr,"Illegal limit for shallow sort\n");
    return 1;
  }
  if(m_MkQsThresh<0 || m_MkQsThresh>cDSSMax_thresh) {
    fprintf(stderr,"Illegal m_MkQsThresh parameter!\n");
    return 1;
  }
  if(m_BlindSortRatio<=0) {
    fprintf(stderr,"blind_sort ratio must be greater than 0!\n");
    return 1;
  }
  return 0;
}


// compute the amount of extra memory required 
// for the input text 
int 
CDSsort::compute_overshoot(void)
{
  return 9+(m_ShallowLimit+Cmp_overshoot);
}

// ----- this function prints any char in a readable form
void 
CDSsort::pretty_putchar(int c)
{
  
  if(c>=32 && c<127)      // printable char
    printf("  %c", c);
  else if(c=='\n')
    printf(" \\n");        // \n
  else if(c=='\t')
    printf(" \\t");        // \t
  else     
    printf(" %02x", c);      // print hex code
}


int 
CDSsort::scmp3(unsigned char *p, unsigned char *q, int *l, int maxl)
{
   int i;
   i = 0;
   while (maxl>0 && *p==*q) {
      p++; q++; i++;
      maxl--;
   }
   *l = i;
   if (maxl>0) return *p-*q;
   return (int)(ptrdiff_t)(q-p);
}








/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   helped.c  
   these routines sort a group of strings which have a common prefix
   using induced sorting (if possible) or a deep-sort routine
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */



// macro to compute the bucket for the suffix
// starting at pos. Note that since pos is evaluated twice
// it should be an expression without side-effects
#define Get_small_bucket(pos) ((m_pText[pos]<<8) + m_pText[pos+1])


/* *****************************************************************
   This procedure sort the strings a[0] ... a[n-1] with the help of an
   anchor. The real sorting is done by the procedure
   anchor_sort(). Here we choose the anchor.  The parameter depth is
   the number of chars that a[0] ... a[n-1] are known to have in
   common (thus a direct comparison among a[i] and a[j] should start
   from position depth) Note that a[] is a subsection of the sa therefore
   a[0] ... a[n-1] are starting position of suffixes
   For every a[i] we look at the anchor a[i]/m_AnchorDist and the one 
   after that. This justifies the definition of m_NumAnchorPts (the size of
   Anchor_ofset[] and m_pAnchorRank[] defined in ds_sort()) as
     m_NumAnchorPts = 2 + (n-1)/m_AnchorDist    
   ***************************************************************** */
void 
CDSsort::helped_sort(INT32 *a, int n, int depth)
{ 
 
  INT32 i, curr_sb, diff, toffset, aoffset;
  INT32 text_pos, anchor_pos, anchor, anchor_rank;
  INT32 min_forw_offset, min_forw_offset_buc, max_back_offset;
  INT32 best_forw_anchor, best_forw_anchor_buc, best_back_anchor; 
  INT32 forw_anchor_index, forw_anchor_index_buc, back_anchor_index;

  m_Calls2HelpedSort++;          // update count
  if(n==1) goto done_sorting;    // simplest case: only one string

  // if there are no anchors use pseudo-anchors or deep_sort
  if(m_AnchorDist==0) {
    pseudo_or_deep_sort(a, n, depth);
    return;
  }

  // compute the current bucket 
  curr_sb = Get_small_bucket(a[0]);

  // init best anchor variables with illegal values
  min_forw_offset = min_forw_offset_buc = INT_MAX;
  max_back_offset = INT_MIN;
  best_forw_anchor = best_forw_anchor_buc = best_back_anchor = -1; 
  forw_anchor_index = forw_anchor_index_buc = back_anchor_index = -1;
  // look at the anchor preceeding each a[i]
  for(i=0;i<n;i++) {
    text_pos = a[i];
    // get anchor preceeding text_pos=a[i]
    anchor = text_pos/m_AnchorDist;
    toffset = text_pos % m_AnchorDist;  // distance of a[i] from anchor
    aoffset = m_pAnchorOffset[anchor];   // distance of sorted suf from anchor 
    if(aoffset<m_AnchorDist) {          // check if it is a "sorted" anchor
      diff = aoffset - toffset;
      assert(diff!=0);
      if(diff>0) {     // anchor <=  a[i] < (sorted suffix)
	if(curr_sb!=Get_small_bucket(text_pos+diff)) {
	  if(diff<min_forw_offset) {
	    min_forw_offset = diff;
	    best_forw_anchor = anchor;
	    forw_anchor_index = i;
	  }
	}
	else {  // the sorted suffix belongs to the same bucket of a[0]..a[n-1]
	  if(diff<min_forw_offset_buc) {
	    min_forw_offset_buc = diff;
	    best_forw_anchor_buc = anchor;
	    forw_anchor_index_buc = i;
	  }
	}
      }
      else {          // diff<0 =>  anchor <= (sorted suffix) < a[i]
	if(diff>max_back_offset) {
	  max_back_offset = diff;
	  best_back_anchor = anchor;
	  back_anchor_index = i;
	}
	// try to find a sorted suffix > a[i] by looking at next anchor
	aoffset = m_pAnchorOffset[++anchor];
	if(aoffset<m_AnchorDist) {
	  diff = m_AnchorDist + aoffset - toffset;
	  assert(diff>0);
	  if(curr_sb!=Get_small_bucket(text_pos+diff)) {
	    if(diff<min_forw_offset) {
	      min_forw_offset = diff;
	      best_forw_anchor = anchor;
	      forw_anchor_index = i;
	    }
	  } else {
	    if(diff<min_forw_offset_buc) {
	      min_forw_offset_buc = diff;
	      best_forw_anchor_buc = anchor;
	      forw_anchor_index_buc = i;
	    }
	  }
	}
      }
    }
  }
  // ------ if forward anchor_sort is possible, do it! --------	    
  if(best_forw_anchor>=0 && min_forw_offset<depth-1) {
    m_Calls2AnchorSortForw++;
    assert(min_forw_offset<2*m_AnchorDist);
    anchor_pos = a[forw_anchor_index] + min_forw_offset;
    anchor_rank = m_pAnchorRank[best_forw_anchor];
    assert(m_pSuffixArray[anchor_rank]==anchor_pos);
    general_anchor_sort(a,n,anchor_pos,anchor_rank,min_forw_offset);
    goto done_sorting;
  }
  // ------ if backward anchor_sort is possible do it! ---------
  if(best_back_anchor>=0) {
    unsigned char *T0, *Ti; int j;

    assert(max_back_offset>-m_AnchorDist && max_back_offset<0);
    // make sure that the offset is legal for all a[i]
    for(i=0;i<n;i++) {
      if(a[i]+max_back_offset<0) 
	goto fail;                    // illegal offset, give up
    }
    // make sure that a[0] .. a[n-1] are preceded by the same substring
    T0 = m_pText + a[0];
    for(i=1;i<n;i++) {
      Ti = m_pText + a[i];
      for(j=max_back_offset; j<= -1; j++)
	if(T0[j]!=Ti[j]) goto fail;   // mismatch, give up
    }
    // backward anchor sorting is possible
    m_Calls2AnchorSortBbackw++;
    anchor_pos = a[back_anchor_index] + max_back_offset;
    anchor_rank = m_pAnchorRank[best_back_anchor];
    assert(m_pSuffixArray[anchor_rank]==anchor_pos);
    general_anchor_sort(a,n,anchor_pos,anchor_rank,max_back_offset);
    goto done_sorting;
  }
 fail:
  // ----- try forward anchor_sort with anchor in the same bucket
  if(best_forw_anchor_buc>=0 && min_forw_offset_buc<depth-1) {
    int equal,lower,upper;

    assert(min_forw_offset_buc<2*m_AnchorDist);
    anchor_pos = a[forw_anchor_index_buc] + min_forw_offset_buc;
    anchor_rank = m_pAnchorRank[best_forw_anchor_buc];
    assert(m_pSuffixArray[anchor_rank]==anchor_pos);

    // establish how many suffixes can be sorted using anchor_sort()
    equal=split_group(a,n,depth,min_forw_offset_buc,
                                forw_anchor_index_buc,&lower);
    if(equal==n) {
      m_Calls2AnchorSortForw++;
      general_anchor_sort(a,n,anchor_pos,anchor_rank,min_forw_offset_buc);
    }
    else {
      //  -- a[0] ... a[n-1] are split into 3 groups: lower, equal, upper
      upper = n-equal-lower;
      assert(upper>=0);
      // printf("Warning! lo=%d eq=%d up=%d a=%pToSort\n",lower,equal,upper,(int)a);
      // sort the equal group 
      m_Calls2AnchorSortForw++;
      if(equal>1)
	general_anchor_sort(a+lower,equal,anchor_pos,anchor_rank,
			    min_forw_offset_buc);

      // sort upper and lower groups using deep_sort
      if(lower>1) pseudo_or_deep_sort(a,lower,depth);
      if(upper>1) pseudo_or_deep_sort(a+lower+equal,upper,depth);
    }       // end if(equal==n) ... else
    goto done_sorting;
  }         // end hard case

  // ---------------------------------------------------------------
  // If we get here it means that everything failed
  // In this case we simply deep_sort a[0] ... a[n-1]
  // ---------------------------------------------------------------
  pseudo_or_deep_sort(a, n, depth);
 done_sorting:
  // -------- update m_pAnchorRank[], m_pAnchorOffset[] ------- 
  if(m_AnchorDist>0) update_anchors(a, n);
}
  


/* *******************************************************************
   try pseudo_anchor sort or deep_sort
   ******************************************************************** */
void 
CDSsort::pseudo_or_deep_sort(INT32 *a, INT32 n, INT32 depth)
{
  INT32 offset, text_pos, sb, pseudo_anchor_pos, max_offset, size;
 
  // ------- search for a useful pseudo-anchor -------------
  if(m_MaxPseudoAnchorOffset>0) {

    max_offset = min(depth-1,m_MaxPseudoAnchorOffset);
    text_pos = a[0];
    for(offset=1;offset<max_offset;offset++) {
      pseudo_anchor_pos = text_pos+offset;
      sb = Get_small_bucket(pseudo_anchor_pos);
      // check if pseudo_anchor is in a sorted bucket
      if(IS_SORTED_BUCKET(sb)) {
	size=BUCKET_SIZE(sb);                     // size of group
	if(size>m_B2gRratio*n) continue;            // discard large groups 
	// sort a[0] ... a[n-1] using pseudo_anchor
	pseudo_anchor_sort(a,n,pseudo_anchor_pos,offset);
	m_Calls2PseudoAnchorSortForw++;        // update count
	return;
      }
    }
  }
  deep_sort(a,n,depth);
}

/* ********************************************************************
   this routine sorts the suffixes a[0] ... a[n-1] using the fact that
   in their common prefix, after offset characters, there is a 
   suffix which is in an already sorted bucket. This suffix is called
   a pseudo anchor since it is used essentially as an anchor, but
   it is not in an anchor position (=position multiple of m_AnchorDist)
   ******************************************************************** */

void 
CDSsort::pseudo_anchor_sort(INT32 *a,INT32 n,INT32 pseudo_anchor_pos, INT32 offset)
{
  INT32 pseudo_anchor_rank;

  // ---------- compute rank ------------
  if(m_UpdateAnchorRanks!=0 && m_AnchorDist>0)
    pseudo_anchor_rank = get_rank_update_anchors(pseudo_anchor_pos);
  else
    pseudo_anchor_rank = get_rank(pseudo_anchor_pos);
  // ---------- check rank --------------
  assert(m_pSuffixArray[pseudo_anchor_rank]==pseudo_anchor_pos);
  // ---------- do the sorting ----------
  general_anchor_sort(a,n,pseudo_anchor_pos,pseudo_anchor_rank,offset);
}


/* ********************************************************
   macros for marking integers: works assuming integers have 
   at least 32 bit and that the 32nd bit is not used
   This simply means that the text size can be at most 2GB
   ********************************************************* */
#define MARKER (1<<31)
#define MARK(i) {                \
  assert(( m_pSuffixArray[i]&MARKER) == 0);  \
  (m_pSuffixArray[i] |= MARKER);             \
}
#define ISMARKED(i) (m_pSuffixArray[i] & MARKER)
#define UNMARK(i) (m_pSuffixArray[i] &= ~MARKER)

/* ********************************************************************
   This routines sorts a[0] ... a[n-1] using the fact that
   in their common prefix, after offset characters, there is a 
   suffix whose rank is known. In this routine we call this suffix anchor
   (and we denote its position and rank with anchor_pos and anchor_rank 
   respectively) but it is not necessarily an anchor (=does not necessarily 
   starts at position multiple of m_AnchorDist) since this function is
   called by pseudo_anchor_sort().
   The routine works by scanning the suffixes before and after the anchor
   in order to find (and mark) those which are suffixes of a[0] ... a[n-1].
   After that, the ordering of a[0] ... a[n-1] is derived with a sigle
   scan of the marked suffixes.
   ******************************************************************** */
void 
CDSsort::general_anchor_sort(INT32 *a, INT32 n, 
                         INT32 anchor_pos, INT32 anchor_rank, INT32 offset)
{
  INT32 sb, lo, hi;
  INT32 curr_lo, curr_hi, to_be_found, i,j;
  INT32 item; 
  void *ris;

  assert(m_pSuffixArray[anchor_rank]==anchor_pos);
  /* ---------- get bucket of anchor ---------- */
  sb = Get_small_bucket(anchor_pos);
  lo = BUCKET_FIRST(sb);
  hi = BUCKET_LAST(sb);
  assert(sb==Get_small_bucket(a[0]+offset));
  // ------ sort pointers a[0] ... a[n-1] as plain integers
  qsort(a,n, sizeof(INT32), integer_cmp);

  // ------------------------------------------------------------------
  // now we scan the bucket containing the anchor in search of suffixes
  // corresponding to the ones we have to sort. When we find one of
  // such suffixes we mark it. We go on untill n sfx's have been marked 
  // ------------------------------------------------------------------
  curr_hi = curr_lo = anchor_rank;

  // the anchor must correspond to a suffix to be sorted
  #if DEBUG
  item = anchor_pos-offset;
  assert(bsearch(&item,a,n,sizeof(INT32), integer_cmp));
  #endif

  MARK(curr_lo);
  // scan suffixes preceeding and following the anchor
  for(to_be_found=n-1;to_be_found>0; ) {
    // invariant: the next positions to check are curr_lo-1 and curr_hi+1
    assert(curr_lo > lo || curr_hi < hi);
    while (curr_lo > lo) {
      item = m_pSuffixArray[--curr_lo]-offset;
      ris = bsearch(&item,a,n,sizeof(INT32), integer_cmp);
      if(ris)	{MARK(curr_lo); to_be_found--;}
      else	break;
    }
    while (curr_hi < hi) {
      item = m_pSuffixArray[++curr_hi]-offset;
      ris = bsearch(&item,a,n,sizeof(INT32), integer_cmp);
      if(ris)	{MARK(curr_hi); to_be_found--;}
      else      break;
    }
  }
  // sort a[] using the marked suffixes
  for(j=0, i=curr_lo;i<=curr_hi;i++) 
    if(ISMARKED(i)) {
      UNMARK(i);
      a[j++] = m_pSuffixArray[i] - offset;
    }
  assert(j==n);  // make sure n items have been sorted
}

/* ********************************************************************
   compute the rank of the suffix starting at pos.
   It is required that the suffix is in an already sorted bucket
   ******************************************************************** */
INT32 
CDSsort::get_rank(INT32 pos)
{
  INT32 sb, lo, hi, j;

  sb = Get_small_bucket(pos);  
  if(!IS_SORTED_BUCKET(sb)) {
    fprintf(stderr,"Illegal call to get_rank! (get_rank1)\n");
    exit(1);
  }
  lo = BUCKET_FIRST(sb);
  hi = BUCKET_LAST(sb);
  for(j=lo;j<=hi;j++) 
    if(m_pSuffixArray[j]==pos) return j;
  fprintf(stderr,"Illegal call to get_rank! (get_rank2)\n");
  exit(1);
  return 1;   // so that the compiler does not complain
}

/* ********************************************************************
   compute the rank of the suffix starting at pos. At the same time
   check if the rank of the suffixes in the bucket containing pos
   can be used to update some entries in m_pAnchorOffset[] and m_pAnchorRank[]
   It is required that the suffix is in an already sorted bucket   
   ******************************************************************** */
INT32 
CDSsort::get_rank_update_anchors(INT32 pos)
{
  INT32 sb, lo, hi, j, toffset, aoffset, anchor, rank;

  assert(m_AnchorDist>0);
  // --- get bucket and verify it is a sorted one
  sb = Get_small_bucket(pos);  
  if(!(IS_SORTED_BUCKET(sb))) {
    fprintf(stderr,"Illegal call to get_rank! (get_rank_update_anchors)\n");
    exit(1);
  }
  // --- if the bucket has been already ranked just compute rank; 
  if(m_BucketRankedArray[sb]) return get_rank(pos);
  // --- rank all the bucket 
  m_BucketRankedArray[sb]=1;
  rank = -1;
  lo = BUCKET_FIRST(sb);
  hi = BUCKET_LAST(sb);
  for(j=lo;j<=hi;j++) {  
    // see if we can update an anchor
    toffset = m_pSuffixArray[j]%m_AnchorDist;
    anchor  = m_pSuffixArray[j]/m_AnchorDist;
    aoffset = m_pAnchorOffset[anchor];  // dist of sorted suf from anchor 
    if(toffset<aoffset) {
      m_pAnchorOffset[anchor] = toffset;
      m_pAnchorRank[anchor] = j;
    }
    // see if we have found the rank of pos, if so store it in rank
    if(m_pSuffixArray[j]==pos) {
      assert(rank==-1); rank=j;
    }
  }
  assert(rank>=0);
  return rank;
}


/* ****************************************************************** 
   comparison function used to sort pointers as if they were integers
   if a pointer does not correspond to an INT32 we must change the
   function accordingly 
   ****************************************************************** */
int 
CDSsort::integer_cmp(const void *a, const void *b)
{
  return *((INT32 *) a) -  *((INT32 *) b); 
}

/* ****************************************************************
   given a SORTED array of suffixes a[0] .. a[n-1]
   updates m_pAnchorRank[] and m_pAnchorOffset[]
   **************************************************************** */
void 
CDSsort::update_anchors(INT32 *a, INT32 n)
{
  INT32 i,anchor,toffset,aoffset,text_pos;

  assert(m_AnchorDist>0);
  for(i=0;i<n;i++) {
    text_pos = a[i];
    // get anchor preceeding text_pos=a[i]
    anchor = text_pos/m_AnchorDist;
    toffset = text_pos % m_AnchorDist;     // distance of a[i] from anchor
    aoffset = m_pAnchorOffset[anchor];  // dist of sorted suf from anchor 
    if(toffset<aoffset) {
      m_pAnchorOffset[anchor] = toffset;
      m_pAnchorRank[anchor] = (INT32)((ptrdiff_t)(a - m_pSuffixArray) + i);
      assert(m_pSuffixArray[m_pAnchorRank[anchor]]==
	     anchor*m_AnchorDist+m_pAnchorOffset[anchor]);
    }
  }
}
 


/* *******************************************************************
   This function takes as input an array a[0] .. a[n-1] of suffixes
   which share the first "depth" chars. "pivot" in an index in 0..n-1
   and offset and integer>0. The function splits a[0] .. a[n-1]
   into 3 groups: first the suffixes which are smaller than a[pivot],
   then those which are equal to a[pivot] and finally those which are 
   greater than a[pivot]. Here, smaller, equal, larger refer to 
   a lexicographic ordering limited to the first depth+offest chars
   (since the first depth chars are equal we only look at the chars
   in position depth, depth+1, ... depth+offset-1).
   The function returns the number "num" of suffixes equal to a[pivot],
   and stores in *first the first of these suffixes. So at the end 
   the smaller suffixes are in a[0] ... a[first-1],
   the equal suffixes in a[first] ... a[first+num-1],
   the larger suffixes in a[first+num] ... a[n-1]
   The splitting is done using a modified mkq()
   ******************************************************************* */
#define swap2(a, b) { t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(i) (*(*(i) + text_depth))
 
INT32 
CDSsort::split_group(INT32 *a, int n, int depth,int offset,INT32 pivot,int *first)
{
  int r, partval;
  INT32 *pa, *pb, *pc, *pd, *pa_old, *pd_old, pivot_pos, t;
  unsigned char *text_depth,*text_limit;

  // --------- initialization ------------------------------------
  pivot_pos = a[pivot];       // starting position in T[] of pivot
  text_depth = m_pText+depth;
  text_limit = text_depth+offset;

  // -------------------------------------------------------------
  // In the following for() loop:
  // [pa ... pd] is the current working region, 
  // pb moves from pa towards pd 
  // pc moves from pd towards pa
  // -------------------------------------------------------------
  pa = a; pd = a + n-1;

  for(  ; pa!=pd && (text_depth<text_limit); text_depth++) {
    assert(pa<pd);
    // ------ the pivot char is m_pText[pivot_pos+depth] where 
    // depth = text_depth-m_pText. This is text_depth[pivot_pos]
    partval = text_depth[pivot_pos];
    // ----- partition ------------ 
    pb = pa_old = pa; 
    pc = pd_old = pd; 
    for (;;) {
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
            if (r == 0) { swap2(pa, pb); pa++; }
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
            if (r == 0) { swap2(pc, pd); pd--; }
            pc--;
        }
        if (pb > pc) break;
        swap2(pb, pc);
        pb++;
        pc--;
    }
    r = (int)min((ptrdiff_t)(pa-pa_old), (ptrdiff_t)(pb-pa)); 
	vecswap2(pa_old,  pb-r, r);
    r = (int)min((ptrdiff_t)(pd-pc), (ptrdiff_t)(pd_old-pd)); 
	vecswap2(pb, pd_old+1-r, r);
    // ------ compute new boundaries ----- 
    pa = pa_old + (pb-pa);     // there are pb-pa chars < partval
    pd = pd_old - (pd-pc);     // there are pd-pc chars > partval

  }
  *first=(int)((ptrdiff_t)(pa-a));        // index in a[] of the first suf. equal to pivot
  assert(pd-pa>=0);
  return (int)((ptrdiff_t)(pd-pa)+1);     // return number of suffixes equal to pivot
}










/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   shallow.c  
   This is the multikey quicksort from bentley-sedgewick modified 
   so that it stops recursion when depth reaches  m_ShallowLimit 
   (that is when two or more suffixes have m_ShallowLimit chars in common).
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */


#define UNROLL 1                   // if !=0 partially unroll shallow_mkq

// ***** entry point for shallow sort routines *****
void 
CDSsort::shallow_sort(INT32 *a, int n) 
{ 
  // init global variables
  m_pShallowTextLimit = m_pText + m_ShallowLimit;
  // call multikey quicksort
  // skip 2 chars since suffixes come from the same bucket 
  switch(m_DsWordSize) {
  case(1): shallow_mkq(a, n, m_pText+2); break;
  case(2): shallow_mkq16(a, n, m_pText+2); break;
  case(4): shallow_mkq32(a, n, m_pText+2); break;
  default:
    fprintf(stderr,
	    "Invalid word size for mkqs (%d) (shallow_sort)\n",m_DsWordSize);
    exit(1);
  }     
}


/* =======================================================
   auxiliary procedures and macro for bentley-sedgewick's
   multikey quicksort
   ======================================================= */
void 
CDSsort::vecswap2(INT32 *a, INT32 *b, int n)
{   while (n-- > 0) {
        INT32 t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

#define swap2(a, b) { t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(i) (*(*(i) + text_depth))

INT32 *
CDSsort::med3func(INT32 *a, INT32 *b, INT32 *c, unsigned char *text_depth)
{   int va, vb, vc;
    if ((va=ptr2char(a)) == (vb=ptr2char(b)))
        return a;
    if ((vc=ptr2char(c)) == va || vc == vb)
        return c;       
    return va < vb ?
          (vb < vc ? b : (va < vc ? c : a ) )
        : (vb > vc ? b : (va < vc ? a : c ) );
}
#define med3(a, b, c) med3func(a, b, c, text_depth)


/* ********************************************************
   recursive multikey quicksort from Bentley-Sedgewick
   stops when text_depth reaches Shallow_depth_limit 
   that is when we have found that the current set of strings
   have m_ShallowLimit chars in common
   ******************************************************** */
void 
CDSsort::shallow_mkq(INT32 *a, int n, unsigned char *text_depth)
{
  int d, r, partval;
  INT32 *pa, *pb, *pc, *pd, *pl, *pn, t;
  INT32 *pm;
  unsigned char *next_depth;

  // ---- On small arrays use insertions sort
  if (n < m_MkQsThresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

  // ----------- choose pivot --------------
 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // -------- partition -----------------
  for (;;) {
    while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
      if (r == 0) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
      if (r == 0) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }

#if UNROLL
  if(pa>pd) {
    // all values were equal to partval: make it simpler
    if( (next_depth = text_depth+1) >= m_pShallowTextLimit) {
      helped_sort(a, n, (int)((ptrdiff_t)(next_depth-m_pText)));
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  // partition a[] into the values smaller, equal, and larger that partval
  pn = a + n;
  r = (int)min((ptrdiff_t)(pa-a), (ptrdiff_t)(pb-pa));    vecswap2(a,  pb-r, r);
  r = (int)min((ptrdiff_t)(pd-pc), ((ptrdiff_t)(pn-pd)-1)); vecswap2(pb, pn-r, r);
  // --- sort smaller strings -------
  if ((r = (int)(ptrdiff_t)(pb-pa)) > 1)
    shallow_mkq(a, r, text_depth);
  // --- sort strings starting with partval -----
  if( (next_depth = text_depth+1) < m_pShallowTextLimit)
    shallow_mkq(a + r, (int)((ptrdiff_t)(pa-pd) + n -1), next_depth);
  else 
    helped_sort(a + r, (int)((ptrdiff_t)(pa-pd)+n-1), (int)(ptrdiff_t)(next_depth-m_pText));
  if ((r = (int)(ptrdiff_t)(pd-pc)) > 1)
    shallow_mkq(a + n-r, r, text_depth);
}



/* ************** 16 *************** */
#define ptr2char16(i) (getword16(*(i) + text_depth))
#define getword16(s) ((unsigned)((*(s) << 8) | *((s)+1)))

void 
CDSsort::shallow_mkq16(INT32 *a, int n, unsigned char *text_depth)
{
  int d, r, partval;
  INT32 *pa, *pb, *pc, *pd, *pl, *pm, *pn;
  INT32 t;
  unsigned char *next_depth;

  // ---- On small arrays use insertions sort
  if (n < m_MkQsThresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

  // ----------- choose pivot --------------
 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char16(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // -------- partition -----------------
  for (;;) {
    while (pb <= pc && (r = ptr2char16(pb)-partval) <= 0) {
      if (r == 0) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char16(pc)-partval) >= 0) {
      if (r == 0) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }
#if UNROLL
  if(pa>pd) {
    // all values were equal to partval: make it simpler
    if( (next_depth = text_depth+2) >= m_pShallowTextLimit) {
      helped_sort(a, n, (int)(next_depth-m_pText));
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  // partition a[] into the values smaller, equal, and larger that partval
  pn = a + n;
  r = (int)min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = (int)min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  // --- sort smaller strings -------
  if ((r = (int)(pb-pa)) > 1)
    shallow_mkq16(a, r, text_depth);
  // --- sort strings starting with partval -----
  if( (next_depth = text_depth+2) < m_pShallowTextLimit)
    shallow_mkq16(a + r, (int)((ptrdiff_t)(pa-pd)+n-1), next_depth);
  else 
    helped_sort(a + r, (int)((ptrdiff_t)(pa-pd)+n-1), (int)(next_depth-m_pText));
  if ((r = (int)(pd-pc)) > 1)
    shallow_mkq16(a + n-r, r, text_depth);
}


/* *************** 32 **************** */
#define ptr2char32(i) (getword32(*(i) + text_depth))
#define getword32(s) ((unsigned)( (*(s) << 24) | ((*((s)+1)) << 16) \
                                  | ((*((s)+2)) << 8) | (*((s)+3)) ))
void 
CDSsort::shallow_mkq32(INT32 *a, int n, unsigned char *text_depth)
{
  UINT32 partval, val;
  INT32 *pa, *pb, *pc, *pd, *pl, *pm, *pn, t, d;
  int r;
  unsigned char *next_depth;

  // ---- On small arrays use insertions sort
  if (n < m_MkQsThresh) {
    shallow_inssort_lcp(a, n, text_depth);
    return;
  }

  // ----------- choose pivot --------------
 repeat:
  pl = a;
  pm = a + (n/2);
  pn = a + (n-1);
  if (n > 30) { // On big arrays, pseudomedian of 9
    d = (n/8);
    pl = med3(pl, pl+d, pl+2*d);
    pm = med3(pm-d, pm, pm+d);
    pn = med3(pn-2*d, pn-d, pn);
  }
  pm = med3(pl, pm, pn);
  swap2(a, pm);
  partval = ptr2char32(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  // -------- partition -----------------
  for (;;) {
    while (pb <= pc &&  (val=ptr2char32(pb)) <=  partval) {
      if (val == partval) { swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (val=ptr2char32(pc)) >= partval) {
      if (val == partval) { swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    swap2(pb, pc);
    pb++;
    pc--;
  }
#if UNROLL
  if(pa>pd) {
    // all values were equal to partval: make it simpler
    if( (next_depth = text_depth+4) >= m_pShallowTextLimit) {
      helped_sort(a, n, (int)(next_depth-m_pText));
      return;
    }
    else {
      text_depth = next_depth;
      goto repeat;
    }
  }
#endif
  // partition a[] into the values smaller, equal, and larger that partval
  pn = a + n;
  r = (int)min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
  r = (int)min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  // --- sort smaller strings -------
  if ((r = (int)(pb-pa)) > 1)
    shallow_mkq32(a, r, text_depth);
  // --- sort strings starting with partval -----
  if( (next_depth = text_depth+4) < m_pShallowTextLimit)
    shallow_mkq32(a + r, (int)((ptrdiff_t)(pa-pd)+n-1), next_depth);
  else 
    helped_sort(a + r, (int)((ptrdiff_t)(pa-pd)+n-1), (int)(next_depth-m_pText));
  if ((r = (int)(pd-pc)) > 1)
    shallow_mkq32(a + n-r, r, text_depth);
}



/* >>>>>>>>>>>>>>>>>>>>>> insertion sort routines >>>>>>>>>>>>>>>>>>>
   This insertion sort routines sorts the suffixes a[0] .. a[n-1]
   which have a common prexif of length text_depth-m_pText.
   The comparisons are done going at most at depth m_ShallowLimit;
   suffixes which have m_ShallowLimit chars in common are sorted using 
   helped_sort().
   This inserion_sort keeps trak of the m_pLCP in order to speed up 
   the sorting.
  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

/* ***********************************************************************
   Function to compare two strings originating from the *b1 and *b2
   The size of the unrolled loop must be at most equal to the costant 
   Cmp_overshoot defined in common.h
   When the function is called m_CmpLeft must contain the maximum number of 
   comparisons the algorithm can do before returning 0 (equal strings)
   At exit m_CmpLeft has been decreased by the # of comparisons done   
   *********************************************************************** */ 

INT32 
CDSsort::cmp_unrolled_shallow_lcp(unsigned char *b1, unsigned char *b2)
{

  unsigned char c1, c2;
  assert(b1 != b2);
  assert(m_CmpLeft > 0);

#ifdef USETHISOLDCODE
  // execute blocks of 16 comparisons until a difference
  // is found or we run out of the string 
  do {
    // 1
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 2
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  1; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 3
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  2; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 4
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  3; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 5
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  4; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 6
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  5; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 7
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  6; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 8
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  7; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 9
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  8; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 10
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -=  9; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 11
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -= 10; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 12
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -= 11; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 13
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -= 12; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 14
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -= 13; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 15
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -= 14; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // 16
    c1 = *b1; c2 = *b2;
    if (c1 != c2) {
      m_CmpLeft -= 15; return ((UINT32)c1 - (UINT32)c2); }
    b1++; b2++; 
    // if we have done enough comparisons the strings are considered equal
    m_CmpLeft -= 16;
    if(m_CmpLeft<=0) return 0;
    // assert( b1<m_pEndOfText && b2<m_pEndOfText);
  } while(1);
  //return (b2-m_pText) - (b1-m_pText);   // we have  b2>b1 <=> *b2<*b1
  return (int)(b2 - b1);

#else
unsigned char *pT1 = b1;
unsigned char *pT2 = b2;
while(m_CmpLeft && ((c1=*pT1++) == (c2=*pT2++)))
	m_CmpLeft--;
if(!m_CmpLeft)
	return(0);
return(c1 > c2 ? 1 : -1);
#endif
} 

/* *****************************************************************
   this is the insertion sort routine called by multikey-quicksort
   for sorting small groups.
   During insertion sort the comparisons are done calling 
   cmp_unrolled_shallow_lcp() and two strings are equal if the coincides 
   for m_ShallowLimit characters.
   After this first phase we sort groups of "equal_string" using 
   helped_sort(). 
   Usage of m_pLCP. 
   For i=1,...n-1 let m_pLCP[i] denote the m_pLCP between a[i] and a[i+1].
   assume a[0] ... a[j-1] are already ordered and that we want to 
   insert a new element ai. If suf(ai) >= suf(a[j-1]) we are done.
   If suf(ai)<suf(a[j-1]) we notice that: if lcpi>m_pLCP[j-2] then
   suf(ai)>suf(a[j-2]) and we can stop since
   j-2 mmmmmmg
   j-1 mmmmmmmmmmmm
   ai  mmmmmmmmmmmf ] lcpi
   so we write a[j-1] in position j and ai in position j-1.
   if lcpi==m_pLCP[j-2] then we need to compare suf(ai) with suf(a[j-2])
   j-2 mmmmmmmmmmm?     we can have either ?<f or ?>f or ?==f  
   j-1 mmmmmmmmmmmm
   j   mmmmmmmmmmmf
   so we move a[j-1] to position j and compare suf(ai) with suf(a[j-2])
   starting from lcpi.
   Finally, if lcpi<m_pLCP[j-2] then
   j-2 mmmmmmmmmmmmmmmmmmg
   j-1 mmmmmmmmmmmmmmmmmmm
   j   mmmmmmmmmmmmmf
   hence we have suf(ai)<suf(a[j-2]) and we consider a[j-3];
   if lcpi<m_pLCP[j-3] we go on look at a[j-4] and go on.
   if m_pLCP[j]>m_pLCP[j-3] we are in the following position:
   j-3 mmmmmmc                   
   j-2 mmmmmmmmmmmmmmmg        
   j-1 mmmmmmmmmmmmmmmm          
   j   mmmmmmmmmmf
   and we know that suf(ai) is larger than suf(a[j-3]). If we find that 
   lcpi==m_pLCP[j-3] then we must compare suf(ai) with suf(a[j-3])
   but starting with position lcpi
   ***************************************************************** */

void 
CDSsort::shallow_inssort_lcp(INT32 *a, INT32 n, unsigned char *text_depth)
{   
  INT32 i, j, j1, lcp_new, r, lcpi;
  INT32 ai;
  INT32 cmp_from_limit;
  unsigned char *text_depth_ai;

  // --------- initialize ----------------
  m_LCPAuxArray[0] = -1;               // set m_pLCP[-1] = -1
  for(i=0;i<n;i++) m_pLCP[i]=0;     // I think this loop is not necessary
  // cmp_from_limit is # of cmp's to be done to reach m_ShallowLimit cmp's
  cmp_from_limit = (int)(m_pShallowTextLimit-text_depth);

  // ----- start insertion sort -----------
  for (i = 1; i< n ; i++) {
    ai = a[i]; lcpi = 0;
    text_depth_ai = ai + text_depth;
    j=i; j1=j-1;                  // j1 is a shorhand for j-1
    while(1) {           

      // ------ compare ai with a[j-1] --------
      m_CmpLeft = cmp_from_limit-lcpi;  
      r = cmp_unrolled_shallow_lcp(lcpi+a[j1]+text_depth,lcpi+text_depth_ai);
      lcp_new = cmp_from_limit - m_CmpLeft;       // m_pLCP between ai and a[j1] 
      assert(r!=0 || lcp_new>= cmp_from_limit);

      if(r<=0) {         // we have a[j-1] <= ai
	m_pLCP[j1]=lcp_new; // ai will be written in a[j]; update m_pLCP[j-1]
	break;
      }

      // --- we have a[j-1]>ai. a[j-1] and maybe other will be moved down 
      // --- use m_pLCP to move down as many elements of a[] as possible
      lcpi = lcp_new;                
      do {
	a[j] = a[j1];               // move down a[j-1]
        m_pLCP[j] = m_pLCP[j1];           // move down m_pLCP[j-1]
	j=j1; j1--;                 // update j and j1=j-1
      } while(lcpi<m_pLCP[j1]);        // recall that m_pLCP[-1]=-1

      if(lcpi>m_pLCP[j1]) break;       // ai will be written in position j

      // if we get here lcpi==m_pLCP[j1]: we will compare them at next iteration

    }     // end for(j=i ...
    a[j]=ai;
    m_pLCP[j]=lcpi;
  }       // end for(i=1 ... 

  // ----- done with insertion sort. now sort groups of equal strings
  for(i=0;i<n-1;i=j+1) {
    for(j=i; j<n ;j++)
      if(m_pLCP[j]<cmp_from_limit) break;
    if(j-i>0) 
      helped_sort(a+i,j-i+1,m_ShallowLimit); 
  }
}











