/*****************************************************************
 * Portions extracted from:
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CShuffle::CShuffle(void)
{
m_pRandomMersenne = new CRandomMersenne(1);		// ensure initially will generate reproducible series of random numbers
}

CShuffle::~CShuffle(void)
{
if(m_pRandomMersenne != NULL)
	delete m_pRandomMersenne;
}


int 
CShuffle::CHOOSE(int Limit) 
{
return(m_pRandomMersenne->IRandom(0,Limit));
};


int
CShuffle::FChoose(int *p, int N)
{
int roll;                   /* random fraction */
int sum;					/* integrated prob */
int   Idx;                    /* counter over the probs */
int Accum;
Accum = 0;
for(Idx=0; Idx < N; Idx++)
	Accum += p[Idx];
roll  = CHOOSE(Accum-1);
sum  = 0;
for (Idx = 0; Idx < N; Idx++)
   {
    sum += p[Idx];
    if (roll < sum) 
		return Idx;
    }
return (CHOOSE(N-1));           /* bulletproof */
}

/* shuffle.c
 * 
 * Routines for randomizing sequences.
 *  
 * All routines are sequence dependent (a,c,g, anf t);
 
 * SeqShuffle()   - shuffled sequence, preserve mono-symbol composition.
 * SeqDPShuffle() - shuffled sequence preserve mono- and di-symbol composition.
 * 
 * SeqMarkov0()   - random sequence, same zeroth order Markov properties.
 * SeqMarkov1()   - random sequence, same first order Markov properties.
 * 
 * SeqReverse()   - simple reversal of string
 * SeqRegionalShuffle() -  mono-symbol shuffled string in regional windows
 *
 * There are also similar routines for shuffling alignments:
 *
 * AlignmentShuffle()   - alignment version of StrShuffle().
 * AlignmentBootstrap() - sample with replacement; a bootstrap dataset.
 * QRNAShuffle()        - shuffle a pairwise alignment, preserving all gap positions.
 * 
 * CVS $Id: shuffle.c,v 1.8 2003/04/14 16:00:16 eddy Exp $
 */


/* Function: StrShuffle()
 * 
 * Purpose:  Returns a shuffled version of pSeq2, in pSeq1.
 *           (pSeq1 and pSeq2 can be identical, to shuffle in place.)
 *  
 * Args:     pSeq1 - allocated space for shuffled sequence.
 *           pSeq2 - sequence to shuffle.
 *           
 * Return:   1 on success.
 */
int
CShuffle::SeqShuffle(int SeqLen, char *pSeq1, char *pSeq2)
{
int  len;
int  pos;
char c;

if (pSeq1 != pSeq2) strcpy(pSeq1, pSeq2);
for (len = SeqLen; len > 1; len--)
{				
    pos       = CHOOSE(len);
    c         = pSeq1[pos];
    pSeq1[pos]   = pSeq1[len-1];
    pSeq1[len-1] = c;
}
return 1;
}

/* Function: SeqDPShuffle()
 * Date:     SRE, Fri Oct 29 09:15:17 1999 [St. Louis]
 *
 * Purpose:  Returns a shuffled version of pSeq2, in pSeq1.
 *           (pSeq1 and pSeq2 may be identical; i.e. a string
 *           may be shuffled in place.) The shuffle is a  
 *           "doublet-preserving" (DP) shuffle. Both
 *           mono- and di-symbol composition are preserved.
 *           
 *           Done by searching for a random Eulerian 
 *           walk on a directed multigraph. 
 *           Reference: S.F. Altschul and B.W. Erickson, Mol. Biol.
 *           Evol. 2:526-538, 1985. Quoted bits in my comments
 *           are from Altschul's outline of the algorithm.
 *
 * Args:     pSeq1   - RETURN: the string after it's been shuffled
 *                    (space for pSeq1 allocated by caller)
 *           pSeq2   - the string to be shuffled
 *
 * Returns:  0 if string can't be shuffled (it's not all [a-zA-z]
 *             alphabetic.
 *           1 on success. 
 */
int
CShuffle::SeqDPShuffle(int SeqLen,unsigned char *pSeq1, unsigned char *pSeq2)
{
int    len;
int    pos;		/* a position in pSeq1 or pSeq2 */
int    x,y;		/* indices of two bases */
unsigned char *EdgeList[4];     /* edge lists: EdgeList[0] is the edge list from vertex A */
int   EdgeListLen[4];    /* lengths of edge lists */
int   iE[4];		/* positions in edge lists */
int    n;		/* tmp: remaining length of an edge list to be shuffled */
unsigned char   sf;		/* last base in pSeq2 */
unsigned char   Z[4];	/* connectivity in last edge graph Z */ 
int    keep_connecting; /* flag used in Z connectivity algorithm */
int    is_eulerian;		/* flag used for when we've got a good Z */

len = SeqLen;
 
  /* "(1) Construct the doublet graph G and edge ordering E
   *      corresponding to S."
   * 
   * Note that these also imply the graph G; and note,
   * for any list x with EdgeListLen[x] = 0, vertex x is not part
   * of G.
   */
   for (x = 0; x < 4; x++)
    {
      EdgeList[x]  = new unsigned char [len-1];
      EdgeListLen[x] = 0; 
    }

  x = pSeq2[0] & 0x03;
  for (pos = 1; pos < len; pos++)
    {
      y = pSeq2[pos]  & 0x03;
      EdgeList[x][EdgeListLen[x]] = y;
      EdgeListLen[x]++;
      x = y;
    }
  
  /* Now we have to find a random Eulerian edge ordering.
   */
  sf = pSeq2[len-1]  & 0x03; 
  is_eulerian = 0;
  while (! is_eulerian)
    {
      /* "(2) For each vertex s in G except s_f, randomly select
       *      one edge from the s edge list of E(S) to be the
       *      last edge of the s list in a new edge ordering."
       *
       * select random edges and move them to the end of each 
       * edge list.
       */
      for (x = 0; x < 4; x++)
	{
	  if (EdgeListLen[x] == 0 || x == sf) continue;
	  
	  pos           = CHOOSE(EdgeListLen[x]);
	  y             = EdgeList[x][pos];		
	  EdgeList[x][pos]     = EdgeList[x][EdgeListLen[x]-1];
	  EdgeList[x][EdgeListLen[x]-1] = y;
	}

      /* "(3) From this last set of edges, construct the last-edge
       *      graph Z and determine whether or not all of its
       *      vertices are connected to s_f."
       * 
       * a probably stupid algorithm for looking at the
       * connectivity in Z: iteratively sweep through the
       * edges in Z, and build up an array (confusing called Z[x])
       * whose elements are 1 if x is connected to sf, else 0.
       */
      for (x = 0; x < 4; x++) Z[x] = 0;
      Z[(int) sf] = keep_connecting = 1;

      while (keep_connecting) {
	keep_connecting = 0;
	for (x = 0; x < 4; x++)
	  {
	    y = EdgeList[x][EdgeListLen[x]-1];            /* xy is an edge in Z */
	    if (Z[x] == 0 && Z[y] == 1)   /* x is connected to sf in Z */
	      {
		Z[x] = 1;
		keep_connecting = 1;
	      }
	  }
      }

      /* if any vertex in Z is tagged with a 0, it's
       * not connected to sf, and we won't have a Eulerian
       * walk.
       */
      is_eulerian = 1;
      for (x = 0; x < 4; x++)
	{
	  if (EdgeListLen[x] == 0 || x == sf) continue;
	  if (Z[x] == 0) {
	    is_eulerian = 0;
	    break;
	  }
	}

      /* "(4) If any vertex is not connected in Z to s_f, the
       *      new edge ordering will not be Eulerian, so return to
       *      (2). If all vertices are connected in Z to s_f, 
       *      the new edge ordering will be Eulerian, so
       *      continue to (5)."
       *      
       * e.g. note infinite loop while is_eulerian is FALSE.
       */
    }

  /* "(5) For each vertex s in G, randomly permute the remaining
   *      edges of the s edge list of E(S) to generate the s
   *      edge list of the new edge ordering E(S')."
   *      
   * Essentially a StrShuffle() on the remaining EdgeListLen[x]-1 elements
   * of each edge list; unfortunately our edge lists are arrays,
   * not strings, so we can't just call out to StrShuffle().
   */
  for (x = 0; x < 4; x++)
    for (n = EdgeListLen[x] - 1; n > 1; n--)
      {
	pos       = CHOOSE(n);
	y         = EdgeList[x][pos];
	EdgeList[x][pos] = EdgeList[x][n-1];
	EdgeList[x][n-1] = y;
      }

  /* "(6) Construct sequence S', a random DP permutation of
   *      S, from E(S') as follows. Start at the s_1 edge list.
   *      At each s_i edge list, add s_i to S', delete the
   *      first edge s_i,s_j of the edge list, and move to
   *      the s_j edge list. Continue this process until
   *      all edge lists are exhausted."
   */ 
  for (x = 0; x < 4; x++) iE[x] = 0; 

  pos = 0; 
  x = pSeq2[0]  & 0x03;
  while (1) 
    {
      pSeq1[pos++] = x;	/* add s_i to S' */
      
      y = EdgeList[x][iE[x]];
      iE[x]++;			/* "delete" s_i,s_j from edge list */
  
      x = y;			/* move to s_j edge list. */

      if (iE[x] == EdgeListLen[x])
	break;			/* the edge list is exhausted. */
    }
  pSeq1[pos++] = sf;

  /* Free and return */
  for(x=0;x<4;x++)
	  delete EdgeList[x]; 

  return 1;
}

  
/* Function: SeqMarkov0()
 * Date:     SRE, Fri Oct 29 11:08:31 1999 [St. Louis]
 *
 * Purpose:  Returns a random string pSeq1 with the same
 *           length and zero-th order Markov properties
 *           as pSeq2. 
 *           
 *           pSeq1 and pSeq2 may be identical, to randomize pSeq2
 *           in place.
 *
 * Args:     pSeq1 - allocated space for random string
 *           pSeq2 - string to base pSeq1's properties on.
 *
 * Returns:  1 on success; 0 if pSeq2 doesn't look alphabetical.
 */
int 
CShuffle::SeqMarkov0(int SeqLen,char *pSeq1, char *pSeq2)
{
  int   len;
  int   pos; 
  int p[4];			/* symbol probabilities */

  len = SeqLen;

  /* Collect zeroth order counts and convert to frequencies.
   */
  p[0]=p[1]=p[2]=p[3]=0;
  for (pos = 0; pos < len; pos++)
	p[pSeq2[pos]] += 1;
  
  /* Generate a random string using those p's.
   */
  for (pos = 0; pos < len; pos++)
    pSeq1[pos] = FChoose(p, 4);
  return 1;
}


/* Function: SeqMarkov1()
 * Date:     SRE, Fri Oct 29 11:22:20 1999 [St. Louis]
 *
 * Purpose:  Returns a random string pSeq1 with the same
 *           length and first order Markov properties
 *           as pSeq2. 
 *           
 *           pSeq1 and pSeq2 may be identical, to randomize pSeq2
 *           in place.
 *
 * Args:     pSeq1 - allocated space for random string
 *           pSeq2 - string to base pSeq1's properties on.
 *
 * Returns:  1 on success; 0 if pSeq2 doesn't look alphabetical.
 */
int 
CShuffle::SeqMarkov1(int SeqLen, char *pSeq1, char *pSeq2)
{
  int   len;
  int   pos; 
  int   x,y;
  int   i;		/* initial symbol */
  int p[4][4];	/* symbol probabilities */

  len = SeqLen;

  /* Collect first order counts and convert to frequencies.
   */
 memset(p,0,sizeof(p));

  i = x = pSeq2[0];
  for (pos = 1; pos < len; pos++)
    {
      y = pSeq2[pos];
      p[x][y] += 1; 
      x = y;
    }

  /* Generate a random string using those p's.
   */
  x = i;
  pSeq1[0] = x;
  for (pos = 1; pos < len; pos++)
    {
      y = FChoose(p[x], 4);
      pSeq1[pos] = y;
      x = y;
    } 
  return 1;
}



/* Function: SeqReverse()
 * Date:     SRE, Thu Nov 20 10:54:52 1997 [St. Louis]
 * 
 * Purpose:  Returns a reversed version of pSeq2, in pSeq1.
 *           (pSeq1 and pSeq2 can be identical, to reverse in place)
 * 
 * Args:     pSeq1 - allocated space for reversed string.
 *           pSeq2 - string to reverse.
 *           
 * Return:   1.
 */                
int
CShuffle::SeqReverse(int SeqLen, char *pSeq1, char *pSeq2)
{
int  len;
int  pos;
char c;
len = SeqLen;
for (pos = 0; pos < len/2; pos++)
	{				/* swap ends */
    c             = pSeq2[len-pos-1];
    pSeq1[len-pos-1] = pSeq2[pos];
    pSeq1[pos]       = c;
	}
if (len%2) /* copy middle residue in odd-len pSeq2 */
	pSeq1[pos] = pSeq2[pos]; 
	 
pSeq1[len] = '\0';
return 1;
}

/* Function: SeqRegionalShuffle()
 * Date:     SRE, Thu Nov 20 11:02:34 1997 [St. Louis]
 * 
 * Purpose:  Returns a regionally shuffled version of pSeq2, in pSeq1.
 *           (pSeq1 and pSeq2 can be identical to regionally 
 *           shuffle in place.) See [Pearson88].
 *           
 * Args:     pSeq1 - allocated space for regionally shuffled string.
 *           pSeq2 - string to regionally shuffle
 *           w  - window size (typically 10 or 20)      
 *           
 * Return:   1.
 */
int
CShuffle::SeqRegionalShuffle(int SeqLen, char *pSeq1, char *pSeq2, int w)
{
int  len;
char c;
int  pos;
int  i, j;

if (pSeq1 != pSeq2) 
	memmove(pSeq1, pSeq2,SeqLen);
len = SeqLen;

#ifdef _WIN32
for (i = 0; i < len; i += w)
for (j = __min(len-1, i+w-1); j > i; j--)
    {
        pos     = i + CHOOSE(j-i);
        c       = pSeq1[pos];
        pSeq1[pos] = pSeq1[j];
        pSeq1[j]   = c;
        }       
#else
for (i = 0; i < len; i += w)
        j = len-1;
        if(j > i+w-1)
                j = i+w-1;
        for (; j > i; j--)
        {
                pos     = i + CHOOSE(j-i);
                c       = pSeq1[pos];
                pSeq1[pos] = pSeq1[j];
                pSeq1[j]   = c;
                }

#endif

return 1;
}


