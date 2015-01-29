//----- Full credit is given to the following authors, this module is very dependent on their code
//FM-index Version 2
//Authors: Paolo Ferragina and Rossano Venturini\n");
//Dipartimento di Informatica, University of Pisa, Italy\n");
//

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "../hdrs/commhdrs.h"
#endif


#include "./fmindexpriv.h"

#define BZ_RUNA 0
#define BZ_RUNB 1
#define BZ_G_SIZE   50
#define BZ_N_ITERS  4
#define BZ_LESSER_ICOST  0
#define BZ_GREATER_ICOST 15


CFMIndex::CFMIndex(void)
{
pmtf_start = NULL;
gmtflen = 0;
m_pIndex = new fm_index;
if(m_pIndex != NULL)
	memset(m_pIndex,0,sizeof(fm_index));
}

CFMIndex::~CFMIndex(void)
{
if(m_pIndex != NULL)
	{
	free_index();
	delete m_pIndex;
	}
}


char *
CFMIndex::GetErrText(int error)
{
return((char *)error_index(error));
}


int
CFMIndex::CreateIndex(char *pszInFile,				// create from contents of this file
					char *pszOutFile,				// write index into this file
					unsigned int *pText_len,		// returned file content length before index created	
					unsigned int *pIndex_len,		// created index length
					int bsl1,						// bucket size level 1 (in 1K increments) must be a multiple of bsl2
					int bsl2,						// bucket size level 2 (in 1K  increments)
					double Freq)					// marker frequency 0.0-1.0
{
unsigned char *text;
UINT32 text_len, index_len;
int error;

if(pText_len != NULL)
	*pText_len = 0;
if(pIndex_len != NULL)
	*pIndex_len = 0;


error = fm_read_file(pszInFile, &text, &text_len);
if (error < 0)
	return(error);

error = fm_build_config (Freq, bsl1, bsl2, 1);
if (error < 0)
	return(error);

error =  fm_build(text, text_len);
if (error < 0)
	return(error);

error = save_index(pszOutFile);
if (error < 0)
	return(error);
		
index_size(&index_len);

	/* How to use int save_index_mem(unsigned char *compress) */
	#if 0
	unsigned char *foo = malloc(index_len * sizeof(unsigned char));
	if(!foo) {
		fprintf(stderr,"Alloc Failed\n");
		exit(1);
		}
	error = save_index_mem( foo);
	IFERROR(error);
	#endif
error = free_index();
if(error < 0)
	return(error);
if(pText_len != NULL)
	*pText_len = text_len;
if(pIndex_len != NULL)
	*pIndex_len = index_len;
return(0);
}


/*
 * Main Search Functions Rossano Venturini 
 */



/*
 * Writes in numocc the number of occurrences of pattern[0..length-1] in
 * index. It also allocates occ (which must be freed by the caller) and
 * writes the locations of the numocc occurrences in occ, in arbitrary
 * order. 
 */
int
CFMIndex::locate (unsigned char * pattern, UINT32 length, UINT32 ** occ,UINT32 * numocc)
{

	multi_count *groups;
	int i, num_groups = 0, state;
	UINT32 *occs = NULL;
	*numocc = 0;
	*occ = NULL;

	if(m_pIndex->smalltext)  //uses Boyer-Moore algorithm
		return fm_boyermoore(pattern, length, occ, numocc);
	
	/* count */
	num_groups = fm_multi_count (pattern, length, &groups);

	if (num_groups <= 0)
		return num_groups;

	for (i = 0; i < num_groups; i++)
		*numocc += groups[i].elements;

	occs = *occ = (UINT32 *)new  UINT32 [sizeof (UINT32) * (*numocc)];
	if (*occ == NULL)
		{
		*numocc = 0;
		return FM_OUTMEM;
		}

	for (i = 0; i < num_groups; i++)
	{
		state = multi_locate (groups[i].first_row, groups[i].elements, occs);
		if (state < 0)
		{
			delete *occ;
			*occ = NULL;
			*numocc = 0;
			return state;
		}
		occs += groups[i].elements;
	}

	free (groups);
	return FM_OK;

}


/*
 * Writes in numocc the number of occurrences of pattern[0..length-1] in
 * index. 
 */
int
CFMIndex::count (unsigned char * pattern, UINT32 length, UINT32 * numocc)
{
	multi_count *groups;
	int i, num_groups = 0;

	*numocc = 0;

	if(m_pIndex->smalltext) { //uses Boyer-Moore algorithm
		UINT32 *occ;
		int error = fm_boyermoore(pattern, length, &occ, numocc);
		if(error < 0)
			*numocc = 0;
		if(*numocc > 0) 
			free(occ);
		return error;		
	}
	
	num_groups = fm_multi_count (pattern, length, &groups);

	if (num_groups <= 0)
		return num_groups;

	for (i = 0; i < num_groups; i++)
		*numocc += groups[i].elements;

	free (groups);
	return FM_OK;

}



#define ADD_LIST(_first_row, _elements) {\
	if ((used+1) == allocated) {\
		allocated += 5000;\
		lista = (multi_count *)realloc(lista, sizeof(multi_count)*allocated);\
		if (lista == NULL) return FM_OUTMEM;\
		}\
	lista[used].first_row = _first_row;\
	lista[used++].elements = _elements;\
}


/*
 * Count ricorsiva per la multilocate. se nel pattern c'e' una occorenza
 * di subchar allora devo dividere la ricerca con due rami distinti. 
 */
int
CFMIndex::count_row_mu (unsigned char * pattern, UINT32 len, UINT32 sp,UINT32 ep)
{
	unsigned char chars_in[ALPHASIZE];
	UINT32 occsp[ALPHASIZE], occep[ALPHASIZE];
	int num_char, i, find = 0;
	unsigned char c;
	UINT32 ssp, sep;
	/*
	 * Versione semplice - possibile fare meglio 
	 */
	while ((sp <= ep) && (len > 0))
	{
		c = pattern[--len];
		find = 0;
		if(sp == 0) 
			num_char = occ_all (0, EOF_shift (ep), occsp, occep, chars_in);
		else 
			num_char = occ_all (EOF_shift (sp - 1), EOF_shift (ep), occsp, occep, chars_in);
        ep = 0;	sp = 1;
		for(i=0; i<num_char; i++) {
			if ((m_pIndex->skip >1) && (chars_in[i] == m_pIndex->specialchar))
			{
				ssp = m_pIndex->bwt_occ[m_pIndex->specialchar] + occsp[m_pIndex->specialchar]; 
				sep = ssp +(occep[m_pIndex->specialchar] - occsp[m_pIndex->specialchar]) - 1;
				assert((ssp<m_pIndex->text_size)&&(sep<m_pIndex->text_size));
				count_row_mu (pattern, len+1, ssp, sep);		
				if (find==1) break;
				find = 1;
			}
			if (chars_in[i] == c) {
				sp = m_pIndex->bwt_occ[c] + occsp[c]; 
				ep = sp + (occep[c] - occsp[c]) - 1;
				assert((sp<m_pIndex->text_size)&&(ep<m_pIndex->text_size));
				if (find==1) break;
				find = 1;
			}
		}
	}
	
		
	/*
	 * return number of occurrences 
	 */
	if (ep < sp)
		return FM_OK;	/* nessuna occorrenza */
	ADD_LIST (sp, ep - sp + 1);	/* Aggiunge gruppo trovato alla lista */

	return FM_OK;
}


/*
 * Count per la ricerca con marcatura per multi locate il problema con
 * questo tipo di marcatura e' che si devono considerare due rami di
 * ricerca ogni volta che si incontra il carattere sostituito nel
 * pattern. Questa divisione peggiora le prestazioni della count. La
 * ricerca sui diversi cammini e' implementata con una chiamata ricorsiva
 * non molto efficiente. 
 */
int
CFMIndex::fm_multi_count (unsigned char * pattern, UINT32 len,multi_count ** list)
{

	UINT32 sp, ep, i, j;
	unsigned char c;

	*list = NULL;
	allocated = len;
	used = 0;

	lista = (multi_count *)malloc (allocated * sizeof (multi_count));
	if (lista == NULL)
		return FM_OUTMEM;

	/*
	 * remap pattern 
	 */
	assert (len > 0);

	for (i = 0; i < len; i++)
	{
		if (m_pIndex->bool_char_map[pattern[i]] == 0)
			{
				for (j = 0; j <= i; j++) 
					pattern[j] = m_pIndex->inv_char_map[pattern[j]];
				free(lista);
				return 0;	/* char not in file */
			}
			
		pattern[i] = m_pIndex->char_map[pattern[i]];	/* remap char */
		if((m_pIndex->skip >1)&&(pattern[i]==m_pIndex->specialchar))
			{
				for (j = 0; j <= i; j++) 
					pattern[j] = m_pIndex->inv_char_map[pattern[j]];
				free(lista);
				return 0;	/* char not in file */
			}
			
	}

	/* get initial sp and ep values */
	c = pattern[len - 1];
	sp = m_pIndex->bwt_occ[c];
	if (c == m_pIndex->alpha_size - 1)
		ep = m_pIndex->text_size - 1;
	else
		ep = m_pIndex->bwt_occ[c + 1]- 1;
	
	count_row_mu (pattern, len - 1, 	sp, ep);	// ricerca per il carattere c

	#if 0
	if (used = 0)
		lista = realloc (lista, sizeof (multi_count) * used);
	else
		free (lista);
	#endif
	
	if(used == 0) {
		free(lista); lista = NULL;
	}
	
	/* inverse remap pattern  */
	for (i = 0; i < len; i++) 
		pattern[i] = m_pIndex->inv_char_map[pattern[i]];
	
	*list = lista;

	return used;
}


/*
 * Restituisce la posizione nel testo memorizzata sul compresso in i-esima 
 * posizione. 
 */
void
CFMIndex::get_pos (UINT32 first_row, UINT32 element, UINT16 step, UINT32 * pos)
{
	int offset, skipbits;
	UINT32 postext, i;
	offset = first_row * m_pIndex->log2textsize;
	fm_init_bit_reader (m_pIndex->start_prologue_occ + (offset >> 3));
	skipbits = offset % 8;	// bits are to be skipped 
	if (skipbits)
	{
		fm_bit_read24 (skipbits, skipbits);
	}

	for (i = 0; i < element; i++, pos++)
	{			// cerca tutto il gruppo
		postext = fm_bit_read (m_pIndex->log2textsize);	// read text pos 
		*pos = postext + step;
	}

}


/*
 * Multilocate 
 */
int
CFMIndex::multi_locate (UINT32 sp, UINT32 element, UINT32 * positions)
{

	if(element == 0) return FM_OK;
	if(m_pIndex->skip == 0) return FM_NOMARKEDCHAR;
	if(m_pIndex->skip == 1) {
		get_pos (sp, element, 0, positions);
		return element;
	}

	UINT32 curr_row, used, recurs, elements, diff;
	UINT32 *elem_array;	// contiene il numero di elementi del sottogruppo
	UINT32 occsp[ALPHASIZE];
	UINT32 occep[ALPHASIZE];
	UINT16 *step_array, step;	// come sopra ma passi
	unsigned char chars_in[ALPHASIZE];
	int j, state, num_char;
	/*
	 * per singola locate 
	 */
	UINT32 occ_sb[ALPHASIZE], occ_b[ALPHASIZE];
	UINT16 localstep;
	unsigned char c, c_sb, cb;

	/* Caso in cui l'ultimo carattere e' m_pIndex->specialchar 
	   ovvero le riga sp e sp+element sono compresse tra  
	   m_pIndex->bwt_occ[m_pIndex->specialchar]  e  
       m_pIndex->bwt_occ[m_pIndex->specialchar+1] */
	if ((sp >= m_pIndex->occcharinf) && (sp+element-1 < m_pIndex->occcharsup)) {
		get_pos (sp - m_pIndex->occcharinf, element, 0, positions);
		return element;
	}
	
	step_array = (UINT16 *)malloc (sizeof (UINT32) * element);
	elem_array = (UINT32 *)malloc (sizeof (UINT32) * element);
	if ((step_array ==NULL) || (elem_array == NULL))
		return FM_OUTMEM;
	
	used = 0;
	recurs = element - 1;
	elem_array[recurs] = element;
	positions[recurs] = sp;
	step_array[recurs] = 0;

	while (recurs < element)
	{
		elements = elem_array[recurs];
		curr_row = positions[recurs];
		step = step_array[recurs++] + 1;

		if ((m_pIndex->bwt_eof_pos >= curr_row) &&
		    (m_pIndex->bwt_eof_pos < curr_row + elements))
			positions[used++] = step - 1;

		num_char =
			occ_all (EOF_shift (curr_row - 1),
				 EOF_shift (curr_row + elements - 1), occsp,occep, chars_in);

		for (j = 0; j < num_char; j++)
		{
			cb = chars_in[j];
			diff = occep[cb] - occsp[cb];
			if (cb == m_pIndex->specialchar)
			{	
				get_pos (occsp[cb], diff, step, positions + used);	
				/* pos of * + step  +1 */
				used += diff;
				continue;
			}

			if (diff == 1)
			{
				curr_row = occsp[cb] + m_pIndex->bwt_occ[cb];
				assert(curr_row<m_pIndex->text_size);
				c_sb = cb;
				localstep = step;
				while (c_sb != m_pIndex->specialchar)
				{
					if (curr_row == m_pIndex->bwt_eof_pos)
					{
						positions[used++] = localstep;
						break;
					}
					state = get_info_sb (EOF_shift
							     (curr_row),
							     occ_sb);
					if (state < 0)
						return FM_SEARCHERR;
					state = get_info_b (NULL_CHAR,
							    EOF_shift
							    (curr_row), occ_b,
							    WHAT_CHAR_IS);
					c = state;
					if (state < 0)
						return FM_SEARCHERR;
					c_sb = m_pIndex->inv_map_sb[c];
					curr_row =
						m_pIndex->bwt_occ[c_sb] +
						occ_sb[c_sb] + occ_b[c] - 1;
					localstep++;
				}
				if (curr_row == m_pIndex->bwt_eof_pos)
					continue;
				get_pos (curr_row - m_pIndex->occcharinf, 1, localstep,
					 positions + used);
				used++;
				continue;
			}

			elem_array[--recurs] = diff;	/* simula ricorsione */
			positions[recurs] = occsp[cb] + m_pIndex->bwt_occ[cb];
			assert(occsp[cb] + m_pIndex->bwt_occ[cb]<m_pIndex->text_size);
			step_array[recurs] = step;
		}
	}

	if (used != element) 
		return FM_GENERR;
	free (elem_array);
	free (step_array);
	
	return FM_OK;
}


/* Boyer-Moore algorithm to support small files */

	
void CFMIndex::preBmBc(unsigned char *x, int m, int bmBc[]) {
   int i;
 
   for (i = 0; i < ALPHASIZE; ++i)
      bmBc[i] = m;
   for (i = 0; i < m - 1; ++i)
      bmBc[x[i]] = m - i - 1;
}
 
 
void CFMIndex::suffixes(unsigned char *x, int m, int *suff) 
{
   int f, g, i;
   f = 0;	
   suff[m - 1] = m;
   g = m - 1;
   for (i = m - 2; i >= 0; --i) {
      if (i > g && suff[i + m - 1 - f] < i - g)
         suff[i] = suff[i + m - 1 - f];
      else {
         if (i < g)
            g = i;
         f = i;
         while (g >= 0 && x[g] == x[g + m - 1 - f])
            --g;
         suff[i] = f - g;
      }
   }
}
 
void CFMIndex::preBmGs(unsigned char *x, int m, int bmGs[]) 
{
   int i, j, *suff;
   suff = (int *)malloc(m * sizeof(int));
   suffixes(x, m, suff);
 
   for (i = 0; i < m; ++i)
      bmGs[i] = m;
   j = 0;
   for (i = m - 1; i >= -1; --i)
      if (i == -1 || suff[i] == i + 1)
         for (; j < m - 1 - i; ++j)
            if (bmGs[j] == m)
               bmGs[j] = m - 1 - i;
   for (i = 0; i <= m - 2; ++i)
      bmGs[m - 1 - suff[i]] = m - 1 - i;
   free(suff);
}

int CFMIndex::fm_boyermoore(unsigned char * pattern, UINT32 length, UINT32 ** occ, UINT32 * numocc) 
{
   
   UINT32 j;
   int i, *bmGs, bmBc[ALPHASIZE];
   UINT32 alloc = 10;
   *numocc = 0;
	if(m_pIndex->text_size < length)
		{
		*occ = NULL;
		return(FM_OK);
		}

   *occ = (UINT32 *)malloc(sizeof(UINT32)*alloc);
   if(*occ == NULL) 
	return FM_OUTMEM;   

   bmGs = (int *)malloc(sizeof(int)*length);
   if(bmGs == NULL)
   {
    free(*occ);
	*occ = NULL;
	return FM_OUTMEM;   
   }
	
   /* Preprocessing */
   preBmGs(pattern, (int) length, bmGs);
   preBmBc(pattern, (int) length, bmBc);
 
   /* Searching */
   j = 0;
   while (j <=  m_pIndex->text_size - length) {
      for (i = length - 1; i >= 0 && pattern[i] == m_pIndex->text[i + j]; --i);
      if (i < 0) {
		 (*numocc)++;
		 if(*numocc == alloc) {
				alloc = MIN(alloc*2,m_pIndex->text_size);
			 	*occ = (UINT32 *)realloc(*occ, sizeof(UINT32)*alloc);
   				if(*occ == NULL)
					{
					free(bmGs);
					return FM_OUTMEM;
					}
			}
		 (*occ)[*numocc-1] = j;
         j += bmGs[0];
      }
      else
         j += MAX((UINT32) bmGs[i], bmBc[m_pIndex->text[i + j]] - length + 1 + i);
   }
   if(*numocc>0) 
	   *occ = (UINT32 *)realloc(*occ, sizeof(UINT32)*(*numocc));
   else
		{
		free(*occ);
		*occ = NULL;
		}
   free(bmGs);
   return FM_OK;
}





/*
 * Read/Write on index Functions Rossano Venturini 
 */

/*
 * Loads index from file called filename.fmi MODIFIED: index that is
 * allocated fmindex.compress fmindex.compress_size
 * 
 * with read_basic_prologue() fmindex.type_compression fmindex.text_size
 * fmindex.bwt_eof_pos fmindex.bucket_size_lev1 fmindex.bucket_size_lev2
 * fmindex.num_bucs_lev1 fmindex.alpha_size fmindex.specialchar
 * fmindex.skip fmindex.start_prologue_occ fmindex.start_prologue_info_sb
 * fmindex.start_prologue_info_b fmindex.subchar fmindex.bool_char_map[]
 * fmindex.char_map[] fmindex.inv_char_map[] fmindex.bits_x_occ
 * 
 * 
 * REQUIRES: filename: index file void *index point to NULL 
 */

int
CFMIndex::load_index (char * filename)
{
	int error;
		
	m_pIndex->compress_owner = 1;
	m_pIndex->owner = 0;
	
	/*
	 * Load index file 
	 */
	error =
		open_file (filename, &(m_pIndex->compress),
			   &(m_pIndex->compress_size));
	if (error)
		return error;

	error = fm_read_basic_prologue();
	if (error)
		return error;

	if(m_pIndex->smalltext) {
		m_pIndex->skip = 0;
		if(m_pIndex->text_size<SMALLSMALLFILESIZE) {
			m_pIndex->text = m_pIndex->compress+4;
			return FM_OK;
		}
		m_pIndex->owner = 1;
		m_pIndex->smalltext = 2;
		error = fm_bwt_uncompress();
		if (error < 0) return error;
		return FM_OK;
	}
	
	/*
	 * init some var 
	 */
	m_pIndex->int_dec_bits =
		int_log2 (int_log2
			  (m_pIndex->bucket_size_lev1 -
			   m_pIndex->bucket_size_lev2));
	
	
	if(m_pIndex->skip >1) {
	m_pIndex->occcharinf = m_pIndex->bwt_occ[m_pIndex->specialchar];
	if(m_pIndex->specialchar==m_pIndex->alpha_size-1)
		m_pIndex->occcharsup = m_pIndex->text_size-1;
	else 
		m_pIndex->occcharsup = m_pIndex->bwt_occ[m_pIndex->specialchar+1];
	
	m_pIndex->num_marked_rows = m_pIndex->occcharsup-m_pIndex->occcharinf;
	} else  m_pIndex->num_marked_rows = 0;
	
	m_pIndex->mtf_seq =
		(unsigned char *) malloc (m_pIndex->bucket_size_lev2 * sizeof (unsigned char));
	if (m_pIndex->mtf_seq == NULL)
		return FM_OUTMEM;
	
	m_pIndex->var_byte_rappr = ((m_pIndex->log2textsize + 7) / 8)*8;
	
	return FM_OK;
}

int
CFMIndex::load_index_mem(unsigned char *compress, UINT32 size)
{

	int error;

	m_pIndex->compress = compress;
	m_pIndex->compress_size = size;
	m_pIndex->compress_owner = 0;
	m_pIndex->text = NULL;
	m_pIndex->lf = NULL;
	m_pIndex->bwt = NULL;
	
	error = fm_read_basic_prologue();
	if (error)
		return error;

	if(m_pIndex->smalltext) {
		m_pIndex->skip = 0;
		if(m_pIndex->text_size<SMALLSMALLFILESIZE) {
			m_pIndex->text = m_pIndex->compress+4;
			return FM_OK;
		}
		m_pIndex->smalltext = 2;
		error = fm_bwt_uncompress();
		if (error < 0) return error;
		return FM_OK;
	}
	
	/*
	 * init some var 
	 */
	m_pIndex->int_dec_bits =
		int_log2 (int_log2
			  (m_pIndex->bucket_size_lev1 -
			   m_pIndex->bucket_size_lev2));
	
	
	if(m_pIndex->skip >1) {
	m_pIndex->occcharinf = m_pIndex->bwt_occ[m_pIndex->specialchar];
	if(m_pIndex->specialchar==m_pIndex->alpha_size-1)
		m_pIndex->occcharsup = m_pIndex->text_size-1;
	else 
		m_pIndex->occcharsup = m_pIndex->bwt_occ[m_pIndex->specialchar+1];
	
	m_pIndex->num_marked_rows = m_pIndex->occcharsup-m_pIndex->occcharinf;
	} else  m_pIndex->num_marked_rows = 0;

	m_pIndex->mtf_seq =
		(unsigned char *) malloc (m_pIndex->bucket_size_lev2 * sizeof (unsigned char));
	if (m_pIndex->mtf_seq == NULL)
		return FM_OUTMEM;
	
	m_pIndex->var_byte_rappr = ((m_pIndex->log2textsize + 7) / 8)*8;
	return FM_OK;
}

/*
 * Open and Read .fmi file (whitout mmap()) 
 */
int
CFMIndex::open_file (char * filename, unsigned char ** file, UINT32 * size)
{

	char *outfilename;
	FILE *outfile = NULL;

	outfilename =
		(char *) malloc ((strlen (filename) + strlen (EXT) + 1) *
				 sizeof (char));
	if (outfilename == NULL)
		return FM_OUTMEM;

	outfilename = strcpy (outfilename, filename);
	outfilename = strcat (outfilename, EXT);

	outfile = fopen (outfilename, "rb");	// b is for binary: required by
	// DOS
	if (outfile == NULL)
		return FM_READERR;

	/*
	 * store input file length 
	 */
	if (fseek (outfile, 0, SEEK_END) != 0)
		return FM_READERR;
	*size = ftell (outfile);

	if (*size < 1)
		return FM_READERR;
	rewind (outfile);

	/*
	 * alloc memory for text 
	 */
	*file = (unsigned char *)malloc ((*size) * sizeof (unsigned char));
	if ((*file) == NULL)
		return FM_OUTMEM;

	UINT32 t =
		(UINT32) fread (*file, sizeof (unsigned char), (size_t) * size,
			       outfile);
	if (t != *size)
		return FM_READERR;

	fclose (outfile);
	free (outfilename);
	return FM_OK;
}


/*
 * read basic prologue from compress 
 */
int
CFMIndex::fm_read_basic_prologue (void)
{

	int i;
	UINT32 size;

	fm_init_bit_reader (m_pIndex->compress);
	m_pIndex->text_size = fm_uint_read ();
	if(m_pIndex->text_size< SMALLFILESIZE){
			m_pIndex->smalltext=1; 
			return FM_OK;
		}
	m_pIndex->smalltext = 0;
	m_pIndex->type_compression = fm_bit_read (8);
	m_pIndex->log2textsize = int_log2 (m_pIndex->text_size - 1);
	m_pIndex->bwt_eof_pos = fm_uint_read ();
	if (m_pIndex->bwt_eof_pos > m_pIndex->text_size)
		return FM_COMPNOTCORR;

	m_pIndex->bucket_size_lev1 = fm_bit_read (16) << 10;
	m_pIndex->bucket_size_lev2 = fm_bit_read (16);

	if (m_pIndex->bucket_size_lev1 % m_pIndex->bucket_size_lev2)
		return FM_COMPNOTCORR;

	m_pIndex->num_bucs_lev1 =
		(m_pIndex->text_size + m_pIndex->bucket_size_lev1 - 1) / m_pIndex->bucket_size_lev1;
	m_pIndex->num_bucs_lev2 =
		(m_pIndex->text_size + m_pIndex->bucket_size_lev2 - 1) / m_pIndex->bucket_size_lev2;

	/* mtf & alphabet information */
	m_pIndex->alpha_size = fm_bit_read (8) + 1;


	/* read Mark mode & starting position of occ list */
	m_pIndex->specialchar = (unsigned char) fm_bit_read (8);
	m_pIndex->skip =  fm_bit_read (32);
	unsigned int start = fm_uint_read();

	m_pIndex->start_prologue_occ = m_pIndex->compress + start;
	m_pIndex->start_prologue_info_sb = fm_uint_read ();
	m_pIndex->subchar = (unsigned char) fm_bit_read (8);	/* remapped cmpress alphabet */

	/* some information for the user */
	#if 0
	fprintf (stdout, "Compression type %d\n", m_pIndex->type_compression);
	fprintf (stdout, "Text Size %lu\n", m_pIndex->text_size);
	fprintf (stdout, "Bwt EOF %lu\n",m_pIndex->bwt_eof_pos);
	fprintf (stdout, "alphasize %d\n",m_pIndex->alpha_size);
	fprintf(stdout, "start prologue %lu\n", m_pIndex->start_prologue_occ);
	fprintf (stdout, "Compression method: ");
	switch (m_pIndex->type_compression)
		{
		case MULTIH:
			fprintf (stdout, "Huffman with multiple tables.\n");
			break;
	
		default:
		return FM_COMPNOTSUP;
	}
	#endif

	/* alphabet info and inverse char maps */
	for (i = 0; i < ALPHASIZE; i++)
		m_pIndex->bool_char_map[i] = fm_bit_read (1);

	for (i = 0, size = 0; i < ALPHASIZE; i++)
		if (m_pIndex->bool_char_map[i])
		{
			m_pIndex->char_map[i] = (unsigned char)size;
			m_pIndex->inv_char_map[size++] = (unsigned char) i;
		}
	assert (size == m_pIndex->alpha_size);

	/* prefix summed char-occ info momorizzate con m_pIndex->log2textsize bits */

	for (i = 0; i < m_pIndex->alpha_size; i++)
	{			// legge somme occorrenze
		// caratteri
		m_pIndex->bwt_occ[i] = fm_bit_read (m_pIndex->log2textsize);
	}

	/*
	 * calcola le occorrenze di ogni carattere nel testo 
	 */
	for (i = 1; i < m_pIndex->alpha_size; i++)
		m_pIndex->char_occ[i - 1] = (m_pIndex->bwt_occ[i]) - (m_pIndex->bwt_occ[i - 1]);

	m_pIndex->char_occ[(m_pIndex->alpha_size) - 1] =
		(m_pIndex->text_size) - (m_pIndex->bwt_occ[(m_pIndex->alpha_size) - 1]);

	/*
	 * Calcolo posizione inizio info buckets 
	 */
	m_pIndex->sb_bitmap_size = (m_pIndex->alpha_size+7)/8;
	
	m_pIndex->start_prologue_info_b = m_pIndex->start_prologue_info_sb + 
		(m_pIndex->sb_bitmap_size*m_pIndex->num_bucs_lev1)
		+ (m_pIndex->alpha_size * sizeof(UINT32) * (m_pIndex->num_bucs_lev1 - 1));

	return FM_OK;
}


/*
 * Frees the memory occupied by index. 
 */
int
CFMIndex::free_index (void)
{
	if (m_pIndex == NULL) 
		return FM_OK;
		
	if (m_pIndex->compress_owner == 2 ) 
		{ // Arrivo dalla build
		free(m_pIndex->compress);
		m_pIndex->compress = NULL;
		if (m_pIndex->smalltext==0)
			{
			free(m_pIndex->mtf_seq); // libero mtf del bucket
			m_pIndex->mtf_seq = NULL;
			}
		m_pIndex->text = m_pIndex->oldtext;
		}			
	
	if (m_pIndex->compress_owner == 1 ) { // Arrivo lettura da file
		if (!m_pIndex->smalltext)
			{
			free(m_pIndex->mtf_seq);
			m_pIndex->mtf_seq = NULL;
			}
		free(m_pIndex->compress);
		m_pIndex->compress = NULL;
		if(m_pIndex->smalltext==2)
			{
			free(m_pIndex->text);
			m_pIndex->text = NULL;
			}
		}
	
	if (m_pIndex->compress_owner == 0) { // Arrivo da lettura di memoria 
		if (!m_pIndex->smalltext)
			{
			free(m_pIndex->mtf_seq);
			m_pIndex->mtf_seq = NULL;
			}
		if (m_pIndex->smalltext == 2)
			{
			free(m_pIndex->text);
			m_pIndex->text = NULL;
			}
		}

	memset(m_pIndex,0,sizeof(fm_index));
	return FM_OK;
} 

/*
 * Obtains the length of the text represented by index. 
 */
int
CFMIndex::get_length (UINT32 * length)
{
*length = m_pIndex->text_size - m_pIndex->num_marked_rows;
return FM_OK;
}

/*
	Writes in size the size, in bytes, of index. This should be the size
	needed by the index to perform any of the operations it implements.
*/
int 
CFMIndex::index_size(UINT32 *size) {

	*size = m_pIndex->compress_size;
	return FM_OK;

}




/*
 * Funzioni per il calcolo delle occorrenze dei caratteri nei bucket e nei 
 * superbucket, e per la decompressione dei bucket.
 * 
 * 
 */

/*
 * Occ all Attenzione richiede che m_pIndex->mtf_seq sia riempita anche se nel
 * bucket e' presente un solo carattere 
 */
int 
CFMIndex::occ_all (UINT32 sp, UINT32 ep, UINT32 * occsp, UINT32 * occep,unsigned char * char_in)
{

	int i, state, diff, mod, b2end, remap;
	UINT32 occ_sb[ALPHASIZE], occ_b[ALPHASIZE];
	UINT32 occ_sb2[ALPHASIZE], occ_b2[ALPHASIZE];
	UINT32 num_buc_ep = ep / m_pIndex->bucket_size_lev2;
	int char_present = 0;	// numero caratteri distinti presenti nel bucket 
	unsigned char *c, d;

	/*
	 * conta il numero di occorrenze fino al subperbucket che contiene sp 
	 */
	state = get_info_sb (sp, occ_sb);	// get occ of all chars in prev
	// superbuckets 
	if (state < 0)
		return state;

	/*
	 * Controlla se sp ed ep stanno nella stesso bucket. se si decodifica 
	 * solo un bucket 
	 */

	if (num_buc_ep == (sp / m_pIndex->bucket_size_lev2))
	{	// stesso bucket
		unsigned char char_map[ALPHASIZE];
	
		if ((num_buc_ep%2 == 0)	&& ( num_buc_ep%(m_pIndex->bucket_size_lev1 / m_pIndex->bucket_size_lev2)!=0) 
			&& (num_buc_ep != m_pIndex->num_bucs_lev2-1) ) 
			{        // bucket dispari
			state = get_info_b ('\0', sp, occ_b, WHAT_CHAR_IS);
			mod =  m_pIndex->bucket_size_lev2 - (ep % m_pIndex->bucket_size_lev2) - 1;
			diff = ep - sp;

		
			for (i = 0; i < m_pIndex->alpha_size_sb; i++)
			{
				occ_b2[i] = occ_b[i];
				char_map[i] = 0;
			}
			
			c = m_pIndex->mtf_seq + mod;
			int i;
			for (i=0; i < diff; i++,c++){
				occ_b2[*c]++;
				if (char_map[*c] == 0)
				{
					char_map[*c] = 1;
					char_in[char_present++] = *c;
				}
			}

			for (i = 0; i < char_present; i++)
			{
				d = m_pIndex->inv_map_sb[char_in[i]];
				occsp[d] = occ_sb[d] + occ_b[char_in[i]];
				occep[d] = occ_sb[d] + occ_b2[char_in[i]];
				char_in[i] = d;
			}
			return char_present;
		} else {
			
		mod = sp % m_pIndex->bucket_size_lev2;
		diff = ep - sp;
		b2end = mod + diff;	// posizione di ep nel bucket > 1023 => bucket diverso

		state = get_info_b ('\0', ep, occ_b2, WHAT_CHAR_IS);
		for (i = 0; i < m_pIndex->alpha_size_sb; i++)
		{
			occ_b[i] = occ_b2[i];
			char_map[i] = 0;
		}

		mod++;
		c = m_pIndex->mtf_seq + mod;

		for (; mod <= b2end; mod++, c++)
		{
			occ_b[*c]--;	// togli occorrenze b2
			if (char_map[*c] == 0)
			{
				char_map[*c] = 1;
				char_in[char_present++] = *c;
			}
		}

		for (i = 0; i < char_present; i++)
		{
			d = m_pIndex->inv_map_sb[char_in[i]];
			occsp[d] = occ_sb[d] + occ_b[char_in[i]];
			occep[d] = occ_sb[d] + occ_b2[char_in[i]];
			char_in[i] = d;
		}

		return char_present;
	}
	}
	// calcola occorrenze fino a k e a k2
	state = get_info_b ('c', sp, occ_b, WHAT_CHAR_IS);
	if (state < 0)
		return state;	// if error return code

	/* add occ of All chars  */
	remap = 0;		/* rimappa i caratteri dall'alfabeto
				 * locale al superbucket a quello del
				 * testo */
	for (i = 0; i < m_pIndex->alpha_size; i++)
	{
		if (m_pIndex->bool_map_sb[i])
			occsp[i] = occ_sb[i] + occ_b[remap++];
		else
			occsp[i] = occ_sb[i];	/* non occorre nel sb corrente ma
						 			   devi dirgli le occorrenze fin qui !!! */
	}
	// calcola occorrenze fino a k e a k2
	state = get_info_sb (ep, occ_sb2);	// get occ of all chars in 
	// prev superbuckets 
	if (state < 0)
		return state;

	state = get_info_b ('c', ep, occ_b2, WHAT_CHAR_IS);
	if (state < 0)
		return state;	// if error return code

	/*
	 * add occ of All chars 
	 */
	remap = 0;		/* rimappa i caratteri dall'alfabeto
				 * locale al superbucket a quello del
				 * testo */
	for (i = 0; i < m_pIndex->alpha_size; i++)
	{
		if (m_pIndex->bool_map_sb[i])
			occep[i] = occ_sb2[i] + occ_b2[remap++];
		else
			occep[i] = occ_sb2[i];	/* non occorre nel sb corrente ma
						 * devi dirgli le occorrenze fin
						 * qui !!! */
		if (occep[i] != occsp[i])
			char_in[char_present++] = i;

	}

	return char_present;

}


/*
 * Read informations from the header of the superbucket. "pos" is a
 * position in the last column (that is the bwt). We we are interested in
 * the information for the superbucket containing "pos". Initializes the
 * data structures: m_pIndex->inv_map_sb, m_pIndex->bool_map_sb, m_pIndex->alpha_size_sb
 * Returns initialized the array occ[] containing the number of
 * occurrences of the chars in the previous superbuckets. All
 * sb-occurences in the index are stored with log2textsize bits. 
 */
int
CFMIndex::get_info_sb (UINT32 pos, UINT32 * occ)
{

	UINT32 size, sb, *occpoint = occ, offset, i;

	if (pos >= m_pIndex->text_size)
		return FM_SEARCHERR;	// Invalid pos

	// ----- initialize the data structures for I/O-reading
	sb = pos / m_pIndex->bucket_size_lev1;	// superbucket containing pos

	offset = m_pIndex->start_prologue_info_sb;


	// --------- go to the superbucket header 
	if (sb > 0)
	{
		// skip bool map in previous superbuckets
		// skip occ in previous superbuckets (except the first one)
		offset += (sb - 1) * (m_pIndex->alpha_size * sizeof(UINT32)) + (sb * m_pIndex->sb_bitmap_size);
	} 

	fm_init_bit_reader (m_pIndex->compress + offset);	// position of sb header

	/* get bool_map_sb[] */
	for (i = 0; i < m_pIndex->alpha_size; i++)
	{
		fm_bit_read24 (1, m_pIndex->bool_map_sb[i]);
	}

	/* compute alphabet size */
	m_pIndex->alpha_size_sb = 0;
	for (i = 0; i < m_pIndex->alpha_size; i++)
		if (m_pIndex->bool_map_sb[i])
			m_pIndex->alpha_size_sb++;

	/* Invert the char-map for this superbucket */
	for (i = 0, size = 0; i < m_pIndex->alpha_size; i++)
		if (m_pIndex->bool_map_sb[i])
			m_pIndex->inv_map_sb[size++] = (unsigned char) i;

	assert (size == m_pIndex->alpha_size_sb);

	/* for the first sb there are no previous occurrences */
	if (sb == 0)
	{
		for (i = 0; i < m_pIndex->alpha_size; i++, occpoint++)
			*occpoint = 0;
	}
	else
	{
		/* otherwise copy # occ_map in previous superbuckets */
		memcpy(occ, m_pIndex->compress + offset + m_pIndex->sb_bitmap_size, m_pIndex->alpha_size*sizeof(UINT32));

	}

	return FM_OK; 
}


/*
 * Read informations from the header of the bucket and decompresses the
 * bucket if needed. Initializes the data structures: m_pIndex->inv_map_b,
 * m_pIndex->bool_map_b, m_pIndex->alpha_size_b Returns initialized the array
 * occ[] containing the number of occurrences of all the chars since the
 * beginning of the superbucket Explicitely returns the character
 * (remapped in the alphabet of the superbucket) occupying the absolute
 * position pos. The decompression of the bucket when ch does not occur in 
 * the bucket (because of m_pIndex->bool_map_b[ch]=0) is not always carried out. 
 * The parameter "flag" setted to COUNT_CHAR_OCC indicates that we want
 * to count the occurreces of ch; in this case when m_pIndex->bool_map_b[ch]==0
 * the bucket is not decompressed. When the flag is setted to
 * WHAT_CHAR_IS, then ch is not significant and we wish to retrieve the
 * character in position k. In this case the bucket is always
 * decompressed. 
 */

int
CFMIndex::get_info_b (unsigned char ch, UINT32 pos, UINT32 *occ, int flag)
{

	UINT32 buc_start_pos, buc, size, nextbuc = 0;
	int i, offset, is_odd = 0, isnotfirst;
	unsigned char ch_in_pos = 0;

	buc = pos / m_pIndex->bucket_size_lev2;	// bucket containing pos
	assert (buc < (m_pIndex->text_size + m_pIndex->bucket_size_lev2 - 1)
		/ m_pIndex->bucket_size_lev2);
	isnotfirst = buc % (m_pIndex->bucket_size_lev1 / m_pIndex->bucket_size_lev2);

	/* read bucket starting position */
	offset = m_pIndex->start_prologue_info_b + m_pIndex->var_byte_rappr/8 * buc;
	fm_init_bit_reader ((m_pIndex->compress) + offset);
	buc_start_pos = fm_bit_read (m_pIndex->var_byte_rappr);

	if((buc%2 == 0) && (isnotfirst)  && (buc != m_pIndex->num_bucs_lev2-1)) {
		is_odd = 1; // bucket per il quale non sono memorizzate le occorrenze
		nextbuc = fm_bit_read (m_pIndex->var_byte_rappr);	
	}
	
	/* move to the beginning of the bucket */
	/* Se e' un bucket senza occ leggo quelle del successivo */
	if (is_odd) offset = nextbuc;
	else offset = buc_start_pos;
	fm_init_bit_reader ((m_pIndex->compress) + offset);

	/* Initialize properly the occ array */
	if (isnotfirst == 0)
	{
		for (i = 0; i < m_pIndex->alpha_size_sb; i++)
			occ[i] = 0;
	}
	else 
		{
		for (i = 0; i < m_pIndex->alpha_size_sb; i++)
		{	
			occ[i] = fm_integer_decode (m_pIndex->int_dec_bits);
		}
	}
	
	if (is_odd) { // se sono uno senza occ mi posiziono su quello corrente 
	fm_init_bit_reader ((m_pIndex->compress) + buc_start_pos);
	}

	/* get bool char map */
	for (i = 0; i < m_pIndex->alpha_size_sb; i++)
	{	
		fm_bit_read24 (1, m_pIndex->bool_map_b[i]);
	}

	/* get bucket alphabet size and the code of ch in this bucket */
	m_pIndex->alpha_size_b = 0;
	for (i = 0; i < m_pIndex->alpha_size_sb; i++)
		if (m_pIndex->bool_map_b[i])
			m_pIndex->alpha_size_b++;	// alphabet size in the bucket

	/* if no occ of this char in the bucket then skip everything */
	if ((flag == COUNT_CHAR_OCC) && (m_pIndex->bool_map_b[ch] == 0))
		return ((unsigned char) 0);	// dummy return

	/* Invert the char-map for this bucket */
	for (i = 0, size = 0; i < m_pIndex->alpha_size_sb; i++)
		if (m_pIndex->bool_map_b[i])
			m_pIndex->inv_map_b[size++] = (unsigned char) i;

	assert (size == m_pIndex->alpha_size_b);
		
	/* decompress and count CH occurrences on-the-fly */
//	switch (m_pIndex->type_compression)
	//{

	//case MULTIH:		/* multihuffman compression of */
		ch_in_pos = get_b_multihuf(pos, occ,is_odd);
		//break;

	//default:
		//return FM_COMPNOTSUP;
	//}
	return (int) (ch_in_pos); /* char represented in [0..m_pIndex->alpha_size.sb-1] */
}

/*
 * Multi-table-Huffman compressed bucket. Update the array occ[] summing
 * up all occurrencs of the chars in its prefix preceding the absolute
 * position k. Note that ch is a bucket-remapped char. 
 */
unsigned char
CFMIndex::get_b_multihuf(UINT32 k, UINT32 * occ, int is_odd)
{
	UINT32 bpos, i, j, mtf_seq_len;
	unsigned char char_returned;

	bpos = k % m_pIndex->bucket_size_lev2;

	if (is_odd) bpos = m_pIndex->bucket_size_lev2 - bpos - 1;

	if (m_pIndex->alpha_size_b == 1)
	{			/* special case bucket with only one char */
		char_returned = m_pIndex->inv_map_b[0];
		if(is_odd) {
			for (j=0; j <= bpos; j++)
			{
				m_pIndex->mtf_seq[j] = char_returned;
				occ[char_returned]--;
				
			}
			occ[char_returned]++;
	    } else {	
			for (i = 0; i <= bpos; i++)
			{
				m_pIndex->mtf_seq[i] = char_returned;
				occ[char_returned]++;
				
			}
		}
		return char_returned;
	}

	/* Initialize Mtf start */
	for (i = 0; i < m_pIndex->alpha_size_b; i++)
			m_pIndex->mtf[i] = (unsigned char)i;
	
	mtf_seq_len =
			fm_multihuf_decompr (m_pIndex->mtf_seq, m_pIndex->alpha_size_b, bpos+1);
	
	assert (mtf_seq_len > bpos);
	assert (mtf_seq_len <= m_pIndex->bucket_size_lev2);

	/* The chars in the unmtf_bucket are already un-mapped */
	unmtf_unmap (m_pIndex->mtf_seq, bpos + 1);

	/* returning char at bwt-position k --> Inv[] not necessary */
	char_returned = m_pIndex->mtf_seq[bpos];
	assert (char_returned < m_pIndex->alpha_size_sb);
	
	if (is_odd) {
		for (i=0; i < bpos; i++)
			occ[m_pIndex->mtf_seq[i]]--;
	} else {	
	/* update occ[]array */
		for (i = 0; i <= bpos; i++)
			occ[m_pIndex->mtf_seq[i]]++;	/* unmtf_bucket contains chars un-mapped */
	}
	return char_returned;
}


/*
 * Receives in input a bucket in the MTF form, having length len_mtf;
 * returns the original bucket where MTF-ranks have been explicitely
 * resolved. The characters obtained from m_pIndex->mtf[] are UNmapped according
 * to the ones which actually occur into the superbucket. Therefore, the
 * array m_pIndex->inv_map_b[] is necessary to unmap those chars from
 * m_pIndex->alpha_size_b to m_pIndex->alpha_size_sb. 
 */
void
CFMIndex::unmtf_unmap (unsigned char * mtf_seq, int len_mtf)
{
	int i, j, rank;
	unsigned char next;

	/* decode "inplace" mtf_seq */
	for (j = 0; j < len_mtf; j++, mtf_seq++)
	{
		rank = *mtf_seq;
		assert (rank < m_pIndex->alpha_size_b);
		next = m_pIndex->mtf[rank];			/* decode mtf rank */
		*mtf_seq = m_pIndex->inv_map_b[next];	/* apply invamp	*/
	    assert(*mtf_seq < m_pIndex->alpha_size_sb);
	
		for (i = rank; i > 0; i--)		/* update mtf list */
			m_pIndex->mtf[i] = m_pIndex->mtf[i - 1];
		m_pIndex->mtf[0] = next;				/* update mtf[0] */
	}
}


/* Functions to compress buckets */

/* 
   compress and write to file a bucket of length "len" starting at in[0].
   the compression is done as follows:
   first the charatcters are remapped (we expect only a few distinct chars
   in a single bucket) then we use mtf and we compress. 
*/ 
int 
CFMIndex::compress_bucket(unsigned char *in, UINT32 len, UINT16 alphasize) {
	
  UINT16 local_alpha_size, j;
  unsigned char c, local_bool_map[256], local_map[256]; 
 
  /* ---------- compute and write local boolean map ------ */
  for(j=0; j<alphasize; j++){     
    local_bool_map[j]=0;
	local_map[j]=0;
  }
  local_alpha_size=0;
  
  for(j=0;j<len;j++) {             // compute local boolean map
    c = in[j];                     // remapped char
    assert(c<alphasize);                              
    local_bool_map[c] = 1;     
  }

  for(j=0; j<alphasize; j++)      // compute local map
    if(local_bool_map[j])
      local_map[j] = (unsigned char)local_alpha_size++; 
	
  for(j=0;j<alphasize;j++)     // write bool char map to file 
   		if(local_bool_map[j]) {fm_bit_write24(1,1);}
   		else {fm_bit_write24(1,0);} 
  
  for(j=0;j<len;j++)             // remap bucket
    in[j]=local_map[in[j]];
  
  int error = 0;
  switch ( m_pIndex->type_compression ) {
	  case ( MULTIH ):
		   if (local_alpha_size == 1) { 
					fm_bit_flush(); 
					return FM_OK;
					}
		    /* compute mtf picture */
  			mtf_string(in, m_pIndex->mtf_seq, len, local_alpha_size);
			error = fm_multihuf_compr(m_pIndex->mtf_seq, len, local_alpha_size);
			if ( error < 0 ) return error;
			fm_bit_flush(); 
			break;
	  default: 
	  		return FM_COMPNOTSUP;
  	}
 
  return FM_OK;

}

/* Compute Move to Front for string */

void 
CFMIndex::mtf_string(unsigned char *in, unsigned char *out, UINT32 len, UINT16 mtflen)
{
	
  UINT32 i,o,m;
  UINT16 j,h;
  unsigned char c;

  if(pmtf_start == NULL || gmtflen < mtflen)
	{
	if(pmtf_start != NULL)
		{
		free(pmtf_start);
		pmtf_start = NULL;
		}
	gmtflen = 0;
	}
  if(pmtf_start == NULL)
	{
	pmtf_start = (unsigned char *)malloc(sizeof(unsigned char) * (mtflen + 1000));
	gmtflen = mtflen + 1000;
	}
 
  for(j=0; j<mtflen; j++) 
	  pmtf_start[j] = (unsigned char)j;
  
  m = mtflen;   // # of chars in the mtf list
  o = 0;   // # of char in the output string

  for(i=0; i<len; i++) {
    c = in[i];
    /* search c in mtf */
  	for(h=0; h<m; h++)
       	if(pmtf_start[h]==c) break;
	  
    /* c found in mtf[h] */
    out[o++] = (unsigned char) h;   
    for(j=h; j>0; j--)   // update mtf[]
      	 pmtf_start[j] = pmtf_start[j-1];
    pmtf_start[0] = c;
  }
  assert(o<=2*len);
}

int 
CFMIndex::fm_bwt_compress(void) 
{

	int error;
	
	m_pIndex->compress = (unsigned char *)malloc((size_t)((400+floor(1.1*m_pIndex->text_size))*sizeof(unsigned char)));
	if (m_pIndex->compress == NULL) 
			return FM_OUTMEM;

	m_pIndex->compress_size = 0;
	fm_init_bit_writer(m_pIndex->compress, &m_pIndex->compress_size);
	fm_uint_write(m_pIndex->text_size);
	fm_uint_write(m_pIndex->bwt_eof_pos);
	
	/* MTF */
	m_pIndex->mtf_seq = (unsigned char *) malloc (m_pIndex->text_size * sizeof (unsigned char));
	if (m_pIndex->mtf_seq == NULL)
		return FM_OUTMEM;
	
	mtf_string(m_pIndex->bwt, m_pIndex->mtf_seq, m_pIndex->text_size, ALPHASIZE);
	
	free(m_pIndex->bwt);
	m_pIndex->bwt = NULL;
	/* compress rle + huffman */
	error = fm_multihuf_compr(m_pIndex->mtf_seq, m_pIndex->text_size, ALPHASIZE);
	if ( error < 0 ) return error;

	fm_bit_flush(); 
	
	free(m_pIndex->mtf_seq);
	m_pIndex->mtf_seq = NULL;
	m_pIndex->compress = (unsigned char *)realloc(m_pIndex->compress, sizeof(unsigned char)*m_pIndex->compress_size);

	return FM_OK;
}

void 
CFMIndex::fm_unmtf(unsigned char *in, unsigned char *out, int length)
{
  int i,k, pos;
  unsigned char mtf_list[256];

  for (i=0; i<256; i++)  // Initialize the MTF list
    mtf_list[i]= (unsigned char) i;
  for (i=0; i < length; i++) {
	pos = (int) in[i];
	out[i] = mtf_list[pos]; // MTF unencode the next symbol
	for (k=pos; k>0; k--) // move-to-front bwt[i]
		mtf_list[k]=mtf_list[k-1];
	mtf_list[0]=out[i];
  }
}

int 
CFMIndex::fm_bwt_uncompress(void) 
{
	UINT32 i;
	int error;
	m_pIndex->bwt_eof_pos = fm_uint_read();
	
	m_pIndex->mtf_seq = (unsigned char *)malloc(sizeof(unsigned char) * m_pIndex->text_size);
	if (m_pIndex->mtf_seq == NULL) return FM_OUTMEM;
	
	fm_multihuf_decompr (m_pIndex->mtf_seq, ALPHASIZE, m_pIndex->text_size);
	m_pIndex->bwt = (unsigned char *)malloc(sizeof(unsigned char) * m_pIndex->text_size);
	if (m_pIndex->bwt == NULL) return FM_OUTMEM;
		
	/* The chars in the unmtf_bucket are already un-mapped */
	fm_unmtf(m_pIndex->mtf_seq, m_pIndex->bwt, m_pIndex->text_size );
	free(m_pIndex->mtf_seq);
	m_pIndex->mtf_seq = NULL;
	
	/* compute bwt occ */
  	for(i=0;i<ALPHASIZE;i++) m_pIndex->bwt_occ[i]=0;
  	for(i=0;i<m_pIndex->text_size;i++) m_pIndex->bwt_occ[m_pIndex->bwt[i]]++;
  	for(i=1;i<ALPHASIZE;i++) m_pIndex->bwt_occ[i] += m_pIndex->bwt_occ[i-1];
  	assert(m_pIndex->bwt_occ[ALPHASIZE-1]==m_pIndex->text_size);
  	for(i=ALPHASIZE-1;i>0;i--) m_pIndex->bwt_occ[i] = m_pIndex->bwt_occ[i-1];
	m_pIndex->bwt_occ[0] = 0;

	error = fm_compute_lf();
	if (error < 0) return error;

	error = fm_invert_bwt();
	if (error < 0) return error;

	return FM_OK;
}


/*
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 bwt-based
 * compression and indexing
 
 hufbzip.c
 compression and decompression
 * using multiple huffman tables
 as in the bzip2 compressors.
 
 P.
 * Ferragina & G. Manzini, 10 June 2000
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
 */ 

	
	 /* ********************************************************************
	 rle+compression of a string using Huffman with multiple tables 
	 input
	 int   len         size of mtf sequence
	 unsigned char *in         input mtf sequence
	 int   alpha_size   size of the alphabet
	 output
	 the compressed string is written in the output file 
	 ******************************************************************* */ 
int 
CFMIndex::fm_multihuf_compr (unsigned char * in, int len, int alpha_size) 
{
int v, t, i, j, gs, ge, totc, bt, bc, iter;
int nSelectors, minLen, maxLen, new_len;
int nGroups;
UINT16 cost[BZ_N_GROUPS];
int fave[BZ_N_GROUPS];
UINT16 * mtfv;
unsigned char * selector;

mtfv = (UINT16 *) malloc ((len + 1) * sizeof (UINT16));
if (mtfv == NULL)
	return FM_OUTMEM;
	
		// encode sequences of 0's using 1-2 coding
new_len = 0;

	{
	int c, z = 0;
	for (i = 0; i < len; i++)
		{
		c = in[i];
		assert (c < alpha_size);
		if (c == 0)
			z++;
		else
			{
				
					/*
					 * ----- check if there are pending zeores ---- 
					 */ 
			if (z > 0)
				{	// 1-2 encoding 
				z++;	// write z+1 in binary least sign bit
					// first 
				while (z > 1)
					{
					mtfv[new_len++] =	(z & 1) ? BZ_RUNB :	BZ_RUNA;
					z = z >> 1;
					}
				z = 0;
				}
			mtfv[new_len++] = (UINT16) c + 1;
			}
		}
		
			// ---- there could be some pending zeroes
	if (z > 0)
		{
			
		z++;	// write z+1 in binary least sign bit
			// first 
		while (z > 1)
			{
			mtfv[new_len++] = (z & 1) ? BZ_RUNB : BZ_RUNA;
			z = z >> 1;
			}
		}
	}
	
mtfv[new_len++] = alpha_size + 1;	// end of block
alpha_size += 2;	// 2 new symbols have been used
if (__Fm_Verbose > 2)
	fprintf (stderr, "block size after MTF & 1-2 coding: %d\n",
			  new_len);
	

		// init mtf_freq[]
for (i = 0; i < alpha_size; i++)
	mtf_freq[i] = 0;
	
for (i = 0; i < new_len; i++)
	mtf_freq[mtfv[i]]++;
	
		// init huf_len[][]
for (t = 0; t < BZ_N_GROUPS; t++)
	for (v = 0; v < alpha_size; v++)
		huf_len[t][v] = BZ_GREATER_ICOST;
	
		// alloc selector[]
selector =	(unsigned char *) malloc ((1 + new_len / BZ_G_SIZE) * sizeof (unsigned char));
if (selector == NULL)
	return FM_OUTMEM;
	

   /*--- Decide how many coding tables to use ---*/ 
assert (new_len > 0);
if (new_len < 200)
	nGroups = 2;
else
	if (new_len < 600)
		nGroups = 3;
	else
		if (new_len < 1200)
			nGroups = 4;
		else
			if (new_len < 2400)
				nGroups = 5;
			else
				nGroups = 6;
	

   /*--- Generate an initial set of coding tables ---*/ 
		// each table uses BZ_LESSER_ICOST for a group of consecutive 
		// chars (gs to ge) and BZ_GREATER_ICOST for the others chars
	{
	int nPart, remF, tFreq, aFreq;
	nPart = nGroups;
	remF = new_len;
	gs = 0;
	while (nPart > 0)
		{
		tFreq = remF / nPart;
		ge = gs - 1;
		aFreq = 0;
		while (aFreq < tFreq && ge < alpha_size - 1)
			{
			ge++;
			aFreq += mtf_freq[ge];
			}
		if (ge > gs 
			     &&nPart != nGroups && nPart != 1 
			     &&((nGroups - nPart) % 2 == 1))
			{
			aFreq -= mtf_freq[ge];
			ge--;
			}
		if (__Fm_Verbose > 2)
			fprintf (stderr,
					  "      initial group %d, [%d .. %d], has %d syms (%4.1f%%)\n", 
						nPart, gs, ge, aFreq,(100.0 * (float) aFreq)/(float) (new_len));
		for (v = 0; v < alpha_size; v++)
			if (v >= gs && v <= ge)
				huf_len[nPart - 1][v] =	BZ_LESSER_ICOST;
			else
				huf_len[nPart - 1][v] =	BZ_GREATER_ICOST;
		nPart--;
		gs = ge + 1;
		remF -= aFreq;
		}
	}
	

   /*--- 
      Iterate up to BZ_N_ITERS times to improve the tables.
   ---*/ 
for (iter = 0; iter < BZ_N_ITERS; iter++)
	{
	for (t = 0; t < nGroups; t++)
		fave[t] = 0;
	for (t = 0; t < nGroups; t++)
		for (v = 0; v < alpha_size; v++)
			rfreq[t][v] = 0;
		
	nSelectors = 0;
	totc = 0;
	gs = 0;
	while (true)
		{
			
				/*
				 * Set group start & end marks. -- 
				 */ 
		if (gs >= new_len)
			break;
		ge = gs + BZ_G_SIZE - 1;	// size is at most BZ_G_SIZE
		if (ge >= new_len)
			ge = new_len - 1;
			
	 /*-- 
            Calculate the cost of this group as coded
            by each of the coding tables.
         --*/ 
		for (t = 0; t < nGroups; t++)
			cost[t] = 0;
		if (nGroups == 6)
			{
			register UINT16 cost0,cost1, cost2, cost3, cost4, cost5;
			cost0 = cost1 = cost2 = cost3 = cost4 =	cost5 = 0;
			for (i = gs; i <= ge; i++)
				{
				UINT16 icv = mtfv[i];
				cost0 += huf_len[0][icv];
				cost1 += huf_len[1][icv];
				cost2 += huf_len[2][icv];
				cost3 += huf_len[3][icv];
				cost4 += huf_len[4][icv];
				cost5 += huf_len[5][icv];
				}
			cost[0] = cost0;
			cost[1] = cost1;
			cost[2] = cost2;
			cost[3] = cost3;
			cost[4] = cost4;
			cost[5] = cost5;
			}
		else
			{
			for (i = gs; i <= ge; i++)
				{
				UINT16 icv = mtfv[i];
				for (t = 0; t < nGroups; t++)
					cost[t] += huf_len[t][icv];
				}
			}
			
	 /*-- 
            Find the coding table which is best for this group,
            and record its identity in the selector table.
         --*/ 
		bc = 999999999;
		bt = -1;
		for (t = 0; t < nGroups; t++)
		if (cost[t] < bc)
			{
			bc = cost[t];
			bt = t;
			};
		totc += bc;
		fave[bt]++;
		selector[nSelectors++] = bt;
			
	 /*-- 
            Increment the symbol frequencies for the selected table.
          --*/ 
		for (i = gs; i <= ge; i++)
			rfreq[bt][mtfv[i]]++;
		gs = ge + 1;	// consider next group 
		}
	if (__Fm_Verbose > 2)
		{
			
fprintf (stderr,
				  "      pass %d: size is %d, grp uses are ",
				  
iter + 1, totc / 8);
			
for (t = 0; t < nGroups; t++)
				
fprintf (stderr, "%d ", fave[t]);
			
fprintf (stderr, "\n");
		
}
		
      /*--
        Recompute the tables based on the accumulated frequencies.
      --*/ 
			for (t = 0; t < nGroups; t++)
			
hbMakeCodeLengths (&(huf_len[t][0]), (UINT32 *)&(rfreq[t][0]),
					    alpha_size, 20);
	
}
	

   /*--- Assign actual codes for the tables. --*/ 
		for (t = 0; t < nGroups; t++)
	{
		
minLen = 32;
		
maxLen = 0;
		
for (i = 0; i < alpha_size; i++)
		{
			
if (huf_len[t][i] > maxLen)
				maxLen = huf_len[t][i];
			
if (huf_len[t][i] < minLen)
				minLen = huf_len[t][i];
		
}
		
assert (!(maxLen > 20));
		
assert (!(minLen < 1));
		
hbAssignCodes ((UINT32 *)&(huf_code[t][0]), &(huf_len[t][0]), 
minLen,
				maxLen, alpha_size);
	
}
	

		// write coding tables (i.e codeword length).
		assert (nGroups < 8);
	
fm_bit_write (3, nGroups);
	
for (t = 0; t < nGroups; t++)
	{
		
int curr = huf_len[t][0];
		
fm_bit_write (5, curr);
		
for (i = 0; i < alpha_size; i++)
		{
			
while (curr < huf_len[t][i])
			{
				fm_bit_write (2, 2);
				curr++;	/* 10 */
			};
			
while (curr > huf_len[t][i])
			{
				fm_bit_write (2, 3);
				curr--;	/* 11 */
			};
			
fm_bit_write (1, 0);
		
}
	
}
	

   /*--- write selectors and compressed data ---*/ 
	{
		
int sel = 0;
		
unsigned char pos[BZ_N_GROUPS], ll_i, tmp2, tmp;
		
for (i = 0; i < nGroups; i++)
			pos[i] = i;
		
gs = 0;
		
while (true)
		{
			
if (gs >= new_len)
				break;
			
ge = gs + BZ_G_SIZE - 1;
			
if (ge >= new_len)
				ge = new_len - 1;	// establish group boundaries
			assert (selector[sel] < nGroups);
			
			{
				
ll_i = selector[sel];	// get mtf rank for selector
				j = 0;
				
tmp = pos[j];
				
while (ll_i != tmp)
				{
					
j++;
					
tmp2 = tmp;
					tmp = pos[j];
					pos[j] = tmp2;
				
};
				
pos[0] = tmp;
				
fm_bit_write (j + 1, 1);	// write selector mtf rank in
				// unary 
			}
			
for (i = gs; i <= ge; i++)
			{
				
assert (mtfv[i] < alpha_size);
				
fm_bit_write (huf_len[selector[sel]][mtfv[i]],
					    
huf_code[selector[sel]][mtfv
								     [i]]);
			
}
			
gs = ge + 1;
			
sel++;
		
}
		
assert (sel == nSelectors);
	
}
	
free (selector);
	
free (mtfv);
	
return FM_OK;

}



	/*
	 *********************************************************
	 decode a unary code read from Infile. the output is the
	 # of zeroes we see before we see a 1 
	 1 --> 0
	 01 --> 1
	 001 --> 2  
	 etc.
	 ********************************************************* */ 
int
CFMIndex::decode_unary (void) 
{
	
int t, i = 0;
	

	do
	{
		
fm_bit_read24 (1, t);
		
if (t != 0)
			break;
		
i++;
	
}
	while (1);
	
return i;

}




	/*
	 ********************************************************************
	 this procedures reads a bucket from file decodes it and writes 
	 it to dest[] (which should be of the appropriate size).
	 The decoding stops when an EOB is encountered or >= limit bytes
	 have been decoded. that is, when >= limit chars have been written
	 to dest the decompression terminates and the procedure returns.
	 Note that more than limit chars can be
	 decoded, so dest() should be large enough to contain the 
	 complete bucket. the procedure returns the number of chars written to
	 dest (which can be less than limit (if a EOB is encountered)) 
	 ******************************************************************** */ 
int
CFMIndex::fm_multihuf_decompr (unsigned char * dest, int alpha_size, int limit) 
{
int t, i, j, minLen, maxLen, len, nGroups;
	

alpha_size += 2;	// we temporarily use a larger alphabet
	
		// get number of groups
		fm_bit_read24 (3, nGroups);
	
  /*--- get the coding tables ---*/ 
	{
		
int curr, uc;
		

for (t = 0; t < nGroups; t++)
		{
			
fm_bit_read24 (5, curr);
			
for (i = 0; i < alpha_size; i++)
			{
				
while (true)
				{
					
if (curr < 1 || curr > 20)
						
return FM_DECERR;
					
fm_bit_read24 (1, uc);
					
if (uc == 0)
						break;
					
fm_bit_read24 (1, uc);
					
if (uc == 0)
						curr++;
					else
						curr--;
				
}
				
huf_len[t][i] = curr;
			
}
		
}
		

    /*--- Create the Huffman decoding tables ---*/ 
			for (t = 0; t < nGroups; t++)
		{
			
minLen = 32;
			
maxLen = 0;
			
for (i = 0; i < alpha_size; i++)
			{
				
if (huf_len[t][i] > maxLen)
					maxLen = huf_len[t][i];
				
if (huf_len[t][i] < minLen)
					minLen = huf_len[t][i];
			
}
			
hbCreateDecodeTables ((UINT32 *)&(huf_limit[t][0]),(UINT32 *)&(huf_base[t][0]),(UINT32 *)&(huf_perm[t][0]),&(huf_len[t][0]),minLen,maxLen, alpha_size);
			
huf_minLens[t] = minLen;
		
}
	
}
	

   /*------- uncompress data -------*/ 
	{
		
int rle_sofar, run, next, rank, gSel, to_be_read;
		
int zn, zj, zvec, *gLimit, *gPerm, *gBase;
		
unsigned char pos[BZ_N_GROUPS], gMinlen = 0;
		

gLimit = gPerm = gBase = NULL;	// to avoid annoying
		// compiler warnings
		for (i = 0; i < nGroups; i++)
			pos[i] = i;
		
len = 0;
		rle_sofar = 0;
		
to_be_read = 0;
		
while (true)
		{
			
if (to_be_read == 0)
			{
				
to_be_read = BZ_G_SIZE;
				
rank = decode_unary ();	// get mtf rank of new group
				assert (rank < nGroups);
				
gSel = pos[rank];
				
for (j = rank; j > 0; j--)
					pos[j] = pos[j - 1];
				
pos[0] = (unsigned char) gSel;
				
					// get tables for this group
					gMinlen = huf_minLens[gSel];
				
gLimit = &(huf_limit[gSel][0]);
				
gPerm = &(huf_perm[gSel][0]);
				
gBase = &(huf_base[gSel][0]);
			
}
			
to_be_read--;
			
				// get next huffman encoded char
				zn = gMinlen;
			
				// zvec = bit_read(zn);
				fm_bit_read24 (zn, zvec);
			
while (zvec > gLimit[zn])
			{
				
zn++;
				
					// zj=bit_read(1);
					fm_bit_read24 (1, zj);
				
zvec = (zvec << 1) | zj;
			
};
			
next = gPerm[zvec - gBase[zn]];
			
				// decode next
				assert (next < alpha_size);
			
if (next == alpha_size - 1)
				break;	// end of bucket
			if (next == BZ_RUNA)
			{	// 0 of a 1-2 encoding
				run = 1 << rle_sofar;
				
for (j = 0; j < run; j++)
					dest[len++] = 0;
				
rle_sofar++;
			
}
			
			else if (next == BZ_RUNB)
			{	// 1 of a 1-2 encoding 
				run = 2 << rle_sofar;
				
for (j = 0; j < run; j++)
					dest[len++] = 0;
				
rle_sofar++;
			
}
			
			else
			{
				
dest[len++] = next - 1;
				
rle_sofar = 0;
			
}
			
if (len >= limit)
				return len;	// only line added to stop when >= limit
		}		// chars have been decoded
	}
	
return len;

}




	/*
	 ************************************************************
	 uncompress the bucket which starts at the current position
	 of infile. the bucket is "len" bytes long and should
	 be written in array dest[]
	 ************************************************************ */ 
int
CFMIndex::fm_uncompress_bucket_multihuf (unsigned char * dest, int len, int alpha_size)
{
int k, j, i, aux_len, rank, local_alpha_size, next;
	
unsigned char mtf[256], inv_local_map[256];
	

		/*
		 * ---------- read local boolean map and compute inverse map
		 * ------ 
		 */ 
		local_alpha_size = 0;
	
for (k = 0; k < alpha_size; k++)
		
if (fm_bit_read (1))
			
inv_local_map[local_alpha_size++] = k;
if (local_alpha_size == 1)
	{			// use this when you have only one 
		// char in bucket
		for (i = 0; i < len; i++)
			dest[i] = inv_local_map[0];
		
return FM_OK;
}
	

		/*
		 * Set initial MTF 
		 */ 
		for (k = 0; k < local_alpha_size; k++)
	{
		
mtf[k] = k;
	
}
	

		/*
		 * ------- decode multiple huffman codes ----- 
		 */ 
		aux_len = fm_multihuf_decompr (dest, local_alpha_size, len);
	
if (aux_len != len)
		
return FM_DECERR;
	

		/*
		 * ------ decode *inplace* mtf_seq -------------------- 
		 */ 
		for (j = 0; j < len; j++)
	{
		
rank = dest[j];
		
assert (rank < local_alpha_size);
		
next = mtf[rank];	// decode mtf rank 
		for (i = rank; i > 0; i--)	// update mtf list
			mtf[i] = mtf[i - 1];
		
mtf[0] = next;	// update mtf[0]
		dest[j] = inv_local_map[next];	// apply invamp
		assert (dest[j] < alpha_size);
	
}
	
return FM_OK;

}





/*
 * Funzioni per la scrittura di bit in memoria 
 */

/*
 * Mem e' il puntatore al primo byte da scrivere e pos_mem riporta il
 * numero di byte scritti. Non Inizializza Buffer. 
 */

void
CFMIndex::fm_init_bit_writer (unsigned char * mem, UINT32 * pos_mem)
{
	__MemAddress = mem;
	__Num_Bytes = pos_mem;
	fm_init_bit_buffer ();
}

void
CFMIndex::fm_init_bit_reader (unsigned char * mem)
{
	__MemAddress = mem;
	__pos_read = 0;
	__Num_Bytes = &__pos_read;
	fm_init_bit_buffer ();
}

/*
 * -----------------------------------------------------------------------------
 * Funzioni per leggere/scrivere piu' di 24 bits usando le funzioni fondamentali
 * ----------------------------------------------------------------------------- 
 */

// ****** Write in Bit_buffer n bits taken from vv (possibly n > 24) 
void
CFMIndex::fm_bit_write (int n, UINT32 vv)
{
	UINT32 v = (UINT32) vv;

	assert (n <= 32);
	if (n > 24)
	{
		fm_bit_write24 ((n - 24), (v >> 24 & 0xffL));
		fm_bit_write24 (24, (v & 0xffffffL));
	}
	else
	{
		fm_bit_write24 (n, v);
	}
}

// ****** Read n bits from Bit_buffer 
int
CFMIndex::fm_bit_read (int n)
{
	UINT32 u = 0;
	int i;
	assert (n <= 32);
	if (n > 24)
	{
		fm_bit_read24 ((n - 24), i);
		u = i << 24;
		fm_bit_read24 (24, i);
		u |= i;
		return ((int) u);
	}
	else
	{
		fm_bit_read24 (n, i);
		return i;
	}
}

void CFMIndex::fm_bit_write24(int bits, UINT32 num) {
					
  	assert(__Bit_buffer_size<8); 
  	assert(bits>0 && bits<=24);
  	assert( num < 1u << bits );	
	__Bit_buffer_size += (bits); 
  	__Bit_buffer |= (num << (32 - __Bit_buffer_size));
	while (__Bit_buffer_size >= 8) {
		*__MemAddress = (unsigned char)(__Bit_buffer>>24);
		__MemAddress++;
    	(*__Num_Bytes)++;
    	__Bit_buffer <<= 8;
    	__Bit_buffer_size -= 8;
		} 
}

/*
 * -----------------------------------------------------------------------------
 * Funzioni per leggere/scrivere piu' di 4 Bytes usando le funzioni fondamentali
 * ----------------------------------------------------------------------------- 
 */

// ****** Write in Bit_buffer four bytes 
void
CFMIndex::fm_uint_write (UINT32 uu)
{
	UINT32 u = (UINT32) uu;
	fm_bit_write24 (8, ((u >> 24) & 0xffL));
	fm_bit_write24 (8, ((u >> 16) & 0xffL));
	fm_bit_write24 (8, ((u >> 8) & 0xffL));
	fm_bit_write24 (8, (u & 0xffL));

}

// ***** Return 32 bits taken from Bit_buffer
UINT32
CFMIndex::fm_uint_read (void)
{

	UINT32 u;
	int i;
	fm_bit_read24 (8, i);
	u = i << 24;
	fm_bit_read24 (8, i);
	u |= i << 16;
	fm_bit_read24 (8, i);
	u |= i << 8;
	fm_bit_read24 (8, i);
	u |= i;
	return ((int) u);
}

/*
 * -----------------------------------------------------------------------------
 * Funzione Bit Flush.
 * ----------------------------------------------------------------------------- 
 */

// ***** Complete with zeroes the first byte of Bit_buffer 
// ***** This way, the content of Bit_buffer is entirely flushed out

void
CFMIndex::fm_bit_flush (void)
{
	if (__Bit_buffer_size != 0)
		fm_bit_write24 ((8 - (__Bit_buffer_size % 8)), 0);	// pad with zero !
}


UINT32
CFMIndex::fm_integer_decode (unsigned short int headbits)
{
	int k, i;
	fm_bit_read24 (headbits, i);
	/*
	 * casi base speciali 
	 */

	if (i == 0)
		return 0;
	if (i == 1)
		return 1;
	if (i == 2)
	{
		fm_bit_read24 (1, i);
		k = 2 + i;
		return k;
	}
	/*
	 * i indica il numero di bits che seguono num = 2^i + k 
	 */
	i--;
	/*
	 * bit_read24 e' corretta fino a superbuckets di dimensione inferiore
	 * a 16384*1024. E' comunque meglio utilizzare questa perche molto
	 * piu' veloce 
	 */
	fm_bit_read24 (i, k);
	return ((1 << i) + k);
}



/*-------------------------------------------------------------*/
/*--- Huffman coding low-level stuff                        ---*/
/*---                                             huffman.c ---*/
/*-------------------------------------------------------------*/

/*--
  This file is a part of bzip2 and/or libbzip2, a program and
  library for lossless, block-sorting data compression.

  Copyright (C) 1996-1999 Julian R Seward.  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

  2. The origin of this software must not be misrepresented; you must 
     not claim that you wrote the original software.  If you use this 
     software in a product, an acknowledgment in the product 
     documentation would be appreciated but is not required.

  3. Altered source versions must be plainly marked as such, and must
     not be misrepresented as being the original software.

  4. The name of the author may not be used to endorse or promote 
     products derived from this software without specific prior written 
     permission.

  THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
  OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Julian Seward, Cambridge, UK.
  jseward@acm.org
  bzip2/libbzip2 version 0.9.5 of 24 May 1999

  This program is based on (at least) the work of:
     Mike Burrows
     David Wheeler
     Peter Fenwick
     Alistair Moffat
     Radford Neal
     Ian H. Witten
     Robert Sedgewick
     Jon L. Bentley

  For more information on these sources, see the manual.
--*/



/*---------------------------------------------------*/
#define WEIGHTOF(zz0)  ((zz0) & 0xffffff00)
#define DEPTHOF(zz1)   ((zz1) & 0x000000ff)
#define MYMAX(zz2,zz3) ((zz2) > (zz3) ? (zz2) : (zz3))

#define ADDWEIGHTS(zw1,zw2)                           \
   (WEIGHTOF(zw1)+WEIGHTOF(zw2)) |                    \
   (1 + MYMAX(DEPTHOF(zw1),DEPTHOF(zw2)))

#define UPHEAP(z)                                     \
{                                                     \
   UINT32 zz, tmp;                                     \
   zz = z; tmp = heap[zz];                            \
   while (weight[tmp] < weight[heap[zz >> 1]]) {      \
      heap[zz] = heap[zz >> 1];                       \
      zz >>= 1;                                       \
   }                                                  \
   heap[zz] = tmp;                                    \
}

#define DOWNHEAP(z)                                   \
{                                                     \
   UINT32 zz, yy, tmp;                                 \
   zz = z; tmp = heap[zz];                            \
   while (true) {                                     \
      yy = zz << 1;                                   \
      if (yy > nHeap) break;                          \
      if (yy < nHeap &&                               \
          weight[heap[yy+1]] < weight[heap[yy]])      \
         yy++;                                        \
      if (weight[tmp] < weight[heap[yy]]) break;      \
      heap[zz] = heap[yy];                            \
      zz = yy;                                        \
   }                                                  \
   heap[zz] = tmp;                                    \
}


/*---------------------------------------------------*/
void 
CFMIndex::hbMakeCodeLengths ( unsigned char *len, 
                         UINT32 *freq,
                         UINT32 alphaSize,
                         UINT32 maxLen )
{
   /*--
      Nodes and heap entries run from 1.  Entry 0
      for both the heap and nodes is a sentinel.
   --*/
   UINT32 nNodes, nHeap, n1, n2, i, j, k;
   bool  tooLong;

   INT32 heap   [ BZ_MAX_ALPHA_SIZE + 2 ];
   INT32 weight [ BZ_MAX_ALPHA_SIZE * 2 ];
   INT32 parent [ BZ_MAX_ALPHA_SIZE * 2 ]; 

   for (i = 0; i < alphaSize; i++)
      weight[i+1] = (freq[i] == 0 ? 1 : freq[i]) << 8;

   while (true) {

      nNodes = alphaSize;
      nHeap = 0;

      heap[0] = 0;
      weight[0] = 0;
      parent[0] = -2;

      for (i = 1; i <= alphaSize; i++) {
         parent[i] = -1;
         nHeap++;
         heap[nHeap] = i;
         UPHEAP(nHeap);
      }

      assert( nHeap < (BZ_MAX_ALPHA_SIZE+2));
   
      while (nHeap > 1) {
         n1 = heap[1]; heap[1] = heap[nHeap]; nHeap--; DOWNHEAP(1);
         n2 = heap[1]; heap[1] = heap[nHeap]; nHeap--; DOWNHEAP(1);
         nNodes++;
         parent[n1] = parent[n2] = nNodes;
         weight[nNodes] = ADDWEIGHTS(weight[n1], weight[n2]);
         parent[nNodes] = -1;
         nHeap++;
         heap[nHeap] = nNodes;
         UPHEAP(nHeap);
      }

      assert( nNodes < (BZ_MAX_ALPHA_SIZE * 2));

      tooLong = false;
      for (i = 1; i <= alphaSize; i++) {
         j = 0;
         k = i;
         while (parent[k] >= 0) { k = parent[k]; j++; }
         len[i-1] = j;
         if (j > maxLen) tooLong = true;
      }
      
      if (! tooLong) break;

      for (i = 1; i < alphaSize; i++) {
         j = weight[i] >> 8;
         j = 1 + (j / 2);
         weight[i] = j << 8;
      }
   }
}


/*---------------------------------------------------*/
void 
CFMIndex::hbAssignCodes ( UINT32 *code,
                     unsigned char *length,
                     UINT32 minLen,
                     UINT32 maxLen,
                     UINT32 alphaSize )
{
   UINT32 n, vec, i;

   vec = 0;
   for (n = minLen; n <= maxLen; n++) {
      for (i = 0; i < alphaSize; i++)
         if (length[i] == n) { code[i] = vec; vec++; };
      vec <<= 1;
   }
}


/*---------------------------------------------------*/
void 
CFMIndex::hbCreateDecodeTables ( UINT32 *limit,
                            UINT32 *base,
                            UINT32 *perm,
                            unsigned char *length,
                            UINT32 minLen,
                            UINT32 maxLen,
                            UINT32 alphaSize )
{
   UINT32 pp, i, j, vec;

   pp = 0;
   for (i = minLen; i <= maxLen; i++)
      for (j = 0; j < alphaSize; j++)
         if (length[j] == i) { perm[pp] = j; pp++; };

   for (i = 0; i < BZ_MAX_CODE_LEN; i++) base[i] = 0;
   for (i = 0; i < alphaSize; i++) base[length[i]+1]++;

   for (i = 1; i < BZ_MAX_CODE_LEN; i++) base[i] += base[i-1];

   for (i = 0; i < BZ_MAX_CODE_LEN; i++) limit[i] = 0;
   vec = 0;

   for (i = minLen; i <= maxLen; i++) {
      vec += (base[i+1] - base[i]);
      limit[i] = vec-1;
      vec <<= 1;
   }
   for (i = minLen + 1; i <= maxLen; i++)
      base[i] = ((limit[i-1] + 1) << 1) - base[i];
}


/*-------------------------------------------------------------*/
/*--- end                                         huffman.c ---*/
/*-------------------------------------------------------------*/




/*
 * Allocates snippet (which must be freed by the caller) and writes
 * text[from..to] on it. Returns in snippet_length the length of the
 * snippet actually extracted (that could be less than to-from+1 if to is
 * larger than the text size). 
 */
#define NEARCLOSE (200) // valore da testare

int 
CFMIndex::extract(UINT32 from, UINT32 to, unsigned char **dest, 
			UINT32 *snippet_length) {

	UINT32 written, numchar;
	UINT32 row = 0; /* numero di riga corrispondente all'ultima position */
	UINT32 scarto = 0;  	/* lo scarto tra la posizione richiesta e la posizione 
			    	       successiva divisibile per m_pIndex->skip */
	UINT32 i,j;
	UINT32 pos_text = 0; /* last readen position */
	unsigned char * text;

	*dest = NULL; 
	*snippet_length = 0;
	
	if ((from >= m_pIndex->text_size) || (from >= to))
			return FM_OK; // Invalid Position
		
	to = MIN(to, m_pIndex->text_size-1);
	
	if(m_pIndex->smalltext) { //uses Boyer-Moore algorithm
		*snippet_length = to-from+1;
		*dest = (unsigned char *)new unsigned char [sizeof(unsigned char)*(*snippet_length)];
		if (*dest==NULL) 
			{
			*snippet_length = 0;
			return FM_OUTMEM;
			}
		memcpy(*dest, m_pIndex->text+from, (*snippet_length)*sizeof(unsigned char));
		return FM_OK;
	}
	
	UINT32 real_text_size;
	if(m_pIndex->skip>1) real_text_size = m_pIndex->text_size-m_pIndex->num_marked_rows; 
	else real_text_size = m_pIndex->text_size;
		
	if ((from == 0) && (to == real_text_size-1)) { // potrebbe essere conveniente anche se inferiore
			int error = fm_unbuild(dest, snippet_length);
			return error;
	}
	
	if (m_pIndex->skip == 0)  
			return FM_NOMARKEDCHAR;
		
	if (to >= m_pIndex->text_size) 
			to = m_pIndex->text_size-1;
				
	numchar = to - from + 1;

	
	/* Legge sequenzalimente le posizioni marcate. E' garantita la presenza 
	   di tutte le posizioni divisibili per m_pIndex->skip */
	
	scarto = m_pIndex->skip - (to%m_pIndex->skip);
	//if (scarto == m_pIndex->skip) scarto = 0;
	
	UINT32 to_new = to + scarto-1;
	if (to_new >= real_text_size-m_pIndex->skip-1) { // vicina alla fine del testo la riga e' la 0
			row = m_pIndex->bwt_eof_pos;
			scarto = real_text_size - to-1;
	}
	/* Inizia a leggere le posizioni finche non trova to_new */
	fm_init_bit_reader(m_pIndex->start_prologue_occ);
	for(i=0; i<m_pIndex->num_marked_rows; i++) { 
		pos_text = fm_bit_read(m_pIndex->log2textsize);

		if((pos_text >= to_new) &&(pos_text < to_new+NEARCLOSE)){
							scarto += pos_text-to_new;
							row = i + m_pIndex->occcharinf;
				  			break;
		}
	}		

	/* 
		Prendi il testo andando all'indietro da row per il numero di caratteri 
	   	che servono per arrivare alla posizione richiesta 
	 */
	text = (unsigned char *)malloc((numchar+scarto)*sizeof(unsigned char));
	if (text == NULL) 
			return FM_OUTMEM;

	/* conosco il primo carattere e' special_char della colonna F in quanto marcato! */
	if (row!=m_pIndex->bwt_eof_pos) {
			text[0] = m_pIndex->inv_char_map[m_pIndex->subchar]; // specifico del tipo di marcamento
			written = go_back(row, scarto + numchar - 1, text+1);
			written++;
	} else written = go_back(row, scarto + numchar, text);
	//fprintf(stderr, "Bucket dec %lu Pos cercata %lu posizione originale %lu scarto %lu\n",written, to_new, to, scarto);

	numchar = MIN(written, numchar); 
	unsigned char *desti = (unsigned char *)new unsigned char [numchar * sizeof(unsigned char)];

	if (desti == NULL) 
			return FM_OUTMEM;
	
	/* reverse the string from the end of text */
	for (j = 0, i = written-1; j<numchar; i--, j++) 
		desti[j] = text[i];
		
	free(text);	
	*snippet_length = numchar;
	*dest = desti;
	
	return FM_OK;
}

/*
   read len chars before the position given by row using the LF mapping
   stop if the beginning of the file is encountered
   return the number of chars actually read 
*/
UINT32 
CFMIndex::go_back(UINT32 row, UINT32 len, unsigned char *dest) {
	
  UINT32 written, curr_row, n, occ_sb[256], occ_b[256];
  unsigned char c, c_sb, cs;
 
  if (row != m_pIndex->bwt_eof_pos) curr_row = EOF_shift(row);
  else curr_row = 0;
  

  for( written=0; written < len; ) {
  
    // fetches info from the header of the superbucket
    get_info_sb(curr_row, occ_sb);  
    // fetches occ into occ_b properly remapped and returns
    // the remapped code for occ_b of the char in the  specified position
    c = get_info_b(NULL_CHAR, curr_row,occ_b, WHAT_CHAR_IS);  
    assert(c < m_pIndex->alpha_size_sb);
  
    c_sb = m_pIndex->inv_map_sb[c];
    assert(c_sb < m_pIndex->alpha_size);
	cs = c_sb;
  
		if ((m_pIndex->skip <= 1) || (cs != m_pIndex->specialchar)) { //skip special char		
    		dest[written++] = m_pIndex->inv_char_map[cs]; 	// store char    
	}
    n = occ_sb[c_sb] + occ_b[c];         	    // # of occ before curr_row

	curr_row = m_pIndex->bwt_occ[c_sb] + n - 1; // get next row
    if(curr_row == m_pIndex->bwt_eof_pos) break;    
    curr_row = EOF_shift(curr_row);        
	
  }

  return written;
}


/*
   write in dest[] the first "len" chars of "row" (unless
   the EOF is encountered).
   return the number of chars actually read 
*/
UINT32 
CFMIndex::go_forw(UINT32 row, UINT32 len, unsigned char *dest) {
	

  UINT32 written;
  unsigned char c,cs;

  for(written=0;written<len; ) {
  
    c = get_firstcolumn_char(row);
    assert(c < m_pIndex->alpha_size);
	cs = c;
	if ((m_pIndex->skip <= 1) || (cs != m_pIndex->specialchar)) { // skip special char 
    	dest[written++] = m_pIndex->inv_char_map[cs];
	}
    // compute the first to last mapping
    row = fl_map(row,c);
    // adjust row to take the EOF symbol into account
	if(row == 0) break; // row = -1
    if(row <= m_pIndex->bwt_eof_pos) row -= 1;
     }
 
  return written;
}

/*
	compute the first-to-last map using binary search
*/
UINT32 
CFMIndex::fl_map(UINT32 row, unsigned char ch) {

  UINT32 i, n, rank, first, last, middle;
  UINT32 occ_sb[ALPHASIZE], occ_b[ALPHASIZE];
  unsigned char c_b,c_sb;              // char returned by get_info
  unsigned char ch_b;

  // rank of c in first column 
  rank = 1 + row - m_pIndex->bwt_occ[ch];
  // get position in the last column using binary search
  first = 0; last = m_pIndex->text_size;
  // invariant: the desired position is within first and last
  while(first<last) {
    middle = (first+last)/2;
    /* -------------------------------------------------------------
       get the char in position middle. As a byproduct, occ_sb
       and occ_b are initialized
       ------------------------------------------------------------- */
    get_info_sb(middle, occ_sb);   // init occ_sb[]  
    c_b = get_info_b(NULL_CHAR, middle, occ_b, WHAT_CHAR_IS); // init occ_b[] 
    assert(c_b < m_pIndex->alpha_size_sb);
    c_sb = m_pIndex->inv_map_sb[c_b];           
    assert(c_sb < m_pIndex->alpha_size);  // c_sb is the char in position middle 
  
    /* --------------------------------------------------------------
       count the # of occ of ch in [0,middle]
       -------------------------------------------------------------- */
    if(m_pIndex->bool_map_sb[ch]==0)
     n=occ_sb[ch];          // no occ of ch in this superbucket
    else {
      ch_b=0;                        // get remapped code for ch
      for(i=0;i<ch;i++)                  
	if(m_pIndex->bool_map_sb[i]) ch_b++;
      assert(ch_b<m_pIndex->alpha_size_sb);
      n = occ_sb[ch] + occ_b[ch_b];  // # of occ of ch in [0,middle]
    }
    /* --- update first or last ------------- */
    if(n>rank)
      last=middle;
    else if (n<rank)
      first=middle+1;
    else {                      // there are exactly "rank" c's in [0,middle]
      if(c_sb==ch) 
	first=last=middle;      // found!
      else 
	last=middle;            // restrict to [first,middle)
    }
  }
  // first is the desired row in the last column
  return first;
}


/*
   return the first character of a given row.
   This routine can be improved using binary search!
*/
unsigned char 
CFMIndex::get_firstcolumn_char(UINT32 row)
{
  int i;

  for(i=1;i<m_pIndex->alpha_size;i++)
    if(row<m_pIndex->bwt_occ[i])
      return i-1;
  return m_pIndex->alpha_size-1;
}


int 
CFMIndex::display(unsigned char *pattern, UINT32 length, UINT32 nums, UINT32 *numocc, 
			unsigned char **snippet_text, UINT32 **snippet_len) 
{

	multi_count *groups;
	int i, num_groups = 0, error;
	unsigned char *snippets;
	UINT32 *snip_len, j, h, len;

	*numocc = 0;
	*snippet_text = NULL;
	*snippet_len = 0;

	len = length + 2*nums;

	if(m_pIndex->smalltext) { //uses Boyer-Moore algorithm
		UINT32 *occ, to, numch, k;
		int error = fm_boyermoore(pattern, length, &occ, numocc);
		if(error<0 || *numocc <= 0) 
			return error;
		snip_len = (UINT32 *) malloc (sizeof (UINT32) * (*numocc));
		if (snip_len == NULL)
			{
			*numocc = 0;
			return FM_OUTMEM;
			}
		
		*snippet_text = (unsigned char *) new unsigned char [sizeof (unsigned char) * len *(*numocc)];
		snippets = *snippet_text;
		
		if (snippets == NULL)
			{
			*numocc = 0;
			return FM_OUTMEM;
			}

		for(k=0;k<*numocc;k++) {
			if(occ[k]<nums) 
				to = 0;
			else 
				to = occ[k]-nums;
				
			if(occ[k]+nums+length-1>m_pIndex->text_size-1) 
				numch =  m_pIndex->text_size - to;
			else 
				numch = occ[k]+nums+length-to;
				
			memcpy(snippets, m_pIndex->text+to, numch);
			snip_len[k] = numch;
			snippets += numch;
		}
		
		free(occ);
		*snippet_len = snip_len;
		return FM_OK;
	}
		

	/* count */
	num_groups = fm_multi_count (pattern, length, &groups);

	if (num_groups <= 0)
		return num_groups;

	for (i = 0; i < num_groups; i++)
		*numocc += groups[i].elements;

	snip_len = (UINT32 *) new UINT32 [sizeof (UINT32) * (*numocc)];
	if (snip_len == NULL)
		{
		*numocc = 0;
		return FM_OUTMEM;
		}

	snippets = (unsigned char *) new unsigned char [sizeof (unsigned char) * len *(*numocc)];
	if (snippets == NULL)
		{
		*numocc = 0;
		return FM_OUTMEM;
		}

	h = 0;
	for (i = 0; i < num_groups; i++)
	{
		
		for(j=0; j<groups[i].elements; j++) {
	
			error =	fm_snippet(groups[i].first_row + j, length, nums,snippets + h*len, &(snip_len[h]));
		
			if (error < 0)
				{
				delete snippets;
				delete snip_len;
				*numocc = 0;
				return error;
				}
			h++;	
		}
	}
	
	*snippet_text = snippets;
	*snippet_len = snip_len;
	free (groups);
	return FM_OK;

}

/*
   display the text sourronding a given pattern. Is is assumed that 
   the pattern starts at "row" (we do not know its position in 
   the input text), and we get the clen chars preceeding and 
   following it.
*/
   
int 
CFMIndex::fm_snippet(UINT32 row, UINT32 plen, UINT32 clen, unsigned char *dest, 
			   UINT32 *snippet_length) {

  UINT32 back, forw, i;
  unsigned char * temptext;
			   
  temptext = (unsigned char *)malloc(sizeof(unsigned char) * clen);
  if (temptext==NULL) return FM_OUTMEM;
			   	
  /* --- get clen chars preceding the current position --- */
  back = go_back(row, clen, temptext);
  assert(back <= clen);

  for(i=0; i<back; i++) // reverse temptext
	  	dest[i] = temptext[back-(i+1)];
  free(temptext);
  
  /* --- get plen+clen chars from the current position --- */
  forw = go_forw(row, clen+plen, dest+back);
  assert(forw <= clen+plen);
  if(forw<plen) return FM_GENERR;
  
  *snippet_length = back+forw;
  return FM_OK;
}




/* UNBUILD */

/* 
   read prologue of a .bwi file
   Output
     m_pIndex->numbucs, m_pIndex->buclist
*/
int 
CFMIndex::read_prologue(void)
{
  bucket_lev1 *sb;  
  UINT32 i, k, offset;
	  
  /* alloc superbuckets */
  m_pIndex->buclist_lev1 = (bucket_lev1 *) malloc(m_pIndex->num_bucs_lev1 * sizeof(bucket_lev1));
  if(m_pIndex->buclist_lev1==NULL) return FM_OUTMEM; 

  /* alloc aux array for each superbucket */
  for(i=0; i< m_pIndex->num_bucs_lev1; i++){
    sb = &(m_pIndex->buclist_lev1[i]);

    /* allocate space for array of occurrences */
    sb->occ = (UINT32 *) malloc((m_pIndex->alpha_size)* sizeof(UINT32));
	if(sb->occ==NULL) return FM_OUTMEM;

    /* allocate space for array of boolean char map */
    sb->bool_char_map = (unsigned char *)malloc((m_pIndex->alpha_size)*sizeof(unsigned char));
	if(sb->bool_char_map == NULL) return FM_OUTMEM;
  }  
	
  offset = m_pIndex->start_prologue_info_sb;
   
  for(i=0; i<m_pIndex->num_bucs_lev1; i++) {
    sb = &(m_pIndex->buclist_lev1[i]);
 	fm_init_bit_reader(m_pIndex->compress + offset);
  
    for(k=0; k<m_pIndex->alpha_size; k++)     /* boolean char_map */
      sb->bool_char_map[k] = fm_bit_read(1);

   	if(i>0)   {                         // read prefix-occ 
      	for(k=0;k<m_pIndex->alpha_size;k++) 
	  		memcpy(sb->occ, m_pIndex->compress + offset + m_pIndex->sb_bitmap_size, m_pIndex->alpha_size*sizeof(UINT32));
		offset += (m_pIndex->alpha_size * sizeof(UINT32) + m_pIndex->sb_bitmap_size);
	} else offset += m_pIndex->sb_bitmap_size;
	  
  }

  /* alloc array for the starting positions of the buckets */
  m_pIndex->start_lev2 =  (UINT32 *) malloc((m_pIndex->num_bucs_lev2)* sizeof(UINT32));
  if(m_pIndex->start_lev2 == NULL) return FM_OUTMEM;
 
  fm_init_bit_reader(m_pIndex->compress + m_pIndex->start_prologue_info_b);
  
  /* read the start positions of the buckets */
  for(i=0;i<m_pIndex->num_bucs_lev2;i++) 
    m_pIndex->start_lev2[i] = fm_bit_read(m_pIndex->var_byte_rappr);
  
  return FM_OK;
}  


/* 
   	retrieve the bwt by uncompressing the data in the input file
   	Output
     	m_pIndex->bwt  
*/ 
int 
CFMIndex::uncompress_data(void)
{
  UINT32 i;
  int error;

  m_pIndex->bwt = (unsigned char *) malloc(m_pIndex->text_size);
  if(m_pIndex->bwt == NULL)
    	return FM_OUTMEM;
 
  for(i=0; i < m_pIndex->num_bucs_lev1; i++){
    	error = uncompress_superbucket( i, m_pIndex->bwt+i*m_pIndex->bucket_size_lev1);
  		if(error < 0) return error;
  }
  
  return FM_OK;
} 


/* 
   expand the superbucket num. the uncompressed data is written
   in the array out[] which should be of the appropriate size  
   (i.e. m_pIndex->bucket_size_lev1 unless num is the last superbucket) 
*/
int 
CFMIndex::uncompress_superbucket(UINT32 numsb, unsigned char *out)
{
  bucket_lev1 sb;  
  unsigned char *dest, c;
  UINT32 sb_start, sb_end, start,  b2, temp_occ[ALPHASIZE];
  int i, k, error, len, is_odd, temp_len;
 
  assert(numsb<m_pIndex->num_bucs_lev1);
  sb = m_pIndex->buclist_lev1[numsb];    	/* current superbucket */
  sb_start = numsb*m_pIndex->bucket_size_lev1;/* starting position of superbucket */
  sb_end = MIN(sb_start+m_pIndex->bucket_size_lev1, m_pIndex->text_size);    
  b2 = sb_start/m_pIndex->bucket_size_lev2; /* initial level 2 bucket */
	
  m_pIndex->alpha_size_sb = 0;                /* build inverse char map for superbucket */
  for(k=0; k<m_pIndex->alpha_size; k++) 
    if(sb.bool_char_map[k]) 
      m_pIndex->inv_map_sb[m_pIndex->alpha_size_sb++] = k;

  for(start=sb_start; start < sb_end; start += m_pIndex->bucket_size_lev2, b2++) {
   
	len = MIN(m_pIndex->bucket_size_lev2, sb_end-start); // length of bucket
    dest = out + (start - sb_start);                       
	
	fm_init_bit_reader((m_pIndex->compress) + m_pIndex->start_lev2[b2]); // go to start of bucket
  	is_odd = 0;
  	if((b2%2 == 0) && (start != sb_start) && (b2 != m_pIndex->num_bucs_lev2-1))
		is_odd = 1; 
	
    if((start != sb_start) && (!is_odd)) // if not the first bucket and not odd skip occ
      for(k=0; k<m_pIndex->alpha_size_sb; k++) 
         error = fm_integer_decode(m_pIndex->int_dec_bits); /* non servono se non mtf2 */
	
	/* Compute bucket inv map */
  	m_pIndex->alpha_size_b = 0;
    for(i=0; i< m_pIndex->alpha_size_sb;i++)     
      if( fm_bit_read(1) ) 
		  m_pIndex->inv_map_b[m_pIndex->alpha_size_b++] = i;
	  
	assert(m_pIndex->alpha_size_sb >= m_pIndex->alpha_size_b);

    /* Applies the proper decompression routine */
    switch (m_pIndex->type_compression) 
	{
      case MULTIH: /* Bzip compression of mtf-ranks */
    	/* Initialize Mtf start */
		for (i = 0; i < m_pIndex->alpha_size_b; i++)
			m_pIndex->mtf[i] = i;
		if(is_odd)
			temp_len = 0;
		else 
			temp_len = len-1;
		
		get_b_multihuf(temp_len, temp_occ,  is_odd); /* temp_occ is not needed 
	  											       get_b_m modify m_pIndex->mtf_seq */
	 	break;

  	  default:  
			return FM_COMPNOTSUP;
      }

    /* remap the bucket according to the superbucket m_pIndex->inv_map_sb */
    for(i=0; i<len; i++) { 
      assert(m_pIndex->mtf_seq[i] < m_pIndex->alpha_size_sb);
	  c = m_pIndex->inv_map_sb[m_pIndex->mtf_seq[i]]; /* compute remapped char */
      assert(c < m_pIndex->alpha_size);          
      if (is_odd) dest[len-(i+1)] = c;	/* reverse this is an odd bucket */
	  else dest[i] = c;                      
    } 
	}
	return FM_OK;
}


/*
	compute the lf mapping (see paper)
   	Input
    	 m_pIndex->bwt, m_pIndex->bwt_occ[], m_pIndex->text_size, m_pIndex->alpha_size
   	Output
    	 m_pIndex->lf   
*/  
int 
CFMIndex::fm_compute_lf(void)
{
  UINT32 i, occ_tmp[ALPHASIZE];

  /* alloc memory */
  m_pIndex->lf = (UINT32 *) malloc(m_pIndex->text_size*sizeof(UINT32));
  if(m_pIndex->lf == NULL)
    return FM_OUTMEM;

  /* copy bwt_occ */
  for(i=0;i<ALPHASIZE;i++)
    occ_tmp[i] = m_pIndex->bwt_occ[i];

  /* now computes lf mapping */
  for(i=0;i<m_pIndex->text_size;i++)
    m_pIndex->lf[i] = occ_tmp[m_pIndex->bwt[i]]++;
  
  return FM_OK;
}


/* 
   compute the inverse bwt using the lf mapping       
   Input
     m_pIndex->bwt, m_pIndex->bwt_eof_pos, m_pIndex->text_size, m_pIndex->lf
   Output
     m_pIndex->text
*/    
int 
CFMIndex::fm_invert_bwt(void)
{
  UINT32 j;
  UINT32 i, real_text_size;
	
  if(m_pIndex->skip>1) real_text_size = m_pIndex->text_size-m_pIndex->num_marked_rows;
  else  
	  real_text_size = m_pIndex->text_size; 
	  
  /* alloc memory */
  m_pIndex->text = (unsigned char *) malloc(real_text_size*sizeof(unsigned char));
  if(m_pIndex->text == NULL)
    return FM_OUTMEM;

  for(j=0, i=real_text_size-1; i>0; i--) {
	if((m_pIndex->skip<=1) || (m_pIndex->bwt[j] != m_pIndex->specialchar)) 
	  m_pIndex->text[i] = m_pIndex->bwt[j];
	else i++;
	
    j = m_pIndex->lf[j];              // No account for EOF

    assert(j < m_pIndex->text_size);

    if(j < m_pIndex->bwt_eof_pos) j++; // EOF is not accounted in c[] and thus lf[]
                              // reflects the matrix without the first row.
                              // Since EOF is not coded, the lf[] is correct
                              // after bwt_eof_pos but it is -1 before.
                              // The ++ takes care of this situation.
	
  }
  
  //fprintf(stderr, "mark %lu realsiz %lu i %lu j %lu eof %lu size %lu\n",m_pIndex->num_marked_rows, real_text_size, i, m_pIndex->lf[j], m_pIndex->bwt_eof_pos, real_text_size); 
  /* i == 0 */
  m_pIndex->text[i] = m_pIndex->bwt[j];
  assert(j<m_pIndex->text_size);
  j = m_pIndex->lf[j];              // No account for EOF
  if(j<m_pIndex->bwt_eof_pos) j++;

  assert(j==m_pIndex->bwt_eof_pos);
  free(m_pIndex->lf); m_pIndex->lf = NULL;
  free(m_pIndex->bwt); m_pIndex->bwt = NULL;
  return FM_OK;
}

void 
CFMIndex::free_unbuild_mem(void) { 
	UINT32 i;
	bucket_lev1 *sb;
	
	free(m_pIndex->start_lev2);
	m_pIndex->start_lev2 = NULL;

	for(i=0; i< m_pIndex->num_bucs_lev1; i++) {
		sb = &(m_pIndex->buclist_lev1[i]);
		free(sb->occ);
		sb->occ = NULL;
		free(sb->bool_char_map);
		sb->bool_char_map = NULL;
	}
	free(m_pIndex->buclist_lev1);
	m_pIndex->buclist_lev1 = NULL;
}


int 
CFMIndex::fm_unbuild(unsigned char ** text, UINT32 *length) {
	
	int error;
	UINT32  i;
	if ((error = read_prologue()) < 0 ) {
			free_unbuild_mem();
			return error;
	}
	if ((error = uncompress_data()) < 0 ) {
			free_unbuild_mem();
			return error;
	}
	if ((error = fm_compute_lf()) < 0 ) {
			free_unbuild_mem();
			return error;
	}

	if ((error = fm_invert_bwt()) < 0 ) {
			free_unbuild_mem();
			return error;
	}

	UINT32 real_text_size;
    if(m_pIndex->skip>1) real_text_size = m_pIndex->text_size-m_pIndex->num_marked_rows;
  	else real_text_size = m_pIndex->text_size;
		
	/* remap text */	
	for(i=0; i<real_text_size; i++) {
		if(m_pIndex->text[i] == m_pIndex->specialchar) m_pIndex->text[i] = m_pIndex->subchar;
    	m_pIndex->text[i] = m_pIndex->inv_char_map[m_pIndex->text[i]];
	}
	*text = m_pIndex->text;
	m_pIndex->text = NULL;
	free_unbuild_mem(); /* libera memoria allocata */

	*length = real_text_size;
	return FM_OK;
}



/*
 * Sends to stderr an error message corresponding to error code e 
 */
const char* 
CFMIndex::error_index (int e)
{
if(e >= 0)
	return("No Errors");

switch (e)
	{
	case FM_GENERR:
		return "Appeared general error";
		break;
	case FM_OUTMEM:
		return "Malloc failed. Not unable to allocate memory";
		break;
	case FM_CONFERR:
		return "Some parameters of configuration are not correct";
		break;
	case FM_COMPNOTSUP:
		return "Compression Algorithm is not supported";
		break;
	case FM_DECERR:
		return "General error on decompression";
		break;
	case FM_COMPNOTCORR:
		return "Compressed file is not correct";
		break;
	case FM_SEARCHERR:
		return "Error during search";
		break;
	case FM_NOMARKEDCHAR:
		return "Impossible to make search/extract without marked chars";
		break;
	case FM_READERR:
		return "Can't open or read the file";
		break;
	case FM_NOTIMPL:
		return "Function not implemented";
		break;
	default:
		return "general error\n";
		break;
	}
}


/*
 * Compute # bits to represent u. Per calcolare log2(u) intero superiore
 * devo passargli (u-1) 
 */
int
CFMIndex::int_log2 (int u)
{
	/*
	 * codifica con if i casi piu frequenti un if costa meno di una
	 * moltiplicazione 
	 */

	if (u < 2)
		return 1;
	if (u < 4)
		return 2;
	if (u < 8)
		return 3;	// da 4 a 7
	if (u < 16)
		return 4;
	if (u < 32)
		return 5;
	if (u < 64)
		return 6;
	if (u < 128)
		return 7;
	if (u < 256)
		return 8;
	if (u < 512)
		return 9;
	if (u < 1024)
		return 10;
	if (u < 2048)
		return 11;
	if (u < 4096)
		return 12;
	if (u < 8192)
		return 13;
	if (u < 16384)
		return 14;
	if (u < 32768)
		return 15;
	if (u < 65536)
		return 16;
	if (u < 131072)
		return 17;
	if (u < 262144)
		return 18;
	if (u < 524288)
		return 19;
	if (u < 1048576)
		return 20;
	int i = 20;
	int r = 1048575;

	while (r < u)
	{
		r = 2 * r + 1;
		i = i + 1;
	}
	return i;
}

int
CFMIndex::int_pow2 (int u)		// compute 2^u
{
	int i, val;

	for (i = 0, val = 1; i < u; i++)
		val *= 2;

	assert (i == u);
	return val;

}



/* 
 * Main Build Functions 
 */


#define POINTPROLOGUE (19);

#if TESTINFO
typedef struct { // Solo per i test: memorizza spazio occupato
	UINT32 bucket_compr; 
	UINT32 bucket_occ;
	UINT32 bucket_alphasize;
	UINT32 bucket_pointer;		
	UINT32 sbucket_occ;       	
	UINT32 sbucket_alphasize; 	
	UINT32 prologue_size;
	UINT32 sbucket_bitmap;
	UINT32 bucket_bitmap;
	UINT32 marked_pos;
	UINT32 temp;
} measures;

measures Test;
#endif

int 
CFMIndex::fm_build_config (double freq, 
							UINT32 bsl1, 
							UINT32 bsl2, 
							UINT16 owner)	// 0==caller retains ownership of text to prcess, 1 == CFMIndex becomes owner and will delete this memory 
{
// instead of returning errors on paparameters force parameters to be reasonable values
// freq must be between 0.0 and 1.0 inclusive
if(freq < 0.0)
	freq = 0.0;
else
	if(freq > 1.0)
		freq = 1.0;

// clamp bsl1,bsl2 to within their min/max limits
if(bsl2 < cMinL2BlockSize)
	bsl2 = cMinL2BlockSize;
else
	if(bsl2 > cMaxL2BlockSize)
		bsl2 = cMaxL2BlockSize;

if(bsl1 < cMinL1BlockSize)
	bsl1 = cMinL1BlockSize;
else
	if(bsl1 > cMaxL1BlockSize)
		bsl1 = cMaxL1BlockSize;

// bsl1 must always be at least 8*bsl2
if(bsl1 < (bsl2*8))
	bsl1 = bsl2*8;
else
	{
	// force bsl1 to be an exact multiple of bsl2
	if(bsl1 % bsl2)
		bsl1 = bsl2 * ((bsl1 + bsl2 - 1)/ bsl2);
	}

m_pIndex->owner = owner;
m_pIndex->type_compression = MULTIH;
m_pIndex->bucket_size_lev1 = bsl1  << 10;
m_pIndex->bucket_size_lev2 = bsl2  << 10;					
if (freq >= 0.5) { m_pIndex->skip = 1; return FM_OK;}
if(freq == 0) { m_pIndex->skip = 0; return FM_OK;}
m_pIndex->skip = (UINT16) (1.0/freq);	/* 1/Mark_freq 2% text size */
return FM_OK;
}


int 
CFMIndex::parse_options(char *optionz) {
	
	int i, numtoken = 1;
  	
	/* default */
	UINT32 bsl1 = cDfltL1BlockSize;
	UINT32 bsl2 = cDfltL2BlockSize;
	UINT16 owner = 1;
	double freq = cDfltMarkerFreq;	

	if (optionz == NULL) 
			return fm_build_config (freq, bsl1, bsl2, owner);
	
	char *options = (char *)malloc(sizeof(char)*(strlen(optionz)+1));
	if (options == NULL) 
			return FM_OUTMEM;
	memcpy(options,optionz,sizeof(char)*(strlen(optionz)+1));
	
	i = 0;
	while (options[i] != '\0' )
		if(options[i++] == ' ')	 numtoken++;

	int j = 0;
	char **pArray;

	pArray = (char **)malloc(numtoken * sizeof(char *));

	pArray[0] = options;
	for(i = 0; i<numtoken-1; i++) {
		while(options[j] != ' ') j++;	
		options[j++] = '\0';
		pArray[i+1] = options+j;	
	}

	i = 0;
	while(i<numtoken) {
		if(!strcmp(pArray[i], "-B")) {
			if (i+1<numtoken ) {
				bsl1 = atoi(pArray[++i]); 	
				i++;
				continue;
				} 
			else 
				{
				free(pArray);
				return FM_CONFERR;
				}
			}
		if(!strcmp(pArray[i], "-b")) {
			if (i+1<numtoken ) {
				bsl2 = atoi(pArray[++i]); 
				i++;
				continue;
			} else
				{
				free(pArray);
				return FM_CONFERR; 
				}
			}
		if(!strcmp(pArray[i], "-a")) {
			if (i+1<numtoken ) {
				owner = atoi(pArray[++i]); 
				i++;
				continue;
			} else 
			{
			free(pArray);
			return FM_CONFERR;
			}
			}
		if(!strcmp(pArray[i], "-F")) {
			if (i+1<numtoken ) {
				freq = atof(pArray[++i]); 
				i++;
				continue;
			} else 
			{
			free(pArray);
			return FM_CONFERR;
			}
			}
			free(pArray);
      		return FM_CONFERR;
	}
	free(options);
	free(pArray);

	return fm_build_config(freq, bsl1, bsl2, owner);
}
/* 
	Creates index from text[0..length-1]. Note that the index is an 
	opaque data type. Any build option must be passed in string 
	build_options, whose syntax depends on the index. The index must 
	always work with some default parameters if build_options is NULL. 
	The returned index is ready to be queried. 
*/
int 
CFMIndex::build_index(unsigned char *text, UINT32 length, char *build_options) {
	
	int error;

	memset(m_pIndex,0,sizeof(fm_index));

	error = parse_options(build_options);
	if (error < 0) return error;
		
	error = fm_build(text, length);
	if (error < 0) {
		return error;
	}
	
	return FM_OK;
}

/* Build */
int 
CFMIndex::fm_build(unsigned char *text, UINT32 length) {

	int error;	
	m_pIndex->compress = NULL;
	m_pIndex->lf = NULL;
	m_pIndex->loc_occ = NULL;
	m_pIndex->bwt = NULL;
	m_pIndex->buclist_lev1 = NULL;
	m_pIndex->start_lev2 = NULL;

	m_pIndex->text = text;
	m_pIndex->text_size = length;
	m_pIndex->compress_size = 0;


	#if TESTINFO
	Test.bucket_compr = 0; 
	Test.bucket_occ = 0;
	Test.bucket_alphasize = 0;
	Test.sbucket_alphasize = 0; 	
	#endif
	if(m_pIndex->text_size < SMALLSMALLFILESIZE) m_pIndex->skip = 0;
	m_pIndex->compress_owner = 2;
	if((m_pIndex->owner != 0) && (m_pIndex->skip <= 1)) {
					
		unsigned char * texts = (unsigned char *)malloc(sizeof(unsigned char)*(m_pIndex->text_size));
		if ( texts == NULL) return FM_OUTMEM;
		memcpy(texts, m_pIndex->text, (m_pIndex->text_size)*sizeof(unsigned char));
		if (m_pIndex->owner == 1)
			free(m_pIndex->text);
		m_pIndex->text = texts;
		}
	
	if(m_pIndex->text_size < SMALLSMALLFILESIZE) {
		
		m_pIndex->smalltext = 1;
		m_pIndex->compress_size = 0;
		m_pIndex->compress = (unsigned char *)malloc((4+m_pIndex->text_size)*sizeof(unsigned char));
		if (m_pIndex->compress == NULL) 
			return FM_OUTMEM;
		m_pIndex->compress_size = 0;
		fm_init_bit_writer(m_pIndex->compress, &m_pIndex->compress_size);
		fm_uint_write(m_pIndex->text_size);
		memcpy(m_pIndex->compress+4,m_pIndex->text, sizeof(unsigned char)*m_pIndex->text_size);
		m_pIndex->compress_size += m_pIndex->text_size;	
		
		return FM_OK;
	}	


	if(m_pIndex->text_size < SMALLFILESIZE) { // writes plain text

		m_pIndex->smalltext = 2;
		error = build_sa();
		if (error) return error;
		error = build_bwt();
		if (error) return error;
		free(m_pIndex->lf);
		m_pIndex->lf = NULL;
		error = fm_bwt_compress();
		if(m_pIndex->owner == 1) 
			{
			free(m_pIndex->text); 
			m_pIndex->text = NULL;
			}
		return(error);
	}
	
	m_pIndex->smalltext = 0;
	
	error = select_subchar();
	if (error < 0) return errore( error);
	if(TESTINFO) fprintf(stderr, "Select char done\n");
		
	m_pIndex->log2textsize = int_log2 (m_pIndex->text_size - 1);
	m_pIndex->int_dec_bits = int_log2 (int_log2(m_pIndex->bucket_size_lev1 - m_pIndex->bucket_size_lev2));
	m_pIndex->mtf_seq = (unsigned char *) malloc (m_pIndex->bucket_size_lev2 * sizeof (unsigned char));
	
	if (m_pIndex->mtf_seq == NULL)
		return FM_OUTMEM;	
	
	/* Build suffix array */  
   	 error = build_sa();		
    	if ( error < 0 ) 
			return errore( error);
	if(TESTINFO) fprintf(stderr,"Build suffix array done\n");
	
	count_occ();  /* conta occorrenze */

	/* make bwt */
	error = build_bwt(); 
	if (error < 0 )
			return errore( error);
	if(TESTINFO) fprintf(stderr,"Build BWT done\n");
	
	if(m_pIndex->skip>1) { free(m_pIndex->text); m_pIndex->text = NULL;}
	else if(m_pIndex->owner == 1) {
		free(m_pIndex->text); 
		m_pIndex->text = NULL;
		}
	
	UINT32 i;
	/* Remap bwt */
	for(i=0; i<m_pIndex->text_size; i++) 
   		m_pIndex->bwt[i] = m_pIndex->char_map[m_pIndex->bwt[i]];
	
	/* Mark position */
	error = compute_locations(); 
	if (error < 0) 
			return errore( error);
	if(TESTINFO) fprintf(stderr,"Compute locations done\n");
	
	error = compute_info_superbuckets();
	if (error < 0 )
			return errore( error);	
			
	/* Compute various infos for each bucket */ 
	error = compute_info_buckets();
	if (error < 0 ) 
			return errore( error);
	
 	/* alloc memory for compress */
	if(m_pIndex->skip != 1) {
			free(m_pIndex->lf);
			m_pIndex->lf = NULL;
		}

	UINT32 stima_len_compress = (UINT32)(100 + m_pIndex->text_size*1.5);
	
	stima_len_compress += m_pIndex->num_marked_rows*sizeof(UINT32); 
	
	m_pIndex->compress = (unsigned char *)malloc(stima_len_compress*sizeof(unsigned char));
	if (m_pIndex->compress == NULL) 
			return errore( FM_OUTMEM);
	
	m_pIndex->compress_size = 0;

	/* Write Prologue on mem */
	write_prologue();	
	if(TESTINFO) fprintf(stderr,"Write prologue done\n");
			
	/* Compress each bucket in all superbuckets */	
	for(i=0; i<m_pIndex->num_bucs_lev1; i++) { /* comprimi ogni bucket */
             error = compress_superbucket( i);
			 if (error < 0 ) 
					return errore( error);	
		 }
	if(TESTINFO) fprintf(stderr,"Compress buckets done\n");
		
    m_pIndex->start_positions = m_pIndex->compress_size;
	
	/* write the starting position of buckets */
  	write_susp_infos(); 
	if(TESTINFO) fprintf(stderr,"Write susp info done\n");

	dealloc_bucketinfo();
	 
	/* write marked chars */
	if (m_pIndex->skip>0) 
			write_locations();
   	if(TESTINFO) fprintf(stderr,"Write locations done\n");
		
	if(m_pIndex->skip == 1) {
			free(m_pIndex->lf);
			m_pIndex->lf = NULL;
	}
	
	if(m_pIndex->skip > 1) {
			free(m_pIndex->loc_occ);
			m_pIndex->loc_occ = NULL;
	}
	
	m_pIndex->compress = (unsigned char *)realloc(m_pIndex->compress, m_pIndex->compress_size);
	dealloc();

	#if TESTINFO
	fprintf(stderr, "<indexSize>%lu</indexSize>/n<textSize>%lu</textSize>\n", m_pIndex->compress_size, m_pIndex->text_size);
	fprintf(stderr, "<textAlphasize>%d</textAlphasize>\n", m_pIndex->alpha_size);
	fprintf(stderr, "<prologueSize>%lu</prologueSize>\n", Test.prologue_size);
	fprintf(stderr, "<numBuckets2>%lu</numBuckets2>\n", m_pIndex->num_bucs_lev2);
	fprintf(stderr, "<numBuckets1>%lu</numBuckets1>\n", m_pIndex->num_bucs_lev1);
	fprintf(stderr, "Bucket pointers size %lu Kb bits x ponter %u\n", Test.bucket_pointer/1024, m_pIndex->var_byte_rappr);	
    fprintf(stderr, "Buckets compressed size %lu Kb\n", Test.bucket_compr/1024); 
	fprintf(stderr, "Buckets prefix chars occ size %lu Kb\n", Test.bucket_occ/1024);
	fprintf(stderr, "Averange buckets alpha size %lu\n", Test.bucket_alphasize/m_pIndex->num_bucs_lev2);
	fprintf(stderr, "Bucket bitmaps size %lu\n\n", Test.bucket_bitmap/1024);

	fprintf(stderr, "Superbuckets size %lu # Superbuckets %lu\n", m_pIndex->bucket_size_lev1, m_pIndex->num_bucs_lev1);
	fprintf(stderr, "Superbuckets prefix chars occ size %lu Kb\n", Test.sbucket_occ/1024);
	fprintf(stderr, "Average superbuckets alpha size %lu\n", Test.sbucket_alphasize/m_pIndex->num_bucs_lev1);
	fprintf(stderr, "Superbucket bitmaps size %lu\n\n", Test.sbucket_bitmap/1024);
	fprintf(stderr, "# marked positions %lu\nMarked positions size %lu Kb\n\n", m_pIndex->num_marked_rows, Test.marked_pos/1024);
	#endif

	return FM_OK;	
}
	 

/* 
	Cerca un carattere non presente e lo inserisce ogni m_pIndex->skip posizioni nel testo originale
*/
int
CFMIndex::select_subchar(void) {
	UINT32 i, newtextsize, pos;
	UINT16 mappa[ALPHASIZE];
	m_pIndex->subchar = 0; // inutile in questa versione
	if (m_pIndex->text_size <= m_pIndex->skip) m_pIndex->skip=1;
	m_pIndex->oldtext = NULL;
	
	if (m_pIndex->skip == 1) {
	    m_pIndex->num_marked_rows  = m_pIndex->text_size;
		return FM_OK;
		}
	if (m_pIndex->skip == 0)  {
		m_pIndex->num_marked_rows = 0;
		return FM_OK;
		}

	/* Controlla se c'e' un carattere non presente nel testo */	
	for(i=0;i<ALPHASIZE;i++) 
		mappa[i] = 0;
	
	for(i=0;i<m_pIndex->text_size;i++) 
		mappa[m_pIndex->text[i]] = 1;

	for(i=0;i<ALPHASIZE;i++) 
		if (!mappa[i]) {
			m_pIndex->specialchar = (unsigned char)i; /* Nuovo carattere aggiunto */
			break;
	}
	
	if(i==ALPHASIZE) {
			fprintf(stderr,"256 chars I cannot insert a new chars. I mark all positions\n");
			m_pIndex->num_marked_rows  = m_pIndex->text_size; m_pIndex->skip = 1;
			return FM_OK;
	}
	
	m_pIndex->num_marked_rows = m_pIndex->text_size/m_pIndex->skip;
	newtextsize = m_pIndex->text_size + m_pIndex->num_marked_rows;
		
	unsigned char * text = (unsigned char *)malloc(newtextsize);
	if (!text) return FM_OUTMEM;
	
	for(i=0,pos = 0;i<m_pIndex->num_marked_rows; i++, pos += m_pIndex->skip) {
		if (pos) 
			text[pos++] = m_pIndex->specialchar;
		memcpy(text+pos, m_pIndex->text+pos-i, m_pIndex->skip);		
		}

	UINT32 offset = m_pIndex->text_size - (pos - i)-1;
	if(offset){
		if (pos) text[pos++] = m_pIndex->specialchar;
		memcpy(text+pos, m_pIndex->text+pos-i, offset);
	}
	
	if (m_pIndex->owner == 1)  
		{ 
		free(m_pIndex->text); 
		m_pIndex->text = NULL;
		m_pIndex->oldtext = NULL;
		} 
	else 
		m_pIndex->oldtext = m_pIndex->text;

	m_pIndex->text = text;
	m_pIndex->text_size = newtextsize;
	return FM_OK;
	
}


/* 
   build suffix array for the text in m_pIndex->text[] not remapped
 */

int 
CFMIndex::build_sa(void) { 
	
  m_pIndex->lf = (UINT32 *)malloc(m_pIndex->text_size * sizeof(UINT32));
  if (m_pIndex->lf == NULL) 
	  	return FM_OUTMEM;

  /* compute Suffix Array with library */  
  SAIS.sais(m_pIndex->text, (INT32 *)m_pIndex->lf, m_pIndex->text_size);
  return FM_OK;

}

/* 	
	This procedure count occurences in the text.
	Init m_pIndex->alpha_size, m_pIndex->bool_char_map, m_pIndex->pfx_char_occ
*/

void 
CFMIndex::count_occ(void) {
	
	/* count occurences */
	UINT32 i;
	unsigned char curchar;

	m_pIndex->alpha_size = 0;
	for ( i = 0; i < ALPHASIZE; i++) { // init occ array
		 m_pIndex->bool_char_map[i] = 0;
		 m_pIndex->pfx_char_occ[i] = 0;		
	}
   
	/* calcolo numero occorrenze caratteri, alpha_size e bool_char_map */
	for (i = 0; i < m_pIndex->text_size; i++) {
		curchar = m_pIndex->text[i];
		if ( m_pIndex->bool_char_map[curchar] == 0 ) {
				m_pIndex->alpha_size++;
				m_pIndex->bool_char_map[curchar] = 1;
		}
		m_pIndex->pfx_char_occ[curchar]++;
	}
	
	/* remap dell'alfabeto */
	UINT16 cfree = 0; /* primo posto libero in m_pIndex->char_map */
	for (i =0; i < ALPHASIZE; i++) {
		if ( m_pIndex->bool_char_map[i] == 1 ) {
			m_pIndex->char_map[i] = (unsigned char)cfree;
			m_pIndex->pfx_char_occ[cfree] = m_pIndex->pfx_char_occ[i];
			cfree++;
		}
	}
	assert(cfree == m_pIndex->alpha_size);
	
	/* Compute prefix sum of char occ */
	UINT32 temp, sum = 0;
    for(i=0; i<m_pIndex->alpha_size; i++) {
    	temp = m_pIndex->pfx_char_occ[i];
    	m_pIndex->pfx_char_occ[i] = sum; 
    	sum += temp;
  	}	

}

/*
   compute BWT of input text.
   Remap alphabet when bwt is computed.
   La rimappatura potrebbe essere fatta in contemporanea con 
   la creazione del bwt ( con m_pIndex->char_map[m_pIndex->text[m_pIndex->lf[i]-1]] )
   ma e' piu lenta di + di 1 secondo su 9.
   Input:
     m_pIndex->text, m_pIndex->text_size, m_pIndex->lf
   Output
     m_pIndex->bwt, m_pIndex->bwt_eof_pos
*/ 

int 
CFMIndex::build_bwt(void) {
  
  UINT32 i;
  
  /* alloc memory for bwt */
  m_pIndex->bwt = (unsigned char *)malloc(m_pIndex->text_size*sizeof(unsigned char));
  if(m_pIndex->bwt==NULL)
  			return FM_OUTMEM;
  
  m_pIndex->bwt[0] = m_pIndex->text[m_pIndex->text_size-1];  //L[0] = Text[n-1]

  /* finche non incontro l'EOF il bwt e' piu' avanti di un elemento 
  	 dopo gli indici di sa e bwt ritornano uguali. 
  */
  
  UINT32 *sa = m_pIndex->lf;		 	 // punta al primo elemento 
  unsigned char *bwt = &(m_pIndex->bwt[1]);  	 // punta al secondo 
  for(i=0; i<m_pIndex->text_size; i++) {  // non leggibile ma piu' performante 
    if(*sa !=0 ){		     	 // al posto di *bwt c'era m_pIndex->bwt[j++] con j=1
      *bwt = m_pIndex->text[(*sa)-1]; 	 // al posto di *safix m_pIndex->lf[i]
	  bwt++; 
	} else
      m_pIndex->bwt_eof_pos = i; // EOF is not written but its position remembered !	
	sa++;
	}

   return FM_OK;	
}

/* 
  	Si ricorda delle posizioni. Funzionante per multilocate 
*/
int 
CFMIndex::compute_locations(void) { 

	UINT32 i,j; /* numero posizioni marcate */
	unsigned char spchar = m_pIndex->char_map[m_pIndex->specialchar];
	UINT32 firstrow = m_pIndex->pfx_char_occ[spchar];

  	if( (m_pIndex->skip==0)|| (m_pIndex->skip == 1)) 
    	return FM_OK;
	
	/* alloc m_pIndex->loc_occ */
  	m_pIndex->loc_occ = (UINT32 *)malloc(sizeof(UINT32) * m_pIndex->num_marked_rows);
	if (m_pIndex->loc_occ == NULL) return FM_OUTMEM;
		
	for(i=firstrow,j=firstrow; i<(m_pIndex->num_marked_rows+firstrow); i++,j++) {
					//if(m_pIndex->lf[i] == 0) { PLEASE FIX ME
					//	j++;
					//	continue;
					//	}
					m_pIndex->loc_occ[i-firstrow] = m_pIndex->lf[j];	
				}

  return FM_OK;
   
}

/* 
   Init m_pIndex->buclist_lev1 (list of superbuckets)
   For each superbuckets init the fields:  
   	occ[], alpha_size, bool_char_map[]
*/             

int 
CFMIndex::compute_info_superbuckets(void)
{
  UINT32 b, temp, occ, k, i;
  bucket_lev1 * sb;
  
  /* compute number of superbuckets  */
  m_pIndex->num_bucs_lev1 = (m_pIndex->text_size + m_pIndex->bucket_size_lev1 - 1) / m_pIndex->bucket_size_lev1;

  /* alloc superbuckets */
  m_pIndex->buclist_lev1 = (bucket_lev1 *) malloc(m_pIndex->num_bucs_lev1 * sizeof(bucket_lev1));
  if(m_pIndex->buclist_lev1==NULL) return FM_OUTMEM; 
                
  /* alloc aux array for each superbucket */
  for(i=0; i< m_pIndex->num_bucs_lev1; i++) {
    sb = &(m_pIndex->buclist_lev1[i]);   // sb points to current superbucket

    /* Allocate space for data structures */
    sb->occ = (UINT32 *)malloc((m_pIndex->alpha_size)* sizeof(UINT32));
    if(sb->occ == NULL) 
			return FM_OUTMEM;
	
    sb->bool_char_map = (unsigned char *)malloc((m_pIndex->alpha_size)*sizeof(UINT16));
    if(sb->bool_char_map == NULL) 
			return FM_OUTMEM;

    /* Initialize data structures */
    sb->alpha_size = 0;    
    for(k=0; k<m_pIndex->alpha_size; k++){
      sb->bool_char_map[k] = 0;
	  sb->occ[k] = 0;}
  	
  }  
  
  /* scan bwt and gather information */
  
  UINT32 currentBuck = 0;  // indice superbuckets
  UINT32 dim = m_pIndex->bucket_size_lev1; // per non dividere
  
  sb =  &m_pIndex->buclist_lev1[0];
  for(i=0; i<m_pIndex->text_size; i++) {
	
	  if ( i == dim ) {  // NON DIVIDERE sono in ordine 
		         currentBuck++;
		         dim += m_pIndex->bucket_size_lev1;
	  			 sb = &(m_pIndex->buclist_lev1[currentBuck]);
		         }
		          
    k = m_pIndex->bwt[i];                 // current char
    		 
    if(sb->bool_char_map[k]==0) { // build char_map of current sb
      sb->bool_char_map[k] = 1; 
	  sb->occ[k] = 0;
      sb->alpha_size++;           // compute alphabet size
    }
    sb->occ[k]++;
  }

  /* compute occ in previous buckets */
  for(k=0; k<m_pIndex->alpha_size;k++) {
    occ = 0;
    for(b=0; b<m_pIndex->num_bucs_lev1; b++) {   //prefix sum on OCC
      temp = m_pIndex->buclist_lev1[b].occ[k];
      m_pIndex->buclist_lev1[b].occ[k]=occ;
	  occ += temp;
	  
	  #if TESTINFO
  	  Test.sbucket_alphasize += m_pIndex->buclist_lev1[b].alpha_size;
	  Test.bucket_bitmap += (m_pIndex->buclist_lev1[b].alpha_size*(m_pIndex->bucket_size_lev1/m_pIndex->bucket_size_lev2));
  	  #endif
    }
  }

  return FM_OK;
}

/* 
   Init m_pIndex->start_lev2 which is the starting position of each bucket
   in the compressed file
*/             
int 
CFMIndex::compute_info_buckets(void)
{
  /* compute number of buckets */
  assert((m_pIndex->bucket_size_lev1 % m_pIndex->bucket_size_lev2) == 0);
  m_pIndex->num_bucs_lev2 = (m_pIndex->text_size + m_pIndex->bucket_size_lev2 - 1) / m_pIndex->bucket_size_lev2;

  /* alloc array for buckets starting positions */
  m_pIndex->start_lev2 =  (UINT32 *)malloc((m_pIndex->num_bucs_lev2)* sizeof(UINT32));
  if(m_pIndex->start_lev2==NULL) 
	  	return FM_OUTMEM;
  return FM_OK;	  
}


/* ---------------------------------------------------------
   The current format of the prologue is the following:
	8 bits	type of compression (4=Multi Table Huff 5=gamma 6=mtf2)
   32 bits 	size of input file (text)
   32 bits  position of eof in m_pIndex->bwt
   16 bits  size of a super bucket divided by 1024
   16 bits  size of a bucket divided (divides the previous)
	8 bits  size-1 of the compacted alphabet of the text
	8 bits 	m_pIndex->specialchar remapped eith m_pIndex->char_map
   32 bits	m_pIndex->skip expected # of chars between 2 marked positions 
   32 bits  occ-explicit list
   32 bits  m_pIndex->start_prologue_info_sb
    8 bits  m_pIndex->char_map[m_pIndex->subchar]
  256 bits  m_pIndex->bool_char_map: bit i is 1 iff char i occurs in the text

   	  		For each c that occurs in ter text, writes the m_pIndex->pfx_char_occ[c]
      		on m_pIndex->log2textsize bits
	x bits  bit_flush() in order to be aligned to byte 
		
  			For each Superbucket i 
			- m_pIndex->alpha_size bits to store the Sb bitmap
			- for each char c that occurs in te text 
			  stores the sb.occ[c] (prefix sum of character occurrences) 
			  with m_pIndex->log2textsize (except for the first sbucket)
	x bits  bit_flush() in order to be aligned to byte 
	
			Stores starting position of each bucket in the compressed file 
			using (m_pIndex->log2textsize+7)/8 * 8 bits. This is byte aligned.
			
*/
  
void 
CFMIndex::write_prologue(void) {

  bucket_lev1 sb;
  UINT32 i,k;

  /* write file and bucket size */
  fm_init_bit_writer(m_pIndex->compress, &m_pIndex->compress_size);
  fm_uint_write(m_pIndex->text_size);
  fm_bit_write24(8, m_pIndex->type_compression);
  fm_uint_write(m_pIndex->bwt_eof_pos);

  assert(m_pIndex->bucket_size_lev1 >> 10 < 65536);
  assert((m_pIndex->bucket_size_lev1 & 0x3ff) == 0);
  fm_bit_write24(16, m_pIndex->bucket_size_lev1 >> 10);

  assert(m_pIndex->bucket_size_lev2 >= 256);
  assert(m_pIndex->bucket_size_lev2 <= 4096);
  fm_bit_write24(16, m_pIndex->bucket_size_lev2);

  /* alphabet information */
  assert(m_pIndex->alpha_size>0 && m_pIndex->alpha_size<=ALPHASIZE);
  fm_bit_write24(8, m_pIndex->alpha_size-1);   
  
  /* write starting byte of occ-list */
  fm_bit_write24(8, m_pIndex->char_map[m_pIndex->specialchar]); /* carattere speciale */
  fm_bit_write(32, m_pIndex->skip);
  fm_uint_write(0); /* spazio libero per occ-explicit list 20esimo byte--*/
  fm_uint_write(0); /* spazio libero per mettere inizio sb info 25esimo byte */
  fm_bit_write24(8, m_pIndex->char_map[m_pIndex->subchar]); /* carattere sostituito */
  
  /* boolean alphabet char map */
  for(i=0; i<ALPHASIZE; i++)
   if (m_pIndex->bool_char_map[i]) {fm_bit_write24(1,1);}
   else {fm_bit_write24(1,0); }

  /* write prefix sum of char occ  Da valutare vantaggi in compressione vs tempo */
  for(i=0; i < m_pIndex->alpha_size; i++) {
   	fm_bit_write(m_pIndex->log2textsize, m_pIndex->pfx_char_occ[i]);   
	}
  fm_bit_flush(); 
 
  m_pIndex->start_prologue_info_sb = m_pIndex->compress_size;

  /* Process superbuckets */
  for(i=0;i<m_pIndex->num_bucs_lev1;i++) {
	sb = m_pIndex->buclist_lev1[i];

	for(k=0; k<m_pIndex->alpha_size; k++)  // boolean char_map
     	if(sb.bool_char_map[k]){ fm_bit_write24(1,1);}
     	else { fm_bit_write24(1,0);}
	fm_bit_flush();
  	if(i>0) // write prefix-occ
		{ 
		unsigned char *dest = m_pIndex->compress+m_pIndex->compress_size;
		dest = (unsigned char *)memcpy(dest, sb.occ, sizeof(UINT32) * m_pIndex->alpha_size);
	  	m_pIndex->compress_size += sizeof(UINT32) * m_pIndex->alpha_size;
	  	fm_init_bit_writer(m_pIndex->compress+m_pIndex->compress_size, &m_pIndex->compress_size);
	 }
  }

  /* leave space for storing the start positions of buckets */
  m_pIndex->var_byte_rappr = ((m_pIndex->log2textsize+7)/8)*8;   // it's byte-aligned
 
  for(i=0;i<m_pIndex->num_bucs_lev2;i++)
  	fm_bit_write(m_pIndex->var_byte_rappr,0);
  
  #if TESTINFO
  Test.bucket_pointer = (m_pIndex->var_byte_rappr)*m_pIndex->num_bucs_lev2/8;
  Test.sbucket_occ = (m_pIndex->num_bucs_lev1-1)*(m_pIndex->log2textsize*m_pIndex->alpha_size)/8;
  Test.sbucket_bitmap = m_pIndex->alpha_size*m_pIndex->num_bucs_lev1/8;
  Test.bucket_bitmap = Test.bucket_bitmap/8;
  Test.prologue_size = m_pIndex->compress_size;
  #endif

}


/* 
   compress a superbucket

   NOTE the content of m_pIndex->bwt is changed (remapped) !!!! 

   The number of occurrences of each char is represented using
   integer_encode.
*/
int 
CFMIndex::compress_superbucket(UINT32 num)
{
  
  bucket_lev1 sb;  
  unsigned char *in, c, char_map[ALPHASIZE];
  UINT32 k, temp, sb_start, sb_end, start,len, bocc[ALPHASIZE], b2,i,j;
  int is_odd;

  assert(num<m_pIndex->num_bucs_lev1);

  sb = m_pIndex->buclist_lev1[num];             // current superbucket
  sb_start = num * m_pIndex->bucket_size_lev1; // starting position of superbucket 
  sb_end = MIN(sb_start+m_pIndex->bucket_size_lev1, m_pIndex->text_size);    
  b2 = sb_start/m_pIndex->bucket_size_lev2;    // initial level 2 bucket il numero del bucket in assoluto

  temp=0;                               // build char map for superbucket 
  for(k=0;k<m_pIndex->alpha_size;k++) 
    if(sb.bool_char_map[k]) 
      char_map[k] = (unsigned char)temp++;
	
  assert(temp==sb.alpha_size);

  for(i=0; i<sb.alpha_size;i++)              // init occ
    bocc[i]=0;

  for(start=sb_start; start<sb_end; start+=m_pIndex->bucket_size_lev2, b2++) {

	m_pIndex->start_lev2[b2] = m_pIndex->compress_size; // start of bucket in compr file

    len = MIN(m_pIndex->bucket_size_lev2, sb_end-start);    // length of bucket
    in = m_pIndex->bwt + start;  // start of bucket
	
	is_odd = 0;
	if((b2%2 ==0) && (start != sb_start) && (b2 != m_pIndex->num_bucs_lev2-1))
		is_odd = 1; // decide quali bucket rovesciare e non mettere occ
	
    if(start != sb_start)       // if not the first bucket write occ to file
			if(!is_odd) 
				for(k=0; k<sb.alpha_size; k++) {
				#if TESTINFO	
				Test.temp = m_pIndex->compress_size;
			 	#endif
			
				fm_integer_encode(bocc[k], m_pIndex->int_dec_bits); 
			
				#if TESTINFO	
				Test.bucket_occ+=(m_pIndex->compress_size-Test.temp);
			 	#endif
				}
	
    for(j=0; j<len; j++) {          // update occ[] and remap
      assert(in[j]< m_pIndex->alpha_size); 
      c = char_map[in[j]];          // compute remapped char
      assert(c < sb.alpha_size);          
      in[j]=c;                      // remap
      bocc[c]++;                    // bocc is used in the next bucket
    }    
	
	// ROVESCIA i bucket dispari
	if(is_odd) {
			unsigned char *intemp;
			intemp = (unsigned char *)malloc(len * sizeof(unsigned char));
			for(j=0; j<len; j++) 
				{
				intemp[len-j-1] = in[j];
				}
			for(j=0; j<len; j++) 
				in[j] = intemp[j];
			free(intemp);
			}
	#if TESTINFO	
	Test.temp = m_pIndex->compress_size;
	#endif

	int error = compress_bucket( in, len, sb.alpha_size);
	
	#if TESTINFO	
	Test.bucket_compr+=(m_pIndex->compress_size-Test.temp-(sb.alpha_size/8));
	#endif
	if (error < 0) return error;
}

  return FM_OK;
}

/*
   write the starting position (in the output file) of each one 
   of the m_pIndex->num_bucs_lev2 buckets. These values are written at the
   end of the prologue just before the beginning of the first bucket.
   It writes also the starting position of the occurrence list
*/ 
void 
CFMIndex::write_susp_infos(void) {
  	 
  	UINT32 i, offset, unuseful = 0;
	unsigned char *write;

  	/* write starting position of points to positions  marked and sb occ list */
  	// warning! the constant POINTPROLOGUE depends on the structure of the prologue!!!
	write = m_pIndex->compress + POINTPROLOGUE;
  	fm_init_bit_writer(write, &unuseful); // vado a scrivere dove avevo lasciato spazio
  	fm_uint_write(m_pIndex->start_positions);
	fm_uint_write(m_pIndex->start_prologue_info_sb); // inizio sbuckets info
		
  	// Warning: the offset heavily depends on the structure of prologue.
  	//          The value of start_level2[0] has been initialized in
  	//          the procedure compress_superbucket()

  	offset = m_pIndex->start_lev2[0] - m_pIndex->num_bucs_lev2*(m_pIndex->var_byte_rappr/8);
	write = m_pIndex->compress + offset;
 	fm_init_bit_writer(write, &unuseful);
  	for(i=0;i<m_pIndex->num_bucs_lev2;i++)
          fm_bit_write(m_pIndex->var_byte_rappr, m_pIndex->start_lev2[i]);

	fm_bit_flush();
	write = m_pIndex->compress + m_pIndex->start_positions;
	fm_init_bit_writer(write, &m_pIndex->compress_size); // ritorno a scrivere in coda alla memoria
}


/*
	Se m_pIndex->skip>1 scrive nel compresso tutte le posizioni contenute in 
	m_pIndex->loc_occ utilizzando m_pIndex->log2textsize bits per ciascuna posizione.
    write the position of the marked chars
*/
void 
CFMIndex::write_locations(void) {

 	UINT32 i;
	fm_init_bit_buffer();
	if(m_pIndex->skip == 1) { // marcamento completo file 	
		for(i=0; i<m_pIndex->text_size; i++) 
			fm_bit_write(m_pIndex->log2textsize, (UINT32)m_pIndex->lf[i]);
		fm_bit_flush();
		#if TESTINFO
		Test.marked_pos = m_pIndex->log2textsize*m_pIndex->num_marked_rows/8;
		#endif
		return;
	}
	
	/* Marcamento sostituisci char bit log2 text_size*/
	for(i=0; i<m_pIndex->num_marked_rows;i++) 
		fm_bit_write(m_pIndex->log2textsize, (UINT32)(m_pIndex->loc_occ[i]-(m_pIndex->loc_occ[i]/(m_pIndex->skip+1)))-1);
		
	fm_bit_flush();

	#if TESTINFO
	Test.marked_pos = m_pIndex->log2textsize*m_pIndex->num_marked_rows/8;
	#endif
}

/* 
   Libera la memoria allocata.  Viene usata in caso di errore 
   o alla fine della compressione.
*/
void 
CFMIndex::dealloc(void) {
	free(m_pIndex->lf);
	m_pIndex->lf = NULL;
	dealloc_bucketinfo();
	free(m_pIndex->start_lev2);
	m_pIndex->start_lev2 = NULL;
	free(m_pIndex->loc_occ);
	m_pIndex->loc_occ = NULL;
	if ((m_pIndex->owner == 0) && (m_pIndex->skip !=0)) {
		if (m_pIndex->subchar != m_pIndex->specialchar) { 
			int i;
			for(i=0;i<m_pIndex->alpha_size;i++)
				if(m_pIndex->text[i] == m_pIndex->subchar) m_pIndex->text[i] = m_pIndex->specialchar;
			}				
	} else {free(m_pIndex->text);m_pIndex->text=NULL;}

}

void 
CFMIndex::dealloc_bucketinfo(void)
{
	UINT32 i;
	bucket_lev1 *sb;

	free(m_pIndex->bwt);
	m_pIndex->bwt = NULL;
	if(m_pIndex->buclist_lev1 != NULL)
			for(i=0; i< m_pIndex->num_bucs_lev1; i++) {
				sb = &(m_pIndex->buclist_lev1[i]); 
    			free(sb->occ);
    			free(sb->bool_char_map);
	}
	free(m_pIndex->buclist_lev1);
	m_pIndex->buclist_lev1 = NULL;
	
}	

int 
CFMIndex::errore(int error) {
	if(m_pIndex->compress != NULL)
		{
		free(m_pIndex->compress);
		m_pIndex->compress = NULL;
		}
	if(m_pIndex->mtf_seq != NULL)
		{
		free(m_pIndex->mtf_seq);
		m_pIndex->mtf_seq = NULL;
		}
	dealloc();
	return error;
}

/*
	Saves index on disk by using single or multiple files, having 
	proper extensions. 
*/
int 
CFMIndex::save_index(char *filename) {

	FILE *outfile;
	char *outfile_name, *ext = EXT;
	
	outfile_name = (char *) malloc(strlen(filename)+strlen(ext)+1);
	if (outfile_name == NULL) return FM_OUTMEM;
	outfile_name = strcpy(outfile_name, filename);
	outfile_name = strcat(outfile_name, ext);
	
	outfile = fopen(outfile_name, "wb");
	if(outfile == NULL) { 
		free(outfile_name);
		return FM_FILEERR;
	}

    	free(outfile_name);

	UINT32 t = (UINT32)fwrite(m_pIndex->compress, sizeof(unsigned char), m_pIndex->compress_size, outfile);
	if(t != m_pIndex->compress_size) {
		fclose(outfile);
		return FM_FILEERR;
	}
    
    fclose(outfile);
	return FM_OK;
}

int 
CFMIndex::save_index_mem(unsigned char *compress)
{
memcpy(compress, m_pIndex->compress, sizeof(unsigned char)*m_pIndex->compress_size);
return FM_OK;
}

/* 
  Opens and reads a text file

*/
int CFMIndex::fm_read_file(char *filename,		// file to read from 
				 unsigned char **textt,			// returned buffer allocated using malloc() containing contents of file + space for suffix sorting 
				 UINT32 *length)		// returned content length, excludes suffix sort length
{

  unsigned char *text;
  UINT32 t;
  FILE *infile;
  
  infile = fopen(filename, "rb"); // b is for binary: required by DOS
  if(infile == NULL) return FM_FILEERR;
  
  /* store input file length */
  if(fseek(infile,0,SEEK_END) !=0 ) return FM_FILEERR;
  *length = ftell(infile);
  
  /* alloc memory for text  */
  text = (unsigned char *)malloc(((*length))*sizeof(*text)); 
  if(text == NULL) return FM_OUTMEM;  
  
  /* read text in one sweep */
  rewind(infile);
  t = (UINT32)fread(text, sizeof(*text), (size_t) *length, infile);
  if(t!=*length) return FM_READERR;
  *textt = text;
  fclose(infile);
  return FM_OK;
}


