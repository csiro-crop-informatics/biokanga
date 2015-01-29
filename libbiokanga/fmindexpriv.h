//----- Full credit is given to the following authors, this application is very dependent on their code
//FM-index Version 2
//Authors: Paolo Ferragina and Rossano Venturini\n");
//Dipartimento di Informatica, University of Pisa, Italy\n");
//
#pragma once

#include <stdlib.h>
#ifndef __FM_COMMON
#define __FM_COMMON
#define TESTINFO (0)

/* Some useful macro */
#define EOF_shift(n) (n < m_pIndex->bwt_eof_pos) ? n+1 :  n
#define MIN(a, b) ((a)<=(b) ? (a) : (b))
#define MAX(a, b) ((a)<=(b) ? (a) : (b))

#endif


/* Errors type */
#ifndef __FM_ERRORS
#define __FM_ERRORS
/* define errors */
#define FM_OK 		( 0)
#define FM_GENERR 	(-1)	   	// Errore generale
#define FM_OUTMEM 	(-2)	   	// Errore allocazione memoria
#define FM_CONFERR 	(-11)      	// Errore nella configurazione 

/* read/write errors */
#define FM_READERR  (-12) 		// Errore nella lettura file 
#define FM_FILEERR 	(-14) 		// Problem with some file

/* decompress errors */
#define FM_FORMATERR 	(-3)  	// Formato file da decomprimere errato
#define FM_COMPNOTSUP 	(-4) 	// Metodo di compressione o decompressione non supportato
#define FM_DECERR 		(-5)    // Errore in decompressione 
#define FM_COMPNOTCORR 	(-8) 	// File compresso non corretto

/* locate/count/extract errors */
#define FM_SEARCHERR 	(-6)    // Errore nella ricerca causato da posizioni errate
#define FM_NOMARKEDCHAR (-7) 	// Nessun carattere marcato impossibile search/extract
#define FM_MARKMODENOTKNOW (-9) // Tipo di marcamento sconosciuto

#define FM_NOTIMPL (-13) 		// Function non implemented
#endif

#ifndef INTERNALDATATYPE 

/* Tipo di compressione */
#define MULTIH 	 (4) 	/* MTF + Wheeler 1/2 + Multi Huffman 	*/

/* constants used for searching in fmi files (see search_main()) */
#define NULL_CHAR      (0)	/* used to skip the count of char-occs */
#define WHAT_CHAR_IS   (1)	/* used to retrieve the char in a given pos */
#define COUNT_CHAR_OCC (2)	/* used to count char-occs before a given pos */

#define INTERNALDATATYPE
#endif

#ifndef _BZLIB_PRIVATE_H
#define _BZLIB_PRIVATE_H

#include <stdlib.h>

#ifndef BZ_NO_STDIO
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#endif

/*-- General stuff. --*/

#define BZ_VERSION  "0.9.5a"

#ifndef __GNUC__
#define __inline__  /* */
#endif 

/*-- Constants for the back end. --*/

#define SMALLFILESIZE (51201)
// 50 Kbytes
#define SMALLSMALLFILESIZE (1025)

#define BZ_RUNA 0
#define BZ_RUNB 1


#define BZ_G_SIZE   50
#define BZ_N_ITERS  4

#define BZ_MAX_SELECTORS (2 + (900000 / BZ_G_SIZE))



/*-- Stuff for randomising repetitive blocks. --*/

extern UINT32 rNums[512];

#define BZ_RAND_DECLS                          \
   UINT32 rNToGo;                               \
   UINT32 rTPos                                 \

#define BZ_RAND_INIT_MASK                      \
   s->rNToGo = 0;                              \
   s->rTPos  = 0                               \

#define BZ_RAND_MASK ((s->rNToGo == 1) ? 1 : 0)

#define BZ_RAND_UPD_MASK                       \
   if (s->rNToGo == 0) {                       \
      s->rNToGo = rNums[s->rTPos];             \
      s->rTPos++;                              \
      if (s->rTPos == 512) s->rTPos = 0;       \
   }                                           \
   s->rNToGo--;



/*-- Stuff for doing CRCs. --*/

extern unsigned int crc32Table[256];

#define BZ_INITIALISE_CRC(crcVar)              \
{                                              \
   crcVar = 0xffffffffL;                       \
}

#define BZ_FINALISE_CRC(crcVar)                \
{                                              \
   crcVar = ~(crcVar);                         \
}

#define BZ_UPDATE_CRC(crcVar,cha)              \
{                                              \
   crcVar = (crcVar << 8) ^                    \
            crc32Table[(crcVar >> 24) ^        \
                       ((unsigned char)cha)];          \
}



/*-- States and modes for compression. --*/

#define BZ_M_IDLE      1
#define BZ_M_RUNNING   2
#define BZ_M_FLUSHING  3
#define BZ_M_FINISHING 4

#define BZ_S_OUTPUT    1
#define BZ_S_INPUT     2

#define BZ_N_RADIX 2
#define BZ_N_QSORT 12
#define BZ_N_SHELL 18
#define BZ_N_OVERSHOOT (BZ_N_RADIX + BZ_N_QSORT + BZ_N_SHELL + 2)

/*-- externs for compression. --*/
extern void 
hbAssignCodes ( UINT32*, unsigned char*, UINT32, UINT32, UINT32 );

extern void 
hbMakeCodeLengths ( unsigned char*, UINT32*, UINT32, UINT32 );



/*-- states for decompression. --*/

#define BZ_X_IDLE        1
#define BZ_X_OUTPUT      2

#define BZ_X_MAGIC_1     10
#define BZ_X_MAGIC_2     11
#define BZ_X_MAGIC_3     12
#define BZ_X_MAGIC_4     13
#define BZ_X_BLKHDR_1    14
#define BZ_X_BLKHDR_2    15
#define BZ_X_BLKHDR_3    16
#define BZ_X_BLKHDR_4    17
#define BZ_X_BLKHDR_5    18
#define BZ_X_BLKHDR_6    19
#define BZ_X_BCRC_1      20
#define BZ_X_BCRC_2      21
#define BZ_X_BCRC_3      22
#define BZ_X_BCRC_4      23
#define BZ_X_RANDBIT     24
#define BZ_X_ORIGPTR_1   25
#define BZ_X_ORIGPTR_2   26
#define BZ_X_ORIGPTR_3   27
#define BZ_X_MAPPING_1   28
#define BZ_X_MAPPING_2   29
#define BZ_X_SELECTOR_1  30
#define BZ_X_SELECTOR_2  31
#define BZ_X_SELECTOR_3  32
#define BZ_X_CODING_1    33
#define BZ_X_CODING_2    34
#define BZ_X_CODING_3    35
#define BZ_X_MTF_1       36
#define BZ_X_MTF_2       37
#define BZ_X_MTF_3       38
#define BZ_X_MTF_4       39
#define BZ_X_MTF_5       40
#define BZ_X_MTF_6       41
#define BZ_X_ENDHDR_2    42
#define BZ_X_ENDHDR_3    43
#define BZ_X_ENDHDR_4    44
#define BZ_X_ENDHDR_5    45
#define BZ_X_ENDHDR_6    46
#define BZ_X_CCRC_1      47
#define BZ_X_CCRC_2      48
#define BZ_X_CCRC_3      49
#define BZ_X_CCRC_4      50



/*-- Constants for the fast MTF decoder. --*/

#define MTFA_SIZE 4096
#define MTFL_SIZE 16




/*-- Macros for decompression. --*/

#define BZ_GET_FAST(cccc)                     \
    s->tPos = s->tt[s->tPos];                 \
    cccc = (unsigned char)(s->tPos & 0xff);           \
    s->tPos >>= 8;

#define BZ_GET_FAST_C(cccc)                   \
    c_tPos = c_tt[c_tPos];                    \
    cccc = (unsigned char)(c_tPos & 0xff);            \
    c_tPos >>= 8;

#define SET_LL4(i,n)                                          \
   { if (((i) & 0x1) == 0)                                    \
        s->ll4[(i) >> 1] = (s->ll4[(i) >> 1] & 0xf0) | (n); else    \
        s->ll4[(i) >> 1] = (s->ll4[(i) >> 1] & 0x0f) | ((n) << 4);  \
   }

#define GET_LL4(i)                             \
   ((((unsigned int)(s->ll4[(i) >> 1])) >> (((i) << 2) & 0x4)) & 0xF)

#define SET_LL(i,n)                          \
   { s->ll16[i] = (UINT16)(n & 0x0000ffff);  \
     SET_LL4(i, n >> 16);                    \
   }

#define GET_LL(i) \
   (((unsigned int)s->ll16[i]) | (GET_LL4(i) << 16))

#define BZ_GET_SMALL(cccc)                        \
      cccc = indexIntoF ( s->tPos, s->cftab );    \
      s->tPos = GET_LL(s->tPos);


/*-- externs for decompression. --*/

extern UINT32 
indexIntoF ( UINT32, UINT32* );


extern void 
hbCreateDecodeTables ( UINT32*, UINT32*, UINT32*, unsigned char*,
                       UINT32,  UINT32, UINT32 );

#endif


/*-- BZ_NO_STDIO seems to make NULL disappear on some platforms. --*/

#ifdef BZ_NO_STDIO
#ifndef NULL
#define NULL 0
#endif
#endif


/*-------------------------------------------------------------*/
/*--- end                                   bzlib_private.h ---*/
/*-------------------------------------------------------------*/



#ifndef __FM_MNGBITS
#define __FM_MNGBITS
/* M A C R O  */

#define fm_init_bit_buffer() {__Bit_buffer = 0; __Bit_buffer_size = 0;}

/* -----------------------------------------------------------------------------
   Funzioni per leggere/scrivere meno di 24 bits.
   L'accesso fuori da memoria allocata non e' gestito a questo livello ma e' 
   compito del chiamante. 
   Implementate come macro per migliorare le prestazione.
   ----------------------------------------------------------------------------- */

/* Uso Corretto write:
   BitsDaWrite
   if ( PosizioneinMemori + (log2(BitsDaWrite)+1) < BytesDisponibili ) 
   				ALLORA REALLOCA e RIFAI init_Byte_write;
   bit_write(BitsDaWrite,IntDaCodificare);
*/
							   
/* (numero bits da leggere, int che conterra' il risultato) */
#define fm_bit_read24(__n,__result)	{					\
	UINT32 __t,__u;										\
	assert(__Bit_buffer_size<8);						\
	assert(__n>0 && __n<=24);							\
	/* --- read groups of 8 bits until size>= n --- */ 	\
  	while(__Bit_buffer_size < __n) {					\
		__t = (UINT32) *__MemAddress;					\
		__MemAddress++; (*__Num_Bytes)++;				\
    	__Bit_buffer |= (__t << (24-__Bit_buffer_size));\
    	__Bit_buffer_size += 8;}						\
  	/* ---- write n top bits in u ---- */				\
  	__u = __Bit_buffer >> (32-__n);						\
  	/* ---- update buffer ---- */						\
  	__Bit_buffer <<= __n;								\
  	__Bit_buffer_size -= __n;							\
  	__result = ((int)__u);}					   

	
#define BIT_READ(__bitz, __result) {					\
  UINT32 __ur = 0;										\
  int __ir;												\
  assert(__bitz <= 32);									\
  if (__bitz > 24){										\
	fm_bit_read24((__bitz-24), __ir);					\
    __ur =  __ir<<24;									\
	fm_bit_read24(24,__ir);								\
    __ur |= __ir;										\
    __result = __ur;									\
  } else {												\
    fm_bit_read24(__bitz,__ir);							\
	__result = __ir;									\
  }}
	
									
/* -----------------------------------------------------------------------------
   Codifica di Interi con una parte fissa ed una variabile. 
   Questa codifica non so se esisteva ma ha prestazioni buone rispetto a write7x8
   Anche queste come macro perche' molto utilizzate in ogni parte di FM-Index
   ------------------------------------------------------------------------------
   */

#define fm_integer_encode(__numero, __log2log2maxvalue) {{		\
	UINT32  __k;										\
	UINT32 __i = 0;									\
	assert(__log2log2maxvalue<6);							\
	switch ( __numero ){									\
		/* primi 4 casi speciali */							\
		case 0: fm_bit_write24(__log2log2maxvalue,0); 			\
				break;										\
		case 1: fm_bit_write24(__log2log2maxvalue,1); 			\
				break;										\
		case 2: fm_bit_write24((__log2log2maxvalue+1),4);		\
				break;										\
		case 3: fm_bit_write24((__log2log2maxvalue+1),5); 		\
				break;										\
	/*	Altri casi. Si calcola i = log2(occ) __numero		\
		di bit per rappresentare occ. si calcola  k=		\
		occ-2^i cioe' la distanza dalla potenza di due		\ 
		inferiore. Scrive i su log2log2BuckSize bits		\
		e k su i-1 bits indica __numero in piu' alla		\
	 	potenza di due precedente 						*/	\
		default: {											\
			UINT32 __pow = 1;									\
			/* 	calcola distanza dalla potenza di due 		\	
			 	precedente i indica l'esponente di 			\
				tale potenza */ 							\
			__i = 3;										\
			__pow = 4;										\
			while(__pow <= __numero){						\
    			__pow= __pow<<1;							\
    			__i = __i+1;								\
  			}												\
			__i--; __pow = __pow>>1;						\
			__k = __numero - __pow; 						\
			fm_bit_write24(__log2log2maxvalue,__i);			\
			fm_bit_write24(__i-1,__k);							\
			}}}}
			
UINT32 fm_integer_decode(unsigned short int);
#endif


