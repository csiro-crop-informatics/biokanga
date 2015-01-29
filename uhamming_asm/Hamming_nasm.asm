; on entry then expect parameters passed as
;		Windows Linux
; P1	RCX		RDI		pHDs			// where to return Hamming differentials for each subsequence
; P2	RDX		RSI		SubSeqLen		// generate Hammings edit distances for subsequences of this length
; p3	R8		RDX		SSofs			// offset between subsequences for current pass
; p4	R9		RCX		pGenomeSeq		// genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
; p5	RSP+40	R8		GenomeLen		// total genome length including chrom separators and genome terminator eBaseEOG

; within main iteration loop processing
; AL holds RelBase
; AH holds RefBase
; BL holds PPHD
; BH holds 
; CH holds PPEOScnt
; DL holds 
; RSI holds pRefSeq    (allows loads into AH/L...DH/L unlike the RN registers)
; R9 holds pRefSeq[RelSeqOfs]
; R10 holds pRefSeq[SubSeqLen1]
; R12 holds pHPlus[SSofs]
; R13 holds pRefSeq[SSofs]
; R15 holds pHPlus

	global	GHamDistWatson		
	SECTION .text

GHamDistWatson: 				; GHamDistWatson, COMDAT
	push	rbp
 	mov		rbp,rsp
	push    rbx
	push	r12
	push	r13
	push	r14
	push	r15

; reorder register parameter assignments to be same as for windows
; this is to allow the initialisation and main Hanmming generation loop to be
; consistent between the linux and windows implementation 
	mov     RAX,R8		;GenomeLen
	mov		R8,RDX		;SSofs
	mov		RDX,RSI		;SubSeqLen
	mov     R9,RCX		;pGenomeSeq
	mov		RCX,RDI		;pHDs
	
; initialisation here...
; P1 RCX = pHDs			// where to return Hamming differentials for each subsequence
; P2 RDX = SubSeqLen	// generate Hammings edit distances for subsequences of this length
; p3 R8 = SSofs			// offset between subsequences for current pass
; p4 R9 = pGenomeSeq	// genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
; p5 EAX = GenomeLen   
;
;pRefSeq = pGenomeSeq;
;SubSeqLen1 = SubSeqLen-1;
;RelSeqOfs = SSofs+SubSeqLen1;
;pHPlus = pHDs;
	mov	RBX, R8			; BX now holds SSofs
	mov RSI,R9			; RSI now holds pRefSeq
	lea R13,[R9+R8]		; R13 now holds pRefSeq[SSofs]
	mov R15,RCX			; R15 now holds pHPlus or pHDs
	lea R12,[R15+R8]	; R12 now holds pHPlus[SSofs]
	lea R9,[RSI+RBX-1] 
	add R9,RDX			; R9 now holds pRefSeq[RelSeqOfs]
	dec RDX				; RDX now holds  SubSeqLen1
	lea R10,[RSI+RDX]	; R10 now holds pRefSeq[SubSeqLen1]
; 
;generate Huffman for 1st SubSeqLen-1 bases
	xor rcx,rcx
	xor rbx,rbx
	push RSI			;pRefSeq
	push R13			;pRefSeq[SSofs]

initloop:				;for(Idx = 0; Idx < SubSeqLen-1; Idx++)
	mov AH,[RSI]		;   RefBase = pRefSeq
	inc RSI
	mov AL,[R13]		;	RelBase = pRefSeq[SSofs];
	inc R13
	cmp AH,07			;	if(RefBase == eBaseEOS)
	jne	skipinc
	add CH,1				;		PPEOScnt += 1;
skipinc:
	cmp AL,07			;	if(RelBase == eBaseEOS)
	jne skipinc1
	add CH,1				;		PPEOScnt += 1;
skipinc1:
	cmp AL,AH		;	if(RefBase != RelBase)
	je skipinc2
	inc BL				;		PPHD += 1;
skipinc2:
	dec DL
	jnz initloop
	pop R13			;pRefSeq[SSofs]
	pop RSI			;pRefSeq
	
	
;Huffman initialisation for 1st SubSeqLen-1 bases completed, now for main loop processing -		
mainloop:
	mov AL,  [R10]	; RefBase = pRefSeq[SubSeqLen1]
	mov AH,AL
	inc R10					; pRefSeq[SubSeqLen1] ready for next iteration
	mov AL,  [R9]	; RelBase = pRefSeq[RelSeqOfs]
	inc R9					; pRefSeq[RelSeqOfs] ready for next iteration
	cmp AL, 07				; 
	jl noinc1				; if(RelBase < eBaseEOS)
	jg bailout				;if(RelBase > eBaseEOS) processing completed
	add CH,1					; if(RelBase == eBaseEOS) PPEOScnt += 1
noinc1:
	cmp AH,07				; if(RefBase == eBaseEOS)
	jne noinc3
	add CH,1				; PPEOScnt += 1;
noinc3:
	cmp	AL,AH				; if(RefBase != RelBase)				
	je  noinc2
	add BL,1				; PPHD += 1;
noinc2:
	test CH,CH				; if(!PPEOScnt)
	jne	nocmp1
	cmp BL, [R15]	; if(*pHPlus > PPHD)
	jge xx1
	mov  [R15],BL	; *pHPlus = PPHD
xx1: 
	cmp BL,[R12]	; if(pHPlus[SSofs] > PPHD)
	jge nocmp1
	mov [R12],BL			; pHPlus[SSofs] = PPHD
nocmp1:
	inc R12						; pHPlus[SSofs] incremented
	inc R15						; pHPlus += 1;

	mov AH, [RSI]	; RefBase = *pRefSeq;
	inc RSI					; pRefSeq++
	mov AL, [R13]	; RelBase = pRefSeq[SSofs];
	inc R13					; pRefSeq[SSofs] incremented ready for next iteration
	
	cmp AH, 07				;if(RefBase == eBaseEOS)
	jne nodec1
	add CH,-1					; PPEOScnt -= 1;
nodec1:
	cmp AH,AL				;if(RefBase != RelBase)		
	je zzzzz
	add BL,-1					;PPHD -= 1;
zzzzz:
	cmp AL, 07				;if(RelBase == eBaseEOS)
	jne mainloop
	add CH,-1					; PPEOScnt -= 1;
	jmp mainloop
		
bailout:	
	pop r15
	pop r14
	pop r13
	pop r12
	pop rbx
	pop rbp
	xor rax,rax
	ret

	global	GHamDistCrick		
	SECTION .text

GHamDistCrick: 				; GHamDistCrick, COMDAT
	push	rbp
 	mov		rbp,rsp
	push    rbx
	push	r12
	push	r13
	push	r14
	push	r15

; reorder register parameter assignments to be same as for windows
; this is to allow the initialisation and main Hanmming generation loop to be
; consistent between the linux and windows implementation 
	mov     RAX,R8		;GenomeLen
	mov		R8,RDX		;SSofs
	mov		RDX,RSI		;SubSeqLen
	mov     R9,RCX		;pGenomeSeq
	mov		RCX,RDI		;pHDs

; initialisation here...
; P1 RCX = pHDs			// where to return Hamming differentials for each subsequence
; P2 RDX = SubSeqLen	// generate Hammings edit distances for subsequences of this length
; p3 R8 = SSofs			// offset between subsequences for current pass
; p4 R9 = pGenomeSeq	// genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
; p5 EAX = GenomeLen   
;
; Initialise such that:
;CL:  SubSeqLen
;R8:  LastIdx = GenomeLen - SSofs - SubSeqLen;
;R12: pHSS1 = pHDs;
;R13: pHSS2 = &pHDs[SSofs];
;R14: pSS1 = pGenomeSeq;
;R15: pSS2 = &pGenomeSeq[SSofs+SubSeqLen-1];
	mov R12,RCX			; R12 now pHSS1 = pHDs
	lea R13,[R12+R8]	; R13 now holds pHSS2 = &pHDs[SSofs]
	mov	R14,R9			; R14  now holds pSS1 = pGenomeSeq
	lea	R15,[R9+R8-1]
	add R15,RDX			; R15 now holds pSS2 = &pGenomeSeq[SSofs+SubSeqLen-1]
	sub	RAX,R8
	sub	RAX,RDX
	mov R8,RAX			; R8 now holds LastIdx = GenomeLen - SSofs - SubSeqLen
	mov CL,DL			; subseqlen now in CL
; Register usage
;R8 has number of remaining main loop iterations
;SubSeqLen in CL
;MaxPMH in BL
;BIdx in BH
;PPHD2 in CH
;PMHD1 in DL
;PPHD1 in DH
;pS1 in RDI
;pS2 in RSI
;pHSS1 in R12
;pHSS2 in R13
;pSS1 in R14
;pSS2 in R15
mainloop:
	mov		al,[R12]			;if((PPHD1 = *pHSS1) > SubSeqLen) continue	
	cmp		al,cl
	ja		IterMainLoop
	mov		dh,al						; dh to hold current PPHD1
	mov		al, [R13]			;if((PPHD2 = *pHSS2) > SubSeqLen) continue
	cmp		al,cl
	ja		IterMainLoop
	mov		ch,al						; ch to hold current PPHD2
	mov		bl,al						; assume PPHD2 will be the max(PPHD2,PPHD1)	
	cmp		bl,dh						
	ja		chkmax0
	mov		bl,dh						; PPHD1 is the max(PPHD2,PPHD1)
chkmax0:
	test	bl,bl						;if the max(PPHD2,PPHD1) is already down at 0 then can't do much better
	jz		IterMainLoop				;next iteration if MaxPMH == 0
; need to now gen the Hamming between pSS1 and the revcpl of pSS2
	mov		RDI,R14						;pS1 = pSS1;
	mov		RSI,R15						;pS2 = pSS2;
	mov		dl,0						;initialise PMHD1 = 0
	mov		BH,CL						;for(BIdx = 0; BIdx < SubSeqLen; BIdx++)
NxtSeqBase:
	mov		ah,[RSI]			; al = *pS2--
	dec		RSI
	mov		al,[RDI]			; al = *pS1++
	inc		RDI
	cmp		ah,03				; MapCpl(*pS2)
	ja		cmpNucs
	xor		ah,03
cmpNucs:
	cmp		al,ah				;if(*pS1 != MapCpl[*pS2] && (++PMHD1 > MaxPMH)) break
	jz		IterNxtSeqBase
	add		DL,1				; ++PMHD1
	cmp		DL,BL				; if PMHD1 > MaxPMH then break
	jg		SeqCmpd
IterNxtSeqBase:
	dec		bh
	jnz		NxtSeqBase
	
SeqCmpd:	
	cmp	    dl,dh				;if(PMHD1 < PPHD1) *pHSS1 = PMHD1;
	jae		chk
	mov		[R12],dl
chk:	
	cmp		dl,ch				;if(PMHD1 < PPHD2) *pHSS2 = PMHD1;
	jae		IterMainLoop
	mov		[R13],dl
IterMainLoop:
	inc		R12
	inc		R13
	inc		R14
	inc		R15
	sub		R8,1
	jnz		mainloop
	
bailout:	
	pop r15
	pop r14
	pop r13
	pop r12
	pop rbx
	pop rbp
	xor rax,rax
	ret