	.file	"cHamming.cpp"
	.text
	.p2align 4,,15
.globl _Z8GHamDistjiPhS_j
	.type	_Z8GHamDistjiPhS_j, @function
_Z8GHamDistjiPhS_j:
.LFB187:
	.cfi_startproc
	.cfi_personality 0x3,__gxx_personality_v0
	pushq	%r15
	.cfi_def_cfa_offset 16
	mov	%edi, %eax
	subl	%edi, %r8d
	cmpl	$1, %esi
	movl	%esi, %r15d
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	pushq	%r13
	.cfi_def_cfa_offset 32
	pushq	%r12
	.cfi_def_cfa_offset 40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	leaq	(%rcx,%rax), %rbp
	.cfi_offset 6, -48
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	pushq	%rbx
	.cfi_def_cfa_offset 56
	movq	%rax, -32(%rsp)
	mov	%r8d, %eax
	movq	%rdx, -40(%rsp)
	leaq	(%rcx,%rax), %rax
	movl	%r8d, -44(%rsp)
	movq	%rax, -24(%rsp)
	jle	.L29
	.cfi_offset 3, -56
	leal	-2(%r15), %r8d
	movq	-24(%rsp), %r14
	movq	-24(%rsp), %rbx
	movq	%rbp, %r10
	movq	%rcx, %r9
	xorl	%eax, %eax
	movq	%r8, -16(%rsp)
	notq	%r8
	xorl	%edx, %edx
	movq	%r8, -8(%rsp)
	addq	%r8, %r14
	xorl	%esi, %esi
	xorl	%edi, %edi
	jmp	.L9
	.p2align 4,,10
	.p2align 3
.L4:
	cmpb	$7, %r12b
	sete	%r13b
	addl	%r13d, %edx
	cmpb	$7, %r11b
	movzbl	%r11b, %r11d
	sete	%r13b
	addl	%r13d, %eax
	cmpb	%r12b, %r8b
	setne	%r12b
	addl	%r12d, %edi
	cmpb	%r8b, MapCpl(%r11)
	setne	%r8b
	addl	%r8d, %esi
	cmpq	%r14, %rbx
	je	.L30
.L9:
	movzbl	(%r9), %r8d
	movzbl	(%r10), %r12d
	addq	$1, %r9
	addq	$1, %r10
	movzbl	(%rbx), %r11d
	subq	$1, %rbx
	cmpb	$7, %r8b
	jne	.L4
	addl	$1, %edx
	addl	$1, %eax
	jmp	.L4
	.p2align 4,,10
	.p2align 3
.L30:
	movq	-16(%rsp), %rbx
	movq	-24(%rsp), %r9
	addq	-8(%rsp), %r9
	addq	$1, %rbx
	leaq	(%rcx,%rbx), %r13
	leaq	(%rbp,%rbx), %rbx
	movq	%r9, -16(%rsp)
.L3:
	movl	-44(%rsp), %r8d
	movq	-40(%rsp), %r14
	movq	%rcx, %r12
	addq	-32(%rsp), %r14
	movq	-40(%rsp), %r11
	subl	%r15d, %r8d
	movl	%r8d, %r15d
	xorl	%r8d, %r8d
	addl	$1, %r15d
	addq	-40(%rsp), %r15
	movq	%r14, -32(%rsp)
	movq	%r15, -8(%rsp)
	jmp	.L10
	.p2align 4,,10
	.p2align 3
.L13:
	cmpb	%r10b, %r9b
	movzbl	%r15b, %r15d
	setne	%r10b
	addl	%r10d, %edi
	cmpb	%r9b, MapCpl(%r15)
	setne	%r9b
	addl	%r9d, %esi
	testb	%dl, %dl
	jne	.L16
	cmpb	(%r11), %dil
	jae	.L17
	movb	%dil, (%r11)
.L17:
	movq	-32(%rsp), %r9
	cmpb	(%r9), %dil
	jae	.L16
	movb	%dil, (%r9)
.L16:
	testb	%al, %al
	jne	.L18
	cmpb	(%r11), %sil
	jae	.L19
	movb	%sil, (%r11)
.L19:
	movq	-8(%rsp), %r14
	cmpb	(%r14,%r8), %sil
	jae	.L18
	movb	%sil, (%r14,%r8)
.L18:
	cmpb	-49(%rsp), %cl
	setne	%r9b
	subb	%r9b, %dil
	movzbl	-40(%rsp), %r9d
	cmpb	%cl, MapCpl(%r9)
	setne	%r9b
	subb	%r9b, %sil
	cmpb	$7, %cl
	je	.L31
.L22:
	cmpb	$7, -49(%rsp)
	sete	%cl
	subb	%cl, %dl
	cmpb	$7, -40(%rsp)
	sete	%cl
	addq	$1, -32(%rsp)
	addq	$1, %r13
	subb	%cl, %al
	addq	$1, %r12
	addq	$1, %rbp
	addq	$1, %r11
	subq	$1, %r8
.L10:
	movzbl	(%rbx), %r10d
	addq	$1, %rbx
	cmpb	$7, %r10b
	ja	.L32
	movzbl	(%rbp), %r14d
	movq	-16(%rsp), %r9
	movzbl	(%r12), %ecx
	movzbl	(%r9,%r8), %r15d
	movzbl	(%r13), %r9d
	movb	%r14b, -49(%rsp)
	movq	-24(%rsp), %r14
	movzbl	(%r14,%r8), %r14d
	movb	%r14b, -40(%rsp)
	sete	%r14b
	addl	%r14d, %edx
	cmpb	$7, %r15b
	sete	%r14b
	addl	%r14d, %eax
	cmpb	$7, %r9b
	jne	.L13
	addl	$1, %edx
	addl	$1, %eax
	jmp	.L13
	.p2align 4,,10
	.p2align 3
.L31:
	subl	$1, %edx
	subl	$1, %eax
	jmp	.L22
	.p2align 4,,10
	.p2align 3
.L32:
	popq	%rbx
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
.L29:
	movq	%rax, -16(%rsp)
	movq	%rbp, %rbx
	movq	%rcx, %r13
	xorl	%eax, %eax
	xorl	%edx, %edx
	xorl	%esi, %esi
	xorl	%edi, %edi
	jmp	.L3
	.cfi_endproc
.LFE187:
	.size	_Z8GHamDistjiPhS_j, .-_Z8GHamDistjiPhS_j
.globl MapCpl
	.data
	.type	MapCpl, @object
	.size	MapCpl, 9
MapCpl:
	.byte	3
	.byte	2
	.byte	1
	.byte	0
	.byte	4
	.byte	5
	.byte	6
	.byte	7
	.byte	8
	.ident	"GCC: (Ubuntu 4.4.1-4ubuntu8) 4.4.1"
	.section	.note.GNU-stack,"",@progbits
