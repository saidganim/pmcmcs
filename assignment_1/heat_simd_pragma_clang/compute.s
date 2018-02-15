	.file	"compute.c"
	.local	N
	.comm	N,4,4
	.local	M
	.comm	M,4,4
	.text
	.globl	do_compute
	.type	do_compute, @function
do_compute:
.LFB2:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$200, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	%rdi, -200(%rbp)
	movq	%rsi, -208(%rbp)
	movq	-200(%rbp), %rsi
	movq	8(%rsi), %rsi
	subq	$1, %rsi
	movq	%rsi, -112(%rbp)
	movq	-200(%rbp), %rsi
	movq	8(%rsi), %rsi
	movq	%rsi, %r8
	movl	$0, %r9d
	movq	-200(%rbp), %rsi
	movq	8(%rsi), %rsi
	leaq	0(,%rsi,8), %rdi
	movq	%rdi, -216(%rbp)
	movq	-200(%rbp), %rsi
	movq	(%rsi), %rsi
	subq	$1, %rsi
	movq	%rsi, -120(%rbp)
	movq	-200(%rbp), %rsi
	movq	8(%rsi), %rsi
	movq	%rsi, %rcx
	movl	$0, %ebx
	movq	-200(%rbp), %rsi
	movq	(%rsi), %rsi
	movq	%rsi, %rax
	movl	$0, %edx
	movq	%rbx, %rdi
	imulq	%rax, %rdi
	movq	%rdx, %rsi
	imulq	%rcx, %rsi
	addq	%rdi, %rsi
	mulq	%rcx
	leaq	(%rsi,%rdx), %rcx
	movq	%rcx, %rdx
	movq	-200(%rbp), %rax
	movq	(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	8(%rax), %rax
	imulq	%rdx, %rax
	salq	$3, %rax
	movq	%rax, %rdi
	call	malloc
	movq	%rax, -128(%rbp)
	movq	-200(%rbp), %rax
	movq	8(%rax), %rax
	addq	$1, %rax
	movq	%rax, -136(%rbp)
	movq	-200(%rbp), %rax
	movq	8(%rax), %rax
	addq	$2, %rax
	movq	%rax, -240(%rbp)
	movq	$0, -232(%rbp)
	movq	-200(%rbp), %rax
	movq	8(%rax), %rax
	addq	$2, %rax
	leaq	0(,%rax,8), %rbx
	movq	-200(%rbp), %rax
	movq	(%rax), %rax
	addq	$1, %rax
	movq	%rax, -144(%rbp)
	movq	-200(%rbp), %rax
	movq	8(%rax), %rax
	addq	$2, %rax
	movq	%rax, %r14
	movl	$0, %r15d
	movq	-200(%rbp), %rax
	movq	(%rax), %rax
	addq	$2, %rax
	movq	%rax, %r12
	movl	$0, %r13d
	movq	%r15, %rdx
	imulq	%r12, %rdx
	movq	%r13, %rax
	imulq	%r14, %rax
	leaq	(%rdx,%rax), %rcx
	movq	%r14, %rax
	mulq	%r12
	addq	%rdx, %rcx
	movq	%rcx, %rdx
	movq	-200(%rbp), %rax
	movq	(%rax), %rax
	leaq	2(%rax), %rdx
	movq	-200(%rbp), %rax
	movq	8(%rax), %rax
	addq	$2, %rax
	imulq	%rdx, %rax
	salq	$3, %rax
	movq	%rax, %rdi
	call	malloc
	movq	%rax, -152(%rbp)
	movq	-200(%rbp), %rax
	movsd	40(%rax), %xmm1
	movsd	.LC0(%rip), %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -56(%rbp)
	pxor	%xmm0, %xmm0
	movsd	%xmm0, -64(%rbp)
	movl	$0, -68(%rbp)
	movq	-200(%rbp), %rax
	movq	(%rax), %rax
	movl	%eax, N(%rip)
	movq	-200(%rbp), %rax
	movq	8(%rax), %rax
	movl	%eax, M(%rip)
	movl	$0, -72(%rbp)
	jmp	.L2
.L3:
	movl	M(%rip), %eax
	movl	%eax, %eax
	leaq	0(,%rax,8), %rdx
	movq	-200(%rbp), %rax
	movq	48(%rax), %rax
	movl	M(%rip), %esi
	movl	-72(%rbp), %ecx
	imull	%esi, %ecx
	movl	%ecx, %ecx
	salq	$3, %rcx
	addq	%rax, %rcx
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-72(%rbp), %eax
	cltq
	imulq	%rsi, %rax
	addq	$1, %rax
	leaq	0(,%rax,8), %rsi
	movq	-152(%rbp), %rax
	addq	%rsi, %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	memcpy
	addl	$1, -72(%rbp)
.L2:
	movl	-72(%rbp), %edx
	movl	N(%rip), %eax
	cmpl	%eax, %edx
	jb	.L3
	movl	$0, -76(%rbp)
	jmp	.L4
.L5:
	movq	%rbx, %rcx
	shrq	$3, %rcx
	movq	-152(%rbp), %rax
	movl	-76(%rbp), %edx
	movslq	%edx, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-152(%rbp), %rax
	movl	-76(%rbp), %edx
	movslq	%edx, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	N(%rip), %eax
	leal	1(%rax), %r9d
	movq	%rbx, %rdi
	shrq	$3, %rdi
	movl	N(%rip), %r8d
	movq	-152(%rbp), %rax
	movl	-76(%rbp), %edx
	movslq	%edx, %rcx
	movl	%r8d, %edx
	imulq	%rdi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-152(%rbp), %rax
	movl	-76(%rbp), %edx
	movslq	%edx, %rcx
	movl	%r9d, %edx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	addl	$1, -76(%rbp)
.L4:
	movl	-76(%rbp), %edx
	movl	M(%rip), %eax
	cmpl	%eax, %edx
	jb	.L5
	movl	M(%rip), %edx
	movq	-152(%rbp), %rax
	movl	%edx, %edx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-152(%rbp), %rax
	movsd	%xmm0, (%rax)
	movl	M(%rip), %eax
	leal	1(%rax), %edx
	movq	-152(%rbp), %rax
	movsd	8(%rax), %xmm0
	movq	-152(%rbp), %rax
	movl	%edx, %edx
	movsd	%xmm0, (%rax,%rdx,8)
	movq	%rbx, %rcx
	shrq	$3, %rcx
	movl	N(%rip), %eax
	leal	1(%rax), %r8d
	movq	%rbx, %rdx
	shrq	$3, %rdx
	movl	N(%rip), %esi
	movl	M(%rip), %edi
	movq	-152(%rbp), %rax
	movl	%edi, %edi
	movl	%esi, %esi
	imulq	%rsi, %rdx
	addq	%rdi, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-152(%rbp), %rax
	movl	%r8d, %edx
	imulq	%rcx, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	movq	%rbx, %rcx
	shrq	$3, %rcx
	movl	N(%rip), %eax
	leal	1(%rax), %edi
	movl	M(%rip), %eax
	leal	1(%rax), %esi
	movq	%rbx, %rdx
	shrq	$3, %rdx
	movl	N(%rip), %eax
	leal	1(%rax), %r8d
	movq	-152(%rbp), %rax
	movl	%r8d, %r8d
	imulq	%r8, %rdx
	addq	$1, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-152(%rbp), %rax
	movl	%esi, %esi
	movl	%edi, %edx
	imulq	%rcx, %rdx
	addq	%rsi, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	leaq	-176(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	gettimeofday
	jmp	.L6
.L18:
	pxor	%xmm0, %xmm0
	movsd	%xmm0, -56(%rbp)
	movl	$0, -80(%rbp)
	jmp	.L7
.L8:
	movq	%rbx, %rcx
	shrq	$3, %rcx
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	M(%rip), %edx
	movq	-152(%rbp), %rax
	movl	%edx, %edi
	movl	-80(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rdi, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-152(%rbp), %rax
	movl	-80(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rcx, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	movq	%rbx, %rcx
	shrq	$3, %rcx
	movl	M(%rip), %eax
	leal	1(%rax), %edi
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-80(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	$1, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-152(%rbp), %rax
	movl	%edi, %esi
	movl	-80(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rcx, %rdx
	addq	%rsi, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	addl	$1, -80(%rbp)
.L7:
	movl	-80(%rbp), %edx
	movl	N(%rip), %eax
	cmpl	%eax, %edx
	jb	.L8
	movl	$1, -84(%rbp)
	jmp	.L9
.L14:
	movl	$1, -88(%rbp)
	jmp	.L10
.L13:
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edi
	movq	-152(%rbp), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rcx
	movslq	%edi, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm1
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-88(%rbp), %eax
	leal	-1(%rax), %edx
	movq	-152(%rbp), %rax
	movslq	%edx, %rcx
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	addsd	%xmm1, %xmm0
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	1(%rax), %edi
	movq	-152(%rbp), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rcx
	movslq	%edi, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm1
	addsd	%xmm1, %xmm0
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-88(%rbp), %eax
	leal	1(%rax), %edx
	movq	-152(%rbp), %rax
	movslq	%edx, %rcx
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm1
	addsd	%xmm1, %xmm0
	movsd	.LC2(%rip), %xmm1
	mulsd	%xmm0, %xmm1
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edx
	movl	-88(%rbp), %eax
	leal	-1(%rax), %ecx
	movq	-152(%rbp), %rax
	movslq	%ecx, %rcx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm2
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	1(%rax), %edx
	movl	-88(%rbp), %eax
	leal	-1(%rax), %ecx
	movq	-152(%rbp), %rax
	movslq	%ecx, %rcx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	addsd	%xmm2, %xmm0
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edx
	movl	-88(%rbp), %eax
	leal	1(%rax), %ecx
	movq	-152(%rbp), %rax
	movslq	%ecx, %rcx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm2
	addsd	%xmm2, %xmm0
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	1(%rax), %edx
	movl	-88(%rbp), %eax
	leal	1(%rax), %ecx
	movq	-152(%rbp), %rax
	movslq	%ecx, %rcx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm2
	addsd	%xmm2, %xmm0
	movsd	.LC3(%rip), %xmm2
	mulsd	%xmm2, %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -160(%rbp)
	movq	-200(%rbp), %rax
	movq	56(%rax), %rax
	movl	-84(%rbp), %edx
	subl	$1, %edx
	movl	%edx, %ecx
	movl	M(%rip), %edx
	imull	%edx, %ecx
	movl	-88(%rbp), %edx
	addl	%ecx, %edx
	subl	$1, %edx
	movl	%edx, %edx
	salq	$3, %rdx
	addq	%rdx, %rax
	movsd	(%rax), %xmm1
	movsd	.LC0(%rip), %xmm0
	subsd	%xmm1, %xmm0
	movsd	.LC4(%rip), %xmm1
	divsd	%xmm1, %xmm0
	movsd	-160(%rbp), %xmm1
	mulsd	%xmm1, %xmm0
	movsd	%xmm0, -160(%rbp)
	movq	-216(%rbp), %r10
	movq	%r10, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edi
	movl	-88(%rbp), %eax
	leal	-1(%rax), %r8d
	movq	%rbx, %r9
	shrq	$3, %r9
	movq	-152(%rbp), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rcx
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%r9, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm1
	movq	-200(%rbp), %rax
	movq	56(%rax), %rax
	movl	-84(%rbp), %edx
	subl	$1, %edx
	movl	%edx, %ecx
	movl	M(%rip), %edx
	imull	%edx, %ecx
	movl	-88(%rbp), %edx
	addl	%ecx, %edx
	subl	$1, %edx
	movl	%edx, %edx
	salq	$3, %rdx
	addq	%rdx, %rax
	movsd	(%rax), %xmm0
	mulsd	%xmm1, %xmm0
	movq	-128(%rbp), %rax
	movslq	%r8d, %rcx
	movslq	%edi, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	movq	%r10, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edi
	movl	-88(%rbp), %eax
	leal	-1(%rax), %r8d
	movq	%r10, %r9
	shrq	$3, %r9
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edx
	movl	-88(%rbp), %eax
	leal	-1(%rax), %ecx
	movq	-128(%rbp), %rax
	movslq	%ecx, %rcx
	movslq	%edx, %rdx
	imulq	%r9, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	addsd	-160(%rbp), %xmm0
	movq	-128(%rbp), %rax
	movslq	%r8d, %rcx
	movslq	%edi, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	%xmm0, (%rax,%rdx,8)
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rcx
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	%r10, %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edx
	movl	-88(%rbp), %eax
	leal	-1(%rax), %ecx
	movq	-128(%rbp), %rax
	movslq	%ecx, %rcx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm1
	subsd	%xmm1, %xmm0
	movsd	.LC5(%rip), %xmm1
	andpd	%xmm1, %xmm0
	ucomisd	-56(%rbp), %xmm0
	jbe	.L11
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-88(%rbp), %edx
	movslq	%edx, %rcx
	movl	-84(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-216(%rbp), %rsi
	shrq	$3, %rsi
	movl	-84(%rbp), %eax
	leal	-1(%rax), %edx
	movl	-88(%rbp), %eax
	leal	-1(%rax), %ecx
	movq	-128(%rbp), %rax
	movslq	%ecx, %rcx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm1
	subsd	%xmm1, %xmm0
	movsd	.LC5(%rip), %xmm1
	andpd	%xmm1, %xmm0
	movsd	%xmm0, -56(%rbp)
.L11:
	addl	$1, -88(%rbp)
.L10:
	movl	-88(%rbp), %edx
	movl	M(%rip), %eax
	cmpl	%eax, %edx
	jbe	.L13
	addl	$1, -84(%rbp)
.L9:
	movl	-84(%rbp), %edx
	movl	N(%rip), %eax
	cmpl	%eax, %edx
	setbe	%al
	testb	%al, %al
	jne	.L14
	movl	$0, -92(%rbp)
	jmp	.L15
.L16:
	movl	M(%rip), %eax
	movl	%eax, %eax
	leaq	0(,%rax,8), %rdx
	movq	%rbx, %rcx
	shrq	$3, %rcx
	movl	-92(%rbp), %eax
	addl	$1, %eax
	cltq
	imulq	%rcx, %rax
	addq	$1, %rax
	leaq	0(,%rax,8), %rcx
	movq	-152(%rbp), %rax
	addq	%rax, %rcx
	movq	-216(%rbp), %rsi
	shrq	$3, %rsi
	movl	-92(%rbp), %eax
	cltq
	imulq	%rsi, %rax
	addq	$1, %rax
	leaq	0(,%rax,8), %rsi
	movq	-128(%rbp), %rax
	addq	%rsi, %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	memcpy
	addl	$1, -92(%rbp)
.L15:
	movl	-92(%rbp), %edx
	movl	N(%rip), %eax
	cmpl	%eax, %edx
	jb	.L16
.L6:
	movl	-68(%rbp), %eax
	leal	1(%rax), %edx
	movl	%edx, -68(%rbp)
	movl	%eax, %edx
	movq	-200(%rbp), %rax
	movq	16(%rax), %rax
	cmpq	%rax, %rdx
	jnb	.L17
	movq	-200(%rbp), %rax
	movsd	40(%rax), %xmm1
	movsd	-56(%rbp), %xmm0
	ucomisd	%xmm1, %xmm0
	ja	.L18
.L17:
	leaq	-192(%rbp), %rax
	movl	$0, %esi
	movq	%rax, %rdi
	call	gettimeofday
	movq	%rbx, %rdx
	shrq	$3, %rdx
	movq	-152(%rbp), %rax
	addq	$1, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-208(%rbp), %rax
	movsd	%xmm0, 16(%rax)
	movq	-208(%rbp), %rax
	movsd	16(%rax), %xmm0
	movq	-208(%rbp), %rax
	movsd	%xmm0, 8(%rax)
	movl	$1, -96(%rbp)
	jmp	.L19
.L26:
	movl	$1, -100(%rbp)
	jmp	.L20
.L25:
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-100(%rbp), %edx
	movslq	%edx, %rcx
	movl	-96(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-208(%rbp), %rax
	movsd	16(%rax), %xmm1
	ucomisd	%xmm1, %xmm0
	jbe	.L21
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-100(%rbp), %edx
	movslq	%edx, %rcx
	movl	-96(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-208(%rbp), %rax
	movsd	%xmm0, 16(%rax)
.L21:
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-100(%rbp), %edx
	movslq	%edx, %rcx
	movl	-96(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm1
	movq	-208(%rbp), %rax
	movsd	8(%rax), %xmm0
	ucomisd	%xmm1, %xmm0
	jbe	.L23
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-100(%rbp), %edx
	movslq	%edx, %rcx
	movl	-96(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movq	-208(%rbp), %rax
	movsd	%xmm0, 8(%rax)
.L23:
	movq	%rbx, %rsi
	shrq	$3, %rsi
	movq	-152(%rbp), %rax
	movl	-100(%rbp), %edx
	movslq	%edx, %rcx
	movl	-96(%rbp), %edx
	movslq	%edx, %rdx
	imulq	%rsi, %rdx
	addq	%rcx, %rdx
	movsd	(%rax,%rdx,8), %xmm0
	movsd	-64(%rbp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, -64(%rbp)
	addl	$1, -100(%rbp)
.L20:
	movl	-100(%rbp), %edx
	movl	M(%rip), %eax
	cmpl	%eax, %edx
	jbe	.L25
	addl	$1, -96(%rbp)
.L19:
	movl	-96(%rbp), %edx
	movl	N(%rip), %eax
	cmpl	%eax, %edx
	jbe	.L26
	movq	-192(%rbp), %rax
	pxor	%xmm1, %xmm1
	cvtsi2sdq	%rax, %xmm1
	movq	-184(%rbp), %rax
	pxor	%xmm0, %xmm0
	cvtsi2sdq	%rax, %xmm0
	movsd	.LC6(%rip), %xmm2
	divsd	%xmm2, %xmm0
	addsd	%xmm1, %xmm0
	movq	-176(%rbp), %rax
	pxor	%xmm2, %xmm2
	cvtsi2sdq	%rax, %xmm2
	movq	-168(%rbp), %rax
	pxor	%xmm1, %xmm1
	cvtsi2sdq	%rax, %xmm1
	movsd	.LC6(%rip), %xmm3
	divsd	%xmm3, %xmm1
	addsd	%xmm2, %xmm1
	subsd	%xmm1, %xmm0
	movq	-208(%rbp), %rax
	movsd	%xmm0, 40(%rax)
	movl	-68(%rbp), %eax
	subl	$1, %eax
	movl	%eax, %edx
	movq	-208(%rbp), %rax
	movq	%rdx, (%rax)
	movl	N(%rip), %edx
	movl	M(%rip), %eax
	imull	%edx, %eax
	movl	%eax, %eax
	testq	%rax, %rax
	js	.L27
	pxor	%xmm0, %xmm0
	cvtsi2sdq	%rax, %xmm0
	jmp	.L28
.L27:
	movq	%rax, %rdx
	shrq	%rdx
	andl	$1, %eax
	orq	%rax, %rdx
	pxor	%xmm0, %xmm0
	cvtsi2sdq	%rdx, %xmm0
	addsd	%xmm0, %xmm0
.L28:
	movsd	-64(%rbp), %xmm1
	divsd	%xmm0, %xmm1
	movapd	%xmm1, %xmm0
	movq	-208(%rbp), %rax
	movsd	%xmm0, 32(%rax)
	movq	-208(%rbp), %rax
	movsd	-56(%rbp), %xmm0
	movsd	%xmm0, 24(%rax)
	nop
	addq	$200, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE2:
	.size	do_compute, .-do_compute
	.section	.rodata
	.align 8
.LC0:
	.long	0
	.long	1072693248
	.align 8
.LC2:
	.long	855738472
	.long	1071824579
	.align 8
.LC3:
	.long	2583490355
	.long	1071284857
	.align 8
.LC4:
	.long	0
	.long	1074790400
	.align 16
.LC5:
	.long	4294967295
	.long	2147483647
	.long	0
	.long	0
	.align 8
.LC6:
	.long	0
	.long	1093567616
	.ident	"GCC: (GNU) 6.4.1 20170727 (Red Hat 6.4.1-1)"
	.section	.note.GNU-stack,"",@progbits
