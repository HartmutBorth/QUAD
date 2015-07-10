# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 15-Apr-2015
# ============================================
# call ifft2s(a  ,c  ,trigs,n  ,lot)
#            (rdi,rsi,rdx  ,rcx,r8 )
#
# rdi: a(n,lot)  : real(8) source array
# rsi: c(n,lot)  : real(8) target array
# rdx: trigs(n)  : sine/cosine values
# rcx: n         : 1st. dimension
# r8 : lot       : 2nd. dimension
# --------------------------------------------

   .globl _ifft2s_
_ifft2s_:

   pushq   %r13
   pushq   %r12
   movl    (%r8 ), %r8d              # r8    = lot
   movl    (%rcx), %ecx              # rcx   = n
   movq    $0x3fe0000000000000, %rax # rax   = 0.5
   movd    %rax, %xmm8               # xmm8  = 0.5
   movq    %rcx, %r9                 # r9    = n
   shrq    $2, %r9                   # r9    = n / 4
   dec     %r9                       # r9    = (n / 4) - 1
   movq    %r9 , %r10                # r10   = loop count
   movq    %rdx, %r11                # r11   = &trigs
   movq    %rdi, %r13                # r13   = &a

LOOP_LOT_2S:
   movapd  %xmm8 , %xmm0             # xmm0  = 0.5
   mulsd   (%rdi), %xmm0             # xmm0  = 0.5 * a(0)
   movsd   %xmm0 , (%rsi)            # c(0)  = xmm0
   movsd   %xmm0 ,8(%rsi)            # c(1)  = xmm0
   leaq    -16(%rdi,%rcx,8), %r12    # r12   = &a(n-2) : ib
   addq    $16, %rdi                 # rdi   = &a(2)   : ia
   addq    $16, %rsi                 # rsi   = &c(2)   : j

LOOP_IFFT2:
   addq    $16, %rdx                 # &trigs
   movsd   (%rdi), %xmm7             # xmm7  = a(ia)
   movsd   (%r12), %xmm3             # xmm3  = a(ib)
   movsd   8(%r12), %xmm6            # xmm6  = a(ib+1)
   movapd  %xmm7, %xmm1              # xmm1  = a(ia)
   movsd   8(%rdi), %xmm2            # xmm2  = a(ia+1)
   subsd   %xmm3, %xmm1              # xmm1  = a(ia) - a(ib)
   addsd   %xmm7, %xmm3              # xmm3  = a(ia) + a(ib)
   movsd   (%rdx), %xmm4             # xmm4  = trigs(ia)
   movapd  %xmm2, %xmm0              # xmm0  = a(ia+1)
   subsd   %xmm6, %xmm2              # xmm2  = a(ia+1) - a(ib+1)
   movsd   8(%rdx), %xmm5            # xmm5  = trigs(ia+1)
   addsd   %xmm6, %xmm0              # xmm0  = a(ia+1) + a(ib+1)
   movsd   %xmm3, (%rsi)             # c(j)  = xmm3
   movapd  %xmm5, %xmm3              # xmm3  = s1 = trigs(ia+1)
   movsd   %xmm2, 16(%rsi)           # c(j+2) = 
   movapd  %xmm4, %xmm2
   mulsd   %xmm1, %xmm2
   mulsd   %xmm0, %xmm3
   mulsd   %xmm5, %xmm1
   mulsd   %xmm4, %xmm0
   subsd   %xmm3, %xmm2
   addsd   %xmm1, %xmm0
   movsd   %xmm2,  8(%rsi)           # c(j+1)
   movsd   %xmm0, 24(%rsi)           # c(j+3)
   addq    $32, %rsi                 # j  = j  + 4
   addq    $16, %rdi                 # ia = ia + 2
   subq    $16, %r12                 # ib = ib - 2
   dec     %r9
   jnz     LOOP_IFFT2

   movsd   (%rdi), %xmm0             # xmm0   = a(ia)
   movsd   %xmm0, (%rsi)             # c(n-2) = a(ia)
   xorpd   %xmm0, %xmm0              # xmm0   = 0
   subsd   8(%rdi), %xmm0            # xmm0   = -a(ia+1)
   movsd   %xmm0, 8(%rsi)            # c(n-1) = -a(ia+1)

   leaq    (%r13,%rcx,8), %rdi       # advance &a + n
   movq    %rdi, %r13                # save new base
   addq    $16, %rsi                 # advance &c + 2
   movq    %r11, %rdx                # restore &trigs
   movq    %r10, %r9                 # count for inner loop
   dec     %r8                       # lot
   jnz     LOOP_LOT_2S

   popq    %r12
   popq    %r13
   ret


# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 08-Jun-2015
# ============================================
# call ifft4s(a  ,c  ,trigs,n  ,lot)
#            (rdi,rsi,rdx  ,rcx, r8)
#
# rdi: a(n,lot)  : real(8) source array
# rsi: c(n,lot)  : real(8) target array
# rdx: trigs(n)  : sine/cosine values
# rcx: n         : 1st. dimension
# r8 : lot       : 2nd. dimension
# --------------------------------------------

   .globl   _ifft4s_
_ifft4s_:

   pushq   %rbp
   movq    %rsp, %rbp
   call    mcount
   pushq   %r15
   pushq   %r14
   pushq   %r13
   pushq   %r12
   pushq   %rbx
   movl    (%r8 ), %r15d             # r15  = lot
   movl    (%rcx), %ecx              # ecx  = n
   movq    %rcx, %rax                # rax  = n
   sarq    %rax                      # rax  = n / 2
   movq    %rcx, %r14                # r14  = n
   sarq    $3, %r14                  # r14  = n / 8
   dec     %r14                      # r14  = n / 8 - 1
   shlq    $3, %rcx                  # rcx  = n * 8

LOOP_O_IFFT4S:

   leaq     16(%rdi)       , %r10    # &a(i0)  i0 = 2
   leaq     16(%rdi,%rax,8), %r11    # &a(i1)  i1 = 2 + n/2
   leaq    -16(%rdi,%rcx)  , %r12    # &a(i2)  i2 = n - 2
   leaq    -32(%r11)       , %r13    # &a(i3)  i3 = i2 - n/2

   movq   $0x3fe0000000000000, %r8   # r8    = 0.5
   movd        %r8    , %xmm6        # xmm6  = 0.5
   movsd      (%rdi)  , %xmm0        # xmm0  = a(0)
   movsd   -16(%r11)  , %xmm4        # xmm4  = a(i1)
   movsd    -8(%r11)  , %xmm5        # xmm5  = a(i1+1)

   mulsd   %xmm6 , %xmm0             # xmm0 = a(0) * 0.5
   movapd  %xmm0 , %xmm1             # xmm1 = a(0) * 0.5
   movapd  %xmm0 , %xmm2             # xmm2 = a(0) * 0.5
   movapd  %xmm0 , %xmm3             # xmm3 = a(0) * 0.5

   addsd   %xmm4 , %xmm0             # xmm0 = a(0) + a(i1)
   subsd   %xmm5 , %xmm1             # xmm1 = a(0) - a(i1+1)
   subsd   %xmm4 , %xmm2             # xmm2 = a(0) - a(i1)
   addsd   %xmm5 , %xmm3             # xmm3 = a(0) + a(i1+1)

   movsd   %xmm0,   (%rsi)           # c(0)   = xmm0
   movsd   %xmm1,  8(%rsi)           # c(1)   = xmm1
   movsd   %xmm2, 16(%rsi)           # c(2)   = xmm2
   movsd   %xmm3, 24(%rsi)           # c(3)   = xmm3

   movq    %rdx, %r8                 # r8   = trigs
   movq    %rdx, %r9                 # r9   = trigs

   pushq   %rdx                      # save trigs
   pushq   %r14                      # save inner loop count

LOOP_I_IFFT4S:                                 

   movsd   (%r10), %xmm2             # xmm2  = a(i0)
   movsd   (%r11), %xmm3             # xmm3  = a(i1)
   movsd   (%r12), %xmm4             # xmm4  = a(i2)
   movsd   (%r13), %xmm5             # xmm5  = a(i3)

   movsd  8(%r10), %xmm8             # xmm8  = a(i0+1)
   movsd  8(%r11), %xmm9             # xmm9  = a(i1+1)
   movsd  8(%r12), %xmm10            # xmm10 = a(i2+1)
   movsd  8(%r13), %xmm11            # xmm11 = a(i3+1)

   addq    $16, %rdx                 # rdx   = &trigs(2)
   addq    $32, %r8                  # r8    = &trigs(4)
   addq    $48, %r9                  # r9    = &trigs(6)
   addq    $64, %rsi                 # rsi   = &c

   movapd  %xmm2 , %xmm0             # xmm0  = a(i0)
   movapd  %xmm3 , %xmm1             # xmm1  = a(i1)
   addsd   %xmm4 , %xmm0             # xmm0  = a(i0) + a(i2)
   addsd   %xmm5 , %xmm1             # xmm1  = a(i1) + a(i3)
   subsd   %xmm4 , %xmm2             # xmm2  = a(i0) - a(i2)
   subsd   %xmm5 , %xmm3             # xmm3  = a(i1) - a(i3)

   movapd  %xmm8 , %xmm6             # xmm6  = a(i0+1)
   movapd  %xmm9 , %xmm7             # xmm7  = a(i1+1)
   addsd   %xmm10, %xmm6             # xmm6  = a(i0+1) + a(i2+1)
   addsd   %xmm11, %xmm7             # xmm7  = a(i1+1) + a(i3+1)
   subsd   %xmm10, %xmm8             # xmm8  = a(i0+1) - a(i2+1)
   subsd   %xmm11, %xmm9             # xmm9  = a(i1+1) - a(i3+1)

   movapd  %xmm0 , %xmm4             # xmm4  = a0p2
   addsd   %xmm1 , %xmm0             # xmm0  = a0p2 + a1p3
   subsd   %xmm1 , %xmm4             # xmm4  = a0p2 - a1p3

   movapd  %xmm8 , %xmm5             # xmm5  = a4m6
   addsd   %xmm9 , %xmm8             # xmm8  = a4m6 + a5m7
   subsd   %xmm9 , %xmm5             # xmm5  = a4m6 - a5m7

   movsd   %xmm0 , -32(%rsi)         # c(0)
   movsd   %xmm8 ,    (%rsi)         # c(4)

   movapd  %xmm2 , %xmm0             # xmm0  = a0m2
   addsd   %xmm7 , %xmm2             # xmm2  = a0m2 + a5p7
   subsd   %xmm7 , %xmm0             # xmm0  = a0m2 - a5p7

   movapd  %xmm6 , %xmm1             # xmm1  = a4p6
   addsd   %xmm3 , %xmm1             # xmm1  = a4p6 + a1m3
   subsd   %xmm3 , %xmm6             # xmm6  = a4p6 - a1m3

   movsd    (%r8), %xmm7             # xmm7  = c2
   movsd    (%r9), %xmm15            # xmm15 = c3
   movsd   8(%r8), %xmm9             # xmm9  = s2
   movsd   8(%r9), %xmm12            # xmm12 = s3

   movapd  %xmm4 , %xmm3             # a0p2m1p3
   mulsd   %xmm7 , %xmm4             # (c2,c3) *
   mulsd   %xmm9 , %xmm3             # (s2,s3) *

   movapd  %xmm2 , %xmm10            # a0m2p5p7
   mulsd   %xmm15, %xmm2             # (c2,c3) *
   mulsd   %xmm12, %xmm10            # (s2,s3) *

   movapd  %xmm5 , %xmm11            # a4m6m5m7
   mulsd   %xmm7 , %xmm5             # (c2,c3) *
   mulsd   %xmm9 , %xmm11            # (s2,s3) *

   movapd  %xmm6 , %xmm14            # a4p6m1m3
   mulsd   %xmm15, %xmm6             # (c2,c3) *
   mulsd   %xmm12, %xmm14            # (s2,s3) *

   subsd   %xmm11, %xmm4             # xmm4  = c(2)
   addsd   %xmm5 , %xmm3             # xmm3  = c(6)

   subsd   %xmm14, %xmm2             # xmm4  = c(3)
   addsd   %xmm6 , %xmm10            # xmm3  = c(7)

   movsd   %xmm4 ,-16(%rsi)          # c(2)
   movsd   %xmm3 , 16(%rsi)          # c(6)
   movsd   %xmm2 , -8(%rsi)          # c(3)
   movsd   %xmm10, 24(%rsi)          # c(7)

   movsd    (%rdx), %xmm7            # xmm7  = c1
   movsd   8(%rdx), %xmm9            # xmm9  = s1

   movapd  %xmm1 , %xmm2             # xmm2  = a4p6p1m3
   mulsd   %xmm7 , %xmm1             # (c1) *
   mulsd   %xmm9 , %xmm2             # (s1) *
   movapd  %xmm0 , %xmm3             # xmm3  = a0m2m5p7
   mulsd   %xmm7 , %xmm0             # (c1) *
   mulsd   %xmm9 , %xmm3             # (s1) *

   subsd   %xmm2 , %xmm0             # xmm0  = c(1)
   addsd   %xmm3 , %xmm1             # xmm1  = c(5)   

   movsd   %xmm0 , -24(%rsi)         # c(1)
   movsd   %xmm1 ,   8(%rsi)         # c(5)

   addq    $16, %r10                 # i0 = i0 + 2
   addq    $16, %r11                 # i1 = i1 + 2
   subq    $16, %r12                 # i2 = i2 - 2
   subq    $16, %r13                 # i3 = i3 - 2
   dec     %r14
   jnz     LOOP_I_IFFT4S

   movsd   (%r10), %xmm5             # xmm5   = a(i0)
   movsd   (%r11), %xmm8             # xmm8   = a(i1)
   movsd  8(%r10), %xmm3             # xmm3   = a(i0+1)
   movsd  8(%r11), %xmm11            # xmm11  = a(i1+1)

   movapd  %xmm5 , %xmm4             # xmm4  = a(i0)
   addsd   %xmm8 , %xmm4             # xmm4  = a(i0) + a(i1)
   subsd   %xmm8 , %xmm5             # xmm5  = a(i0) - a(i1)

   movapd  %xmm3 , %xmm9             # xmm6  = a(i0+1)
   addsd   %xmm11, %xmm9             # xmm6  = a(i0+1) + a(i1+1)
   subsd   %xmm3, %xmm11             # xmm11 = a(i1+1) - a(i0+1)

   movapd  %xmm5 , %xmm8             # xmm8  = a(i0) - a(i1)
   movsd   %xmm4 , 32(%rsi)          # c(n-3)
   subsd   %xmm9, %xmm8              # xmm8  = 

   movq  $0x3FE6A09E667F3BCD, %rdx   # rdx   = sin45 = qrt(0.5)
   movd    %rdx , %xmm4              # xmm4  = sin45
   addsd   %xmm9, %xmm5              # xmm9  = 
   mulsd   %xmm4, %xmm8              # xmm8  = c(n-2)
   mulsd   %xmm4, %xmm5              # xmm5  = c(n)
   movsd   %xmm8 , 40(%rsi)          # c(n-2)

   movq  $0x8000000000000000, %rdx   # rdx   = sign bit
   movd    %rdx , %xmm8              # xmm8  = sign bit

   movsd   %xmm11, 48(%rsi)          # c(n-1)
   xorpd   %xmm8 , %xmm5             # negate
   movsd   %xmm5 , 56(%rsi)          # c(n)

   addq    %rcx, %rdi                # advance &a
   addq    $64 , %rsi                # advance &c
   popq    %r14                      # restore loop count
   popq    %rdx                      # restore trigs
   movq    %rdx, %r8                 # r8   = trigs
   movq    %rdx, %r9                 # r9   = trigs

   dec     %r15
   jnz     LOOP_O_IFFT4S

   popq    %rbx
   popq    %r12
   popq    %r13
   popq    %r14
   popq    %r15
   movq    %rbp, %rsp
   popq    %rbp
   ret


# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 01-Jun-2015
# ============================================
# call ifft4m(a  ,c  ,trigs,n  ,la ,lot)
#            (rdi,rsi,rdx  ,rcx,r8 ,r9 )
#
# rdi: a(n,lot)  : real(8) source array
# rsi: c(n,lot)  : real(8) target array
# rdx: trigs(n)  : sine/cosine values
# rcx: n         : 1st. dimension
# r8 : la        : factor
# r9 : lot       : # of lines
# --------------------------------------------

   .globl _ifft4m_
_ifft4m_:

   pushq   %r15
   pushq   %r14
   pushq   %r13
   pushq   %r12
   pushq   %rbp
   pushq   %rbx

   subq   $64, %rsp                  # reserve memory
   movl   (%r9), %r9d                # r9d = lot
   movq   %r9 , (%rsp)               # store lot
   movq   %rdx, %r13                 # &trigs       PERMANENT
   movl   (%rcx), %r11d              # r11 = n      PERMANENT
   movq   %r11, %rax                 # rax = n
   movq   %r11, %r15                 # r15 = n
   shrq   %r15                       # r15 = i5 = n / 2
   movl   (%r8), %r8d                # r8  = la     PERMANENT
   leaq   (,%r8,8), %r12             # r12 = la * 8 PERMANENT
   leaq   (%r8,%r8,2), %r14          # r14 = la * 3
   movq   %r8 , %rbx                 # rbx = la
   shrq   %rbx                       # rbx = la / 2
   movq   %r11, %r10                 # r10 = n
   subq   %r8 , %r10                 # r10 = i2 = n - la
   movq   %r15, %r9                  # r9  = n / 2
   subq   %r8 , %r9                  # r9  = i1 = n / 2 - la
   movq   %r9 ,  8(%rsp)             # i1
   movq   %r10, 16(%rsp)             # i2
   movq   %r15, 24(%rsp)             # i5
   movq   %r14, 32(%rsp)             # la * 3
   movq   %rdi, 40(%rsp)             # &a
   movq   %rsi, 48(%rsp)             # &c
   movq   %rbx, 56(%rsp)             # loop count la

LOOP_LOT_IFFT4M:
   movq    8(%rsp), %r9              # i1
   movq   16(%rsp), %r10             # i2
   movq   24(%rsp), %r15             # i5
   movq   32(%rsp), %r14             # la * 3
   movq   40(%rsp), %rdi             # &a
   movq   48(%rsp), %rsi             # &c
   movq   56(%rsp), %rbx             # loop count la

LOOP_1_I4FFTM:
   movapd  (%rdi)       , %xmm0      # xmm0  = a(i0)
   movapd  (%rdi,%r9 ,8), %xmm1      # xmm1  = a(i1)
   movapd  (%rdi,%r10,8), %xmm2      # xmm2  = a(i2)
   movapd  (%rdi,%r15,8), %xmm5      # xmm2  = a(i5)
   movapd  %xmm0, %xmm9              # xmm9  = a(i0)
   addpd   %xmm2, %xmm9              # xmm9  = a(i0) + a(i2)
   subpd   %xmm2, %xmm0              # xmm0  = a(i0) - a(i2)
   movapd  %xmm9, %xmm4              # xmm4  = a(i0) + a(i2)
   movapd  %xmm0, %xmm8              # xmm8  = a(i0) - a(i2)
   addpd   %xmm1, %xmm4              # xmm4  = a(i0) + a(i2) + a(i1)
   movapd  %xmm4, (%rsi)             # store   c(i)
   subpd   %xmm1, %xmm9              # xmm9  = a(i0) + a(i2) - a(i1)
   movapd  %xmm9, (%rsi,%r12,2)      # store   c(la*2+i)
   addpd   %xmm5, %xmm8              # xmm8  = a(i0) - a(i2) + a(i5)
   movapd  %xmm8, (%rsi,%r14,8)      # store   c(la*3+i)
   subpd   %xmm5, %xmm0              # xmm0  = a(i0) - a(i2) - a(i5)
   movapd  %xmm0, (%rsi,%r12)        # store   c(la+i)
   addq    $16,%rdi                  # advance &a
   addq    $16,%rsi                  # advance &c
   dec     %rbx
   jnz     LOOP_1_I4FFTM             # loop

   movq   %r11, %rbx                 # rbx = n
   leaq   (,%r8,4), %rax             # rax = la * 4
   subq   %rax, %rbx                 # rbx = i2 = n - la * 4
   movq   %rbx, %rdx                 # rdx = i2
   subq   %r15, %rdx                 # rdx = i3 = i2 - n / 2

   leaq   (%rdi,%r15,8), %r9         # r9    = &a(i1) : i1 = i0 + n / 2
   leaq   (%rdi,%rbx,8), %r10        # r10   = &a(i2) : i2 = n - 3 * la
   leaq   (%rdi,%rdx,8), %r15        # r15   = &a(i3) : i3 = n - 3 * la - n / 2
   leaq   (%r8 ,%r8 ,2), %rax        # rax   = la * 3 : rsi = &c(la)
   leaq   (%rsi,%rax,8), %rsi        # rsi   = &c(j0) : j0 = la * 4
   leaq   (%r8 ,%r8), %rbx           # rbx   = la * 2

   movq   %r11, %rax                 # rax = n
   shrq   $3, %rax                   # rax = n / 8
   movl   $0, %edx                   # rdx = 0
   divl   %r8d                       # rax = n / 8 / la
   dec    %rax                       # rax = n / 8 / la - 1
   movq   %rax, %rcx                 # store loop count
   jz     L0_I4FFTM                  # jump on zero loop

LOOP_O_I4FFTM:
   movddup  (%r13,%rbx,8), %xmm5     # xmm5  = c1 = trigs(kb)
   movddup 8(%r13,%rbx,8), %xmm6     # xmm6  = s1

   leaq     (%rbx,%rbx), %rax        # rax   = kc
   movddup  (%r13,%rax,8), %xmm7     # xmm7  = c2
   movddup 8(%r13,%rax,8), %xmm8     # xmm8  = s2

   leaq     (%rax,%rbx), %rax        # rax   = kd
   movddup  (%r13,%rax,8), %xmm9     # xmm9  = c3
   movddup 8(%r13,%rax,8), %xmm10    # xmm10 = s3

   movq   %r8 , %rax                 # rax   = loop count = la

LOOP_I_IFFT4M:
   movapd  (%rdi), %xmm1             # a(i0)
   movapd  (%r10), %xmm2             # a(i2)
   movapd  %xmm1 , %xmm15
   addpd   %xmm2 , %xmm15            # a(i0) + a(i2)
   subpd   %xmm2 , %xmm1             # a(i0) - a(i2)

   movapd  (%r9 ), %xmm11            # a(i1)
   movapd  (%r15), %xmm0             # a(i3)
   movapd  %xmm11, %xmm14
   addpd   %xmm0 , %xmm14            # a(i1) + a(i3)
   subpd   %xmm0 , %xmm11            # a(i1) - a(i3)

   movapd  (%rdi,%r12), %xmm13       # a(i4)
   movapd  (%r10,%r12), %xmm2        # a(i6)
   movapd  %xmm13, %xmm0
   addpd   %xmm2 , %xmm0             # a(i4) + a(i6)
   subpd   %xmm2 , %xmm13            # a(i4) - a(i6)

   movapd  (%r9 ,%r12), %xmm12       # a(i5)
   movapd  (%r15,%r12), %xmm3        # a(i7)
   movapd  %xmm12, %xmm2
   addpd   %xmm3 , %xmm2             # a(i5) + a(i7)
   subpd   %xmm3 , %xmm12            # a(i5) - a(i7)

   movapd  %xmm15, %xmm4             # xmm4  = a0p2
   subpd   %xmm14, %xmm4             # xmm4  = a0p2 - a1p3
   movapd  %xmm13, %xmm3             # xmm3  = a4m6
   subpd   %xmm12, %xmm3             # xmm3  = a4m6 - a5m7
   addpd   %xmm15, %xmm14            # xmm14 = a0p2 + a1p3

   movapd  %xmm14, (%rsi)            # c(j0) = a0p2 + a1p3

   addpd   %xmm13, %xmm12

   movapd  %xmm12, (%rsi,%r12,4)     # c(j4) = a4m6 - a5m7

   movapd  %xmm7, %xmm12
   mulpd   %xmm4, %xmm12
   movapd  %xmm8, %xmm13
   mulpd   %xmm3, %xmm13
   subpd   %xmm13, %xmm12

   movapd  %xmm12, (%rsi,%r12,2)     # c(j2) = c2 * a0p2m1p3 - s2 * a4m6m5m7

   mulpd   %xmm8, %xmm4
   mulpd   %xmm7, %xmm3
   addpd   %xmm4, %xmm3

   leaq    (%rsi,%r12,4), %rdx
   movapd  %xmm3, (%rdx,%r12,2)      # c(j6) = s2 * a0p2m1p3 + c2 * a4m6m5m7

   movapd  %xmm1, %xmm4
   subpd   %xmm2, %xmm4
   movapd  %xmm0, %xmm3
   addpd   %xmm11, %xmm3
   movapd  %xmm4, %xmm12
   mulpd   %xmm5, %xmm12
   movapd  %xmm3, %xmm13
   mulpd   %xmm6, %xmm13
   subpd   %xmm13, %xmm12
   movapd  %xmm12, (%rsi,%r12)       # c(j1)

   mulpd   %xmm6, %xmm4
   mulpd   %xmm5, %xmm3
   addpd   %xmm4, %xmm3

   movapd   %xmm3, (%rdx,%r12)        # c(j5)

   addpd   %xmm2, %xmm1
   subpd   %xmm11, %xmm0
   movapd  %xmm1, %xmm2
   mulpd   %xmm9, %xmm2
   movapd  %xmm0, %xmm3
   mulpd   %xmm10, %xmm3
   subpd   %xmm3, %xmm2

   leaq    (%rsi,%r12,2), %rdx
   movapd  %xmm2, (%rdx,%r12)        # c(j3)

   mulpd   %xmm10, %xmm1
   mulpd   %xmm9, %xmm0
   addpd   %xmm1, %xmm0

   leaq    (%rdx,%r12,4), %rdx
   movapd  %xmm0, (%rdx,%r12)        # c(j7)

   addq   $16, %rdi                  # &a(i0)++
   addq   $16, %r10                  # &a(i2)++
   addq   $16, %r9                   # &a(i1)++
   addq   $16, %r15                  # &a(i3)++
   addq   $16, %rsi                  # &c(j0)++
   subl   $2,%eax
   jnz   LOOP_I_IFFT4M

   leaq   (%rdi,%r12), %rdi          # adjust &a(i0)
   leaq   (%r9 ,%r12), %r9           # adjust &a(i1)
   leaq   (%r8 ,%r8 ,2), %rax        # rax = la * 3
   shlq   $3, %rax                   # rax * 8
   subq   %rax, %r10                 # adjust &a(i2)
   subq   %rax, %r15                 # adjust &a(i3)
   subq   %r12, %rsi
   leaq   (%rsi,%r12,8), %rsi        # adjust &c(j0)
   leaq   (%rbx,%r8,2), %rbx         # rbx += la * 2
   dec    %rcx                       # loop counter
   jnz   LOOP_O_I4FFTM

L0_I4FFTM:
   movq  $0x3FE6A09E667F3BCD, %rax   # rax   = sin45 = qrt(0.5)
   movd    %rax , %xmm5              # xmm5  = (sin45,   0 )
   movddup %xmm5, %xmm5              # xmm5  = (sin45,sin45)
   movq  $0x8000000000000000, %rax   # rax   = sign bit
   movd    %rax , %xmm6              # xmm6  = sign bit
   movddup %xmm6, %xmm6              # xmm6  = sign bits
   movq    %r8 , %rax                # loop count = la

LOOP_2_I4FFTM:
   movapd  (%rdi), %xmm1             # a(i0)
   movapd  (%r9 ), %xmm2             # a(i1)
   movapd  %xmm1, %xmm0
   addpd   %xmm2, %xmm0
   movapd  %xmm0, (%rsi)             # c(j0)
   subpd   %xmm2, %xmm1
   movapd  (%rdi,%r12), %xmm4        # a(i0+la)
   movapd  (%r9 ,%r12), %xmm2        # a(i1+la)
   movapd  %xmm4, %xmm0
   addpd   %xmm2, %xmm0
   movapd  %xmm1, %xmm3
   subpd   %xmm0, %xmm3
   mulpd   %xmm5, %xmm3
   movapd  %xmm3, (%rsi,%r12)        # c(j1)
   subpd   %xmm4, %xmm2
   movapd  %xmm2, (%rsi,%r12,2)      # c(j2)
   addpd   %xmm1, %xmm0
   mulpd   %xmm5, %xmm0
   xorpd   %xmm6, %xmm0
   leaq    (%rsi,%r12,2), %r14
   movapd  %xmm0, (%r14,%r12)        # c(j3)
   addq   $16, %rdi                  # &a(i0)++
   addq   $16, %r9                   # &a(i1)++
   addq   $16, %rsi                  # &c(j0)++
   subl   $2, %eax
   jnz    LOOP_2_I4FFTM

   movq   40(%rsp), %rdi             # &a
   movq   48(%rsp), %rsi             # &c
   leaq   (%rdi,%r11,8), %rdi        # &a += n * 8
   leaq   (%rsi,%r11,8), %rsi        # &c += n * 8
   movq   %rdi, 40(%rsp)             # &a
   movq   %rsi, 48(%rsp)             # &c

   decl    (%rsp)                     # --lot
   jnz    LOOP_LOT_IFFT4M

   addq   $64, %rsp                   # restore %rsp
   popq   %rbx
   popq   %rbp
   popq   %r12
   popq   %r13
   popq   %r14
   popq   %r15
   ret

# Fast Double Precision Fourier Transformation
# ============================================
# E. Kirk - 19-Mar-2015
# ============================================
# call ifft8(a  ,b  ,n  ,lot)
#           (rdi,rsi,rdx,rcx)
#
# rdi: a(n,lot)  : real(8) array
# rsi: n         : 1st. dimension
# rdx: lot       : 2nd. dimension
# ---------------------------------------


.globl _ifft8e_                      # MAC OSX names
_ifft8e_:

.globl ifft8e_                       # Linux names
ifft8e_:

    push    %rbp                     # save frame pointer
    movq    %rsp, %rbp
    call    mcount                   # enable profiling
    push    %r15
    push    %r14
    push    %r13
    push    %r12

    movl    (%rcx), %ecx             # lot
    movl    (%rdx), %edx             # n
    movl    %edx, %eax               # n
    sarl    $3, %eax                 # n / 8
    cltq
    movq    %rax, %r14               # i1 = la
    movq    %r14, %r13
    addq    %r14, %r13               # i2
    movq    %r13, %r12
    addq    %r14, %r12               # i3
    movq    %r12, %r11
    addq    %r14, %r11               # i4
    movq    %r11, %r10
    addq    %r14, %r10               # i5
    movq    %r10, %r9 
    addq    %r14, %r9                # i6
    movq    %r9 , %r8 
    addq    %r14, %r8                # i7
    movabsq $4609047870845172684, %rax
    movd    %rax, %xmm12             # SQRT(2.0)
    movlhps %xmm12, %xmm12           # expand

LOOP_O_IFFT8:
    movq    %r14, %r15               # loop count = la

LOOP_I_IFFT8:

# a0p7 = a(i0) + a(i7)
# a1m5 = a(i1) - a(i5)
#---------------------
    movupd  (%rdi)       , %xmm0     # a(i0)
    movapd  %xmm0, %xmm8             # a(i0)
    movupd  (%rdi,%r8 ,8), %xmm7     # a(i7)
    addpd    %xmm7, %xmm0            # a(i0) + a(i7)
    subpd    %xmm7, %xmm8            # a(i0) - a(i7)

# a1p5 = a(i1) + a(i5)
# a1m5 = a(i1) - a(i5)
#---------------------
    movupd  (%rdi,%r14,8), %xmm1     # a(i1)
    movapd  %xmm1, %xmm9             # a(i1)
    movupd  (%rdi,%r10,8), %xmm5     # a(i5)
    addpd   %xmm5, %xmm1             # a(i1) + a(i5)
    subpd   %xmm5, %xmm9             # a(i1) - a(i5)

# a2p6 = a(i2) + a(i6)
# a2m6 = a(i2) - a(i6)
#---------------------
    movupd  (%rdi,%r13,8), %xmm2     # a(i2)
    movapd  %xmm2, %xmm10            # a(i2)
    movupd  (%rdi,%r9 ,8), %xmm6     # a(i6)
    addpd   %xmm6, %xmm2             # a(i2) + a(i6)
    subpd   %xmm6, %xmm10            # a(i2) - a(i6)

# a0p7p3 = a0p7 + a(i3)
# a0p7m3 = a0p7 - a(i3)
#----------------------
    movupd  (%rdi,%r12,8), %xmm3     # a(i3)
    movapd  %xmm0, %xmm6             # a0p7
    movupd  (%rdi,%r11,8), %xmm4     # a(i4)
    addpd   %xmm3, %xmm0             # a0p7 + a(i3)
    subpd   %xmm3, %xmm6             # a0p7 - a(i3)

# a0m7p4   = 2.0 * (a0m7 + a(i4))
# a0m7m4   = 2.0 * (a0m7 - a(i4))
#--------------------------------
    movapd  %xmm8 , %xmm11           # a0m7
    addpd   %xmm4 , %xmm8            # a0m7 + a(i4)
    subpd   %xmm4 , %xmm11           # a0m7 - a(i4)
    addpd   %xmm8 , %xmm8            # * 2
    addpd   %xmm11, %xmm11           # * 2
        
# a1m5p2p6 = SQRT2 * (a1m5 + a2p6)
# a1m5m2p6 = SQRT2 * (a1m5 - a2p6)
#---------------------------------
    movapd  %xmm9 , %xmm14            # a1m5
    addpd   %xmm2 , %xmm9             # a1m5 + a2p6
    subpd   %xmm2 , %xmm14            # a1m5 - a2p6
    mulpd   %xmm12, %xmm9             # * SQRT(2.0)
    mulpd   %xmm12, %xmm14            # * SQRT(2.0)

# a(i0)  = 2.0 * (a0p7p3 + a1p5)
# a(i4)  = 2.0 * (a0p7p3 - a1p5)
#-------------------------------
    movapd  %xmm0, %xmm4             # a0p7p3
    addpd   %xmm1, %xmm0             # a0p7p3 + a1p5
    subpd   %xmm1, %xmm4             # a0p7p3 - a1p5
    addpd   %xmm0, %xmm0             # * 2
    addpd   %xmm4, %xmm4             # * 2
    movupd  %xmm0, (%rsi)            # store a(i0)
    movupd  %xmm4, (%rsi,%r11,8)     # store a(i4)

# a(i6)  = 2.0 * (a0p7m3 + a2m6)
# a(i2)  = 2.0 * (a0p7m3 - a2m6)
#-------------------------------
    movapd  %xmm6 , %xmm2            # a0p7m3
    addpd   %xmm10, %xmm6            # a0p7m3 + a2m6
    subpd   %xmm10, %xmm2            # a0p7m3 - a2m6
    addpd   %xmm6 , %xmm6            # * 2
    addpd   %xmm2 , %xmm2            # * 2
    movupd  %xmm6 , (%rsi,%r9 ,8)    # store a(i6)
    movupd  %xmm2 , (%rsi,%r13,8)    # store a(i2)

# a(i1)  = a0m7m4 + a1m5m2p6
# a(i5)  = a0m7m4 - a1m5m2p6
#---------------------------
    movapd  %xmm11, %xmm5            # a0m7m4
    addpd   %xmm14, %xmm11           # a0m7m4 + a1m5m2p6
    subpd   %xmm14, %xmm5            # a0m7m4 - a1m5m2p6
    movupd  %xmm11, (%rsi,%r14,8)    # store a(i1)
    movupd  %xmm5 , (%rsi,%r10,8)    # store a(i5)

# a(i3)  = a0m7p4 - a1m5p2p6
# a(i7)  = a0m7p4 + a1m5p2p6
#---------------------------
    movapd  %xmm8, %xmm3             # a0m7p4
    addpd   %xmm9, %xmm8             # a0m7p4 - a1m5p2p6
    subpd   %xmm9, %xmm3             # a0m7p4 - a1m5p2p6
    movupd  %xmm8, (%rsi,%r8 ,8)     # store a(i7)
    movupd  %xmm3, (%rsi,%r12,8)     # store a(i3)

    addq    $16, %rdi                # shift base of a
    addq    $16, %rsi                # shift base of b
    dec     %r15
    dec     %r15
    jnz     LOOP_I_IFFT8

    movq    %rdx, %rax               # n
    subq    %r14, %rax               # n - la
    leaq    (%rdi,%rax,8), %rdi      # next line
    leaq    (%rsi,%rax,8), %rsi      # next line
    dec     %ecx                     # --lot
    jnz     LOOP_O_IFFT8

    pop     %r12
    pop     %r13
    pop     %r14
    pop     %r15
    leave
    ret                              # finito


# Fast Double Precision Matrix Transposition
# ==========================================
# E. Kirk - 26-Jun-2015
# ==========================================
# call fast_mtp(a  ,n  )
#              (rdi,rsi)
#
# rdi: a(n,n)    : real(8) matrix
# rsi: n         : dimension
# ------------------------------------------

   .globl _fast_mtp
_fast_mtp:
   .globl _fast_mtp_
_fast_mtp_:

   movl    (%rsi), %esi              # rsi = n

   movq    $8, %rax                  # (1,0) : 8
   leaq    (,%rsi,8), %r8            # (0,1) : 8 * n
   leaq    -1(%rsi), %rdx            # outer loop : n-1

LOOP_O_FAST_MTP:
   movq    %rax, %r10                # (i,j)
   movq    %r8 , %r11                # (j,i)
   movq    %rdx, %rcx                # inner loop count

LOOP_I_FAST_MTP:
   movsd   (%rdi,%r10), %xmm0        # a(i,j)
   movsd   (%rdi,%r11), %xmm1        # a(j,i)
   movsd   %xmm0, (%rdi,%r11)
   movsd   %xmm1, (%rdi,%r10)

   addq    $8, %r10                  # += 1
   leaq    (%r11,%rsi,8), %r11       # += n

   dec     %rcx                      # inner loop counter
   jnz     LOOP_I_FAST_MTP

   leaq   8(%rax,%rsi,8), %rax       # rax += n + 1
   leaq   8(%r8 ,%rsi,8), %r8        # r8  += n + 1
   
   dec     %rdx                      # outer loop counter
   jnz     LOOP_O_FAST_MTP

   ret


# Fast Double Precision Fourier Transposition
# ==========================================
# E. Kirk - 25-Jun-2015
# ==========================================
# call fast_ftp(a  ,b  ,n  )
#              (rdi,rsi,rdx)
#
# rdi: a(*,*)    : complex fourier coefficients
# rsi: b(n,n)    : square matrix
# rdx: n         : dimension
# ------------------------------------------
# rax: multi purpose
# rbx: multi purpose, loop length
# rcx: loop counter
# rdx: 3rd. parameter, division rest
# rsp: only chnaged by push / pop
# rbp: unused
# rdi: array pointer a
# rsi: array pointer b
# r8 : byte count for one row of matrix a
# r9 : index register for a(i,k)
# r10: column address for a(k,i)
# r11: n/2 - k (# of fill positions at end of row)
# r12: n/2
# r13: index for accessing a(k,i)
# r14: byte count for one row of matrix b (n * 8)
# r15: loop index rows
# ------------------------------------------
# xmm0 : 0.0
# xmm5 : 0.5 , 0.5

# source array example for T2
#             0           1           2
# -------------------------------------
#  0   R00    0    R10  I10    R20  I20
#  1   R01  I01    R11  I11    R21  I21
#  2   R02  I02    R12  I12    R22  I22
# -2   P02  Q02    P12  Q12    P22  Q22
# -1   P01  Q01    P11  Q11    P21  Q21
#
# target array example for N8
#        0    1      2    3      4    5      6    7
# -------------------------------------------------
#  0   X00    0    X10  Y10    X20  Y20      0    0
#  1     0    0      0    0      0    0      0    0
#  2   X02    0    X12  Y12    X22  Y22      0    0
#  3   X03    0    X13  Y13    X23  Y23      0    0
#  4   X04    0    X14  Y14    X24  Y24      0    0
#  5   X05    0    X15  Y15    X25  Y25      0    0
#  6     0    0      0    0      0    0      0    0
#  7     0    0      0    0      0    0      0    0

# (X00,  0) = (R00,  0)
# (X10,Y10) = (R01,I01)   (X20,Y20) = (R02,I02)
# (X02,  0) = (R10,  0)   (X03,  0) = (I10,  0)   ...
# (X12,Y12) = 0.5 * (R11,I11) +- (P11,Q11))
# (X13,Y13) = 0.5 * (P11,Q11) -+ (R11,I11))
# (X22,Y22) = 0.5 * (R12,I12) +- (P12,Q12))
# (X23,Y23) = 0.5 * (P12,Q12) -+ (R12,I12))


   .globl _fast_ftp
_fast_ftp:
   .globl _fast_ftp_
_fast_ftp_:


   push    %rbx
   push    %r12
   push    %r13
   push    %r14
   push    %r15

   xorpd   %xmm0, %xmm0              # xmm0 = 0
   movq    $0x3fe0000000000000, %rax # rax  = 0.5
   movd    %rax, %xmm5               # xmm5 = 0.5
   movddup %xmm5, %xmm5              # xmm5 = 0.5 , 0.5
   movl    (%rdx), %eax              # rax  = n
   movq    %rax, %r11                # r11  = n
   movq    %rax, %r14                # r14  = n
   shl     $3, %r14                  # r14  = n * 8 (byte count)
   movl    $3, %ebx                  # rbx  = 3
   xorq    %rdx, %rdx                # rdx  = 0
   divl    %ebx                      # rax  = n/3
   inc     %rax                      # rax  = k = (n/3) + 1
   sarq    %r11                      # r11  = n/2
   subq    %rax, %r11                # r11  = n/2 - k
   movq    %rax, %r8                 # r8   = k
   movq    %rax, %rbx                # rbx  = k
   dec     %rbx                      # wavenumber excluding 0
   shlq    $4, %r8                   # r8  = row increment [bytes]

# the imaginary part of mode (0,0) is zero by definition

   movsd   (%rdi), %xmm1             # xmm1 = a(0,0)
   movapd  %xmm1, (%rsi)             # b(0,0) = xmm1
   movapd  %xmm0, (%rsi,%r14)        # b(0,1) = 0.0
   mov     %r8, %r12                 # r12 = row
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   movq    %rax, %rcx                # rcx = k
   dec     %rcx

# transpose row 0 and set row 1 to zero

LOOP_1_FAST_FTP:
   movapd  (%rdi,%r12), %xmm1        # xmm1 = a(0,i)
   movapd  %xmm1, (%rsi)             # b(i,0) = xmm1
   movapd  %xmm0, (%rsi,%r14)        # b(i,1) = 0.0
   addq    %r8, %r12                 # r12+= row
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_1_FAST_FTP

# clear rest of rows

   movq    %r11, %rcx                # rcx = n/2 - k

LOOP_2_FAST_FTP:
   movapd  %xmm0, (%rsi)             # b(i,0) = 0.0
   movapd  %xmm0, (%rsi,%r14)        # b(i,1) = 0.0
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_2_FAST_FTP
   addq    %r14, %rsi                # advance one row

# compute column 0

   movq    $16, %r10                 # r10 = index j
   movq    %r8, %r13                 # r13 = rowsize
   sal     $1, %rax                  # rax = k * 2
   subl    $2, %eax                  # rax = k * 2 - 3
   mul     %r13d                     # rax = rowsize * (k-1)
   leaq    (%r10, %rax), %r9         # r9  = index k

   movq    %rbx, %r15                # loop count rows

LOOP_O_FAST_FTP:
   movq    %r10, %r12                # j
   movq    %r9 , %r13                # k
   movsd   (%rdi,%r12), %xmm1        # xmm1 = a(i  ,0)
   movsd   8(%rdi,%r12), %xmm2       # xmm1 = a(i+1,0)
   movapd  %xmm1, (%rsi)             # b(0,i  ) = xmm1
   movapd  %xmm2, (%rsi,%r14)        # b(0,i+1) = xmm2
   addq    $16, %rsi                 # &b += 16 bytes (one complex)

# compute fourier elemnt for wavenumber i,j (4 values)

   movq    %rbx, %rcx                # loop count = k
   addq    %r8, %r12                 # r12+= row
   

LOOP_I_FAST_FTP:
   movapd  (%rdi,%r12), %xmm1        # xmm1 = a(j,i)
   movapd  %xmm1, %xmm7              # xmm7 = a(j,i)
   movapd  (%rdi,%r13), %xmm6        # xmm6 = a(j,k)
   xorpd   %xmm2, %xmm2              # xmm2 = 0
   subpd   %xmm6, %xmm2              # xmm2 = -a(j,k)
   addsubpd %xmm2, %xmm1             # xmm1 = a(j,i) -+ a(j,k)
   mulpd   %xmm5, %xmm1              # xmm1 * 0.5
   movapd  %xmm1, (%rsi)             # b(i,j)
   addsubpd %xmm7, %xmm6             # xmm6 = a(j,k) -+ a(j,i)
   mulpd   %xmm5, %xmm6              # xmm6 * 0.5
   shufpd  $1, %xmm6, %xmm6          # swap quad words
   movapd  %xmm6, (%rsi,%r14)        # b(i,j+n)

   addq    %r8, %r12                 # r12 += row
   subq    %r8, %r13                 # r13 -= row
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_I_FAST_FTP

# fill rest of row

   movq    %r11, %rcx                # rcx = n/2 - k

LOOP_3_FAST_FTP:
   movapd  %xmm0, (%rsi)             # b(i,j  ) = 0
   movapd  %xmm0, (%rsi,%r14)        # b(i,j+1) = 0
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rcx
   jnz     LOOP_3_FAST_FTP

   addq    %r14,  %rsi               # advance row
   addq    $16, %r10                 # advance column
   addq    $16, %r9

   dec     %r15
   jnz     LOOP_O_FAST_FTP

# fill rest of target array

   movq    %r11, %rax                # rax = n/2 - k
   mul     %r14                      # rax = (n*8) * (n/2 - k)
   sar     $3, %rax                  # rax = n * (n/2 - k)
   
LOOP_4_FAST_FTP:
   movapd  %xmm0, (%rsi)             # b(i,j  ) = 0
   addq    $16, %rsi                 # &b += 16 bytes (one complex)
   dec     %rax
   jnz     LOOP_4_FAST_FTP

   pop     %r15
   pop     %r14
   pop     %r13
   pop     %r12
   pop     %rbx

   ret

