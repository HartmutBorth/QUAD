! =============
! MODULE FFTMOD
! =============

module fftmod
integer, parameter :: NRES = 9
integer :: nallowed(NRES) = [16,32,64,128,256,512,1024,2048,4096]

! T3    - N16   : 8-2
! T10   - N32   : 8-2-2
! T15   - N48   : 8-3-2
! T21   - N64   : 8-4-2
! T31   - N96   : 8-4-3
! T42   - N128  : 8-4-4
! T85   - N256  : 8-4-4-2
! T127  - N384  : 8-4-4-3
! T170  - N512  : 8-4-4-4
! T341  - N1024 : 8-4-4-4-2
! T682  - N2048 : 8-4-4-4-4
! T1365 - N4096 : 8-4-4-4-4-2

integer :: lastn = 0 ! last used value for n
integer :: lsize = 0 ! last used size for w
integer :: nbit  = 0 ! bit # set for value of n
real(8),allocatable :: trigs(:)
real(8),allocatable :: w(:)
end module fftmod

! *******************
! * fourier_to_grid *
! *******************

subroutine fourier_to_grid(fc,gp)
use quadmod
implicit none
real(8),intent(in ) :: fc(*)
real(8),intent(out) :: gp(*)

call fast_ftp(fc,gp,ngx)
call fc2gp(gp,ngy,nfx+1)
call fast_mtp(gp,ngx)
call fc2gp(gp,ngx,ngy)

return
end


! *******************
! * grid_to_fourier *
! *******************

subroutine grid_to_fourier(gp,fc)
use quadmod
implicit none
real(8),intent(in ) :: gp(ngx,ngy)
real(8),intent(out) :: fc(0:nfx,0:nfy)
real(8)             :: zz(ngx,ngy)
integer             :: j

zz(:,:) = gp(:,:)
call gp2fc(zz,ngx,ngy)
zz(:,:) = transpose(zz)
call gp2fc(zz,ngy,nfx+2)

! transpose, reorder and compute real to complex coefficients

fc(:,0) = zz(1,1:nfx+1)

fc(0,  1:nky     ) =  zz(3:nfy+2:2,1)
fc(0,nfy:nky+1:-1) =  zz(3:nfy+2:2,1)
fc(1,  1:nky     ) =  zz(4:nfy+2:2,1)
fc(1,nfy:nky+1:-1) = -zz(4:nfy+2:2,1)

do j = 2 , nfx , 2
   fc(j  ,  1:nky     ) = zz(3:nfy+2:2,j+1) - zz(4:nfy+2:2,j+2)
   fc(j+1,  1:nky     ) = zz(4:nfy+2:2,j+1) + zz(3:nfy+2:2,j+2)
   fc(j  ,nfy:nky+1:-1) = zz(3:nfy+2:2,j+1) + zz(4:nfy+2:2,j+2)
   fc(j+1,nfy:nky+1:-1) = zz(3:nfy+2:2,j+2) - zz(4:nfy+2:2,j+1)
enddo
 
return
end


! ================
! SUBROUTINE GP2FC
! ================

subroutine gp2fc(a,n,lot)
use fftmod
implicit none
real(8) :: a(n,lot)
integer :: n
integer :: lot

integer :: l,la

if (n /= lastn) then
   if (allocated(trigs)) deallocate(trigs)
   allocate(trigs(n))
   lastn = n
   call fftini(n,0)
endif

call dfft8(a,n,lot)
la = n / 8
do while (la >= 4)
   call dfft4(a,trigs,n,lot,la)
enddo

if (la == 3) then
   do l = 1 , lot
      call dfft3(a(1,l),trigs,n)
   enddo
endif

if (la == 2) then
   do l = 1 , lot
      call dfft2(a(1,l),trigs,n)
   enddo
endif
return
end subroutine gp2fc

! ================
! SUBROUTINE FC2GP
! ================

subroutine fc2gp(a,n,lot)
use fftmod
implicit none
real(8) :: a(n*lot)
integer :: n,lot,j,la
logical :: ld

if (n /= lastn .or. n * lot > lsize) call fftini(n,lot)

la = 2 + 2 * mod(nbit,2) ! start with la=2 or la=4
ld = .true.              ! start ifft4m with w -> a

if (la == 2) then
   call ifft2s(a,w,trigs,n,lot)
else
   call ifft4s(a,w,trigs,n,lot)
endif

do j = 1 , nbit / 2
   if (ld) then
      call ifft4m(w,a,trigs,n,la,lot)
   else
      call ifft4m(a,w,trigs,n,la,lot)
   endif
   la = la * 4
   ld = .not. ld
enddo

! make sure, that the final result is stored in array a

if (ld) then
   call ifft8e(w,a,n,lot)
else
   call ifft8e(a,a,n,lot) ! ifft8e allows target == source 
endif

return
end subroutine fc2gp

! =================
! SUBROUTINE FFTINI
! =================

subroutine fftini(n,lot)
use fftmod
implicit none
integer :: n,lot
logical labort

integer :: j,k,ibit
real(8) :: angle,del

! check for allowed values of n

labort = .true.
do j = 1 , NRES
   if (n == nallowed(j)) labort = .false.
enddo

if (labort) then
   write (*,*) '*** FFT does not support n = ',n,' ***'
   write (*,*) 'Following resolutions may be used:'
   write (*,*) '----------------------------------'
   do j = 1 , NRES
      write (*,1000) nallowed(j), nallowed(j)/2, nallowed(j)/3
   enddo
   stop
endif
 1000 format(' NLON=',I5,'  NLAT=',I5,'  NTRU=',I5)

! compute the bit # that is set in the integer n (power of 2)

ibit = n / 32 ! factors 8 and 4 set alreay
nbit = 0
do while (ibit > 0)
   ibit = ibit / 2
   nbit = nbit + 1
enddo  

! allocate trigs for new length n

if (allocated(trigs)) deallocate(trigs)
allocate(trigs(n))
lastn = n

! allocate w for new length n * lot

if (lot > 0) then
   lsize = n * lot
   if (allocated(w)) deallocate(w)
   allocate(w(lsize))
endif

lastn = n

! compute values for trigs

del = 4.0 * asin(1.0) / n
do k=0,n/2-1
  angle = k * del
  trigs(2*k+1) = cos(angle)
  trigs(2*k+2) = sin(angle)
enddo
return
end subroutine fftini

! ================
! SUBROUTINE DFFT2
! ================

subroutine dfft2(a,trigs,n)
implicit none
real(8) :: a(n)
real(8) :: c(n)
real(8) :: trigs(n)
integer :: n,ja,jb,i
real(8) :: c1,s1
real(8) :: a1p3,a3m1

c(1) = a(1) + a(2)
c(2) = 0.0

ja = 3
jb = n - 1

do i=3,n-5,4
   c1 = trigs(ja  )
   s1 = trigs(ja+1)
   a1p3 = c1 * a(i+1) + s1 * a(i+3)
   a3m1 = c1 * a(i+3) - s1 * a(i+1)
   c(ja  ) = a(i) + a1p3
   c(jb  ) = a(i) - a1p3
   c(ja+1) = a3m1 + a(i+2)
   c(jb+1) = a3m1 - a(i+2)
   ja = ja + 2
   jb = jb - 2
enddo

c(ja  ) =  a(n-1)
c(ja+1) = -a(n  )

a = c
return
end subroutine dfft2

! ================
! SUBROUTINE DFFT3
! ================

subroutine dfft3(a,trigs,n)
implicit none
real(8), parameter :: SIN60 = 0.866025403784438D0
real(8)            :: a(n),c(n),trigs(n)
integer :: n
integer :: i
integer :: ja,jb,jc
real(8) :: a1,a2,a3
real(8) :: b1,b2,b3
real(8) :: c1,c2
real(8) :: s1,s2

ja = 1              !  1
jb = 2 * (n/3)  + 1 ! 65
jc = jb             ! 65

c(ja  ) = a(1) + a(2) + a(3)
c(ja+1) = 0.0
c(jb  ) = a(1) - 0.5 * (a(2) + a(3))
c(jb+1) =      SIN60 * (a(3) - a(2))

ja = 3         !  3, 5, 7, ... ,31
jb = jb + 2    ! 67,69,71, ... ,95
jc = jc - 2    ! 63,61,59, ... ,35

do i = 4 , n-8 , 6 ! 88
   c1 = trigs(ja  )
   s1 = trigs(ja+1)
   c2 = trigs(ja+ja-1)
   s2 = trigs(ja+ja  )
   a1 = (c1*a(i+1)+s1*a(i+4))+(c2*a(i+2)+s2*a(i+5))
   b1 = (c1*a(i+4)-s1*a(i+1))+(c2*a(i+5)-s2*a(i+2))
   a2 = a(i  ) - 0.5 * a1
   b2 = a(i+3) - 0.5 * b1
   a3 = SIN60*((c1*a(i+1)+s1*a(i+4))-(c2*a(i+2)+s2*a(i+5)))
   b3 = SIN60*((c1*a(i+4)-s1*A(i+1))-(c2*a(i+5)-s2*a(i+2)))
   c(ja  ) = a(i  ) + a1
   c(ja+1) = a(i+3) + b1
   c(jb  ) = a2 + b3
   c(jb+1) = b2 - a3
   c(jc  ) = a2 - b3
   c(jc+1) =-b2 - a3
   ja = ja + 2
   jb = jb + 2
   jc = jc - 2
enddo

if (ja <= jc) then ! ja=33  jc=33
   c(ja  ) = a(n-2) + 0.5 * (a(n-1) - a(n)) ! 33
   c(ja+1) =       -SIN60 * (a(n-1) + a(n)) ! 34
endif
a(:) = c(:)
return
end subroutine dfft3

! ================
! SUBROUTINE DFFT4
! ================

subroutine dfft4(a,trigs,n,lot,la)
implicit none
real(8) :: a(n,lot)
real(8) :: c(n,lot)
real(8) :: trigs(n)
integer :: n,lot,la

integer :: i,j,k,l
integer :: ibase,jink,kb,kc,kd
integer :: i1,i2,i3,i4,i5,i6,i7
integer :: j0,j1,j2,j3,j4,j5,j6,j7
real(8) :: a0,a1,a2,a3
real(8) :: b0,b1,b2,b3
real(8) :: a0p2,a1p3,a1m3
real(8) :: a1p5,a2p6,a3p7
real(8) :: a5m1,a6m2,a7m3
real(8) :: c1,c2,c3
real(8) :: s1,s2,s3
real(8) :: sin45

la = la / 4

i1 = la
i2 = la + i1
i3 = la + i2
i4 = la + i3
i5 = la + i4
i6 = la + i5
i7 = la + i6

j1 = n/2 - la
j2 = n - la
j3 = j1
j5 = j1 + la

do i=1,la
   do l=1,lot
   a0p2 = a(i   ,l) + a(i2+i,l)
   a1p3 = a(i1+i,l) + a(i3+i,l)
   c(   i,l) = a0p2 + a1p3
   c(j2+i,l) = a0p2 - a1p3
   c(j1+i,l) = a(   i,l) - a(i2+i,l)
   c(j5+i,l) = a(i3+i,l) - a(i1+i,l)
   enddo
enddo

jink = 2 * la
j0 = la
j1 = j1 + jink
j2 = j2 - jink
j3 = j3 - jink
j4 = j0 + la
j5 = j1 + la
j6 = j2 + la
j7 = j3 + la

ibase=4*la

do 450 k=la,(n-4)/8,la
   kb=k+k
   kc=kb+kb
   kd=kc+kb
   c1=trigs(kb+1)
   s1=trigs(kb+2)
   c2=trigs(kc+1)
   s2=trigs(kc+2)
   c3=trigs(kd+1)
   s3=trigs(kd+2)

   i=ibase+1
   do j=1,la
      do l=1,lot
      a1p5 = c1 * a(i1+i,l) + s1 * a(i5+i,l)
      a2p6 = c2 * a(i2+i,l) + s2 * a(i6+i,l)
      a3p7 = c3 * a(i3+i,l) + s3 * a(i7+i,l)
      a5m1 = c1 * a(i5+i,l) - s1 * a(i1+i,l)
      a6m2 = c2 * a(i6+i,l) - s2 * a(i2+i,l)
      a7m3 = c3 * a(i7+i,l) - s3 * a(i3+i,l)
      a0 = a(i,l) + a2p6
      a2 = a(i,l) - a2p6
      a1 = a1p5 + a3p7
      a3 = a3p7 - a1p5
      b0 = a(i4+i,l) + a6m2
      b2 = a(i4+i,l) - a6m2
      b1 = a5m1 + a7m3
      b3 = a5m1 - a7m3
      c(j0+j,l) = a0+a1
      c(j2+j,l) = a0-a1
      c(j4+j,l) = b0+b1
      c(j6+j,l) = b1-b0
      c(j1+j,l) = a2+b3
      c(j3+j,l) = a2-b3
      c(j5+j,l) = a3+b2
      c(j7+j,l) = a3-b2
      enddo
      i=i+1
   enddo

  ibase=ibase+8*la
  j0 = j0 + jink
  j1 = j1 + jink
  j2 = j2 - jink
  j3 = j3 - jink
  j4 = j0 + la
  j5 = j1 + la
  j6 = j2 + la
  j7 = j3 + la
450   continue
if (j1 <= j2) then
   sin45=sqrt(0.5)
   i=ibase+1
   do j=1,la
      do l=1,lot
      a1p3 = sin45 * (a(i1+i,l) + a(i3+i,l))
      a1m3 = sin45 * (a(i1+i,l) - a(i3+i,l))
      c(j0+j,l) =  a(   i,l) + a1m3
      c(j1+j,l) =  a(   i,l) - a1m3
      c(j4+j,l) = -a(i2+i,l) - a1p3
      c(j5+j,l) =  a(i2+i,l) - a1p3
      enddo
      i=i+1
   enddo
endif
if (la == 1) then
   do l=1,lot
   a(1,l) = c(1,l)
   a(2,l) = 0.0
   a(3:n,l) = c(2:n-1,l)
   enddo
else
   a = c
endif
return
end subroutine dfft4

! ================
! SUBROUTINE DFFT8
! ================

subroutine dfft8(a,n,lot)
implicit none
real(8) :: a(n*lot)
integer :: n
integer :: lot

integer :: i,k,la
integer :: i0,i1,i2,i3,i4,i5,i6,i7

real(8) :: z,zsin45
real(8) :: a0p4,a1p5,a2p6,a3p7,a5m1,a7m3,a0m4,a6m2
real(8) :: a0p4p2p6,a1p5p3p7,a7m3p5m1,a7m3m5m1

la = n / 8
z  = 1.0 / n
zsin45 = z * sqrt(0.5_8)

do i=0,la*lot-1
   i0 = (i/la) * n + mod(i,la) + 1
   i1 = i0 + la
   i2 = i1 + la
   i3 = i2 + la
   i4 = i3 + la
   i5 = i4 + la
   i6 = i5 + la
   i7 = i6 + la

   a0p4 =  a(i0) + a(i4)
   a1p5 =  a(i1) + a(i5)
   a2p6 =  a(i2) + a(i6)
   a3p7 =  a(i3) + a(i7)
   a5m1 =  a(i5) - a(i1)
   a7m3 =  a(i7) - a(i3)
   a0m4 = (a(i0) - a(i4)) * z
   a6m2 = (a(i6) - a(i2)) * z

   a0p4p2p6 = a0p4 + a2p6
   a1p5p3p7 = a1p5 + a3p7
   a7m3p5m1 = (a7m3 + a5m1) * zsin45
   a7m3m5m1 = (a7m3 - a5m1) * zsin45

   a(i0) = z * (a0p4p2p6 + a1p5p3p7)
   a(i7) = z * (a0p4p2p6 - a1p5p3p7)
   a(i3) = z * (a0p4 - a2p6)
   a(i4) = z * (a3p7 - a1p5)
   a(i1) = a0m4 + a7m3m5m1
   a(i5) = a0m4 - a7m3m5m1
   a(i2) = a7m3p5m1 + a6m2
   a(i6) = a7m3p5m1 - a6m2
enddo
return
end subroutine dfft8

