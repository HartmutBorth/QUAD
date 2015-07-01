module quadmod

!===========================================!
!                                           !
!       Fluid Simulator (FluSim)            !
!                                           !
!      Version: 1.0    March 2015           !
!                                           !
! (Based on a code of Annalisa Bracco)      !
!                                           !
!===========================================!
! Theoretical Meteorology - KlimaCampus     !
!       University of Hambuurg              !
!===========================================!
! Valerio Lucarini                          !
!                                           !
! Edilbert Kirk - Hartmut Borth             !
!===========================================!
! The latest version of the code and the    !
! the User's Guide can be downloaded from   !
!                                           !
! github   .....                            !
!===========================================!


!===========================================!
! For more details on parameters, variables !
! and default values, see the User's Guide. !
!===========================================!

!=====================!
! Basic control flags !
!=====================!


!===========================================!
! Basic physical and mathematical constants !
!===========================================!
real(8), parameter :: pi    = 4.0d0 * atan(1.0d0)
real(8), parameter :: twopi = pi+pi

complex(8), parameter :: ci = (0.0,1.0)    ! complex unit

!===============================!
! Basic input/output parameters !
!===============================!

!--- i/o units
integer, parameter :: nunl     = 13  ! i/o-unit for namelist  
integer, parameter :: nuini    = 10  
integer, parameter :: nuphys   = 20
integer, parameter :: nunserv  = 30
integer, parameter :: nunstat  = 31

!--- i/o file names
character (256) :: quad_namelist       = "quad_namelist"

!--- header for output of service format
integer :: ihead(8)


!===================================!
! Basic parameters of model physics !
!===================================!

!--- size of fluid domain
real(8) :: lx = twopi   ! length [m]
real(8) :: ly = twopi   ! width [m]

!--- parameters of evolution equations
real(8) :: alpha = 0.0      ! modification parameter for finite Rossby radius 
                            ! alpha =  1/L_R**2 [1/m**2]
real(8) :: beta  = 0.0      ! strength of ambient vorticity gradient [1/m*s]

!--- Laplacian viscosity (dissipation on small scales)
real(8) :: sig   =  5.d-6    ! time-scale of hyperfriction [1/s]
real(8) :: klow  =  64       ! low cutoff wave number  [1/m]       
real(8) :: psig  =  2.d0     ! power of Laplacian modelling dissipation

!--- Laplacian friction (dissipation on large scales)
real(8) :: lam   =  0.65d-3  ! time-scale of hypofriction [1/s]
real(8) :: kmax  =  2        ! maximum resolved wave number [1/m]
real(8) :: plam  = -4.d0     ! power of Laplacian modelling dissipation
                             ! plam <= 0


! parafor.h   Forcing characteristics
! kk0 is the forcing shell, (forcing on all K with modulus ]kk0-dk,kk0+dk[ )
! The amplitude of the forcing is fixed (AMPFORCING is the energy of the forcing).
! The phases are random, and are described
! as a Markov process, whose decorrelation time is TF (adimensional units).

real(8) :: ampforc = 0.015   ! amplitude of forcing
real(8) :: tforc   = 0.001   ! memory time scale of forcing [s]

! commfor.h   Common della fase del forcing

integer :: idum
integer :: in(1000,2)
integer :: nk
integer :: itau

real(8) :: ampcoeff

complex(8) :: phi

integer :: nforc  = 0       ! forcing switch 0/1/2  no/Markov/white noise     
real(8) :: kfmin  = 4.0+1.0 ! min radius = sqrt(kfmin) of spectr. forcing
real(8) :: kfmax  = 4.0-1.0 ! max radius = sqrt(kfmax) of spectr. forcing
!============================================


!======================================!
! Basic parameters of numerical scheme !
!======================================!

!--- grid in physical space
integer :: ngx = 64    ! number of grid points (x-direction)
integer :: ngy = 64    ! number of grid points (y-direction)

integer :: nxy  = 4096                      ! ngx*ngy
real(8) :: rnx  = 1.56250000000000000E-002  ! 1/ngx
real(8) :: rny  = 1.56250000000000000E-002  ! 1/ngy
real(8) :: rnxy = 2.44140625000000000E-004  ! 1/nxy


!--- grid in Fourier space
integer :: nkx  = 21  ! ngx/3 (max. wave number in x-direction)
integer :: nky  = 21  ! ngy/3 (max. wave number in y-direction)
integer :: nfx  = 43  ! 2*nkx+1 (x-dimension in fourier domain)
integer :: nfy  = 42  ! 2*nky   (y-dimension in fourier domain)


!--- time integration
integer :: nstop = 10000  ! last time step of integration
integer :: nout  = 100    ! time steps between data output (-1 = no output)
integer :: ndiag = 10     ! time steps between diagnostic output (-1 no output)
integer :: ncfl  = 10     ! time steps between cfl-check (-1 no cfl check)
 
real(8) :: dt    = 1.d-3  ! length of time step [s]

integer :: nstep = 0      ! first time step counter


!--- variables in physical (grid) space
real(8), allocatable :: gvo(:,:)  ! vorticity            [1/s]
real(8), allocatable :: gpsi(:,:) ! stream function      [m^2/s]
real(8), allocatable :: gu(:,:)   ! velocity x-direction [m/s]
real(8), allocatable :: gv(:,:)   ! velocity y-direction [m/s]

real(8), allocatable :: gp1(:,:) ! u or u*zeta
real(8), allocatable :: gp2(:,:) ! v or v*zeta


!--- variables in spectral (Fourier) space 
real(8), allocatable :: fpsi(:,:) ! stream function
real(8), allocatable :: fu(:,:)   ! velocity x-direction
real(8), allocatable :: fv(:,:)   ! velocity y-direction


real(8), allocatable :: fp1(:,:) ! u or u*zeta
real(8), allocatable :: fp2(:,:) ! v or v*zeta

complex(8), allocatable :: cvo(:,:)  ! vorticity
complex(8), allocatable :: cnl0(:,:) ! Jacobian at time 0
complex(8), allocatable :: cnl1(:,:) ! Jacobian at time -1
complex(8), allocatable :: cnl2(:,:) ! Jaconbian at time -2


!============================================
! barotropic model with finite Rossby radius f = (L/R)^2
!real(8) :: f = 0.0

! strength of ambient vorticity gradient
!real(8) :: beta = 0.0

! parafor.h   Forcing characteristics
! kk0 is the forcing shell, (forcing on all K with modulus ]kk0-dk,kk0+dk[ )
! The amplitude of the forcing is fixed (AMPFORCING is the energy of the forcing).
! The phases are random, and are described
! as a Markov process, whose decorrelation time is TF (adimensional units).

!real(8) :: ampforc = 0.015   ! amplitude of forcing
!real(8) :: tforc   = 0.001   !  

! commfor.h   Common della fase del forcing

!integer :: idum
!integer :: in(1000,2)
!integer :: nk
!integer :: itau

!real(8) :: ampcoeff

!complex(8) :: phi

!integer :: nforc  = 0       ! forcing switch 0/1/2  no/Markov/white noise     
!real(8) :: kfmin  = 4.0+1.0 ! min radius = sqrt(kfmin) of spectr. forcing
!real(8) :: kfmax  = 4.0-1.0 ! max radius = sqrt(kfmax) of spectr. forcing
!============================================


! fft parameters
integer   , allocatable :: ki(:)
integer   , allocatable :: kj(:)
real(8)   , allocatable :: kx  (:), ky  (:)
real(8)   , allocatable :: kx2 (:), ky2 (:)
real(8)   , allocatable :: pkx2(:), pky2(:)
complex(8), allocatable :: cikx(:), ciky(:)
complex(8), allocatable :: cli(:,:)

end module quadmod

! ******************************
! * SUBROUTINE ALLOCATE_ARRAYS *
! ******************************

subroutine allocate_arrays
use quadmod
implicit none

integer :: j

nxy  = ngx * ngy       ! # of gridpoints
print *, "nxy = ", nxy
nkx  = ngx / 3       ! max wavenumber in x
print *, "nkx = ", nkx
nky  = ngy / 3       ! max wavenumber in y
print *, "nky = ", nky
nfx  = nkx * 2 + 1   ! x dimension in fourier domain
print *, "nfx = ", nfx
nfy  = nky * 2       ! y dimension in fourier domain
print *, "nfy = ", nfy
rnx  = 1.0d0 / ngx   ! reverse of ngx
print *, "rnx = ", rnx
rny  = 1.0d0 / ngy   ! reverse of ngy
print *, "rny = ", rny
rnxy = 1.0d0 / nxy   ! reverse of gridpoint number
print *, "rnxy = ", rnxy


allocate(ki(0:nkx)) ; ki(:) = [(j,j=0,nkx)]
allocate(kj(0:nfy)) ; kj(:) = [(j,j=0,nky),(j,j=-nky,-1)]

allocate(gp1(1:ngx,1:ngy)) ; gp1(:,:) = 0.0 ! u or u*zeta in grid space
allocate(gp2(1:ngx,1:ngy)) ; gp2(:,:) = 0.0 ! v or v*zeta in grid space
allocate(fp1(0:nfx,0:nfy)) ; fp1(:,:) = 0.0 ! u or u*zeta in Fourier space 
allocate(fp2(0:nfx,0:nfy)) ; fp2(:,:) = 0.0 ! v or v*zeta in Fourier space
allocate(kx(0:nkx))        ; kx(:)    = 0.0
allocate(ky(0:nfy))        ; ky(:)    = 0.0
allocate(kx2(0:nkx))       ; kx2(:)   = 0.0
allocate(ky2(0:nfy))       ; ky2(:)   = 0.0
allocate(pkx2(0:nkx))      ; pkx2(:)  = 0.0
allocate(pky2(0:nfy))      ; pky2(:)  = 0.0

allocate(cikx(0:nkx))      ; cikx(:)   = (0.0,0.0)
allocate(ciky(0:nfy))      ; ciky(:)   = (0.0,0.0)
allocate(cli(0:nkx,0:nfy)) ; cli(:,:)  = (0.0,0.0)
allocate(cvo(0:nkx,0:nfy)) ; cvo(:,:)  = (0.0,0.0) ! vorticity in fourier domain
allocate(cnl0(0:nkx,0:nfy)); cnl0(:,:) = (0.0,0.0) ! nonlinear time derivative
allocate(cnl1(0:nkx,0:nfy)); cnl1(:,:) = (0.0,0.0) ! nonlinear previous time deriv.
allocate(cnl2(0:nkx,0:nfy)); cnl2(:,:) = (0.0,0.0) ! nonlinear preprevious time deriv.

return
end subroutine allocate_arrays


! ********
! * QUAD *
! ********

program quad
use quadmod
implicit none

real(8) :: time
real(8) :: dt1
real    :: cpu_start
real    :: cpu_stop
real    :: runtime
real    :: rundnq
real    :: rundnt
integer :: ispm


call cpu_time(cpu_start)

write (*,*) '********'
write (*,*) '* Quad *'
write (*,*) '********'

call iniquad  ! get the initial vorticity field

! open output file in service format

open(nunserv,file='quad.srv',form='unformatted')

! initialize the time variables.

time = 0.d0
nstep = 0

if (nforc .ge. 1) call iniforcing

! >>>>>>>>>>  first steps  <<<<<<<<<<<<<<<<<<<<<<<<<<<

! spectral: z(t)

dt1=dt/2.

! 1: leap-frog step
! do half a step by forward euler

call jacob(cnl1)
cnl2(:,:) = cnl1(:,:)
call teuler(cnl1,dt1)

! do one step by leap-frog without forcing

call jacob(cnl0)    
call teuler(cnl0,dt)
nstep=nstep+1
time=time+dt

! spectral: z(2dt),cnl0(dt),cnl1(dt),cnl2(0)
! 2: 2nd. order adam-bashforth

call jacob(cnl0)
call tadba2(cnl0,cnl1)
cnl1(:,:) = cnl0(:,:)
nstep=nstep+1
time=time+dt

! 3: 3rd. order adam-bashforth
! >>>>>>>>>>  time evolution <<<<<<<<<<<<<<<<<<<<<<<<

do while (nstep <= nstop)
   call jacob(cnl0)

   ! time integration part

   call tadba3 (cnl0,cnl1,cnl2)
   cnl2(:,:) = cnl1(:,:)
   cnl1(:,:) = cnl0(:,:)
  
   if (mod(nstep,1000) == 0) write(*,*)' time step ',nstep

   nstep = nstep + 1
   time = time + dt
enddo
close (nunserv)

call cpu_time(cpu_stop)
runtime = cpu_stop - cpu_start
rundnq  = runtime / nstop * rnxy *1.0e9
rundnt  = rundnq  * rnx
ispm    = nint(60.0 * nstop / runtime)
open(nunstat,file='quad.rt',position='append')
write(nunstat,'(i4," x ",i4,i10," [steps / minute]",2f10.4)') &
      ngx,ngx,ispm,rundnq,rundnt
close(nunstat)

stop
end program quad


! ========
! DEFKBETA
! ========

subroutine defkbeta
use quadmod
implicit none

real(8) :: k2
complex(8) :: aniso,scra
integer :: i,j

! termine per l'integrazione
!
! sig: viscosita' ultravioletta
! lam: viscosita infrarossa
! beta: strength of coriolis force

! {linear operator} =
! exp[-dt*(sig*(k**2)**psig + lam*(k**2)**plam) * (k**2)/(k**2+alpha)] *
! (dt*beta*kx/(k**2+alpha))

do j = 0, nfy   
   do i = 0, nkx
      k2 = (pkx2(i)+pky2(j))
      scra = (sig*k2**psig + lam*k2**plam)*k2/(k2+alpha)
      aniso = cmplx(0.0,beta * ki(i) / (k2+alpha))
      cli(i,j) = exp(-dt*scra)*exp(dt*aniso)
   enddo
enddo
!write (88,*) "L"
!write (88,'(4e12.4)') cli
cli(0,0) = (0.0,0.0)

return
end


!*********! 
! iniquad !
!*********! 
subroutine iniquad
use quadmod
implicit none

integer :: i,j,ios
real(8) :: z_r,z_i

namelist /quad_nl/ nstop,nout,ndiag,ncfl,dt,alpha,beta,lx,ly,sig,psig,lam, &
                   plam,nforc,kfmin,kfmax,ampforc,tforc

! read relative vorticity gridpoint array

open (nuini,file='init.srv',form='unformatted')
read (nuini) ihead
ngx = ihead(5)
ngy = ihead(6)
allocate(gvo(ngx,ngy))
read (nuini) gvo(:,:)
close(nuini)

call allocate_arrays


open (nunl,file="quad_namelist",status='old',iostat=ios)
if (ios == 0) then
   read (nunl,quad_nl)
   close(nunl)
endif

write (*,quad_nl)

! inizializza i numeri onda e quello che serve per le fft

call init_fft

call grid_to_fourier(gvo,cvo)
cvo(0,0) = (0.0,0.0) ! force mean of vorticity to zero

! inizializza quanto serve per le viscosita uv, ir e beta

call defkbeta

return
end


! integrazione secondo euler sulla parte non lineare ed
! esatta su quella lineare.

! z=cvo(t)
! cnl=nlt(t)
! znew=cvo(t+1)
! dimensioni: (0:ngx+1,0:ngy+1)   versione ngx, ngy

! ======
! TEULER
! ======

subroutine teuler(cnl,pdt)
use quadmod
implicit none

real(8), intent(in) :: pdt
complex(8), intent(in) :: cnl(0:nkx,0:nfy)

integer :: i
integer :: j

do j = 0, nfy
   do i = 0, nkx
      cvo(i,j)   = (cvo(i,j)+pdt*cnl(i,j))*cli(i,j) 
   enddo
enddo

cvo(0,0)= (0.0,0.0)

return
end


! calcolo del termine non lineare j=d(psi,z)/d(x,y)
! in input: z
! in output: cnl = - d(psi,z)/d(x,y) 
!            v (=v) e um(=-u)
! dimensioni: (n,n) - spazio di fourier (codificato)
! versione ngx, ngy

! =====
! JACOB
! =====

subroutine jacob(cnl)
use quadmod
implicit none

complex(8) :: cnl(0:nkx,0:nfy)

call deriv

! trasforma in reale

call fourier_to_grid(cvo,gvo)
call fourier_to_grid(fp1,gp1)
call fourier_to_grid(fp2,gp2)

if(nstep.gt.0.and.mod(nstep,nout).eq.0) call ecrire(nunserv,gvo)

! check whether stability is obtained.

if(mod(nstep,ncfl).eq.0) call cfl_test

! calcola i prodotti psi1 = v*z e psi2 = um*z

gp1(:,:) = gp1(:,:) * gvo(:,:)
gp2(:,:) = gp2(:,:) * gvo(:,:)

! trasforma in fourier

call grid_to_fourier(gp1,fp1)
call grid_to_fourier(gp2,fp2)

! calcola cnl= - jacobiano

call nonlinear_terms(cnl)
return
end


! calcola (u,-v) a partire da z nello spazio di fourier 
! fp1 = ki * z / (-k**2) =  v
! fp2 = kj * z / (-k**2) = -u

! =====
! DERIV
! =====

subroutine deriv
use quadmod
implicit none

integer :: i,j
real(8) :: k2

do j = 0, nfy
   do i = 0, nkx
      k2 = 1.0d0 / (pkx2(i) + pky2(j) + alpha)
      fp1(i+i  ,j) =  ki(i) * aimag(cvo(i,j)) * k2
      fp1(i+i+1,j) = -ki(i) * real (cvo(i,j)) * k2
      fp2(i+i  ,j) =  kj(j) * aimag(cvo(i,j)) * k2
      fp2(i+i+1,j) = -kj(j) * real (cvo(i,j)) * k2
   enddo
enddo

return
end


! =======
! FORCING
! =======

subroutine forcing
use quadmod
implicit none
integer :: i,j,ic,ifk
real(8) :: k2
complex(8) :: psif
real(8) :: eni,enf,ran3,ran4

! creates a forcing streamfunction (psifr,psifi)
! the phase phi is uniformly distributed
! between 0 and 2 pi, kept constant for the time tforc
! and then changed.

eni=0.d0
enf=0.d0

select case (nforc)
  case (1)
    ic=nstep+1
    if (mod(ic,itau).eq.0.) phi=ran3(idum)*twopi*ci
    ! print*, phi,itau,ic
    psif=ampcoeff*exp(phi)
    ! add the forcing to the spectral vorticity z -> z+ k2 *psif
    do ifk=1,nk
      i = in(ifk,1)
      j = in(ifk,2)
      k2 = pkx2(i)+pky2(j)+alpha
      eni=eni+(cvo(i,j)*cvo(i,j))/k2
      cvo(i,j) = cvo(i,j)- k2*psif*rnxy
      enf=enf+(cvo(i,j)*cvo(i,j))/k2
    enddo
  case (2)
    do ifk=1,nk
      i = in(ifk,1)
      j = in(ifk,2)
      call random_number(ran4)
      psif = ampcoeff*exp(ci*twopi*ran4)
      k2 = pkx2(i)+pky2(j)+alpha
      eni=eni+(cvo(i,j)*cvo(i,j))/k2
      cvo(i,j) = cvo(i,j)-k2*psif*rnxy
      enf=enf+(cvo(i,j)*cvo(i,j))/k2
    enddo
end select

write (32,*) (enf-eni)/dt

return
end



subroutine iniforcing
use quadmod
implicit none

real(8) :: k2,ent,ran3
integer :: i,j

! initialize fortran random number generator
call random_seed()

! initialize internal random number generator
ent=0.d0
do j = 0, nfy
  do i = 0, nkx
    k2 = pkx2(i)+pky2(j)
    if(j.ne.0.or.i.ne.0)then
      if (sqrt(k2).gt.kfmin.and.sqrt(k2).lt.kfmax) then
        nk = nk + 1
        in(nk,1) = i
        in(nk,2) = j
        ent=ent+k2+alpha
      endif
    endif
  enddo
enddo


ampcoeff=sqrt(ampforc*dt/ent)*nxy
if (tforc.le.dt) then 
  itau=1
  else
    itau=tforc/dt
endif
idum=-12895673
phi=ran3(idum)*twopi
print*, phi 
print *, "the forcing ring contains ",nk," wavevectors."
return
end

real(8) function ran3(idum)
save
parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
dimension ma(55)
data iff /0/
if (idum.lt.0.or.iff.eq.0)then
   iff=1
   mj=mseed-iabs(idum)
   mj=mod(mj,mbig)
   ma(55)=mj
   mk=1
   do i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if (mk.lt.mz)mk=mk+mbig
      mj=ma(ii)
   enddo
   do k=1,4
      do i=1,55
         ma(i)=ma(i)-ma(1+mod(i+30,55))
         if(ma(i).lt.mz)ma(i)=ma(i)+mbig
      enddo
   enddo
   inext=0
   inextp=31
   idum=1
endif
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.mz)mj=mj+mbig
ma(inext)=mj
ran3=mj*fac
return
end


! integrazione adams-bashforth sulla parte non lineare ed
! esatta su quella lineare.
! schema di legras
! z=z(t)
! cnl=nlt(t)
! cnlold=nlt(t-1)
!       cnlold2=nlt(t-2)
! znew=z(t+1)
! dimensioni: (0:ngx/2,0:ngy+1) 
! versione ngx, ngy

! ======
! TADBA2
! ======

subroutine tadba2 (cnl,cnlold)
use quadmod
implicit none

complex(8) :: cnl(0:nkx,0:nfy),cnlold(0:nkx,0:nfy)
integer :: i,j
real(8) :: dt2

dt2=dt*0.5

do j = 0, nfy
   do i = 0, nkx
      cvo(i,j)   = cvo(i,j)*cli(i,j)                          &
    &                    + dt2*(3.0*cnl(i,j)*cli(i,j) -       &
    &                      cnlold(i,j)*(cli(i,j)*cli(i,j)))
   enddo
enddo

cvo(0,0) = (0.0,0.0)

if (nforc .ge. 1) call forcing

return
end
! ------------------------------------------------------------
! ------------------------------------------------------------

subroutine tadba3 (cnl,cnlold,cnlold2)
use quadmod
implicit none

complex(8) :: cnl(0:nkx,0:nfy),cnlold(0:nkx,0:nfy)
complex(8) :: cnlold2(0:nkx,0:nfy)
integer :: i,j
real(8) :: c

c = dt/12.0

do j = 0, nfy
   do i = 0, nkx
      cvo(i,j) = cvo(i,j)*cli(i,j)                         + &
    &       c*(23.0*cnl(i,j)*cli(i,j)                      - &
    &       16.0*cnlold(i,j)*(cli(i,j)*cli(i,j))           + &
    &       5.0*cnlold2(i,j)*(cli(i,j)*cli(i,j)*cli(i,j)))
   enddo
enddo

cvo(0,0) = (0.0,0.0)

if (nforc .ge. 1) call forcing
       
return
end


! scrive n real*4 nel file binario ku dal vettore u. complementare ad lire

! ======
! ECRIRE
! ======

subroutine ecrire(ku,pz)
use quadmod
implicit none

integer :: ku
real(8) :: pz(ngx,ngy)
integer :: yy,mo,dd,hh,mm,ii

! Build a header for service format

mm = mod(nstep,60)
ii = nstep / 60
hh = mod(ii,24)
ii = ii / 24
dd = mod(ii,30) + 1
ii = ii / 30
mo = mod(ii,12) + 1
yy = ii / 12

ihead(1) =    138 ! relative vorticity
ihead(2) =      0 ! no level
ihead(3) = dd + 100 * mo + 10000 * yy
ihead(4) = mm + 100 * hh
ihead(5) =    ngx ! 1st. dim
ihead(6) =    ngy ! 2nd. dim
ihead(7) =      0
ihead(8) =      0

write (ku) ihead
write (ku) pz(:,:)

return
end


! calcolo del termine non lineare.
! in input: zp1 = z*d1psi = z*v
!    zp2 = z*d2psi =  z*um (=-z*u)
! in output: cnl = - d(psi,z)/d(x,y) 
! dimensioni: (n,n) - spazio di fourier (codificato)
! versione ngx, ngy

! ===============
! NONLINEAR_TERMS
! ===============

subroutine nonlinear_terms(pnl)
use quadmod
implicit none

complex(8) :: pnl(0:nkx,0:nfy) 
integer i,j

do j = 0, nfy
   do i = 0, nkx
      pnl(i,j) = cmplx(kj(j) * fp1(i+i+1,j) - ki(i) * fp2(i+i+1,j) &
                      ,ki(i) * fp2(i+i  ,j) - kj(j) * fp1(i+i  ,j))
   enddo
enddo

return
end


! write energy and enstrophy
! z in spettrale codificato

! ========
! CFL_TEST
! ========

subroutine cfl_test
use quadmod
implicit none

real(8) :: cfl
real(8) :: maxu,maxv
real(8) :: ee,zz

zz = 0.5 * rnxy *  sum(gvo(:,:) * gvo(:,:))
ee = 0.5 * rnxy * (sum(gp1(:,:) * gp1(:,:)) + sum(gp2(:,:) * gp2(:,:)))

! checking whether the numerical stability is obtained, using the
! courant number (cfl criteria). this number should be less than 0.5.

maxu = maxval(abs(gp2(:,:))) * dt * ngx / lx
maxv = maxval(abs(gp1(:,:))) * dt * ngy / ly

if (maxu > maxv) then
   cfl = maxu
else
   cfl = maxv
endif

write(nuphys,*) real(nstep*dt),ee,zz,real(cfl)

return
end


!******************************************************************************|
! fft.f90, the fft package for diablo.                             version 0.3e
!
! this file isolates all calls to the fftw package (available at: www.fftw.org)
! these wrapper routines were written by t. bewley (spring 2001).
!******************************************************************************|


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! the arrangement of the significant real numbers in the arrays (denoted by +)
! in physical space, in fourier space, and in fourier space after packing are
! shown below for the 2d (x-z) plane.  the third direction (y) is handled in
! an identical matter as the z direction shown here.
!
!       oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
!       oooooooooooooooooo         oooooooooooooooooo         oooooooooooooooooo
! ngy-1 ++++++++++++++++oo     -1  ++++++++++++oooooo         oooooooooooooooooo
!       ++++++++++++++++oo     -2  ++++++++++++oooooo         oooooooooooooooooo
!       ++++++++++++++++oo     -3  ++++++++++++oooooo         oooooooooooooooooo
!       ++++++++++++++++oo         ++++++++++++oooooo         oooooooooooooooooo
!       ++++++++++++++++oo    -nky ++++++++++++oooooo         oooooooooooooooooo
!       ++++++++++++++++oo         oooooooooooooooooo     -1  ++++++++++++oooooo
!       ++++++++++++++++oo         oooooooooooooooooo     -2  ++++++++++++oooooo
!       ++++++++++++++++oo         oooooooooooooooooo     -3  ++++++++++++oooooo
!       ++++++++++++++++oo         oooooooooooooooooo         ++++++++++++oooooo
!       ++++++++++++++++oo         oooooooooooooooooo    -nky ++++++++++++oooooo
!       ++++++++++++++++oo     nky ++++++++++++oooooo     nky ++++++++++++oooooo
!       ++++++++++++++++oo         ++++++++++++oooooo         ++++++++++++oooooo
!    3  ++++++++++++++++oo      3  ++++++++++++oooooo      3  ++++++++++++oooooo
!    2  ++++++++++++++++oo      2  ++++++++++++oooooo      2  ++++++++++++oooooo
!    1  ++++++++++++++++oo      1  ++++++++++++oooooo      1  ++++++++++++oooooo
!    0  ++++++++++++++++oo      0  +o++++++++++oooooo      0  +o++++++++++oooooo
!       ^^^^           ^           ^ ^ ^     ^                ^ ^ ^     ^
!       0123           ngx-1        0 1 2     nkx              0 1 2     nkx
!
!       physical space              fourier space         fourier space (packed)
!
! after the real->fourier transform, the significant coefficients are put next
! to each other in the array, so a loop such as
!
!        do j=0,nfy           [where nfy = 2*nky = 2*(ngy/3) ]
!          do i=0,nkx          [where  nkx = ngx/3             ]
!            cp(i,j)= ...
!          end do
!        end do
!
! includes all the fourier coefficients of interest.  the subsequent loops in
! fourier space just work on these coefficients in the matrix.
!  
! before a fourier->real transform, the significant coefficients are unpacked
! and the higher wavenumbers are set to zero before the inverse transform.
! this has the effect of doing the required dealiasing.
!
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine init_fft
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      use quadmod
      implicit none
      integer :: i,j

        do i=0,nkx
          kx(i)=i*pi
          kx2(i)=kx(i)*kx(i)
          pkx2(i)=kx2(i)/(pi*pi)
          cikx(i)=ci*kx(i)
        end do

        pkx2(0) = 1.0 ! avoid divison by zero

        do j=0,nky
          ky(j)=j*pi
        end do
        do j=1,nky
          ky(nfy+1-j)=-j*pi
        end do
        do j=0,nfy
          ky2(j)=ky(j)*ky(j)
          ciky(j)=ci*ky(j)
          pky2(j)=ky2(j)/(pi*pi)
        end do

!write (88,*) "pkx2"
!write (88,'(4e16.8)') pkx2
!write (88,*) "pky2"
!write (88,'(4e16.8)') pky2

      return
      end

! **************
! * check_diff *
! **************

subroutine check_diff(pold,pnew,ytext)
use quadmod
real(8) :: pold(0:ngx+1,0:ngy+1)
real(8) :: pnew(0:nfx,0:nfy)
character(*) :: ytext

write (88,*) ytext
write (88,*) "max old  = ",maxval(pold)
write (88,*) "max new  = ",maxval(pnew)
write (88,*) "max diff = ",maxval(pnew-pold(0:nfx,0:nfy))
write (88,*) "min old  = ",minval(pold)
write (88,*) "min new  = ",minval(pnew)
write (88,*) "min diff = ",minval(pnew-pold(0:nfx,0:nfy))
return
end
