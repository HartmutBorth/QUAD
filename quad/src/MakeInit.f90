program makeinit
implicit none

integer :: jx
integer :: jy
integer :: nx
integer :: ny

real(8) :: pi
real(8) :: dx
real(8) :: dy 
real(8), allocatable :: z(:,:)

character(80) :: arg

call get_command_argument(1,arg)
read(arg,*) nx
ny = nx

write (*,'("Allocate z(",i4,",",i4,")")') nx,ny

allocate(z(nx,ny))

pi = 4.0 * atan(1.0d0)
dx = 2.0d0 * pi / nx
dy = dx * 2.0d0

do jy = 1 , ny
   do jx = 1 , nx
!      z(jx,jy) = &
!               + 1.1 * cos((jx-1) * dx) &
!               + 0.1 * sin((jy-1) * dx)

!     z(jx,jy) = 0.0d0 

     z(jx,jy) =  cos((jx-1) * dx +(jy-1) * dx)
     z(jx,jy) =  1.0 * cos((jx-1)*dx)   +  3.0 * cos((jx-1)*dx*2.0) &
               + 10.0 * sin((jy-1) * dy) + 20.0 * sin((jy-1) * dy * 2.0)
   enddo
enddo
!    z(nx/2,ny/2) = 2.0

write(50) 138,0,10101,0,nx,ny,0,0
write(50) z
stop
end
      


