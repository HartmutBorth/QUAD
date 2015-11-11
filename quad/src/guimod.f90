! ========
! GUISTART
! ========

subroutine guistart
implicit none

integer :: ngui    =  1 ! set by caller
integer :: model   =  4 ! QUAD
integer :: nguidbg =  0 ! Debug mode
integer :: nlat    = 32 ! dummy value
integer :: mrpid   =  0 ! # process ID
integer :: mrnum   =  1 ! # of instances

character (80) :: yplanet = "QUAD" ! dummy

if (ngui == 0) return

call initgui(model,nguidbg,nlat,mrpid,mrnum,trim(yplanet)//char(0))

return
end subroutine guistart

! =======
! GUISTOP
! =======

subroutine guistop
implicit none

integer :: ngui = 1 ! set by caller

if (ngui > 0) call guiclose

return
end subroutine guistop

! ============
! GUISTEP_QUAD
! ============

subroutine guistep_quad
implicit none

integer  :: nshutdown     ! user action
integer  :: ndatim(6) = 0 ! date & time info
real (4) :: parc(5)       ! for timeseries display


interface 
   integer(kind=4) function iguistep(parc,idatim)
      real   (kind=4), intent(inout) :: parc(*)
      integer(kind=4), intent(in)    :: idatim(6)
   end function iguistep
end interface

integer (kind=4) idatim(6)


parc(1)   = 1.0
parc(2)   = 2.0
parc(3)   = 3.0
parc(4)   = 4.0
parc(5)   = 5.0
idatim(:) = ndatim(:)
nshutdown = iguistep(parc,idatim)    ! GUI event handler
   
return
end subroutine guistep_quad
