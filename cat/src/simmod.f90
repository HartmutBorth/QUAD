module simmod
use catmod

!--- flags and switches
logical :: lsimnl = .false.    ! true if <sim_namelist> exists
                               ! and cat is set-up for a predefined
                               ! simulation

integer :: ios_simnl = 1       ! 0 if <sim_namelist> is readable

!--- i/o units
integer, parameter :: nusimnl    = 110  ! sim namelist

!--- i/o file names
character (256) :: sim_namelist = "sim_namelist"

end module simmod


! ***********************
! * SUBROUTINE SIMSTART *
! ***********************
subroutine simstart
use simmod
implicit none


!--- check if sim_namelist is present
inquire(file=sim_namelist,exist=lsimnl)
if (lsimnl) then
  open(nusimnl,file=sim_namelist,iostat=ios_simnl)
else
  return
endif


return
end subroutine simstart


! **********************
! * SUBROUTINE SIMSTEP *
! **********************
subroutine simstep
use simmod
implicit none


return
end subroutine simstep


! **********************
! * SUBROUTINE SIMSTOP *
! **********************
subroutine simstop
use simmod
implicit none


return
end subroutine simstop
