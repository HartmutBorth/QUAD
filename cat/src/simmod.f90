module simmod
use catmod

!--- flags and switches
logical :: lsimnl = .false.    ! .true. if sim_namelist is present in run 
                               ! directory

integer :: ios_simnl = 1       ! 0 if <sim_namelist> is readable

!--- i/o units
integer, parameter :: nusimnl    = 110  ! sim namelist

!--- i/o file names
character (256) :: sim_namelist = "sim_namelist"

!--- predefined experimets
character (256) :: sim = "dec01" ! type of predefined simulation

!--- parameters of dec01 (decaying turbulence)

end module simmod


! ***********************
! * SUBROUTINE SIMSTART *
! ***********************
subroutine simstart
use simmod
implicit none

!--- define sim_namelist
namelist /sim_nl/ sim

!--- check if sim_namelist is present
inquire(file=sim_namelist,exist=lsimnl)


!--- read sim_nl
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
