! **********
! * RSTMOD *
! **********
module rstmod
use quadmod, only: nurstini,quad_rstini, nurstfin, nudiag

integer, parameter :: nrstdim  = 200   ! Max number of records
integer            :: nexcheck =   1   ! Extended checks
integer            :: nrstnum  =   0   ! Actual number of records
integer            :: nlastrec =   0   ! Last read record
                                          
character (len=16) :: yrstnam(nrstdim) ! Array of record names

end module rstmod


! ************************
! * SUBROUTINE CHECK_RES *
! ************************
subroutine check_rst
use rstmod

character (len=16) :: yn ! variable name


write(nudiag, &
 '(/," *************************************************")')
write(nudiag,'("  Checking for restart file <quad_rstini>")')
write(nudiag, &
 '(" *************************************************",/)')



do
   read (nurstini,IOSTAT=iostat) yn
   if (iostat /= 0) exit
   nrstnum = nrstnum + 1
   yrstnam(nrstnum) = yn
   read (nurstini,IOSTAT=iostat)
   if (iostat /= 0) exit
   if (nrstnum >= nrstdim) then
      write(nudiag,*) 'Too many variables in restart file'
      write(nudiag,*) 'Increase NRESDIM in module rstmod'
      write(nudiag,*) '*** Error Stop ***'
      stop
   endif
enddo

write(nudiag,'(a,i4,3a/)') 'Found ',nrstnum, &
      ' variables in file <',trim(quad_rstini),'>'
do j = 1 , nrstnum
   write(nudiag,'(i4," : ",8x,1x,a)') j,yrstnam(j)
enddo
nlastrec = nrstnum

return
end subroutine check_rst


! **********************************
! * SUBROUTINE GET_RESTART_INTEGER *
! **********************************
subroutine get_restart_integer(yn,kv)
use rstmod

character (len=*) :: yn
integer :: kv

do j = 1 , nrstnum
   if (trim(yn) == trim(yrstnam(j))) then
      call fileseek(yn,j)
      read (nurstini) kv
      nlastrec = nlastrec + 1
      return
   endif
enddo
if (nexcheck == 1) then
   write(nudiag,*) '*** Error in get_restart_integer ***'
   write(nudiag,*) 'Requested integer {',yn,'} was not found'
   stop
endif
return
end subroutine get_restart_integer


! ********************************
! * SUBROUTINE GET_RESTART_ARRAY *
! ********************************
subroutine get_restart_array(yn,pa,k1,k2,k3)
use rstmod

character (len=*) :: yn
real :: pa(k2,k3)

do j = 1 , nrstnum
   if (trim(yn) == trim(yrstnam(j))) then
      call fileseek(yn,j)
      read (nurstini) pa(1:k1,:)
      nlastrec = nlastrec + 1
      return
   endif
enddo
if (nexcheck == 1) then
   write(nudiag,*) '*** Error in get_restart_array ***'
   write(nudiag,*) 'Requested array {',yn,'} was not found'
   stop
endif

return
end subroutine get_restart_array


! ********************************** 
! * SUBROUTINE PUT_RESTART_INTEGER *
! ********************************** 
subroutine put_restart_integer(yn,kv)
use rstmod

      character (len=*)  :: yn
      character (len=16) :: yy
      integer :: kv

      yy = yn
      write(nurstfin) yy
      write(nurstfin) kv
      return
      end subroutine put_restart_integer


! ********************************
! * SUBROUTINE PUT_RESTART_ARRAY *
! ********************************
subroutine put_restart_array(yn,pa,k1,k2,k3)
use rstmod

character (len=*)  :: yn
character (len=16) :: yy
integer :: k1,k2,k3
real :: pa(k2,k3)

yy = yn
write(nurstfin) yy
write(nurstfin) pa(1:k1,1:k3)

return
end subroutine put_restart_array


! ***********************
! * SUBROUTINE FILESEEK *
! ***********************

subroutine fileseek(yn,k)
use rstmod

character (len=*)  :: yn
character (len=16) :: yy

if (k <= nlastrec) then
   rewind nurstini
   nlastrec = 0
endif

do
   read (nurstini,iostat=iostat) yy
   if (iostat /= 0) exit
   if (trim(yn) == trim(yy)) return ! success
   read (nurstini,iostat=iostat)    ! skip data
   if (iostat /= 0) exit
   nlastrec = nlastrec + 1
enddo
write(nudiag,*) 'Variable <',trim(yn),'> not in restart file'

return
end


! *****************************
! * SUBROUTINE CHECK_EQUALITY *
! *****************************
subroutine check_equality(yn,pa,pb,k1,k2)

character (len=*) :: yn
real :: pa(k1,k2)
real :: pb(k1,k2)

do j2 = 1 , k2
   do j1 = 1 , k1
      if (pa(j1,j2) /= pb(j1,j2)) then
         write(nudiag,*) 'No Equality on ',yn,'(',j1,',',j2,')',pa(j1,j2),pb(j1,j2)
         return
      endif
   enddo
enddo
 write(nudiag,*) 'Array {',yn,'} is OK'

return
end


! **********************
! * SUBROUTINE VARSEEK *
! **********************
subroutine varseek(yn,knum)
use rstmod

character (len=*)  :: yn
character (len=16) :: ytmp
integer :: k, knum

knum = 0
do k = 1,nrstdim
   ytmp = yrstnam(k)
   if (trim(yn) == trim(ytmp)) then 
      knum = k
   endif
enddo

return
end
