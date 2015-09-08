! **********
! * RESMOD *
! **********
module resmod
use quadmod, only: nuresini,quad_resini, nuresfin, nudiag

integer, parameter :: nresdim  = 200   ! Max number of records
integer            :: nexcheck =   1   ! Extended checks
integer            :: nresnum  =   0   ! Actual number of records
integer            :: nlastrec =   0   ! Last read record
                                          
character (len=16) :: yresnam(nresdim) ! Array of record names

end module resmod


! ************************
! * SUBROUTINE CHECK_RES *
! ************************
subroutine check_res
use resmod

character (len=16) :: yn ! variable name

do
   read (nuresini,IOSTAT=iostat) yn
   if (iostat /= 0) exit
   nresnum = nresnum + 1
   yresnam(nresnum) = yn
   read (nuresini,IOSTAT=iostat)
   if (iostat /= 0) exit
   if (nresnum >= nresdim) then
      write(nudiag,*) 'Too many variables in restart file'
      write(nudiag,*) 'Increase NRESDIM in module resmod'
      write(nudiag,*) '*** Error Stop ***'
      stop
   endif
enddo

write(nudiag,'(a,i4,3a/)') 'Found ',nresnum, &
      ' variables in file <',trim(quad_resini),'>'
do j = 1 , nresnum
   write(nudiag,'(i4," : ",8x,1x,a)') j,yresnam(j)
enddo
nlastrec = nresnum

return
end subroutine check_res


! **********************************
! * SUBROUTINE GET_RESTART_INTEGER *
! **********************************
subroutine get_restart_integer(yn,kv)
use resmod

character (len=*) :: yn
integer :: kv

do j = 1 , nresnum
   if (trim(yn) == trim(yresnam(j))) then
      call fileseek(yn,j)
      read (nuresini) kv
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
use resmod

character (len=*) :: yn
real :: pa(k2,k3)

do j = 1 , nresnum
   if (trim(yn) == trim(yresnam(j))) then
      call fileseek(yn,j)
      read (nuresini) pa(1:k1,:)
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
use resmod

      character (len=*)  :: yn
      character (len=16) :: yy
      integer :: kv

      yy = yn
      write(nuresfin) yy
      write(nuresfin) kv
      return
      end subroutine put_restart_integer


! ********************************
! * SUBROUTINE PUT_RESTART_ARRAY *
! ********************************
subroutine put_restart_array(yn,pa,k1,k2,k3)
use resmod

character (len=*)  :: yn
character (len=16) :: yy
integer :: k1,k2,k3
real :: pa(k2,k3)

yy = yn
write(nuresfin) yy
write(nuresfin) pa(1:k1,1:k3)

return
end subroutine put_restart_array


! ***********************
! * SUBROUTINE FILESEEK *
! ***********************

subroutine fileseek(yn,k)
use resmod

character (len=*)  :: yn
character (len=16) :: yy

if (k <= nlastrec) then
   rewind nuresini
   nlastrec = 0
endif

do
   read (nuresini,iostat=iostat) yy
   if (iostat /= 0) exit
   if (trim(yn) == trim(yy)) return ! success
   read (nuresini,iostat=iostat)    ! skip data
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
use resmod

character (len=*)  :: yn
character (len=16) :: ytmp
integer :: k, knum

knum = 0
do k = 1,nresdim
   ytmp = yresnam(k)
   if (trim(yn) == trim(ytmp)) then 
      knum = k
   endif
enddo

return
end
