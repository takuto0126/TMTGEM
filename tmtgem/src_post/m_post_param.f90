!# Modified on 2018.01.19
!# bl,br,bt,bb became array with size of nslice
module post_param
use param ! 2020.12.10

type param_post
 integer(4)    :: itype_process ! 1:m0(avrage) + dm, 2: m1 & m2 2018.10.24
 character(70) :: mshfile
 !#
 integer(4)    :: ifileflag   !1:condfile, 2:modelfile 2018.10.24
 character(70) :: m0condfile
 character(70) :: dmcondfile
 !#
 character(70) :: connectfile ! 2018.10.24
 character(70) :: m0modelfile ! 2018.10.24
 character(70) :: dmmodelfile ! 2018.10.24
 !#
 character(70) :: outputfile(5) ! 2018.10.24
 real(8)       :: xyprof(2)
 integer(4)    :: nlayer
 real(8),allocatable,dimension(:) :: zlayer ! [nlayer + 1]
end type

type profile
 integer(4)   :: nlayer
 real(8),allocatable,dimension(:,:) :: bottoplayer ![2,nlayer] 1:botz,2topz
 real(8),allocatable,dimension(:)   :: value       ![nlayer]
end type

type param_slice ! 2017.12.20
 character(70) :: mshfile
 ! 0:condfile, 1:modelfile 2018.01.05
 ! 2: dcond,   3:dmodel    2018.03.19
 integer(4)    :: ifileflag
  !# model file version ( ifileflag = 0 or 2 ) ===== ! 2018.03.12
 integer(4)  :: ncondfile                           ! 2018.01.05
 character(70),allocatable,dimension(:) :: condfile ! 2018.01.05
 !# model file version ( ifileflag = 1 ) ========== ! 2018.01.05
 character(70) :: connectfile
 integer(4)    :: nmodelfile=0 ! 2017.12.21
 character(70),allocatable,dimension(:) :: modelfile ! 2017.12.21
 !#
 character(70),allocatable,dimension(:) :: outpolygonhead
 !# upper and lower bound 2018.04.17
 integer(4)  :: ibound ! 0:nothing, 1:lower bound, 2: upper bound
 real(8)     :: lbound,ubound  ! 2018.04.17
 !# Ax +By +Cz +D = 0
 !# C=0 is assumed on 2017.12.20
 integer(4)  :: nslice ! 2017.12.21
 real(8),allocatable,dimension(:)  :: A  ! A should be equally larger than 0
 real(8),allocatable,dimension(:)  :: B
 real(8),allocatable,dimension(:)  :: D
 real(8),allocatable,dimension(:)  :: bl,br,bb,bt ! 2018.01.19
end type


contains

!========================================================
!# 2017.12.20
subroutine sliceparamread(g_paramslice)
implicit none
type(param_slice),intent(out) :: g_paramslice
integer(4) :: input=1,i
integer(4) :: ifileflag ! 2018.01.05
integer(4),       parameter :: n = 1000 ! 2020.12.10
character(200),dimension(n) :: lines    ! 2020.12.10
character(100)              :: paramfile ! 2020.12.10

write(*,*) ""
write(*,*) "<Please input the forward parameter file>" ! 2020.12.10
read(*,'(a)') paramfile                 ! 2020.12.10
call readcontrolfile(paramfile,n,lines,ikeep=1) ! 2020.12.10 see m_param.ctl

open(input,file="tmp.ctl") ! 2020.12.10
!# read parameter
read(input,'(20x,a)') g_paramslice%mshfile
write(*,*) "Please input 0: condfile, 1: modelfile, 2:dcond,3:dcond" ! 2018.01.05
read(input,'(20x,i10)') g_paramslice%ifileflag        ! 2018.01.05
ifileflag = g_paramslice%ifileflag                    ! 2018.01.05
   if ( ifileflag .eq. 0 .or. ifileflag .eq. 2 ) then !==  read  cond file
      read(input,'(20x,i10)') g_paramslice%ncondfile      ! 2018.01.06
      write(*,'(i10)')        g_paramslice%ncondfile      ! 2020.12.10
      allocate(g_paramslice%condfile(g_paramslice%ncondfile)) ! 2018.01.06
      allocate(g_paramslice%outpolygonhead(g_paramslice%ncondfile)) ! 2018.01.06
      do i=1,g_paramslice%ncondfile                        ! 2018.01.06
         read(input,'(20x,a)') g_paramslice%outpolygonhead(i)! 2017.01.06
         read(input,'(20x,a)') g_paramslice%condfile(i)      ! 2018.01.06
         write(*,'(a)') g_paramslice%outpolygonhead(i)       ! 2020.12.10
         write(*,'(a)') g_paramslice%condfile(i)             ! 2020.12.10
      end do ! 2017.12.21

   else if ( ifileflag .eq. 1 .or. ifileflag .eq. 3) then ! read model file
      read(input,'(20x,a)') g_paramslice%connectfile
      write(*,'(a)') g_paramslice%connectfile              ! 2020.12.10
      read(input,'(20x,i10)') g_paramslice%nmodelfile      ! 2017.12.21
      write(*,'(i10)') g_paramslice%nmodelfile             ! 2020.12.10
      allocate(g_paramslice%modelfile(g_paramslice%nmodelfile)) ! 2017.12.21
      allocate(g_paramslice%outpolygonhead(g_paramslice%nmodelfile)) ! 2017.12.21
      do i=1,g_paramslice%nmodelfile                       ! 2017.12.21
         read(input,'(20x,a)') g_paramslice%outpolygonhead(i)! 2017.21.21
         write(*,*) g_paramslice%outpolygonhead(i)           ! 2020.12.10
         read(input,'(20x,a)') g_paramslice%modelfile(i)     ! 2017.21.21
         write(*,*) g_paramslice%modelfile(i)                 ! 2020.12.10
      end do ! 2017.12.21

   else    ! ============================================! 2018.01.06
      write(*,*) "GEGEGE ifileflag",ifileflag,"should be 0 or 1"
      stop   ! 2018.01.06
  end if  ! ============================================! 2018.01.06

! upper or lower bound 2018.04.17
write(*,*) "ibound 0: no bound, 1:lower, 2: upper bound" ! 2018.04.17
read(input,'(20x,i10)') g_paramslice%ibound    ! 2018.04.17
write(*,*) "ibound",g_paramslice%ibound        ! 2020.12.10
if      ( g_paramslice%ibound .eq. 1 ) then
 read(input,'(20x,g15.7)') g_paramslice%lbound ! 2018.04.17
else if ( g_paramslice%ibound .eq. 2 ) then
 read(input,'(20x,g15.7)') g_paramslice%ubound ! 2018.04.17
end if
read(input,'(20x,i10)') g_paramslice%nslice ! 2017.12.21
allocate(g_paramslice%A(g_paramslice%nslice))! 2017.12.21
allocate(g_paramslice%B(g_paramslice%nslice))! 2017.12.21
allocate(g_paramslice%D(g_paramslice%nslice))! 2017.12.21
allocate(g_paramslice%bl(g_paramslice%nslice))! 2018.01.19
allocate(g_paramslice%br(g_paramslice%nslice))! 2018.01.19
allocate(g_paramslice%bb(g_paramslice%nslice))! 2018.01.19
allocate(g_paramslice%bt(g_paramslice%nslice))! 2018.01.19
do i=1,g_paramslice%nslice ! 2017.12.21
write(*,*) "Please input A (should be >= 0)"
read(input,'(20x,g15.7)') g_paramslice%A(i) ! 2017.12.21
read(input,'(20x,g15.7)') g_paramslice%B(i) ! 2017.12.21
read(input,'(20x,g15.7)') g_paramslice%D(i) ! 2017.12.21
read(input,'(20x,g15.7)') g_paramslice%bl(i) ! boundary left   2018.01.19
read(input,'(20x,g15.7)') g_paramslice%br(i) ! boundary right  2018.01.19
read(input,'(20x,g15.7)') g_paramslice%bb(i) ! boundary bottom 2018.01.19
read(input,'(20x,g15.7)') g_paramslice%bt(i) ! boundary top    2018.01.19
end do ! 2017.12.21

write(*,'(a)') "    1st slice | 2nd slice |"
write(*,10)    " A ",g_paramslice%A
write(*,10)    " B ",g_paramslice%B
write(*,10)    " D ",g_paramslice%D
write(*,10)    "lb ",g_paramslice%bl
write(*,10)    "rb ",g_paramslice%br
write(*,10)    "bb ",g_paramslice%bb
write(*,10)    "tb ",g_paramslice%bt

write(*,*) "g_paramslice%nmodelfile",g_paramslice%nmodelfile
!do i=1,g_paramslice%nmodelfile
!write(*,*) "polygonhead i",i,g_paramslice%outpolygonhead(i)
!end do

do i=1,g_paramslice%nslice
if (g_paramslice%A(i) .lt. 0.d0) then
 write(*,*) " GEGEGE! g_paramslice%A=",g_paramslice%A(i)
 write(*,*) "should be >= 0.d0 !!"
stop
end if
end do

close(input) ! 2020.12.10

return
10 format(a,3f9.3)

end subroutine

!========================================================
subroutine postparamread(g_parampost)
implicit none
type(param_post),intent(out) :: g_parampost
integer(4)  :: input=5
integer(4)  :: i

!# header
read(input,*)

!# file data
read(input, '(20x,i10)')g_parampost%itype_process !1:m0,dm, 2: m1,m2
read(input, '(20x,a)') g_parampost%mshfile
write(*,*) "g_parampost%mshfile",g_parampost%mshfile
read(input, '(20x,i10)') g_parampost%ifileflag ! 2018.10.24
if (g_parampost%ifileflag .eq. 1 ) then ! condfile 2018.10.24
 read(input, '(20x,a)') g_parampost%m0condfile
 read(input, '(20x,a)') g_parampost%dmcondfile
else if (g_parampost%ifileflag .eq. 2) then ! modelfile 2018.10.24
 read(input, '(20x,a)') g_parampost%connectfile ! 2018.10.24
 read(input, '(20x,a)') g_parampost%m0modelfile  ! 2018.10.24
 read(input, '(20x,a)') g_parampost%dmmodelfile  ! 2018.10.24
else                   ! 2018.10.24
 write(*,*) "GEGEGE"   ! 2018.10.24
 stop                  ! 2018.10.24
end if                 ! 2018.10.24
do i=1,5 ! 2018.10.24
 read(input,'(20x,a)') g_parampost%outputfile(i)
end do

!# layer data
read(input,'(20x,2f12.5)') g_parampost%xyprof
read(input,'(20x,i10)')  g_parampost%nlayer
write(*,*) "nlayer=",g_parampost%nlayer
allocate( g_parampost%zlayer(g_parampost%nlayer + 1))
 do i=1,g_parampost%nlayer + 1
  read(input,'(20x,g15.7)') g_parampost%zlayer(i)
 end do
end subroutine

end module
