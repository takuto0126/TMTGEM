!# coded on 2022.10.24
module cond_model_type
use param
use mesh_type
use cond_type
implicit none

type info_cond_model
 integer(4)     :: iflag_cond ! 0: condfile, 1:model file
 integer(4)     :: iflag_ele
 character(100) :: inputcondfile
 character(100) :: modelfile
 character(100) :: connectfile
 character(70)  :: mshfile
 character(70)  :: topo_2d_mshfile ! 2022.10.24
 character(70)  :: outcondfile
 integer(4)     :: ncuboid
 real(8)        :: rhohomo
 type(cuboid),allocatable,dimension(:) :: g_cuboid
end type

contains
!##################################################### SETCOND
!# coded on 2018.10.04
!# case where rhohomo is given
subroutine setcond(h_cond,g_mesh,rhohomo)
 use mesh_type
 use param     ! 2018.10.05
 implicit none
 type(mesh),       intent(in)    :: g_mesh
 real(8),          intent(in)    :: rhohomo
 type(param_cond), intent(inout) :: h_cond
 integer(4)                      :: nphys2,ntet,i
 integer(4),allocatable,dimension(:,:) :: n4flag
 real(8),   allocatable,dimension(:)   :: rho, sigma

 !# set
 ntet   = g_mesh%ntet
 allocate( n4flag(ntet,2) )
 n4flag = g_mesh%n4flag

 !# cal culate nphys2
 nphys2 = 0
 do i=1,ntet
  if ( n4flag(i,2) .ge. 3 ) nphys2 = nphys2 + 1 ! for TMTGEM 2024.02.01
 end do
 write(*,*) "nphys2=",nphys2

 !# set output
 allocate(h_cond%sigma(nphys2))
 allocate(h_cond%rho(  nphys2))
 allocate(h_cond%index(nphys2))
 h_cond%nphys2   = nphys2
 h_cond%rho(:)   = rhohomo
 h_cond%sigma(:) = 1.d0/rhohomo

 return
 end
!##################################################################################
subroutine READGENCOND(g_info_cond_model,h_mesh,g_cond)
 implicit none
 type(info_cond_model),intent(in)  :: g_info_cond_model
 type(mesh),           intent(in)  :: h_mesh
 type(param_cond),     intent(out) :: g_cond
 integer(4)     :: iflag_cond   ! 0: condfile, 1:model file
 character(100) :: connectfile,modelfile
 real(8)        :: rhohomo

  !#[1]## set
    iflag_cond = g_info_cond_model%iflag_cond ! 0: condfile, 1:model file
    !# flag = 0 requires:
     rhohomo = g_info_cond_model%rhohomo            ! 2018.10.04
    !# flag = 1 requires:
     g_cond%condfile = g_info_cond_model%inputcondfile      ! 2017.07.19
    !# flag = 2 requires:
     connectfile = g_info_cond_model%connectfile
     modelfile   = g_info_cond_model%modelfile

  !#[2]## generate cond from either of cond or model
   if     ( iflag_cond .eq. 0 ) then                ! 0 for homogeneous
     rhohomo = g_info_cond_model%rhohomo   
     CALL SETCOND(g_cond,h_mesh,rhohomo)     

   elseif ( iflag_cond .eq. 1 ) then                ! 1 for read condfile
     CALL READCOND(g_cond)                          ! 2017.07.19

   elseif ( iflag_cond .eq. 2 ) then                ! 2 for read modelfile
     CALL READMODEL2COND(g_cond,connectfile,modelfile)

   else                                                ! 2018.10.04
     write(*,*) "GEGEGE! iflag_cond",iflag_cond   ! 2018.10.04
     stop                                               ! 2018.10.04

    end if

  !#[3]## set nphys1
    CALL SETNPHYS1INDEX2COND(h_mesh,g_cond)            ! set nphys1

return
end

!############################################## subroutine setNPHYS1INDEX2COND
!Coded on 2017.05.12
subroutine SETNPHYS1INDEX2COND(g_mesh,g_cond)
use param
use mesh_type
implicit none
type(mesh),      intent(in)    :: g_mesh
type(param_cond),intent(inout) :: g_cond
integer(4) :: nphys1,nphys2,i,ntet

  nphys2        = g_cond%nphys2 ! # of elements in 2nd physical volume (land)
  nphys1        = g_mesh%ntet - nphys2
  ntet          = g_mesh%ntet
  g_cond%nphys1 = nphys1
  g_cond%ntet   = ntet
  write(*,*) "nphys1=",nphys1,"nphys2=",nphys2,"ntet=",g_mesh%ntet
  do i=1,nphys2
   g_cond%index(i) = nphys1 +i ! element id for whole element space
  end do

return
end

!##################################################### READMODEL2COND
!# coded on 2018.06.21
subroutine readmodel2cond(r_cond,connectfile,modelfile)
! use modelpart
 use param
 implicit none
 type(param_cond), intent(inout)     :: r_cond
 character(70),    intent(in)        :: connectfile
 character(70),    intent(in)        :: modelfile
 integer(4)                          :: nphys2,nmodel,nmodel2
 integer(4)                          :: i
 integer(4),allocatable,dimension(:) :: id,ele2model
 real(8),   allocatable,dimension(:) :: rho

 !#[1]## read connectfile
  open(1,file=connectfile,status='old',err=90)
  read(1,*) nphys2,nmodel
  allocate(id(nphys2),ele2model(nphys2))
  do i=1,nphys2
   read(1,'(2i10)') id(i),ele2model(i)
  end do
  close(1)

 !#[2]## read model
 allocate(rho(nmodel))
 open(1,file=modelfile,status='old',err=80) ! 20127.12.21
  read(1,*,err=81) nmodel2
  if ( nmodel .ne. nmodel2 ) then
   write(*,*) "GEGEGE nmodel",nmodel,"nmodel2",nmodel2
   stop
  end if
  do i=1,nmodel
   read(1,*) rho(i)
  end do
 close(1)

 !#[3]## gen r_cond
 allocate(r_cond%rho(  nphys2) )
 allocate(r_cond%sigma(nphys2) )
 allocate(r_cond%index(nphys2)) ! 2018.03.20
 r_cond%nphys2 = nphys2
 r_cond%sigma  = -9999
 do i=1,nphys2
  r_cond%rho(i)=rho(ele2model(i))
  if (abs(r_cond%rho(i)) .gt. 1.d-10 ) r_cond%sigma(i) = 1.d0/r_cond%rho(i)
 end do

 return
  90 continue
  write(*,*) "File is not exist",connectfile
  stop
  80 continue
  write(*,*) "File is not exist",modelfile
  stop
  81 continue
  write(*,*) "File is not exist",modelfile,"line",i
  stop
 end




end module
