! modified so that the outmost models are combined 2017.07.20
! coded on May 10, 2017
module modelpart
use param
use matrix

! parameters for partition
type modelpara
 character(50) :: mshfile
 ! 0: no combining
 ! 1:combining outside models and not fix the resistivity ,2017.09.20
 ! 2:combining outside and fix the resistivity by the value of reference 2018.03.16
 integer(4)    :: icombine
 integer(4)    :: nxdiv,nydiv,nzdiv
 real(8),allocatable,dimension(:) :: xdiv
 real(8),allocatable,dimension(:) :: ydiv
 real(8),allocatable,dimension(:) :: zdiv
end type

type model
 integer(4)                          :: nmodel    ! # of models for inversion
 integer(4)                          :: nphys1    ! # of elements within air
 integer(4)                          :: nphys2    ! # of elements within land
 integer(4)                          :: nmodelactive!# of model for inversion 2018.03.18
 integer(4),allocatable,dimension(:) :: iactive   ! [nmodel] 1:active,0:innactive 20180318
 integer(4),allocatable,dimension(:) :: index     ! [nphys2] added on May 13, 2017
 integer(4),allocatable,dimension(:) :: ele2model ! [nphys2] only nphysinv is for inv
 real(8),allocatable,dimension(:)    :: rho_model ! [nmodel]
 real(8),allocatable,dimension(:)    :: logrho_model ! [nmodel]
 type(real_crs_matrix)               :: model2ele    ! id for global element space
 ! only for inv_mdm 2017.10.31
 integer(4),allocatable,dimension(:) :: dm2modelptr  ! [1:nmodel] id for mother model
 !# for icombine=2
 integer(4) :: icombine   ! 2018.03.17 0:no combine,1:combine outside 2:combine and fix
end type

contains

!###################################################
! coded on 2017.05.10
subroutine readmodelpara(idev,g_modelpara)
implicit none
type(modelpara),intent(out) :: g_modelpara
integer(4),intent(in) :: idev
integer(4) :: i,j,nxdiv,nydiv,nzdiv

!open(idev,file=ifile)
 write(*,*) "" !2020.09.29
 write(*,40) " < Input icombine >"
 write(*,40) " < icombine = 0: not combine the outside model >"
 write(*,40) " < icombine = 1: combine the outside models and make it variable through inversion >"! 2020.09.29
 write(*,40) " < icombine = 2: combine the outside model and fix it to the initial model >"    ! 2020.09.29
 read(idev,'(20x,i5)') g_modelpara%icombine              ! 2017.09.20
 write(*,42) " icombine =", g_modelpara%icombine ! 2020.09.29
!# xdiv
 write(*,*) "" !2020.09.29
 write(*,*) "< Input nxdiv for partition of model space in x direction>" ! 2020.09.29
 read(idev,*) nxdiv
 g_modelpara%nxdiv = nxdiv
 write(*,42) " nxdiv =",g_modelpara%nxdiv  ! 2020.09.29
 allocate(g_modelpara%xdiv(nxdiv))
 write(*,*) "" !2020.09.29
 write(*,*) "< Input nxdiv x coordinates>" ! 2020.09.29
 read(idev,'(g15.7)') (g_modelpara%xdiv(i),i=1,nxdiv)
 write(*,43) (g_modelpara%xdiv(i),i=1,nxdiv)

! write(*,'(g15.7)')   (g_modelpara%xdiv(i),i=1,nxdiv) commented out 2017.12.25
!# ydiv
 write(*,*) "" !2020.09.29
 write(*,*) "< Input nxdiv for partition of model space in y direction>" ! 2020.09.29
 read(idev,*) nydiv
 g_modelpara%nydiv = nydiv
 write(*,42) " nydiv =",g_modelpara%nydiv   ! 2020.09.29
 allocate(g_modelpara%ydiv(nydiv))
 write(*,*) "" !2020.09.29
 write(*,*) "< Input nydiv y coordinates>" ! 2020.09.29
 read(idev,'(g15.7)') (g_modelpara%ydiv(i),i=1,nydiv)
 write(*,43) (g_modelpara%ydiv(i),i=1,nydiv)
! write(*,'(g15.7)')   (g_modelpara%ydiv(i),i=1,nydiv) commented out 2017.12.25
!# zdiv
 write(*,*) "" !2020.09.29
 write(*,*) "< Input nzdiv for partition of model space in z direction>" ! 2020.09.29
 read(idev,*) nzdiv
 g_modelpara%nzdiv = nzdiv
 write(*,42) " nzdiv =",g_modelpara%nzdiv   ! 2020.09.29
 allocate(g_modelpara%zdiv(nzdiv))
 write(*,*) "" !2020.09.29
 write(*,*) "< Input nzdiv z coordinates>" ! 2020.09.29
 read(idev,'(g15.7)') (g_modelpara%zdiv(i),i=1,nzdiv)
 write(*,43) (g_modelpara%zdiv(i),i=1,nzdiv)
! write(*,'(g15.7)')   (g_modelpara%zdiv(i),i=1,nzdiv) commented out 2017.12.25
!close(idev)
40 format(a)! 2020.09.29
41 format(a,a)
42 format(a,i3)
43 format(10f9.4)
return
end subroutine

!#################################################### genmodelspace
! coded on 2017.05.10
subroutine genmodelspace(g_mesh,g_modelpara,g_model,g_param,g_cond)
use mesh_type
use param
implicit none
type(mesh),             intent(in)     :: g_mesh
type(modelpara),        intent(in)     :: g_modelpara
type(model),            intent(out)    :: g_model
type(param_forward),    intent(in)     :: g_param
type(param_cond),       intent(in)     :: g_cond
real(8),   allocatable, dimension(:)   :: xdiv,ydiv,zdiv
real(8),   allocatable, dimension(:,:) :: xyz
integer(4),allocatable, dimension(:,:) :: n4
integer(4),allocatable, dimension(:)   :: ele2model,elecount_model,model2model
integer(4),allocatable, dimension(:)   :: elecount_model2 ! 2017.05.15
integer(4),allocatable, dimension(:)   :: stack,item      ! 2017.05.15
real(8)    :: xyzminmax(6),xyzcen(3)
integer(4) :: i,j,k,ie, i1,j1,j2,k1,iele,ii
integer(4) :: nmodel,nxdiv,nydiv,nzdiv
integer(4) :: nphys2,ntet,node,icount,imodel,nphys1
integer(4) :: icombine                   ! 2017.09.20
integer(4),allocatable,dimension(:)    :: index
type(real_crs_matrix) :: crsout

!#[1]## set input
 xyzminmax = g_param%xyzminmax
 nxdiv     = g_modelpara%nxdiv
 nydiv     = g_modelpara%nydiv
 nzdiv     = g_modelpara%nzdiv
 icombine  = g_modelpara%icombine      ! 2017.09.20
 write(*,*) "check1"
 write(*,*) "nxdiv,nydiv,nzdiv=",nxdiv,nydiv,nzdiv
 allocate(xdiv(nxdiv + 2),ydiv(nydiv + 2),zdiv(nzdiv + 2))
 xdiv(2:nxdiv+1) = g_modelpara%xdiv
 ydiv(2:nydiv+1) = g_modelpara%ydiv
 zdiv(2:nzdiv+1) = g_modelpara%zdiv
 xdiv(1) = xyzminmax(1) ;  xdiv(nxdiv+2) = xyzminmax(2)
 ydiv(1) = xyzminmax(3) ;  ydiv(nydiv+2) = xyzminmax(4)
 zdiv(1) = xyzminmax(5) ;  zdiv(nzdiv+2) = xyzminmax(6)
 if (      icombine .eq. 0 ) then       ! no combine 2017.09.20
  nmodel = (nxdiv + 1)*(nydiv +1 )*(nzdiv + 1)  ! commented out on 2017.09.20
 else if ( icombine .eq. 1 .or. icombine .eq. 2 ) then ! combine outside model 2018.03.16
  nmodel = (nxdiv - 1)*(nydiv - 1 )*(nzdiv - 1) + 1   ! 2017.09.20
 else
  write(*,*) "GEGEGE! stop!"                    ! 2017.09.20
 end if
 write(*,*) "# of model parameters: nmodel=",nmodel
!
 ntet = g_mesh%ntet
 node = g_mesh%node
 allocate(xyz(3,node),n4(ntet,4))
 xyz = g_mesh%xyz
 n4  = g_mesh%n4
 nphys2 = g_cond%nphys2
 nphys1 = g_cond%nphys1
 allocate(ele2model(nphys2),index(nphys2))
 index = g_cond%index
 ele2model(:)=-1 ! default with not assigned to any model  2018.03.16

!#[2]## make partition
do ie=1,nphys2
 xyzcen = 0.d0
 iele = index(ie)
 do i=1,4
  xyzcen(1:3) = xyzcen(1:3) + xyz(1:3,n4(iele,i))/4.d0
 end do
 do k=1,nzdiv+1
  if ( zdiv(k) .lt. xyzcen(3) .and. xyzcen(3) .le. zdiv(k+1)) then
  do j=1,nydiv+1
   if ( ydiv(j) .lt. xyzcen(2) .and. xyzcen(2) .le. ydiv(j+1)) then
   do i=1,nxdiv+1
    if ( xdiv(i) .lt. xyzcen(1) .and. xyzcen(1) .le. xdiv(i+1)) then
     k1 = k ; j1 = j ; i1 = i
     goto 100
    end if
   end do
   end if
  end do
  end if
 end do
 write(*,*) "GEGEGE iele=",iele,"xyzcen=",xyzcen(1:3)
 stop
 100 continue
 if ( icombine .eq. 1 .or. icombine .eq. 2 ) then    ! 2018.03.16
  if ( k1 .eq. 1 .or. k1 .eq. nzdiv+1 .or. j1 .eq. 1 .or. j1 .eq. nydiv+1 .or.&
   &    i1 .eq. 1 .or. i1 .eq. nxdiv+1 ) then        ! 2017.07.20
   ele2model(ie) = (nydiv-1)*(nxdiv-1)*(nzdiv-1) + 1
  else   ! 2017.07.20
    ele2model(ie) = (k1-2)*(nydiv-1)*(nxdiv-1) + (j1-2)*(nxdiv-1) + i1-1! 2017.07.20
  end if
 else if (icombine .eq. 0 ) then ! 2017.09.20
  ele2model(ie) = (k1-1)*(nydiv+1)*(nxdiv+1) + (j1-1)*(nxdiv+1) + i1  ! 2017.09.20
 end if ! 2017.09.20
! write(*,*) "ie=",ie,"index(ie)=",index(ie),"ele2model(ie)=",ele2model(ie)
! write(*,*) "xyzcen=",xyzcen
end do

!#[2]## exclude model with no elements
 allocate(elecount_model(nmodel),model2model(nmodel))
 allocate(elecount_model2(nmodel))
 elecount_model   = 0
 do i=1,nphys2
  imodel = ele2model(i)
  elecount_model(imodel) =  elecount_model(imodel) + 1
 end do
 icount=0
 do i=1,nmodel
  if (elecount_model(i) .ne. 0) then
   icount=icount+1
   model2model(i) = icount
   elecount_model2(icount) = elecount_model(i) ! added on 2017.05.15
!   write(*,*) "icount=",icount,"# of ele=",elecount_model2(icount)
  end if
 end do
 do i=1,nphys2
  ele2model(i)=model2model(ele2model(i))
 end do
 write(*,*) "nmodel",nmodel,"->",icount
 nmodel = icount

!#[3]## model2ele (crsmatrix) ! added on 2017.05.15
 allocate(stack(0:nmodel),item(nphys2))
 stack(0)=0
 do i=1,nmodel
  stack(i) = stack(i-1) + elecount_model2(i)
 end do
 if (stack(nmodel) .ne. nphys2 ) then
  write(*,*) "GEGEGE stack(nmodel)=",stack(nmodel),"nphys2=",nphys2
  stop
 end if
 elecount_model = 0
 do i=1,nphys2
   ii = ele2model(i)
   elecount_model(ii) = elecount_model(ii) + 1
   item(stack(ii-1)+elecount_model(ii)) = index(i) ! element id for whole element space
 end do
!# cehck
! do i=1,nphys2
!  if (ele2model(i) .eq. 5) write(*,*) i,"ele2model=",ele2model(i)
! end do
! do i=5,5 ! nmodel
!  write(*,*) i,"elecount_model2=",elecount_model2(i)
! end do
! do i=5,5 ! nmodel
!  write(*,*) "model",i,"elements:",(item(j),j=stack(i-1)+1,stack(i))
! end do
! stop

!#[4]## output
 g_model%nmodel           = nmodel
 allocate(g_model%ele2model(nphys2))
 g_model%ele2model        = ele2model
 g_model%nphys2           = nphys2
 g_model%nphys1           = nphys1
 allocate(g_model%index(nphys2))
 g_model%index            = index
 allocate(g_model%model2ele%stack(0:nmodel) )
 allocate(g_model%model2ele%item(nphys2))
 g_model%model2ele%nrow   = nmodel       ! 2017.06.05
 g_model%model2ele%ncolm  = ntet         ! 2017.06.05 max(item) = ntet
 g_model%model2ele%ntot   = nphys2
 g_model%model2ele%stack  = stack
 g_model%model2ele%item   = item         ! id for element in whole element space
 g_model%icombine         = icombine     ! 2018.03.19
 !#
 allocate(g_model%model2ele%val(nphys2)) ! 2018.03.16
 g_model%model2ele%val(:) = 1.d0         ! 2017.06.05
 !# assume some model is not used in inversion 2018.03.18
 allocate(g_model%iactive(nmodel))       ! 2018.03.18
 g_model%iactive      = 1                ! 2018.03.18
 g_model%nmodelactive = nmodel           ! 2018.03.18
 !# when icombine = 2, remain nmodel-th model as initial value 2018.06.25
 if (icombine .eq. 2) then               ! 2018.03.18
  g_model%iactive(nmodel) = 0            ! 2018.03.18
  g_model%nmodelactive    = nmodel - 1   ! 2018.03.18
 end if

if (.false.) then
 crsout = g_model%model2ele
 write(*,*) "crs%nrow=",crsout%nrow
 write(*,*) "crs%ncolm=",crsout%ncolm
 write(*,*) "crs%ntot=",crsout%ntot
 do i=1,crsout%nrow
  if ( i .eq. 3120) then
 if (crsout%stack(i)-crsout%stack(i-1) .ne. 0) then
!  write(*,*) i,"# of content",crsout%stack(i)-crsout%stack(i-1)!,"ele2model=",ele2model(i)
  j1=crsout%stack(i-1)+1;j2=crsout%stack(i)
  write(*,'(5g15.7)') (crsout%item(j),j=j1,j2)
! write(*,'(5g15.7)') (crsout%val(j),j=j1,j2)
  write(*,*) ""
  end if
 end if
 end do
end if

return
end subroutine
!##################################################### modelparam
! coded on 2017.08.28
! to show model parameter space
subroutine modelparam(g_model)
implicit none
type(model),intent(inout) :: g_model
integer(4) :: nmodel, i
real(8)    :: z

!#[1]## set
nmodel = g_model%nmodel

!#[2]# output
if ( .not. allocated(g_model%logrho_model) ) allocate(g_model%logrho_model(nmodel))
call random_number(g_model%logrho_model(1:nmodel))
do i=1,nmodel
 !# 0 - 1 -> 1 - 5
 g_model%logrho_model(i) = g_model%logrho_model(i) * 4. + 1.
end do

return
end subroutine
!##################################################### assignmodelrho
! coded on 2017.05.10
! averaged rho is assigned to each model region
subroutine assignmodelrho(g_cond,g_model)
use mesh_type
use param
implicit none
type(model),            intent(inout) :: g_model
type(param_cond),       intent(in)    :: g_cond
real(8),   allocatable, dimension(:)  :: rho_model,logrho_model
integer(4),allocatable, dimension(:)  :: icount_model,ele2model
real(8),   allocatable, dimension(:)  :: rho
integer(4),allocatable, dimension(:)  :: index
integer(4) :: i, imodel,nphys2,nmodel

!#[1]##
 nmodel    = g_model%nmodel
 nphys2    = g_cond%nphys2
 allocate(ele2model(nphys2),rho(nphys2))
 allocate(index(nphys2))
 ele2model = g_model%ele2model
 rho       = g_cond%rho
 index     = g_cond%nphys2

!#[2]##
 allocate(rho_model(nmodel),logrho_model(nmodel),icount_model(nmodel))
 icount_model(:)=0 ; rho_model=0.d0 ; logrho_model=0.d0
 do i=1,nphys2
  imodel= ele2model(i)
  icount_model(imodel) =  icount_model(imodel) + 1
  logrho_model(imodel) = logrho_model(imodel) + log10(rho(i))
 end do

!#[3]## average
  do imodel=1,nmodel
   if (icount_model(imodel) .ne. 0) then
    logrho_model(imodel)    =  logrho_model(imodel)/dble(icount_model(imodel))
    rho_model(imodel) = 10**(logrho_model(imodel))
!    write(*,*) "i=",imodel,"rho_model=",rho_model(imodel)
   end if
  end do

!#[4]## output
   if ( allocated(g_model%rho_model)    ) deallocate(g_model%rho_model   )!2017.07.19
   if ( allocated(g_model%logrho_model) ) deallocate(g_model%logrho_model)!2017.07.19
   allocate(g_model%rho_model(   nmodel) )   !
   allocate(g_model%logrho_model(nmodel) )
   g_model%rho_model    = rho_model
   g_model%logrho_model = logrho_model

return
end subroutine

!######################################################  model2cond
! modified on 2017.05.18
subroutine model2cond(g_model,g_cond,iflag) ! 2018.11.08
implicit none
type(model),           intent(in)    :: g_model
type(param_cond),      intent(inout) :: g_cond
integer(4)                           :: nphys1,nphys2, ntet, nmodel,i
integer(4),allocatable,dimension(:)  :: ele2model
real(8),   allocatable,dimension(:)  :: logrho_model
real(8)                              :: rho,sigma_air
character(50)                        :: condfile
type(param_cond)                     :: h_cond
!   0: normal, 1:g_cond%rho=g_model%logrho_model
integer(4),intent(in),optional       :: iflag

!#[1]## check and deallocate
if (allocated(g_cond%rho) ) call deallocatecond(g_cond) ! see m_param.f90

!#[2]## set
nmodel           = g_model%nmodel
nphys1           = g_model%nphys1
nphys2           = g_model%nphys2
ntet             = nphys1 + nphys2
allocate(ele2model(nphys2),logrho_model(nmodel))
ele2model        = g_model%ele2model
logrho_model     = g_model%logrho_model   ! logrho is input for g_cond
sigma_air        = g_cond%sigma_air       ! added on 2017.05.18
condfile         = g_cond%condfile        ! added on 2017.05.18

!#[3]## gen new g_cond
allocate(h_cond%index(nphys2),h_cond%rho(nphys2),h_cond%sigma(nphys2))
h_cond%nphys1    = nphys1
h_cond%nphys2    = nphys2
h_cond%ntet      = ntet
h_cond%index     = g_model%index
h_cond%condflag  = 1         ! 2017.05.18 (not homogeneous solid earth)
h_cond%condfile  = condfile  ! 2017.05.18
h_cond%sigma_air = sigma_air ! 2017.05.18
write(*,*) size(logrho_model)
do i=1,nphys2
! write(*,*) "i",i,ele2model(i),ele2model(i),logrho_model(ele2model(i))
 rho = 10.d0**logrho_model(ele2model(i))
 h_cond%rho(i) = rho
 h_cond%sigma(i) = 1.d0/rho
 if (present(iflag) .and. iflag .eq. 1 ) then         ! 2018.11.08
  rho = logrho_model(ele2model(i)) ! 2018.11.08
  h_cond%rho(i)   = rho            ! 2018.11.08
  h_cond%sigma(i) = 0.d0           ! 2018.11.08
 end if                            ! 2018.11.08
end do

!#[4]## set output
  g_cond = h_cond

return
end subroutine

end module
