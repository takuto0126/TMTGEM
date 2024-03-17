!# adopted "precision.inc" for len_file on 2016.11.15
!# coded on 2016.09.14
module oceanFvxyz
use mesh_type
use param      ! include "precision.inc"
use param_mesh
use FROMCOMCOT
use constants
use solution_ana ! for OUTANA, see m_solution_ana.f90
use matrix     ! 2017.11.02
use triangle   ! 2017.11.02
implicit none

! ocean data is only for extrude mesh
type ocean_data
 integer(4) :: node   ! # of all the nodes
 integer(4) :: nodes  ! # of ocean nodes
 integer(4) :: IPL    ! # of ocean surface nodes
 integer(4),allocatable,dimension(:)   :: nz
 integer(4),allocatable,dimension(:)   :: surfptr ! [IPL] node id on ocean surface in 3-D mesh space
 integer(4),allocatable,dimension(:)   :: ocean2emptr ! [nodes] node id in 3-D em mesh
 integer(4),allocatable,dimension(:)   :: em2oceanptr ! [node]  node id in 3-D ocean mesh
 real(8),   allocatable,dimension(:,:) :: Fxyz ! (3,npoints_in_ocean) 2021.06.24
 real(8),   allocatable,dimension(:,:) :: vxyz ! (3,npoints_in_ocean) 2021.06.24
 real(8),   allocatable,dimension(:)   :: z_tri! z_tri(IPL) 2021.12.08 wave height
 !# for div error 2017.11.02
 integer(4)            :: ntets            ! number of ocean elements         2017.11.02
 integer(4)            :: ntris            ! number of surface triangles      2017.11.02
 type(real_crs_matrix) :: oceanele2column  ! nrow: # of triangles, ntot: nele 2017.11.02
 real(8),allocatable,dimension(:,:) :: xy_tri ! [2,ntris] center coord of triangles 2017.11.02
end type ocean_data

type comcot_data
 integer(4)          :: iflag_comcot_spherical ! comcot is conducted by 0:cartesian, 1:spherical 2021.12.08
 integer(4)          :: nlayer ! # of layers of nested grid 2021.12.02
 integer(4)          :: IPL
 integer(4)          :: nodes  ! # of ocean nodes in 3-D mesh
 !# time step info             !      2021.12.02
 integer(4)          :: it_vxyz1, it_vxyz2
 real(8)             :: t_vxyz1,   t_vxyz2
 real(8),   allocatable, dimension(:,:) :: vxyz1,vxyz2   ! layer integrated velocity for nodes
 real(8),   allocatable, dimension(:)   :: z_tri1,z_tri2 ! layer integrated z for IPL 2021.12.08
 type(layer_data),allocatable,dimension(:)   :: layer  ! # of comcot layers 2021.12.06
end type comcot_data

type layer_data                ! 2021.12.06
 integer(4)         :: nx,ny
 type(grid_data_2D) :: h_grd !(nlayer) see m_FROMCOMCOT.f90 2021.12.02
 integer(4),allocatable, dimension(:,:) :: ivxyzflagkl  !inserted on Oct. 25, 2015
 real(8),   allocatable, dimension(:,:) :: vxyzcoef     !inserted on Oct. 25, 2015
end type layer_data

contains


!########################################## prepoceanele2column
!# coded on 2017.11.02
subroutine prepoceanele2column(g_param,em_mesh,h_ocean)
implicit none
type(param_forward),   intent(in)     :: g_param
type(mesh),            intent(in)     :: em_mesh
type(ocean_data),      intent(inout)  :: h_ocean
type(mesh)                            :: ki_mesh
integer(4)                            :: ntets
integer(4)                            :: nodek,ntrik
integer(4)                            :: node, ntet
real(8),   allocatable,dimension(:,:) :: xyz
real(8),   allocatable,dimension(:,:) :: xyzk
integer(4),allocatable,dimension(:,:) :: n3k
integer(4),allocatable,dimension(:,:) :: n4
integer(4),allocatable,dimension(:,:) :: n4flag
integer(4)                            :: ndat
type(real_crs_matrix)                 :: oceanele2column
integer(4)                            :: nx,ny
real(8)                               :: xyzminmax(6),a3(3),xyobs(2),ele_xyz(3,4)
type(grid_list_type)                  :: glist
integer(4)                            :: i,j,iele,ii, ishift, icount_max
integer(4),allocatable,dimension(:,:) :: band_ind
integer(4),allocatable,dimension(:)   :: icount
real(8),   allocatable,dimension(:,:) :: xy_tri

!#[0]## read ki_mesh and set input
   xyzminmax = g_param%xyzminmax
   node      = em_mesh%node
   allocate( xyz(3,node) )
   xyz       = em_mesh%xyz
   ntet      = em_mesh%ntet
   allocate( n4(ntet,4),n4flag(ntet,2) )
   n4        = em_mesh%n4
   n4flag    = em_mesh%n4flag
   !
   CALL READMESH_TOTAL(ki_mesh, g_param%pokimsh) ! ki_mesh includes both land and ocean 2d mesh
   nodek     = ki_mesh%node
   ntrik     = ki_mesh%ntri
   allocate( xyzk(3,nodek), n3k(ntrik,3)   )
   xyzk      = ki_mesh%xyz
   n3k       = ki_mesh%n3

!#[2]## Prepare list
  nx=300;ny=300
  CALL allocate_2Dgrid_list(nx,ny,ntrik,glist)   ! see m_mesh_type.f90
  CALL gen2Dgridforlist(xyzminmax,glist) ! see m_mesh_type.f90
  CALL classifytri2grd(ki_mesh,glist)   ! classify ele to glist,see

!#[3]## 
  allocate( band_ind(50,ntrik) )
  allocate( icount(ntrik)      )
  icount(:)  = 0
  ntets      = 0
  do i=1,ntet
   if ( n4flag(i,2) .ne. 2 ) cycle ! 1:air, 2:ocean, 3:earth
   ntets = ntets + 1
   do j=1,4
    ele_xyz(1:3,j) = xyz(1:3,n4(i,j))
   end do
   xyobs(1) = sum(ele_xyz(1,:))/4.d0 ! x of center gravity
   xyobs(2) = sum(ele_xyz(2,:))/4.d0 ! y of center gravity
   call findtriwithgrid(ki_mesh,glist,xyobs(1:2),iele,a3) ! iele is element id in ki_mesh
   !
   icount(iele) = icount(iele) + 1
!   write(*,'(a,i10,a,i10,a,i10,a,2f15.7)') "i",i,"iele",iele,"icount(iele)",icount(iele),"xy",xyobs(1:2)
!   write(*,'(2g15.7)') (xyzk(1:2,n3k(iele,j)),j=1,3)
   band_ind(icount(iele),iele) = i
  end do
  write(*,*) "ntet =",ntet
  write(*,*) "ntets=",ntets


!#[2]## sort ocean element to column
  oceanele2column%nrow = ntrik
  oceanele2column%ntot = ntets
  allocate( oceanele2column%stack(0:ntrik) )
  allocate( oceanele2column%item( ntets  ) )
  allocate( oceanele2column%val(  ntets  ) )
  oceanele2column%stack(:) = 0
  icount_max = 0
  allocate( xy_tri(2,ntrik) )
  xy_tri(:,:) = 0.d0
  do i=1,ntrik
   ishift = oceanele2column%stack(i-1)
   oceanele2column%stack(i) = ishift + icount(i)
   do j=1,icount(i)
    oceanele2column%item(ishift + j) = band_ind(j,i)
    oceanele2column%val( ishift + j) = band_ind(j,i)*1.d0
   end do
   icount_max = max(icount_max,icount(i))
   ! xy_tri
   do j=1,3
    xy_tri(1,i) = xy_tri(1,i) + xyzk(1,n3k(i,j))/3.d0
    xy_tri(2,i) = xy_tri(2,i) + xyzk(2,n3k(i,j))/3.d0
   end do
  end do
  oceanele2column%ncolm = icount_max

!#[3]## set output
  h_ocean%ntets           = ntets
  h_ocean%ntris           = ntrik
  h_ocean%oceanele2column = oceanele2column
  allocate(h_ocean%xy_tri(2,ntrik))
  h_ocean%xy_tri          = xy_tri

 write(*,*) "### PREPOCEANELE2COLUMN END!! ###"
 return

!# Error
 99 continue
 write(*,*) "GEGEGE! ndat=",ndat,"is not equal to nodek",nodek
 stop

end
!########################################## prepareFxyz
subroutine prepareFxyz(h_ocean,em_mesh,g_param,g_meshpara)
implicit none
type(mesh),         intent(in) :: em_mesh
type(param_forward),intent(in) :: g_param
type(meshpara),     intent(in) :: g_meshpara
type(ocean_data),intent(inout) :: h_ocean
integer(4)                     :: nodes,iflag_geomag
integer(4)                     :: i,j,node,inod,unit ! 2019.01.25

!#[0]## set input
!# iflag_geomag =
!# 0 : IGRF
!# 1 : file
!# 2 : constant
iflag_geomag  = g_param%iflag_geomag

!#[1]## for Fxyz

 !#[1-0]##  iflag_geomag = 0 for IGRF value
 if      ( iflag_geomag .eq. 0  ) then                ! 2018.11.14
  call geomagigrf(h_ocean,em_mesh,g_param,g_meshpara) ! 2018.11.14

 !#[1-1]##  iflag_geomag = 1 for file value
 else if      (iflag_geomag .eq. 1 ) then
  call geomagfile(h_ocean,em_mesh,g_param,g_meshpara)

 !#[1-2]## iflag_geomag = 2 for homogeneous value
 else if (iflag_geomag .eq. 2 ) then
  call geomaghomo(h_ocean,g_param)

 end if

!#[2]## OUTPUT prepared geomag field 2019.01.25
 open(file=trim(g_param%outbxyzfolder)//"geomag.out",newunit=unit) ! 2019.01.25
  do i=1,h_ocean%IPL
    j=h_ocean%surfptr(i)             ! 2019.01.25
    write(unit,'(5g15.7)') em_mesh%xyz(1:2,j),h_ocean%Fxyz(1:3,j) ! 2019.01.25
  end do
 close(unit)


return
end

!########################################## preparevxyz
subroutine preparevxyz(h_ocean,em_mesh,g_param,g_meshpara,it,dt,h_comcot)
implicit none
type(mesh),         intent(in)        :: em_mesh
type(param_forward),intent(in)        :: g_param
type(meshpara),     intent(in)        :: g_meshpara
integer(4),         intent(in)        :: it
real(8),            intent(in)        :: dt
type(ocean_data),   intent(inout)     :: h_ocean  ! include vxyz,z_tri ! z_tri added 2021.12.08
type(comcot_data),  intent(inout)     :: h_comcot
type(grid_data_2D),allocatable,dimension(:)  :: vxyh_grd ! 2021.12.02
type(grid_data_2D),allocatable,dimension(:)  ::    z_grd ! 2021.12.02
integer(4)                            :: iflag_velocity, nodes
integer(4),allocatable,dimension(:)   :: ocean2emptr
real(8),   allocatable,dimension(:,:) :: vxyz1, vxyz2
real(8),   allocatable,dimension(:)   :: z_tri1,z_tri2 ! sea surface elevation 2021.12.08
real(8)    :: t  ! [sec]
real(8)    :: dt_vxyz ![sec]
real(8)    :: r,r1
integer(4) :: it_vxyz1, it_vxyz2, i,IPL,ii ! 2021.12.08 IPL is added
integer(4) :: nlayer,nx,ny                    ! 2021.12.02
real(8)    :: t_vxyz1, t_vxyz2

!#[0]## set input
  iflag_velocity = g_param%iflag_velocity
  nodes          = h_ocean%nodes
  IPL            = h_ocean%IPL
  allocate( ocean2emptr(nodes))
  allocate( vxyz1(3,nodes),vxyz2(3,nodes))
  allocate( z_tri1(IPL), z_tri2(IPL) )           ! 2021.12.08
  ocean2emptr    = h_ocean%ocean2emptr
  nlayer         = g_param%nlayer  ! 2021.12.02
  allocate( vxyh_grd(nlayer) )     ! 2021.12.02
  allocate(    z_grd(nlayer) )     ! 2021.12.02

!### for vxyz

 !#[1]## iflag_velocity = 1 for COMCOT results
 if      (iflag_velocity .eq. 1 ) then

  !#[1-1]## coefficient preparation for COMCOT files (before it iterations)
  if      ( it .eq. 0 ) then
     call gencomcotcoef(h_ocean,em_mesh,g_param,g_meshpara,h_comcot) ! see below

     !# initialize h_comcot
     h_comcot%t_vxyz1  = 0.d0
     h_comcot%t_vxyz2  = 0.d0
     h_comcot%it_vxyz1 = 0      ! it_vxyz1 start with 0
     h_comcot%it_vxyz2 = 0      ! it_vxyz2 start with 0
     h_comcot%vxyz1    = 0.d0   ! comcot layer integrated velocity for nodes
     h_comcot%vxyz2    = 0.d0   ! comcot layer integrated velocity for nodes
     h_comcot%z_tri1   = 0.d0   ! layer integrated surface displacement for IPL 2021.12.08
     h_comcot%z_tri2   = 0.d0   ! layer integrated surface displacement for IPL 2021.12.08

  !#[1-2]## calculation of vxyz from COMCOT files
  else if ( it .ge. 1 ) then

     !#[1-2-1]## set input from h_comcot
      t        = it*dt
      t_vxyz1  = h_comcot%t_vxyz1
      t_vxyz2  = h_comcot%t_vxyz2
     it_vxyz1  = h_comcot%it_vxyz1
     it_vxyz2  = h_comcot%it_vxyz2
     dt_vxyz   = g_param%dt_comcot
        vxyz1  = h_comcot%vxyz1
        vxyz2  = h_comcot%vxyz2
        z_tri1 = h_comcot%z_tri1 ! 2021.12.08
        z_tri2 = h_comcot%z_tri2 ! 2021.12.08

     do i=1,nlayer ! 2021.12.02 layer loop
      CALL ALLOCATEGRDDATA(vxyh_grd(i),h_comcot%layer(i)%h_grd%xygrid,2)! z_grd is added
      CALL ALLOCATEGRDDATA(   z_grd(i),h_comcot%layer(i)%h_grd%xygrid,1)! z_grd is added
     end do

     !#[1-2-2]## cal vxyz
     if ( t_vxyz2 .lt. t) then
      10 continue

      !# move comcot time forward
      it_vxyz1  = it_vxyz2
      it_vxyz2  = it_vxyz2 + 1
       t_vxyz1  = dt_vxyz*it_vxyz1
       t_vxyz2  = dt_vxyz*it_vxyz2

      if (  t_vxyz2 .lt. t ) goto 10

      !# after it gets t_vxyz1 < t <= t_vxyz2
      
	 vxyz1(:,:) = vxyz2(:,:)                   ! renew vxyz1
      z_tri1(:)  = z_tri2(:)                    ! renew z      2021.12.08
      CALL READMN(it_vxyz2,vxyh_grd,z_grd,g_param)!get vxyh_grd,z_grd, see m_FROMCOMCOT.f90 2021.12.06

!       write(*,*) "size(vxyh_grd(1)%dat(1,:,:))",size(vxyh_grd(1)%dat(1,:,:))
!       nx=vxyh_grd(1)%xygrid%nx
!       ny=vxyh_grd(1)%xygrid%ny
!       write(*,*) "nx,ny",nx,ny
!       call countnonzero3( vxyh_grd(1)%dat(1:2,:,:),2,nx,ny,ii)
!       write(*,*) "vxyz_grd(1)%dat(1:2,:,:) count",ii
!       call countnonzero3(    z_grd(1)%dat,1,nx,ny,ii)
!       write(*,*) "z_grd(1)%dat(1:1,:,:) count",ii

      CALL CALOCEANVXYZ(vxyh_grd,z_grd,h_ocean,h_comcot,em_mesh,g_meshpara,vxyz2,z_tri2) !2021.12.08 z_tri added
!      call countnonzero1( z_tri2,IPL,ii)
!      write(*,*) "z_tri2(1:IPL) count",ii
!      call countnonzero2( vxyz2,3,nodes,ii)
!      write(*,*) "vxyz2(1:3,nodes) count",ii

     end if
     r  = ( t - t_vxyz1)/(t_vxyz2 - t_vxyz1)
     r1 = 1.d0 -r
     do i=1,nodes
      h_ocean%vxyz(1:3,i) = r1*vxyz1(1:3,i) + r * vxyz2(1:3,i)
     end do
     do i=1,IPL  ! 2021.12.08
      h_ocean%z_tri(   i) = r1*z_tri1(i)    + r * z_tri2(   i) ! 2021.12.08
     end do      ! 2021.12.08

     !#[1-2-3]## select components of seawater velocity, 2017.07.04
     if (g_param%iflag_vF .eq. 1) then ! only vh
      write(*,*) "Only vh is used."
      h_ocean%vxyz(3,1:nodes)   = 0.d0
     else if (g_param%iflag_vF .eq. 2) then ! only vz
      write(*,*) "Only vz is used."
      h_ocean%vxyz(1:2,1:nodes) = 0.d0
     else if (g_param%iflag_vF .eq. 3) then ! both
      write(*,*) "Both vh and vz are used."
     else
      write(*,*) "GEGEGE g_param%iflag_vF should be 1, 2, or 3.",g_param%iflag_vF
	stop
     end if

     !#[1-2-4]## set output
     h_comcot%t_vxyz1  = t_vxyz1
     h_comcot%t_vxyz2  = t_vxyz2
     h_comcot%it_vxyz1 = it_vxyz1
     h_comcot%it_vxyz2 = it_vxyz2
     h_comcot%vxyz1    = vxyz1
     h_comcot%vxyz2    = vxyz2
     h_comcot%z_tri1   = z_tri1 ! 2021.12.08
     h_comcot%z_tri2   = z_tri2 ! 2021.12.08

  end if

 !#[2]## iflag_velocity = 2 for analytical velocity
 else if (iflag_velocity .eq. 2 ) then
  if ( it .eq. 0 ) then
   CALL OUTANA(g_param) ! see m_solution_ana.f90
  else if ( it .ge. 1 ) then
   t = it * dt
   call velocityana(h_ocean,em_mesh,g_param,t)
  end if
 end if

return
end
!################################################
! 2021.12.08
subroutine countnonzeroint2(array,n,m,ii)
implicit none
integer(4), intent(in)  :: n,m
integer(4), intent(in)  :: array(n,m)
integer(4), intent(out) :: ii
integer(4)              :: i,j,k
ii=0
do j=1,m;do i=1,n
if ( array(i,j) .ne. 0 ) ii=ii+1
end do;end do
return
end

!################################################
! 2021.12.08
subroutine countnonzero3(array,idim,n,m,ii)
implicit none
integer(4), intent(in)  :: n,m,idim
real(8),    intent(in)  :: array(idim,n,m)
integer(4), intent(out) :: ii
integer(4)              :: i,j,k
real(8)                 :: threshold=1.d-5
ii=0
do j=1,m;do i=1,n;do k=1,idim
if ( abs(array(k,i,j)) .gt. threshold ) ii=ii+1
end do;end do;end do
return
end
!################################################
! 2021.12.08
subroutine countnonzero2(array,n,m,ii)
implicit none
integer(4), intent(in)  :: n,m
real(8),    intent(in)  :: array(n,m)
integer(4), intent(out) :: ii
integer(4)              :: i,j,k
real(8)                 :: threshold=1.d-5
ii=0
do j=1,m;do i=1,n
if ( abs(array(i,j)) .gt. threshold ) ii=ii+1
end do;end do
return
end
!################################################
! 2021.12.08
subroutine countnonzero1(array,m,ii)
implicit none
integer(4), intent(in)  :: m
real(8),    intent(in)  :: array(m)
integer(4), intent(out) :: ii
integer(4)              :: i,j,k
real(8)                 :: threshold=1.d-5
ii=0
do j=1,m
if ( abs(array(j)) .gt. threshold ) ii=ii+1
end do
return
end
!########################################## calcoeffcomcot
subroutine gencomcotcoef(h_ocean,em_mesh,g_param,g_meshpara,h_comcot)
implicit none
type(mesh),         intent(in)    :: em_mesh
type(param_forward),intent(in)    :: g_param
type(meshpara),     intent(in)    :: g_meshpara
type(ocean_data),   intent(in)    :: h_ocean
type(comcot_data),  intent(out)   :: h_comcot
type(grid_xy),      allocatable,dimension(:) :: xygrid ! 2021.12.02
type(grid_data_2D), allocatable,dimension(:) :: vxyh_grd,h_grd ! 2021.12.02
character(70),      allocatable,dimension(:) :: xgrdfile, ygrdfile,hfile
integer(4)         :: nodes,IPL! IPL is the numver of surface ocean nodes
integer(4)         :: nlayer,i,j ! # of COMCOT layer 2021.12.02
integer(4),allocatable, dimension(:,:,:) :: ivxyzflagkl  !2021.12.08
real(8),   allocatable, dimension(:,:,:) :: vxyzcoef     !2021.12.08

!#[0]## set input
 IPL       = em_mesh%npoi
 nodes     = h_ocean%nodes
 nlayer    = g_param%nlayer         ! 2021.12.02
 allocate( vxyh_grd(nlayer) )       ! 2021.12.08
 allocate(    h_grd(nlayer) )       ! 2021.12.08
 allocate( xgrdfile(nlayer) )       ! 2021.12.02
 allocate( ygrdfile(nlayer) )       ! 2021.12.02
 allocate( hfile(   nlayer) )       ! 2021.12.02
 allocate( xygrid(  nlayer) )       ! 2021.12.08
 xgrdfile(1:nlayer) = g_param%xgrdfile_v(1:nlayer)
 ygrdfile(1:nlayer) = g_param%ygrdfile_v(1:nlayer)
 hfile(   1:nlayer) = g_param%hgrdfile(  1:nlayer)
 allocate( ivxyzflagkl(IPL,3,nlayer)) ! 2021.12.08
 allocate( vxyzcoef(IPL,2,nlayer))    ! 2021.12.08

 do i =1,nlayer  ! 2021.12.02
  !#[0]## read COMCOT grid data and Fx, Fy, Fz, added on Oct. 25, 2015
  CALL COUNTCOMCOT(  xygrid(i), xgrdfile(i), ygrdfile(i)) ! get nx, ny  < FROMCOMCOT.f90 2021.12.08
  CALL READCOMCOTGRD(xygrid(i), xgrdfile(i), ygrdfile(i),g_param,g_meshpara) ! 2021.12.08

  !#[2]## allocate vxyh_grd and h_grd
  CALL ALLOCATEGRDDATA(   vxyh_grd(i), xygrid(i),2)  ! 2021.12.02     ! see FROMCOMCOT.f90
  CALL ALLOCATEGRDDATA(      h_grd(i), xygrid(i),1)  ! 2021.12.02
  CALL READGRDDATA(hfile(i), h_grd(i), 1)  ! 2021.12.02 see FROMCOMCOT.f90

  !#[3]##
  CALL MKVXYZCOEF(xygrid(i),h_ocean,em_mesh,ivxyzflagkl(:,:,i),vxyzcoef(:,:,i)) ! see FROMCOMCOT.f90
 end do

!#[4]## set h_comcot
 !# general info 2021.12.06
 h_comcot%iflag_comcot_spherical = g_param%iflag_comcot_spherical ! 2021.12.08
 h_comcot%IPL            = IPL     ! # of nodes at the sea surface
 h_comcot%nodes          = nodes
 allocate(h_comcot%vxyz1(3,nodes))
 allocate(h_comcot%vxyz2(3,nodes))
 allocate(h_comcot%z_tri1(IPL)   ) ! 2021.12.08
 allocate(h_comcot%z_tri2(IPL)   ) ! 2021.12.08
 h_comcot%nlayer         = nlayer          ! 2021.12.02
 !# layer info
 allocate(h_comcot%layer(nlayer))          ! 2021.12.02
 do i=1,nlayer
  h_comcot%layer(i)%nx             = xygrid(i)%nx
  h_comcot%layer(i)%ny             = xygrid(i)%ny
  CALL ALLOCATEGRDDATA(h_comcot%layer(i)%h_grd,xygrid(i),1) ! 2021.12.02
  h_comcot%layer(i)%h_grd          = h_grd(i) ! 2021.12.02
  allocate(h_comcot%layer(i)%ivxyzflagkl(IPL,3))
  allocate(h_comcot%layer(i)%vxyzcoef(   IPL,2))
  h_comcot%layer(i)%ivxyzflagkl(:,:) = ivxyzflagkl(:,:,i)
  h_comcot%layer(i)%vxyzcoef(:,:)    = vxyzcoef(:,:,i)
 end do

 write(*,*) "### GENCOMCOTCOEF END!! ###"
return
end
!########################################## geomagigrf
!# coded 2018.11.14
 subroutine geomagigrf(h_ocean,em_mesh,g_param,g_meshpara)
 use IGRF
 implicit none
 type(mesh),         intent(in)    :: em_mesh
 type(param_forward),intent(in)    :: g_param
 type(meshpara),     intent(in)    :: g_meshpara
 type(ocean_data),   intent(inout) :: h_ocean
 type(grid_xy)      :: xygrid
 type(grid_data_2D) :: h_grd,Fxyz_grd
 character(70)      :: xgrdfile, ygrdfile,magfile
 integer(4)         :: ivxyzflagkl(h_ocean%IPL,3)
 real(8)            :: vxyzcoef(h_ocean%IPL,2)

!#[1]## set input
 xgrdfile = g_param%xgrdfile_f
 ygrdfile = g_param%ygrdfile_f
 magfile  = g_param%geomagfile

!#[2]## read grid
 CALL COUNTCOMCOT(xygrid, xgrdfile, ygrdfile)   ! get nx, ny  <
 CALL READCOMCOTGRD(xygrid,xgrdfile,ygrdfile,g_param,g_meshpara) !2021.12.08

!#[3]## allocate grid data
 CALL ALLOCATEGRDDATA(Fxyz_grd, xygrid, 3)
 CALL ALLOCATEGRDDATA(   h_grd, xygrid, 3)
 CALL CALGRDIGRF(h_grd,g_param) ! see m_IGRF.f90
 Fxyz_grd%dat(1,:,:) =   h_grd%dat(2,:,:) ! eastward  [nT]
 Fxyz_grd%dat(2,:,:) =   h_grd%dat(1,:,:) ! northward [nT]
 Fxyz_grd%dat(3,:,:) = - h_grd%dat(3,:,:) ! upward    [nT]

!#[4]## calculate coefficient of ivxyzflagkl and vxyzcoef
 CALL MKVXYZCOEF(xygrid,h_ocean,em_mesh,ivxyzflagkl,vxyzcoef) ! see FROMCOMCOT.f90

!#[5]## FXYZ is created
 CALL CALOCEANFXYZ(Fxyz_grd,h_ocean,ivxyzflagkl,vxyzcoef) ! cal h_ocean%Fxyz

return
end
!########################################## geomagfile
subroutine geomagfile(h_ocean,em_mesh,g_param,g_meshpara)
implicit none
type(mesh),         intent(in)    :: em_mesh
type(param_forward),intent(in)    :: g_param
type(meshpara),     intent(in)    :: g_meshpara
type(ocean_data),   intent(inout) :: h_ocean
type(grid_xy)      :: xygrid
type(grid_data_2D) :: h_grd,Fxyz_grd
character(70)      :: xgrdfile, ygrdfile,magfile
integer(4)         :: ivxyzflagkl(h_ocean%IPL,3)
real(8)            :: vxyzcoef(h_ocean%IPL,2)

!#[1]## set input
 xgrdfile = g_param%xgrdfile_f
 ygrdfile = g_param%ygrdfile_f
 magfile  = g_param%geomagfile

!#[2]## read grid
 CALL COUNTCOMCOT(xygrid, xgrdfile, ygrdfile)   ! get nx, ny  <
CALL READCOMCOTGRD(xygrid, xgrdfile, ygrdfile,g_param,g_meshpara) !2021.12.08

!#[3]## allocate grid data
 CALL ALLOCATEGRDDATA(Fxyz_grd, xygrid, 3)
 CALL ALLOCATEGRDDATA(   h_grd, xygrid, 3)
 CALL READGRDDATA(magfile,h_grd,3)
 Fxyz_grd%dat(1,:,:) =   h_grd%dat(2,:,:) ! eastward  [nT]
 Fxyz_grd%dat(2,:,:) =   h_grd%dat(1,:,:) ! northward [nT]
 Fxyz_grd%dat(3,:,:) = - h_grd%dat(3,:,:) ! upward    [nT]

!#[4]## calculate coefficient of ivxyzflagkl and vxyzcoef
 CALL MKVXYZCOEF(xygrid,h_ocean,em_mesh,ivxyzflagkl,vxyzcoef) ! see below

!#[5]## FXYZ is created
 CALL CALOCEANFXYZ(Fxyz_grd,h_ocean,ivxyzflagkl,vxyzcoef) ! cal h_ocean%Fxyz

return
end
!########################################   MKVXYZCOEF
subroutine MKVXYZCOEF(xygrid,h_ocean,em_mesh,ivxyzflagkl, vxyzcoef)
implicit none
type(ocean_data),intent(in)    :: h_ocean
type(mesh),      intent(in)    :: em_mesh
type(grid_xy),   intent(in)    :: xygrid       ! 2021.12.06
real(8),         intent(out)   :: vxyzcoef(h_ocean%IPL,2)
integer(4),      intent(out)   :: ivxyzflagkl(h_ocean%IPL,3)
!#   vxycoef(nodes,1:2) =[ A , B ]
!#   ivxyzflagkl(nodes,1) =0 : v=0, 1 : v can be provided
!#   ivxyzflagkl(nodes,2) =kset  ;  xgrd(kset) < x < xgrd(kset+1)
!#   ivxyzflagkl(nodes,3) =lset  ;  ygrd(lset) <  y  < ygrd(lset+1)
integer(4)                             :: IPL
integer(4),allocatable,dimension(:)    :: surfptr
real(8),   allocatable,dimension(:,:)  :: xyzg  ! coordinates of em3d.msh
real(8),   allocatable,dimension(:)    :: xx, yy
real(8)                                :: x1,y1, A, B
integer(4)                             :: i,j, npoint, k, l, kset, lset, nx,ny

!#[0]## set input
IPL       = h_ocean%IPL    ! nuber of ocean surface nodes
allocate(surfptr(IPL),xyzg(3,em_mesh%node))
surfptr   = h_ocean%surfptr
xyzg      = em_mesh%xyz
nx        = xygrid%nx
ny        = xygrid%ny
allocate( xx(nx),yy(ny)) ! 2021.12.08
xx        = xygrid%xgrd    ! 2021.12.08
yy        = xygrid%ygrd    ! 2021.12.08
!do i=1,nx
!write(*,*) "nx,i,xx,yy",i,xx(i),yy(i)
!end do

!#[1]# surface node loop start
ivxyzflagkl(1:IPL,1:3)=0
vxyzcoef(1:IPL,1:2)=0.d0
do i=1,IPL
!  write(*,*) "i=",i,"IPL=",IPL
!  write(*,*) "i=",i,"IPL=",IPL
  npoint=surfptr(i)
  x1=xyzg(1,npoint)
  y1=xyzg(2,npoint)
!  write(*,*) "x1,y1=",x1,y1
!  write(*,*) "x1=",x1,"y1=",y1
  !#[1-1]# find location of (x1, y1)
  do k=1,nx-1
    if ( xx(k) .lt. x1 .and. x1 .le. xx(k+1) ) then
      kset=k
	 A=(xx(k+1) - x1)/ ( xx(k+1) - xx(k) )  ! coeff for values at x(kset)
	 goto 10
    end if
  end do
  goto 200  ! kset not found
  10 continue  ! lset, A, A1, kset, B, B1
!  write(*,*) "kset=",kset, "x(kset),x(kset+1)=",xx(kset),xx(kset+1)
  do l=1,ny-1
    if ( yy(l) .lt. y1 .and. y1 .le. yy(l+1) ) then
      lset=l
	 B=(yy(l+1) - y1)/ (yy(l+1) - yy(l) ) ! coeff for values at y(lset)
	 goto 20
    end if
  end do
  goto 200  ! lset not found
  20 continue
!  write(*,*) "lset=",lset,"y(lset),y(lset+1)=",yy(lset),yy(lset+1)
  ivxyzflagkl(i,1:3)=(/ 1, kset, lset/) ! vxyz can be calculated
  vxyzcoef(i,1:2)=(/A, B/)
!  write(*,*) "found!!", "kset,lset=",kset,lset, "nx,ny=",nx,ny
  200 continue ! ivxyzflagkl is not set about npoint-th node column
end do

write(*,*) "### MKVXYZCOEF END!! ###"
return
end
!###################################################### CALOCEANVXYZ
subroutine CALOCEANVXYZ(vxyh_grd,z_grd,h_ocean,h_comcot,em_mesh,g_meshpara,vxyz,z_tri)!2021.12.08
implicit none
type(ocean_data),         intent(in)  :: h_ocean
type(comcot_data),        intent(in)  :: h_comcot
type(mesh),               intent(in)  :: em_mesh
type(meshpara),           intent(in)  :: g_meshpara
type(grid_data_2D),       intent(in)  :: vxyh_grd(h_comcot%nlayer) ! 2021.12.06
type(grid_data_2D),       intent(in)  :: z_grd(   h_comcot%nlayer) ! 2021.12.06
real(8),                  intent(out) :: z_tri(h_ocean%IPL)   ! 2021.12.08
real(8),                  intent(out) :: vxyz(3,h_ocean%nodes)! vx, vy, vz at nodes in ocean
type(grid_data_2D),allocatable,dimension(:) :: h_grd               ! 2021.12.06
integer(4),allocatable,dimension(:,:) :: ivxyzflagkl
integer(4),allocatable,dimension(:)   :: surfptr, nz
real(8),   allocatable,dimension(:,:) :: vxyzcoef
real(8),   allocatable,dimension(:)   :: xgrd, ygrd ! vx,vy,vz on comcot grid
real(8),   allocatable,dimension(:,:) :: vx, vy, vz, vxh, vyh, h,z ! 2021.12.08 z added
integer(4) :: IPL, nodes,nlayer,ilayer ! 2021.12.08 nlayer is added
real(8)    :: latorigin, lonorigin
real(8),allocatable,dimension(:,:) :: xyz
real(8)    ::  vx1, vy1, vz1, A,B, A1, B1, h0, dx,dy,z1
integer(4) :: i,j, icount, k,l, ii
integer(4) :: nx,ny

!#[0]# set input
IPL         = h_comcot%IPL
allocate(ivxyzflagkl(IPL,3),surfptr(IPL), nz(IPL),vxyzcoef(IPL,2))
allocate(xyz(3,em_mesh%node))
latorigin   = g_meshpara%latorigin
lonorigin   = g_meshpara%lonorigin
surfptr     = h_ocean%surfptr
nz          = h_ocean%nz
nodes       = h_ocean%nodes
xyz         = em_mesh%xyz
nlayer      = h_comcot%nlayer          ! 2021.12.06
vxyz  =0.d0  ! 2021.12.09
z_tri=0.d0 ! 2021.12.09

do ilayer = 1, nlayer ! 2021.12.07

 ivxyzflagkl = h_comcot%layer(ilayer)%ivxyzflagkl
 vxyzcoef    = h_comcot%layer(ilayer)%vxyzcoef
 nx          = h_comcot%layer(ilayer)%nx
 ny          = h_comcot%layer(ilayer)%ny
 allocate(xgrd(nx), ygrd(ny),vx(nx,ny), vy(nx,ny), vz(nx,ny))
 allocate(z(nx,ny)) ! 2021.12.08
 allocate(vxh(nx,ny), vyh(nx,ny), h(nx,ny))
 xgrd        = h_comcot%layer(ilayer)%h_grd%xygrid%xgrd(:)
 ygrd        = h_comcot%layer(ilayer)%h_grd%xygrid%ygrd(:)
 vxh         = vxyh_grd(ilayer)%dat(1,:,:)
 vyh         = vxyh_grd(ilayer)%dat(2,:,:)
 h           = h_comcot%layer(ilayer)%h_grd%dat(1,:,:)    ! [m]
 z           = z_grd(ilayer)%dat(1,:,:) ! 2021.12.08
 vx(:,:)=0.d0
 vy(:,:)=0.d0
 vz(:,:)=0.d0

!write(*,*) "CALOCEANVXYZ check 1"
!#[1]## cal vx, vy, vz (z=0)
!write(*,*) "check1"
dy=(ygrd(2) - ygrd(1))*d2r*earthrad
!open(30,file="tmp.dat")
do j=2, ny
  do i=2,nx
    if ( h(i,j) .le. 0.d0 ) goto 110 ! 2018.11.14 neglect the water velocity on land
    dx=(xgrd(i) - xgrd(i-1))*d2r*earthrad*dcos(latorigin*d2r)  ! vxh, vyh [ m/s * m ]
    ! Note in COMCOT that M(I,J)=P(I,J)=P_{i+1/2,j}, Q(I,J)=N(I,J)=Q_{i,j+1/2} 2021.12.08 see l.450 in moment.f90
    ! free surface displacement z(i,j)=eta_{i,j}, and ocean depth h(i,j) = H_{i,j} ! 2021.12.08
    !write(30,*) vxh(i,j), vxh(i-1,j), z(i,j), h(i,j), i, j, nx, ny
    if (abs(z(i,j)+h(i,j)) .lt. 1.d-4) go to 110  !1.d-4[km]  2024.01.23 TT
    vx(i,j)=  (vxh(i,j)+vxh(i-1,j))/2.d0/((z(i,j)+h(i,j))*1.d-3)    ![m/s *m /km] = [mm/s] ! z added 2021.12.08
    vy(i,j)=  (vyh(i,j)+vyh(i,j-1))/2.d0/((z(i,j)+h(i,j))*1.d-3)    ![m/s *m /km] = [mm/s] ! z added 2021.12.08
    !# vz = - int_-h^z (dvx/dx + dvy/dy) dz
    !#    = - (dvx/dx + dvy/dy)*(z+h)
    !#    = - (dvxh/dx + dvyh/dy) *(1 + z/h) [m/sec]
    !# [m/s * m * 1/km] = [mm/s]
    vz(i,j)= - ((vxh(i,j)-vxh(i-1,j))/dx + (vyh(i,j)-vyh(i,j-1))/dy) ! value at the surface
    110 continue ! 2018.11.14
  end do
end do
!close(30)
!write(*,*) "nx,ny",nx,ny
!call countnonzero2( vx,nx,ny,ii) ; write(*,*) "vx count",ii
!call countnonzero2( vy,nx,ny,ii) ; write(*,*) "vy count",ii
!call countnonzero2( vz,nx,ny,ii) ; write(*,*) "vz count",ii
!write(*,*) "IPL",IPL
!call countnonzeroint2(ivxyzflagkl,IPL,3,ii); write(*,*) "ivxyzflagkl count",ii
!call countnonzero2(    vxyzcoef,  IPL,2,ii); write(*,*) " vxyzcoef count",ii

!#[2]## cal vx,vy,vz for nodes of tetrahedral mesh
!write(*,*) "check2"
!icount=0

!write(*,*) "h_ocean%nodes =", h_ocean%nodes
!write(*,*) "h_comcot%nodes =", h_ocean%nodes
!write(*,*) "surfptr(IPL),nz(IPL)=",surfptr(IPL),nz(IPL)
!open(25,file="iklAB.dat")

do i=1, IPL
!    write(*,*) "i=",i,"IPL=",IPL
!    write(*,*) "surfptr(i)=",surfptr(i),"nz(i)=",nz(i)
  if ( ivxyzflagkl(i,1) .eq. 1 ) then
    k=ivxyzflagkl(i,2) ;  A=vxyzcoef(i,1) ; A1=1.d0 - A   ! A, A1 are dimensionless
    l=ivxyzflagkl(i,3) ;  B=vxyzcoef(i,2) ; B1=1.d0 - B    ! B, B1 are dimensionless
!    write(25,*) "i k l A B",i,k,l,A,B
    vx1=vx(k,l)*A*B + vx(k+1,l)*A1*B + vx(k,l+1)*A*B1 + vx(k+1,l+1)*A1*B1
    vy1=vy(k,l)*A*B + vy(k+1,l)*A1*B + vy(k,l+1)*A*B1 + vy(k+1,l+1)*A1*B1
    vz1=vz(k,l)*A*B + vz(k+1,l)*A1*B + vz(k,l+1)*A*B1 + vz(k+1,l+1)*A1*B1
     z1= z(k,l)*A*B +  z(k+1,l)*A1*B +  z(k,l+1)*A*B1 +  z(k+1,l+1)*A1*B1 ! 2021.12.08
    h0= - xyz(3, surfptr(i)+nz(i)-1 )  ! [km]
    if ( h0 .lt. 0.001d0 ) goto 100 ! if h is less than 1 m
    z_tri(i) = z1 ! 2021.12.08
    do j=1,nz(i)
	ii=surfptr(i)+(j-1)
!	write(*,*)  "ii=",ii,"nodes=",nodes
      vxyz(1,ii)=vx1
      vxyz(2,ii)=vy1
      vxyz(3,ii)=vz1*(1.d0 + xyz(3, ii)/h0 ) ! commented out on 2016.11.19
    end do
  end if
  100 continue
end do
!close(25)
!write(*,*) "nodes",nodes
!call countnonzero2( vxyz,3,nodes,ii) ; write(*,*) "vxyz count",ii
!call countnonzero1( z_tri,IPL,ii) ; write(*,*) "z_tri count",ii

deallocate(xgrd, ygrd,vx, vy, vz, vxh, vyh, h,z) ! 2021.12.08

end do ! layer loop 2021.12.07

!call countnonzero2( vxyz,3,nodes,ii) ; write(*,*) "vxyz count",ii
!call countnonzero1( z_tri,IPL,ii) ; write(*,*) "z_tri count",ii

!#[3]## check
!if (icount .ne. nodes ) then
!write(*,*) "GEGEGE! icount", icount, "is not equal to nodes", nodes
!write(*,*) "### SOMETHING WRONG in CALOCEANVXYZ !! ###"
!end if

write(*,*) "### CALOCEANVXYZ END!! ###"
return
end

!####################################################
! only for extrude mesh
subroutine CALOCEANFXYZ(Fxyz_grd_2D,h_ocean,ivxyzflagkl,vxyzcoef)
implicit none
type(ocean_data),intent(inout) :: h_ocean
type(grid_data_2D), intent(in) :: Fxyz_grd_2D ! x: eastward, y: northward, z: upward
real(8),            intent(in) :: vxyzcoef(h_ocean%IPL,2)
integer(4),         intent(in) :: ivxyzflagkl(h_ocean%IPL,3)
integer(4) :: nodes, node
real(8),   allocatable, dimension(:,:,:) :: Fxyz_grd
real(8),   allocatable, dimension(:,:)   :: Fxyz
integer(4),allocatable, dimension(:)     :: nz,surfptr, ocean2emptr
integer(4) :: i, j, k, l, ii,nx,ny,IPL
real(8)    :: A, A1, B, B1, Fx, Fy, Fz

!#[0]## set input
nx       = Fxyz_grd_2D%xygrid%nx
ny       = Fxyz_grd_2D%xygrid%ny
IPL      = h_ocean%IPL
allocate(nz(IPL), surfptr(IPL))
nz       = h_ocean%nz
surfptr  = h_ocean%surfptr
node     = h_ocean%node
allocate(Fxyz_grd(3,nx,ny))
allocate(Fxyz(3, node))
nodes    = h_ocean%nodes
allocate(ocean2emptr(nodes))
Fxyz_grd = Fxyz_grd_2D%dat
ocean2emptr = h_ocean%ocean2emptr

!#[1]## generate Fxyz(1:3, 1:nodes)
Fxyz     = 0.d0
do i=1,IPL
  if ( ivxyzflagkl(i,1) .eq. 1) then ! if i-th node column is in tsunami simulation area
     k=ivxyzflagkl(i,2)
     l=ivxyzflagkl(i,3)
     A=vxyzcoef(i,1) ; A1=1.d0 - A
     B=vxyzcoef(i,2) ; B1=1.d0 - B
     Fx=Fxyz_grd(1,k,l)*A*B+Fxyz_grd(1,k+1,l)*A1*B+Fxyz_grd(1,k,l+1)*A*B1+Fxyz_grd(1,k+1,l+1)*A1*B1
     Fy=Fxyz_grd(2,k,l)*A*B+Fxyz_grd(2,k+1,l)*A1*B+Fxyz_grd(2,k,l+1)*A*B1+Fxyz_grd(2,k+1,l+1)*A1*B1
     Fz=Fxyz_grd(3,k,l)*A*B+Fxyz_grd(3,k+1,l)*A1*B+Fxyz_grd(3,k,l+1)*A*B1+Fxyz_grd(3,k+1,l+1)*A1*B1
  else
    Fx=0.d0 ; Fy=0.d0; Fz=0.d0
  end if
  do j=1, nz(i)
     ii=surfptr(i)+j-1
     Fxyz(1,ii)=Fx
     Fxyz(2,ii)=Fy
     Fxyz(3,ii)=Fz
  end do
!  write(*,*) "i=",i,"IPL=",IPL,"nz(i)=",nz(i),"surfptr(i)=",surfptr(i)
!  write(*,*) "Fxyz(1:3,ii)=",Fxyz(1:3,ii)
end do

!#[2]## set output
 h_ocean%Fxyz(1:3,:) = 0.d0
 do i=1,nodes
  h_ocean%Fxyz(1:3,i) = Fxyz(1:3,i)
 end do

!if ( icount .ne. nodes ) then
!write(*,*) "GEGEGE! icount", icount, "is not equal to nodes", nodes
!write(*,*) "### SOMETHING WRONG in CALOCEANFXYZ !! ###"
!end if

write(*,*) "### CALOCEANFXYZ END!! ###"
return
end

!########################################## velocityana
subroutine velocityana(h_ocean,em_mesh,g_param,t)
implicit none
type(ocean_data),   intent(inout) :: h_ocean
type(mesh),         intent(in)    :: em_mesh
type(param_forward),intent(in)    :: g_param
real(8),            intent(in)    :: t       ! [sec]
real(8)    :: elm_xyz(3)
integer(4) :: i,inod

  h_ocean%vxyz(1:3,:) = 0.d0
  do i=1,h_ocean%nodes
   elm_xyz = em_mesh%xyz(1:3,i)
   call calvxyz(elm_xyz,h_ocean%vxyz(1:3,i),t,g_param)
  end do

return
end
!########################################## geomaghomo
subroutine geomaghomo(h_ocean,g_param)
implicit none
type(ocean_data),   intent(inout) :: h_ocean
type(param_forward),intent(in)    :: g_param
integer(4) :: i,inod

  h_ocean%Fxyz(1:3,:) = 0.d0
  do i=1,h_ocean%nodes
   h_ocean%Fxyz(1:3,i) = g_param%geomagvector(1:3)
  end do

return
end

!########################################## genoceanptr
subroutine genoceanptr(h_ocean,em_mesh,g_param)
implicit none
type(mesh),         intent(in)   :: em_mesh
type(param_forward),intent(in)   :: g_param
type(ocean_data),   intent(out)  :: h_ocean
integer(4),         parameter    :: maxnele=70
type(mesh)                       :: ocean_mesh
integer(4),allocatable,dimension(:,:) :: node2ele,n4
! nele(i) is # of element including i-th node
integer(4)                          :: ntet, node, npoi ! 2017.11.02
integer(4)                          :: n1, i, j, icount
integer(4),allocatable,dimension(:) :: em2oceanptr, ocean2emptr,nele
integer(4),allocatable,dimension(:) :: nz, surfptr
integer(4),allocatable,dimension(:) :: igroup ! 2017.11.02
integer(4)                          :: nodes, IPL

!#[0]## set input
ntet    = em_mesh%ntet ! 2017.11.02
npoi    = em_mesh%npoi
node    = em_mesh%node
IPL     = em_mesh%npoi ! # of ocean surface nodes in extrude mesh
allocate(node2ele(maxnele,node),nele(node))
allocate(em2oceanptr(node),n4(ntet,4))  ! 2017.11.02
allocate(ocean2emptr(node),nz(npoi), surfptr(npoi))
allocate( igroup(ntet) )                ! 2017.11.02
n4      = em_mesh%n4
surfptr = em_mesh%n1(:,1)
igroup  = em_mesh%n4flag(:,2) ! 1:air, 2:ocean, 3:earth 2017.11.02

goto 90
!#[1]## gen node-element connectivity data
nele(:)=0
do i=1,ntet ! 2017.11.02
 do j=1,4
  n1=n4(i,j)
  if (maxnele .eq. nele(n1)) goto 998
  nele(n1)=nele(n1)+1     ! increase nele(n1) by 1
  node2ele(nele(n1),n1)=i ! set i-th element to element list of n1 node
 end do
end do
write(*,*) "node-element connectivity end!"

!#[2]## calculate em2oceanptr, ocean2emptr, and nodes
em2oceanptr(:)=0
ocean2emptr(:)=0
icount = 0
do i=1,node
 do j=1,nele(i)
  if ( igroup(node2ele(j,i)) .eq. 2 ) then ! if node belongs to ocean element
   icount = icount + 1
   em2oceanptr(i     ) = icount  ! em -> ocean
   ocean2emptr(icount) = i
   goto 100
  end if
 end do
 100 continue
end do
nodes = icount ! if there are ocean depth = 0 coastline, nodes obtained above is different

90 continue
CALL READMESH_TOTAL(ocean_mesh,g_param%oceanfile)
nodes = ocean_mesh%node

!#[3]## set IPL, surfptr, nz
call MAKENZ3(nz,IPL,surfptr,nodes)

!#[4]## set em2oceanptr and ocean2emptr
 call allocateoceandata(h_ocean,nodes,node,IPL)  ! allocate Fxyz,vxyz,z_tri,em2oceanptr, ocean2emptr 2021.12.08 IPL is added
 h_ocean%em2oceanptr(1:node)  = em2oceanptr(1:node)
 h_ocean%ocean2emptr(1:nodes) = ocean2emptr(1:nodes)
 allocate(h_ocean%surfptr(IPL), h_ocean%nz(IPL))
 h_ocean%IPL     = IPL
 h_ocean%nz      = nz
 h_ocean%surfptr = surfptr

return

998 continue
write(*,*) "GEGEGE! # maxnele=",maxnele,"too small for node-element connectivity data"
stop

end
!##############################################   MAKENZ3
subroutine MAKENZ3(nz,IPL,surfptr,nodes)
implicit none
integer(4) :: IPL,i,nodes
integer(4),dimension(IPL),intent(in) ::surfptr
integer(4),dimension(IPL),intent(out) ::nz
do i=1,IPL-1
nz(i)=surfptr(i+1)-surfptr(i)
end do
nz(IPL)=nodes-surfptr(IPL)+1
write(*,*) "IPL=",IPL
write(*,*) "nz(IPL)=",nz(IPL)
write(*,*) "### MAKENZ3 END ###"
return
END

!########################################## allocateoceandata
subroutine allocateoceandata(h_ocean,nodes,node,IPL) ! 2021.12.08 IPL is added
implicit none
integer(4),intent(in) :: nodes,node,IPL
type(ocean_data),intent(out) :: h_ocean

h_ocean%node  = node
h_ocean%nodes = nodes
allocate(h_ocean%ocean2emptr(nodes))
allocate(h_ocean%em2oceanptr(node))
allocate(h_ocean%Fxyz(3,node))
allocate(h_ocean%vxyz(3,node))
allocate(h_ocean%z_tri(IPL)) ! 2021.12.08

return
end

!##################################################
subroutine OUTVXYZ(h_ocean, em_mesh, it,iflag,g_param) ! 2018.11.14
implicit none
type(ocean_data),        intent(in)   :: h_ocean
type(mesh),              intent(in)   :: em_mesh
integer(4),              intent(in)   :: it
integer(4),              intent(in)   :: iflag   ! 0 for coordinate file, 1 for nothing
type(param_forward),     intent(in)   :: g_param ! 2018.11.14
integer(4)                            :: node, IPL
integer(4),allocatable,dimension(:)   :: surfptr
real(8),   allocatable,dimension(:,:) :: vxyz, xyz
real(8),   allocatable,dimension(:)   :: z_tri    ! 2021.12.08
real(8)       :: vh,ph
character(6)  :: num     ! 2019.01.23 4 -> 6
character(70) :: vhfile, vzfile, zfile ! 2021.12.08 zfile added
integer(4)    :: i, ii,nodes ! 2021.12.08 nodes is added
character(70) :: outfldr  ! 2018.11.14
integer(4)    :: lenout   ! 2018.11.14

!#[1]## set input
node    = h_ocean%node
nodes   = h_ocean%nodes ! 2021.12.08 correction from node to nodes
IPL     = h_ocean%IPL
allocate(surfptr(IPL),vxyz(3,nodes), xyz(3,node)) ! 2021.12.08 vxyz(node) to vxyz(nodes)
allocate(z_tri(IPL))     ! 2021.12.08
surfptr = h_ocean%surfptr
vxyz    = h_ocean%vxyz
 xyz    = em_mesh%xyz
z_tri   = h_ocean%z_tri  ! 2021.12.08
outfldr = g_param%outvxyzfolder ! 2018.11.14
lenout  = len_trim(outfldr)     ! 2018.11.14

!#[1]## coordinate
if (iflag .eq. 0) then
 open(1,file=outfldr(1:lenout)//"coord_xy2D.dat") ! 2018.11.14
 write(1,*) IPL
 do i=1,IPL
  write(1,'(3g15.7)') (xyz(ii,surfptr(i)),ii=1,3)
 end do
end if

!#[2]## output
write(num,'(i6.6)') it ! 2019.01.23
vhfile=outfldr(1:lenout)//"vh"//num(1:6)//".dat" ! 2019.01.19
vzfile=outfldr(1:lenout)//"vz"//num(1:6)//".dat" ! 2019.01.19
 zfile=outfldr(1:lenout)//"z"//num(1:6)//".dat"  ! 2021.12.08
!open(1,file=vhfile)
!open(2,file=vzfile)
open(3,file= zfile) ! 2021.12.08
do i=1,IPL
!  ii=surfptr(i)
!  vh=dsqrt(vxyz(1,ii)**2.d0 + vxyz(2,ii)**2.d0)
!  ph =  datan2(vxyz(2,ii), vxyz(1,ii))*r2d
!  write(1,*) vh
!  write(2,*) vxyz(3,ii)
  write(3,*) z_tri(i)        ! 2021.12.08
end do
!close(1)
!close(2)
close(3)! 2021.12.08
return
end

!####################################################### calvxyzout
! Coded on June 23, 2016, only for this folder, ana_comp/emtet4_edge/
subroutine calvxyzout(xyzminmax,it,dt,g_param)
use param
implicit none
type(param_forward),intent(in) :: g_param
integer(4),intent(in) :: it
real(8),intent(in)    :: dt
real(8),intent(in)   :: xyzminmax(6)
! (1,(1,2,3)) for edge basis fun (x,y,z), (2,(1,2,3)) for face basis fun
integer(4),parameter :: nx=300,ny=300
character(70) :: vxyzfile
character(6) :: num
integer(4) :: i,j
real(8)    :: x3(3),v3(3),xmin,xmax,ymin,ymax,dy,dz,t
xmin=xyzminmax(1)
xmax=xyzminmax(2)
ymin=xyzminmax(3)
ymax=xyzminmax(4)

x3(3)=-3.d0 ! [km]
write(num,'(i6.6)') it
vxyzfile="./vxyz/vx"//num(1:6)//".dat"

open(1,file=vxyzfile)

do i=1,nx
 x3(1)=(xmax-xmin)/nx * (i-1) + xmin
do j=1,ny
 x3(2)=(ymax-ymin)/ny * (j-1) + ymin
 t = dt*it
 call calvxyz(x3,v3,t,g_param) ! see forward_rhs.f90
 write(1,'(3g15.7)') x3(1),x3(2),v3(1) ! v3 [mm/s] only x component
end do
end do

close(1)
write(*,*) "### calvxyzout end!! ###"
return
end


!###########################################
subroutine calvxyz(x3,v3,t,g_param)
implicit none
type(param_forward),intent(in)  :: g_param
real(8),            intent(in)  :: x3(3) ! coordinate [km]
real(8),            intent(in)  :: t    ! [sec]
real(8),            intent(out) :: v3(3)
real(8) :: chi,tlarge,depth,g,c,sigma,x0,w,w1,eta
real(8) :: x1,x2,y1,y2,k,omega
real(8) :: lambda, xwavelength,yextent
!#[1]## set input
tlarge = g_param%h_anapara%tlarge ! [sec]
depth  = g_param%h_anapara%depth  ! ocean depth[km]
lambda = g_param%h_anapara%lambda ! wavelength [km]
xwavelength = g_param%h_anapara%xwavelength
yextent     = g_param%h_anapara%yextent
eta = g_param%h_anapara%eta

!#[2]## set parameters
g=9.8d0
k=2.*pi/lambda ! [1/km]
c=sqrt(g*depth*1.d3) ! [m/s]
omega=k*c/1.d3 ! [rad/s]
y2=yextent
x2=xwavelength/2.*lambda
x1=x2 - 20.d0
y1=y2 - 20.d0

if ( - depth - 0.01 .le. x3(3) .and. x3(3) .le. 0.d0 ) then
 !t=it*dt  ! this is time for velocity when bz_n at (it-1)*dt and bz_n+1 at it*dt
else
 write(*,*) "GEGEGE! x3(1:3)=",x3(1:3),"-depth=",-depth
 stop
end if

!x0=- 200.d0 + 2.5/2.*lambda + c*t/1.d3   ! [km] initial x coordinate of wave peak
x0=0.d0

!# time dependence
if ( t .lt. tlarge) then
 chi=t/tlarge
else
 chi=1.d0
end if

!# eta
!eta=1.d0*exp(-(x3(1) - c*t/1.d3 -x0)**2./2./sigma**2.) ! [m]
eta=eta*cos(k*(x3(1) -c*t/1.d3 - x0)  ) ! [m]

!# dependence on y
if ( y2 .lt. x3(2) ) w=0.d0
if ( y1 .le. x3(2) .and. x3(2) .le. y2 ) w=-1./(y2-y1)*(x3(2)-y2)
if (-y1 .lt. x3(2) .and. x3(2) .lt. y1 ) w=1.d0
if (-y2 .le. x3(2) .and. x3(2) .le.-y1 ) w= 1./(y2-y1)*(x3(2)+y2)
if ( x3(2)  .lt. -y2 ) w=0.d0

!# dependence on x
if ( x2+x0 .lt. x3(1) ) w1=0.d0
if ( x1+x0 .le. x3(1) .and. x3(1) .le. x2+x0 ) w1=-1./(x2-x1)*(x3(1)-(x2+x0))
if (-x1+x0 .lt. x3(1) .and. x3(1) .lt. x1+x0 ) w1=1.d0
if (-x2+x0 .le. x3(1) .and. x3(1) .le.-x1+x0 ) w1= 1./(x2-x1)*(x3(1)-(-x2+x0))
if ( x3(1)  .lt. -x2+x0 ) w1=0.d0

!# calculate v3
v3(1:3)=0.d0
v3(1)=chi*w*w1*(c*eta/depth) ! [mm/s]

return
end

end module oceanFvxyz
