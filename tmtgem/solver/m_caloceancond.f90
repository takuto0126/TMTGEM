! coded on 2018.11.13
module caloceancond
use global_ocean_cond_map ! m_global_ocean_cond.f90
use mesh_type             ! m_mesh_type.f90
use param                 ! m_param.f90
use spherical             ! m_spherical.f90
use oceanFvxyz            ! m_oceanFvxyz.f90 2019.01.28
implicit none

contains
!################################################################### setoceancond
!# coded on 2018.11.13
subroutine setoceancond(g_mesh,g_param,h_ocean) ! 2019.01.28
implicit none
type(mesh),            intent(inout)  :: g_mesh  ! m_mesh_type.f90
type(param_forward),   intent(in)     :: g_param ! m_param.f90
type(ocean_data),      intent(in)     :: h_ocean ! m_oceanFvxyz.f90 2019.01.28
type(ocean_cond)                      :: g_cond ! 2019.01.27
integer(4)                            :: ntet
integer(4)                            :: node
real(8),   allocatable,dimension(:,:) :: xyz
integer(4),allocatable,dimension(:,:) :: n4
integer(4),allocatable,dimension(:)   :: igroup
real(8),   allocatable,dimension(:)   :: lon,lat,depth,cond
real(8),   allocatable,dimension(:)   :: d_top,d_bot   ! 2019.01.27
real(8),   allocatable,dimension(:)   :: cmodel
character(50)                         :: woafile
integer(4)                            :: ntet_ocean, nshift, i,j,k,ii
real(8)                               :: lonorigin,latorigin,lonlatalt(3)
real(8),   allocatable,dimension(:,:) :: xc            ! 2018.11.15
integer(4)                            :: iflag_oceancond
integer(4),allocatable,dimension(:)   :: ocean2ntetptr ! 2018.11.14
real(8)                               :: oceancond_default,a,b ! 2018.11.15
character(50)                         :: fna      ! 2019.01.28
type(ocean_cond) :: m_cond

!# set
 ntet              = g_mesh%ntet
 node              = g_mesh%node
 allocate( igroup(ntet),n4(ntet,4),xyz(3,node) )
 xyz               = g_mesh%xyz
 n4                = g_mesh%n4
 igroup            = g_mesh%n4flag(:,2)
 woafile           = g_param%woafile
 lonorigin         = g_param%lonorigin
 latorigin         = g_param%latorigin
 iflag_oceancond   = g_param%iflag_oceancond ! 1 for constant, 2 for woafile
 oceancond_default = g_param%cond(2)         ! 2018.11.15

!#[0]## gen default cmodel
 allocate(cmodel(ntet))
 do i=1,ntet
  if ( igroup(i) .eq. 1 ) cmodel(i) = g_param%cond(1) ! air conductivity
  if ( igroup(i) .eq. 2 ) cmodel(i) = g_param%cond(2) ! default ocean cond [S/m]
  if ( igroup(i) .eq. 3 ) cmodel(i) = g_param%cond(3) ! solid earth conductivity
 end do

 !#[1]## set ocean2ntetptr
  ntet_ocean = 0
  do i=1,ntet
   if ( igroup(i) .eq. 2 ) ntet_ocean = ntet_ocean + 1
  end do
  write(*,*) "ntet_ocean=",ntet_ocean ! 2019.01.25
  allocate(ocean2ntetptr(ntet_ocean))
  j=0
  do i=1,ntet
   if ( igroup(i) .eq. 2 ) then
    j=j+1
    ocean2ntetptr(j) = i
   end if
  end do

!#========== iflag_oceancond = 2 ( woa file is given) ========================  woafile !
if (iflag_oceancond .eq. 2 ) then ! when woa file is given

 !#[2]## get global ocean cod
  call read_global_cond(g_cond,woafile) ! see m_globalocean_cond.f90
  !# replace fillval by oceancond_default 2019.01.28
  do i=1,g_cond%nlon ; do j=1,g_cond%nlat ; do k=1,g_cond%ndepth
    if ( g_cond%val(i,j,k) .lt. 0.d0 ) g_cond%val(i,j,k) = oceancond_default
  end do ; end do; end do

 !#[3]## prepare lon, lat, depth
  allocate( lon(ntet_ocean),lat(ntet_ocean),depth(ntet_ocean),xc(3,ntet_ocean))
  allocate( d_top(ntet_ocean), d_bot(ntet_ocean)) ! 2019.01.27
  xc(:,:)  = 0.d0 ;  d_top(:) = 0.d0 ; d_bot = 0.d0 ! [m] downward positive 2019.01.27
  do i=1,ntet_ocean
   ii = ocean2ntetptr(i) ! 2019.01.27
   do k=1,4
    xc(1:3,i)=xc(1:3,i) + xyz(1:3,n4(ii,k))/4.d0
   end do
   call xyz2lonlatalt(xc(:,i),lonorigin,latorigin,lonlatalt)
   lon(i)   =   lonlatalt(1)
   lat(i)   =   lonlatalt(2)
   depth(i) = - lonlatalt(3)*1000. ! upward positive [km] -> downward positive [m]
   d_top(i) = -1000.* max(xyz(3,n4(ii,1)),xyz(3,n4(ii,2)),xyz(3,n4(ii,3)),xyz(3,n4(ii,4))) ! 2019.01.27
   d_bot(i) = -1000.* min(xyz(3,n4(ii,1)),xyz(3,n4(ii,2)),xyz(3,n4(ii,3)),xyz(3,n4(ii,4))) ! 2019.01.27
  end do

 !#[4]## cal ocean cond
  allocate( cond(ntet_ocean))

  !#[4-1]## output cmodel without complementation of woacond
   if ( .true. ) then
    call calvolumetricmeanofoceancond(ntet_ocean,lon,lat,d_top,d_bot,cond,g_cond) ! 2019.01.27
    do i=1,ntet_ocean
     cmodel(ocean2ntetptr(i)) = cond(i)
    end do
    fna="surfcond_wocomp.out"
    call outsurfcond(g_param,g_mesh,cmodel,ntet,ntet_ocean,ocean2ntetptr,fna)!2019.01.28
    fna="bottcond_wocomp.out"
    call outbottcond(g_param,g_mesh,cmodel,ntet,ntet_ocean,ocean2ntetptr,fna,h_ocean)!2019.01.28
   end if
   m_cond = g_cond

  !#[4-2]## cal cmodel with g_cond complementation
  ! call calpointoceancond(ntet_ocean,lon,lat,depth,cond,g_cond) ! commented out on 2019.01.27
   call complementemptybinofwoacond(g_param,g_cond) ! 2019.01.27 see below
   call calvolumetricmeanofoceancond(ntet_ocean,lon,lat,d_top,d_bot,cond,g_cond) ! 2019.01.27
   do i=1,ntet_ocean
    cmodel(ocean2ntetptr(i)) = cond(i)
   end do

 end if
!#========== iflag_oceancond = 2 ( woa file is given) ======================== end woafile case!
   if (.false. ) then ! 2019.01.28
    a=142.8151 ;   b=40.566 ! shoice of lon lat where you want to check vertical profile
    call search(m_cond%nlon,m_cond%nlat,m_cond%lonbnds,m_cond%latbnds,a,b,i,j)
    do k=1,g_cond%ndepth
     write(*,*) k,m_cond%val(i,j,k),"->",g_cond%val(i,j,k)
    end do
   end if

!# output final surfcond.out
   fna="surfcond.out"
   call outsurfcond(g_param,g_mesh,cmodel,ntet,ntet_ocean,ocean2ntetptr,fna)!2019.01.28
   fna="bottcond.out"
   call outbottcond(g_param,g_mesh,cmodel,ntet,ntet_ocean,ocean2ntetptr,fna,h_ocean)!2019.01.28

!# generate output
  ! cmodel is allocated by "allocate (h_mesh%cmodel(h_mesh%ntet))" in mesh/m_mesh_type.f90
  g_mesh%cmodel = cmodel

return
end
!######################################################  outsurfcond
 subroutine outsurfcond(g_param,g_mesh,cmodel,ntet,ntet_ocean,ocean2ntetptr,fna)
 implicit none
 type(param_forward),intent(in) :: g_param
 type(mesh),         intent(in) :: g_mesh
 integer(4),         intent(in) :: ntet, ntet_ocean
 integer(4),         intent(in) :: ocean2ntetptr(ntet_ocean)
 real(8),            intent(in) :: cmodel(ntet)
 character(50),      intent(in) :: fna
 real(8)                        :: xyz_ele(3,4)
 integer(4)                     :: i,j,icount
 real(8)                        :: xc(2)

 !# output cmodel ! 2019.01.25
  open(1,file=trim(g_param%outbxyzfolder)//trim(fna)) ! 2019.01.28
   do i=1,ntet_ocean
    icount = 0 ! 2019.01.27
    do j=1,4
     xyz_ele(1:3,j)=g_mesh%xyz(1:3,g_mesh%n4(ocean2ntetptr(i),j))
     if ( xyz_ele(3,j) .gt. -0.001 ) icount = icount + 1 ! 2019.01.27
    end do
     if ( icount .eq. 3 ) then ! one triangle is on the ocean 2019.01.27
     xc(1) = sum(xyz_ele(1,:))/4.d0
     xc(2) = sum(xyz_ele(2,:))/4.d0
    write(1,'(3f12.5)') xc(1:2),cmodel(ocean2ntetptr(i))
    end if
   end do
  close(1)

 return
 end

!######################################################  outbottcond
!# coded on 2019.01.28
 subroutine outbottcond(g_param,g_mesh,cmodel,ntet,ntet_ocean,ocean2ntetptr,fna,h_ocean)
 implicit none
 type(param_forward),intent(in) :: g_param
 type(ocean_data),   intent(in) :: h_ocean ! m_oceanFvxyz.f90
 type(mesh),         intent(in) :: g_mesh
 integer(4),         intent(in) :: ntet, ntet_ocean
 integer(4),         intent(in) :: ocean2ntetptr(ntet_ocean)
 real(8),            intent(in) :: cmodel(ntet)
 character(50),      intent(in) :: fna
 real(8)                        :: xyz_ele(3,4)
 integer(4)                     :: i,j,icount,node,IPL,ii
 integer(4),allocatable,dimension(:) :: surfptr,nz
 logical,   allocatable,dimension(:) :: iflag_bot
 real(8)                        :: xc(2)

 !#[1]## set
 IPL     = h_ocean%IPL
 allocate(surfptr(IPL),nz(IPL))
 surfptr = h_ocean%surfptr
 nz      = h_ocean%nz
 node    = g_mesh%node
 allocate( iflag_bot(node) )
 iflag_bot = .false.

 !#[2]## bottom node label
  do i=1,IPL
   iflag_bot(surfptr(i) + nz(i) - 1) = .true.
  end do

 !#[3]##
  open(1,file=trim(g_param%outbxyzfolder)//trim(fna)) ! 2019.01.28
   do i=1,ntet_ocean
    icount = 0 ! 2019.01.27
    do j=1,4
     ii = g_mesh%n4(ocean2ntetptr(i),j)
     xyz_ele(1:3,j)=g_mesh%xyz(1:3,ii)
     if ( iflag_bot(ii) ) icount = icount + 1 ! 2019.01.27
    end do
     if ( icount .eq. 3 ) then ! one triangle is on the ocean 2019.01.27
     xc(1) = sum(xyz_ele(1,:))/4.d0
     xc(2) = sum(xyz_ele(2,:))/4.d0
     write(1,'(3f12.5)') xc(1:2),cmodel(ocean2ntetptr(i))
    end if
   end do
  close(1)

 return
 end

!##################################################### complementemptybinofwoacond
!# coded on 2019.01.27
  subroutine complementemptybinofwoacond(g_param,g_cond) ! 2019.01.27
  implicit none
  type(param_forward),intent(in)    :: g_param
  type(ocean_cond),   intent(inout) ::  g_cond
  integer(4) :: nlon,nlat,ndepth,ilon,jlat,idepth
  real(8),  allocatable,dimension(:,:) :: lonbnds,depthbnds,latbnds
  real(8)    :: xyzminmax(6),lonb(2),latb(2),altb(2)
  real(8)    :: xc1(3),xc2(3)
  real(8)    :: lonorigin,latorigin,lonlatalt1(3),lonlatalt2(3)
  type(ocean_cond) :: h_cond
  integer(4) :: icount,i,j,k
  real(8)    :: cond

  !#[1]## set
   nlon      = g_cond%nlon
   nlat      = g_cond%nlat
   ndepth    = g_cond%ndepth
   allocate(lonbnds(2,nlon),latbnds(2,nlat),depthbnds(2,ndepth))
   lonbnds   = g_cond%lonbnds
   latbnds   = g_cond%latbnds
   depthbnds = g_cond%depthbnds
   xyzminmax = g_param%xyzminmax
   lonorigin         = g_param%lonorigin
   latorigin         = g_param%latorigin

  !#[2]## lonlataltminmax for em_mesh
   xc1(1:3)=(/xyzminmax(1),xyzminmax(3),xyzminmax(5)/)
   xc2(1:3)=(/xyzminmax(2),xyzminmax(4),xyzminmax(6)/)
   call xyz2lonlatalt(xc1,lonorigin,latorigin,lonlatalt1) ! min coordinate
   call xyz2lonlatalt(xc2,lonorigin,latorigin,lonlatalt2) ! max coordinate
   lonb=(/lonlatalt1(1),lonlatalt2(1)/)
   latb=(/lonlatalt1(2),lonlatalt2(2)/)

  !#[4]## complement with horizontally surrounding bins
   h_cond    = g_cond
   do i=1,nlon
    if ( lonb(1) .le. lonbnds(2,i) .and.  lonbnds(1,i) .le. lonb(2) ) then
    do j=1,nlat
     if ( latb(1) .le. latbnds(2,j) .and.  latbnds(1,j) .le. latb(2) ) then
     do k=1,ndepth
      if ( .not. g_cond%exist(i,j,k) ) then ! when not exist
       icount = 0 ; cond = 0
       ilon = i - 1 ; jlat = j      ! left　bin
       call addcond(cond,icount,ilon,jlat,k,g_cond)
       ilon = i + 1 ; jlat = j      ! right bin
       call addcond(cond,icount,ilon,jlat,k,g_cond)
       ilon = i     ; jlat = j - 1  ! down bin
       call addcond(cond,icount,ilon,jlat,k,g_cond)
       ilon = i     ; jlat = j + 1  ! up
       call addcond(cond,icount,ilon,jlat,k,g_cond)
       if ( icount .eq. 0 ) then
        ilon = i - 1 ; jlat = j - 1 ! left　down
        call addcond(cond,icount,ilon,jlat,k,g_cond)
        ilon = i - 1 ; jlat = j + 1 ! left up
        call addcond(cond,icount,ilon,jlat,k,g_cond)
        ilon = i + 1 ; jlat = j - 1 ! right down
        call addcond(cond,icount,ilon,jlat,k,g_cond)
        ilon = i + 1 ; jlat = j + 1 ! right up
        call addcond(cond,icount,ilon,jlat,k,g_cond)
       end if
       if ( icount .ne. 0 ) then ! when a bin with value exist surrounding the column
	  cond = cond / (1.*icount)
	  h_cond%exist(i,j,k) = .true.
	  h_cond%val(i,j,k)   = cond
	 end if
      end if
     end do
     end if
    end do
    end if
   end do
   g_cond = h_cond

  !#[3]## complement downward when the surface cond exists 2019.01.28
   do i=1,nlon
    if ( lonb(1) .le. lonbnds(2,i) .and.  lonbnds(1,i) .le. lonb(2) ) then
    do j=1,nlat
     if ( latb(1) .le. latbnds(2,j) .and.  latbnds(1,j) .le. latb(2) ) then
     if ( g_cond%exist(i,j,1) ) then ! when the top cond exist
      !# depth loop
      do k=2,ndepth
       if ( .not. g_cond%exist(i,j,k) ) then ! when not exist
	  g_cond%exist(i,j,k:ndepth) = .true.  ! generate
        g_cond%val(  i,j,k:ndepth) = g_cond%val(i,j,k-1) ! set homogeneous cond below kth bin
	  goto 110
       end if
      end do
	!# depth loop end !!
      110 continue
     end if
     end if
    end do ! lat loop end
    end if
   end do

  !#[5]## output
  !g_cond%exist = h_cond%exist
  !g_cond%val   = h_cond%val

  write(*,*) "### complementemptybinofwoacond END!! ###"
  return
  end
!#####
  subroutine addcond(cond,icount,ilon,jlat,k,g_cond)
  implicit none
  integer(4),      intent(in)    :: ilon,jlat,k
  type(ocean_cond),intent(in)    :: g_cond
  real(8),         intent(inout) :: cond
  integer(4),      intent(inout) :: icount

       if ( g_cond%exist(ilon,jlat,k) ) then
        icount = icount + 1
        cond   = cond + g_cond%val(ilon,jlat,k)
       end if

  return
  end

!####################################################################  caloceancond
!# coded on Aug 31, 2018
 subroutine calpointoceancond(npoint,lon,lat,depth,cond,g_cond) ! 2019.01.27
 implicit none
 integer(4),                 intent(in)   :: npoint
 real(8),                    intent(in)   :: lon(npoint),lat(npoint),depth(npoint)
 real(8),                    intent(out)  :: cond(npoint)          ! 2019.01.25
 type(ocean_cond),           intent(in)   :: g_cond ! 2019.01.25
 integer(4)                               :: i,j,k,ii,ilon,jlat,kdepth
 integer(4)                               :: nlon,nlat,ndepth
 real(8),   allocatable,dimension(:,:)    :: lonbnds  ![deg] lon bin
 real(8),   allocatable,dimension(:,:)    :: latbnds  ![deg] lat bin
 real(8),   allocatable,dimension(:,:)    :: depthbnds! [m] depth bin(down positive)
 real(8),   allocatable,dimension(:,:,:)  :: val      ! [S/m] global cond grid data
 real(8)                                  :: dd       ! 2018.11.14

!#[2]# set ##
 nlon      = g_cond%nlon
 nlat      = g_cond%nlat
 ndepth    = g_cond%ndepth
 write(*,*) "nlon",nlon,"nlat",nlat,"ndepth",ndepth             ! 2019.01.25
 allocate( lonbnds(2,nlon),latbnds(2,nlat),depthbnds(2,ndepth)) ! 2019.01.25
 lonbnds   = g_cond%lonbnds
 latbnds   = g_cond%latbnds
 depthbnds = g_cond%depthbnds

!#[3]## cal cond(1:npoint)
  do ii=1,npoint

   !# depth adjust 2018.11.14
   dd = depth(ii)
   if ( depth(ii) .gt. depthbnds(2,ndepth) ) then ! depth is downward positive
    dd = depthbnds(2,ndepth)
   end if

   do i=1,nlon
    if ( lonbnds(1,i) .le. lon(ii) .and. lon(ii) .le. lonbnds(2,i)) then
    do j=1,nlat
     if ( latbnds(1,j) .le. lat(ii) .and. lat(ii) .le. latbnds(2,j)) then
     do k=1,ndepth
     if ( depthbnds(1,k) .le. dd .and. dd .le. depthbnds(2,k)) then
      ilon   = i
      jlat   = j
      kdepth = k
      goto 100
     end if
     end do
     end if
    end do
    end if
   end do
   write(*,*) "GEGEGE not found lon lat depth=",lon(ii),lat(ii),dd
   stop
100 continue
    cond(ii) = g_cond%val(ilon,jlat,kdepth)
  end do

return
end

!######################################################## calvolumetricmeanofoceancond
!# This subroutine doesn't consider whether the value exist in the bin or not 2019.01.28
!# coded on Jan 27, 2019
 subroutine calvolumetricmeanofoceancond(npoint,lon,lat,d_top,d_bot,cond,g_cond) ! 2019.01.27
 implicit none
 integer(4),                 intent(in)   :: npoint
 real(8),                    intent(in)   :: lon(npoint),lat(npoint)
 real(8),                    intent(in)   :: d_top(npoint),d_bot(npoint) ! 2019.01.27
 real(8),                    intent(out)  :: cond(npoint)                ! 2019.01.25
 type(ocean_cond),           intent(in)   :: g_cond                      ! 2019.01.27
 integer(4)                               :: i,j,k,k1,k2,ii,ilon,jlat,kdepth
 integer(4)                               :: nlon,nlat,ndepth
 real(8),   allocatable,dimension(:,:)    :: lonbnds   ![deg] lon bin
 real(8),   allocatable,dimension(:,:)    :: latbnds   ![deg] lat bin
 real(8),   allocatable,dimension(:,:)    :: depthbnds ! [m] depth bin(down positive)
 real(8),   allocatable,dimension(:,:,:)  :: val       ! [S/m] global cond grid data
 real(8)                                  :: d1,d2,dlen,len! 2019.01.27

!#[2]# set ##
 nlon      = g_cond%nlon
 nlat      = g_cond%nlat
 ndepth    = g_cond%ndepth
 write(*,*) "nlon",nlon,"nlat",nlat,"ndepth",ndepth             ! 2019.01.25
 allocate( lonbnds(2,nlon),latbnds(2,nlat),depthbnds(2,ndepth)) ! 2019.01.25
 lonbnds   = g_cond%lonbnds
 latbnds   = g_cond%latbnds
 depthbnds = g_cond%depthbnds

!#[3]## cal cond(1:npoint)
  do ii=1,npoint

   !#[3-1]## depth adjust 2018.11.14
   d1 = d_top(ii) ! 2019.01.27
   d2 = d_bot(ii) ! 2019.01.27 d is downward positive
   if ( d1 .gt. depthbnds(2,ndepth) ) d1 = depthbnds(2,ndepth)
   if ( d2 .gt. depthbnds(2,ndepth) ) d2 = depthbnds(2,ndepth)

   do i=1,nlon
    if ( lonbnds(1,i) .le. lon(ii) .and. lon(ii) .le. lonbnds(2,i)) then
    ilon   = i ! 2019.01.27
    do j=1,nlat
     if ( latbnds(1,j) .le. lat(ii) .and. lat(ii) .le. latbnds(2,j)) then
     jlat   = j ; k1 = 0 ; k2 = 0 ! 2019.01.27
     do k=1,ndepth
      if ( depthbnds(1,k) .le. d1 .and. d1 .le. depthbnds(2,k)) k1 = k ! 2019.01.27
      if ( depthbnds(1,k) .le. d2 .and. d2 .le. depthbnds(2,k)) k2 = k ! 2019.01.27
      if ( k1 .ne. 0 .and. k2 .ne. 0 ) goto 100 ! 2019.01.27
     end do !
     end if
    end do ! nlat loop
    end if
   end do ! nlon loop
   write(*,*) "GEGEGE not found lon lat depth=",lon(ii),lat(ii),d1,d2
   stop
100 continue

  !#[3-2]## calculate volumetric mean
   if ( k1 .eq. k2 ) then ! in the same bin
    cond(ii) = g_cond%val(ilon,jlat,k1)
   else if ( k1 .lt. k2 ) then
    !# top bin
    dlen     = depthbnds(2,k1) - d1
    len      = dlen
    cond(ii) = g_cond%val(ilon,jlat,k1) * dlen ! conductance [S/m * m]
    !# intermediate bins
    do k=k1+1,k2-1
     dlen = depthbnds(2,k) - depthbnds(1,k)
     len  = len + dlen
     cond(ii) = cond(ii) + g_cond%val(ilon,jlat,k) * dlen
    end do
    !# bottom bin
    dlen     = d2 - depthbnds(1,k2)
    len      = len + dlen
    cond(ii) = cond(ii) + g_cond%val(ilon,jlat,k2) * dlen ! conductance [S/m * m]
    !# conductance -> conductivity
    cond(ii) = cond(ii) / len
   else
     write(*,*) "GEGEGE stop! "
   end if

  end do ! npoint loop

return
end

end module


