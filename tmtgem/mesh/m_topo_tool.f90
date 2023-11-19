!# coded on 2018.08.27
module topo_tool
implicit none

contains

!##########################################################   countnode1
!# moved from coastline.f90 on 2019.02.27
subroutine countnode1(gebcofile,node)
implicit none
integer(4) :: node
character(70) :: gebcofile
node=0
open(1,file=gebcofile)
do while (node .ge. 0)
  read(1,*,end=99)
  node=node+1
end do
99 continue
close(1)
write(*,'(a,a,a,i10)') " # of poins of ", trim(gebcofile)," is", node  ! inserted on Oct. 24, 2015
write(*,*) "### countnode1 end! ###"
return
end subroutine countnode1

!##########################################################   readgebco2
subroutine readgebco2(gebcofile,gebco_grd)
use coastline_data
implicit none
character(70),  intent(in)    :: gebcofile
type(grid_data),intent(inout) :: gebco_grd
integer(4)                    :: i
real(8)                       :: alt

!#[1]## read grd data
open(2,file=gebcofile)
 do i=1, gebco_grd%node ! node is set
  read(2,*) gebco_grd%lonlatalt(1:2,i),alt ! alt [m] (downward positive) 2019.03.06
  gebco_grd%lonlatalt(3,i) = - alt ! -> upward positive 2019.03.06
  if ( alt .eq. 0.d0 ) gebco_grd%lonlatalt(3,i) = 1.d0 ! inserted by T.M. on Oct. 24, 2015
 end do
close(2)

!#[2]## calculate neast and nsouth
do i=1,gebco_grd%node
 if (gebco_grd%lonlatalt(2,i) .ne. gebco_grd%lonlatalt(2,i+1)) goto 99    ! latitude
end do
99 continue
gebco_grd%neast=i
gebco_grd%nsouth=gebco_grd%node/gebco_grd%neast ! neast and nsout are set for h_grd

write(*,*) "### readgebco2 end! ###"
return
end subroutine readgebco2

!########################################################### lonlattoxy4
! moved from coastline.f90 on 2019.02.27
subroutine lonlattoxy4(g_grd,g_meshpara)
use constants
use coastline_data
use param_mesh
implicit none
type(meshpara),    intent(in)      :: g_meshpara
type(grid_data),   intent(inout)   :: g_grd
real(8)                            :: lonorigin,latorigin
integer(4)                         :: i,node
real(8),allocatable,dimension(:,:) :: lonlatalt ! 2018.08.28

!#[0]## set
 node      = g_grd%node           ! 2108.08.28
 allocate( lonlatalt(3,node))     ! 2018.08.28
 lonlatalt = g_grd%lonlatalt      ! lonlatalt(3,g_grd%node) 2018.08.28
 lonorigin = g_meshpara%lonorigin
 latorigin = g_meshpara%latorigin

!#[1]## calculate xyz
do i=1,g_grd%node
  g_grd%xyz(1,i) = earthrad*dcos(latorigin*d2r)*(lonlatalt(1,i)-lonorigin)*d2r ! eastward [km]
  g_grd%xyz(2,i) = earthrad*(lonlatalt(2,i)-latorigin)*d2r      ! northward [km]
  g_grd%xyz(3,i) = lonlatalt(3,i)
end do

write(*,*) "### lonlattoxy4 end!###"
return
end
!########################################################## neasatsouth3
 subroutine TRUNCATEGRD3(h_grd,g_grd,g_meshpara)
 use constants
 use coastline_data
 use param_mesh
 ! ne1 is # of node read from west, ns1 is that from north
 implicit none
 type(meshpara), intent(in)  :: g_meshpara
 type(grid_data),intent(in)  :: h_grd ! whole grid
 type(grid_data),intent(out) :: g_grd ! truncated grid
 integer(4) :: node, ne1, ne2, ns1, ns2
 integer(4) :: i,j,ii,jj
 real(8)    :: wlon,elon,slat,nlat,lenout,lonorigin,latorigin
 real(8)    :: wlon0,elon0,slat0,nlat0
 real(8),allocatable,dimension(:,:) :: lonlatalt ! 2018.08.28
 logical    :: iflag_topright_land               ! 2019.02.25
 integer(4) :: neast ! 2019.02.25

!#[0]## set parameter
 allocate( lonlatalt(3,h_grd%node))
 lonlatalt = h_grd%lonlatalt  ! lonlatalt(3,h_grd%node)
 wlon      = g_meshpara%wlon
 elon      = g_meshpara%elon
 slat      = g_meshpara%slat
 nlat      = g_meshpara%nlat
 lenout    = g_meshpara%lenout
 lonorigin = g_meshpara%lonorigin
 latorigin = g_meshpara%latorigin
 neast     = h_grd%neast      ! 2019.02.26
 write(*,*) "neast",neast

!#[1]## calculate wlon0,elon0,slat0,nlat0
 wlon0 = wlon-lenout/(earthrad*d2r*dcos(latorigin*d2r))
 elon0 = elon+lenout/(earthrad*d2r*dcos(latorigin*d2r))
 slat0 = slat-lenout/(earthrad*d2r)
 nlat0 = nlat+lenout/(earthrad*d2r)

!#[2]## check boundary of topo file
 if ( wlon0 .lt. lonlatalt(1,1)) then
  write(*,*) "GEGEGE! western boundary of topofile is too east wlon0",wlon0,"topo west",lonlatalt(1,1)
  stop
 end if
 if ( elon0 .gt. lonlatalt(1,neast) ) then
  write(*,*) "GEGEGE! eastern boundary of topofile is too west elon0",elon0,"topo east",lonlatalt(1,neast)
  stop
 end if
 if ( nlat0 .gt. lonlatalt(2,1)) then
  write(*,*) "GEGEGE! north boundary of topofile is too south nlat0",nlat0,"topo north",lonlatalt(2,1)
  stop
 end if
 if ( slat0 .lt. lonlatalt(2,h_grd%node)) then
  write(*,*) "GEGEGE! south boundary of topofile is too north slat0",slat0,"topo south",lonlatalt(2,h_grd%node)
  stop
 end if

!#[2]## create ne1,ne2
!#      wlon0           elon0
!#        |               |
!#  |   | <---- ne2 ------> |            |
!#  1  ne1              (ne1-1)+ne2    neast
!#      1                   ne2
ne1=0;ne2=0
do i=1,h_grd%neast
 if (lonlatalt(1,i) .lt. wlon0) ne1=i              ! start longitude
 if (lonlatalt(1,i) .lt. elon0) ne2=i+1-(ne1-1)  ! end longitude
end do

!#[3]## calculate ns1 and ns2
ns1=0;ns2=0
do i=1,h_grd%nsouth
 if (lonlatalt(2,(i-1)*h_grd%neast +1) .gt. nlat0 ) ns1=i
 if (lonlatalt(2,(i-1)*h_grd%neast +1) .gt. slat0 ) ns2=i+1-(ns1-1)
! write(*,*) "i",i,"ns1",ns1,"ns2",ns2
end do
!stop

!#[4]## generate new grd
 node=ne2*ns2
 CALL ALLOCATEGRD(g_grd,node)
 g_grd%neast=ne2
 g_grd%nsouth=ns2
 jj=0
 do j=ns1, (ns1-1)+ns2
  do i=ne1,(ne1-1)+ne2
   ii= (j-1)*h_grd%neast + i
   jj=jj+1
   g_grd%lonlatalt(1:3,jj)=h_grd%lonlatalt(1:3,ii)
!   g_grd%xyz(1:3,jj)      =h_grd%xyz(1:3,ii)      ! commented out 2019.02.25
  end do
 end do

!#[5]## iflag_topright_land ### 2019.02.25
 iflag_topright_land=.false.
 if ( g_grd%lonlatalt(3,ne2) .ge. 0.d0 ) iflag_topright_land=.true. ! 2019.02.25
 g_grd%iflag_topright_land = iflag_topright_land              ! 2019.02.25

write(*,*) "ns1=",ns1,"ne1=",ne1
write(*,*) "ns2=",ns2,"ne2=",ne2,"truncated node=",g_grd%node
!write(*,*) "grd%xyz(3,neast)=",g_grd%xyz(3,g_grd%neast)
write(*,'(a,l3)') "iflag_toplight_land", g_grd%iflag_topright_land

write(*,*) "### TRUNCATEGRD3 end!###"

return
end subroutine TRUNCATEGRD3
end module
