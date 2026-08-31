!# adopt "precision.inc" for len_file on 2016.11.15
!# 
!# grid_data_2D is added on May 16, 2016
!#
!# Coded on Oct. 25, 2015
!# to calculate (vx,vy,vz) from COMCOT results
!# Contents :
!#   subroutine  COUNTCOMCOT
!#   subroutine  READCOMCOTGRD
!#   subroutine  READGRD
!#   subroutine  READMN
!#   subroutine CALOCEANVXYZ
!#   subroutine CALOCEANFXYZ
!#  subroutine GENLVEC3COM(Fxyz, nodes, g_node, nodtot, Fxyzip, ip)
!########################################   COUNTCOMCOT
module FROMCOMCOT
use param
use param_mesh ! 2021.12.08
implicit none

type grid_xy
 integer(4) :: nx
 integer(4) :: ny
 real(8),allocatable,dimension(:)   :: xgrd,ygrd ! Cartesian coordinate [km]
 real(8),allocatable,dimension(:)   :: xgrd_sphe,ygrd_sphe ! Spherical coordinate [deg]2021.12.08
end type grid_xy

type grid_data_2D
 type(grid_xy) :: xygrid
 real(8),allocatable,dimension(:,:,:) :: dat
 ! dat(i,j,k) is
 ! i-th component at coordinate (j,k) in nx*ny grid
end type grid_data_2D

!===
contains
!########################################   COUNTCOMCOT
subroutine COUNTCOMCOT(xygrid, xgrdfile, ygrdfile)
implicit none
character(70), intent(in)  :: xgrdfile, ygrdfile
type(grid_xy), intent(out) :: xygrid
integer(4) :: nx,ny
real(8) :: tmp
nx=0 ; ny=0
open(1,file=xgrdfile)
do
  read(1,*,end=10) tmp
  nx=nx+1
end do
10 write(*,'(a,i8,a,a)') "nx=", nx, "file=",trim(xgrdfile) ! 2021.12.08
close(1)
open(2,file=ygrdfile)
do
  read(2,*,end=20) tmp
  ny=ny+1
end do
20 write(*,'(a,i8,a,a)') "ny=", ny, "file=",trim(ygrdfile) ! 2021.12.08
close(2)
xygrid%nx=nx
xygrid%ny=ny
write(*,*) "### COUNTCOMCOT END!! ###"
return
end
!########################################   READCOMCOTGRD
subroutine READCOMCOTGRD(xygrid, xgrdfile, ygrdfile,g_param,g_meshpara) ! 2021.12.08 g_param,g_meshparam is added
use constants
implicit none
type(param_forward),intent(in)    :: g_param    ! 2021.12.08
type(meshpara),     intent(in)    :: g_meshpara ! 2021.12.08
character(70),      intent(in)    :: xgrdfile, ygrdfile
type(grid_xy),      intent(inout) :: xygrid
integer(4)                        :: i,j
real(8)                           :: lonorigin, latorigin,RX ! 2021.12.08

lonorigin = g_meshpara%lonorigin ! 2021.12.08
latorigin = g_meshpara%latorigin ! 2021.12.08

!#[1]## allocate xgrd(1:nx) read
allocate(xygrid%xgrd(xygrid%nx)) ! [km]
allocate(xygrid%ygrd(xygrid%ny)) ! [km]
allocate(xygrid%xgrd_sphe(xygrid%nx)) ! [deg] 2021.12.08
allocate(xygrid%ygrd_sphe(xygrid%ny)) ! [deg] 2021.12.08

!#[2]## read grid
open(1,file=xgrdfile)
do i=1,xygrid%nx
  read(1,*) xygrid%xgrd(i)
end do
close(1)
! ygrd(1:ny) read
open(2,file=ygrdfile)
do i=1,xygrid%ny
  read(2,*) xygrid%ygrd(i)
end do
close(2)

!#[3]# grid conversion from COMCOT coordinate to TMTGEM coordinate
!---------------------------------------------------------------------- comcot is spherical
if ( g_param%iflag_comcot_spherical .eq. 1) then ! 2021.12.08
 ! 0: cartesian (grid is given in meter), 1: spherical (grid is given in degree)! 2021.12.08
 xygrid%xgrd_sphe = xygrid%xgrd ! 2021.12.08 [deg]
 xygrid%ygrd_sphe = xygrid%ygrd ! 2021.12.08 [deg]
 do j=1,xygrid%nx   ! 2021.12.08
  xygrid%xgrd(j)=(xygrid%xgrd(j)-lonorigin)*d2r*earthrad*dcos(latorigin*d2r) ! [km] 2021.12.08
 end do                ! 2021.12.08
 do j=1,xygrid%ny   ! 2021.12.08
  xygrid%ygrd(j)=(xygrid%ygrd(j)-latorigin)*d2r*earthrad ! [km] 2021.12.08
 end do      ! 2021.12.08
!---------------------------------------------------------------------- comcot is cartesian
else if ( g_param%iflag_comcot_spherical .eq. 0 ) then  ! 2021.12.08
 xygrid%xgrd(:)=xygrid%xgrd(:)*1.d-3 ! [m] -> [km]  ! 2021.12.08
 xygrid%ygrd(:)=xygrid%ygrd(:)*1.d-3 ! [m] -> [km]  ! 2021.12.08
 RX=1./d2r/earthrad/dcos(latorigin*d2r)
 do j=1,xygrid%nx   ! 2021.12.08
  xygrid%xgrd_sphe(j) = xygrid%xgrd(j)*RX + lonorigin! 2021.12.08 [deg]
 end do
 RX=1./d2r/earthrad    ! 2021.12.08
 do j=1,xygrid%ny
  xygrid%ygrd_sphe(j) = xygrid%ygrd(j)*RX +latorigin ! 2021.12.08 [deg]
 end do
else       ! 2021.12.08
 write(*,*) "GEGEGE iflag_comcot_spherical should be 0 or 1" ! 2021.12.08
 stop      ! 2021.12.08
end if     ! 2021.12.08

write(*,*) "### READCOMCOTGRD END!! ###"
return
end

!######################################## ALLOCATEGRDDATA
subroutine ALLOCATEGRDDATA(h_grd,xygrid,idim)
implicit none
integer(4),        intent(in)    :: idim
type(grid_xy),     intent(in)    :: xygrid
type(grid_data_2D),intent(inout) :: h_grd
integer(4) :: nx,ny ! 2021.12.02
nx=xygrid%nx        ! 2021.12.02
ny=xygrid%ny        ! 2021.12.02

h_grd%xygrid = xygrid ! deep copy nx,ny,xgrd,ygrd,xgrd_sphe,ygrd_sphe ! 2021.12.08
allocate(h_grd%dat(idim,nx,ny))

return
end
!######################################## ALLOCATEGRDDATA
subroutine DEALLOCATEGRDDATA(h_grd)
implicit none
type(grid_data_2D),intent(inout) :: h_grd

deallocate(h_grd%xygrid%xgrd)
deallocate(h_grd%xygrid%ygrd)
deallocate(h_grd%dat)

return
end
!
!########################################   READFXYZ
subroutine READGRDDATA(infile,h_grd,ncomp)
implicit none
type(grid_data_2D),intent(inout) :: h_grd
integer(4),        intent(in)    :: ncomp
character(70),     intent(in)    :: infile
integer(4) :: i,j,k
real(8) :: f
!#[1]## Fxyz read
write(*,*) "infile=",infile
open(1,file=infile)
do j=1,h_grd%xygrid%ny
  do i=1,h_grd%xygrid%nx
     read(1,*) ( h_grd%dat(k,i,j),k=1,ncomp ) !,f
!    read(1,*) Fxyz_grd%dat(2,i,j), Fxyz_grd%dat(1,i,j), Fxyz_grd%dat(3,i,j),f
!    write(*,*) Fxyz_grd%dat(2,i,j), Fxyz_grd%dat(1,i,j), Fxyz_grd%dat(3,i,j),f
  end do
end do
close(1)
!Fxyz_grd%dat(3,:,:)=-Fxyz_grd%dat(3,:,:)

write(*,*) "### READGRDDATA END!! ###"
return
end
!##########################################################
! nlayer is added and z_grd is added on 2021.12.02
subroutine READMN(it, vxyh_grd, z_grd, g_param)  ! inserted on Oct. 25, 2015
implicit none
integer(4),         intent(in)    :: it
type(param_forward),intent(in)    :: g_param
type(grid_data_2D), intent(inout) :: vxyh_grd(g_param%nlayer) ! 2021.12.02
type(grid_data_2D), intent(inout) :: z_grd(   g_param%nlayer) ! 2021.12.02
character(70)                     :: mpre, npre, zpre   ! zpre is added 2021.12.02
character(70)                     :: mfile, nfile,zfile ! zfile added 2021.12.02
integer(4)                        :: itime_comcot, i, j,nx,ny
integer(4)                        :: nlayer,ilayer      ! 2021.12.02
character(6)                      :: num
real(8)                           :: dt_comcot          ! 2018.11.15

nlayer = g_param%nlayer ! 2021.12.02
do ilayer=1,nlayer      ! 2021.12.02

!#[0]## set input
  mpre      = g_param%mpre(ilayer) ! 2021.12.02
  npre      = g_param%npre(ilayer) ! 2021.12.02
  zpre      = g_param%zpre(ilayer) ! 2021.12.02
  dt_comcot = g_param%dt_comcot    ! 2018.11.15

!#[1]## generate file name
  itime_comcot = it*NINT(dt_comcot/g_param%dt_sim_comcot)  ! 2024.01.19 TT
  !itime_comcot= it*dt_comcot  ! 2018.11.15
  nx           = vxyh_grd(ilayer)%xygrid%nx
  ny           = vxyh_grd(ilayer)%xygrid%ny
  write(num,'(i6.6)') itime_comcot
  mfile = mpre(1:len_trim(mpre))//num(1:6)//".dat"
  nfile = npre(1:len_trim(npre))//num(1:6)//".dat"
  zfile = zpre(1:len_trim(zpre))//num(1:6)//".dat"    ! 2021.12.02
  write(*,'(a,a)') " mfile=",mfile(1:len_trim(mfile)) ! 2019.01.23
  write(*,'(a,a)') " nfile=",nfile(1:len_trim(nfile)) ! 2019.01.23
  write(*,'(a,a)') " zfile=",zfile(1:len_trim(zfile)) ! 2021.12.02

!#[2]## read vxh, vyh
open(1,file=mfile)
open(2,file=nfile)
open(3,file=zfile) ! 2021.12.02
do j=1,ny
!  write(*,*) "j=",j,"ny=",ny
  read(1,'(15f9.3)') (vxyh_grd(ilayer)%dat(1,i,j),i=1,nx) ! 2021.12.02
  read(2,'(15f9.3)') (vxyh_grd(ilayer)%dat(2,i,j),i=1,nx) ! 2021.12.02
  read(3,'(15f9.3)') (z_grd(   ilayer)%dat(1,i,j),i=1,nx) ! 2021.12.02
end do
close(1)
close(2)
close(3) !2021.12.02

end do   ! 2021.12.02

write(*,*) "### READMN END!! ###"
return
end
!
!
!#################################################
subroutine GENLVEC3COM(Fxyz, nodes, g_node, nodtot, Fxyzip, ip)
implicit none
integer(4),intent(in) :: nodes, ip
real(8),intent(in) :: Fxyz(3, nodes)
integer(4),intent(in) :: nodtot, g_node(nodtot)
real(8), intent(out) :: Fxyzip(3, nodtot)
integer(4) :: i
do i=1,nodtot
  if ( g_node(i) .le. nodes) then
    Fxyzip(1:3,i) = Fxyz(1:3, g_node(i))
  else
    Fxyzip(1:3,i)=0.d0
  end if
end do
write(*,*) "### GENLVEC3COM END!! ### ip=",ip
return
end

end module FROMCOMCOT
