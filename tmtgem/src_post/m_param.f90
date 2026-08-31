module param
implicit none

type param_forward
 ! mesh info
 integer(4)    :: itopoflag    ! 0: no topo, 1: topo are from topofiles 2017.09.29
 integer(4)    :: nfile        ! 2017.09.27
 character(50),allocatable,dimension(:)   :: topofile ! 2017.09.27
 real(8),      allocatable,dimension(:,:) :: lonlatshift  ! 2017.09.27
 character(50) :: g_meshfile   ! global mesh file
 integer(4)    :: surface_id_ground=6 ! surface consisting of triangles for ground
 character(50) :: z_meshfile    ! 2d triangle mesh file
 character(50) :: g_lineinfofile ! added on April 21, 2016
 character(50) :: outputfolder
 character(50) :: header2d
 character(50) :: header3d
 ! frequency info
 integer(4)    :: nfreq ! # of frequency for this
 real(8),allocatable,dimension(:)   :: freq ! frequencies
 ! observatory info
 integer(4)    :: nobs  ! # of observatories
 integer(4)    :: lonlatflag !# 1: lonlat [deg], 2: xy [km]
 real(8)       :: wlon,elon,slat,nlat ! [deg]
 real(8)       :: lonorigin,latorigin ! [deg]
 real(8)       :: lenout    ! [km]
 real(8)       :: upzin     ! [km]
 real(8)       :: downzin   ! [km]
 real(8)       :: zmax      ! [km]
 real(8)       :: zmin      ! [km]
 real(8)       :: sizein    ! [km]
 real(8)       :: sizebo    ! [km]
 real(8)       :: sigma_obs ! [km] ! converted to dimension
 real(8)       :: A_obs     ! [km] ! on2016.10.20
 real(8)       :: dlen_source ! [km]  2017.10.12
 real(8)       :: sigma_src ! [km]
 real(8)       :: A_src     ! [km]
 integer(4)    :: nobsr
 real(8),      allocatable,dimension(:,:) :: xyz_r
 real(8),      allocatable,dimension(:)   :: A_r,sigma_r
 real(8)                                  :: xbound(4),ybound(4),zbound(4)
 character(3)                             :: UTM ! UTM zone
 real(8),      allocatable,dimension(:,:) :: lonlataltobs!lon[deg],lat[deg],alt[km]
 real(8),      allocatable,dimension(:,:) :: xyzobs ! x eastward[km],y northward[km],z[km]
 character(50),allocatable,dimension(:)   :: obsname
 !# output spatial distribution 2017.10.11
 integer(4)    :: ixyflag      ! 0: nothing, 1: (nx,ny) surface values, 2: triangle surface
 character(50) :: xyfilehead   ! 2017.10.11 (ixyflag = 1 or 2)
 integer(4)    :: nx,ny        ! 2017.10.11 (ixyflag = 1     )
 !#
 real(8)       :: xyzminmax(6) ! added on 2017.02.21
 real(8)       :: zorigin=0.d0 ! added on 2017.02.21
 !# conductivity structure
 character(50) :: condfile
 integer(4)    :: condflag
end type

type param_source
 integer(4) :: lonlatflag !# 1:lonlat [deg], 2: xy [km]
 integer(4) :: nsource    ! # of sources  added on 2017.07.11
 character(50),allocatable,dimension(:) :: sourcename ! 2017.07.11
 real(8),allocatable,dimension(:,:) :: xs1  ! [km] xyz coordinate of start electrode
 real(8),allocatable,dimension(:,:) :: xs2  ! [km] xyz coordinate of end electrode
 real(8),allocatable,dimension(:,:) :: lonlats1 ! [deg] lonlat of start electrode
 real(8),allocatable,dimension(:,:) :: lonlats2 ! [deg] xyz coordinate of end electrode
 real(8)                            :: I ! source current through wire [A]
end type

type param_cond
 integer(4)    :: condflag ! 0: (partitioned) homegeneous earth, 1: file is given
 character(50) :: condfile
 integer(4)    :: ntet
 integer(4)    :: nphys1 ! # of elements of air
 integer(4)    :: nphys2 ! # of elements in land
 real(8)       :: sigma_air=1.d-8 ! [S/m]
 integer(4),allocatable,dimension(:) :: index ! element id for nphys 2, 2017.05.10
 real(8),   allocatable,dimension(:) :: sigma ! [S/m]
 real(8),   allocatable,dimension(:) :: rho   ! [Ohm.m]
 integer(4)    :: nvolume     ! # of region for land  ( condflag = 0 )2017.09.28
 real(8),   allocatable,dimension(:) :: sigma_land  ! ( condflag = 0 )2017.09.28
end type

!type param_bell ! inserted on 2016.10.12 for Kusatsushirane
! character(50) :: bellgeofile
! integer(4) :: ilonlatflag ! 1 for lonlat, 2 for xy [km]
! real(8)  :: lon_bell, lat_bell
! real(8)  :: radius![km]
! real(8)  :: width ![km]
! real(8)  :: ztop  ! < 0 [km] how deep is the top of bell
! real(8)  :: zbot  ! < 0 [km] how deep is the bottom of bell
! real(8)  :: zthic ! [km] thickness of top of the bell
! real(8)  :: reso_bell ! [km] mesh resolution for inside bell
!--
! real(8)  :: xyz_bell(3)
!end type

contains
!################################################# duplicatescalarcond
!# Coded on 2017.05.14
subroutine duplicatescalarcond(g_cond,h_cond)
implicit none
type(param_cond),intent(in)  :: g_cond
type(param_cond),intent(out) :: h_cond

h_cond%condflag   = g_cond%condflag   ! 0: homegeneous earth, 1: file is given
h_cond%condfile   = g_cond%condfile
h_cond%nvolume    = g_cond%nvolume    !       2017.09.28
allocate(h_cond%sigma_land(h_cond%nvolume)) ! 2017.09.28
h_cond%sigma_land = g_cond%sigma_land ! (used [S/m] (only when condflag = 0 )
h_cond%ntet       = g_cond%ntet
h_cond%nphys1     = g_cond%nphys1     ! # of elements of air
h_cond%nphys2     = g_cond%nphys2     ! # of elements in land
h_cond%sigma_air  = g_cond%sigma_air  ! [S/m] this will be also used in inversion 2018.10.04


return
end subroutine

!################################################## subroutine readparam
!# change the
subroutine readparam(c_param,sparam,g_cond)
implicit none
type(param_forward),          intent(out) :: c_param
type(param_source),           intent(out) :: sparam
type(param_cond),   optional, intent(out) :: g_cond
integer(4)                  :: i,j,nobs,input=2,nsource
character(50)               :: site
real(8)                     :: lonorigin,latorigin,a
character(100)              :: paramfile
integer(4),    parameter    :: n = 1000 ! 2020.09.28
character(200),dimension(n) :: lines    ! 2020.09.28
!#
write(*,*) ""
write(*,*) "<Please input the forward parameter file>" ! 2020.09.28
read(*,'(a)') paramfile           ! 2020.09.28
call readcontrolfile(paramfile,n,lines) ! 2020.09.28

open(input,file="tmp.ctl")
!#[2]# read mesh_param file
write(*,*) ""
write(*,*) "<input 0 for no topo, 1 for topofiles (itopoflag)>" ! 2018.10.02
read(input,*)    c_param%itopoflag      ! 11->* 2021.09.02
write(*,'(a,i3)') " itopoflag =", c_param%itopoflag    ! 2020.09.17

if ( c_param%itopoflag .eq. 0 ) then ! 2017.09.29
 goto 101 ! 2017.09.29
else if (c_param%itopoflag .eq. 1 ) then ! 2017.09.29
  read(input,*) c_param%nfile            ! 12->* 2021.09.02
  allocate(c_param%topofile(c_param%nfile) )            ! 2017.09.27
  allocate(c_param%lonlatshift(2,c_param%nfile) )       ! 2017.09.27
 do i=1,c_param%nfile                                  ! 2017.09.27
  read(input,10) c_param%topofile(i)
  read(input,*) c_param%lonlatshift(1:2,i)             ! 2021.09.28
 end do                                                ! 2017.09.27
else         ! 2017.09.29
 write(*,*) "GEGEGE! itopoflag should be 0 or 1: itopoflag",c_param%itopoflag!2017.09.29
 stop        ! 2017.09.29
end if       ! 2017.09.29
101 continue ! 2017.09.29
write(*,*) ""
write(*,   40) " <Please enter global meshfile>"
read(input,10) c_param%g_meshfile
write(*,   41) " global mesh file : ",trim(c_param%g_meshfile)     ! 2020.09.29
read(input,10) c_param%z_meshfile
write(*,   41) " 2dz    mesh file : ",trim(c_param%z_meshfile)     ! 2020.09.29
read(input,10) c_param%g_lineinfofile
write(*,   41) " line info   file : ",trim(c_param%g_lineinfofile) ! 2020.09.29
read(input,10) c_param%outputfolder
write(*,   41) " output folder    : ",trim(c_param%outputfolder)   ! 2020.09.29
read(input,10) c_param%header2d
write(*,   41) "header2d          : ", trim(c_param%header2d)      ! 2020.10.31
read(input,10) c_param%header3d
write(*,   41) "header3d          : ", trim(c_param%header3d)      ! 2020.10.31

!# frequency info
read(input,*) c_param%nfreq          ! 12->* 2021.09.02
write(*,'(a,i3)') " # of frequency : ",c_param%nfreq
allocate(c_param%freq(c_param%nfreq))
do i=1, c_param%nfreq
 read(input,*) c_param%freq(i)        ! 12->* 2021.09.02
 write(*,'(f9.4,a)')     c_param%freq(i)," [Hz]" ! 2020.09.29
end do

!# xy plane boundary by xy coordinate [km]
 write(*,*) ""
 write(*,*) "< input xbound and ybound >"    ! 2020.09.29
 read(input,*) c_param%xbound(2)             ! 12->* 2021.09.02
 read(input,*) c_param%xbound(3)             ! 12->* 2021.09.02
 read(input,*) c_param%ybound(2)             ! 12->* 2021.09.02
 read(input,*) c_param%ybound(3)             ! 12->* 2021.09.02
 write(*,43) " xbound(2)=",c_param%xbound(2) ! 2020.09.29
 write(*,43) " xbound(3)=",c_param%xbound(3) ! 2020.09.29
 write(*,43) " ybound(2)=",c_param%ybound(2) ! 2020.09.29
 write(*,43) " ybound(3)=",c_param%ybound(3) ! 2020.09.29

!# mesh info
 write(*,*) ""
 write(*,*) "< input meshinfo >" ! 2020.09.29
 read(input,*) c_param%lenout    ! 12->* 2021.09.02
 read(input,*) c_param%upzin     ! 12->* 2021.09.02
 read(input,*) c_param%downzin   ! 12->* 2021.09.02
 read(input,*) c_param%zmax      ! 12->* 2021.09.02
 read(input,*) c_param%zmin      ! 12->* 2021.09.02
 read(input,*) c_param%sizein    ! 12->* 2021.09.02
 read(input,*) c_param%sizebo    ! 12->* 2021.09.02
 read(input,*) c_param%sigma_obs ! commented out on 2016.10.20 ! 12->* 2021.09.02
 read(input,*) c_param%A_obs     ! commented out on 2016.10.20 ! 12->* 2021.09.02
 read(input,*) c_param%dlen_source ! [km] 2017.10.12 ! 12->* 2021.09.02
 read(input,*) c_param%sigma_src ! 12->* 2021.09.02
 read(input,*) c_param%A_src     ! 12->* 2021.09.02
 write(*,43) " lenout      =",c_param%lenout ! 2020.09.29
 write(*,43) " upzin       =",c_param%upzin
 write(*,43) " downzin     =",c_param%downzin
 write(*,43) " zmax        =",c_param%zmax
 write(*,43) " zmin        =",c_param%zmin
 write(*,43) " sizein      =",c_param%sizein
 write(*,43) " sizebo      =",c_param%sizebo
 write(*,43) " sigma_obs   =",c_param%sigma_obs ! commented out on 2016.10.20
 write(*,43) " A_obs       =",c_param%A_obs     ! commented out on 2016.10.20
 write(*,43) " dlen_source =",c_param%dlen_source ! [km] 2017.10.12
 write(*,43) " sigma_src   =",c_param%sigma_src
 write(*,43) " A_src       =",c_param%A_src

 !# calculate xbound(1:4), ybound(1:4), zbound(1:4)
 c_param%zbound(1:4) = (/c_param%zmin,c_param%downzin,c_param%upzin,c_param%zmax/)
 c_param%xbound(1) =   c_param%xbound(2) - c_param%lenout
 c_param%xbound(4) = - c_param%xbound(1)
 c_param%ybound(1) =   c_param%ybound(2) - c_param%lenout
 c_param%ybound(4) = - c_param%ybound(1)
 write(*,*) "" ! 2020.09.29
 write(*,'(a,4f9.4)') " xbound =",c_param%xbound(1:4) ! 2021.09.02
 write(*,'(a,4f9.4)') " ybound =",c_param%ybound(1:4) ! 2021.09.02
 write(*,'(a,4f9.4)') " zbound =",c_param%zbound(1:4) ! 2021.09.02
 write(*,*) "" ! 2021.09.02
!# observatory info
 read(input,*) c_param%nobs ! 11->* 2021.09.02
 write(*,'(a,i3)') " # of observatories : ",c_param%nobs ! 2020.09.17
 allocate(c_param%lonlataltobs(3,c_param%nobs))
 allocate(c_param%xyzobs      (3,c_param%nobs))
 allocate(c_param%obsname(c_param%nobs))
! allocate(c_param%A_obs(c_param%nobs))       ! added on 2016.10.20
! allocate(c_param%sigma_obs(c_param%nobs))   ! added on 2016.10.20

!# lonlatflag
 read(input,*)  c_param%lonlatflag    ! 1: lonlatalt, 2:xyz ! 11->* 2021.09.02
 if ( c_param%lonlatflag .eq. 1) then ! lonlat
  read(input,*) c_param%UTM           ! 10->* 2021.09.02
  write(*,*) "UTM zone : ",c_param%UTM ! 2020.09.29
  read(input,*) c_param%lonorigin,c_param%latorigin ! 12->* 2021.09.02
  write(*,*) "lonorigin =",c_param%lonorigin ! 2020.09.17
  write(*,*) "latorigin =",c_param%latorigin  ! 2020.09.17
 end if
 lonorigin = c_param%lonorigin
 latorigin = c_param%latorigin
 nobs = c_param%nobs

!# read site information
if (c_param%lonlatflag .eq. 1) write(*,'(a)') "< conversion of Lon Lat to UTM x (east), y (north) [km]>"!2021.09.29
write(*,*) "" !2021.09.29

 do i=1, nobs
   read(input,10) c_param%obsname(i)
   site=c_param%obsname(i) ! 2021.09.02

   ! lonlatalt
   if (c_param%lonlatflag .eq. 1 ) then
    read(input,*) (c_param%lonlataltobs(j,i),j=1,3)
    write(*,'(1x,a,a,3f15.7)') trim(site)," :",c_param%lonlataltobs(1:3,i) !2021.09.02
    call UTMXY(c_param%lonlataltobs(1:2,i),&
        & lonorigin,latorigin,c_param%xyzobs(1:2,i),c_param%UTM)
    c_param%xyzobs(3,i) = c_param%lonlataltobs(3,i)
    write(*,'(1x,a,2f15.7,a)') " UTM>",c_param%xyzobs(1:2,i)," [km]" ! 2021.09.29
    write(*,*) ""  ! 2021.09.29

   ! xyz
   else if (c_param%lonlatflag .eq. 2 ) then ! xyz
    read(input,*) (c_param%xyzobs(j,i),j=1,3)

   end if
  end do

!#[2]# xy reading  2017.10.11
 write(*,*) ""
 write(*,*) "< input ixyflag: 0 for nothing, 1 for xy plan view >"
 read(input,*) c_param%ixyflag    ! 12->* 2021.09.02
 write(*,'(a,i3)') " ixyflag =",c_param%ixyflag ! 2020.09.17
 if ( c_param%ixyflag .eq. 1) then ! 2017.10.11
  read(input,*) c_param%xyfilehead     ! 2021.09.02
  read(input,*) c_param%nx,c_param%ny  ! 2021.09.02
  write(*,*) "nx,ny =",c_param%nx,c_param%ny
  write(*,*) "xyfilehead = ",c_param%xyfilehead
 end if                                ! 2017.10.11
 if ( c_param%ixyflag .eq. 2) then     ! 2017.10.11
  read(input,*) c_param%xyfilehead     ! 2021.09.02
 end if                                ! 2017.10.11

!#[3]# read source param
 sparam%lonlatflag = c_param%lonlatflag

 ! # of sources 2017.07.11
 write(*,*) ""
 write(*,*) "< input # of source wires > "!2017.07.11
 read(input,*) sparam%nsource       !2017.07.11, ! 11->* 2021.09.02
 nsource = sparam%nsource
 write(*,'(a,i3)') " # of source wires (nsource) =",nsource ! 2020.09.17
 allocate(sparam%xs1(3,nsource),    sparam%xs2(3,nsource))!2017.07.11
 allocate(sparam%lonlats1(2,nsource),sparam%lonlats2(2,nsource))!2017.07.11
 allocate(sparam%sourcename(nsource))

 ! lonlat input
 do i=1,nsource ! 2017.07.11
  read(input,*) sparam%sourcename(i) ! 10->* 2021.09.02
  write(*,41) " source name: ",sparam%sourcename(i) ! 2020.09.29
  if (sparam%lonlatflag .eq. 1) then ! lonlat
   read(input,*) sparam%lonlats1(1:2,i),sparam%xs1(3,i)
   read(input,*) sparam%lonlats2(1:2,i),sparam%xs2(3,i)
   write(*,44) " Start point (lon,lat,z)=",sparam%xs1(1:3,i)
   write(*,44) " End   point (lon,lat,z)=",sparam%xs2(1:3,i)
   call UTMXY(sparam%lonlats1(1:2,i),lonorigin,latorigin,sparam%xs1(1:2,i),c_param%UTM)
   call UTMXY(sparam%lonlats2(1:2,i),lonorigin,latorigin,sparam%xs2(1:2,i),c_param%UTM)
   write(*,*) " converted with UTM to :"
   write(*,'(a,3f15.7)') " xs1 (x,y,z)=",sparam%xs1(1:3,i) !2020.09.29
   write(*,'(a,3f15.7)') " xs2 (x,y,z)=",sparam%xs2(1:3,i) !2020.09.29
 !# xy input
 else if (sparam%lonlatflag .eq. 2) then ! xyz
  read(input,*) sparam%xs1(1:3,i)
  read(input,*) sparam%xs2(1:3,i)

 else
  goto 99
 end if
 end do ! 2017.07.11

 read(input,*) sparam%I ! [A] ! 11->* 2021.09.02
 write(*,'(a,g15.7,a)') " Electric source current :",sparam%I," [A]" ! 2020.09.29
!#[4]## lonlatflag : 1 -> 2 because xyzobs is laready set
 c_param%lonlatflag = 2


!#[5]## read conductivity information
 if ( present(g_cond) ) then     ! if the optional argument, g_cond, is present
  read(input,*) g_cond%sigma_air ! 11->* 2021.09.02
  write(*,'(a,g15.7,a)') " sigma_air =",g_cond%sigma_air," [S/m]"
  read(input,*) g_cond%condflag  ! 0:homogeneous, 1:file
  if (g_cond%condflag .eq. 0 ) then                         ! 2017.09.29
   write(*,*) "" ! 2020.09.29
   write(*,*) "<Input # of physical volumes in land region>"
   read(input,*) g_cond%nvolume ! # of physical volume in land, ! 11->* 2021.09.02
   write(*,*) "nvolume=",g_cond%nvolume
   allocate( g_cond%sigma_land(g_cond%nvolume) )                  ! 2017.09.28
   write(*,*) "" ! 2020.09.29
   write(*,*) "<Inuput land sigma [S/m] for each physical volume>"  ! 2017.09.28
   do i=1,g_cond%nvolume
    read(input,*) g_cond%sigma_land(i) ! conductivity in land region 2017.09.28, ! 12->* 2021.09.02
    write(*,*) i,"sigma_land=",g_cond%sigma_land(i),"[S/m]"
   end do
  else if (g_cond%condflag .eq. 1) then ! file
   read(input,'(a50)') g_cond%condfile  ! 2021.10.04
   write(*,*) "cond file is",g_cond%condfile
   CALL READCOND(g_cond)          ! read conductivity structure
  else
   write(*,*) "GEGEGE condflag should be 0 or 1 : condflag=",g_cond%condflag
   stop
  end if
 end if

 close(input)
 write(*,*) "### READ FORWAR PARAM END!! ###"

return

99 continue
write(*,*) "GEGEGE!"
write(*,*) "c_param%lonlatflag=",c_param%lonlatflag
write(*,*) "sparam%lonlatflag=",sparam%lonlatflag
stop
! in 10 to 32, "20x" is omitted, initial 20 characters are omitted in control file in readcontrol 2021.09.02
10 format(a)      ! 2021.09.02
40 format(a)
41 format(a,a)
43 format(a,f9.4)
44 format(a,3f9.4)
end subroutine readparam

!----------------------------------------- READCOND
subroutine READCOND(g_cond)
implicit none
type(param_cond),     intent(inout)  :: g_cond
real(8),   allocatable,dimension(:)  :: rho
integer(4),allocatable,dimension(:)  :: index
integer(4) :: nphys2
integer(4) :: i,input=11

!#[1]## read condfile
 write(*,*) "condfile : ",g_cond%condfile ! 2020.07.19
 open(input,file=g_cond%condfile)
!# header
 do i=1,11      ! 2017.09.11 changed from 8 to 11
  read(input,*)
 end do
  read(input,*) nphys2
 allocate(rho(nphys2),index(nphys2)) !!! Important !!!!
 do i=1,nphys2
  read(input,*) index(i),rho(i)
 end do
 close(input)

!#[2]## set g_cond
 allocate(g_cond%rho(  nphys2))
 allocate(g_cond%sigma(nphys2))
 allocate(g_cond%index(nphys2)) ! added on 2017.05.10
 g_cond%rho    = rho
 g_cond%nphys2 = nphys2
 g_cond%sigma  = -9999. ! 2017.11.06
 do i=1,nphys2
  if (abs(rho(i)) .gt. 1.d-10 ) g_cond%sigma(i) = 1.d0/rho(i) ! 2017.11.06
!  write(*,*) i,"g_cond%sigma(i)=",g_cond%sigma(i)
 end do
write(*,*) "g_cond%nphys2 =",g_cond%nphys2 ! 2021.01.08
write(*,*) "### READCOND  END!! ###"       ! 2020.09.29
return
end subroutine
!----------------------------------------- deallocatecond
! Coded on 2017.05.14
subroutine deallocatecond(g_cond)
implicit none
type(param_cond),intent(inout) :: g_cond

g_cond%ntet   = 0
g_cond%nphys1 = 0! # of elements of air
g_cond%nphys2 = 0! # of elements in land
if (allocated(g_cond%index)) deallocate(g_cond%index) ! element id for nphys 2, 2018.10.05
if (allocated(g_cond%sigma)) deallocate(g_cond%sigma) ! [S/m]   2018.10.05
if (allocated(g_cond%rho)  ) deallocate(g_cond%rho  ) ! [Ohm.m] 2018.10.05

return
end subroutine



!#################  subroutine calobsr ################################################# calobsr
! Coded on 2016.10.14
! generate virtual observatories for mesh refinement along source dipole 2021.09.29
subroutine calobsr(s_param,g_param)
!use param
implicit none
type(param_forward),intent(inout)       :: g_param
type(param_source), intent(in)          :: s_param
type source !========================= type source definition 2017.07.11
 integer(4)                             :: nobsr_per_src ! 2021.09.29
 real(8),    allocatable,dimension(:,:) :: xyz_src
end type
!===================================== type source definition end 2021.09.29
type(source),allocatable,dimension(:)   :: ss ! 2017.07.11
real(8),     allocatable,dimension(:,:) :: xs1,xs2
real(8),     allocatable,dimension(:,:) :: xyz
real(8),     allocatable,dimension(:,:) :: xyz_r
real(8),     allocatable,dimension(:)   :: sigma_r,A_r
real(8)                                 :: len_source,xs12(3)
real(8)                                 :: dlen_source
real(8)                                 :: sigma_obs,A_obs
real(8)                                 :: sigma_src,A_src
integer(4)                              :: nobsr,nobs,i,j,icount
integer(4)                              :: nsr ! # of source wires 2017.07.11

write(*,'(a)') " ### CALOBSR START!! ###" ! 2021.09.29

!#[1]## set
 nobs        = g_param%nobs
 nsr         = s_param%nsource   ! 2017.07.11 # of dipole sources
 allocate(xs1(3,nsr),xs2(3,nsr)) ! 2017.07.11
 allocate(xyz(3,nobs))           ! 2017.07.11
 allocate(ss(nsr))               ! 2017.07.11
 xs1         = s_param%xs1
 xs2         = s_param%xs2
 dlen_source = g_param%dlen_source ! 2017.10.12
 xyz         = g_param%xyzobs
 sigma_obs   = g_param%sigma_obs
 A_obs       = g_param%A_obs
 sigma_src   = g_param%sigma_src
 A_src       = g_param%A_src


 !#[2]##
 do j=1,nsr ! 2017.07.11
  xs12(1:3) = xs2(1:3,j) - xs1(1:3,j)
  len_source = sqrt(xs12(1)**2. + xs12(2)**2. + xs12(3)**2. )
  write(*,*) "len_source =",len_source,"[km]"
  ss(j)%nobsr_per_src = int(len_source/dlen_source) +1  ! 2021.09.29
!  ss(j)%nsource = 1     ! commented out on 2017.10.12
  write(*,*) "nsource    =",ss(j)%nobsr_per_src,"for source",j ! 2021.09.29
  allocate( ss(j)%xyz_src(3,ss(j)%nobsr_per_src) ) !2021.0.29
  if (ss(j)%nobsr_per_src .gt. 1) then ! 2021.09.29
  do i=1, ss(j)%nobsr_per_src ! 2021.09.29
   ss(j)%xyz_src(1:3,i) = xs1(1:3,j) + float(i-1)/float(ss(j)%nobsr_per_src - 1)*xs12(1:3)!2021.09.29
  end do
  else if ( ss(j)%nobsr_per_src .eq. 1) then ! 2021.09.29
   ss(j)%xyz_src(1:3,1) = xs1(1:3,j) + xs12(1:3)/2.
  end if
 end do ! 2017.07.11

 !#[3]## calculate
  nobsr = nobs
 do j=1,nsr ! 2017.07.11
  nobsr = nobsr + ss(j)%nobsr_per_src ! 2021.09.29
 end do
 allocate(xyz_r(3,nobsr),A_r(nobsr),sigma_r(nobsr))
 xyz_r(:,1:nobs)       = xyz(:,1:nobs)
 icount=0                     ! 2017.07.11
 do j=1,nsr                   ! 2017.07.11
  do i=1,ss(j)%nobsr_per_src  ! 2021.09.29
   icount = icount + 1        ! 2017.07.11
   xyz_r(:,nobs+icount) = ss(j)%xyz_src(:,i) ! 2017.07.11
  end do                      ! 2017.07.11
 end do                       ! 2017.07.11
 A_r(     1:nobs)      = A_obs
 A_r(nobs+1:nobsr)     = A_src
 sigma_r(     1:nobs ) = sigma_obs
 sigma_r(nobs+1:nobsr) = sigma_src

!#[4]## set output
allocate(g_param%xyz_r(3,nobsr),g_param%A_r(nobsr),g_param%sigma_r(nobsr))
g_param%nobsr   = nobsr
g_param%xyz_r   = xyz_r
g_param%A_r     = A_r
g_param%sigma_r = sigma_r
write(*,*) "" ! 2021.09.29
write(*,*) "< observatories for receiver points >"!2021.09.29
write(*,*) " - Note that z is relative to the ground surface - " ! 2021.09.29
do i=1,nobs
 write(*,'(i3,1x,a3,a,3f14.7,a)') i,g_param%obsname(i)," (x,y,z)=",g_param%xyz_r(1:3,i)," [km]"!sigma_r(i),A_r(i) 2021.09.21
end do
write(*,*) "" ! 2021.09.30
write(*,*) "< virtual observatories along source wire for mesh refinement >"!2021.09.29
do i=nobs+1,nobsr
 write(*,'(i3,a,3f14.7,a)') i,"     (x,y,z)=",g_param%xyz_r(1:3,i)," [km]"!sigma_r(i),A_r(i) 2021.09.21
end do
write(*,*) "" !2021.09.29

write(*,'(a)') " ### CALOBSR END!! ###"

return
end subroutine
!######################################################### UTMXY
!# Coded on 2016.10.12 by T.MINAMI
subroutine UTMXY(lonlat,lonorigin,latorigin,xyout,zone)
implicit none
real(8),intent(in) :: lonlat(2),lonorigin,latorigin
character(3),intent(in) :: zone
real(8) ,intent(out) :: xyout(2)
real(8) :: xorigin,yorigin,x,y

call UTMGMT(lonlat(1),lonlat(2), x,       y,      zone,0)
call UTMGMT(lonorigin,latorigin, xorigin, yorigin,zone,0)
!write(*,*) "lon=",lon,"lonorigin=",lonorigin,"lat=",lat,"latorigin=",latorigin
!write(*,*) "x=",x,"xorigin=",xorigin,"y=",y,"yorigin=",yorigin

xyout(1) = (x - xorigin)/1.d3 ! [km]
xyout(2) = (y - yorigin)/1.d3 ! [km]

!write(*,*) "xout=",xout,"yout=",yout

return
end subroutine
!######################################################### UTMGMT
subroutine UTMGMT(xin,yin,xout,yout,zone,iflag)
implicit none
real(8),intent(in) :: xin,yin
real(8),intent(out) :: xout,yout
character(3),intent(in) :: zone
integer(4),intent(in) :: iflag ! 0: LONLAT2UTM, 1:UTM2LONLAT

integer(4) ::  izone          ! 20200728
read(zone(1:2),*) izone       ! 20200728

!#[1]## prepare values
!write(values,'(g18.10,1x,g18.10)') xin,yin ! commented out 20200728

!#[2]## use mapproject
if (iflag .eq. 0 ) then ! LONLAT2UTM
! CALL system("echo "//values//" | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C > tmp.dat")! commented out 20200728
 CALL utm_geo_tm(xin,yin,xout,yout,izone,iflag) ! 20200728
else if (iflag .eq. 1 ) then ! UTM2LONLAT
! CALL system("echo "//values//" | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C -I > tmp.dat")! commented out 20200728
 CALL utm_geo_tm(xout,yout,xin,yin,izone,iflag) ! 20200728
else
 write(*,*) "GEGEGE! iflag should be 0 (LONLAT2UTM) or 1(UTM2LONLAT), iflag=",iflag
 stop
end if

!#[3]## read xout and yout ! commented out 20200728
! open(12,file="tmp.dat")
!  read(12,*) xout,yout
! close(12)
! call system("rm tmp.dat")

return
end subroutine

!######################################################### UTMGMT_N
subroutine UTMGMT_N(n,xin,yin,xout,yout,zone,iflag)
implicit none
integer(4),intent(in) :: n
real(8),   intent(in) :: xin(n),yin(n)
real(8),   intent(out) :: xout(n),yout(n)
character(3),intent(in) :: zone
integer(4),intent(in) :: iflag ! 0: LONLAT2UTM, 1:UTM2LONLAT
integer(4) :: i

integer(4) ::  izone          ! 20200728
read(zone(1:2),*) izone       ! 20200728
write(*,*) "izone",izone

!#[1]## prepare input file ! commented out 20200728
!open(11,file="in.dat")
! write(11,'(2g18.10)') (xin(i),yin(i),i=1,n)
!close(11)

!#[2]## use mapproject
if (iflag .eq. 0 ) then ! LONLAT2UTM
 ! CALL system("cat in.dat | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C > tmp.dat")! commented out 20200728
 do i=1,n ! 20200728
  CALL utm_geo_tm(xin(i),yin(i),xout(i),yout(i),izone,iflag) ! 20200728
 end do   ! 20200728
else if (iflag .eq. 1 ) then ! UTM2LONLAT
 !CALL system("cat in.dat | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C -I > tmp.dat") ! commented out 20200728
 do i=1,n ! 20200728
  CALL utm_geo_tm(xout(i),yout(i),xin(i),yin(i),izone,iflag) ! 20200728
 end do   ! 20200728
else
 write(*,*) "GEGEGE! iflag should be 0 (LONLAT2UTM) or 1(UTM2LONLAT), iflag=",iflag
 stop
end if

!#[3]## read xout and yout   commented out 20200728
! open(12,file="tmp.dat")
!  do i=1,n
!   read(12,*) xout(i),yout(i)
!  end do
! close(12)
! call system("rm tmp.dat")

return
end subroutine

!######################################################### showcond
! Coded on 2017.05.18
subroutine showcond(g_cond,n)
implicit none
integer(4),intent(in) :: n
type(param_cond),intent(in) :: g_cond
integer(4) :: i

write(*,*) "condflag=",g_cond%condflag
write(*,*) "condfile=",g_cond%condfile
write(*,*) "sigma_land",g_cond%sigma_land
write(*,*) "ntet=",g_cond%ntet
write(*,*) "nphys1=",g_cond%nphys1
write(*,*) "nphys2=",g_cond%nphys2
write(*,*) "sigma_air=",g_cond%sigma_air
do i=1,n!g_cond%nphys2
 write(*,10) i,"index=",g_cond%index(i),&
 &             " sigma",g_cond%sigma(i)," rho=",g_cond%rho(i)
end do

return
10 format(i3,a,i7,a,g15.7,a,g15.7)
end subroutine

!##################################################################### utm_geo
!  Added by Takuto Minami on July 28, 2020
!================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

  subroutine utm_geo_tm(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway) ! _tm is added 20200728 Takuto Minami

! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

  implicit none

! include "constants.h"   ! commented out 20200728 Takuto Minami
 logical :: suppress_utm_projection=.false. ! inserted on Dec. 17,2015
                        !refer page 22 of specfem3d-manual.pdf

! flag for projection from latitude/longitude to UTM, and back ! 20200728 Takuto Minami from
  integer, parameter :: ILONGLAT2UTM = 0, IUTM2LONGLAT = 1     ! 20200728 Takuto Minami
  double precision, parameter :: PI = 3.141592653589793d0      ! 20200728 Takuto Minami

!
!-----CAMx v2.03
!
!     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!
!     This is a Fortran version of the BASIC program "Transverse Mercator
!     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!     Based on algorithm taken from "Map Projections Used by the USGS"
!     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!
!     Input/Output arguments:
!
!        rlon                  Longitude (deg, negative for West)
!        rlat                  Latitude (deg)
!        rx                    UTM easting (m)
!        ry                    UTM northing (m)
!        UTM_PROJECTION_ZONE  UTM zone
!        iway                  Conversion type
!                              ILONGLAT2UTM = geodetic to UTM
!                              IUTM2LONGLAT = UTM to geodetic
!

  integer utm_projection_zone,iway
  double precision rx,ry,rlon,rlat

  double precision, parameter :: degrad=pi/180., raddeg=180./pi
  double precision, parameter :: semimaj=6378206.4d0, semimin=6356583.8d0
  double precision, parameter :: scfa=.9996d0
  double precision, parameter :: north=0.d0, east=500000.d0

  double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
  double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision rx_save,ry_save,rlon_save,rlat_save

  if(suppress_utm_projection) then
    if (iway == ilonglat2utm) then
      rx = rlon
      ry = rlat
    else
      rlon = rx
      rlat = ry
    endif
    return
  endif

! save original parameters
  rlon_save = rlon
  rlat_save = rlat
  rx_save = rx
  ry_save = ry

! define parameters of reference ellipsoid
  e2=1.0-(semimin/semimaj)**2.0
  e4=e2*e2
  e6=e2*e4
  ep2=e2/(1.-e2)

  if (iway == iutm2longlat) then
    xx = rx
    yy = ry
  else
    dlon = rlon
    dlat = rlat
  endif
!
!----- Set Zone parameters
!
  zone = dble(utm_projection_zone)
  cm = zone*6.0 - 183.0
  cmr = cm*degrad
!
!---- Lat/Lon to UTM conversion
!
  if (iway == ilonglat2utm) then

  rlon = degrad*dlon
  rlat = degrad*dlat

  delam = dlon - cm
  if (delam < -180.) delam = delam + 360.
  if (delam > 180.) delam = delam - 360.
  delam = delam*degrad

  f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat
  f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.
  f2 = f2*sin(2.*rlat)
  f3 = 15.*e4/256.*45.*e6/1024.
  f3 = f3*sin(4.*rlat)
  f4 = 35.*e6/3072.
  f4 = f4*sin(6.*rlat)
  rm = semimaj*(f1 - f2 + f3 - f4)
  if (dlat == 90. .or. dlat == -90.) then
    xx = 0.
    yy = scfa*rm
  else
    rn = semimaj/sqrt(1. - e2*sin(rlat)**2)
    t = tan(rlat)**2
    c = ep2*cos(rlat)**2
    a = cos(rlat)*delam

    f1 = (1. - t + c)*a**3/6.
    f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2
    f2 = f2*a**5/120.
    xx = scfa*rn*(a + f1 + f2)
    f1 = a**2/2.
    f2 = 5. - t + 9.*c + 4.*c**2
    f2 = f2*a**4/24.
    f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2
    f3 = f3*a**6/720.
    yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
  endif
  xx = xx + east
  yy = yy + north

!
!---- UTM to Lat/Lon conversion
!
  else

  xx = xx - east
  yy = yy - north
  e1 = sqrt(1. - e2)
  e1 = (1. - e1)/(1. + e1)
  rm = yy/scfa
  u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256.
  u = rm/(semimaj*u)

  f1 = 3.*e1/2. - 27.*e1**3./32.
  f1 = f1*sin(2.*u)
  f2 = 21.*e1**2/16. - 55.*e1**4/32.
  f2 = f2*sin(4.*u)
  f3 = 151.*e1**3./96.
  f3 = f3*sin(6.*u)
  rlat1 = u + f1 + f2 + f3
  dlat1 = rlat1*raddeg
  if (dlat1 >= 90. .or. dlat1 <= -90.) then
    dlat1 = dmin1(dlat1,dble(90.) )
    dlat1 = dmax1(dlat1,dble(-90.) )
    dlon = cm
  else
    c1 = ep2*cos(rlat1)**2.
    t1 = tan(rlat1)**2.
    f1 = 1. - e2*sin(rlat1)**2.
    rn1 = semimaj/sqrt(f1)
    r1 = semimaj*(1. - e2)/sqrt(f1**3)
    d = xx/(rn1*scfa)

    f1 = rn1*tan(rlat1)/r1
    f2 = d**2/2.
    f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2
    f3 = f3*d**2*d**2/24.
    f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
    f4 = f4*(d**2)**3./720.
    rlat = rlat1 - f1*(f2 - f3 + f4)
    dlat = rlat*raddeg

    f1 = 1. + 2.*t1 + c1
    f1 = f1*d**2*d/6.
    f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2.
    f2 = f2*(d**2)**2*d/120.
    rlon = cmr + (d - f1 + f2)/cos(rlat1)
    dlon = rlon*raddeg
    if (dlon < -180.) dlon = dlon + 360.
    if (dlon > 180.) dlon = dlon - 360.
  endif
  endif

  if (iway == iutm2longlat) then
    rlon = dlon
    rlat = dlat
    rx = rx_save
    ry = ry_save
  else
    rx = xx
    ry = yy
    rlon = rlon_save
    rlat = rlat_save
  endif

  end subroutine utm_geo_tm

!############################################ whoelread ! 2020.09.28
subroutine readcontrolfile(filename,n,lines,ikeep) ! 2021.10.04 ikeep is added
implicit none
character(100),intent(in)  :: filename
integer(4),    intent(in)  :: n
character(200),intent(out) :: lines(n)
integer(4) :: i
integer(4),optional,intent(in) ::  ikeep ! 1 means keep 20 columns 2021.10.04
i=0
write(*,'(a,a)') " file is ",trim(filename)
lines(:)=" " ! initialize 2020.12.10

open(1,file=filename)
do
i=i+1
read(1,'(a)',end=99) lines(i)
end do
99 continue

close(1)

!do i=1,20
!write(*,'(a)') trim(lines(i))
!end do

call subtractcommentout(n,lines)

open(1,file="tmp.ctl")
do i=1,n
if ( present(ikeep) .and. ikeep .eq. 1 ) then !2021.10.04
 write(1,'(a)')  lines(i)(1:len_trim(lines(i))) ! cut initial 20 characters 2021.10.04
else   ! 2021.10.04
 write(1,'(a)')  lines(i)(21:len_trim(lines(i))) ! cut initial 20 characters 2021.09.02
end if ! 2021.10.04
! write(1,'(a)') trim(lines(i))
end do
close(1)

write(*,'(a)') " ### readcontrolfile END!! ###"

return
end subroutine
!############################################# subtractcommentout ! 2020.09.28
subroutine subtractcommentout(n,lines)
implicit none
integer(4),    intent(in)    :: n
character(200),intent(inout) :: lines(n)
character(200)               :: lines_out(n)
integer(4) :: i,j
j=0
do i=1,n
if ( lines(i)(1:2) .ne. "##" ) then
  j=j+1
  lines_out(j) = lines(i)
!  write(*,'(a)')trim(lines_out(j))
 end if
end do

lines = lines_out

return
end subroutine

!######################################################### UTMGMT_N2 2021.09.29
! moved from meshgen2_bell.f90 on 2021.09.29 Name is changed from UTMGMT_N to UTMGMT_N2
subroutine UTMGMT_N2(n,xin,yin,xout,yout,zone,iflag)
implicit none
integer(4),  intent(in)  :: n
real(8),     intent(in)  :: xin(n),yin(n)
real(8),     intent(out) :: xout(n),yout(n)
character(3),intent(in)  :: zone
integer(4),intent(in)    :: iflag ! 0: LONLAT2UTM, 1:UTM2LONLAT
integer(4)               :: i
!#[1]## prepare input file
open(11,file="in.dat")
 write(11,'(2g18.10)') (xin(i),yin(i),i=1,n)
close(11)

!#[2]## use mapproject
if (iflag .eq. 0 ) then ! LONLAT2UTM
 CALL system("cat in.dat | mapproject -Ju"//zone(1:3)//"/1.0 -F -C > tmp.dat")
else if (iflag .eq. 1 ) then ! UTM2LONLAT
 CALL system("cat in.dat | mapproject -Ju"//zone(1:3)//"/1.0 -F -C -I > tmp.dat")
else
 write(*,*) "GEGEGE! iflag should be 0 (LONLAT2UTM) or 1(UTM2LONLAT), iflag=",iflag
 stop
end if

!#[3]## read xout and yout
 open(12,file="tmp.dat")
  do i=1,n
   read(12,*) xout(i),yout(i)
  end do
 close(12)
 ! call system("rm tmp.dat")
return
end subroutine

end module param
