! Coded on 2017.02.22
module cond_type
implicit none

type cuboid
 real(8) :: xminmax(2) ! [km]
 real(8) :: yminmax(2) ! [km]
 real(8) :: zminmax(2) ! [km]
 real(8) :: rho
 real(8) :: ratio      ! rho_n = 10^(ratio*log(rho) + (1-ratio)*log(rho_ori))
end type

type param_cmod
 character(50)  :: mshfile
 character(50)  :: outcondfile
 character(100) :: inmodelfile,outmodelfile,connectfile ! 2021.01.05
 real(8)        :: ratio ! rho_n = rho*ratio
 integer(4)     :: ncuboid
 type(cuboid),allocatable,dimension(:) :: g_cuboid
end type

contains

!################################# readcmodparam
subroutine readcmodparam(m_param)
implicit none
type(param_cmod) :: m_param
integer(4)       :: input=5,i

!#[1]# skip header
 read(input,*)

!#[2]# read
 read(input,10) m_param%mshfile
 read(input,10) m_param%outcondfile

!#[3]# ratio
 read(input,11) m_param%ratio

!#[4]# ratio
 read(input,12) m_param%ncuboid
 allocate(m_param%g_cuboid(m_param%ncuboid))

!#[5]#
 do i=1,m_param%ncuboid
   read(input,13) m_param%g_cuboid(i)%xminmax(1:2)
   read(input,13) m_param%g_cuboid(i)%yminmax(1:2)
   read(input,13) m_param%g_cuboid(i)%zminmax(1:2)
   read(input,11) m_param%g_cuboid(i)%rho
   read(input,11) m_param%g_cuboid(i)%ratio
 end do

10 format(20x,a)
11 format(20x,g15.7)
12 format(20x,i10)
13 format(20x,2g15.7)
return
end subroutine

end module cond_type
