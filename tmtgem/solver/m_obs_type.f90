! coded on May 20, 2016
!
module obs_type
use matrix
implicit none

type obs_info
 integer(4) :: nobs
 character(4) :: name
 real(8),pointer,dimension(:,:)  :: xyz_obs
 real(8),pointer,dimension(:,:)  :: xyzspherical ! 2016.11.20
 real(8),pointer,dimension(:,:)  :: lonlatalt    ! 2016.11.20
 type (real_crs_matrix)          :: coeff(2,3)
 type (real_crs_matrix)          :: coeff_vF
 integer(4),pointer,dimension(:) :: devnum_b(:)  ! 2021.07.26
 integer(4),pointer,dimension(:) :: devnum_e(:)  ! 2021.07.26
end type

contains

!################################################## OPENOBSFILE
! output for electric field is added on 2021.07.26
subroutine OPENOBSFILE(obs_sites,g_param)
use param
implicit none
type(param_forward),intent(in)    :: g_param
type(obs_info),     intent(inout) :: obs_sites
character(70)                     :: filename1,filename2 ! 2021.07.26
character(70)                     :: obsname
integer(4)                        :: inum_b,inum_e ! 2021.07.26
integer(4)                        :: i       ! 2021.07.26
character(70)                     :: outbxyz,outexyz ! 2021.07.27
integer(4)                        :: lenbxyz,lenexyz ! 2021.07.27

allocate(obs_sites%devnum_b(obs_sites%nobs)) ! 2021.07.26
allocate(obs_sites%devnum_e(obs_sites%nobs)) ! 2021.07.26
outbxyz = g_param%outbxyzfolder ! 2018.11.14
outexyz = g_param%outexyzfolder ! 2021.07.27
lenbxyz = len_trim(outbxyz)     ! 2018.11.14
lenexyz = len_trim(outexyz)     ! 2018.11.14

do i=1,obs_sites%nobs
 inum_b=30+i                  ! 2021.07.26
 inum_e=60+i                  ! 2021.07.26
 obs_sites%devnum_b(i)=inum_b ! 2021.07.26
 obs_sites%devnum_e(i)=inum_e ! 2021.07.26
 obsname = g_param%obsname(i)
 !# 2018.11.14 output folder is modified from bxyz/ to outputfolder
 filename1=outbxyz(1:lenbxyz)//obsname(1:len_trim(obsname))//"_bxyz_ts.dat"
 filename2=outexyz(1:lenexyz)//obsname(1:len_trim(obsname))//"_exyz_ts.dat"! 2021.07.26
 open(inum_b,file=filename1) ! 2021.07.26
 open(inum_e,file=filename2) ! 2021.07.26
end do

write(*,*) "### OPENOBSFILE END!! ###"

return
end

!################################################## CLOSEOBSFILE
subroutine CLOSEOBSFILE(obs_sites)
implicit none
type(obs_info),intent(in) :: obs_sites
integer(4) :: i

do i=1,obs_sites%nobs
 close(obs_sites%devnum_b(i)) ! 2021.07.26
 close(obs_sites%devnum_e(i)) ! 2021.07.26
end do

write(*,*) "### CLOSEOBSFILE END!! ###"

return
end

!################################################## BXYZOUT_TS
! Name is changed to BXYZOUT_TS from BEXYZOUT_TS 2021.07.26
subroutine BXYZOUT_TS(obs_sites,bxyz_sites,it,dt)
implicit none
type(obs_info),intent(in) :: obs_sites
integer(4),    intent(in) :: it
real(8),       intent(in) :: dt
real(8),       intent(in) :: bxyz_sites(obs_sites%nobs,3) ! 2021.07.26
integer(4) :: i
real(8) :: t,bxyz_ana(3)

t=it*dt
do i=1, obs_sites%nobs
! call calbxyzana(obs_sites%xyz_obs(1:3,i),it,dt,bxyz_ana)
! write(obs_sites%devnum_b(i),'(7g15.7)') t,bxyz_sites(i,1:3),bxyz_ana(1:3)
 write(obs_sites%devnum_b(i),'(4g15.7)') t,bxyz_sites(i,1:3) !
end do

write(*,*) "### BXYZOUT_TS END!! ###"

return
end
!################################################## EXYZOUT_TS
! Coded on 2021.07.26
subroutine EXYZOUT_TS(obs_sites,exyz_sites,vfxyz_sites,exyz_total_sites,it,dt)
implicit none
type(obs_info),intent(in) :: obs_sites
integer(4),    intent(in) :: it
real(8),       intent(in) :: dt
real(8),       intent(in) :: exyz_sites(obs_sites%nobs,3) ! 2021.07.26
real(8),       intent(in) :: vfxyz_sites(obs_sites%nobs,3) ! 2021.07.26
real(8),       intent(in) :: exyz_total_sites(obs_sites%nobs,3) ! 2021.07.26
integer(4) :: i
real(8) :: t,bxyz_ana(3)

t=it*dt
do i=1, obs_sites%nobs
! call calbxyzana(obs_sites%xyz_obs(1:3,i),it,dt,bxyz_ana)
! write(obs_sites%devnum_b(i),'(7g15.7)') t,bxyz_sites(i,1:3),bxyz_ana(1:3)
 write(obs_sites%devnum_e(i),'(10g15.7)') t,exyz_sites(i,1:3),vfxyz_sites(i,1:3),exyz_total_sites(i,1:3)
end do

write(*,*) "### EXYZOUT_TS END!! ###"

return
end

!##################################################
subroutine calbxyzana(x3,it,dt,bxyz)
implicit none
real(8),intent(in) :: x3(3),dt
real(8),intent(out) :: bxyz(3)
integer(4),intent(in) :: it
real(8) :: v3(3),Fz=-30000.d0 ! [nT]
real(8),parameter :: h=8000.d0, g=9.8d0
real(8),parameter :: sigma=3.3d0, pi=4.d0*atan(1.d0)
real(8),parameter :: mu=4.*pi*1.d-7
complex(8) :: iunit=(0.d0,1.d0)
complex(8) :: CT,eta,bz
real(8) :: t,k,c,ch
t=(it+0.5d0)*dt

c=sqrt(g*h) ! [m/sec]
ch=2./mu/sigma/h
k=2.*pi/80.d0 ! [rad/km]

eta=cdexp(iunit*(k*(x3(1)-(c*1.d-3)*t))) ! [m]
bz=c/(c+iunit*ch)*eta/h*Fz
bxyz(3)=real(bz)
bxyz(1)=real(iunit*bz)
bxyz(2)=0.d0
return
end

end module obs_type
