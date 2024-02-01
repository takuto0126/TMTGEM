! Modified on 2020.12.10
! Coded on 2017.12.20
! to calculate slice data for GMT
program slice
use mesh_type
use post_param
use param      ! 2018.01.06
use slicesub   ! 2020.12.04
implicit none
type(param_slice) :: g_paramslice
type(mesh)        :: g_mesh
integer(4)        :: nphys2,nmodel,nmodel2,i,j
integer(4),allocatable,dimension(:) :: id,ele2model
real(8),   allocatable,dimension(:,:) :: rho
integer(4)        :: nmodelfile
integer(4)        :: ncondfile ! 2018.01.06
type(param_cond),allocatable,dimension(:) :: g_cond ! 2018.01.06
integer(4)        :: ifileflag ! 2018.01.06 0: condfile, 1:modelfile
integer(4),allocatable,dimension(:) :: iexist ! 0:no model, 1:success

!#[0]## read param
 call sliceparamread(g_paramslice) ! see m_post_param.f90
 ifileflag = g_paramslice%ifileflag ! 2018.01.06

!#[1]## read mesh
 call READMESH_TOTAL(g_mesh,g_paramslice%mshfile)

!#[2]## modelfile case
 !#[2-1]## read connection
 if ( ifileflag .eq. 1 .or. ifileflag .eq. 3 ) then ! modelfile case 2018.03.19
  open(1,file=g_paramslice%connectfile)
  read(1,*) nphys2,nmodel
  allocate(id(nphys2),ele2model(nphys2))
  do i=1,nphys2
   read(1,'(2i10)') id(i),ele2model(i)
  end do
  close(1)
  write(*,*) "### READ CONNECTIONFILE END!! ###"

 !#[2-2]## read model
 nmodelfile = g_paramslice%nmodelfile
 allocate(rho(nmodel,nmodelfile),iexist(nmodelfile))
 do j=1,nmodelfile    ! 2017.12.21
 open(1,file=g_paramslice%modelfile(j),status='old',err=90) ! 20127.12.21
  read(1,*,err=80) nmodel2
  do i=1,nmodel
   read(1,*) rho(i,j)
  end do
 close(1)
 iexist(j)=1 ! 2018.01.17
 write(*,*) "### READ MODELFILE END!! ###",j
 goto 100
  80 continue
  close(1)
  iexist(j)=0
  write(*,*) "Cannot read",j
  goto 100
  90 continue
  write(*,*) "File is not exist",j,g_paramslice%modelfile(j)
  iexist(j)=0 ! 2018.01.17
 100 continue ! 2018.01.17
 end do ! 2017.12.21
 end if

!#[3]## condfile case 2018.01.06
 if ( ifileflag .eq. 0 .or. ifileflag .eq. 2 ) then ! condfile,dcondfile case 2018.03.12
  nmodelfile = g_paramslice%ncondfile
  allocate( g_cond(nmodelfile) )
  do i=1,nmodelfile
   g_cond(i)%condfile = g_paramslice%condfile(i)
   call READCOND(g_cond(i))        ! 2018.01.06
  end do
  nphys2 = g_cond(1)%nphys2        ! 2018.01.06
  allocate( ele2model(nphys2) )
  do i=1,nphys2                    ! 2018.01.06
   ele2model(i) = i                ! 2018.01.06
  end do                           ! 2018.01.06
  nmodel = nphys2                  ! 2018.01.06
  allocate(rho(nmodel,nmodelfile)) ! 2018.01.06
  do i=1,nmodelfile                ! 2018.01.06
   rho(1:nmodel,i) = g_cond(i)%rho(1:nphys2) ! 2018.01.06
  end do                           ! 2018.01.06
 end if                            ! 2018.01.06

!#[4]## output slice polygon data
 ! 2018.01.17 iexist is added
 call OUTSLICE(g_paramslice,g_mesh,rho,ele2model,nmodel,nmodelfile,nphys2,ifileflag,iexist) !2018.03.12

end program

