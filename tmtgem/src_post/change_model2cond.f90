! Takuto Minami
! coded on 2022.10.24
program change_model2cond
use param
use mesh_type
use cond_model_type ! m_cond_model_type.f90
implicit none
type(mesh)            :: h_mesh
type(param_forward)   :: g_param
type(param_cond)      :: g_cond
type(param_cond)      :: h_cond ! see ../common/m_param.f90
type(info_cond_model) :: g_info_cond_model
integer(4) :: i

!#[1]##
 call readctl(g_info_cond_model) ! see below

!#[2]## read 3d msh tetrahedral mesh
 CALL READMESH_TOTAL(h_mesh,g_info_cond_model%mshfile) ! read 3dmesh
 CALL GENXYZMINMAX(h_mesh,g_param)

!#[3]## read homo, cond, or model file to generate g_cond
 CALL READGENCOND(g_info_cond_model,h_mesh,g_cond) ! see m_cond_model_type.f90

write(*,*) "nphys2",g_cond%nphys2
do i=1,g_cond%nphys2
!write(*,*) "i",i,"rho",g_cond%rho(i)
end do
!#[4]## change the model
 CALL ASSIGNRHO(h_mesh,g_cond,h_cond,g_info_cond_model) ! gen h_cond

!#[5]# output new cond (h_cond)
 h_cond%condfile = g_info_cond_model%outcondfile
 CALL OUTCOND(h_cond,h_mesh) ! 2022.10.22

end program

!##################################################################################
subroutine readctl(g_info_cond_model) ! see m_cond_model_type.f90
use cond_model_type
implicit none
type(info_cond_model),intent(out) :: g_info_cond_model ! m_cond_type.f90
integer(4)       :: input=5,i,iflag_cond,iflag_ele

!#[1]# skip header
 read(input,*)

!#[2]# read conductivity strucuture info
 read(input,10) g_info_cond_model%mshfile
 read(input,12) iflag_cond
 write(*,*) "iflag",iflag_cond
 g_info_cond_model%iflag_cond = iflag_cond

 if ( iflag_cond == 0 ) then !read homogeneous info
   read(input,11) g_info_cond_model%rhohomo
   write(*,*) "rho homo [Ohm.m]",g_info_cond_model%rhohomo

 else if ( iflag_cond == 1 ) then !read cond info
   read(input,10) g_info_cond_model%inputcondfile

 else if ( iflag_cond == 2 ) then ! read model info
   read(input,10) g_info_cond_model%connectfile
   write(*,*) "connectfile ",trim(g_info_cond_model%connectfile)
   read(input,10) g_info_cond_model%modelfile
   write(*,*) "model file ",trim(g_info_cond_model%modelfile)

 end if

!#[3]##
 read(input,10) g_info_cond_model%outcondfile
    write(*,*) "outcond file ",trim(g_info_cond_model%outcondfile)

!#[3]## elevation of depth from surf: 0: ele, 1: depth from the surface
 read(input,12) iflag_ele
 write(*,*) "iflag_ele",iflag_ele
 g_info_cond_model%iflag_ele = iflag_ele

if (iflag_ele .eq. 1 ) read(input,10) g_info_cond_model%topo_2d_mshfile

!#[4]# ratio
 read(input,12) g_info_cond_model%ncuboid
 allocate(g_info_cond_model%g_cuboid(g_info_cond_model%ncuboid))

!#[5]#
 do i=1,g_info_cond_model%ncuboid
   read(input,13) g_info_cond_model%g_cuboid(i)%xminmax(1:2)
   read(input,13) g_info_cond_model%g_cuboid(i)%yminmax(1:2)
   read(input,13) g_info_cond_model%g_cuboid(i)%zminmax(1:2)
   read(input,11) g_info_cond_model%g_cuboid(i)%rho
   write(*,*) "cuboid",i
   write(*,*) "xminmax",g_info_cond_model%g_cuboid(i)%xminmax(1:2)
   write(*,*) "yminmax",g_info_cond_model%g_cuboid(i)%yminmax(1:2)
   write(*,*) "zminmax",g_info_cond_model%g_cuboid(i)%zminmax(1:2)
   write(*,*) "rho",g_info_cond_model%g_cuboid(i)%rho,"[Ohm.m]"
 end do

10 format(20x,a)
11 format(20x,g15.7)
12 format(20x,i10)
13 format(20x,2g15.7)
return
end

!###################################################################################  change cond
!# 2022.10.24
subroutine ASSIGNRHO(h_mesh,g_cond,h_cond,g_info_cond_model)
use mesh_type
use cond_model_type
use param
use triangle
implicit none
type(info_cond_model),intent(in)      :: g_info_cond_model
type(mesh),           intent(in)      :: h_mesh
type(param_cond),     intent(in)      :: g_cond
type(param_cond),     intent(out)     :: h_cond
type(mesh)                            :: topo2d_mesh
integer(4),allocatable,dimension(:,:) :: n4,n3k
integer(4),allocatable,dimension(:)   :: n4flag
real(8),   allocatable,dimension(:,:) :: xyz,cxyz,xyzk
type(grid_list_type)                  :: glist         ! see ../common/m_mesh_type.f90
integer(4)                            :: node,ntet,ntri,nphys2
integer(4)                            :: i,j,l,nphys1,n1,n2,n3,nx,ny,iele
integer(4)                            :: iflag_ele
real(8)                               :: x2(2),y2(2),z2(2),z_topo,a3(3),xyzminmax(6)
!
!#[1]##
ntet       = h_mesh%ntet
node       = h_mesh%node
allocate(n4flag(ntet), cxyz(3,ntet))
n4         = h_mesh%n4     ! allocate and fill (ntet,4)
n4flag(:)  = h_mesh%n4flag(:,1)
xyz        = h_mesh%xyz    ! allocatae and fill (3,node)
nphys2     = g_cond%nphys2 ! 2022.10.24
nphys1     = g_cond%nphys1
iflag_ele  = g_info_cond_model%iflag_ele
h_cond     = g_cond        ! output default is the same as input
xyzminmax  = h_mesh%xyzminmax
!write(*,*) "xyzminmax",xyzminmax 

!#[2]## assign center of gravity
 do i=1,nphys2
   ![2-1]# center of gravity of this tetrahedron
   l = nphys1 + i
   do j=1,3
     cxyz(j,l)=(xyz(j,n4(l,1))+xyz(j,n4(l,2))+xyz(j,n4(l,3))+xyz(j,n4(l,4)))/4.d0
   end do
 end do
 !write(*,*) "check1"

!#[2]## cal z for nobsr
 if ( iflag_ele == 1) then! depth from surface
   CALL READMESH_TOTAL(topo2d_mesh,g_info_cond_model%topo_2d_mshfile)
   ntri       = topo2d_mesh%ntri
   n3k        = topo2d_mesh%n3
   xyzk       = topo2d_mesh%xyz
   write(*,*) "ntri",ntri
   nx=300;ny=300
   CALL allocate_2Dgrid_list(nx,ny,ntri,glist)  ! see m_mesh_type.f90
   CALL gen2Dgridforlist(xyzminmax,glist)       ! see m_mesh_type.f90
   CALL classifytri2grd(topo2d_mesh,glist)      ! classify ele to glist,see
   write(*,*) "classify2grd end!!"
 end if

!#[3]## change model

 do i=1,g_info_cond_model%ncuboid
   x2=g_info_cond_model%g_cuboid(i)%xminmax
   y2=g_info_cond_model%g_cuboid(i)%yminmax
   z2=g_info_cond_model%g_cuboid(i)%zminmax
   
   do j=1,nphys2 ! tetrahedron loop
     l=nphys1 + j

     if (x2(1) .lt. cxyz(1,l) .and. cxyz(1,l) .lt. x2(2) ) then
       if (y2(1) .lt. cxyz(2,l) .and. cxyz(2,l) .lt. y2(2) ) then
          if ( iflag_ele == 0 ) then ! normal cuboid
            if (z2(1) .lt. cxyz(3,l) .and. cxyz(3,l) .lt. z2(2) ) then
!              write(*,*) "In cxyz",cxyz(1:3,l)
              h_cond%rho(j) = g_info_cond_model%g_cuboid(i)%rho
!              write(*,*) "rho is change from",g_cond%rho(j)," ->",h_cond%rho(j)
            end if
          else if ( iflag_ele == 1) then ! depth from surface
!             write(*,*) "cxyz(1:2,l)",cxyz(1:2,l)
             call findtriwithgrid(topo2d_mesh,glist,cxyz(1:2,l),iele,a3) ! 2017.07.14
!             write(*,*) "iele=",iele
             n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3) ! 3 node for iele triangle
             z_topo = a3(1)*xyzk(3,n1)+a3(2)*xyzk(3,n2)+a3(3)*xyzk(3,n3)
             if (z_topo + z2(1) .lt. cxyz(3,l) .and. cxyz(3,l) .lt. z_topo + z2(2) ) then
               h_cond%rho(j) = g_info_cond_model%g_cuboid(i)%rho
             end if
          end if
       end if
      end if

   end do ! nphys2 loop

  end do ! cuboid loop

  do j=1,nphys2
   h_cond%sigma(i)=1.d0/h_cond%rho(j)
  end do

write(*,*) "### ASSIGNRHO END!! ###"
return
end
!###################################################################
! copied from n_ebfem_bxyz.f90 on 2017.05.10
subroutine GENXYZMINMAX(em_mesh,g_param)
use param ! 2016.11.20
use mesh_type
implicit none
type(mesh),            intent(inout) :: em_mesh
type(param_forward),   intent(inout) :: g_param
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
real(8) :: xyz(3,em_mesh%node),xyzminmax(6)
integer(4) :: i
xyz = em_mesh%xyz ! normal
xmin=xyz(1,1) ; xmax=xyz(1,1)
ymin=xyz(2,1) ; ymax=xyz(2,1)
zmin=xyz(3,1) ; zmax=xyz(3,1)

do i=1,em_mesh%node
 xmin=min(xmin,xyz(1,i))
 xmax=max(xmax,xyz(1,i))
 ymin=min(ymin,xyz(2,i))
 ymax=max(ymax,xyz(2,i))
 zmin=min(zmin,xyz(3,i))
 zmax=max(zmax,xyz(3,i))
end do

xyzminmax(1:6)=(/xmin,xmax,ymin,ymax,zmin,zmax/)

!# set output
g_param%xyzminmax = xyzminmax
em_mesh%xyzminmax = xyzminmax    ! 2021.12.27

write(*,*) "### GENXYZMINMAX END!! ###"
return
end
!########################################### OUTCOND
! Output folder is changed on 2017.09.04
! Coded on 2017.05.18
subroutine OUTCOND(g_cond,g_mesh) ! 2017.09.11
use mesh_type
use param
implicit none
type(param_cond), intent(in) :: g_cond
type(mesh),       intent(in) :: g_mesh
integer(4)    :: j, nphys2,npoi,nlin,ntri,ishift,nphys1
real(8),   allocatable,dimension(:)    :: rho
integer(4),allocatable,dimension(:) :: index
character(100) :: outcondfile

!#[0]## set
outcondfile  = g_cond%condfile
nphys1       = g_cond%nphys1
nphys2       = g_cond%nphys2
allocate(index(nphys2),rho(nphys2))
index        = g_cond%index
rho          = g_cond%rho
npoi         = g_mesh%npoi
nlin         = g_mesh%nlin
ntri         = g_mesh%ntri
write(*,*) "nphys1=",nphys1
write(*,*) "nphys2=",nphys2

!#[1]## output rho
 open(1,file=outcondfile)

 !# standard info
! CALL MESHOUT(1,g_mesh)

 write(1,'(a)') "$MeshFormat"     ! 2017.09.13
 write(1,'(a)') "2.2 0 8"        ! 2017.09.13
 write(1,'(a)') "$EndMeshFormat"  ! 2017.09.13
 write(1,'(a)') "$ElementData"
 write(1,'(a)') "1"
 write(1,'(a)') '"A rho model view"'
 write(1,'(a)') "1"
 write(1,'(a)') "0.0"
 write(1,'(a)') "3"
 write(1,'(a)') "0"
 write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
 write(1,'(i10)') nphys2
 ishift = npoi + nlin + ntri
 do j=1,nphys2
!  write(*,*) "j=",j,"nphys2=",nphys2,"ele2model(j)=",ele2model(j),"nmodel=",nmodel
  write(1,*) ishift+index(j),rho(j)
 end do
 write(1,'(a)') "$EndElementData"
close(1)

write(*,*) "### OUTCOND END!! ###"
return
end
