module slicesub ! 2020.12.04
implicit none
contains
!===================================================================
!# 2018.03.12 ifileflag = 2 for d_cond.msh is added
!# When A=0 and B=0, C=1 is set and horizontal slice can be generated
subroutine OUTSLICE(g_paramslice,g_mesh,rho_total,ele2model,nmodel,nmodelfile,nphys2,ifileflag,iexist)
use mesh_type
use post_param
use outerinnerproduct
implicit none
type(mesh),         intent(in) :: g_mesh
type(param_slice),  intent(in) :: g_paramslice
integer(4),         intent(in) :: nmodel,nphys2,nmodelfile
integer(4),         intent(in) :: ele2model(nphys2)
real(8),            intent(in) :: rho_total(nmodel,nmodelfile)
integer(4),         intent(in) :: ifileflag ! 2018.03.12
integer(4),optional,intent(in) :: iexist(nmodelfile)! 2018.01.17
real(8),   allocatable,dimension(:,:) :: xyz
integer(4),allocatable,dimension(:,:) :: n4
integer(4)    :: node,ntet,nphys1
real(8)       :: bl,br,bb,bt,A,B,C,D
integer(4),allocatable,dimension(:)  :: ind
logical,   allocatable,dimension(:)  :: plus,hout1,hout2
integer(4)    :: n(4),icount_plus,icount_minus,ind_plus(4),ind_minus(4)
integer(4)    :: ii,i,j,k
real(8)       :: xh,h1,h2,h3,h4
real(8),allocatable,dimension(:) :: SC ! 2018.04.17
real(8)       :: x1(3),x2(3),xt(3),x(3)
integer(4)    :: isig,input
integer(4)    :: imodel,nslice,islice ! 2017.12.21
integer(4)    :: ioutdiv(nmodelfile)
character(70) :: filename,head    ! 2017.12.21
character(2)  :: num ! 2017.12.21
integer(4)    :: ihorflag ! 0:verticalslice, 1:horizontal slice 2018.01.19
real(8),dimension(3) :: xn,zn,xdir,ydir,xo
real(8) :: CC(4,4)

 !#[0]## set
 node = g_mesh%node
 ntet = g_mesh%ntet
 nphys1 = ntet - nphys2
 write(*,*) "ntet",ntet
 write(*,*) "nphys1",nphys1,"nphys2",nphys2
 allocate(xyz(3,node),n4(ntet,4))
 xyz  = g_mesh%xyz
 n4   = g_mesh%n4
 nslice = g_paramslice%nslice
 allocate( plus(node),hout1(node),hout2(node))
 allocate( SC(nmodelfile) ) ! 2018.04.17

 !#[1]## slice loop start
 do islice=1,nslice          ! 2017.12.21 slice loop start

 !#[3]## open polygonfiles for each model
 do imodel=1,nmodelfile ! 2017.12.21
  if( present(iexist) ) then ! 2018.01.18
   if ( iexist(imodel) .eq. 0 ) goto 80 ! 2018.01.17
  end if                     ! 2018.01.18
  ioutdiv(imodel) = 10 + imodel
  head=g_paramslice%outpolygonhead(imodel)
  write(num,'(i2.2)') islice
  filename=head(1:len_trim(head))//"_"//num(1:2)//".dat"
  open(ioutdiv(imodel),file=filename)
  80 continue                ! 2018.01.18
 end do

 !# Ax + By + D = 0
 A   = g_paramslice%A(islice) ! 2017.12.21
 B   = g_paramslice%B(islice) ! 2017.12.21
 D   = g_paramslice%D(islice) ! 2017.12.21
 bl  = g_paramslice%bl(islice) ! 2018.01.19
 br  = g_paramslice%br(islice) ! 2018.01.19
 bb  = g_paramslice%bb(islice) ! 2018.01.19
 bt  = g_paramslice%bt(islice) ! 2018.01.19

 !# Horisonal slice can be made when A=0 and B=0 ! 2018.01.19
 C=0.d0 ! 2018.01.19
 ihorflag = 0
 if ( A .lt. 1.d-5 .and. B .lt. 1.d-5 ) then
  C=1.d0 ! 2018.01.19
  ihorflag = 1 ! horizontal slice
 end if

 !#[2]## separate nodes to two region
  hout1(1:node) = .false.
  hout2(1:node) = .false.
  plus(1:node)  = .false.
  zn = (/ 0.d0,0.d0,1.d0/)
  xn = (/A,B,C/)/sqrt(A**2. + B**2. + C**2.)
  if ( A .lt. 1.d-5  .and. B .gt. 1.d-5 ) xn = (/0.d0,-1.d0,0.d0/)
  if ( ihorflag .eq. 1) zn = (/0.d0,1.d0,0.d0/)
!  write(*,*) "zn",zn
!  write(*,*) "xn",xn
  xdir = outer(zn,xn)
  if ( B .ne. 0. .and.  -A/B .lt. 0. ) xdir = outer(xn,zn) ! 2018.03.28
!  write(*,*) "xdir",xdir
  xdir(1:3) = xdir(1:3) / sqrt(xdir(1)**2. + xdir(2)**2. + xdir(3)**2.d0)
  ydir = outer(xn,xdir)
  if ( B .ne. 0. .and.  -A/B .lt. 0. ) ydir = outer(xdir,xn) ! 2018.03.28
!  write(*,*) "xdir",xdir
!  write(*,*) "ydir",ydir
  !# xo : vector from origin to the origin on the plane for drawing
  xo(1:3) = -D*(/A,B,C/)/(A**2. + B**2. + C**2.)

  CC(1,1:4)=(/xdir(1:3),-bl/) ! xmin
  CC(2,1:4)=(/xdir(1:3),-br/) ! xmax
  CC(3,1:4)=(/ydir(1:3),-bb/) ! ymin
  CC(4,1:4)=(/ydir(1:3),-bt/) ! ymax

 do i=1,node
  x(1:3)=xyz(1:3,i)
  if ( A*x(1)+B*x(2)+C*x(3)+D .gt. 0.d0 ) plus(i)=.true.
  h1=inner(CC(1,1:3),x) + CC(1,4)  ! Ax + By + Cz + D
  h2=inner(CC(2,1:3),x) + CC(2,4)
  if ( h1 .lt. 0.d0 .or. 0.d0 .lt. h2 ) hout1(i)=.true.
  h3=inner(CC(3,1:3),x) + CC(3,4)
  h4=inner(CC(4,1:3),x) + CC(4,4)
  if ( h3 .lt. 0.d0 .or. 0.d0 .lt. h4 ) hout2(i)=.true.
 end do
 write(*,*) "### categolize nodes end !! ###"

 !#[3]##
 do i=1,nphys2
  if (mod(i,100000) .eq. 0) write(*,*) "i",i,"nphys2",nphys2
  n(1:4) = n4(i+nphys1,1:4) ! node id for whole node space
  if ( hout1(n(1)) .and. hout1(n(2)) .and. hout1(n(3)) .and. hout1(n(4)) ) cycle
  if ( hout2(n(1)) .and. hout2(n(2)) .and. hout2(n(3)) .and. hout2(n(4)) ) cycle
  icount_plus  = 0
  icount_minus = 0
  ind_minus = 0
  ind_plus  = 0
  do j=1,4
   if ( plus(n(j)) ) then
    icount_plus = icount_plus + 1
    ind_plus(icount_plus) = j
   else
    icount_minus = icount_minus + 1
    ind_minus(icount_minus) = j
   end if
  end do
  if ( icount_plus .eq. 0 .or. icount_plus .eq. 4 ) cycle
  !# color
  do imodel=1,nmodelfile ! treat with many model files simultaneously
   if( present(iexist) ) then            ! 2018.01.17
    if ( iexist(imodel) .eq. 0 ) goto 70 ! 2018.01.17
   end if                                ! 2018.01.17
   if (      ifileflag .eq. 0 .or. ifileflag .eq. 1 ) then ! cond or model
    if ( rho_total(ele2model(i),imodel) .lt. 1.d-10 ) then ! 2022.07.05
      SC(imodel) = -10.0 ! 2022.07.05
    else ! 2022.07.05
      SC(imodel) = log10(rho_total(ele2model(i),imodel))
    end if ! 2022.07.05
   else if ( ifileflag .eq. 2 .or. ifileflag .eq. 3 ) then ! d_cond or d_model 2018.03.19
    SC(imodel) = rho_total(ele2model(i),imodel)
   else
    write(*,*) "GEGEGE"
    stop
   end if
   if ( g_paramslice%ibound .eq. 1 .and. SC(imodel) .lt. g_paramslice%lbound ) cycle ! 2018.04.17
   if ( g_paramslice%ibound .eq. 2 .and. SC(imodel) .gt. g_paramslice%ubound ) cycle ! 2018.04.17
   write(ioutdiv(imodel),'(a4,f0.5)')"> -Z",SC(imodel)
   70 continue ! 2018.01.17
  end do
  !# when 2 plus and 2 minus
  if ( icount_plus .eq. 2 ) then
    ii=ind_plus(1)
    x1(1:3) = xyz(1:3,n(ii))
    do k=1,2
     x2(1:3) = xyz(1:3,n(ind_minus(k)))
     if ( present(iexist) ) then
      call pointwrite(x1,x2,xdir,ydir,xo,A,B,C,D,ioutdiv,nmodelfile,SC,g_paramslice,iexist)
     else
      call pointwrite(x1,x2,xdir,ydir,xo,A,B,C,D,ioutdiv,nmodelfile,SC,g_paramslice)
     end if
    end do
    ii=ind_plus(2)
    x1(1:3) = xyz(1:3,n(ii))
    do k=2,1,-1
     x2(1:3) = xyz(1:3,n(ind_minus(k)))
     if ( present(iexist) ) then ! 2018.01.17
      call pointwrite(x1,x2,xdir,ydir,xo,A,B,C,D,ioutdiv,nmodelfile,SC,g_paramslice,iexist)
     else
      call pointwrite(x1,x2,xdir,ydir,xo,A,B,C,D,ioutdiv,nmodelfile,SC,g_paramslice)
     end if
    end do
  else
   !# only 1 or 3 plus
   if ( icount_plus .eq. 1 ) ii=ind_plus(1)
   if ( icount_plus .eq. 3 ) ii=ind_minus(1)
   do j=1,4
     if ( j .eq. ii ) cycle
     x1(1:3)  = xyz(1:3,n(ii))
     x2(1:3)  = xyz(1:3,n(j))
     if ( present(iexist) ) then
      call pointwrite(x1,x2,xdir,ydir,xo,A,B,C,D,ioutdiv,nmodelfile,SC,g_paramslice,iexist)
     else
      call pointwrite(x1,x2,xdir,ydir,xo,A,B,C,D,ioutdiv,nmodelfile,SC,g_paramslice)
     end if
   end do
  end if
 end do
 100 continue ! 2018.01.19

 !#[4]##
 do imodel=1,nmodelfile ! 2017.12.21
  if( present(iexist) ) then ! 2018.01.17
    if ( iexist(imodel) .eq. 0 ) goto 60 ! 2018.01.17
   end if                     ! 2018.01.17
  close(ioutdiv(imodel))! 2017.12.21
  60 continue
 end do

 end do ! 2017.12.21 slice loop

return
end subroutine OUTSLICE ! 2018.01.19

!#########################################
!# modified on 2018.01.20
!# xo   : vector representing the origin at the plane
!# xdir : unit vector in xaxis in the plane
!# ydir : unit vector in yaxis in the plane
subroutine pointwrite(x1,x2,xdir,ydir,xo,A,B,C,D,ioutdiv,nmodelfile,SC,g_paramslice,iexist) ! 2018.01.17
use outerinnerproduct
use post_param ! 2018.04.17
implicit none
real(8),            intent(in) :: x1(3),x2(3),xo(3),xdir(3),ydir(3)
integer(4),         intent(in) :: nmodelfile          ! 2017.12.21
integer(4),         intent(in) :: ioutdiv(nmodelfile) ! 2017.12.21
real(8),            intent(in) :: A,B,C,D
integer(4),optional,intent(in) :: iexist(nmodelfile)  ! 2018.01.17
type(param_slice),  intent(in) :: g_paramslice        ! 2018.04.17
real(8),            intent(in) :: SC(nmodelfile)      ! 2018.04.17
real(8)    :: x12(3),s,xt(3),xh,x,y
integer(4) :: imodel
real(8)    :: xa(3)

!#[1]## calculate the intersecting point
x12 = x2 - x1
xa(1:3)=(/A,B,C/)
s=(-D-inner(xa,x1))/inner(xa,x12)     ! 2018.01.19
if ( s .lt. 0.d0 .or. 1.d0 .lt. s ) goto 100
xt(1:3)  = x1(:) + s*x12(:)

!#[2]## coordinate in the plane of Ax + By + Cz + D = 0
x=inner(xt-xo,xdir)
y=inner(xt-xo,ydir)

do imodel=1,nmodelfile ! 2017.12.21
 if ( present(iexist) ) then           ! 2018.01.17
  if ( iexist(imodel) .eq. 0) goto 60  ! 2018.01.17
 end if
 if ( g_paramslice%ibound .eq. 1 .and. SC(imodel) .lt. g_paramslice%lbound ) cycle ! 2018.04.17
 if ( g_paramslice%ibound .eq. 2 .and. SC(imodel) .gt. g_paramslice%ubound ) cycle ! 2018.04.17
 write(ioutdiv(imodel),'(2g15.7)') x,y ! 2018.01.20
 60 continue                           ! 2018.01.17
end do ! 2017.12.21
!write(*,    '(2g15.7)') xh,xt(3)

return

100 continue
write(*,*) "GEGEGE s should be [0:1] s=",s
stop
end subroutine

end module slicesub ! 2020.12.04

