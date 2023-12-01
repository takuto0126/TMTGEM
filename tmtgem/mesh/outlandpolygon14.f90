!########################################################### outlandpolygon14
!# coded on 2019.03.07
!# land polygon for GMT plotting
subroutine outlandpolygon14(h_poly,ilflg,xcorner,ycorner)
use coastline_data
use param_mesh ! see m_param_mesh.f90, added on 2017.06.28
implicit none
type(poly_data),       intent(in)     :: h_poly ! after nclose polygons are combined
real(8),               intent(in)     :: xcorner(0:3),ycorner(0:3)  ! 2019.02.20 0 for right top
integer(4),            intent(in)     :: ilflg(0:3)
character(70)                         :: polygonfile
integer(4)                            :: ncmax,lpmax,nclose
logical                               :: iflag_topright_land ! 2021.05.29
integer(4)                            :: lpoly0,i,j
integer(4),allocatable,dimension(:)   :: lpoly
real(8),   allocatable,dimension(:,:) :: xpoly,ypoly
real(8),   allocatable,dimension(:,:) :: loc

!#[0]## set input
 ncmax        = h_poly%ncmax                   ! 2018.08.28
 lpoly0       = h_poly%lpoly0
 nclose       = h_poly%nclose
 lpmax        = h_poly%lpmax
 allocate(      xpoly(lpmax,ncmax)          )  ! 2018.08.28
 allocate(      ypoly(lpmax,ncmax)          )  ! 2018.08.28
 allocate(      lpoly(lpmax)                 )
 lpoly        = h_poly%lpoly                    ! 2018.08.28
 xpoly        = h_poly%xypoly(1,:,:)   ! east
 ypoly        = h_poly%xypoly(2,:,:)   ! north
 iflag_topright_land = h_poly%iflag_topright_land    ! 2019.02.21
 allocate( loc(nclose,2) )    ! 2019.02.21
 loc          = h_poly%loc    ! 2019.02.21

!#[0]## starting
 write(*,*) "start outlandpolygon14"
 write(*,'(a,l3)') "iflag_topright_land",iflag_topright_land

!#[1]## set and open file
 polygonfile="landpolygon.gmt"
 write(*,*) "polygonfile",polygonfile
 open(1,file=polygonfile)

 do i=1,lpoly0
   write(1,'(a4)') "> -Z"
   do j=1,lpoly(i)
      write(1,'(2g15.7)') xpoly(i,j),ypoly(i,j)
   end do

  !# add corners
  if ( i .eq. 1 .and. iflag_topright_land ) then
    do j=3,0,-1
      if ( j .lt. loc(1,2) ) write(1,'(2g15.7)') ycorner(j),xcorner(j)
    end do
    do j=3,1,-1
      if ( loc(1,1) .lt. j) write(1,'(2g15.7)') ycorner(j),xcorner(j)
    end do
  else
    do j=3,1,-1
      if ( ilflg(j) .eq. i) write(1,'(2g15.7)') ycorner(j),xcorner(j)
    end do
  end if

 end do
 close(1)

 write(*,*) "### outlandpolygon14 END!! ###"

return
end
