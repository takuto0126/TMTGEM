!# Modified to permit land at top right corner on 2019.02.20
!# May 22, 2014, zlabel is introduced to identify zero-depth nodes.

module addcorner11m
use coastline_data
implicit none
contains

!########################################################### addcorner11
! ilflg(i) = j : i-th corner is included by j-th land polygon, j=0 means ocean
!  i=[0:3]     : 0 means top right corner 2019.02.20
!
! iflg(i) = k  : i-th corner node is between kth and k+1 land polygons
!  i=[0:3]     : default = 999
!
! corner numbering  : i=0 -> top right      2019.02.20
!                   : i=1 -> bottom right, i=2 -> bottom left, i=3 -> top left
!
! loc(i,1) : value for the start boundary of i-th unclosed polygon,
! loc(i,2) : value for the end   boundary of i-th unclosed polygon
!  i=[1:nclose]
!
! boundary value for loc: 0.5 -> right, 1.5 -> bottom, 2.5 -> left, 3.5 -> top
!
! lpoly(i)=j : j is total # of coastline nodes for i-th polygon
!
!!!NOTE that this program assumes top right corner is in the ocean!!!
!
! ind2(i,j,1) =  k : k-th grid node for the starting node (j=1) or end node (j=2) of ith polygon
! ind2(i,j,2) =  1 : gebco node is on the left side
! ind2(i,j,2) =  2 : gebco node is on the upside
!  i=[1:nclose]
!  j=[1:2]
!
!          ind2(i,j,1)     <-   ind2(i,j,2)=1
!            ^
!         ind2(i,j,2)=2
!
! zlabel(i)   : 0 -> depth is zero, 1 -> depth is not zero
! i=[1:ncmax]
!
subroutine addcorner11(h_poly,i_poly,g_grd,g_bound,iflg,ilflg,xcorner,ycorner)
implicit none
type(grid_data),       intent(in)        :: g_grd  ! see mesh/m_coastline_data.f90
type(poly_data),       intent(inout)     :: h_poly ! 1st land polygon is changed when iflag_topright_land=.true.
type(poly_data),       intent(  out)     :: i_poly ! ocean polygon is included
type(bound_data),      intent(inout)     :: g_bound
integer(4),            intent(inout)     :: ilflg(0:3),iflg(0:3)      ! 2019.02.20
real(8),               intent(out)       :: xcorner(0:3),ycorner(0:3) ! 2019.02.25
integer(4)                               :: ncmax, lpmax ! 2018.08.28
integer(4), allocatable,dimension(:)     :: zlabel,ibelong
integer(4), allocatable,dimension(:)     :: lpoly
integer(4), allocatable,dimension(:,:,:) :: ind2
real(8)   , allocatable,dimension(:,:)   :: loc
real(8),    allocatable,dimension(:,:)   :: xpoly,ypoly ! 4 for corner nodes
real(8),    allocatable,dimension(:)     :: xpoly2,ypoly2
! ncoast includes real # of coast points + 4 for corners
real(8),    allocatable,dimension(:)     :: x,y
integer(4)                               :: neast,node,nclose,lpoly0
integer(4)                               :: i,j,ii,jj,lpoly2
logical                                  :: iflag_topright_land  ! 2020.05.29
integer(4)                               :: is,ie,istart_polygon ! 2019.02.21

!#[0]## set input
 node   = g_grd%node
 neast  = g_grd%neast
 ncmax  = g_bound%ncmax     ! 2018.08.28
 lpmax  = h_poly%lpmax      ! 2018.08.28
 nclose = h_poly%nclose     ! # of unclosed polygons
 allocate(x(node),y(node))  ! 2018.08.28
 allocate(zlabel(ncmax),ibelong(ncmax)) ! 2018.08.28
 allocate(xpoly(lpmax,ncmax), ypoly(lpmax,ncmax)) ! 2018.08.28
 allocate(xpoly2(ncmax),      ypoly2(ncmax)  )    ! 2018.08.28
 allocate(lpoly(lpmax))              ! 2018.08.28
 if ( nclose .ge. 1) then            ! 2019.03.06
  allocate(loc(nclose,2)   )         ! 2018.08.28
  allocate(ind2(nclose,2,2))         ! 2018.08.28
  loc    = h_poly%loc                ! 2019.02.25
  ind2(:,:,:) = h_poly%ind2(:,:,:)   ! 2019.02.21
 end if
 x(:)   = g_grd%xyz(2,:)    ! northward [km] 2018.08.28
 y(:)   = g_grd%xyz(1,:)    ! eastward [km]  2018.08.28
 iflag_topright_land = g_grd%iflag_topright_land ! 2019.02.25
 !
 lpoly0      = h_poly%lpoly0
 lpoly(:)    = h_poly%lpoly(:)      ! 2019.02.21
 ypoly(:,:)  = h_poly%xypoly(1,:,:) ! east
 xpoly(:,:)  = h_poly%xypoly(2,:,:) ! north

!do i=1,lpoly0
! write(*,*) "i",i,"lpoly",lpoly(i)
!end do

 if ( iflag_topright_land ) then                               ! 2019.02.20
  write(*,'(a)') " top right corner is in LAND !!"             ! 2019.03.05
  iflg(0:3)  = 999 ! corner is included by land
  ilflg(0:3) = 1   ! corner is included by 1st land
 else                                                          ! 2019.03.05
  write(*,'(a)') " top right corner is in OCEAN !!"            ! 2019.02.20
  iflg(0:3)  = 0   ! all corners between 0 and 1 land polygon  ! 2019.02.21
  ilflg(0:3) = 0   ! all corners are in ocean                  ! 2019.02.21
 end if                                                        ! 2019.02.20

 if ( ncmax .ne. h_poly%ncmax ) goto 998                       ! 2019.02.20
 ii=0


!# [2] ### set iflg(1:3) and ilflg(0:3)
  !# iflg(i) = k  -> [ kth polygon - ith corner - k+1 polygon ]
  !#[2-1]# for the first land polygon
   if ( iflag_topright_land ) then ! 2019.02.21
    do i=1,3
     if ( loc(1,2) .lt. i .and. i .lt. loc(1,1) ) then
      ilflg(i) = 0  ! in ocean
      iflg(i)  = 1  ! between first and second polygon
     end if
    end do
    istart_polygon=2
   else ! top right is in ocean
    istart_polygon=1
   end if

  !# [2-2] ####  istart_polygon to nclose polygon, here is skipped when nclose = 0 2019.03.06
    do jj=istart_polygon,nclose !! investigate if corner is present between jj-th and jj+1-th boundary polygons
     do i=1,3
      if ( loc(jj,1) .lt. i .and. i .lt. loc(jj,2) ) then ! included in land
       ilflg(i)=jj
       iflg(i)=999
	end if
      if ( loc(jj,2) .lt. i ) then ! i th corner is after the end of jj polygon
       ilflg(i)=0
       iflg(i)=jj
      end if
     end do
    end do

  write(*,*) "iflg(0:3)=",(iflg(i),i=0,3)   ! 2019.02.25
  write(*,*) "ilflg(0:3)=",(ilflg(i),i=0,3) ! 2019.02.25

!# [3] #### set corner nodes
  xcorner(0)=x(neast)        ; ycorner(0)=y(neast)        ! top right 2019.02.20
  xcorner(1)=x(node)         ; ycorner(1)=y(node)         ! bottom right corner
  xcorner(2)=x(node-neast+1) ; ycorner(2)=y(node-neast+1) ! bottom left corner
  xcorner(3)=x(1)            ; ycorner(3)=y(1)                               ! top left corner

!#[4]## integrate nclose polygons to one ocean polygon
!# [4-1] #### set top right starting node(s)
 if ( iflag_topright_land ) then ! 2019.02.20
   lpoly2=0
 else
  xpoly2(1)=xcorner(0)          ! 2019.02.02
  ypoly2(1)=ycorner(0)          ! 2019.02.20
  lpoly2=1
  zlabel(1)=1  ! zlabel     = 1 indicates depth here is not zero (in ocean)
  ibelong(1)=1 ! ibelong(i) = 1 indicates boundary line starting from i-th node is not on coastlines
  do j=1,3  ! corner loop 2019.02.21
   if (iflg(j) .eq. 0) then
    xpoly2(lpoly2+1) = xcorner(j)
    ypoly2(lpoly2+1) = ycorner(j)
    lpoly2 = lpoly2 + 1
    zlabel(  lpoly2 ) = 1 ! at added corner, the ocean depth is not zero
    ibelong( lpoly2 ) = 1 ! On lines starting from corner nodes, it is necessary to calculate depths
   end if
  end do
 end if

!# [4-2] #### after the first to the last polygon
 do i=1,nclose !# this loop is skipped when nclose = 0   2019.03.06
  is=lpoly2 + 1
  ie=lpoly2 + lpoly(i)
  xpoly2(is:ie)   = xpoly(i,1:lpoly(i))
  ypoly2(is:ie)   = ypoly(i,1:lpoly(i))
  zlabel(is:ie)   = 0 ! All of the lpoly(i) nodes are on coastlines.
  ibelong(is:ie-1)= 0 ! the lines starting from the first and inner nodes are on coastlines
  ibelong(ie)     = 1 ! the last node is on coast line, but also on the boundaries.
  lpoly2 = ie         ! Then on lines starting from last nodes, it is necessary to calculate ocean depths
  do j=1,3  ! corner loop 2019.02.21
   if (iflg(j) .eq. i) then
    xpoly2(  lpoly2+1 ) = xcorner(j)
    ypoly2(  lpoly2+1 ) = ycorner(j)
    lpoly2 = lpoly2 + 1
    zlabel(  lpoly2 ) = 1 ! at added corner, the ocean depth is not zero
    ibelong( lpoly2 ) = 1 ! On lines starting from corner nodes, it is necessary to calculate depths
   end if
  end do
 end do

!# [5] #### deal with [nclose + 1 : lpoly0] polygons
 write(*,*) "addcorner [5] lpoly0",lpoly0,"nclose",nclose
 if ( lpoly0 - nclose + 1 .ge. 2 ) then ! here is skipped when lpoly0 = nclose = 0, 2019.03.06
  write(*,*) "lpoly0",lpoly0,"nclose",nclose,"lpmax",lpmax,"ncmax",ncmax
  write(*,*) "size(xpoly(:,1))",size(xpoly(:,1))
!  do i=1,lpoly0-nclose
  xpoly(1+1: 1+ lpoly0 - nclose,:) = xpoly(nclose+1:lpoly0,:) ! 20200805
  ypoly(1+1: 1+ lpoly0 - nclose,:) = ypoly(nclose+1:lpoly0,:) ! 20200805
  lpoly(1+1: 1+ lpoly0 - nclose  ) = lpoly(nclose+1:lpoly0  ) ! 20200805
!   xpoly(1+i,:)=xpoly(nclose+i,:) ! 2021.05.24
!   ypoly(1+i,:)=ypoly(nclose+i,:) ! 2021.05.24
!   lpoly(1+i  )=lpoly(nclose+i  ) ! 2021.05.24
!  end do
  if ( lpoly0 - nclose +2 .le. h_poly%lpmax ) then ! set zero for nul space
   xpoly(lpoly0-nclose+2:h_poly%lpmax,:) = 0.d0
   ypoly(lpoly0-nclose+2:h_poly%lpmax,:) = 0.d0
   lpoly(lpoly0-nclose+2:h_poly%lpmax)   = 0
  end if
 end if

 xpoly(1,:)=xpoly2(:) ! 1st ocean polygon
 ypoly(1,:)=ypoly2(:)
 lpoly(1)=lpoly2
 lpoly0 = lpoly0 - nclose + 1    ! nclose polygons has been combined to 1 polygon
 zlabel( lpoly2 + 1 : ncmax)=0   ! after the nclose polygon, all the nodes are on coastlines
 ibelong(lpoly2 + 1 : ncmax)=0   ! after the nclose polygon, all boundary lines are on coastlines
! write(*,*) "lpoly0=",lpoly0,"lpoly(1)=",lpoly(1),"nclose=",nclose,"h_poly%lpoly0=",h_poly%lpoly0

!#[6]## set output
 call allocatepoly(i_poly,lpoly0,ncmax,0) ! nclose for i_poly = 0
 i_poly%lpoly0          = lpoly0
 i_poly%lpoly(1:lpoly0) = lpoly(1:lpoly0)
 i_poly%xypoly(1,1:lpoly0,1:ncmax)   = ypoly(1:lpoly0,1:ncmax)
 i_poly%xypoly(2,1:lpoly0,1:ncmax)   = xpoly(1:lpoly0,1:ncmax)
 i_poly%nclose          = 0
 g_bound%ibelong        = ibelong
 g_bound%zlabel         = zlabel
 i_poly%iflag_topright_land = iflag_topright_land ! 2019.02.21

 write(*,*) "### addcorner11 end!###"

return

998 continue
  write(*,*) "GEGEGE! ncmax=",ncmax,"h_poly%ncmax=",h_poly%ncmax
  stop

end subroutine addcorner11
end module addcorner11m
