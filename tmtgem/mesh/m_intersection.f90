!# coded on 2019.02.27
module intersection
use coastline_data
implicit none

contains

!########################################  avoid intersections_interpolygons
!# coded on 2019.02.27
subroutine avoidintersections_interpolygons(h_poly)
implicit none
type(poly_data),       intent(inout) :: h_poly
integer(4)                           :: i,j,k,l,jj,i1,i2,itype
real(8)                              :: xy1(2),xy2(2),xy3(2),xy4(2),xy5(2),xy6(2)
logical                              :: exist,exist1,exist2,exist3
real(8),allocatable,dimension(:,:,:) :: xypoly
integer(4)                           :: lpoly0,ncmax,lpmax,nclose
integer(4),allocatable,dimension(:)  :: lpoly
real(8),   allocatable,dimension(:,:) :: xypj,xypk
logical                              :: iline_available(5)
integer(4)                           :: j_start,j_end,k_start,k_end
real(8),              dimension(2,4) :: xyk ! 2019.03.04
real(8),              dimension(2,3) :: xyj ! 2019.03.04

!#[0]## set
 lpmax = h_poly%lpmax
 ncmax = h_poly%ncmax
 allocate(xypoly(2,lpmax,ncmax))
 xypoly = h_poly%xypoly
 lpoly0 = h_poly%lpoly0
 allocate( lpoly(lpoly0))
 lpoly  = h_poly%lpoly
 nclose = h_poly%nclose ! 2019.02.28

!#[1]##
 do i1=2,lpoly0          ! mother polygon
  j_start = 1 ; j_end = lpoly(i1)
  allocate( xypj(2,lpoly(i1)))
  xypj = xypoly(1:2,i1,1:lpoly(i1))
  do j=1,lpoly(i1)-1     ! mother line
   xy4(1:2) = xypj(1:2,j  )
   xy5(1:2)=  xypj(1:2,j+1)

   do i2=1,i1-1           ! target polygon
    allocate( xypk(2,lpoly(i1)))
    xypk = xypoly(1:2,i2,1:lpoly(i2))
    k_start=1 ; k_end = lpoly(i2)

    do k=1,lpoly(i2) - 1 ! target line
	xy1(1:2) = xypk(1:2,k)
      xy2(1:2) = xypk(1:2,k+1)
      call intersectionexist(xy1,xy2,xy4,xy5,exist)
      if ( exist ) then
       write(*,*) "intersection exist between polygons",i1,i2
       write(*,'(4(a,i5,1x))') " lines: j/lply(i1)",j,"/",lpoly(i1),"k/lpoly(i2)",k,"/",lpoly(i2)!2019.03.05
       call setsevenpoints(j,xypj,j_start,j_end,k,xypk,k_start,k_end,xyj,xyk,iline_available) ! 2019.03.04
       call avoidintersection7points(xyk,xyj,iline_available,itype,0,0) ! 20200805
       if ( itype .ge. 1) then
        xypk(1:2,k  ) = xyk(:,2) ! 2019.03.04
        xypk(1:2,k+1) = xyk(:,3) ! 2019.03.04
        xypj(1:2,j+1) = xyj(:,2) ! 2019.03.04
       end if
	end if
    end do ! k loop end

    xypoly(1:2,i2,1:lpoly(i2)) = xypk(1:2,1:lpoly(i2))
    deallocate( xypk )
   end do

  end do

  xypoly(1:2,i1,1:lpoly(i1)) = xypj(1:2,1:lpoly(i1))
  deallocate(xypj)
 end do

!#[2]## output 2019.02.27
 h_poly%xypoly = xypoly

write(*,*) "### avoidintersections_interpolygons END!! ###"

return
end
!########################################  avoid intersections
!# coded on 2019.02.27
subroutine avoidintersections(h_poly)
implicit none
type(poly_data),       intent(inout) :: h_poly
integer(4)                           :: i,j,k,jj
real(8)                              :: xy1(2),xy2(2),xy3(2),xy4(2),xy5(2),xy6(2)
logical                              :: exist,exist1
real(8),allocatable,dimension(:,:,:) :: xypoly
integer(4)                           :: lpoly0,ncmax,lpmax
integer(4),allocatable,dimension(:)  :: lpoly
integer(4)                           :: nclose,k_start,j_start,j_end,itype
real(8),  allocatable,dimension(:,:) :: xyp
logical,              dimension(5)   :: iline_available ! 2019.03.04
real(8),              dimension(2,4) :: xyk ! 2019.03.04
real(8),              dimension(2,3) :: xyj ! 2019.03.04

!#[0]## set
 lpmax = h_poly%lpmax
 ncmax = h_poly%ncmax
 allocate(xypoly(2,lpmax,ncmax))
 xypoly = h_poly%xypoly
 lpoly0 = h_poly%lpoly0
 allocate( lpoly(lpoly0))
 lpoly  = h_poly%lpoly
 nclose = h_poly%nclose

!#[1]##
do i=1,lpoly0
 if ( i .le. nclose ) then ! unclosed polygons
  allocate( xyp(2,1:lpoly(i))) ! 2019
  j_start=1
 else                      ! closed polygons
  allocate( xyp(2,0:lpoly(i)) )                  ! 2019.02.28
  j_start=0
  xyp(1:2,0)=xypoly(1:2,i,lpoly(i)) ! 2019.02.28
 end if
 j_end = lpoly(i)
 xyp(1:2,1:j_end)=xypoly(1:2,i,1:j_end)

 do j = j_start,lpoly(i)-1 ! j_start= 1 for unclosed polygon, 0 for closed polygon
!  write(*,*) "i,j",i,j,"lpoly(i)",lpoly(i)
  xy4(1:2) = xyp(1:2,j  )
  xy5(1:2)=  xyp(1:2,j+1)
  k_start = 1
  if ( i .gt. nclose .and. .not. j .eq. lpoly(i) - 1 ) k_start=0 ! for closed polygon
  do k = k_start,j-2
   xy1(1:2) = xyp(1:2,k)
   xy2(1:2) = xyp(1:2,k+1)
   call intersectionexist(xy1,xy2,xy4,xy5,exist)
!   if ( .false. ) then ! i .eq. ** can be used for debugging 20200805
!    write(*,*) "--"
!    write(*,*) "j,k,",j,k,"exist",exist
!    write(*,*) "j  ",j,"xy",  xyp(1:2,j  )
!    write(*,*) "j+1",j,"xy",  xyp(1:2,j+1)
!    write(*,*) "k  ",k,"xy",  xyp(1:2,k  )
!    write(*,*) "k+1",k+1,"xy",  xyp(1:2,k+1)
!    write(*,*) "--"
!   end if

   if ( exist ) then
    ! set xyj, xyk, iline_available
    write(*,*) "intersection exist within polygon",i
    write(*,*) "lines: j,k / lpoly(i)",j,k,"/",lpoly(i)
!    write(*,*) "xy1",xy1
!    write(*,*) "xy2",xy2
!    write(*,*) "xy4",xy4
!    write(*,*) "xy5",xy5
    call setsevenpoints(j,xyp,j_start,j_end,k,xyp,j_start,j_end,xyj,xyk,iline_available)

    !#[1]## whether twisted case or no
    call lineintersectseven(xyk,xyj,iline_available,2,5,exist1) ! 2019.03.04
    if ( exist1 .eqv. .false. ) then  ! 2020.05.29 "eq"->"eqv" line 2 and line 5 are not intersected
    if (  j - k .eq. 2 .or.  (j_end + k) - j .eq. 2 .and. i .gt. nclose  ) then! avoid same lines
     write(*,*) "j,k",j,k,"intersection (type twist) exist!!"
     write(*,*) "xyk2",xyk(1:2,2)
     write(*,*) "xyk3",xyk(1:2,3)
     write(*,*) "xyj1",xyj(1:2,1)
     write(*,*) "xyj2",xyj(1:2,2)
     xyp(1:2,k+1) = xyj(1:2,1)
     xyp(1:2,j  ) = xyk(1:2,3)
     write(*,*) "->"
     write(*,*) "xyk3",xyj(1:2,1) ! 2019.03.04
     write(*,*) "xyj1",xyk(1:2,3) ! 2019.03.04
     goto 100
    end if
    end if

    !#[2]## two point intersections
    write(*,*) "j,k",j,k,"intersection (not twist) exist!!"
    call avoidintersection7points(xyk,xyj,iline_available,itype,k,j) ! 2019.03.04
    if ( itype .ge. 1) then
     xyp(1:2,k  ) = xyk(:,2) ! 2019.03.04
     xyp(1:2,k+1) = xyk(:,3) ! 2019.03.04
     xyp(1:2,j+1) = xyj(:,2) ! 2019.03.04
    end if
    100 continue

   end if ! case for exist = .true.
  end do ! k loop
 end do  ! j loop

 h_poly%xypoly(1:2,i,1:lpoly(i)) = xyp(1:2,1:lpoly(i))
 deallocate(xyp)
end do

write(*,*) "### avoidintersections end!! ###"

return
end

!######################################## setsevenpoints
!# k-1  k  k+1 k+2   j  j+1 j+2
!#  |---|---|---|    |---|---|
!#    1   2   3        4   5   <- line id
!#  ^   ^   ^   ^    ^   ^   ^
!#  1   2   3   4    1   2   3
!#  xyk(1:2,1:4)     xyj(1:2,1:3)
!#
!# assume line 2 and 4 intersect each other
!# coded on 2019.03.04
 subroutine setsevenpoints(j,xypj,j_start,j_end,k,xypk,k_start,k_end,xyj,xyk,iline_available)
 implicit none
 integer(4),intent(in)  :: j,k,j_start,j_end,k_start,k_end
 real(8),   intent(in)  :: xypj(2,j_start:j_end)
 real(8),   intent(in)  :: xypk(2,k_start:k_end)
 real(8),   intent(out) :: xyj(2,3),xyk(2,4)
 logical,   intent(out) :: iline_available(5)

!    write(*,*) "## setsevenpoints start!! ###"
    !# default line 2 and line 4
     xyk=0.d0
     xyj=0.d0
     iline_available(:)=.false.
     iline_available(2)=.true.
     iline_available(4)=.true.
     xyk(1:2,2)=xypk(1:2,k)
     xyk(1:2,3)=xypk(1:2,k+1)
     xyj(1:2,1)=xypj(1:2,j)
     xyj(1:2,2)=xypj(1:2,j+1)

     ! line 1
     if ( k-1 .ge. k_start ) then
      iline_available(1) = .true.
      xyk(1:2,1)=xypk(1:2,k-1)
     end if

     ! line 3
     if ( k+2 .le. k_end) then
      iline_available(3) = .true.
      xyk(1:2,4)=xypk(1:2,k+2)
     end if

     ! line 5
     if ( j + 2 .le. j_end ) then
      iline_available(5) = .true.
      xyj(1:2,3)=xypj(1:2,j+2)
     end if

     if ( .false. ) then
      write(*,*) "xyk1",xyk(1:2,1)
      write(*,*) "xyk2",xyk(1:2,2)
      write(*,*) "xyk3",xyk(1:2,3)
      write(*,*) "xyk4",xyk(1:2,4)
      write(*,*) "xyj1",xyj(1:2,1)
      write(*,*) "xyj2",xyj(1:2,2)
      write(*,*) "xyj3",xyj(1:2,3)
      write(*,'(a,5l3)') "iline_available",iline_available(1:5)
      write(*,*) "## set sevenpoints end!! ##"
     end if

     return
     end

!####################################### lineintersectseven
subroutine lineintersectseven(xyk,xyj,iline_available,l1,l2,exist)
implicit none
integer(4),intent(in)  :: l1 ! 1-3
integer(4),intent(in)  :: l2 ! 4-5 ! line id see setsevenpoints
logical,   intent(in)  :: iline_available(5)
real(8),   intent(in)  :: xyj(2,3),xyk(2,4)
logical,   intent(out) :: exist
real(8)                :: xy1(2),xy2(2),xy3(2),xy4(2)
integer(4)             :: i1,i2

if ( .not. (1 .le. l1 .and. l1 .le. 3 )) goto 100
if ( .not. (4 .le. l2 .and. l2 .le. 5 )) goto 101

exist=.false.

if ( iline_available(l1) ) then
 xy1(1:2) = xyk(1:2,l1  )
 xy2(1:2) = xyk(1:2,l1+1)
else
 return
end if

if ( iline_available(l2) ) then
 xy3(1:2) = xyj(1:2,l2-3)
 xy4(1:2) = xyj(1:2,l2-2)
else
 return
end if

call intersectionexist(xy1,xy2,xy3,xy4,exist)


return

100 continue
write(*,*) "GEGEGE l1 shold be [1:3]; l1=",l1
stop
101 continue
write(*,*) "GEGEGE l2 shold be [4:5]; l2=",l2
stop
end

!########################################  avoidintersection7points
!# 2019.03.04
!# assume xy1 -> xy2 intersect xy4 -> xy5
subroutine avoidintersection7points(xyk,xyj,iline_available,itype,k,j) ! 20200805
implicit none
real(8),dimension(2),intent(inout) :: xyk(2,4),xyj(2,3)
logical,             intent(in)    :: iline_available(5)
integer(4),          intent(out)   :: itype
logical                            :: exist,exist1,exist2,exist3,exist4
integer(4),          intent(in)    :: k,j

!#[0]## check

     if ( .true. ) then
      write(*,*) "## avoidintersection7points start!! ##"
      write(*,*) "xyk1",xyk(1:2,1)
      write(*,*) "xyk2",xyk(1:2,2)
      write(*,*) "xyk3",xyk(1:2,3)
      write(*,*) "xyk4",xyk(1:2,4)
      write(*,*) "xyj1",xyj(1:2,1)
      write(*,*) "xyj2",xyj(1:2,2)
      write(*,*) "xyj3",xyj(1:2,3)
      write(*,'(a,5l3)') " iline_available",iline_available(1:5)
     end if


!#[1]## check intersection type
      call lineintersectseven(xyk,xyj,iline_available,2,4,exist)
      call lineintersectseven(xyk,xyj,iline_available,3,4,exist1) ! type 1
      call lineintersectseven(xyk,xyj,iline_available,1,5,exist2) ! type 2
      call lineintersectseven(xyk,xyj,iline_available,2,5,exist3) ! type 3
      call lineintersectseven(xyk,xyj,iline_available,3,5,exist4) ! type 4
      if ( j-k .eq. 2 ) exist1 = .false. ! 20200805
      write(*,'(4(a,l3))') "exist",exist," exist1",exist1," exist2",exist2," exist3",exist3," exist4",exist4

!#[2]##
      itype = 0

      if ( .not. exist ) then
       write(*,*) "GEGEGE 2,4 should be intersected!"
       write(*,'(a,l3)') "exist",exist
       stop
      end if

      if (     exist1 ) then
       itype = 1
       write(*,'(a)') " type 1"
       write(*,*) "k+1,xyk(1:2,3)",xyk(1:2,3)
       xyk(1:2,3) = xyk(1:2,2)/3. + xyk(1:2,3)/3. + xyk(1:2,4)/3. ! k + 1
       write(*,*) "-> k+1 ",xyk(1:2,3)

      elseif ( exist2 ) then
       itype = 2
       write(*,'(a)') " type 2"
       write(*,*) "k,   xyk(1:2,2)",xyk(1:2,2)
       write(*,*) "j+1, xyj(1:2,2)",xyj(1:2,2)
       xyk(1:2,2) = xyk(1:2,1)/3. + xyk(1:2,2)/3. + xyk(1:2,3)/3. ! k
       xyj(1:2,2) = xyj(1:2,1)/3. + xyj(1:2,2)/3. + xyj(1:2,3)/3. ! j + 1
       write(*,*) "-> k  ",xyk(1:2,2)
       write(*,*) "-> j+1",xyj(1:2,2)

	elseif ( exist3 ) then
       write(*,'(a)') " type 3"
	 itype = 3
       write(*,*) "j+1,xyj(1:2,2)",xyj(1:2,2)
       xyj(1:2,2) = xyj(1:2,1)/3. + xyj(1:2,2)/3. + xyj(1:2,3)/3. ! j + 1
       write(*,*) " -> j+1 ",xyj(1:2,2)

	elseif ( exist4 ) then
	 itype = 4
       write(*,'(a)') " ## type 4 ##"
       write(*,*) "k+1, xyk(1:2,2)",xyk(1:2,3)
       write(*,*) "j+1, xyj(1:2,2)",xyj(1:2,2)
       xyk(1:2,3) = xyk(1:2,2)/3. + xyk(1:2,3)/3. + xyk(1:2,4)/3. ! k + 1
       xyj(1:2,2) = xyj(1:2,1)/3. + xyj(1:2,2)/3. + xyj(1:2,3)/3. ! j + 1
       write(*,*) "-> k+1",xyk(1:2,3)
       write(*,*) "-> j+1",xyj(1:2,2)

      else
       write(*,'(4(a,l3))') " exist1",exist1," exist2",exist2," exist3",exist3,"exist4",exist4
       write(*,'(a)') " GEGEGE! not followed"
      end if

      write(*,*) "## avoid intersection7points end!! ##"
return
end
!########################################  intersectionEX
!# check intersection between line xy1 -> xy2 and line xy3 -> xy4
!# http://www5d.biglobe.ne.jp/~tomoya03/shtml/algorithm/IntersectionEX.htm
!# NOTE that to share a point by two lines are recognized as intersection
subroutine intersectionexist(xy1,xy2,xy3,xy4,exist)
implicit none
real(8), intent(in)  :: xy1(2),xy2(2),xy3(2),xy4(2)
logical, intent(out) :: exist
integer(4)           :: i

do i=1,2
if (xy1(i) >= xy2(i)) then
  if ( (xy1(i) < xy3(i) .and. xy1(i) < xy4(i)) .or. (xy2(i) > xy3(i) .and. xy2(i) > xy4(i)) ) then
   exist=.false.
   return
  end if
else
  if ( (xy2(i) < xy3(i) .and. xy2(i) < xy4(i)) .or. (xy1(i) > xy3(i) .and. xy1(i) > xy4(i))) then
   exist = .false.
   return
  end if
end if
end do

if( ((xy1(1) - xy2(1))*(xy3(2)-xy1(2))+(xy1(2)-xy2(2))*(xy1(1)-xy3(1)))* &
&   ((xy1(1) - xy2(1))*(xy4(2)-xy1(2))+(xy1(2)-xy2(2))*(xy1(1)-xy4(1))) > 0. ) then
   exist=.false.
   return
end if
if( ((xy3(1) - xy4(1))*(xy1(2)-xy3(2))+(xy3(2)-xy4(2))*(xy3(1)-xy1(1)))* &
&   ((xy3(1) - xy4(1))*(xy2(2)-xy3(2))+(xy3(2)-xy4(2))*(xy3(1)-xy2(1))) > 0. ) then
   exist=.false.
   return
end if
exist=.true.

return
end

end module intersection
