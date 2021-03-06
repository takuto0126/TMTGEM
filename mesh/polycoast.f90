! Coded by T. MINAMI on Jan 02, 2014
! input file "coast.dat" is generated by coastline.f90
! format of coast.dat
! x [km], y [km], ind(1) , ind(2)
! x : northward coordinate
! y : eastward coordinate
! ind(1) : grid node # of original gebco file on left/up side
! ind(2) : 1 -> ind(1) indicates left node, 2 -> ind(2) is top side node
program polycoast
implicit none
integer(4) :: i,j,n,l,node
real(8),allocatable,dimension(:) :: x,y
integer(4),allocatable,dimension(:) :: label
integer(4),allocatable,dimension(:,:) :: ind
parameter ( node=8069761 )
real(8),dimension(node) :: gx,gy
!label : polygon number, default = 0, meaning not belonging to any polygon
!
!#### gebco file read ###################
open(2,file="gebco_1min_120_20_176_60.xyz")
read(2,*) (gx(i),gy(i),i=1,node)
close(2)
write(*,*) "gebco read end!"
!#### gebco file read end ###############
!
!####  coastline nodes read #############
j=1
open(1,file="coast.dat")
do while ( j > 0 )
read(1,*,end=99)
j=j+1
end do
99 continue
n=j-1
write(*,*) "# of coastline nodes is",n
allocate (x(n),y(n))
allocate (ind(n,2))
rewind(1)
read(1,*) (x(i),y(i),ind(i,1),ind(i,2),i=1,n)
close(1)
write(*,*) "coast nodes read end!"
!####  coastline nodes read end !!#######
!
!### make polygon #######
allocate (label(n))
l=0
do i=1,n
if (label(i) .eq. 0 ) then
 l=l+1

 if ( l=label()) then
 write(*,*) "polygon", l, "finished!"
 end if
end if
end do
end program polycoast
