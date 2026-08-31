module matrix
implicit none

type real_crs_matrix ! compressed row matrix
 integer(4) :: nrow
 integer(4) :: ncolm ! added on 2017.05.11 (required for mat-mat multiplication)
 integer(4) :: ntot
 integer(4),dimension(:),allocatable :: item !item(1:ntot) :colmun id for every entry for
 integer(4),dimension(:),allocatable :: stack!stack(0:nrow) accumulated # of row members
 real(8),dimension(:),allocatable :: val  !val(1:ntot) : values
end type

type complex_crs_matrix ! compressed row matrix (2017.05.15)
 integer(4) :: nrow
 integer(4) :: ncolm ! added on 2017.05.11 (required for mat-mat multiplication)
 integer(4) :: ntot
 integer(4),dimension(:),allocatable :: item !item(1:ntot) :colmun id for every entry for
 integer(4),dimension(:),allocatable :: stack!stack(0:nrow) accumulated # of row members
 complex(8),dimension(:),allocatable :: val  !val(1:ntot) : values
end type

type real_ccs_matrix ! compressed colmun matrix, added on 2017.05.10
 integer(4) :: nrow ! added on 2017.05.11 (required for mat-mat multiplication)
 integer(4) :: ncolm
 integer(4) :: ntot
 integer(4),dimension(:),allocatable :: item !item(1:ntot) :colmun id for every entry for
 integer(4),dimension(:),allocatable :: stack!stack(0:nrow) accumulated # of row members
 real(8),dimension(:),allocatable :: val  !val(1:ntot) : values
end type

type complex_ccs_matrix ! compressed colmun matrix, added on 2017.05.15
 integer(4) :: nrow ! added on 2017.05.11 (required for mat-mat multiplication)
 integer(4) :: ncolm
 integer(4) :: ntot
 integer(4),dimension(:),allocatable :: item !item(1:ntot) :colmun id for every entry for
 integer(4),dimension(:),allocatable :: stack!stack(0:nrow) accumulated # of row members
 complex(8),dimension(:),allocatable :: val  !val(1:ntot) : values
end type

contains
!#############################################  allocate_real_crs_with_steady_ncolm
!# coded on May 20, 2016
!# when the # of non zero members in every colmun
!# set mcrs%nrow, mcrs%ntot, and allocate mcrs%stack, mcrs%item, mcrs%val
!# and set mcrs%stack
subroutine allocate_real_crs_with_steady_ncolm(mcrs,nrow,ncolm)
implicit none
integer(4),intent(in) :: nrow,ncolm
type(real_crs_matrix),intent(inout) :: mcrs
integer(4) :: i
 mcrs%nrow=nrow ! number of observation points
 mcrs%ntot=nrow * ncolm ! for lines of tetrahedral mesh
 allocate(mcrs%stack(0:mcrs%nrow))
 allocate(mcrs%item(   mcrs%ntot))
 allocate(mcrs%val(    mcrs%ntot))
 mcrs%stack(0)=0
 do i=1,nrow
  mcrs%stack(i)=i*ncolm
 end do
return
end subroutine

!#############################################  mul_matcrs_rv
subroutine mul_matcrs_rv(crsmat,rv,doftot,vout)
implicit none
type(real_crs_matrix),intent(in) :: crsmat
integer(4),intent(in)  :: doftot
real(8),   intent(in)  :: rv(doftot)
real(8),   intent(out) :: vout(crsmat%nrow)
integer(4) :: irow,jtot
vout(:)=0.d0
do irow=1,crsmat%nrow
! write(*,*) "irow=",irow,"crsmat%nrow"
 do jtot=crsmat%stack(irow-1)+1,crsmat%stack(irow)
!  write(*,*) "jtot=",jtot,"crsmat%item(jtot)=",crsmat%item(jtot),"doftot=",doftot
  vout(irow)=vout(irow) + crsmat%val(jtot)*rv(crsmat%item(jtot))
 end do
end do
return
end subroutine

!#############################################  mul_matcrs_cv
subroutine mul_matcrs_cv(crsmat,cv,doftot,vout)
implicit none
type(real_crs_matrix),intent(in) :: crsmat
integer(4),intent(in)   :: doftot
complex(8),  intent(in) :: cv(doftot)
complex(8),  intent(out):: vout(crsmat%nrow)
integer(4) :: irow,jtot
!#[1] decompose cv to rcv and icv
vout(:)=(0.d0,0.d0)
do irow=1,crsmat%nrow
! write(*,*) "irow=",irow,"crsmat%nrow=",crsmat%nrow
 do jtot=crsmat%stack(irow-1)+1,crsmat%stack(irow)
! write(*,*) "jtot=",jtot,"crsmat%item(jtot)=",crsmat%item(jtot),"cv(crsmat%item(jtot))=",cv(crsmat%item(jtot))
  vout(irow)=vout(irow) + crsmat%val(jtot)*cv(crsmat%item(jtot))
 end do
end do
return
end subroutine
!#############################################  mul_matcrscomp_cv
!# OMP is applied on 2017.07.21
subroutine mul_matcrscomp_cv(crsmat,cv,doftot,vout)
implicit none
type(complex_crs_matrix),intent(in)  :: crsmat
integer(4),              intent(in)  :: doftot
complex(8),              intent(in)  :: cv(doftot)
complex(8),              intent(out) :: vout(crsmat%nrow)
integer(4) :: irow,jtot
integer(4) :: chunk    ! 2017.07.21
!#[1] decompose cv to rcv and icv
vout(:)=(0.d0,0.d0)

chunk = 1      ! 2017.07.21
!$OMP PARALLEL SHARED(crsmat,vout,cv,CHUNK) PRIVATE(irow,jtot) ! 2017.07.21
!!  TID = OMP_GET_THREAD_NUM()                                 ! 2017.07.21
!$OMP DO SCHEDULE (STATIC,chunk)                               ! 2017.07.21

do irow=1,crsmat%nrow
! write(*,*) "irow=",irow,"crsmat%nrow=",crsmat%nrow
 do jtot=crsmat%stack(irow-1)+1,crsmat%stack(irow)
! write(*,*) "jtot=",jtot,"crsmat%item(jtot)=",crsmat%item(jtot),"cv(crsmat%item(jtot))=",cv(crsmat%item(jtot))
  vout(irow)=vout(irow) + crsmat%val(jtot)*cv(crsmat%item(jtot))
 end do
end do
!$OMP END PARALLEL ! 2017.07.21

return
end subroutine

!###############################################
! Coded on 2017.05.12
subroutine conv_full2crs(full,nrow,ncolm,crs,threshold)
implicit none
real(8),   intent(in)  :: threshold
real(8),   intent(in)  :: full(nrow,ncolm)
integer(4),intent(in)  :: nrow,ncolm
type(real_crs_matrix),   intent(out) :: crs
integer(4) :: i,j,row_count(nrow),ii
integer(4) :: item_work(nrow,ncolm),ntot

!#[1]## count
row_count = 0 ; ntot=0
do i=1,nrow
 do j=1,ncolm
  if (abs(full(i,j)) .gt. threshold) then
   row_count(i) = row_count(i) + 1
   ii = row_count(i)
   item_work(i,ii) = j
  end if
 end do
 ntot = ntot + row_count(i)
end do

!#[2]## generate crs
 crs%ntot = ntot
 allocate(crs%val(ntot),crs%item(ntot),crs%stack(0:nrow))
 crs%stack(0)=0
 do i=1,nrow
  crs%stack(i) = crs%stack(i-1) + row_count(i)
  do j=1,row_count(i)
  crs%item(crs%stack(i-1) + j) = item_work(i,j)
  crs%val(crs%stack(i-1) + j)  = full(i,item_work(i,j))
  end do
 end do

!#[3]## set out
 crs%nrow  = nrow
 crs%ncolm = ncolm

return
end subroutine

!###############################################
! Coded on 2017.06.16
subroutine conv_full2ccs(full,nrow,ncolm,ccs,threshold)
implicit none
real(8),   intent(in)  :: threshold
real(8),   intent(in)  :: full(nrow,ncolm)
integer(4),intent(in)  :: nrow,ncolm
type(real_ccs_matrix),   intent(out) :: ccs
integer(4) :: i,j,colm_count(ncolm),ii
integer(4) :: item_work(ncolm,nrow),ntot

!#[1]## count
colm_count = 0 ; ntot=0
do i=1,ncolm
 do j=1,nrow
  if (abs(full(j,i)) .gt. threshold) then
   colm_count(i) = colm_count(i) + 1
   ii = colm_count(i)
   item_work(i,ii) = j
!   write(*,*) "i,j",i,j,"colm_count=",colm_count(i)
  end if
 end do
 ntot = ntot + colm_count(i)
end do

!#[2]## generate crs
 ccs%ntot = ntot
 allocate(ccs%val(ntot),ccs%item(ntot),ccs%stack(0:ncolm))
 ccs%stack(0)=0
 do i=1,ncolm
  ccs%stack(i) = ccs%stack(i-1) + colm_count(i)
  do j=1,colm_count(i)
  ccs%item(ccs%stack(i-1) + j) = item_work(i,j)
  ccs%val(ccs%stack(i-1) + j)  = full(item_work(i,j),i)
  end do
 end do

!#[3]## set out
 ccs%nrow  = nrow
 ccs%ncolm = ncolm

return
end subroutine

!############################################### add_crs_crs_crs
! Coded on 2017.05.17
subroutine add_crs_crs_crs(crs1,crs2,crs3)
implicit none
type(real_crs_matrix),intent(in)  :: crs1,crs2
type(real_crs_matrix),intent(out) :: crs3
integer(4) :: nrow1,ncolm1,ntot1
integer(4) :: nrow2,ncolm2,ntot2
integer(4),allocatable,dimension(:) :: stack1,stack2,item1,item2
real(8),   allocatable,dimension(:) :: val1,val2
integer(4) :: ntot3w
integer(4),allocatable,dimension(:) :: itemw,stackw
real(8),   allocatable,dimension(:) :: valw
integer(4) :: i,j1,j2,jcol1,jcol2,icount

!#[0]## set
nrow1  = crs1%nrow
ncolm1 = crs1%ncolm
ntot1  = crs1%ntot
nrow2  = crs2%nrow
ncolm2 = crs2%ncolm
ntot2  = crs2%ntot
allocate(stack1(0:nrow1),item1(ntot1),val1(ntot1))
allocate(stack2(0:nrow2),item2(ntot2),val2(ntot2))
stack1 = crs1%stack ; item1 = crs1%item ; val1 = crs1%val
stack2 = crs2%stack ; item2 = crs2%item ; val2 = crs2%val

!#[1]## check
if (nrow1  .ne. nrow2 ) goto 99
if (ncolm1 .ne. ncolm2) goto 99

!#[2]## count ntot3
 ntot3w = ntot1 + ntot2
 allocate(itemw(ntot3w),valw(ntot3w),stackw(0:nrow1))
 icount = 0 ; stackw(0)=0
 do i=1,nrow1
!  write(*,*) "i=",i
  j1 = 1; j2 = 1
  100 continue
  if ( j1 .gt. stack1(i)-stack1(i-1) .and. &
  &    j2 .gt. stack2(i)-stack2(i-1) ) cycle
!  write(*,*) "j1=",j1,"j2=",j2
  jcol1 = item1(stack1(i-1)+j1)
  jcol2 = item2(stack2(i-1)+j2)
  if ( j1 .gt. stack1(i)-stack1(i-1) ) jcol1 = ncolm1 +1
  if ( j2 .gt. stack2(i)-stack2(i-1) ) jcol2 = ncolm1 +1
  !# case 1
  if ( jcol1 .eq. jcol2 ) then
   icount = icount + 1 ; itemw(icount) = jcol1 ; stackw(i)=icount
   valw(icount) = val1(stack1(i-1)+j1) + val2(stack2(i-1)+j2)
   j1 = j1 + 1; j2 = j2 + 1 ; goto 100
  end if
  !# case 2
  if ( jcol1 .lt. jcol2 ) then
   icount = icount + 1 ; itemw(icount) = jcol1 ; stackw(i)=icount
   valw(icount) = val1(stack1(i-1)+j1)
   j1 = j1 + 1;  goto 100
  end if
  !# case 3
  if ( jcol1 .gt. jcol2 ) then
   icount = icount + 1 ; itemw(icount) = jcol2 ; stackw(i)=icount
   valw(icount) = val2(stack2(i-1)+j2)
   j2 = j2 + 1;  goto 100
  end if
 end do

!#[3]## set crs3
 crs3%nrow  = nrow1
 crs3%ncolm = ncolm1
 crs3%ntot  = icount
 allocate(crs3%stack(0:nrow1),crs3%item(icount),crs3%val(icount))
 crs3%stack          = stackw
 crs3%item(1:icount) = itemw(1:icount)
 crs3%val(1:icount)  = valw(1:icount)

return
99 continue
 write(*,*) "GEGEGE rank mismatch"
 write(*,*) "nrow1=",nrow1,"ncolm1=",ncolm1
 write(*,*) "nrow2=",nrow2,"ncolm2=",ncolm2
 stop
end subroutine
!###############################################
! Coded on 2017.05.15
subroutine convcomp_full2crs(full,nrow,ncolm,crs,threshold)
implicit none
complex(8),                 intent(in)  :: full(nrow,ncolm)
integer(4),                 intent(in)  :: nrow,ncolm
type(complex_crs_matrix),   intent(out) :: crs
real(8),                    intent(in)  :: threshold
integer(4),allocatable,dimension(:)   :: row_count
integer(4),allocatable,dimension(:,:) :: item_work
integer(4) :: i,j,ii,ntot

!write(*,*) "nrow=",nrow
!write(*,*) "ncolm=",ncolm
allocate(row_count(nrow),item_work(nrow,ncolm))

!#[1]## count
row_count(:) = 0 ; ntot=0
do i=1,nrow
 do j=1,ncolm
  if (abs(full(i,j)) .gt. threshold) then
   row_count(i) = row_count(i) + 1
   ii = row_count(i)
   item_work(i,ii) = j
  end if
 end do
 ntot = ntot + row_count(i)
end do

!#[2]## generate crs
 crs%ntot = ntot
 allocate(crs%val(ntot),crs%item(ntot),crs%stack(0:nrow))
 crs%stack(0)=0
 do i=1,nrow
  crs%stack(i) = crs%stack(i-1) + row_count(i)
  do j=1,row_count(i)
  crs%item(crs%stack(i-1) + j) = item_work(i,j)
  crs%val(crs%stack(i-1) + j)  = full(i,item_work(i,j))
  end do
 end do

!#[3]## set out
 crs%nrow  = nrow
 crs%ncolm = ncolm

return
end subroutine
!######################################################### lateralcomb_real_crs_mat
!# 2017.11.01
subroutine lateralcomb_real_crs_mat(crs1,crs2,crs3)
implicit none
type(real_crs_matrix),intent(in)  :: crs1
type(real_crs_matrix),intent(in)  :: crs2
type(real_crs_matrix),intent(out) :: crs3
integer(4)  :: ntot1,ntot2,ntot3
integer(4)  :: ncolm1,ncolm2,ncolm3
integer(4)  :: nrow1,nrow2,nrow3
integer(4)  :: istack1,istack2,istack3
integer(4)  :: nstack1,nstack2,nstack3
integer(4)  :: i

!#[0] set
ntot1  = crs1%ntot
ntot2  = crs2%ntot
ntot3  = ntot1 + ntot2
nrow1  = crs1%nrow
nrow2  = crs2%nrow
nrow3  = nrow1
ncolm1 = crs1%ncolm
ncolm2 = crs2%ncolm
ncolm3 = ncolm1 + ncolm2
if ( nrow1 .ne. nrow2 ) then
 write(*,*) "GEGEGE nrow1",nrow1," .ne. nrow2",nrow2
 stop
end if

!#[1]## allocate crs3
crs3%ntot  = ntot3
crs3%nrow  = nrow3
crs3%ncolm = ncolm3
allocate(crs3%stack(0:nrow3),crs3%item(ntot3),crs3%val(ntot3))

!#[2]## cal crs3
 crs3%stack(0) = 0
 do i=1,nrow3
  istack1 = crs1%stack(i-1)
  istack2 = crs2%stack(i-1)
  istack3 = crs3%stack(i-1)
  nstack1 = crs1%stack(i) - istack1
  nstack2 = crs2%stack(i) - istack2
  nstack3 = nstack1 + nstack2
  crs3%stack(i) = nstack3 + istack3
  crs3%item(istack3+1:istack3+nstack1)=crs1%item(istack1+1:istack1+nstack1)
  crs3%val( istack3+1:istack3+nstack1)=crs1%val( istack1+1:istack1+nstack1)
  crs3%item(istack3+nstack1+1:istack3+nstack3)=crs2%item(istack2+1:istack2+nstack2)+ncolm1
  crs3%val( istack3+nstack1+1:istack3+nstack3)=crs2%val( istack2+1:istack2+nstack2)
 end do

return
end subroutine
!######################################################### add_real_crs_mat
!# ntet of rmat should be set before this subrouinte
subroutine combine_real_crs_mat(rmat,pmat,irow)
implicit none
integer(4) :: irow ! starting row index
type(real_crs_matrix),intent(in)    :: pmat
type(real_crs_matrix),intent(inout) :: rmat
integer(4) :: nrowp,ntotp,ntotb

!#[1]## check
 if ( rmat%ncolm .ne. pmat%ncolm ) then
  write(*,*) "GEGEGE! rmat%ncolm=",rmat%ncolm,"pmat%ncolm",pmat%ncolm
  stop
 else if ( rmat%nrow - (irow -1) .lt. pmat%nrow ) then
  write(*,*) "GEGEGE! rmat%nrow=",rmat%nrow,"pmat%nrow",pmat%nrow,"irow",irow
  stop
 else if ( rmat%ntot - rmat%stack(irow-1) .lt. pmat%ntot ) then
  write(*,*) "GEGEGE! rmat%ntot=",rmat%ntot,"rmat%stack(irow-1)",&
		& rmat%stack(irow-1),"pmat%ntot",pmat%ntot,"irow",irow
 end if

!#[2]## set pmat to rmat from irow-th row
 nrowp = pmat%nrow
 ntotp = pmat%ntot
 ntotb = rmat%stack(irow-1)
 rmat%stack(irow:irow-1+nrowp) =  ntotb + pmat%stack(1:nrowp)
 rmat%item(ntotb+1:ntotb+ntotp) = pmat%item(1:ntotp)
 rmat%val(ntotb+1:ntotb+ntotp)  = pmat%val(1:ntotp)

return
end subroutine

!######################################################### combine_real_ccs_mat
!# ntot of rmat should be set before this subrouinte
subroutine combine_real_ccs_mat(rmat,pmat,icolm)
implicit none
integer(4) :: icolm ! starting row index
type(real_ccs_matrix),intent(in)    :: pmat
type(real_ccs_matrix),intent(inout) :: rmat
integer(4) :: ncolmp,ntotp,ntotb

!#[1]## check
 if ( rmat%nrow .ne. pmat%nrow ) then
  write(*,*) "GEGEGE! rmat%nrow=",rmat%nrow,"pmat%nrow",pmat%nrow
  stop
 else if ( rmat%ncolm - (icolm -1) .lt. pmat%ncolm ) then
  write(*,*) "GEGEGE! rmat%ncolm=",rmat%ncolm,"pmat%ncolm",pmat%ncolm,"icolm",icolm
  stop
 else if ( rmat%ntot - rmat%stack(icolm-1) .lt. pmat%ntot ) then
  write(*,*) "GEGEGE! rmat%ntot=",rmat%ntot,"rmat%stack(icolm-1)",&
		& rmat%stack(icolm-1),"pmat%ntot",pmat%ntot,"icolm",icolm
 end if

!#[2]## set pmat to rmat from icolm-th colm
 ncolmp = pmat%ncolm
 ntotp  = pmat%ntot
 ntotb  = rmat%stack(icolm-1)
 rmat%stack(icolm:icolm-1+ncolmp) = ntotb + pmat%stack(1:ncolmp)
 rmat%item(ntotb+1:ntotb+ntotp)   = pmat%item(1:ntotp)
 rmat%val(ntotb+1:ntotb+ntotp)    = pmat%val(1:ntotp)

return
end subroutine

!###############################################
! Modified on 2017.12.25
! Coded on 2017.05.11
! Please comment out include line and set iflag = 0
! multiply crs by ccs to gen crs
! --------- When T is given ---------- 2017.12.25
! If T = 'N' or 'n', then C := A*B
! If T = 'T' or 't' or 'C' or 'c', then C := AT*B.
!------------------------------------- 2017.12.25
subroutine mulreal_crs_crs_crs(crsmat1,crsmat2,crsout,T) ! T is added on 2017.12.25
implicit none
!include "mkl.fi" ! 2017.09.06
type(real_crs_matrix),intent(in)  :: crsmat1,crsmat2
type(real_crs_matrix),intent(out) :: crsout
character(1),optional,intent(in)  :: T
type(real_ccs_matrix) :: ccsmat1,ccsmat2
!===== below are for mkl_dcsrmulcsr ! 2017.09.06
real(8),   allocatable,dimension(:) :: a,b,c
integer(4),allocatable,dimension(:) :: ia,ib,ic ! 2017.09.06
integer(4),allocatable,dimension(:) :: ja,jb,jc ! 2017.09.06
integer(4)   :: iflag = 0 ! 0 for normal, 1 for when mkl library available
character(1) :: trans
integer(4)   :: sort,request,m1,n1,k1,m2,n2,k2,m3,n3,nzmax,info
!================================================

if ( iflag .eq. 0 ) then    ! 2017.09.06
 call conv_crs2ccs(crsmat2,ccsmat1)
! write(*,*) "conv_crs2ccs end!"
 call mulreal_crs_ccs_ccs(crsmat1,ccsmat1,ccsmat2)
! write(*,*) "mulreal_crs_ccs_ccs end!"
 call conv_ccs2crs(ccsmat2,crsout)
! write(*,*) "conv_ccs2crs end!"
else if ( iflag .eq. 1) then ! mkl is available 2017.09.06
 sort=8 ! reordering only jc and c
 m1 = crsmat1%nrow
 n1 = crsmat1%ncolm
 k1 = crsmat1%ntot
 m2 = crsmat2%nrow
 n2 = crsmat2%ncolm
 k2 = crsmat2%ntot
 allocate(ia(m1+1),ib(m2+1) )
 allocate(ja(k1  ),jb(k2  ) )
 allocate( a(k1  ), b(k2  ) )
 ia(1:m1+1) = crsmat1%stack(0:m1)+1 ! 2017.12.26
 ja         = crsmat1%item          ! 2017.12.26
  a         = crsmat1%val           ! 2017.12.26
 ib(1:m2+1) = crsmat2%stack(0:m2)+1 ! 2017.12.26
 jb         = crsmat2%item          ! 2017.12.26
  b         = crsmat2%val           ! 2017.12.26
 request=1 ! calculate only the values of ic
 trans="n"
 if ( present(T) ) trans = T ! csrmat1 is used as it is 2017.12.25
 if ( trans .eq. "T" .or. trans .eq. "t" ) then
   m3 = n1 ; n3 = n2  ! 2017.12.26 crsout%nrow = n1
 else
   m3 = m1 ; n3 = n2  ! 2017.12.26 crsout%nrow = m1
 end if
 allocate( ic(m3+1) ) ! 2017.12.26
 ! write(*,*) "start!"
 !call mkl_dcsrmultcsr(trans,request,sort,m1,n1,n2,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,info) ! commented out on 2020.12.04 for gfortran compilation
 !write(*,*) "first end!!"
 request=0
 allocate(c(ic(m3+1)-1),jc(ic(m3+1)-1)) ! 2017.12.26
 nzmax  =   ic(m3+1)-1                  ! 2017.12.26
 !write(*,*) "nzmax=",nzmax
!call mkl_dcsrmultcsr(trans,request,sort,m1,n1,n2,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,info) ! commented out on 2020.12.04 for gfortran compilation
 !write(*,*) "end!!"
 ! output
 crsout%ntot  = ic(m3+1)-1 ! 2017.12.26
 crsout%nrow  = m3 ! 2017.12.26
 crsout%ncolm = n3 ! 2017.12.26
 allocate(crsout%stack(0:m3),crsout%item(crsout%ntot),crsout%val(crsout%ntot))
 crsout%stack(0:m3) = ic(1:m3+1) - 1 ! 2017.12.26
 crsout%item        = jc
 crsout%val         = c
end if

return
end subroutine
!###############################################
! Coded on 2017.05.11
! multiply crs by ccs to gen crs
subroutine mulreal_crs_ccs_crs(crsmat,ccsmat,crsout)
implicit none
type(real_crs_matrix),intent(in)  :: crsmat
type(real_ccs_matrix),intent(in)  :: ccsmat
type(real_crs_matrix),intent(out) :: crsout
type(real_ccs_matrix) :: ccsout

call mulreal_crs_ccs_ccs(crsmat,ccsmat,ccsout)
call conv_ccs2crs(ccsout,crsout)

return
end subroutine
!###############################################
! Coded on 2017.05.11
! multiply crs by ccs to gen ccs
subroutine mulreal_crs_ccs_ccs(crsmat,ccsmat,ccsout)
implicit none
type(real_crs_matrix),intent(in)  :: crsmat
type(real_ccs_matrix),intent(in)  :: ccsmat
type(real_ccs_matrix),intent(out) :: ccsout
integer(4) :: nrowcrs,ncolmcrs,ntotcrs,nrowccs,ncolmccs,ntotccs
integer(4) :: i,j,k,icount,l,TID, OMP_GET_THREAD_NUM
real(8),allocatable,dimension(:)  :: v,vout
integer(4),allocatable,dimension(:) :: stackccs, itemccs, stackcrs,itemcrs
integer(4),allocatable,dimension(:) :: itemwork,stackccsout
real(8),   allocatable,dimension(:) :: valccs,valwork,valcrs
real(8),parameter,dimension(8)  :: maxscale = (/0.01,0.05,0.1,0.2,0.3,0.5,0.7,1.0/)
integer(4) :: ntotoutmax,ntotout
integer(4) :: chunk ! 2017.06.18

!#[1]## set
!# crs
 nrowcrs     = crsmat%nrow
 ncolmcrs    = crsmat%ncolm
 ntotcrs  = crsmat%ntot
 allocate(stackcrs(0:nrowcrs), itemcrs(ntotcrs),valcrs(ntotcrs))
 stackcrs = crsmat%stack
 itemcrs  = crsmat%item
 valcrs   = crsmat%val
!# ccs
 nrowccs  = ccsmat%nrow
 ncolmccs = ccsmat%ncolm
 ntotccs  = ccsmat%ntot
 allocate(stackccs(0:ncolmccs),itemccs(ntotccs),valccs(ntotccs))
 stackccs = ccsmat%stack
 itemccs  = ccsmat%item
 valccs   = ccsmat%val
 if ( ncolmcrs .ne. nrowccs) then
  write(*,*) "GEGEGE! ncolmcrs",ncolmcrs,".ne. nrowccs",nrowccs
  stop
 end if

!#[2]## multiply
 allocate(v(nrowccs),vout(nrowcrs))
 allocate(stackccsout(0:ncolmccs))
 stackccsout = 0
 k=2
 100 continue
 k=k+1
  if ( k .gt. 8 ) then
   write(*,*) "k=",k,"GEGEGE"
   stop
  end if
  ! not that ntotumax should less than 2147483647
  ntotoutmax = (nrowcrs * maxscale(k)) * ncolmccs
!  write(*,10) "k,nrowcrs,ncolmccs,maxxcale,ntotmax",nrowcrs,ncolmccs,maxscale(k),ntotoutmax
!10 format(a,2i8,g15.7,i10)
  if (ntotoutmax .eq. 0 ) goto 100
  allocate(valwork(ntotoutmax),itemwork(ntotoutmax))
  ntotout = 0
 do i=1,ncolmccs
!  if ( mod(i,1000) .eq. 0 ) write(*,*) "i=",i,"ncolmccs",ncolmccs
  v = 0.d0 ; vout = 0.d0
!  write(*,*) "i=",i,"ncolm=",ncolm
  if (stackccs(i-1) .lt. stackccs(i)) then
  do j=stackccs(i-1)+1,stackccs(i)
   v(itemccs(j)) = valccs(j)
  end do
chunk = 1000
!$OMP PARALLEL SHARED(nrowcrs,crsmat,vout,v,CHUNK) PRIVATE(j,l,TID) ! 2017.06.18
!!  TID = OMP_GET_THREAD_NUM()
!$OMP DO SCHEDULE (STATIC,chunk)
  do j=1,nrowcrs
!   write(*,*) "TID",TID,"j",j
   do l=crsmat%stack(j-1)+1,crsmat%stack(j)
     vout(j) = vout(j) + crsmat%val(l)*v(crsmat%item(l))
   end do
  end do
!$OMP END PARALLEL

!  call mul_matcrs_rv(crsmat,v,nrowccs,vout) ! v [nrowccs], vout[nrowcrs]
  do j=1,nrowcrs
   if (abs(vout(j)) .gt. 1.d-15) then
    ntotout = ntotout + 1
!   write(*,*) "i,j=",i,j,"vout=",vout(j),"ntotout=",ntotout
    if (ntotout .gt. ntotoutmax) then
     deallocate(valwork,itemwork)
     goto 100
    end if
    valwork(ntotout)  = vout(j)
    itemwork(ntotout) = j
   end if
  end do
  end if
  stackccsout(i) = ntotout
 end do

!#[3]## output
 allocate(ccsout%val(ntotout),ccsout%item(ntotout),ccsout%stack(0:ncolmccs))
 ccsout%ntot  = ntotout
 ccsout%nrow  = nrowcrs
 ccsout%ncolm = ncolmccs
 ccsout%stack = stackccsout
 ccsout%val   = valwork(1:ntotout)
 ccsout%item  = itemwork(1:ntotout)

return
end subroutine
!###############################################
! Coded on 2017.05.11
! convert compressed row matrix to compressed column matrix
subroutine conv_ccs2crs(ccsmat,crsmat)
implicit none
type(real_ccs_matrix),intent(in)  :: ccsmat
type(real_crs_matrix),intent(out) :: crsmat
integer(4) :: nrow, ncolm, ntot
integer(4) :: i,j,ii,irow,itot
integer(4),allocatable,dimension(:) :: row,colm,stackcrs,itemcrs
integer(4),allocatable,dimension(:) :: stackccs,itemccs,count_colm
real(8),   allocatable,dimension(:) :: valcrs,valccs

!#[0]## set
ncolm = ccsmat%ncolm
nrow  = ccsmat%nrow
ntot  = ccsmat%ntot
allocate(stackccs(0:ncolm),valccs(ntot),itemccs(ntot))
stackccs = ccsmat%stack
itemccs  = ccsmat%item
valccs   = ccsmat%val

!#[1]## set row and colm
allocate(row(ntot),colm(ntot))
ii=0
do i=1,ncolm
 do j=stackccs(i-1)+1,stackccs(i)
  colm(j)   = i
  row(j)    = itemccs(j)
 end do
end do

!#[2]## gen crs matrix
allocate(stackcrs(0:nrow),count_colm(nrow),itemcrs(ntot),valcrs(ntot))
count_colm = 0
do ii=1,ntot
count_colm(row(ii)) = count_colm(row(ii)) + 1
end do
!#[2-1]## stackcrs
stackcrs(0:nrow)=0
do i=1,nrow
 stackcrs(i)=stackcrs(i-1)+count_colm(i)
end do
!#[2-2]## itemcrs, valcrs
count_colm(:)=0
do ii=1,ntot
 irow = row(ii)
 count_colm(irow) = count_colm(irow) + 1
 itot = stackcrs(irow-1)+count_colm(irow)
 itemcrs(itot) = colm(ii)
 valcrs(itot)  = valccs(ii)
end do

!#[3]## set ccs matrix
allocate(crsmat%val(ntot),crsmat%stack(0:nrow),crsmat%item(ntot))
crsmat%stack = stackcrs
crsmat%val   = valcrs
crsmat%item  = itemcrs
crsmat%ntot  = ntot
crsmat%nrow  = nrow
crsmat%ncolm = ncolm

return
end subroutine
!###############################################
! Coded on 2017.05.11
! convert compressed row matrix to compressed column matrix
subroutine convcomp_ccs2crs(ccsmat,crsmat)
implicit none
type(complex_ccs_matrix),intent(in) :: ccsmat
type(complex_crs_matrix),intent(out) :: crsmat
integer(4) :: nrow, ncolm, ntot
integer(4) :: i,j,ii,irow,itot
integer(4),allocatable,dimension(:) :: row,colm,stackcrs,itemcrs
integer(4),allocatable,dimension(:) :: stackccs,itemccs,count_colm
complex(8),allocatable,dimension(:) :: valcrs,valccs

!#[0]## set
ncolm = ccsmat%ncolm
nrow  = ccsmat%nrow
ntot  = ccsmat%ntot
allocate(stackccs(0:ncolm),valccs(ntot),itemccs(ntot))
stackccs = ccsmat%stack
itemccs  = ccsmat%item
valccs   = ccsmat%val

!#[1]## set row and colm
allocate(row(ntot),colm(ntot))
ii=0
do i=1,ncolm
 do j=stackccs(i-1)+1,stackccs(i)
  colm(j)   = i
  row(j)    = itemccs(j)
 end do
end do

!#[2]## gen crs matrix
allocate(stackcrs(0:nrow),count_colm(nrow),itemcrs(ntot),valcrs(ntot))
count_colm = 0
do ii=1,ntot
count_colm(row(ii)) = count_colm(row(ii)) + 1
end do
!#[2-1]## stackcrs
stackcrs(0:nrow)=0
do i=1,nrow
 stackcrs(i)=stackcrs(i-1)+count_colm(i)
end do
!#[2-2]## itemcrs, valcrs
count_colm(:)=0
do ii=1,ntot
 irow = row(ii)
 count_colm(irow) = count_colm(irow) + 1
 itot = stackcrs(irow-1)+count_colm(irow)
 itemcrs(itot) = colm(ii)
 valcrs(itot)  = valccs(ii)
end do

!#[3]## set ccs matrix
allocate(crsmat%val(ntot),crsmat%stack(0:nrow),crsmat%item(ntot))
crsmat%stack = stackcrs
crsmat%val   = valcrs
crsmat%item  = itemcrs
crsmat%ntot  = ntot
crsmat%nrow  = nrow
crsmat%ncolm = ncolm

return
end subroutine

!###############################################
! Coded on 2017.05.11
! convert compressed row matrix to compressed column matrix
subroutine conv_crs2ccs(crsmat,ccsmat)
implicit none
type(real_crs_matrix),intent(in) :: crsmat
type(real_ccs_matrix),intent(out) :: ccsmat
integer(4) :: nrow, ncolm, ntot
integer(4) :: i,j,ii,icolm,itot
integer(4),allocatable,dimension(:) :: row,colm,stackcrs,itemcrs
integer(4),allocatable,dimension(:) :: stackccs,itemccs,count_row
real(8),allocatable,dimension(:)    :: valcrs,valccs

!#[0]## set
 nrow  = crsmat%nrow
 ncolm = crsmat%ncolm
 ntot = crsmat%ntot
 allocate(stackcrs(0:nrow),valcrs(ntot),itemcrs(ntot))
 stackcrs = crsmat%stack
 itemcrs  = crsmat%item
 valcrs   = crsmat%val

!#[1]## set row and colm
 allocate(row(ntot),colm(ntot))
 ii=0
 do i=1,nrow
  do j=stackcrs(i-1)+1,stackcrs(i)
   row(j)  = i
   colm(j) = itemcrs(j)
  end do
 end do

!#[2]## gen ccs matrix
 allocate(stackccs(0:ncolm),count_row(ncolm),itemccs(ntot),valccs(ntot))
 count_row = 0
 do ii=1,ntot
  count_row(colm(ii)) = count_row(colm(ii)) + 1
 end do
 !#[2-1]## stackccs
 stackccs(0:ncolm)=0
 do i=1,ncolm
  stackccs(i)=stackccs(i-1)+count_row(i)
 end do
 !#[2-2]## itemccs, valccs
 count_row(:)=0
 do ii=1,ntot
  icolm = colm(ii)
  count_row(icolm) = count_row(icolm) + 1
  itot = stackccs(icolm-1)+count_row(icolm)
  itemccs(itot) = row(ii)
  valccs(itot)  = valcrs(ii)
 end do

!#[3]## set ccs matrix
 allocate(ccsmat%val(ntot),ccsmat%stack(0:ncolm),ccsmat%item(ntot))
 ccsmat%stack = stackccs
 ccsmat%val   = valccs
 ccsmat%item  = itemccs
 ccsmat%ntot  = ntot
 ccsmat%ncolm = ncolm
 ccsmat%nrow  = nrow

return
end subroutine
!#############################################  transpose crs -> crs
! Coded on 2018.01.23
!# transposed crs1 is provided as crs2
subroutine trans_crs2crs(crs1,crs2)
implicit none
type(real_crs_matrix),intent(in)  :: crs1
type(real_crs_matrix),intent(out) :: crs2
type(real_ccs_matrix) :: ccs1

call trans_crs2ccs(crs1,ccs1)
call conv_ccs2crs(ccs1,crs2)

return
end subroutine
!#############################################  transpose
! debugged on 2017.05.17
! Coded    on 2017.05.12
subroutine trans_crs2ccs(crs,ccs)
implicit none
type(real_crs_matrix),intent(in) :: crs
type(real_ccs_matrix),intent(out) :: ccs
integer(4) :: nrow,ncolm,ntot

!#[1]## set
nrow  = crs%nrow
ncolm = crs%ncolm
ntot  = crs%ntot

!#[2]## gen transverse ccs
ccs%nrow  = ncolm
ccs%ncolm = nrow
ccs%ntot  = ntot
allocate(ccs%stack(0:ccs%ncolm),ccs%item(ntot),ccs%val(ntot))
ccs%stack = crs%stack
ccs%item  = crs%item
ccs%val   = crs%val

return
end subroutine
!#############################################  transpose
! Coded on 2017.05.15
subroutine transcomp_crs2ccs(crs,ccs)
implicit none
type(complex_crs_matrix),intent(in) :: crs
type(complex_ccs_matrix),intent(out) :: ccs
integer(4) :: nrow,ncolm,ntot

!#[1]## set
nrow  = crs%nrow
ncolm = crs%ncolm
ntot  = crs%ntot

!#[2]## gen transverse ccs
ccs%nrow  = ncolm
ccs%ncolm = nrow
ccs%ntot  = ntot
allocate(ccs%stack(0:ccs%ncolm),ccs%item(ntot),ccs%val(ntot))
ccs%stack = crs%stack
ccs%item  = crs%item
ccs%val   = crs%val

return
end subroutine
!############################################# duplicate_crsmat
!# coded on 2017.05.15
subroutine duplicate_crsmat(crsin,crsout)
implicit none
type(real_crs_matrix),intent(in)  :: crsin
type(real_crs_matrix),intent(out) :: crsout
integer(4) :: nrow, ncolm,ntot

nrow  = crsin%nrow
ncolm = crsin%ncolm
ntot  = crsin%ntot
allocate(crsout%stack(0:nrow),crsout%item(ntot),crsout%val(ntot))
crsout%nrow  = nrow
crsout%ncolm = ncolm
crsout%ntot  = ntot
crsout%stack = crsin%stack
crsout%item  = crsin%item
crsout%val   = crsin%val

return
end subroutine
!#############################################  mul_atb
subroutine  mul_atb( a, b, c,  aj_leng, bi_leng, cj_leng  )
!    multiply matrix 'c' = 'a't * 'b'
!      c(aj_leng,cj_leng) = a(bi_leng,aj_leng) * b(bi_leng,cj_leng)
implicit    none
integer(4) , intent(in)  ::  aj_leng !  2-nd element length of array 'a'
integer(4) , intent(in)  ::  bi_leng  ! 1-st element length of array 'b'
integer(4) , intent(in)  ::  cj_leng  !  2-nd element length of array 'c'
real (8), dimension(bi_leng, aj_leng), intent(in)  ::  a! input matrix 'a'
real (8), dimension(bi_leng, cj_leng), intent(in)  ::  b !   input matrix 'b'
real (8), dimension(aj_leng, cj_leng), intent(out) ::  c   !  input matrix 'c'
integer(4 )  ::  i, j, k
do i=1,aj_leng
  do j=1,cj_leng
    c(i,j) = 0.d0
  do k=1,bi_leng
  c(i,j) = c(i,j) + a(k,i)*b(k,j)
enddo
enddo
enddo
return
end subroutine  mul_atb!

!###############################################  mul_ab
subroutine  mul_ab ( a, b, c,  ai_leng, bi_leng, cj_leng  )
!    multiply matrix 'c' = 'a' * 'b'
!      c(ai_leng,cj_leng) = a(ai_leng,bi_leng) * b(bi_leng,cj_leng)
implicit    none
integer(4), intent(in)  ::  ai_leng ! 1-st element length of array 'a'
integer(4), intent(in)  ::  bi_leng ! 1-st element length of array 'b'
integer(4), intent(in)  ::  cj_leng !  2-nd element length of array 'c'
real (8), dimension(ai_leng, bi_leng), intent(in)  ::  a !   input matrix 'a'
real (8), dimension(bi_leng, cj_leng), intent(in)  ::  b  !  input matrix 'b'
real (8), dimension(ai_leng, cj_leng), intent(out) ::  c ! input matrix 'c'
integer(4 )  ::  i, j, k
do i=1,ai_leng
  do j=1,cj_leng
    c(i,j) = 0.d0
    do k=1,bi_leng
      c(i,j) = c(i,j) + a(i,k)*b(k,j)
    enddo
  enddo
enddo
return
end subroutine  mul_ab
!############################################### mul_abt
subroutine  mul_abt( a, b, c, ai_leng, bj_leng, cj_leng  )
!    multiply matrix 'c' = 'a' * 'b't
!      c(ai_leng,cj_leng) = a(ai_leng,bj_leng) * b(cj_leng,bj_leng)
implicit    none
integer(4), intent(in)  ::  ai_leng !  1-st element length of array 'a'
integer(4), intent(in)  ::  bj_leng !  2-nd element length of array 'b'
integer(4), intent(in)  ::  cj_leng !  2-nd element length of array 'c'
real   (8), dimension(ai_leng, bj_leng), intent(in)  ::  a!   input matrix 'a'
real   (8), dimension(cj_leng, bj_leng), intent(in)  ::  b !  input matrix 'b'
real   (8), dimension(ai_leng, cj_leng), intent(out) ::  c!  input matrix 'c'
integer(4 )  ::  i, j, k
do i=1,ai_leng
  do j=1,cj_leng
    c(i,j) = 0.d0
    do k=1,bj_leng
      c(i,j) = c(i,j) + a(i,k)*b(j,k)
    enddo
  enddo
enddo
return
end subroutine  mul_abt
!############################################### realcrsout
!# coded on 2017.05.17
subroutine realcrsout(crsout)
implicit none
type(real_crs_matrix),intent(in) :: crsout

write(*,*) "size(item)=",size(crsout%item)
write(*,*) "size(stack)=",size(crsout%stack)
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
write(*,*) "crs%stack=",crsout%stack
write(*,*) "crs%item=",crsout%item
write(*,*) "crs%val=",crsout%val

return
end subroutine
!############################################### compcrsout
!# coded on 2017.05.17
subroutine compcrsout(crsout)
implicit none
type(complex_crs_matrix),intent(in) :: crsout

write(*,*) "size(item)=",size(crsout%item)
write(*,*) "size(stack)=",size(crsout%stack)
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
write(*,*) "crs%stack=",crsout%stack
write(*,*) "crs%item=",crsout%item
write(*,*) "crs%val=",crsout%val

return
end subroutine

!############################################### dewallocate real_crsmatrix
!# 2017.06.12
subroutine deallocate_real_crsmat(crsmat)
implicit none
type(real_crs_matrix),intent(inout) :: crsmat

crsmat%ntot = 0.d0
crsmat%nrow = 0.d0
crsmat%ncolm = 0.d0
if (allocated(crsmat%stack)) deallocate(crsmat%stack)
if (allocated(crsmat%item) ) deallocate(crsmat%item)
if (allocated(crsmat%val)  ) deallocate(crsmat%val)

return
end subroutine

end module matrix
