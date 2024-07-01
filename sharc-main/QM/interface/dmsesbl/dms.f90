!===================================================================================
!Module for NN DMS
module nndms
  use nn_class
  implicit none
  integer :: natoms
  type(nn), allocatable :: ANN(:) !Feed-Forward Neural Network
  integer :: na ! number of NNs
  integer :: ncoord   !total number of internal coords
  integer :: nstates  !number of electronic states
  character(len=99), dimension(:),  allocatable ::  wfile
  real*8 :: bonddiff
end module nndms
!===================================================================================
subroutine dplcoord(natom,ncoords,cgeom,igeom)
  implicit none
  integer, intent(in) :: natom, ncoords
  real*8, intent(in) :: cgeom(3,natom)
  real*8, intent(out) :: igeom(ncoords)
  real*8 :: dx(3),dcoord(12),x(3,4)
  integer :: i,j,k,info

  x=cgeom(1:3,1:4)
  call order(x,info)

  k=0
  do i=1,natom-1
    do j=i+1,natom
      k=k+1
      dx=cgeom(:,i)-cgeom(:,j)
      igeom(k)=1.d0/sqrt(dot_product(dx,dx))
    end do
  end do

  call OOP0(cgeom(:,1),cgeom(:,2),cgeom(:,3),cgeom(:,4),igeom(ncoords),dcoord,0)

  return
end
!===============================================================================
!initialize NN DMS
subroutine dmsinit
  use nndms
  implicit none
  integer :: i

  natoms=4
  nstates=2
  ncoord=7
  na=10
  bonddiff=5.d-2
  call initPotential()

  allocate(ANN(na),wfile(na))

  wfile(1)='WB01.txt'
  wfile(2)='WB02.txt'
  wfile(3)='WB03.txt'
  wfile(4)='WB04.txt'
  wfile(5)='WB05.txt'
  wfile(6)='WB06.txt'
  wfile(7)='WB07.txt'
  wfile(8)='WB08.txt'
  wfile(9)='WB09.txt'
  wfile(10)='WB10.txt'

  do i=1,na
    ANN(i)%structfl='struct'
    call ANN(i)%init()
    call ANN(i)%savenn(wfile(i),0)
  end do

  return
end subroutine dmsinit
!==================================================================================
subroutine EvaluateDpl0(i,igeom,dpl)
  use nndms
  implicit none
  integer, intent(in) :: i
  real*8, intent(in) :: igeom(ncoord)
  real*8, intent(out) :: dpl(3,nstates,nstates)
  integer :: R

  R=ANN(i)%RX
  call ANN(i)%output(igeom(1:R))

  !(1,1)
  dpl(1,1,1)=ANN(i)%y(1)
  dpl(2,1,1)=ANN(i)%y(2)
  dpl(3,1,1)=ANN(i)%y(3)*igeom(ncoord)
  !(1,2)
  dpl(1,1,2)=ANN(i)%y(4)*igeom(ncoord)
  dpl(2,1,2)=ANN(i)%y(5)*igeom(ncoord)
  dpl(3,1,2)=ANN(i)%y(6)
  !(2,1)
  dpl(:,2,1)=dpl(:,1,2)
  !(2,2)
  dpl(1,2,2)=ANN(i)%y(7)
  dpl(2,2,2)=ANN(i)%y(8)
  dpl(3,2,2)=ANN(i)%y(9)*igeom(ncoord)

  return
end subroutine EvaluateDpl0
!==================================================================================
subroutine EvaluateDpl(igeom,dpl)
  use nndms
  implicit none
  real*8, intent(in) :: igeom(ncoord)
  real*8, intent(out) :: dpl(3,nstates,nstates)
  real*8 :: dpl0(3,nstates,nstates)
  integer :: i

  dpl=0.d0
  do i=1,na
    call EvaluateDpl0(i,igeom,dpl0)
    dpl=dpl+dpl0
  end do
  dpl=dpl/dble(na)

  return
end subroutine EvaluateDpl
!==================================================================================
subroutine EvaluateDipole(cgeom,dpl)
  use nndms
  implicit none
  integer, parameter :: lwork=99
  real*8, intent(in) :: cgeom(3,natoms)
  real*8, intent(out) :: dpl(3,nstates,nstates)
  real*8, allocatable :: gmek(:,:),coord(:),work(:)
  real*8, dimension(:), allocatable :: e
  real*8, dimension(:,:), allocatable :: ckl
  real*8, dimension(:,:,:), allocatable :: cg, dcg
  real*8 :: T(3,3)
  integer :: i,j,info

  allocate(gmek(3,natoms),coord(ncoord),work(lwork))
  allocate(ckl(nstates,nstates))
  allocate(e(nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  !standard orientation
  gmek=cgeom
  call order(gmek,info)
  call orientation(gmek,T,info)

  !AtD transformation
  call EvaluateSurfgen(gmek,e,cg,ckl,dcg)
  call dsyev('V','U',nstates,ckl,nstates,e,work,lwork,info)
  if(info.ne.0) stop 'Failed to solve eigenvectors!'

  !diabatic dipole
  call Evaluatediadpl(gmek,dpl)

  !adiabatic dipole
  do i=1,3
    dpl(i,:,:)=matmul(transpose(ckl),matmul(dpl(i,:,:),ckl))
  end do

  !rotation
  do i=1,nstates
    do j=1,nstates
      dpl(:,i,j)=matmul(transpose(T),dpl(:,i,j))
    end do
  end do

  return
end
!=================================================================================
subroutine checkbond(cgeom,info)
  use nndms, only: natoms,bonddiff
  implicit none
  real*8, intent(in) :: cgeom(3,natoms)
  integer, intent(out) :: info
  real*8, allocatable :: bond(:),gmek(:,:)
  integer :: nbond
  real*8 :: d1,d2

  nbond=natoms*(natoms-1)/2
  allocate(bond(nbond),gmek(3,natoms))

  gmek=cgeom
  call order(gmek,info)
  call cart2dist(natoms,gmek,bond)

  d1=bond(2)-bond(1)
  d2=bond(3)-bond(2)

  if(d1.ge.bonddiff .and. d2.ge.bonddiff) then
    info=1
  else
    info=0
  end if

  return
end
!================================================================================
subroutine Evaluatediadpl(cgeom,dpl)
  use nndms
  implicit none
  integer, parameter :: nmax=999
  real*8, intent(in) :: cgeom(3,natoms) !already at standard orientation
  real*8, intent(out) :: dpl(3,nstates,nstates)
  real*8 :: T(3,3)
  real*8 :: e1(3),e2(3),e3(3),fval(9)
  real*8, allocatable :: gmek(:,:),coord(:)
  integer :: i,j,k,o,p,info,id
  !shepard interpolation
  real*8 :: dr !bond disp
  integer :: n,nq,nr,nw
  integer, allocatable :: lcell(:,:,:),lnext(:)
  real*8, allocatable :: x(:),y(:),z(:),f(:,:),a(:,:),rsq(:)
  real*8 :: xyzmin(3),xyzdel(3),rmax
  real*8, external :: qs3val

  dr=bonddiff
  allocate(gmek(3,natoms),coord(ncoord))
  allocate(x(nmax),y(nmax),z(nmax),f(nmax,9),lnext(nmax),a(9,nmax),rsq(nmax))

  !check bond distances
  call checkbond(cgeom,info)

  if(info .eq. 1) then
    !get diabatic dipole directly from surface
    call dplcoord(natoms,ncoord,cgeom,coord)
    call EvaluateDpl(coord,dpl)
  else !get diabatic dipole through shepard interpolation
    !unit vectors
    e1=cgeom(:,2)-cgeom(:,1)
    e2=cgeom(:,3)-cgeom(:,1)
    e3=cgeom(:,4)-cgeom(:,1)
    e1=e1/sqrt(dot_product(e1,e1))
    e2=e2/sqrt(dot_product(e2,e2))
    e3=e3/sqrt(dot_product(e3,e3))
    n=0
    do i=-2,2
      do j=-2,2
        do k=-2,2
          !disp
          gmek(:,1)=cgeom(:,1)
          gmek(:,2)=cgeom(:,2)+dble(i)*dr*e1
          gmek(:,3)=cgeom(:,3)+dble(j)*dr*e2
          gmek(:,4)=cgeom(:,4)+dble(k)*dr*e3
          call checkbond(gmek,info)
          if(info.eq.1) then
            n=n+1
            x(n)=dble(i)
            y(n)=dble(j)
            z(n)=dble(k)
            call order(gmek,id)
            call orientation(gmek,T,info)
            call dplcoord(natoms,ncoord,gmek,coord)
            call EvaluateDpl(coord,dpl)
            do o=1,nstates
              do p=1,nstates
                dpl(:,o,p)=matmul(transpose(T),dpl(:,o,p))
              end do
            end do
            call invpermdiadpl(nstates,id,dpl)
            f(n,1:3)=dpl(1:3,1,1)
            f(n,4:6)=dpl(1:3,1,2)
            f(n,7:9)=dpl(1:3,2,2)
          end if
        end do
      end do
    end do
    !shepard interpolation
    if(n.gt.nmax .or. n .lt. 10) stop 'n.gt.nmax .or. n .lt. 10!'
    nq=min(n-1,17)
    nw=min(n-1,32)
    nr=int((dble(n)/3.d0)**(1.d0/3.d0))+1
    allocate(lcell(nr,nr,nr))
    do i=1,9
      call qshep3(n,x(1:n),y(1:n),z(1:n),f(1:n,i),nq,nw,nr,lcell,lnext(1:n),xyzmin,&
                  xyzdel,rmax,rsq(1:n),a,info)
      if(info.ne.0) stop 'shepard interpolation failed!'
      fval(i)=qs3val(0.d0,0.d0,0.d0,n,x(1:n),y(1:n),z(1:n),f(1:n,i),nr,lcell,&
                     lnext,xyzmin,xyzdel,rmax,rsq(1:n),a)
    end do
    dpl(1:3,1,1)=fval(1:3)
    dpl(1:3,1,2)=fval(4:6)
    dpl(1:3,2,1)=fval(4:6)
    dpl(1:3,2,2)=fval(7:9)
  end if

  return
end
!===================================================================================
