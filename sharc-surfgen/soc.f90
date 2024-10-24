!===================================================================================
!Module for NN DMS
module nnsoc
  use nn_class
  implicit none
  integer :: natoms
  type(nn), allocatable :: ANN(:) !Feed-Forward Neural Network
  integer :: na ! number of NNs
  integer :: ncoord   !total number of internal coords
  integer :: nstates  !number of electronic states
  character(len=99), dimension(:),  allocatable ::  wfile
  real*8 :: bonddiff
end module nnsoc
!===============================================================================
!initialize NN DMS
subroutine socinit
  use nnsoc
  implicit none
  integer :: i

  natoms=4
  nstates=2
  ncoord=7
  na=1
  bonddiff=5.d-2
  call initPotential()

  allocate(ANN(na),wfile(na))

  wfile(1)='WBsoc.txt'

  do i=1,na
    ANN(i)%structfl='struct-soc'
    call ANN(i)%init()
    call ANN(i)%savenn(wfile(i),0)
  end do

  return
end subroutine socinit
!==================================================================================
subroutine Evaluatedsoc(i,igeom,soc)
  use nnsoc
  implicit none
  integer, intent(in) :: i
  real*8, intent(in) :: igeom(ncoord)
  real*8, intent(out) :: soc(3,nstates)
  integer :: R
  real*8  :: z

  R=ANN(i)%RX
  call ANN(i)%output(igeom(1:R))
  !(1,1)
  soc(1,1)=ANN(i)%y(1)
  soc(2,1)=ANN(i)%y(2)
  soc(3,1)=ANN(i)%y(3)*igeom(ncoord)
  !(1,2)
  soc(1,2)=ANN(i)%y(4)*igeom(ncoord)
  soc(2,2)=ANN(i)%y(5)*igeom(ncoord)
  soc(3,2)=ANN(i)%y(6)
  z=10d0
  soc=soc*z

  return
end subroutine Evaluatedsoc
!==================================================================================
subroutine Evaluateasoc(cgeom,asoc)
  use nnsoc
  implicit none
  integer, parameter :: lwork=99
  real*8, intent(in) :: cgeom(3,natoms)
  real*8, intent(out) :: asoc(2,3,nstates)
  real*8, allocatable :: gmek(:,:),coord(:),work(:)
  real*8, dimension(:), allocatable :: e
  real*8, dimension(:,:), allocatable :: ckl
  real*8, dimension(:,:,:), allocatable :: cg, dcg
  real*8 :: T(3,3),dsoc(3,2),dsoctmp(3,2),theta,trial(2),diff1,diff2,h11(2),h22(2),err
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

  !diabatic soc
  call dplcoord(natoms,ncoord,gmek,coord)
  call Evaluatedsoc(1,coord,dsoc)
  !rotation
  do i=1,nstates
    dsoc(:,i)=matmul(transpose(T),dsoc(:,i))
  end do
  !get diabatization rotation angle theta=1/2atan(2h12/(h22-h11))
  call EvaluateSurfgen(gmek,e,cg,ckl,dcg)
  trial(1)=0.5d0*atan(-2*ckl(1,2)/(ckl(2,2)-ckl(1,1)))
  trial(2)=trial(1)+0.5d0*3.14159265358979
!select trn dipole sign and k
  h11(1)=cos(trial(1))**2*e(1)+sin(trial(1))**2*e(2)
  h11(2)=cos(trial(2))**2*e(1)+sin(trial(2))**2*e(2)
  h22(1)=sin(trial(1))**2*e(1)+cos(trial(1))**2*e(2)
  h22(2)=sin(trial(2))**2*e(1)+cos(trial(2))**2*e(2)
  diff1=abs(h11(1)-ckl(1,1))+abs(h22(1)-ckl(2,2))
  diff2=abs(h11(2)-ckl(1,1))+abs(h22(2)-ckl(2,2))
  if (diff1.lt.diff2) then
    err=diff1
    theta=trial(1)
  else
    err=diff2
    theta=trial(2)
  endif
  if(err.gt.1d-4) then
    print *,"error determine theta stop program"
    stop
  endif

  !diabatic rotation
  dsoctmp=dsoc
  dsoc(:,1)=dsoctmp(:,1)*cos(theta)+dsoctmp(:,2)*sin(theta)
  dsoc(:,2)=-dsoctmp(:,1)*sin(theta)+dsoctmp(:,2)*cos(theta)

  !fill adiabatic soc matrix elements
  !asoc(x,y,z) x=1,2 is real/imag part
  !y=1,2,3 for lx, ly, lz
  !z=1,2 for s0/t1 and s1/t1 soc
  !see Fedorov, D. G.; Gordon, M. S. (2002) Chapter 13, pp276–297.
  do i=1,2
    asoc(1,1,i)=dsoc(2,i)/1.4143
    asoc(2,1,i)=-dsoc(1,i)/1.4143
    asoc(1,2,i)=0
    asoc(2,2,i)=dsoc(3,i)/1.4143
    asoc(1,3,i)=asoc(1,1,i)
    asoc(2,3,i)=-asoc(2,1,i)
  enddo
  return
end subroutine evaluateasoc
