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
end module nnsoc
!===================================================================================
!initialize NN DMS
subroutine socinit
  use nnsoc
  implicit none
  integer :: i

  natoms=4
  nstates=2
  ncoord=7
  na=1

  allocate(ANN(na),wfile(na))

  wfile(1)='WBsoc.txt'

  do i=1,na
    ANN(i)%structfl='struct.soc'
    call ANN(i)%init()
    call ANN(i)%savenn(wfile(i),0)
  end do

  return
end subroutine socinit
!==================================================================================
subroutine EvaluateSoc(i,igeom,soc)
  use nndms
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
end subroutine EvaluateSoc
!==================================================================================
subroutine Evaluatediasoc(cgeom,soc)
  use nndms
  implicit none
  integer, parameter :: lwork=99
  real*8, intent(in) :: cgeom(3,natoms)
  real*8, intent(out) :: soc(3,nstates)
  real*8, allocatable :: gmek(:,:),coord(:),work(:)
  real*8 :: T(3,3)
  integer :: i,j,info

  allocate(gmek(3,natoms),coord(ncoord),work(lwork))

  gmek=cgeom
  call order(gmek,info)
  call orientation(gmek,T,info)
  !diabatic dipole
  call dplcoord(natoms,ncoord,gmek,coord)
  call EvaluateSoc(1,coord,soc)

  !rotation
  do i=1,nstates
    soc(:,i)=matmul(transpose(T),soc(:,i))
  end do

  return
end
!=================================================================================
