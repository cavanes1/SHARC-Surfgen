!============================================================================
program main
  use nndms
  implicit none
  real*8, parameter :: bohr2angs=0.52917721067d0
  real*8, parameter :: au2cm=219474.63
  real*8, parameter :: degree=acos(-1.d0)/180.d0
  real*8, allocatable :: geom(:),veloc(:),e(:),cg(:,:,:),h(:,:),dcg(:,:,:)
  real*8, allocatable :: dip(:,:,:),et(:)
  integer :: i,j,trajid,ndim
  character(len=3) :: nm
  integer  :: nstatestol
  real*8  :: z,cgt(12,1,1),soc(2,3,2),diasoc(3,2)

  call dmsinit
  natoms=4
  nstates=2
  nstatestol=5
  allocate(geom(3*natoms),veloc(3*natoms))
  allocate(e(nstates),et(nstates),h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))
  allocate(dip(3,nstates,nstates))

  open(100,file='QM.in',position="rewind")
  open(101,file='dip.out')
!read QM.in
  read(100,*)natoms
  read(100,*)trajid
  do i=1,natoms
    read(100,*)nm,geom((3*i-2):(3*i)),veloc((3*i-2):(3*i))
  enddo
!note the QM.in is in angstrom, need au for surfgen
  geom=geom/bohr2angs
  call EvaluateDipole(geom,dip)
  write(101,"(3f18.8)")dip(:,1,1)
  write(101,"(3f18.8)")dip(:,2,2)
  write(101,"(3f18.8)")dip(:,1,2)

end
!============================================================================

