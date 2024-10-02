!==========================================================================
program main
  implicit none
  real*8, parameter :: au2cm=219474d0
  real*8, parameter :: bohr2angs=0.52917721067d0
  character(3) :: sym
  real*8 :: anums
  real*8, dimension(:), allocatable :: e,geom,veloc
  real*8, dimension(:,:), allocatable :: h,gref
  real*8, dimension(:,:,:), allocatable :: cg, dcg
  real*8 :: r,phi
  real*8 :: eshift
  integer :: i,j,k,total,ios,idx,natoms,nstates,trajid
  call initPotential()
  call getinfo(natoms, nstates)
  total=500
  allocate(e(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))
  allocate(geom(3*natoms),veloc(3*natoms))

  open(100,file='QM.in',position="rewind")
  !read QM.in
  read(100,*)natoms
  read(100,*)trajid
  do i=1,natoms
    read(100,*)sym,geom((3*i-2):(3*i)),veloc((3*i-2):(3*i))
  enddo
  !QM.in is in angstrom, need au for surfgen
  geom=geom/bohr2angs

  open(200,file='e.out')
  call EvaluateSurfgen(geom,e,cg,h,dcg)
  write(200,"(2f18.8)") e(:)
  write(200,"(12f18.8)") cg(:,1,1)
  write(200,"(12f18.8)") cg(:,2,2)
  write(200,"(12f18.8)") cg(:,1,2)
  write(200,"(12f18.8)") cg(:,2,1)

  stop
end
!==========================================================================
