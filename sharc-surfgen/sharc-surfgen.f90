!============================================================================
program main
  use nndms
  implicit none
  real*8, parameter :: bohr2angs=0.52917721067d0
  real*8, parameter :: degree=acos(-1.d0)/180.d0
  real*8, allocatable :: geom(:),veloc(:),e(:),cg(:,:,:),h(:,:),dcg(:,:,:)
  real*8, allocatable :: dip(:,:,:)
  integer :: i,j,trajid,ndim
  character(len=3) :: nm
  real*8  :: z,cavity
  real*8  :: ht(1),dht(12,1,1)
!  call pesinit
  call dmsinit
  natoms=4
  nstates=2
  allocate(geom(3*natoms),veloc(3*natoms))
  allocate(e(nstates),h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))
  allocate(dip(3,nstates,nstates))

  open(100,file='QM.in',position="rewind")
  open(101,file='QM.out')
!read QM.in
  read(100,*)natoms
  read(100,*)trajid
  do i=1,natoms
    read(100,*)nm,geom((3*i-2):(3*i)),veloc((3*i-2):(3*i))
  enddo
!note the QM.in is in angstrom, need au for surfgen
  geom=geom/bohr2angs
!  print *,geom
!call surface
!  call Evaluatehf(geom,e,cg,h,dcg)
  call Evaluatesurfgen(geom,e,cg,h,dcg)
  call EvaluateDipole(geom,dip)
  cavity=0.1*dip(3,1,2)
!Hamiltonian
  write(101,"(A)")"! 1 Hamiltonian Matrix (2x2, complex)"
  write(101,"(2i2)")nstates,nstates
  z=0d0
  write(101,"(4e14.6)")e(1),z,z,z
!  write(101,"(4e14.6)")e(1),z,cavity,z
!  write(101,"(4e14.6)")cavity,z,e(2),z
  write(101,"(4e14.6)")z,z,e(2),z
  write(101,"(A)")
!dipole
  ndim=3
  write(101,"(A)")"! 2 Dipole Moment Matrices (3x2x2, complex)"
  write(101,"(2i2)")nstates,nstates
  write(101,"(4e14.6)")dip(1,1,1),z,dip(1,1,2),z
  write(101,"(4e14.6)")dip(1,2,1),z,dip(1,2,2),z
  write(101,"(2i2)")nstates,nstates
  write(101,"(4e14.6)")dip(2,1,1),z,dip(2,1,2),z
  write(101,"(4e14.6)")dip(2,2,1),z,dip(2,2,2),z
  write(101,"(2i2)")nstates,nstates
  write(101,"(4e14.6)")dip(3,1,1),z,dip(3,1,2),z
  write(101,"(4e14.6)")dip(3,2,1),z,dip(3,2,2),z
  write(101,"(A)")
!Gradient
  write(101,"(A)")"! 3 Gradient Vectors (2x4x3, real)"
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),2,2)
  enddo
  write(101,"(A)")
!Nonadiabatic coupling
  write(101,"(A)")"! 5 nacdr Vectors (2x2x4x3, real)"
!1
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0   m2 1 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0   m2 1 s2 2 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),1,2)
  enddo
!2 
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0   m2 1 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),2,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0   m2 1 s2 2 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo

end
!============================================================================
