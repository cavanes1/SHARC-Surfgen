!============================================================================
program main
  use nnsoc
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

  call socinit
  natoms=4
  nstates=2
  nstatestol=5
  allocate(geom(3*natoms),veloc(3*natoms))
  allocate(e(nstates),et(nstates),h(nstates,nstates))
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
!read the surfgen generated surfaces
  open(201,file="../s0_s1.out")
  open(202,file="../t1.out")
  open(203,file="../dip.out")
!singlet, energy, gradient, derivative coupling
  read(201,*)e(:)
  read(201,*)cg(:,1,1)
  read(201,*)cg(:,2,2)
  read(201,*)cg(:,1,2)
  read(201,*)cg(:,2,1)
!triplet
  read(202,*)et(:)
  read(202,*)cg(:,1,1)
!dipole
  read(203,*)dip(:,1,1)
  read(203,*)dip(:,2,2)
  read(203,*)dip(:,1,2)
  dip(:,2,1)=-dip(:,1,2)
!call surface
  call Evaluateasoc(geom,soc)
!Hamiltonian
  write(101,"(A)")"! 1 Hamiltonian Matrix (5x5, complex)"
  write(101,"(2i2)")nstatestol,nstatestol
  z=0d0
!property surfaces sometimes blow up at region that are rarely visited
!set the bound to be 50cm-1
  do i=1,2
    do j=1,3
      if(soc(1,j,i).gt.50d0 .or. soc(2,j,i).gt.50d0) then
        soc(:,j,i)=50d0
      endif
    enddo
  enddo
  soc=soc/au2cm
  !Hermitian Hamiltonian
  write(101,"(10e14.6)")e(1),z,z,z,soc(:,1,1),soc(:,2,1),soc(:,3,1)
  write(101,"(10e14.6)")z,z,e(2),z,soc(:,1,2),soc(:,2,2),soc(:,3,2)
  write(101,"(10e14.6)")soc(1,1,1),-soc(2,1,1),soc(1,1,2),-soc(2,1,2),et(1),z,z,z,z,z
  write(101,"(10e14.6)")soc(1,2,1),-soc(2,2,1),soc(1,2,2),-soc(2,2,2),z,z,et(1),z,z,z
  write(101,"(10e14.6)")soc(1,3,1),-soc(2,3,1),soc(1,3,2),-soc(2,3,2),z,z,z,z,et(1),z
  write(101,"(A)")
!dipole
  ndim=3
  write(101,"(A)")"! 2 Dipole Moment Matrices (3x5x5, complex)"
  write(101,"(2i2)")nstatestol,nstatestol
  write(101,"(10e14.6)")dip(1,1,1),z,dip(1,1,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(2i2)")nstatestol,nstatestol
  write(101,"(10e14.6)")dip(2,1,1),z,dip(2,1,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(2,2,1),z,dip(2,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(2i2)")nstatestol,nstatestol
  write(101,"(10e14.6)")dip(3,1,1),z,dip(3,1,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(3,2,1),z,dip(3,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(10e14.6)")dip(1,2,1),z,dip(1,2,2),z,z,z,z,z,z,z
  write(101,"(A)")
!Gradient
  write(101,"(A)")"! 3 Gradient Vectors (5x4x3, real)"
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),2,2)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 0"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 1"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(A)")
!Nonadiabatic coupling
  write(101,"(A)")"! 5 nacdr Vectors (5x5x4x3, real)"
!1
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0   m2 1 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0   m2 1 s2 2 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),1,2)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0   m2 3 s2 1 ms2 -1"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0   m2 3 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 1 ms1 0   m2 3 s2 1 ms2 1"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
!
!2
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0   m2 1 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),2,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0   m2 1 s2 2 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cg((3*i-2):(3*i),2,2)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0   m2 3 s2 1 ms2 -1"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0   m2 3 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 1 s1 2 ms1 0   m2 3 s2 1 ms2 1"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
!
!3
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 1 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 1 s2 2 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 -1"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 1"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
!
!4
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 1 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 1 s2 2 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 -1"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 1"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
!
!5
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 1 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 1 s2 2 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 -1"
  do i=1,natoms
    write(101,"(3e20.12)")cgt((3*i-2):(3*i),1,1)
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 0"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
  write(101,"(2i2,x,A)")natoms,ndim,"! m1 3 s1 1 ms1 -1   m2 3 s2 1 ms2 1"
  do i=1,natoms
    write(101,"(3e20.12)")z,z,z
  enddo
!

end
!============================================================================

