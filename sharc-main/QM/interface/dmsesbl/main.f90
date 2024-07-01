!============================================================================
program main
  use nndms
  implicit none
  real*8, allocatable :: gmek(:,:)
  real*8, allocatable :: adiadpl(:,:,:)
  integer :: i,j,k,m,o,total,ios,idx,info,id
  character(3) :: sym
  real*8 :: anums,dipole(3)

  !init DMS
  call dmsinit

  allocate(gmek(3,natoms))
  allocate(adiadpl(3,nstates,nstates))

  open(unit=100,file='xyz.txt',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(unit=101,file='dipole1.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(200,file='dipole1.txt')
  do i=1,1000
    read(100,*,iostat=ios) anums
    read(100,*,iostat=ios) anums
    do j=1,natoms
      read(100,*,iostat=ios) sym,gmek(:,j)
    end do
    read(101,*,iostat=ios) dipole
    if(ios.ne.0) then
      total=i-1
      exit
    end if
    call EvaluateDipole(gmek,adiadpl)
    write(200,"(7f15.9)") dble(i+79),dipole,adiadpl(:,1,1)
  end do

  print*, total, 'valid points.'

  stop
end
!============================================================================
