!==============================================================================
subroutine permutation(geom,id)
  !atom list: N H H H
  implicit none
  real*8, intent(inout) :: geom(3,4)
  integer, intent(in) :: id
  real*8 :: x(3,4)

  x=geom

  select case(id)
  case(1) !E
    geom=x
  case(2) !(123)
    geom(:,2)=x(:,3)
    geom(:,3)=x(:,4)
    geom(:,4)=x(:,2)
  case(3) !(132)
    geom(:,2)=x(:,4)
    geom(:,3)=x(:,2)
    geom(:,4)=x(:,3)
  case(4) !(12)
    geom(:,2)=x(:,3)
    geom(:,3)=x(:,2)
  case(5) !(23)
    geom(:,3)=x(:,4)
    geom(:,4)=x(:,3)
  case(6) !(13)
    geom(:,2)=x(:,4)
    geom(:,4)=x(:,2)
  case default
    print*, 'Wrong permutation number. No permutation performed!'
    geom=x
  end select

  return
end
!==============================================================================
subroutine order(geom,id)
  implicit none
  real*8, intent(inout) :: geom(3,4)
  integer, intent(out) :: id
  real*8 :: x(3,4),r(6)
  integer :: i

  do i=1,6
    x=geom
    call permutation(x,i)
    call cart2dist(4,x,r)
    if(r(1).le.r(2) .and. r(1).le.r(3) .and. r(2).le.r(3)) then
      id=i
      exit
    end if
  end do

  call permutation(geom,id)

  return
end
!==============================================================================
subroutine setorigin(geom)
  implicit none
  real*8, intent(inout) :: geom(3,4)
  real*8 :: dx(3)
  integer :: i
  dx=geom(:,1)
  do i=1,4
    geom(:,i)=geom(:,i)-dx
  end do
  return
end
!==============================================================================
subroutine rotatex(geom,R,info)
  implicit none
  real*8, intent(inout) :: geom(3,4)
  real*8, intent(out) :: R(3,3)
  integer, intent(out) :: info
  real*8 :: a(3),b(3),n(3),t,norm

  !normalized vectors
  a=geom(:,2)-geom(:,1)
  t=sqrt(dot_product(a,a))
  a=a/t
  b(1)=1.d0;b(2:3)=0.d0

  !rotation angle
  t=acos(dot_product(a,b))

  !rotation axis
  call cross_product(a,b,n)
  norm=sqrt(dot_product(n,n))
  if(norm .lt. 1.d-3) then
    write(*,"('Warning! Small norm of vector in rotatex: ',e15.6)") norm
    info=-1
    return
  end if

  !rotation matrix
  call rotation0(n,t,R)

  !rotation
  geom=matmul(R,geom)
  info=0

  return
end
!===============================================================================
subroutine rotatez(geom,R,info)
  implicit none
  real*8, intent(inout) :: geom(3,4)
  real*8, intent(out) :: R(3,3)
  integer, intent(out) :: info
  real*8 :: t,d,a(3),b(3),c(3)
  integer :: i

  !new z axis
  a(1)=1.d0;a(2:3)=0.d0
  b=geom(:,3)-geom(:,1)
  call cross_product(a,b,c)
  d=sqrt(dot_product(c,c))
  b=c/d
  if(d .lt. 1.d-3) then
    write(*,"('Small norm of vector in rotatez determining new z axis: ',e15.6)") d
    print*, 'Warning: atoms are colinear or nearly colinear.'
    info=-1
    return
  end if
  !old z axis
  a(1:2)=0.d0;a(3)=1.d0

  !rotation angle
  t=-acos(dot_product(a,b))

  !rotation axis
  call cross_product(a,b,c)
  d=sqrt(dot_product(c,c))

  if(d .lt. 1.d-3) then
  write(*,"('Small norm of vector in rotatez determining rotation axis: ',e15.6)") d
  info=-1
  return
  end if

  !rotation matrix
  call rotation0(c,t,R)

  !rotation
  geom=matmul(R,geom)
  info=0

  return
end
!===============================================================================
subroutine rotation0(n,t,R)
  implicit none
  real*8, intent(inout) :: n(3)
  real*8, intent(in) :: t
  real*8, intent(out) :: R(3,3)
  real*8 :: d,sint,cost

  !normalized rotation axis
  d=sqrt(dot_product(n,n))
  n=n/d

  !rotation matrix
  sint=sin(t)
  cost=cos(t)
  R(1,1)=cost+n(1)**2*(1.d0-cost)
  R(1,2)=n(1)*n(2)*(1.d0-cost)-n(3)*sint
  R(1,3)=n(1)*n(3)*(1.d0-cost)+n(2)*sint
  R(2,1)=n(2)*n(1)*(1.d0-cost)+n(3)*sint
  R(2,2)=cost+n(2)**2*(1.d0-cost)
  R(2,3)=n(2)*n(3)*(1.d0-cost)-n(1)*sint
  R(3,1)=n(3)*n(1)*(1.d0-cost)-n(2)*sint
  R(3,2)=n(3)*n(2)*(1.d0-cost)+n(1)*sint
  R(3,3)=cost+n(3)**2*(1.d0-cost)

  return
end
!==============================================================================
subroutine orientation(geom,R,info)
  implicit none
  real*8, intent(inout) :: geom(3,4)
  real*8, intent(out) :: R(3,3)
  integer, intent(out) :: info

  call setorigin(geom)
  call bfxyz(geom(:,1),geom(:,2),geom(:,3),R,info)
  R=transpose(R)
  geom=matmul(R,geom)

  return
end
!==============================================================================
subroutine randomx(x1,x2,x)
  implicit none
  real*8, intent(in) :: x1,x2
  real*8, intent(out) :: x
  call random_number(x)
  x=x1+(x2-x1)*x
  return
end
!==============================================================================
subroutine randomspherepoint(n)
  implicit none
  real*8, intent(out) :: n(3)
  real*8, parameter :: pi=acos(-1.d0)
  real*8 :: theta, u
  call randomx(0.d0,2.d0*pi,theta)
  call randomx(-1.d0,1.d0,u)
  n(1)=sqrt(abs(1.d0-u**2))*cos(theta)
  n(2)=sqrt(abs(1.d0-u**2))*sin(theta)
  n(3)=u
  return
end
!==============================================================================
subroutine cross_product(a,b,c)
  implicit none
  real*8, intent(in) :: a(3),b(3)
  real*8, intent(out) :: c(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=-a(1)*b(3)+a(3)*b(1)
  c(3)=a(1)*b(2)-a(2)*b(1)
  return
end
!==============================================================================
subroutine cart2dist(natoms,x,r)
  implicit none
  integer, intent(in) :: natoms
  real*8, intent(in) :: x(3,natoms)
  real*8, intent(out) :: r(natoms*(natoms-1)/2)
  real*8 :: dx(3)
  integer :: i,j,k
  k=0
  do i=1,natoms-1
    do j=i+1,natoms
      k=k+1
      dx=x(:,i)-x(:,j)
      r(k)=sqrt(dot_product(dx,dx))
    end do
  end do
  return
end
!==============================================================================
subroutine bfxyz(a1,a2,a3,xyz,info)
  implicit none
  real*8, intent(in) :: a1(3),a2(3),a3(3)
  real*8, intent(out) :: xyz(3,3)
  integer, intent(out) :: info
  real*8 :: dx(3),norm

  info=0
  !x
  dx=a2-a1
  norm=sqrt(dot_product(dx,dx))
  if(norm .lt. 1.d-6) then
    write(*,"('Small norm of vector in determining new x axis: ',e15.6)") norm
    info=-1
  end if
  xyz(:,1)=dx/norm

  !z
  call cross_product(a2-a1,a3-a1,dx)
  norm=sqrt(dot_product(dx,dx))
  if(norm .lt. 1.d-6) then
    write(*,"('Small norm of vector in determining new z axis: ',e15.6)") norm
    print*, 'Warning: atoms may be colinear or nearly colinear.'
    info=-2
  end if
  xyz(:,3)=dx/norm

  !y
  call cross_product(xyz(:,3),xyz(:,1),xyz(:,2))

  return
end
!==============================================================================
