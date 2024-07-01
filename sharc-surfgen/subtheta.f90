subroutine Evaluatek
  implicit none
  real*8:: gmek(3,4)
  real*8:: the
  real*8:: error
  integer:: k
  real*8  :: hd(2,2),h12(2)
  real*8  :: dplmek(3,2,2)
  real*8  ::e(2),cg(12),dcg(12)
  integer :: i,j,info,idx(4)
  real*8  :: adplmek(3,2,2),dis(4)
  real*8  :: T(3,3),theta(4),dplt(3,2,2,4),adpl(3,2,2,4),s
  parameter pi=3.14159265

  gmek=geom(:,:,i)
  call EvaluateSurfgen(gmek,e,cg,hd,dcg)
  print *,hd
  theta(1)=0.5d0*atan(-2*hd(1,2)/(hd(2,2)-hd(1,1)))
  theta(2)=theta(1)+0.5d0*pi
!select trn dipole sign and k
  h12(1)=-0.5d0*sin(2d0*theta(1))*(e(1)-e(2))-hd(1,2)
  h12(2)=-0.5d0*sin(2d0*theta(2))*(e(1)-e(2))-hd(1,2)
  if (abs(h12(1)).lt.abs(h12(2))) then
    k=0
    the=theta(1)
    error=h12(1)
  else
    k=1
    the=theta(2)
    error=h12(2)
  endif
  write(21,"(i5,2f14.8)")k,the,error
  return
end subroutine

subroutine overlap(n,x,y,t)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(n),y(n)
  real*8, intent(out) :: t
  real*8 :: dx,dy
  dx=sqrt(dot_product(x(1:n),x(1:n)))
  dy=sqrt(dot_product(y(1:n),y(1:n)))
  dx=max(dx,1.d-20)
  dy=max(dy,1.d-20)
  t=dot_product(x(1:n),y(1:n))/(dx*dy)
  return
end subroutine

subroutine calcdpl(adiadpl,theta,diadpl)
  implicit none
  real*8, intent(in) :: adiadpl(3,2,2),theta
  real*8, intent(out) :: diadpl(3,2,2)
  real*8 :: D(2,2),dpl(3,2,2)
  integer :: i

  D(1,1)=dcos(theta)
  D(1,2)=-dsin(theta)
  D(2,1)=-D(1,2)
  D(2,2)=D(1,1)

  dpl=adiadpl
!  dpl(:,1,2)=-adiadpl(:,1,2)
!  dpl(:,2,1)=dpl(:,1,2)
  do i=1,3
    diadpl(i,:,:)=matmul(matmul(D,dpl(i,:,:)),transpose(D))
  end do

  return
end subroutine

subroutine sort0(n,a,idx)
  use ifport
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: a(n)
  integer, intent(out) :: idx(n)
  integer(SIZEOF_SIZE_T) :: array_len,array_size
  integer :: i
  array_len=n
  array_size=4
  forall(i=1:n) idx(i)=i
  call qsort(idx,array_len,array_size,compar)
  contains
  function compar(i,j)
    integer*2 compar
    integer :: i,j
    if(a(i).lt.a(j)) then
      compar=-1
    else if(a(i).gt.a(j)) then
      compar=1
    else
      compar=0
    end if
    return
  end function
end
