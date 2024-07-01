!=============================================================================
function triple_product(a,b,c)
  implicit none
  real*8 :: triple_product
  real*8, intent(in) :: a(3),b(3),c(3)
  triple_product = a(1)*b(2)*c(3) - a(1)*b(3)*c(2) - a(2)*b(1)*c(3) + &
                   a(2)*b(3)*c(1) + a(3)*b(1)*c(2) - a(3)*b(2)*c(1)
  return
end
!=============================================================================
subroutine OOP0(a,b,c,d,coord,dcoord,id)
  implicit none
  real*8, intent(in) :: a(3),b(3),c(3),d(3)
  integer, intent(in) :: id
  real*8, intent(out) :: coord, dcoord(12)
  real*8, external :: triple_product
  real*8 :: ab(3), ac(3), ad(3), rab, rac, rad
  real*8 :: abx(12),acx(12),adx(12),dtp(12)

  ab=a-b
  ac=a-c
  ad=a-d
  rab=sqrt(dot_product(ab,ab))
  rac=sqrt(dot_product(ac,ac))
  rad=sqrt(dot_product(ad,ad))
  coord=triple_product(ab,ac,ad)/(rab*rac*rad)

  if(id.ne.0) then
    abx=0.d0
    abx(1:3)=ab/rab
    abx(4:6)=-abx(1:3)
    acx=0.d0
    acx(1:3)=ac/rac
    acx(7:9)=-acx(1:3)
    adx=0.d0
    adx(1:3)=ad/rad
    adx(10:12)=-adx(1:3)
    !tp = ab(1)*ac(2)*ad(3) - ab(1)*ac(3)*ad(2) - ab(2)*ac(1)*ad(3) + &
    !     ab(2)*ac(3)*ad(1) + ab(3)*ac(1)*ad(2) - ab(3)*ac(2)*ad(1)
    dtp(1)=ac(2)*ad(3)-ac(3)*ad(2)-ab(2)*ad(3)+ab(2)*ac(3)+ab(3)*ad(2)-ab(3)*ac(2)
    dtp(2)=ab(1)*ad(3)-ab(1)*ac(3)-ac(1)*ad(3)+ac(3)*ad(1)+ab(3)*ac(1)-ab(3)*ad(1)
    dtp(3)=ab(1)*ac(2)-ab(1)*ad(2)-ab(2)*ac(1)+ab(2)*ad(1)+ac(1)*ad(2)-ac(2)*ad(1)
    dtp(4)=-ac(2)*ad(3)+ac(3)*ad(2)
    dtp(5)=ac(1)*ad(3)-ac(3)*ad(1)
    dtp(6)=-ac(1)*ad(2)+ac(2)*ad(1)
    dtp(7)=ab(2)*ad(3)-ab(3)*ad(2)
    dtp(8)=-ab(1)*ad(3)+ab(3)*ad(1)
    dtp(9)=ab(1)*ad(2)-ab(2)*ad(1)
    dtp(10)=-ab(2)*ac(3)+ab(3)*ac(2)
    dtp(11)=ab(1)*ac(3)-ab(3)*ac(1)
    dtp(12)=-ab(1)*ac(2)+ab(2)*ac(1)
    dcoord=dtp/(rab*rac*rad)-coord*(abx/rab+acx/rac+adx/rad)
  end if

  return
end
!================================================================================
subroutine permdiadpl(nstates,id,diadpl)
  integer, intent(in) :: nstates,id
  real*8, intent(inout) :: diadpl(3,nstates,nstates)

  select case(id)
  case(1) !E
  case(2) !(123)
  case(3) !(132)
  case(4) !(12)
    diadpl(:,1,2)=-diadpl(:,1,2)
    diadpl(:,2,1)=diadpl(:,1,2)
  case(5) !(23)
    diadpl(:,1,2)=-diadpl(:,1,2)
    diadpl(:,2,1)=diadpl(:,1,2)
  case(6) !(13)
    diadpl(:,1,2)=-diadpl(:,1,2)
    diadpl(:,2,1)=diadpl(:,1,2)
  case default
    stop 'Wrong permutation number!'
  end select

  return
end
!===============================================================================
subroutine invpermdiadpl(nstates,id,diadpl)
  integer, intent(in) :: nstates,id
  real*8, intent(inout) :: diadpl(3,nstates,nstates)

  select case(id)
  case(1) !E <-> E
  case(2) !(123) <-> (132)
  case(3) !(132) <-> (123)
  case(4) !(12) <-> (12)
    diadpl(:,1,2)=-diadpl(:,1,2)
    diadpl(:,2,1)=diadpl(:,1,2)
  case(5) !(23) <-> (23)
    diadpl(:,1,2)=-diadpl(:,1,2)
    diadpl(:,2,1)=diadpl(:,1,2)
  case(6) !(13) <-> (13)
    diadpl(:,1,2)=-diadpl(:,1,2)
    diadpl(:,2,1)=diadpl(:,1,2)
  case default
    stop 'Wrong permutation number!'
  end select

  return
end
!===============================================================================
