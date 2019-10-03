module initialize_update
  !------------------------------------!
  !Use system parameters to initialize
  !the positions and update the positions.
  !------------------------------------!
  implicit none

  contains


! subroutine Initialize_position
!   !--------------------------------------!
!   !Initialize linear chains.
!   !
!   !Input
!   !  pos
!   !Output
!   !  pos
!   !External Variables
!   !  
!   !Routine Referenced:
!   !   rij_and_rr, period_condition_rij
!   !Reference:
!   !   Spherical coordinate on Wiki to generate uniform distribution
!   !   on the sphere surface.
!   !--------------------------------------!
!   use global_variables
!   implicit none
!   integer :: i, j, k, l, m, n, x, y, p
!   real*8 :: theta, rnd1, rnd2, rnd3, rsqr
!   real*8, dimension(3) :: rij


!   do i=1, Ngl
!     l = (i-1)*Nml + 1
!     x=(i-1)/nint(sqrt(Ngl*1.))+1
!     y=mod(i-1,nint(sqrt(Ngl*1.)))+1
!     pos(l,1)=Lx/nint(sqrt(Ngl*1.))*(x-0.5)-Lx/2
!     pos(l,2)=Ly/nint(sqrt(Ngl*1.))*(y-0.5)-Ly/2
!     pos(l,3)=Lz/nint(sqrt(Ngl*1.))*0.5-Lz/2
!     do k=2, Nml
!       l=l+1
!       m=1
!       p=0
!       do while (m==1)
!         m=0
!         call random_number(rnd1)
!         if (p<10) then
!           rnd1=rnd1/2
!         else
!           rnd1=rnd1**2
!         end if
!         call random_number(rnd2)
!         pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
!         pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
!         pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
!         !periodic condition
!         call periodic_condition(pos(l,1:3))
!         !
!         !Judge whether the particle is close the former paritcle
!         !too much.
!         do n=1,l-1
!           call rij_and_rr(rij,rsqr,n,l)
!           if (rsqr<0.81 .or. pos(l,3)<(-Lz/2+1) ) then
!             m=1
!             p=p+1
!             cycle
!           end if
!         end do
!       end do
!     end do
!   end do

! end subroutine Initialize_position


subroutine Initialize_position
  !------------------------------------!
  !Initialize position
  !   This program is used to initialize the position of
  !   Polyelectrolytes and ions, and parameters, energy of
  !   the potential.
  !Input
  !   pos, random_or_uniform
  !Output
  !   pos
  !External Variables
  !   pos, random_or_uniform
  !Routine Referenced:
  !1.subroutine random_grafted
  !   initialize chains by randomly grafting on the plate
  !2.subroutine uniform_grafted
  !   initialize chains by uniformly grafting on the plate
  !3.subroutine initialize_ions
  !   initialize ions in the system
  !------------------------------------!
  use global_variables
  implicit none
  integer :: i,j,k,l,m
  real*8 :: rnd1, rnd2, rnd3, rr
  real*8, dimension(3) :: rij

  pos=0

  do i = 1, Ngl
    !
    ! the start monomer of each polymer
    m = 1
    k = (i-1)*Nml+1
    do while (m==1)
      m = 0
      call random_number(rnd1)
      call random_number(rnd2)
      call random_number(rnd3)
      pos(k,1) = rnd1*Lx-Lx/2
      pos(k,2) = rnd2*Ly-Ly/2
      pos(k,3) = rnd3*Lz-Lz/2
      !
      !Jugde whether the particle is close the former paritcle
      !too much.
      do l = 1, k-1
        call rij_and_rr(rij,rr,l,k)
        if (rr<0.8) then
          m=1
          cycle
        end if
      end do
    end do
    !
    !the remaining monomers of that polymer
    do j = 2, Nml
      k = k + 1
      m = 1
      do while (m == 1)
        m = 0
        call random_number(rnd1)
        call random_number(rnd2)
        pos(k,1) = pos(k-1,1) + R_bond*cos(2*pi*rnd2)*sin(pi*rnd1)
        pos(k,2) = pos(k-1,2) + R_bond*sin(2*pi*rnd2)*sin(pi*rnd1)
        pos(k,3) = pos(k-1,3) + R_bond*cos(pi*rnd1)
        !periodic condition
        call periodic_condition(pos(k,1:3))
        !
        !Jugde whether the particle is close the former paritcle
        !too much.
        do l = 1, k-1
          call rij_and_rr(rij,rr,l,k)
          if (rr<0.81) then
            m=1
            cycle
          end if
        end do
      end do
    end do 

  end do

end subroutine Initialize_position


subroutine CBMC_Move( EE, DeltaE )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out)   :: DeltaE
  integer :: m, n
  real*8 :: EE1, EE2, wo, wn
  logical :: new_conf

  n = 0
  do while ( n <= NN )
  !     call total_energy(EE1)

    call retrace( wo )

    call regrow( wn )

    call Move_or_not( wo, wn )

    n = n + ib
    !
    !test EE2-EE1 = DeltaE
  !     call total_energy(EE2)
  !     write(*,*) EE2 - EE1, DeltaE, EE2, EE1  
  end do

end subroutine CBMC_Carlo_Move


subroutine retrace( w )
  !------------------------------------!
  !   
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(out) :: w
  integer :: i, j, k, l, n
  real*8 :: rnd1, rnd2, ib, ic, sumw, rnd(3)
  real*8 :: eni, xt(k_try,3), xn(Nml,3)

  ib = ibnewconf
  ic = ic_newconf

  do i = 1, Nml
    j = ic*Nml - ib + i
    xn(i,:) =  pos(j,:)
  end do

  w = 1
  do i = 1, ib
    if ( i == Nml ) then
      xn(1,:) = pos(ic*Nml-ib+1,:)
      xt(1,:) = xn(1,:)
      call enerex(ic, ib, xt(1,:), ib-i+1, eni)
      w = k_try*exp(-Beta*eni)
    else
      sumw = 0
      do j = 1, k
        if ( j==1 ) then
          xt(1) = pos(ic*Nml-i+1,:)
        else
          call next_ci(ib, xt(j,:), xn, ib-i+1)
        end if
        call enerex(ic, ib, xt(j,:), ib-i+1, eni)
        wt(j) = exp( -Beta * eni )
        sumw = sumw + wt(j)
      end do
      w = w * sumw
    end if
  end do

end subroutine retrace


subroutine regrow( w )
  !------------------------------------!
  !   
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(out) :: w
  integer :: i, j, k, l, n
  real*8 :: rnd1, rnd2, ib, ic, sumw, rnd(3)
  real*8 :: eni, xt(k_try,3), xn(Nml,3)

  call random_number(rnd1)
  call random_number(rnd2)
  ic = int( rnd1 * Ngl ) + 1
  ib = int( rnd2 * Nml ) + 1
  ic_newconf = ic
  ib_newconf = ib
  deallocate( pos_new )
  allocate( pos_new(ib_newconf, 3) )

  do i = 1, ib - 1
    j = ic*Nml - ib + i
    xn(i,:) =  pos(j,:)
  end do

  w = 1
  do i = 1, ib
    if ( ib == Nml .and. i == 1 ) then 
      call random_number(rnd)
      xn(1, 1) = rnd(1) * Lx
      xn(1, 2) = rnd(2) * Ly
      xn(1, 3) = rnd(3) * Lz
      xt(1, :) = xn(1, :)
      call enerex(ic, ib, xt(1,:), i, eni)
      w = k_try*exp(-Beta*eni)
    else
      sumw = 0
      do j = 1, k
        call next_ci(ib, xt(j,:), xn, i)
        call enerex(ic, ib, xt(j,:), i, eni)
        wt(j) = exp( -Beta * eni )
        sumw = sumw + wt(j)
      end do
      w = w * sumw
      call select( wt, sumw, n )
      xn(i, :) = xt(n, :)
    end if
  end do
  pos_new = xn( Nml-ib+1:Nml, : )

end subroutine regrow


subroutine next_ci( ib, xt, xn, i )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: ib
  integer, intent(in) :: i
  real*8, dimension(3),    intent(out) :: xt
  real*8, dimension(Nml,3), intent(in) :: xn
  integer :: j

  j = Nml - ib + i
  if ( j == 2 ) then
    call next_c2(xn, xt)
  elseif( j == 3 ) then
    call next_c3(xn, xt, j)
  else
    call next_cn(xn, xt, j)
  end if

end subroutine next_ci


subroutine next_c2(xn, xt)
  use global_variables
  implicit none 
  real*8, dimension(Nml,3), intent(in) :: xn
  real*8, dimension(3), intent(out) :: xt
  real*8 :: l, b(3)

  call bond_l(l) 
  call ran_or(b)
  xt = xn(1,:) + l*b

end subroutine next_c2


subroutine next_c3(xn, xt)
  use global_variables
  implicit none
  real*8, dimension(Nml,3), intent(in) :: xn
  real*8, dimension(3), intent(out) :: xt
  integer, intent(in) :: j
  real*8 :: l, b(3)

  call bandl(l)
  call bond_a(xn, b, j)
  xt = xn(2,:) + l*b  

end subroutine next_c3  


subroutine next_cn(xn, xt, j)
  use global_variables
  implicit none
  real*8, dimension(Nml,3), intent(in) :: xn
  real*8, dimension(3), intent(out) :: xt
  integer, intent(in) :: j
  real*8 :: l, b(3)

  call bondl(l)
  call tors_bonda(xn, b, j)
  xt = xn(j-1) + l*b

end subroutine next_cn


subroutine bondl(l)
  use global_variables
  implicit none
  real*8, intent(out) :: l
  real*8 :: sigma, l0, R0, kFENE, Kvib, a, rnd, Us
  logical :: FENE, ready

  FENE = .false.
  ready = .false.
  if (FENE) then
    R0 = 1.5
    kFENE = 30
    do while ( .not. ready )
      call random_number(rnd)
      l = rnd * R0
      Us = - 0.5 * R0 * R0 * log( 1 - l**2 / R0**2 )
      call random_number(rnd)
      if ( rnd < exp(-beta*Us) ) then
        ready = .true.
      end if
    end do
  else
    sigma = sqrt(1/(beta*Kvib))
    l0 = 1
    a = ( l0 + 3*simga )**2
    do while (.not. ready)
      call gauss(sigma, l0, l)
      call random_number(rnd)
      if ( rnd < l*l/a ) then
        ready = .true.
      end if
    end do
  end if

end subroutine bondl


subroutine ran_or(b)
  use global_variables
  implicit none
  real*8, dimension(3), intent(out) :: b
  real*8 :: ransq, ran1, ran2, ranh

  ransq = 2
  do while( ransq > 1)
    call random_number(ran1)
    call random_number(ran2)
    ran1 = 1 - 2 * ran1
    ran2 = 1 - 2 * ran2
    ransq = ran1 * ran1 + ran2 * ran2
  end do
  ranh = 2 * sqrt( 1 - ransq )
  b(1) = ran1 * ranh
  b(2) = ran2 * ranh
  b(3) = 1 - 2 * ransq

end subroutine ran_or


subroutine bond_a(xn, b, j)
  use global_variables
  implicit none
  real*8, dimension(Nml,3), intent(in) :: xn
  real*8, dimension(3), intent(out) :: b
  integer, intent(in) :: j
  logical :: ready
  real*8  :: dx1x2(3), b(3), rnd, ubb
  real*8  :: k_phi, phi0, phi, rr

  ready = .false.
  do while ( .not. ready )
    call ran_or(b)
    dx1x2 = xn(j-1,:) - xn(j-2,:)
    rr = sqrt( dot_product(dx1x2, dx1x2) )
    dx1x2 = dx1x2 / rr
    phi = acos( dot_product(dx1x2, b) )
    ubb = k_phi * ( phi - phi0 )**2
    call random_number(rnd)
    if ( rnd < exp(-beta*ubb) ) then
      ready = .true.
    end if
  end do

end subroutine bond_a


subroutine tors_bonda(xn, b, j)
  use global_variables
  implicit none
  real*8, dimension(Nml,3), intent(in) :: xn
  real*8, dimension(3), intent(out) :: b
  integer :: j
  logical :: ready
  real*8 :: dx1x2(3), dx2x3(3), xx1(3), xx2(3), r12, r23
  real*8 :: phi, phi0, k_phi, theta, theta0, k_theta
  real*8 :: rnd, ubb, utors, usum

  ready = .false.
  do while( .not. ready )
    call ran_or(b)

    dx1x2 = xn(i-1,:) - xn(i-2,:)
    r12 = sqrt( dot_product(dx1x2, dx1x2) )
    dx1x2 = dx1x2 / r12

    dx2x3 = xn(i-2,:) - xn(i-3,:)
    r23 = sqrt( dot_product(dx2x3, dx2x3) )
    dx2x3 = dx2x3 / r23

    phi = acos(dot_product(b,dx1x2))
    ubb = k_phi * ( phi - phi0 )**2

    call cross_product(b, dx1x2, xx1)
    call cross_product(dx1x2, dx2x3, xx2)
    theta = acos(dot_product(xx1,xx2))
    utors = k_theta * (theta - theta0)**2

    usum = ubb + utors
    call random_number(rnd)
    if ( rnd < exp(-beta*usum) ) then
      ready = .true.
    end if
  end do

end subroutine tors_bonda


subroutine Move_or_not(wo, wn, m)
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   pos, pos_ip0, pos_ip1, ip, Beta
  !Routine Referenced:
  !1.
  !Reference:
  !
  !--------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8,  intent(in) :: wo
  real*8,  intent(in) :: wn
  integer, intent(in) :: m
  real*8 :: rnd
  integer :: i, j
  !
  !Judge whether move or not
  call random_number( rnd )
  if (rnd < (wn/wo) ) then
    i = ic_newconf * Nml - ib_newconf + 1
    j = ic_newconf * Nml
    pos(i:j,:) = pos_new(:,:)
    accpt_num = accpt_num + m
  end if
  total_num = total_num + m
end subroutine Move_or_not


end module initialize_update



