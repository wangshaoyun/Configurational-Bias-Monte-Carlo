module initialize_update
  !------------------------------------!
  !Use system parameters to initialize
  !the positions and update the positions.
  !------------------------------------!
  implicit none

  contains


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
  real*8 :: rnd1, rnd2, rnd3, rr, rr1
  real*8, dimension(3) :: rij, rnd(3)

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
      pos(k,1) = rnd1*Lx-Lx/2
      pos(k,2) = rnd2*Ly-Ly/2
      pos(k,3) = -Lz/2 + 1
      !
      !Jugde whether the particle is close the former paritcle
      !too much.
      do l = 1, i-1
        call rij_and_rr(rij,rr,(l-1)*Nml+1,k)
        rr1=pos(k,1)*pos(k,1)+pos(k,2)*pos(k,2)
        if (rr<0.8 .or. rr1<r_cy*r_cy) then
          m=1
          cycle
        end if
      end do
    end do
    if (mod(1,man)==0) then
      pos(k,4) = qq
    end if
    pos(k,5) = 1
    pos(k,6) = 1
    pos(k,7) = 1
    !
    !the remaining monomers of that polymer
    do j = 2, Nml
      k = k + 1
      pos(k,1) = pos(k-1,1)
      pos(k,2) = pos(k-1,2)
      pos(k,3) = pos(k-1,3) + R_bond
      if (mod(j,2)==0) then
        pos(k,4) = qq
        pos(k,5) = 1
      end if
      pos(k,6) = j
      pos(k,7) = j
    end do  
  end do

  !
  !aions of polymer
  do i=1, Nq_pe
    k = k + 1
    m = 1
    do while(m==1)
      m = 0
      call random_number(rnd)
      pos(k,1) = rnd(1)*Lx-Lx/2
      pos(k,2) = rnd(2)*Ly-Ly/2
      pos(k,3) = rnd(3)*Lz-Lz/2
      do l = 1, k-1
        call rij_and_rr(rij,rr,k,l)
        rr1 = pos(k,1)*pos(k,1)+pos(k,2)*pos(k,2)
        if (rr<0.8 .or. rr1<r_cy*r_cy) then
          m=1
          cycle
        end if
      end do
    end do
    pos(k,4) = -qq
    pos(k,5) = 1
  end do

  !
  !salt
  do i=1,Nq_salt_ions
    k = k + 1
    m = 1
    do while(m==1)
      m = 0
      call random_number(rnd)
      pos(k,1) = rnd(1)*Lx-Lx/2
      pos(k,2) = rnd(2)*Ly-Ly/2
      pos(k,3) = rnd(3)*Lz-Lz/2
      do l = 1, k-1
        call rij_and_rr(rij,rr,k,l)
        rr1 = pos(k,1)*pos(k,1)+pos(k,2)*pos(k,2)
        if (rr<0.8 .or. rr1<r_cy*r_cy) then
          m=1
          cycle
        end if
      end do
    end do
    pos(k,4) = qqi
    pos(k,5) = 1
  end do 

  !
  !aions of salt
  do i = 1, Nq_salt_ions*abs(qqi)
    k = k + 1
    m = 1
    do while(m==1)
      m = 0
      call random_number(rnd)
      pos(k,1) = rnd(1)*Lx-Lx/2
      pos(k,2) = rnd(2)*Ly-Ly/2
      pos(k,3) = rnd(3)*Lz-Lz/2
      do l = 1, k-1
        call rij_and_rr(rij,rr,k,l)
        rr1 = pos(k,1)*pos(k,1)+pos(k,2)*pos(k,2)
        if (rr<0.8 .or. rr1<r_cy*r_cy) then
          m=1
          cycle
        end if
      end do
    end do
    pos(k,4) = - qqi/abs(qqi)
    pos(k,5) = 1
  end do 

  !
  !aions of cylinder
  do i = 1, Nqc
    k = k + 1
    m = 1
    do while(m==1)
      m = 0
      call random_number(rnd)
      pos(k,1) = rnd(1)*Lx-Lx/2
      pos(k,2) = rnd(2)*Ly-Ly/2
      pos(k,3) = rnd(3)*Lz-Lz/2
      do l = 1, k-1
        call rij_and_rr(rij,rr,k,l)
        rr1 = pos(k,1)*pos(k,1)+pos(k,2)*pos(k,2)
        if (rr<0.8 .or. rr1<r_cy*r_cy) then
          m=1
          cycle
        end if
      end do
    end do
    pos(k,4) = qq
    pos(k,5) = 1
  end do 

  !
  !cylinder
  do i = 1, Npc
    k = k + 1
    m = 1
    do while(m==1)
      m = 0
      call random_number(rnd)
      pos(k,1) = rnd(1)*Lx-Lx/2
      pos(k,2) = rnd(2)*Ly-Ly/2
      pos(k,3) = rnd(3)*Lz-Lz/2
      do l = 1, k-1
        call rij_and_rr(rij,rr,k,l)
        rr1 = pos(k,1)*pos(k,1)+pos(k,2)*pos(k,2)
        if (rr<0.8 .or. rr1<r_cy*r_cy) then
          m=1
          cycle
        end if
      end do
    end do
    pos(k,4) = qqc
    pos(k,5) = 1
  end do   

end subroutine Initialize_position



subroutine CBMC_Move( EE, DeltaE, n0 )
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
  integer, intent(inout) :: n0
  integer :: i, j, k, l, m, n
  real*8 :: EE1, EE2, wo, wn

  i = n0
  do while ( i <= NN )

    call choose_chains

    call retrace( wo, DeltaE )

    call regrow( wn, DeltaE )

    call CBMC_Move_or_not( wo, wn, DeltaE, EE )

    j= nint(1.*num_newconf*Nq_salt_net/Npe) 

    do k = 1, j

      if (mod(k+i,Deltastep)==0) then

        call choose_particle_pH

        call Delta_energy_pH(DeltaE)

        call pH_move_or_not(DeltaE,EE)

      end if

      call choose_ions

      call delta_energy_ions_move(DeltaE)

      call ions_move_or_not(DeltaE,EE)

    end do

    i = i + num_newconf + j

  end do

  n0 = mod(i, NN) 

end subroutine CBMC_Move


subroutine choose_chains
  use global_variables
  implicit none
  real*8 :: rnd
  integer :: i, base

  call random_number(rnd)
  ic_newconf = int(rnd*Ngl) + 1        !chains

  call random_number(rnd)
  num_newconf = int(rnd*Nml) + 1       !regrow particles

  call random_number(rnd)
  cas = int(4*rnd) + 1                 !case 

  base = (ic_newconf-1)*Nml
  ib_newconf = base + Nml - num_newconf + 1  !starting particle

  !
  ! retrace from Nml to ib_newconf+1
  ! regrow from ib_newconf+1 to Nml
  select case ( cas )

    case (1) ! Remove from end and add to end

      pos_old = pos(base+1:base+Nml,:)
      pos_new = pos_old

    case (2) ! Remove from end and add to start

      pos_old = pos(base+1:base+Nml,:)
      pos_new(1:Nml-num_newconf,1:6) = pos(ib_newconf-1:base+1:-1,1:6)
      pos_new(Nml-num_newconf+1:Nml,5:6) = pos(base+Nml:ib_newconf:-1,5:6)
      do i = 1, Nml
        pos_new(pos_new(i,5),7) = i
      end do

    case (3) ! Remove from start and add to end

      pos_old = pos(base+Nml:base+1:-1,:)
      pos_new(1:Nml-num_newconf,1:6) = pos(base+num_newconf+1:base+Nml,1:6)
      pos_new(Nml-num_newconf+1:Nml,5:6) = pos(base+1:base+num_newconf,5:6)
      do i = 1, Nml
        pos_old(pos_old(i,5),7) = i
        pos_new(pos_new(i,5),7) = i
      end do

    case (4) ! Remove from start and add to start

      pos_old = pos(base+Nml:base+1:-1,:)
      pos_new(1:Nml-num_newconf,1:6) = pos(base+Nml:base+num_newconf+1:-1,1:6)
      pos_new(Nml-num_newconf+1:Nml,5:6) = pos(base+num_newconf:base+1:-1,5:6)
      do i = 1, Nml
        pos_old(pos_old(i,5),7) = i
        pos_new(pos_new(i,5),7) = i
      end do 

  end select

end subroutine choose_chains


subroutine choose_ions
  use global_variables
  implicit none
  real*8 :: rnd
  real*8 :: rnd1(3)

  call random_number(rnd)
  ip = int(rnd*(NN-Npe))+1+Npe
  do while(pos(ip,4)==0) 
    call random_number(rnd)
    ip = int(rnd*(NN-Npe))+1+Npe
  end do    

  pos_ip0=pos(ip,1:4)
  call random_number(rnd1)
  pos_ip1(1:3)=pos_ip0+rnd1*dr
  call periodic_condition(pos_ip1(1:3))

end subroutine choose_ions


subroutine choose_particle_pH
  use global_variables
  implicit none
  real*8 :: rnd
  real*8 :: rnd1(3)

  call random_number(rnd)
  ip = int(rnd*Npe)+1
  do while(pos(ip,5)==0)
    call random_number(rnd)
    ip = int(rnd*Npe)+1
  end do

  call random_number(rnd)
  ip1=Npe+int(rnd*Nq_pe)+1
  if (pos(ip,4)/=0)
    do while(pos(ip1,4)==0) then
      call random_number(rnd)
      ip1=Npe+int(rnd*Nq_pe)+1
    end do
  else
    do while(pos(ip1,4)/=0) then
      call random_number(rnd)
      ip1=Npe+int(rnd*Nq_pe)+1
    end do
  end if

  pos_ip0=pos(ip,1:4)
  pos_ip1=pos_ip0
  pos_ipi0=pos(ip1,1:4)
  if (pos_ip0(4)/=0) then
    pos_ip1(4)=0
    call random_number(rnd1)
    pos_ipi1(1)=rnd(1)*Lx-Lx/2
    pos_ipi1(2)=rnd(2)*Ly-Ly/2
    pos_ipi1(3)=rnd(3)*Lz-Lz/2
    pos_ipi1(4)=0
  else
    pos_ip1(4)=qq
    call random_number(rnd1)
    pos_ipi1(1)=rnd(1)*Lx-Lx/2
    pos_ipi1(2)=rnd(2)*Ly-Ly/2
    pos_ipi1(3)=rnd(3)*Lz-Lz/2
    pos_ipi1(4)=-qq
  end if

end subroutine choose_particle_pH


subroutine retrace( w, DeltaE )
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
  real*8, intent(out) :: DeltaE
  integer :: i, j, n
  real*8 :: eni, xt(3), wt, sumw

  w = 1
  DeltaE = 0
  do i = 1, num_newconf
    n = Nml - num_newconf + i
    call retrace_list(pos_old(n,7)+base,pos_old(n,1:4))
    if ( i == Nml ) then           ! When retrace only 1 monomer
      xt(:) = pos_old(1,1:4)
      call enerex_short(xt(:), eni)
      w = k_try*exp(-Beta*eni)
      DeltaE = DeltaE - eni
    else
      sumw = 0
      do j = 1, k_try
        if ( j==1 ) then
          xt(1,:) = pos_old(n,1:3)
        else
          call next_ci(n, xt(j,:))
        end if
        call enerex_short(xt(j,:), eni)
        wt(j) = exp( -Beta * eni )
        sumw = sumw + wt(j)
        if (j==1) then
          DeltaE = DeltaE - eni
        end if
      end do
      w = w * sumw
    end if
  end do

end subroutine retrace


subroutine regrow( w , DeltaE )
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
  real*8, intent(out) :: DeltaE
  integer :: i, j, m, n
  real*8  :: eni(k_try), xt(k_try,3), wt(k_try), sumw
  real*8  :: rnd(3)

  w = 1
  DeltaE = 0
  do i = 1, num_newconf
    n = Nml - num_newconf + i 
    if ( num_newconf == Nml .and. i == 1 ) then 
      call random_number(rnd)
      xt(1,1) = rnd(1) * Lx
      xt(1,2) = rnd(2) * Ly
      xt(1,3) = rnd(3) * Lz
      do while(xt(1,1)*xt(1,1)+xt(1,2)*xt(1,2)<r_cy) then
        xt(1,1) = rnd(1) * Lx
        xt(1,2) = rnd(2) * Ly
        xt(1,3) = rnd(3) * Lz
      end do
      call enerex_short(xt(1,:), eni(1))
      w = k_try*exp(-Beta*eni)
      DeltaE = DeltaE + eni
    else
      sumw = 0
      do j = 1, k_try
        call next_ci(n, xt(j,:))
        call enerex_short(xt(j,:), eni(j))
        wt(j) = exp( -Beta * eni )
        sumw = sumw + wt(j)
      end do
      w = w * sumw
      call select_move( wt, sumw, m )
      pos_new(n, :) = xt(m, :)
      DeltaE = DeltaE + eni(n)
    end if
    call regrow_list(pos_old(n,7)+base,pos_new(n,1:4))
  end do

end subroutine regrow


subroutine next_ci( n, xt )
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
  integer, intent(in) :: n
  real*8, dimension(3), intent(out) :: xt
  integer :: j

  if ( n == 2 ) then
    call next_c2(n, xt)
  elseif( n == 3 ) then
    call next_c3(n, xt)
  else
    call next_cn(n, xt)
  end if

end subroutine next_ci


subroutine next_c2(n, xt)
  use global_variables
  implicit none 
  integer, intent(in) :: n
  real*8, dimension(3), intent(out) :: xt
  real*8 :: l, b(3)

  call bond_l(l) 
  call ran_or(b)
  xt = pos(base+1,1:3) + l*b

end subroutine next_c2


subroutine next_c3(n, xt)
  use global_variables
  implicit none
  integer, intent(in) :: n
  real*8, dimension(3), intent(out) :: xt
  real*8 :: l, b(3)

  call bandl(l)
  call bond_a(n,b)
  xt = pos(base+2,:) + l*b  

end subroutine next_c3  


subroutine next_cn(n, xt)
  use global_variables
  implicit none
  real*8,intent(in) :: n
  real*8, dimension(3), intent(out) :: xt
  real*8 :: l, b(3)

  call bondl(l)
  call tors_bonda(n,b)
  xt = pos(base+n-1,1:3) + l*b

end subroutine next_cn


subroutine bondl(l)
  use global_variables
  implicit none
  real*8, intent(out) :: l
  real*8 :: sigma, l0, R0, kFENE, Kvib, a, rnd, Us
  logical :: FENE, ready

  FENE = .true.
  ready = .false.
  if (FENE) then
    R0 = 1.5
    R0_2 = R0*R0
    kFENE = 30
    do while ( .not. ready )
      call random_number(rnd)
      l = rnd * R0
      Us = - 0.5 * R0_2 * log( 1 - l*l / R0_2 )
      call random_number(rnd)
      if ( rnd < exp(-beta*Us) ) then
        ready = .true.
      end if
    end do
  else
    kvib = 400
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


subroutine bond_a(n, b)
  use global_variables
  implicit none
  real*8, dimension(Nml,3), intent(in) :: xn
  real*8, dimension(3), intent(out) :: b
  logical :: ready
  real*8  :: dx1x2(3), b(3), rnd, ubb
  real*8  :: k_phi, phi0, phi, rr

  k_phi = 30
  phi0 = 0
  ready = .false.
  do while ( .not. ready )
    call ran_or(b)
    dx1x2 = pos(base+n-1,:) - pos(base+n-2,:)
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


subroutine tors_bonda(n, b)
  use global_variables
  implicit none
  real*8, intent(in) :: n
  real*8, dimension(3), intent(out) :: b
  logical :: ready
  real*8 :: dx1x2(3), dx2x3(3), xx1(3), xx2(3), r12, r23
  real*8 :: phi, phi0, k_phi, theta, theta0, k_theta
  real*8 :: rnd, ubb, utors, usum

  ready = .false.
  do while( .not. ready )
    call ran_or(b)

    dx1x2 = pos(base+n-1,:) - pos(base+n-2,:)
    r12 = sqrt( dot_product(dx1x2, dx1x2) )
    dx1x2 = dx1x2 / r12

    dx2x3 = pos(base+n-2,:) - pos(base+n-3,:)
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


subroutine select_move(wt, sumw, n)
  use global_variables
  implicit none
  real*8, dimension(k_try), intent(in) :: wt
  real*8, intent(in) :: sumw
  integer, intent(out) :: n
  real*8 :: rnd, ws, cumw

  call random_number(rnd)
  ws = rnd * sumw
  cumw = wt(1)
  n = 1
  do while ( cumw < ws )
    n = n + 1
    cumw = cumw + wt(n)
  end do

end subroutine select_move


subroutine CBMC_Move_or_not(wo, wn, EE, DeltaE)
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
  implicit none
  real*8, intent(in) :: wo
  real*8, intent(in) :: wn
  real*8, intent(inout) :: EE
  real*8, intent(in) :: DeltaE
  real*8 :: rnd, Del_E
  integer :: i, j
  !
  !Long range energy difference
  call energy_long(Del_E)
  !
  !Judge whether move or not
  call random_number( rnd )
  if (rnd < (wn/wo) ) then
    pos((ic_newconf-1)*Nml+1:ic_newconf*Nml,:) = pos_new(:,:)
    accpt_num = accpt_num + num_newconf
    EE = EE + De
    call grow_list(.true.)
  else
    call grow_list(.false.)
  end if
  total_num = total_num + num_newconf

end subroutine CBMC_Move_or_not


subroutine ions_move_or_not(De,EE)
  use global_variables
  implicit none
  real*8, intent(in) :: De
  real*8, intent(inout) :: EE
  real*8 :: rnd

  if (De<0) then
    pos(ip,1:4) = pos_ip1
    call update_cell_list
    EE = EE + De
  else
    call random_number(rnd)
    if (rnd<exp(-De*beta)) then
      pos(ip,1:4) = pos_ip1
      call update_cell_list
      EE = EE + De
    end if
  end if
 
end subroutine ions_move_or_not


subroutine pH_move_or_not(De)
  use global_variables
  implicit none
  real*8, intent(in) :: De
  real*8, intent(inout) :: EE
  real*8 :: rnd

  if (De<0) then
    pos(ip,1:4) = pos_ip1
    pos(ip1,1:4) = pos_ipi1
    call update_cell_list
    EE = EE + De
  else
    call random_number(rnd)
    if (rnd<exp(-De*beta)) then
      pos(ip,1:4) = pos_ip1
      pos(ip1,1:4) = pos_ipi1
      call update_cell_list
      EE = EE + De
    end if
  end if

end subroutine pH_move_or_not


subroutine grow_list(new_conf)
  use compute_energy
  implicit none
  logical, intent(in) :: new_conf

  call add_rho_k_CBMC(new_conf)

  call add_cell_list_r1_CBMC(new_conf)

  call add_cell_list_r2_CBMC(new_conf)

  call add_cell_list_lj_CBMC(new_conf)

end subroutine grow_list


subroutine retrace_list(np, rr)
  use global_variables
  use compute_energy
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr

  if (rr(4)/=0) then

    call delete_cell_list_real_CBMC(np,rr)

  end if

  call delete_cell_list_lj_CBMC(np,rr)

end subroutine retrace_list


subroutine regrow_list(np, rr)
  use global_variables
  use compute_energy
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr

  if (pos(np,4)/=0) then

    call add_cell_list_real_CBMC(np,rr)

  end if

  call add_cell_list_lj_CBMC(np,rr)

end subroutine regrow_list


end module initialize_update



