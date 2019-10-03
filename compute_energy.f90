module compute_energy
  !--------------------------------------!
  !Input:
  ! pos, pos_ip0, pos_ip1, ip
  ! and the system parameters
  !Output: 
  ! EE, DeltaE 
  !--------------------------------------!
  implicit none

  save

!############coefficient in potential function#############!
!
!lj potential
  real*8,  private :: epsilon     !Energy unit epsilon in lj potential
  real*8,  private :: sigma       !Distance sigma in lj potential
  real*8,  private :: rc_lj       !Cut off radius of LJ potential
  real*8,  private :: rv_lj       !Verlet list radius of LJ potential
  real*8,  private :: rsk_lj      !Skin between cut off sphere and verlet list 
                                  !sphere
  integer, private :: npair1      !number of pairs in the lj verlet sphere
!
!fene potential
  real*8,  private :: R0_2        !Max bond length R0 square in FENE potential
  real*8,  private :: Kvib        !FENE spring constant k
  integer, private :: N_bond      !Number of Chemical bond of polymers
!##########end coefficient in potential function###########!


!##########################arrays##########################!
  integer, allocatable, dimension( : ), private :: lj_pair_list  
                                  !LJ potential verlet list
  integer, allocatable, dimension( : ), private :: lj_point
                                  !the particles near i are from
                                  !lj_pair_list(lj_point(i-1)) to 
                                  !lj_pair_list(lj_point(i))
!########################end arrays########################!


contains


subroutine initialize_energy_parameters
  !--------------------------------------!
  !Initial parameters are not inputted from file and compute
  !the total energy of the system.
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !Reference:
  !The computation of alpha, rc_real et al are refered to
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.304-306.
  !--------------------------------------!
  use global_variables
  implicit none
  !
  !read energy parameters from file
  call read_energy_parameters

  if (rc_lj<Lx/20) then
    !
    !Initialize lj parameters and array allocate.
    call initialize_lj_parameters
    !
    !build lj_pair_list and lj_point
    call build_lj_verlet_list
  end if

end subroutine initialize_energy_parameters


subroutine total_energy (EE)
  !--------------------------------------!
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
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: EE

  EE=0

  call LJ_Energy(EE)

end subroutine total_energy


subroutine enerex(ic, ib, xt, i, eni)
  use global_variables
  implicit none
  integer, intent(in) :: ic
  integer, intent(in) :: ib
  integer, intent(in) :: i
  real*8, intent(out) :: eni
  real*8, dimension(3), intent(in) :: xt

  

end subroutine enerex


subroutine LJ_energy (EE)
  !--------------------------------------!
  !Compute total LJ potential energy,
  !including LJ energy of wall.
  !   
  !Input
  !   EE
  !Output
  !   EE
  !External Variables
  !   lj_point, lj_pair_list, pos, 
  !   epsilon, sigma, rc_lj, Lz
  !Routine Referenced:
  !1. rij_and_rr( rij, rr, i, j )
  !Reference:
  !1.In fact, the cut-off radius in good solvent is 2^(1/6), 
  !  which was first obatained by JOHN D. WEEKS, DAVID CHANDLER
  !  and HANS C. ANDERSEN. So it is called WCA potential.
  !  JOHN D. WEEKS, DAVID CHANDLER and HANS C. ANDERSEN, 'Role of 
  !  Repulsive Forces in Determining the Equilibrium Structure of
  !  Simple Liquids', THE JOURNAL OF CHEMICAL PHYSICS, Vol. 54, 
  !  pp.5237-5247, (1971).
  !2.The potential of particle and wall are 9-3 LJ potential which is
  !  cut off at 0.4^(1/6) = 0.86.
  !  Yu-Fan Ho, et al, 'Structure of Polyelectrolyte Brushes Subject
  !  to Normal Electric Fields', Langmuir, 19, pp.2359-2370, (2013).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: EE
  integer :: i, j, k, l, m
  real*8  :: rr, rij(3), inv_rr2, inv_rr6

  if (rc_lj>Lx/10) then
    do i = 1, NN-1
      do j = i+1, NN
        call rij_and_rr( rij, rr, i, j )
        inv_rr2  = sigma*sigma/rr
        inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
        EE = EE + 4 * epsilon * ( inv_rr6 * inv_rr6 - inv_rr6 + 0.25D0)
      end do 
    end do
  else
    do i = 1, NN
      if ( i == 1) then
        k = 1
        l = lj_point(1)
      else
        k = lj_point(i-1)+1
        l = lj_point(i)
      end if
      do m = k, l
        j = lj_pair_list(m)
        call rij_and_rr( rij, rr, i, j )
        if ( rr < rc_lj * rc_lj ) then
          inv_rr2  = sigma*sigma/rr
          inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
          EE = EE + 4 * epsilon * ( inv_rr6 * inv_rr6 - inv_rr6 + 0.25D0) / 2
          ! ! ! must divided by 2 because of the repeating cycle
        end if
      end do
    end do
  end if

end subroutine LJ_energy


subroutine Delta_Energy(DeltaE)
  !--------------------------------------!
  !Compute change of energy.
  !   
  !Input
  !   
  !Output
  !   DeltaE
  !External Variables
  !   pos_ip0, pos_ip1, ip
  !   inv_charge, DeltaE, EF
  !Routine Referenced:
  !1.Delta_LJ_Energy(DeltaE)
  !2.Delta_FENE_Energy(DeltaE)
  !3.Delta_real_Energy(DeltaE)
  !4.Delta_Reciprocal_Energy(DeltaE)
  !--------------------------------------!
  use global_variables
  implicit none
	real*8,  intent(out) :: DeltaE

  DeltaE = 0
  !
  !Compute energy of LJ potential
  call Delta_LJ_Energy(DeltaE)
  !
  !Compute Delta energy of FENE potential
  call Delta_FENE_Energy(DeltaE)

end subroutine Delta_Energy


subroutine Delta_lj_Energy(DeltaE)
  !--------------------------------------!
  !Compute change of LJ potential Energy.
  !   
  !Input
  !   DeltaE
  !Output
  !   DeltaE
  !External Variables
  !   pos, lj_pair_list, lj_point
  !   pos_ip0, pos_ip1, ip
  !   Lx, Ly, Lz, sigma, epsilon, rc_lj
  !Routine Referenced:
  !
  !Reference:
  !In fact, the cut-off radius in good solvent is 2^(1/6), 
  !which was first obatained by JOHN D. WEEKS, DAVID CHANDLER
  !and HANS C. ANDERSEN. So it is called WCA potential.
  !JOHN D. WEEKS, DAVID CHANDLER and HANS C. ANDERSEN, 'Role of 
  !Repulsive Forces in Determining the Equilibrium Structure of
  !Simple Liquids', THE JOURNAL OF CHEMICAL PHYSICS, Vol. 54, 
  !pp.5237-5247, (1971).
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8  :: EE, sigma2, rc_lj2
  real*8  :: rij(3), rr, inv_rr2, inv_rr6, inv_rr12
  integer :: i, j, k, l

  EE     = 0
  sigma2 = sigma * sigma
  rc_lj2 = rc_lj * rc_lj

  if (rc_lj>Lx/20) then
    do i = 1, NN
      if ( i == ip ) cycle
      !
      !Energy of old configuration
      !
      rij = pos(i, 1:3) - pos_ip0(1:3)
      !
      !periodic condition
      call periodic_condition(rij)
      !
      !lj energy
      rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
      inv_rr2  = sigma2 / rr
      inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
      inv_rr12 = inv_rr6 * inv_rr6
      EE       = EE + inv_rr6 - inv_rr12 - 0.25D0
      !
      !Energy of new configuration
      !
      rij = pos(i, 1:3) - pos_ip1(1:3)
      !
      !periodic condition
      call periodic_condition(rij)
      !
      !lj energy
      rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
      inv_rr2  = sigma2 / rr
      inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
      inv_rr12 = inv_rr6 * inv_rr6
      EE       = EE + inv_rr12 - inv_rr6 + 0.25D0
    end do
  else
    if (ip==1) then
      k = 1
      l = lj_point( ip )
    else
      k = lj_point( ip-1 ) + 1
      l = lj_point( ip )
    end if

    do j= k, l
      i = lj_pair_list(j)
      !
      !Energy of old configuration
      !
      rij = pos(i, 1:3) - pos_ip0(1:3)
      !
      !periodic condition
      call periodic_condition(rij)
      !
      !lj energy
      rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
      if ( rr < rc_lj2 ) then
        inv_rr2  = sigma2 / rr
        inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
        inv_rr12 = inv_rr6 * inv_rr6
        EE       = EE + inv_rr6 - inv_rr12 - 0.25D0
      end if
      !
      !Energy of new configuration
      !
      rij = pos(i, 1:3) - pos_ip1(1:3)
      !
      !periodic condition
      call periodic_condition(rij)
      !
      !lj energy
      rr = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
      if ( rr < rc_lj2 ) then
        inv_rr2  = sigma2 / rr
        inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
        inv_rr12 = inv_rr6 * inv_rr6
        EE       = EE + inv_rr12 - inv_rr6 + 0.25D0
      end if
    end do
  end if

  DeltaE = DeltaE + 4 * epsilon * EE

end subroutine Delta_lj_Energy


subroutine read_energy_parameters
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  open(unit=100, file='energy_data.txt')
    read(100,*) epsilon
    read(100,*) sigma
    read(100,*) rc_lj
    read(100,*) rv_lj
    read(100,*) rsk_lj
    read(100,*) R0_2
    read(100,*) Kvib
  close(100)

  if (rc_lj>Lx/20) then
    rc_lj = Lx/2
    write(*,*) 'Does not use verlet list!'
    write(*,*) 'rc_lj=', rc_lj
  end if

end subroutine read_energy_parameters


subroutine compute_pressure (pressure)
  !----------------------------------------!
  !input:
  !  pos
  !output:
  !  pressure
  !External Variables:
  !  Ngl, Nml, Npe, NN,
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.52, Eq. (3.4.1).
  !----------------------------------------!
  use global_variables
  implicit none
  real*8, intent(out) :: pressure
  integer i,j,k
  real*8 :: rr, vir, inv_r2, inv_r6, rc_lj2
  real*8, dimension(3) :: rij, fij

  vir = 0
  rc_lj2 = rc_lj * rc_lj
  do i = 1, NN
    do j = i+1, NN
      call rij_and_rr(rij, rr, i, j)
      if (rr<rc_lj2) then
        inv_r2 = sigma*sigma / rr
        inv_r6 = inv_r2*inv_r2*inv_r2
        fij = 48 * epsilon * inv_r2 * inv_r6 * (inv_r6-0.5) * rij
        vir = vir + dot_product(fij,rij)/3
      end if
    end do 
    if ( mod(i,Nml)==1 ) then
      call rij_and_rr(rij, rr, i, i+1)
      fij = Kvib * ( 1 - sqrt(R0_2/rr) ) * rij
      vir = vir + dot_product(fij,rij)/3/2
    elseif ( mod(i,Nml)==0 ) then
      call rij_and_rr(rij, rr, i, i-1)
      fij = Kvib * ( 1 - sqrt(R0_2/rr) ) * rij
      vir = vir + dot_product(fij,rij)/3/2
    else
      call rij_and_rr(rij, rr, i, i+1)
      fij = Kvib * ( 1 - sqrt(R0_2/rr) ) * rij
      vir = vir + dot_product(fij,rij)/3/2
      call rij_and_rr(rij, rr, i, i-1)
      fij = Kvib * ( 1 - sqrt(R0_2/rr) ) * rij
      vir = vir + dot_product(fij,rij)/3/2
    end if
  end do
  pressure = rho / Beta + vir / (Lx*Ly*Lz)

end subroutine compute_pressure


end module compute_energy


















