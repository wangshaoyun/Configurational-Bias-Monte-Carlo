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
  !coulomb
  real*8,  private :: lb          !Bjerrum length
  real*8,  private :: tau_rf      !time ratio of real space and fourier space
  real*8,  private :: alpha       !Ewald screening parameter alpha
  real*8,  private :: alpha2      !alpha2=alpha*alpha
  real*8,  private :: tol
  !
  !real space, cut off at small radius
  real*8,  private :: rcc0         !Cut off radius of real space
  real*8,  private :: rcc02        !rcc2=rcc*rcc  !
  real*8,  private :: clx0         !length of cell in x direction
  real*8,  private :: cly0         !length of cell in y direction
  real*8,  private :: clz0         !length of cell in z direction
  integer, private :: nclx0        !number of cell in x direction
  integer, private :: ncly0        !number of cell in y direction
  integer, private :: nclz0        !number of cell in z direction 
  !
  !real space, cut off at large radius
  real*8,  private :: rcc1         !Cut off radius of real space
  real*8,  private :: rcc12        !rcc2=rcc*rcc  !
  real*8,  private :: clx1         !length of cell in x direction
  real*8,  private :: cly1         !length of cell in y direction
  real*8,  private :: clz1         !length of cell in z direction
  integer, private :: nclx1        !number of cell in x direction
  integer, private :: ncly1        !number of cell in y direction
  integer, private :: nclz1        !number of cell in z direction 
  !
  !reciprocal space
  integer, private :: Kmax1       !max wave number of x direction
  integer, private :: Kmax2       !max wave number of y direction 
  integer, private :: Kmax3       !max wave number of z direction 
  integer, private :: K_total     !Total wave number in reciprocal space
  !
  !lj
  real*8,  private :: epsilon     !Energy unit epsilon in lj potential
  real*8,  private :: sigma       !Distance sigma in lj potential
  real*8,  private :: rcl         !Cut off radius of LJ potential
  real*8,  private :: rcl2        !rcc2=rcc*rcc  !
  real*8,  private :: clx         !length of cell in x direction
  real*8,  private :: cly         !length of cell in y direction
  real*8,  private :: clz         !length of cell in z direction
  integer, private :: nclx        !number of cell in x direction
  integer, private :: ncly        !number of cell in y direction
  integer, private :: nclz        !number of cell in z direction 


!##########################arrays##########################!
  !
  !charge number to monomer number        
  integer, allocatable, dimension( : )          :: charge
  !  
  !monomer number to charge number
  integer, allocatable, dimension( : ), private :: inv_charge
  !
  !cell list of charge
  integer, allocatable, dimension( : ), private :: cell_list_q
  !
  !inverse cell list of charge
  integer, allocatable, dimension( : ), private :: inv_cell_list_q
  !
  !neighbor cells of the center cell
  integer, allocatable, dimension(:,:,:), private :: cell_near_list_lj
  !
  !cell list in real space
  integer, allocatable, dimension( : ), private :: cell_list_lj
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list_lj
  !
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:), private :: hoc_lj  
  !
  ! head of chains, inverse cell list
  integer, allocatable, dimension(:,:,:), private :: inv_hoc_lj
  !
  !neighbor cells of the center cell
  integer, allocatable, dimension(:,:,:), private :: cell_near_list_r0
  !
  !cell list in real space
  integer, allocatable, dimension( : ), private :: cell_list_r0
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list_r0
  !
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:), private :: hoc_r0     
  !
  ! head of chains, inverse cell list
  integer, allocatable, dimension(:,:,:), private :: inv_hoc_r0
  !
  !neighbor cells of the center cell
  integer, allocatable, dimension(:,:,:), private :: cell_near_list_r1
  !
  !cell list in real space
  integer, allocatable, dimension( : ), private :: cell_list_r1
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list_r1
  !
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:), private :: hoc_r1     
  !
  ! head of chains, inverse cell list
  integer, allocatable, dimension(:,:,:), private :: inv_hoc_r1
  !
  !Coulomb energy of i,j in real space
  real,  allocatable, dimension(:), private :: real_ij 
  !
  !coefficients in Fourier space
  real*8,  allocatable, dimension( : ), private :: exp_ksqr
  !
  !structure factor
  complex(kind=8), allocatable, dimension( : ), private :: rho_k
  !
  !difference of structure factor
  complex(kind=8), allocatable, dimension( : ), private :: delta_rhok
  !
  !difference of structure factor
  complex(kind=8), allocatable, dimension( : ), private :: delta_rhok1
  !
  !difference of structure factor
  complex(kind=8), allocatable, dimension( : ), private :: delta_cosk
  !
  !wave vector ordinal number
  integer, allocatable, dimension(:,:), private :: totk_vectk
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

  !
  !
  call initialize_lj_parameters 
  !
  !Initialize ewald parameters and array allocate.
  call Initialize_ewald_parameters
  !
  !Construct the array totk_vectk(K_total,3), and allocate
  !rho_k(K_total), delta_rhok(K_total).
  call build_totk_vectk
  !
  !Construct the coefficients vector in Fourier space
  call build_exp_ksqr
  !
  !
  call pre_calculate_real_space

end subroutine initialize_energy_parameters


subroutine Initialize_energy_arrays
  use global_variables
  implicit none

end subroutine Initialize_energy_arrays


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
    read(100,*) rcl
    read(100,*) lb
    read(100,*) tol
    read(100,*) tau_rf
  close(100)

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


















