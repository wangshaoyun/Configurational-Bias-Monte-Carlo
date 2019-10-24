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
  !real space, cut off at large radius
  real*8,  private :: rcc         !Cut off radius of real space
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
  integer, allocatable, dimension(:,:,:), private :: cell_near_list_r
  !
  !cell list in real space
  integer, allocatable, dimension( : ), private :: cell_list_r
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list_r
  !
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:), private :: hoc_r
  !
  ! head of chains, inverse cell list
  integer, allocatable, dimension(:,:,:), private :: inv_hoc_r
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
  if ( qq /= 0 ) then
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
  end if

end subroutine initialize_energy_parameters


subroutine Initialize_energy_arrays
  use global_variables
  implicit none

  if ( qq /= 0 ) then
    !
    !Initialize charge with lind list. From this subroutine, pos array is needed.
    call Build_Charge_Ewald
    !
    !Initialize cell list of charge
    call Initialize_cell_list_q_Ewald
    !
    !Initialize real cell list
    call Initialize_real_cell_list_Ewald
    !
    !Construct the structure factor rho_k
    call build_rho_k
  end if

  call write_energy_parameters_Ewald

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


subroutine enerex_short(xt, eni)
  use global_variables
  implicit none
  real*8, intent(out) :: eni
  real*8, dimension(4), intent(in) :: xt
  integer :: icelx, icely, icelz, ncel
  integer :: i, j, k
  real*8 :: rij(3), rr, EE1, EE2
  real*8 :: inv_rr2, inv_rr6, inv_rr12

  EE1 = 0
  if (xt(4)/=0) then
    icelx = int(rr(1)/clx1) + 1
    icely = int(rr(2)/cly1) + 1
    icelz = int(rr(3)/clz1) + 1
    ncel = (icelx-1)*ncly1*nclz1+(icely-1)*nclz1+icelz
    do i = cell_near_list_r(ncel,28,1)
      icelx = cell_near_list_r(ncel,i,1)
      icely = cell_near_list_r(ncel,i,2)
      icelz = cell_near_list_r(ncel,i,3)
      j = hoc_r(icelx,icely,icelz) 
      do while (j/=0) 
        rij = xt(1:3) - pos(k,1:3)
        call periodic_condition(rij)
        rr = sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
        if (rr<rcc) then
          EE1=EE1+pos(j,4)*erfc(alpha * rr) / rr
        end if
        j = cell_list_r(j)
      end do
    end do
    EE1 = EE1 * xt(4)
  end if

  EE2 = 0
  icelx = int(rr(1)/clx1) + 1
  icely = int(rr(2)/cly1) + 1
  icelz = int(rr(3)/clz1) + 1  
  ncel = (icelx-1)*ncly1*nclz1+(icely-1)*nclz1+icelz
  do i = cell_near_list_lj(ncel,28,1)
    icelx = cell_near_list_lj(ncel,i,1)
    icely = cell_near_list_lj(ncel,i,2)
    icelz = cell_near_list_lj(ncel,i,3)
    j = hoc_lj(icelx,icely,icelz) 
    do while (j/=0) 
      rij = xt(1:3) - pos(k,1:3)
      call periodic_condition(rij)
      rr = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
      if ( rr < rc_lj2 ) then
        inv_rr2  = sigma2 / rr
        inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
        inv_rr12 = inv_rr6 * inv_rr6
        EE2      = EE2 + inv_rr12 - inv_rr6 + 0.25D0
      end if
      j = cell_list_r(j)
    end do
  end do
  EE2 = EE2 * 4 * epsilon * EE2

  eni = EE1 + EE2

  rr = xt(1)*xt(1) + xt(2)*xt(2)
  if (rr<r_cy) then
    eni = eni + 1e10
  end if

end subroutine enerex_short


subroutine energy_long(del_E)
  use global_variables
  implicit none
  real*8, intent(out) :: del_E
  real*8  :: Del_Recip_erg
  complex(kind=8) :: eikx0(num_newconf, -Kmax1:Kmax1 )
  complex(kind=8) :: eikx1(num_newconf, -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0(num_newconf, -Kmax2:Kmax2 )
  complex(kind=8) :: eiky1(num_newconf, -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0(num_newconf, -Kmax3:Kmax3 )
  complex(kind=8) :: eikz1(num_newconf, -Kmax3:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: c1, c2, c3
  integer :: ord(3), i, m, p, q, r

  c1 = 2*pi/Lx
  c2 = 2*pi/Ly
  c3 = 2*pi/Lz/Z_empty

  do i = 1, num_newconf

    m = Nml - num_newconf + i 

    eikx0(i,0)  = (1,0)
    eiky0(i,0)  = (1,0)
    eikz0(i,0)  = (1,0)

    eikx0(i,1)  = cmplx( cos(c1*pos_old(m,1)), sin(c1*pos_old(m,1)), 8 )
    eiky0(i,1)  = cmplx( cos(c2*pos_old(m,2)), sin(c2*pos_old(m,2)), 8 )
    eikz0(i,1)  = cmplx( cos(c3*pos_old(m,3)), sin(c3*pos_old(m,3)), 8 )

    eikx0(i,-1) = conjg(eikx0(i,1))
    eiky0(i,-1) = conjg(eiky0(i,1))
    eikz0(i,-1) = conjg(eikz0(i,1))

    eikx1(i,0)  = (1,0)
    eiky1(i,0)  = (1,0)
    eikz1(i,0)  = (1,0)

    eikx1(i,1)  = cmplx( cos(c1*pos_new(m,1)), sin(c1*pos_new(m,1)), 8 )
    eiky1(i,1)  = cmplx( cos(c2*pos_new(m,2)), sin(c2*pos_new(m,2)), 8 )
    eikz1(i,1)  = cmplx( cos(c3*pos_new(m,3)), sin(c3*pos_new(m,3)), 8 )

    eikx1(i,-1) = conjg(eikx1(i,1))
    eiky1(i,-1) = conjg(eiky1(i,1))
    eikz1(i,-1) = conjg(eikz1(i,1))

  end do    

  do p=2, Kmax1
    do m=1, num_newconf
      eikx0(m,p)=eikx0(m,p-1)*eikx0(m,1)
      eikx0(m,-p)=conjg(eikx0(m,p))
      eikx1(m,p)=eikx1(m,p-1)*eikx1(m,1)
      eikx1(m,-p)=conjg(eikx1(m,p))
    end do
  end do
  do q=2, Kmax2
    do m=1, num_newconf
      eiky0(m,q)=eiky0(m,q-1)*eiky0(m,1)
      eiky0(m,-q)=conjg(eiky0(m,q))
      eiky1(m,q)=eiky1(m,q-1)*eiky1(m,1)
      eiky1(m,-q)=conjg(eiky1(m,q))
    end do
  end do
  do r=2, Kmax3
    do m=1, num_newconf
      eikz0(m,r)=eikz0(m,r-1)*eikz0(m,1)
      eikz0(m,r)=conjg(eikz0(m,q))
      eikz1(m,r)=eikz1(m,r-1)*eikz1(m,1)
      eikz1(m,r)=conjg(eikz1(m,q))
    end do
  end do

  do i = 1, K_total
    ord = totk_vectk(i,:)
    do m = 1, num_newconf
      rho_k(i) = rho_k(i) + &
                 zq(m) * eikx(m,ord(1)) * eiky(m,ord(2)) * eikz(m,ord(3))
    end do
  end do

end subroutine energy_long



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
    read(100,*) rcc0
    read(100,*) lb
    read(100,*) tol
    read(100,*) tau_rf
  close(100)

end subroutine read_energy_parameters


subroutine initialize_lj_parameters
  use global_variables
  implicit none

  rcl2 = rcl*rcl
  nclx = int(Lx/(rcl+1))
  ncly = int(Ly/(rcl+1))
  nclz = int(Lz/(rcl+1))

  clx = Lx/nclx
  cly = Ly/ncly
  clz = Lz/nclz

end subroutine initialize_lj_parameters


subroutine Initialize_ewald_parameters
  use global_variables
  implicit none

  alpha    = ( tau_rf * pi**3 * Nq / (Lx*Ly*Lz*Z_empty)**2 ) ** (1.D0/6)
  alpha2   = alpha * alpha
  rcc1  = tol / alpha
  rcc12 = rcc1 * rcc1
  if (rcc1<min(Lx/3,Lz/3)) then
    !
    !use verlet list in real space
    Kmax1 = ceiling(tol*Lx*alpha/pi)
    Kmax2 = ceiling(tol*Ly*alpha/pi)
    Kmax3 = ceiling(tol*Lz*alpha/pi)
  else
    rcc1 = min(Lx/3,Lz/3)
    Kmax1    = ceiling(tol*tol/pi*Lx/rcc)
    Kmax2    = ceiling(tol*tol/pi*Ly/rcc)
    Kmax3    = ceiling(tol*tol/pi*Lz/rcc)
    alpha    = tol / rcc1
    alpha2   = alpha * alpha
    rcc12 = rcc1 * rcc1
  end if
  !
  !Cell list parameters
  nclx1 = int(Lx/rcc1)     !cell numbers in x direction
  ncly1 = int(Ly/rcc1)
  nclz1 = int(Lz/rcc1)
  clx1 = Lx/nclx1         !cell length    
  cly1 = Ly/ncly1
  clz1 = Lz/nclz1

  rcc02 = rcc0*rcc0

  nclx0 = int(Lx/rcc0)    !cell numbers in x direction
  ncly0 = int(Ly/rcc0)
  nclz0 = int(Lz/rcc0)
  clx0 = Lx/nclx1         !cell length    
  cly0 = Ly/ncly1
  clz0 = Lz/nclz1

end subroutine Initialize_ewald_parameters


subroutine pre_calculate_real_space
  use global_variables
  implicit none
  integer :: n=200000, i
  real*8 :: del_r, rr

  if (allocated(real_ij0)) deallocate(real_ij0)
  if (allocated(real_ij1)) deallocate(real_ij1)
  allocate( real_ij0(n) )
  allocate( real_ij1(n) )

  del_r = rcc1/n
  do i = 1, n
    rr = del_r*i
    real_ij1(i) = erfc(alpha*rr)/rr 
  end do

  del_r = rcc0/n
  do i = 1, n
    rr = del_r*i
    real_ij0(i) = erfc(alpha*rr)/rr 
  end do

end subroutine pre_calculate_real_space


subroutine build_totk_vectk
  !--------------------------------------!
  !exp_ksqr, rho_k, delta_rhok are all vectors with size of
  !K_total. For i = 1 to K_total, we often need to know 
  !corresponding wave number kx,ky,kz. This progam build a 
  !array totk_vectk(1:K_total,3) to store kx,ky,kz.
  !What's more, rho_k and delta_rhok are allocated here.
  !   
  !Input
  !   
  !Output
  !   totk_vectk
  !   K_total
  !External Variables
  !   K_total, Kmax1, Kmax2, Kmax3
  !   totk_vectk
  !Routine Referenced:
  !   
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l
  real*8  :: ksqr, k1, k2, k3, factor, kcut

  K_total=0
  do k = 0, Kmax3
    do i = -Kmax1, Kmax1
      do j = -Kmax2, Kmax2
        kcut = (1.D0*i/Kmax1) * (1.D0*i/Kmax1) &
             + (1.D0*j/Kmax2) * (1.D0*j/Kmax2) &
             + (1.D0*k/Kmax3) * (1.D0*k/Kmax3)
        if ( kcut>1 .or. kcut==0 ) cycle
        K_total = K_total + 1
      end do
    end do
  end do

  if ( allocated(totk_vectk) ) deallocate(totk_vectk)
  if ( allocated(rho_k)      ) deallocate(rho_k)
  if ( allocated(delta_rhok) ) deallocate(delta_rhok)
  if ( allocated(delta_rhok1) ) deallocate(delta_rhok1)
  if ( allocated(delta_cosk) ) deallocate(delta_cosk)
  allocate( totk_vectk( K_total, 3 ) )
  allocate( rho_k( K_total )         )
  allocate( delta_rhok( K_total )    )
  allocate( delta_rhok1( K_total )   )
  allocate( delta_cosk( K_total )    )
  totk_vectk = 0
  rho_k      = 0
  delta_rhok = 0
  delta_rhok1 = 0
  delta_cosk = 0

  l=0
  do k = 0, Kmax3
    do i = -Kmax1, Kmax1
      do j = -Kmax2, Kmax2
        kcut = (1.D0*i/Kmax1) * (1.D0*i/Kmax1) &
             + (1.D0*j/Kmax2) * (1.D0*j/Kmax2) &
             + (1.D0*k/Kmax3) * (1.D0*k/Kmax3)
        if ( kcut>1 .or. kcut==0 ) cycle
        l = l + 1
        totk_vectk( l, 1 ) = i
        totk_vectk( l, 2 ) = j
        totk_vectk( l, 3 ) = k
      end do
    end do
  end do

end subroutine build_totk_vectk


subroutine build_exp_ksqr
  !--------------------------------------!
  !Reciprocal energy is divided to three parts: 
  !1.structrure factor is referred to rho_k.
  !2.difference of structure factor between new and old
  !position is referred to delta_rhok.
  !3.the other which includes exp(k^2/4/alpha) is referred 
  !to exp_ksqr.
  !This program is used to bulid the third part.
  !
  !Input
  !   
  !Output
  !   exp_ksqr
  !External Variables
  !   K_total
  !   Kmax1, Kmax2, Kmax3
  !   alpha2, lb
  !   Lx, Ly, Lz, Z_empty
  !   Beta
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.300(12.1.25),
  !however his alpha is alpha^2 in this program.b 
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, l, ord(3)
  real*8  :: ksqr, k1, k2, k3, factor

  if ( allocated(exp_ksqr) ) deallocate(exp_ksqr)
  allocate( exp_ksqr(K_total) )
  exp_ksqr = 0

  l = 0
  do i = 1, K_total
    ord = totk_vectk(i,:)
    if ( ord(3) == 0 ) then
      factor = 1
    else
      factor = 2
    end if
    k1   = 2*pi*ord(1) / Lx
    k2   = 2*pi*ord(2) / Ly
    k3   = 2*pi*ord(3) / Lz
    ksqr = k1*k1 + k2*k2 + k3*k3 
    exp_ksqr(i) = factor * 4*pi / (Lx*Ly*Lz) *  &
                  exp(-ksqr/4/alpha2) / ksqr * lb / Beta     
  end do

end subroutine build_exp_ksqr


subroutine build_rho_k
  !--------------------------------------!
  !Calculate the structure factor array.
  !   
  !Input
  !   
  !Output
  !   rho_k
  !External Variables
  !   pos, charge
  !   Nq, Lx, Ly, Lz, Z_empty, K_total
  !Routine Referenced:
  !1.
  !--------------------------------------!
  use global_variables
  implicit none
  complex(kind=8) :: eikx(1:Nq, -Kmax1:Kmax1)
  complex(kind=8) :: eiky(1:Nq, -Kmax2:Kmax2)
  complex(kind=8) :: eikz(1:Nq, 0:Kmax3)
  integer i,j,l,m,n,p,q,r,ord(3)
  real*8 :: c1, c2, c3
  real*8 :: zq(Nq)
  rho_k = 0
  zq = 0
  eikx = 0
  eiky = 0
  eikz = 0

  m = cell_list_q(Nq+1)
  do while( cell_list_q(m)/=0 )
    zq(m) = pos(charge(m),4)
    m = cell_list_q(m)
  end do

  c1 = 2*pi/Lx
  c2 = 2*pi/Ly
  c3 = 2*pi/Lz
  m = cell_list_q(Nq+1)
  do while(m /= 0)
    i = charge(m)
    eikx(m,0)  = (1,0)
    eiky(m,0)  = (1,0)
    eikz(m,0)  = (1,0)

    eikx(m,1)  = cmplx( cos(c1*pos(i,1) ), sin( c1*pos(i,1) ), 8 )
    eiky(m,1)  = cmplx( cos(c2*pos(i,2) ), sin( c2*pos(i,2) ), 8 )
    eikz(m,1)  = cmplx( cos(c3*pos(i,3) ), sin( c3*pos(i,3) ), 8 )

    eikx(m,-1) = conjg(eikx(m,1))
    eiky(m,-1) = conjg(eiky(m,1))
    m = cell_list_q(m)
  end do

  do p=2, Kmax1
    m = cell_list_q(Nq+1)
    do while(m/=0)
      eikx(m,p)=eikx(m,p-1)*eikx(m,1)
      eikx(m,-p)=conjg(eikx(m,p))
      m = cell_list_q(m)
    end do
  end do
  do q=2, Kmax2
    m = cell_list_q(Nq+1)
    do while(m/=0)
      eiky(m,q)=eiky(m,q-1)*eiky(m,1)
      eiky(m,-q)=conjg(eiky(m,q))
      m = cell_list_q(m)
    end do
  end do
  do r=2, Kmax3
    m = cell_list_q(Nq+1)
    do while(m/=0)
      eikz(m,r)=eikz(m,r-1)*eikz(m,1)
      m = cell_list_q(m)
    end do
  end do

  do i = 1, K_total
    ord = totk_vectk(i,:)
    m = cell_list_q(Nq+1)
    do while(m/=0)
      rho_k(i) = rho_k(i) + &
                 zq(m) * eikx(m,ord(1)) * eiky(m,ord(2)) * eikz(m,ord(3))
      m = cell_list_q(m)
    end do
  end do

end subroutine build_rho_k


subroutine delete_cell_list_real(np,rr)
  use global_variables
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr
  integer :: icelx,icely,icelz
  integer :: nti,bfi  

  icelx = int(rr(1)/clx)+1
  icely = int(rr(2)/cly)+1
  icelz = int(rr(3)/clz)+1     

  bfi = cell_list_r(np)
  nti = inv_cell_list_r(np)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list_r(nti) = bfi
    inv_cell_list_r(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list_r(nti) = bfi
    inv_hoc_r(icelx,icely,icelz) = nti
  elseif ( bfi/=0 .and. nti==0 ) then    !the last one
    hoc_r(icelx,icely,icelz) = bfi
    inv_cell_list_r(bfi) = nti
  else                                   !only one
    hoc_r(icelx,icely,icelz) = nti
    inv_hoc_r(icelx,icely,icelz) = bfi
  end if

end subroutine delete_cell_list_real


subroutine add_cell_list_r(np,rr)
  use global_variables
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr
  integer :: icelx,icely,icelz

  icelx = int(rr(1)/clx)+1
  icely = int(rr(2)/cly)+1
  icelz = int(rr(3)/clz)+1 

  inv_cell_list_r(np) = 0
  if ( inv_hoc_r(icelx,icely,icelz) /=0 ) then
    inv_cell_list_r( hoc_r(icelx,icely,icelz) ) = np
  else
    inv_hoc_r(icelx,icely,icelz) = np
  end if

  cell_list_r(ii) = hoc_r(icelx,icely,icelz)
  hoc_r(icelx,icely,icelz) = np

end subroutine add_cell_list_r


subroutine delete_cell_list_lj(np,rr)
  use global_variables
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr
  integer :: icelx,icely,icelz

  icelx = int(rr(1)/clx)+1
  icely = int(rr(2)/cly)+1
  icelz = int(rr(3)/clz)+1 

  bfi = cell_list_lj(np)
  nti = inv_cell_list_lj(np)

  if ( bfi/=0 .and. nti/=0 ) then        !middle
    cell_list_lj(nti) = bfi
    inv_cell_list_lj(bfi) = nti
  elseif ( bfi==0 .and. nti/=0 ) then    !the first one
    cell_list_lj(nti) = bfi
    inv_hoc_lj(icelx,icely,icelz) = nti
  elseif ( bfi/=0 .and. nti==0 ) then    !the last one
    hoc_lj(icelx,icely,icelz) = bfi
    inv_cell_list_lj(bfi) = nti
  else                                   !only one
    hoc_lj(icelx,icely,icelz) = nti
    inv_hoc_lj(icelx,icely,icelz) = bfi
  end if

end subroutine delete_cell_list_lj


subroutine add_cell_list_lj(np,rr)
  use global_variables
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr
  integer :: icelx,icely,icelz

  icelx = int(rr(1)/clx)+1
  icely = int(rr(2)/cly)+1
  icelz = int(rr(3)/clz)+1 

  inv_cell_list_lj(np) = 0
  if ( inv_hoc_lj(icelx,icely,icelz) /=0 ) then
    inv_cell_list_lj( hoc_r(icelx,icely,icelz) ) = np
  else
    inv_hoc_r(icelx,icely,icelz) = np
  end if

  cell_list_r(ii) = hoc_r(icelx,icely,icelz)
  hoc_r(icelx,icely,icelz) = np

end subroutine add_cell_list_lj




end module compute_energy


















