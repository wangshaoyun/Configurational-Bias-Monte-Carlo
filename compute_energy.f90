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
  real*8,  private :: rcc          !Cut off radius of real space
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
  integer, private :: nclx        !number of cell in x direction
  integer, private :: ncly        !number of cell in y direction
  integer, private :: nclz        !number of cell in z direction 
  real*8 :: clx                   !length of cell in x direction
  real*8 :: cly                   !length of cell in y direction
  real*8 :: clz                   !length of cell in z direction


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
  integer, allocatable, dimension( : ) :: cell_list_lj
  !
  !inverse cell list in real space
  integer, allocatable, dimension( : ), private :: inv_cell_list_lj
  !
  ! head of chains, cell list
  integer, allocatable, dimension(:,:,:) :: hoc_lj  
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
  !
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
    call Initialize_lj_cell_list
    !
    !Initialize real cell list
    call Initialize_real_cell_list
    !
    !
    call Initialize_cell_list_q
    !
    !Construct the structure factor rho_k
    call build_rho_k
  end if

  call write_energy_parameters_Ewald

end subroutine Initialize_energy_arrays


subroutine error_analysis(EE1)
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: EE0,tol1, rmse
  real*8, intent(out) :: EE1


  tol1 = tol
  tol = 5                
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
  !Initialize real cell list
  call Initialize_real_cell_list
  !
  !Construct the structure factor rho_k
  call build_rho_k
  !
  !
  call total_energy(EE0)

  tol=tol1
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
  !Initialize real cell list
  call Initialize_real_cell_list
  !
  !Construct the structure factor rho_k
  call build_rho_k
  !
  !
  call total_energy(EE1)

  rmse = abs(EE1-EE0)/EE0

  write(*,*) rmse

end subroutine error_analysis


subroutine total_energy(EE)
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

  call Coulomb_Energy(EE)

end subroutine total_energy


subroutine LJ_Energy (EE)
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: EE
  integer :: i, j, k
  integer :: icelx, icely, icelz, ncel
  real*8  :: rr, rij(3), inv_rr2, inv_rr6, inv_rr12, sigma2

  EE = 0
  do i = 1, NN
    icelx = int(pos(i,1)/clx) + 1
    icely = int(pos(i,2)/cly) + 1
    icelz = int(pos(i,3)/clz) + 1  
    ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz 
    do j = 1, 27
      icelx = cell_near_list_lj(ncel,i,1)
      icely = cell_near_list_lj(ncel,i,2)
      icelz = cell_near_list_lj(ncel,i,3)
      k = hoc_lj(icelx,icely,icelz)    
      do while(k/=0)
        if (k==i) cycle
        call rij_and_rr(rij,rr,k,i)
        if ( rr < rcl2 ) then
          inv_rr2  = sigma2 / rr
          inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
          inv_rr12 = inv_rr6 * inv_rr6
          EE       = EE + inv_rr12 - inv_rr6 + 0.25D0
        end if
        k = cell_list_lj(k)
      end do
    end do   
  end do

  EE = EE * 4 * epsilon

end subroutine LJ_energy


subroutine Coulomb_Energy(EE)
  use global_variables
  implicit none
  real*8, intent(out) :: EE 
  integer :: i, j, k
  integer :: icelx, icely, icelz, ncel
  real*8  :: rr, rij(3), EE1

  do i = 1, NN
    if (pos(i,4)/=0) then
      icelx = int(pos(i,1)/clx) + 1
      icely = int(pos(i,2)/cly) + 1
      icelz = int(pos(i,3)/clz) + 1  
      ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz 
      do j = 1, 27
        icelx = cell_near_list_r(ncel,i,1)
        icely = cell_near_list_r(ncel,i,2)
        icelz = cell_near_list_r(ncel,i,3)
        k = hoc_r(icelx,icely,icelz)    
        do while(k/=0)
          if (k==i) cycle
          call rij_and_rr(rij,rr,k,i)
          rr = sqrt(rr)
          if (rr<rcc) then
            EE1 = EE1 + pos(k,4) * erfc(alpha * rr) / rr
          end if
          k = cell_list_r(k)
        end do
      end do 
    end if
    EE = EE + EE1 * pos(i,4)
  end do
  !
  !fourier space
  EE = EE + EE/Beta*lb + sum( exp_ksqr * real( conjg(rho_k) * rho_k ) )/2.D0 

end subroutine Coulomb_Energy


subroutine energy_short(xt, eni)
  use global_variables
  implicit none
  real*8, intent(out) :: eni
  real*8, dimension(4), intent(in) :: xt
  integer :: icelx, icely, icelz, ncel
  integer :: i, j, k
  real*8 :: rij(3), rr, EE1, EE2, sigma2
  real*8 :: inv_rr2, inv_rr6, inv_rr12

  EE1 = 0
  if (xt(4)/=0) then
    icelx = int(xt(1)/clx1) + 1
    icely = int(xt(2)/cly1) + 1
    icelz = int(xt(3)/clz1) + 1
    ncel = (icelx-1)*ncly1*nclz1+(icely-1)*nclz1+icelz
    do i = 1, 27
      icelx = cell_near_list_r(ncel,i,1)
      icely = cell_near_list_r(ncel,i,2)
      icelz = cell_near_list_r(ncel,i,3)
      j = hoc_r(icelx,icely,icelz) 
      do while (j/=0)
        call Lagrange_to_Euler(j,k) 
        rij = xt(1:3) - pos(k,1:3)
        call periodic_condition(rij)
        rr = sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
        if (rr<rcc) then
          EE1 = EE1 + pos(j,4)*erfc(alpha * rr) / rr
        end if
        j = cell_list_r(j)
      end do
    end do
    EE1 = EE1 * xt(4)
  end if

  EE2 = 0
  icelx = int(xt(1)/clx) + 1
  icely = int(xt(2)/cly) + 1
  icelz = int(xt(3)/clz) + 1  
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, 27
    icelx = cell_near_list_lj(ncel,i,1)
    icely = cell_near_list_lj(ncel,i,2)
    icelz = cell_near_list_lj(ncel,i,3)
    j = hoc_lj(icelx,icely,icelz) 
    do while (j/=0) 
      call Lagrange_to_Euler(j,k)
      rij = xt(1:3) - pos(k,1:3)
      call periodic_condition(rij)
      rr = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
      if ( rr < rcl2 ) then
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

end subroutine energy_short


subroutine energy_long(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  complex(kind=8) :: eikx0(num_newconf, -Kmax1:Kmax1 )
  complex(kind=8) :: eikx1(num_newconf, -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0(num_newconf, -Kmax2:Kmax2 )
  complex(kind=8) :: eiky1(num_newconf, -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0(num_newconf, 0:Kmax3 )
  complex(kind=8) :: eikz1(num_newconf, 0:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: c1, c2, c3
  integer :: ord(3), i, m, n1, n2, p, q, r
  integer, dimension(:), allocatable :: charge1
  integer, dimension(:), allocatable :: charge2

  n1 = 0
  n2 = 0
  do i = 1, num_newconf
    m = Nml - num_newconf + i
    if (pos_old(m,4)/=0) then
      n1 = n1 + 1
    end if
    if (pos_new(m,4)/=0) then
      n2 = n2 + 1
    end if
  end do
  allocate(charge1(n1))
  allocate(charge2(n2))
  n1 = 0
  n2 = 0
  do i = 1, num_newconf
    m = Nml - num_newconf + i
    if (pos_old(m,4)/=0) then
      n1 = n1 + 1
      charge1(n1) = m
    end if
    if (pos_new(m,4)/=0) then
      n2 = n2 + 1
      charge2(n2) = m
    end if
  end do

  c1 = 2*pi/Lx
  c2 = 2*pi/Ly
  c3 = 2*pi/Lz

  do i = 1, n1

    m = charge1(i)

    eikx0(i,0)  = (1,0)
    eiky0(i,0)  = (1,0)
    eikz0(i,0)  = (1,0)

    eikx0(i,1)  = cmplx( cos(c1*pos_old(m,1)), sin(c1*pos_old(m,1)), 8 )
    eiky0(i,1)  = cmplx( cos(c2*pos_old(m,2)), sin(c2*pos_old(m,2)), 8 )
    eikz0(i,1)  = cmplx( cos(c3*pos_old(m,3)), sin(c3*pos_old(m,3)), 8 )

    eikx0(i,-1) = conjg(eikx0(i,1))
    eiky0(i,-1) = conjg(eiky0(i,1))

    m = charge2(i)

    eikx1(i,0)  = (1,0)
    eiky1(i,0)  = (1,0)
    eikz1(i,0)  = (1,0)

    eikx1(i,1)  = cmplx( cos(c1*pos_new(m,1)), sin(c1*pos_new(m,1)), 8 )
    eiky1(i,1)  = cmplx( cos(c2*pos_new(m,2)), sin(c2*pos_new(m,2)), 8 )
    eikz1(i,1)  = cmplx( cos(c3*pos_new(m,3)), sin(c3*pos_new(m,3)), 8 )

    eikx1(i,-1) = conjg(eikx1(i,1))
    eiky1(i,-1) = conjg(eiky1(i,1))

  end do

  do p=2, Kmax1
    do m=1, n1
      eikx0(m,p)=eikx0(m,p-1)*eikx0(m,1)
      eikx0(m,-p)=conjg(eikx0(m,p))
      eikx1(m,p)=eikx1(m,p-1)*eikx1(m,1)
      eikx1(m,-p)=conjg(eikx1(m,p))
    end do
  end do
  do q=2, Kmax2
    do m=1, n1
      eiky0(m,q)=eiky0(m,q-1)*eiky0(m,1)
      eiky0(m,-q)=conjg(eiky0(m,q))
      eiky1(m,q)=eiky1(m,q-1)*eiky1(m,1)
      eiky1(m,-q)=conjg(eiky1(m,q))
    end do
  end do
  do r=2, Kmax3
    do m=1, n1
      eikz0(m,r)=eikz0(m,r-1)*eikz0(m,1)
      eikz1(m,r)=eikz1(m,r-1)*eikz1(m,1)
    end do
  end do

  delta_rhok = 0
  do i = 1, K_total
    ord = totk_vectk(i,:)
    do m = 1, n1
      delta_rhok(i) = delta_rhok(i)  &
      +(pos_new(charge2(m),4)*eikx1(m,ord(1))*eiky1(m,ord(2))*eikz1(m,ord(3)) &
      - pos_old(charge1(m),4)*eikx0(m,ord(1))*eiky0(m,ord(2))*eikz0(m,ord(3)))
    end do
  end do

  DeltaE = sum( exp_ksqr * ( 2*Real(rho_k*delta_rhok) &
                      + conjg(delta_rhok)*delta_rhok ) )

  deallocate(charge1)
  deallocate(charge2)

end subroutine energy_long


subroutine delta_energy_ions_move(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  real*8  :: Del_Recip_erg
  complex(kind=8) :: eikx0( -Kmax1:Kmax1 ), eikx1( -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0( -Kmax2:Kmax2 ), eiky1( -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0( -Kmax3:Kmax3 ), eikz1( -Kmax3:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: rij(3), rr
  real*8  :: inv_rr2, inv_rr6, inv_rr12
  real*8  :: c1, c2, c3
  integer :: i, j, k, ord(3), p, q, r
  integer :: icelx,icely,icelz,ncel
  
  DeltaE = 0

  call Delta_lj_Energy(.true., ip, pos_ip1, DeltaE)

  call Delta_lj_Energy(.false., ip, pos_ip0, DeltaE)

  call Delta_real_energy(.true., ip, pos_ip1, DeltaE)

  call Delta_real_energy(.false., ip, pos_ip0, DeltaE)

  call Delta_Reciprocal_Energy(pos_ip0, pos_ip1, DeltaE)

end subroutine delta_energy_ions_move


subroutine Delta_energy_pH(DeltaE)
  use global_variables
  implicit none
  real*8, intent(out) :: DeltaE
  real*8 :: EE1, EE2

  DeltaE = 0
  if ( pos_ip0(4) == 0 ) then    !add

    call Delta_lj_energy(.true., ipi, pos_ipi1, DeltaE)

    call Delta_real_energy(.true., ipi, pos_ipi1, DeltaE)

    call Delta_real_energy(.true., ip, pos_ip1, DeltaE)

    call Delta_Reciprocal_Energy_pH(.true., pos_ip1, pos_ipi1, DeltaE)

  else        !delete

    call Delta_lj_energy(.false., ipi, pos_ipi0, DeltaE)

    call Delta_real_energy(.false., ipi, pos_ipi0, DeltaE)

    call Delta_real_energy(.false., ip, pos_ip1, DeltaE)

    call Delta_Reciprocal_Energy_pH(.false., pos_ip0, pos_ipi0, DeltaE)

  end if

end subroutine Delta_energy_pH


subroutine Delta_lj_Energy(lg, np, ri, DeltaE)
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8, dimension(4), intent(in) :: ri
  integer, intent(in) :: np
  logical, intent(in) :: lg
  real*8  :: EE, sigma2
  real*8  :: rij(3), rr, inv_rr2, inv_rr6, inv_rr12
  integer :: i, j, k
  integer :: icelx, icely, icelz, ncel

  EE     = 0
  sigma2 = sigma * sigma

  icelx = int(ri(1)/clx) + 1
  icely = int(ri(2)/cly) + 1
  icelz = int(ri(3)/clz) + 1  
  ncel = (icelx-1)*ncly*nclz+(icely-1)*nclz+icelz
  do i = 1, 27
    icelx = cell_near_list_lj(ncel,i,1)
    icely = cell_near_list_lj(ncel,i,2)
    icelz = cell_near_list_lj(ncel,i,3)
    j = hoc_lj(icelx,icely,icelz) 
    do while (j/=0) 
      call Lagrange_to_Euler(j,k)
      rij = ri(1:3) - pos(k,1:3)
      call periodic_condition(rij)
      rr = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
      if ( rr < rcl2 ) then
        inv_rr2  = sigma2 / rr
        inv_rr6  = inv_rr2 * inv_rr2 * inv_rr2
        inv_rr12 = inv_rr6 * inv_rr6
        EE       = EE + inv_rr12 - inv_rr6 + 0.25D0
      end if
      j = cell_list_r(j)
    end do
  end do

  if (lg) then
    DeltaE = DeltaE + 4 * epsilon * EE
  else
    DeltaE = DeltaE - 4 * epsilon * EE
  end if

end subroutine Delta_lj_Energy


subroutine Delta_real_energy(lg, np, ri, DeltaE)
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8, dimension(4), intent(in) :: ri
  integer, intent(in) :: np
  logical, intent(in) :: lg  
  real*8  :: rij(3), rr, EE
  integer :: i, j, k
  integer :: icelx, icely, icelz, ncel

  EE = 0
  icelx = int(ri(1)/clx1) + 1
  icely = int(ri(2)/cly1) + 1
  icelz = int(ri(3)/clz1) + 1
  ncel = (icelx-1)*ncly1*nclz1+(icely-1)*nclz1+icelz
  do i = 1, 27
    icelx = cell_near_list_r(ncel,i,1)
    icely = cell_near_list_r(ncel,i,2)
    icelz = cell_near_list_r(ncel,i,3)
    j = hoc_r(icelx,icely,icelz) 
    do while (j/=0) 
      call Lagrange_to_Euler(j,k)
      rij = ri(1:3) - pos(k,1:3)
      call periodic_condition(rij)
      rr = sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))
      if (rr<rcc) then
        EE = EE + pos(j,4)*erfc(alpha * rr) / rr
      end if
      j = cell_list_r(j)
    end do
  end do
  EE = EE * ri(4)

  if (lg) then
    DeltaE = DeltaE + EE + ri(4)*ri(4)*lb/beta
  else 
    DeltaE = DeltaE - EE - ri(4)*ri(4)*lb/beta
  end if

end subroutine Delta_real_energy


subroutine Delta_Reciprocal_Energy(r1, r2, DeltaE)
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8, dimension(4), intent(in) :: r1
  real*8, dimension(4), intent(in) :: r2
  real*8  :: Del_Recip_erg
  complex(kind=8) :: eikx0( -Kmax1:Kmax1 ), eikx1( -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0( -Kmax2:Kmax2 ), eiky1( -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0( -Kmax3:Kmax3 ), eikz1( -Kmax3:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: c1, c2, c3
  integer :: ord(3), p, q, r, i

  c1 = 2*pi / Lx
  c2 = 2*pi / Ly
  c3 = 2*pi / Lz

  eikx0(0)  = ( 1,0 )
  eiky0(0)  = ( 1,0 )
  eikz0(0)  = ( 1,0 )
  eikx0(1)  = cmplx( cos( c1 * r1(1) ), sin( -c1 * r1(1) ), 8 )
  eiky0(1)  = cmplx( cos( c2 * r1(2) ), sin( -c2 * r1(2) ), 8 )
  eikz0(1)  = cmplx( cos( c3 * r1(3) ), sin( -c3 * r1(3) ), 8 )
  eikx0(-1) = conjg( eikx0(1) )
  eiky0(-1) = conjg( eiky0(1) )
  eiky0(-1) = conjg( eiky0(1) )

  do p = 2, Kmax1
    eikx0(p)  = eikx0(p-1) * eikx0(1)
    eikx0(-p) = conjg( eikx0(p) )
  end do
  do q = 2, Kmax2
    eiky0(q)  = eiky0(q-1) * eiky0(1)
    eiky0(-q) = conjg(eiky0(q))
  end do
  do r = 2, Kmax3
    eikz0(r)  = eikz0(r-1) * eikz0(1)
  end do

  eikx1(0)  = ( 1,0 )
  eiky1(0)  = ( 1,0 )
  eikz1(0)  = ( 1,0 )
  eikx1(1)  = cmplx( cos( c1 * r2(1) ), sin( -c1 * r2(1) ), 8 )
  eiky1(1)  = cmplx( cos( c2 * r2(2) ), sin( -c2 * r2(2) ), 8 )
  eikz1(1)  = cmplx( cos( c3 * r2(3) ), sin( -c3 * r2(3) ), 8 )
  eikx1(-1) = conjg( eikx1(1) )
  eiky1(-1) = conjg( eiky1(1) )
  eikz1(-1) = conjg( eikz1(1) )

  do p=2, Kmax1
    eikx1(p)  = eikx1(p-1) * eikx1(1)
    eikx1(-p) = conjg( eikx1(p) )
  end do
  do q=2, Kmax2
    eiky1(q)  = eiky1(q-1) * eiky1(1)
    eiky1(-q) = conjg(eiky1(q))
  end do
  do r=2, Kmax3
    eikz1(r)  = eikz1(r-1) * eikz1(1)
  end do

  do i=1, K_total
    ord = totk_vectk(i,1:3)
    eikr0 = eikx0(ord(1)) * eiky0(ord(2)) * eikz0(ord(3))
    eikr1 = eikx1(ord(1)) * eiky1(ord(2)) * eikz1(ord(3))
    delta_rhok(i) = eikr1 - eikr0
    delta_cosk(i) = 1 - real( conjg(eikr1) * eikr0 )
  end do

  delta_rhok = delta_rhok * r1(4)

  delta_cosk = delta_cosk * ( r1(4) * r1(4) )

  Del_Recip_erg = sum( exp_ksqr * ( Real( rho_k * delta_rhok ) + delta_cosk ) )

  DeltaE = DeltaE + Del_Recip_erg

end subroutine Delta_Reciprocal_Energy


subroutine Delta_Reciprocal_Energy_pH(lg, r1, r2, DeltaE)
  use global_variables
  implicit none
  real*8, intent(inout) :: DeltaE
  real*8, dimension(4), intent(in) :: r1
  real*8, dimension(4), intent(in) :: r2
  logical, intent(in) :: lg   
  real*8  :: Del_Recip_erg
  complex(kind=8) :: eikx0( -Kmax1:Kmax1 ), eikx1( -Kmax1:Kmax1 )
  complex(kind=8) :: eiky0( -Kmax2:Kmax2 ), eiky1( -Kmax2:Kmax2 )
  complex(kind=8) :: eikz0( -Kmax3:Kmax3 ), eikz1( -Kmax3:Kmax3 )
  complex(kind=8) :: eikr0, eikr1
  real*8  :: c1, c2, c3
  integer :: ord(3), p, q, r, i

  c1 = 2*pi / Lx
  c2 = 2*pi / Ly
  c3 = 2*pi / Lz

  eikx0(0)  = ( 1,0 )
  eiky0(0)  = ( 1,0 )
  eikz0(0)  = ( 1,0 )
  eikx0(1)  = cmplx( cos( c1 * r1(1) ), sin( -c1 * r1(1) ), 8 )
  eiky0(1)  = cmplx( cos( c2 * r1(2) ), sin( -c2 * r1(2) ), 8 )
  eikz0(1)  = cmplx( cos( c3 * r1(3) ), sin( -c3 * r1(3) ), 8 )
  eikx0(-1) = conjg( eikx0(1) )
  eiky0(-1) = conjg( eiky0(1) )
  eiky0(-1) = conjg( eiky0(1) )

  do p = 2, Kmax1
    eikx0(p)  = eikx0(p-1) * eikx0(1)
    eikx0(-p) = conjg( eikx0(p) )
  end do
  do q = 2, Kmax2
    eiky0(q)  = eiky0(q-1) * eiky0(1)
    eiky0(-q) = conjg(eiky0(q))
  end do
  do r = 2, Kmax3
    eikz0(r)  = eikz0(r-1) * eikz0(1)
  end do

  eikx1(0)  = ( 1,0 )
  eiky1(0)  = ( 1,0 )
  eikz1(0)  = ( 1,0 )
  eikx1(1)  = cmplx( cos( c1 * r2(1) ), sin( -c1 * r2(1) ), 8 )
  eiky1(1)  = cmplx( cos( c2 * r2(2) ), sin( -c2 * r2(2) ), 8 )
  eikz1(1)  = cmplx( cos( c3 * r2(3) ), sin( -c3 * r2(3) ), 8 )
  eikx1(-1) = conjg( eikx1(1) )
  eiky1(-1) = conjg( eiky1(1) )
  eikz1(-1) = conjg( eikz1(1) )

  do p=2, Kmax1
    eikx1(p)  = eikx1(p-1) * eikx1(1)
    eikx1(-p) = conjg( eikx1(p) )
  end do
  do q=2, Kmax2
    eiky1(q)  = eiky1(q-1) * eiky1(1)
    eiky1(-q) = conjg(eiky1(q))
  end do
  do r=2, Kmax3
    eikz1(r)  = eikz1(r-1) * eikz1(1)
  end do

  if (lg) then
    do i=1, K_total
      ord = totk_vectk(i,1:3)
      eikr0 = eikx0(ord(1)) * eiky0(ord(2)) * eikz0(ord(3))
      eikr1 = eikx1(ord(1)) * eiky1(ord(2)) * eikz1(ord(3))
      delta_rhok(i) = eikr1 - eikr0
      delta_cosk(i) = 1 - real( conjg(eikr1) * eikr0 )
    end do
  else
    do i=1, K_total
      ord = totk_vectk(i,1:3)
      eikr0 = eikx0(ord(1)) * eiky0(ord(2)) * eikz0(ord(3))
      eikr1 = eikx1(ord(1)) * eiky1(ord(2)) * eikz1(ord(3))
      delta_rhok(i) = eikr0 - eikr1
      delta_cosk(i) = 1 - real( conjg(eikr1) * eikr0 )
    end do    
  end if

  delta_rhok = delta_rhok * r1(4)

  delta_cosk = delta_cosk * ( r1(4) * r1(4) )

  Del_Recip_erg = sum( exp_ksqr * ( Real( rho_k * delta_rhok ) + delta_cosk ) )

  DeltaE = DeltaE + Del_Recip_erg

end subroutine Delta_Reciprocal_Energy_pH


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


subroutine initialize_lj_parameters
  use global_variables
  implicit none

  rcl2 = rcl*rcl
  nclx = int(Lx/rcl)
  ncly = int(Ly/rcl)
  nclz = int(Lz/rcl)

  clx = Lx/nclx
  cly = Ly/ncly
  clz = Lz/nclz

end subroutine initialize_lj_parameters


subroutine Initialize_ewald_parameters
  use global_variables
  implicit none

  alpha    = ( tau_rf * pi**3 * Nq / (Lx*Ly*Lz)**2 ) ** (1.D0/6)
  alpha2   = alpha * alpha
  rcc  = tol / alpha
  if (rcc<min(Lx/3,Lz/3)) then
    !
    !use verlet list in real space
    Kmax1 = ceiling(tol*Lx*alpha/pi)
    Kmax2 = ceiling(tol*Ly*alpha/pi)
    Kmax3 = ceiling(tol*Lz*alpha/pi)
  else
    rcc = min(Lx/3,Lz/3)
    Kmax1    = ceiling(tol*tol/pi*Lx/rcc)
    Kmax2    = ceiling(tol*tol/pi*Ly/rcc)
    Kmax3    = ceiling(tol*tol/pi*Lz/rcc)
    alpha    = tol / rcc
    alpha2   = alpha * alpha
  end if
  !
  !Cell list parameters
  nclx1 = int(Lx/rcc)     !cell numbers in x direction
  ncly1 = int(Ly/rcc)
  nclz1 = int(Lz/rcc)
  clx1 = Lx/nclx1         !cell length    
  cly1 = Ly/ncly1
  clz1 = Lz/nclz1

end subroutine Initialize_ewald_parameters


subroutine Build_Charge_Ewald
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j

  if (allocated(charge)) deallocate(charge)
  allocate(charge(Nq))
  if (allocated(inv_charge)) deallocate(inv_charge)
  allocate(inv_charge(NN))

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      charge(j) = i
    end if
  end do

  j = 0
  do i = 1, NN
    if ( pos(i,4) /= 0 ) then
      j = j + 1
      inv_charge(i) = j
    else
      inv_charge(i) = 0
    end if
  end do

end subroutine Build_Charge_Ewald


subroutine build_totk_vectk
  !--------------------------------------!
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
  !
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
  !
  !--------------------------------------!
  use global_variables
  implicit none
  complex(kind=8) :: eikx(1:Nq, -Kmax1:Kmax1)
  complex(kind=8) :: eiky(1:Nq, -Kmax2:Kmax2)
  complex(kind=8) :: eikz(1:Nq, 0:Kmax3)
  integer i,j,l,m,n,p,q,r,ord(3)
  real*8 :: c1, c2, c3
  rho_k = 0
  eikx = 0
  eiky = 0
  eikz = 0

  c1 = 2*pi/Lx
  c2 = 2*pi/Ly
  c3 = 2*pi/Lz
  i = cell_list_q(NN+1)
  do while(i /= 0)
    m = inv_charge(i)
    eikx(m,0)  = (1,0)
    eiky(m,0)  = (1,0)
    eikz(m,0)  = (1,0)

    eikx(m,1)  = cmplx( cos(c1*pos(i,1) ), sin( c1*pos(i,1) ), 8 )
    eiky(m,1)  = cmplx( cos(c2*pos(i,2) ), sin( c2*pos(i,2) ), 8 )
    eikz(m,1)  = cmplx( cos(c3*pos(i,3) ), sin( c3*pos(i,3) ), 8 )

    eikx(m,-1) = conjg(eikx(m,1))
    eiky(m,-1) = conjg(eiky(m,1))
    i = cell_list_q(i)
  end do

  do p=2, Kmax1
    i = cell_list_q(NN+1)
    do while(i/=0)
      m = inv_charge(i)
      eikx(m,p)=eikx(m,p-1)*eikx(m,1)
      eikx(m,-p)=conjg(eikx(m,p))
      i = cell_list_q(i)
    end do
  end do
  do q=2, Kmax2
    i = cell_list_q(NN+1)
    do while(i/=0)
      m = inv_charge(i)
      eiky(m,q)=eiky(m,q-1)*eiky(m,1)
      eiky(m,-q)=conjg(eiky(m,q))
      i = cell_list_q(i)
    end do
  end do
  do r=2, Kmax3
    i = cell_list_q(NN+1)
    do while(i/=0)
      m = inv_charge(i)
      eikz(m,r)=eikz(m,r-1)*eikz(m,1)
      i = cell_list_q(i)
    end do
  end do

  do i = 1, K_total
    ord = totk_vectk(i,:)
    m = cell_list_q(NN+1)
    do while(m/=0)
      j = inv_charge(m)
      rho_k(i) = rho_k(i) + &
                 pos(m,4) * eikx(j,ord(1)) * eiky(j,ord(2)) * eikz(j,ord(3))
      m = cell_list_q(m)
    end do
  end do

end subroutine build_rho_k


subroutine Initialize_real_cell_list
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n, p, q, r, x, y, z
  integer :: icelx, icely, icelz

  !
  ! maxium situation, (125,125,100), 6.2Mb RAM is needed.
  if (allocated(hoc_r)) deallocate(hoc_r)
  if (allocated(inv_hoc_r)) deallocate(inv_hoc_r)
  allocate(hoc_r(nclx1,ncly1,nclz1))
  allocate(inv_hoc_r(nclx1,ncly1,nclz1))
  hoc_r = 0
  inv_hoc_r = 0

  if (allocated(cell_list_r)) deallocate(cell_list_r)
  if (allocated(inv_cell_list_r)) deallocate(inv_cell_list_r)
  allocate(cell_list_r(NN))
  allocate(inv_cell_list_r(NN))
  cell_list_r = 0
  inv_cell_list_r = 0

  do i = 1, NN
    if (pos(i,4)/=0) then
      icelx = int((pos(i,1)-1)/clx1) + 1
      icely = int((pos(i,2)-1)/cly1) + 1
      icelz = int((pos(i,3)-1)/clz1) + 1
      cell_list_r(i) = hoc_r(icelx,icely,icelz)
      hoc_r(icelx,icely,icelz) = i
    end if
  end do

  do i = NN, 1, -1
    if (pos(i,4)/=0) then
      icelx = int((pos(i,1)-1)/clx1) + 1
      icely = int((pos(i,2)-1)/cly1) + 1
      icelz = int((pos(i,3)-1)/clz1) + 1
      inv_cell_list_r(i) = inv_hoc_r(icelx,icely,icelz)
      inv_hoc_r(icelx,icely,icelz) = i
    end if
  end do

  !
  ! maxium situation, (125*125*100,28,3), 500Mb RAM is needed.
  if(allocated(cell_near_list_r)) deallocate(cell_near_list_r)
  allocate(cell_near_list_r(nclx1*ncly1*nclz1,27,3))
  cell_near_list_r = 0
  m = 0
  do i = 1, nclx1
    do j = 1, ncly1
      do k = 1, nclz1
        m = m + 1
        n = 0
        do p = -1, 1
          do q = -1, 1
            do r = -1, 1
              x = i + p
              y = j + q
              z = k + r
              n = n + 1
              if (x>nclx1) then
                x = x - nclx1
              elseif (x<=0) then
                x = x + nclx1
              end if
              if (y>ncly1) then
                y = y - ncly1
              elseif (y<=0) then
                y = y + ncly1
              end if
              if (z>nclz1) then
                z = z - nclz1
              elseif (z<=0) then
                z = z + nclz1
              end if
              cell_near_list_r(m,n,1) = x
              cell_near_list_r(m,n,2) = y
              cell_near_list_r(m,n,3) = z
            end do
          end do
        end do
      end do
    end do
  end do

  open(113,file='./data/cell_list_r.txt')
    do i = 1, NN
      write(113,*) i, cell_list_r(i), inv_cell_list_r(i)
    end do
  close(113)

end subroutine Initialize_real_cell_list


subroutine Initialize_cell_list_q
  use global_variables
  implicit none
  integer :: i, j, k 

  if (allocated(cell_list_q)) deallocate(cell_list_q)
  if (allocated(inv_cell_list_q)) deallocate(inv_cell_list_q)
  allocate(cell_list_q(NN+1))
  allocate(inv_cell_list_q(NN+1))
  cell_list_q = 0
  inv_cell_list_q = 0

  do i = 1, NN
    if (pos(i,4) /= 0) then
      cell_list_q(i) = cell_list_q(NN+1)
      cell_list_q(NN+1) = i
    end if
  end do

  do i = NN, 1, -1
    if (pos(i,4) /= 0) then
      inv_cell_list_q(i) = inv_cell_list_q(NN+1)
      inv_cell_list_q(NN+1) = i
    end if
  end do

end subroutine Initialize_cell_list_q


subroutine Initialize_lj_cell_list
  use global_variables
  implicit none
  integer :: i, j, k, l, m, n, p, q, r, x, y, z
  integer :: icelx, icely, icelz

  !
  ! maxium situation, (125,125,100), 6.2Mb RAM is needed.
  if (allocated(hoc_lj)) deallocate(hoc_lj)
  if (allocated(inv_hoc_lj)) deallocate(inv_hoc_lj)
  allocate(hoc_lj(nclx,ncly,nclz))
  allocate(inv_hoc_lj(nclx,ncly,nclz))
  hoc_lj = 0
  inv_hoc_lj = 0

  if (allocated(cell_list_lj)) deallocate(cell_list_lj)
  if (allocated(inv_cell_list_lj)) deallocate(inv_cell_list_lj)
  allocate(cell_list_r(NN-Npc))
  allocate(inv_cell_list_r(NN-Npc))
  cell_list_lj = 0
  inv_cell_list_lj = 0

  do i = 1, NN-Npc
    icelx = int((pos(i,1)-1)/clx) + 1
    icely = int((pos(i,2)-1)/cly) + 1
    icelz = int((pos(i,3)-1)/clz) + 1
    cell_list_lj(i) = hoc_lj(icelx,icely,icelz)
    hoc_lj(icelx,icely,icelz) = i
  end do

  do i = NN-Npc, 1, -1
    icelx = int((pos(i,1)-1)/clx) + 1
    icely = int((pos(i,2)-1)/cly) + 1
    icelz = int((pos(i,3)-1)/clz) + 1
    inv_cell_list_lj(i) = inv_hoc_lj(icelx,icely,icelz)
    inv_hoc_lj(icelx,icely,icelz) = i
  end do

  !
  ! maxium situation, (125*125*100,28,3), 500Mb RAM is needed.
  if(allocated(cell_near_list_lj)) deallocate(cell_near_list_lj)
  allocate(cell_near_list_lj(nclx*ncly*nclz,27,3))
  cell_near_list_lj = 0
  m = 0
  do i = 1, nclx
    do j = 1, ncly
      do k = 1, nclz
        m = m + 1
        n = 0
        do p = -1, 1
          do q = -1, 1
            do r = -1, 1
              x = i + p
              y = j + q
              z = k + r
              n = n + 1
              if (x>nclx) then
                x = x - nclx
              elseif (x<=0) then
                x = x + nclx
              end if
              if (y>ncly) then
                y = y - ncly
              elseif (y<=0) then
                y = y + ncly
              end if
              if (x>nclz) then
                z = z - nclz
              elseif (x<=0) then
                z = z + nclz
              end if
              cell_near_list_lj(m,n,1) = x
              cell_near_list_lj(m,n,2) = y
              cell_near_list_lj(m,n,3) = z
            end do
          end do
        end do
      end do
    end do
  end do

  open(113,file='./data/cell_list_lj.txt')
    do i = 1, NN-Npc
      write(113,*) i, cell_list_lj(i), inv_cell_list_lj(i)
    end do
  close(113)

end subroutine Initialize_lj_cell_list


subroutine delete_cell_list_r(np,rr)
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

end subroutine delete_cell_list_r


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

  cell_list_r(np) = hoc_r(icelx,icely,icelz)
  hoc_r(icelx,icely,icelz) = np

end subroutine add_cell_list_r


subroutine delete_cell_list_lj(np,rr)
  use global_variables
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr
  integer :: icelx,icely,icelz
  integer :: bfi, nti

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

  cell_list_r(np) = hoc_r(icelx,icely,icelz)
  hoc_r(icelx,icely,icelz) = np

end subroutine add_cell_list_lj


subroutine update_rhok
  use global_variables
  implicit none

  rho_k = rho_k + Conjg( delta_rhok )

end subroutine update_rhok


subroutine update_cell_list_ion_move
  use global_variables
  implicit none

  !
  !lj
  call delete_cell_list_lj(ip,pos_ip0)

  call add_cell_list_lj(ip,pos_ip1)

  !
  !real
  call delete_cell_list_r(ip,pos_ip0)

  call add_cell_list_r(ip,pos_ip1)


end subroutine update_cell_list_ion_move


subroutine update_cell_list_pH(lg)
  use global_variables
  implicit none
  logical, intent(in) :: lg

  if (lg) then

    call add_cell_list_r(ip, pos_ip1)

    call add_cell_list_lj(ipi, pos_ipi1)

    call add_cell_list_r(ip, pos_ip1)

    call add_cell_list_r(ipi, pos_ipi1)

  else

    call delete_cell_list_lj(ip, pos_ip1)

    call delete_cell_list_lj(ipi, pos_ipi1)

    call delete_cell_list_r(ip, pos_ip1)

    call delete_cell_list_r(ipi, pos_ipi1)

  end if

end subroutine update_cell_list_pH


subroutine grow_list(new_conf)
  use global_variables
  implicit none
  logical, intent(in) :: new_conf
  integer :: i, np1, np2, n

  if ( .not. new_conf ) then
    do i = 1, num_newconf
      n = Nml + 1 - i
      if (pos_new(n,4)/=0) then
        np1 = pos_new(n,6) + base
        call delete_cell_list_r(np1, pos_new(n,1:4))
      end if
      np2 = pos_new(n,6) + base
      call delete_cell_list_lj(np2, pos_new(n,1:4))
    end do

    do i = 1, num_newconf
      n = Nml + 1 - i
      if (pos_old(n,4)/=0) then
        np1 = pos_old(n,6) + base
        call add_cell_list_r(np1, pos_old(n,1:4))
      end if
      np2 = pos_old(n,6) + base 
      call add_cell_list_lj(np2, pos_old(n,1:4))
    end do
  end if

end subroutine grow_list


subroutine retrace_list(np, rr)
  use global_variables
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr

  if (rr(4)/=0) then

    call delete_cell_list_r(np,rr)

  end if

  call delete_cell_list_lj(np,rr)

end subroutine retrace_list


subroutine regrow_list(np, rr)
  use global_variables
  implicit none
  integer, intent(in) :: np
  real*8, dimension(4), intent(in) :: rr

  if (pos(np,4)/=0) then

    call add_cell_list_r(np,rr)

  end if

  call add_cell_list_lj(np,rr)

end subroutine regrow_list


subroutine write_energy_parameters_Ewald
  !--------------------------------------!
  !
  !--------------------------------------!
  use global_variables
  implicit none

  write(*,*) '******************  Potential  *********************'
  write(*,*) 'Kmax1      :', Kmax1
  write(*,*) 'Kmax2      :', Kmax2
  write(*,*) 'Kmax3      :', Kmax3
  write(*,*) 'K_total    :', K_total
  write(*,*) 'alpha      :', alpha
  write(*,*) 'tol        :', tol
  write(*,*) 'tau_rf     :', tau_rf
  write(*,*) 'rcc        :', rcc
  write(*,*) 'nclx1      :', nclx1
  write(*,*) 'ncly1      :', ncly1
  write(*,*) 'nclz1      :', nclz1
  write(*,*) 'clx1       :', clx1
  write(*,*) 'cly1       :', cly1
  write(*,*) 'clz1       :', clz1
  write(*,*) 'nclx       :', nclx
  write(*,*) 'ncly       :', ncly
  write(*,*) 'nclz       :', nclz
  write(*,*) 'clx        :', clx
  write(*,*) 'cly        :', cly
  write(*,*) 'clz        :', clz
  write(*,*) '****************************************************'

end subroutine write_energy_parameters_Ewald


end module compute_energy


















