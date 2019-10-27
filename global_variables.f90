module global_variables

  implicit none
  save

!########################constants#########################!
  real*8, parameter :: pi=3.141592653589793D0     !Circumference ratio pi
  real*8, parameter :: gamma=.5772156649015329D0  !Euler constants gamma
!########################constants#########################!

!####################systems coefficient###################!
  integer :: Ngl      !Number of linear chains
  integer :: Nml      !Number of monomers in each chain
  integer :: NN       !Total particles in the system
  integer :: Nq       !Total charge in the system
  integer :: Nqc      !charges on cylinder
  integer :: Npc      !particles on cylinder
  integer :: Npe      !Total monomers in Polyelectrolytes(PE)
  integer :: Nq_PE    !charges on PE
  integer :: Nq_pe_net    ! net charges on PE
  integer :: Nq_salt_ions ! salt ions
  integer :: Nq_ions_net  ! net ions
  integer :: Nq_net   ! net charges in system 
  integer :: man      !Manning effect, each man particle have one charge
  real*8  :: r_cy     !radius of cylinder
  real*8  :: rho      !Polymer density
  real*8  :: rho_c    !line charge density on cylinder 
  real*8  :: qq       !Charge of charged monomers
  real*8  :: qqi      !charge of salt
  real*8  :: qqc      !charge of particles on cylinder
  real*8  :: ion_ratio!ratio of salt
  real*8  :: Lx       !Length of cell in x direction
  real*8  :: Ly       !Length of cell in y direction
  real*8  :: Lz       !Length of cell in z direction
  real*8  :: ratio_xz !Rotio of length x and height z of the box
  real*8  :: R_bond   !Initial bond length of polymers
  real*8  :: Beta     !Beta=1/(kB*T), T is temperature, 
                      !kB is Boltzmann constant
  real*8  :: pH_pKa   !pH - pKa
!##################end systems coefficient#################!


!##################running and Histogram###################!
  integer :: restart_or_continue  !restart or continue after breaking off 
  integer :: multistep_or_not     !0 : not use multistep, 1: use multistep
  integer :: StepNum0             !steps of preheating
  integer :: StepNum              !steps of running
  integer :: DeltaStep            !steps of pH titeration
  integer :: DeltaStep1           !step inteval, physical quantities
  integer :: DeltaStep2           !step histogram
  integer :: DeltaStep3           !step write data
  integer :: step                 !steps of calculate the physical quantities
  integer :: multistep            !multistep
  real*8  :: dr                   !length of each moving
  real*8  :: total_num = 0        !Total choose number
  real*8  :: accept_num = 0       !accepted number
  real*8  :: accept_ratio_p       !accepted ratio of polymer
  real*8  :: accept_ratio_ion     !accepted ratio of ions
  real*8  :: accept_ratio_pH      !accepted ratio of pH
  real*8  :: best_accept_ratio    !best accepted ratio
  real*8  :: delta_dr             !adjust move distance
  !
  !timing
  real*8  :: started    = 0       !time at starting
  real*8  :: finished   = 0       !time at finishing
  real*8  :: total_time = 0       !total time of the simulation
  !
  !histogram
  integer :: SizeHist             !number of histogram which is equally divided
!################end running and Histogram#################!

!##########################arrays##########################!
 !pos(i,6): position i now is occupied by particle pos(i,6)
 !pos(i,7): particle i now is in position pos(i,7)
 real*8, allocatable, dimension(:,:) :: pos     !old position array
 real*8, allocatable, dimension(:,:) :: pos_old !new position of part of chains
 real*8, allocatable, dimension(:,:) :: pos_new !new position of part of chains
 integer :: ic_newconf                          !The chain that is choosed
 integer :: num_newconf                         !numbers of monomers regrowed
 integer :: base                                !(ic_newconf)*Nml
 integer :: k_try                               !try numbers
 integer :: ip
 integer :: ipi
 real*8  :: pos_ip0(4)
 real*8  :: pos_ip1(4)
 real*8  :: pos_ipi0(4)
 real*8  :: pos_ipi1(4)
!########################end arrays########################!

contains 

subroutine periodic_condition(rr)
  !--------------------------------------!
  !3D Peridodic condition of position vector
  !   
  !Input
  !   rr
  !Output
  !   rr
  !External Variables
  !   Lx, Ly, Lz
  !Routine Referenced:
  !1.
  !--------------------------------------!
  implicit none
  real*8, intent(inout) :: rr(3)

  if ( rr(1) > Lx/2 ) then
    rr(1) = rr(1) - Lx
  elseif( rr(1) <= -Lx/2 ) then
    rr(1) = rr(1) + Lx
  end if
  if ( rr(2) > Ly/2 ) then
    rr(2) = rr(2) - Ly
  elseif( rr(2) <= -Ly/2 ) then
    rr(2) = rr(2) + Ly
  end if
  if ( rr(3) > Lz/2 ) then
    rr(3) = rr(3) - Lz
  elseif( rr(3) <= -Lz/2 ) then
    rr(3) = rr(3) + Lz
  end if

end subroutine periodic_condition


subroutine rij_and_rr(rij, rsqr, i, j)
  !-----------------------------------------!
  !compute displacement vector and displacement of two particles
  !input:
  !  i, j(particle number) 
  !output:
  !  rij(displacement vecter), rr(square of displacement)
  !External Variant:
  !  Lx,Ly,Lz(used in period condition)
  !  pos
  !note:
  !  including period condition
  !-----------------------------------------!
  implicit none
  real*8, dimension(3), intent(out) :: rij
  real*8, intent(out) :: rsqr
  integer, intent(in) :: i
  integer, intent(in) :: j

  rij = pos(i,1:3) - pos(j,1:3)

  ! Periodic Condition
  if ( rij(1) > Lx/2 ) then
    rij(1) = rij(1) - Lx
  elseif( rij(1) <= -Lx/2 ) then
    rij(1) = rij(1) + Lx
  end if
  if ( rij(2) > Ly/2 ) then
    rij(2) = rij(2) - Ly
  elseif( rij(2) <= -Ly/2 ) then
    rij(2) = rij(2) + Ly
  end if
  if ( rij(3) > Lz/2 ) then
    rij(3) = rij(3) - Lz
  elseif( rij(3) <= -Lz/2 ) then
    rij(3) = rij(3) + Lz
  end if

  rsqr = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

end subroutine rij_and_rr


subroutine gauss(sigma, mu, x)
  !--------------------------------------!
  !Gaussian distribution function
  !   
  !Input
  !   mu, sigma
  !Output
  !   x
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !--------------------------------------!
  implicit none
  real*8, intent(in)  :: simga
  real*8, intent(in)  :: mu
  real*8, intent(out) :: x
  real*8 :: r, v1, v2, rnd1, rnd2

  r = 2
  do while ( r > 1 ) 
    call random_number(rnd1)
    call random_number(rnd2)
    v1 = 2 * rnd1 - 1
    v2 = 2 * rnd2 - 1
    r = v1 * v1 + v2 * v2
  end do
  x = v1 * sqrt( -2*log(r)/r )
  x = mu + simga * x

end subroutine gauss


subroutine cross_product(x, y, z)
  implicit none
  real*8, dimension(3), intent(in)  :: x
  real*8, dimension(3), intent(in)  :: y
  real*8, dimension(3), intent(out) :: z

  z(1) = x(2)*y(3) - x(3)*y(2)
  z(2) = x(3)*y(1) - x(1)*y(3)
  z(3) = x(1)*y(2) - x(2)*y(1)

end subroutine cross_product

end module global_variables

