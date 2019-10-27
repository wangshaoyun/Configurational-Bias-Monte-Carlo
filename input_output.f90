module input_output
  implicit none
  
  save
  real*8, allocatable, dimension(:,:), private :: rdf
  real*8, allocatable, dimension(:,:), private :: gr_p
  real*8, allocatable, dimension(:,:), private :: gr_p_a
  real*8, allocatable, dimension(:,:), private :: gr_s
  real*8, allocatable, dimension(:,:), private :: gr_s_a
  real*8, allocatable, dimension(:,:), private :: gr_c_a


  contains

subroutine initialize_parameters
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none
  logical alive 
  !
  !Input parameters
  call read_data
  !
  ! data operation
  call data_operation
  !
  !Write data
  call write_data_to_screen
  !
  !Allocate arrays and initialize them
  call data_allocate

end subroutine initialize_parameters

subroutine data_operation
  use global_variables
  implicit none
  integer :: charge_ions

  Npe = Nml * Ngl
  if ( abs(qq) == 0 ) then
    Nq_PE = 0
  else
    if (man/=0) then
      Nq_PE = Nml/man*Nga
    else
      Nq_PE = 0
    end if
  end if
  Nq_pe_net = Nq_PE
  !
  ! number of salt particles
  if ( abs(qqi) ==0 ) then
    Nq_salt_ions = 0
  else
    charge_ions = nint( ion_ratio * Nq_PE * abs(qq) )
    Nq_salt_ions = nint( charge_ions / abs(qqi) )
  end if
  !
  ! charges on cylinder
  Nqc = nint( Lz * rho_c )
  qqc = Npc / Nqc * (-qq)/abs(qq)
  !
  ! total charges and particles
  Nq = Nq_PE * ( abs(qq)+1 ) + Nq_salt_ions * ( abs(qqi) + 1 ) + Npc + Nqc
  NN = Npe + Nq - Nq_PE
  !
  !
  Nq_ions_net = NN - Npc - Npe
  Nq_net = Nq_pe_net + Nq_ions_net
  !
  !System size
  Lz = Nml*1.2*R_bond
  Lx = sqrt(Ngl/rho/Lz)
  Ly = Lx
  ratio_xz = Lx / Lz
  !
  !Judge whether restart or continue
  if (restart_or_continue==1) then
    Inquire(file='start_time.txt',exist=alive)
    if (alive) then
      open(11,file='./start_time.txt')
        read(11,*) restart_or_continue
      close(11)
    else
      restart_or_continue=0
    end if
  end if

end subroutine data_operation


subroutine read_data
  use global_variables
  implicit none

  open(unit=100, file='system_data.txt')
    read(100,*) restart_or_continue
    read(100,*) multistep_or_not
    read(100,*) rho
    read(100,*) Beta
    read(100,*) Nml
    read(100,*) Ngl
    read(100,*) man
    read(100,*) qq
    read(100,*) qqi
    read(100,*) rho_c
    read(100,*) ion_ratio
    read(100,*) R_bond
    read(100,*) r_cy
    read(100,*) k_try
    read(100,*) Npc
    read(100,*) StepNum0
    read(100,*) StepNum
    read(100,*) DeltaStep
    read(100,*) multistep
    read(100,*) DeltaStep1
    read(100,*) DeltaStep2
    read(100,*) DeltaStep3
    read(100,*) dr
    read(100,*) best_accept_ratio
    read(100,*) delta_dr
    read(100,*) pH_pKa
  close(100)

end subroutine read_data


subroutine write_data_to_screen
  use global_variables
  implicit none

  write(*,*)
  write(*,*)
  write(*,*) '******************system_data***********************'
  write(*,*) 'Total chains,                          Ngl:', Ngl
  write(*,*) 'Particles of each chain,               Nml:', Nml
  write(*,*) 'Total particles,                        NN:', NN
  write(*,*) 'total charged particles,                Nq:', Nq
  write(*,*) 'total charged particles in polymer   Nq_PE:', Nq_pe
  write(*,*) 'total charged salt particles  Nq_salt_ions:', Nq_salt_ions
  write(*,*) 'total charged particles on cylind      Npc:', Npc
  write(*,*) 'Bond length of polymer,             R_bond:', R_bond
  write(*,*) 'Length of the box,                      Lx:', Lx
  write(*,*) 'Width of the box,                       Ly:', Ly
  write(*,*) 'Height of the box,                      Lz:', Lz
  write(*,*) 'radius of cylinder,                   r_cy:', r_cy
  write(*,*) 'line charge density on cylinder,     rho_c:', rho_c
  write(*,*) 'Beta=1/kT,                            Beta:', Beta
  write(*,*) 'Charges of polymer,                     qq:', qq
  write(*,*) 'Charges of salt ions,                  qqi:', qqi
  write(*,*) 'Each man_s monomers with one charge, man_s:', man
  write(*,*) 'pH-pKa,                             pH-pKa:', pH_pKa
  write(*,*) '****************************************************'

  write(*,*)
  write(*,*) '******************running_steps*********************'
  write(*,*) 'restart (0), continue (0), restart_continue:',restart_or_continue
  write(*,*) 'multistep                                  :',multistep_or_not
  write(*,*) 'Preheating steps                           :',StepNum0
  write(*,*) 'Running steps                              :',StepNum
  write(*,*) 'Total steps                                :',(StepNum0+StepNum)
  write(*,*) 'DeltaStep                                  :',DeltaStep
  write(*,*) 'DeltaStep1                                 :',DeltaStep1
  write(*,*) 'DeltaStep2                                 :',DeltaStep2
  write(*,*) 'DeltaStep3                                 :',DeltaStep3
  write(*,*) 'multistep                                  :',multistep
  write(*,*) 'dr                                         :',dr
  write(*,*) 'best_accept_ratio                          :',best_accept_ratio
  write(*,*) '****************************************************'
  write(*,*)
  write(*,*)

end subroutine write_data_to_screen


subroutine data_allocate
  use global_variables
  implicit none

  ! pos(i,6): position i now is occupied by particle pos(i,6)
  ! pos(i,7): particle i now is in position pos(i,7)
  allocate( pos(NN, 7)      )
  allocate( pos_new(Nml,7)  )
  allocate( pos_old(Nml,7)  )
  allocate( rdf(SizeHist,2) )
  allocate( gr_p(SizeHist,2) )
  allocate( gr_p_a(SizeHist,2) )
  allocate( gr_s(SizeHist,2) )
  allocate( gr_s_a(SizeHist,2) )
  allocate( gr_c_a(SizeHist,2) )

  rdf = 0
  gr_p = 0
  gr_p_a = 0
  gr_s = 0
  gr_s_a = 0
  gr_c_a = 0

end subroutine data_allocate


subroutine continue_read_data(l)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none
  integer, intent(out) :: l
  integer :: i, j 

  open(20,file='./data/pos1.txt')
    read(20,*) ((pos(i,j),j=1,7),i=1,NN)
  close(20)
  open(19,file='./start_time.txt')
    read(19,*)
    read(19,*) l
    read(19,*) dr
    read(19,*) total_time
  close(19)

  open(21, file = './data/rdf.txt')
  open(22, file = './data/gr.txt' )
    read(21,*) ((rdf(i,j),j=1,2),i=1,SizeHist)
    read(22,*) ((gr(i,j),j=1,2), i=1,SizeHist)
  close(22)
  close(21)

  Nq_pe_net = 0
  do i = 1, Npe
    if (pos(i,4)/=0) then
      Nq_pe_net = Nq_pe_net + 1
    end if
  end do

  Nq_ions_net
  do i = Npe+1, NN-Npc
    if (pos(i,4)/=0) then
      Nq_ions_net = Nq_ions_net + 1
    end if
  end do
  
  Nq_net = Nq_pe_net + Nq_ions_net

end subroutine continue_read_data


subroutine compute_physical_quantities
  !----------------------------------------!
  ! compute pressure
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
  use compute_energy
  implicit none
  integer i, j, k, m, n
  real*8 :: rr, Rg, Rgi
  real*8, dimension(3) :: rij
  
  Rg = 0
  do i = 1, Ngl
    Rgi = 0
    do j = 1, Nml-1
      do k = j+1, Nml
        m = (i-1)*Nml + j
        n = (i-1)*Nml + k
          call rij_and_rr(rij, rr, m, n)
          Rgi  = Rgi + rr
      end do
    end do
    Rgi = Rgi / Nml / (Nml-1)
    Rg = Rg + Rgi
  end do
  Rg = Rg / Ngl
  !
  !Output Pressure
  open(37, position='append', file='./data/Rg.txt')
    write(37,370) 1.*step, Rg
    370 format(3F15.6)
  close(37)
  
end subroutine compute_physical_quantities


subroutine compute_radial_distribution_function
  use global_variables
  implicit none
  integer :: i,j,k
  real*8  :: rr, del_r
  real*8, dimension(3) :: rij

  del_r = Lx/2/SizeHist

  do i = 1, Npe-1
    do j = i+1, Npe
      call rij_and_rr(rij,rr,i,j)
      if (sqrt(rr)<Lx/2) then
        k = int( sqrt(rr)/del_r ) + 1
        rdf(k,2) = rdf(k,2) + 2 
      end if
    end do 
  end do 

  do i = 1, Npe
    rr = sqrt(pos(i,1)*pos(i,1) + pos(i,2)*pos(i,2))
    if (rr<Lx/2) then
      k = int( sqrt(rr)/del_r ) + 1
      gr_p(k,2) = gr_p(k,2) + 1
    end if
  end do

  do i = Npe+1, Npe+Nq_PE
    rr = sqrt(pos(i,1)*pos(i,1) + pos(i,2)*pos(i,2))
    if (rr<Lx/2) then
      k = int( sqrt(rr)/del_r ) + 1
      gr_p_a(k,2) = gr_p_a(k,2) + 1
    end if
  end do

  do i = Npe+Nq_PE+1, Npe+Nq_PE+Nq_salt_ions
    rr = sqrt(pos(i,1)*pos(i,1) + pos(i,2)*pos(i,2))
    if (rr<Lx/2) then
      k = int( sqrt(rr)/del_r ) + 1
      gr_s(k,2) = gr_s(k,2) + 1
    end if
  end do

  do i = Npe+Nq_PE+Nq_salt_ions+1, Npe+Nq_PE+Nq_salt_ions*(1+abs(qqi))
    rr = sqrt(pos(i,1)*pos(i,1) + pos(i,2)*pos(i,2))
    if (rr<Lx/2) then
      k = int( sqrt(rr)/del_r ) + 1
      gr_s_a(k,2) = gr_s_a(k,2) + 1
    end if
  end do

  do i = NN-Npc-Nqc+1, NN-Npc
    rr = sqrt(pos(i,1)*pos(i,1) + pos(i,2)*pos(i,2))
    if (rr<Lx/2) then
      k = int( sqrt(rr)/del_r ) + 1
      gr_c_a(k,2) = gr_c_a(k,2) + 1
    end if
  end do  

  open(31,file='./data/rdf.txt')
    do i=1, SizeHist
      write(31,310) del_r*i, rdf(i,2), gr_p(i,2), gr_p_a(i,2), gr_s(i,2), &
                    gr_s_a(i,2), gr_c_a(i,2)
      310 format(7F20.6)
    end do
  close(31)

end subroutine compute_radial_distribution_function


subroutine write_pos
  !----------------------------------------!
  !write position to pos.txt
  !input:
  !  pos
  !External Variants:
  !  NN
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: i

  open(30,file='./data/pos.txt')
    do i=1, NN
      write(30,300) pos(i,1), pos(i,2), pos(i,3), pos(i,4), &
                    pos(i,5), pos(i,6), pos(i,7)
      300 format(7F15.6)
    end do
  close(30)

end subroutine write_pos


subroutine write_pos1(j)
  !----------------------------------------!
  !write position to pos1.txt
  !input:
  !  pos
  !External Variants:
  !  NN
  !----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: j
  integer :: i

  open(30,file='./data/pos1.txt')
    do i=1, NN
      write(30,300) pos(i,1), pos(i,2), pos(i,3), pos(i,4), &
                    pos(i,5), pos(i,6), pos(i,7)
      300 format(7F15.6)
    end do
  close(30)

  open(109,file='./start_time.txt')
    write(109,*) 1
    write(109,*) j
    write(109,*) dr
    call cpu_time(finished)
    total_time=total_time+finished-started
    call cpu_time(started)
    write(109,*) total_time
    write(109,*) 'time:(minutes):', real(total_time/60)
    write(109,*) 'time:(hours)  :', real(total_time/3600)
    write(109,*) 'time:(days)   :', real(total_time/86400)
    write(109,*) 'accept_ratio   :', accept_ratio
    write(109,*) 'Lx            :', Lx
    write(109,*) 'Ly            :', Ly
    write(109,*) 'Lz            :', Lz
    write(109,*) 'Nq_net        :', Nq_net
    write(109,*) 'Nq_pe_net     :', Nq_pe_net
    write(109,*) 'Nq_ions_net   :', Nq_ions_net
    write(109,*) 'Ngl           :', Ngl
    write(109,*) 'Nml           :', Nml
    write(109,*) 'Nq            :', Nq
    write(109,*) 'NN            :', NN
    write(109,*) 'Nq_PE         :', Nq_PE
    write(109,*) 'Npe           :', Npe
    write(109,*) 'Nq_salt_ions  :', Nq_salt_ions
    write(109,*) 'rho           :', rho
    write(109,*) 'rho_c         :', rho_c
    write(109,*) 'Beta          :', Beta
    write(109,*) 'qq            :', qq
    write(109,*) 'qqi           :', qqi
    write(109,*) 'man           :', man
    write(109,*) 'pH-pKa        :', pH_pKa
    write(109,*) 'restart_continue:',restart_or_continue
    write(109,*) 'StepNum0      :', StepNum0
    write(109,*) 'StepNum       :', StepNum
    write(109,*) 'StepNum0+StepNum:', (StepNum0+StepNum)
    write(109,*) 'DeltaStep     :', DeltaStep
    write(109,*) 'DeltaStep1    :', DeltaStep1
    write(109,*) 'DeltaStep2    :', DeltaStep2
    write(109,*) 'DeltaStep3    :', DeltaStep3
    write(109,*) 'multistep     :', multistep
  close(109)

end subroutine write_pos1


subroutine write_physical_quantities(j, EE, EE1, DeltaE)
  !----------------------------------------!
  !
  !----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: j
  real*8,  intent(in) :: EE
  real*8,  intent(in) :: EE1
  real*8,  intent(in) :: DeltaE

  open(37,position='append', file='./data/energy_and_time.txt')
    write(37,370) 1.*j, EE, EE1, DeltaE, accept_ratio_p, accept_ratio_ion, &
                  accept_ratio_pH
    370 format(7F15.6)
  close(37)

end subroutine write_physical_quantities



end module input_output
