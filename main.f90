program main
use global_variables
use input_output
use initialize_update
use compute_energy
implicit none

!##########Data Dictionary############!
  integer :: i, j, k
  real*8  :: EE        ! Total Energy before move
  real*8  :: EE1       ! Total Energy after move
  real*8  :: DeltaE    ! Energy difference
!#####################################!

!#############Initialize##############!
  call cpu_time(started)
  call random_seed()
  !
  !input and initialize system, timing and histogram parameters.
  call initialize_parameters
  !
  !initialize energy and parameters of potential
  call initialize_energy_parameters
  !
  !
  if (restart_or_continue /= 1 ) then
    !
    !initialize position
    call Initialize_position
    !
    !initialize energy arrays
    call Initialize_energy_arrays
    !
    !output data
    call write_pos
    call write_pos1(1)
    !
    !Compute total energy
    call error_analysis(EE)
    i=1
  else
    !
    !read position and histogram data
    call continue_read_data(i)
    !
    !initialize energy arrays
    call initialize_energy_arrays
    !
    !Compute total energy
    call error_analysis(EE)
  end if
!#####################################!

  k = 0
!##############Preheation#############!
  if ( i <= StepNum0 ) then
    do step = i, StepNum0
      call CBMC_Move( EE, DeltaE, k )
      if ( mod(step,DeltaStep1) == 0 ) then
        call compute_physical_quantities
        call total_energy(EE1)
        call write_physical_quantities( step, EE, EE1, DeltaE )
      end if
      if ( mod(step,DeltaStep2) == 0 ) then
        call write_pos1(step)
      end if
    end do
    i = step
  end if
!#####################################!

  call total_energy(EE)
!###############Running###############!
  do step=i, StepNum+StepNum0
    call CBMC_Move( EE, DeltaE, k )
    if ( mod(step,DeltaStep1) == 0 ) then 
      call compute_physical_quantities
      call total_energy(EE1)
      call compute_radial_distribution_function
      call write_physical_quantities( step, EE, EE1, DeltaE )
    end if
    if ( mod(step, DeltaStep2) == 0 ) then
      call write_pos1(step)
    end if
  end do
!#####################################!

!###############Finished##############!
  call cpu_time(finished)
  total_time=finished-started+total_time
  call write_pos1(step)
  write(*,*) 'finished!'
!#####################################!

end program main








