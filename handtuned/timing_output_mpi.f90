module timing_output_mpi_mod
  use, intrinsic :: iso_fortran_env, only : real64
  use mpi_f08
  use timing_mod, only : model_timers, region_total, region_initialization, region_total_excluding_init, &
    region_propagation, region_forcing, region_source_terms, region_statistics
  implicit none

  private
  public :: write_timers_json_mpi

contains

  subroutine write_timers_json_mpi(timers, implementation, comm)
    type(model_timers), intent(in) :: timers
    character(len=*), intent(in) :: implementation
    type(MPI_Comm), intent(in) :: comm
    real(real64), allocatable :: all_wall(:, :)
    real(real64), allocatable :: all_cpu(:, :)
    real(real64), allocatable :: step_wall(:)
    real(real64), allocatable :: step_cpu(:)
    real(real64), allocatable :: all_step_wall(:, :)
    real(real64), allocatable :: all_step_cpu(:, :)
    real(real64), allocatable :: step_region_wall(:, :)
    real(real64), allocatable :: step_region_cpu(:, :)
    real(real64), allocatable :: all_step_region_wall(:, :, :)
    real(real64), allocatable :: all_step_region_cpu(:, :, :)
    integer, allocatable :: all_calls(:, :)
    integer, allocatable :: step_region_calls(:, :)
    integer, allocatable :: all_step_region_calls(:, :, :)
    integer, allocatable :: step_number(:)
    integer, allocatable :: all_step_number(:, :)
    integer, allocatable :: all_step_counts(:)
    integer :: rank
    integer :: comm_size
    integer :: max_steps
    integer :: step
    integer :: ierr
    integer :: region_count

    region_count = size(timers%wall_seconds)
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, comm_size, ierr)
    max_steps = size(timers%steps)
    allocate(all_wall(region_count, comm_size))
    allocate(all_cpu(region_count, comm_size))
    allocate(all_calls(region_count, comm_size))
    allocate(all_step_counts(comm_size))
    allocate(step_wall(max_steps))
    allocate(step_cpu(max_steps))
    allocate(step_number(max_steps))
    allocate(all_step_wall(max_steps, comm_size))
    allocate(all_step_cpu(max_steps, comm_size))
    allocate(all_step_number(max_steps, comm_size))
    allocate(step_region_wall(region_count, max_steps))
    allocate(step_region_cpu(region_count, max_steps))
    allocate(step_region_calls(region_count, max_steps))
    allocate(all_step_region_wall(region_count, max_steps, comm_size))
    allocate(all_step_region_cpu(region_count, max_steps, comm_size))
    allocate(all_step_region_calls(region_count, max_steps, comm_size))

    step_wall = 0.0_real64
    step_cpu = 0.0_real64
    step_number = 0
    step_region_wall = 0.0_real64
    step_region_cpu = 0.0_real64
    step_region_calls = 0
    do step = 1, timers%step_count
      step_wall(step) = timers%steps(step)%wall_seconds
      step_cpu(step) = timers%steps(step)%cpu_seconds
      step_number(step) = timers%steps(step)%step
      step_region_wall(:, step) = timers%steps(step)%region_wall_seconds
      step_region_cpu(:, step) = timers%steps(step)%region_cpu_seconds
      step_region_calls(:, step) = timers%steps(step)%region_calls
    end do

    call MPI_Gather(timers%wall_seconds, region_count, MPI_DOUBLE_PRECISION, all_wall, region_count, &
      MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Gather(timers%cpu_seconds, region_count, MPI_DOUBLE_PRECISION, all_cpu, region_count, &
      MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Gather(timers%calls, region_count, MPI_INTEGER, all_calls, region_count, MPI_INTEGER, 0, comm, ierr)
    call MPI_Gather(timers%step_count, 1, MPI_INTEGER, all_step_counts, 1, MPI_INTEGER, 0, comm, ierr)
    call MPI_Gather(step_wall, max_steps, MPI_DOUBLE_PRECISION, all_step_wall, max_steps, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Gather(step_cpu, max_steps, MPI_DOUBLE_PRECISION, all_step_cpu, max_steps, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Gather(step_number, max_steps, MPI_INTEGER, all_step_number, max_steps, MPI_INTEGER, 0, comm, ierr)
    call MPI_Gather(step_region_wall, region_count * max_steps, MPI_DOUBLE_PRECISION, all_step_region_wall, &
      region_count * max_steps, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Gather(step_region_cpu, region_count * max_steps, MPI_DOUBLE_PRECISION, all_step_region_cpu, &
      region_count * max_steps, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Gather(step_region_calls, region_count * max_steps, MPI_INTEGER, all_step_region_calls, &
      region_count * max_steps, MPI_INTEGER, 0, comm, ierr)

    if (rank == 0) then
      call write_gathered_timers_json(timers, implementation, all_wall, all_cpu, all_calls, all_step_counts, &
        all_step_wall, all_step_cpu, all_step_number, all_step_region_wall, all_step_region_cpu, all_step_region_calls)
    end if

  end subroutine write_timers_json_mpi

  subroutine write_gathered_timers_json(root_timers, implementation, all_wall, all_cpu, all_calls, all_step_counts, &
      all_step_wall, all_step_cpu, all_step_number, all_step_region_wall, all_step_region_cpu, all_step_region_calls)
    type(model_timers), intent(in) :: root_timers
    character(len=*), intent(in) :: implementation
    real(real64), intent(in) :: all_wall(:, :)
    real(real64), intent(in) :: all_cpu(:, :)
    integer, intent(in) :: all_calls(:, :)
    integer, intent(in) :: all_step_counts(:)
    real(real64), intent(in) :: all_step_wall(:, :)
    real(real64), intent(in) :: all_step_cpu(:, :)
    integer, intent(in) :: all_step_number(:, :)
    real(real64), intent(in) :: all_step_region_wall(:, :, :)
    real(real64), intent(in) :: all_step_region_cpu(:, :, :)
    integer, intent(in) :: all_step_region_calls(:, :, :)
    integer :: unit
    integer :: ios
    integer :: rank_index
    integer :: region
    integer :: step
    integer :: rank_count
    integer :: region_count

    rank_count = size(all_step_counts)
    region_count = size(all_wall, 1)
    open(newunit=unit, file='timers.json', status='replace', action='write', iostat=ios)
    if (ios /= 0) return

    write(unit, '(a)') '{'
    write(unit, '(a,a,a)') '  "implementation": "', trim(implementation), '",'
    write(unit, '(a)') '  "ranks": ['
    do rank_index = 1, rank_count
      write(unit, '(a)') '    {'
      write(unit, '(a,i0,a)') '      "rank": ', rank_index - 1, ','
      write(unit, '(a)') '      "regions": {'
      do region = 1, region_count
        call write_region(unit, region, all_wall(region, rank_index), all_cpu(region, rank_index), &
          all_calls(region, rank_index), region < region_count)
      end do
      write(unit, '(a)') '      },'
      write(unit, '(a)') '      "steps": ['
      do step = 1, all_step_counts(rank_index)
        call write_gathered_step(unit, root_timers, step, all_step_number(step, rank_index), &
          all_step_wall(step, rank_index), all_step_cpu(step, rank_index), &
          all_step_region_wall(:, step, rank_index), all_step_region_cpu(:, step, rank_index), &
          all_step_region_calls(:, step, rank_index), step < all_step_counts(rank_index))
      end do
      write(unit, '(a)') '      ]'
      if (rank_index < rank_count) then
        write(unit, '(a)') '    },'
      else
        write(unit, '(a)') '    }'
      end if
    end do
    write(unit, '(a)') '  ],'
    write(unit, '(a)') '  "schema": "ecwam-timers-v2",'
    write(unit, '(a)') '  "summary": {'
    write(unit, '(a)') '    "regions": {'
    do region = 1, region_count
      call write_summary_region_from_ranks(unit, region, all_wall(region, :), all_cpu(region, :), all_calls(region, :), &
        region < region_count)
    end do
    write(unit, '(a)') '    }'
    write(unit, '(a)') '  }'
    write(unit, '(a)') '}'
    close(unit)
  end subroutine write_gathered_timers_json

  subroutine write_gathered_step(unit, root_timers, step_index, step_number, wall_seconds, cpu_seconds, &
      region_wall_seconds, region_cpu_seconds, region_calls, has_comma)
    integer, intent(in) :: unit
    type(model_timers), intent(in) :: root_timers
    integer, intent(in) :: step_index
    integer, intent(in) :: step_number
    real(real64), intent(in) :: wall_seconds
    real(real64), intent(in) :: cpu_seconds
    real(real64), intent(in) :: region_wall_seconds(:)
    real(real64), intent(in) :: region_cpu_seconds(:)
    integer, intent(in) :: region_calls(:)
    logical, intent(in) :: has_comma
    character(len=14) :: model_time
    character(len=14) :: forcing_time
    character(len=1) :: comma
    integer :: region

    model_time = ''
    forcing_time = ''
    if (step_index >= 1 .and. step_index <= root_timers%step_count) then
      model_time = root_timers%steps(step_index)%model_time
      forcing_time = root_timers%steps(step_index)%forcing_time
    end if
    comma = merge(',', ' ', has_comma)
    write(unit, '(a,es24.16,a,a,a,a,a,i0,a,es24.16,a)') '        {"cpu_seconds": ', cpu_seconds, &
      ', "forcing_time": "', trim(forcing_time), '", "model_time": "', trim(model_time), &
      '", "step": ', step_number, ', "wall_seconds": ', wall_seconds, ','
    write(unit, '(a)') '         "regions": {'
    do region = region_total_excluding_init, region_statistics
      call write_step_region(unit, region, region_wall_seconds(region), region_cpu_seconds(region), &
        region_calls(region), region < region_statistics)
    end do
    write(unit, '(a)') '         }'
    write(unit, '(a,a)') '        }', comma
  end subroutine write_gathered_step

  subroutine write_step_region(unit, region, wall_seconds, cpu_seconds, calls, has_comma)
    integer, intent(in) :: unit
    integer, intent(in) :: region
    real(real64), intent(in) :: wall_seconds
    real(real64), intent(in) :: cpu_seconds
    integer, intent(in) :: calls
    logical, intent(in) :: has_comma
    character(len=1) :: comma

    comma = merge(',', ' ', has_comma)
    write(unit, '(a,a,a,es24.16,a,es24.16,a,i0,a,a)') '           "', trim(region_name(region)), &
      '": {"cpu_seconds": ', cpu_seconds, ', "wall_seconds": ', wall_seconds, ', "calls": ', calls, '}', comma
  end subroutine write_step_region

  subroutine write_summary_region_from_ranks(unit, region, wall_seconds, cpu_seconds, calls, has_comma)
    integer, intent(in) :: unit
    integer, intent(in) :: region
    real(real64), intent(in) :: wall_seconds(:)
    real(real64), intent(in) :: cpu_seconds(:)
    integer, intent(in) :: calls(:)
    logical, intent(in) :: has_comma
    character(len=1) :: comma
    real(real64) :: rank_count

    rank_count = real(size(wall_seconds), real64)
    comma = merge(',', ' ', has_comma)
    write(unit, '(a,a,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,i0,a,a)') &
      '        "', trim(region_name(region)), '": {"cpu_seconds_max": ', maxval(cpu_seconds), &
      ', "cpu_seconds_mean": ', sum(cpu_seconds) / rank_count, ', "cpu_seconds_sum": ', sum(cpu_seconds), &
      ', "wall_seconds_max": ', maxval(wall_seconds), ', "wall_seconds_mean": ', sum(wall_seconds) / rank_count, &
      ', "wall_seconds_sum": ', sum(wall_seconds), ', "calls_sum": ', sum(calls), '}', comma
  end subroutine write_summary_region_from_ranks

  subroutine write_region(unit, region, wall_seconds, cpu_seconds, calls, has_comma)
    integer, intent(in) :: unit
    integer, intent(in) :: region
    real(real64), intent(in) :: wall_seconds
    real(real64), intent(in) :: cpu_seconds
    integer, intent(in) :: calls
    logical, intent(in) :: has_comma
    character(len=1) :: comma

    comma = merge(',', ' ', has_comma)
    write(unit, '(a,a,a,es24.16,a,es24.16,a,i0,a,a)') '        "', trim(region_name(region)), &
      '": {"cpu_seconds": ', cpu_seconds, ', "wall_seconds": ', wall_seconds, ', "calls": ', calls, '}', comma
  end subroutine write_region

  function region_name(region) result(name)
    integer, intent(in) :: region
    character(len=32) :: name

    select case (region)
    case (region_total)
      name = 'total'
    case (region_initialization)
      name = 'initialization'
    case (region_total_excluding_init)
      name = 'total_excluding_init'
    case (region_propagation)
      name = 'propagation'
    case (region_forcing)
      name = 'forcing'
    case (region_source_terms)
      name = 'source_terms'
    case (region_statistics)
      name = 'statistics'
    case default
      name = 'unknown'
    end select
  end function region_name

end module timing_output_mpi_mod
