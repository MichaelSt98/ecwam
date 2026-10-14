! Lightweight model timing with no synchronization or I/O in hot paths.
!
! Region totals are local-rank measurements; timing_output_mpi aggregates them.
! Step records retain model and forcing times for post-initialization analysis.
module timing_mod
  use, intrinsic :: iso_fortran_env, only : int64, real64
  implicit none

  private
  public :: model_timers, timer_mark
  public :: original_timers
  public :: init_timers, start_timer, add_elapsed, record_step
  public :: region_total, region_initialization, region_propagation
  public :: region_forcing, region_source_terms, region_statistics
  public :: region_total_excluding_init

  integer, parameter :: region_total = 1
  integer, parameter :: region_initialization = 2
  integer, parameter :: region_total_excluding_init = 3
  integer, parameter :: region_propagation = 4
  integer, parameter :: region_forcing = 5
  integer, parameter :: region_source_terms = 6
  integer, parameter :: region_statistics = 7
  integer, parameter :: region_count = 7

  type timer_mark
    integer(int64) :: wall_count = 0_int64
    integer(int64) :: wall_rate = 1_int64
    real(real64) :: cpu_seconds = 0.0_real64
  contains
    procedure, pass(mark) :: start => start_timer
  end type timer_mark

  type step_timing
    integer :: step = 0
    character(len=14) :: model_time = ''
    character(len=14) :: forcing_time = ''
    real(real64) :: wall_seconds = 0.0_real64
    real(real64) :: cpu_seconds = 0.0_real64
    real(real64) :: region_wall_seconds(region_count) = 0.0_real64
    real(real64) :: region_cpu_seconds(region_count) = 0.0_real64
    integer :: region_calls(region_count) = 0
  end type step_timing

  type model_timers
    real(real64) :: wall_seconds(region_count) = 0.0_real64
    real(real64) :: cpu_seconds(region_count) = 0.0_real64
    integer :: calls(region_count) = 0
    integer :: step_count = 0
    real(real64) :: recorded_region_wall_seconds(region_count) = 0.0_real64
    real(real64) :: recorded_region_cpu_seconds(region_count) = 0.0_real64
    integer :: recorded_region_calls(region_count) = 0
    type(step_timing), allocatable :: steps(:)
  contains
    procedure, pass(timers) :: initialize => init_timers
    procedure, pass(timers) :: add_elapsed
    procedure, pass(timers) :: record_step
  end type model_timers

  type(model_timers), save :: original_timers

contains

  subroutine init_timers(timers, max_steps)
    class(model_timers), intent(out) :: timers
    integer, intent(in) :: max_steps

    allocate(timers%steps(max_steps))
  end subroutine init_timers

  subroutine start_timer(mark)
    class(timer_mark), intent(out) :: mark

    call system_clock(count=mark%wall_count, count_rate=mark%wall_rate)
    if (mark%wall_rate <= 0_int64) mark%wall_rate = 1_int64
    call cpu_time(mark%cpu_seconds)
  end subroutine start_timer

  subroutine elapsed_since(mark, wall_seconds, cpu_seconds)
    type(timer_mark), intent(in) :: mark
    real(real64), intent(out) :: wall_seconds
    real(real64), intent(out) :: cpu_seconds
    integer(int64) :: wall_now
    real(real64) :: cpu_now

    call system_clock(count=wall_now)
    call cpu_time(cpu_now)
    wall_seconds = real(wall_now - mark%wall_count, real64) / real(mark%wall_rate, real64)
    cpu_seconds = cpu_now - mark%cpu_seconds
  end subroutine elapsed_since

  subroutine add_elapsed(timers, region, mark)
    class(model_timers), intent(inout) :: timers
    integer, intent(in) :: region
    type(timer_mark), intent(in) :: mark
    real(real64) :: wall_seconds
    real(real64) :: cpu_seconds

    call elapsed_since(mark, wall_seconds, cpu_seconds)
    call add_region_time(timers, region, wall_seconds, cpu_seconds)
  end subroutine add_elapsed

  subroutine add_region_time(timers, region, wall_seconds, cpu_seconds)
    type(model_timers), intent(inout) :: timers
    integer, intent(in) :: region
    real(real64), intent(in) :: wall_seconds
    real(real64), intent(in) :: cpu_seconds

    timers%wall_seconds(region) = timers%wall_seconds(region) + wall_seconds
    timers%cpu_seconds(region) = timers%cpu_seconds(region) + cpu_seconds
    timers%calls(region) = timers%calls(region) + 1
  end subroutine add_region_time

  subroutine record_step(timers, step, model_time, forcing_time, mark, wall_seconds, cpu_seconds)
! Capacity is fixed at initialization to keep integration allocation-free. If
! records overflow, aggregate timing remains complete and extra rows are omitted.
    class(model_timers), intent(inout) :: timers
    integer, intent(in) :: step
    character(len=*), intent(in) :: model_time
    character(len=*), intent(in) :: forcing_time
    type(timer_mark), intent(in) :: mark
    real(real64), intent(out) :: wall_seconds
    real(real64), intent(out) :: cpu_seconds
    real(real64) :: region_wall_seconds(region_count)
    real(real64) :: region_cpu_seconds(region_count)
    integer :: region_calls(region_count)

    call elapsed_since(mark, wall_seconds, cpu_seconds)
    call add_region_time(timers, region_total_excluding_init, wall_seconds, cpu_seconds)
    region_wall_seconds = timers%wall_seconds - timers%recorded_region_wall_seconds
    region_cpu_seconds = timers%cpu_seconds - timers%recorded_region_cpu_seconds
    region_calls = timers%calls - timers%recorded_region_calls
    timers%recorded_region_wall_seconds = timers%wall_seconds
    timers%recorded_region_cpu_seconds = timers%cpu_seconds
    timers%recorded_region_calls = timers%calls
    if (timers%step_count < size(timers%steps)) then
      timers%step_count = timers%step_count + 1
      timers%steps(timers%step_count)%step = step
      timers%steps(timers%step_count)%model_time = model_time
      timers%steps(timers%step_count)%forcing_time = forcing_time
      timers%steps(timers%step_count)%wall_seconds = wall_seconds
      timers%steps(timers%step_count)%cpu_seconds = cpu_seconds
      timers%steps(timers%step_count)%region_wall_seconds = region_wall_seconds
      timers%steps(timers%step_count)%region_cpu_seconds = region_cpu_seconds
      timers%steps(timers%step_count)%region_calls = region_calls
    end if
  end subroutine record_step

end module timing_mod
