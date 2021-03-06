module diag_mod

  use ieee_arithmetic
  use parallel_mod
  use mesh_mod
  use data_mod
  use types_mod
  use log_mod
  use params_mod

  implicit none

  private

  public diag_init
  public diag_run
  public diag_final
  public diag_total_energy
  public diag_type
  public diag

  type diag_type
    real total_mass
    real total_energy
    real total_enstrophy
    real, allocatable :: vor(:,:)
    real, allocatable :: div(:,:)
  end type diag_type

  type(diag_type) diag

contains

  subroutine diag_init()

    if (.not. allocated(diag%vor)) call parallel_allocate(diag%vor, half_lon=.true., half_lat=.true.)
    if (.not. allocated(diag%div)) call parallel_allocate(diag%div)

    call log_notice('Diag module is initialized.')

  end subroutine diag_init

  subroutine diag_run(state)

    type(state_type), intent(in) :: state

    real vm1, vp1, um1, up1
    integer i, j
    real hd_cornor(parallel%half_lon_start_idx: parallel%half_lon_end_idx, &
                   parallel%half_lat_start_idx: parallel%half_lat_end_idx )

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        um1 = state%u(i-1,j)
        up1 = state%u(i,j)
        vm1 = state%v(i,j-1) * mesh%half_cos_lat(j-1)
        vp1 = state%v(i,j) * mesh%half_cos_lat(j)
        diag%div(i,j) = (up1 - um1) / coef%full_dlon(j) + (vp1 - vm1) / coef%full_dlat(j)
      end do
    end do
    call parallel_fill_halo(diag%div, all_halo=.true.)

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        um1 = state%u(i,j) * mesh%full_cos_lat(j)
        up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
        vm1 = state%v(i,j)
        vp1 = state%v(i+1,j)
        diag%vor(i,j) = (vp1 - vm1) / coef%half_dlon(j) - (up1 - um1) / coef%half_dlat(j)
      end do
    end do
    call parallel_fill_halo(diag%vor, all_halo=.true.)

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        hd_cornor(i,j) = 0.25 * (state%gd(i,j) + state%gd(i+1,j) + state%gd(i,j+1) + state%gd(i+1,j+1)) / g
      end do
    end do 
    diag%total_enstrophy = 0.0
     do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        diag%total_enstrophy = diag%total_enstrophy + 0.5 *(diag%vor(i,j) + coef%half_f(j))**2 / hd_cornor(i,j) * mesh%dlon * mesh%dlat * mesh%half_cos_lat(j)
      end do
    end do 
    diag%total_enstrophy = diag%total_enstrophy * radius**2

    if (ieee_is_nan(diag%total_enstrophy)) then
      call log_error('Total potential total_enstrophy is NaN!')
    end if


    diag%total_mass = 0.0
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        diag%total_mass = diag%total_mass + mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat * state%gd(i,j)
      end do
    end do
    diag%total_mass = diag%total_mass * radius**2

    if (ieee_is_nan(diag%total_mass)) then
      call log_error('Total mass is NaN!')
    end if

    diag%total_energy = diag_total_energy(state)

    if (ieee_is_nan(diag%total_energy)) then
      call log_error('Total energy is NaN!')
    end if

  end subroutine diag_run

  subroutine diag_final()

    if (allocated(diag%vor)) deallocate(diag%vor)
    if (allocated(diag%div)) deallocate(diag%div)

  end subroutine diag_final

  real function diag_total_energy(state) result(res)

    type(state_type), intent(in) :: state

    integer i, j

    res = 0.0
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + state%iap%u(i,j)**2 * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + state%iap%v(i,j)**2 * mesh%half_cos_lat(j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + (state%gd(i,j) + static%ghs(i,j))**2 * mesh%full_cos_lat(j)
      end do
    end do
    res = res * 0.5 * radius**2
  end function diag_total_energy

end module diag_mod
