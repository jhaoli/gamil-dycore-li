module diffusion_mod

  use parallel_mod
  use params_mod
  use mesh_mod
  use types_mod
  use data_mod
  use diag_mod
  use filter_mod
  use log_mod

  implicit none

  private

  public diffusion_init
  public divergence_diffusion
  public ordinary_diffusion_direct
  public ordinary_diffusion_split
  public ordinary_diffusion_limiter
  public ordinary_diffusion_expand
  public ordinary_diffusion_nonlinear
  public ordinary_diffusion_nonlinear2
  public shapiro_filter
  public diffusion_final

  real, allocatable :: ud(:,:)
  real, allocatable :: vd(:,:)
  real, allocatable :: gdd(:,:)

  real, allocatable :: u(:,:)
  real, allocatable :: v(:,:)
  real, allocatable :: gd(:,:)

  real, allocatable :: diffusion_coef_gd(:)
  real, allocatable :: diffusion_coef_u(:)
  real, allocatable :: diffusion_coef_v(:)

  real, allocatable :: diffusion_coef_gd_x(:)
  real, allocatable :: diffusion_coef_u_x(:)
  real, allocatable :: diffusion_coef_v_x(:)

  real, allocatable :: diffusion_coef_gd_y(:)
  real, allocatable :: diffusion_coef_u_y(:)
  real, allocatable :: diffusion_coef_v_y(:)

  real, allocatable :: gd_x(:,:), gd_y(:,:)
  real, allocatable :: gdd_x(:,:), gdd_y(:,:)

  real, allocatable :: u_x(:,:), u_y(:,:)
  real, allocatable :: ud_x(:,:), ud_y(:,:)

  real, allocatable :: v_x(:,:), v_y(:,:)
  real, allocatable :: vd_x(:,:), vd_y(:,:)

  real, allocatable :: gd_tmp(:,:), u_tmp(:,:), v_tmp(:,:)
  real, allocatable :: half_grad_x(:,:), half_grad_y(:,:)
  real, allocatable :: half_h_grad_x(:,:), half_h_grad_y(:,:)
  ! for nonlinear diffusion
  real, allocatable :: tension_strain_full(:,:), shear_strain_corner(:,:), horizontal_viscosity_full(:,:), horizontal_viscosity_corner(:,:)
contains

  subroutine diffusion_init()

    if (.not. allocated(ud))  call parallel_allocate(ud, half_lon=.true.)
    if (.not. allocated(vd))  call parallel_allocate(vd, half_lat=.true.)
    if (.not. allocated(gdd)) call parallel_allocate(gdd)
    if (.not. allocated(u))   call parallel_allocate(u,  half_lon=.true.)
    if (.not. allocated(v))   call parallel_allocate(v,  half_lat=.true.)
    if (.not. allocated(gd))  call parallel_allocate(gd)

    if (.not. allocated(diffusion_coef_gd)) call parallel_allocate(diffusion_coef_gd)
    if (.not. allocated(diffusion_coef_u))  call parallel_allocate(diffusion_coef_u, half_lon=.true.)
    if (.not. allocated(diffusion_coef_v))  call parallel_allocate(diffusion_coef_v, half_lat=.true.)
    if (.not. allocated(diffusion_coef_gd_x)) call parallel_allocate(diffusion_coef_gd_x)
    if (.not. allocated(diffusion_coef_gd_y)) call parallel_allocate(diffusion_coef_gd_y)
    if (.not. allocated(diffusion_coef_u_x))  call parallel_allocate(diffusion_coef_u_x, half_lon=.true.)
    if (.not. allocated(diffusion_coef_u_y))  call parallel_allocate(diffusion_coef_u_y, half_lon=.true.)
    if (.not. allocated(diffusion_coef_v_x))  call parallel_allocate(diffusion_coef_v_x, half_lat=.true.)
    if (.not. allocated(diffusion_coef_v_y))  call parallel_allocate(diffusion_coef_v_y, half_lat=.true.)
    if (.not. allocated(gd_x)) call parallel_allocate(gd_x)
    if (.not. allocated(gd_y)) call parallel_allocate(gd_y)
    if (.not. allocated(gdd_x)) call parallel_allocate(gdd_x)
    if (.not. allocated(gdd_y)) call parallel_allocate(gdd_y)

    if (.not. allocated(u_x)) call parallel_allocate(u_x)
    if (.not. allocated(u_y)) call parallel_allocate(u_y)
    if (.not. allocated(ud_x)) call parallel_allocate(ud_x)
    if (.not. allocated(ud_y)) call parallel_allocate(ud_y)
    if (.not. allocated(v_x)) call parallel_allocate(v_x)
    if (.not. allocated(v_y)) call parallel_allocate(v_y)
    if (.not. allocated(vd_x)) call parallel_allocate(vd_x)
    if (.not. allocated(vd_y)) call parallel_allocate(vd_y)
    if (.not. allocated(gd_tmp)) call parallel_allocate(gd_tmp)
    if (.not. allocated(u_tmp)) call parallel_allocate(u_tmp)
    if (.not. allocated(v_tmp)) call parallel_allocate(v_tmp)
    if (.not. allocated(half_grad_x)) call parallel_allocate(half_grad_x)
    if (.not. allocated(half_grad_y)) call parallel_allocate(half_grad_y)
    if (.not. allocated(half_h_grad_x)) call parallel_allocate(half_h_grad_x)
    if (.not. allocated(half_h_grad_y)) call parallel_allocate(half_h_grad_y)
    
    ! for nonlinear diffusion
    if (.not. allocated(tension_strain_full)) call parallel_allocate(tension_strain_full)
    if (.not. allocated(shear_strain_corner)) call parallel_allocate(shear_strain_corner)
    if (.not. allocated(horizontal_viscosity_full)) call parallel_allocate(horizontal_viscosity_full)
    if (.not. allocated(horizontal_viscosity_corner)) call parallel_allocate(horizontal_viscosity_corner)
    
    call log_notice('diffusion module is initialized.')

  end subroutine diffusion_init

  subroutine divergence_diffusion(dt, diag, state)
    real, intent(in) :: dt
    type(diag_type), intent(in) :: diag
    type(state_type), intent(inout) :: state

    real, parameter :: vd = 1.0e5
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + dt * vd * (diag%div(i+1,j) - diag%div(i,j)) / coef%full_dlon(j)
      end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u(i,j) = state%u(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j))
      end do
    end do
    do j = parallel%half_lon_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + dt * vd * (diag%div(i,j+1) - diag%div(i,j)) / radius / mesh%dlat
      end do
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v(i,j) = state%v(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1))
      end do
    end do

    call parallel_fill_halo(state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v(:,:), all_halo=.true.)
    call parallel_fill_halo(state%u(:,:),     all_halo=.true.)
    call parallel_fill_halo(state%v(:,:),     all_halo=.true.)
  end subroutine divergence_diffusion

  subroutine ordinary_diffusion_direct(dt, state)
    real, intent(in) :: dt
    type(state_type), intent(inout) :: state

    real sp, np
    integer i, j, order, sign

    u(:,:) = state%iap%u(:,:)
    v(:,:) = state%iap%v(:,:)
    gd(:,:) = state%iap%gd(:,:)


    ! Scalar diffusion:
    !
    ! 2nd order:
    !
    !   ∂ F         1      ∂² F        1      ∂          ∂ F
    !   --- = 𝞶 ---------- ---- + 𝞶 --------- ---(cos(φ) ---)
    !   ∂ t     a² cos²(φ) ∂ λ²     a² cos(φ) ∂ φ        ∂ φ
    !
    !
    ! Vector diffusion:
    !
    ! 2nd order:
    !
    ! ∂² u                 u        2 sin𝞿   ∂ v
    ! ----           - --------- - --------- ---
    ! ∂ λ²             a² cos²𝞿    a² cos²𝞿  ∂ 𝞴
    !
    ! ∂          ∂ v       v        2 sin𝞿   ∂ u
    ! --- cos(φ) --- - --------- + --------- ---
    ! ∂ φ        ∂ φ   a² cos²𝞿    a² cos²𝞿  ∂ 𝞴
    

    do order = 1, diffusion_order / 2
      ! H
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          gdd(i,j) = (gd(i+1,j) - 2 * gd(i,j) + gd(i-1,j)) / coef%full_dlon(j)**2 + &
                     ((gd(i,j+1) - gd(i,j  )) * mesh%half_cos_lat(j  ) - &
                      (gd(i,j  ) - gd(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
        end do
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        sp = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          sp = sp + gd(i,j+1) - gd(i,j)
        end do
        call parallel_zonal_sum(sp)
        gdd(:,j) = sp * mesh%half_cos_lat(j) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      end if
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        np = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          np = np - (gd(i,j) - gd(i,j-1))
        end do
        call parallel_zonal_sum(np)
        gdd(:,j) = np * mesh%half_cos_lat(j-1) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      end if
      ! U
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          ud(i,j) = (u(i+1,j) - 2 * u(i,j) + u(i-1,j)) / coef%full_dlon(j)**2 + &
                    ((u(i,j+1) - u(i,j  )) * mesh%half_cos_lat(j  ) - &
                     (u(i,j  ) - u(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) !- &
                    ! (u(i,j) - mesh%full_sin_lat(j) * (v(i+1,j-1) + v(i+1,j) - v(i,j-1) - v(i,j)) / mesh%dlon) / &
                    !  radius**2 / mesh%full_cos_lat(j)**2
        end do
      end do
      ! V
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = (v(i+1,j) - 2 * v(i,j) + v(i-1,j)) / coef%half_dlon(j)**2
        end do
      end do
      do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = vd(i,j) + ((v(i,j+1) - v(i,j  )) * mesh%full_cos_lat(j+1) - &
                               (v(i,j  ) - v(i,j-1)) * mesh%full_cos_lat(j  )) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j) !- &
                              ! (v(i,j) + mesh%half_sin_lat(j) * (u(i,j) + u(i,j+1) - u(i-1,j) - u(i-1,j+1)) / mesh%dlon) / &
                              !  radius**2 / mesh%half_cos_lat(j)**2
        end do
      end do
      if (parallel%has_south_pole) then
        j = parallel%half_lat_south_pole_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = vd(i,j) + (v(i,j+1) - v(i,j)) * mesh%full_cos_lat(j+1) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j) !- &
                              ! (v(i,j) + mesh%half_sin_lat(j) * (u(i,j+1) - u(i-1,j+1)) / mesh%dlon) / radius**2 / mesh%half_cos_lat(j)**2
        end do
      end if
      if (parallel%has_north_pole) then
        j = parallel%half_lat_north_pole_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = vd(i,j) - (v(i,j) - v(i,j-1)) * mesh%full_cos_lat(j) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j) !- &
                              ! (v(i,j) + mesh%half_sin_lat(j) * (u(i,j) - u(i-1,j)) / mesh%dlon) / radius**2 / mesh%half_cos_lat(j)**2
        end do
      end if
      if (order /= diffusion_order / 2) then
        call parallel_fill_halo(gdd, all_halo=.true.)
        call parallel_fill_halo(ud,  all_halo=.true.)
        call parallel_fill_halo(vd,  all_halo=.true.)
        gd(:,:) = gdd(:,:)
        u(:,:)  = ud(:,:)
        v(:,:)  = vd(:,:)
      end if
    end do

    ! Do FFT filter.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (filter_full_zonal_tend(j)) then
        call filter_array_at_full_lat(j, gdd(:,j))
        call filter_array_at_full_lat(j, ud(:,j))
      end if
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (filter_half_zonal_tend(j)) then
        call filter_array_at_half_lat(j, vd(:,j))
      end if
    end do
    ! diffusion coefficient
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      diffusion_coef_gd(j) = 1.0E05 !  / dt * ( coef%full_dlon(j) / 2 )**diffusion_order 
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      diffusion_coef_v(j) = 1.0E05 ! / dt * ( coef%half_dlon(j) / 2 )**diffusion_order
    end do

    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%gd(i,j) = state%iap%gd(i,j) + sign * dt * diffusion_coef_gd(j) * gdd(i,j)
      end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u(i,j) = state%iap%u(i,j) + sign * dt * diffusion_coef_gd(j) * ud(i,j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v(i,j) = state%iap%v(i,j) + sign * dt * diffusion_coef_v(j) * vd(i,j)
      end do
    end do

    ! Transform from IAP to normal state.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%gd(i,j) = state%iap%gd(i,j)**2
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%iap%u(i,j) * 2.0 / (state%iap%gd(i,j) + state%iap%gd(i+1,j))
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%iap%v(i,j) * 2.0 / (state%iap%gd(i,j) + state%iap%gd(i,j+1))
      end do
    end do
    call parallel_fill_halo(state%gd, all_halo=.true.)
    call parallel_fill_halo(state%u,  all_halo=.true.)
    call parallel_fill_halo(state%v,  all_halo=.true.)

    ! call iap_transform(state)
  end subroutine ordinary_diffusion_direct

  subroutine ordinary_diffusion_split(dt, state)
    real, intent(in) :: dt
    type(state_type), intent(inout) :: state

    real sp, np
    integer i, j, order, sign
    integer i0

    ! u(:,:) = state%iap%u(:,:)
    ! v(:,:) = state%iap%v(:,:)
    ! gd(:,:) = state%iap%gd(:,:)
    u(:,:)  = state%u(:,:)
    v(:,:)  = state%v(:,:)
    gd(:,:) = state%gd(:,:)

    gd_x(:,:) = gd(:,:)
    gd_y(:,:) = gd(:,:)
    u_x(:,:) = u(:,:)
    u_y(:,:) = u(:,:)
    v_x(:,:) = v(:,:)
    v_y(:,:) = v(:,:)
      ! H on longitude
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do order = 1, diffusion_order / 2
          do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
            gdd_x(i,j) = (gd_x(i+1,j) - 2 * gd_x(i,j) + gd_x(i-1,j)) / coef%full_dlon(j)**2          
          end do
          if (order /= diffusion_order / 2) then
            call parallel_fill_halo(gdd_x, all_halo=.true.)
            gd_x(:,j) = gdd_x(:,j)
          end if
        end do 
      end do

      ! H on latitude
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        i0 = i + mesh%num_full_lon / 2 
        if (i0 > mesh%num_full_lon / 2) then 
          i0 = i0 - mesh%num_full_lon / 2
        end if 
        do order = 1, diffusion_order / 2
          do j = parallel%full_lat_start_idx , parallel%full_lat_end_idx 
            if(j == parallel%full_lat_start_idx) then
              gd_y(i,j-1) = gd_y(i0,j+1)
              mesh%half_cos_lat(j-1) = mesh%half_cos_lat(j+1)
            elseif(j == parallel%full_lat_end_idx) then
              gd_y(i,j+1) = gd_y(i0,j-1)
            end if 
            gdd_y(i,j) = ((gd_y(i,j+1) - gd_y(i,j  )) * mesh%half_cos_lat(j  ) - &
                          (gd_y(i,j  ) - gd_y(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
          end do
          if (order /= diffusion_order / 2) then
            call parallel_fill_halo(gdd_y, all_halo=.true.)
            gd_y(i,:) = gdd_y(i,:)
          end if
        end do 
      end do 

      ! U on longitude
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do order = 1, diffusion_order / 2
          do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
            ud_x(i,j) = (u_x(i+1,j) - 2 * u_x(i,j) + u_x(i-1,j)) / coef%full_dlon(j)**2 
          end do
          if (order /= diffusion_order / 2) then
            call parallel_fill_halo(ud_x, all_halo=.true.)
            u_x(:,j) = ud_x(:,j)
          end if
        end do 
      end do
      ! U on latitude
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        do order = 1, diffusion_order / 2
          do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
            ud_y(i,j) = ((u_y(i,j+1) - u_y(i,j  )) * mesh%half_cos_lat(j  ) - &
                         (u_y(i,j  ) - u_y(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)  
          end do
          if (order /= diffusion_order / 2) then
            call parallel_fill_halo(ud_y, all_halo=.true.)
            u_y(i,:) = ud_y(i,:)
          end if
        end do 
      end do 
     
      ! V on longitude
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do order = 1, diffusion_order / 2
          do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
            vd_x(i,j) = (v_x(i+1,j) - 2 * v_x(i,j) + v_x(i-1,j)) / coef%full_dlon(j)**2
          end do
          if (order /= diffusion_order / 2) then
            call parallel_fill_halo(vd_x, all_halo=.true.)
            v_x(:,j) = vd_x(:,j)
          end if
        end do 
      end do
      ! V on latitude
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        do order = 1, diffusion_order / 2
          do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
            i0 = i + mesh%num_full_lon / 2
            if (i0 > mesh%num_full_lon / 2) then
              i0 = i0 - mesh%num_full_lon / 2
            end if
            if (j == parallel%half_lat_start_idx) then
              v_y(i,j-1) = v_y(i0,j)
            else if (j == parallel%half_lat_end_idx) then
              v_y(i,j+1) = v_y(i0,j)
              mesh%full_cos_lat(j+1) = mesh%full_cos_lat(j-1)
            end if 
            vd_y(i,j) = ((v_y(i,j+1) - v_y(i,j  )) * mesh%full_cos_lat(j+1) - &
                         (v_y(i,j  ) - v_y(i,j-1)) * mesh%full_cos_lat(j)) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j)
          end do
          if (order /= diffusion_order / 2) then
            call parallel_fill_halo(vd_y, all_halo=.true.)
            v_y(i,:) = vd_y(i,:)
          end if
        end do 
      end do
        
    ! diffusion coefficient
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      diffusion_coef_gd_x(j) = 1.0 / dt * (coef%full_dlon(j) / 2)**diffusion_order
      diffusion_coef_gd_y(j) = radius**4 / dt / ((1+mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2 * (4 / mesh%dlat**2) +&
                                16 / mesh%dlat**4 ) 

      diffusion_coef_u_x(j) = diffusion_coef_gd_x(j)
      diffusion_coef_u_y(j) = diffusion_coef_gd_y(j)
    end do 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      diffusion_coef_v_x(j) = 1.0 / dt * (coef%half_dlon(j) / 2)**diffusion_order
      diffusion_coef_v_y(j) = radius**4 / dt / ((1+mesh%half_cos_lat(j)**2) / mesh%half_cos_lat(j)**2 * (4 / mesh%dlat**2) +&
                                16 / mesh%dlat**4 )
    end do 

    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      ! do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
      !   state%iap%gd(i,j) = state%iap%gd(i,j) + sign * dt * (diffusion_coef_gd_x(j) * gdd_x(i,j) + diffusion_coef_gd_y(j) * gdd_y(i,j) )
      ! end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + sign * dt * (diffusion_coef_gd_x(j) * ud_x(i,j) + diffusion_coef_gd_y(j) * ud_y(i,j) )
      end do 
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + sign * dt * (diffusion_coef_v_x(j) * vd_x(i,j) + diffusion_coef_v_y(j) * vd_y(i,j) )
      end do
    end do

    ! Transform from IAP to normal state.
    ! do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
    !   do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
    !     state%gd(i,j) = state%iap%gd(i,j)**2
    !   end do
    ! end do
    ! do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
    !   do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
    !     state%u(i,j) = state%iap%u(i,j) * 2.0 / (state%iap%gd(i,j) + state%iap%gd(i+1,j))
    !   end do
    ! end do
    ! do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
    !   do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
    !     state%v(i,j) = state%iap%v(i,j) * 2.0 / (state%iap%gd(i,j) + state%iap%gd(i,j+1))
    !   end do
    ! end do
    call iap_transform(state)

  end subroutine ordinary_diffusion_split

  subroutine ordinary_diffusion_limiter(dt, state)

    real, intent(in) :: dt
    type(state_type), intent(inout) :: state

    real sp, np
    integer i, j, order, sign

    u(:,:) = state%iap%u(:,:)
    v(:,:) = state%iap%v(:,:)
    gd(:,:) = state%iap%gd(:,:)

    gd_tmp(:,:) = gd(:,:)  
    gd_x(:,:) = gd(:,:)
    gd_y(:,:) = gd(:,:)

    u_tmp(:,:) = u(:,:)
    u_x(:,:) = u(:,:)
    u_y(:,:) = u(:,:)

    v_tmp(:,:) = v(:,:)
    v_x(:,:) = v(:,:)
    v_y(:,:) = v(:,:)
    !====================H==============================================
    !1) calculate even order Laplacian
    ! H on longitude
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do order = 1, diffusion_order / 2 -1
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          gdd_x(i,j) = (gd_x(i+1,j) - 2 * gd_x(i,j) + gd_x(i-1,j)) / coef%full_dlon(j)**2          
        end do
        if (order /= diffusion_order / 2 -1) then
          call parallel_fill_halo(gdd_x, all_halo=.true.)
          gd_x(:,j) = gdd_x(:,j)
        end if
      end do 
    end do
    ! H on latitude
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
      do order = 1, diffusion_order / 2 - 1
        do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
          gdd_y(i,j) = ((gd_y(i,j+1) - gd_y(i,j  )) * mesh%half_cos_lat(j  ) - &
                        (gd_y(i,j  ) - gd_y(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
        end do
        if (parallel%has_south_pole) then
          j = parallel%full_lat_south_pole_idx
          gdd_y(i,j) = (gd_y(i,j+1) - gd_y(i,j)) * mesh%half_cos_lat(j) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) 
        end if 
        if (parallel%has_north_pole) then
          j = parallel%full_lat_north_pole_idx
          gdd_y(i,j) = -(gd_y(i,j) - gd_y(i,j-1)) * mesh%half_cos_lat(j-1) / coef%full_dlat(j-1)**2 * mesh%full_cos_lat(j-1)
        end if 
        if (order /= diffusion_order / 2 - 1) then
          call parallel_fill_halo(gdd_y, all_halo=.true.)
          gd_y(i,:) = gdd_y(i,:)
        end if
      end do         
    end do 

    !2) Calculate gradient of Laplacian above
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
        half_grad_x(i,j) = (gdd_x(i+1,j) - gdd_x(i,j)) / coef%full_dlon(j) 
      end do   
    end do 
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        half_grad_y(i,j) = (gdd_y(i,j+1) - gdd_y(i,j)) / radius / mesh%dlat 
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        half_grad_y(i,j) = (gdd_y(i,j+1) - gdd_y(i,j)) / radius / mesh%dlat 
      end if 
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        half_grad_y(i,j) = (gdd_y(i,j) - gdd_y(i,j-1)) / radius / mesh%dlat 
      end if 
    end do 
    !3) Calculate gradient of the forecast variable 
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
        half_h_grad_x(i,j) = (gd_tmp(i+1,j) - gd_tmp(i,j)) / coef%full_dlon(j) 
      end do   
    end do
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        half_h_grad_y(i,j) = (gd_tmp(i,j+1) - gd_tmp(i,j)) / radius / mesh%dlat 
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        half_h_grad_y(i,j) = (gd_tmp(i,j+1) - gd_tmp(i,j)) / radius / mesh%dlat 
      end if 
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        half_h_grad_y(i,j) = (gd_tmp(i,j) - gd_tmp(i,j-1)) / radius / mesh%dlat 
      end if 
    end do        
    !4) Decide the sign of diffusion flux to be zero or not
    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
        if(sign * half_grad_x(i,j) * half_h_grad_x(i,j) <= 0) then
          half_grad_x(i,j) = 0.0
        end if 
      end do   
    end do
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        if(sign * half_grad_y(i,j) * half_h_grad_y(i,j) <= 0) then
          half_grad_y(i,j) = 0.0
        end if 
      end do
    end do   
    !5) Calculate the time tendency 
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
         gdd_x(i,j) = (half_grad_x(i,j) - half_grad_x(i-1,j)) / coef%full_dlon(j) 
      end do   
    end do
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
         gdd_y(i,j) = (half_grad_y(i,j) - half_grad_y(i,j-1)) / radius / mesh%dlat
      end do
    end do      
    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (half_grad_y(i,j+1) - half_grad_y(i,j))
      end do
      call parallel_zonal_sum(sp)
      gdd_y(:,j) = sp / radius / mesh%dlat / mesh%num_full_lon
    end if
    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np + (half_grad_y(i,j) - half_grad_y(i,j-1))
      end do
      call parallel_zonal_sum(np)
      gdd_y(:,j) = np / radius / mesh%dlat  / mesh%num_full_lon
      end if
    gdd(:,:) = gdd_x(:,:) + gdd_y(:,:)
    !====================U==============================================
    !1) calculate even order Laplacian
    ! U on longitude
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do order = 1, diffusion_order / 2 - 1
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          ud_x(i,j) = (u_x(i+1,j) - 2 * u_x(i,j) + u_x(i-1,j)) / coef%full_dlon(j)**2 
        end do
        if (order /= diffusion_order / 2 - 1) then
          call parallel_fill_halo(ud_x, all_halo=.true.)
          u_x(:,:) = ud_x(:,:)
        end if
      end do 
    end do
    ! U on latitude
    do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
      do order = 1, diffusion_order / 2 - 1
        do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
          ud_y(i,j) = ((u_y(i,j+1) - u_y(i,j  )) * mesh%half_cos_lat(j  ) - &
                       (u_y(i,j  ) - u_y(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)  
        end do
        if(parallel%has_south_pole) then
          j = parallel%full_lat_south_pole_idx
          ud_y(i,j) = (u_y(i,j+1) - u_y(i,j)) * mesh%half_cos_lat(j) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)  
        end if 
        if(parallel%has_north_pole) then
          j = parallel%full_lat_north_pole_idx
          ud_y(i,j) = -(u_y(i,j) - u_y(i,j-1)) * mesh%half_cos_lat(j-1) / coef%full_dlat(j-1)**2 * mesh%full_cos_lat(j-1)
        end if 
        if (order /= diffusion_order / 2 - 1) then
          call parallel_fill_halo(ud_y, all_halo=.true.)
          u_y(i,:) = ud_y(i,:)
        end if
      end do 
    end do 
    !2) Calculate gradient of Laplacian above
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx  
        half_grad_x(i,j) = (ud_x(i+1,j) - ud_x(i,j)) / coef%full_dlon(j) 
      end do   
    end do 
    do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx 
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        half_grad_y(i,j) = (ud_y(i,j+1) - ud_y(i,j)) / radius / mesh%dlat 
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        half_grad_y(i,j) = (ud_y(i,j+1) - ud_y(i,j)) / radius / mesh%dlat 
      end if 
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        half_grad_y(i,j) = (ud_y(i,j) - ud_y(i,j-1)) / radius / mesh%dlat 
      end if 
    end do 
    !3) Calculate gradient of the variable 
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx  
        half_h_grad_x(i,j) = (u_tmp(i+1,j) - u_tmp(i,j)) / coef%full_dlon(j) 
      end do   
    end do
    do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx 
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        half_h_grad_y(i,j) = (u_tmp(i,j+1) - u_tmp(i,j)) / radius / mesh%dlat 
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        half_h_grad_y(i,j) = (u_tmp(i,j+1) - u_tmp(i,j)) / radius / mesh%dlat 
      end if 
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        half_h_grad_y(i,j) = (u_tmp(i,j) - u_tmp(i,j-1)) / radius / mesh%dlat 
      end if 
    end do        
    !4) Decide the sign of diffusion flux to be zero or not
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx  
        if(sign * half_grad_x(i,j) * half_h_grad_x(i,j) <=0 ) then
          half_grad_x(i,j) = 0.0
        end if 
      end do   
    end do
    do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx 
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        if(sign * half_grad_y(i,j) * half_h_grad_y(i,j) <=0 ) then
          half_grad_y(i,j) = 0.0
        end if 
      end do
    end do   
    !5) Calculate the time tendency 
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx  
         ud_x(i,j) = (half_grad_x(i,j) - half_grad_x(i-1,j)) / coef%full_dlon(j) 
      end do   
    end do
    do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx 
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
         ud_y(i,j) = (half_grad_y(i,j) - half_grad_y(i,j-1)) / radius / mesh%dlat
      end do
    end do 
    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (half_grad_y(i,j+1) - half_grad_y(i,j))
      end do
      call parallel_zonal_sum(sp)
      ud_y(:,j) = sp / radius / mesh%dlat / mesh%num_full_lon
    end if
    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np + (half_grad_y(i,j) - half_grad_y(i,j-1))
      end do
      call parallel_zonal_sum(np)
      ud_y(:,j) = np / radius / mesh%dlat / mesh%num_full_lon
    end if
    
    ud(:,:) = ud_x(:,:) + ud_y(:,:)
    !====================V==============================================
    !1) calculate even order Laplacian
    ! V on longitude
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do order = 1, diffusion_order / 2 - 1
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd_x(i,j) = (v_x(i+1,j) - 2 * v_x(i,j) + v_x(i-1,j)) / coef%full_dlon(j)**2
        end do
        if (order /= diffusion_order / 2 - 1) then
          call parallel_fill_halo(vd_x, all_halo=.true.)
          v_x(:,j) = vd_x(:,j)
        end if
      end do 
    end do
    ! V on latitude
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
      do order = 1, diffusion_order / 2 - 1
        do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
          vd_y(i,j) = ((v_y(i,j+1) - v_y(i,j  )) * mesh%full_cos_lat(j+1) - &
                     (v_y(i,j  ) - v_y(i,j-1)) * mesh%full_cos_lat(j)) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j)
        end do
        if(parallel%has_south_pole) then
          j = parallel%full_lat_south_pole_idx
          vd_y(i,j) = (v_y(i,j+1) - v_y(i,j)) * mesh%full_cos_lat(j+1) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j)  
        end if 
        if(parallel%has_north_pole) then
          j = parallel%full_lat_north_pole_idx
          vd_y(i,j) = -(v_y(i,j) - v_y(i,j-1)) * mesh%full_cos_lat(j) / coef%half_dlat(j-1)**2 * mesh%half_cos_lat(j-1)
        end if 
        if (order /= diffusion_order / 2 - 1) then
          call parallel_fill_halo(vd_y, all_halo=.true.)
          v_y(i,:) = vd_y(i,:)
        end if
      end do 
    end do
    !2) Calculate gradient of Laplacian above
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
        half_grad_x(i,j) = (vd_x(i+1,j) - vd_x(i,j)) / coef%half_dlon(j) 
      end do   
    end do 
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
        half_grad_y(i,j) = (vd_y(i,j+1) - vd_y(i,j)) / radius / mesh%dlat 
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        half_grad_y(i,j) = (vd_y(i,j+1) - vd_y(i,j)) / radius / mesh%dlat 
      end if 
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        half_grad_y(i,j) = (vd_y(i,j) - vd_y(i,j-1)) / radius / mesh%dlat 
      end if 
    end do 
    !3) Calculate gradient of the variable 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
        half_h_grad_x(i,j) = (v_tmp(i+1,j) - v_tmp(i,j)) / coef%half_dlon(j) 
      end do   
    end do
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
        half_h_grad_y(i,j) = (v_tmp(i,j+1) - v_tmp(i,j)) / radius / mesh%dlat 
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        half_h_grad_y(i,j) = (v_tmp(i,j+1) - v_tmp(i,j)) / radius / mesh%dlat 
      end if 
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        half_h_grad_y(i,j) = (v_tmp(i,j) - v_tmp(i,j-1)) / radius / mesh%dlat 
      end if 
    end do        
    !4) Decide the sign of diffusion flux to be zero or not
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
        if(sign * half_grad_x(i,j) * half_h_grad_x(i,j) <= 0) then
          half_grad_x(i,j) = 0.0
        end if 
      end do   
    end do
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        if(sign * half_grad_y(i,j) * half_h_grad_y(i,j) <= 0) then
          half_grad_y(i,j) = 0.0
        end if 
      end do
    end do   
    !5) Calculate the time tendency 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
         vd_x(i,j) = (half_grad_x(i,j) - half_grad_x(i-1,j)) / coef%half_dlon(j) 
      end do   
    end do
    do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
      do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
         vd_y(i,j) = (half_grad_y(i,j) - half_grad_y(i,j-1)) / radius / mesh%dlat
      end do
    end do 
    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + (half_grad_y(i,j+1) - half_grad_y(i,j))
      end do
      call parallel_zonal_sum(sp)
      vd_y(:,j) = sp / radius / mesh%dlat / mesh%num_full_lon
    end if
    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np + (half_grad_y(i,j) - half_grad_y(i,j-1))
      end do
      call parallel_zonal_sum(np)
      vd_y(:,j) = np / radius / mesh%dlat / mesh%num_full_lon
    end if         
    vd(:,:) = vd_x(:,:) + vd_y(:,:)
    !====================================================    
    ! Do FFT filter.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (filter_full_zonal_tend(j)) then
        call filter_array_at_full_lat(j, gdd(:,j))
        call filter_array_at_full_lat(j, ud(:,j))
      end if
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (filter_half_zonal_tend(j)) then
        call filter_array_at_half_lat(j, vd(:,j))
      end if
    end do  
    ! deffusion coefficient
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      ! diffusion_coef_gd(j) = 1.0 / dt * ( coef%full_dlon(j) / 2 )**diffusion_order 
      diffusion_coef_gd(j) = 1.0 / dt * (coef%full_dlon(j)**diffusion_order * coef%full_dlat(j)**diffusion_order) / &
                                        (coef%full_dlon(j)**diffusion_order + coef%full_dlat(j)**diffusion_order) * &
                                         0.5**diffusion_order
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      ! diffusion_coef_v(j) = 1.0 / dt * ( coef%half_dlon(j) / 2 )**diffusion_order
      diffusion_coef_v(j) = 1.0 / dt * (coef%half_dlon(j)**diffusion_order * coef%half_dlat(j)**diffusion_order) / &
                                       (coef%half_dlon(j)**diffusion_order + coef%half_dlat(j)**diffusion_order) * &
                                       0.5**diffusion_order
    end do

    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%gd(i,j) = state%iap%gd(i,j) + sign * dt * diffusion_coef_gd(j) * gdd(i,j)
      end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u(i,j) = state%iap%u(i,j) + sign * dt * diffusion_coef_gd(j) * ud(i,j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v(i,j) = state%iap%v(i,j) + sign * dt * diffusion_coef_v(j) * vd(i,j)
      end do
    end do

    ! Transform from IAP to normal state.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%gd(i,j) = state%iap%gd(i,j)**2
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%iap%u(i,j) * 2.0 / (state%iap%gd(i,j) + state%iap%gd(i+1,j))
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%iap%v(i,j) * 2.0 / (state%iap%gd(i,j) + state%iap%gd(i,j+1))
      end do
    end do
  end subroutine ordinary_diffusion_limiter

  subroutine ordinary_diffusion_expand(dt, state)
    real, intent(in) :: dt
    type(state_type), intent(inout) :: state
    real sp, np, beta_x, beta_y
    integer i, j, order, sign, i0

    u(:,:) = state%u(:,:)
    v(:,:) = state%v(:,:)
    gd(:,:) = state%gd(:,:)
   
    ! H
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdd_x(i,j) = (gd(i+2,j) - 4 * gd(i+1,j) + 6 * gd(i,j) - 4 * gd(i-1,j) + gd(i-2,j)) / coef%full_dlon(j)**4
        i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_full_lon / 2) then
          i0 = i0 - mesh%num_full_lon / 2
        end if 

        if(j == parallel%full_lat_start_idx_no_pole) then ! the next-to-last circle on south
          gd(i,j-2) = gd(i0,j)
        elseif(j == parallel%full_lat_south_pole_idx) then ! the south pole 
          gd(i,j-1) = gd(i0,j+1)
          gd(i,j-2) = gd(i0,j+2)
        elseif(j == parallel%full_lat_end_idx_no_pole) then   ! the next-to-last circle on north
          gd(i,j+2) = gd(i0,j)
        elseif(j == parallel%full_lat_north_pole_idx) then ! the north pole 
          gd(i,j+1) = gd(i0,j-1)
          gd(i,j+2) = gd(i0,j-2)
        end if 
        gdd_y(i,j) = -1 / radius**4 * (mesh%full_sin_lat(j) / mesh%full_cos_lat(j)**3 *&
                    (gd(i,j+1) - gd(i,j-1)) / (2 * mesh%dlat) + &
                    (1 + mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2 *&
                    (gd(i,j+1) - 2 * gd(i,j) + gd(i,j-1)) / mesh%dlat**2 +&
                    2 * mesh%full_sin_lat(j) / mesh%full_cos_lat(j) * &
                    (gd(i,j+2) - gd(i,j-2) + 2 * gd(i,j-1) -2 * gd(i,j+1)) / (2 * mesh%dlat**3) -&
                    (gd(i,j+2) - 4 * gd(i,j+1) + 6 * gd(i,j) - 4 * gd(i,j-1) + gd(i,j-2)) / mesh%dlat**4)
      end do
    end do
    !==============================================   
    ! U
    do j = parallel%full_lat_start_idx_no_pole , parallel%full_lat_end_idx_no_pole 
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ud_x(i,j) = (u(i+2,j) - 4 * u(i+1,j) + 6 * u(i,j) - 4 * u(i-1,j) + u(i-2,j)) / coef%full_dlon(j)**4
        i0 = i + mesh%num_half_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_half_lon / 2) then
          i0 = i0 - mesh%num_half_lon / 2
        end if 
        if (j == parallel%full_lat_start_idx_no_pole) then
          u(i,j-2) = u(i0,j)
        elseif(j == parallel%full_lat_end_idx_no_pole) then
          u(i,j+2) = u(i0,j)
        end if 
        ud_y(i,j) = -1 / radius**4 *(mesh%full_sin_lat(j) / mesh%full_cos_lat(j)**3 * &
                    (u(i,j+1) - u(i,j-1)) / (2.0 * mesh%dlat) + &
                    (1.0 + mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2* &
                    (u(i,j+1) - 2.0 * u(i,j) + u(i,j-1)) / mesh%dlat**2 + &
                    2 * mesh%full_sin_lat(j) / mesh%full_cos_lat(j) * &
                    (u(i,j+2) - u(i,j-2) + 2.0 * u(i,j-1) - 2.0 * u(i,j+1)) / (2.0 * mesh%dlat**3) - &
                    (u(i,j+2) - 4.0 * u(i,j+1) + 6.0 * u(i,j) - 4.0 * u(i,j-1) + u(i,j-2)) / mesh%dlat**4)
      end do
    end do   
    ! !==============================================
    !V 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd_x(i,j) = (v(i+2,j) - 4 * v(i+1,j) + 6 * v(i,j) - 4 * v(i-1,j) + v(i-2,j)) / coef%half_dlon(j)**4
        i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_full_lon / 2) then
          i0 = i0 - mesh%num_full_lon / 2
        end if 
        if (j == parallel%half_lat_start_idx_no_pole) then      !! the next-to-last circle on south 4th
          v(i,j-2) = v(i0,j-1)
        elseif(j == parallel%half_lat_start_idx) then
          v(i,j-1) = v(i0,j)
          v(i,j-2) = v(i0,j+1)
        elseif(j == parallel%half_lat_end_idx_no_pole ) then
          v(i,j+2) = v(i0,j+1)
        elseif(j == parallel%half_lat_end_idx) then
          v(i,j+1) = v(i0,j)
          v(i,j+2) = v(i0,j-1)
        end if 
        vd_y(i,j) = -1 / radius**4 *(mesh%half_sin_lat(j) / mesh%half_cos_lat(j)**3 * (v(i,j+1) - v(i,j-1))/(2*mesh%dlat) +&
                    (1+mesh%half_cos_lat(j)**2) / mesh%half_cos_lat(j)**2 * (v(i,j+1) - 2* v(i,j) + v(i,j-1)) / mesh%dlat**2 +&
                    2 * mesh%half_sin_lat(j) / mesh%half_cos_lat(j) * (v(i,j+2) - v(i,j-2) + 2 * v(i,j-1) - 2 * v(i,j+1))/(2*mesh%dlat**3) -&
                    (v(i,j+2) - 4 * v(i,j+1) + 6 * v(i,j) - 4 * v(i,j-1) + v(i,j-2)) / mesh%dlat**4 )
      end do 
    end do 

    !============================================ 
    ! diffusion coefficient
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      ! beta_y = (1+mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2 * 1.0E-03 * 1.0E-01
     !  beta_x = 1 / mesh%full_cos_lat(j) ** 4 * (mesh%full_sin_lat(j)**2 / (1 + mesh%full_sin_lat(j)**2))**4 * 1.0E-07
      diffusion_coef_gd_x(j) = 1.0 / dt * (coef%full_dlon(j) / 2)**diffusion_order !* beta_x
      diffusion_coef_gd_y(j) = radius**4 / dt / ((1+mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2 * (4 / mesh%dlat**2) +&
                                16 / mesh%dlat**4 ) !* beta_y  

      diffusion_coef_u_x(j) = diffusion_coef_gd_x(j)
      diffusion_coef_u_y(j) = diffusion_coef_gd_y(j)
    end do 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      ! beta_y = (1+mesh%half_cos_lat(j)**2) / mesh%half_cos_lat(j)**2 * 1.0E-03 * 1.0E-01
      ! beta_x = 1 / mesh%half_cos_lat(j) ** 4 * (mesh%half_sin_lat(j)**2 / (1 + mesh%half_sin_lat(j)**2))**4 * 1.0E-07
      diffusion_coef_v_x(j) = 1.0 / dt * (coef%half_dlon(j) / 2)**diffusion_order !* beta_x
      diffusion_coef_v_y(j) = radius**4 / dt / ((1+mesh%half_cos_lat(j)**2) / mesh%half_cos_lat(j)**2 * (4 / mesh%dlat**2) +&
                                16 / mesh%dlat**4 ) !* beta_y 
    end do 
    !============================================ 

    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
       ! do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
       !   state%gd(i,j) = state%gd(i,j) + sign * dt * (diffusion_coef_gd_x(j) * gdd_x(i,j) + diffusion_coef_gd_y(j) * gdd_y(i,j))
       ! end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + sign * dt * (diffusion_coef_u_x(j) * ud_x(i,j) + diffusion_coef_u_y(j) * ud_y(i,j))
        ! state%u(i,j) = state%u(i,j) + sign * dt * diffusion_coef_u_y(j) * ud_y(i,j)
      end do 
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + sign * dt * (diffusion_coef_v_x(j) * vd_x(i,j) + diffusion_coef_v_y(j) * vd_y(i,j))
        ! state%v(i,j) = state%v(i,j) + sign * dt * diffusion_coef_v_y(j) * vd_y(i,j)
      end do
    end do

    call iap_transform(state)

  end subroutine ordinary_diffusion_expand

  subroutine ordinary_diffusion_nonlinear(dt, state)
    ! Smagorinsky(1963) Nonlinear second order horizontal diffusion
    real, intent(in) :: dt
    type(state_type), intent(inout) :: state
    real :: length_scale, shear_strain_tmp, tension_strain_tmp, deformation_flow
    real :: sp, np
    integer :: i, j, order, i0, sign
    real :: um1, up1, vm1, vp1
  
    u(:,:) = state%u(:,:)
    v(:,:) = state%v(:,:)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_full_lon / 2) then
          i0 = i0 - mesh%num_full_lon / 2
        end if 
        um1 = state%u(i-1,j)
        up1 = state%u(i,j)
        vm1 = state%v(i,j-1) * mesh%half_cos_lat(j-1)
        vp1 = state%v(i,j) * mesh%half_cos_lat(j)
        if (j == parallel%full_lat_start_idx) then
!           vm1 = state%v(i0,j) * mesh%half_cos_lat(j)
          vm1 = -state%v(i0,j) * cos(mesh%half_lat(j) - mesh%dlat)
        elseif (j == parallel%full_lat_end_idx) then
!           vp1 = state%v(i0,j-1) * mesh%half_cos_lat(j-1)
          vp1 = - state%v(i0,j-1) * cos(mesh%half_lat(j-1) + mesh%dlat)
        end if 
        tension_strain_full(i,j) = (up1 - um1) / coef%full_dlon(j) - (vp1 - vm1) / coef%full_dlat(j)
      end do 
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        um1 = state%u(i,j) * mesh%full_cos_lat(j)
        up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
        vm1 = state%v(i,j)
        vp1 = state%v(i+1,j)
        shear_strain_corner(i,j) = (vp1 - vm1) / coef%half_dlon(j) + (up1 - um1) / coef%half_dlat(j)
      end do
    end do
    call parallel_fill_halo(tension_strain_full, all_halo=.true.)
    call parallel_fill_halo(shear_strain_corner, all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        if (j == parallel%full_lat_start_idx) then
          shear_strain_tmp = (shear_strain_corner(i,j)**2 + shear_strain_corner(i-1,j)**2) * 0.5
        elseif (j == parallel%full_lat_end_idx) then
          shear_strain_tmp = (shear_strain_corner(i,j-1)**2 + shear_strain_corner(i-1,j-1)**2) * 0.5
        else 
          shear_strain_tmp = (shear_strain_corner(i-1,j)**2 + shear_strain_corner(i-1,j-1)**2 + shear_strain_corner(i,j)**2 + shear_strain_corner(i,j-1)**2) * 0.25 
        end if 
        length_scale = radius * mesh%full_cos_lat(j) * mesh%dlon     
        deformation_flow = sqrt( tension_strain_full(i,j)**2 + shear_strain_tmp )
        horizontal_viscosity_full(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
      end do 
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        length_scale = radius * mesh%half_cos_lat(j) * mesh%dlon
        tension_strain_tmp = (tension_strain_full(i,j)**2 + tension_strain_full(i+1,j)**2 + tension_strain_full(i,j+1)**2 + tension_strain_full(i+1,j+1)**2) * 0.25
        deformation_flow = sqrt( tension_strain_tmp + shear_strain_corner(i,j)**2 )
        horizontal_viscosity_corner(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
      end do 
    end do 
    call parallel_fill_halo(horizontal_viscosity_full,   all_halo=.true.)
    call parallel_fill_halo(horizontal_viscosity_corner, all_halo=.True.)
    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ud(i,j) = (horizontal_viscosity_full(i+1,j) * tension_strain_full(i+1,j) - horizontal_viscosity_full(i,j) * tension_strain_full(i,j)) / coef%full_dlon(j) + &
                  (horizontal_viscosity_corner(i,j) * shear_strain_corner(i,j) - horizontal_viscosity_corner(i,j-1) * shear_strain_corner(i,j-1)) / radius / mesh%dlat
      end do
    end do   
    !V 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx        
        vd(i,j) = (horizontal_viscosity_corner(i,j) * shear_strain_corner(i,j) - horizontal_viscosity_corner(i-1,j) * shear_strain_corner(i-1,j) ) / coef%half_dlon(j) -&
                  (horizontal_viscosity_full(i,j+1) * tension_strain_full(i,j+1) - horizontal_viscosity_full(i,j) * tension_strain_full(i,j)) / radius / mesh%dlat
      end do 
    end do   

    call parallel_fill_halo(ud,  all_halo=.true.)
    call parallel_fill_halo(vd,  all_halo=.true.)

    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + sign * dt * ud(i,j) 
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + sign * dt * vd(i,j) 
      end do 
    end do 

    call parallel_fill_halo(state%u,  all_halo=.true.)
    call parallel_fill_halo(state%v,  all_halo=.true.)

    call iap_transform(state)
  end subroutine ordinary_diffusion_nonlinear
  
  subroutine ordinary_diffusion_nonlinear2(dt, state)
    !  David L. Williamson. 1978.
    ! The Relative Importance of Resolution, Accuracy and Diffusion in 
    ! Short-Range Forecasts with the NCAR Global Circulation Model
    real, intent(in) :: dt
    type(state_type), intent(inout) :: state
    real :: length_scale, shear_strain_tmp, tension_strain_tmp, deformation_flow
    real :: sp, np
    integer :: i, j, order, i0, sign
    real :: um1, up1, vm1, vp1

    u(:,:) = state%u(:,:)
    v(:,:) = state%v(:,:)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_full_lon / 2) then
          i0 = i0 - mesh%num_full_lon / 2
        end if 
        um1 = state%u(i-1,j)
        up1 = state%u(i,j)
        vm1 = state%v(i,j-1) * mesh%half_cos_lat(j-1)
        vp1 = state%v(i,j) * mesh%half_cos_lat(j)
        if (j == parallel%full_lat_start_idx) then
          vm1 = state%v(i0,j) * mesh%half_cos_lat(j)
        elseif (j == parallel%full_lat_end_idx) then
          vp1 = state%v(i0,j-1) * mesh%half_cos_lat(j-1)
        end if 
        tension_strain_full(i,j) = (up1 - um1) / coef%full_dlon(j) - (vp1 - vm1) / coef%full_dlat(j)
      end do 
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        um1 = state%u(i,j) * mesh%full_cos_lat(j)
        up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
        vm1 = state%v(i,j)
        vp1 = state%v(i+1,j)
        shear_strain_corner(i,j) = (vp1 - vm1) / coef%half_dlon(j) + (up1 - um1) / coef%half_dlat(j)
      end do
    end do
    call parallel_fill_halo(tension_strain_full, all_halo=.true.)
    call parallel_fill_halo(shear_strain_corner, all_halo=.true.)
    
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_full_lon / 2) then
          i0 = i0 - mesh%num_full_lon / 2
        end if 
        if (j == parallel%full_lat_start_idx) then
          shear_strain_tmp = (shear_strain_corner(i,j)**2 + shear_strain_corner(i-1,j)**2) * 0.5
        elseif (j == parallel%full_lat_end_idx) then
          shear_strain_tmp = (shear_strain_corner(i,j-1)**2 + shear_strain_corner(i-1,j-1)**2) * 0.5
        else 
          shear_strain_tmp = (shear_strain_corner(i-1,j)**2 + shear_strain_corner(i-1,j-1)**2 + shear_strain_corner(i,j)**2 + shear_strain_corner(i,j-1)**2) * 0.25 
        end if 
        length_scale = radius * mesh%full_cos_lat(j) * mesh%dlon     
        deformation_flow = sqrt( tension_strain_full(i,j)**2 + shear_strain_tmp )
        horizontal_viscosity_full(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
      end do 
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        length_scale = radius * mesh%half_cos_lat(j) * mesh%dlon
        tension_strain_tmp = (tension_strain_full(i,j)**2 + tension_strain_full(i+1,j)**2 + tension_strain_full(i,j+1)**2 + tension_strain_full(i+1,j+1)**2) * 0.25
        deformation_flow = sqrt( tension_strain_tmp + shear_strain_corner(i,j)**2 )
        horizontal_viscosity_corner(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
      end do 
    end do 

    call parallel_fill_halo(horizontal_viscosity_full, all_halo=.true.)
    call parallel_fill_halo(horizontal_viscosity_corner, all_halo=.true.)
    ! U
    do j = parallel%full_lat_start_idx_no_pole+2 , parallel%full_lat_end_idx_no_pole-2
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx        
        ud(i,j) = (horizontal_viscosity_full(i+1,j) * (u(i+1,j) - u(i,j)) - horizontal_viscosity_full(i,j) * (u(i,j) -u(i-1,j))) / coef%full_dlon(j)**2 +&
                  ( mesh%half_cos_lat(j) * horizontal_viscosity_corner(i,j) * (u(i,j+1) - u(i,j)) - mesh%half_cos_lat(j-1) * horizontal_viscosity_corner(i,j-1) *(u(i,j) - u(i,j-1)) ) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
      end do
    end do 
    !V 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd(i,j) = ( horizontal_viscosity_corner(i,j) * (v(i+1,j) - v(i,j)) - horizontal_viscosity_corner(i,j) * (v(i+1,j) - v(i,j)) ) / coef%half_dlon(j) +&
                  ( mesh%full_cos_lat(j+1) * horizontal_viscosity_full(i,j+1) * (v(i,j+1) - v(i,j)) - mesh%full_cos_lat(j) * horizontal_viscosity_full(i,j) * (v(i,j) - v(i,j-1)) ) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j)
      end do 
    end do  
        
    call parallel_fill_halo(ud, all_halo=.true.)
    call parallel_fill_halo(vd, all_halo=.true.)

    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%u(i,j) = state%u(i,j) + sign * dt * ud(i,j)
        state%v(i,j) = state%v(i,j) + sign * dt * vd(i,j)
      end do
    end do

    call iap_transform(state)
  end subroutine ordinary_diffusion_nonlinear2

  subroutine shapiro_filter(state)
    !Shapiro, 1975, "Linear Filtering," Mathematics of Computation,
    !                   Volume 29, #132, pp. 1094-1097
    !Table 2: Stencils for different-order filters:
    ! filter order         x(i)        x(i+/-1)    x(i+/-2)    x(i+/-3)    x(i+/-4)
    !       1   1/2^2        2           1
    !       2   1/2^4        10          4           -1
    !       3   1/2^6        44          15          -6          1
    !       4   1/2^8        186         56          -28         8           -1

    type(state_type), intent(inout) :: state
    integer :: i, j, i0
    real :: uxm1, uxm2, uxm3, uxp1, uxp2, uxp3, &
            vxm1, vxm2, vxm3, vxp1, vxp2, vxp3
    real :: uym1, uym2, uym3, uyp1, uyp2, uyp3, &
            vym1, vym2, vym3, vyp1, vyp2, vyp3          
    real, parameter :: lat0 = -90.0

    u(:,:) = state%u(:,:)
    v(:,:) = state%v(:,:)
    call parallel_fill_halo(u,  all_halo=.true.)
    call parallel_fill_halo(v,  all_halo=.true.)

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        uym1 = u(i,j-1)
        uym2 = u(i,j-2)
        uym3 = u(i,j-3)
        uyp1 = u(i,j+1)
        uyp2 = u(i,j+2)
        uyp3 = u(i,j+3)
        i0 = i + mesh%num_half_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_half_lon / 2) then
          i0 = i0 - mesh%num_half_lon / 2
        end if 
        if (j == parallel%full_lat_start_idx_no_pole) then
          uym2 = -u(i0,j)
          uym3 = -u(i0,j+1)
        elseif (j == parallel%full_lat_start_idx_no_pole + 1) then
          uym3 = -u(i0,j-1)
        elseif (j == parallel%full_lat_end_idx_no_pole) then
          uyp2 = -u(i0,j)
          uyp3 = -u(i0,j-1)
        elseif (j == parallel%full_lat_end_idx_no_pole - 1) then
          uyp3 = -u(i0,j+1)
        end if
        uxm1 = u(i-1,j)
        uxm2 = u(i-2,j)
        uxm3 = u(i-3,j)
        uxp1 = u(i+1,j)
        uxp2 = u(i+2,j)
        uxp3 = u(i+3,j)
        if (i == parallel%half_lon_start_idx) then
          uxm1 = u(mesh%num_half_lon,j)
          uxm2 = u(mesh%num_half_lon-1,j)
          uxm3 = u(mesh%num_half_lon-2,j)
        elseif (i == parallel%half_lon_end_idx) then
          uxp1 = u(1,j)
          uxp2 = u(2,j)
          uxp3 = u(3,j)
        endif 
        if (shapiro_filter_order == 2) then 
          u_x(i,j) = ( 10.0 * u(i,j) + 4 * (uxm1 + uxp1) - 1.0 * (uxm2 + uxp2) ) / 16.0
          u_y(i,j) = ( 10.0 * u(i,j) + 4 * (uym1 + uyp1) - 1.0 * (uym2 + uyp2) ) / 16.0
        elseif(shapiro_filter_order == 3) then
          u_x(i,j) = ( 44.0 * u(i,j) + 15.0 * (uxm1 + uxp1) - 6.0 * (uxm2 + uxp2) + 1.0 * (uxm3 + uxp3) ) / 64.0
          u_y(i,j) = ( 44.0 * u(i,j) + 15.0 * (uym1 + uyp1) - 6.0 * (uym2 + uyp2) + 1.0 * (uym3 + uyp3) ) / 64.0
        end if 
      end do 
    end do 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vym1 = v(i,j-1)
        vym2 = v(i,j-2)
        vym3 = v(i,j-3)
        vyp1 = v(i,j+1)
        vyp2 = v(i,j+2)
        vyp3 = v(i,j+3)
        i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
        if (i0 > mesh%num_full_lon / 2) then
          i0 = i0 - mesh%num_full_lon / 2
        end if 
        if (j == parallel%half_lat_start_idx) then
          vym1 = -v(i0,j)
          vym2 = -v(i0,j+1)
          vym3 = -v(i0,j+2)
        elseif (j == parallel%half_lat_start_idx + 1) then    
          vym2 = -v(i0,j-1)
          vym3 = -v(i0,j)
        elseif (j == parallel%half_lat_start_idx + 2) then    
          vym3 = -v(i0,j-2)
        elseif (j == parallel%half_lat_end_idx) then
          vyp1 = -v(i0,j)
          vyp2 = -v(i0,j-1)
          vyp3 = -v(i0,j-2)
        elseif (j == parallel%half_lat_end_idx - 1) then
          vyp2 = -v(i0,j+1)
          vyp3 = -v(i0,j)
        elseif (j == parallel%half_lat_end_idx - 2) then
          vyp3 = -v(i0,j+2)
        end if 
        vxm1 = v(i-1,j)
        vxm2 = v(i-2,j)
        vxm3 = v(i-3,j)
        vxp1 = v(i+1,j)
        vxp2 = v(i+2,j)
        vxp3 = v(i+3,j)
        if (i == parallel%full_lon_start_idx) then
          vxm1 = v(mesh%num_full_lon,j)
          vxm2 = v(mesh%num_full_lon-1,j)
          vxm3 = v(mesh%num_full_lon-2,j)
        elseif (i == parallel%half_lon_end_idx) then
          vxp1 = v(1,j)
          vxp2 = v(2,j)
          vxp3 = v(3,j)
        endif 
        if (shapiro_filter_order == 2) then 
          v_x(i,j) = ( 10.0 * v(i,j) + 4 * (vxm1 + vxp1) - 1.0 * (vxm2 + vxp2) ) / 16.0
          v_y(i,j) = ( 10.0 * v(i,j) + 4 * (vym1 + vyp1) - 1.0 * (vym2 + vyp2) ) / 16.0     
        else if (shapiro_filter_order == 3) then
          v_x(i,j) = ( 44.0 * v(i,j) + 15.0 * (vxm1 + vxp1) - 6.0 * (vxm2 + vxp2) + 1.0 * (vxm3 + vxp3) ) / 64.0
          v_y(i,j) = ( 44.0 * v(i,j) + 15.0 * (vym1 + vyp1) - 6.0 * (vym2 + vyp2) + 1.0 * (vym3 + vyp3) ) / 64.0
        end if 
      end do 
    end do 
    state%u(:,:) = u_y(:,:) !(u_x(:,:) + u_y(:,:)) * 0.5
    state%v(:,:) = v_y(:,:) !(v_x(:,:) + v_y(:,:)) * 0.5
    call parallel_fill_halo(state%u,  all_halo=.true.)
    call parallel_fill_halo(state%v,  all_halo=.true.)

    call iap_transform(state)
  end subroutine shapiro_filter

  subroutine diffusion_final()

    if (allocated(ud))  deallocate(ud)
    if (allocated(vd))  deallocate(vd)
    if (allocated(gdd)) deallocate(gdd)
    if (allocated(u))   deallocate(u)
    if (allocated(v))   deallocate(v)
    if (allocated(gd))  deallocate(gd)

    if (allocated(diffusion_coef_gd)) deallocate(diffusion_coef_gd)
    if (allocated(diffusion_coef_u))  deallocate(diffusion_coef_u)
    if (allocated(diffusion_coef_v))  deallocate(diffusion_coef_v)

    if (allocated(diffusion_coef_gd_x)) deallocate(diffusion_coef_gd_x)
    if (allocated(diffusion_coef_u_x))  deallocate(diffusion_coef_u_x)
    if (allocated(diffusion_coef_v_x))  deallocate(diffusion_coef_v_x)

    if (allocated(diffusion_coef_gd_y)) deallocate(diffusion_coef_gd_y)
    if (allocated(diffusion_coef_u_y))  deallocate(diffusion_coef_u_y)
    if (allocated(diffusion_coef_v_y))  deallocate(diffusion_coef_v_y)
    if (allocated(gd_x)) deallocate(gd_x)
    if (allocated(gd_y)) deallocate(gd_y)
    if (allocated(u_x)) deallocate(u_x)
    if (allocated(u_y)) deallocate(u_y)
    if (allocated(ud_x)) deallocate(ud_x)
    if (allocated(ud_y)) deallocate(ud_y)
    if (allocated(u_tmp)) deallocate(u_tmp)

    if (allocated(v_x)) deallocate(v_x)
    if (allocated(v_y)) deallocate(v_y)
    if (allocated(vd_x)) deallocate(vd_x)
    if (allocated(vd_y)) deallocate(vd_y)  
    if (allocated(v_tmp)) deallocate(v_tmp) 
    ! for nonlinear diffsuion
    if (allocated(tension_strain_full)) deallocate(tension_strain_full)
    if (allocated(shear_strain_corner)) deallocate(shear_strain_corner)
    if (allocated(horizontal_viscosity_full)) deallocate(horizontal_viscosity_full) 
    if (allocated(horizontal_viscosity_corner)) deallocate(horizontal_viscosity_corner) 
 
  end subroutine diffusion_final

end module diffusion_mod
