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
  ! public divergence_diffusion
  public ordinary_diffusion_limiter
  public ordinary_diffusion_direct
  ! public ordinary_diffusion_nonlinear
  ! public ordinary_diffusion_nonlinear2
  public shapiro_filter
  public diffusion_final


  real, allocatable :: u(:,:)
  real, allocatable :: v(:,:)
  real, allocatable :: gd(:,:)
  ! 4th diffusion coefficient
  real, allocatable :: cx_full_lat(:)
  real, allocatable :: cy_full_lat(:)
  real, allocatable :: cx_half_lat(:)
  real, allocatable :: cy_half_lat(:)

  ! for limiter
  real, allocatable :: gx_u(:,:), gy_u(:,:)
  real, allocatable :: dudx(:,:), dudy(:,:)
  real, allocatable :: gx_v(:,:), gy_v(:,:)
  real, allocatable :: dvdx(:,:), dvdy(:,:)
  ! for direct 
  real, allocatable :: lgd_x(:,:), lgd_y(:,:)
  real, allocatable :: lu_x(:,:), lu_y(:,:)
  real, allocatable :: lv_x(:,:), lv_y(:,:)
  ! for nonlinear diffusion
  real, allocatable :: tension_strain_full(:,:), shear_strain_corner(:,:), horizontal_viscosity_full(:,:), horizontal_viscosity_corner(:,:)
  ! for shapiro_filter
  real, allocatable :: u_x(:,:), u_y(:,:), v_x(:,:), v_y(:,:)
contains

  subroutine diffusion_init()

    integer j

    call diffusion_final()

    if (.not. allocated(u))   call parallel_allocate(u,  half_lon=.true.)
    if (.not. allocated(v))   call parallel_allocate(v,  half_lat=.true.)
    if (.not. allocated(gd))  call parallel_allocate(gd)

    allocate(cx_full_lat(parallel%full_lat_start_idx:parallel%full_lat_end_idx))
    allocate(cy_full_lat(parallel%full_lat_start_idx:parallel%full_lat_end_idx))
    allocate(cx_half_lat(parallel%half_lat_start_idx:parallel%half_lat_end_idx))
    allocate(cy_half_lat(parallel%half_lat_start_idx:parallel%half_lat_end_idx))   

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      cx_full_lat(j) = (mesh%full_cos_lat(j) * mesh%dlon * 0.5)**4  !&
                          ! exp(-100.0 * (pi05 - abs(global_mesh%full_lat(j)))**2)
      cy_full_lat(j) = 1.0 / ((1.0 + mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2 * 4.0 / &
                        mesh%dlat**2 + 16.0 / mesh%dlat**4)  
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      cx_half_lat(j) = (mesh%half_cos_lat(j) * mesh%dlon * 0.5)**4 !* &
                      ! exp(-100.0 * (pi05 - abs(global_mesh%half_lat(j)))**2)
      cy_half_lat(j) = 1.0 / ((1.0 + mesh%half_cos_lat(j)**2) / mesh%half_cos_lat(j)**2 * 4.0 / &
                        mesh%dlat**2 + 16.0 / mesh%dlat**4) * 0.1
    end do

    select case(diffusion_method)
    case('limiter') ! Xue(2000)
      if (.not. allocated(gx_u)) call parallel_allocate(gx_u                                  )
      if (.not. allocated(gy_u)) call parallel_allocate(gy_u, half_lon=.true., half_lat=.true.)
      if (.not. allocated(dudx)) call parallel_allocate(dudx                                  )
      if (.not. allocated(dudy)) call parallel_allocate(dudy, half_lon=.true., half_lat=.true.)
  
      if (.not. allocated(gx_v)) call parallel_allocate(gx_v, half_lon=.true., half_lat=.true.)
      if (.not. allocated(gy_v)) call parallel_allocate(gy_v                                  )
      if (.not. allocated(dvdx)) call parallel_allocate(dvdx, half_lon=.true., half_lat=.true.)
      if (.not. allocated(dvdy)) call parallel_allocate(dvdy                                  )
    case('direct') ! ordinary artificial diffusion 
      if (.not. allocated(lgd_x)) call parallel_allocate(lgd_x)
      if (.not. allocated(lgd_y)) call parallel_allocate(lgd_y)
      if (.not. allocated(lu_x))  call parallel_allocate(lu_x, half_lon=.true.)
      if (.not. allocated(lu_y))  call parallel_allocate(lu_y, half_lon=.true.)
      if (.not. allocated(lv_x))  call parallel_allocate(lv_x, half_lat=.true.)
      if (.not. allocated(lv_y))  call parallel_allocate(lv_y, half_lat=.true.)
    case('nonlinear') ! nonlinear diffusion
      if (.not. allocated(tension_strain_full)) call parallel_allocate(tension_strain_full)
      if (.not. allocated(shear_strain_corner)) call parallel_allocate(shear_strain_corner)
      if (.not. allocated(horizontal_viscosity_full)) call parallel_allocate(horizontal_viscosity_full)
      if (.not. allocated(horizontal_viscosity_corner)) call parallel_allocate(horizontal_viscosity_corner)
    case('shapiro_filter') ! shapiro_filter
      if (.not. allocated(u_x)) call parallel_allocate(u_x, half_lon=.true.)
      if (.not. allocated(u_y)) call parallel_allocate(u_y, half_lat=.true.)
      if (.not. allocated(v_x)) call parallel_allocate(v_x, half_lat=.true.)
      if (.not. allocated(v_y)) call parallel_allocate(v_y, half_lat=.true.)
    case default
      write(6, *) '[Error]: Unkown diffusion method: '// trim(diffusion_method) // '!'
    end select 

    call log_notice('diffusion module is initialized.')

  end subroutine diffusion_init

  ! subroutine divergence_diffusion(dt, diag, state)
  !   real, intent(in) :: dt
  !   type(diag_type), intent(in) :: diag
  !   type(state_type), intent(inout) :: state

  !   real, parameter :: vd = 1.0e5
  !   integer i, j

  !   do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       state%u(i,j) = state%u(i,j) + dt * vd * (diag%div(i+1,j) - diag%div(i,j)) / coef%full_dlon(j)
  !     end do
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       state%iap%u(i,j) = state%u(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j))
  !     end do
  !   end do
  !   do j = parallel%half_lon_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       state%v(i,j) = state%v(i,j) + dt * vd * (diag%div(i,j+1) - diag%div(i,j)) / radius / mesh%dlat
  !     end do
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       state%iap%v(i,j) = state%v(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1))
  !     end do
  !   end do

  !   call parallel_fill_halo(state%iap%u(:,:), all_halo=.true.)
  !   call parallel_fill_halo(state%iap%v(:,:), all_halo=.true.)
  !   call parallel_fill_halo(state%u(:,:),     all_halo=.true.)
  !   call parallel_fill_halo(state%v(:,:),     all_halo=.true.)
  ! end subroutine divergence_diffusion


  subroutine ordinary_diffusion_limiter(dt, state)

    real, intent(in) :: dt
    type(state_type), intent(inout) :: state

    integer i, j, ip

    ! u(:,:) = state%iap%u(:,:)
    ! v(:,:) = state%iap%v(:,:)
    ! gd(:,:) = state%iap%gd(:,:)
    u(:,:)  = state%u(:,:)
    v(:,:)  = state%v(:,:)
    gd(:,:) = state%gd(:,:)
    !====================H==============================================
    !1) calculate damping flux and the gradient of forecast variable at interfaces.
    ! H on longitude
  
    !====================U==============================================
    !1) calculate damping flux and the gradient of forecast variable at interfaces.
    ! U on longitude
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gx_u(i,j) = (-u(i-2,j) + 3 * u(i-1,j) - 3 * u(i,j) + u(i+1,j)) / (mesh%full_cos_lat(j) * mesh%dlon)**3
        dudx(i,j) = u(i,j) - u(i-1,j)
      end do
    end do
    ! U on latitude
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ip = i + mesh%num_half_lon / 2
        if (ip > mesh%num_half_lon / 2) then
          ip = ip - mesh%num_half_lon / 2
        end if
        if (j == parallel%half_lat_start_idx) then
          u(i,j-1) = u(ip,j+1)
        else if (j == parallel%half_lat_end_idx) then
          u(i,j+2) = u(ip,j)
        end if
        gy_u(i,j) = -(u(i,j+1) - u(i,j)) / mesh%dlat / mesh%half_cos_lat(j)**2 -    &
                  mesh%half_sin_lat(j) / mesh%half_cos_lat(j) *                     &
                  (u(i,j+2) - u(i,j+1) - u(i,j) + u(i,j-1)) / (2 * mesh%dlat**2) +  &
                  (-u(i,j-1) + 3 * u(i,j) - 3 * u(i,j+1) + u(i,j+2)) / mesh%dlat**3  
        dudy(i,j) = u(i,j+1) - u(i,j)
      end do 
    end do 
    !2) Limit damping flux to avoid upgradient (Xue 2000)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gx_u(i,j) = gx_u(i,j) * max(0.0, sign(1.0, -gx_u(i,j) * dudx(i,j)))
      end do
    end do
    call parallel_fill_halo(gx_u, all_halo=.true.)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        gy_u(i,j) = gy_u(i,j) * max(0.0, sign(1.0, -gy_u(i,j) * dudy(i,j)))
      end do
    end do
    call parallel_fill_halo(gy_u, all_halo=.true.)
    !3) Damp physical variable at last.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) - (                                                  &
                 (gx_u(i+1,j) - gx_u(i,j)) / mesh%full_cos_lat(j) / mesh%dlon * cx_full_lat(j) + &
                 (gy_u(i,j) - gy_u(i,j-1)) / mesh%dlat * cy_full_lat(j)                          &
                          )
      end do
    end do
    !====================V==============================================
    !1) calculate damping flux and the gradient of forecast variable at interfaces.
    ! V on longitude
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        gx_v(i,j) = (-v(i-1,j) + 3 * v(i,j) - 3 * v(i+1,j) + v(i+2,j)) / (mesh%half_cos_lat(j) * mesh%dlon)**3
        dvdx(i,j) = v(i+1,j) - v(i,j)
      end do 
    end do
    ! V on latitude
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip = i + mesh%num_full_lon / 2
        if (ip > mesh%num_full_lon / 2) then
          ip = ip - mesh%num_full_lon / 2
        end if  
        if (j == parallel%full_lat_start_idx_no_pole) then
          v(i,j-2) = v(ip,j-1)
        else if (j == parallel%full_lat_end_idx_no_pole) then
          v(i,j+1) = v(ip,j)
        end if
        gy_v(i,j) = -(v(i,j) - v(i,j-1)) / mesh%full_cos_lat(j)**2 / mesh%dlat -       &
                      mesh%full_sin_lat(j) / mesh%full_cos_lat(j) *                    &
                      (v(i,j+1) - v(i,j) - v(i,j-1) + v(i,j-2)) / (2 * mesh%dlat**2) + &
                     (-v(i,j-2) + 3 * v(i,j-1) - 3 * v(i,j) + v(i,j+1)) / mesh%dlat**3 
        dvdy(i,j) = v(i,j) - v(i,j-1)
      end do 
    end do
    !2) Limit damping flux to avoid upgradient (Xue 2000)
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        gx_v(i,j) = gx_v(i,j) * max(0.0, sign(1.0, -gx_v(i,j) * dvdx(i,j)))
      end do
    end do
    call parallel_fill_halo(gx_v, all_halo=.true.)
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gy_v(i,j) = gy_v(i,j) * max(0.0, sign(1.0, -gy_v(i,j) * dvdy(i,j)))
      end do
    end do
    call parallel_fill_halo(gy_v, all_halo=.true.)
    ! Damp physical variable at last.
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) - (                                           &
          (gx_v(i,j) - gx_v(i-1,j)) / mesh%half_cos_lat(j) / mesh%dlon * cx_half_lat(j) + &
          (gy_v(i,j+1) - gy_v(i,j)) / mesh%dlat * cy_half_lat(j)                          &
                                      )
      end do
    end do

    ! call parallel_fill_halo(state%gd, all_halo=.true.)
    call parallel_fill_halo(state%u,  all_halo=.true.)
    call parallel_fill_halo(state%v,  all_halo=.true.)

    call iap_transform(state)

    !====================================================    
    ! Do FFT filter.
    ! do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
    !   if (filter_full_zonal_tend(j)) then
    !     call filter_array_at_full_lat(j, gdd(:,j))
    !     call filter_array_at_full_lat(j, ud(:,j))
    !   end if
    ! end do
    ! do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
    !   if (filter_half_zonal_tend(j)) then
    !     call filter_array_at_half_lat(j, vd(:,j))
    !   end if
    ! end do  
    !====================================================

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
    ! call parallel_fill_halo(state%gd, all_halo=.true.)
    ! call parallel_fill_halo(state%u,  all_halo=.true.)
    ! call parallel_fill_halo(state%v,  all_halo=.true.)

  end subroutine ordinary_diffusion_limiter

  subroutine ordinary_diffusion_direct(dt, state)
    real, intent(in) :: dt
    type(state_type), intent(inout) :: state
    integer i, j, ip

    u(:,:) = state%u(:,:)
    v(:,:) = state%v(:,:)
    gd(:,:) = state%gd(:,:)
   
    ! H
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        lgd_x(i,j) = (gd(i+2,j) - 4 * gd(i+1,j) + 6 * gd(i,j) - 4 * gd(i-1,j) + gd(i-2,j)) /&
                     (mesh%full_cos_lat(j) * mesh%dlon)**4
        ip = i + mesh%num_full_lon / 2 ! ip is the opposite grid of i 
        if (ip > mesh%num_full_lon / 2) then
          ip = ip - mesh%num_full_lon / 2
        end if 

        if(j == parallel%full_lat_start_idx_no_pole) then ! the next-to-last circle on south
          gd(i,j-2) = gd(ip,j)
        elseif(j == parallel%full_lat_south_pole_idx) then ! the south pole 
          gd(i,j-1) = gd(ip,j+1)
          gd(i,j-2) = gd(ip,j+2)
        elseif(j == parallel%full_lat_end_idx_no_pole) then   ! the next-to-last circle on north
          gd(i,j+2) = gd(ip,j)
        elseif(j == parallel%full_lat_north_pole_idx) then ! the north pole 
          gd(i,j+1) = gd(ip,j-1)
          gd(i,j+2) = gd(ip,j-2)
        end if 
        lgd_y(i,j) = -mesh%full_sin_lat(j) / mesh%full_cos_lat(j)**3 *&
                    (gd(i,j+1) - gd(i,j-1)) / (2 * mesh%dlat) - &
                    (1 + mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2 *&
                    (gd(i,j+1) - 2 * gd(i,j) + gd(i,j-1)) / mesh%dlat**2 -&
                    2 * mesh%full_sin_lat(j) / mesh%full_cos_lat(j) * &
                    (gd(i,j+2) - gd(i,j-2) + 2 * gd(i,j-1) -2 * gd(i,j+1)) / (2 * mesh%dlat**3) +&
                    (gd(i,j+2) - 4 * gd(i,j+1) + 6 * gd(i,j) - 4 * gd(i,j-1) + gd(i,j-2)) / mesh%dlat**4
      end do
    end do
    !==============================================   
    ! U
    do j = parallel%full_lat_start_idx_no_pole , parallel%full_lat_end_idx_no_pole 
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        lu_x(i,j) = (u(i+2,j) - 4 * u(i+1,j) + 6 * u(i,j) - 4 * u(i-1,j) + u(i-2,j)) / (mesh%full_cos_lat(j) * mesh%dlon)**4
        ip = i + mesh%num_half_lon / 2 ! i0 is the opposite grid of i 
        if (ip > mesh%num_half_lon / 2) then
          ip = ip - mesh%num_half_lon / 2
        end if 
        if (j == parallel%full_lat_start_idx_no_pole) then
          u(i,j-2) = u(ip,j)
        elseif(j == parallel%full_lat_end_idx_no_pole) then
          u(i,j+2) = u(ip,j)
        end if 
        lu_y(i,j) = -mesh%full_sin_lat(j) / mesh%full_cos_lat(j)**3 * &
                    (u(i,j+1) - u(i,j-1)) / (2.0 * mesh%dlat) - &
                    (1.0 + mesh%full_cos_lat(j)**2) / mesh%full_cos_lat(j)**2* &
                    (u(i,j+1) - 2.0 * u(i,j) + u(i,j-1)) / mesh%dlat**2 - &
                    2 * mesh%full_sin_lat(j) / mesh%full_cos_lat(j) * &
                    (u(i,j+2) - u(i,j-2) + 2.0 * u(i,j-1) - 2.0 * u(i,j+1)) / (2.0 * mesh%dlat**3) + &
                    (u(i,j+2) - 4.0 * u(i,j+1) + 6.0 * u(i,j) - 4.0 * u(i,j-1) + u(i,j-2)) / mesh%dlat**4
      end do
    end do   
    !==============================================
    !V 
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        lv_x(i,j) = (v(i+2,j) - 4 * v(i+1,j) + 6 * v(i,j) - 4 * v(i-1,j) + v(i-2,j)) / (mesh%half_cos_lat(j) * mesh%dlon)**4
        ip = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
        if (ip > mesh%num_full_lon / 2) then
          ip = ip - mesh%num_full_lon / 2
        end if 
        if (j == parallel%half_lat_start_idx_no_pole) then      !! the next-to-last circle on south 4th
          v(i,j-2) = v(ip,j-1)
        elseif(j == parallel%half_lat_start_idx) then
          v(i,j-1) = v(ip,j)
          v(i,j-2) = v(ip,j+1)
        elseif(j == parallel%half_lat_end_idx_no_pole ) then
          v(i,j+2) = v(ip,j+1)
        elseif(j == parallel%half_lat_end_idx) then
          v(i,j+1) = v(ip,j)
          v(i,j+2) = v(ip,j-1)
        end if 
        lv_y(i,j) = -mesh%half_sin_lat(j) / mesh%half_cos_lat(j)**3 * &
                    (v(i,j+1) - v(i,j-1))/(2*mesh%dlat) - &
                    (1+mesh%half_cos_lat(j)**2) / mesh%half_cos_lat(j)**2 * &
                    (v(i,j+1) - 2* v(i,j) + v(i,j-1)) / mesh%dlat**2 - &
                    2 * mesh%half_sin_lat(j) / mesh%half_cos_lat(j) * &
                    (v(i,j+2) - v(i,j-2) + 2 * v(i,j-1) - 2 * v(i,j+1)) / (2 * mesh%dlat**3) + &
                    (v(i,j+2) - 4 * v(i,j+1) + 6 * v(i,j) - 4 * v(i,j-1) + v(i,j-2)) / mesh%dlat**4
      end do 
    end do 

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
       ! do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
         ! state%gd(i,j) = state%gd(i,j) - (cx_full_lat(j) * lgd_x(i,j) + cy_full_lat(j) * lgd_y(i,j))
       ! end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ! state%u(i,j) = state%u(i,j) - (cx_full_lat(j) * lu_x(i,j) + cy_full_lat(j) * lu_y(i,j))
        state%u(i,j) = state%u(i,j) - cy_full_lat(j) * lu_y(i,j)
      end do 
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ! state%v(i,j) = state%v(i,j) - cx_half_lat(j) * lv_x(i,j) + cy_half_lat(j) * lv_y(i,j))
        state%v(i,j) = state%v(i,j) - cy_half_lat(j) * lv_y(i,j)
      end do
    end do

    call iap_transform(state)

  end subroutine ordinary_diffusion_direct

  ! subroutine ordinary_diffusion_nonlinear(dt, state)
  !   ! Smagorinsky(1963) Nonlinear second order horizontal diffusion
  !   real, intent(in) :: dt
  !   type(state_type), intent(inout) :: state
  !   real :: length_scale, shear_strain_tmp, tension_strain_tmp, deformation_flow
  !   real :: sp, np
  !   integer :: i, j, order, i0, sign
  !   real :: um1, up1, vm1, vp1
  
  !   u(:,:) = state%u(:,:)
  !   v(:,:) = state%v(:,:)

  !   do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
  !       if (i0 > mesh%num_full_lon / 2) then
  !         i0 = i0 - mesh%num_full_lon / 2
  !       end if 
  !       um1 = state%u(i-1,j)
  !       up1 = state%u(i,j)
  !       vm1 = state%v(i,j-1) * mesh%half_cos_lat(j-1)
  !       vp1 = state%v(i,j) * mesh%half_cos_lat(j)
  !       if (j == parallel%full_lat_start_idx) then
  !         ! vm1 = state%v(i0,j) * mesh%half_cos_lat(j)
  !         vm1 = -state%v(i0,j) * cos(mesh%half_lat(j) - mesh%dlat)
  !       elseif (j == parallel%full_lat_end_idx) then
  !         ! vp1 = state%v(i0,j-1) * mesh%half_cos_lat(j-1)
  !         vp1 = - state%v(i0,j-1) * cos(mesh%half_lat(j-1) + mesh%dlat)
  !       end if 
  !       tension_strain_full(i,j) = (up1 - um1) / coef%full_dlon(j) - (vp1 - vm1) / coef%full_dlat(j)
  !     end do 
  !   end do 

  !   do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       um1 = state%u(i,j) * mesh%full_cos_lat(j)
  !       up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
  !       vm1 = state%v(i,j)
  !       vp1 = state%v(i+1,j)
  !       shear_strain_corner(i,j) = (vp1 - vm1) / coef%half_dlon(j) + (up1 - um1) / coef%half_dlat(j)
  !     end do
  !   end do
  !   call parallel_fill_halo(tension_strain_full, all_halo=.true.)
  !   call parallel_fill_halo(shear_strain_corner, all_halo=.true.)

  !   do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       if (j == parallel%full_lat_start_idx) then
  !         shear_strain_tmp = (shear_strain_corner(i,j)**2 + shear_strain_corner(i-1,j)**2) * 0.5
  !       elseif (j == parallel%full_lat_end_idx) then
  !         shear_strain_tmp = (shear_strain_corner(i,j-1)**2 + shear_strain_corner(i-1,j-1)**2) * 0.5
  !       else 
  !         shear_strain_tmp = (shear_strain_corner(i-1,j)**2 + shear_strain_corner(i-1,j-1)**2 + shear_strain_corner(i,j)**2 + shear_strain_corner(i,j-1)**2) * 0.25 
  !       end if 
  !       length_scale = radius * mesh%full_cos_lat(j) * mesh%dlon     
  !       deformation_flow = sqrt( tension_strain_full(i,j)**2 + shear_strain_tmp )
  !       horizontal_viscosity_full(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
  !     end do 
  !   end do 

  !   do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       length_scale = radius * mesh%half_cos_lat(j) * mesh%dlon
  !       tension_strain_tmp = (tension_strain_full(i,j)**2 + tension_strain_full(i+1,j)**2 + tension_strain_full(i,j+1)**2 + tension_strain_full(i+1,j+1)**2) * 0.25
  !       deformation_flow = sqrt( tension_strain_tmp + shear_strain_corner(i,j)**2 )
  !       horizontal_viscosity_corner(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
  !     end do 
  !   end do 
  !   call parallel_fill_halo(horizontal_viscosity_full,   all_halo=.true.)
  !   call parallel_fill_halo(horizontal_viscosity_corner, all_halo=.True.)
  !   ! U
  !   do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       ud(i,j) = (horizontal_viscosity_full(i+1,j) * tension_strain_full(i+1,j) - horizontal_viscosity_full(i,j) * tension_strain_full(i,j)) / coef%full_dlon(j) + &
  !                 (horizontal_viscosity_corner(i,j) * shear_strain_corner(i,j) - horizontal_viscosity_corner(i,j-1) * shear_strain_corner(i,j-1)) / radius / mesh%dlat
  !     end do
  !   end do   
  !   !V 
  !   do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx        
  !       vd(i,j) = (horizontal_viscosity_corner(i,j) * shear_strain_corner(i,j) - horizontal_viscosity_corner(i-1,j) * shear_strain_corner(i-1,j) ) / coef%half_dlon(j) -&
  !                 (horizontal_viscosity_full(i,j+1) * tension_strain_full(i,j+1) - horizontal_viscosity_full(i,j) * tension_strain_full(i,j)) / radius / mesh%dlat
  !     end do 
  !   end do   

  !   call parallel_fill_halo(ud,  all_halo=.true.)
  !   call parallel_fill_halo(vd,  all_halo=.true.)

  !   sign = (-1)**(diffusion_order / 2 + 1)

  !   do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       state%u(i,j) = state%u(i,j) + sign * dt * ud(i,j) 
  !     end do
  !   end do

  !   do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       state%v(i,j) = state%v(i,j) + sign * dt * vd(i,j) 
  !     end do 
  !   end do 

  !   call parallel_fill_halo(state%u,  all_halo=.true.)
  !   call parallel_fill_halo(state%v,  all_halo=.true.)

  !   call iap_transform(state)
  ! end subroutine ordinary_diffusion_nonlinear
  
  ! subroutine ordinary_diffusion_nonlinear2(dt, state)
  !   !  David L. Williamson. 1978.
  !   ! The Relative Importance of Resolution, Accuracy and Diffusion in 
  !   ! Short-Range Forecasts with the NCAR Global Circulation Model
  !   real, intent(in) :: dt
  !   type(state_type), intent(inout) :: state
  !   real :: length_scale, shear_strain_tmp, tension_strain_tmp, deformation_flow
  !   real :: sp, np
  !   integer :: i, j, order, i0, sign
  !   real :: um1, up1, vm1, vp1

  !   u(:,:) = state%u(:,:)
  !   v(:,:) = state%v(:,:)

  !   do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
  !       if (i0 > mesh%num_full_lon / 2) then
  !         i0 = i0 - mesh%num_full_lon / 2
  !       end if 
  !       um1 = state%u(i-1,j)
  !       up1 = state%u(i,j)
  !       vm1 = state%v(i,j-1) * mesh%half_cos_lat(j-1)
  !       vp1 = state%v(i,j) * mesh%half_cos_lat(j)
  !       if (j == parallel%full_lat_start_idx) then
  !         vm1 = state%v(i0,j) * mesh%half_cos_lat(j)
  !       elseif (j == parallel%full_lat_end_idx) then
  !         vp1 = state%v(i0,j-1) * mesh%half_cos_lat(j-1)
  !       end if 
  !       tension_strain_full(i,j) = (up1 - um1) / coef%full_dlon(j) - (vp1 - vm1) / coef%full_dlat(j)
  !     end do 
  !   end do 

  !   do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       um1 = state%u(i,j) * mesh%full_cos_lat(j)
  !       up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
  !       vm1 = state%v(i,j)
  !       vp1 = state%v(i+1,j)
  !       shear_strain_corner(i,j) = (vp1 - vm1) / coef%half_dlon(j) + (up1 - um1) / coef%half_dlat(j)
  !     end do
  !   end do
  !   call parallel_fill_halo(tension_strain_full, all_halo=.true.)
  !   call parallel_fill_halo(shear_strain_corner, all_halo=.true.)
    
  !   do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       i0 = i + mesh%num_full_lon / 2 ! i0 is the opposite grid of i 
  !       if (i0 > mesh%num_full_lon / 2) then
  !         i0 = i0 - mesh%num_full_lon / 2
  !       end if 
  !       if (j == parallel%full_lat_start_idx) then
  !         shear_strain_tmp = (shear_strain_corner(i,j)**2 + shear_strain_corner(i-1,j)**2) * 0.5
  !       elseif (j == parallel%full_lat_end_idx) then
  !         shear_strain_tmp = (shear_strain_corner(i,j-1)**2 + shear_strain_corner(i-1,j-1)**2) * 0.5
  !       else 
  !         shear_strain_tmp = (shear_strain_corner(i-1,j)**2 + shear_strain_corner(i-1,j-1)**2 + shear_strain_corner(i,j)**2 + shear_strain_corner(i,j-1)**2) * 0.25 
  !       end if 
  !       length_scale = radius * mesh%full_cos_lat(j) * mesh%dlon     
  !       deformation_flow = sqrt( tension_strain_full(i,j)**2 + shear_strain_tmp )
  !       horizontal_viscosity_full(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
  !     end do 
  !   end do 

  !   do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
  !       length_scale = radius * mesh%half_cos_lat(j) * mesh%dlon
  !       tension_strain_tmp = (tension_strain_full(i,j)**2 + tension_strain_full(i+1,j)**2 + tension_strain_full(i,j+1)**2 + tension_strain_full(i+1,j+1)**2) * 0.25
  !       deformation_flow = sqrt( tension_strain_tmp + shear_strain_corner(i,j)**2 )
  !       horizontal_viscosity_corner(i,j) = 2 * diffusion_coef**2 * length_scale**2 * deformation_flow
  !     end do 
  !   end do 

  !   call parallel_fill_halo(horizontal_viscosity_full, all_halo=.true.)
  !   call parallel_fill_halo(horizontal_viscosity_corner, all_halo=.true.)
  !   ! U
  !   do j = parallel%full_lat_start_idx_no_pole+2 , parallel%full_lat_end_idx_no_pole-2
  !     do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx        
  !       ud(i,j) = (horizontal_viscosity_full(i+1,j) * (u(i+1,j) - u(i,j)) - horizontal_viscosity_full(i,j) * (u(i,j) -u(i-1,j))) / coef%full_dlon(j)**2 +&
  !                 ( mesh%half_cos_lat(j) * horizontal_viscosity_corner(i,j) * (u(i,j+1) - u(i,j)) - mesh%half_cos_lat(j-1) * horizontal_viscosity_corner(i,j-1) *(u(i,j) - u(i,j-1)) ) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
  !     end do
  !   end do 
  !   !V 
  !   do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       vd(i,j) = ( horizontal_viscosity_corner(i,j) * (v(i+1,j) - v(i,j)) - horizontal_viscosity_corner(i,j) * (v(i+1,j) - v(i,j)) ) / coef%half_dlon(j) +&
  !                 ( mesh%full_cos_lat(j+1) * horizontal_viscosity_full(i,j+1) * (v(i,j+1) - v(i,j)) - mesh%full_cos_lat(j) * horizontal_viscosity_full(i,j) * (v(i,j) - v(i,j-1)) ) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j)
  !     end do 
  !   end do  
        
  !   call parallel_fill_halo(ud, all_halo=.true.)
  !   call parallel_fill_halo(vd, all_halo=.true.)

  !   sign = (-1)**(diffusion_order / 2 + 1)

  !   do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
  !     do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
  !       state%u(i,j) = state%u(i,j) + sign * dt * ud(i,j)
  !       state%v(i,j) = state%v(i,j) + sign * dt * vd(i,j)
  !     end do
  !   end do

  !   call iap_transform(state)
  ! end subroutine ordinary_diffusion_nonlinear2

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
          uym2 = u(i0,j)
          uym3 = u(i0,j+1)
        elseif (j == parallel%full_lat_start_idx_no_pole + 1) then
          uym3 = u(i0,j-1)
        elseif (j == parallel%full_lat_end_idx_no_pole) then
          uyp2 = u(i0,j)
          uyp3 = u(i0,j-1)
        elseif (j == parallel%full_lat_end_idx_no_pole - 1) then
          uyp3 = u(i0,j+1)
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
          vym1 = v(i0,j)
          vym2 = v(i0,j+1)
          vym3 = v(i0,j+2)
        elseif (j == parallel%half_lat_start_idx + 1) then    
          vym2 = v(i0,j-1)
          vym3 = v(i0,j)
        elseif (j == parallel%half_lat_start_idx + 2) then    
          vym3 = v(i0,j-2)
        elseif (j == parallel%half_lat_end_idx) then
          vyp1 = v(i0,j)
          vyp2 = v(i0,j-1)
          vyp3 = v(i0,j-2)
        elseif (j == parallel%half_lat_end_idx - 1) then
          vyp2 = v(i0,j+1)
          vyp3 = v(i0,j)
        elseif (j == parallel%half_lat_end_idx - 2) then
          vyp3 = v(i0,j+2)
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

    if (allocated(u))   deallocate(u)
    if (allocated(v))   deallocate(v)
    if (allocated(gd))  deallocate(gd)

    if (allocated(cx_full_lat)) deallocate(cx_full_lat)
    if (allocated(cy_full_lat)) deallocate(cy_full_lat)
    if (allocated(cx_half_lat)) deallocate(cx_half_lat)
    if (allocated(cy_half_lat)) deallocate(cy_half_lat)


    if (allocated(gx_u)) deallocate(gx_u)
    if (allocated(gy_u)) deallocate(gy_u)
    if (allocated(gx_v)) deallocate(gx_v)
    if (allocated(gy_v)) deallocate(gy_v)

    if (allocated(dudx)) deallocate(dudx)
    if (allocated(dudy)) deallocate(dudy)
    if (allocated(dvdx)) deallocate(dvdx)
    if (allocated(dvdy)) deallocate(dvdy)


    ! for nonlinear diffsuion
    if (allocated(tension_strain_full)) deallocate(tension_strain_full)
    if (allocated(shear_strain_corner)) deallocate(shear_strain_corner)
    if (allocated(horizontal_viscosity_full)) deallocate(horizontal_viscosity_full) 
    if (allocated(horizontal_viscosity_corner)) deallocate(horizontal_viscosity_corner) 
 
  end subroutine diffusion_final

end module diffusion_mod
