! Martin Karp 13/3-2023
module user
  use neko
  use fast3d
  implicit none

  ! Global user variables
  real(kind=rp), allocatable :: f(:,:,:), f_2d_x(:,:), f_2d_z(:,:) ! work array to store field in a block
  real(kind=rp), allocatable :: fb(:,:,:) ! work array to including the ending boundary value
  real(kind=rp), allocatable :: ident(:,:) ! work array for the identity matrix
  integer :: kx_cutoff = 2
  integer :: kz_cutoff = 2
  integer :: nelvx = 8
  integer :: nelvy = 8
  integer :: nelvz = 8
  integer :: lx, ly, lz ! GLL points in an element
  integer :: nx, ny, nz ! no-overlapping points in the domain
  type(matrix_t) :: wt, wtt, wt_inv, wt_invt ! weighting matrix for mapping between different collocation poits
  ! operator to link field on Fourier collocation space and state space.
  complex, allocatable :: phi_x(:,:), phi_x_inv(:,:), phi_z(:,:), phi_z_inv(:,:)
  real(kind=rp), allocatable :: filter_amp_x(:,:), filter_amp_z(:,:), filter_x(:,:), filter_z(:,:)

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_ic => user_ic
    u%user_mesh_setup => user_mesh_scale
    u%user_check => user_field_filtering
    u%user_init_modules => user_initialize
    u%user_finalize_modules => user_finalize
  end subroutine user_setup

  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, nvert
    real(kind=rp) :: y, beta

    beta = 1.4_rp

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = pi*msh%points(i)%x(1)

       y = msh%points(i)%x(2) + 1 ! store it from (-1,1) to (0,2)
       ! distort and map into (-1,1)
       msh%points(i)%x(2) = tanh(beta* (y - 1.0_rp))/tanh(beta)

       msh%points(i)%x(3) = pi*msh%points(i)%x(3) 
    end do
    
  end subroutine user_mesh_scale


  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i
    real(kind=rp) :: uvw(3)

    do i = 1, u%dof%size()
       uvw = channel_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
  end subroutine user_ic

  ! Kind of brute force with rather large initial disturbances 
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, uy, uz, Re_tau, yp, Re_b
    real(kind=rp) :: eps1, alpha1, beta1, kx1, kz1
    real(kind=rp) :: eps2, alpha2, beta2, kx2, kz2
    real(kind=rp) :: C, k

      Re_tau = 180.0
      C      = 5.17
      k      = 0.41
      Re_b   = 2860

      yp = (1-y)*Re_tau
      if (y.lt.0) yp = (1+y)*Re_tau
      
      ! Reichardt function
      ux  = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
            (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
      ux  = ux * Re_tau/Re_b

      eps1 = 3e-2
      kx1  = 5
      kz1  = 2

      alpha1 = kx1 * 2*PI/(pi*2)
      beta1  = kz1 * 2*PI/(pi)

      eps2 = 5e-3
      kx2  = 5
      kz2  = 5

      alpha2 = kx2 * 2*PI/(pi*2)
      beta2  = kz2 * 2*PI/(pi)

      ! ! add perturbation to trigger turbulence 
      uvw(1)  = ux  + eps1*beta1  * sin(alpha1*x)*cos(beta1*z) + eps2*beta2  * sin(alpha2*x)*cos(beta2*z)
      uvw(2)  =       eps1        * sin(alpha1*x)*sin(beta1*z) + eps2        * sin(alpha2*x)*sin(beta2*z)
      uvw(3)  =     - eps1*alpha1 * cos(alpha1*x)*sin(beta1*z) - eps2*alpha2 * cos(alpha2*x)*sin(beta2*z)

      ! uvw(1) = (z-pi/2)**2    ! allocate(f_local(nelvz_local*(u%xh%lz-1),nelvy*(u%xh%ly-1),nelvx*(u%xh%lx-1)))
      
  end function channel_ic

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    lx = u%xh%lx
    ly = u%xh%ly
    lz = u%xh%lz
    nx = nelvx * (lx - 1)
    ny = nelvy * ly
    nz = nelvz * (lz - 1)

    allocate(f(nx,ny,nz))
    allocate(fb(nx+1,ny,nz+1))
    allocate(f_2d_z(nx*ny,nz))
    allocate(f_2d_x(nx,ny*nz))

    call map_wt_init(u)

    allocate(phi_x(nx,nx))
    allocate(phi_x_inv(nx,nx))
    allocate(filter_amp_x(nx,nx))
    allocate(filter_x(nx,nx))

    allocate(phi_z(nz,nz))
    allocate(phi_z_inv(nz,nz))
    allocate(filter_amp_z(nz,nz))
    allocate(filter_z(nz,nz))

    call Fourier_init(phi_x, phi_x_inv, filter_amp_x, filter_x, kx_cutoff, nx)
    call Fourier_init(phi_z, phi_z_inv, filter_amp_z, filter_z, kz_cutoff, nz)
    
    call user_field_filtering(t, 0, u, v, w, p, coef, params)

  end subroutine user_initialize

  ! User-defined routine called at the end of every time step
  subroutine user_field_filtering(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    
    integer :: i

    ! Interpolate to equidistant points on homogeneous directions
    call map_collocation_pts(u, "z", wt%x, wtt%x, ident, lx)
    ! Filtering over homogeneous directions
    call field_global_filtering_z(u, coef)
   !  call field_global_filtering_z(v, coef)
   !  call field_global_filtering_z(w, coef)
   !  call field_global_filtering_z(p, coef)

    ! Interpolate back to GLL points on homogeneous directions
    call map_collocation_pts(u, "z", wt_inv%x, wt_invt%x, ident, lx)
  
  end subroutine user_field_filtering

  subroutine map_wt_init(u)
    type(field_t), intent(in) :: u
    real(kind=rp) :: x_equid
    integer :: i

    call wt%init(lx, lx)
    call wtt%init(lx, lx)
    call wt_inv%init(lx, lx)
    call wt_invt%init(lx, lx)
    allocate(ident(lx,lx))

    ident = 0.0_rp
    do i = 1, lx
      x_equid = -1.0_rp + (i-1) * 2.0_rp/(lx-1)
      call fd_weights_full(x_equid, u%xh%zg(:,1), lx-1, 0, wtt%x(:,i))
      wt%x(i,:) = wtt%x(:,i)
      ident(i,i) = 1.0_rp
    end do

    wt_inv%x = wt%x
    wt_invt%x = wtt%x
    call wt_inv%inverse()
    call wt_invt%inverse()

  end subroutine map_wt_init

  subroutine Fourier_init(phi, phi_inv, filter_amp, filter, k_cutoff, n)
    complex, intent(inout) :: phi(n,n), phi_inv(n,n)
    real(kind=rp), intent(inout) :: filter_amp(n,n), filter(n,n)
    integer :: k_cutoff
    integer, intent(in) :: n
    integer :: i, j, k(n)
    real(kind=rp) :: ksi(n)
    complex :: w1(n,n) ! work array
    
    if (mod(n,2).ne.0) then
       call neko_error("Fourier global filtering for odd number grid points has not been implemented")
    end if

    ! wave number array kx and kz
    do i = 1, int(n/2)
       k(i) = i - 1
    end do
    do i = int(n/2)+1, n
       k(i) = -1 + i - n
    end do
    ! local coordinates between 0 and 2pi
    do i = 1, n
       ksi(i) = 2*pi*(i-1)/n
    end do
    
    !!! Step1: Construct phi and phi_inv such that phi f_hat(:,i,j) = f(:,i,j)
    do i = 1, n ! row
       do j = 1, n ! column
          phi_inv(i,j) = exp(cmplx(0.0_rp, -2.0_rp*pi*(i-1)*(j-1)/n))/n
          phi(i,j) = exp(cmplx(0.0_rp, k(j)*ksi(i)))
       end do
    end do
    !!! Step2: Perform global modal cutoff
    filter_amp = 0.0_rp
    filter_amp(1,1) = 1.0_rp
    do i = 1, k_cutoff
       filter_amp(i+1,i+1) = 1.0_rp
       filter_amp(nx-i+1,nx-i+1) = 1.0_rp
    end do
    
    !!! Step3: Formulate the filtering operator
    w1 = cmplx(filter_amp, 0.0_rp)
    w1 = matmul(w1,phi_inv)
    w1 = matmul(phi,w1)
    filter = real(w1)

  end subroutine Fourier_init

  ! 1d mapping between GLL collocation points and Fourier collocation points
  subroutine map_collocation_pts(u, hom_dir, wt, wtt, ident, lx)
    type(field_t), intent(inout) :: u
    character(len=1), intent(in) :: hom_dir
    integer :: lx
    real(kind=rp) :: wt(lx,lx), wtt(lx,lx), ident(lx,lx)

    ! interpolate the field
    if (trim(hom_dir) .eq. 'x') then
       call tnsr1_3d(u%x, lx, lx, wt, ident, ident, u%msh%nelv)
    else if (trim(hom_dir) .eq. 'z') then
       call tnsr1_3d(u%x, lx, lx, ident, ident, wtt, u%msh%nelv)
    else
       call neko_error("homogeneous direction is wrongly set")
    end if

  end subroutine map_collocation_pts
  
  ! cross-element filtering designated for structured mesh
  subroutine field_global_filtering_z(u, coef)
    type(field_t), intent(inout) :: u
    type(coef_t), intent(inout) :: coef
    integer :: i, j, k, m, jk(2)
    integer :: i_start, i_step, i_step_pre, ierr
    real(kind=rp), allocatable :: f_2d(:,:)
    real(kind=rp) :: ref(28)    
    
    !!! Step0: Formulate a global array of the field
    call rzero(f,size(f))
    call assemble_global_field(f, nx, ny, nz, u, nelvx, nelvy, nelvz, lx, ly, lz)
    
    !!! Step3: Perform operator parallelly
    !! Step 3.1: form a 2d field array for tensor product
    if (mod(nx*ny,pe_size) .ne. 0) then
       if (pe_rank .ne. pe_size - 1) then
          i_step = int(nx*ny/pe_size)
          i_step_pre = i_step
       else
          i_step = mod(nx*ny,pe_size) + int(nx*ny/pe_size)
          i_step_pre = int(nx*ny/pe_size)
       end if
    else
       i_step = int(nx*ny/pe_size)
       i_step_pre = i_step
    end if
    i_start = pe_rank*i_step_pre + 1
    do i = i_start, i_start + i_step - 1
       jk = index_element_1d_to_2d(i, nx)     
       j = jk(1)
       k = jk(2)
       f_2d_z(i,:) = f(j,k,:)
    end do
    call rzero(f,size(f)) ! clear work array f
    
    !! Step 3.2: filtering
    f_2d_z(i_start: i_start + i_step - 1, :) = &
       transpose(matmul(filter_z, &
       transpose(f_2d_z(i_start: i_start + i_step - 1, :))))
    
    !! Step 3.2: reshape the array
    do i = i_start, i_start + i_step - 1
       jk = index_element_1d_to_2d(i, nx)
       j = jk(1)
       k = jk(2)
       f(j,k,:) = f_2d_z(i,:)
    end do

    !! Step 3.3: collect information from all ranks
    call MPI_Allreduce(MPI_IN_PLACE, f, size(f), &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

    
    !!! Step4: disassemble the global array into u%x
    fb(1:nx,1:ny,1:nz) = f
    !!! resume the periodicity boundary condition
    fb(:,:,nz+1) = fb(:,:,1)
    fb(nx+1,:,:) = fb(1,:,:)
    call disassemble_global_field(fb, nx+1, ny, nz+1, u, nelvx, nelvy, nelvz, lx, ly, lz)
    
  end subroutine field_global_filtering_z

  function index_element_1d_to_3d(i_1d, n1, n2) result(ijk)
    integer, intent(in) :: i_1d, n1, n2
    integer :: remainder, i, j, k
    integer :: ijk(3)

    k = int((i_1d - 1) / (n1 * n2)) + 1
    remainder = i_1d - (k - 1) * n1 * n2
    j= int((remainder - 1) / n1) + 1
    i = remainder - (j - 1) * n1

    ijk(1) = i
    ijk(2) = j
    ijk(3) = k
  end function index_element_1d_to_3d

  function index_element_1d_to_2d(i_1d, n1) result(ij)
    integer, intent(in) :: i_1d, n1
    integer :: remainder, i, j
    integer :: ij(2)

    j = int((i_1d - 1) / n1) + 1
    i = i_1d - (j - 1) * n1

    ij(1) = i
    ij(2) = j
  end function index_element_1d_to_2d

  subroutine assemble_global_field(field_glb, nx, ny, nz, u, nelvx, nelvy, nelvz, lx, ly, lz)
    integer, intent(in) :: nx, ny, nz
    real(kind=rp), intent(inout) :: field_glb(nx,ny,nz)
    type(field_t), intent(in) :: u
    integer, intent(in) :: nelvx, nelvy, nelvz, lx, ly, lz

    integer :: index_global_element_xyz(3), glb_index_element
    integer :: index_global_element_x, index_global_element_y, index_global_element_z
    integer :: index_global_pts_x, index_global_pts_y, index_global_pts_z
    integer :: i, ierr

    do i = 1, u%msh%nelv
       glb_index_element = i + u%msh%offset_el
       ! break down the 1d global index into 3d for box mesh
       index_global_element_xyz = index_element_1d_to_3d(glb_index_element, nelvx, nelvy)
       index_global_element_x = index_global_element_xyz(1)
       index_global_element_y = index_global_element_xyz(2)
       index_global_element_z = index_global_element_xyz(3)
       ! form the global field locally
       index_global_pts_x = (index_global_element_x-1) * (lx-1) + 1
       index_global_pts_y = (index_global_element_y-1) * ly + 1
       index_global_pts_z = (index_global_element_z-1) * (lz-1) + 1

       field_glb(index_global_pts_x:index_global_pts_x + (lx-2), &
         index_global_pts_y:index_global_pts_y + (ly-1), &
         index_global_pts_z:index_global_pts_z + (lz-2)) = u%x(1:lx-1, :, 1:lz-1, i)
    end do
    ! Assemble the global field
    call MPI_Allreduce(MPI_IN_PLACE, f, size(f), &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

  end subroutine assemble_global_field

  subroutine disassemble_global_field(fb, nxb, nyb, nzb, u, nelvx, nelvy, nelvz, lx, ly, lz)
    integer, intent(in) :: nxb, nyb, nzb
    real(kind=rp), intent(inout) :: fb(nxb,nyb,nzb)
    type(field_t), intent(inout) :: u
    integer, intent(in) :: nelvx, nelvy, nelvz, lx, ly, lz

    integer :: index_global_element_xyz(3), glb_index_element
    integer :: index_global_element_x, index_global_element_y, index_global_element_z
    integer :: index_global_pts_x, index_global_pts_y, index_global_pts_z
    integer :: i, ierr

    do i = 1, u%msh%nelv
       glb_index_element = i + u%msh%offset_el
       ! break down the 1d global index into 3d for box mesh
       index_global_element_xyz = index_element_1d_to_3d(glb_index_element, nelvx, nelvy)
       index_global_element_x = index_global_element_xyz(1)
       index_global_element_y = index_global_element_xyz(2)
       index_global_element_z = index_global_element_xyz(3)
       ! form the global field locally
       index_global_pts_x = (index_global_element_x-1) * (lx-1) + 1
       index_global_pts_y = (index_global_element_y-1) * ly + 1
       index_global_pts_z = (index_global_element_z-1) * (lz-1) + 1

       u%x(1:lx, 1:ly, 1:lz, i) = &
          fb(index_global_pts_x:index_global_pts_x + (lx-1), &
                    index_global_pts_y:index_global_pts_y + (ly-1), &
                    index_global_pts_z:index_global_pts_z + (lz-1))
    end do

  end subroutine disassemble_global_field

  ! User-defined finalization routine called at the end of the simulation
  subroutine user_finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    deallocate(f)
    deallocate(fb)
    deallocate(ident)
    call wt%free()
    call wtt%free()
    call wt_inv%free()
    call wt_invt%free()
    deallocate(phi_z)
    deallocate(phi_z_inv)
    deallocate(phi_x)
    deallocate(phi_x_inv)
    deallocate(filter_amp_x)
    deallocate(filter_amp_z)
    deallocate(filter_x)
    deallocate(filter_z)
    deallocate(f_2d_x)
    deallocate(f_2d_z)

  end subroutine user_finalize
end module user
