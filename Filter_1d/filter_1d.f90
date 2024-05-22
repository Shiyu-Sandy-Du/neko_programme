program filter_1d
  use neko
  use elementwise_filter, only : elementwise_filter_t
  implicit none

  type(file_t) :: output_file
  character(len=64) :: path
  type(vector_t) :: vec_out 
  type(space_t) :: Xh
  real(kind=rp), allocatable :: x(:), f(:), f_filtered_b(:), f_filtered_nb(:)
  real(kind=rp) :: xmin, xmax, x_elem_min, len_elem
  integer :: lx, nelem, i, j, unit
  type(elementwise_filter_t) :: test_filter_b, test_filter_nb

  call neko_init  
  
  ! output_file = file_init("output.csv")

  lx = 8
  nelem = 1
  xmin = -1.0_rp
  xmax = 1.0_rp
  path = "/scratch/shiyud/nekoexamples/test/"

  call Xh%init(GLL, lx, lx, lx)

  ! Initialize the filter
  call test_filter_b%init(lx, "Boyd")
  test_filter_b%trnsfr(test_filter_b%nx-0) = 0.0   
  test_filter_b%trnsfr(test_filter_b%nx-1) = 0.05
  test_filter_b%trnsfr(test_filter_b%nx-2) = 0.50
  test_filter_b%trnsfr(test_filter_b%nx-3) = 0.95  
  call test_filter_b%build_1d()

  call test_filter_nb%init(lx, "nonBoyd")
  test_filter_nb%trnsfr(test_filter_nb%nx-0) = 0.0   
  test_filter_nb%trnsfr(test_filter_nb%nx-1) = 0.05
  test_filter_nb%trnsfr(test_filter_nb%nx-2) = 0.50
  test_filter_nb%trnsfr(test_filter_nb%nx-3) = 0.95  
  call test_filter_nb%build_1d()

  ! function setup, elements are assumed to be of equal size
  allocate(f(lx*nelem))
  allocate(x(lx*nelem))
  len_elem = (xmax - xmin)/nelem
  do i = 1, nelem
     x_elem_min = xmin + (i-1)*len_elem
     do j = 1, lx
        x((i-1)*lx + j) =  x_elem_min + (Xh%zg(j,1)+1.0_rp)*len_elem/2.0_rp
     end do
  end do
  
  do i = 1, lx*nelem
     f(i) = - 2.0_rp*sin(pi*x(i)) &
            - sin(3.0_rp*pi*x(i))
  end do

  ! filtering the function
  allocate(f_filtered_b(lx*nelem))
  allocate(f_filtered_nb(lx*nelem))
  
  do i = 1, nelem
     call mxm(test_filter_b%fh, lx, f((i-1)*lx+1:i*lx), lx, f_filtered_b((i-1)*lx+1:i*lx), 1)
     call mxm(test_filter_b%fh, lx, f((i-1)*lx+1:i*lx), lx, f_filtered_nb((i-1)*lx+1:i*lx), 1)
  end do
  
  ! output
  call vec_out%init(lx*nelem)

  output_file = file_init(trim(path)//trim("x.csv"))
  vec_out%x = x
  call output_file%write(vec_out)
  
  output_file = file_init(trim(path)//trim("f.csv"))
  vec_out%x = f
  call output_file%write(vec_out)
  
  output_file = file_init(trim(path)//trim("f_filtered_b.csv"))
  vec_out%x = f_filtered_b
  call output_file%write(vec_out)
  
  output_file = file_init(trim(path)//trim("f_filtered_nb.csv"))
  vec_out%x = f_filtered_nb
  call output_file%write(vec_out)
  
  call file_free(output_file)

  if (pe_rank .eq. 0) write(*,*) 'Done'
  
  call neko_finalize

end program filter_1d
