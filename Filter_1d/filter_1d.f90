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
     call mxm(test_filter_b%fh, lx, f((i-1)*lx+1:i*lx), lx, f_filtered_b((i-1)*lx+1:i*lx), lx)
     call mxm(test_filter_b%fh, lx, f((i-1)*lx+1:i*lx), lx, f_filtered_nb((i-1)*lx+1:i*lx), lx)
  end do

  open(unit = 10, file = "data.csv") 
  do i = 1, lx*nelem
        do j = 1, 4
            if (j == 4) then
                write(unit, '(I0)', advance='no') (/x(i), f(i), f_filtered_b(i), f_filtered_nb(i)/)
            else
                write(unit, '(I0, ",")', advance='no') (/x(i), f(i), f_filtered_b(i), f_filtered_nb(i)/)
            end if
        end do
        write(unit, *) ! New line after each row
  end do
  close(unit)


  if (pe_rank .eq. 0) write(*,*) 'Done'
  
  call neko_finalize

end program filter_1d
