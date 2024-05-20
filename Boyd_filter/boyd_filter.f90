program boyd_filter
  use neko
  use elementwise_filter, only : elementwise_filter_t
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, field_fname, hom_dir, output_fname
  type(file_t) :: field_file, output_file
  type(fld_file_data_t) :: field_data
  type(space_t) :: Xh
  type(vector_ptr_t), allocatable :: fields(:)
  integer :: argc, i, lx, j, file_precision
  logical :: dp_precision
  type(elementwise_filter_t) :: test_filter
  real(kind=rp) :: t1, t2

  argc = command_argument_count()

  if ((argc .lt. 3) .or. (argc .gt. 3)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./boyd_filter field.fld outfield.fld precision'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, fmt='(A)') field_fname
  call get_command_argument(2, inputchar)
  read(inputchar, fmt='(A)') output_fname
  call get_command_argument(3, inputchar)
  read(inputchar, *) dp_precision

  if (dp_precision) then
     file_precision = dp
  else
     file_precision = sp
  end if

  field_file = file_t(trim(field_fname),precision=file_precision)

  call field_data%init()

  if (pe_rank .eq. 0) write(*,*) 'Reading file:', 1
  call field_file%read(field_data)

  lx = field_data%lx

  call Xh%init(GLL, field_data%lx, field_data%ly, field_data%lz)

  ! Initialize the filter
  call test_filter%init(lx, "Boyd")
  test_filter%trnsfr(test_filter%nx-0) = 0.0   
  test_filter%trnsfr(test_filter%nx-1) = 0.0
  test_filter%trnsfr(test_filter%nx-2) = 0.0
  test_filter%trnsfr(test_filter%nx-3) = 0.0   
  test_filter%trnsfr(test_filter%nx-4) = 0.05
  test_filter%trnsfr(test_filter%nx-5) = 0.50   
  test_filter%trnsfr(test_filter%nx-6) = 0.95
  test_filter%trnsfr(test_filter%nx-7) = 1.00
  call test_filter%build_1d()

  ! filtering the field at t=0
  allocate(fields(field_data%size()))
  call field_data%get_list(fields,field_data%size())
  if (pe_rank .eq. 0) write(*,*) 'Filtering:'
  do i = 1, field_data%size()
     call test_filter%filter_3d(fields(i)%ptr%x, fields(i)%ptr%x, field_data%nelv)
  end do
  if (pe_rank .eq. 0) write(*,*) 'Writing file:'
  ! output at t=0
  output_file = file_t(trim(output_fname),precision=file_precision)
  call output_file%write(field_data, field_data%time)

  ! filtering field for t>0
  do i = 1, field_data%meta_nsamples-1

     t1 = MPI_Wtime()
     if (pe_rank .eq. 0) write(*,*) " "
     if (pe_rank .eq. 0) write(*,*) 'Reading file:', i+1
     call field_file%read(field_data)
     t2 = MPI_Wtime()
     if (pe_rank .eq. 0) write(*,*) 'Reading Time:', t2 - t1

     call field_data%get_list(fields,field_data%size())
     t1 = MPI_Wtime()
     if (pe_rank .eq. 0) write(*,*) 'get_list Time:', t1 - t2

     do j = 1, field_data%size()
        call test_filter%filter_3d(fields(j)%ptr%x, fields(j)%ptr%x, field_data%nelv)
     end do
     t2 = MPI_Wtime()
     if (pe_rank .eq. 0) write(*,*) 'Filtering Time:', t2 - t1

     ! output for t>0
     call output_file%write(field_data, field_data%time)
     t1 = MPI_Wtime()
     if (pe_rank .eq. 0) write(*,*) 'Writing file Time:', t1 - t2

  end do


  if (pe_rank .eq. 0) write(*,*) 'Done'

  call neko_finalize

end program boyd_filter
