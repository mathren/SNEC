subroutine opacity(rho_x,temp_x,kappa_x,kappa_table_x,dkappadt_x)

  use blmod, only: logT, logR_op, opacity_floor
  use parameters
  use physical_constants
  implicit none

  !input:
  real*8 rho_x(imax)
  real*8 temp_x(imax)

  !output:
  real*8 kappa_x(imax)
  real*8 kappa_table_x(imax)
  real*8 dkappadt_x(imax)

  !local:
  integer i
  real*8 t6
  real*8 R_op
  real*8 log_kappa

  integer no_tabular_value

  !---------------------------------------------------------------------------
  ! the envelope metallicity is assumed to be the nominal metallicity
  ! of the OPAL Type II tables; see blmod.F90

  do i=1, imax - 1

     !definitions of R_op and t6 are given in the codata files (opacity tables)
     t6 = temp_x(i)*1.0d-6

     logT(i) = log10(temp_x(i))

     R_op = rho_x(i)/(t6**3)

     logR_op(i) = log10(R_op)

     call bicubic_interpolation(i,logR_op(i),logT(i),log_kappa,no_tabular_value)

     if (no_tabular_value.eq.0) then

         kappa_x(i) = 10.0d0**log_kappa

     else

        if(i.eq.1) then
           kappa_x(i) = opacity_floor(i)
        else
           kappa_x(i) = kappa_x(i-1)
        end if

     end if

     !derivative is not used in the current version of the code
     dkappadt_x(i) = 0.0d0

  end do

  kappa_table_x(1:imax-1) = kappa_x(1:imax-1)

  do i=1, imax-1
    if( kappa_x(i) .lt. opacity_floor(i)) then
        kappa_x(i)   =  opacity_floor(i)
    end if
  end do

  ! Since the temperature temp(imax) is evaluated in the evolution, but
  ! set by the boundary condition, it doesn't make sense to find
  ! kappa(imax) from the table, so we assume:
  kappa_x(imax) = kappa_x(imax-1)
  kappa_table_x(imax) = kappa_table_x(imax-1)
  dkappadt_x(imax) = dkappadt_x(imax-1)

end subroutine opacity

!################ COMPOSING OPACITY TABLES FROM OPAL TABLES ###################
!####### this is done once at the beginning (called in problem.F90) ###########

subroutine compose_opacity_tables_OPAL

  use blmod, only: logR_array, logT_array, opacity_tables, rpoints, tpoints, &
      op_rows, op_cols, envelope_metallicity, xxo, xxc, H_number, C_number, &
      O_number, comp
  use parameters
  use physical_constants
  implicit none

  integer :: i,j,k
  real*8 :: R_op, t6

  !for OPAL interpolation routine
  real*4 opact,dopact,dopacr,dopactd
  common/e/ opact,dopact,dopacr,dopactd

  ! --- Cache ---
  character(len=256) :: cache_dir
  integer :: mkdir_stat
  character(len=256) :: cache_file
  logical :: cache_exists
  integer :: io_unit = 99, io_stat
  ! --- Cache fingerprint ---
  integer :: cached_imax
  real*8  :: cached_metallicity
  real*8, allocatable :: cached_comp(:,:)

  !---------------------------------------------------------------------------

  cache_dir = trim(adjustl(outdir))//'/cache'
  cache_file = trim(adjustl(cache_dir))//'/opacity_cache.bin'

  inquire(file=cache_file, exist=cache_exists)

  if (cache_exists) then
    write(*,*) "Loading opacity tables from cache..."
    open(unit=io_unit, file=cache_file, form='unformatted', status='old', &
         action='read', iostat=io_stat)
    if (io_stat /= 0) goto 10  ! fall through to recompute if read fails

    ! --- read and validate fingerprint before loading tables ---
    allocate(cached_comp(imax, size(comp, 2)))
    read(io_unit, iostat=io_stat) cached_imax, cached_metallicity, cached_comp
    if (io_stat /= 0 .or. &
        cached_imax /= imax .or. &
        abs(cached_metallicity - envelope_metallicity) > 1.0d-12 .or. &
        any(abs(cached_comp - comp) > 1.0d-12)) then
      write(*,*) "Cache fingerprint mismatch or unreadable — recomputing."
      close(io_unit)
      deallocate(cached_comp)
      goto 10
    end if
    deallocate(cached_comp)

    read(io_unit, iostat=io_stat) tpoints, rpoints
    if (io_stat /= 0) then
       close(io_unit)
       goto 10
    end if

    allocate(logT_array(tpoints))
    allocate(logR_array(rpoints))
    allocate(opacity_tables(imax, tpoints, rpoints))
    allocate(op_rows(tpoints, 2))
    allocate(op_cols(rpoints, 2))

    read(io_unit, iostat=io_stat) logT_array, logR_array, &
                                   op_rows, op_cols, opacity_tables
    close(io_unit)

    if (io_stat /= 0) then
      write(*,*) "Cache read error — recomputing..."
      deallocate(logT_array, logR_array, opacity_tables, op_rows, op_cols)
      goto 10
    end if

    write(*,*) "Opacity tables loaded from cache."
    return
  end if

10 continue  ! <-- recompute label

!------------------------------------------------------------------------------

  write(*,*) "composing opacity tables from OPAL (may take a few minutes)..."

  tpoints = 134     !maximum size of the table in logT direction
  rpoints = 19      !maximum size of the table in logR direction

  allocate(opacity_tables(imax,tpoints,rpoints)) !log kappa(gridpoint:logT:logR)

!------------- Assigning values to the log T and log R arrays -----------------

  allocate(logT_array(tpoints)) !values of log T
  allocate(logR_array(rpoints)) !values of log R

  logT_array(1) = 2.7d0
  do i=2, 5
    logT_array(i) = logT_array(i-1)+0.05d0
  end do
  logT_array(6) = logT_array(5)+0.01d0
  do i=7, 11
    logT_array(i) = logT_array(i-1)+0.02d0
  end do
  do i=12, 60
    logT_array(i) = logT_array(i-1)+0.01d0
  end do
  do i=61, 110
    logT_array(i) = logT_array(i-1)+0.05d0
  end do
  do i=111, 131
    logT_array(i) = logT_array(i-1)+0.1d0
  end do
  do i=132,134
    logT_array(i) = logT_array(i-1)+0.2d0
  end do

  do j=1, rpoints
    logR_array(j) = -8.0d0 + (j-1)*0.5d0
  end do

!------------------------ Structure of the tables -----------------------------

  ! opacity tables may miss values in the corner regions, so
  ! op_rows and op_cols designate the region of tables where the
  ! meaningful data is present

  !op_rows(i,1) is the index of the first meaningful value of the i-th row
  !op_rows(i,2) is the index of the last  meaningful value of the i-th row

  !op_cols(j,1) is the index of the first meaningful value of the j-th column
  !op_cols(j,2) is the index of the last  meaningful value of the j-th column

  allocate(op_rows(tpoints,2))
  allocate(op_cols(rpoints,2))

  do i=1, tpoints
    op_rows(i,1) = 1
  end do
  do i=1, tpoints-13
    op_rows(i,2) = 19
  end do
  op_rows(tpoints-12,2) = 18
  op_rows(tpoints-11,2) = 17
  op_rows(tpoints-10,2) = 17
  op_rows(tpoints- 9,2) = 16
  op_rows(tpoints- 8,2) = 16
  op_rows(tpoints- 7,2) = 16
  op_rows(tpoints- 6,2) = 16
  op_rows(tpoints- 5,2) = 15
  op_rows(tpoints- 4,2) = 15
  op_rows(tpoints- 3,2) = 15
  op_rows(tpoints- 2,2) = 15
  op_rows(tpoints- 1,2) = 15
  op_rows(tpoints,2)    = 14

  do j=1, rpoints
    op_cols(j,1) = 1
  end do
  do j=1, 14
    op_cols(j,2) = 134
  end do
  op_cols(15,2) = 133
  op_cols(16,2) = 128
  op_cols(17,2) = 124
  op_cols(18,2) = 122
  op_cols(19,2) = 121

  do k=1, imax

    !definitions of xxc and xxo are given in opal_opacity.F90
    xxc(k) = MAX( 0.0d0, comp(k,C_number) - envelope_metallicity*C_frac_sol )

    xxo(k) = MAX( 0.0d0, comp(k,O_number) - envelope_metallicity*O_frac_sol )

    do i = 1, tpoints
        do j = 1, rpoints
            if (i.lt.op_cols(j,1) .or. i.gt.op_cols(j,2) &
                .or. j.lt.op_rows(i,1) .or. j.gt.op_rows(i,2)) then

                !out of the meaningful region
                opacity_tables(k,i,j) = 0.0d0

            else

                !if the values of log T and log R are at the very boundaries
                !of meaningful regions, shift them inwards on 1.0d-5, so that
                !OPAL interpolator does not return an error
                if (i.eq.op_cols(j,1)) then
                    t6 = (10.0d0**(logT_array(i)+1.0d-5))*1.0d-6
                else if (i.eq.op_cols(j,2)) then
                    t6 = (10.0d0**(logT_array(i)-1.0d-5))*1.0d-6
                else
                    t6 = (10.0d0**logT_array(i))*1.0d-6
                end if
                if (j.eq.op_rows(i,1)) then
                    R_op = 10.0d0**(logR_array(j)+1.0d-5)
                else if (j.eq.op_rows(i,2)) then
                    R_op = 10.0d0**(logR_array(j)-1.0d-5)
                else
                    R_op = 10.0d0**logR_array(j)
                end if

                call opac( real(envelope_metallicity), real(comp(k,H_number)), &
                    real(xxc(k)), real(xxo(k)), real(t6), real(R_op) )

                opacity_tables(k,i,j) = dble(opact)

            end if
        end do
    end do

  end do

  write(*,*) "finished composing opacity tables"

   ! --- Write cache ---
  write(*,*) "Saving opacity tables to cache..."
  call execute_command_line('mkdir -p '//trim(cache_dir), exitstat=mkdir_stat)
  if (mkdir_stat /= 0) write(*,*) "Warning: could not create cache directory: ", trim(cache_dir)

  open(unit=io_unit, file=cache_file, form='unformatted', status='replace', &
       action='write', iostat=io_stat)
  if (io_stat == 0) then
     ! --- write fingerprint first ---
    write(io_unit) imax, envelope_metallicity, comp
    write(io_unit) tpoints, rpoints
    write(io_unit) logT_array, logR_array, op_rows, op_cols, opacity_tables
    close(io_unit)
    write(*,*) "Cache saved to: ", cache_file
  else
    write(*,*) "Warning: could not write opacity cache."
  end if

end subroutine compose_opacity_tables_OPAL
