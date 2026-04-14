subroutine output_all(modeflag)

  use blmod
  use parameters
  use physical_constants
  use hdf5
  implicit none

  character(len=1024) :: filename
  integer             :: modeflag
  integer(HID_T)      :: file_id
  integer             :: hdferr
  logical             :: file_exists

  !---------------------------------------------------------------------------
  filename = trim(adjustl(outdir))//"/output.h5"

  call h5open_f(hdferr)

  inquire(file=trim(filename), exist=file_exists)
  if (file_exists) then
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdferr)
  else
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdferr)
  end if

  ! create groups
  call hdf5_ensure_group(file_id, "fields")
  call hdf5_ensure_group(file_id, "scalars")

  ! Always stamp the current time at the root
  call hdf5_write_scalar_r8(file_id, "time", time)

  !---------------------------------------------------------------------------
  if (modeflag == 1) then

    ! ---- mass-coordinate fields -------------------------------------------
    call hdf5_write_1d(file_id, "fields/mass",     mass(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/vel",      vel(1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/rho",      rho(1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/ye",       ye(1:imax),   imax)
    call hdf5_write_1d(file_id, "fields/press",    p(1:imax),    imax)
    call hdf5_write_1d(file_id, "fields/cs2",      cs2(1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/Q",        Q(1:imax),    imax)
    call hdf5_write_1d(file_id, "fields/eps",      eps(1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/temp",     temp(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/lum",      lum(1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/tau",      tau(1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/delta_time", delta_time(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/radius",   r(1:imax),    imax)
    call hdf5_write_1d(file_id, "fields/kappa",    kappa(1:imax),imax)
    call hdf5_write_1d(file_id, "fields/kappa_table", kappa_table(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/logR_op",  logR_op(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/logT",     logT(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/p_rad",    p_rad(1:imax),imax)
    call hdf5_write_1d(file_id, "fields/free_electron_frac", &
                                 free_electron_frac(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/E_shell",  E_shell(1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/time_diff",time_diff(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/time_exp", time_exp(1:imax), imax)
    call hdf5_write_1d(file_id, "fields/photosphere_tracer", &
                                 photosphere_tracer(1:imax), imax)

    ! Ion fractions
    call hdf5_write_1d(file_id, "fields/He_1", &
                        ion_fractions(He_number,1,1:imax), imax)
    call hdf5_write_1d(file_id, "fields/He_2", &
                        ion_fractions(He_number,2,1:imax), imax)
    call hdf5_write_1d(file_id, "fields/He_3", &
                        ion_fractions(He_number,3,1:imax), imax)
    call hdf5_write_1d(file_id, "fields/H_1",  &
                        ion_fractions(H_number,1,1:imax),  imax)
    call hdf5_write_1d(file_id, "fields/H_2",  &
                        ion_fractions(H_number,2,1:imax),  imax)

    ! Optional Ni56 deposit
    if (Ni_mass > 0) then
      call hdf5_write_1d(file_id, "fields/Ni_deposit_function", &
                          Ni_deposit_function(1:imax), imax)
    end if

  !---------------------------------------------------------------------------
  else if (modeflag == 2) then

    ! ---- scalars (one value appended per call) -----------------------------
    call hdf5_append_scalar(file_id, "scalars/T_eff",              T_eff)
    call hdf5_append_scalar(file_id, "scalars/Ni_total_luminosity", Ni_total_luminosity)
    call hdf5_append_scalar(file_id, "scalars/lum_observed",        lum_observed)
    call hdf5_append_scalar(file_id, "scalars/lum_photo",           lum_photo)
    call hdf5_append_scalar(file_id, "scalars/mass_photo",          mass_photo)
    call hdf5_append_scalar(file_id, "scalars/vel_photo",           vel_photo)
    call hdf5_append_scalar(file_id, "scalars/rad_photo",           rad_photo)
    call hdf5_append_scalar(file_id, "scalars/mass_lumshell",       mass_lumshell)
    call hdf5_append_scalar(file_id, "scalars/time",                time)

    call hdf5_append_int(file_id, "scalars/index_photo",       index_photo)
    call hdf5_append_int(file_id, "scalars/opacity_corrupted", opacity_corrupted)
    call hdf5_append_int(file_id, "scalars/index_lumshell",    index_lumshell)

    if (innerBC == "inflow") then
      call hdf5_append_int(file_id, "scalars/inner_boundary", iBC)
    end if

  end if

  call h5fclose_f(file_id, hdferr)
  call h5close_f(hdferr)

end subroutine output_all


!=============================================================================
! Write/overwrite a 1-D real*8 dataset (used for field snapshots)
!=============================================================================
subroutine hdf5_write_1d(file_id, dsetname, data, n)
  use hdf5
  implicit none
  integer(HID_T), intent(in) :: file_id
  character(*),   intent(in) :: dsetname
  integer,        intent(in) :: n
  real(8),        intent(in) :: data(n)

  integer(HID_T)   :: dset_id, dspace_id
  integer(HSIZE_T) :: dims(1)
  integer          :: hdferr
  logical          :: exists

  dims(1) = n

  call h5lexists_f(file_id, dsetname, exists, hdferr)
  if (exists) call h5ldelete_f(file_id, dsetname, hdferr)

  call h5screate_simple_f(1, dims, dspace_id, hdferr)
  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                   dset_id, hdferr)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5sclose_f(dspace_id, hdferr)

end subroutine hdf5_write_1d


!=============================================================================
! Write/overwrite a scalar real*8 dataset
!=============================================================================
subroutine hdf5_write_scalar_r8(file_id, dsetname, val)
  use hdf5
  implicit none
  integer(HID_T), intent(in) :: file_id
  character(*),   intent(in) :: dsetname
  real(8),        intent(in) :: val

  integer(HID_T)   :: dset_id, dspace_id
  integer(HSIZE_T) :: dims(1)
  integer          :: hdferr
  logical          :: exists

  dims(1) = 1

  call h5lexists_f(file_id, dsetname, exists, hdferr)
  if (exists) call h5ldelete_f(file_id, dsetname, hdferr)

  call h5screate_simple_f(1, dims, dspace_id, hdferr)
  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                   dset_id, hdferr)
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, val, dims, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5sclose_f(dspace_id, hdferr)

end subroutine hdf5_write_scalar_r8


!=============================================================================
! Append one real*8 value to a chunked unlimited-dimension dataset
!=============================================================================
subroutine hdf5_append_scalar(file_id, dsetname, val)
  use hdf5
  implicit none
  integer(HID_T), intent(in) :: file_id
  character(*),   intent(in) :: dsetname
  real(8),        intent(in) :: val

  integer(HID_T)   :: dset_id, dspace_id, memspace_id, plist_id
  integer(HSIZE_T) :: cur_dims(1), max_dims(1), new_dims(1), one(1), offset(1)
  integer          :: hdferr
  logical          :: exists

  one(1)      = 1
  max_dims(1) = H5S_UNLIMITED_F

  call h5lexists_f(file_id, dsetname, exists, hdferr)

  if (.not. exists) then
    ! Create chunked extendable dataset with first value
    cur_dims(1) = 1
    call h5screate_simple_f(1, cur_dims, dspace_id, hdferr, max_dims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, hdferr)
    call h5pset_chunk_f(plist_id, 1, one, hdferr)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, &
                     dset_id, hdferr, plist_id)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, val, one, hdferr)
    call h5pclose_f(plist_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
  else
    ! Extend by one and write into the new slot
    call h5dopen_f(file_id, dsetname, dset_id, hdferr)
    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sget_simple_extent_dims_f(dspace_id, cur_dims, max_dims, hdferr)
    call h5sclose_f(dspace_id, hdferr)

    new_dims(1) = cur_dims(1) + 1
    call h5dset_extent_f(dset_id, new_dims, hdferr)

    offset(1) = cur_dims(1)
    call h5screate_simple_f(1, one, memspace_id, hdferr)
    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, one, hdferr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, val, one, hdferr, &
                    memspace_id, dspace_id)
    call h5sclose_f(memspace_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
  end if

  call h5dclose_f(dset_id, hdferr)

end subroutine hdf5_append_scalar


!=============================================================================
! Append one integer value to a chunked unlimited-dimension dataset
!=============================================================================
subroutine hdf5_append_int(file_id, dsetname, val)
  use hdf5
  implicit none
  integer(HID_T), intent(in) :: file_id
  character(*),   intent(in) :: dsetname
  integer,        intent(in) :: val

  integer(HID_T)   :: dset_id, dspace_id, memspace_id, plist_id
  integer(HSIZE_T) :: cur_dims(1), max_dims(1), new_dims(1), one(1), offset(1)
  integer          :: hdferr
  logical          :: exists

  one(1)      = 1
  max_dims(1) = H5S_UNLIMITED_F

  call h5lexists_f(file_id, dsetname, exists, hdferr)

  if (.not. exists) then
    cur_dims(1) = 1
    call h5screate_simple_f(1, cur_dims, dspace_id, hdferr, max_dims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, hdferr)
    call h5pset_chunk_f(plist_id, 1, one, hdferr)
    call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
                     dset_id, hdferr, plist_id)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, val, one, hdferr)
    call h5pclose_f(plist_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
  else
    call h5dopen_f(file_id, dsetname, dset_id, hdferr)
    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sget_simple_extent_dims_f(dspace_id, cur_dims, max_dims, hdferr)
    call h5sclose_f(dspace_id, hdferr)

    new_dims(1) = cur_dims(1) + 1
    call h5dset_extent_f(dset_id, new_dims, hdferr)

    offset(1) = cur_dims(1)
    call h5screate_simple_f(1, one, memspace_id, hdferr)
    call h5dget_space_f(dset_id, dspace_id, hdferr)
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, one, hdferr)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, val, one, hdferr, &
                    memspace_id, dspace_id)
    call h5sclose_f(memspace_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
  end if

  call h5dclose_f(dset_id, hdferr)

end subroutine hdf5_append_int


!******************************************************************************
!output of the variable versus grid point number for a given moment of time
!used to output the initial values of some variables
subroutine output_screenshot(var,filename,imaximum)

  implicit none
  real*8 var(*)
  character(len=100) filename
  integer nt
  integer i
  integer imaximum

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
       form='formatted',position="append")

  do i=1, imaximum
      write(666,"(I5.4, E25.16E3)") i,var(i)
  enddo

  close(666)

end subroutine output_screenshot


!=============================================================================
! Create a group if it doesn't already exist
!=============================================================================
subroutine hdf5_ensure_group(file_id, groupname)
  use hdf5
  implicit none
  integer(HID_T), intent(in) :: file_id
  character(*),   intent(in) :: groupname

  integer(HID_T) :: grp_id
  integer        :: hdferr
  logical        :: exists

  call h5lexists_f(file_id, groupname, exists, hdferr)
  if (.not. exists) then
    call h5gcreate_f(file_id, groupname, grp_id, hdferr)
    call h5gclose_f(grp_id, hdferr)
  end if

end subroutine hdf5_ensure_group
