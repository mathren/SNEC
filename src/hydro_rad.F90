subroutine hydro_rad

  use blmod
  use parameters
  use physical_constants
  implicit none

  integer :: i,j,k,n
  integer :: keytemp,keyerr
  integer :: iBC_new
  real*8 :: dtv

  ! Variables for inversion of Jacobian
  external :: dgbsv
  integer :: info
  integer, parameter :: kl = 1
  integer, parameter :: ku = 1
  integer, parameter :: ldab=2*kl+ku+1
  real*8 :: ab(ldab,imax-1), b(imax-1)
  integer :: ipiv(imax-1)

  real*8 :: delta_max
  integer :: location_max

  real*8 :: Aarray(imax-1), Barray(imax-1), Carray(imax-1), Darray(imax-1)
  real*8 :: p_temp(imax), eps_temp(imax), lum_temp(imax), temp_temp(imax)
  real*8 :: lambda_temp(imax)

  real*8, parameter :: EPSTOL = 1.0d-7
  integer, parameter :: ITMAX = 300


  !------------------------------------------------------------------------------

  ! follow dt convention of MB93:
  dtv = 0.5d0 * (dtime + dtime_p)

  location_max = iBC

  ! copy over data into _p arrays
  rho_p(:) = rho(:)
  vel_p(:) = vel(:)
  r_p(:)   = r(:)
  cr_p(:)  = cr(:)
  eps_p(:) = eps(:)
  p_p(:) = p(:)
  temp_p(:) = temp(:)
  kappa_p(:) = kappa(:)

  !---------------------------- update velocities -------------------------------
  if(do_piston) then
     if (time.ge.piston_tstart .and. time.le.piston_tend) then
        vel(1) = piston_vel
        vel(2) = piston_vel
     else
        vel(1) = 0.0d0
     end if
  else if(do_bomb) then
     vel(1) = 0.0d0
  else if(inject_BE) then
     vel(1) = 0.0d0
  end if

  do i=iBC+1, imax
     vel(i) = vel_p(i) &
          ! gravity
          - dtv * ggrav*mass(i) / r(i)**2 *gravity_switch   &
          ! pressure
          - dtv * 4.0d0*pi*r(i)**2 * (p(i) - p(i-1)) / delta_cmass(i-1)   &
          ! artificial viscosity
          - dtv * 4.0d0*pi * (cr(i)**2 * Q(i) - cr(i-1)**2 * Q(i-1))/delta_cmass(i-1)
  end do

  ! inflow-only (no backreaction) inner boundary
  if (innerBC == "inflow") then
     vel(iBC) = min(0.0d0, vel(iBC+1))
  else
     vel(iBC) = 0.0d0
  end if

  !----------------------- update the radial coordinates-------------------------
  do i=iBC, imax
     r(i) = r_p(i) + dtime * vel(i)
  end do

  ! now we have updated radii and velocities: do we need to move the
  ! inner boundary index (iBC)? Everything that reached r < rBC_initial is
  ! excised. Take the OUTERMOST face that fell inside, so that zone
  ! crossings during fallback are absorbed instead of aborting the run.
  if (innerBC == "inflow") then
     iBC_new = iBC
     do i=iBC, imax-1
        if (r(i) <= rBC_initial) iBC_new = i
     end do

     if (iBC_new.ge.imax-1) then
        write(*,*) "hydro_rad: inner boundary reached the outer zones, stopping"
        stop
     end if

     iBC = iBC_new
     ! the excision radius is fixed: pull the new boundary face back onto it
     r(iBC)   = max(r(iBC), rBC_initial)
     vel(iBC) = min(0.0d0, vel(iBC+1))

     if (iBC>1) then
        ! flatten everything inside the inner boundary to prevent pressure,
        ! temperature and internal energy gradients that could push back
        vel(1:iBC-1) = 0.0d0
        ! keep the dead radii monotonically ordered and below r(iBC): several
        ! routines (optical_depth, nickel_heating/map_find_index, analysis)
        ! assume a monotonic r(1:imax)
        do i=1, iBC-1
           r(i) = 0.95d0*r(iBC)*dble(i-1)/dble(iBC-1)
        end do
        ! zero-gradient extrapolation into the excised region
        eps(1:iBC-1)  = eps(iBC)
        p(1:iBC-1)    = p(iBC)
        temp(1:iBC-1) = temp(iBC)
        rho(1:iBC-1)  = rho(iBC)
     end if
  end if

  do i=iBC+1, imax ! check radial ordering outside inner boundary
     if ((r(i).lt.r(i-1))) then
        write(*,*) 'radius of a gridpoint', i, 'is less than preceding'
        write(*,*) 'boundary at cell', iBC, ' - retrying with a smaller dt'
        scratch_step = .true.
        return
     end if
  end do


  !------------------------- update the zone densities --------------------------

  do i=iBC,imax-1
     rho(i) = delta_mass(i) / (4.0d0*pi * (r(i+1)**3 - r(i)**3)/3.0d0)
  end do
  rho(imax) = 0.0d0 !passive boundary condition

  !------------------------- update zone center radius --------------------------
  do i=iBC,imax-1
     cr(i) = ( ( r(i)**3 + r(i+1)**3 ) / 2.0d0 )**(1.0d0/3.0d0)
  end do
  cr(imax) = r(imax) + (r(imax) - cr(imax-1))
  !passive boundary condition, used in the expression for the velocity update,
  !but multiplied by the artificial viscosity, which is zero at the last point

  ! update the artificial viscosity
  call artificial_viscosity

  !----------- update the temperature, pressure and internal energy -------------

  !calculate heating term due to Ni
  if(time.ge.time_Ni) then
     time_Ni = time_Ni + Ni_period
     call nickel_heating
  end if

  !calculate heating term due to bomb
  if(do_bomb .and. time.ge.bomb_tstart .and. time.le.bomb_tend) then
     call bomb_pattern
  else if (inject_be .and. nt == 0) then
     call inject_progenitor_binding_energy
  else
     bomb_heating(:) = 0.0d0
  end if

  ! no heating inside the excised region
  if (iBC>1) then
     bomb_heating(1:iBC-1) = 0.0d0
     Ni_heating(1:iBC-1)   = 0.0d0
  end if

  !Initial guess for the quantities at the next time step
  p_temp(1:imax) = p(1:imax)
  eps_temp(1:imax) = eps(1:imax)
  temp_temp(1:imax) = temp(1:imax)

  ! number of unknowns: the active zones iBC ... imax-1
  n = imax - iBC

  ! NR iterations start
  do k=1, ITMAX

     ! re-impose the zero-gradient condition on the excised zones so that the
     ! ghost state stays slaved to the current iterate of zone iBC
     if (iBC>1) then
        temp_temp(1:iBC-1) = temp_temp(iBC)
        p_temp(1:iBC-1)    = p_temp(iBC)
        eps_temp(1:iBC-1)  = eps_temp(iBC)
     end if

     keytemp = 1
     call eos(rho(1:imax-1),temp_temp(1:imax-1),ye(1:imax-1), &
          abar(1:imax-1),p_temp(1:imax-1),eps_temp(1:imax-1), &
          cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), entropy(1:imax-1), &
          p_rad(1:imax-1),keyerr,keytemp,eoskey)


     call luminosity(r(:), temp_temp(:), kappa_p(:), &
          lambda_temp(:), inv_kappa(:), lum_temp(:))

     !calculate the coefficients of the equation:
     ! A(i)*\delta T(i+1) + B(i)*\delta T(i) + C(i)*\delta T(i-1) = D(i)
     call matrix_arrays(temp_temp(:), lambda_temp(:), inv_kappa(:), &
          eps_temp(:), p_temp(:), lum_temp(:), &
          Aarray(:), Barray(:), Carray(:), Darray(:))

     ! initialize to zero at every iteration
     ab(:,:)=0.0d0
     !assemble the matrix in the form used by lapack, for the active
     !zones only: local index j = i - iBC + 1, j = 1 ... n
     ab(2,2:n)   = Aarray(iBC:imax-2)
     ab(3,1:n)   = Barray(iBC:imax-1)
     ab(4,1:n-1) = Carray(iBC+1:imax-1)

     b(1:n) = Darray(iBC:imax-1)

     !invert the matrix
     !if the inversion fails, the whole matrix is written in 'failed_matrix.dat'
     info = 0
     call dgbsv(n,kl,ku,1,ab,ldab,ipiv,b,n,info)

     if(info.ne.0) then
        open(unit=666, &
             file=trim(adjustl(trim(adjustl(outdir))//"/failed_matrix.dat")), &
             status="unknown",form='formatted',position="append")
        do j=1,n
           write(666,*) ab(2,j), ab(3,j), ab(4,j), Darray(iBC+j-1)
        end do
        close(666)
        stop "problem in the matrix inversion (see Data/failed_matrix.dat)"
     end if

     !check if the iteration procedure converged
     delta_max = 0.0d0

     do j=1,n
        i = iBC + j - 1
        if(abs(b(j)/temp_temp(i)).gt.delta_max) then
           delta_max = abs(b(j)/temp_temp(i))
           location_max = i
        end if
     end do

     if((delta_max.le.EPSTOL)) goto 101

     !add the increment to the temperature
     do j=1,n
        i = iBC + j - 1
        temp_temp(i) = temp_temp(i) + b(j)
        if(temp_temp(i).lt.0.0d0) then
           goto 100
        end if
     end do

  end do

100 continue

  write(6,*) "EOS problem", delta_max, location_max, iBC
  scratch_step = .true.

101 continue


  eps(iBC:imax-1) = eps_temp(iBC:imax-1)
  p(iBC:imax-1)   = p_temp(iBC:imax-1)
  temp(iBC:imax-1)  = temp_temp(iBC:imax-1)

  if (iBC>1) then
     ! zero-gradient extrapolation into the excised region
     eps(1:iBC-1)  = eps(iBC)
     p(1:iBC-1)    = p(iBC)
     temp(1:iBC-1) = temp(iBC)
  end if

  !passive boundary conditions, do not participate in the evolution
  temp(imax) = 0.0d0
  eps(imax) = 0.0d0

  !active boundary condition, used in the velocity update
  p(imax) = 0.0d0

  ! NOTE: opacity() and luminosity() declare their dummy arguments with an
  ! explicit size imax and loop over 1:imax(-1) internally, skipping the
  ! excised region through the module variable iBC. They must therefore be
  ! called with WHOLE arrays -- passing rho(iBC:imax) shifts every index by
  ! iBC-1 and writes past the end of the output arrays.
  call opacity(rho(:), temp_temp(:), kappa(:), kappa_table(:), dkappadt(:))

  call luminosity(r(:), temp(:), kappa(:), lambda(:), inv_kappa(:), lum(:))

end subroutine hydro_rad
