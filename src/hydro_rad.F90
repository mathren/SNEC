subroutine hydro_rad

  use blmod
  use parameters
  use physical_constants
  implicit none

  integer :: i,k,j
  integer :: keytemp,keyerr
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

  do i=iBC+1,imax
     vel(i) = vel_p(i) &
          ! gravity
          - dtv * ggrav*mass(i) / r(i)**2 *gravity_switch   &
          ! pressure
          - dtv * 4.0d0*pi*r(i)**2 * (p(i) - p(i-1)) / delta_cmass(i-1)   &
          ! artificial viscosity
          - dtv * 4.0d0*pi * (cr(i)**2 * Q(i) - cr(i-1)**2 * Q(i-1))/delta_cmass(i-1)
  end do


  !----------------------- update the radial coordinates-------------------------

  do i=iBC+1,imax
     r(i) = r_p(i) + dtime * vel(i)
     if((i>iBC) .and. &
          (r(i).lt.r(i-1)) .and. &
          innerBC /= "inflow") then
        write(*,*) 'radius of a gridpoint', i, 'is less than preceding'
        write(*,*) 'boundary at cell', iBC
        stop
     end if
  end do

  ! now we have updated radii and velocities: do we need to move the
  ! inner boundary index (iBC)? cut everything reaching r smaller thana
  ! the initial smaller radius, a fixed boundary
  if (innerBC == "inflow") then
     i = imax
     do while (i > iBC)
        if (r(i) <= max(rBC_initial, r(iBC))) then
           ! ----------------------!
           ! update inner boundary !
           ! update iBC            !
           ! ----------------------!
           iBC = i+1
           if (iBC > imax) stop "inner boundary == outer cell"
           r(iBC) = max(rBC_initial, r(iBC))
           exit
        end if
        i = i - 1 ! loop inward
     end do

     if (iBC>1) then
        ! wipe velocities below
        vel(1:iBC-1) = 0.0d0
        ! fix radii below
        ! linearly spaced grid between 0 and 95% of r(iBC)
        do i=1,iBC-1,1
           r(i) = (r(iBC)*0.95d0)*real(i-1)/real(iBC-1)
        end do
     end if

     ! hack inner boundary to be inflow no backreaction
     vel(iBC) = min(0.0d0, vel(iBC+1)) !vel(iBC+1) !
  end if

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

  !Initial guess for the quantities at the next time step
  p_temp(1:imax) = p(1:imax)
  eps_temp(1:imax) = eps(1:imax)
  temp_temp(1:imax) = temp(1:imax)

  do k=1, ITMAX

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


     !assemble the matrix in the form used by lapack
     ab(2,2:imax-1) = Aarray(1:imax-2)
     ab(3,1:imax-1) = Barray(1:imax-1)
     ab(4,1:imax-2) = Carray(2:imax-1)

     b(1:imax-1) = Darray(1:imax-1)

     !invert the matrix
     !if the inversion fails, the whole matrix is written in 'failed_matrix.dat'
     info = 0
     call dgbsv(imax-1,kl,ku,1,ab,ldab,ipiv,b,imax-1,info)

     if(info.ne.0) then
        open(unit=666, &
             file=trim(adjustl(trim(adjustl(outdir))//"/failed_matrix.dat")), &
             status="unknown",form='formatted',position="append")
        do i=1,imax-1
           write(666,*) ab(2,i), ab(3,i), ab(4,i), Darray(i)
        end do
        close(666)
        stop "problem in the matrix inversion (see Data/failed_matrix.dat)"
     end if

     !check if the iteration procedure converged
     delta_max = 0.0d0

     j = iBC

     if (iBC>1) then
        ! flatten everything inside inner boundary
        ! prevent pressure, temperature and internal
        ! energy gradients which could cause backreaction
        eps(1:iBC) = eps(iBC+1)
        p(1:iBC) = p(iBC+1)
        temp(1:iBC) = temp(iBC+1)
        temp_temp(1:iBC) = temp_temp(iBC+1)
        if (vel(iBC)<0) then
           ! do not check change in T at inner boundary
           j = iBC+1
        end if
     end if

     do i=j,imax-1 ! loop from j set above
        if(abs(b(i)/temp_temp(i)).gt.delta_max) then
           delta_max = abs(b(i)/temp_temp(i))
           location_max = i
        end if
     end do

     if((delta_max.le.EPSTOL)) goto 101

     !add the increment to the temperature
     do i=j, imax-1
        temp_temp(i) = temp_temp(i) + b(i)
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

  !passive boundary conditions, do not participate in the evolution
  temp(imax) = 0.0d0
  eps(imax) = 0.0d0

  !active boundary condition, used in the velocity update
  p(imax) = 0.0d0

  call opacity(rho(iBC:imax),temp_temp(iBC:imax),kappa(iBC:imax),kappa_table(iBC:imax),dkappadt(iBC:imax))

  call luminosity(r(iBC:imax),temp(iBC:imax),kappa(iBC:imax),lambda(iBC:imax),inv_kappa(iBC:imax),lum(iBC:imax))

end subroutine hydro_rad
