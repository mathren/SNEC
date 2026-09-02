subroutine matrix_arrays(temp_x, lambda_x, inv_kappa_x, eps_x, p_x, lum_x, &
    Aarray_x, Barray_x, Carray_x, Darray_x)

  use blmod, only: delta_mass, delta_cmass, dtime, theta, &
        eps_p, p_p, lum, dedt, dpdt, rho, rho_p, r, Qterm, &
        bomb_heating, Ni_heating, iBC
  use parameters
  use physical_constants
  implicit none

!input
  real*8 :: temp_x(imax)
  real*8 :: lambda_x(imax)
  real*8 :: inv_kappa_x(imax)
  real*8 :: eps_x(imax), p_x(imax), lum_x(imax)

!output:
  real*8 :: Aarray_x(imax-1)
  real*8 :: Barray_x(imax-1)
  real*8 :: Carray_x(imax-1)
  real*8 :: Darray_x(imax-1)
  real*8 :: const
  real*8 :: inv_delta_mass
  real*8 :: inv_delta_cmass
  real*8 :: inv_delta_cmass_iM1
  real*8 :: r_iP1_4
  real*8 :: r_4
  real*8 :: lambda_inv_kappa
  real*8 :: temp_3
  real*8 :: inv_rho
  real*8 :: inv_rho_p

!local:
  integer :: i

!------------------------------------------------------------------------------
  !** common constant used in all terms
  const = -theta * dtime * 64.0d0 * pi * pi * a_rad * clite / 3.0d0

  if (iBC.ge.imax-1) stop "matrix_arrays: inner boundary reached outermost zone"

  !** zones inside the moving inner boundary are not part of the system.
  !** Give them a trivial identity row (delta T = 0) so that the banded
  !** solve stays regular even if it is called on the full 1:imax-1 range.
  Aarray_x(:) = 0.0d0
  Barray_x(:) = 1.0d0
  Carray_x(:) = 0.0d0
  Darray_x(:) = 0.0d0

  do i=iBC, imax-2
     inv_delta_mass = 1.0d0 / delta_mass(i)
     inv_delta_cmass = 1.0d0 / delta_cmass(i)
     r_iP1_4 = r(i+1)**4
     r_4 = r(i)**4
     lambda_inv_kappa = lambda_x(i) * inv_kappa_x(i)
     temp_3 = temp_x(i)**3
     inv_rho = 1.0d0 / rho(i)
     inv_rho_p = 1.0d0 / rho_p(i)

  !************* A array **************** coefficients of \delta T(i+1)
     Aarray_x(i) = const*inv_delta_mass * &
          r_iP1_4 * &
          (lambda_x(i+1)*inv_kappa_x(i+1) * &
          temp_x(i+1)**3*inv_delta_cmass )

     if (i.gt.iBC) then
  !--------------------------- interior zones ---------------------------------
        inv_delta_cmass_iM1 = 1.0d0 / delta_cmass(i-1)

  !************* B array **************** coefficients of \delta T(i)
        Barray_x(i) = dedt(i) + 0.5d0*dpdt(i)* &
             (inv_rho-inv_rho_p) &
             - const*inv_delta_mass * &
             r_iP1_4  * &
             (lambda_x(i+1)*inv_kappa_x(i+1)* &
             temp_3*inv_delta_cmass ) &
             - const*inv_delta_mass * &
             r_4 * &
             (lambda_inv_kappa* &
             temp_3*inv_delta_cmass_iM1 )

  !************* C array **************** coefficients of \delta T(i-1)
        Carray_x(i) = const*inv_delta_mass * r_4 * &
             (lambda_inv_kappa*temp_x(i-1)**3*inv_delta_cmass_iM1 )

  !************* D array **************** free terms
        Darray_x(i) = &
             eps_p(i) - 0.5d0*p_p(i)*(inv_rho-inv_rho_p) &
             - (1-theta)*(lum(i+1)-lum(i))*dtime*inv_delta_mass &
             + Qterm(i) &
             - ( eps_x(i) + 0.5d0*p_x(i)*(inv_rho-inv_rho_p) &
             + theta*dtime*inv_delta_mass * ( lum_x(i+1)-lum_x(i) ) ) &
             + bomb_heating(i)*dtime + Ni_heating(i)*dtime

     else
  !------------------- innermost active zone: i == iBC -------------------------
  !** zero-gradient (zero-flux) inner boundary condition: L(iBC) = 0 and
  !** dL(iBC)/dT = 0, so the inner-face terms are dropped from B and D, and
  !** C is identically zero. This is the generalisation of the i=1 rows of
  !** the original code to a moving inner boundary. Without it the Jacobian
  !** retains an inner-face flux derivative (evaluated with the placeholder
  !** lambda=1, inv_kappa=1 that luminosity() writes into 1:iBC) that the
  !** residual D does not contain, and row iBC stays coupled to the frozen
  !** zone iBC-1 through C.
        Barray_x(i) = dedt(i) + 0.5d0*dpdt(i)*(inv_rho-inv_rho_p) &
             - const*inv_delta_mass * &
             r_iP1_4 * &
             (lambda_x(i+1)*inv_kappa_x(i+1)* &
             temp_3*inv_delta_cmass )

        Carray_x(i) = 0.0d0

        Darray_x(i) = &
             eps_p(i) - 0.5d0*p_p(i)*(inv_rho-inv_rho_p) &
             - (1-theta)*lum(i+1)*dtime*inv_delta_mass &
             + Qterm(i) &
             - ( eps_x(i) + 0.5d0*p_x(i)*(inv_rho-inv_rho_p) &
             + theta*dtime*inv_delta_mass * lum_x(i+1) ) &
             + bomb_heating(i)*dtime + Ni_heating(i)*dtime
     end if
  end do

  !********** Outermost zone: L(imax) = L(imax-1) ******************************
  Aarray_x(imax-1) = 0.0d0
  Carray_x(imax-1) = 0.0d0

  Barray_x(imax-1) = &
       dedt(imax-1) + 0.5d0*dpdt(imax-1)* &
       (1.0d0/rho(imax-1)-1.0d0/rho_p(imax-1))

  Darray_x(imax-1) = &
       eps_p(imax-1) - 0.5d0*p_p(imax-1)*&
       (1.0d0/rho(imax-1)-1.0d0/rho_p(imax-1)) &
       + Qterm(imax-1) &
       - (eps_x(imax-1)+0.5d0*p_x(imax-1)* &
       (1.0d0/rho(imax-1)-1.0d0/rho_p(imax-1)))&
       + bomb_heating(imax-1)*dtime + Ni_heating(imax-1)*dtime

end subroutine matrix_arrays
