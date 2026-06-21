module checkpoint

   implicit none

contains

   subroutine save_checkpoint(filename)

      use blmod

      implicit none

      character(len=*), intent(in) :: filename

      integer :: unit

      open(newunit=unit, &
         file=filename, &
         form='unformatted', &
         access='stream', &
         status='replace')

      !! Write consistently with read_checkpoint
      close(unit)

   end subroutine save_checkpoint



   subroutine read_checkpoint(filename)

      use blmod, only: mass, cmass, vel, rho, temp, ncomps, ye, abar, comp_details,&
         eps, p, cs2, dedt, dpdt, entropy, zav, p_rad
      use parameters
      use physical_constants
      use eosmodule, only: init_ionpot
      implicit none

      integer :: profile_zones
      integer :: i,l
      integer :: ibuffer
      integer :: keytemp, keyerr
      integer :: unit

      real*8,allocatable :: pmass(:), pradius(:), ptemp(:), prho(:), pvel(:)

      character(len=*), intent(in) :: filename

      open(newunit=unit, &
         file=filename, &
         form='unformatted', &
         access='stream', &
         status='old')

      !! * Write consistently with read_checkpoint (to be written too)
      !! * to produce the same output as read_profile.F90 read_profile subroutine
      close(unit)

   end subroutine read_checkpoint


   subroutine get_ncomps_from_checkpoint(restart_file,xncomps)

      implicit none

      !Input:
      character(*)         :: restart_file

      !Output:
      integer              :: xncomps

      !! * Write consistently with get_ncomps_from_profile in read_profile.F90

   end subroutine get_ncomps_from_checkpoint

   subroutine get_inner_outer_mass_from_checkpoint(restart_file, inner_mass, outer_mass)
      implicit none

      !Input:
      character(*)         :: restart_file

      !Output:
      real*8               :: inner_mass
      real*8               :: outer_mass

      !Local:
      integer              :: profile_zones
      integer              :: i,ibuffer
      real*8               :: buffer
      real*8, allocatable  :: pmass(:)

      !! * Write consistently with get_inner_outer_mass_from_profile in read_profile.F90

   end subroutine get_inner_outer_mass_from_checkpoint


   subroutine get_inner_outer_radius_from_checkpoint(restart_file,inner_mass, inner_radius, outer_radius)
      implicit none

      !Input:

      !Input:
      character(*)         :: restart_file
      real*8               :: inner_mass

      !Output:
      real*8               :: inner_radius
      real*8               :: outer_radius

      !Local:
      integer              :: profile_zones
      integer              :: i,ibuffer
      real*8               :: buffer
      real*8, allocatable  :: pradius(:)
      real*8, allocatable  :: pmass(:)

      !! * Write consistently with get_inner_outer_mass_from_profile in read_profile.F90

   end subroutine get_inner_outer_radius_from_checkpoint
end module checkpoint
