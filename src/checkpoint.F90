module checkpoint

   implicit none

   !! Simple binary format shared by save_checkpoint/read_checkpoint and the
   !! get_*_from_checkpoint helpers below. magic+version let every entry
   !! point fail loudly (instead of silently misreading) if ever pointed at
   !! a file that isn't a SNEC checkpoint, or one written by an incompatible
   !! version of this module.
   !!
   !! On-disk layout (all unformatted/stream, written/read in this order):
   !!    header: magic(4), version, nzones, ncomps, iBC, rBC_initial
   !!    body  : mass(nzones), r(nzones), vel(nzones), rho(nzones),
   !!            temp(nzones), metallicity(nzones),
   !!            eps(nzones), p(nzones), cs2(nzones), dedt(nzones),
   !!            dpdt(nzones), entropy(nzones), p_rad(nzones),
   !!            [if ncomps>0] comp_details(ncomps,2), comp(nzones,ncomps)
   !!                        zav(ncomps,nzones),
   !!            H_number, He_number, C_number, O_number, Ni_number
   character(len=7), parameter :: ckpt_magic   = "SNEC+fb"
   integer,          parameter :: ckpt_version = 1

contains

   subroutine save_checkpoint()
      use blmod
      use parameters, only: imax, outdir, restart_file

      implicit none

      integer :: unit

      open(newunit=unit, &
           file=trim(adjustl(outdir))//"/"//trim(adjustl(restart_file)), &
           form='unformatted', &
           access='stream', &
           action='write', &
           status='replace')

      !! Write consistently with read_checkpoint

      !----------------------------- header ----------------------------------

      write(unit) ckpt_magic
      write(unit) ckpt_version
      write(unit) imax    ! number of zones in *this* checkpoint
      write(unit) ncomps  ! number of isotopes in *this* checkpoint
      write(unit) iBC
      write(unit) rBC_initial

      !--------------------------- grid coordinates ---------------------------
      write(unit) mass(1:imax)
      write(unit) r(1:imax)

      !----------------------------- hydro state -------------------------------
      write(unit) vel(1:imax)
      write(unit) rho(1:imax)
      write(unit) temp(1:imax)
      ! Note: ye and abar are removed here because they are calculated explicitly from composition
      write(unit) metallicity(1:imax)
      write(unit) eps(1:imax)
      write(unit) p(1:imax)
      write(unit) cs2(1:imax)
      write(unit) dedt(1:imax)
      write(unit) dpdt(1:imax)
      write(unit) entropy(1:imax)
      write(unit) p_rad(1:imax)

      !----------------------------- composition --------------------------------
      if (ncomps.gt.0) then
         write(unit) comp_details(1:ncomps,1:2)  !(A,Z) of each isotope
         write(unit) comp(1:imax,1:ncomps)        !mass fractions
         write(unit) zav(1:ncomps,1:imax)         !average ionic charge
      endif

      write(unit) H_number, He_number, C_number, O_number, Ni_number

      close(unit)

      write(*,*) "Wrote checkpoint with ", imax, " zones to ", trim(adjustl(outdir))//"/"//trim(adjustl(restart_file))

   end subroutine save_checkpoint


   subroutine read_checkpoint()
      use blmod, only: mass, cmass, vel, rho, temp, ncomps, ye, abar, comp_details,&
         eps, p, cs2, dedt, dpdt, entropy, zav, p_rad, &
         comp, metallicity, H_number, He_number, C_number, O_number, Ni_number, iBC, rBC_initial
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
      real*8,allocatable :: pmetallicity(:)
      real*8,allocatable :: peps(:), pp(:), pcs2(:), pdedt(:), pdpdt(:)
      real*8,allocatable :: pentropy(:), pp_rad(:)
      real*8,allocatable :: pcomp(:,:), pzav(:,:)

      character(len=256) :: filename
      character(len=7) :: magic
      integer :: version
      integer :: ncomps_chk

      filename = trim(adjustl(outdir))//"/"//trim(adjustl(restart_file))

      open(newunit=unit, &
           file=trim(filename), &
           form='unformatted', &
           access='stream', &
           action='read', &
           status='old')

      !---------------------------- header ----------------------------------
      read(unit) magic
      read(unit) version

      if (magic.ne.ckpt_magic) then
         write(*,*) "Error: ", trim(filename)," does not look like a SNEC checkpoint file."
         stop
      endif
      if (version.ne.ckpt_version) then
         write(*,*) "Error: checkpoint file ", trim(filename), &
            " has version ",version," but this code reads version ",ckpt_version
         stop
      endif

      read(unit) profile_zones
      read(unit) ncomps_chk
      read(unit) iBC
      read(unit) rBC_initial

      if (ncomps_chk.ne.ncomps) then
         write(*,*) "Error: checkpoint file ", trim(filename)," has ncomps=",ncomps_chk, &
            " but this run was set up with ncomps=",ncomps
         write(*,*) "(ncomps should have already been set from get_ncomps_from_checkpoint)"
         stop
      endif

      write(*,*) "We have ",profile_zones, "checkpoint zones."

      !------------------------- read the checkpoint arrays --------------------
      allocate(pmass(profile_zones))
      allocate(pradius(profile_zones))
      allocate(pvel(profile_zones))
      allocate(prho(profile_zones))
      allocate(ptemp(profile_zones))
      allocate(pmetallicity(profile_zones))
      allocate(peps(profile_zones))
      allocate(pp(profile_zones))
      allocate(pcs2(profile_zones))
      allocate(pdedt(profile_zones))
      allocate(pdpdt(profile_zones))
      allocate(pentropy(profile_zones))
      allocate(pp_rad(profile_zones))

      read(unit) pmass(1:profile_zones)
      read(unit) pradius(1:profile_zones)
      read(unit) pvel(1:profile_zones)
      read(unit) prho(1:profile_zones)
      read(unit) ptemp(1:profile_zones)
      read(unit) pmetallicity(1:profile_zones)
      read(unit) peps(1:profile_zones)
      read(unit) pp(1:profile_zones)
      read(unit) pcs2(1:profile_zones)
      read(unit) pdedt(1:profile_zones)
      read(unit) pdpdt(1:profile_zones)
      read(unit) pentropy(1:profile_zones)
      read(unit) pp_rad(1:profile_zones)

      if(ncomps.gt.0) then
         allocate(pcomp(profile_zones,ncomps))
         allocate(pzav(ncomps,profile_zones))

         read(unit) comp_details(1:ncomps,1:2)  !not grid-dependent, copy directly
         read(unit) pcomp(1:profile_zones,1:ncomps)
         read(unit) pzav(1:ncomps,1:profile_zones)
      endif

      read(unit) H_number, He_number, C_number, O_number, Ni_number

      close(unit)

      !----------------- map the checkpoint onto the current grid --------------
      do i=iBC,imax-1
         rho(i) = prho(i)
         vel(i) = pvel(i)
         temp(i) = ptemp(i)
         metallicity(i) = pmetallicity(i)
         eps(i) = peps(i)
         p(i) = pp(i)
         cs2(i) = pcs2(i)
         dedt(i) = pdedt(i)
         dpdt(i) = pdpdt(i)
         entropy(i) = pentropy(i)
         p_rad(i) = pp_rad(i)
      enddo

      rho(imax)         = 0.0d0
      temp(imax)        = 0.0d0
      eps(imax)         = 0.0d0
      p(imax)           = 0.0d0
      metallicity(imax) = metallicity(imax-1)
      cs2(imax)         = cs2(imax-1)
      dedt(imax)        = dedt(imax-1)
      dpdt(imax)        = dpdt(imax-1)
      entropy(imax)     = entropy(imax-1)
      p_rad(imax)       = p_rad(imax-1)

      if (iBC>1) then
         vel(1:iBC-1)         = 0.0d0
         rho(1:iBC-1)         = rho(iBC)
         temp(1:iBC-1)        = temp(iBC)
         eps(1:iBC-1)         = eps(iBC)
         p(1:iBC-1)           = p(iBC)
         metallicity(1:iBC-1) = metallicity(iBC)
         cs2(1:iBC-1)         = cs2(iBC)
         dedt(1:iBC-1)        = dedt(iBC)
         dpdt(1:iBC-1)        = dpdt(iBC)
         entropy(1:iBC-1)     = entropy(iBC)
         p_rad(1:iBC-1)       = p_rad(iBC)
         ye(1:iBC-1)=ye(iBC)
         abar(1:iBC-1) = 1
      end if

      ! Map compositions
      do l=1,ncomps
         comp(1:imax,l) = pcomp(1:imax,l)
         zav(l,1:imax)  = pzav(l,1:imax)
      enddo

      ! recalculate ye and abar
      do i=iBC,imax
         ye(i) = sum(comp(i,1:ncomps) * &
             (comp_details(1:ncomps,2)/comp_details(1:ncomps,1)))
         abar(i) = 1.0d0/sum(comp(i,1:ncomps)/comp_details(1:ncomps,1))
      enddo


      deallocate(pmass)
      deallocate(pradius)
      deallocate(pvel)
      deallocate(prho)
      deallocate(ptemp)
      deallocate(pmetallicity)
      deallocate(peps)
      deallocate(pp)
      deallocate(pcs2)
      deallocate(pdedt)
      deallocate(pdpdt)
      deallocate(pentropy)
      deallocate(pp_rad)
      if(ncomps.gt.0) then
         deallocate(pcomp)
         deallocate(pzav)
      endif

      if(eoskey.eq.2) then
         call init_ionpot
      endif

      write(*,*) "Restarted from checkpoint: ", trim(adjustl(outdir))//"/"//trim(adjustl(restart_file))

   end subroutine read_checkpoint

   subroutine get_ncomps_from_checkpoint(xncomps)
      use parameters, only: outdir, restart_file

      implicit none

      !Output:
      integer              :: xncomps

      !Local:
      integer              :: unit
      character(len=7)     :: magic
      integer              :: version
      integer              :: profile_zones

      character(len=254)   :: filename

      filename = trim(adjustl(outdir))//"/"//trim(adjustl(restart_file))

      open(newunit=unit, &
         file=trim(filename), &
         form='unformatted', &
         access='stream', &
         action='read', &
         status='old')

      read(unit) magic
      read(unit) version

      if (magic.ne.ckpt_magic) then
         write(*,*) "Error: ",trim(filename)," does not look like a SNEC checkpoint file."
         stop
      endif
      if (version.ne.ckpt_version) then
         write(*,*) "Error: checkpoint file ",filename, &
            " has version ",version," but this code reads version ",ckpt_version
         stop
      endif

      read(unit) profile_zones
      read(unit) xncomps

      close(unit)

   end subroutine get_ncomps_from_checkpoint


   subroutine get_inner_outer_mass_from_checkpoint(inner_mass, outer_mass)
      use parameters, only: outdir, restart_file
      use blmod, only: iBC, rBC_initial
      implicit none

      !Output:
      real*8               :: inner_mass
      real*8               :: outer_mass

      !Local:
      integer              :: profile_zones
      integer              :: ncomps_chk
      integer              :: unit
      character(len=7)     :: magic
      integer              :: version
      real*8, allocatable  :: pmass(:)

      !! * Write consistently with get_inner_outer_mass_from_profile in read_profile.F90

      open(newunit=unit, &
         file=trim(adjustl(outdir))//"/"//trim(adjustl(restart_file)), &
         form='unformatted', &
         access='stream', &
         action='read', &
         status='old')

      read(unit) magic
      read(unit) version

      if (magic.ne.ckpt_magic) then
         write(*,*) "Error: ",trim(adjustl(restart_file))," does not look like a SNEC checkpoint file."
         stop
      endif
      if (version.ne.ckpt_version) then
         write(*,*) "Error: checkpoint file ",trim(adjustl(restart_file)), &
            " has version ",version," but this code reads version ",ckpt_version
         stop
      endif

      read(unit) profile_zones
      read(unit) ncomps_chk
      read(unit) iBC
      read(unit) rBC_initial

      allocate(pmass(profile_zones))
      read(unit) pmass(1:profile_zones)
      close(unit)

      outer_mass = pmass(profile_zones)
      inner_mass = pmass(1)

      deallocate(pmass)

   end subroutine get_inner_outer_mass_from_checkpoint

   subroutine get_inner_outer_radius_from_checkpoint(inner_mass, inner_radius, outer_radius)
      use parameters, only: outdir, restart_file
      use blmod, only: iBC, rBC_initial
      implicit none

      !Input:

      !Input:
      real*8               :: inner_mass

      !Output:
      real*8               :: inner_radius
      real*8               :: outer_radius

      !Local:
      integer              :: profile_zones
      integer              :: ncomps_chk
      integer              :: unit
      character(len=7)     :: magic
      integer              :: version
      real*8, allocatable  :: pradius(:)
      real*8, allocatable  :: pmass(:)

      !! * Write consistently with get_inner_outer_mass_from_profile in read_profile.F90

      open(newunit=unit, &
         file=trim(adjustl(outdir))//"/"//trim(adjustl(restart_file)), &
         form='unformatted', &
         access='stream', &
         action='read', &
         status='old')

      read(unit) magic
      read(unit) version

      if (magic.ne.ckpt_magic) then
         write(*,*) "Error: ",trim(adjustl(restart_file))," does not look like a SNEC checkpoint file."
         stop
      endif
      if (version.ne.ckpt_version) then
         write(*,*) "Error: checkpoint file ",trim(adjustl(restart_file)), &
            " has version ",version," but this code reads version ",ckpt_version
         stop
      endif

      read(unit) profile_zones
      read(unit) ncomps_chk
      read(unit) iBC
      read(unit) rBC_initial

      allocate(pmass(profile_zones))
      allocate(pradius(profile_zones))

      read(unit) pmass(1:profile_zones)
      read(unit) pradius(1:profile_zones)
      close(unit)

      inner_radius = rBC_initial
      outer_radius = pradius(profile_zones)
      print *, "outer_radius", outer_radius
      deallocate(pradius)
      deallocate(pmass)

   end subroutine get_inner_outer_radius_from_checkpoint

end module checkpoint
