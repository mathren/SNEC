! mock_modules.f90
! Mocking the dependencies to make the checkpoint module testable standalone
! ==============================================================================
! MAIN TEST DRIVER
! ==============================================================================
program test_checkpoint
   use parameters
   use blmod
   use checkpoint
   implicit none

   integer :: test_ncomps
   real(8) :: test_m_in, test_m_out
   real(8) :: test_r_in, test_r_out
   integer :: i

   print *, "=== Starting Checkpoint Unit Test ==="

   ! 1. Set the file path BEFORE anything else touches I/O
   outdir                 = "."
   profile_name           = "../profiles/15Msol_RSG.short"
   composition_profile_name      = "../profiles/15Msol_RSG.iso.dat"
   grid_pattern_name      = "../profiles/Grid_pattern_15RSG.dat"
   initial_data           =  "Thermal_Bomb"
   bomb_mode              = 1
   final_energy           = 0.0d0
   bomb_tstart            = 0.0d0
   bomb_tend              = 1.0d-2
   bomb_mass_spread       = 0.1
   bomb_start_point       = 1
   imax                   = 3208
   gridding               = "from_file_by_mass"
   mass_excision          = 1
   mass_excised           = 1.4
   radiation              = 1
   eoskey                 = 2
   Ni_by_hand             = 0
   Ni_switch              = 0
   Ni_mass                = 0.0d0
   Ni_boundary_mass       = 2.9d0
   Ni_period              = 1.0d4
   saha_ncomps            = 3
   boxcar_smoothing       = 1
   ntmax                  = 10
   tend                   = 100
   max_t_dense_out        = 100.0d0
   dtout                  = 100.0d0
   dtout_scalar           = 100.0d0
   dtout_check            = 3600.0d0

   ntout                  = -1
   ntout_scalar           = -1
   ntout_check            = 10000

   ntinfo                 = 1000

   dtmin                  = 1.0d-10
   dtmax                  = 1.0d2
   sedov                  = 0
   innerBC                = "inflow"
   restart                = 0
   restart_file           = "test.checkpoint"


   ! 3. Populate live memory arrays with distinct mock values
   do i = 1, imax
      mass(i)        = i * 1.0d0
      r(i)           = i * 10.0d0
      vel(i)         = i * 0.1d0
      rho(i)         = i * 5.0d0
      temp(i)        = i * 100.0d0
      metallicity(i) = 0.02d0
      eps(i)         = 1.1d0
      p(i)           = 2.2d0
      cs2(i)         = 3.3d0
      dedt(i)        = 4.4d0
      dpdt(i)        = 5.5d0
      entropy(i)     = 6.6d0
      p_rad(i)       = 7.7d0
      comp(i, 1)     = 0.75d0
      comp(i, 2)     = 0.25d0
      zav(1, i)      = 1.0d0
      zav(2, i)      = 2.0d0
   enddo

   comp_details(1, 1) = 1.0d0 ! A for H
   comp_details(1, 2) = 1.0d0 ! Z for H
   comp_details(2, 1) = 4.0d0 ! A for He
   comp_details(2, 2) = 2.0d0 ! Z for He

   ! 4. WRITE THE CHECKPOINT FILE FIRST
   print *, "--> Step 1: Writing checkpoint file..."
   call save_checkpoint()

   ! 5. TEST PARTIAL-READ HELPER ROUTINES
   print *, "--> Step 2: Testing partial-read helper functions..."

   call get_ncomps_from_checkpoint(test_ncomps)
   if (test_ncomps /= ncomps) then
      print *, "❌ FAIL: ncomps mismatch! Expected:", ncomps, "Got:", test_ncomps
      stop 1
   endif

   call get_inner_outer_mass_from_checkpoint(test_m_in, test_m_out)
   if (abs(test_m_in - mass(1)) > 1.0d-5 .or. abs(test_m_out - mass(imax)) > 1.0d-5) then
      print *, "❌ FAIL: Mass bounds mismatch!"
      stop 1
   endif

   call get_inner_outer_radius_from_checkpoint(test_m_in, test_r_in, test_r_out)
   if (abs(test_r_in - rBC_initial) > 1.0d-5 .or. abs(test_r_out - r(imax)) > 1.0d-5) then
      print *, "❌ FAIL: Radius bounds mismatch!"
      stop 1
   endif

   ! 6. CLEAR LIVE MEMORY
   ! We blank out the arrays to verify that read_checkpoint actually pulls data from the file
   rho = 0.0d0
   temp = 0.0d0
   comp = 0.0d0

   ! 7. READ CHECKPOINT BACK IN FULL
   print *, "--> Step 3: Testing full read_checkpoint..."
   call read_checkpoint()

   ! 8. VALIDATION
   if (abs(rho(5) - 25.0d0) > 1.0d-5 .or. abs(temp(5) - 500.0d0) > 1.0d-5) then
      print *, "❌ FAIL: Full state restore failed data validation!"
      stop 1
   endif

   print *, "✅ SUCCESS: All checkpoint routines working flawlessly!"

end program test_checkpoint
