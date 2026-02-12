subroutine read_BolCorr

  use blmod, only: bol_corr, temp_bol_corr, nlines_bol_corr
  use parameters

  implicit none

  character(len=256) :: snec_dir



  integer :: nlines, i

  !------------------------------------------------------------------------------

  call get_environment_variable('SNEC_DIR', snec_dir)

  open(666,file=trim(snec_dir) // trim("tables/BolCorr.dat"),status='unknown',&
       form='formatted',action='read')

  read(666,*)
  nlines = 0
  do
     read(666,*,end=15)
     nlines = nlines + 1
  end do
15 close(666)

  nlines_bol_corr = nlines

  allocate(temp_bol_corr(nlines))
  allocate(bol_corr(nlines,11))

  open(666,file=trim(snec_dir) // trim("tables/BolCorr.dat"),status='unknown',&
       form='formatted',action='read')
  read(666,*)
  do i=1,nlines
     read(666,*) temp_bol_corr(i), bol_corr(i,1:11)
  enddo
  close(666)

end subroutine read_BolCorr
