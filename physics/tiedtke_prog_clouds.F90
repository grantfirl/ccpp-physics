!> \file tiedtke_prog_clouds.F90
!!  Contains the Tiedtke prognostic cloud scheme

!> This module contains the CCPP-compliant Tiedtke prognostic cloud scheme.
      module tiedtke_prog_clouds

      contains

!> \section arg_table_tiedtke_prog_clouds_init Argument Table
!! \htmlinclude tiedtke_prog_clouds_init.html
!!
      subroutine tiedtke_prog_clouds_init ()
         use machine, only : kind_phys
         
         implicit none

         character(len=*),     intent(out) :: errmsg
         integer,              intent(out) :: errflg
         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

      end subroutine tiedtke_prog_clouds_init

!> \section arg_table_tiedtke_prog_clouds_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_run.html
!!
      subroutine tiedtke_prog_clouds_run (errmsg,errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg


! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0


      end subroutine tiedtke_prog_clouds_run


    end module tiedtke_prog_clouds
