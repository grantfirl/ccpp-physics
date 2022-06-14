!> \file tiedtke_prog_clouds_post.F90
!!  Contains code to execute after the Tiedtke prognostic cloud scheme

!> This module contains code to execute after the Tiedtke prognostic cloud scheme.
    module tiedtke_prog_clouds_post

    contains

!> \section arg_table_tiedtke_prog_clouds_post_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_post_run.html
!!
    subroutine tiedtke_prog_clouds_post_run (idim, kdim, errmsg, errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      !GJF: need to at least call lscloud_driver.F90/update_fields_and_tendencies

    end subroutine tiedtke_prog_clouds_post_run
    
    end module tiedtke_prog_clouds_post