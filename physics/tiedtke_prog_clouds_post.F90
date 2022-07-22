!> \file tiedtke_prog_clouds_post.F90
!!  Contains code to execute after the Tiedtke prognostic cloud scheme

!> This module contains code to execute after the Tiedtke prognostic cloud scheme.
    module tiedtke_prog_clouds_post

    contains

!> \section arg_table_tiedtke_prog_clouds_post_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_post_run.html
!!
    subroutine tiedtke_prog_clouds_post_run (idim, kdim, do_liq_num, do_ice_num, ntqv, ntcw, ntiw, ntclamt, ntlnc, ntinc, ST, SQ, SL, SI, SA, SN, SNI, gt0, gq0, errmsg, errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim, ntqv, ntcw, ntiw, ntclamt, ntlnc, ntinc
      logical, intent(in) :: do_liq_num, do_ice_num
      real(kind=kind_phys), intent(in) :: ST(:,:), SQ(:,:), SL(:,:), SI(:,:), SA(:,:), SN(:,:), SNI(:,:)
      
      real(kind=kind_phys), intent(inout) :: gt0(:,:)
      real(kind=kind_phys), intent(inout) :: gq0(:,:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
            
      gq0(:,:,ntcw)    = gq0(:,:,ntcw)    + SL
      gq0(:,:,ntiw)    = gq0(:,:,ntiw)    + SI
      gq0(:,:,ntclamt) = gq0(:,:,ntclamt) + SA
      if (do_liq_num) then
        gq0(:,:,ntlnc) = gq0(:,:,ntlnc)   + SN
      end if
      if (do_ice_num) then
        gq0(:,:,ntinc) = gq0(:,:,ntinc)   + SNi
      end if
      
      gt0              = gt0              + ST
      gq0(:,:,ntqv)    = gq0(:,:,ntqv)    + SQ

    end subroutine tiedtke_prog_clouds_post_run
    
    end module tiedtke_prog_clouds_post