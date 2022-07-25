!> \file tiedtke_prog_clouds_post.F90
!!  Contains code to execute after the Tiedtke prognostic cloud scheme

!> This module contains code to execute after the Tiedtke prognostic cloud scheme.
    module tiedtke_prog_clouds_post

    contains

!> \section arg_table_tiedtke_prog_clouds_post_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_post_run.html
!!
    subroutine tiedtke_prog_clouds_post_run (idim, kdim, i_macro, i_temp, do_liq_num, do_ice_num, ldiag3d, qdiag3d, ntqv, ntcw, ntiw, ntclamt, ntlnc, ntinc, ST, SQ, SL, SI, SA, SN, SNI, gt0, gq0, dtend, dtidx, errmsg, errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim, i_macro, i_temp, ntqv, ntcw, ntiw, ntclamt, ntlnc, ntinc
      logical, intent(in) :: do_liq_num, do_ice_num, ldiag3d, qdiag3d
      real(kind=kind_phys), intent(in) :: ST(:,:), SQ(:,:), SL(:,:), SI(:,:), SA(:,:), SN(:,:), SNI(:,:)
      
      real(kind=kind_phys), intent(inout) :: gt0(:,:)
      real(kind=kind_phys), intent(inout) :: gq0(:,:,:)
      
      integer, intent(in) :: dtidx(:,:)
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:) 
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      integer :: idtend
            
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
      
      if (ldiag3d) then
        idtend = dtidx(i_temp,i_macro)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + ST
        endif
        
        if (qdiag3d) then
          idtend = dtidx(100+ntqv,i_macro)
          if(idtend>0) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + SQ
          endif
          
          idtend = dtidx(100+ntcw,i_macro)
          if(idtend>0) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + SL
          endif
          
          idtend = dtidx(100+ntiw,i_macro)
          if(idtend>0) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + SI
          endif
          
          idtend = dtidx(100+ntclamt,i_macro)
          if(idtend>0) then
             dtend(:,:,idtend) = dtend(:,:,idtend) + SA
          endif
          
          if (do_liq_num) then
            idtend = dtidx(100+ntlnc,i_macro)
            if(idtend>0) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + SN
            endif
          end if
          
          if (do_ice_num) then
            idtend = dtidx(100+ntinc,i_macro)
            if(idtend>0) then
               dtend(:,:,idtend) = dtend(:,:,idtend) + SNi
            endif
          end if
          
        end if
      end if

    end subroutine tiedtke_prog_clouds_post_run
    
    end module tiedtke_prog_clouds_post