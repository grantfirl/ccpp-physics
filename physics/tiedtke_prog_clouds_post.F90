!> \file tiedtke_prog_clouds_post.F90
!!  Contains code to execute after the Tiedtke prognostic cloud scheme

!> This module contains code to execute after the Tiedtke prognostic cloud scheme.
    module tiedtke_prog_clouds_post

    contains

!> \section arg_table_tiedtke_prog_clouds_post_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_post_run.html
!!
    subroutine tiedtke_prog_clouds_post_run (idim, kdim, i_macro, i_temp, do_liq_num, do_ice_num, ldiag3d, qdiag3d, ntqv, ntcw, ntiw, ntclamt, ntlnc, ntinc, dt, qmin, tfreeze, ST, SQ, SL, SI, SA, SN, SNI, dcond_ls_l, dcond_ls_i, d_eros, qa_upd, ql_upd, qi_upd, gt0, gq0, dtend, dtidx, d_eros_l, d_eros_i, nerosc, nerosi, dqcdt, dqidt, errmsg, errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim, i_macro, i_temp, ntqv, ntcw, ntiw, ntclamt, ntlnc, ntinc
      logical, intent(in) :: do_liq_num, do_ice_num, ldiag3d, qdiag3d
      real(kind=kind_phys), intent(in) :: dt, qmin, tfreeze
      real(kind=kind_phys), intent(in) :: ST(:,:), SQ(:,:), SL(:,:), SI(:,:), SA(:,:), SN(:,:), SNI(:,:)
      real(kind=kind_phys), intent(in) :: dcond_ls_l(:,:), dcond_ls_i(:,:), d_eros(:,:), qa_upd(:,:), ql_upd(:,:), qi_upd(:,:)
      
      real(kind=kind_phys), intent(inout) :: gt0(:,:)
      real(kind=kind_phys), intent(inout) :: gq0(:,:,:)
      
      integer, intent(in) :: dtidx(:,:)
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:) 
      
      real(kind=kind_phys), intent(out) :: d_eros_l(:,:), d_eros_i(:,:), nerosc(:,:), nerosi(:,:), dqcdt(:,:), dqidt(:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      integer :: i, k, idtend
      
      real(kind=kind_phys) :: dt_inv
      
      real(kind=kind_phys), dimension(idim,kdim) :: dcond_ls_tot
      
      dt_inv = 1.0/dt
      dcond_ls_tot = 0.0
          
      !gq0(:,:,ntcw)    = gq0(:,:,ntcw)    + SL
      !gq0(:,:,ntiw)    = gq0(:,:,ntiw)    + SI
      gq0(:,:,ntclamt) = gq0(:,:,ntclamt) + SA
      !if (do_liq_num) then
      !  gq0(:,:,ntlnc) = gq0(:,:,ntlnc)   + SN
      !end if
      !if (do_ice_num) then
      !  gq0(:,:,ntinc) = gq0(:,:,ntinc)   + SNi
      !end if
      !gt0              = gt0              + ST
      !gq0(:,:,ntqv)    = gq0(:,:,ntqv)    + SQ
      
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
      
      
      !GJF: In AM4, the following code is executed prior to calling microphysics (specifically MG) in ls_cloud_microphysics
      do k=1,kdim
        do i=1,idim
          dcond_ls_tot(i,k) = dcond_ls_l(i,k) + dcond_ls_i(i,k) 
          D_eros_i(i,k) = -qi_upd(i,k)*D_eros(i,k)*dt_inv
          D_eros_l(i,k) = -ql_upd(i,k)*D_eros(i,k)*dt_inv
          if (ql_upd(i,k) >= qmin .and. do_liq_num) then
            !GJF: can use gq0(:,:,ntlnc) as calculated above, since it should be equivalent to qn_upd
            !     I don't understand why we need to use qa_upd here (as in AM4), instead of gq0(:,:,ntclamt)
            !write(*,*) k, D_eros_l(i,k), gq0(i,k,ntlnc), qa_upd(i,k), D_eros(i,k)
            nerosc(i,k) = D_eros_l(i,k)/ql_upd(i,k)*gq0(i,k,ntlnc)/MAX(0.0001, qa_upd(i,k))
          else
            nerosc(i,k) = 0.
          endif
          if (qi_upd(i,k) >= qmin .and. do_ice_num) then
            !GJF: can use gq0(:,:,ntinc) as calculated above, since it should be equivalent to qni_upd
            nerosi(i,k) = D_eros_i(i,k)/qi_upd(i,k)*gq0(i,k,ntinc)/MAX(0.0001, qa_upd(i,k))
          else
            nerosi(i,k) = 0.
          endif
          !write(*,*) k, nerosc(i,k), nerosi(i,k)
          if (dcond_ls_tot(i,k) > 0.) then
            if (gt0(i,k) <= (tfreeze - 40.) ) then
              dqcdt (i,k) = 0.
              dqidt(i,k) = dcond_ls_tot(i,k)*dt_inv
            else
              dqidt (i,k) = 0.
              dqcdt(i,k) = dcond_ls_tot(i,k)*dt_inv
            endif
          else
            if (gt0(i,k) <= tfreeze) then
              dqcdt(i,k) = MAX(dcond_ls_tot(i,k),-ql_upd(i,k))
              dqidt(i,k) = MAX(dcond_ls_tot(i,k) - dqcdt(i,k), -qi_upd(i,k))
              dqcdt(i,k) = dqcdt(i,k)*dt_inv
              dqidt(i,k) = dqidt(i,k)*dt_inv
            else
              dqidt(i,k) = 0.
              dqcdt(i,k) = MAX(dcond_ls_tot(i,k),-ql_upd(i,k))*dt_inv
            endif
          endif
        end do   
      end do   

    end subroutine tiedtke_prog_clouds_post_run
    
    end module tiedtke_prog_clouds_post