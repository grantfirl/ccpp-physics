!> \file tiedtke_prog_clouds_post.F90
!!  Contains code to execute after the Tiedtke prognostic cloud scheme

!> This module contains code to execute after the Tiedtke prognostic cloud scheme.
    module tiedtke_prog_clouds_post

    contains

!> \section arg_table_tiedtke_prog_clouds_post_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_post_run.html
!!
    subroutine tiedtke_prog_clouds_post_run (idim, kdim, i_macro, i_temp, do_liq_num, do_ice_num, do_mynnedmf, &
      ldiag3d, qdiag3d, ntqv, ntcw, ntiw, ntrw, ntsw, ntgl, ntclamt, ntlnc, ntinc, con_hvap, con_hfus, con_cp, dt, qmin, tfreeze, &
      ST, SQ, SL, SI, SA, SN, SNI, dcond_ls_l, dcond_ls_i, d_eros, qa_upd, ql_upd, qi_upd, gt0, gq0, &
      dtend, dtidx, d_eros_l, d_eros_i, nerosc, nerosi, dqcdt, dqidt, dqadt_pbl, dqcdt_pbl, ovhd_cldcov, &
      ap, ap_cld, ap_clr, errmsg, errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim, i_macro, i_temp, ntqv, ntcw, ntiw, ntrw, ntsw, ntgl, ntclamt, ntlnc, ntinc
      logical, intent(in) :: do_liq_num, do_ice_num, ldiag3d, qdiag3d, do_mynnedmf
      real(kind=kind_phys), intent(in) :: dt, qmin, tfreeze, con_hvap, con_hfus, con_cp
      real(kind=kind_phys), intent(in) :: SA(:,:), SN(:,:), SNI(:,:)
      real(kind=kind_phys), intent(inout) :: ST(:,:), SQ(:,:), SL(:,:), SI(:,:)
      real(kind=kind_phys), intent(in) :: dcond_ls_l(:,:), dcond_ls_i(:,:), d_eros(:,:), qa_upd(:,:), ql_upd(:,:), qi_upd(:,:), dqadt_pbl(:,:), dqcdt_pbl(:,:)
      
      real(kind=kind_phys), intent(inout) :: gt0(:,:)
      real(kind=kind_phys), intent(inout) :: gq0(:,:,:)
      
      integer, intent(in) :: dtidx(:,:)
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:) 
      
      real(kind=kind_phys), intent(out) :: d_eros_l(:,:), d_eros_i(:,:), nerosc(:,:), nerosi(:,:), dqcdt(:,:), dqidt(:,:), ovhd_cldcov(:,:), ap(:,:), ap_cld(:,:), ap_clr(:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      integer :: i, k, idtend
      
      real(kind=kind_phys), parameter :: delta = 1.0E-6
      real(kind=kind_phys) :: dt_inv
      
      real(kind=kind_phys), dimension(idim,kdim) :: dcond_ls_tot
      
      logical :: precip_found
      real(kind=kind_phys), dimension(idim,kdim) :: gq0_precip
      real(kind=kind_phys) :: del_ovhd_cldcov, del_ap_cld2clr, del_ap_clr2cld
      
      dt_inv = 1.0/dt
      dcond_ls_tot = 0.0
      
      !apply tendencies from realizability only
      
      ! gq0(:,:,ntcw)    = gq0(:,:,ntcw)    + SL
      ! gq0(:,:,ntiw)    = gq0(:,:,ntiw)    + SI
      ! if (do_mynnedmf) then
      !    gq0(:,:,ntclamt) = gq0(:,:,ntclamt) + SA + dqadt_pbl*dt
      !    !print*,dqadt_pbl
      !    !test
      !    !do k = 1,kdim
      !    !  do i = 1,idim
      !    !    gq0(i,k,ntclamt) = gq0(i,k,ntclamt) + SA(i,k) + dqadt_pbl(i,k)*dt
      !    !  enddo
      !    !enddo
      ! else
      ! gq0(:,:,ntclamt) = gq0(:,:,ntclamt) + SA
      ! endif
      ! if (do_liq_num) then
      !  gq0(:,:,ntlnc) = gq0(:,:,ntlnc)   + SN
      ! end if
      ! if (do_ice_num) then
      !  gq0(:,:,ntinc) = gq0(:,:,ntinc)   + SNi
      ! end if
      ! gt0(:,:)         = gt0(:,:)         + ST
      ! gq0(:,:,ntqv)    = gq0(:,:,ntqv)    + SQ
      !write(*,*) 'SL', SL
      !write(*,*) 'SI', SI
      !write(*,*) 'ST', ST, gt0
      !write(*,*) 'SQ', SQ, gq0(:,:,ntqv)
      
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
              dqcdt(i,k) = MAX(dcond_ls_tot(i,k)*dt_inv,-ql_upd(i,k)*dt_inv)
            endif
          endif
        end do   
      end do
      
      do k=1,kdim
        do i=1,idim
          SL(i,k) = SL(i,k) + (dqcdt(i,k) + D_eros_l(i,k))*dt
          SI(i,k) = SI(i,k) + (dqidt(i,k) + D_eros_i(i,k))*dt
          SQ(i,k) = SQ(i,k) - SL(i,k) - SI(i,k)
          ST(i,k) = ST(i,k) + ((con_hvap + con_hfus)/con_cp)*SI(i,k) + &
            (con_hvap/con_cp)*SL(i,k)
        end do
      end do
      gq0(:,:,ntcw)    = gq0(:,:,ntcw)    + SL
      gq0(:,:,ntiw)    = gq0(:,:,ntiw)    + SI
      
      gq0(:,:,ntclamt) = gq0(:,:,ntclamt) + SA
      gq0(:,:,ntqv)    = gq0(:,:,ntqv)    + SQ
      gt0(:,:)         = gt0(:,:)         + ST
      
      
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
            
      !This expression “gives maximum overlap for clouds in adjacent levels with cloud fraction monotonically increasing or decreasing with height, 
      !and random overlap for clouds either separated by clear levels of for levels of changing sign in the vertical gradient of cloud fraction” (Jakob and Klein, 2000)
      do i=1,idim
        ovhd_cldcov(i,kdim) = gq0(i,kdim,ntclamt)
        do k=kdim-1,1,-1
          ovhd_cldcov(i,k) = 1.0 - (1.0 - ovhd_cldcov(i,k+1))*&
            (1.0 - MAX(gq0(i,k,ntclamt),gq0(i,k+1,ntclamt)))/(1.0 - MIN(gq0(i,k+1,ntclamt),1.0 - delta))          
        end do
      end do
      
      !Chosson_et_al (2014) Appendix A algorithm for precip fraction
      ap(:,:) = 0.0
      ap_cld(:,:) = 0.0
      ap_clr(:,:) = 0.0
      !calculate total precipitation mixing ratio (depends on MP scheme, but just for Thompson MP for now)
      gq0_precip(:,:) = gq0(:,:,ntrw) + gq0(:,:,ntsw) + gq0(:,:,ntgl)
      do i=1,idim
        precip_found = .false.
        do k=kdim,1,-1
          if (gq0_precip(i,k) > 0.0) then !should this be a different threshold?
            if (.not. precip_found) then
              ap_cld(i,k) = gq0(i,k,ntclamt)
              precip_found = .true.
            else
              del_ovhd_cldcov = ovhd_cldcov(i,k) - ovhd_cldcov(i,k+1)
              del_ap_cld2clr = ap_cld(i,k+1) - MIN(gq0(i,k,ntclamt) - del_ovhd_cldcov,ap_cld(i,k+1))
              del_ap_clr2cld = MAX(0.0,MIN(ap_clr(i,k+1),gq0(i,k,ntclamt) - del_ovhd_cldcov - gq0(i,k+1,ntclamt)))
              ap_cld(i,k) = gq0(i,k+1,ntclamt) + del_ap_clr2cld - del_ap_cld2clr !first term is different than Firl(2009), which has ap_cld(i,k+1)
              ap_clr(i,k) = ap_clr(i,k+1) - del_ap_clr2cld + del_ap_cld2clr
            end if
            ap(i,k) = ap_cld(i,k) + ap_clr(i,k)
          end if
        end do
      end do

    end subroutine tiedtke_prog_clouds_post_run
    
    end module tiedtke_prog_clouds_post