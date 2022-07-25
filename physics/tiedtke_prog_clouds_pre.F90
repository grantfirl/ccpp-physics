!> \file tiedtke_prog_clouds_pre.F90
!!  Contains the preparation code for Tiedtke prognostic cloud scheme

!> This module contains the CCPP-compliant Tiedtke prognostic cloud scheme.
    module tiedtke_prog_clouds_pre

    contains

!> \section arg_table_tiedtke_prog_clouds_pre_init Argument Table
!! \htmlinclude tiedtke_prog_clouds_pre_init.html
!!    
    subroutine tiedtke_prog_clouds_pre_init(hlv, hlf, rvgas, cp, errmsg, errflg)
      use machine  , only : kind_phys
      use ccpp_saturation, only: get_qs_init
      
      implicit none
      
      real(kind=kind_phys), intent(in)    :: hlv, hlf, rvgas, cp
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      call get_qs_init(hlv, hlf, rvgas, cp) 
      
    end subroutine tiedtke_prog_clouds_pre_init

!> \section arg_table_tiedtke_prog_clouds_pre_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_pre_run.html
!!
    subroutine tiedtke_prog_clouds_pre_run (idim, kdim, pdf_org, pdfcld, super_ice_opt, do_liq_num, do_ice_num, do_aero_eros, hlv, hlf, rvgas, cp, rdgas, qmin, tin, qin, pfull, cnvc, ql_in, qi_in, qa_in, qn_in, qni_in, dtdt_rad_turb_gwd, sa, sl, si, sq, st, sn, sni, drop1, qs, qsl, qsi, dqsdT, gamma, radturbten, ql_upd, qi_upd, qa_upd, qn_upd, qni_upd, U_ca, U01, airdens, errmsg, errflg)

      use machine  , only : kind_phys
      use ccpp_saturation, only: get_qs
      
      implicit none
      
      integer,              intent(in)    :: idim, kdim

      logical,              intent(in)    :: pdf_org, pdfcld, do_liq_num, do_ice_num, do_aero_eros
      integer,              intent(in)    :: super_ice_opt
      
      real(kind=kind_phys), intent(in)    :: hlv, hlf, rvgas, cp, rdgas, qmin
      real(kind=kind_phys), intent(in)    :: tin(:,:), qin(:,:), pfull(:,:), cnvc(:,:), ql_in(:,:), qi_in(:,:), qa_in(:,:), qn_in(:,:), qni_in(:,:)
      real(kind=kind_phys), intent(in)    :: dtdt_rad_turb_gwd(:,:)
      real(kind=kind_phys), intent(inout) :: sa(:,:), sl(:,:), si(:,:), sq(:,:), st(:,:), sn(:,:), sni(:,:), drop1(:,:)
      real(kind=kind_phys), intent(out)   :: qs(:,:), qsl(:,:), qsi(:,:), dqsdT(:,:), gamma(:,:), radturbten(:,:)
      real(kind=kind_phys), intent(out)   :: ql_upd(:,:), qi_upd(:,:), qa_upd(:,:), qn_upd(:,:), qni_upd(:,:)
      real(kind=kind_phys), intent(out)   :: U_ca(:,:), U01(:,:), airdens(:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      !temporary until set up in GFS_typedefs
      
      integer :: i,k
      integer :: alg_choice, alg_flatau_92
      
      real(kind=kind_phys) :: conv_frac_max, tfreeze, hls
      real(kind=kind_phys) :: es(idim,kdim), esl(idim,kdim), esi(idim,kdim)
      real(kind=kind_phys) :: rh_wtd_conv_area(idim,kdim), convective_area(idim,kdim)
      real(kind=kind_phys) :: qrf(idim,kdim), env_qv(idim,kdim), env_fraction(idim,kdim), humidity_ratio(idim,kdim)
      
      logical :: ql_too_small(idim,kdim), qi_too_small(idim,kdim)
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      alg_choice = 1
      alg_flatau_92 = 1
      conv_frac_max = 1.0
      tfreeze = 273.16 !should be a physical constant coming in from host
      hls = hlv + hlf
            
      call get_qs(tin, pfull, alg_choice, alg_flatau_92, qs, qsl, qsi, es, esl, esi, dqsdT, gamma)
      
      !GFJ: calc rh_wtd_conv_area (relative humidity-weighted conv area; 
      !this is useful if a deep convective scheme assumes a RH other than 100% for convective cloud area fraction)
      !from AM4 convection_driver.F90/compute_convective_area
      do k=1, kdim
        do i=1, idim
          !assume 100% RH in convective cloud cover area
          rh_wtd_conv_area(i,k) = cnvc(i,k)
        end do
      end do
      
      !from AM4 convection_driver.F90/compute_convective_area
      !calculate convective_humidity_ratio to prepare to calculate U_ca
      do k=1, kdim
        do i=1, idim
          convective_area(i,k) = min(cnvc(i,k), conv_frac_max)
          env_fraction(i,k) = 1.0 - convective_area(i,k)
          qrf(i,k) = max(qin(i,k), 0.0)
          env_qv(i,k) = qrf(i,k) - qs(i,k)*rh_wtd_conv_area(i,k)
          if (qrf(i,k) /= 0.0 .and. env_qv(i,k) > 0.0) then
            if (env_fraction(i,k) > 0.0) then
              humidity_ratio(i,k) = max(qrf(i,k)*env_fraction(i,k)/env_qv(i,k), 1.0)
            else
              humidity_ratio(i,k) = -10.0
            end if
          else
            humidity_ratio(i,k) = 1.0
          end if
        end do
      end do
      
      !calc U_ca from compute_qs_a
      if (super_ice_opt < 1) then
        where (humidity_ratio > 0)
          U_ca = min(max(0.0,qin/(humidity_ratio*qs)),1.0)
        elsewhere
          U_ca = 0.0
        end where
        U01 = U_ca
      else
        do k=1, kdim
          do i=1, idim
            if (humidity_ratio(i,k) > 0) then
              U_ca(i,k) = max(0.0, qin(i,k)/(humidity_ratio(i,k)*qs(i,k)))
              if (tin(i,k) < tfreeze) then
                U01(i,k) = max(0.0, qin(i,k)/(humidity_ratio(i,k)*qsi(i,k)))
              else
                U01(i,k) = max(0.0, qin(i,k)/(humidity_ratio(i,k)*qsl(i,k)))
              end if
            else
              U_ca(i,k) = 0.0
              U01(i,k) = 0.0
            end if
          end do
        end do
      end if
      
      !GJF: from lscloud_driver.F90/impose_realizability
      !-----------------------------------------------------------------------
!    for the non-pdf scheme,  assure that cloud fraction is greater than 
!    qmin.  if it is not, set it to 0.0 and save the tendency and updated
!    value.
!-----------------------------------------------------------------------
      if (.not. pdfcld) then
        where (qa_in .le. qmin)
          SA = SA - qa_in
          qa_upd = 0.
        elsewhere
          qa_upd = qa_in     
        end where
        
!------------------------------------------------------------------------ 
!    define the max cloud area permitted (U01), while maintaining the 
!    grid-box-mean RH, under the assumption that the cloudy air is satur-
!    ated and the temperature inside and outside of the cloud are ~ the 
!    same. 
!------------------------------------------------------------------------ 
        U01 = min(U01, 1.)
        where (qa_upd .gt. U01)
          SA = SA + U01 - qa_upd
          qa_upd = U01      
        end where
      endif
      
      !-------------------------------------------------------------------------
!    define the conditions under which liquid and ice filling must be done,
!    for both the pdf and non-pdf cases.
!    for the non-pdf scheme, the filling requires that cloud liquid, cloud
!    ice, and cloud fraction are greater than qmin.
!    for pdf clouds, cloud fraction need not be considered since it is 
!    diagnosed later from the PDF cloud field.
!-------------------------------------------------------------------------
      if (.not. pdfcld) then
        if (do_liq_num) then
          ql_too_small = (ql_in .le. qmin .or.   &
                          qa_in .le. qmin .or.   &
                          qn_in .le. qmin)
        else
          ql_too_small = (ql_in .le. qmin .or.   &
                          qa_in .le. qmin)
        endif
        if (do_ice_num) then
          qi_too_small = (qi_in .le.  qmin .or.   &
                          qa_in .le. qmin .or.   &
                          qni_in .le. qmin)
        else
          qi_too_small = (qi_in .le.  qmin .or.   &
                          qa_in .le. qmin)
        endif
      else
        if (do_liq_num) then
          ql_too_small = (ql_in .le.  qmin  .or.   &
                          qn_in .le.  qmin)
        else
          ql_too_small = (ql_in .le.  qmin)
        endif
        if (do_ice_num) then
          qi_too_small = (qi_in .le.  qmin .or.   &
                          qni_in .le. qmin)
        else
          qi_too_small = (qi_in .le.  qmin )
        endif
      endif
      
!------------------------------------------------------------------------
!    call subroutine adjust_condensate to conservatively fill ql if needed.
!------------------------------------------------------------------------
      call adjust_condensate (ql_too_small, sl, sq, st, ql_in, hlv, cp, ql_upd)
      
!------------------------------------------------------------------------
!    adjust cloud droplet numbers as needed when those fields are being 
!    predicted. if droplet number is not being predicted, values were 
!    set at allocation.  
!------------------------------------------------------------------------
      if (do_liq_num) then
        call adjust_particle_number (ql_too_small, sn, qn_in, qn_upd)
      endif
      
!------------------------------------------------------------------------
!    call subroutine adjust_condensate to conservatively fill qi if needed.
!------------------------------------------------------------------------
      call adjust_condensate (qi_too_small, SI, SQ, ST, qi_in, HLS, cp, qi_upd)    
      
      !------------------------------------------------------------------------
!    adjust ice particle numbers as needed when those fields
!    are being predicted. 
!------------------------------------------------------------------------
      if (do_ice_num) then
        call adjust_particle_number (qi_too_small, SNi, qni_in, qni_upd)
      endif
      
      !GJF: need to calculate radturbten or use the already calculated process_split_cumulative_tendency_of_air_temperature; the latter also includes dT_dt due to GWD, which is not mentioned in the Tiedtke
      radturbten = dtdt_rad_turb_gwd
      
      !GJF: need to convert mass_number_concentration_of_cloud_liquid_water_particles_in_air_of_new_state (kg-1) to volume_number_concentration_of_cloud_liquid_water_particles_in_air_of_new_state (cm-3)
      airdens = pfull/(rdgas*tin)
      !GJF: convert mass number concentration (# per kg) to volume number concentration (# per cm3)
      if (do_aero_eros) then
        drop1 = qn_upd*airdens*1.0E-6
      endif
      
    end subroutine tiedtke_prog_clouds_pre_run
    
    subroutine adjust_condensate (   &
                mask, delta_cond, delta_q, delta_T, cond_in, lh, cp, cond_out)
      use machine  , only : kind_phys

!------------------------------------------------------------------------
!    subroutine adjust_condensate eliminates unacceptably small values of
!    condensate in a way which assures mass and enthalpy conservation.
!    tendency terms and output condensate fields are adjusted. 
!------------------------------------------------------------------------

      logical, dimension(:,:),  intent(in)    :: mask
      real(kind=kind_phys), dimension(:,:),     intent(in)    :: cond_in
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: delta_cond, delta_q, delta_T
      real(kind=kind_phys), dimension(:,:),     intent(out)   :: cond_out
      real(kind=kind_phys),                     intent(in)    :: lh, cp

!---------------------------------------------------------------------
!   local variables:

      integer :: idim, kdim    
      integer :: i, k

      idim = size(mask,1)
      kdim = size(mask,2)

!---------------------------------------------------------------------
      do k=1,kdim
        do i=1,idim
          if (mask(i,k)) then
            delta_cond(i,k) = delta_cond(i,k) -   cond_in(i,k)
            delta_q(i,k) = delta_q(i,k) + cond_in(i,k)
            delta_T(i,k) = delta_T(i,k) - lh*cond_in(i,k)/cp
            cond_out(i,k) = 0.
          else
            cond_out(i,k) = cond_in(i,k)
          endif
        end do
      end do

!------------------------------------------------------------------------

    end subroutine adjust_condensate
    
    subroutine adjust_particle_number (mask, delta_particles,   &
                                             particles_in, particles_out)
      use machine  , only : kind_phys

      logical, dimension(:,:),  intent(in)    :: mask
      real(kind=kind_phys), dimension(:,:),     intent(in)    :: particles_in
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: delta_particles
      real(kind=kind_phys), dimension(:,:),     intent(out)   :: particles_out
      
      !---------------------------------------------------------------------
      !   local variables:
      
            integer :: idim, kdim    
            integer :: i, k
      
            idim = size(mask,1)
            kdim = size(mask,2)
      
      !---------------------------------------------------------------------
            do k=1,kdim
              do i=1,idim
                if (mask(i,k)) then
                  delta_particles(i,k) = delta_particles(i,k) -  &
                                                   particles_in(i,k)
                  particles_out(i,k) = 0.
                else
                  particles_out(i,k) = particles_in(i,k)
                endif
              end do
            end do
      
      !------------------------------------------------------------------------      
    end subroutine adjust_particle_number
      
    end module tiedtke_prog_clouds_pre
