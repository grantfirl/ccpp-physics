!> \file tiedtke_prog_clouds_pre.F90
!!  Contains the preparation code for Tiedtke prognostic cloud scheme

!> This module contains the CCPP-compliant Tiedtke prognostic cloud scheme.
    module tiedtke_prog_clouds_pre

    contains

!> \section arg_table_tiedtke_prog_clouds_pre_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_pre_run.html
!!
    subroutine tiedtke_prog_clouds_pre_run (idim, kdim, pdf_org, pdfcld, super_ice_opt, hlv, hlf, rvgas, cp, tin, qin, pfull, cnvc, ql_in, qi_in, qa_in, qs, qsl, qsi, dqsdT, gamma, ql_upd, qi_upd, qa_upd, U_ca, U01, errmsg, errflg)

      use machine  , only : kind_phys
      use ccpp_saturation, only: get_qs
      
      implicit none
      
      integer,              intent(in)  :: idim, kdim
      
      logical,              intent(in)  :: pdf_org, pdfcld
      integer,              intent(in)  :: super_ice_opt
      
      real(kind=kind_phys), intent(in)  :: hlv, hlf, rvgas, cp
      real(kind=kind_phys), intent(in)  :: tin(:,:), qin(:,:), pfull(:,:), cnvc(:,:), ql_in(:,:), qi_in(:,:), qa_in(:,:)
      real(kind=kind_phys), intent(out) :: qs(:,:), qsl(:,:), qsi(:,:), dqsdT(:,:), gamma(:,:)
      real(kind=kind_phys), intent(out) :: ql_upd(:,:), qi_upd(:,:), qa_upd(:,:)
      real(kind=kind_phys), intent(out) :: U_ca(:,:), U01(:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      !temporary until set up in GFS_typedefs
      
      integer :: i,k
      integer :: alg_choice, alg_flatau_92
      
      real(kind=kind_phys) :: conv_frac_max, tfreeze
      real(kind=kind_phys) :: es(idim,kdim), esl(idim,kdim), esi(idim,kdim)
      real(kind=kind_phys) :: rh_wtd_conv_area(idim,kdim), convective_area(idim,kdim)
      real(kind=kind_phys) :: qrf(idim,kdim), env_qv(idim,kdim), env_fraction(idim,kdim), humidity_ratio(idim,kdim)
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      alg_choice = 1
      alg_flatau_92 = 1
      conv_frac_max = 1.0
      tfreeze = 273.16 !should be a physical constant coming in from host
            
      call get_qs(tin, pfull, alg_choice, alg_flatau_92, hlv, hlf, rvgas, cp, qs, qsl, qsi, es, esl, esi, dqsdT, gamma)
      
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
      
      !GJF: still need impose_realizability? very likely, yes
      ql_upd = ql_in
      qi_upd = qi_in
      qa_upd = qa_in
      
      !GJF: need to calculate radturbten
      !GJF: need to convert mass_number_concentration_of_cloud_liquid_water_particles_in_air_of_new_state (kg-1) to volume_number_concentration_of_cloud_liquid_water_particles_in_air_of_new_state (cm-3)
      
    end subroutine tiedtke_prog_clouds_pre_run
      
      
    end module tiedtke_prog_clouds_pre
