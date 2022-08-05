!> \file tiedtke_prog_clouds.F90
!!  Contains the Tiedtke prognostic cloud scheme

!> This module contains the CCPP-compliant Tiedtke prognostic cloud scheme.
      module tiedtke_prog_clouds

      contains

!> \section arg_table_tiedtke_prog_clouds_init Argument Table
!! \htmlinclude tiedtke_prog_clouds_init.html
!!
      subroutine tiedtke_prog_clouds_init (errmsg, errflg)
         use machine, only : kind_phys
         
         implicit none

         character(len=*),     intent(out) :: errmsg
         integer,              intent(out) :: errflg
         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0
         
         !from lscloud_driver.F90/lscloud_driver_init:
         
         !GJF - might not be needed if we use the same sat_vap_pres used in Thompson (Flatau curve)
         !call polysvp_init -> call CALL sat_vapor_pres_init from sat_vapor_pres_mod (from FMS), requires sat_vapor_pres_k_mod (also from FMS)
         
         !GJF - probably don't need this - can use CCN as other schemes in CCPP
         !call aerosol_cloud_init (aerosol_cloud.F90) -> aer_ccn_act_init -> aer_ccn_act_k_init; ice_nucl_wpdf_init -> ghquad
         
         !call tiedtke_macro_init

      end subroutine tiedtke_prog_clouds_init

!> \section arg_table_tiedtke_prog_clouds_run Argument Table
!! \htmlinclude tiedtke_prog_clouds_run.html
!!
      subroutine tiedtke_prog_clouds_run (idim, kdim, do_pdf_clouds, single_gaussian_pdf, &
          do_aero_eros, u00_profile, add_ahuco, eros_choice, include_neg_mc, super_ice_opt, &
          betaP, limit_conv_cloud_frac, use_qabar, pdf_org, &
          dmin, qmin, qthalfwidth, con_cp, con_g, con_hvap, con_hfus, con_rv, con_t0c, con_pi, &
          ae_ub, ae_lb, ae_N_ub, ae_N_lb, U00, eros_scale, eros_scale_c, eros_scale_t, &
          mc_thresh, diff_thresh, efact, dt, SA, SQ, ST, gamma, qs, dqsdT, U_ca, qsl, qsi, pfull, &
          phalf, mc_full, convective_humidity_area, diff_t, drop1, rh_crit, rh_crit_min, &
          tin, qin, ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, omega, radturbten, airdens, &
          hom, qvg, da_ls, dcond_ls, D_eros, errmsg, errflg)

      use machine, only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim, super_ice_opt, betaP
      logical, intent(in) :: do_pdf_clouds, single_gaussian_pdf, do_aero_eros, u00_profile, add_ahuco, eros_choice, include_neg_mc, limit_conv_cloud_frac, use_qabar, pdf_org
      real(kind=kind_phys), intent(in) :: dt, ae_ub, ae_lb, ae_N_ub, ae_N_lb, U00, eros_scale, eros_scale_c, eros_scale_t, mc_thresh, diff_thresh, efact, dmin, qmin, qthalfwidth
      real(kind=kind_phys), intent(in) :: con_cp, con_g, con_hvap, con_hfus, con_rv, con_t0c, con_pi
      real(kind=kind_phys), intent(in) :: gamma(:,:), qs(:,:), dqsdT(:,:), U_ca(:,:), qsl(:,:), qsi(:,:) !i,k dims
      real(kind=kind_phys), intent(in) :: pfull(:,:), phalf(:,:) !i,k, corresponds to prsl and prsi respectively
      real(kind=kind_phys), intent(in) :: mc_full(:,:) !total net convective mass flux on full levels kg m-2 s-1
      real(kind=kind_phys), intent(in) :: convective_humidity_area(:,:) !grid box area affected by the convective clouds
      real(kind=kind_phys), intent(in) :: diff_t(:,:) !eddy diffusivity for heat
      real(kind=kind_phys), intent(in) :: drop1(:,:) !cloud droplet number concentration (cm-3) (determined in aerosol_cloud.F90/determine_activated_aerosol)
      real(kind=kind_phys), intent(in) :: rh_crit(:,:) !critical relative humidity used in ice nuclei activation (determined in aerosol_cloud.F90/determine_activated_aerosol)
      real(kind=kind_phys), intent(in) :: rh_crit_min(:,:) !minimum critical relative humidity used in ice nuclei activation (determined in aerosol_cloud.F90/determine_activated_aerosol)
      real(kind=kind_phys), intent(in) :: tin(:,:) !i,k dims; this refers to the scheme's input state temperature (can either be process-split or time-split)
      real(kind=kind_phys), intent(in) :: qin(:,:) !i,k dims; this refers to the scheme's input state water vapor (can either be process-split or time-split)
      real(kind=kind_phys), intent(in) :: ql_in(:,:), qi_in(:,:), qa_in(:,:) !these vars are initialized every physics timestep to the advected tracers values of ql, qi in lscloud_driver.F90/lscloud_alloc
      real(kind=kind_phys), intent(in) :: ql_upd(:,:), qi_upd(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_alloc
      real(kind=kind_phys), intent(in) :: omega(:,:) !lagrangian_tendency_of_air_pressure
      real(kind=kind_phys), intent(in) :: radturbten(:,:) !sum of temperature tendency due to radiation and turbulence (K s-1) (interstitial -- reset every timestep)
      real(kind=kind_phys), intent(in) :: airdens(:,:) !air density (kg m-3) calculated in lscloud_driver
      real(kind=kind_phys), intent(in) :: hom(:,:) !serves as a flag representing homogeneous ice nucleation (calculated in ice_nucl.F90/ice_nucl_k)
      real(kind=kind_phys), intent(inout) :: qa_upd(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_alloc
      real(kind=kind_phys), intent(inout) :: qvg(:,:), da_ls(:,:), dcond_ls(:,:), D_eros(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_alloc
      !real(kind=kind_phys), intent(inout) :: dqa_dt_prod(:,:), dqa_dt_loss(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_driver
      
      real(kind=kind_phys), intent(in)    :: SQ(:,:) !this is set to zero in lscloud_driver.F90/lscloud_driver; change in specific humidity due to impose_realiziabilty/adjust_condensate; should be interstitial var
      real(kind=kind_phys), intent(in)    :: ST(:,:) !this is set to zero in lscloud_driver.F90/lscloud_driver; change in temperature due to impose_realiziabilty/adjust_condensate; should be interstitial var
      real(kind=kind_phys), intent(inout) :: SA(:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
!------------------------------------------------------------------------
!----local variables----
!
!       dtinv          inverse of physics timestep
!       mdum           factor used in aero_eros calculation
!       bdum           factor used in aero_eros calculation
!       edum           factor used in aero_eros calculation
!       U00p           critical relative humidity fraction which may be 
!                      a function of pressure 
!       U00pr          gridbox mean rh needed for condensation, after
!                      accounting for convective cloud area
!       erosion_scale  cloud erosion scale
!       i,j,k          do-loop indices
!-------------------------------------------------------------------------
      integer :: i,k
      real(kind=kind_phys) :: dt_inv, hls
      real(kind=kind_phys), dimension(idim,kdim) :: mdum,bdum,edum,U00p,U00pr,erosion_scale,delta_cf

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      !GJF need to call compute_qs_a to calculate qs,gamma,dqsdT,U_ca,qsl,qsi
      
      dt_inv = 1.0/dt
      
      hls = con_hvap + con_hfus
      
!------------------------------------------------------------------------
!    process the non-convective condensation for pdf clouds. 
!                                                                      
!                NON-CONVECTIVE CONDENSATION                          
!                                                                      
!                STATISTICAL CLOUD FRACTION                            
!                                                                      
!------------------------------------------------------------------------
      if (do_pdf_clouds) then
        if (single_gaussian_pdf) then      
          call tiedtke_macro_Single_Gaussian_pdf(idim, kdim, pdf_org, qmin, qthalfwidth, con_pi, & 
            gamma, qs, qin, &
            ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, SQ, qvg, da_ls, delta_cf, dcond_ls, &
            SA)
        else
          call tiedtke_macro_pdf (     &
            idim, kdim, betaP, pdf_org, qmin, qthalfwidth, gamma, qs, qin, ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, &
            SQ, qvg, D_eros, da_ls, delta_cf, dcond_ls, SA)
        endif
      else  !(do_pdf_clouds)
!------------------------------------------------------------------------
!    process the non-convective condensation for non-pdf clouds. some
!    additional calculations are needed.
!
!                 NON-CONVECTIVE CONDENSATION                          
!                                                                     
!                 TIEDTKE (1993) CLOUD FRACTION                    
!                                                                    
!       ANALYTIC INTEGRATION OF SATURATED VOLUME FRACTION EQUATION     
!
!       Do non-convective condensation following Tiedtke, pages 3044-5.
!       In this formulation stratiform clouds are only formed/destroyed 
!       when there is upward or downward motion to support/destroy it. 
!-------------------------------------------------------------------------

!-----------------------------------------------------------------------
!    calculate aerosol erosion term which optionally may be used.
!-----------------------------------------------------------------------
        if (do_aero_eros ) then 
          mdum = (ae_ub - ae_lb)/(ae_N_ub - ae_N_lb)
          bdum = ae_lb - mdum*ae_N_lb
          edum = mdum*drop1 + bdum 
        else
          edum = 1.
        endif

!------------------------------------------------------------------------
!    compute pressure dependent critical relative humidity for onset of
!    condensation (U00p) following ECMWF formula if desired. otherwise,
!    it is given by the nml variable u00.
!------------------------------------------------------------------------
        U00p = U00
        if (u00_profile) then
          DO k=1,kdim
            where (pfull(:,k).gt.   &
                          0.8*phalf(:,KDIM+1)) 
              U00p(:,k) = U00 + (1.- U00)* &
                         (((pfull(:,k) -   &
                            (0.8*phalf(:,KDIM+1)))/  &
                               (0.2*phalf(:,KDIM+1)) )**2.)
            end where
          END DO
        endif       

!------------------------------------------------------------------------
!    modify u00p to account for humidity in convective system.  see 
!    "Tiedtke u00 adjustment" notes, 10/22/02 -- ljd
!    u00p = critical rh for condensation in the grid box, accounting for 
!    the fact that ahuco of the box is saturated (has convective cloud) 
!    better written : ahuco*100% + (1-ahuco)*u00p = gridboxmean rh
!------------------------------------------------------------------------
        IF (add_ahuco) THEN
          u00pr = u00p + (1. - u00p)*convective_humidity_area
          u00p = u00pr
        END IF

!-----------------------------------------------------------------------
!    Theory for eros_scale
!
!    If eros_choice equals false, then a single erosion time scale
!    is used in all conditions (eros_scale).  If eros_choice equals
!    true then it is assumed that the timescale for turbulent 
!    evaporation is a function of the conditions in the grid box.  
!    Specifically, if the flow is highly turbulent then the scale is 
!    short, and eros_scale is large.  Likewise if convection is 
!    occurring, then it is assumed that the erosion term is larger 
!    than backround conditions. 
!
!    Here are the typical values for the timescales and the 
!    switches used (subject to changes via namelist):
!
!         Mixing type      eros_scale (sec-1)          Indicator
!       ----------------   ------------------     --------------------
!
!       Background            1.e-06              always present
!       Convective layers     5.e-06              Mc > Mc_thresh
!       Turbulent  layers     5.e-05              diff_t > diff_thresh
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    define array with background erosion scale.
!-----------------------------------------------------------------------
        erosion_scale =     eros_scale
               
!-----------------------------------------------------------------------
!    Do enhanced erosion in convective or turbulent layers?
!
!                   IMPORTANT NOTE
!                
!    note that convection is considered first, so that if turbulence and 
!    convection occur in the same layer, the erosion rate for turbulence 
!    is selected.                
!-----------------------------------------------------------------------
        if (eros_choice) then
          do k=1,kdim
              do i=1,idim

!-----------------------------------------------------------------------
!    Enhanced erosion in convective layers
!-----------------------------------------------------------------------
                if (include_neg_mc) then
                  if (abs(mc_full(i,k)) .gt. mc_thresh) then 
                    erosion_scale(i,k) = eros_scale_c
                  endif
                else
                  if (mc_full(i,k) .gt. mc_thresh) then
                    erosion_scale(i,k) = eros_scale_c
                  endif
                endif
                if ((diff_t(i,K) .gt. diff_thresh) .or.&
                    (diff_t(i,min(k+1,KDIM)) .gt.  &
                                                       diff_thresh) ) then
                  erosion_scale(i,K) = eros_scale_t
                endif 
              END DO
          END DO
        end if   !(eros_choice)        

!------------------------------------------------------------------------
!    enhance erosion scale if that option has been selected.
!------------------------------------------------------------------------
        IF (do_aero_eros) erosion_scale = edum*erosion_scale

!------------------------------------------------------------------------
!     determine condensation for different supersaturation assumptions.
!------------------------------------------------------------------------
        IF (super_ice_opt .LT. 1 ) THEN

!------------------------------------------------------------------------
!    if supersaturation is not to be allowed in the clear-sky portion
!    of the grid box, call subroutine tiedtke_macro_nopdf_nosuper to 
!    calculate the condensation.
!------------------------------------------------------------------------
          call tiedtke_macro_nopdf_nosuper (    &
               idim, kdim, limit_conv_cloud_frac, use_qabar, mc_full, convective_humidity_area, omega, radturbten, pfull, airdens, dqsdT, U_ca, gamma, qs, &
               qa_upd, ql_upd, qi_upd, SQ, da_ls, D_eros, dcond_ls, delta_cf, &
               edum, SA, U00p, erosion_scale, dt, dmin, qmin, efact, con_cp, con_g)
        else

!------------------------------------------------------------------------
!    if supersaturation is to be allowed in the clear-sky portion of the
!    gridbox, call subroutine tiedtke_macro_nopdf_super to calculate the 
!    condensation.
!------------------------------------------------------------------------
          call tiedtke_macro_nopdf_super  (     &
               idim, kdim, limit_conv_cloud_frac, mc_full, convective_humidity_area, omega, radturbten, pfull, tin, qin, airdens, dqsdT, gamma, qs, qsl, qsi, rh_crit, rh_crit_min, &
               qa_upd, ql_upd, qi_upd, ST, SQ, da_ls, D_eros, dcond_ls, delta_cf, hom,   &
               edum, SA, U00p, erosion_scale, dt, dmin, qmin, efact, con_cp, con_g, con_hvap, hls, con_rv, con_t0c)
        endif
      endif

      end subroutine tiedtke_prog_clouds_run

      SUBROUTINE tiedtke_macro_Single_Gaussian_pdf (idim, kdim, pdf_org, qmin, qthalfwidth, pi, gamma, qs,   &
                              qin, &
                              ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, SQ, qvg, da_ls, delta_cf, dcond_ls, &
                              SA)  
        use machine, only : kind_phys
        implicit none
      !----------------------------------------------------------------------!
      !                                                                      !
      !                     STATISTICAL CLOUD SCHEME                         !
      !                                                                      !
      !-----------------------------------------------------------------------

      integer,                             intent(in)    :: idim, kdim
      logical,              intent(in)                    :: pdf_org
      real(kind=kind_phys), intent(in)                    :: qmin, qthalfwidth, pi
      real(kind=kind_phys), intent(in),    dimension(:,:) :: gamma, qs !i,k dims
      real(kind=kind_phys), intent(in),    dimension(:,:) :: qin !i,k dims
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_in, qi_in, qa_in
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_upd, qi_upd
      real(kind=kind_phys), intent(inout), dimension(:,:) :: qa_upd
      real(kind=kind_phys), intent(inout), dimension(:,:) :: qvg, da_ls, delta_cf, dcond_ls
      real, dimension(idim,kdim),   intent(in)    :: SQ !, GJF: ST is not used
      real, dimension(idim,kdim),   intent(inout) :: SA
      
      !-----------------------------------------------------------------------
      !-----local variables---------------------------------------------------
       
      !-----------------------------------------------------------------------
      !                STATISTICAL CLOUD SCHEME VARIABLES
      !
      !
      !       qcg            total condensate diagnosed        kg condensate /
      !                      from Single-Gaussian PDF          kg air        
      !
      !       qag            cloud fraction diagnosed          dimensionless 
      !                      from Single-Gaussian PDF                                  
      !
      !
      !       qtbar          total water specific humidity     kg water /
      !                      which is equal to the sum of      kg air
      !                      liquid water, ice water, and
      !                      water vapor
      !
      !       s              s = (qt -qs)/(1 + L/Cp * dqsdT)   kg condensate/kg air
      !                     
      !       stddev_s       standard deviation of s           kg condensate/kg air
      !
      !       if ignoring liquid potential temperature variance,  
      !       and correlation between liquid potential temperature
      !       and total water, 
      !
      !       stddev_s = stddev_qt /(1 + L/Cp * dqsdT)
      !                = stddev_qt /(1 + gamma ) 


            real, dimension(idim,kdim)    :: qa1, qa0, qcg, qag
            real, dimension(idim,kdim)    :: qta, qtqsa
       
      !---> h1g, 2015-07-22
           real, dimension(idim,kdim)    ::  s,   stdev_s
           real                               ::  s_mellor_tol, zeta
           integer                            ::  i, j, k
      !-----------------------------------------------------------------------
      !     Compute total water and excess saturation
      !-----------------------------------------------------------------------
            if (pdf_org) then
      !RSH  change to qin ??
      !       qta(:,:,:) = max(qmin, qv_in +  &
              qta(:,:) = max(qmin, qin +  &
                                           ql_in + qi_in)
            else
      !       qta(:,:,:) = max(qmin, qv_in + SQ +  &
              qta(:,:) = max(qmin, qin + SQ +  &
                                           ql_upd + qi_upd)
            endif
            qtqsa(:,:) = qta(:,:) - qs      

            s =  qtqsa(:,:) / (1.0 +  gamma )

      !---> h1g, 2015-07-23
      ! 99.7% data are within mean +- 3 stdev 
            stdev_s =     qthalfwidth * qta * 0.333
            stdev_s = stdev_s / ( 1.0 +  gamma ) 

            s_mellor_tol = 1.e-8
            do k= 1,kdim
              do i = 1,idim
                if( stdev_s(i,k)  > s_mellor_tol ) then
                  zeta = s(i,k)/stdev_s(i,k)
                  qag(i,k)  = 0.5 * ( 1.0 + erf( zeta /sqrt(2.0) ) )
                  qcg(i,k)  = s(i,k) * qag(i,k)  &
                              + stdev_s(i,k)/sqrt(2.0*pi) * exp( -0.5 *zeta*zeta )
                else
                  if ( s(i,k) < 0.0 ) then
                    qag(i,k)  = 0.0
                    qcg(i,k)  = 0.0
                  else
                    qag(i,k)  = 1.0
                    qcg(i,k)  = s(i,k)
                  endif  
                endif

      !---> h1g, 2015-07-23
      !     qtbar:     grid-average total specific humidity
      !     qv,clr:    clear-sky  water vapor specific humidity
      !     qc,cld:    in-cloud  cloud condensate
      !     qs:        satuation vapor specific humidity
      !
      !     qtbar = (1.-CF) * qv,clr + CF * ( qc,cld + qs )
      !     qv,clr = (qtbar - CF * ( qc,cld + qs ) )/(1.-CF)      
      !<--- h1g, 2015-07-23
                qvg(i,k) = (qta(i,k) - qcg(i,k) - qag(i,k)*qs(i,k))  &
                                              /(max( (1.0 - qag(i,k)), qmin) )
              enddo
            enddo        

            !do adjustment of cloud fraction
            qa0 = qa_in
            qa1 = qag

            !set total tendency term and update cloud fraction    
            SA   = (SA   + qa1) - qa0
            qa_upd     = qa1

            !define da_ls and tmp5 needed when do_liq_num = .true. (cjg)
            da_ls = max(qa1-qa0,0.)
            delta_cf = max(qa1-qa0,0.)

            !compute large-scale condensation / evaporation
            dcond_ls = qcg -    &
                                    (ql_upd + qi_upd)


      !--->h1g, the following is inherited from nc_cond_pdf, 2015-07-22
            if ( .not. pdf_org ) then
            !!!! INVESTIGATE!!!
            ! make sure super/subsat is not created here
            ! this is different from the original PDF assumption
            ! (as is saturation adjustment) 

      !RSH    dcond_ls = MAX( ((qv_in + SQ -  &
              dcond_ls = MAX( ((qin      + SQ -  &
                      qs)/(1.+gamma)),    &
                                                     dcond_ls)
        
            end if
      !<--- h1g, 2015-07-22

      !-----------------------------------------------------------------------
      !     Diagnostics
      !-----------------------------------------------------------------------

            !if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
            !    diag_4d(:,:,:,diag_pt%qadt_lsform ) =  max(qa1 - qa0, 0.)*  &
            !                                               inv_dtcloud 
            !end if
            !if (diag_id%qadt_lsdiss + diag_id%qa_lsdiss_col > 0) then
            !    diag_4d(:,:,:,diag_pt%qadt_lsdiss ) =  max(qa0 - qa1, 0.)* &
            !                                               inv_dtcloud
            !end if
      !------------------------------------------------------------------------

      end SUBROUTINE tiedtke_macro_Single_Gaussian_pdf
      
      SUBROUTINE tiedtke_macro_pdf (  &
                   idim, kdim, betaP, pdf_org, qmin, qthalfwidth, gamma, qs, qin, ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, &
                   SQ, qvg, D_eros, da_ls, delta_cf, dcond_ls, SA)
        use machine, only : kind_phys
        implicit none
      !----------------------------------------------------------------------!
      !                                                                      !
      !                     STATISTICAL CLOUD SCHEME                         !
      !                                                                      !
      !-----------------------------------------------------------------------

      integer,                           intent(in )   :: idim, kdim, betaP
      logical,              intent(in)                    :: pdf_org
      real(kind=kind_phys), intent(in)                    :: qmin, qthalfwidth
      real(kind=kind_phys), intent(in),    dimension(:,:) :: gamma, qs
      real(kind=kind_phys), intent(in),    dimension(:,:) :: qin
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_in, qi_in, qa_in
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_upd, qi_upd
      real(kind=kind_phys), intent(inout), dimension(:,:) :: qa_upd
      real(kind=kind_phys), intent(inout), dimension(:,:) :: qvg, D_eros, da_ls, delta_cf, dcond_ls
      !real(kind=kind_phys), intent(inout), dimension(:,:) :: dqa_dt_prod, dqa_dt_loss
      real, dimension(idim,kdim),   intent(in )   :: SQ
      real, dimension(idim,kdim),   intent(inout) :: SA
      
      !-----------------------------------------------------------------------
      !-----local variables
       
      !-----------------------------------------------------------------------
      !                STATISTICAL CLOUD SCHEME VARIABLES
      !
      !
      !       qcg            equilibrium value of cloud      kg condensate /
      !                      condensate that PDF clouds      kg air
      !                      wants           
      !
      !       qag            equilibrium value of cloud      dimensionless 
      !                      fraction for statistical 
      !                      cloud scheme           
      !
      !       qs_norm        the difference between the      dimensionless
      !                      saturation specific humidity    
      !                      and qtmin normalized by deltaQ
      !
      !       icbp           the value of the incomplete     dimensionless
      !                      beta function evaluated with
      !                      x=qs_norm, p=betaP, and q=betaP
      !
      !       icbp1          the value of the incomplete     dimensionless
      !                      beta function evaluated with                  
      !                      x=qs_norm, p=betaP+1, and 
      !                      q=betaP
      !
      !       qtbar          total water specific humidity   kg water /
      !                      which is equal to the sum of    kg air
      !                      liquid water, ice water, and
      !                      water vapor
      !
      !       deltaQ         the width of the total water    kg water /
      !                      subgrid distribution (= qtmax   kg air
      !                      minus qtmin)
      !
      !       qtmin          the minimum value to the total  kg water /
      !                      sub-grid scale distribution     kg air
      !
            real, dimension(idim, kdim)   :: qa1, qa0, qcg, qag
            real, dimension(idim, kdim)   :: qta, qtqsa
            real, dimension(idim)         :: qagtmp, qcgtmp, qvgtmp, qtbar,&
                                                  deltaQ, qtmin, qs_norm     
            real                               :: icbp, icbp1
            integer                            :: i,j,k

       
      !-----------------------------------------------------------------------
      !    Compute total water and excess saturation.
      !-----------------------------------------------------------------------
            if (pdf_org) THEN
              qta(:,:) = max(qmin, qin +  &
                                        ql_in + qi_in)
            else
              qta(:,:) = max(qmin, qin + SQ +  &
                                   ql_upd + qi_upd)
            endif 

            qtqsa(:,:) = qta(:,:) - qs
              

#ifdef SKIP
            qcg(:,:) = 0.
            qvg(:,:) = 0.
            qag(:,:) = 0.
#endif
             
      !-----------------------------------------------------------------------
      !    set Tiedtke erosion term to zero.
      !-----------------------------------------------------------------------
            D_eros = 0.
              
      !-----------------------------------------------------------------------
      !    compute pdf cloud fraction and condensate
      ! 
      !    Note that the SYMMETRIC beta distribution is used here.
      !
      !
      !    Initialize grid-box mean values of cloud fraction (qag),
      !    cloud condensate(qcg), and clear sky water vapor (qvg)
      !-----------------------------------------------------------------------
           ks_loop: DO k=1,kdim
             !! qcg(:,k) = 0.
              ! qvg(:,k) = 0.
              !!qag (:,k)= 0.
              
#ifdef SKIP
              !Create loop over sub-levels within a grid box
         
              sublevel_loop: do ns = 1, nsublevels
              
                   !calculate normalized vertical level
                   ! 0. = top of gridbox
                   ! 1. = bottom of gridbox
              
                   pnorm =  (real(ns) - 0.5 )/real(nsublevels)
              
                   !First step is to calculating the minimum (qtmin)
                   !of the total water distribution and 
                   !the width of the qt distribution (deltaQ)
                   !
                   !For diagnostic variance this is set to (1.-qthalfwidth)*qtbar
                   !and 2*qthalfwidth*qtbar, respectively, where qtbar is the
                   !mean total water in the grid box.        
                   !
                   !

                   qtbar = qta4(2,:,:,k) + pnorm*((qta4(3,:,:,k) -   &
                                      qta4(2,:,:,k)) + qta4(4,:,:,k)*(1-pnorm) )
                   
                   qtbar  = max(qmin ,qtbar )
                   deltaQ  = 2.*qthalfwidth *qtbar
                   qtmin = (1.- qthalfwidth )*qtbar 
#endif
                   qtbar = qta(:,k)
                   qtbar  = max(qmin ,qtbar )
                   deltaQ  = 2.*qthalfwidth *qtbar
                   qtmin = (1.-qthalfwidth )*qtbar
              
                   !From this the variable normalized saturation specific
                   !humidity qs_norm is calculated.
                   !
                   !  qs_norm = (qs(Tl) - qtmin)/(qtmax-qtmin)
                   !
                   !          = 0.5  - (qtbar - qs(Tl))/deltaQ
                   !
                   !Note that if qs_norm > 1., the grid box is fully clear.
                   !If qs_norm < 0., the grid box is fully cloudy.
              
      !            qs_norm = qtqsa4(2,:,:,k)+  &
      !                      pnorm*( (qtqsa4(3,:,:,k)-qtqsa4(2,:,:,k)) + &
      !                      qtqsa4(4,:,:,k)*(1-pnorm) )

                   qs_norm = qtqsa(:,k)
                   qs_norm = 0.5 - ( qs_norm/deltaQ )
                   
                   !Calculation of cloud fraction (qagtmp), cloud condensate 
                   !(qcgtmp), and water vapor in clear air part of the grid 
                   !box (qvgtmp)
                   !
                   !Formulas (from Tompkins, and personal derivations):
                   !
                   !  Define icbp  = incomplete_beta(qs_norm,p,q)
                   !         icbp1 = incomplete_beta(qs_norm,p+1,q)
                   !
                   !  qagtmp = 1. - icbp
                   !
                   !  qcgtmp = aThermo * {  (qtbar-qtmin)*(1.-icbp1) - 
                   !                       qs_norm*deltaQ*(1.-icbp ) }
                   !
                   !
                   !  qvgtmp = qtmin + (p/(p+q))*(icbp1/icbp)*deltaQ
                   !
                   !  
                   ! where aThermo = 1./(1.+(L/cp)*dqsdT)
                   !
                   ! note that in the qvg formula below the factor of 0.5
                   ! is equal to (p/(p+q)).
                   !

                   
                   
                   do i = 1,idim
                    
                    !GJF: incomplete_beta requires INPUT/BetaDistributionTable.txt (which I can't get from FTP!)
                   ! if (qs_norm(i).le.1.) then
                   ! 
                   !     icbp = incomplete_beta(max(0.,qs_norm(i)), &
                   !                          p = betaP    , q =     betaP)
                   !     icbp1= incomplete_beta(max(0.,qs_norm(i)), &
                   !                          p = betaP + 1, q =     betaP)
                   !     qagtmp(i) = 1.-icbp
                   !     qcgtmp(i) = (qtbar(i)-qtmin(i))*(1.-icbp1)&
                   !                   - qs_norm(i)*deltaQ(i)*(1.-icbp)    
                   !     qcgtmp(i) = qcgtmp(i)/   &
                   !                              (1.+gamma(i,k))
                   !     qvgtmp(i) = qtmin(i,j) + &
                   !                   0.5*(icbp1/max(icbp, qmin))*deltaQ(i)
                   ! 
                   !     !bound very very small cloud fractions which may
                   !     !cause negative cloud condensates due to roundoff 
                   !     !errors or similar errors in the beta table lookup.
                   !     if((qagtmp(i).lt.0.).or.(qcgtmp(i).le.0.))then
                   !          qagtmp(i) = 0.
                   !          qcgtmp(i) = 0.
                   !          qvgtmp(i) = qtbar(i)
                   !     end if
                   ! 
                   ! else             
                       qagtmp(i) = 0.
                       qcgtmp(i) = 0.
                       qvgtmp(i) = qtbar(i)             
                   ! end if
                   
                   enddo
                   
                    
#ifdef SKIP
                   !sum vertically
                   !
                   !note special averaging of clear-sky water vapor
                   !this is weighting clear-sky relative humidity by the 
                   !clear-sky fraction
               
                   qag(:,k) = qag(:,k) + qagtmp
                   qcg(:,k) = qcg(:,k) + qcgtmp
                   qvg(:,k) =    &
                           qvg(:,k)+ &
                                     (1.-qagtmp)*min(max(qvgtmp/max(qmin, &
                                         (qtbar+((qs_norm-0.5)*deltaQ))),0.),1.)
                   
              enddo sublevel_loop!for number of sublevels loop

              


              !compute grid-box average cloud fraction, cloud condensate
              !and water vapor
              
              if (nsublevels.gt.1) then
                   qag (:,k)= qag(:,k) / real(nsublevels)
                   qcg(:,k) = qcg(:,k) / real(nsublevels)
                   
                   !note special averaging of clear-sky water vapor
                  
                   
                   do id = 1,idim
                        if ((1.-qag(id,k)).gt. qmin) then
                          qvg(id,k) =  &
                             qvg(id,k)/real(nsublevels)/&
                                                              (1.-qag(id,k))
                          qvg(id,k) =  &
                                    qvg(id,k)*   &
                                                       qs(id,k)
                        else
                          qvg(id,k) =   &
                                                      qs(id,k)
                        end if
                   enddo
                   
                
                   
      !       else
#endif
                   qag(:,k) = qagtmp
                   qcg(:,k) = qcgtmp
                   qvg(:,k) = qvgtmp
      !        end if
           
             
         END DO ks_loop
              
              !do adjustment of cloud fraction
              qa0 = qa_in
              qa1 = qag

              !set total tendency term and update cloud fraction    
              SA   = (SA   + qa1) - qa0
              qa_upd     = qa1


              !define da_ls and tmp5 needed when do_liq_num = .true. (cjg)
              da_ls = max(qa1-qa0,0.)
              delta_cf = max(qa1-qa0,0.)

              !compute large-scale condensation / evaporation
              dcond_ls = qcg -    &
                                    (ql_upd + qi_upd)



         IF ( .NOT. pdf_org ) THEN
         !!!! INVESTIGATE!!!
         !make sure super/subsat is not created here
         ! this is different from the original PDF assumption
         !  (as is saturation adjustment) 

           dcond_ls = MAX( ((qin + SQ -  &
                      qs)/(1.+gamma)),    &
                                                     dcond_ls)
             
            
         endif 

      !-----------------------------------------------------------------------
      !     Diagnostics
      !-----------------------------------------------------------------------

            !diag_dt_prod = max(qa1 - qa0, 0.)*inv_dtcloud
            !diag_dt_loss = max(qa0 - qa1, 0.)*inv_dtcloud
            
      !------------------------------------------------------------------------
       


      end SUBROUTINE tiedtke_macro_pdf
      
      SUBROUTINE tiedtke_macro_nopdf_nosuper (    &
                  idim, kdim, limit_conv_cloud_frac, use_qabar, mc_full, convective_humidity_area, omega, radturbten, pfull, airdens, dqsdT, U_ca, gamma, qs, &
                  qa_upd, ql_upd, qi_upd, SQ, da_ls, D_eros, dcond_ls, delta_cf, &
                  edum, SA, U00p, erosion_scale, dtcloud, dmin, qmin, efact, cp_air, grav)
        use machine, only : kind_phys
        implicit none
      !------------------------------------------------------------------------
      !   subroutine tiedtke_macro_nopdf_nosuper calculates non-convective
      !   condensation while not allowing supersaturation in the clear-sky 
      !   portion of grid boxes (as in original Tiedtke scheme).
      !------------------------------------------------------------------------

      INTEGER,                         INTENT(IN )   :: idim, kdim
      logical,              intent(in) :: limit_conv_cloud_frac, use_qabar
      real(kind=kind_phys), intent(in) :: cp_air, grav, dmin, qmin, efact, dtcloud
      real(kind=kind_phys), intent(in) :: mc_full(:,:) !total net convective mass flux on full levels kg m-2 s-1
      real(kind=kind_phys), intent(in) :: convective_humidity_area(:,:) !grid box area affected by the convective clouds
      real(kind=kind_phys), intent(in) :: omega(:,:) !lagrangian_tendency_of_air_pressure
      real(kind=kind_phys), intent(in) :: radturbten(:,:) !sum of temperature tendency due to radiation and turbulence (K s-1) (interstitial -- reset every timestep)
      real(kind=kind_phys), intent(in) :: pfull(:,:)
      real(kind=kind_phys), intent(in) :: airdens(:,:) !air density (kg m-3) calculated in lscloud_driver
      real(kind=kind_phys), intent(in) :: dqsdT(:,:)
      real(kind=kind_phys), intent(in) :: gamma(:,:), qs(:,:)
      real(kind=kind_phys), intent(in) :: U_ca(:,:) !grid box relative humidity (fraction)
      real(kind=kind_phys), intent(inout) :: qa_upd(:,:)
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_upd, qi_upd
      real(kind=kind_phys), intent(inout) :: da_ls(:,:)
      real(kind=kind_phys), intent(inout) :: D_eros(:,:)
      real(kind=kind_phys), intent(inout) :: dcond_ls(:,:)
      real(kind=kind_phys), intent(inout) :: delta_cf(:,:)
      REAL, dimension(idim,kdim), INTENT(IN)    :: SQ, edum, U00p,  &
                                                        erosion_scale
      REAL, dimension(idim,kdim), INTENT(INOUT) :: SA
      
      !-------------------------------------------------------------------------
      !----local variables----

      !-------------------------------------------------------------------------
      !       dqs_ls         change in qs due to large 
      !                      scale processes
      !       A_dt           product of A and dtcloud in     dimensionless in 
      !                      in the analytic integration     qa integration
      !                      of the qa equation, or C and
      !                      dtcloud in the analytic         kg condensate/
      !                      integration of the ql and qi    kg air in ql or 
      !                      equations.                      qi integration
      !
      !       B_dt           product of B and dtcloud in     dimensionless in
      !                      in the analytic integration     qa, ql, and qi
      !                      of the qa equation, or D and    integration
      !                      dtcloud in the analytic         
      !                      integration of the ql and qi    
      !                      equations.                      
      !
      !       qa0            value of cloud fraction or      dimensionless or
      !                      cloud condensate at the         kg condensate /
      !                      initial time                    kg air
      !        
      !       qa1            value of cloud fraction or      dimensionless or
      !                      cloud condensate at the final   kg condensate /
      !                      time                            kg air
      !
      !       qaeq           equilibrium value of cloud      dimensionless or
      !                      fraction or cloud condensate    kg condensate /
      !                      that the analytic integration   kg air
      !                      approaches                   
      !       qabar          mean value of cloud fraction    dimensionless or
      !                      or cloud condensate over the    kg condensate /
      !                      t0 to t0 + dtcloud interval     kg air
      !
      !------------------------------------------------------------------------
           
            real, dimension(idim, kdim)   :: dqs_ls
            real, dimension(idim, kdim)   :: A_dt, B_dt
            real, dimension(idim, kdim)   :: qa1, qa0, qabar
            real, dimension(idim, kdim)   :: qaeq
            real, dimension(idim, kdim)   :: tmp1
            REAL                               :: eslt, qvs, qs_d, qvi, esit, ul
            REAL                               :: ttmp, qtmp, qs_t, dqsdT1,  &
                                                  qvmax, esat0, gamma1,  &
                                                  tmp1s, qs_l, qs_i
            INTEGER                            :: ns, id
            INTEGER                            :: i,k

      !-----------------------------------------------------------------------
      !
      !       The first step is to compute the change in qs due to large-
      !       scale processes, dqs_ls.   In Tiedtke, it has contributions from 
      !       large-scale uplift, convection induced compensating subsidence,
      !       turbulence cooling and radiative cooling.  dqs_ls has the form:
      !
      !               (((omega+ grav*Mc)/airdens/cp)+radturbten)*dqsdT*dtcloud
      !   (6) dqs_ls= --------------------------------------------------------
      !                  1.  +   ( qa +  (da_ls/2.) ) * gamma
      !
      !       Here da_ls is the increase in cloud fraction due to non-
      !       convective processes.  Because this increase is also a function
      !       of dqs_ls, a quadratic equation must be solved for dqs_ls in
      !       the case that da_ls is not equal to zero.
      !
      !------------------------------------------------------------------------
            do k=1,kdim
              do i=1,idim
                dqs_ls(i,k) = ((                    &
                  (omega(i,k) + grav*mc_full(i,k))/  &
                                    airdens(i,k)/cp_air)  +   &
                             radturbten(i,k))* &
                                             dtcloud*dqsdT(i,k)
                if (dqs_ls(i,k) .le. 0. .and.    &
                           U_ca(i,k) .ge. U00p(i,k) .and.   &
                                    qa_upd(i,k) .lt. 1.)  then
                  tmp1(i,k) = sqrt( (1. + qa_upd(i,k)*  &
                                           gamma(i,k))**2. -  &
                                      (1. - qa_upd(i,k))*  &
                                         (1. - qa_upd(i,k))*&
                                     gamma(i,k)*dqs_ls(i,k)/&
                                              qs(i,k)/      &
                        max(1.-U_ca(i,k), qmin) ) -   &
                            (1. + qa_upd(i,k)*  &
                                                   gamma(i,k))
                  tmp1(i,k) = -1.*tmp1(i,k)/((1. -   &
                                           qa_upd(i,k))*   &
                                   (1. - qa_upd(i,k))*  &
                                                gamma(i,k)/&
                                    qs(i,k)/   &
                                max(1. - U_ca(i,k), qmin)/2.)
                  dqs_ls(i,k) = min(tmp1(i,k),dqs_ls(i,k)/(1. +   &
                                  0.5*(1. + qa_upd(i,k))*  &
                                                 gamma(i,k)))
                else
                  dqs_ls(i,k) = dqs_ls(i,k)/(1. +   &
                       qa_upd(i,k)*gamma(i,k))
                endif
              end do
            end do

      !------------------------------------------------------------------------
      !       The next step is to compute the change in saturated volume
      !       fraction due to non-convective condensation, da_ls.   This 
      !       occurs in two conditions:
      !
      !       (a)  dqs_ls < 0. and U00 < U < 1., where U00 is the threshold
      !            relative humidity for non-convective condensation. Note 
      !            that if U is greater than or equal to 1., ideally qa = 1,
      !            and da_ls = 0.  However this may not be the case for 
      !            numerical reasons so this must be assured after analytic 
      !            integration of the qa equation.
      !
      !            For these cases the change in saturated volume fraction is:
      !
      !   (7)      da_ls = - (1.-qa)*(1.-qa)*dqs_ls/2./qs/(1.-U)
      !
      !            This formula arises from the assumption that vapor is uni-
      !            formly distributed in the range [qv_clr - (qs - qv_clr),qs]
      !            where qv_clr is the amount of vapor in the unsaturated 
      !            volume and is given from the following equation:
      !
      !   (8)      qv  =   qa * qs      +   (1.-qa) * qv_clr
      !          
      !            Implicit in equation (7) is the following assumption:
      !            As qsat changes, the distribution of qv+ql+qi 
      !            remains constant.  That is as qsat rises, portions where
      !            qv+ql+qi > qsat+dqsat remain saturated.  This can only
      !            occur if it is assumed that ql+qi evaporate-sublimate or
      !            condense-deposit to keep qv = qsat. 
      !
      !       (b)  dqs_ls > 0.  Ideally some portion of the cloud should
      !            evaporate however this is not accounted for at present.
      !            
      !------------------------------------------------------------------------
            do k=1,kdim
              do i=1,idim
                if ((dqs_ls(i,k) .le. 0. .and.    &
                   U_ca(i,k) .ge. U00p(i,k)) .and. &
                    (qa_upd(i,k) +  &
                      convective_humidity_area(i,k) .le. 1.))   then
                  da_ls(i,k) =    &
                          -0.5*(1. - qa_upd(i,k) -  &
                                   convective_humidity_area(i,k))*  &
                                (1. - qa_upd(i,k) -    &
                                  convective_humidity_area(i,k) )* &
                                       dqs_ls(i,k)/qs(i,k) /  &
                                      max(1.-U_ca(i,k), qmin)
                else
                  da_ls(i,k) = 0.
                endif
              end do
            end do

      !------------------------------------------------------------------------
      !       Turbulent erosion of clouds
      !
      !       As in Tiedtke (1993) this is calculated using the eros_scale
      !       parameter as:
      !
      !   (9) dql/dt    =  - qa * eros_scale * (qs - qv) * (ql/ ql+qi )
      !
      !  (10) dqi/dt    =  - qa * eros_scale * (qs - qv) * (qi/ ql+qi )
      !
      !  (11) dqa/dt    =  - qa * eros_scale * (qs - qv) * (qa/ ql+qi )
      !
      !       for which the erosion sink term (B in equation 13) is
      !
      !  (12) B = qa * eros_scale * (qs - qv) / (ql + qi)  
      !
      !------------------------------------------------------------------------
            do k=1,kdim
              do i=1,idim
                if (ql_upd(i,k) .gt. qmin .or.   &
                            qi_upd (i,k).gt. qmin) then
                  D_eros(i,k) =   &
                        qa_upd(i,k) * erosion_scale(i,k) *   &
                                   dtcloud * qs(i,k) *      &
                          (1.-U_ca(i,k)) /   &
                       (qi_upd(i,k) + ql_upd(i,k))

                  if (pfull(i,k) .gt. 400.e02) then
                    D_eros(i,k) =  &
                         D_eros(i,k) + efact*  &
                                          D_eros(i,k)* &
                              ((pfull(i,kdim) -   &
                                                   pfull(i,k))/  &
                                    (pfull(i,kdim) - 400.e02))
                  else
                    D_eros(i,k)=  &
                          D_eros(i,k) +  &
                                         efact*D_eros(i,k)

                  endif
                else
                  D_eros(i,k) = 0.
                endif
              END DO        
            END DO    

      !------------------------------------------------------------------------
      !       The next step is to analytically integrate the saturated volume
      !       fraction equation.  This follows the Tiedtke approach
      !
      !       The qa equation is written in the form:
      !
      !  (13) dqa/dt    =   (1.-qa) * A   -  qa * B 
      !
      !       Note that over the physics time step, A, B are assumed to be 
      !       constants.
      !
      !       Defining qa(t) = qa0 and qa(t+dtcloud) = qa1, the analytic
      !       solution of the above equation is:
      !
      !  (14) qa1 = qaeq -  (qaeq - qa0) * exp (-(A+B)*dtcloud)
      ! 
      !       where qaeq is the equilibrium cloud fraction that is approached
      !       with an time scale of 1/(A+B),
      !
      !  (15) qaeq  =  A/(A+B)
      !
      !
      !       To diagnose the magnitude of each of the right hand terms of
      !       (13) integrated over the time step, define the average cloud
      !       fraction in the interval t to t + dtcloud qabar as:
      !
      !  (16) qabar  = qaeq - [ (qa1-qa0) / ( dtcloud * (A+B) ) ]
      ! 
      !       from which the magnitudes of the A and B terms integrated
      !       over the time step are:
      !
      !       A * (1-qabar)    and    -B * (qabar)
      !
      !       Additional notes on this analytic integration:
      !
      !       1.   For large-scale cloud formation or destruction from 
      !            the dqs_ls term the contributions to A or B are defined
      !            from:
      !
      !  (19)      A_ls * (1. - qa) = da_ls / dtcloud      if da_ls >= 0.
      ! 
      !  (20)      B_ls * qa        = da_ls / dtcloud      if da_ls < 0.
      !
      !
      !       3.   Qa goes to zero only in the case of ql and qi less than or
      !            equal to     qmin; see 'cloud destruction code' near the end of 
      !            this loop over levels.
      !------------------------------------------------------------------------
              
      !------------------------------------------------------------------------
      !    compute A_dt; This is assigned to the large-scale source term
      !    following (18). Reset B_dt.
      !------------------------------------------------------------------------
            do k=1,kdim
              do i=1,idim
                A_dt(i,k) = da_ls(i,k)/   &
                                  max((1.-qa_upd(i,k)), qmin)
                B_dt(i,k) = D_eros(i,k)
        
      !------------------------------------------------------------------------
      !    do analytic integration.      
      !------------------------------------------------------------------------
                if ( (A_dt(i,k) .gt. Dmin) .or.   &
                     (B_dt(i,k) .gt. Dmin) )  then
                  qa0(i,k) = qa_upd(i,k)
                  qaeq(i,k) = A_dt(i,k)/(A_dt(i,k)  + B_dt(i,k))
                  qa1(i,k) = qaeq(i,k) - (qaeq(i,k) - qa0(i,k)) * &
                                          exp ( -1.*(A_dt(i,k)+B_dt(i,k)) )
                  qabar(i,k) = qaeq(i,k) - ((qa1(i,k) - qa0(i,k))/  &
                                                   (A_dt(i,k) + B_dt(i,k)))
                else
                  qa0(i,k)   = qa_upd(i,k)
                  qaeq(i,k)  = qa0(i,k)
                  qa1(i,k)   = qa0(i,k)   
                  qabar(i,k) = qa0(i,k)  
                endif 
              END DO        
            END DO    

      !       do k=1,kdim
      !         do i=1,idim
      ! !------------------------------------------------------------------------
      ! !    save some diagnostics.
      ! !------------------------------------------------------------------------
      !           if ( (A_dt(i,k) .gt.     Dmin) .and.   &
      !                (B_dt(i,k) .gt.     Dmin) )  then
      !             if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
      !               diag_4d(i,k,diag_pt%qadt_lsform) =  A_dt(i,k)*  &
      !                                           (1.-qabar(i,k))*inv_dtcloud
      !             end if
      !             if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
      !               diag_4d(i,k,diag_pt%qadt_eros)  =    &
      !                    ((qa1(i,k) - qa0(i,k))*inv_dtcloud )- &
      !                                    diag_4d(i,k,diag_pt%qadt_lsform)   
      !             end if
      !           else if (A_dt(i,k) .gt. Dmin) then
      !             if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
      !               diag_4d(i,k,diag_pt%qadt_lsform) =   &
      !                         (qa1(i,k) - qa0(i,k))*inv_dtcloud
      !             end if 
      !           else if (B_dt(i,k) .gt. Dmin)  then
      !             if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
      !               diag_4d(i,k,diag_pt%qadt_eros)  =   &
      !                         (qa1(i,k) - qa0(i,k))*inv_dtcloud
      !             end if 
      !           endif
      !         END DO    
      !       END DO    

      !-----------------------------------------------------------------------
      !   define value needed for diagnostic calculation.
      ! !-----------------------------------------------------------------------
      !       if (diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0) then
      !         diag_4d(:,:,diag_pt%qadt_ahuco) = qa1(:,:)
      !       end if

      !------------------------------------------------------------------------
      !    limit cloud area to be no more than that which is not being
      !    taken by convective clouds.
      !------------------------------------------------------------------------
            if (limit_conv_cloud_frac) then
              qa1 = MIN(qa1, 1.0 -convective_humidity_area)
            endif
                       
      !------------------------------------------------------------------------
      !    complete calculation of diagnostic.
      !------------------------------------------------------------------------
            ! if (diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0) then
            !   diag_4d(:,:,diag_pt%qadt_ahuco) =    &
            !             (qa1(:,:) - diag_4d(:,:,diag_pt%qadt_ahuco))* & 
            !                                                         inv_dtcloud 
            ! end if

      !------------------------------------------------------------------------
      !    set total tendency term and update cloud fraction    
      !------------------------------------------------------------------------
            SA = (SA + qa1) - qa0
            qa_upd = qa1
              
      !------------------------------------------------------------------------
      !       The next step is to calculate the change in condensate
      !       due to non-convective condensation, dcond_ls. Note that this is
      !       not the final change but is used only to apportion condensate
      !       change between phases. According to Tiedtke 1993 this takes the
      !       form:
      !
      !  (21) dcond_ls = -1. * (qa +  0.5*da_ls) * dqs_ls
      !
      !       Here the 0.5*da_ls represents using a midpoint cloud fraction.
      !       This is accomplished by using the variable qabar.
      !
      !       Also, save term needed when predicted droplets is active (tmp5).
      !------------------------------------------------------------------------
            IF (use_qabar) THEN
              dcond_ls = -1.*qabar*dqs_ls
              delta_cf = A_dt*(1.-qabar)
            ELSE
              dcond_ls = -1.*qa_upd*dqs_ls   
              delta_cf = A_dt*(1.-qa_upd)
            END IF      

      !-----------------------------------------------------------------------


      end subroutine tiedtke_macro_nopdf_nosuper
      
      SUBROUTINE tiedtke_macro_nopdf_super (   &
                         idim, kdim, limit_conv_cloud_frac, mc_full, convective_humidity_area, omega, radturbten, pfull, tin, qin, airdens, dqsdT, gamma, qs, qsl, qsi, rh_crit, rh_crit_min, &
                         qa_upd, ql_upd, qi_upd, ST, SQ, da_ls, D_eros, dcond_ls, delta_cf, hom,  &
                         edum, SA, U00p, erosion_scale, dtcloud, dmin, qmin, efact, cp_air, grav, hlv, hls, rvgas, tfreeze) 
        use machine, only : kind_phys
        use ccpp_saturation, only: get_qs
        implicit none
      !------------------------------------------------------------------------
      !   subroutine tiedtke_macro_nopdf_super calculates non-convective
      !   condensation while allowing ice supersaturation in the clear-sky 
      !   portion of grid boxes, as is needed in order to activate ice crystals 
      !   (see Tompkins et al 2007, Salzmann et al , 2010).
      !------------------------------------------------------------------------

      INTEGER,                            INTENT(IN )    :: idim, kdim
      logical,              intent(in) :: limit_conv_cloud_frac
      real(kind=kind_phys), intent(in) :: dtcloud, dmin, qmin, efact, cp_air, grav, hlv, hls, rvgas, tfreeze
      real(kind=kind_phys), intent(in) :: mc_full(:,:) !total net convective mass flux on full levels kg m-2 s-1
      real(kind=kind_phys), intent(in) :: convective_humidity_area(:,:) !grid box area affected by the convective clouds
      real(kind=kind_phys), intent(in) :: omega(:,:) !lagrangian_tendency_of_air_pressure
      real(kind=kind_phys), intent(in) :: radturbten(:,:) !sum of temperature tendency due to radiation and turbulence (K s-1) (interstitial -- reset every timestep)
      real(kind=kind_phys), intent(in) :: pfull(:,:)
      real(kind=kind_phys), intent(in) :: tin(:,:) !i,k dims; this refers to the scheme's input state temperature (can either be process-split or time-split)
      real(kind=kind_phys), intent(in) :: qin(:,:) !i,k dims; this refers to the scheme's input state water vapor (can either be process-split or time-split)
      real(kind=kind_phys), intent(in) :: airdens(:,:) !air density (kg m-3) calculated in lscloud_driver
      real(kind=kind_phys), intent(in) :: dqsdT(:,:)
      real(kind=kind_phys), intent(in) :: gamma(:,:), qs(:,:), qsl(:,:), qsi(:,:)
      real(kind=kind_phys), intent(in) :: rh_crit(:,:), rh_crit_min(:,:)
      real(kind=kind_phys), intent(inout) :: qa_upd(:,:)
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_upd, qi_upd
      real(kind=kind_phys), intent(inout) :: da_ls(:,:)
      real(kind=kind_phys), intent(inout) :: D_eros(:,:)
      real(kind=kind_phys), intent(inout) :: dcond_ls(:,:)
      real(kind=kind_phys), intent(inout) :: delta_cf(:,:)
      real(kind=kind_phys), intent(in) :: hom(:,:)
      REAL, dimension(idim,kdim),    INTENT(IN )    :: ST, SQ, edum, U00p,&
                                                            erosion_scale
                                                            
      REAL, dimension(idim,kdim),    INTENT(INout ) :: SA

      !-------------------------------------------------------------------------
      !-----local variables-----
      ! 
      !       A_dt           product of A and dtcloud in     dimensionless in 
      !                      in the analytic integration     qa integration
      !                      of the qa equation, or C and
      !                      dtcloud in the analytic         kg condensate/
      !                      integration of the ql and qi    kg air in ql or 
      !                      equations.                      qi integration
      !
      !       B_dt           product of B and dtcloud in     dimensionless in
      !                      in the analytic integration     qa, ql, and qi
      !                      of the qa equation, or D and    integration
      !                      dtcloud in the analytic         
      !                      integration of the ql and qi    
      !                      equations.                      
      !
      !       qa0            value of cloud fraction or      dimensionless or
      !                      cloud condensate at the         kg condensate /
      !                      initial time                    kg air
      !        
      !       qa1            value of cloud fraction or      dimensionless or
      !                      cloud condensate at the final   kg condensate /
      !                      time                            kg air
      !
      !       qaeq           equilibrium value of cloud      dimensionless or
      !                      fraction or cloud condensate    kg condensate /
      !                      that the analytic integration   kg air
      !                      approaches                   
      !       qabar          mean value of cloud fraction    dimensionless or
      !                      or cloud condensate over the    kg condensate /
      !                      t0 to t0 + dtcloud interval     kg air
      !       qve            Tiedtke environmenmtal 
      !                      humidity
      !
      !------------------------------------------------------------------------
            real, dimension(idim, kdim) :: deltpg, dqs_ls, drhcqs_ls, dum,&
                                                qve 
            real, dimension(idim, kdim) :: A_dt, B_dt, qa1, qa0,  &
                                                qabar, qaeq, qa_t, tmp0, tmp1, &
                                                drhcqsdT, beta, ttmp, qtmp,  &
                                                qs_t, qs_l, qs_i, dqsdT1, gamma1, &
                                                es_t, es_l, es_i
            real, dimension(idim)       :: qagtmp,qcgtmp,qvgtmp
            REAL                             :: eslt, qvs, qs_d, qvi, esit, ul
            REAL                             :: qvmax, esat0, tmp1s            
            INTEGER                          :: ns, id
            INTEGER                          :: i,k
            
            integer :: alg_choice, alg_flatau_92
            alg_choice = 1
            alg_flatau_92 = 1

      !------------------------------------------------------------------------
      !    calculate the derivative of the threshold relative humidity over water!    (for temps > 233) and over ice (temps < 233) with respect to 
      !    temperature, and the L/Cp drhqsdT term (beta) in the
      !    temperature equation to be used in this calculation. for
      !    the water case (T >= 233), no supersaturation is allowed. For the ice
      !    case, superstaurations up to rh_crit are permitted.
      !------------------------------------------------------------------------
            do k=1,kdim
              do i=1,idim
                if (tin(i,k) .GE. TFREEZE - 40. ) then 
                  drhcqsdT(i,k) =  HLV*qsl(i,k)/  &
                                            (RVGAS*tin(i,k)**2)
                  beta(i,k) = drhcqsdT(i,k)*HLV/CP_AIR
                else
                  drhcqsdT(i,k) = rh_crit(i,k)*   &
                              HLS*qsi(i,k)/  &
                                           (RVGAS*tin(i,k)**2) +  &
                               hom(i,k)*qsi(i,k)*  &
                           (2.*0.0073*(tin(i,k) - TFREEZE) + 1.466)
                  beta(i,k) = drhcqsdT(i,k)*HLS/CP_AIR
                endif 
              end do
            end do

            do k=1,kdim
              do i=1,idim

      !------------------------------------------------------------------------
      !    calculate environmental specific humidity (qve) based on input cloud 
      !    fraction. take into account convective cloud presence.
      !------------------------------------------------------------------------
                qa_t(i,k) =    &
                  MAX(MIN(qa_upd(i,k) +   &
                              convective_humidity_area(i,k), 1.),0.)
                qve(i,k) =  (qin(i,k) - qa_t(i,k)*  &
                                qs(i,k) ) /   &
                                           MAX((1. - qa_t(i,k) ), qmin)

      !------------------------------------------------------------------------
      !       The first step is to compute the change in qs due to large-
      !       scale processes, dqs_ls.   In Tiedtke, it has contributions from 
      !       large-scale uplift, convection induced compensating subsidence,
      !       turbulence cooling and radiative cooling.  dqs_ls has the form:
      !
      !               (((omega+ grav*Mc)/airdens/cp)+radturbten)*dqsdT*dtcloud
      !   (6) dqs_ls= --------------------------------------------------------
      !                  1.  +   ( qa +  (da_ls/2.) ) * gamma
      !
      !       Here da_ls is the increase in cloud fraction due to non-
      !       convective processes.  Because this increase is also a function
      !       of dqs_ls, a quadratic equation must be solved for dqs_ls in
      !       the case that da_ls is not equal to zero.
      !------------------------------------------------------------------------
                dum(i,k) = (((omega(i,k) + grav*  &
                                 mc_full(i,k))/ &
                                     airdens(i,k)/cp_air) + &
                                       radturbten(i,k))*dtcloud
                dqs_ls(i,k) = dum(i,k) * dqsdT(i,k)
                drhcqs_ls(i,k) = dum(i,k) *drhcqsdT(i,k)
              end do
            end do

       !------------------------------------------------------------------------
       !       The next step is to compute the change in saturated volume
       !       fraction due to non-convective condensation, da_ls.   This 
       !       occurs in two conditions:
       !
       !       (a)  dqs_ls < 0. and U00 < U < 1., where U00 is the threshold
       !            relative humidity for non-convective condensation. Note 
       !            that if U is greater than or equal to 1., ideally qa = 1,
       !            and da_ls = 0.  However this may not be the case for 
       !            numerical reasons so this must be assured after analytic 
       !            integration of the qa equation.
       !
       !            For these cases the change in saturated volume fraction is:
       !
       !   (7)      da_ls = - (1.-qa)*(1.-qa)*dqs_ls/2./qs/(1.-U)
       !
       !            This formula arises from the assumption that vapor is uni-
       !            formly distributed in the range [qv_clr - (qs - qv_clr),qs]
       !            where qv_clr is the amount of vapor in the unsaturated 
       !            volume and is given from the following equation:
       !
       !   (8)      qv  =   qa * qs      +   (1.-qa) * qv_clr
       !          
       !            Implicit in equation (7) is the following assumption:
       !            As qsat changes, the distribution of qv+ql+qi 
       !            remains constant.  That is as qsat rises, portions where
       !            qv+ql+qi > qsat+dqsat remain saturated.  This can only
       !            occur if it is assumed that ql+qi evaporate-sublimate or
       !            condense-deposit to keep qv = qsat. 
       !
       !       (b)  dqs_ls > 0.  Ideally some portion of the cloud should
       !            evaporate however this is not accounted for at present.
       !----------------------------------------------------------------------- 
            where (dqs_ls .le. 0. .and.    &
                        qve .ge. U00p*rh_crit_min*qs &
                                                             .and. qa_t .lt. 1.)
              tmp0 =  (1 + gamma*qa_t)*drhcqs_ls - beta*qa_t* &
                                                                          dqs_ls
              tmp1 = SQRT( (1.+qa_t*gamma)**2. - (1.-qa_t)* &
                               beta * tmp0 /MAX( (  rh_crit *  &
                                qs - qve ), qmin ) ) -  &
                                                  (1.+qa_t*gamma) 
              tmp1 = -1.*tmp1/((1. - qa_t)*beta/(2.*MAX(   &
                         (rh_crit*qs - qve), qmin )) )
              drhcqs_ls = min(tmp1, tmp0 / MAX(1.+ gamma*qa_t +&
                                               beta*0.5*(1. - qa_t), qmin))
              da_ls = -0.5*(1. - qa_t)*drhcqs_ls/   &
                   MAX( rh_crit*qs - qve, qmin)
              dqs_ls = (dqs_ls - gamma*0.5*  &
                        da_ls*drhcqs_ls)/   &
                                                (1. + qa_t*gamma)
            elsewhere
              drhcqs_ls = 0.
              da_ls = 0.
              dqs_ls = dqs_ls/(1. + qa_t*gamma)
            endwhere

      !------------------------------------------------------------------------
      !       Turbulent erosion of clouds
      !
      !       As in Tiedtke (1993) this is calculated using the eros_scale
      !       parameter as:
      !
      !   (9) dql/dt    =  - qa * eros_scale * (qs - qv) * (ql/ ql+qi )
      !
      !  (10) dqi/dt    =  - qa * eros_scale * (qs - qv) * (qi/ ql+qi )
      !
      !  (11) dqa/dt    =  - qa * eros_scale * (qs - qv) * (qa/ ql+qi )
      !
      !       for which the erosion sink term (B in equation 13) is
      !
      !  (12) B = qa * eros_scale * (qs - qv) / (ql + qi)  
      !
      !------------------------------------------------------------------------
            DO k=1,kdim
              DO i=1,idim
                if (ql_upd(i,k) .gt. qmin .or. &
                               qi_upd (i,k).gt. qmin) then
                  D_eros(i,k) =    &
                            qa_upd(i,k)*erosion_scale(i,k)*   &
                                                               dtcloud*    &
                       (MAX(qs(i,k) - qin(i,k),&
                             0.)/(1. + gamma(i,k)) ) /      &
                        (qi_upd(i,k) + ql_upd(i,k))

                  if (pfull(i,k) .gt. 400.e02) then
                    D_eros(i,k) =   &
                         D_eros(i,k) + efact*  &
                                    D_eros(i,k)*  &
                        ((pfull(i,kdim) - pfull(i,k))/  &
                                    (pfull(i,kdim) - 400.e02))

                  else
                    D_eros(i,k) =  &
                         D_eros(i,k) + efact*  &
                                              D_eros(i,k)
                  endif
                else
                  D_eros(i,k) = 0.
                endif
              END DO        
            END DO    

      !------------------------------------------------------------------------ 
      !       The next step is to analytically integrate the saturated volume
      !       fraction equation.  This follows the Tiedtke approach
      !
      !       The qa equation is written in the form:
      !
      !  (13) dqa/dt    =   (1.-qa) * A   -  qa * B 
      !
      !       Note that over the physics time step, A, B are assumed to be 
      !       constants.
      !
      !       Defining qa(t) = qa0 and qa(t+dtcloud) = qa1, the analytic
      !       solution of the above equation is:
      !
      !  (14) qa1 = qaeq -  (qaeq - qa0) * exp (-(A+B)*dtcloud)
      ! 
      !       where qaeq is the equilibrium cloud fraction that is approached
      !       with an time scale of 1/(A+B),
      !
      !  (15) qaeq  =  A/(A+B)
      !
      !
      !       To diagnose the magnitude of each of the right hand terms of
      !       (13) integrated over the time step, define the average cloud
      !       fraction in the interval t to t + dtcloud qabar as:
      !
      !  (16) qabar  = qaeq - [ (qa1-qa0) / ( dtcloud * (A+B) ) ]
      ! 
      !       from which the magnitudes of the A and B terms integrated
      !       over the time step are:
      !
      !       A * (1-qabar)    and    -B * (qabar)
      !
      !       Additional notes on this analytic integration:
      !
      !       1.   For large-scale cloud formation or destruction from 
      !            the dqs_ls term the contributions to A or B are defined
      !            from:
      !
      !  (19)      A_ls * (1. - qa) = da_ls / dtcloud      if da_ls >= 0.
      ! 
      !  (20)      B_ls * qa        = da_ls / dtcloud      if da_ls < 0.
      !
      !
      !       3.   Qa goes to zero only in the case of ql and qi less than or
      !            equal to     qmin; see 'cloud destruction code' near the end of 
      !            this loop over levels.
      !------------------------------------------------------------------------
              
      !-------------------------------------------------------------------------
      !    compute A_dt; This is assigned to the large-scale source term
      !    following (18). Reset B_dt.
      !-------------------------------------------------------------------------
         
            do k=1,kdim
              do i=1,idim
                A_dt(i,k) = da_ls(i,k)/   &
                                                max((1.-qa_t(i,k)), qmin)
                B_dt(i,k) = D_eros(i,k)
                if ( (A_dt(i,k) .gt. Dmin) .or.   &
                                          (B_dt(i,k) .gt. Dmin) ) then 
                  qa0(i,k)   = qa_upd(i,k)
                  qaeq(i,k)  =                                     &
                                     A_dt(i,k)/(A_dt(i,k) + B_dt(i,k))
                  qa1(i,k)  = qaeq(i,k) - (qaeq(i,k) - qa0(i,k))* &
                                         exp(-1.*(A_dt(i,k) + B_dt(i,k)) )
                  qabar(i,k) = qaeq(i,k) - ((qa1(i,k) - qa0(i,k))/  &
                                             (A_dt(i,k) + B_dt(i,k)))
                else
                  qa0(i,k)   = qa_upd(i,k)
                  qaeq(i,k)  = qa0(i,k)   
                  qa1(i,k)   = qa0(i,k)   
                  qabar(i,k) = qa0(i,k)  
                endif
              end do
            end do

      !-------------------------------------------------------------------------
      !    output some diagnostics.
      !-------------------------------------------------------------------------
            ! do k=1,kdim
            !   do i=1,idim
            !     if ( (A_dt(i,k) .gt. Dmin) .and.   &
            !          (B_dt(i,k) .gt. Dmin) )  then
            !       if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
            !         diag_4d(i,k,diag_pt%qadt_lsform) =  A_dt(i,k)*  &
            !                       (1. - qabar(i,k))*inv_dtcloud
            !       end if
            !       if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
            !         diag_4d(i,k,diag_pt%qadt_eros)  =   &
            !              ((qa1(i,k) - qa0(i,k))*inv_dtcloud ) - &
            !                                diag_4d(i,k,diag_pt%qadt_lsform)  
            !       end if
            !     else if (A_dt(i,k) .gt. Dmin) then
            !       if (diag_id%qadt_lsform + diag_id%qa_lsform_col > 0) then
            !         diag_4d(i,k,diag_pt%qadt_lsform) =    &
            !                   (qa1(i,k) - qa0(i,k))*inv_dtcloud
            !       end if 
                ! else if (B_dt(i,k) .gt. Dmin)  then
            !       if ( diag_id%qadt_eros + diag_id%qa_eros_col > 0 ) then
            !         diag_4d(i,k,diag_pt%qadt_eros)  =    &
            !                    (qa1(i,k) - qa0(i,k))*inv_dtcloud
            !       end if 
            !     endif
            !   end do
            ! end do

      !-----------------------------------------------------------------------
      !   define value needed for diagnostic calculation.
      !-----------------------------------------------------------------------
            !if ( diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0 ) then
            !  diag_4d(:,:,:,diag_pt%qadt_ahuco)  = qa1
            !end if

      !-------------------------------------------------------------------------
      !    limit cloud area to be no more than that which is not being
      !    taken by convective clouds.
      !-------------------------------------------------------------------------
            if (limit_conv_cloud_frac) then
              qa1 = MIN(qa1, 1.0 - convective_humidity_area)
            endif

      !------------------------------------------------------------------------
      !    complete calculation of diagnostic.
      !------------------------------------------------------------------------
            ! if ( diag_id%qadt_ahuco + diag_id%qa_ahuco_col > 0 ) then
            !   diag_4d(:,:,diag_pt%qadt_ahuco)  =   &
            !        (qa1(:,:) - diag_4d(:,:,diag_pt%qadt_ahuco))*inv_dtcloud 
            ! end if
                       
      !-------------------------------------------------------------------------
      !    set total tendency term and update cloud fraction.    
      !-------------------------------------------------------------------------
            SA = (SA + qa1) - qa0
            qa_upd = qa1
            delta_cf = MAX(qa1 - qa0 , 0.)

      !-------------------------------------------------------------------------
      !       The next step is to calculate the change in condensate
      !       due to non-convective condensation, dcond_ls. Note that this is
      !       not the final change but is used only to apportion condensate
      !       change between phases. According to Tiedtke 1993 this takes the
      !       form:
      !
      !  (21) dcond_ls = -1. * (qa +  0.5*da_ls) * dqs_ls
      !
      !       Here the 0.5*da_ls represents using a midpoint cloud fraction.
      !       This is accomplished by using the variable qabar.

      !cms but qabar is not limited to area outside  conv. clouds ...
      !!      dcond_ls = -1. * MIN(qabar, 1.- ahuco)  * dqs_ls
      !-----------------------------------------------------------------------
            dcond_ls = -1.*qa_t*dqs_ls -    &
                 0.5*MAX(MIN(da_ls, 1. - qa_t),0.)*drhcqs_ls 

      !-----------------------------------------------------------------------
      !     compute qs and condensation using updated temps and vapor. 
      !     see Tompkins et al., 2007
      !-----------------------------------------------------------------------
            ttmp=tin + ST        
            qtmp= qin + SQ         
            !CALL compute_qs_x1 (idim, kdim, ttmp, pfull, &
            !                                 qs_t, qs_l, qs_i, dqsdT1, gamma1 )
            
            call get_qs(ttmp, pfull, alg_choice, alg_flatau_92, qs_t, qs_l, qs_i, es_t, es_l, es_i, dqsdT1, gamma1)
            
            DO k=1,kdim
              DO i=1, idim
                qvmax = qs_t(i,k)*(qa_t(i,k) + (1. -  qa_t(i,k))*  &
                                                 rh_crit(i,k) )
                tmp1s =  max(0., (qtmp(i,k) - qvmax) )/(1. + gamma1(i,k))
                ! limit 
                IF (dcond_ls(i,k) .GT. 0. .OR.    &
                                                        tmp1s .GT. 0. )  THEN
                  dcond_ls(i,k) =   &
                              MAX(dcond_ls(i,k), tmp1s) 
                ELSE
                !don't evaporate into saturated grid cells
                !CHECK
                  if ( qve(i,k) .lt. qs_t(i,k)) then
                    tmp1s =  max(0.,( qs_t(i,k) - qtmp(i,k) ))/   &
                                   (1. + gamma1(i,k)) ! evap .le. qs - qtmp
                    dcond_ls(i,k) = MAX( -1.*tmp1s,  &
                                             dcond_ls(i,k))
                  end if
                END IF
              END DO
            END DO

      !-------------------------------------------------------------------  


      end SUBROUTINE tiedtke_macro_nopdf_super
    end module tiedtke_prog_clouds
