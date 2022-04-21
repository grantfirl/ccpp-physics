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
      subroutine tiedtke_prog_clouds_run (idim,kdim,do_pdf_clouds,single_gaussian_pdf,do_aero_eros,ae_ub,ae_lb,ae_N_ub,ae_N_lb,dt,SA,SQ,gamma,qs,qin,ql_in,qi_in,qa_in,ql_upd,qi_upd,qa_upd,qvg,da_ls,delta_cf,dcond_ls,errmsg,errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim
      logical, intent(in) :: do_pdf_clouds, single_gaussian_pdf, do_aero_eros
      real(kind=kind_phys), intent(in) :: dt, ae_ub, ae_lb, ae_N_ub, ae_N_lb
      real(kind=kind_phys), intent(in) :: gamma(:,:), qs(:,:) !i,k dims
      real(kind=kind_phys), intent(in) :: drop1(:,:) !cloud droplet number concentration (cm-3) (determined in aerosol_cloud.F90/determine_activated_aerosol)
      real(kind=kind_phys), intent(in) :: qin(:,:) !i,k dims; this refers to the scheme's input state water vapor
      real(kind=kind_phys), intent(in) :: ql_in(:,:), qi_in(:,:), qa_in(:,:) !these vars are initialized every physics timestep to the advected tracers values of ql, qi in lscloud_driver.F90/lscloud_alloc
      real(kind=kind_phys), intent(in) :: ql_upd(:,:), qi_upd(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_alloc
      real(kind=kind_phys), intent(inout) :: qa_upd(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_alloc
      real(kind=kind_phys), intent(inout) :: qvg(:,:), da_ls(:,:), delta_cf(:,:), dcond_ls(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_alloc
      real(kind=kind_phys), intent(inout) :: dqa_dt_prod(:,:), dqa_dt_loss(:,:) !these vars are initialized every physics timestep to 0 in lscloud_driver.F90/lscloud_driver
      
      real(kind=kind_phys), intent(in)    :: SQ(:,:) !this is set to zero in lscloud_driver.F90/lscloud_driver
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
!       i,j,k          do-loop indices
!-------------------------------------------------------------------------
      real(kind=kind_phys) :: dt_inv
      real(kind=kind_phys), dimension(idim,kdim) :: mdum,bdum,edum

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      !GJF need to call compute_qs_a to calculate qs,gamma
      
      dt_inv = 1.0/dt
      
!------------------------------------------------------------------------
!    process the non-convective condensation for pdf clouds. 
!                                                                      
!                NON-CONVECTIVE CONDENSATION                          
!                                                                      
!                STATISTICAL CLOUD FRACTION                            
!                                                                      
!------------------------------------------------------------------------
      if (do_pdf_clouds) then
        if (single_gaussion_pdf) then      
          call tiedtke_macro_Single_Gaussian_pdf(idim, kdim,    & 
            gamma, qs, qin, &
            ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, SQ, qvg, da_ls, delta_cf, dcond_ls, &
            SA)
        else
          call tiedtke_macro_pdf (     &
            idim, kdim, gamma, qs, qin, ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, &
            SQ, qvg, D_eros, da_ls, delta_cf, dcond_ls, dqa_dt_prod, dqa_dt_loss, SA)
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
          edum = mdum*Particles%drop1 + bdum 
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
            where (Input_mp%pfull(:,:,k).gt.   &
                          0.8*Input_mp%phalf(:,:,KDIM+1)) 
              U00p(:,:,k) = U00 + (1.- U00)* &
                         (((Input_mp%pfull(:,:,k) -   &
                            (0.8*Input_mp%phalf(:,:,KDIM+1)))/  &
                               (0.2*Input_mp%phalf(:,:,KDIM+1)) )**2.)
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
          u00pr = u00p + (1. - u00p)*C2ls_mp%convective_humidity_area
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
            do j=1,jdim
              do i=1,idim

!-----------------------------------------------------------------------
!    Enhanced erosion in convective layers
!-----------------------------------------------------------------------
                if (include_neg_mc) then
                  if (abs(C2ls_mp%mc_full(i,j,k)) .gt. mc_thresh) then 
                    erosion_scale(i,j,k) = eros_scale_c
                  endif
                else
                  if (C2ls_mp%mc_full(i,j,k) .gt. mc_thresh) then
                    erosion_scale(i,j,k) = eros_scale_c
                  endif
                endif
                if ((Input_mp%diff_t(i,j,K) .gt. diff_thresh) .or.&
                    (Input_mp%diff_t(i,j,min(k+1,KDIM)) .gt.  &
                                                       diff_thresh) ) then
                  erosion_scale(i,j,K) = eros_scale_t
                endif 
              END DO
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
               idim, jdim, kdim, C2ls_mp, Input_mp, Atmos_state,   &
               Cloud_state, ST, SQ, Cloud_processes, Particles,   &
               Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,   &
               Lsdiag_mp_control%diag_id, Lsdiag_mp_control%diag_pt,  &
                                       edum, SA, U00p, erosion_scale)
        else

!------------------------------------------------------------------------
!    if supersaturation is to be allowed in the clear-sky portion of the
!    gridbox, call subroutine tiedtke_macro_nopdf_super to calculate the 
!    condensation.
!------------------------------------------------------------------------
          call tiedtke_macro_nopdf_super  (     &
               idim, jdim, kdim, C2ls_mp, Input_mp, Atmos_state,   &
               Cloud_state, ST, SQ, Cloud_processes, Particles,   &
               Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,   &
               Lsdiag_mp_control%diag_id, Lsdiag_mp_control%diag_pt,  &
                                            edum, SA, U00p, erosion_scale)
        endif
      endif

      end subroutine tiedtke_prog_clouds_run

      SUBROUTINE tiedtke_macro_Single_Gaussian_pdf (idim, kdim, gamma, qs,   &
                              qin, &
                              ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, SQ, qvg, da_ls, delta_cf, dcond_ls, &
                              diag_4d, diag_id, diag_pt,       SA)  

      !----------------------------------------------------------------------!
      !                                                                      !
      !                     STATISTICAL CLOUD SCHEME                         !
      !                                                                      !
      !-----------------------------------------------------------------------

      integer,                             intent(in)    :: idim, kdim     
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
              qta(:,:,:) = max(qmin, qin +  &
                                           ql_in + qi_in)
            else
      !       qta(:,:,:) = max(qmin, qv_in + SQ +  &
              qta(:,:,:) = max(qmin, qin + SQ +  &
                                           ql_upd + qi_upd)
            endif
            qtqsa(:,:,:) = qta(:,:,:) - qs      

            s =  qtqsa(:,:,:) / (1.0 +  gamma )

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
                   idim, kdim, gamma, qs, qin, ql_in, qi_in, qa_in, ql_upd, qi_upd, qa_upd, &
                   SQ, qvg, D_eros, da_ls, delta_cf, dcond_ls, dqa_dt_prod, dqa_dt_loss, SA)

      !----------------------------------------------------------------------!
      !                                                                      !
      !                     STATISTICAL CLOUD SCHEME                         !
      !                                                                      !
      !-----------------------------------------------------------------------

      integer,                           intent(in )   :: idim, kdim     
      real(kind=kind_phys), intent(in),    dimension(:,:) :: gamma, qs
      real(kind=kind_phys), intent(in),    dimension(:,:) :: qin
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_in, qi_in, qa_in
      real(kind=kind_phys), intent(in),    dimension(:,:) :: ql_upd, qi_upd
      real(kind=kind_phys), intent(inout), dimension(:,:) :: qa_upd
      real(kind=kind_phys), intent(inout), dimension(:,:) :: qvg, D_eros, da_ls, delta_cf, dcond_ls
      real(kind=kind_phys), intent(inout), dimension(:,:) :: dqa_dt_prod, dqa_dt_loss
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
              
                   if (qs_norm(i).le.1.) then
                       
                       icbp = incomplete_beta(max(0.,qs_norm(i)), &
                                            p = betaP    , q =     betaP)
                       icbp1= incomplete_beta(max(0.,qs_norm(i)), &
                                            p = betaP + 1, q =     betaP)
                       qagtmp(i) = 1.-icbp
                       qcgtmp(i) = (qtbar(i)-qtmin(i))*(1.-icbp1)&
                                     - qs_norm(i)*deltaQ(i)*(1.-icbp)    
                       qcgtmp(i) = qcgtmp(i)/   &
                                                (1.+gamma(i,k))
                       qvgtmp(i) = qtmin(i,j) + &
                                     0.5*(icbp1/max(icbp, qmin))*deltaQ(i)
                   
                       !bound very very small cloud fractions which may
                       !cause negative cloud condensates due to roundoff 
                       !errors or similar errors in the beta table lookup.
                       if((qagtmp(i).lt.0.).or.(qcgtmp(i).le.0.))then
                            qagtmp(i) = 0.
                            qcgtmp(i) = 0.
                            qvgtmp(i) = qtbar(i)
                       end if
                       
                   else             
                       qagtmp(i) = 0.
                       qcgtmp(i) = 0.
                       qvgtmp(i) = qtbar(i)             
                   end if
                   
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

            diag_dt_prod = max(qa1 - qa0, 0.)*inv_dtcloud
            diag_dt_loss = max(qa0 - qa1, 0.)*inv_dtcloud
            
      !------------------------------------------------------------------------
       


      end SUBROUTINE tiedtke_macro_pdf
    end module tiedtke_prog_clouds
