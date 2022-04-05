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
      subroutine tiedtke_prog_clouds_run (idim,kdim,do_pdf_clouds,single_gaussian_pdf,dt,SA,errmsg,errflg)

      use machine  , only : kind_phys
      
      implicit none
      
      integer, intent(in) :: idim, kdim
      logical, intent(in) :: do_pdf_clouds, single_gaussian_pdf
      real(kind=kind_phys), intent(in) :: dt
      real(kind=kind_phys), intent(inout) :: SA(:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
!local variables
      real(kind=kind_phys) :: dt_inv

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      !GJF need to call compute_qs_a to calculate atmos_state%qs,gamma
      
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
      !   Constants, Atmos_state, &
            Atmos_state, Input_mp, &
            Cloud_state, SQ, Cloud_processes, &
            SA)
        else
          call tiedtke_macro_pdf (     &
            idim, jdim, kdim, Atmos_state, Input_mp, Cloud_state, ST,  &
            SQ, Cloud_processes, Particles, Lsdiag_mp_control%n_diag_4d,&
            Lsdiag_mp%diag_4d, Lsdiag_mp_control%diag_id,   &
            Lsdiag_mp_control%diag_pt, SA)
        endif
      else  !(do_pdf_clouds)
      endif

      end subroutine tiedtke_prog_clouds_run

      SUBROUTINE tiedtke_macro_Single_Gaussian_pdf (idim, kdim,  Atmos_state,   &
                              Input_mp, &
                              Cloud_state, SQ, Cloud_processes, Particles, &
                              n_diag_4d, diag_4d, diag_id, diag_pt,       SA)  

      !----------------------------------------------------------------------!
      !                                                                      !
      !                     STATISTICAL CLOUD SCHEME                         !
      !                                                                      !
      !-----------------------------------------------------------------------

      integer,                           intent(in)    :: idim, kdim     
      !type(strat_nml_type),              intent(in)    :: Nml
      type(atmos_state_type),            intent(inout) :: Atmos_state !GJF: only need gamma, qs
      type(mp_input_type),               intent(inout) :: Input_mp    !GJF: only need qin
      type(cloud_state_type),            intent(inout) :: Cloud_state !GJF: only need ql_in, qi_in, ql_upd, qi_upd, qa_in, qa_upd
      type(cloud_processes_type),        intent(inout) :: Cloud_processes !GJF: only need qvg, da_ls, delta_cf, dcond_ls
      !type(particles_type),              intent(inout) :: Particles  !GJF: not used
      !type(strat_constants_type),        intent(inout) :: Constants  
      real, dimension(idim,kdim),   intent(in)    :: SQ !, GJF: ST is not used
      real, dimension(idim,kdim),   intent(inout) :: SA
      !integer,                           intent(in)    :: n_diag_4d
      !real, dimension(idim,kdim,0:n_diag_4d),         &
      !                                   intent(inout) :: diag_4d
      !type(diag_id_type),                intent(in)    :: diag_id
      !type(diag_pt_type),                intent(in)    :: diag_pt
      !integer,                           intent(in)    :: otun

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
      !                = stddev_qt /(1 + Atmos_state%gamma ) 


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
      !RSH  change to Input_mp%qin ??
      !       qta(:,:,:) = max(qmin, Atmos_state%qv_in +  &
              qta(:,:,:) = max(qmin, Input_mp%qin +  &
                                           Cloud_state%ql_in + Cloud_state%qi_in)
            else
      !       qta(:,:,:) = max(qmin, Atmos_State%qv_in + SQ +  &
              qta(:,:,:) = max(qmin, Input_mp%qin + SQ +  &
                                           Cloud_state%ql_upd + Cloud_state%qi_upd)
            endif
            qtqsa(:,:,:) = qta(:,:,:) - Atmos_state%qs      

            s =  qtqsa(:,:,:) / (1.0 +  Atmos_state%gamma )

      !---> h1g, 2015-07-23
      ! 99.7% data are within mean +- 3 stdev 
            stdev_s =     qthalfwidth * qta * 0.333
            stdev_s = stdev_s / ( 1.0 +  Atmos_state%gamma ) 

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
                Cloud_processes%qvg(i,k) = (qta(i,k) - qcg(i,k) - qag(i,k)*Atmos_state%qs(i,k))  &
                                              /(max( (1.0 - qag(i,k)), qmin) )
              enddo
            enddo        

            !do adjustment of cloud fraction
            qa0 = Cloud_state%qa_in
            qa1 = qag

            !set total tendency term and update cloud fraction    
            SA   = (SA   + qa1) - qa0
            Cloud_state%qa_upd     = qa1

            !define da_ls and tmp5 needed when do_liq_num = .true. (cjg)
            Cloud_processes%da_ls = max(qa1-qa0,0.)
            Cloud_processes%delta_cf = max(qa1-qa0,0.)

            !compute large-scale condensation / evaporation
            Cloud_processes%dcond_ls = qcg -    &
                                    (Cloud_state%ql_upd + Cloud_state%qi_upd)


      !--->h1g, the following is inherited from nc_cond_pdf, 2015-07-22
            if ( .not. pdf_org ) then
            !!!! INVESTIGATE!!!
            ! make sure super/subsat is not created here
            ! this is different from the original PDF assumption
            ! (as is saturation adjustment) 

      !RSH    Cloud_processes%dcond_ls = MAX( ((Atmos_state%qv_in + SQ -  &
              Cloud_processes%dcond_ls = MAX( ((Input_mp%qin      + SQ -  &
                      Atmos_State%qs)/(1.+Atmos_state%gamma)),    &
                                                     Cloud_processes%dcond_ls)
        
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
    end module tiedtke_prog_clouds
