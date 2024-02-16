!>\file UFS_rad_time_vary.fv3.F90
!!  Contains code related to GFS radiation suite setup (radiation part of time_vary_step)
   module UFS_rad_time_vary

      implicit none

      private
      
      !> new data input control variables (set/reset in subroutine radupdate):
      integer :: month0 = 0
      integer :: iyear0 = 0
      integer :: monthd = 0
      
      !> control flag for the first time of reading climatological ozone data
      !! (set/reset in subroutines radinit/radupdate, it is used only if the
      !! control parameter ntoz=0)
      logical :: loz1st = .true.
      
      public UFS_rad_time_vary_timestep_init

      contains

!>\defgroup mod_UFS_rad_time_vary GFS Radiation Time Update
!! This module contains code related to GFS radiation setup.
!> @{
!> \section arg_table_UFS_rad_time_vary_timestep_init Argument Table
!! \htmlinclude UFS_rad_time_vary_timestep_init.html
!!
      subroutine UFS_rad_time_vary_timestep_init (me, idate, jdate, ictm, isol, iaermdl, ico2, ntoz, deltsw, deltim, con_pi, aeros_file, co2dat_file, co2glb_file, ozphys, lrseeds, rseeds,                     &
              lslwr, lsswr, isubc_lw, isubc_sw, icsdsw, icsdlw, cnx, cny, isc, jsc,    &
              imap, jmap, sec, kdt, imp_physics, imp_physics_zhao_carr, ipsd0, ipsdlim,&
              ps_2delt, ps_1delt, t_2delt, t_1delt, qv_2delt, qv_1delt, t, qv, ps,     &
              slag, sdec, cdec, solcon, errmsg, errflg)

         use module_radiation_astronomy, only : sol_update
         use module_radiation_aerosols,  only : aer_update
         use module_radiation_gases,     only : gas_update
         use module_ozphys,              only : ty_ozphys
         use mersenne_twister,           only : random_setseed, random_index, random_stat
         use machine,                    only : kind_phys
         use radcons,                    only : qmin, con_100

         implicit none

         ! Interface variables
         integer,                intent(in)    :: idate(:)
         integer,                intent(in)    :: jdate(:)
         logical,                intent(in)    :: lrseeds
         integer,                intent(in)    :: rseeds(:,:)
         integer,                intent(in)    :: me, ictm, isol, iaermdl, ico2, ntoz
         integer,                intent(in)    :: isubc_lw, isubc_sw, cnx, cny, isc, jsc, kdt
         integer,                intent(in)    :: imp_physics, imp_physics_zhao_carr, ipsd0, ipsdlim
         logical,                intent(in)    :: lslwr, lsswr
         integer,                intent(inout) :: icsdsw(:), icsdlw(:)
         integer,                intent(in)    :: imap(:), jmap(:)
         character(len=26),      intent(in)    :: aeros_file, co2dat_file, co2gbl_file
         real(kind_phys),        intent(in)    :: deltsw, deltim, con_pi
         real(kind_phys),        intent(in)    :: sec
         real(kind_phys),        intent(inout) :: ps_2delt(:)
         real(kind_phys),        intent(inout) :: ps_1delt(:)
         real(kind_phys),        intent(inout) :: t_2delt(:,:)
         real(kind_phys),        intent(inout) :: t_1delt(:,:)
         real(kind_phys),        intent(inout) :: qv_2delt(:,:)
         real(kind_phys),        intent(inout) :: qv_1delt(:,:)
         type(ty_ozphys),        intent(inout) :: ozphys
         real(kind_phys),        intent(in)    :: t(:,:), qv(:,:), ps(:)
         real(kind_phys),        intent(out)   :: slag
         real(kind_phys),        intent(out)   :: sdec
         real(kind_phys),        intent(out)   :: cdec
         real(kind_phys),        intent(out)   :: solcon
         character(len=*),       intent(out)   :: errmsg
         integer,                intent(out)   :: errflg

         ! Local variables
         integer :: iyear, imon, iday, ihour
         integer :: kyear, kmon, kday, khour
         logical :: lmon_chg       ! month change flag
         logical :: lco2_chg       ! cntrl flag for updating co2 data
         logical :: lsol_chg       ! cntrl flag for updating solar constant
         
         type (random_stat) :: stat
         integer :: ix, j, i, ipseed
         integer :: numrdm(cnx*cny*2)

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0
         
         ! Set up time stamp at fcst time and that for green house gases
         iyear = jdate(1)
         imon  = jdate(2)
         iday  = jdate(3)
         ihour = jdate(5)
         
         ! Set up time stamp used for green house gases (** currently co2 only)
         ! get external data at initial condition time 
         if ( ictm==0 .or. ictm==-2 ) then 
            kyear = idate(1)
            kmon  = idate(2)
            kday  = idate(3)
            khour = idate(5)
         ! get external data at fcst or specified time 
         else
            kyear = iyear
            kmon  = imon
            kday  = iday
            khour = ihour
         endif       
         
         if ( month0 /= imon ) then
           lmon_chg = .true.
           month0   = imon
         else
           lmon_chg = .false.
         endif
         
         !> -# Call module_radiation_astronomy::sol_update(), yearly update, no
         !! time interpolation.
         if (lsswr) then
           if ( isol == 0 .or. isol == 10 ) then
             lsol_chg = .false.
           elseif ( iyear0 /= iyear ) then
             lsol_chg = .true.
           else
             lsol_chg = ( isol==4 .and. lmon_chg )
           endif
           iyear0 = iyear

           call sol_update(jdate,kyear,deltsw,deltim,lsol_chg,me,    & !inputs
              slag,sdec,cdec,solcon,con_pi,errmsg,errflg)              !outputs
         endif 
         
         !> -# Call module_radiation_aerosols::aer_update(), monthly update, no
         !! time interpolation
         if ( lmon_chg ) then
           call aer_update (iyear, imon, me, iaermdl, aeros_file, errflg, errmsg)
         endif

         !> -# Call co2 and other gases update routine:
         !! module_radiation_gases::gas_update()
         if ( monthd /= kmon ) then
           monthd = kmon
           lco2_chg = .true.
         else
           lco2_chg = .false.
         endif

         call gas_update ( kyear,kmon,kday,khour,lco2_chg, me, co2dat_file, &
            co2gbl_file, ictm, ico2, errflg, errmsg )
         if (ntoz == 0) then
           call ozphys%update_o3clim(kmon, kday, khour, loz1st)
         endif

         if ( loz1st ) loz1st = .false.
         
         if (lsswr .or. lslwr) then

           !--- call to GFS_radupdate_timestep_init is now in GFS_rrtmg_setup_timestep_init

           !--- set up random seed index in a reproducible way for entire cubed-sphere face (lat-lon grid)
           if ((isubc_lw==2) .or. (isubc_sw==2)) then
             !NRL If random seeds supplied by NEPTUNE
             if(lrseeds) then
               do ix=1,size(jmap)
                 icsdsw(ix) = rseeds(ix,1)
                 icsdlw(ix) = rseeds(ix,2)
               enddo
             else
               ipseed = mod(nint(con_100*sqrt(sec)), ipsdlim) + 1 + ipsd0
               call random_setseed (ipseed, stat)
               call random_index (ipsdlim, numrdm, stat)

               do ix=1,size(jmap)
                 j = jmap(ix)
                 i = imap(ix)
                 !--- for testing purposes, replace numrdm with '100'
                 icsdsw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx)
                 icsdlw(ix) = numrdm(i+isc-1 + (j+jsc-2)*cnx + cnx*cny)
               enddo
             end if !lrseeds
           endif  ! isubc_lw and isubc_sw

           if (imp_physics == imp_physics_zhao_carr) then
             if (kdt == 1) then
               t_2delt  = t
               t_1delt  = t
               qv_2delt = qv
               qv_1delt = qv
               ps_2delt = ps
               ps_1delt = ps
             endif
           endif

         endif

      end subroutine UFS_rad_time_vary_timestep_init

!> @}

   end module UFS_rad_time_vary
