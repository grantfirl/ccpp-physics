!> \file clouds.F90
!!  Contains the Tiedtke prognostic cloud scheme

!> This module contains the CCPP-compliant Tiedtke prognostic cloud scheme.
      module clouds

      contains

!> \section arg_table_clouds_init Argument Table
!! \htmlinclude clouds_init.html
!!
      subroutine clouds_init(errmsg, errflg)
         use machine, only : kind_phys
         
         implicit none

         character(len=*),     intent(out) :: errmsg
         integer,              intent(out) :: errflg
         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

      end subroutine clouds_init

!> \section arg_table_clouds_run Argument Table
!! \htmlinclude clouds_run.html
!!
      subroutine clouds_run (kdim, dtp, dt_inner, tgrs, prsl, omega, spechum, qc, qa, dtdt_rad_turb_gwd, dcond_ls, errmsg, errflg)

      use machine, only : kind_phys
      
      implicit none
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Input
      integer, intent(in) :: kdim
      real(kind_phys), intent(in) :: dtp, dt_inner
      real(kind_phys), intent(in) :: prsl(:,:), omega(:,:), dtdt_rad_turb_gwd(:,:)
      
      ! Input / Output
      real(kind_phys), intent(inout) :: spechum(:,:), qc(:,:), qa(:,:), tgrs(:,:), dcond_ls(:,:)

      ! Local
      integer :: i, k
      integer :: idim = 1
      real, allocatable :: dqs_ls(:,:), qv(:,:), rho(:,:), qsl(:,:), mc_full(:,:), dqsdT(:,:), U(:,:), U00(:,:), tmp1(:,:), gamm(:,:), da_ls(:,:)
      real, allocatable :: A_dt(:,:), B_dt(:,:), qa0(:,:), qa1(:,:), qabar(:,:), qaeq(:,:), D_eros(:,:), lvap(:,:), qa_tend(:,:)
      real, allocatable :: ocp(:,:)
      real :: con_eps = 0.608
      real :: con_rd = 287.
      real :: grav = 9.8
      real :: cp_air = 1004.
      real :: rvgas = 461.50
      real :: hlv = 2.5e6
      logical, allocatable :: do_subgrid_clouds(:,:)
      
      allocate(dqs_ls(1,kdim))
      allocate(qv(1,kdim))
      allocate(rho(1,kdim))
      allocate(qsl(1,kdim))
      allocate(qa_tend(1,kdim))
      allocate(mc_full(1,kdim))
      allocate(dqsdT(1,kdim))
      allocate(U(1,kdim))
      allocate(U00(1,kdim))
      allocate(tmp1(1,kdim))
      allocate(gamm(1,kdim))
      allocate(da_ls(1,kdim))
      allocate(A_dt(1,kdim))
      allocate(B_dt(1,kdim))
      allocate(qa0(1,kdim))
      allocate(qa1(1,kdim))
      allocate(qabar(1,kdim))
      allocate(qaeq(1,kdim))
      allocate(D_eros(1,kdim))
      allocate(ocp(1,kdim))
      allocate(lvap(1,kdim))
      allocate(do_subgrid_clouds(1,kdim))
      
! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      write(*,*) 'IN THE CLOUDS', kdim, dtp, dt_inner
      
      
      do k = 1, kdim
         do i = 1, idim
            ! Process rates
            dcond_ls(i,k) = 0.
            qa_tend(i,k) = 0.

            ! Env vars
            qv(i,k) = spechum(i,k) / (1.0_kind_phys - spechum(i,k))
            rho(i,k) = con_eps * prsl(i,k) / (con_rd*tgrs(i,k)*(qv(i,k)+con_eps))
            lvap(i,k) = hlv + (2106.0 - 4218.0)*(tgrs(i,k)-273.15)
            ocp(i,k) = 1./(cp_air*(1.+0.887*qv(i,k)))

            ! More env vars
            qsl(i,k) = get_qsl(tgrs(i,k), prsl(i,k))
            dqsdT(i,k) = lvap(i,k) * qsl(i,k) / (rvgas*tgrs(i,k)**2)      
            gamm(i,k) = dqsdT(i,k) * lvap(i,k) * ocp(i,k)
            mc_full(i,k) = 0.0
            U(i,k) = qv(i,k) / qsl(i,k)
            U00(i,k) = 0.8

            ! Flag for subgrid
            do_subgrid_clouds(i,k) = .false.

            qa(i,k) = min(max(qa(i,k),0.),1.)
            
            ! Make/destroy subgrid clouds where RH < 100%, temperature > -20C, cloud fraction < 1
            ! Otherwise, let microphysics handle things
            if(U(i,k).le.1..and.tgrs(i,k).gt.253.15) then
               do_subgrid_clouds(i,k) = .true.
            endif

!            if(U(i,k).gt.1.002.and.qc(i,k).gt.1.e-4) then
!               qa(i,k) = 1.
!            endif
            
            if (qc(i,k).le.1.e-12) then
               qa(i,k) = 0.
            endif

         enddo
      enddo

      do k = 1, kdim
         do i = 1, idim
            if (do_subgrid_clouds(i,k)) then
               
               ! check on qc
               if(qc(i,k).le.1.e-12) then
                  qv(i,k) = qv(i,k) + qc(i,k)
                  tgrs(i,k) = tgrs(i,k) - lvap(i,k)*ocp(i,k)*qc(i,k)
                  qc(i,k) = 0.
                  qa(i,k) = 0.
               else
                  if(qa(i,k).lt.0.001) then
                     qa(i,k) = 0.001
                  endif
               endif

               ! check on qv
               if(qv(i,k).lt.1.e-10) then
                  qv(i,k) = 1.e-10
               endif

               ! can we assume that in-cloud value are coming in (change here to grid values)
!               qc(i,k) = qc(i,k)*qa(i,k)

               ! Calculation
               dqs_ls(i,k) = (((omega(i,k) + grav*mc_full(i,k))/rho(i,k)*ocp(i,k))+dtdt_rad_turb_gwd(i,k))*dt_inner*dqsdT(i,k)
               if (dqs_ls(i,k).lt.0..and.U(i,k).ge.U00(i,k).and.qa(i,k).lt.1.)  then
                  tmp1(i,k) = sqrt( (1. + qa(i,k)*gamm(i,k))**2. - (1. - qa(i,k))*(1. - qa(i,k))*gamm(i,k)*dqs_ls(i,k)/qsl(i,k)/max(1.-U(i,k), 1.e-12)) - (1. + qa(i,k)* gamm(i,k))
                  tmp1(i,k) = -1.*tmp1(i,k)/((1. - qa(i,k))*(1. - qa(i,k))*gamm(i,k)/qsl(i,k)/max(1. - U(i,k), 1.e-12)/2.)
                  dqs_ls(i,k) = min(tmp1(i,k),dqs_ls(i,k)/(1. +0.5*(1. + qa(i,k))* gamm(i,k)))
               else
                  dqs_ls(i,k) = dqs_ls(i,k)/(1. + qa(i,k)*gamm(i,k))
               endif
               
               if ((dqs_ls(i,k).le.0..and.U(i,k).ge.U00(i,k)).and.(qa(i,k).le.1.)) then
                  da_ls(i,k) = -0.5*(1. - qa(i,k))*(1. - qa(i,k))* dqs_ls(i,k)/qsl(i,k) / max(1.-U(i,k), 1.e-12)
               else
                  da_ls(i,k) = 0.
               endif

               ! Erosion
               D_eros(i,k) = 0.
               if (qc(i,k).gt.1.e-12) then
                  D_eros(i,k) = -1.e-6*qa(i,k)*(qsl(i,k)-qv(i,k))*qa(i,k) / qc(i,k)
               endif
               
               A_dt(i,k) = da_ls(i,k) / max((1.-qa(i,k)), 1.e-12)
               B_dt(i,k) = 1.e-5
!!!               B_dt(i,k) = 0.
               
               !------------------------------------------------------------------------
               !    do analytic integration.      
               !------------------------------------------------------------------------
               if ( (A_dt(i,k) .gt. 1.e-12) .or. (B_dt(i,k) .gt.1.e-12 ) )  then
                  qa0(i,k) = qa(i,k)
                  qaeq(i,k) = A_dt(i,k)/(A_dt(i,k)  + B_dt(i,k))
                  qa1(i,k) = qaeq(i,k) - (qaeq(i,k) - qa0(i,k)) * exp ( -1.*(A_dt(i,k)+B_dt(i,k)) )
                  qabar(i,k) = qaeq(i,k) - ((qa1(i,k) - qa0(i,k))/ (A_dt(i,k) + B_dt(i,k)))
               else
                  qa0(i,k)   = qa(i,k)
                  qaeq(i,k)  = qa0(i,k)
                  qa1(i,k)   = qa0(i,k)   
                  qabar(i,k) = qa0(i,k)  
               endif

               ! Add tendencies where doing subgrid clouds
               qa_tend(i,k) = (qa1(i,k)-qa0(i,k)) / dt_inner
               dcond_ls(i,k) = -1. * qabar(i,k) * dqs_ls(i,k)

               ! Limit loss of condensation based on available cloud water
               if ((-1.*dcond_ls(i,k)*dt_inner).gt.qc(i,k)) then
                  dcond_ls(i,k) = -qc(i,k)/dt_inner
               endif

               ! Likely needed cloud fraction destruction
!               if (dcond_ls(i,k).lt.0..and.qc(i,k).gt.1.e-12) then
!                  qa_tend(i,k) = dcond_ls(i,k)/qc(i,k)*qa(i,k)
!               endif
               
               qa(i,k) = qa(i,k) + qa_tend(i,k)*dt_inner
               qc(i,k) = qc(i,k) + dcond_ls(i,k)*dt_inner
               qv(i,k) = qv(i,k) - dcond_ls(i,k)*dt_inner
               tgrs(i,k) = tgrs(i,k) + lvap(i,k)*ocp(i,k)*dcond_ls(i,k)*dt_inner

               ! Send in-cloud values out
!               if (qa(i,k).gt.0.001) then
!                  qc(i,k) = qc(i,k)/qa(i,k)
!               else
!                  qa(i,k) = 0.
!               endif

               if (qc(i,k).le.1.e-12) then
                  qa(i,k) = 0.
               endif

            else
               if (qc(i,k).le.1.e-12) then
                  qa(i,k) = 0.
               endif
            endif
         enddo
      enddo

      !-----------------------------------------------------------------------
    end subroutine clouds_run

    real function get_qsl(temp, pres)
      implicit none
      
      real, intent(in)  :: temp, pres
      real            :: x, esl
      real, parameter :: c0= .611583699E03
      real, parameter :: c1= .444606896E02
      real, parameter :: c2= .143177157E01
      real, parameter :: c3= .264224321E-1
      real, parameter :: c4= .299291081E-3
      real, parameter :: c5= .203154182E-5
      real, parameter :: c6= .702620698E-8
      real, parameter :: c7= .379534310E-11
      real, parameter :: c8=-.321582393E-13
      
      x   = MAX(-80.,temp-273.16)
      esl = c0 + x*(c1 + x* (c2 + x*(c3 + x*(c4 + x*(c5 + x*(c6 + x*(c7 + x*c8)))))))
      esl = MIN(esl, pres*0.15)        ! Even with P=1050mb and T=55C, the sat. vap. pres only contributes to ~15% of total pres.                                                                                                            
      get_qsl = .622*esl/max(1.e-4,(pres - esl))
      
    end function get_qsl
  
  end module clouds
