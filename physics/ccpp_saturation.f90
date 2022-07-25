module ccpp_saturation

  use machine, only: kind_phys
  
  implicit none
  save
  private
  
  real(kind=kind_phys), parameter :: min_blend_temp = 253.16
  real(kind=kind_phys), parameter :: max_blend_temp = 273.16
  real(kind=kind_phys) :: hlv, hlf, rvgas, cp_air
    
  public get_qs_init, get_qs 
  
  contains
    
  subroutine get_qs_init(hlv_in, hlf_in, rvgas_in, cp_air_in)
    real(kind=kind_phys), intent(in)  :: hlv_in, hlf_in, rvgas_in, cp_air_in
    
    hlv = hlv_in
    hlf = hlf_in
    rvgas = rvgas_in
    cp_air = cp_air_in
  end subroutine get_qs_init
    
  elemental subroutine get_qs(temp, pres, alg_choice, alg_flatau_92, qs, qsl, qsi, es, esl, esi, dqsdT, gamma)
     real(kind=kind_phys), intent(in)  :: temp, pres
     integer             , intent(in)  :: alg_choice, alg_flatau_92
     real(kind=kind_phys), intent(out) :: qs, qsl, qsi, es, esl, esi, dqsdT, gamma
     
     real(kind=kind_phys) :: hls
     real(kind=kind_phys) :: lifrac, dqsdT_l, dqsdT_i, gamma_l, gamma_i
      
     hls = hlv + hlf
     
     if (alg_choice == alg_flatau_92) then
        call get_qs_flatau_92_over_liquid(temp, pres, esl, qsl)
        call get_qs_flatau_92_over_ice(temp, pres, esi, qsi)
     else
       qs = -999.9
       RETURN
     end if
     
     dqsdT_l = hlv*qsl/(rvgas*temp**2)
     dqsdT_i = hls*qsi/(rvgas*temp**2)
     gamma_l = dqsdT_l*hlv/cp_air
     gamma_i = dqsdT_i*hls/cp_air
     
     !blend qsl, qsi => qs; esl, esi => es, etc.
     if ( temp <= min_blend_temp ) then
       qs = qsi
       es = esi
       dqsdT = dqsdT_i
       gamma = gamma_i
     else if ( temp >= max_blend_temp ) then
       qs = qsl
       es = esl
       dqsdT = dqsdT_l
       gamma = gamma_l
     else
       lifrac = (temp - min_blend_temp) / (max_blend_temp - min_blend_temp)
       qs = qsi + lifrac*(qsl - qsi)
       es = esi + lifrac*(esl - esi)
       dqsdT = dqsdT_i + lifrac*(dqsdT_l - dqsdT_i)
       gamma = gamma_i + lifrac*(gamma_l - gamma_i)
     end if
    
  end subroutine get_qs
  
  elemental subroutine get_qs_flatau_92_over_liquid(temp, pres, esl, qsl)
    real(kind=kind_phys), intent(in)  :: temp, pres
    real(kind=kind_phys), intent(out) :: esl, qsl
    
    !code obtained from module_mp_thompson.F90 in ccpp-physics
    
    real(kind=kind_phys)            :: x
    real(kind=kind_phys), parameter :: c0= .611583699E03
    real(kind=kind_phys), parameter :: c1= .444606896E02
    real(kind=kind_phys), parameter :: c2= .143177157E01
    real(kind=kind_phys), parameter :: c3= .264224321E-1
    real(kind=kind_phys), parameter :: c4= .299291081E-3
    real(kind=kind_phys), parameter :: c5= .203154182E-5
    real(kind=kind_phys), parameter :: c6= .702620698E-8
    real(kind=kind_phys), parameter :: c7= .379534310E-11
    real(kind=kind_phys), parameter :: c8=-.321582393E-13
    
    x   = MAX(-80.,temp-273.16)
    esl = c0 + x*(c1 + x* (c2 + x*(c3 + x*(c4 + x*(c5 + x*(c6 + x*(c7 + x*c8)))))))
    esl = MIN(esl, pres*0.15)        ! Even with P=1050mb and T=55C, the sat. vap. pres only contributes to ~15% of total pres.
    qsl = .622*esl/max(1.e-4,(pres - esl))
    
  end subroutine get_qs_flatau_92_over_liquid
  
  elemental subroutine get_qs_flatau_92_over_ice(temp, pres, esi, qsi)
    real(kind=kind_phys), intent(in)  :: temp, pres
    real(kind=kind_phys), intent(out) :: esi, qsi
    
    !code obtained from module_mp_thompson.F90 in ccpp-physics
    
    real(kind=kind_phys)            :: x
    real(kind=kind_phys), parameter :: c0= .609868993E03
    real(kind=kind_phys), parameter :: c1= .499320233E02
    real(kind=kind_phys), parameter :: c2= .184672631E01
    real(kind=kind_phys), parameter :: c3= .402737184E-1
    real(kind=kind_phys), parameter :: c4= .565392987E-3
    real(kind=kind_phys), parameter :: c5= .521693933E-5
    real(kind=kind_phys), parameter :: c6= .307839583E-7
    real(kind=kind_phys), parameter :: c7= .105785160E-9
    real(kind=kind_phys), parameter :: c8= .161444444E-12
    
    x   = MAX(-80.,temp-273.16)
    esi = c0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + x*(c6 + x*(c7 + x*c8)))))))
    esi = MIN(esi, pres*0.15)
    qsi =.622*esi/max(1.e-4,(pres-esi))
    
  end subroutine get_qs_flatau_92_over_ice

end module ccpp_saturation
