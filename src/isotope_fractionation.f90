      subroutine iso_frac
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Compute liquid-vapor equilibrium and kinetic isotopic fractionation
!!    during soil-surface evaporation, then update the soil layer-1 delta
!!    value and compute the delta of the evaporation vapor.
!!
!!    Ports Isotope_fractionation (Watson, Birkel, Kralisch 2023):
!!
!!    1. Equilibrium fractionation factor alpha_plus (Horita & Wesolowski 1994)
!!       for delta-18O:
!!         alpha = exp(1/1000 * (1158.8*T^3/1e9 - 1620.1*T^2/1e6 +
!!                               794.84*T/1e3 - 161.04 + 2.9992e9/T^3))
!!       where T is temperature in Kelvin.
!!
!!    2. Kinetic fractionation (Merlivat 1978 / Gibson et al. 2016):
!!         eps_k = 0.9755 * (1 - 0.9755) * 1000 * (1 - rh)
!!
!!    3. Enrichment slope m and limiting composition delta* (Gibson 2016):
!!         m      = (rh - 1e-3*(eps_k + eps_mas/alpha)) /
!!                  (1 - rh + 1e-3*eps_k)
!!         delta_A = (delta_P - k*eps_mas) / (1 + eps_mas*1e-3)
!!         delta*  = (rh*dA + eps_k + eps_mas/alpha) /
!!                   (rh - 1e-3*(eps_k + eps_mas/alpha))
!!
!!    4. Residual soil-water delta after evaporation (Craig & Gordon 1965):
!!         delta_S = delta_S0 - delta* * (1-x)^m + delta*
!!
!!    5. Evaporation vapor delta:
!!         delta_E = ((delta_S - eps_mas)/alpha - rh*dA - eps_k) /
!!                   (1 - rh + 1e-3*eps_k)
!!
!!    Only layer 1 is updated (surface evaporation).  Transpiration is
!!    assumed non-fractionating (Evaristo et al. 2015).
!!
!!    NOTE: alpha formula above is for delta-18O.  For delta-D a different
!!    set of coefficients applies; extend num_iso and this subroutine.

      use hru_module,    only : hru, ihru, es_day
      use climate_module, only : w
      use isotope_module

      implicit none

      integer :: j    = 0   !none   |HRU number
      real :: tempC   = 0.  !degC   |daily average air temperature
      real :: tempK   = 0.  !K      |temperature in Kelvin
      real :: rh      = 0.  !frac   |relative humidity (0-1)
      real :: alphamas = 0. !none   |liquid-vapor equilibrium fractionation
      real :: epsimas  = 0. !permil |equilibrium epsilon
      real :: epsk_H   = 0. !permil |kinetic fractionation epsilon
      real :: enr_slp  = 0. !none   |enrichment slope m
      real :: d_atm    = 0. !permil |atmospheric vapor delta (delta_A)
      real :: dstar    = 0. !permil |limiting isotopic composition (delta*)
      real :: d_soil   = 0. !permil |residual soil water delta
      real :: d_evap   = 0. !permil |evaporation vapor delta
      real :: d_precip = 0. !permil |precipitation delta (proxy for d_A calc)
      real :: es_frac  = 0. !frac   |fraction of soil evaporation from layer 1
      real :: denom    = 0. !permil |denominator in delta* calculation

      j = ihru

      !! skip if fractionation is not meaningful
      if (iso_on == 0 .or. es_day <= 0.) return

      tempC = w%tave
      rh    = w%rhum      ! SWAT+ stores rhum as a fraction (0-1)

      !! guard against extreme conditions
      rh    = max(0.01, min(0.99, rh))
      tempC = max(-40., min(100., tempC))
      tempK = tempC + 273.15

      !! --- equilibrium fractionation (Horita & Wesolowski 1994, Eq.2) ---
      !! valid for delta-18O; T in Kelvin
      alphamas = exp((1. / 1000.) *                                       &
                     (1158.8 * tempK**3 / 1.e9                            &
                    - 1620.1 * tempK**2 / 1.e6                            &
                    +  794.84 * tempK   / 1.e3                            &
                    -  161.04                                              &
                    + 2.9992e9 / tempK**3))

      !! epsilon (‰ notation)
      epsimas = (alphamas - 1.) * 1000.

      !! --- kinetic fractionation (Merlivat 1978) ---
      !! diffusivity ratio D(H2-18O)/D(H2-16O) ≈ 0.9755 (laminar-turbulent)
      epsk_H = 0.9755 * (1. - 0.9755) * 1000. * (1. - rh)

      !! --- enrichment slope m (Gibson et al. 2016) ---
      denom = 1. - rh + 1.e-3 * epsk_H
      if (abs(denom) < 1.e-10) denom = 1.e-10
      enr_slp = (rh - 1.e-3 * (epsk_H + epsimas / alphamas)) / denom

      !! --- atmospheric vapor delta from precip (Gibson et al. 2008) ---
      !! delta_A = (delta_P - k * eps_mas) / (1 + eps_mas * 1e-3)
      d_precip = iso_precip(ob_wst_now(j), 1)  ! delta-18O of current precip
      d_atm    = (d_precip - iso_k * epsimas) / (1. + epsimas * 1.e-3)

      !! --- limiting isotopic composition delta* (Gonfiantini 1986) ---
      denom = rh - 1.e-3 * (epsk_H + epsimas / alphamas)
      if (abs(denom) < 1.e-10) denom = 1.e-10
      dstar = (rh * d_atm + epsk_H + epsimas / alphamas) / denom

      !! --- residual soil-water delta after evaporation ---
      d_soil = iso_soil(j, 1, 1) - dstar * (1. - iso_x)**enr_slp + dstar

      !! --- evaporation vapor delta (Craig & Gordon 1965) ---
      denom = 1. - rh + 1.e-3 * epsk_H
      if (abs(denom) < 1.e-10) denom = 1.e-10
      d_evap = ((d_soil - epsimas) / alphamas - rh * d_atm - epsk_H) / denom

      !! update soil layer 1 and output variables
      iso_soil(j, 1, 1)  = d_soil
      iso_d_evap(j, 1)   = d_evap
      iso_alphamas(j)    = alphamas
      iso_epsimas(j)     = epsimas
      iso_epsk(j)        = epsk_H
      iso_dstar(j)       = dstar
      iso_enr_slope(j)   = enr_slp

      return

      contains

      !! helper: get weather station index for HRU j
      integer function ob_wst_now(jhru)
        use hru_module, only : hru
        use hydrograph_module, only : ob
        integer, intent(in) :: jhru
        ob_wst_now = ob(hru(jhru)%obj_no)%wst
      end function ob_wst_now

      end subroutine iso_frac
