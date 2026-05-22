      subroutine iso_hydsep
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Perform isotope-based hydrograph separation for the current HRU.
!!    Implements both:
!!
!!    (1) Binary 2-component separation (Watson et al. 2022):
!!        f_rain = (delta_stream - delta_gw)  / (delta_rain - delta_gw)
!!        f_gw   = (delta_stream - delta_rain) / (delta_gw   - delta_rain)
!!
!!        The corrected (standard) mixing equations are used here.
!!        The Java Binary_isotope_mixing code contains a precedence ambiguity;
!!        this implementation follows the textbook 2-component IHS formula
!!        (Sklash & Farvolden 1979; Pinder & Jones 1969).
!!
!!    (2) Tertiary 3-component proportioning (Watson et al. 2022):
!!        Computes the "measured" and "simulated" isotope flux loads.
!!
!!        comp_A = delta_stream * Q_total
!!                 (isotope load from observed stream, used for validation)
!!
!!        When precipitation delta is available (iso_d_rain != -99):
!!          comp_B = delta_rain * Q_surface
!!                 + delta_gw   * Q_baseflow
!!                 + delta_sw   * Q_lateral
!!        When no rain isotope data:
!!          comp_B = delta_gw   * Q_baseflow
!!                 + delta_sw   * Q_lateral
!!
!!        (Note: the Java Tertiary_Isotope_mixing has an inverted if-branch;
!!        this implementation uses the logically correct interpretation.)
!!
!!    Inputs used:
!!      iso_d_surq(j,1)  – simulated delta in surface runoff (≈ delta_rain proxy)
!!      iso_aqu(j,1)     – simulated delta in aquifer (≈ delta_gw)
!!      iso_soil(j,*,1)  – simulated delta in soil water (≈ monthly delta_sw)
!!      iso_d_rain(j,1)  – precipitation delta for current time step
!!      surfq(j), latq(j), sepbtm(j) – flow components (mm)
!!
!!    Results written to iso_comp_rain(j), iso_comp_gw(j), etc.

      use hru_module,    only : hru, ihru, surfq, latq, sepbtm, qdr
      use soil_module
      use isotope_module

      implicit none

      integer :: j      = 0    !none   |HRU number
      integer :: nly    = 0    !none   |number of soil layers
      real :: d_stream  = 0.   !permil |simulated stream isotope delta
      real :: d_rain    = 0.   !permil |precipitation delta
      real :: d_gw      = 0.   !permil |aquifer (groundwater) delta
      real :: d_sw      = 0.   !permil |soil-water delta (volume-weighted profile)
      real :: q_surf    = 0.   !mm     |surface runoff
      real :: q_lat     = 0.   !mm     |lateral / interflow
      real :: q_base    = 0.   !mm     |baseflow (approx from sepbtm)
      real :: q_total   = 0.   !mm     |total HRU water yield
      real :: f_rain    = 0.   !frac   |binary rain fraction
      real :: f_gw      = 0.   !frac   |binary GW fraction
      real :: frac_avail = 0.  !frac   |1 - comp_SW (available non-interflow)
      real :: denom     = 0.   !permil |denominator of binary equations
      real :: vol_sw    = 0.   !mm     |total soil water
      real :: wt_sum    = 0.   !mm     |weight sum for soil average
      integer :: jj     = 0    !none   |soil layer counter

      j   = ihru
      nly = soil(j)%nly

      !! --- get flow components ---
      q_surf  = surfq(j)
      q_lat   = latq(j)
      q_base  = sepbtm(j)          ! percolation as proxy for baseflow
      q_total = q_surf + q_lat + q_base

      if (q_total <= 0.) then
        iso_comp_rain(j)   = 0.
        iso_comp_gw(j)     = 0.
        iso_comp_rain_n(j) = 0.
        iso_comp_gw_n(j)   = 0.
        iso_comp_a(j)      = 0.
        iso_comp_b(j)      = 0.
        return
      end if

      !! --- get isotope end-member deltas ---
      d_rain = iso_d_rain(j, 1)      ! today's precipitation delta
      d_gw   = iso_aqu(j, 1)        ! aquifer delta

      !! volume-weighted average soil water delta across all layers
      vol_sw  = 0.
      wt_sum  = 0.
      do jj = 1, nly
        vol_sw = soil(j)%phys(jj)%st
        wt_sum = wt_sum + iso_soil(j, jj, 1) * vol_sw
        vol_sw = vol_sw + vol_sw   ! accumulate total
      end do
      vol_sw = 0.
      do jj = 1, nly
        vol_sw = vol_sw + soil(j)%phys(jj)%st
      end do
      if (vol_sw > 0.) then
        d_sw = wt_sum / vol_sw
      else
        d_sw = iso_soil(j, 1, 1)
      end if

      !! simulated stream delta: volume-weighted mix of all flow components
      if (q_total > 0.) then
        d_stream = (iso_d_surq(j,1) * q_surf                             &
                  + iso_d_latq(j,1) * q_lat                              &
                  + d_gw            * q_base) / q_total
      end if

      !! ================================================================
      !! (1) Binary hydrograph separation
      !!     Standard 2-component IHS (Sklash & Farvolden 1979)
      !! ================================================================
      denom = d_rain - d_gw
      if (abs(denom) > 1.e-6) then
        f_rain = (d_stream - d_gw)  / denom
        f_gw   = (d_stream - d_rain) / (-denom)
      else
        f_rain = 0.
        f_gw   = 1.
      end if

      !! apply calibration floors
      if (f_rain <= 0.) f_rain = iso_min_comp_rain
      if (f_gw   <= 0.) f_gw   = iso_min_comp_gw

      !! normalized fractions (account for simulated soil-water component)
      frac_avail = 1. - (q_lat / max(q_total, 1.e-10))
      iso_comp_rain(j)   = f_rain
      iso_comp_gw(j)     = f_gw
      iso_comp_rain_n(j) = f_rain * frac_avail
      iso_comp_gw_n(j)   = f_gw  * frac_avail

      !! ================================================================
      !! (2) Tertiary proportioning
      !!     comp_A = simulated isotope load (stream × total Q)
      !!     comp_B = sum of component contributions
      !! ================================================================
      iso_comp_a(j) = d_stream * q_total

      !! rain data missing flag: use -99 to indicate no isotope in precip
      if (d_rain > -98.) then
        !! full 3-component: surface (rain) + baseflow (GW) + lateral (soil)
        iso_comp_b(j) = d_rain * q_surf                                   &
                      + d_gw   * q_base                                   &
                      + d_sw   * q_lat
      else
        !! 2-component only (no rain isotope data available)
        iso_comp_b(j) = d_gw * q_base + d_sw * q_lat
      end if

      return
      end subroutine iso_hydsep
