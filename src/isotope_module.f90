      module isotope_module
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Data structures and global variables for stable water isotope tracking.
!!    Ports J2000 isotope components (Watson et al.) to SWAT+.
!!    Supports delta-18O and delta-D (deuterium) tracking through the
!!    soil-aquifer-channel continuum, including evaporative fractionation
!!    and hydrograph separation diagnostics.
!!
!!    Reference:
!!    Watson A, Vystavna Y, Kralisch S, Helmschrot J, van Rooyen J, Miller J.
!!    Development of an isotope-enabled rainfall-runoff model.
!!    Horita J, Wesolowski DJ (1994) Liquid-vapor fractionation of O and H
!!    isotopes of water from the freezing to the critical temperature.
!!    Geochim Cosmochim Acta 58:3425-3437.

      implicit none

      !! --- simulation control ---
      integer :: iso_on = 0     !! none  |0=off; 1=simulate stable water isotopes
      integer :: num_iso = 1    !! none  |number of isotope tracers (1=d18O; 2=d18O+dD)

      !! --- soil isotope delta values ---
      !! iso_soil(ihru, layer) stores the volume-weighted delta value (permil)
      !! in the soil water of each layer.  Maximum 20 soil layers.
      integer, parameter :: iso_mxlyr = 20
      real, dimension(:,:,:), allocatable :: iso_soil
      !! permil |delta value in soil water (ihru, layer, iso_type)

      !! --- aquifer isotope delta values ---
      !! per-HRU shallow aquifer delta (ihru, iso_type).
      !! Simplified: tracks percolation-weighted average, updated each time
      !! step from incoming percolation and depleted by baseflow discharge.
      real, dimension(:,:), allocatable :: iso_aqu
      !! permil |delta value in shallow aquifer water (ihru, iso_type)

      !! --- precipitation isotope input ---
      !! Daily delta in precipitation per weather station (iwst, iso_type).
      !! Loaded each time step from the precip_isotope monthly table.
      real, dimension(:,:), allocatable :: iso_precip
      !! permil |current-day delta in precipitation (iwst, iso_type)

      !! Monthly precipitation delta per weather station (iwst, month, iso_type).
      !! Read once from precip.iso file during initialization.
      real, dimension(:,:,:), allocatable :: iso_precip_mo
      !! permil |monthly-average delta in precipitation (iwst, 12, iso_type)

      !! --- fractionation parameters (global) ---
      real :: iso_k = 1.0    !! none  |seasonality factor for atmospheric vapor
      real :: iso_x = 0.9    !! none  |soil-water exchange/turnover factor (0-1)

      !! --- calibration floors for binary hydrograph separation ---
      real :: iso_min_comp_rain = 0.  !! frac |minimum rainfall fraction (replaces 0)
      real :: iso_min_comp_gw   = 0.  !! frac |minimum groundwater fraction

      !! --- HRU-level daily isotope flux delta values ---
      !! Volume-weighted delta of each outflow from the HRU soil.
      real, dimension(:,:), allocatable :: iso_d_surq  !! permil (ihru, iso)
      real, dimension(:,:), allocatable :: iso_d_latq  !! permil (ihru, iso)
      real, dimension(:,:), allocatable :: iso_d_tile  !! permil (ihru, iso)
      real, dimension(:,:), allocatable :: iso_d_perc  !! permil (ihru, iso) bottom perc
      real, dimension(:,:), allocatable :: iso_d_evap  !! permil (ihru, iso) evap vapor
      real, dimension(:,:), allocatable :: iso_d_rain  !! permil (ihru, iso) input precip

      !! --- fractionation intermediate variables (per HRU, single isotope) ---
      real, dimension(:), allocatable :: iso_alphamas  !! none  |liquid-vapor eq. frac. factor
      real, dimension(:), allocatable :: iso_epsimas   !! permil|eq. fractionation epsilon
      real, dimension(:), allocatable :: iso_epsk      !! permil|kinetic fractionation epsilon
      real, dimension(:), allocatable :: iso_dstar     !! permil|limiting isotopic composition
      real, dimension(:), allocatable :: iso_enr_slope !! none  |enrichment slope m

      !! --- hydrograph separation diagnostics (per HRU) ---
      !! Binary 2-component separation (rain vs. groundwater):
      real, dimension(:), allocatable :: iso_comp_rain    !! frac  |rain fraction in streamflow
      real, dimension(:), allocatable :: iso_comp_gw      !! frac  |GW fraction in streamflow
      real, dimension(:), allocatable :: iso_comp_rain_n  !! frac  |normalized rain fraction
      real, dimension(:), allocatable :: iso_comp_gw_n    !! frac  |normalized GW fraction
      !! Tertiary 3-component separation (rain / soil-water / groundwater):
      real, dimension(:), allocatable :: iso_comp_a  !! permil*mm |stream isotope × Q_total
      real, dimension(:), allocatable :: iso_comp_b  !! permil*mm |weighted component sum

      !! --- flags / file names ---
      character(len=1) :: iso_atmo = "n"  !! y/n  |read isotope data from precip.iso

      end module isotope_module
