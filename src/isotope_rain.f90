      subroutine iso_rain
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Update the current-day precipitation delta value from the monthly
!!    table, then mix incoming precipitation into soil layer 1 using
!!    volume-weighted (mass-flux) mixing.
!!
!!    This ports the IsotopeMixer_ (Watson et al.) application to the
!!    precipitation-to-soil-layer-1 pathway in SWAT+.
!!
!!    Called once per HRU per time step, after precipitation is known
!!    (after sq_snom / inflpcp are set) and before iso_lch.

      use hru_module,    only : hru, ihru
      use soil_module
      use climate_module, only : wst, w
      use hydrograph_module, only : ob
      use time_module
      use isotope_module

      implicit none

      external :: iso_mix_binary

      integer :: j    = 0   !none   |HRU number
      integer :: iwst = 0   !none   |weather station index
      integer :: iob  = 0   !none   |object number
      integer :: iiso = 0   !none   |isotope counter
      integer :: imo  = 0   !none   |calendar month (1-12)
      real    :: vol_soil_1 = 0.  !mm  |soil water in layer 1
      real    :: vol_precip = 0.  !mm  |precipitation depth
      real    :: d_precip   = 0.  !permil |isotope delta of precipitation
      real    :: d_soil_old = 0.  !permil |previous soil layer 1 delta
      real    :: iso_mix_binary   !function return value

      j   = ihru
      iob = hru(j)%obj_no
      iwst = ob(iob)%wst
      imo  = time%mo

      !! depth of precipitation on this day
      vol_precip = w%precip

      !! no isotope update if no precipitation
      if (vol_precip <= 0.) return

      !! update current-day precipitation delta from monthly table
      do iiso = 1, num_iso
        iso_precip(iwst, iiso) = iso_precip_mo(iwst, imo, iiso)
        iso_d_rain(j, iiso) = iso_precip(iwst, iiso)
      end do

      !! mix precipitation into soil layer 1
      vol_soil_1 = soil(j)%phys(1)%st

      do iiso = 1, num_iso
        d_precip   = iso_precip(iwst, iiso)
        d_soil_old = iso_soil(j, 1, iiso)
        iso_soil(j, 1, iiso) = iso_mix_binary(d_soil_old, vol_soil_1, &
                                               d_precip,   vol_precip)
      end do

      return
      end subroutine iso_rain
