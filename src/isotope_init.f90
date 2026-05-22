      subroutine iso_init
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Allocate and initialize all arrays required by the isotope tracking
!!    module.  Read initial delta values from precip.iso (if present).
!!
!!    Called once during the model initialization sequence, after HRU,
!!    soil, and aquifer arrays have been allocated.
!!
!!    Input file format  (precip.iso):
!!    Line 1 : title (ignored)
!!    Line 2 : "iso_on num_iso k x min_comp_rain min_comp_gw"
!!    Lines 3+: "wst_name  d18O_jan d18O_feb ... d18O_dec  [dD_jan ... dD_dec]"
!!    Soil and aquifer initial delta = annual-average precip delta of
!!    the nearest weather station (simplification; can be overridden).

      use hydrograph_module, only : sp_ob
      use soil_module
      use climate_module,    only : wst
      use maximum_data_module, only : db_mx
      use isotope_module
      use input_file_module

      implicit none

      integer :: ihru   = 0
      integer :: iaqu   = 0
      integer :: iwst   = 0
      integer :: imo    = 0
      integer :: iiso   = 0
      integer :: ly     = 0
      integer :: eof    = 0
      integer :: mhru   = 0
      integer :: maqu   = 0
      integer :: mwst   = 0
      real    :: d18O_mo(12) = 0.
      real    :: dD_mo(12)   = 0.
      real    :: d18O_aa     = 0.
      real    :: dD_aa       = 0.
      character(len=25) :: wst_nm = " "
      character(len=80) :: titldum = " "
      logical :: i_exist

      !! nothing to do if isotopes are off
      if (iso_on == 0) return

      mhru = sp_ob%hru
      maqu = mhru               !! one per-HRU aquifer proxy
      mwst = db_mx%wst

      !! allocate soil isotope arrays
      allocate (iso_soil(mhru, iso_mxlyr, num_iso), source = 0.)
      allocate (iso_aqu(maqu, num_iso),              source = 0.)
      allocate (iso_precip(mwst, num_iso),           source = 0.)
      allocate (iso_precip_mo(mwst, 12, num_iso),    source = 0.)

      !! allocate HRU flux arrays
      allocate (iso_d_surq(mhru, num_iso), source = 0.)
      allocate (iso_d_latq(mhru, num_iso), source = 0.)
      allocate (iso_d_tile(mhru, num_iso), source = 0.)
      allocate (iso_d_perc(mhru, num_iso), source = 0.)
      allocate (iso_d_evap(mhru, num_iso), source = 0.)
      allocate (iso_d_rain(mhru, num_iso), source = 0.)

      !! fractionation intermediates (single-isotope scalars per HRU)
      allocate (iso_alphamas(mhru),  source = 1.)
      allocate (iso_epsimas(mhru),   source = 0.)
      allocate (iso_epsk(mhru),      source = 0.)
      allocate (iso_dstar(mhru),     source = 0.)
      allocate (iso_enr_slope(mhru), source = 1.)

      !! hydrograph separation arrays
      allocate (iso_comp_rain(mhru),   source = 0.)
      allocate (iso_comp_gw(mhru),     source = 0.)
      allocate (iso_comp_rain_n(mhru), source = 0.)
      allocate (iso_comp_gw_n(mhru),   source = 0.)
      allocate (iso_comp_a(mhru),      source = 0.)
      allocate (iso_comp_b(mhru),      source = 0.)

      !! read precip.iso if it exists
      inquire (file="precip.iso", exist=i_exist)
      if (.not. i_exist) return

      open (unit=119, file="precip.iso", status="old")
      read (119,*,iostat=eof) titldum
      if (eof < 0) then
        close (119); return
      endif
      !! second line: global parameters
      read (119,*,iostat=eof) iso_on, num_iso, iso_k, iso_x,  &
                               iso_min_comp_rain, iso_min_comp_gw
      if (eof < 0) then
        close (119); return
      endif

      !! read monthly delta values per weather station
      do iwst = 1, mwst
        if (num_iso >= 1) then
          read (119,*,iostat=eof) wst_nm, (d18O_mo(imo), imo=1,12)
          if (eof < 0) exit
          do imo = 1, 12
            iso_precip_mo(iwst, imo, 1) = d18O_mo(imo)
          end do
          d18O_aa = sum(d18O_mo) / 12.
        end if
        if (num_iso >= 2) then
          read (119,*,iostat=eof) wst_nm, (dD_mo(imo), imo=1,12)
          if (eof < 0) exit
          do imo = 1, 12
            iso_precip_mo(iwst, imo, 2) = dD_mo(imo)
          end do
          dD_aa = sum(dD_mo) / 12.
        end if

        !! initialize current-day precip delta to January value
        do iiso = 1, num_iso
          iso_precip(iwst, iiso) = iso_precip_mo(iwst, 1, iiso)
        end do
      end do
      close (119)

      !! initialize soil and aquifer delta to annual-average precip delta
      !! of weather station 1 (simple default; user can provide a separate
      !! init file later)
      if (mwst > 0) then
        do iiso = 1, num_iso
          d18O_aa = sum(iso_precip_mo(1, :, iiso)) / 12.
          do ihru = 1, mhru
            do ly = 1, iso_mxlyr
              iso_soil(ihru, ly, iiso) = d18O_aa
            end do
            iso_aqu(ihru, iiso) = d18O_aa
          end do
        end do
      end if

      return
      end subroutine iso_init
