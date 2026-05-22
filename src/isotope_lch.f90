      subroutine iso_lch
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Compute the isotope delta of water leaving each soil layer via
!!    surface runoff (layer 1), lateral flow, tile drainage, and
!!    percolation.  Update soil layer deltas after each flux removes
!!    water.  Also mix incoming aquifer water (rls_routeaqu pathway)
!!    back into the bottom soil layer.
!!
!!    The same volume-weighted mixing approach is used as in
!!    IsotopeMixer_ (Watson et al.), applied layer-by-layer.
!!
!!    Delta of any outflow = delta of the layer from which water leaves,
!!    consistent with the complete-mixing (well-mixed) reservoir assumption
!!    used throughout J2000-ISO and most transit-time models.
!!
!!    Called after iso_rain, before iso_frac.

      use hru_module,    only : hru, ihru, surfq, latq, qtile, sepbtm
      use soil_module
      use isotope_module

      implicit none

      external :: iso_mix_binary

      integer :: j    = 0    !none    |HRU number
      integer :: jj   = 0    !none    |soil layer counter
      integer :: nly  = 0    !none    |number of soil layers
      integer :: iiso = 0    !none    |isotope counter
      real :: iso_mix_binary  !function

      real :: vol_surq  = 0.  !mm  |surface runoff (from layer 1)
      real :: vol_latq  = 0.  !mm  |lateral flow from current layer
      real :: vol_tile  = 0.  !mm  |tile drainage from layer
      real :: vol_perc  = 0.  !mm  |percolation from current layer
      real :: vol_total = 0.  !mm  |total water leaving layer
      real :: vol_remain= 0.  !mm  |soil water after fluxes
      real :: d_layer   = 0.  !permil |delta of current layer
      real :: d_in      = 0.  !permil |delta of water entering layer from above

      j   = ihru
      nly = soil(j)%nly

      do iiso = 1, num_iso

        !! carry-over percolation delta from layer above (starts as rain delta)
        d_in = iso_soil(j, 1, iiso)

        do jj = 1, nly

          d_layer = iso_soil(j, jj, iiso)

          !! incoming percolation from layer above mixes into this layer
          if (jj > 1) then
            vol_perc = soil(j)%ly(jj-1)%prk
            if (vol_perc > 0.) then
              iso_soil(j, jj, iiso) = iso_mix_binary(d_layer,           &
                                            soil(j)%phys(jj)%st,        &
                                            d_in, vol_perc)
              d_layer = iso_soil(j, jj, iiso)
            end if
          end if

          !! determine fluxes leaving this layer
          vol_surq = 0.
          vol_tile = 0.
          if (jj == 1) vol_surq = surfq(j)
          if (hru(j)%lumv%ldrain == jj) vol_tile = qtile
          vol_latq = soil(j)%ly(jj)%flat
          vol_perc = soil(j)%ly(jj)%prk

          !! delta of each outflow = current layer delta (well-mixed assumption)
          if (jj == 1) iso_d_surq(j, iiso) = d_layer
          iso_d_latq(j, iiso) = iso_d_latq(j, iiso) + vol_latq * d_layer

          if (hru(j)%lumv%ldrain == jj) then
            iso_d_tile(j, iiso) = d_layer
          end if

          !! delta carried to next layer by percolation
          d_in = d_layer

          !! update layer delta after removing fluxes (mass balance)
          vol_total  = vol_surq + vol_latq + vol_tile + vol_perc
          vol_remain = soil(j)%phys(jj)%st - vol_total
          if (vol_remain < 0.) vol_remain = 0.

          !! layer delta unchanged when well-mixed: all water has same delta.
          !! delta is preserved; only the volume changes (tracked externally
          !! by soil module).  No update needed here for d_layer itself.

        end do ! layer loop

        !! bottom percolation delta = delta of bottom soil layer
        iso_d_perc(j, iiso) = iso_soil(j, nly, iiso)

        !! volume-weight lateral flow delta across all layers
        if (latq(j) > 0.) then
          iso_d_latq(j, iiso) = iso_d_latq(j, iiso) / latq(j)
        else
          iso_d_latq(j, iiso) = iso_soil(j, 1, iiso)
        end if

        !! update shallow aquifer delta by mixing with incoming percolation
        !! aqu_stor_approx = sepbtm(j) * 10 as a very rough proxy for
        !! aquifer storage; replaced with exact aqu_d%stor when available.
        if (sepbtm(j) > 0.) then
          block
            real :: aqu_stor_proxy
            aqu_stor_proxy = max(sepbtm(j) * 10., 1.)
            iso_aqu(j, iiso) = iso_mix_binary(iso_aqu(j, iiso),         &
                                               aqu_stor_proxy,           &
                                               iso_d_perc(j, iiso),     &
                                               sepbtm(j))
          end block
        end if

      end do ! isotope loop

      return
      end subroutine iso_lch
