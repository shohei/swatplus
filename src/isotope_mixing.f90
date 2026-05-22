!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    Utility functions for stable-isotope mixing.
!!    Ports IsotopeMixer_ (binary sequential mixing), IsotopeMixer
!!    (array bidirectional mixing with proportion), and IsotopeMultiMixer
!!    (one-to-many proportional mixing) from J2000-ISO (Watson et al.).
!!
!!    All functions use volume-weighted (mass-flux) mixing:
!!      delta_mixed = (delta_A * vol_A + delta_B * vol_B) / (vol_A + vol_B)
!!
!!    Reference:
!!    Watson A, Birkel C, Kralisch S (2023).  IsotopeMixer components
!!    for J2000-ISO, FSU Jena / Stellenbosch University.


      !! -----------------------------------------------------------------
      !! iso_mix_binary
      !!   Binary sequential mixing of two water volumes.
      !!   Returns the mass-flux weighted delta of the combined water.
      !!   Equivalent to IsotopeMixer_ (Watson et al.).
      !!
      !!   delta_a  [permil] – delta of existing storage (vol_a)
      !!   vol_a    [mm]     – volume of existing storage
      !!   delta_b  [permil] – delta of incoming water (vol_b)
      !!   vol_b    [mm]     – volume of incoming water
      !! -----------------------------------------------------------------
      real function iso_mix_binary (delta_a, vol_a, delta_b, vol_b)

      implicit none
      real, intent(in) :: delta_a, vol_a, delta_b, vol_b
      real :: total_vol

      total_vol = vol_a + vol_b
      if (total_vol > 0.) then
        iso_mix_binary = (delta_a * vol_a + delta_b * vol_b) / total_vol
      else
        iso_mix_binary = delta_a
      end if

      return
      end function iso_mix_binary


      !! -----------------------------------------------------------------
      !! iso_mix_array
      !!   Bidirectional array mixing with a mixing proportion.
      !!   Equivalent to IsotopeMixer (Watson, Birkel, Kralisch 2023).
      !!
      !!   For each element i:
      !!     vol_eff  = vol_b(i) / prop          (effective volume of B)
      !!     x        = (dA(i)*vol_A(i) + dB(i)*vol_eff) / (vol_A(i)+vol_eff)
      !!     if bidir: dA(i) = x
      !!               dB(i) = x
      !!
      !!   n     – array length
      !!   vol_a – volumes of compartment A (mm)  [intent in]
      !!   conc_a– delta of compartment A (permil)[intent inout if bidir]
      !!   vol_b – volumes of compartment B (mm)  [intent in]
      !!   conc_b– delta of compartment B (permil)[intent inout]
      !!   bidir – .true. = update both A and B; .false. = update B only
      !!   prop  – mixing proportion (0 < prop <= 1)
      !! -----------------------------------------------------------------
      subroutine iso_mix_array (n, vol_a, conc_a, vol_b, conc_b, bidir, prop)

      implicit none
      integer,  intent(in)    :: n
      real,     intent(in)    :: vol_a(n), vol_b(n), prop
      real,     intent(inout) :: conc_a(n), conc_b(n)
      logical,  intent(in)    :: bidir

      integer :: i
      real    :: vol_eff, x_mix, total_vol

      if (prop <= 0.) return

      do i = 1, n
        vol_eff   = vol_b(i) / prop
        total_vol = vol_a(i) + vol_eff
        if (total_vol > 0.) then
          x_mix = (conc_a(i) * vol_a(i) + conc_b(i) * vol_eff) / total_vol
        else
          x_mix = conc_a(i)
        end if
        if (bidir) conc_a(i) = x_mix
        conc_b(i) = x_mix
      end do

      return
      end subroutine iso_mix_array


      !! -----------------------------------------------------------------
      !! iso_mix_multi
      !!   One-to-many mixing: distribute one source (A) into N
      !!   destinations (B(1..n)), weighted by the volume fraction of each
      !!   B relative to the total B volume.
      !!   Equivalent to IsotopeMultiMixer (Kralisch 2023).
      !!
      !!   Algorithm:
      !!     1. weight(i)  = vol_b(i) / sum(vol_b)
      !!     2. vol_a_w(i) = vol_a * weight(i)
      !!     3. x(i) = (dA * vol_a_w(i) + dB(i)*vol_b(i)) / (vol_a_w(i)+vol_b(i))
      !!     4. dA   = sum(x(i) * weight(i))   (volume-weighted average)
      !!     5. dB(i)= x(i)
      !! -----------------------------------------------------------------
      subroutine iso_mix_multi (vol_a, conc_a, n, vol_b, conc_b)

      implicit none
      real,     intent(in)    :: vol_a
      real,     intent(inout) :: conc_a
      integer,  intent(in)    :: n
      real,     intent(in)    :: vol_b(n)
      real,     intent(inout) :: conc_b(n)

      integer :: i
      real    :: total_vol, weight, vol_a_w, x_mix, conc_a_new

      total_vol = 0.
      do i = 1, n
        total_vol = total_vol + vol_b(i)
      end do
      if (total_vol <= 0.) return

      conc_a_new = 0.
      do i = 1, n
        weight  = vol_b(i) / total_vol
        vol_a_w = vol_a * weight
        if (vol_a_w + vol_b(i) > 0.) then
          x_mix = (conc_a * vol_a_w + conc_b(i) * vol_b(i)) /  &
                  (vol_a_w + vol_b(i))
        else
          x_mix = conc_a
        end if
        conc_a_new = conc_a_new + x_mix * weight
        conc_b(i)  = x_mix
      end do
      conc_a = conc_a_new

      return
      end subroutine iso_mix_multi
