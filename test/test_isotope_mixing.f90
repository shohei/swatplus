!! test_isotope_mixing.f90
!!
!! Standalone Fortran unit tests for iso_mix_binary, iso_mix_array, iso_mix_multi.
!! These three routines have no SWAT+ module dependencies and can be compiled
!! directly with the source file:
!!
!!   gfortran test/test_isotope_mixing.f90 src/isotope_mixing.f90 \
!!            -o build/test_isotope_mixing && ./build/test_isotope_mixing
!!
!! Exit code 0 = all tests passed.  Exit code 1 = at least one failure.

      program test_isotope_mixing

      implicit none

      !! ---- declarations -------------------------------------------------------
      integer  :: n_pass  = 0
      integer  :: n_fail  = 0
      real     :: got, tol

      !! mixing function (external)
      real, external :: iso_mix_binary

      !! arrays for iso_mix_array and iso_mix_multi tests
      real :: va(3), ca(3), vb(3), cb(3)
      real :: vol_a_s, conc_a_s, vol_b_s(3), conc_b_s(3)

      tol = 1.e-5

      write(*,'(a)') "=== isotope mixing unit tests ==="

      !! =========================================================================
      !! iso_mix_binary
      !! =========================================================================

      !! TC-1: equal volumes → average of the two deltas
      got = iso_mix_binary(-5.0, 10.0, -3.0, 10.0)
      call check("binary equal vols     ", got, -4.0, tol, n_pass, n_fail)

      !! TC-2: all water in A, none in B → result = delta_A
      got = iso_mix_binary(-8.0, 50.0, -3.0, 0.0)
      call check("binary zero B vol     ", got, -8.0, tol, n_pass, n_fail)

      !! TC-3: all water in B → result = delta_B
      got = iso_mix_binary(-8.0, 0.0, -3.0, 50.0)
      call check("binary zero A vol     ", got, -3.0, tol, n_pass, n_fail)

      !! TC-4: both volumes zero → fall back to delta_A
      got = iso_mix_binary(-7.0, 0.0, -2.0, 0.0)
      call check("binary both zero      ", got, -7.0, tol, n_pass, n_fail)

      !! TC-5: mass conservation; vol_A=30, vol_B=10 → weighted mean
      !! expected = (-6*30 + -2*10) / 40 = (-180 - 20)/40 = -5.0
      got = iso_mix_binary(-6.0, 30.0, -2.0, 10.0)
      call check("binary mass conserv   ", got, -5.0, tol, n_pass, n_fail)

      !! TC-6: same delta → result is unchanged regardless of volumes
      got = iso_mix_binary(-4.5, 100.0, -4.5, 1.0)
      call check("binary same delta     ", got, -4.5, tol, n_pass, n_fail)

      !! TC-7: positive delta values
      !! expected = (2*20 + 8*5) / 25 = (40+40)/25 = 3.2
      got = iso_mix_binary(2.0, 20.0, 8.0, 5.0)
      call check("binary positive delta ", got, 3.2, tol, n_pass, n_fail)


      !! =========================================================================
      !! iso_mix_array (bidirectional, prop=1)
      !! =========================================================================

      !! TC-8: n=2, bidir=true, prop=1, equal vols → both updated to average
      !! A = [-6, -4], B = [-2, 0], vol_A = vol_B = [10, 10]
      !! expected x = [(-6*10 + -2*10)/20, (-4*10 + 0*10)/20] = [-4, -2]
      va  = (/ 10., 10., 0. /)
      ca  = (/ -6., -4., 0. /)
      vb  = (/ 10., 10., 0. /)
      cb  = (/ -2., 0., 0.  /)
      call iso_mix_array(2, va, ca, vb, cb, .true., 1.0)
      call check("array bidir ca(1)     ", ca(1), -4.0, tol, n_pass, n_fail)
      call check("array bidir ca(2)     ", ca(2), -2.0, tol, n_pass, n_fail)
      call check("array bidir cb(1)     ", cb(1), -4.0, tol, n_pass, n_fail)
      call check("array bidir cb(2)     ", cb(2), -2.0, tol, n_pass, n_fail)

      !! TC-9: bidir=false → only B is updated
      va  = (/ 20., 20., 0. /)
      ca  = (/ -6., -4., 0. /)
      vb  = (/ 10., 10., 0. /)
      cb  = (/ -2., 0., 0.  /)
      call iso_mix_array(2, va, ca, vb, cb, .false., 1.0)
      !! ca should be unchanged
      call check("array unidir ca(1)    ", ca(1), -6.0, tol, n_pass, n_fail)
      call check("array unidir ca(2)    ", ca(2), -4.0, tol, n_pass, n_fail)
      !! cb: x = (ca*va + cb*vb/1)/(va+vb) = (−6*20 + −2*10)/30 = −160/30 ≈ −4.667
      call check("array unidir cb(1)    ", cb(1), -14.0/3.0, tol, n_pass, n_fail)

      !! TC-10: prop=0 → skip (guard branch); arrays unchanged
      va  = (/ 10., 10., 0. /)
      ca  = (/ -5., -5., 0. /)
      vb  = (/ 10., 10., 0. /)
      cb  = (/ -3., -3., 0. /)
      call iso_mix_array(2, va, ca, vb, cb, .true., 0.0)
      call check("array prop=0 guard    ", cb(1), -3.0, tol, n_pass, n_fail)

      !! TC-11: prop=0.5 → vol_eff = vol_b / 0.5 = 2*vol_b
      !! n=1, va=10, ca=-6; vb=10, prop=0.5 → vol_eff=20
      !! x = (-6*10 + -2*20)/(10+20) = (-60-40)/30 = -3.333
      va(1) = 10.; ca(1) = -6.
      vb(1) = 10.; cb(1) = -2.
      call iso_mix_array(1, va, ca, vb, cb, .true., 0.5)
      call check("array prop=0.5        ", cb(1), -10.0/3.0, tol, n_pass, n_fail)


      !! =========================================================================
      !! iso_mix_multi
      !! =========================================================================

      !! TC-12: 3 equal-volume destinations → each gets 1/3 of source
      !! source: vol_a=30, conc_a=-6
      !! dest: vol_b=[10,10,10], conc_b=[0,0,0]
      !! weight_i = 1/3; vol_a_w = 10 for each
      !! x_i = (-6*10 + 0*10)/(10+10) = -3 for each
      !! conc_a_new = -3 (weighted avg of identical x_i)
      vol_a_s   = 30.;  conc_a_s  = -6.
      vol_b_s   = (/ 10., 10., 10. /)
      conc_b_s  = (/  0.,  0.,  0. /)
      call iso_mix_multi(vol_a_s, conc_a_s, 3, vol_b_s, conc_b_s)
      call check("multi eq dest cb(1)   ", conc_b_s(1), -3.0, tol, n_pass, n_fail)
      call check("multi eq dest cb(2)   ", conc_b_s(2), -3.0, tol, n_pass, n_fail)
      call check("multi eq dest cb(3)   ", conc_b_s(3), -3.0, tol, n_pass, n_fail)
      call check("multi eq dest ca new  ", conc_a_s,    -3.0, tol, n_pass, n_fail)

      !! TC-13: all dest vols zero → source unchanged (guard branch)
      vol_a_s   = 10.;  conc_a_s  = -5.
      vol_b_s   = (/ 0., 0., 0. /)
      conc_b_s  = (/ -3., -3., -3. /)
      call iso_mix_multi(vol_a_s, conc_a_s, 3, vol_b_s, conc_b_s)
      call check("multi zero dest guard ", conc_a_s, -5.0, tol, n_pass, n_fail)

      !! TC-14: two unequal destinations; mass conservation
      !! source: vol_a=20, conc_a=-8
      !! dest: vol_b=[10,30], conc_b=[-2,-4]
      !! total_b = 40; w1=0.25, w2=0.75
      !! vol_a_w1 = 20*0.25=5, vol_a_w2 = 20*0.75=15
      !! x1 = (-8*5 + -2*10)/(5+10)  = (-40-20)/15 = -4.0
      !! x2 = (-8*15 + -4*30)/(15+30) = (-120-120)/45 = -5.333
      !! conc_a_new = x1*0.25 + x2*0.75 = -4*0.25 + -5.333*0.75 = -1 + -4 = -5.0
      vol_a_s   = 20.;  conc_a_s  = -8.
      vol_b_s   = (/ 10., 30., 0. /)
      conc_b_s  = (/ -2., -4., 0. /)
      call iso_mix_multi(vol_a_s, conc_a_s, 2, vol_b_s, conc_b_s)
      call check("multi unequal cb(1)   ", conc_b_s(1), -4.0,       tol, n_pass, n_fail)
      call check("multi unequal cb(2)   ", conc_b_s(2), -16.0/3.0,  tol, n_pass, n_fail)
      call check("multi unequal ca new  ", conc_a_s,    -5.0,        tol, n_pass, n_fail)

      !! =========================================================================
      !! Summary
      !! =========================================================================
      write(*,'(a)') repeat("-", 44)
      write(*,'(a,i3,a,i3,a)') "RESULT: ", n_pass, " passed, ", n_fail, " failed"

      if (n_fail > 0) then
        stop 1
      end if

      contains

      subroutine check(label, got, expected, tol, n_pass, n_fail)
        character(len=*), intent(in) :: label
        real,    intent(in)  :: got, expected, tol
        integer, intent(inout) :: n_pass, n_fail
        real :: diff
        diff = abs(got - expected)
        if (diff <= tol) then
          write(*,'(a,a,a)') "  PASS  ", label
          n_pass = n_pass + 1
        else
          write(*,'(a,a,a,f12.6,a,f12.6,a,e10.3)') &
            "  FAIL  ", label, "  got=", got, "  exp=", expected, "  diff=", diff
          n_fail = n_fail + 1
        end if
      end subroutine check

      end program test_isotope_mixing
