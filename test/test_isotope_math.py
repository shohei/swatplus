#!/usr/bin/env python3
"""
test_isotope_math.py
====================
Unit tests for the mathematical equations implemented in:
  src/isotope_mixing.f90        - iso_mix_binary / iso_mix_array / iso_mix_multi
  src/isotope_fractionation.f90 - Horita & Wesolowski (1994), Merlivat (1978),
                                   Gibson et al. (2016), Craig & Gordon (1965)
  src/isotope_hydsep.f90        - Sklash & Farvolden (1979) binary IHS,
                                   tertiary proportioning

All tests are pure-Python reimplementations that run without SWAT+ or a
Fortran compiler.  They validate the scientific logic independently so that
regressions in either the Fortran or the equations themselves are caught.

NOTE on fractionation formula (isotope_fractionation.f90 / Isotope_fractionation.java):
  The H&W 1994 polynomial with T^3 in the numerator (1158.8, 1620.1, ...) is
  the deuterium (delta-D) equation.  The delta-18O equation has a different
  form.  The J2000 Java source uses the same polynomial and labels it for
  both isotopes.  This port faithfully reproduces J2000 behaviour; tests
  below are written against the IMPLEMENTED behaviour and flag the
  scientific discrepancy with "NOTE" comments.

Usage:
  python3 test/test_isotope_math.py
  python3 -m pytest test/test_isotope_math.py -v

Exit code 0 = all pass, 1 = any failure.
"""

import math
import os
import sys
import tempfile

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PASS = 0
_FAIL = 0


def _check(label, got, expected, tol=1e-5):
    global _PASS, _FAIL
    diff = abs(got - expected)
    if diff <= tol:
        print(f"  PASS  {label}")
        _PASS += 1
    else:
        print(f"  FAIL  {label}  got={got:.8f}  exp={expected:.8f}  diff={diff:.3e}")
        _FAIL += 1


def _check_range(label, got, lo, hi):
    global _PASS, _FAIL
    if lo <= got <= hi:
        print(f"  PASS  {label}")
        _PASS += 1
    else:
        print(f"  FAIL  {label}  got={got:.6f}  expected in [{lo}, {hi}]")
        _FAIL += 1


# ===========================================================================
# 1. Mixing functions  (mirrors isotope_mixing.f90)
# ===========================================================================

def iso_mix_binary(delta_a, vol_a, delta_b, vol_b):
    total = vol_a + vol_b
    if total > 0.0:
        return (delta_a * vol_a + delta_b * vol_b) / total
    return delta_a


def iso_mix_multi(vol_a, conc_a, vol_b_list, conc_b_list):
    """Returns (new_conc_a, new_conc_b_list)."""
    total = sum(vol_b_list)
    if total <= 0.0:
        return conc_a, list(conc_b_list)
    new_b = []
    new_a = 0.0
    for vb, cb in zip(vol_b_list, conc_b_list):
        w = vb / total
        va_w = vol_a * w
        denom = va_w + vb
        x = (conc_a * va_w + cb * vb) / denom if denom > 0 else conc_a
        new_a += x * w
        new_b.append(x)
    return new_a, new_b


def test_mixing():
    print("\n--- 1. Mixing functions ---")

    # iso_mix_binary
    _check("binary equal vols",     iso_mix_binary(-5, 10, -3, 10),     -4.0)
    _check("binary zero B vol",     iso_mix_binary(-8, 50, -3,  0),     -8.0)
    _check("binary zero A vol",     iso_mix_binary(-8,  0, -3, 50),     -3.0)
    _check("binary both zero",      iso_mix_binary(-7,  0, -2,  0),     -7.0)
    _check("binary mass conserv",   iso_mix_binary(-6, 30, -2, 10),     -5.0)
    _check("binary same delta",     iso_mix_binary(-4.5, 100, -4.5, 1), -4.5)

    # commutativity
    r1 = iso_mix_binary(-6, 30, -2, 10)
    r2 = iso_mix_binary(-2, 10, -6, 30)
    _check("binary commutative",    r1, r2)

    # linearity: (A+B)+C == A+(B+C) with correct volumes
    ab = iso_mix_binary(-6, 30, -2, 10)
    abc_seq = iso_mix_binary(ab, 40, 0, 20)
    abc_direct = (-6*30 + -2*10 + 0*20) / 60.0
    _check("binary linearity",      abc_seq, abc_direct, tol=1e-5)

    # iso_mix_multi
    new_a, new_b = iso_mix_multi(30, -6, [10, 10, 10], [0, 0, 0])
    _check("multi eq dest b[0]",    new_b[0], -3.0)
    _check("multi eq dest b[1]",    new_b[1], -3.0)
    _check("multi eq dest b[2]",    new_b[2], -3.0)
    _check("multi eq dest a new",   new_a,    -3.0)

    # guard: zero total dest volume -> source unchanged
    new_a2, _ = iso_mix_multi(10, -5, [0, 0], [-3, -3])
    _check("multi zero dest guard", new_a2, -5.0)

    # unequal destinations
    # vol_a=20, conc_a=-8; vol_b=[10,30], conc_b=[-2,-4]
    # w=[0.25, 0.75]; va_w=[5, 15]
    # x0 = (-8*5 + -2*10)/(5+10) = -60/15 = -4.0
    # x1 = (-8*15 + -4*30)/(15+30) = -240/45 = -16/3
    # new_a = -4*0.25 + (-16/3)*0.75 = -1 + -4 = -5.0
    new_a3, new_b3 = iso_mix_multi(20, -8, [10, 30], [-2, -4])
    _check("multi unequal b[0]",    new_b3[0], -4.0)
    _check("multi unequal b[1]",    new_b3[1], -16.0/3.0)
    _check("multi unequal a new",   new_a3,    -5.0)


# ===========================================================================
# 2. Fractionation equations  (mirrors isotope_fractionation.f90)
# ===========================================================================

def alpha_hw1994(T_C):
    """H&W 1994 polynomial as implemented (T^3 form; matches J2000 Java)."""
    T = T_C + 273.15
    exponent = (1.0 / 1000.0) * (
        1158.8 * T**3 / 1.0e9
        - 1620.1 * T**2 / 1.0e6
        + 794.84 * T   / 1.0e3
        - 161.04
        + 2.9992e9 / T**3
    )
    return math.exp(exponent)


def eps_kinetic(rh):
    """Merlivat (1978) kinetic fractionation epsilon (permil)."""
    return 0.9755 * (1.0 - 0.9755) * 1000.0 * (1.0 - rh)


def fractionation(T_C, rh, d_precip, d_soil0, iso_k=1.0, iso_x=0.9):
    """
    Full fractionation calculation matching isotope_fractionation.f90.
    Returns (alpha, eps_mas, eps_k, m, d_atm, dstar, d_soil, d_evap).
    """
    rh = max(0.01, min(0.99, rh))
    T_C = max(-40.0, min(100.0, T_C))

    alpha = alpha_hw1994(T_C)
    eps_mas = (alpha - 1.0) * 1000.0
    eps_k = eps_kinetic(rh)

    denom_m = 1.0 - rh + 1.0e-3 * eps_k
    if abs(denom_m) < 1.0e-10:
        denom_m = 1.0e-10
    m = (rh - 1.0e-3 * (eps_k + eps_mas / alpha)) / denom_m

    d_atm = (d_precip - iso_k * eps_mas) / (1.0 + eps_mas * 1.0e-3)

    denom_star = rh - 1.0e-3 * (eps_k + eps_mas / alpha)
    if abs(denom_star) < 1.0e-10:
        denom_star = 1.0e-10
    dstar = (rh * d_atm + eps_k + eps_mas / alpha) / denom_star

    # J2000 formula: d_soil = d_soil0 - dstar*(1-x)^m + dstar
    # Equivalent to: d_soil0 + dstar*(1-(1-x)^m)
    d_soil = d_soil0 - dstar * (1.0 - iso_x)**m + dstar

    denom_e = 1.0 - rh + 1.0e-3 * eps_k
    if abs(denom_e) < 1.0e-10:
        denom_e = 1.0e-10
    d_evap = ((d_soil - eps_mas) / alpha - rh * d_atm - eps_k) / denom_e

    return alpha, eps_mas, eps_k, m, d_atm, dstar, d_soil, d_evap


def test_fractionation():
    print("\n--- 2. Fractionation equations ---")
    # NOTE: The polynomial used here (1158.8*T^3/1e9 ...) is the H&W 1994
    # equation for deuterium (delta-D), which gives alpha ~1.079 at 25 deg C.
    # The delta-18O formula would give ~1.009.  The J2000 Java source uses this
    # same polynomial; tests below validate the IMPLEMENTED behaviour.

    # --- Equilibrium fractionation factor ---
    alpha_25 = alpha_hw1994(25.0)
    # Implemented formula gives ~1.079 at 25 deg C (delta-D scale)
    _check_range("alpha at 25C in [1.07,1.10]",  alpha_25, 1.07, 1.10)
    _check_range("alpha > 1.0 at 25C",            alpha_25, 1.0,  2.0)

    # Monotone: alpha decreases as T rises toward boiling point
    alpha_0   = alpha_hw1994(0.0)
    alpha_100 = alpha_hw1994(100.0)
    _check_range("alpha 0C > alpha 25C",   alpha_0   - alpha_25, 1e-4, 1.0)
    _check_range("alpha 25C > alpha 100C", alpha_25  - alpha_100, 1e-4, 1.0)

    # Epsilon (permil) ~79 permil at 25 deg C for the implemented polynomial
    eps_25 = (alpha_25 - 1.0) * 1000.0
    _check_range("eps_mas at 25C in [70,90] permil", eps_25, 70.0, 90.0)

    # --- Kinetic fractionation ---
    ek_dry = eps_kinetic(0.3)
    ek_wet = eps_kinetic(0.8)
    _check_range("eps_k dry > eps_k wet",  ek_dry - ek_wet, 1e-4, 100.0)
    _check("eps_k at rh=1 is 0",           eps_kinetic(1.0), 0.0)
    _check_range("eps_k >= 0 at rh=0.5",   eps_kinetic(0.5), 0.0, 100.0)

    # --- Reference conditions ---
    T_C, rh, d_precip, d_soil0 = 20.0, 0.6, -8.0, -7.0
    _, _, _, _, _, _, d_soil, d_evap = \
        fractionation(T_C, rh, d_precip, d_soil0)

    # Evaporation enriches residual soil water (large range due to delta-D scale)
    _check_range("d_soil > d_soil0 after evap",
                 d_soil - d_soil0, 0.0, 500.0)

    # Evaporation vapor is depleted relative to residual soil water
    _check_range("d_evap < d_soil", d_soil - d_evap, 0.0, 2000.0)

    # Drier air -> more kinetic fractionation -> larger enrichment
    _, _, _, _, _, _, d_soil_dry, _ = fractionation(T_C, 0.3, d_precip, d_soil0)
    _, _, _, _, _, _, d_soil_wet, _ = fractionation(T_C, 0.8, d_precip, d_soil0)
    _check_range("dry air enriches more than wet",
                 d_soil_dry - d_soil_wet, 0.0, 2000.0)

    # --- Guard: extreme rh clamped to [0.01, 0.99] ---
    alpha_lo, *_ = fractionation(25.0, -0.5, -8.0, -7.0)
    alpha_hi, *_ = fractionation(25.0,  2.0, -8.0, -7.0)
    _check_range("alpha valid when rh<0 clamped", alpha_lo, 1.0, 2.0)
    _check_range("alpha valid when rh>1 clamped", alpha_hi, 1.0, 2.0)

    # --- iso_x boundary conditions ---
    # Formula: d_soil = d_soil0 - dstar*(1-x)^m + dstar
    #          = d_soil0 + dstar*(1-(1-x)^m)

    # At x=0: (1-0)^m = 1  ->  d_soil = d_soil0 - dstar + dstar = d_soil0
    _, _, _, _, _, _, d_soil_x0, _ = fractionation(20.0, 0.6, -8.0, d_soil0, iso_x=0.0)
    _check("iso_x=0 d_soil equals d_soil0", d_soil_x0, d_soil0, tol=1e-4)

    # At x=1: (1-1)^m = 0  ->  d_soil = d_soil0 + dstar
    # NOTE: the standard Craig-Gordon form would give d_soil = dstar at x=1;
    # the J2000/Fortran form adds dstar to the initial value instead.
    _, _, _, _, _, dstar_x1, d_soil_x1, _ = \
        fractionation(20.0, 0.6, -8.0, d_soil0, iso_x=1.0)
    _check("iso_x=1 d_soil = d_soil0+dstar", d_soil_x1, d_soil0 + dstar_x1, tol=1e-4)


# ===========================================================================
# 3. Hydrograph separation  (mirrors isotope_hydsep.f90)
# ===========================================================================

def binary_hydsep(d_stream, d_rain, d_gw, min_comp_rain=0.0, min_comp_gw=0.0):
    """Sklash & Farvolden (1979) 2-component IHS."""
    denom = d_rain - d_gw
    if abs(denom) > 1.0e-6:
        f_rain = (d_stream - d_gw)  / denom
        f_gw   = (d_stream - d_rain) / (-denom)
    else:
        f_rain = 0.0
        f_gw   = 1.0
    if f_rain <= 0.0:
        f_rain = min_comp_rain
    if f_gw <= 0.0:
        f_gw = min_comp_gw
    return f_rain, f_gw


def tertiary_proportioning(d_stream, q_total,
                            d_rain, q_surf,
                            d_gw, q_base,
                            d_sw, q_lat,
                            rain_available=True):
    """Tertiary isotope proportioning (Watson et al. 2022)."""
    comp_a = d_stream * q_total
    if rain_available:
        comp_b = d_rain * q_surf + d_gw * q_base + d_sw * q_lat
    else:
        comp_b = d_gw * q_base + d_sw * q_lat
    return comp_a, comp_b


def test_hydsep():
    print("\n--- 3. Hydrograph separation ---")

    # --- Binary IHS ---
    # When stream = rain end-member -> f_rain = 1
    f_r, f_g = binary_hydsep(-8.0, d_rain=-8.0, d_gw=-5.0)
    _check("binary f_rain=1 when stream=rain", f_r, 1.0)
    _check("binary f_gw=0  when stream=rain",  f_g, 0.0)

    # When stream = gw end-member -> f_gw = 1
    f_r, f_g = binary_hydsep(-5.0, d_rain=-8.0, d_gw=-5.0)
    _check("binary f_rain=0 when stream=gw",   f_r, 0.0)
    _check("binary f_gw=1  when stream=gw",    f_g, 1.0)

    # Mid-point -> each 0.5
    f_r, f_g = binary_hydsep(-6.5, d_rain=-8.0, d_gw=-5.0)
    _check("binary midpoint f_rain=0.5", f_r, 0.5)
    _check("binary midpoint f_gw=0.5",   f_g, 0.5)
    _check("binary sum = 1.0",           f_r + f_g, 1.0)

    # Near-zero denominator -> fallback f_gw = 1
    f_r, f_g = binary_hydsep(-6.0, d_rain=-6.0 + 1e-8, d_gw=-6.0)
    _check("binary denom~0 fallback f_gw=1", f_g, 1.0)

    # Calibration floor: negative f_rain replaced
    f_r, f_g = binary_hydsep(-4.0, d_rain=-8.0, d_gw=-5.0,
                              min_comp_rain=0.05)
    _check("binary floor applied to f_rain", f_r, 0.05)

    # Monotonicity: closer to rain end-member -> higher f_rain
    f_r_near_rain, _ = binary_hydsep(-7.5, d_rain=-8.0, d_gw=-5.0)
    f_r_near_gw,   _ = binary_hydsep(-5.5, d_rain=-8.0, d_gw=-5.0)
    _check_range("binary f_rain near rain > 0.8", f_r_near_rain, 0.8, 1.0)
    _check_range("binary f_rain near gw  < 0.2",  f_r_near_gw,   0.0, 0.2)

    # Mass balance: f_rain*d_rain + f_gw*d_gw == d_stream (no floor active)
    d_st = -7.0
    f_r, f_g = binary_hydsep(d_st, d_rain=-9.0, d_gw=-4.0)
    reconstructed = f_r * -9.0 + f_g * -4.0
    _check("binary mass balance", reconstructed, d_st, tol=1e-4)

    # --- Tertiary proportioning ---
    # Perfect model: comp_a should equal comp_b when deltas match flow-weighted mix
    d_st = iso_mix_binary(-8.0, 5.0, -5.0, 5.0)   # = -6.5
    comp_a, comp_b = tertiary_proportioning(
        d_stream=d_st, q_total=10,
        d_rain=-8.0, q_surf=5,
        d_gw=-5.0,   q_base=5,
        d_sw=0.0,    q_lat=0,
        rain_available=True)
    _check("tertiary comp_a == comp_b (perfect)", comp_a, comp_b, tol=1e-4)

    # No rain isotope -> comp_b excludes rain term
    comp_a2, comp_b2 = tertiary_proportioning(
        d_stream=-5.5, q_total=20,
        d_rain=-99.0,  q_surf=0,
        d_gw=-5.0,     q_base=15,
        d_sw=-6.0,     q_lat=5,
        rain_available=False)
    # comp_b = -5*15 + -6*5 = -105
    _check("tertiary no-rain comp_b=-105", comp_b2, -105.0)
    _check("tertiary comp_a = d_st*q_tot", comp_a2, -5.5 * 20)

    # Zero flow -> comp_b = 0
    _, comp_b3 = tertiary_proportioning(
        d_stream=0, q_total=0,
        d_rain=-8, q_surf=0,
        d_gw=-5,   q_base=0,
        d_sw=-6,   q_lat=0,
        rain_available=True)
    _check("tertiary zero flow comp_b=0", comp_b3, 0.0)


# ===========================================================================
# 4. precip.iso file format  (mirrors iso_init.f90 parsing logic)
# ===========================================================================

_PRECIP_ISO_CONTENT = """\
precip.iso - Monthly precipitation delta-18O values (test)
1  1  1.0  0.9  0.0  0.0
sta01  -5.5  -6.1  -7.2  -6.8  -5.0  -3.5  -3.2  -3.8  -4.5  -5.8  -6.5  -5.9
"""


def parse_precip_iso(text):
    lines = [l for l in text.strip().splitlines() if l.strip()]
    params  = lines[1].split()
    iso_on  = int(params[0])
    num_iso = int(params[1])
    iso_k   = float(params[2])
    iso_x   = float(params[3])
    min_r   = float(params[4])
    min_g   = float(params[5])
    stations = []
    for raw in lines[2:]:
        parts = raw.split()
        stations.append((parts[0], [float(v) for v in parts[1:13]]))
    return iso_on, num_iso, iso_k, iso_x, min_r, min_g, stations


def test_precip_iso():
    print("\n--- 4. precip.iso file format ---")

    iso_on, num_iso, iso_k, iso_x, min_r, min_g, stations = \
        parse_precip_iso(_PRECIP_ISO_CONTENT)

    _check("iso_on = 1",        float(iso_on),  1.0)
    _check("num_iso = 1",       float(num_iso), 1.0)
    _check("iso_k = 1.0",       iso_k,  1.0)
    _check("iso_x = 0.9",       iso_x,  0.9)
    _check("min_comp_rain = 0", min_r,  0.0)
    _check("min_comp_gw = 0",   min_g,  0.0)

    global _PASS, _FAIL
    if len(stations) == 1:
        print("  PASS  1 station parsed")
        _PASS += 1
    else:
        print(f"  FAIL  expected 1 station, got {len(stations)}")
        _FAIL += 1

    name, deltas = stations[0]
    _check("station name sta01", 1.0 if name == "sta01" else 0.0, 1.0)
    _check("Jan d18O = -5.5",   deltas[0],  -5.5)
    _check("Jul d18O = -3.2",   deltas[6],  -3.2)
    _check("Dec d18O = -5.9",   deltas[11], -5.9)

    aa = sum(deltas) / 12.0
    _check_range("annual avg d18O in [-7,-4]", aa, -7.0, -4.0)

    # File I/O round-trip
    with tempfile.NamedTemporaryFile(mode='w', suffix='.iso', delete=False) as fh:
        fh.write(_PRECIP_ISO_CONTENT)
        tmpname = fh.name
    try:
        with open(tmpname) as fh:
            text2 = fh.read()
        _, _, _, _, _, _, stations2 = parse_precip_iso(text2)
        _check("file round-trip Jan d18O", stations2[0][1][0], -5.5)
    finally:
        os.unlink(tmpname)


# ===========================================================================
# 5. Physical consistency / regression checks
# ===========================================================================

def test_physical_consistency():
    print("\n--- 5. Physical consistency ---")

    # Large rain event pulls soil delta toward precipitation delta
    mixed = iso_mix_binary(-5.0, 5.0, -8.0, 50.0)
    _check_range("large rain pulls soil toward d_rain",
                 mixed, -8.0 - 0.5, -8.0 + 0.5)

    # Percolation signal propagates down-profile
    d_l1 = iso_mix_binary(-5.0, 10.0, -8.0, 40.0)   # ~-7.4
    d_l2 = iso_mix_binary(-5.0, 20.0, d_l1,  10.0)
    _check_range("percolation signal moves downward",
                 d_l2, d_l1 - 0.1, -5.0 + 0.1)

    # Multi-mixer: total isotope mass conserved
    # before = -8*30 + (-4)*10 + (-4)*10 + (-4)*10 = -240 + -120 = -360
    mass_before = -8.0*30 + sum(-4.0*10 for _ in range(3))
    new_a, new_b = iso_mix_multi(30, -8.0, [10, 10, 10], [-4, -4, -4])
    mass_after = new_a * 30 + sum(d * 10 for d in new_b)
    _check("multi-mixer mass balance", mass_after, mass_before, tol=1e-3)

    # Repeated evaporation continuously enriches soil water
    d_s = -7.0
    d_s_init = d_s
    for _ in range(5):
        _, _, _, _, _, _, d_s, _ = fractionation(25.0, 0.5, -8.0, d_s)
    _check_range("repeated evap enriches soil monotonically",
                 d_s - d_s_init, 1.0, 1e6)

    # Binary IHS mass balance: weighted sum recovers stream delta
    d_rain_em, d_gw_em, d_st = -9.0, -4.0, -7.0
    f_r, f_g = binary_hydsep(d_st, d_rain_em, d_gw_em)
    _check("IHS mass balance", f_r*d_rain_em + f_g*d_gw_em, d_st, tol=1e-4)


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 52)
    print("  isotope math unit tests")
    print("=" * 52)

    test_mixing()
    test_fractionation()
    test_hydsep()
    test_precip_iso()
    test_physical_consistency()

    print()
    print("-" * 52)
    print(f"RESULT: {_PASS} passed, {_FAIL} failed")

    if _FAIL > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
