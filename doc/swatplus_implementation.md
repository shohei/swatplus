# SWAT+ Fortran Implementation Guide

## Table of Contents

1. [Overview](#1-overview)
2. [Repository Structure](#2-repository-structure)
3. [Build System](#3-build-system)
4. [Module Architecture](#4-module-architecture)
5. [Initialization Sequence](#5-initialization-sequence)
6. [Time Control and Daily Simulation](#6-time-control-and-daily-simulation)
7. [Routing Network and Command Loop](#7-routing-network-and-command-loop)
8. [Hydrologic Response Units (HRU)](#8-hydrologic-response-units-hru)
9. [Spatial Objects and Data Structures](#9-spatial-objects-and-data-structures)
10. [Key Process Modules](#10-key-process-modules)
11. [Constituent Transport Framework](#11-constituent-transport-framework)
12. [Configuration Codes](#12-configuration-codes)
13. [Input/Output Files](#13-inputoutput-files)
14. [Extending SWAT+](#14-extending-swat)

---

## 1. Overview

SWAT+ (Soil & Water Assessment Tool Plus) is a modular, daily-timestep watershed hydrology and water quality model written in Fortran. It is the successor to SWAT2012 and adopts a flexible spatial object routing architecture that allows arbitrary watershed topologies.

The model simulates:
- Surface runoff, lateral flow, percolation, and tile drainage
- Evapotranspiration (Priestley-Taylor, Penman-Monteith, or Hargreaves)
- Snow accumulation and melt
- Sediment detachment and transport
- Plant growth and nutrient cycling (N, P, C)
- Channel routing (variable storage or Muskingum)
- Reservoir and wetland operations
- Groundwater (shallow aquifer with optional 2D MODFLOW-style)
- Constituent transport: salts, boron, selenium, pathogens, pesticides
- Stable water isotopes (δ¹⁸O, δD) — added in this branch

---

## 2. Repository Structure

```
swatplus/
├── src/               # All Fortran source files (~200 .f90 files)
├── doc/               # Documentation
├── j2k_iso/           # J2000 Java isotope sources (reference)
├── CMakeLists.txt     # Top-level CMake build definition
├── CMakePresets.json  # Build presets (gfortran_release_macbook, etc.)
└── build/             # CMake build output (generated)
    └── release/       # Release build directory
        └── swatplus   # Compiled executable
```

All simulation logic lives in `src/`. The files follow a loose naming convention:

| Prefix | Domain |
|--------|--------|
| `proc_*` | Top-level initialization subroutines |
| `hru_*` | HRU-level operations |
| `cs_*` | Constituent (solute) transport |
| `iso_*` | Isotope tracking |
| `cha_*` / `rte_*` | Channel routing |
| `res_*` | Reservoir operations |
| `aqu_*` / `gwflow_*` | Aquifer / groundwater |
| `wet_*` | Wetland |
| `plant_*` / `grow_*` | Plant growth |
| `nu_*` | Nutrient cycling |
| `*_module.f90` | Data structure definitions (no executable code) |
| `*_init.f90` | One-time initialization routines |
| `*_read.f90` | Input file readers |

---

## 3. Build System

SWAT+ uses CMake with Fortran support. The recommended workflow:

```bash
# Configure (re-run when source files are added)
cmake --preset gfortran_release_macbook

# Compile (parallel)
cmake --build build/release -j$(sysctl -n hw.logicalcpu)
```

The `gfortran_release_macbook` preset (defined in `CMakePresets.json`) selects gfortran, sets optimization flags, and targets the `build/release/` output directory.

**Important**: CMake's Fortran module dependency scanner automatically determines compilation order from `USE` statements. When you add a new module file that other files `USE`, re-run `cmake --preset ...` before building so CMake rescans dependencies. Failing to do so causes `Fatal Error: Cannot open module file '*.mod'`.

---

## 4. Module Architecture

### 4.1 Dependency Hierarchy

SWAT+ uses Fortran `module` / `use` for all shared state. The dependency chain from most-fundamental to most-derived:

```
time_module              ← simulation calendar (time%, nbyr, etc.)
maximum_data_module      ← array size limits (db_mx%)
    ↓
basin_module             ← watershed parameters (bsn%, bsn_cc%)
soil_module              ← soil physical/chemical arrays (soil())
hru_module               ← HRU attributes (hru(), wst())
climate_module           ← weather types (weather_daily)
    ↓
hydrograph_module        ← routing network (ob(), sp_ob, hyd_output)
    ↓
constituent_mass_module  ← solute state (constituent_mass type)
cs_module                ← solute balance types (cs_balance)
    ↓
isotope_module           ← isotope state (iso_soil, iso_aqu, etc.)
    ↓
process subroutines      ← hru_control, cs_lch, iso_lch, etc.
```

Each `*_module.f90` file declares a Fortran module containing:
- Derived type definitions (structs)
- Global allocatable arrays
- Module-level scalars

Subroutines access global state by `use module_name` rather than passing large argument lists.

### 4.2 Key Module Files

| File | Contents |
|------|----------|
| `time_module.f90` | `time_t` type, `time%`, `nbyr`, `iyr`, `jday` |
| `maximum_data_module.f90` | `db_mx%` (array bounds: `wst`, `hru`, `sol`, etc.) |
| `basin_module.f90` | `basin_control_codes`, `bsn_cc%`, `bsn%` |
| `hydrograph_module.f90` | `spatial_objects` (`sp_ob`), `hyd_output`, `ob()` array |
| `hru_module.f90` | `hydrologic_response_unit`, `hru()` array |
| `soil_module.f90` | `soil_phys`, `soil_chem`, `soil()` array |
| `climate_module.f90` | `weather_daily`, `wst()`, `w` (current weather) |
| `constituent_mass_module.f90` | `constituent_mass`, `cs_soil()`, `cs_db` |

---

## 5. Initialization Sequence

`src/main.f90` is the program entry point. It calls ~40 initialization subroutines before entering the time loop:

```fortran
program main
  ! --- Basin / time / database ---
  call proc_bsn          ! read basin parameters (file.cio, bsn.cio)
  call proc_date_time    ! parse simulation dates
  call proc_db           ! read all database files (plants, soils, land use)
  call proc_read         ! read spatial input tables (.hru, .sub, etc.)

  ! --- Connectivity ---
  call hyd_connect       ! build the routing network (ob() linked list)
  call recalldb_read     ! read recall hydrographs
  call exco_db_read      ! read export coefficient data
  call dr_db_read        ! read drainage routing data

  ! --- Build cmd_next linked list ---
  allocate(cmd_next(sp_ob%objs))
  icmd = sp_ob1%objs
  do while (icmd /= 0)
    cmd_next(iob) = icmd
    icmd = ob(icmd)%cmd_next
  end do

  ! --- Climate / output setup ---
  call cli_lapse         ! compute lapse-rate climate adjustments
  call object_read_output  ! configure output files

  ! --- Water management ---
  call om_water_init
  call pest_cha_res_read
  call salt_cha_read
  call cs_cha_read

  ! --- HRU initialization (most critical block) ---
  call lsu_read_elements  ! landscape unit → HRU mapping
  call proc_hru           ! allocate and initialize all HRU state
  call proc_cha           ! initialize channels
  call proc_aqu           ! initialize aquifers

  ! --- Decision tables / management ---
  call dtbl_lum_read
  call hru_lte_read
  call proc_cond

  ! --- Reservoir / wetland ---
  call res_read_weir
  call dtbl_res_read
  call proc_res
  call wet_read_hyd, wet_read, wet_read_salt_cs
  if (db_mx%wet_dat > 0) call wet_all_initial
  call wet_fp_init

  ! --- Soil nutrient initialization ---
  do ihru = 1, sp_ob%hru
    call soil_nutcarb_init(isol)
  end do

  ! --- Calibration / output ---
  call proc_cal
  call proc_open
  call unit_hyd_ru_hru   ! unit hydrograph for subdaily runoff
  call dr_ru             ! drainage area ratios
  call hyd_connect_out   ! output connectivity

  ! --- Enter time loop ---
  call time_control
end program
```

### 5.1 proc_hru Sequence

`src/proc_hru.f90` handles all HRU-level initialization:

```fortran
subroutine proc_hru
  call hru_allo            ! allocate hru(), soil(), etc.
  call hru_read            ! read .hru input files
  call hrudb_init          ! look up database entries
  call hru_lum_init_all    ! land use management setup
  call topohyd_init        ! topographic/hydrologic parameters
  call hru_output_allo     ! allocate output arrays
  call carbon_read         ! initial C pools
  call soils_init          ! soil layer initialization
  call structure_init      ! tile drain, crack structures
  call plant_all_init      ! initial plant state
  call cn2_init_all        ! curve number adjustment
  call hydro_init          ! hydrologic parameters
  call pesticide_init      ! (if num_pests > 0)
  call pathogen_init       ! (if num_paths > 0)
  call salt_hru_init       ! (if num_salts > 0)
  call cs_hru_init         ! (if num_cs > 0)
  call iso_init            ! stable water isotope initialization
end subroutine
```

---

## 6. Time Control and Daily Simulation

`src/time_control.f90` implements the main simulation loop:

```fortran
subroutine time_control
  do curyr = 1, nbyr               ! annual loop
    do julian_day = day_start, day_end_yr  ! daily loop
      call climate_control         ! read/interpolate weather
      call command                 ! process all spatial objects for this day
      ! end-of-day bookkeeping
    end do
    ! end-of-year outputs
  end do
end subroutine
```

Key time variables (from `time_module`):

| Variable | Meaning |
|----------|---------|
| `time%yrc` | Current calendar year |
| `time%day` | Julian day of year (1–365/366) |
| `time%mo` | Month (1–12) |
| `time%day_mo` | Day of month |
| `iyr` | Internal year counter (1 to nbyr) |
| `nbyr` | Total simulation years |
| `time%step` | Timestep: 0=daily, <0=annual export coeff |

Leap years are handled by adjusting `day_end_yr` (365 or 366).

---

## 7. Routing Network and Command Loop

### 7.1 Architecture

SWAT+ replaces SWAT2012's fixed HRU→subbasin→reach hierarchy with a flexible linked-list routing network. Every watershed element (HRU, channel, reservoir, aquifer, etc.) is a **spatial object** (`ob()` array).

```fortran
! In hydrograph_module.f90
type :: spatial_object
  integer :: typ          ! object type (see below)
  integer :: num          ! index within that type's array
  integer :: cmd_next     ! index of next ob() to process
  integer :: rcv_tot      ! number of upstream inflow sources
  type(hyd_output) :: hin ! inflow hydrograph
  type(hyd_output) :: hout! outflow hydrograph
  type(hyd_output), allocatable :: hin_sur(:)  ! surface inflows
  type(hyd_output), allocatable :: hin_til(:)  ! tile/lat inflows
  ! ...
end type
type(spatial_object), allocatable :: ob(:)
```

`sp_ob` (type `spatial_objects`) stores the count of each object type:

```fortran
type :: spatial_objects
  integer :: objs    ! total spatial objects
  integer :: hru     ! number of HRUs
  integer :: aqu     ! aquifers
  integer :: chan    ! channels
  integer :: res     ! reservoirs
  integer :: outlet  ! outlet objects
  ! ... (ru, gwflow, recall, exco, dr, canal, pump, wet, chandeg, aqu2d)
end type
```

### 7.2 Command Loop (command.f90)

Each day, `command` traverses the linked list in watershed order:

```fortran
subroutine command
  icmd = sp_ob1%objs           ! start object (set by hyd_connect)
  do while (icmd /= 0)
    ! 1. Zero hydrograph accumulator
    ob(icmd)%hin = hz           ! hz = zero hydrograph

    ! 2. Accumulate all upstream inflows
    do ircv = 1, ob(icmd)%rcv_tot
      ob(icmd)%hin = ob(icmd)%hin + upstream_contribution
    end do

    ! 3. Dispatch to object-type process routine
    select case (ob(icmd)%typ)
      case (hru_typ)    ; call hru_control
      case (aqu_typ)    ; call aqu_control
      case (chan_typ)   ; call cha_control
      case (res_typ)    ; call res_control
      case (recall_typ) ; call recall_control
      case (outlet_typ) ; call outlet_control
      ! ... (ru, gwflow, exco, dr, canal, wet, chandeg, aqu2d)
    end select

    ! 4. Advance to next object in linked list
    icmd = ob(icmd)%cmd_next
  end do
end subroutine
```

### 7.3 Hydrograph Structure

All water fluxes use the `hyd_output` type:

```fortran
type :: hyd_output
  real :: flo        ! total flow (mm or m³/s)
  real :: sed        ! sediment
  real :: orgn, sedp ! organic N, sediment P
  real :: no3, solp  ! nitrate, soluble P
  real :: chla       ! chlorophyll-a
  real :: nh3, no2   ! ammonia, nitrite
  real :: cbod, dox  ! carbonaceous BOD, dissolved O₂
  real :: san, sil, cla, sag, lag, grv  ! particle sizes
end type
```

Hydrograph component indices (`hd(1:5)` for disaggregated tracking):

| Index | Component |
|-------|-----------|
| `hd(1)` | Total flow |
| `hd(2)` | Recharge / groundwater |
| `hd(3)` | Surface runoff |
| `hd(4)` | Lateral flow |
| `hd(5)` | Tile drainage |

---

## 8. Hydrologic Response Units (HRU)

The HRU is the fundamental land unit. All soil/plant/runoff processes execute at the HRU level once per day, called from `command` via `hru_control`.

### 8.1 HRU Data Structure

Defined in `hru_module.f90`:

```fortran
type :: hydrologic_response_unit
  character(len=16) :: name
  type(hru_databases) :: dbs     ! pointers to DB entries (soil, plant, wst...)
  type(topography)   :: topo     ! slope, area, latitude, elevation
  type(hydrology_parms) :: hyd   ! CN, esco, ov_n, perco, latq_co...
  type(hru_land_use) :: lumv     ! USLE factors, cover management
  integer :: wst                 ! weather station index
  integer :: tiledrain           ! tile drain structure index (0=none)
  ! ... (plant, management, septic, etc.)
end type
type(hydrologic_response_unit), allocatable :: hru(:)  ! hru(1:sp_ob%hru)
```

Soil state lives in `soil_module.f90`:

```fortran
type :: soil_phys_layer
  real :: st    ! soil water content (mm)
  real :: prk   ! percolation out of layer (mm/day)
  real :: flat  ! lateral flow (mm/day)
  real :: bd, por, awc, k, ...
end type
type :: soil_profile
  integer :: nly       ! number of layers
  type(soil_phys_layer), allocatable :: phys(:)  ! phys(1:nly)
  type(soil_chem_layer), allocatable :: ly(:)    ! ly(1:nly)
  ! ...
end type
type(soil_profile), allocatable :: soil(:)  ! soil(1:sp_ob%hru)
```

### 8.2 hru_control Daily Sequence

`src/hru_control.f90` orchestrates the daily HRU simulation. The call order within one day:

```
snow melt / accumulation          (snow_mlt)
canopy interception               (canopy)
irrigation check                  (mgt_ops_hru)
soil crack redistribution         (soil_crack)
surface storage / routing         (surface)
  ├─ CN runoff                    (surq_cnno / surq_cngr)
  └─ green-ampt infiltration      (grnampt_cn)
potential ET                      (et_pot)
actual ET                         (etact)
  ├─ plant transpiration
  ├─ soil evaporation
  └─ canopy evaporation
lateral flow                      (surq_rchg)
percolation                       (percmain)
tile drainage                     (tiledrain_perc)
groundwater contribution          (gwflow)
plant growth                      (plantmod → grow_parms)
management operations             (mgt_ops_hru)
sediment yield                    (ero_musle / ero_usle)
organic matter decomposition      (decomp → carbon_new)
nutrient transformations          (nutcarb_lch, no3_nloss...)
pesticide transport               (pest_lch)
pathogen transport                (path_lch)
salt transport                    (salt_lch)
constituent transport             (cs_rain, cs_lch)   ← cs_module
isotope tracking                  (iso_rain, iso_lch, iso_frac, iso_hydsep)
HRU output aggregation
```

---

## 9. Spatial Objects and Data Structures

### 9.1 Object Types

| Type constant | Spatial object | Process routine |
|---------------|---------------|-----------------|
| `hru_typ` | HRU | `hru_control` |
| `hru_lte_typ` | Long-term export HRU | `hru_lte_control` |
| `ru_typ` | Routing unit (HRU aggregate) | `ru_control` |
| `gwflow_typ` | 2D MODFLOW aquifer cell | `gwflow_simulate` |
| `aqu_typ` | Lumped shallow aquifer | `aqu_control` |
| `chan_typ` | Channel reach | `cha_control` |
| `chandeg_typ` | Channel degradation | `chandeg_simulate` |
| `res_typ` | Reservoir | `res_control` |
| `recall_typ` | Point source recall | `recall_control` |
| `exco_typ` | Export coefficient object | `exco_control` |
| `dr_typ` | Drainage area object | `dr_control` |
| `canal_typ` | Irrigation canal | `canal_control` |
| `pump_typ` | Pump | `pump_control` |
| `wet_typ` | Wetland | `wet_control` |
| `aqu2d_typ` | 2D aquifer (alternative) | — |
| `outlet_typ` | Watershed outlet | `outlet_control` |

### 9.2 Weather Data

Current-day weather for the HRU being processed is in the global variable `w` (type `weather_daily`, from `climate_module.f90`):

```fortran
type :: weather_daily
  real :: precip  ! precipitation (mm)
  real :: tmax    ! max temperature (°C)
  real :: tmin    ! min temperature (°C)
  real :: tave    ! average temperature (°C)
  real :: rhum    ! relative humidity (fraction, 0–1)
  real :: wndspd  ! wind speed (m/s)
  real :: solrad  ! solar radiation (MJ/m²)
  real :: snotmp  ! snowfall temperature threshold
end type
weather_daily :: w  ! current HRU's weather
```

Weather stations are indexed by `hru(ihru)%wst`. Multiple stations are defined in `wst(1:db_mx%wst)`.

---

## 10. Key Process Modules

### 10.1 Runoff Generation

- **Curve Number method** (`surq_cnno.f90`, `surq_cngr.f90`): Standard SCS-CN with daily antecedent moisture correction. CN2 is adjusted for slope, plant cover, and soil moisture.
- **Green-Ampt infiltration** (`grnampt_cn.f90`): Optional sub-daily alternative.

### 10.2 Evapotranspiration

Controlled by `bsn_cc%pet`:

| Code | Method | File |
|------|--------|-------|
| 0 | Priestley-Taylor | `et_pot.f90` |
| 1 | Penman-Monteith | `et_pot.f90` |
| 2 | Hargreaves | `et_pot.f90` |
| 3 | User-defined | — |

### 10.3 Soil Water Routing

- **Percolation** (`percmain.f90`): Travel-time based, layer by layer. `soil(j)%ly(jj)%prk` = mm percolated.
- **Lateral flow** (`surq_rchg.f90`): Kinematic storage model. `soil(j)%ly(jj)%flat` = mm lateral flow.
- **Tile drainage** (`tiledrain_perc.f90`): Hooghoudt's equation.

### 10.4 Channel Routing

Controlled by `bsn_cc%rte`:

| Code | Method |
|------|--------|
| 0 | Variable storage (SWAT default) |
| 1 | Muskingum |

Files: `cha_control.f90`, `rte_route.f90`, `musk_route.f90`.

### 10.5 Plant Growth

`plantmod.f90` → `grow_parms.f90`:
- Radiation use efficiency (RUE) based biomass accumulation
- Water, temperature, N, P stress factors
- LAI development and senescence
- Harvest: `mgt_ops_hru.f90` reads management schedules

### 10.6 Nutrient Cycling

- **Carbon**: CENTURY (cswat=2), C-FARM (cswat=1), or static (cswat=0)
- **Nitrogen**: nitrification, denitrification, mineralization (`nutcarb_lch.f90`, `no3_nloss.f90`)
- **Phosphorus**: sorption, mineralization (`phosmin.f90`)

### 10.7 Aquifer (Shallow Groundwater)

`aqu_control.f90`:
- Single lumped storage per `aqu_typ` object
- `aqu_d%rchrg` = recharge from percolation out of soil bottom
- `aqu_d%flo` = baseflow to channel (Boussinesq recession)
- `aqu_d%stor` = current storage (mm)
- Percolation into deep aquifer: `sepbtm(j)` in soil layer context

`rls_routeaqu.f90` routes aquifer discharge back to its receiving channel.

---

## 11. Constituent Transport Framework

SWAT+ has a general constituent transport framework (`cs_module`, `constituent_mass_module`) for tracking conservative and reactive solutes (boron, selenium, custom constituents). This serves as the template for all new solute tracking.

### 11.1 Data Structures

```fortran
! cs_module.f90
type :: cs_balance
  real, allocatable :: soil(:)  ! (nly) — in each soil layer
  real :: surq    ! in surface runoff
  real :: latq    ! in lateral flow
  real :: tile    ! in tile drainage
  real :: perc    ! percolating to aquifer
  real :: rain    ! from precipitation
  real :: dryd    ! dry deposition
  real, allocatable :: conc(:) ! (nly) — concentration in soil water
end type

! constituent_mass_module.f90
type :: constituent_mass
  type(cs_balance), allocatable :: cs(:)   ! (num_cs) per constituent
  type(cs_balance), allocatable :: csc(:)  ! cumulative
end type
type(constituent_mass), allocatable :: cs_soil(:)  ! cs_soil(0:sp_ob%hru)
```

### 11.2 Extension Pattern

To add a new constituent (e.g., isotopes), the pattern is:

1. Create `*_module.f90` with new state arrays
2. Create `*_init.f90` to allocate and initialize
3. Create `*_rain.f90` for atmospheric input
4. Create `*_lch.f90` for soil layer transport
5. Add `use new_module` in `proc_hru.f90` and `hru_control.f90`
6. Call `init` from `proc_hru`, call process routines from `hru_control`
7. Re-run CMake to scan new `USE` dependencies

The isotope module (`isotope_module.f90` etc.) follows this pattern exactly.

---

## 12. Configuration Codes

`bsn_cc` (type `basin_control_codes`, read from `basins.bsn`) controls major process switches:

| Field | Values | Effect |
|-------|--------|--------|
| `pet` | 0/1/2 | ET method: Priestley-Taylor / Penman-Monteith / Hargreaves |
| `rte` | 0/1 | Routing: variable storage / Muskingum |
| `sed_det` | 0/1 | Sediment detachment: MUSLE / USLE |
| `cswat` | 0/1/2 | Carbon model: static / C-FARM / Century |
| `tdrn` | 0/1 | Tile drainage on/off |
| `wtdn` | 0/1 | Water table depth model |
| `crk` | 0/1 | Soil crack flow |
| `gwflow` | 0/1 | 2D MODFLOW groundwater |
| `gampt` | — | GAMPT plant parameter set |
| `qual2e` | 0/1 | QUAL2E in-stream water quality |
| `swift_out` | 0/1 | Write SWIFT input files |

---

## 13. Input/Output Files

### 13.1 Primary Input Files

| File | Contents |
|------|----------|
| `file.cio` | Master file list (all other inputs listed here) |
| `time.sim` | Simulation period (start/end dates) |
| `print.prt` | Output variable selection |
| `basins.bsn` | Basin parameters and control codes |
| `*.hru` | HRU physical attributes |
| `*.sol` | Soil profiles |
| `*.pcp`, `*.tmp`, `*.hmd`, `*.wnd`, `*.slr` | Climate inputs by weather station |
| `*.mgt` | Management schedules |
| `hru-data.hru` | HRU-to-database linkage |
| `topography.hyd` | Slope, area, connectivity |

### 13.2 Key Output Files

| File | Contents |
|------|----------|
| `simulation.out` | Run log and status |
| `success.fin` | Created on successful completion |
| `hru_wb.txt` | HRU water balance |
| `hru_nb.txt` | HRU nutrient balance |
| `basin_wb.txt` | Basin-wide water balance |
| `cha_wb.txt` | Channel water balance |
| `aqu_wb.txt` | Aquifer water balance |
| `checker.out` | Initialization parameter check |
| `erosion.out` | Sediment by HRU |

### 13.3 Isotope Input (this branch)

| File | Contents |
|------|----------|
| `precip.iso` | Monthly δ¹⁸O per weather station; `iso_on` flag |

Format:
```
! title comment
iso_on  num_iso  k  x  min_comp_rain  min_comp_gw
! station 1 monthly delta-18O (Jan-Dec)
-5.2  -4.8  -4.1  -3.5  -3.0  -2.8  -2.9  -3.4  -4.0  -4.6  -5.0  -5.3
! station 2 ...
```

---

## 14. Extending SWAT+

### 14.1 Adding a New Process

1. **Module**: Create `src/myprocess_module.f90` with state arrays and parameters.
2. **Init**: Create `src/myprocess_init.f90`; allocate arrays; read input if needed.
3. **Process**: Create `src/myprocess_calc.f90` with the daily calculation.
4. **Hook into proc_hru**: In `src/proc_hru.f90`, add `use myprocess_module` and `external :: myprocess_init`, then `call myprocess_init`.
5. **Hook into hru_control**: In `src/hru_control.f90`, add `use myprocess_module`, `external :: myprocess_calc`, and `call myprocess_calc` at the appropriate point in the daily sequence.
6. **Reconfigure**: `cmake --preset gfortran_release_macbook` (forces CMake to rescan module dependencies).
7. **Build**: `cmake --build build/release -j$(sysctl -n hw.logicalcpu)`.

### 14.2 Adding a New Spatial Object Type

1. Add a new type constant in `hydrograph_module.f90`.
2. Add the count field in `spatial_objects`.
3. Write a `*_read.f90` to populate `ob()` entries.
4. Add a `case` branch in `command.f90`'s `select case (ob(icmd)%typ)`.
5. Write the `*_control.f90` simulation routine.

### 14.3 Accessing Upstream/Downstream State

Within a process called from `command.f90`:
- `ob(icmd)%hin` — aggregated inflow from all upstream objects (already accumulated by the command loop)
- `ob(icmd)%hout` — set this to the object's outflow (read by downstream objects)
- `ob(icmd)%num` — index into the type-specific array (e.g., `hru(ob(icmd)%num)`)

The command loop guarantees that all upstream objects have completed before the current object is processed, because `hyd_connect` orders `cmd_next` in topological (upstream-first) order.

---

## References

- Arnold, J.G. et al. (2012). SWAT: Model use, calibration, and validation. *Transactions of the ASABE*, 55(4): 1491–1508.
- Bieger, K. et al. (2017). Introduction to SWAT+, a completely restructured version of the Soil and Water Assessment Tool. *JAWRA*, 53(1): 115–130.
- Neitsch, S.L. et al. (2011). *Soil and Water Assessment Tool Theoretical Documentation Version 2009*. Texas Water Resources Institute.
- SWAT+ source repository: https://github.com/swat-model/swatplus
