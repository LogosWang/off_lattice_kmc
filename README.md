# Off-Lattice KMC Simulator

An off-lattice Kinetic Monte Carlo (KMC) simulator for studying irradiation-assisted oxidation in austenitic stainless steel. The model captures the coupled kinetics of solute radiation induced segregation (RIS) by implementing effective chemical potential gradients, and oxygen absorption, diffusion and oxidation at a grain boundary (GB).

---

## Physical Model

The simulation box is a 3D cube (default 1 mm edge length) containing a tilted grain boundary plane (`A·x + B·y + C·z + D = 0`). Four particle species are tracked:

| Species | Type Code | Role |
|---------|-----------|------|
| Si | 0 | Solute, irradiation-driven depletion at GB |
| Cr | 1 | Solute, irradiation-driven enrichment at GB |
| Ni | 3 | Solute, irradiation-driven depletion at GB |
| O  | 2 | Oxygen, adsorbs on GB surface and diffuses in-plane |

### Events

| Event | Description |
|-------|-------------|
| `SILI_DIFF` | Si random walk, modulated by local stress and irradiation chemical potential |
| `CROM_DIFF` | Cr random walk, modulated by local stress and irradiation chemical potential |
| `NIK_DIFF`  | Ni random walk, modulated by local stress and irradiation chemical potential |
| `OXYGEN_ABS` | New O atom adsorbs at the GB–surface intersection line |
| `OXYGEN_DIFF` | O atom diffuses within the GB plane |

### Rate Expression

Metal atom hop rates follow a stress- and irradiation-modified Arrhenius form:

```
rate = v0 * exp(-(ΔE + tr(σ)·Ω - μ_eff) / kT) × propensity_scalar
```

- `tr(σ)·Ω`: local hydrostatic stress (interpolated from `stress.csv`) times activation volume
- `μ_eff = kT · dose · 0.15 · ln(...) · exp(-d_GB / 5 μm)`: irradiation-induced effective chemical potential, decaying exponentially away from the GB
- `propensity_scalar`: spatially adaptive resolution factor (see below)

### Spatially Adaptive Resolution

The simulator implements adaptive spatial resolution near the GB: finer resolution is applied close to the GB where oxidation chemistry is active, and coarser resolution is used far away. This is achieved by jointly scaling the jump distance and hop rate while keeping the diffusivity D invariant.

A position-dependent scalar `s(d_GB)` is defined via a sigmoid function, ranging from **0.1** at the GB to **1.0** far away:

```
s(d) = 10^( (1 - 2) / (1 + exp((d - 30 μm) / length_scale)) )
```

The actual jump distance and propensity used at each step are:

```
jump_actual      = jump_base × s          →  smaller near GB (finer resolution)
propensity_actual = propensity_base / s²  →  larger near GB (more frequent steps)
```

Because diffusivity D ∝ rate × jump², the product is conserved:

```
D ∝ (rate / s²) × (jump × s)² = rate × jump²  ✓
```

This ensures the physical diffusivity is unchanged across the simulation domain, while the GB region is sampled at up to 10× finer spatial and temporal resolution.

Oxidation occurs stochastically when a solute atom and an O atom come within `OXIDATION_THRESH` (250 nm), with a probability based on the oxidation Gibbs free energy of each species (SiO₂ > Cr₂O₃ > NiO thermodynamic driving force order), and Cr oxide count to capture the protective effect of Cr oxide.

### Key Physical Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Temperature | 673 K (400 °C) | |
| Migration barrier | 1.0 eV | |
| Pre-factor v0 | 1×10¹³ s⁻¹ | |
| Activation volume Ω | 1×10⁻²⁸ m³ | |
| Oxidation threshold | 250 nm | Coarse-grain reaction radius |
| Max O atoms | 6000 | Adsorption cap |

---

## Build

Requires a C++17-capable compiler (tested with `clang++`).

```bash
make
```

This produces the `kmc_sim` executable.

```bash
make clean   # remove compiled objects and executable
```

---

## Usage

```bash
./kmc_sim -i <input_file.i>
```

If `-i` is omitted, the simulator runs with hardcoded default parameters.

### Input File Format (`.i`)

Plain text, one keyword per line. Lines starting with `#` are treated as comments.

```
dose            5.0       # irradiation dose (dpa), scales irradiation chemical potential
GB_normal_stress  500.0   # GB normal stress (MPa), scales O adsorption rate
output_dir      50_0500   # directory for all output files (created automatically)
```

Several example input files are provided:

| File | Dose (dpa) | GB Normal Stress (MPa) |
|------|-----------|------------------------|
| `00.i` | 0.0 | 0 |
| `05.i` | 0.5 | 0 |
| `10.i` | 1.0 | 0 |
| `30.i` | 3.0 | 0 |
| `50.i` | 5.0 | 0 |
| `50_0000.i` | 5.0 | 0 |
| `50_0500.i` | 5.0 | 500 |
| `50_1000.i` | 5.0 | 1000 |
| `50_1500.i` | 5.0 | 1500 |
| `50_2000.i` | 5.0 | 2000 |

---

## Required Input Data Files

These files must be present in the working directory at runtime:

### `sites3type_more.csv`
Initial particle positions. Format (with header):
```
id,x,y,z,type
0,0.000123,0.000456,0.000789,1
...
```

### `stress.csv`
Static stress field from external FEM/MD calculation. Format (with header, 7 columns):
```
x,y,z,trace,grad_x,grad_y,grad_z
...
```
- `trace`: hydrostatic stress trace tr(σ) in Pa
- `grad_x/y/z`: stress gradient vector (used in TRACKDIFF mode, currently unused)

---

## Output Files

All outputs are written to `output_dir/`.

| File | Description |
|------|-------------|
| `sites_<step>.xyz` | Particle positions in XYZ format at step N (VMD-compatible) |
| `sitesresult<dir>.csv` | Final positions of all non-oxidized particles |
| `oxidegrow<dir>.csv` | Time series of oxide counts: `time, Si_oxide, Cr_oxide, Ni_oxide, Total_oxide` |
| `propensity<time>.csv` | Snapshot of all event propensities at a given time |

Trajectory snapshots (`sites_*.xyz`) are written every 10⁶ steps. Final state is written at simulation end.

---

## Simulation Parameters (in `main.cpp`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `num_sites` | 1000 | Number of metal atoms |
| `box_size` | 1×10⁻³ m | Simulation box edge length |
| `max_sim_time` | 2.0 s | Wall-clock equivalent simulation time |
| `max_sim_steps` | 3×10¹⁰ | Maximum KMC steps |
| `jump_distance` | 1×10⁻⁶ m | Base hop distance |
| `CSVflag` | 1 | 1 = load sites from CSV, 0 = random initialization |
| `Diff_mod` | 0 | 0 = `RANDDIFF` (random walk), 1 = `TRACKDIFF` (deprecated) |

---

## Visualization

Use OVITO to load trajectory files:

```
vmd -e load_sites_trajectory.tcl
```

Jupyter notebooks for analysis:

| Notebook | Purpose |
|----------|---------|
| `oxidationvstime.ipynb` | Plot oxide count vs. time |
| `sitesdistribution.ipynb` | Spatial distribution of sites |
| `sitesgeneration.ipynb` | Generate initial site configurations |
| `stressfield.ipynb` | Extract stress field |
| `stress_grad.ipynb` | Compute stress gradients |
| `kmc_realtime_visualization_STRICT.ipynb` | Real-time solute concentration visualization |

---

## Code Structure

| File | Description |
|------|-------------|
| `main.cpp` | Entry point, parameter definitions |
| `KMC_Simulator.h/cpp` | Core KMC engine: event management, Gillespie algorithm, I/O |
| `Site.h/cpp` | Particle data structure |
| `Event.h/cpp` | Event data structure (type, propensity, site_id) |
| `GB.h/cpp` | Grain boundary plane geometry, O adsorption site generation, 2D PBC |
| `StressPoint.h/cpp` | Stress field data point structure |
