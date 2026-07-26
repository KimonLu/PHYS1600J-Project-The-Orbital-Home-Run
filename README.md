<p align="center">
  <img src="assets/background.png" alt="An astronaut playing baseball on the Moon" width="52%">
</p>

<h1 align="center">The Orbital Home Run</h1>
<p align="center"><strong>PHYS1600J Project</strong></p>

<p align="center">
  <strong>English</strong> | <a href="README.zh-CN.md">中文</a>
</p>

Could a baseball struck on the Moon enter a low lunar orbit, travel around the
Moon, and return to the batter from behind? This repository presents our
answer to PHYS1600J Problem B through analytic orbital mechanics, terrain-aware
high-fidelity numerical modelling, a formal project report, and a bilingual
interactive trajectory solver.

- **Website:** [Orbital Home Run Solver](https://kimonlu.github.io/PHYS1600J-Project-The-Orbital-Home-Run/)
- **Report:** [Physics1600J_Project.pdf](Project/docs/Physics1600J_Project.pdf)
- **Repository:** [KimonLu/PHYS1600J-Project-The-Orbital-Home-Run](https://github.com/KimonLu/PHYS1600J-Project-The-Orbital-Home-Run)

## Project Information

- **Course:** PHYS1600J — Honors Physics
- **Problem:** Problem B — The Orbital Home Run
- **Full report title:** *The Orbital Home Run — From a Measure-Zero Ideal
  Orbit to a Terrain-Aware Lunar Return*
- **Team:** Team 09
- **Team members:** Kemeng Lu, Jiahao Yin, Yuxin Wu
- **Supervisor:** Prof. Zijie Qu
- **Affiliation:** Global College, Shanghai Jiao Tong University
- **Main deliverables:** scientific report, reproducible Python model,
  versioned derived data and figures, and a bilingual browser solver
- **Reference models:** LOLA LDEM64 terrain and GRAIL GRGM1200B gravity
- **Implementation:** Python, TypeScript, Three.js, Vite, LaTeX
- **License:** [MIT](LICENSE)

## Project Background

On Earth, a home run is primarily a short atmospheric flight. On the Moon,
where there is no appreciable atmosphere, a sufficiently fast baseball must
instead be treated as an orbiting point mass. The familiar “first cosmic
velocity” is therefore only the beginning of the problem.

A self-hitting lunar home run must satisfy four separate conditions:

1. the trajectory remains gravitationally bound to the Moon;
2. it completes a revolution without intersecting the lunar surface or
   terrain;
3. it returns to the rotating launch site, rather than merely to the original
   point in an inertial frame; and
4. the conclusion remains meaningful after realistic gravity, terrain,
   numerical error, and uncertainty are considered.

The project develops a model hierarchy from an exact spherical two-body
analysis to a three-dimensional terrain-aware propagation model. The final
reference calculation uses the full LOLA LDEM64 raster, GRGM1200B gravity
truncated at a convergence-tested degree and order 600, uniform lunar
rotation, and fixed short-arc Earth and solar differential tides.

## Project Results

- **Ideal speed scale.** At the mean lunar radius, the circular speed is
  **1.680 km/s** and the escape speed is **2.376 km/s**.
- **Surface-launch condition.** On a smooth, nonrotating spherical Moon, a
  bound launch from the reference surface is non-penetrating only for an
  exactly horizontal launch with $1 \leq v/v_c < \sqrt{2}$. The feasible
  speed-angle set is therefore measure zero and structurally fragile.
- **Lunar rotation.** Rotation does not prohibit a return, but it changes the
  required encounter time. An equatorial prograde circular orbit revisits the
  moving field after the synodic time $2\pi/(n-\Omega)$.
- **Concrete terrain-aware return.** The degree-600 boundary-value solution
  launches from $5.4296875^\circ\text{N},\,201.3671875^\circ\text{E}$, with
  the ball centre **30.000 km above the reference sphere** and **19.243 km
  above the local LDEM64 terrain**.
- **Solved initial state.** The surface-relative launch speed is
  **1661.4559 m/s**, at **0.531611° elevation** and **88.383063° azimuth
  clockwise from north**. The scheduled return time is **111.439 min**.
- **Clearance and numerical residual.** The minimum LDEM64 clearance is
  **10.279 km**, and the fixed-time boundary-value residual is **4.690 mm**.
  The latter is an optimizer residual inside the deterministic model, not a
  claim of millimetre-scale physical accuracy.
- **Physical interpretation.** The return speed corresponds to approximately
  **200 kJ** of kinetic energy. The calculation rules out an unaided human
  hit and does not establish feasible launch or capture hardware.
- **Browser cross-check.** Against direct degree-600 Python trajectories, the
  browser gravity atlas differs by **6.06 m at 15 km** altitude and **4.04 m
  at 30 km** altitude over the validation circuit.

The project establishes a reproducible numerical initial state within a
clearly stated model. It does not claim a navigation-grade or physically
certified lunar experiment; such a claim would additionally require dated
SPICE geometry, local terrain surveying, gravity and launch covariance, solar
radiation pressure, and a robust constrained correction.

## Solver Website

The [bilingual static solver](https://kimonlu.github.io/PHYS1600J-Project-The-Orbital-Home-Run/)
runs entirely in the browser. Users specify a body-fixed lunar launch site,
height above terrain, surface-relative speed, elevation, azimuth, propagation
time, return tolerance, assumed position-error bound, and optional fixed
Earth/Sun tide terms.

The solver provides:

- English and Chinese interfaces that can be switched without clearing the
  current result;
- verified examples for return, uncertain return, terrain impact, escape, and
  no-return classifications;
- an interactive inertial-frame 3D trajectory with synchronized lunar
  rotation and 1×–1000× playback;
- an equirectangular body-fixed ground track;
- LDEM64 terrain clearance, first-impact coordinates, and the closest
  post-launch approach to the moving launch site; and
- model and domain warnings, including the numerical step-doubling diagnostic.

Terrain is streamed as lossless LDEM64 tiles. Gravity combines an analytic
central term with an altitude-layered acceleration-correction atlas generated
from GRGM1200B at degree and order 600. Rotating-frame dynamics include
Coriolis and centrifugal terms and are integrated with fourth-order
Runge–Kutta step doubling. The website is an exploratory research tool, not a
navigation, landing, or safety product.

## Project File Tree

### `Project/` — scientific model and report

```text
Project/
├── README.md
├── docs/
│   ├── main.tex                       # report source
│   ├── references.bib                 # bibliography
│   ├── setting.cls                    # report style
│   ├── project_statement.pdf          # original assignment
│   └── Physics1600J_Project.pdf       # final report
├── scripts/
│   ├── 01_ideal_models.py             # analytic two-body orbit family
│   ├── 02_sensitivity.py              # angle, height, and Monte Carlo tests
│   ├── 03_rotation_and_resonance.py   # lunar rotation and return resonance
│   ├── 04_realistic_perturbations.py  # low-order and third-body hierarchy
│   ├── 05_terrain_envelope.py         # conservative terrain envelope
│   ├── 06_validation.py               # central-model convergence
│   ├── 07_gravity_convergence.py      # GRGM truncation selection
│   ├── 08_prepare_web_data.py         # browser terrain/gravity products
│   ├── 09_validate_web_solver.py      # browser/direct trajectory checks
│   ├── 10_high_fidelity_case.py       # degree-600 boundary-value solution
│   ├── 11_case_sensitivity.py         # Jacobian convergence and sensitivity
│   ├── 12_validate_web_terrain.py     # lossless terrain-tile validation
│   ├── 13_terrain_visualizations.py   # terrain and trajectory figures
│   ├── 14_surface_feasibility.py      # bounded great-circle terrain study
│   ├── 15_high_fidelity_validation.py # integrator and gravity-order checks
│   ├── 16_height_continuation.py      # bounded release-height continuation
│   ├── general_solver.py              # public arbitrary-input solver
│   ├── lunar_gravity.py               # spherical-harmonic gravity model
│   ├── lunar_terrain.py               # LDEM64 access and interpolation
│   ├── orbital_home_run.py            # shared orbital dynamics
│   ├── plotting.py                    # shared plotting utilities
│   ├── generate_all.py                # reproduction profile runner
│   ├── validate_results.py            # packaged-result validation
│   ├── download_science_data.py       # verified NASA data downloader
│   └── requirements.txt
├── data/
│   ├── input/                         # constants and source-data manifest
│   ├── external/                      # downloaded NASA data; not versioned
│   └── output/                        # derived CSV/JSON/NPZ results
└── figures/                           # publication-ready figures
```

### `Web/` — bilingual browser solver

```text
Web/
├── index.html
├── package.json
├── package-lock.json
├── tsconfig.json
├── README.md
├── src/
│   ├── main.ts                        # interface and result presentation
│   ├── solver.worker.ts               # trajectory integration worker
│   ├── moon-view.ts                   # interactive 3D lunar view
│   ├── ground-track-view.ts           # body-fixed ground-track view
│   ├── data-store.ts                  # streamed scientific-data loader
│   ├── examples.ts                    # verified example trajectories
│   ├── i18n.ts                        # English/Chinese text
│   ├── types.ts                       # shared TypeScript types
│   └── style.css
└── public/
    ├── assets/                        # lunar visual assets
    └── data/
        ├── terrain/                   # lossless LDEM64 tiles and manifest
        └── gravity/                   # degree-600 atlas tiles and manifest
```

## Reproduction

### 1. Clone the repository and retrieve Git LFS data

Git, [Git LFS](https://git-lfs.com/), Python 3.10 or newer, and Node.js 22 are
recommended. TeX Live 2025 is only required to rebuild the report.

```bash
git clone https://github.com/KimonLu/PHYS1600J-Project-The-Orbital-Home-Run.git
cd PHYS1600J-Project-The-Orbital-Home-Run
git lfs install
git lfs pull
```

### 2. Reproduce the `Project/` calculations

Install the Python dependencies:

```bash
cd Project
python -m pip install -r scripts/requirements.txt
```

The lightweight analytic, sensitivity, rotation, perturbation, and central
validation workflow requires no large external download:

```bash
python scripts/generate_all.py --profile core
```

For the high-fidelity calculation, first download the authoritative LOLA
LDEM64 and GRAIL GRGM1200B products. The downloader verifies the expected
sizes and SHA-256 hashes.

```bash
python scripts/download_science_data.py
python scripts/generate_all.py --profile science
```

To reproduce every calculation, publication figure, and browser data product:

```bash
python scripts/generate_all.py --profile full
```

High-degree gravity and global-terrain stages are computationally intensive.
To check the committed machine-readable results without regenerating them:

```bash
python scripts/generate_all.py --profile validate
```

The nominal return is recorded in
[`data/output/high_fidelity_case_summary.json`](Project/data/output/high_fidelity_case_summary.json)
and its sampled trajectory in
[`data/output/high_fidelity_case_trajectory.csv`](Project/data/output/high_fidelity_case_trajectory.csv).

### 3. Reproduce the `Project/` visualizations and report

The `core` and `science` profiles generate their corresponding plots under
`Project/figures/`. After downloading the scientific data, the terrain
visualizations can also be regenerated directly:

```bash
python scripts/13_terrain_visualizations.py
```

To rebuild the report:

```bash
cd docs
latexmk -pdf main.tex
```

### 4. Run the `Web/` solver locally

The committed Git LFS data is sufficient to run the browser solver. From the
repository root:

```bash
cd Web
npm ci
npm run dev
```

Vite prints the local URL in the terminal. To reproduce the production build
and preview it locally:

```bash
npm run build
npm run preview
```

To regenerate the browser terrain and gravity products from the authoritative
NASA sources rather than use the committed LFS objects:

```bash
cd ../Project
python scripts/download_science_data.py
python scripts/generate_all.py --profile web
```
