# Orbital Home Run — Web solver

An English, static, client-side solver for the PHYS1600J “The Orbital Home
Run” project. It accepts a body-fixed lunar launch site and a surface-relative
velocity vector, then returns:

- an inertial-frame 3D trajectory reconstructed from the rotating-frame solution;
- an equirectangular ground track;
- LDEM64 terrain clearance and first impact coordinates;
- the closest post-launch approach to the moving launch site;
- a return / uncertain return / impact / escape / no-return classification.

## Interface

- Switch between English and Chinese without rerunning or clearing a result.
- Select verified examples for all five physical trajectory classifications.
- Focus the 3D canvas and use WASD to pan; drag to orbit and scroll to zoom.
- Play the ball and lunar rotation on one physical simulation clock from 1× to
  1000× real time. The inertial trajectory and the Moon use the same sidereal
  rotation rate; the marker sizes are visual aids rather than physical scale.
- Scroll and drag the equirectangular ground-track map to inspect close passes
  and terrain intersections. Double-click, or press 0 while focused, to reset
  the map.
- Open the model badges for the LOLA LDEM64 archive, the GRAIL GRGM source, and
  RK4 step-doubling documentation.

## Scientific model

- **Terrain:** global NASA LOLA LDEM64, 64 pixels per degree (about 474 m at
  the equator). The original 531 MB PDS image is transformed into reversible,
  horizontally delta-coded 10° gzip tiles. No elevation samples are discarded.
- **Gravity:** NASA GRAIL GRGM1200B truncated at degree and order 600. Degree
  600 was selected by comparison with degree 1200 at 12, 15, and 30 km
  reference altitudes. The browser streams a 0.125° × 0.125°,
  altitude-layered acceleration-correction atlas and evaluates the central
  term analytically.
- **Dynamics:** body-fixed equations include lunar gravity, Coriolis and
  centrifugal terms, plus optional differential Earth and solar tides. The
  trajectory uses fixed short-arc third-body geometry.
- **Integrator:** fourth-order Runge–Kutta with step doubling. The accumulated
  local difference is reported as a conservative numerical diagnostic.

The 15 km and 30 km cross-checks differ from direct degree-600 propagation by
6.06 m and 4.04 m after one circuit. A 12 km reference trajectory falls below
the conservative Brillouin sphere and differs by 24.0 m; the UI explicitly
flags this domain. These are numerical model comparisons, not physical
confidence intervals.

The reversible terrain tiles are also checked at random points, at every
10-degree tile seam, and across the 0/360-degree meridian. Their bilinear
queries reproduce direct LDEM64 values exactly in the validation sample.

## Run locally

```text
npm install
npm run dev
```

Build the static site with:

```text
npm run build
```

The Vite base is relative, so `dist/` works both locally and under a GitHub
Pages repository subpath.

## GitHub Pages deployment

The repository workflow builds `Web/` and deploys `Web/dist`. Before its first
run, a repository administrator must enable the Pages site once:

```text
Settings > Pages > Build and deployment > Source > GitHub Actions
```

This cannot be bootstrapped by the workflow's default `GITHUB_TOKEN`.
Subsequent pushes to `main`, or a manual `workflow_dispatch`, use the normal
`configure-pages`, `upload-pages-artifact`, and `deploy-pages` pipeline.

## Data provenance and limitations

The lunar color and display-height textures are from the
[NASA Scientific Visualization Studio CGI Moon Kit](https://svs.gsfc.nasa.gov/4720).
The collision surface is the scientific LDEM64 binary product, not the display
texture. GRGM1200B is provided by NASA GSFC's Planetary Geodynamics Data
Archive.

The solver is a course-project research tool, not a navigation, landing, or
safety product. Rock-scale hazards, gravity covariance, launch-state
covariance, solar radiation pressure, and a complete SPICE frame/ephemeris
treatment are outside the browser model.
