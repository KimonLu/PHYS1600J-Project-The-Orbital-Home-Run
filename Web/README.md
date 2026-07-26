# Orbital Home Run Web Solver

This static client-side solver accepts a body-fixed lunar launch site and a
surface-relative velocity vector, then returns:

- an inertial-frame 3D trajectory reconstructed from the rotating-frame result;
- an equirectangular ground track;
- LDEM64 terrain clearance and first-impact coordinates;
- the closest post-launch approach to the moving launch site;
- a return, uncertain-return, impact, escape, or no-return classification.

## Interface

- Switch between English and Chinese without clearing the current result.
- Select verified examples for all five trajectory classifications.
- Focus the 3D canvas and use WASD to pan; drag to orbit and scroll to zoom.
- Play the ball and lunar rotation on one physical simulation clock from 1x to
  1000x real time.
- Scroll and drag the ground-track map; double-click, or press 0 while focused,
  to reset it.
- Open the model badges for LOLA LDEM64, GRAIL GRGM1200B, RK4 step-doubling,
  and the NASA Scientific Visualization Studio Moon Kit.

## Scientific model

- Terrain uses the global NASA LOLA LDEM64 grid at 64 pixels per degree. The
  browser tiles preserve every elevation sample.
- Gravity uses a browser acceleration-correction atlas generated from
  GRGM1200B truncated at degree and order 600, with the central term evaluated
  analytically.
- Rotating-frame dynamics include Coriolis and centrifugal terms plus optional
  fixed-geometry differential Earth and solar tides.
- Integration uses fourth-order Runge-Kutta with step doubling. Its result is a
  numerical diagnostic, not a rigorous physical error bound.

The assumed total position-error bound is an external input, not a quantity
calculated by the solver. Classification uses the larger of that assumed bound
and the numerical diagnostic. `RETURN` is therefore conditional on the stated
input assumptions; a return-radius interval crossing the threshold is
`RETURN_UNCERTAIN`.

The current Earth and Sun controls only enable or disable fixed short-arc
geometry. They do not select a direction or a physical epoch. A real
epoch-dependent trajectory would require a complete SPICE ephemeris, frame,
attitude, and libration chain.

## Local build

```text
npm install
npm run dev
npm run build
```

The Vite base is relative, so `dist/` works locally and under a GitHub Pages
repository subpath. `.github/workflows/deploy-pages.yml` builds and deploys
`Web/dist` after GitHub Pages is configured to use GitHub Actions.

## Limitations

This is a course-project research tool, not a navigation, landing, or safety
product. Rock-scale hazards, gravity covariance, launch-state covariance,
solar radiation pressure, and full SPICE geometry are outside the browser
model.
