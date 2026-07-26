# PHYS1600J: The Orbital Home Run

This repository contains the scientific model and static trajectory solver for
PHYS1600J Problem B.

- [`Project/`](Project/) contains the Python analysis, versioned derived data,
  figures, and an archived paper build.
- [`Web/`](Web/) contains the bilingual client-side solver deployed by GitHub
  Pages.
- [`.github/workflows/deploy-pages.yml`](.github/workflows/deploy-pages.yml)
  builds `Web/` and deploys the static artifact.

The current realistic workflow uses LOLA LDEM64 terrain, GRGM1200B gravity
truncated at a convergence-tested degree, fixed short-arc Earth and solar
tides, and explicit numerical cross-checks. Large authoritative NASA source
files are downloaded and verified by script rather than committed to Git.

See [`Project/README.md`](Project/README.md) for scientific reproduction and
[`Web/README.md`](Web/README.md) for the browser model and local build.
