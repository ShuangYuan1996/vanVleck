# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-11-23

First public release.

### Added
- Reading of CIF (several data blocks) and XYZ (several frames) files; symmetry expansion of the
  asymmetric unit and growth of complete molecules.
- 3D viewer (ball-and-stick, capped sticks, wireframe, spacefill), unit cell, view along a/b/c.
- Multi-file workspace with removal of files.
- Selection of a complex and of its inner coordination sphere (metal + 6 donor atoms).
- Point-group assignment of the ML6 sphere, including the variants D2h^d/D2h^v, C2v^d/C2v^v,
  C2^z/C2^xy, Cs^z/Cs^xy, and drawing of the symmetry elements.
- Symmetry-adapted standard frame and van Vleck coordinates Q1–Q15 with user-defined masses,
  |Q|, and a check of all zero / equal / opposite relations required by the point group.
- OctaDist-type distortion parameters: <d>, zeta, Delta, Sigma, Theta, Theta_min, Theta_max and
  octahedral volume.
- Compare window: consistent frames for molecules of the same material, RMSD to the reference,
  grouped bar chart, export of Q to tab-separated text.
- Synthetic example file and regression tests.
