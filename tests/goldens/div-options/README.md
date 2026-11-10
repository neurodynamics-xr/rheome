# Divergence / curl option goldens (written by nxr; read by Rheome and the Cortical Flow Editor)

Regenerate: `python tests/goldens/div_options_goldens.py`. Format as `tests/goldens/vortex/`: `manifest.json` lists every array (file, dtype, shape, unit, note), the nxr commit and source file:line; `.bin` are little-endian row-major, no header.
Ellipsoid 70/55/45 mm (icosphere-4, 2562 vertices). **Positions: nxr runs on `V_mm/1000` (metres); divergence, curl and H are per METRE and a port reproduces them only with `V_mm/1000`** (with `V_mm` they come out 1000x smaller). Positions are cast to float32 first.
Identity (exact to rounding): `div3d = div_tan + div_curv`, `curl3d = curl_tan + curl_curv`. `div_curv` is div3D(J_n n); `div_curv_2H` = 2 H J_n is its continuum estimate. Counterparts: Editor `flowFieldCpu` `tangential` (true = tan, false = 3d); Rheome `rheome.differential.curl/divergence` 4th argument.
