# Rheome

**A MATLAB toolbox for the multiscale geometry of human cortical dynamics.**

*rheo-* (flow, current) + *-ome* (the complete set, as in connectome). Rheome measures spatiotemporal
dynamics on graphs (cortical meshes, connectomes and sensor arrays) with graph wavelet frames, the
joint (λ, ω) plane, phase geometry and cortical flow from MEG source estimates.

The core is three filterbank classes modelled on MATLAB's `cwtfilterbank`, so the names carry
over (`scales`, `centerFrequencies`, `powerbw`, `qfactor`, `framebounds`, `wt`/`iwt`):

| class | what it is |
|:--|:--|
| `rheome.graphfilterbank` | spectral filterbank on a self-adjoint operator (a mesh Laplacian, a connectome, a sensor graph) |
| `rheome.jointfilterbank` | filterbank on the joint time-vertex plane (λ, ω): spatial scale, rate and speed |
| `rheome.timefilterbank`  | a constant-Q tight frame on the frequency axis, each member at its own rate |

Around them are surface operators (`rheome.operators`, `rheome.eigen`) and kernels (`rheome.filters`).
Then the differential geometry of flow fields (`rheome.differential`, `rheome.detect`) and the
kinematic route (`rheome.dynamics`). Then MEG source mapping (`rheome.forward`, `rheome.inverse`,
`rheome.source`, `rheome.flow`) and sensor arrays (`rheome.sensors`).

**One namespace.** Every function and class lives under `rheome.*`, so nothing collides with
MATLAB's own `+utils`, `+io` or `+filters`, or with another toolbox. `help rheome.flow` lists a
package, and the root `Contents.m` lists them all.

## Install

**Toolbox (recommended).** Download `rheome-1.0.0.mltbx` from the
[release page](https://github.com/neurodynamics-xr/rheome/releases) and double-click it, or:

```matlab
matlab.addons.install('rheome-1.0.0.mltbx');
```

**From source.** Clone the repository and put its root on the path (the `+rheome` folder resolves
from there):

```matlab
addpath('/path/to/rheome');
```

Requires MATLAB R2023b or later. See *Requirements* below for the add-on toolboxes.

## Quick start

A Mexican-hat graph wavelet bank on a sphere, and the energy per scale of a polar cap:

```matlab
[V, F] = rheome.geom.icosphere(4);                              % 2562-vertex unit sphere
[L, M] = rheome.operators.laplace_beltrami(V, F, 'galerkin');   % cotangent Laplace-Beltrami + mass
gfb = rheome.graphfilterbank.fromOperator(L, 'Mass', M);        % bank on L's spectral range (Chebyshev)
x   = double(V(:,3) > 0.9);                                     % a cap at the north pole
E   = vertexSpectrum(gfb, x)                                    % energy per scale  [M x 1]
```

Each demo plants a feature whose answer is known analytically, prints PASS/FAIL lines, and
returns `out.ok`:

```matlab
out = rheome.demos.filterbank_graph;     % bumps of known width on a sphere, recovered per scale
out.ok
rheome.demos.filterbank_joint            % a rigidly rotating pattern at known speed
rheome.demos.filterbank_cwt              % the same five figures with MATLAB's cwtfilterbank (Wavelet Toolbox)
rheome.demos.sensor_wavelets             % scale, rate and speed planted on a sensor array, then recovered
```

`help rheome.demos` lists the rest. Pass a folder to export a demo's figures:
`rheome.demos.filterbank_graph('figs')`.

⚠ The graph demo recovers **√2·σ**, not σ: `vertexSpectrum` is an *energy* marginal, so the
relevant integral is of the squared gain. `help rheome.demos` states all three calibrations.

## Documentation

The user guide opens in MATLAB's Help browser once the toolbox is installed, under
*Supplemental Software > Rheome* (the Getting Started guide also opens from the Add-On Manager). It contains a tutorial (Getting Started), how-to guides
(importing a Brainstorm study, flow maps, per-participant measures, figures), explanations
(sensors to cortex, the Helmholtz-Hodge split, surface optical flow, scale, the gauge, the
resolution floor) and a function reference generated from the help text. Every example runs on a
synthetic Brainstorm study, so no participant data is needed. The pages are also in
[`doc/`](doc/) in this repository.

## Data and outputs

Two folders, both outside the installed toolbox:

| | holds | default (installed toolbox) | override |
|:--|:--|:--|:--|
| `rheome.load.root()` | cached inputs, regenerable: surfaces, eigenbases, studies (`rheome.import.*` writes, `rheome.load.*` reads) | `<userpath>/rheome/data` | `setenv('RHEOME_DATA', …)` |
| `rheome.load.outroot()` | what analyses found: tables, `.mat` results, figures, reports (`rheome.load.outpath`) | `<userpath>/rheome/results` | `setenv('RHEOME_OUT', …)` |

In a git clone the defaults are `+data/` and `results/` beside the code (both gitignored).
Brainstorm data enter through `rheome.import.*`, e.g.
`rheome.import.dataset(name, studyDir, cortexFile, dataName)`; `rheome.load.list` shows what is cached.

**The HCP-1065 connectome template is a separate download.** `rheome.connectome.resolve` falls back
to a group-average tractography template when a subject has no tractography of its own. The
template is derived from WU-Minn HCP data, which may be redistributed only under the
[HCP Open Access Data Use Terms](https://www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms),
not under the MIT licence. So `hcp_template.mat` is attached to the release on its own. Put it at
`fullfile(rheome.load.root(), 'hcp_template.mat')`, or point `RHEOME_HCP_TEMPLATE` at it.
`rheome.connectome.build_hcp_template` rebuilds it from your own copy.

## Requirements

MATLAB R2023b or later. The filterbank classes, the surface operators and the quick start need
base MATLAB only. Other parts use add-on toolboxes (measured with
`matlab.codetools.requiredFilesAndProducts`):

| toolbox | used by |
|:--|:--|
| Signal Processing | `rheome.flow` (envelopes, vortex tracking), `rheome.scale`, `rheome.forward.simulate`, `rheome.dynamics.pde_fit`, `rheome.filters.firbandpass`, `rheome.detect.spindle`, `rheome.flowpage` |
| Statistics and Machine Learning | `rheome.flow`, `rheome.scale`, `rheome.inverse.resolution`, `rheome.forward`, `rheome.detect.peaks`/`onstate`, `rheome.connectome.resolve`, `rheome.operators.connectome`, `rheome.sensors.tree` |
| Wavelet | `rheome.flowpage`, `rheome.ingest.bank`/`reducepaged`, `rheome.demos.filterbank_cwt` (which skips cleanly without it) |
| Optimization | `rheome.detect.risefall` |
| MATLAB Report Generator | `rheome.report` (optional: results written up as a mini-article) |

Brainstorm itself is not needed: its files are read from disk as plain `.mat`.

## Citation

If you use Rheome, please cite it (see [`CITATION.cff`](CITATION.cff)):

> Basha, D., & Baillet, S. (2026). *Rheome: a MATLAB toolbox for the multiscale geometry of human
> cortical dynamics* (Version 1.0.0) [Computer software]. https://github.com/neurodynamics-xr/rheome

A DOI will be added here when the release is archived on Zenodo.

## For developers

```matlab
addpath(pwd); addpath(fullfile(pwd, 'tests'), fullfile(pwd, 'build'));
runtests('tests')                 % run headless: matlab -nodisplay -batch "..."
build_toolbox                     % -> build/out/rheome-<version>.mltbx, from a clean commit
```

Run the suite headless (`-nodisplay`): several demo tests open five figures each. A few tests
check results on a real cortex and are skipped unless `RHEOME_TEST_SUBJECT` names a subject in the
data cache (see `tests/rheomeTestSubject.m`).

The version lives in the root `Contents.m`. `build/check_install.m` installs a built `.mltbx` into a
clean MATLAB path and runs the quick start. The documentation is rebuilt by `doc/tools/build_docs.m`
and its examples checked by `doc/tools/run_doc_examples.m`. The repository holds what the add-on
leaves out: `tests/`, `build/`, `doc/source/` and `doc/tools/`.

## Licence

MIT; see [LICENSE](LICENSE). The separately distributed HCP-1065 template is under the HCP Open
Access Data Use Terms.

## Authors

Diellor Basha (McGill University; CRCHUM, Université de Montréal) and Sylvain Baillet
(Université de Montréal; CRCHUM; McGill University).
