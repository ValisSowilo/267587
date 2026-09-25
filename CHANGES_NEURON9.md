# NEURON 9 port

The NEURON 9 port of this model (ModelDB #267587, Guet-McCreight et al. 2022) was made by
**Valis Sowilo** ([github.com/ValisSowilo](https://github.com/ValisSowilo)), September 2026.
The model itself, its parameters and the original code are by the authors listed in README.md.

- Based on: ModelDB #267587 as published at github.com/ModelDBRepository/267587 (commit 82cdd91)
- Source of this port: github.com/ValisSowilo/267587
- License: GPL-3.0, like the original (see LICENSE). Every file changed by the port carries a
  "Modified 2026-09-24 by Valis Sowilo" notice, as section 5(a) of the license requires.
- NEURON9_PORT_SHA256SUMS.txt lists SHA-256 checksums of every file changed or added by the
  port, so a copy can be checked against this release with `sha256sum -c NEURON9_PORT_SHA256SUMS.txt`.

## Changes

Mod files (all five mod/ folders): Gfluct.mod, ProbAMPANMDA(2).mod, ProbUDFsyn(2).mod
- Removed "RANGE new_seed" from Gfluct.mod (new_seed is also a PROCEDURE, which NEURON 9 rejects).
- Removed the VERBATIM forward declarations of nrn_random_pick / nrn_random_arg, which conflict
  with NEURON 9's own headers. No change in behaviour.

circuit.py (default_circuit, heterogeneous_circuit, PSP Simulations)
- LFP and dipole recording moved to LFPy 2.1+ probes (RecExtElectrode, CurrentDipoleMoment);
  the saved SPIKES, OUTPUT and DIPOLEMOMENT files keep their previous format.
- SomaAsPointElectrode reproduces the LFPy 2.0 method "soma_as_point" exactly.
- cellnums is set after the TESTING block, and plots are skipped in test runs, so TESTING = True works.

circuit_functions.py and L5Circuit Analyses
- Four-sphere EEG calls updated to the LFPy 2.1+ interface.
- np.trapz -> scipy.integrate.trapezoid; DataFrame.append -> pd.concat; invalid escape sequence fixed.

PSP Simulations/test_syn_Fig7.py
- Network.simulate(electrode=) -> probes=; cell.ymid/zmid -> cell.y/z; np.trapz and DataFrame.append
  as above; the output folder is created if Circuit_output/ does not exist.

Population Analysis AnalyzeResults_FigS4top.py / AnalyzeResults_FigS4bottom.py
- set_xticks([a, b]) with one-element arrays -> set_xticks(numpy.ravel([a, b])): current matplotlib
  requires 1-D tick positions. The ticks are at the same places as before.

Optimizations/L5Pyr_young_AltMorph_Scinet and L5Pyr_old_AltMorph_Scinet
- init_5recordings.py: numpy.int (removed in numpy 1.24) -> int.
- job.sh, job_debug.sh: create $SCRATCH/.ipython before "ipython profile create", which current
  IPython no longer does itself.

circuit.py (all three): ranks without cells
- NEURON 9 aborts with "nrn_calc_fast_imem: Assertion `vec_sav_d' failed" on an MPI rank that owns
  no cells, which happens when there are more ranks than cells (e.g. TESTING = True with more than
  4 ranks). Such a rank now creates one unconnected placeholder section. It belongs to no
  population and is not recorded, and ranks that own cells are unaffected.

Windows (circuit.py in all three folders, PSP Simulations/test_syn_Fig7.py)
- On Windows, np.random.randint without a dtype returns 32-bit integers, so LFPy 2.3.7 fails in
  Network.connect ("high is out of bounds for int32") when it seeds the stochastic synapses'
  random number generators with np.random.randint(0, 2**32 - 1). On Windows only, these scripts
  now make np.random.randint default to 64-bit integers, as on Linux and macOS. This also makes the
  random numbers the same as on Linux, so a seed gives the same connectivity (and, in the
  heterogeneous circuit, the same choice of cell models) on every system.
- The "Mechanisms found" message also recognizes mod/nrnmech.dll (Windows) and mod/arm64/special
  (Apple silicon), instead of printing False there.

Bugs in the original code (not related to NEURON 9)
- PSP Simulations/circuit.py: agegroup = 'o_rescue' referred to a model file that does not exist;
  now 'y' (the options are 'y' and 'o').
- Population Analysis/Morph2/SimulateModel.py loaded HL5PN1.swc, which is only in Morph1; it now
  loads HL5PN2.swc, the morphology in its own folder.
- out_1SimulateModel.py (Scinet folders) used an undefined variable, single_cell_data; it now
  checks target_feature_type == 'Automatic', as init_8plot.py does.
- out_2SimulateHocModel.py (Scinet folders) hard-coded the template name 'interneuron', but these
  folders build 'pyramidal' cells; the name is now read from init_1morphology.py.
- out_3AnalyzeResults.py (Scinet folders) and Population Analysis AnalyzeResults_FigS4*.py: the
  numbered hall-of-fame labels in Quality_SagVsRMP are now drawn with clip_on=True. Before, a label
  far outside the axes (a hall-of-fame model with a large RMP or sag error; Population Analysis fixes
  the x axis at 0-0.9 SD) made savefig(..., bbox_inches='tight', dpi=300) enlarge the image to
  include it, which could need many gigabytes of memory (11 GB in one test). Labels inside the axes
  are unchanged: on normal data the saved figure is pixel-identical to before.
  Note: out_3AnalyzeResults.py still needs at least one model within 2 SD on every feature (as after
  a full optimization); with none it stops with an IndexError, as in the original.
- The output folders that scripts write into (figs*, results, PLOTfiles, work,
  highly_ranked_models, output_readout, ...) are now included as empty folders (with a .gitkeep
  file), so a fresh copy runs without "No such file or directory" errors.

## Installation

- Linux: `pip install -r requirements_neuron9.txt`. macOS also has NEURON packages on PyPI and
  should work the same way, but was not tested.
- Windows: PyPI has no NEURON package for Windows, so NEURON comes from its installer and LFPy is
  installed with --no-deps. The steps are at the top of requirements_neuron9_windows.txt.
- Single-cell code (Optimizations, Population Analysis): also install
  requirements_neuron9_singlecell.txt, as described at the top of that file.

New files: requirements_neuron9.txt, requirements_neuron9_windows.txt,
requirements_neuron9_singlecell.txt, CHANGES_NEURON9.md, NEURON9_PORT_SHA256SUMS.txt, and a
NEURON 9 section in README.md.

## Validation (NEURON 9.0.2, LFPy 2.3.7, Python 3.12)

- All 12 mod folders compile.
- Single-cell F-I curves (L5Pyr_old, L5Pyr_young) are identical to the bundled simdata/HL5PN1_FI.txt.
- PSP simulations (Fig. 7): voltage traces identical to the original code run on NEURON 8.2.7 with
  LFPy 2.0.7, for all six connection/age blocks.
- Circuit (20-cell reduced network): on the same NEURON version, spike times and dipole moments are
  identical to the original code with LFPy 2.0.7 when the same connectivity sampling is used.
- Full-length (7 s) circuit run produces all output files and figures.

For a given random seed, circuit results on NEURON 9 are not identical to the published runs: NEURON 9
removed the default generator of `new Random(seed)` (used for the background noise in net_functions.hoc),
and LFPy 2.1+ samples the random connectivity differently. See README.md.
