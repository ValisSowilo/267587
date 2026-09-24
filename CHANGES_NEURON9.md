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

New files: requirements_neuron9.txt, CHANGES_NEURON9.md, NEURON9_PORT_SHA256SUMS.txt, and a
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
