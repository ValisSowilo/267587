# Human L5 Cortical Circuit Model --- Guet-McCreight-et-al.-2022
==============================================================================
Author: Alexandre Guet-McCreight

This is the readme for the model associated with the paper:

Guet-McCreight A, Chameh HM, Mahallati S, Wishart M, Tripathy SJ, Valiante TA, Hay E (2022) Age-dependent increased sag amplitude in human pyramidal neurons dampens baseline cortical activity. Cerebral Cortex.


Network Simulations:
Simulation code associated with the default L5 circuit used throughout the manuscript is in the /L5Circuit/default_circuit/ directory. Additional tests using heterogeneous circuits can be found in the /L5Circuit/heterogeneous/ directory.

To run simulations, install all of the necessary python modules (see lfpy_env.yml), compile the mod files within the mod folder, and submit the simulations in parallel (e.g., see job.sh). Running younger and older circuits is controlled for by changing the agegroup variable in circuit.py to either 'y' or 'o', e.g.:

agegroup = 'y' or 'o'

In job.sh, the 2nd to last number in the mpiexec command (see below) controls the random seed used for the circuit variance (i.e., connection matrix, synapse placement, etc.), while the last number controls the random seed number used for stimulus variance (i.e. Ornstein Uhlenbeck noise and stimulus presynaptic spike train timing).

mpiexec -n 400 python circuit.py 1234 1

All code used for analyzing the circuit simulation results is found in the /L5Circuit Analyses/ directory. These codes perform analyses of multiple simulations using several circuit random seeds (e.g. fig 3) or stimulus random seeds (e.g. fig 5).


Single Cell Optimizations and Simulations:
All code associated with single cell optimizations is found in the /Single Cell Modelling/Optimizations directory. For usage of this code, we recommend starting from the README.md file and associated code in either the /L5Pyr_young_AltMorph_Scinet/ or /L5Pyr_old_AltMorph_Scinet/ directories.

Analysis code of single cell optimization results can be found in the /Single Cell Modelling/Population Analysis/ folder. Current-step single-cell simulation code can be found in the /Single Cell Modelling/Current-Step Simulations/ folder. Code for simulating PSP summation of different connection types can be found in the /Single Cell Modelling/PSP Simulations/ folder.


Running on NEURON 9 (updated 2026):
NEURON 9 port by Valis Sowilo (github.com/ValisSowilo), 24-25 September 2026. See CHANGES_NEURON9.md for the full list of changes and how they were validated. The model and the port are distributed under the GNU General Public License v3.0 (see LICENSE).

The code has been updated to run on NEURON 9 with current Python packages (tested with Python 3.12, NEURON 9.0.2 and LFPy 2.3.7; see requirements_neuron9.txt). The original NEURON 7.7 / LFPy 2.0 environment is still described in lfpy_env.yml.

pip install -r requirements_neuron9.txt
cd L5Circuit/default_circuit/mod && nrnivmodl && cd ..
mpiexec -n 400 python circuit.py 1234 1

On Windows, install NEURON 9.0.2 with its installer and Microsoft MPI first, then follow the steps at the top of requirements_neuron9_windows.txt (PyPI has no NEURON package for Windows, so LFPy is installed with --no-deps). Before running nrnivmodl on Windows, set MAKEFLAGS=EXTRA_FLAGS=-O2: NEURON 9.0.2 compiles the mod files without optimization on Windows, which made simulations up to about 4 times slower. The single-cell optimization and population analysis code also needs requirements_neuron9_singlecell.txt.

Changes made for NEURON 9 (model behaviour is unchanged):
- Mod files (Gfluct.mod, ProbAMPANMDA(2).mod, ProbUDFsyn(2).mod): removed the "RANGE new_seed" line (new_seed is also a PROCEDURE, which NEURON 9 rejects) and two VERBATIM forward declarations of nrn_random_pick / nrn_random_arg that conflict with NEURON 9's own headers.
- circuit.py (both circuits and PSP Simulations): LFPy 2.1+ records the LFP and the current dipole moment through "probes" (RecExtElectrode and CurrentDipoleMoment) instead of return values of Network.simulate(). The probe data is repackaged into the same OUTPUT and DIPOLEMOMENT formats as before, so the saved .npy files and the L5Circuit Analyses scripts are unchanged. LFPy 2.0's method="soma_as_point" is reproduced exactly by the SomaAsPointElectrode class (LFPy's newer "root_as_point" would only treat the first soma on each MPI rank as a point source).
- circuit.py: cellnums is now set after the TESTING block, so TESTING = True runs with 1 cell per population instead of failing (the plots, which assume the full-length simulation, are skipped in test runs).
- circuit_functions.py and L5Circuit Analyses: updated the four-sphere EEG calls to the LFPy 2.1+ interface (FourSphereVolumeConductor(electrode_positions, radii=..., sigmas=...).get_dipole_potential(p.T, location)).
- Python package updates: np.trapz -> scipy.integrate.trapezoid, DataFrame.append -> pd.concat, cell.ymid/zmid -> cell.y/z, and the 'electrode' argument of Network.simulate() -> 'probes'.

Validation: all 12 mod folders compile on NEURON 9.0.2. The single-cell F-I curves (L5Pyr_old and L5Pyr_young, step_current.hoc) are identical to the bundled simdata/HL5PN1_FI.txt. On the same NEURON version, the updated circuit code with LFPy 2.3.7 reproduces the original code with LFPy 2.0.7 exactly (identical spike times and dipole moments), apart from the one LFPy difference noted below.

Differences from the original simulations for a given random seed:
- NEURON 9 removed the old default generator of new Random(seed) (ACG), which net_functions.hoc uses to seed the Ornstein-Uhlenbeck background noise. The noise has the same statistics, but a given seed produces a different noise realization than on NEURON 7/8, so spike times for a given seed are not identical to the published runs.
- LFPy 2.1+ draws the random connectivity with np.random.binomial instead of np.random.rand < p. The connection probabilities are unchanged, but a given circuit seed produces a different connection matrix than LFPy 2.0.
Results should therefore be compared across seeds, as in the paper, rather than seed by seed.
