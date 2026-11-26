# Code and data for Paper 3: numerical parameters of time-stepping stability analysis

This repository contains the code and the data used in the paper *Amplitude, time step, sampling and field weighting in time-stepping stability analysis with a sharp interface* (A. Okojunu, M. Sussman and U. S. Nair). Every table of the paper can be recomputed from it, and every simulation can be repeated with it.

The simulations were run with the two-phase solver `amrex_implicit_interfaces` of M. Sussman, at a fixed commit, and analyzed with the code in the folder `LSA_Rayleigh_Fortran/`, which is placed inside the solver tree to run them. The runs with the restarting variant also need a change to the solver, given as a patch. Levels 1 and 2 below need neither the solver nor a cluster.

## Contents

| Path | Contents |
|---|---|
| `LSA_Rayleigh_Fortran/src`, `tests`, `Makefile` | the analysis code (Fortran 2008, LAPACK): the input decks, the twin differences, the interface-amplitude fit, the dynamic mode decomposition (DMD) and the dispersion relations, with its test suite |
| `LSA_Rayleigh_Fortran/scripts`, `plotting` | the scripts that run the simulations under Slurm and analyze them |
| `LSA_Rayleigh_Fortran/krylov_test` | the scripts of the restarting variant, among them `krylov_fix.sh`, which runs it with the changes of the patch |
| `LSA_Rayleigh_Fortran/paper3/code` | the scripts written for this paper: the local growth rate in bands of amplitude and the table of the field weighting |
| `LSA_Rayleigh_Fortran/paper3/data` | the runs made for this paper: `2026-11-24_paper3` (low Ohnesorge numbers and fields stored every four steps) and `2026-11-25_sets34` (runs started from the leading mode, and the restarting variant with its changes) |
| `LSA_Rayleigh_Fortran/results` | the stored output of the dissertation runs that the paper uses |
| `LSA_Rayleigh_Fortran/*.txt` | the case files, one line per simulation |
| `LSA_Rayleigh_Fortran/reproduce_paper3.sh` | the commands of the three levels below |
| `solver_patches/lsa_fix_restart.patch` | the change to the restarting driver of the solver (Level 3 only) |

## Level 1. Check the code (about one minute)

On Linux, or WSL on Windows, with gfortran, make and LAPACK:

```bash
git clone https://github.com/abelokoj/lsa-timestepping-parameters.git
cd lsa-timestepping-parameters/LSA_Rayleigh_Fortran
./reproduce_paper3.sh 1
```

The test suite must end with `ALL TEST PROGRAMS PASSED`. The script then prints the unbounded two-fluid relation, the reference of Tables 1 to 3: 0.328689, 0.242064 and 0.174896 at k R0 = 0.7, 0.9 and 0.95 for Oh = 0.02, and 0.341648 at k R0 = 0.7 for Oh = 0.002.

## Level 2. Recompute the tables from the stored data (about one minute)

```bash
python3 -m pip install -r plotting/requirements.txt
./reproduce_paper3.sh 2
```

The output is written to `figures/paper3_tables.txt`. Some of it is labeled with the table numbers of the dissertation on which the paper is based.

| In the paper | Printed by Level 2 as |
|---|---|
| Table 1 and Section 4.1, the amplitude | dissertation Table 6.2, then the smallest amplitude refitted over the window of the largest |
| Section 4.2, the windows | the refitted window of each run started from rest or from the leading mode, beside the stored one |
| Table 2, the local rate | the output of `paper3/code/band_rates.py`, compared with the stored output |
| Table 3 and Section 4.3 | the refitted interface rates of sets 1 and 2, the stored decomposition rates, and the rates at a fixed rank over the strides |
| Table 4 and Section 4.4, the weighting | the output of `paper3/code/weighting_table.py`, and the stored test on the runs of sets 1 to 3 |
| Table 5 and Section 4.5, the budget | the rows at k R0 = 0.7 of dissertation Tables 5.3 and 6.1, and the changes due to the amplitude, grid, time step and wall |
| Section 4.6, the restarting variant | the member rate and the decomposition rates of each run of set 4 and of the runs started from the leading mode |

Level 2 refits every interface growth rate from its amplitude record, recomputes the local rates and the weighting table, and prints the stored output of the decompositions. It does not repeat the decompositions, because their inputs, the twin-difference fields and the members of the restarting variant, are too large for this repository. They are repeated in Level 3. The rates over the strides are printed as rates, with their spread and their largest change from the rate at stride 1, since the paper uses both.

The last column of Table 1 is recomputed by fitting the record of the smallest amplitude over the snapshots of the window of the largest, both ends included. This gives −0.34%, −0.55% and −0.76% on 8, 16 and 32 cells per radius; the value on 8 cells is obtained without the last snapshot of that window.

## Level 3. Run the simulations again

The runs were made on the FSU Research Computing Center cluster with `gnu/13` and `openmpi/4.1.6`, under Slurm. They used commit `dddb827c` of a fork of the solver, whose solver source is identical to commit `0b8df42` of the original repository; the run records give the commit, the number of ranks and the checksum of the executable for every run. The restarting variant and the runs started from the leading mode used the same source with the patch in `solver_patches/`, built in a separate tree so that the executable of the other runs is not changed.

```bash
export AII=$HOME/amrex_implicit_interfaces
git clone https://github.com/msussman42/amrex_implicit_interfaces.git $AII
cd $AII && git checkout 0b8df42
cp -r /path/to/lsa-timestepping-parameters/LSA_Rayleigh_Fortran $AII/NS_fluids_lib/
module purge; module load gnu/13; module load openmpi
cd $AII/NS_fluids_lib && make -j16              # the unmodified solver, amr2d.gnu.MPI.ex
git -C $AII worktree add ~/aii_fix 0b8df42       # a second tree for the patched solver
cd ~/aii_fix && git apply /path/to/lsa-timestepping-parameters/solver_patches/lsa_fix_restart.patch
cd NS_fluids_lib && make -j16
cd $AII/NS_fluids_lib/LSA_Rayleigh_Fortran
ACCOUNT=<your Slurm account> FIX_SOLVER_EXE=$(ls ~/aii_fix/NS_fluids_lib/amr2d*.ex) ./reproduce_paper3.sh 3
```

The patch adds two inputs to the restarting driver, `ns.LSA_fix_restart` and `ns.LSA_levelset_scale`, whose defaults leave the solver unchanged. Use commit `0b8df42` and not the current head: at some later commits a calculation with fluids only stops at start-up with `expecting num_FSI_outer_sweeps>=2` (`LSA_Rayleigh_Fortran/docs/solver_notes.md`).

Level 3 submits the case files of the dissertation runs that Tables 1, 2, 4 and 5 use, the two case files of sets 1 and 2 (`campaign_p3_16.txt`, `campaign_p3_16_str04.txt`), and the seven runs of sets 3 and 4 with `krylov_test/krylov_fix.sh`, with the arguments recorded in their manifests in `paper3/data/2026-11-25_sets34/krylov_test/`. The number of ranks and the wall-clock limit of each case file are set in `scripts/submit_campaign.sh`. Each perturbed run has an unperturbed twin, and the two must run on the same number of ranks; the case files and the scripts ensure this. When every job has ended, `./reproduce_paper3.sh 3-analyse` runs the studies, the decompositions over the strides and the weighting tests on the new runs. The decompositions of the runs started from the leading mode over 0.01 < a ≤ 0.05 R0 (Section 4.2) were made with `run_lsa analyse --from-export` over the windows t = 5.73 to 12.24 (k R0 = 0.9) and 7.89 to 16.83 (k R0 = 0.95). Bitwise reproducibility has not been tested. Another compiler, MPI library or rank count changes the round-off, so the interface rates are expected to agree closely with the stored ones, and the decomposition rates, which come from the difference of two runs, less well.

## Origin

The analysis code is that of the dissertation of A. Okojunu (Florida State University, 2026), published at https://github.com/abelokoj/LSA_Rayleigh_Fortran, with the runs of sets 1 to 4, the changes to the restarting driver and the scripts in `paper3/code` added for this paper.
