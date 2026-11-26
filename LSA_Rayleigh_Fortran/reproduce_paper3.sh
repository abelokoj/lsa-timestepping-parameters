#!/usr/bin/env bash
# Reproduces the results of Paper 3 (amplitude, time step, sampling and field weighting in time-stepping stability analysis) at three levels.
#   ./reproduce_paper3.sh 1          build the analysis code, run its tests and print the reference growth rates
#   ./reproduce_paper3.sh 2          recompute the tables of the paper from the stored data (no solver)
#   ./reproduce_paper3.sh 3          submit the simulations (this folder must be inside amrex_implicit_interfaces/NS_fluids_lib)
#   ./reproduce_paper3.sh 3-analyse  analyse the new runs once every job has ended
# Run it from this folder. Level 2 writes its output to figures/paper3_tables.txt.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
SNAP=results/followup_2026-10-22
S12=paper3/data/2026-11-24_paper3       # sets 1 and 2: low Ohnesorge numbers and fields stored every four steps
S34=paper3/data/2026-11-25_sets34       # sets 3 and 4: start from the leading mode, and the restarting variant with its fixes
OH002="--mu-l 0.002 --mu-g 2.6e-5"

level1 () {
  make -j8
  make test | tail -12
  # the unbounded two-fluid relation, the reference of Tables 1 to 3, at Oh = 0.02 and at Oh = 0.002
  ./bin/run_lsa theory --kr0 0.7,0.9,0.95 | tail -6
  ./bin/run_lsa theory --kr0 0.7 $OH002 | tail -3
}

level2 () {
  [[ -x bin/run_lsa ]] || level1
  "$PY" -c "import numpy, scipy, matplotlib" 2>/dev/null || { echo "install the Python packages first: $PY -m pip install -r plotting/requirements.txt"; exit 1; }
  mkdir -p figures
  "$PY" plotting/thesis/print_thesis_tables.py > figures/dissertation_tables.txt
  {
    echo "== Table 1 (Section 4.1): the interface rate at three amplitudes on three grids"
    sed -n '/^Table 6.2/,/^ *32 /p' figures/dissertation_tables.txt
    echo "-- last column: the smallest amplitude refitted over the snapshots of the window of the largest, both ends included"
    early_window amp08_025 amp08_100 8
    early_window sw070_16 lin100_16 16
    early_window amp32_025 amp32_100 32
    echo
    echo "== Section 4.2: windows of the runs started from rest and from the leading mode, refitted from their records"
    for c in "oh002_16 0.7 $OH002" "lm070oh002_16 0.7 $OH002" "sw090_16 0.9" "lm090_16 0.9" "sw095e1_16 0.95" "lm095_16 0.95"; do
      set -- $c; l=$1; k=$2; shift 2
      printf '%-15s ' "$l"; ./bin/run_lsa record --record $S34/records/$l.dat --kr0 $k "$@" | awk '/window t/{w=$0} /two-fluid visc/{print w " | " $0}' | tr -s ' '
      printf '%-15s ' "  stored"; awk '/window t/{w=$0} /two-fluid visc/{print w " | " $0}' $S34/window/$l.txt | tr -s ' '
    done
    echo "-- Table 2: local rate of the interface amplitude in bands of a/R0 (stored output: $S34/band_rates.txt)"
    "$PY" paper3/code/band_rates.py $S34/records oh002_16=0.341648129 lm070oh002_16=0.341648129 sw090_16=0.242064407 \
        lm090_16=0.242064407 sw095e1_16=0.174895587 lm095_16=0.174895587 > figures/band_rates.txt
    cat figures/band_rates.txt
    diff -q <(tr -d '\r' < $S34/band_rates.txt) figures/band_rates.txt > /dev/null && echo "(identical to the stored output)" || echo "(DIFFERS from the stored output)"
    echo "-- decompositions over 0.01 < a <= 0.05 R0 (stored result):"
    sed -n '/^| lm090_16 | 5.73/,/^| lm095_16 | 7.89/p' $S34/README.md
    echo
    echo "== Section 4.3 and Table 3: the sampling interval and the Ohnesorge number"
    scripts/refit_records.sh $S12 | grep -E "^(case|oh00|str04)"
    echo "-- decomposition rates (median over ranks) and IQR, against the two-fluid relation (stored output of the analysis job)"
    awk -F, 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} $1~/^(sw070_16|str04_16|oh005_16|oh002_16|str04oh002_16)$/ {
        s=$c["visc_two_fluid"]; printf "%-15s velocity %+6.2f%% (IQR %.2f%%)  level set %+6.2f%% (IQR %.2f%%)\n", $1,
        100*($c["dmd_rankscan_median"]/s-1), 100*$c["dmd_rankscan_iqr"]/$c["dmd_rankscan_median"],
        100*($c["dmd_levelset_median"]/s-1), 100*$c["dmd_levelset_iqr"]/$c["dmd_levelset_median"]}' $S12/analysis/summary.csv
    echo "-- velocity and level-set rates at fixed rank over the strides, each over the whole window (stored studies)"
    "$PY" - $S12/analysis <<'PY'
import csv, glob, os, sys
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*", "studies", "stride*_*.csv"))):
    rows = list(csv.DictReader(open(f)))
    for rank in ("4", "8"):
        best = {}   # the whole window of each stride: the row with the most snapshots
        for r in rows:
            if r["rank_requested"].strip() == rank:
                s = int(r["stride"])
                if s not in best or int(r["n_snap"]) > int(best[s]["n_snap"]):
                    best[s] = r
        if best:
            g = [float(best[s]["growth"]) for s in sorted(best)]
            moved = max(abs(x / g[0] - 1) for x in g)   # largest change from the rate at stride 1
            print(f"{os.path.relpath(f, sys.argv[1]):40} rank {rank}: strides {sorted(best)}: rates {' '.join(f'{x:.6f}' for x in g)};"
                  f" spread {100 * (max(g) - min(g)) / min(g):.2f}%, largest change from stride 1 {100 * moved:.2f}%")
PY
    echo
    echo "== Table 4 (Section 4.4): joint decomposition against the weight of the level set"
    (cd paper3/code && "$PY" weighting_table.py)
    echo "-- the same test on the runs of sets 1 to 3 (stored output)"
    cat $S34/weighting_lowoh_2026-11-25.txt
    echo
    echo "== Table 5 (Section 4.5): sensitivity budget at k R0 = 0.7, 16 cells"
    sed -n '/^Table 5.3/,/^Velocity:/p' figures/dissertation_tables.txt | grep -E "x |0\.7 "
    sed -n '/^Table 6.1/,/^Two-fluid, wall/p' figures/dissertation_tables.txt | grep -E "x |0\.7 "
    awk -F, 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next} {r[$1]=$c["rate_interface"]} END{
        printf "amplitude: eps 2.5e-3 against 5e-3      %.3f%%\n", 100*(r["lin050_16"]/r["sw070_16"]-1);
        printf "grid: 16 against 32 cells (eps 5e-3)    %.3f%%\n", 100*(r["meas32"]/r["meas16"]-1);
        printf "time step halved                        %.3f%%\n", 100*(r["dts125_16"]/r["sw070_16"]-1);
        printf "wall at 6 R0 against 4 R0               %.3f%%\n", 100*(r["dom06_16"]/r["sw070_16"]-1)}' <(tr -d ' ' < $SNAP/analysis/summary.csv)
    echo
    echo "== Section 4.6: the restarting variant (stored output of the runs)"
    for d in $S34/krylov_test/*; do
      printf '%-15s member rate %s | level set: %s | velocity: %s\n' "$(basename $d)" "$(tail -1 $d/member_rates.txt | awk '{print $3}')" \
        "$(grep -m1 'median rate' $d/dmd_ls.log | tr -s ' ')" "$(grep -m1 'median rate' $d/dmd_vel.log | tr -s ' ')"
    done
    echo "-- the restarting variant without the fixes (dissertation):"; sed -n '1,30p' results/restarting_members_2026-10-23/README.md
  } 2>&1 | tee figures/paper3_tables.txt
  echo "Output: figures/paper3_tables.txt (see README.md for the table of the paper that each part reproduces)."
}

early_window () {   # small big cells: fit the record of the small amplitude over the window of the large one
  local lo hi
  read -r lo hi <<< "$(awk -F, -v l="$2" '{gsub(/ /,"")} $1==l {print $3 - 1e-4, $4 + 1e-4}' $SNAP/analysis/summary.csv)"
  printf '%2s cells: ' "$3"
  ./bin/run_lsa record --record $SNAP/records/$1.dat --kr0 0.7 --t-min "$lo" --t-max "$hi" | awk '/window t/{w=$0} /two-fluid visc/{print w " | " $0}' | tr -s ' '
}

level3 () {
  [[ -x ../amr2d.gnu.MPI.ex ]] || { echo "build the solver first: ../amr2d.gnu.MPI.ex not found (see README.md, Level 3)"; exit 1; }
  local acc="${ACCOUNT:?set ACCOUNT to your Slurm account}"
  make -j8 && make test | tail -1
  # the dissertation runs that Tables 1, 2, 4 and 5 use, and sets 1 and 2 (unmodified solver)
  ACCOUNT="$acc" scripts/submit_campaign.sh campaign_16_wide.txt campaign_16.txt rerun_meas.txt campaign_next_16.txt \
      campaign_amp08.txt campaign_amp32.txt campaign_check.txt campaign_p3_16.txt campaign_p3_16_str04.txt
  # sets 3 and 4: the restarting driver with the fixes, on the patched solver built in its own tree (README.md, Level 3)
  local exe="${FIX_SOLVER_EXE:?set FIX_SOLVER_EXE to the solver built with solver_patches/lsa_fix_restart.patch}"
  for fix in 4 5 6 7; do
    SOLVER_EXE="$exe" FIX=$fix sbatch -A "$acc" --export=ALL krylov_test/krylov_fix.sh kfx${fix}_070_16 0.7 16 4.0 4.3 160 40 0.25
  done
  SOLVER_EXE="$exe" FIX=7 CONTINUE=1 sbatch -A "$acc" --export=ALL krylov_test/krylov_fix.sh lm070oh002_16 0.7 16 4.0 4.3 160 60 0.25 $OH002
  SOLVER_EXE="$exe" FIX=7 CONTINUE=1 sbatch -A "$acc" --export=ALL krylov_test/krylov_fix.sh lm090_16 0.9 16 4.0 4.3 160 60 0.25
  SOLVER_EXE="$exe" FIX=7 CONTINUE=1 sbatch -A "$acc" --export=ALL krylov_test/krylov_fix.sh lm095_16 0.95 16 4.0 5.2 200 60 0.25
  echo "When every job has ended (squeue --me), run: ./reproduce_paper3.sh 3-analyse"
}

level3_analyse () {
  local acc="${ACCOUNT:?set ACCOUNT to your Slurm account}"
  sbatch -A "$acc" -n 1 -c 4 --mem=16G -t 04:00:00 -o runs/array/studies-%j.out scripts/run_studies.sh all
  sbatch -A "$acc" -n 1 -c 4 --mem=16G -t 04:00:00 -o runs/array/stride52-%j.out scripts/stride52.sh str04oh002_16 str04_16
  "$PY" plotting/weighting_test.py runs/analysis sw040_16 sw050_16 sw060_16 sw070_16 sw080_16 meas16 meas32 conf15_32 visc100r8_16 > runs/analysis/weighting_test_output.txt
  "$PY" plotting/weighting_test.py runs/analysis oh002_16 oh005_16 str04oh002_16 str04_16 lm070oh002_16 > runs/analysis/weighting_lowoh.txt
  "$PY" plotting/thesis/compare_summaries.py runs/analysis/summary.csv
  scripts/refit_records.sh runs
  echo "When the batch jobs have ended, Level 2 can be pointed at the new runs (README.md, Level 3)."
}

case "${1:-}" in
  1) level1 ;;
  2) level2 ;;
  3) level3 ;;
  3-analyse) level3_analyse ;;
  *) sed -n '2,8p' "$0"; exit 1 ;;
esac
