#!/bin/bash
# krylov_fix.sh -- Paper 3, sets 3 and 4: the solver's restarting Krylov driver with the fixes of the dissertation (lsa_fix_restart.patch).
#
# Solver: the dissertation solver, upstream commit 0b8df42 (= fork commit dddb827c, same source), with lsa_fix_restart.patch, so that the
#   only difference from the dissertation runs (campaign and the restarting trial kry070_16b) is the fixes. 0b8df42 has no LSA_activate = 2.
# Set 4 (accuracy of the fixed restarting variant at x = 0.7): the members (LSA_activate = 1), decomposed as in krylov_test/krylov_run.sh, with
#   ns.LSA_fix_restart = FIX (bits: 1 keep the level set outside the band, 2 carry the pressure, 4 one scalar factor per member; 4 is required by 1 and 2)
#   ns.LSA_levelset_scale = LS_SCALE (level-set amplitude of every member in finest cells; 0.04 = 2.5e-3 R0 at 16 cells, the campaign amplitude).
#   With bit 4 every member is scaled by one factor c_k, printed as "LSA_scalar_factor k c_k". The member plotfiles are exported (run_lsa export
#   --krylov), the factors are divided out (rescale_members.py, member_scales.txt), and the level set in the band and the velocity are decomposed
#   from members SKIP+1..N (default SKIP = 20, members 21-40 as in the dissertation). The member growth 1/c_k also gives a power-iteration rate,
#   s_k = -ln(c_k)/tau, written to member_rates.txt.
# Set 3 (start from the leading mode): CONTINUE=1 runs, after the chain, a continuous twin from its last member: the perturbed run restarts from
#   chk<M>LSA<N> (base + the rescaled N-th member), the base twin from chk<M>LSA00000 (the unperturbed member 0), both at step M, with LSA off,
#   into $LSA/runs/<label> and <label>.base, and writes runs/records/<label>.dat and <label>.base.dat, so analyse_campaign.sh treats the pair
#   like any campaign pair. The deck is the campaign deck of the given row, so the window can be compared with the corresponding sw run.
#
# Build, once, in its own worktree of the campaign solver (the campaign executable is never touched):
#   git -C $AII worktree add ~/aii_fix $(git -C $AII rev-parse HEAD)      # HEAD must be 0b8df42 or dddb827c
#   cd ~/aii_fix && git apply /path/to/lsa-timestepping-parameters/solver_patches/lsa_fix_restart.patch && cd NS_fluids_lib && make -j16
#
# Usage (from $LSA):
#   SOLVER_EXE=$(ls -t ~/aii_fix/NS_fluids_lib/amr2d*.ex | head -1) FIX=7 [LS_SCALE=0.04] [SKIP=20] [CONTINUE=0] \
#     sbatch --export=ALL krylov_test/krylov_fix.sh LABEL KR0 CELLS RMAX NPER SNAPS NMEM TSAMPLE [make_inputs options]
#   NPER and SNAPS are the columns of the campaign row (the continuation deck); NMEM members of TSAMPLE (in t_sigma) each.
#SBATCH --job-name=lsa_krylov_fix
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --account=genacc_q
#SBATCH --output=krylov_test/slurm-%x-%j.out
#SBATCH --error=krylov_test/slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
TEST="$PWD/krylov_test"; RUNROOT="$PWD/runs"
[[ $# -ge 8 ]] || { sed -n '2,28p' "$0"; exit 1; }
LABEL="$1"; KR0="$2"; CELLS="$3"; RMAX="$4"; NPER="$5"; SNAPS="$6"; NMEM="$7"; TSAMPLE="$8"; shift 8
FIX="${FIX:?set FIX (4, 5, 6 or 7)}"; LS_SCALE="${LS_SCALE:-0.04}"; SKIP="${SKIP:-20}"; CONTINUE="${CONTINUE:-0}"
EXE="${SOLVER_EXE:?set SOLVER_EXE to the patched executable}"
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
NPROCS="${SLURM_NTASKS:-16}"
SOLVER_OPTS="${SOLVER_OPTS:-mac.mac_abs_tol=1.0e-10 mac.visc_abs_tol=1.0e-10}"
unset SLURM_STEP_ID SLURM_STEPID SLURM_STEP_NUM_TASKS SLURM_STEP_NODELIST || true
set +u
if type module >/dev/null 2>&1; then module purge; module load gnu/13; module load openmpi; fi
set -u
[[ -x "$EXE" ]] || { echo "error: solver $EXE not found" >&2; exit 1; }
grep -q "LSA_fix_restart" "$EXE" || { echo "error: $EXE was not built with lsa_fix_restart.patch" >&2; exit 1; }

DIR="$TEST/runs/$LABEL"; OUT="$TEST/analysis/$LABEL"
mkdir -p "$DIR/decks" "$OUT"
# The campaign deck of the row: same generator call as run_all.sh, the unperturbed (.base) deck, since the driver supplies the perturbation.
./bin/make_inputs single --kr0 "$KR0" --eps 2.5e-3 --cells-per-r0 "$CELLS" --rmax "$RMAX" \
    --n-periods "$NPER" --snapshots "$SNAPS" "$@" --outdir "$DIR/decks" > "$DIR/plan.txt"
base=$(grep -E '/inputs\.k.*\.base$' "$DIR/plan.txt" | awk '{print $NF}' | tail -1)
cp "$base" "$DIR/inputs"
dt=$(awk -F= '/^ns.fixed_dt/{print $2+0}' "$DIR/inputs")
M=$(awk -v t="$TSAMPLE" -v d="$dt" 'BEGIN{m=int(t/d+0.5); if (m<1) m=1; print m}')
TAU=$(awk -v m="$M" -v d="$dt" 'BEGIN{printf "%.10g", m*d}')
LAUNCH=(srun --cpu-bind=none -n "$NPROCS")
COMMON=($SOLVER_OPTS amr.max_grid_size=16 amr.LSA_nsteps_krylov_subspace_method=$NMEM amr.LSA_plot_index=0
        ns.LSA_fix_restart=$FIX ns.LSA_levelset_scale=$LS_SCALE)
{
  echo "--- $(date '+%Y-%m-%d %H:%M:%S') krylov_fix.sh"
  echo "label = $LABEL  kr0 = $KR0  cells = $CELLS  rmax = $RMAX  nper = $NPER  snaps = $SNAPS  extra = $*"
  echo "members = $NMEM  M = $M  tau = $TAU  fix = $FIX  ls_scale = $LS_SCALE  skip = $SKIP  continue = $CONTINUE"
  echo "solver = $EXE  md5 = $(md5sum "$EXE" | cut -c1-32)  commit = $(git -C "$(dirname "$EXE")" rev-parse --short HEAD 2>/dev/null || echo ?)"
} >> "$DIR/manifest.txt"
echo "case $LABEL: M = $M steps per member, tau = $TAU, $NMEM members, fix = $FIX"
cd "$DIR"

# 1. base state (as krylov_upstream.sh)
if [[ ! -f chk00000LSA00000/Header ]]; then
  "${LAUNCH[@]}" "$EXE" inputs "${COMMON[@]}" amr.LSA_activate=0 max_step=1 amr.check_int=1 amr.plot_int=1 \
      < /dev/null > base.out 2> base.err
fi
[[ -f chk00000LSA00000/Header ]] || { echo "error: chk00000LSA00000 not written; see $DIR/base.err" >&2; exit 1; }

# 2. the members (LSA_activate=1)
last=$(printf "chk%05dLSA%05d" "$M" "$NMEM"); zero=$(printf "chk%05dLSA00000" "$M")
if [[ ! -f "$last/Header" ]]; then
  echo "[$(date +%T)] members (LSA_activate=1)"
  "${LAUNCH[@]}" "$EXE" inputs "${COMMON[@]}" amr.restart=chk00000LSA00000 amr.LSA_activate=1 \
      amr.LSA_initial_levelSteps=0 max_step=$M stop_time=1.0e6 amr.plot_int=$M amr.check_int=$M \
      < /dev/null > lsa1.out 2> lsa1.err || echo "WARNING: LSA_activate=1 ended with an error; see $DIR/lsa1.err"
fi
[[ -f "$last/Header" ]] || { echo "error: $last not written; see $DIR/lsa1.err" >&2; exit 1; }

# Factors: the member plotfile k is d_k = A w_{k-1}, with w_k = c_k d_k re-added, so A^k w_0 = d_k / (c_1 ... c_{k-1}); growth over one member is 1/c_k.
cp lsa1.out "$OUT/"
awk '$1=="LSA_scalar_factor" && $2>=1 {c[$2]=$3; if ($2>n) n=$2} END{p=1; for (k=1; k<=n; k++) {printf "%.17e\n", 1/p; p*=c[k]}}' lsa1.out > "$OUT/member_scales.txt"
awk -v tau="$TAU" '$1=="LSA_scalar_factor" && $2>=1 {printf "%d %.17e %.8f\n", $2, $3, -log($3)/tau}' lsa1.out > "$OUT/member_rates.txt"
nf=$(wc -l < "$OUT/member_scales.txt")
[[ "$nf" -eq "$NMEM" ]] || echo "WARNING: $nf scalar factors for $NMEM members (bit 4 of FIX off, or a member failed)"
echo "member rates s_k = -ln(c_k)/tau, last five (k, c_k, s_k):"; tail -5 "$OUT/member_rates.txt"
cd - > /dev/null

# 3. decomposition of the members (as krylov_run.sh), from the export with the factors divided out
# Fluid properties from the deck, as analyse_campaign.sh (props_of), so that the reference rate is that of the run's own Ohnesorge number.
props=$(awk '/^ns\.denconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],d," "); rg=d[2]}
             /^ns\.viscconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],v," "); ml=v[1]; mg=v[2]}
             END{if (ml!="") printf "--mu-l %s ", ml; if (rg!="") printf "--rho-g %s ", rg; if (mg!="") printf "--mu-g %s", mg}' "$DIR/inputs")
for kind in ls vel; do
  if [[ $kind == ls ]]; then flds=L0101; opt=(--band 2); ranks=2:12; else flds=x_velocity,y_velocity; opt=(); ranks=4:16; fi
  if [[ ! -s "$OUT/export_$kind/perturbation.bin" ]]; then
    ./bin/run_lsa export --krylov "$DIR" --fields "$flds" --export "$OUT/export_raw_$kind" > "$OUT/export_$kind.log" 2>&1 \
      && python3 "$TEST/rescale_members.py" "$OUT/export_raw_$kind" "$OUT/member_scales.txt" "$OUT/export_$kind" \
      || { echo "WARNING: $kind export or rescaling failed; see $OUT/export_$kind.log"; continue; }
  fi
  ./bin/run_lsa analyse --from-export "$OUT/export_$kind" --kr0 "$KR0" $props --weight volume ${opt[@]+"${opt[@]}"} --skip "$SKIP" \
      --rank-scan "$ranks" --archive "$OUT/dmd_$kind" > "$OUT/dmd_$kind.log" 2>&1 || true
  echo "$kind (members $((SKIP+1))-$NMEM): $(grep -E 'median rate|median vs' "$OUT/dmd_$kind.log" | tr '\n' ' ')"
done
cd "$DIR"

# 4. continuous twin from the last member (Set 3)
if [[ "$CONTINUE" == 1 ]]; then
  [[ -f "$zero/Header" ]] || { echo "error: $zero (unperturbed member 0 at step $M) missing" >&2; exit 1; }
  mkdir -p "$RUNROOT/records"
  for pair in "$LABEL:$last" "$LABEL.base:$zero"; do
    lab="${pair%%:*}"; start="${pair#*:}"; run="$RUNROOT/$lab"
    rows=$(grep -cve '^[[:space:]]*#' -e '^[[:space:]]*$' "$RUNROOT/records/$lab.dat" 2>/dev/null || true)
    (( ${rows:-0} >= 2 )) && { echo "  $lab: $rows record rows already, skipped"; continue; }
    mkdir -p "$run"; cp "$DIR/inputs" "$run/inputs"
    [[ -d "$run/$start" ]] || cp -r "$start" "$run/"
    # Line 6 of a checkpoint Header is LSA_current_step (AmrCore::checkPoint), and main.cpp aborts a restart with LSA_activate = 0
    # unless it is 0 ("initial_LSA_current_step invalid"). The member checkpoint has N there; the copy (not the chain's original) is set to 0.
    awk 'NR==1 && $0!="CheckPointVersion_1.0" {exit 1} NR==6 && $0!~/^[0-9]+$/ {exit 1}' "$run/$start/Header" \
      || { echo "error: unexpected Header layout in $run/$start" >&2; exit 1; }
    sed -i '6s/^[0-9][0-9]*$/0/' "$run/$start/Header"
    rm -f "$run"/Backtrace.*
    # LSA storage stays on, or the member checkpoint cannot be read; LSA_activate=0 adds no perturbation. Resume from the newest checkpoint.
    rs=$(cd "$run" && for c in chk[0-9]*LSA[0-9]*; do [[ -f "$c/Header" ]] && echo "$c"; done | sort -V | tail -1)
    echo "[$(date +%T)] $lab: continuation from $rs"
    ( cd "$run" && "${LAUNCH[@]}" "$EXE" inputs $SOLVER_OPTS amr.max_grid_size=16 amr.LSA_nsteps_krylov_subspace_method=$NMEM \
        amr.LSA_activate=0 amr.restart="$rs" amr.check_int="${CHECK_INT:-520}" < /dev/null >> solver.out 2>> solver.err ) \
      || echo "WARNING: $lab ended with an error; see $run/solver.err"
    # With LSA storage on the solver names plotfiles nddataPLT<step>LSA00000; the extractor reads nddataPLT<step> only.
    ( cd "$run" && for d in nddataPLT*LSA[0-9]*; do [[ -e "$d" ]] && mv "$d" "${d%LSA*}"; done; true )
    lambda=$(awk '/^ns.yblob/{print $3; exit}' "$run/inputs")
    ( cd "$RUNROOT/.." && ./bin/dump_amplitude "$run" "$lambda" "$RUNROOT/records/$lab.dat" > /dev/null 2>&1 ) \
      || echo "WARNING: dump_amplitude failed for $lab"
    echo "  $lab: $(grep -cve '^[[:space:]]*#' "$RUNROOT/records/$lab.dat" 2>/dev/null || echo 0) record rows"
  done
  echo "continuation done; analyse with the post job (scripts/_post_task.sh) and compare with the corresponding sw run"
fi
echo "[$(date +%T)] done; results in $OUT"
