#!/bin/bash
# stride52.sh -- Paper 3: strides 1, 13, 31 and 32 of a 4-step run over the interface window, as in run_studies.sh stride.
# Every 13th snapshot of a 4-step run is the 52-step interval of sw070_16; every 31st and 32nd (124, 128 steps) bracket 126.
# Reads only the exports (runs/analysis/<label>/export, export_ls) and the record window, so it works after the plotfiles are deleted.
# From $LSA:
#   sbatch -A genacc_q -n 1 -c 4 --mem=16G -t 04:00:00 -o runs/array/stride52-%j.out scripts/stride52.sh str04oh002_16 str04_16
# Output: runs/analysis/<label>/studies/stride52_{vel,ls}.csv and .log
cd "${SLURM_SUBMIT_DIR:-.}"
[ -x bin/run_lsa ] || { echo "bin/run_lsa missing: submit from \$LSA" >&2; exit 1; }
# Fluid properties from the deck (as analyse_campaign.sh), so that err_vs_exact_pct uses each run's own Ohnesorge number.
props_of () {
  local deck="runs/$1/inputs"
  [[ -f "$deck" ]] || return 0
  awk '/^ns\.denconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],d," "); rg=d[2]}
       /^ns\.viscconst[ \t]*=/{sub(/#.*/,""); split($0,a,"="); split(a[2],v," "); ml=v[1]; mg=v[2]}
       END{if (ml!="") printf "--mu-l %s ", ml; if (rg!="") printf "--rho-g %s ", rg; if (mg!="") printf "--mu-g %s", mg}' "$deck"
}
for label in "$@"; do
  rec=runs/analysis/records/$label.csv; out=runs/analysis/$label/studies; mkdir -p "$out"
  [ -f "$rec" ] || { echo "$label: no record $rec"; continue; }
  for kind in vel ls; do
    if [ $kind = vel ]; then exp=runs/analysis/$label/export; opt="--weight volume"
    else exp=runs/analysis/$label/export_ls; opt="--weight volume --band 2"; fi
    [ -f "$exp/perturbation.bin" ] || { echo "$label/$kind: no export"; continue; }
    read lo hi <<< "$(awk -F, 'NR==2{print $3, $4}' "$rec")"
    read sk mx <<< "$(awk -v lo="$lo" -v hi="$hi" '{t=$1+0; if (t < lo-1e-9) s++; else if (t <= hi+1e-9) m++} END{print s+0, m+0}' "$exp/times.txt")"
    ./bin/run_lsa converge --from-export "$exp" --kr0 0.7 $(props_of "$label") $opt --max-snapshots $((sk + mx)) \
      --starts "$sk" --strides 1,13,31,32 --ranks 4,8,16,0 --csv "$out/stride52_$kind.csv" > "$out/stride52_$kind.log" 2>&1
    echo "$label/$kind: window skip $sk, $mx snapshots -> $out/stride52_$kind.csv (exit $?)"
  done
done
