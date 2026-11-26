# Sets 3 and 4 of Paper 3 (runs of 24 and 25 Nov 2026)

Solver: fork commit dddb827c (= upstream 0b8df42, the dissertation solver) with `runs/lsa_fix_restart.patch`; job script `runs/krylov_fix.sh`. Collected on the cluster as `LSA_Rayleigh/analysis and results/2026-11-25_1154_paper3_sets34` (`MANIFEST_cluster.txt`). Reference: the unbounded two-fluid viscous relation (`visc_two_fluid` of summary.csv, the reference of Paper 3's tables): 0.341648129 (x = 0.7, Oh = 0.002), 0.3286887 (x = 0.7, Oh = 0.02), 0.242064407 (x = 0.9), 0.174895587 (x = 0.95).

## Set 4: restarting driver with the fixes, x = 0.7, Oh = 0.02, 40 members of tau = 0.25 t_sigma, a = 2.5e-3 R0

`krylov_test/kfx{4,5,6,7}_070_16`. FIX bits: 1 level set kept outside the band, 2 member pressure carried, 4 one scalar factor per member. `lsa1_key_lines.txt` shows that each bit was active (pressure component nonzero only with bit 2, 0.0044915; level-set maximum over the domain 0.1308 with bit 1, against 0.0025 without), yet the four chains are identical to the last digit:

| Estimate (members 21-40) | rate | vs two-fluid |
|---|---|---|
| member growth, -ln(c_k)/tau (last members) | 0.32009 | -2.62% |
| level-set DMD, median ranks 2-12 (IQR 1.1e-4) | 0.320221 | -2.58% |
| velocity DMD, median ranks 4-16 (IQR 5.9e-5) | 0.321490 | -2.19% |

The same chains of set 3 (60 members, all fixes), members 21-60, one-fluid viscous reference from run_lsa:

| Case | one-fluid exact | level-set DMD | velocity DMD | rate deficit (level set, velocity) |
|---|---|---|---|---|
| x = 0.7, Oh = 0.02 (kfx) | 0.329367 | 0.320221 (-2.78%) | 0.321490 (-2.39%) | 0.0091, 0.0079 |
| x = 0.9 (lm090_16) | 0.242730 | 0.233428 (-3.83%) | 0.233613 (-3.76%) | 0.0093, 0.0091 |
| x = 0.95 (lm095_16) | 0.175494 | 0.165663 (-5.60%) | 0.167062 (-4.81%) | 0.0098, 0.0084 |
| x = 0.7, Oh = 0.002 (lm070oh002_16) | 0.341863 | 0.333856 (-2.34%) | 0.334078 (-2.28%) | 0.0080, 0.0078 |

## Set 3: continuous twin started from the last member (leading mode), a0 = 2.6e-3 R0

Pipeline (window rule of the campaign, summary.csv), departures from the two-fluid relation:

| Pair | window t | samples, e-folds | interface | velocity DMD (IQR) | level-set DMD (IQR) |
|---|---|---|---|---|---|
| oh002_16 (campaign) | 8.85-12.50 | 48, 1.25 | -0.11% | +1.45% (0.47%) | -0.20% (0.08%) |
| lm070oh002_16 | 2.02-10.79 | 114, 3.00 | +0.08% | -0.21% (0.34%) | -0.17% (0.10%) |
| sw090_16 (campaign) | 16.43-17.64 | 12, 0.30 | +1.41% | +8.2% | -1.4% |
| lm090_16 | 14.00-15.11 | 11, 0.27 | +1.47% | +22% | -1.0% |
| sw095e1_16 (campaign, eps 1e-3) | 28.60-29.34 | 6, 0.14 | +3.84% | +16% | -1.9% |
| lm095_16 | 19.81-20.70 | 7, 0.16 | +3.64% | -4.1% | -3.0% |

Local interface rate in amplitude bands: `band_rates.txt` (`code/band_rates.py`). Decomposition over the band 0.01 < a <= 0.05 R0 (run on the cluster 25 Nov 2026, `run_lsa analyse --from-export`, same ranks as the campaign):

| Pair | band t | snapshots | velocity DMD (IQR) | level-set DMD (IQR) |
|---|---|---|---|---|
| lm090_16 | 5.73-12.24 | 60 | 0.242246, +0.08% (1.7e-3) | 0.242779, +0.30% (1.3e-4) |
| lm095_16 | 7.89-16.83 | 61 | 0.176118, +0.70% (3.1e-3) | 0.175664, +0.44% (8.5e-4) |

The pairs lm*_16 have a plotfile at the restart step (t = tau) whose spacing to the next differs from the rest; it was set aside (`nddataPLT...161.restart`, `...175.restart`) before the export. The amplitude records include it.
