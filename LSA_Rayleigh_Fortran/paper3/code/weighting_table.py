"""
weighting_table.py
==================
Rebuilds the field-weighting table (dissertation Table 6.3) from the raw output of plotting/weighting_test.py
(data/weighting_test_output_2026-10-22.txt): for each case, the rate of the joint velocity and level-set decomposition at
each weight lambda of the level-set block, as a departure (%) from the confined two-fluid relation for the run's own wall and
fluids (code/confined_dispersion.py). Writes ../data/weighting_table.csv and prints it.

    python weighting_table.py
"""
import csv, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from confined_dispersion import growth_rate

OH, RHO2 = 0.02, 1.225e-3
# label: (description, x, Rw, mu2)
CASES = {'sw040_16': ('x = 0.4', 0.4, 4.0, 2.6e-4), 'sw050_16': ('x = 0.5', 0.5, 4.0, 2.6e-4), 'sw060_16': ('x = 0.6', 0.6, 4.0, 2.6e-4),
         'sw070_16': ('x = 0.7', 0.7, 4.0, 2.6e-4), 'meas32': ('x = 0.7, 32 cells', 0.7, 4.0, 2.6e-4),
         'conf15_32': ('x = 0.7, R_w = 1.5, 32 cells', 0.7, 1.5, 2.6e-4), 'visc100r8_16': ('x = 0.7, mu2/mu1 = 1, R_w = 8', 0.7, 8.0, 0.02),
         'sw080_16': ('x = 0.8 (rank scatter)', 0.8, 4.0, 2.6e-4), 'meas16': ('x = 0.7, 16 cells (rank scatter)', 0.7, 4.0, 2.6e-4)}
COLS = ['velocity only', 'lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'level set only (band)']
HERE = os.path.dirname(os.path.abspath(__file__))
text = open(os.path.join(HERE, '..', 'data', 'weighting_test_output_2026-10-22.txt'), encoding='utf-8').read()
rows = []
for block in re.split(r'\n(?===)', text):
    m = re.match(r'== (\S+):.*?at lambda = 1: ([\d.]+)', block)
    if not m or m.group(1) not in CASES:
        continue
    lab, ratio = m.group(1), float(m.group(2))
    desc, x, Rw, mu2 = CASES[lab]
    ref = growth_rate(x, OH, RHO2, mu2, Rw)
    rates = dict(re.findall(r'^\s+(velocity only|lambda = [\d.]+|level set only \(band\))\s+rate ([\d.]+)', block, re.M))
    rows.append([lab, desc, f'{ratio:.2f}'] + [f'{100 * (float(rates[c]) / ref - 1):+.2f}' for c in COLS])
assert len(rows) == len(CASES), f'found {len(rows)} of {len(CASES)} cases'
with open(os.path.join(HERE, '..', 'data', 'weighting_table.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['run', 'case', 'norm_ratio'] + COLS); w.writerows(rows)
for r in rows:
    print(' | '.join(r))
