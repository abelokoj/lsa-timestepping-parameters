"""Local growth rate of an interface-amplitude record in amplitude bands (Paper 3, sets 3 and 4).

For each record (runs/records/<label>.dat: time, amplitude, ...), fit ln a against t by least squares within each band
a0 < a <= a1 (a <= 0.1 R0 only), and from the first sample up to a given amplitude, and print the departure from the
two-fluid relation. Usage: python band_rates.py <records dir> label=two_fluid_rate [label=rate ...]
"""
import sys
import numpy as np

BANDS = ((0.0, 0.005), (0.005, 0.01), (0.01, 0.02), (0.02, 0.05), (0.05, 0.1))


def rate(t, a):
    return np.polyfit(t, np.log(a), 1)[0]


def main(folder, cases):
    for lab, s0 in cases:
        d = np.loadtxt(f"{folder}/{lab}.dat")
        t, a = d[:, 0], d[:, 1]
        m = a <= 0.1
        t, a = t[m], a[m]
        print(f"== {lab}: t {t[0]:.2f}-{t[-1]:.2f}, a {a[0]:.4f}-{a[-1]:.4f}, two-fluid {s0:.6f}")
        for lo, hi in BANDS:
            k = (a > lo) & (a <= hi)
            if k.sum() >= 4:
                s = rate(t[k], a[k])
                print(f"   a in ({lo:.3f},{hi:.3f}]: {k.sum():3d} samples, t {t[k][0]:6.2f}-{t[k][-1]:6.2f}, rate {s:.6f} ({100*(s/s0-1):+.2f}%)")
        for amax in (0.01, 0.02, 0.05, 0.1):
            k = a <= amax
            s = rate(t[k], a[k])
            print(f"   first sample to a = {amax}: {k.sum()} samples, {s*(t[k][-1]-t[k][0]):.2f} e-folds, rate {s:.6f} ({100*(s/s0-1):+.2f}%)")


def self_check():
    t = np.linspace(0, 10, 101)
    assert abs(rate(t, 1e-3 * np.exp(0.3 * t)) - 0.3) < 1e-12


if __name__ == "__main__":
    self_check()
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], [(c.split("=")[0], float(c.split("=")[1])) for c in sys.argv[2:]])
