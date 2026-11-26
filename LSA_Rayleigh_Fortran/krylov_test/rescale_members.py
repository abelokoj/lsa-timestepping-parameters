"""Divide the scalar factors out of an export of the restarting driver's members (run_lsa export --krylov), so that the columns are a Krylov sequence again.

Column k (member k, k = 1..N) of the export is d_k = A w_{k-1}, where w_{k-1} = c_{k-1} d_{k-1} is the member re-added to the base state and c_j is the
factor printed as "LSA_scalar_factor j c_j" (lsa_fix_restart.patch, bit 4). Then A^k w_0 = d_k / (c_1 ... c_{k-1}), so column k is multiplied by
f_k = 1 / (c_1 ... c_{k-1}), f_1 = 1; krylov_fix.sh writes the f_k to member_scales.txt.

Usage: python3 rescale_members.py EXPORT_DIR member_scales.txt OUT_DIR
"""
import os, shutil, sys
import numpy as np


def rescale(src, scales, out):
    meta = {k.strip(): v.strip() for k, v in (l.split("=", 1) for l in open(f"{src}/perturbation.txt") if "=" in l)}
    nsp, nsn = int(meta["n_space"]), int(meta["n_snap"])
    f = np.loadtxt(scales, ndmin=1)
    assert f.size == nsn, f"{f.size} factors for {nsn} members"
    X = np.fromfile(f"{src}/perturbation.bin", "<f8").reshape((nsp, nsn), order="F")
    os.makedirs(out, exist_ok=True)
    (X * f[None, :]).reshape(-1, order="F").astype("<f8").tofile(f"{out}/perturbation.bin")
    for name in ("perturbation.txt", "times.txt"):
        shutil.copy(f"{src}/{name}", f"{out}/{name}")


def self_check(tmp):
    # a chain v_k = c_k A v_{k-1} with a single growing mode: after rescaling, consecutive columns grow by the true factor
    lam, c = 1.3, [0.7, 0.8, 0.9]
    w = np.array([1.0, 2.0]); cols = []
    for k in range(len(c) + 1):
        d = lam * w; cols.append(d)          # d_k = A w_{k-1}
        if k < len(c):
            w = c[k] * d                     # w_k = c_k d_k
    X = np.array(cols).T
    os.makedirs(f"{tmp}/in", exist_ok=True)
    X.reshape(-1, order="F").astype("<f8").tofile(f"{tmp}/in/perturbation.bin")
    open(f"{tmp}/in/perturbation.txt", "w").write(f"n_space = 2\nn_snap = {X.shape[1]}\n")
    np.savetxt(f"{tmp}/in/times.txt", np.arange(X.shape[1]))
    f = 1 / np.concatenate(([1.0], np.cumprod(c)))
    np.savetxt(f"{tmp}/s.txt", f)
    rescale(f"{tmp}/in", f"{tmp}/s.txt", f"{tmp}/out")
    Y = np.fromfile(f"{tmp}/out/perturbation.bin", "<f8").reshape((2, X.shape[1]), order="F")
    assert np.allclose(Y[:, 1:] / Y[:, :-1], lam), Y[:, 1:] / Y[:, :-1]


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-check":
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            self_check(t)
        print("self-check passed")
    elif len(sys.argv) == 4:
        rescale(*sys.argv[1:])
    else:
        sys.exit(__doc__)
