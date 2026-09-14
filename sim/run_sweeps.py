"""Sequential cavity sweeps: N_mirror (loading), gmax & W_y (convergence),
plus suspended-vs-on-substrate comparison. Writes JSONs to results/."""
import json, subprocess, sys, time

A = "193.4"
jobs = []
# N_mirror sweep at gmax=2, W_y=4 (Q vs N -> mirror transmission)
for N in (4, 6, 8, 12):
    jobs.append([A, "2.0", "130", "0.14", str(N), "4.0"])
# gmax convergence at N=10
for g in (1.5, 2.5, 3.0):
    jobs.append([A, str(g), "130", "0.14", "10", "4.0"])
# W_y convergence
jobs.append([A, "2.0", "170", "0.14", "10", "5.0"])
for j in jobs:
    t = time.time()
    r = subprocess.run([sys.executable, "run_cavity.py"] + j,
                       capture_output=True, text=True)
    tail = r.stdout.strip().splitlines()[-22:]
    ok = "cavity mode" in r.stdout
    print(" ".join(j), f"[{time.time()-t:.0f}s]", "OK" if ok else "FAIL", flush=True)
    if not ok:
        print(r.stdout[-1500:], r.stderr[-1500:], flush=True)
