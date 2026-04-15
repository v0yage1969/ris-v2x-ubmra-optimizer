"""
Example 1 — Single-point RIS activation advisor
================================================
Given a head-vehicle position, the UB-MRA strategy instantly tells you:
  - How many RIS elements to activate
  - Expected tail-vehicle rate and outage probability

Run from the repo root:
    python examples/example_single_point.py
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))

from ris_ubmra_optimizer import UBMRAOptimizer, DEFAULT

# ── Choose a head-vehicle position ──────────────────────────────────────
x_head = 200   # metres; try different values between 60 and 260

opt = UBMRAOptimizer(DEFAULT)

# Step 1: UB-MRA pre-selection (fast, no Monte Carlo needed)
sel = opt.select_N_act(x_head)
print(f"\n--- UB-MRA recommendation for x_head = {x_head} m ---")
print(f"Activate  : {sel['N_act']} RIS elements")
print(f"Tail UB   : {sel['R_tail_UB']:.3f} bits/s/Hz  (threshold >= {opt.R_UB_thresh})")
print(f"Fairness  : {sel['J_UB']:.4f}               (threshold >= {opt.J_thresh})")

# Step 2: Monte Carlo to verify real performance
perf = opt.simulate(x_head, sel['N_act'], seed=0)
print(f"\n--- Monte Carlo verification ({opt.N_MC} runs) ---")
print(f"Avg platoon rate  : {perf['avg_rate']:.3f} bits/s/Hz")
print(f"Tail P5 rate      : {perf['tail_p5']:.3f} bits/s/Hz")
print(f"Tail outage       : {perf['outage_pct']:.2f}%")
print(f"Fairness (Jain)   : {perf['fairness']:.4f}")

# Comparison with always-on Fixed-256
perf_full = opt.simulate(x_head, 256, seed=0)
ee = opt.norm_ee(perf['avg_rate'], sel['N_act'],
                 ref_rate=perf_full['avg_rate'], ref_N=256)
print(f"\n--- vs Fixed-256 ---")
print(f"Elements saved    : {256 - sel['N_act']} ({(256 - sel['N_act'])/256*100:.1f}%)")
print(f"Rate change       : {(perf['avg_rate'] - perf_full['avg_rate'])/perf_full['avg_rate']*100:+.1f}%")
print(f"Norm. energy eff. : {ee:.3f}x")
