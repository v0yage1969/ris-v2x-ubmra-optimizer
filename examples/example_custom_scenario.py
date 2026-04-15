"""
Example 2 — Custom scenario: change nodes, platoon size, and thresholds
=======================================================================
Shows how to adapt UB-MRA to a different deployment scenario
(e.g. larger platoon, different RIS location, tighter fairness requirement).

Run from the repo root:
    python examples/example_custom_scenario.py
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))

import numpy as np
import matplotlib.pyplot as plt
from ris_ubmra_optimizer import UBMRAOptimizer, DEFAULT

# ── Custom scenario: 6-vehicle platoon, tighter QoS ─────────────────────
custom_params = dict(DEFAULT)
custom_params.update(dict(
    N_cars       = 6,          # 6 vehicles instead of 4
    d_car        = 15,         # tighter spacing
    R_UB_thresh  = 8.0,        # stricter tail-rate threshold
    J_thresh     = 0.98,       # stricter fairness
    R_out_thresh = 6.5,        # stricter outage gate
    N_MC         = 300,
))

opt = UBMRAOptimizer(custom_params)

# ── Sweep head-vehicle positions ─────────────────────────────────────────
x_vec    = np.arange(60, 261, 10, dtype=float)
n_act    = []
avg_rate = []
tail_p5  = []
outage   = []

print("Sweeping platoon positions …")
for x in x_vec:
    sel  = opt.select_N_act(x)
    perf = opt.simulate(x, sel['N_act'], seed=7)
    n_act.append(sel['N_act'])
    avg_rate.append(perf['avg_rate'])
    tail_p5.append(perf['tail_p5'])
    outage.append(perf['outage_pct'])
    print(f"  x={x:.0f}m  N_act={sel['N_act']:3d}  rate={perf['avg_rate']:.2f}  outage={perf['outage_pct']:.1f}%")

# ── Plot ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
fig.suptitle("Custom Scenario: 6-Vehicle Platoon, Tighter QoS", fontsize=13)

axes[0].plot(x_vec, avg_rate, 'g-o', markersize=5)
axes[0].set(title="Avg Platoon Rate", xlabel="x_head (m)",
            ylabel="Rate (bits/s/Hz)")

axes[1].plot(x_vec, tail_p5, 'g-D', markersize=5)
axes[1].axhline(custom_params['R_out_thresh'], color='k', linestyle='--',
                label='Outage threshold')
axes[1].set(title="Tail P5 Rate", xlabel="x_head (m)",
            ylabel="Rate (bits/s/Hz)")
axes[1].legend()

axes[2].step(x_vec, n_act, where='post', color='steelblue', linewidth=2)
axes[2].set_yticks(custom_params['N_cand'])
axes[2].set(title="UB-MRA Selected N_act", xlabel="x_head (m)",
            ylabel="Active RIS elements")

plt.tight_layout()
plt.savefig("examples/custom_scenario_results.png", dpi=130)
print("\nPlot saved to examples/custom_scenario_results.png")
plt.show()
