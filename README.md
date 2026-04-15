# RIS-UBMRA: Upper Bound-guided Minimum RIS Activation for V2X Relay Systems
---

## Overview

This repository provides the official Python implementation of the **UB-MRA** (Upper Bound-guided Minimum RIS Activation) strategy for RIS-assisted vehicle platooning communication systems.

### The core problem

Traditional RIS deployments activate **all elements all the time** (e.g., Fixed-256). This wastes energy when fewer elements would already satisfy the QoS requirements of every vehicle in the platoon.

### Our solution — UB-MRA

For each lead-vehicle position, UB-MRA uses the **large-scale channel upper bound** (via Jensen's inequality) to quickly pre-select the **minimum number of RIS elements** that satisfies two constraints simultaneously:

| Constraint | Meaning |
|---|---|
| Tail-vehicle rate upper bound ≥ 7.5 bits/s/Hz | Weakest vehicle in the platoon must be served |
| Jain fairness index ≥ 0.97 | Rate gap between head and tail must stay small |

Only then does the system activate that number of elements and run the actual link.

### Key results

| Scheme | Avg Active N | Avg Rate | Tail P5 | Outage | Norm. EE |
|---|---|---|---|---|---|
| Fixed-256 | 256 | 8.63 | 6.55 | 1.34% | 1.000 |
| Fixed-128 | 128 | 7.32 | 5.64 | 11.86% | 1.209 |
| **UB-MRA (ours)** | **173.7** | **7.52** | **6.44** | **1.69%** | **1.131** |

- Reduces active RIS elements by **32.1%** vs Fixed-256  
- Tail-user outage improves by **85.8%** vs Fixed-128  
- Normalized energy efficiency improves by **13.0%** vs Fixed-256

---

## System Model

```
 [S] Base Station (0, 0, 15 m)
      |
      |──── Direct path (blocked) ────────────────────────┐
      |                                                    |
      |──→ [RIS] (75, 10, 8 m) ──→ [Car 1] ... [Car 4]  ←──── [R] Relay (150, 0, 8 m)
            N_act elements             ← platoon →
```

- **S → RIS → Vehicle**: RIS-assisted path with ideal phase alignment  
- **S → R → Vehicle**: DF (Decode-and-Forward) relay path, half-duplex  
- **Total rate**: `R = R_DF + R_RIS`
- **Path loss**: 3GPP TR 38.901 UMi Street Canyon LOS  
- **Fading**: Rician, K = 5

---

## Installation

```bash
git clone https://github.com/v0yage1969/ris-v2x-ubmra-optimizer.git
cd ris-v2x-ubmra-optimizer
pip install -r requirements.txt
```

**Requirements:** Python ≥ 3.10, numpy, matplotlib, scipy, pandas, tqdm

---

## Quick Start

### 1 — Single-point advisor (recommended first try)

Get an instant RIS activation recommendation for a specific platoon position:

```bash
python ris_ubmra_optimizer.py --x_head 150
```

Output:
```
  UB-MRA pre-selection result
  Minimum N satisfying constraints : N_act = 64
  Tail-vehicle rate upper bound    : 8.426 bits/s/Hz  (threshold ≥ 7.5)
  Platoon fairness upper bound     : 0.985            (threshold ≥ 0.97)

  Monte Carlo results (200 runs)
  Scheme         N_act  AvgRate   TailP5  Fairness  Outage%
  Fixed-256        256   9.634   10.988    0.9776    0.00%  EE=1.000
  Fixed-128        128   7.958    9.100    0.9765    0.00%  EE=1.175
  UB-MRA            64   6.759    7.074    0.9853    0.50%  EE=1.265

  Recommendation: activate 64 RIS elements (75.0% fewer than full config)
```

### 2 — Full sweep simulation (reproduces all paper figures)

```bash
python ris_ubmra_optimizer.py
```

Runs Monte Carlo over all head-vehicle positions (60 m → 260 m, step 10 m), produces 4 plots and 3 CSV files.

### 3 — Use as a Python library

```python
from ris_ubmra_optimizer import UBMRAOptimizer, DEFAULT

# Create optimizer with default parameters
opt = UBMRAOptimizer(DEFAULT)

# Query the UB-MRA strategy for a given platoon position
result = opt.select_N_act(x_head=180)
print(f"Recommended N_act: {result['N_act']}")
print(f"Tail-rate upper bound: {result['R_tail_UB']:.3f} bits/s/Hz")

# Run Monte Carlo to get real performance metrics
perf = opt.simulate(x_head=180, N_act=result['N_act'], seed=42)
print(f"Average platoon rate : {perf['avg_rate']:.3f} bits/s/Hz")
print(f"Tail P5 rate         : {perf['tail_p5']:.3f} bits/s/Hz")
print(f"Tail outage          : {perf['outage_pct']:.2f}%")
```

---

## Command-line Options

| Argument | Default | Description |
|---|---|---|
| `--x_head` | None | Single-point mode: head vehicle x-position (m). Omit for full sweep. |
| `--n_mc` | 500 | Monte Carlo iterations per position |
| `--x_min` | 60 | Sweep start position (m) |
| `--x_max` | 260 | Sweep end position (m) |
| `--x_step` | 10 | Sweep step size (m) |
| `--n_cars` | 4 | Number of vehicles in platoon |
| `--d_car` | 20 | Inter-vehicle spacing (m) |
| `--fc` | 3.5 | Carrier frequency (GHz) |
| `--R_ub` | 7.5 | Tail-rate upper bound threshold (bits/s/Hz) |
| `--J_min` | 0.97 | Jain fairness index threshold |
| `--no_plot` | False | Skip figure output (CSV only) |
| `--out_dir` | `.` | Output directory for CSV and PNG files |

---

## Output Files

After a full sweep, four files are saved to `--out_dir`:

| File | Content |
|---|---|
| `RIS_UBMRA_results.png` | 4-panel comparison figure |
| `RIS_UBMRA_position_results.csv` | Per-position metrics for all three schemes |
| `RIS_UBMRA_summary_results.csv` | Aggregate statistics across all positions |
| `RIS_UBMRA_evaluation_stats.csv` | Paired t-test results for statistical significance |

---

## Customising Parameters

All system parameters are defined in the `DEFAULT` dictionary at the top of `ris_ubmra_optimizer.py`. You can edit them directly or override via CLI:

```python
DEFAULT = dict(
    fc_GHz   = 3.5,          # Carrier frequency (GHz)
    S_pos    = (0,   0, 15), # Base station 3D position (m)
    R_pos    = (150, 0,  8), # DF relay 3D position (m)
    RIS_pos  = (75, 10,  8), # RIS panel centre (m)
    N_cars   = 4,            # Platoon size
    d_car    = 20,           # Inter-vehicle spacing (m)
    N_cand   = (64, 128, 192, 256),  # Candidate activation sizes
    R_UB_thresh = 7.5,       # Upper-bound pre-selection threshold
    J_thresh    = 0.97,      # Jain fairness threshold
    R_out_thresh = 6.0,      # Outage rate threshold (bits/s/Hz)
    K_rice   = 5,            # Rician K factor
    N_MC     = 500,          # Monte Carlo iterations
    ...
)
```

---

## Algorithm: UB-MRA

```
For each head-vehicle position x_head:
  1. Compute large-scale path gains (3GPP UMi LOS path loss)
  2. For N in {64, 128, 192, 256}:
       Compute Jensen upper bound on tail-vehicle rate:
         γ_RIS_UB  = ρ · η² · N² · PL(S→RIS) · PL(RIS→tail)
         R_tail_UB = log₂(1 + γ_RIS_UB) + ½·log₂(1 + min(γ_SR, γ_RD))
       Compute Jain fairness index on upper-bound rates
       If R_tail_UB ≥ 7.5  AND  J ≥ 0.97:
           N_act = N  ← select this minimum N and stop
  3. Activate exactly N_act elements
  4. (Evaluation) Run 500-round Monte Carlo with Rician fading
```

The upper bound is **intentionally optimistic**: it sets a stricter pre-selection gate (7.5 bits/s/Hz) than the final reliability requirement (6 bits/s/Hz) to compensate for the optimism gap.

---

## Project Structure

```
ris-v2x-ubmra-optimizer/
├── ris_ubmra_optimizer.py   # Main optimizer — run this
├── requirements.txt          # Python dependencies
├── README.md
└── LICENSE
```

---

## References

1. 3GPP TR 38.901, *Study on channel model for frequencies from 0.5 to 100 GHz*
2. 3GPP TR 22.886, *Study on enhancement of 3GPP support for 5G V2X services*
3. Zhang et al., "Active RIS vs. Passive RIS: Which Will Prevail in 6G?", *IEEE Trans. Commun.*, 2023
4. CN 114584587 B, *RIS与中继结合的协同车联网部署方案*

> ⚠️ **Note:** All simulation data are generated from the stochastic channel model described above. They are reproducible simulation results, not field measurements.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{ris_ubmra_2025,
  title   = {RIS-Assisted V2X Relay Systems: UB-MRA Strategy},
  author  = {<Your Names>},
  year    = {2025},
  url     = {https://github.com/v0yage1969/ris-v2x-ubmra-optimizer},
  note    = {SRTP Project, SWJTU}
}
```
