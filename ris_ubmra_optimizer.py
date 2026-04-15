"""
=======================================================================
RIS辅助车联网中继系统：UB-MRA 策略优化程序
RIS-Assisted V2X Relay System: UB-MRA Strategy Optimizer
=======================================================================

核心创新（UB-MRA）：
    对每个头车位置，先用大尺度信道上界快速预判，找出满足
    "队尾车辆QoS + 编队公平性" 的最小RIS激活单元数 N_act，
    避免始终满配带来的不必要能耗。

三种对比方案：
    Fixed-256   : 始终激活 256 个 RIS 元件（性能上限基准）
    Fixed-128   : 始终激活 128 个 RIS 元件（节能基准）
    UB-MRA      : 自适应按需激活（本文提出）

运行方式：
    pip install numpy matplotlib scipy pandas tqdm
    python ris_ubmra_optimizer.py              # 完整仿真 + 出图
    python ris_ubmra_optimizer.py --x_head 150  # 单点分析
    python ris_ubmra_optimizer.py --n_mc 1000   # 增加精度
    python ris_ubmra_optimizer.py --help         # 查看所有选项

参考：
    [1] 3GPP TR 38.901 UMi Street Canyon LOS 路径损耗模型
    [2] CN 114584587 B，混合 RIS + DF 中继车联网方案
=======================================================================
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy import stats
import argparse
import sys
import time
from pathlib import Path

# ── Windows 控制台 UTF-8 输出 ──────────────────────────────────────────
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace', line_buffering=True)
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8',
                                  errors='replace', line_buffering=True)

# ── 进度条（可选）──────────────────────────────────────────────────────
try:
    from tqdm import tqdm
    _TQDM = True
except ImportError:
    _TQDM = False

matplotlib.rcParams.update({
    'font.size': 11,
    'axes.grid': True,
    'grid.alpha': 0.35,
    'lines.linewidth': 2,
    'lines.markersize': 7,
    'figure.dpi': 120,
})

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  默认系统参数（可直接修改此区域，或通过命令行覆盖）
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DEFAULT = dict(
    # 载波与噪声
    fc_GHz   = 3.5,      # 载波频率 (GHz)
    BW       = 20e6,     # 系统带宽 (Hz)
    NF_dB    = 7.0,      # 接收机噪声系数 (dB)

    # 节点三维坐标 [x, y, z]，单位 m
    S_pos    = (0,   0,  15),   # 源节点 / 基站
    R_pos    = (150, 0,   8),   # DF 中继
    RIS_pos  = (75,  10,  8),   # RIS 面板中心

    # 编队
    N_cars   = 4,        # 车辆数量
    d_car    = 20,       # 纵向车间距 (m)
    z_car    = 1.5,      # 车辆天线高度 (m)

    # 头车位置扫描
    x_min    = 60,
    x_max    = 260,
    x_step   = 10,

    # UB-MRA 候选激活规模与阈值
    N_cand         = (64, 128, 192, 256),
    R_UB_thresh    = 7.5,   # 上界预选门限 (bits/s/Hz)，比最终可靠门限严格
    J_thresh       = 0.97,  # Jain 公平性指数门限
    R_out_thresh   = 6.0,   # 队尾中断判据 (bits/s/Hz)

    # 信道
    K_rice   = 5,        # Rician K 因子
    eta      = 1.0,      # RIS 反射幅度（理想值）

    # 蒙特卡洛
    N_MC     = 500,

    # 功耗（工程假设，用于归一化能效相对对比）
    P_s      = 1.0,      # 源发射功率 (W)
    P_r      = 1.0,      # 中继发射功率 (W)
    P_c      = 1.5,      # 电路功耗 (W)
    P_elem   = 0.02,     # 每个 RIS 元件功耗 (W)
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  信道与路径损耗工具函数
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def path_loss_lin(d3d: float, fc_GHz: float) -> float:
    """
    3GPP TR 38.901 UMi Street Canyon LOS 路径损耗（线性）
    PL [dB] = 32.4 + 21·log10(d) + 20·log10(fc_GHz)
    """
    d3d = max(d3d, 1.0)          # 避免 log(0)
    pl_dB = 32.4 + 21 * np.log10(d3d) + 20 * np.log10(fc_GHz)
    return 10 ** (-pl_dB / 10)


def rician_channel(K: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """
    生成 Rician 复信道系数，长度 n，大尺度增益归一化为 1（E[|h|²]=1）。
    实际使用时乘以 sqrt(path_loss) 引入大尺度衰落。
    """
    los  = np.sqrt(K / (K + 1)) * np.exp(1j * 2 * np.pi * rng.random(n))
    nlos = np.sqrt(1 / (K + 1)) * (rng.standard_normal(n) +
                                    1j * rng.standard_normal(n)) / np.sqrt(2)
    return los + nlos


def jain_index(rates: np.ndarray) -> float:
    """Jain 公平性指数 J = (Σr)² / (N·Σr²)"""
    s  = rates.sum()
    ss = (rates ** 2).sum()
    return float(s ** 2 / (len(rates) * ss + 1e-30))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  UB-MRA 核心优化器
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class UBMRAOptimizer:
    """
    UB-MRA 策略优化器。

    给定系统参数后，对任意头车位置可：
      - select_N_act(x_head)  → 快速上界预选最小激活规模
      - simulate(x_head, N_act)→ 蒙特卡洛评估真实性能
    """

    def __init__(self, params: dict):
        p = params
        self.fc_GHz  = p['fc_GHz']
        self.S  = np.array(p['S_pos'],   dtype=float)
        self.R  = np.array(p['R_pos'],   dtype=float)
        self.I  = np.array(p['RIS_pos'], dtype=float)

        self.N_cars = p['N_cars']
        self.d_car  = p['d_car']
        self.z_car  = p['z_car']

        self.N_cand      = sorted(p['N_cand'])
        self.R_UB_thresh = p['R_UB_thresh']
        self.J_thresh    = p['J_thresh']
        self.R_out_thresh = p['R_out_thresh']

        self.K_rice = p['K_rice']
        self.eta    = p['eta']
        self.N_MC   = p['N_MC']

        # 噪声与 SNR
        kB     = 1.38e-23
        T0     = 290
        NF_lin = 10 ** (p['NF_dB'] / 10)
        sigma2 = kB * T0 * p['BW'] * NF_lin
        self.rho_s = p['P_s'] / sigma2
        self.rho_r = p['P_r'] / sigma2

        # 功耗
        self.P_s    = p['P_s']
        self.P_r    = p['P_r']
        self.P_c    = p['P_c']
        self.P_elem = p['P_elem']

        # 预计算固定链路路径损耗
        self.pl_SR = path_loss_lin(np.linalg.norm(self.S - self.R), self.fc_GHz)
        self.pl_SI = path_loss_lin(np.linalg.norm(self.S - self.I), self.fc_GHz)

    # ── 编队位置 ──────────────────────────────────────────────────────
    def platoon_positions(self, x_head: float) -> np.ndarray:
        """返回各车三维坐标，形状 (N_cars, 3)，车辆1=头车"""
        xs = x_head - np.arange(self.N_cars) * self.d_car
        return np.column_stack([xs,
                                np.zeros(self.N_cars),
                                np.full(self.N_cars, self.z_car)])

    def _link_gains(self, x_head: float):
        """计算各车到 RIS 和中继的路径损耗（线性）"""
        cars = self.platoon_positions(x_head)
        pl_ID = np.array([path_loss_lin(np.linalg.norm(self.I - c), self.fc_GHz)
                           for c in cars])
        pl_RD = np.array([path_loss_lin(np.linalg.norm(self.R - c), self.fc_GHz)
                           for c in cars])
        return pl_ID, pl_RD

    # ── 上界预选 ──────────────────────────────────────────────────────
    def select_N_act(self, x_head: float) -> dict:
        """
        UB-MRA 预选算法：
          遍历 N_cand（从小到大），找到第一个满足
          R_tail_UB >= R_UB_thresh 且 J_UB >= J_thresh 的 N，
          返回详细的预选结果。
        """
        pl_ID, pl_RD = self._link_gains(x_head)

        gamma_SR_ls = self.rho_s * self.pl_SR

        for N in self.N_cand:
            # 各车 RIS 速率上界（Jensen 不等式 + 大尺度上界）
            gamma_RIS_UB = self.rho_s * (self.eta ** 2) * (N ** 2) * self.pl_SI * pl_ID
            R_RIS_UB = np.log2(1 + gamma_RIS_UB)

            # 各车 DF 速率（仅大尺度，无小尺度波动）
            R_DF_UB = np.array([
                0.5 * np.log2(1 + min(gamma_SR_ls, self.rho_r * pl_RD[k]))
                for k in range(self.N_cars)
            ])

            R_sum_UB = R_RIS_UB + R_DF_UB
            R_tail_UB = float(R_sum_UB[-1])   # 尾车（最后一辆）
            J_UB = jain_index(R_sum_UB)

            if R_tail_UB >= self.R_UB_thresh and J_UB >= self.J_thresh:
                return dict(N_act=N, R_tail_UB=R_tail_UB, J_UB=J_UB,
                            R_sum_UB=R_sum_UB, satisfied=True)

        # 全部候选都不满足 → 选最大
        N = self.N_cand[-1]
        gamma_RIS_UB = self.rho_s * (self.eta ** 2) * (N ** 2) * self.pl_SI * pl_ID
        R_sum_UB = np.log2(1 + gamma_RIS_UB) + np.array([
            0.5 * np.log2(1 + min(gamma_SR_ls, self.rho_r * pl_RD[k]))
            for k in range(self.N_cars)
        ])
        return dict(N_act=N, R_tail_UB=float(R_sum_UB[-1]),
                    J_UB=jain_index(R_sum_UB), R_sum_UB=R_sum_UB,
                    satisfied=False)

    # ── 蒙特卡洛评估 ─────────────────────────────────────────────────
    def simulate(self, x_head: float, N_act: int,
                 seed: int | None = None) -> dict:
        """
        对给定头车位置和激活规模进行 N_MC 次蒙特卡洛仿真。
        返回：平均编队速率、队尾P5、公平性、中断率。
        """
        rng = np.random.default_rng(seed)
        pl_ID, pl_RD = self._link_gains(x_head)

        rates = np.zeros((self.N_MC, self.N_cars))

        sq_SI = np.sqrt(self.pl_SI)
        sq_SR = np.sqrt(self.pl_SR)
        sq_ID = np.sqrt(pl_ID)
        sq_RD = np.sqrt(pl_RD)

        for m in range(self.N_MC):
            h_SR = sq_SR * rician_channel(self.K_rice, 1, rng)
            gamma_SR = self.rho_s * float(np.abs(h_SR[0]) ** 2)

            h_SI = sq_SI * rician_channel(self.K_rice, N_act, rng)

            for k in range(self.N_cars):
                h_ID = sq_ID[k] * rician_channel(self.K_rice, N_act, rng)
                h_RD = sq_RD[k] * rician_channel(self.K_rice, 1, rng)

                # RIS 速率（理想相位对齐：所有反射路径同相叠加）
                eff_gain  = float(np.sum(np.abs(h_SI) * np.abs(h_ID)))
                gamma_RIS = self.rho_s * (self.eta ** 2) * eff_gain ** 2
                R_RIS     = np.log2(1 + gamma_RIS)

                # DF 中继速率（半双工，两时隙）
                gamma_RD = self.rho_r * float(np.abs(h_RD[0]) ** 2)
                R_DF     = 0.5 * np.log2(1 + min(gamma_SR, gamma_RD))

                rates[m, k] = R_DF + R_RIS

        tail_rates = rates[:, -1]   # 尾车速率序列

        J_vals = np.array([jain_index(rates[m]) for m in range(self.N_MC)])

        return dict(
            avg_rate   = float(rates.mean()),
            tail_p5    = float(np.percentile(tail_rates, 5)),
            fairness   = float(J_vals.mean()),
            outage_pct = float((tail_rates < self.R_out_thresh).mean() * 100),
            rates_all  = rates,
        )

    # ── 归一化能效 ────────────────────────────────────────────────────
    def norm_ee(self, avg_rate: float, N_act: int,
                ref_rate: float, ref_N: int = 256) -> float:
        """EE = AvgRate / P_total；结果归一化到 Fixed-ref_N 基准"""
        P      = self.P_s + self.P_r + self.P_c + N_act    * self.P_elem
        P_ref  = self.P_s + self.P_r + self.P_c + ref_N    * self.P_elem
        return (avg_rate / P) / (ref_rate / P_ref)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  单点分析（命令行 --x_head 模式）
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def single_point_analysis(opt: UBMRAOptimizer, x_head: float):
    """对单个头车位置做完整分析并打印报告"""
    print(f"\n{'='*60}")
    print(f"  单点分析  |  头车位置 x_head = {x_head:.0f} m")
    print(f"{'='*60}")

    # 编队信息
    cars = opt.platoon_positions(x_head)
    print(f"\n  编队车辆位置（共 {opt.N_cars} 辆）：")
    for k, c in enumerate(cars):
        label = "头车" if k == 0 else ("尾车" if k == opt.N_cars - 1 else f"车{k+1}")
        print(f"    [{label}] ({c[0]:.0f}, {c[1]:.0f}, {c[2]:.1f}) m")

    # UB-MRA 预选
    sel = opt.select_N_act(x_head)
    print(f"\n  ── UB-MRA 预选结果 ──────────────────────────")
    print(f"  最小满足约束的激活规模  : N_act = {sel['N_act']}")
    print(f"  队尾速率上界            : {sel['R_tail_UB']:.3f} bits/s/Hz  "
          f"（门限 ≥ {opt.R_UB_thresh}）")
    print(f"  编队公平性上界          : {sel['J_UB']:.4f}  "
          f"（门限 ≥ {opt.J_thresh}）")
    if not sel['satisfied']:
        print("  ⚠ 警告：所有候选规模均未满足约束，已选最大值")

    # 三方案 MC 仿真
    schemes = [('Fixed-256', 256), ('Fixed-128', 128),
               ('UB-MRA',   sel['N_act'])]
    print(f"\n  ── 蒙特卡洛仿真（{opt.N_MC} 次）─────────────────")
    print(f"  {'方案':<14} {'N_act':>6} {'AvgRate':>9} {'TailP5':>9} "
          f"{'Fairness':>10} {'Outage%':>9}")
    print(f"  {'-'*60}")

    results = {}
    ref_rate = None
    for name, N in schemes:
        res = opt.simulate(x_head, N, seed=42)
        results[name] = res
        if name == 'Fixed-256':
            ref_rate = res['avg_rate']
        ee = opt.norm_ee(res['avg_rate'], N, ref_rate)
        print(f"  {name:<14} {N:>6} {res['avg_rate']:>9.4f} "
              f"{res['tail_p5']:>9.4f} {res['fairness']:>10.4f} "
              f"{res['outage_pct']:>9.2f}%   EE={ee:.3f}")

    # 策略建议
    print(f"\n  ── 策略建议 ──────────────────────────────────")
    N_rec = sel['N_act']
    r_ubmra = results['UB-MRA']
    r_f256  = results['Fixed-256']
    delta_n = (N_rec - 256) / 256 * 100
    delta_r = (r_ubmra['avg_rate'] - r_f256['avg_rate']) / r_f256['avg_rate'] * 100
    print(f"  建议激活 {N_rec} 个 RIS 元件（较满配减少 {abs(delta_n):.1f}%）")
    print(f"  预计平均速率变化：{delta_r:+.1f}%")
    print(f"  预计队尾5%分位速率：{r_ubmra['tail_p5']:.3f} bits/s/Hz")
    print(f"  预计队尾中断率：{r_ubmra['outage_pct']:.2f}%")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  完整扫描仿真
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def full_simulation(opt: UBMRAOptimizer, params: dict) -> pd.DataFrame:
    """对所有头车位置运行完整 UB-MRA 仿真，返回位置级 DataFrame"""
    x_vec  = np.arange(params['x_min'], params['x_max'] + 1, params['x_step'])
    n_pos  = len(x_vec)
    schemes = [('Fixed-256', 256), ('Fixed-128', 128), ('UB-MRA', None)]

    records = []
    it = tqdm(range(n_pos), desc='仿真进度') if _TQDM else range(n_pos)

    for i in it:
        x_head = float(x_vec[i])
        sel    = opt.select_N_act(x_head)
        N_ubmra = sel['N_act']

        row = {'x_head': x_head, 'N_act_UBMRA': N_ubmra}

        for name, N in schemes:
            N_use = N_ubmra if N is None else N
            res = opt.simulate(x_head, N_use)
            key = name.replace('-', '')
            row[f'AvgRate_{key}']  = res['avg_rate']
            row[f'TailP5_{key}']   = res['tail_p5']
            row[f'Fairness_{key}'] = res['fairness']
            row[f'Outage_{key}']   = res['outage_pct']

        if not _TQDM:
            pct = (i + 1) / n_pos * 100
            print(f'  [{pct:5.1f}%] x_head={x_head:.0f}m  N_act={N_ubmra}'
                  f'  AvgRate(UB)={row["AvgRate_UBMRA"]:.3f}', end='\r')

        records.append(row)

    if not _TQDM:
        print()
    return pd.DataFrame(records)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  汇总统计与显著性检验
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_summary(df: pd.DataFrame, opt: UBMRAOptimizer) -> dict:
    avg_N = {
        'Fixed256': 256.0,
        'Fixed128': 128.0,
        'UBMRA':    df['N_act_UBMRA'].mean(),
    }
    avg_rate = {k: df[f'AvgRate_{k}'].mean() for k in avg_N}
    tail_p5  = {k: df[f'TailP5_{k}'].mean()  for k in avg_N}
    fairness = {k: df[f'Fairness_{k}'].mean() for k in avg_N}
    outage   = {k: df[f'Outage_{k}'].mean()   for k in avg_N}

    # 归一化能效
    ee_norm = {}
    ref_rate = avg_rate['Fixed256']
    for k, N in avg_N.items():
        ee_norm[k] = opt.norm_ee(avg_rate[k], N, ref_rate)

    # 配对 t 检验
    def paired_t(a, b):
        t, p = stats.ttest_rel(a, b)
        return float(t), float(p)

    tests = {
        'AvgRate: UB-MRA vs Fixed-256': paired_t(
            df['AvgRate_UBMRA'], df['AvgRate_Fixed256']),
        'AvgRate: UB-MRA vs Fixed-128': paired_t(
            df['AvgRate_UBMRA'], df['AvgRate_Fixed128']),
        'Outage:  UB-MRA vs Fixed-128': paired_t(
            df['Outage_UBMRA'],  df['Outage_Fixed128']),
        'TailP5:  UB-MRA vs Fixed-128': paired_t(
            df['TailP5_UBMRA'],  df['TailP5_Fixed128']),
    }

    return dict(avg_N=avg_N, avg_rate=avg_rate, tail_p5=tail_p5,
                fairness=fairness, outage=outage, ee_norm=ee_norm,
                tests=tests)


def print_summary(summ: dict):
    keys  = ['Fixed256', 'Fixed128', 'UBMRA']
    names = {'Fixed256': 'Fixed-256', 'Fixed128': 'Fixed-128',
             'UBMRA': 'Proposed-UBMRA'}

    print(f"\n{'='*70}")
    print("  汇总结果（全位置平均）")
    print(f"{'='*70}")
    hdr = f"  {'方案':<20} {'Avg N':>8} {'AvgRate':>9} {'TailP5':>9} " \
          f"{'Fairness':>10} {'Outage%':>9} {'NormEE':>8}"
    print(hdr)
    print(f"  {'-'*67}")
    for k in keys:
        print(f"  {names[k]:<20} {summ['avg_N'][k]:>8.2f} "
              f"{summ['avg_rate'][k]:>9.4f} {summ['tail_p5'][k]:>9.4f} "
              f"{summ['fairness'][k]:>10.4f} {summ['outage'][k]:>9.4f} "
              f"{summ['ee_norm'][k]:>8.4f}")

    ub  = summ['avg_N']['UBMRA']
    f256 = summ['avg_N']['Fixed256']
    print(f"\n  [关键结论]  较 Fixed-256：")
    dn = (ub - f256) / f256 * 100
    dr = (summ['avg_rate']['UBMRA'] - summ['avg_rate']['Fixed256']) / \
          summ['avg_rate']['Fixed256'] * 100
    dt = (summ['tail_p5']['UBMRA'] - summ['tail_p5']['Fixed256']) / \
          summ['tail_p5']['Fixed256'] * 100
    de = (summ['ee_norm']['UBMRA'] - 1.0) * 100
    print(f"    激活单元数  : {dn:+.1f}%  ({ub:.1f} vs 256)")
    print(f"    平均速率    : {dr:+.1f}%")
    print(f"    队尾P5速率  : {dt:+.1f}%")
    print(f"    归一化能效  : {de:+.1f}%")

    print(f"\n  [统计检验]")
    for name, (t, p) in summ['tests'].items():
        sig = "✓ 显著 (p<0.05)" if p < 0.05 else "— 不显著"
        print(f"    {name:<42}  t={t:+.3f}  p={p:.4f}  {sig}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  可视化
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

COLORS  = {'Fixed256': '#1f77b4', 'Fixed128': '#ff7f0e', 'UBMRA': '#2ca02c'}
MARKERS = {'Fixed256': 'o', 'Fixed128': 's', 'UBMRA': 'D'}
LABELS  = {'Fixed256': 'Fixed-256', 'Fixed128': 'Fixed-128',
           'UBMRA': 'Proposed-UBMRA'}


def plot_results(df: pd.DataFrame, opt: UBMRAOptimizer,
                 out_dir: Path):
    """绘制四张对比图并保存"""
    x = df['x_head'].values
    keys = ['Fixed256', 'Fixed128', 'UBMRA']

    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    fig.suptitle('RIS辅助车联网：UB-MRA 策略优化结果', fontsize=14, fontweight='bold')
    gs = gridspec.GridSpec(2, 2, figure=fig)

    axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(2)]

    # ── 图1：平均编队速率 ──
    ax = axes[0]
    for k in keys:
        ax.plot(x, df[f'AvgRate_{k}'],
                color=COLORS[k], marker=MARKERS[k], markevery=2,
                label=LABELS[k])
    ax.set_xlabel('Head-vehicle position $x_{head}$ (m)')
    ax.set_ylabel('Average platoon rate (bits/s/Hz)')
    ax.set_title('Fig.1  Average Platoon Rate')
    ax.legend(fontsize=9)

    # ── 图2：队尾5%分位速率 ──
    ax = axes[1]
    for k in keys:
        ax.plot(x, df[f'TailP5_{k}'],
                color=COLORS[k], marker=MARKERS[k], markevery=2,
                label=LABELS[k])
    ax.axhline(opt.R_out_thresh, color='k', linestyle='--', linewidth=1.2,
               label=f'Outage threshold ({opt.R_out_thresh} b/s/Hz)')
    ax.set_xlabel('Head-vehicle position $x_{head}$ (m)')
    ax.set_ylabel('Tail-user 5th percentile rate (bits/s/Hz)')
    ax.set_title('Fig.2  Tail-user P5 Rate  ← 最关键图')
    ax.legend(fontsize=9)

    # ── 图3：UB-MRA 激活规模 ──
    ax = axes[2]
    ax.step(x, df['N_act_UBMRA'], where='post',
            color=COLORS['UBMRA'], linewidth=2, label='UB-MRA selected N')
    ax.plot(x, df['N_act_UBMRA'], 'D', color=COLORS['UBMRA'], markersize=6)
    ax.axhline(256, color=COLORS['Fixed256'], linestyle='--',
               linewidth=1, label='Fixed-256')
    ax.axhline(128, color=COLORS['Fixed128'], linestyle='--',
               linewidth=1, label='Fixed-128')
    ax.set_yticks(sorted(opt.N_cand))
    ax.set_ylim(0, max(opt.N_cand) * 1.15)
    ax.set_xlabel('Head-vehicle position $x_{head}$ (m)')
    ax.set_ylabel('Selected active RIS elements')
    ax.set_title('Fig.3  UB-MRA 按需激活规模')
    ax.legend(fontsize=9)

    # ── 图4：队尾中断率 ──
    ax = axes[3]
    for k in keys:
        ax.plot(x, df[f'Outage_{k}'],
                color=COLORS[k], marker=MARKERS[k], markevery=2,
                label=LABELS[k])
    ax.set_xlabel('Head-vehicle position $x_{head}$ (m)')
    ax.set_ylabel('Tail-user outage probability (%)')
    ax.set_title(f'Fig.4  Outage (tail rate < {opt.R_out_thresh} b/s/Hz)')
    ax.legend(fontsize=9)

    save_path = out_dir / 'RIS_UBMRA_results.png'
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\n  图表已保存：{save_path}")
    plt.show()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  CSV 导出
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def export_csv(df: pd.DataFrame, summ: dict, tests: dict,
               out_dir: Path):
    # 位置级结果
    p1 = out_dir / 'RIS_UBMRA_position_results.csv'
    df.to_csv(p1, index=False, float_format='%.6f')

    # 汇总结果
    keys  = ['Fixed256', 'Fixed128', 'UBMRA']
    names = ['Fixed-256', 'Fixed-128', 'Proposed-UBMRA']
    rows  = []
    for k, n in zip(keys, names):
        rows.append({
            'Scheme':    n,
            'Avg_N_act': round(summ['avg_N'][k], 2),
            'Avg_Rate':  round(summ['avg_rate'][k], 4),
            'Tail_P5':   round(summ['tail_p5'][k], 4),
            'Fairness':  round(summ['fairness'][k], 4),
            'Outage_pct': round(summ['outage'][k], 4),
            'Norm_EE':   round(summ['ee_norm'][k], 4),
        })
    p2 = out_dir / 'RIS_UBMRA_summary_results.csv'
    pd.DataFrame(rows).to_csv(p2, index=False)

    # 统计检验
    test_rows = [{'Test': name, 't_stat': round(t, 4), 'p_value': round(p, 4)}
                 for name, (t, p) in summ['tests'].items()]
    p3 = out_dir / 'RIS_UBMRA_evaluation_stats.csv'
    pd.DataFrame(test_rows).to_csv(p3, index=False)

    print(f"  CSV 已保存：\n    {p1}\n    {p2}\n    {p3}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  命令行接口
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description='RIS V2X UB-MRA 策略优化程序',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument('--x_head', type=float, default=None,
                    help='单点分析：头车位置 (m)。不指定则运行完整扫描。')
    ap.add_argument('--n_mc',   type=int,   default=DEFAULT['N_MC'],
                    help='蒙特卡洛仿真次数')
    ap.add_argument('--x_min',  type=float, default=DEFAULT['x_min'],
                    help='头车位置扫描起点 (m)')
    ap.add_argument('--x_max',  type=float, default=DEFAULT['x_max'],
                    help='头车位置扫描终点 (m)')
    ap.add_argument('--x_step', type=float, default=DEFAULT['x_step'],
                    help='头车位置扫描步长 (m)')
    ap.add_argument('--n_cars', type=int,   default=DEFAULT['N_cars'],
                    help='编队车辆数')
    ap.add_argument('--d_car',  type=float, default=DEFAULT['d_car'],
                    help='车间距 (m)')
    ap.add_argument('--fc',     type=float, default=DEFAULT['fc_GHz'],
                    help='载波频率 (GHz)')
    ap.add_argument('--R_ub',   type=float, default=DEFAULT['R_UB_thresh'],
                    help='UB-MRA 队尾速率预选门限 (bits/s/Hz)')
    ap.add_argument('--J_min',  type=float, default=DEFAULT['J_thresh'],
                    help='编队公平性门限（Jain指数）')
    ap.add_argument('--no_plot', action='store_true',
                    help='不显示图形（仅导出CSV）')
    ap.add_argument('--out_dir', type=str, default='.',
                    help='输出文件目录')
    return ap


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#  主程序
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    args = build_parser().parse_args()

    # 将命令行参数合并到参数字典
    params = dict(DEFAULT)
    params['N_MC']        = args.n_mc
    params['N_cars']      = args.n_cars
    params['d_car']       = args.d_car
    params['fc_GHz']      = args.fc
    params['R_UB_thresh'] = args.R_ub
    params['J_thresh']    = args.J_min
    params['x_min']       = args.x_min
    params['x_max']       = args.x_max
    params['x_step']      = args.x_step
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  RIS V2X UB-MRA 策略优化程序")
    print("=" * 60)
    print(f"  载波频率   : {params['fc_GHz']} GHz")
    print(f"  节点位置   : S={params['S_pos']}  R={params['R_pos']}  RIS={params['RIS_pos']}")
    print(f"  编队规模   : {params['N_cars']} 辆 / 间距 {params['d_car']} m")
    print(f"  候选激活N  : {list(params['N_cand'])}")
    print(f"  上界门限   : R_tail ≥ {params['R_UB_thresh']} b/s/Hz，J ≥ {params['J_thresh']}")
    print(f"  MC次数     : {params['N_MC']}")

    opt = UBMRAOptimizer(params)

    # ── 单点分析模式 ──
    if args.x_head is not None:
        single_point_analysis(opt, args.x_head)
        return

    # ── 完整扫描模式 ──
    print(f"\n  头车位置扫描：{params['x_min']:.0f} ~ {params['x_max']:.0f} m，"
          f"步长 {params['x_step']:.0f} m\n")

    t0 = time.time()
    df = full_simulation(opt, params)
    elapsed = time.time() - t0
    print(f"\n  仿真完成，耗时 {elapsed:.1f} 秒")

    summ = compute_summary(df, opt)
    print_summary(summ)

    export_csv(df, summ, summ['tests'], out_dir)

    if not args.no_plot:
        plot_results(df, opt, out_dir)

    print("\n  完成！")


if __name__ == '__main__':
    main()
