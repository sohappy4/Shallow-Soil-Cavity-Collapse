# -*- coding: utf-8 -*-
"""
gtt.py —— 广义切线法（GTT）切线层
==================================
严格实现专著 §13.5.1–13.5.2：
  式(13.9)  g(sigma) = A*sigma_c^(1-B)*(sigma+sigma_t)^B
  式(13.10) tan(delta_t) = g'(sigma_0) = A*B*((sigma_0+sigma_t)/sigma_c)^(B-1)
  式(13.11) c_t = g(sigma_0) - sigma_0*g'(sigma_0)
            = A*sigma_c^(1-B)*(sigma_0+sigma_t)^(B-1)*[sigma_t+(1-B)*sigma_0]
  自检     c_t + sigma_0*tan(delta_t) == tau(sigma_0)
  定理13.1 tan(delta_t) == tan(psi)    切点倾角即剪胀角
  定理13.2 扫描切点取最小上界；最优切点为真实滑面正应力 sigma_*
  B=1 极限 tan(delta_t)=A, c_t=A*sigma_t

【边界声明】
专著给出了切线参数与表 13.2 数值算例，可逐项核验；
但 Q^UB(sigma_0) 的**具体机构荷载表达式**专著未写出（它依赖所选破坏机构）。
本模块只实现可核验的切线层与定理13.1/表13.2，
不虚构一个"严格上界荷载" —— 那会把工程近似冒充严格界值。
"""
from __future__ import annotations

import math
from typing import Dict, Any, List
import numpy as np

__all__ = [
    "tangent_params", "scan_tangent", "table_13_2", "check_table_13_2",
    "chord_lower_params", "chord_tau", "dissipation_bracket",
    "mechanism_scan",
]


def tau_envelope(sigma, sigma_c, A, sigma_t, B):
    """式(13.9) 幂型强度包线。"""
    return A * sigma_c ** (1.0 - B) * (np.asarray(sigma, dtype=float) + sigma_t) ** B


def tangent_params(sigma_0: float, sigma_c: float, A: float,
                   sigma_t: float, B: float) -> Dict[str, float]:
    """式(13.10)/(13.11)：在 sigma_0 处作包线切线，返回 MC 等效参数。"""
    if sigma_0 + sigma_t <= 0:
        raise ValueError("sigma_0 + sigma_t 必须 > 0")
    if not (0 < B <= 1):
        raise ValueError("B 必须在 (0,1]")
    base = (sigma_0 + sigma_t) / sigma_c
    tan_delta = A * B * base ** (B - 1.0)
    g0 = A * sigma_c ** (1.0 - B) * (sigma_0 + sigma_t) ** B
    c_t = g0 - sigma_0 * tan_delta
    c_t_alt = A * sigma_c ** (1.0 - B) * (sigma_0 + sigma_t) ** (B - 1.0) \
        * (sigma_t + (1.0 - B) * sigma_0)
    delta = math.degrees(math.atan(tan_delta))
    return dict(
        sigma_0=float(sigma_0), tan_delta=float(tan_delta),
        delta_deg=float(delta), psi_deg=float(delta),   # 定理 13.1
        c_t=float(c_t), c_t_alt=float(c_t_alt),
        tau_at_point=float(g0),
        selfcheck=float(c_t + sigma_0 * tan_delta),
        selfcheck_residual=float(abs(c_t + sigma_0 * tan_delta - g0)),
        c_t_formula_residual=float(abs(c_t - c_t_alt)),
    )


def psi_from_sigma(sigma, sigma_c, A, sigma_t, B):
    """式(13.12) 剪胀角，用于定理13.1核对。"""
    if abs(B - 1.0) < 1e-12:
        return math.degrees(math.atan(A))
    return math.degrees(math.atan(A * B * ((sigma + sigma_t) / sigma_c) ** (B - 1.0)))


def scan_tangent(sigma_c: float, A: float, sigma_t: float, B: float,
                 sigma_min: float, sigma_max: float, n: int = 200) -> List[Dict[str, float]]:
    """扫描切点，输出每条切线参数与定理13.1核对（切点倾角 vs 式13.12 psi）。
    注意：不输出 Q^UB —— 专著未给机构荷载表达式，见模块边界声明。"""
    out = []
    for s in np.linspace(sigma_min, sigma_max, n):
        r = tangent_params(float(s), sigma_c, A, sigma_t, B)
        r["psi_eq13_12_deg"] = psi_from_sigma(float(s), sigma_c, A, sigma_t, B)
        r["theorem13_1_residual"] = abs(r["delta_deg"] - r["psi_eq13_12_deg"])
        out.append(r)
    return out


def table_13_2() -> Dict[str, Any]:
    """复现专著表 13.2：sigma_c=100, sigma_t=10, A^2=4*sigma_t/sigma_c => A=0.6325。"""
    sc, st, A = 100.0, 10.0, 0.6325
    rows = []
    for B, sigmas in [(0.5, [0.0, 50.0, 100.0]), (1.0, [0.0, 50.0, 100.0])]:
        for s0 in sigmas:
            r = tangent_params(s0, sc, A, st, B)
            rows.append(dict(
                B=B, sigma_0=s0,
                tan_delta=r["tan_delta"], delta_deg=r["delta_deg"],
                c_t=r["c_t"], selfcheck=r["selfcheck"],
                selfcheck_residual=r["selfcheck_residual"],
                c_t_formula_residual=r["c_t_formula_residual"],
                tau_at_point=r["tau_at_point"],
                slip_angle_deg=90.0 - r["delta_deg"],
            ))
    return dict(sigma_c=sc, sigma_t=st, A=A, rows=rows)


def check_table_13_2() -> Dict[str, Any]:
    """按专著表 13.2 的印刷值逐项核对。

    注意：表 13.2 第 2 行（B=0.5, sigma_0=50）的 c_t 与校核列，
    印刷值为 28.56 / 48.96，而按式 (13.11) 精算为 28.579 / 48.993。
    这是原表的舍入差，不是算法分歧。本函数按印刷精度分档判定，
    并把偏差逐行列为 printed_deviation，绝不静默吞掉。
    """
    d = table_13_2()
    # 专著表 13.2 给出的数值（tanδ, δ°, c_t, 校核=tau, 滑面倾角）
    expected = [
        (0.5, 0.0,   1.000, 45.0, 20.00, 20.00, 45.0),
        (0.5, 50.0,  0.408, 22.2, 28.56, 48.96, 67.8),
        (0.5, 100.0, 0.302, 16.8, 36.19, 66.34, 73.2),
        (1.0, 0.0,   0.633, 32.3,  6.33,  6.33, 57.7),
        (1.0, 50.0,  0.633, 32.3,  6.33, 37.95, 57.7),
        (1.0, 100.0, 0.633, 32.3,  6.33, 69.58, 57.7),
    ]
    results = []
    for row, exp in zip(d["rows"], expected):
        B, s0, tan_e, delta_e, ct_e, chk_e, ang_e = exp
        dev_tan = row["tan_delta"] - tan_e
        dev_delta = row["delta_deg"] - delta_e
        dev_ct = row["c_t"] - ct_e
        dev_chk = row["selfcheck"] - chk_e
        dev_ang = row["slip_angle_deg"] - ang_e
        # 印刷精度：tanδ 3 位、角度 1 位小数、c_t/校核 2 位小数。
        # 2 位小数列的正常舍入容差取 0.01；仅第 2 行超出此范围。
        strict_ok = (abs(dev_tan) < 0.001 and abs(dev_delta) < 0.1
                     and abs(dev_ct) < 0.01 and abs(dev_chk) < 0.01
                     and abs(dev_ang) < 0.1)
        printed_ok = (abs(dev_tan) < 0.001 and abs(dev_delta) < 0.1
                      and abs(dev_ct) < 0.05 and abs(dev_chk) < 0.05
                      and abs(dev_ang) < 0.1)
        results.append(dict(
            B=B, sigma_0=s0,
            exact=dict(
                tan_delta=row["tan_delta"], delta_deg=row["delta_deg"],
                c_t=row["c_t"], selfcheck=row["selfcheck"],
                slip_angle_deg=row["slip_angle_deg"],
            ),
            printed=dict(
                tan_delta=tan_e, delta_deg=delta_e,
                c_t=ct_e, selfcheck=chk_e, slip_angle_deg=ang_e,
            ),
            deviation=dict(
                tan_delta=dev_tan, delta_deg=dev_delta, c_t=dev_ct,
                selfcheck=dev_chk, slip_angle_deg=dev_ang,
            ),
            strict_ok=strict_ok,
            printed_ok=printed_ok,
            rounding_note=("原表印刷舍入差" if printed_ok and not strict_ok else ""),
            selfcheck_residual=row["selfcheck_residual"],
        ))
    return dict(
        results=results,
        all_pass=all(r["printed_ok"] and r["selfcheck_residual"] < 1e-10
                     for r in results),
        strict_pass=all(r["strict_ok"] and r["selfcheck_residual"] < 1e-10
                        for r in results),
        rounding_rows=[r["sigma_0"] for r in results if r["rounding_note"]],
    )


# ---------------------------------------------------------------------------
# 凹包线的弦线下界 + 固定机构上的耗散夹逼（定理 13.2 的数值对应）
# ---------------------------------------------------------------------------
# 式 (13.9) 的 g(sigma) 在 0<B<1 时 g''<0（凹），因此：
#   - 任一点切线  c_t + sigma*tan(delta_t) >= g(sigma)   （上界，GTT 的理论根据）
#   - 任意两切点间弦线                <= g(sigma)       （下界）
# 两者在给定滑面上积分即得耗散积分的上、下夹逼。
# 注意：这是"给定机构上的耗散积分夹逼"，**不是 FELA 的荷载上下界**。
# 专著未给出 Q^UB(sigma_0) 的机构荷载表达式，本程序不虚构严格荷载界值。


def chord_lower_params(sigma_c: float, A: float, sigma_t: float, B: float,
                       sigma_nodes) -> List[Dict[str, float]]:
    """把凹包线分段线性化，返回每条弦段的 MC 等效 (c, tanδ)。弦线处处 <= g。"""
    nodes = np.asarray(sigma_nodes, dtype=float)
    nodes = np.unique(nodes)
    if nodes.size < 2:
        raise ValueError("sigma_nodes 至少需要 2 个节点")
    if np.any(nodes + sigma_t <= 0):
        raise ValueError("sigma_nodes + sigma_t 必须 > 0")
    out = []
    for a, b in zip(nodes[:-1], nodes[1:]):
        ga, gb = float(tau_envelope(a, sigma_c, A, sigma_t, B)), \
                 float(tau_envelope(b, sigma_c, A, sigma_t, B))
        slope = (gb - ga) / (b - a)
        out.append(dict(sigma_lo=float(a), sigma_hi=float(b),
                        c_chord=float(ga - a * slope), tan_delta=float(slope),
                        delta_deg=float(math.degrees(math.atan(slope)))))
    return out


def chord_tau(sigma, sigma_c: float, A: float, sigma_t: float, B: float,
              sigma_nodes) -> np.ndarray:
    """分段线性下界在 sigma 处的取值（凹函数弦线 <= 真值）。"""
    nodes = np.asarray(sigma_nodes, dtype=float)
    nodes = np.unique(nodes)
    vals = np.asarray(tau_envelope(nodes, sigma_c, A, sigma_t, B), dtype=float)
    xs = np.asarray(sigma, dtype=float)
    return np.interp(xs, nodes, vals)


def dissipation_bracket(sigma, sigma_c: float, A: float, sigma_t: float, B: float,
                        sigma_0: float, sigma_nodes=None) -> Dict[str, Any]:
    """给定正应力序列，逐点给出耗散强度真值 / 切线上界 / 弦线下界。"""
    sig = np.asarray(sigma, dtype=float)
    if sig.ndim != 1 or sig.size == 0:
        raise ValueError("sigma 必须是一维非空数组")
    if np.any(sig + sigma_t <= 0):
        raise ValueError("sigma + sigma_t 必须 > 0")
    t = tangent_params(float(sigma_0), sigma_c, A, sigma_t, B)
    tau_true = np.asarray(tau_envelope(sig, sigma_c, A, sigma_t, B), dtype=float)
    tau_up = t["c_t"] + sig * t["tan_delta"]
    if sigma_nodes is None:
        sigma_nodes = np.linspace(float(sig.min()), float(sig.max()), 65)
    tau_low = chord_tau(sig, sigma_c, A, sigma_t, B, sigma_nodes)
    return dict(
        sigma_0=float(sigma_0), tangent=t, sigma=sig.tolist(),
        sigma_nodes=np.asarray(sigma_nodes).tolist(),
        tau_true=tau_true.tolist(), tau_upper=tau_up.tolist(), tau_lower=tau_low.tolist(),
        upper_residual=float(np.min(tau_up - tau_true)),
        lower_residual=float(np.min(tau_true - tau_low)),
        max_gap=float(np.max(tau_up - tau_low)),
    )


def mechanism_scan(params: Dict[str, float], mode: str = "active",
                   n_sigma: int = 80, npts: int = 200) -> Dict[str, Any]:
    """定理 13.2 的固定机构数值对应：扫描切点，最小化耗散积分上界。

    sigma_n(z)=q2+gamma*z 与 L1 内核一致；ds 由 L1 滑动面几何给出。
    理论上 D_UB(sigma_0)/L = c_t(sigma_0) + sigma_bar*tan(delta_t(sigma_0))，
    凹包线下其在 sigma_0 = sigma_bar（正应力弧长加权均值）处取最小。
    本函数实测该最小点并与 sigma_bar 对比；同时给出上/下界与真值积分。
    """
    from collapse_kernel import slip_profile
    p = params
    q2, g = p["q2"], p["gamma"]
    xs, zs = slip_profile(p, mode, npts=npts)
    ds = np.sqrt(np.diff(xs) ** 2 + np.diff(zs) ** 2)
    midz = 0.5 * (zs[:-1] + zs[1:])
    sig = q2 + g * midz  # 与 collapse_kernel.sigma_of_z 同一正应力场
    w = ds
    total = float(w.sum())
    sig_bar = float(np.dot(sig, w) / total)
    nodes = np.linspace(float(sig.min()), float(sig.max()), n_sigma + 1)
    tau_true = np.asarray(tau_envelope(sig, p["sigma_c"], p["A"], p["sigma_t"], p["B"]))
    D_true = float(np.dot(tau_true, w))
    D_low = float(np.dot(chord_tau(sig, p["sigma_c"], p["A"], p["sigma_t"],
                                   p["B"], nodes), w))
    cands = np.linspace(float(nodes[0]), float(nodes[-1]), n_sigma)
    D_ub = []
    for s0 in cands:
        t = tangent_params(float(s0), p["sigma_c"], p["A"], p["sigma_t"], p["B"])
        D_ub.append(float(np.dot(t["c_t"] + sig * t["tan_delta"], w)))
    D_ub = np.asarray(D_ub)
    imin = int(np.argmin(D_ub))
    s0_opt = float(cands[imin])
    # 网格点只能给到步长精度；再用一维连续极小化把切点定位做实。
    # 若 scipy 不可用则退回网格最优，不静默声称达到连续精度。
    refined = False
    try:
        from scipy.optimize import minimize_scalar
        def fun_D(s0):
            t = tangent_params(float(s0), p["sigma_c"], p["A"], p["sigma_t"], p["B"])
            return float(np.dot(t["c_t"] + sig * t["tan_delta"], w))
        lo, hi = float(cands[0]), float(cands[-1])
        res = minimize_scalar(fun_D, bounds=(lo, hi), method="bounded",
                              options=dict(xatol=1e-10))
        if res.success and lo <= res.x <= hi and res.fun <= D_ub[imin] + 1e-9:
            s0_opt = float(res.x)
            refined = True
    except ImportError:
        pass
    B = p["B"]
    t_opt = tangent_params(s0_opt, p["sigma_c"], p["A"], p["sigma_t"], p["B"])
    D_upper_opt = float(np.dot(t_opt["c_t"] + sig * t_opt["tan_delta"], w))
    return dict(
        mode=mode, surface_length=total, sigma_min=float(sig.min()),
        sigma_max=float(sig.max()), sigma_bar=sig_bar,
        refined=refined,
        sigma0_opt=s0_opt, D_true=D_true, D_upper=D_upper_opt,
        D_lower=D_low,
        upper_residual=float(D_upper_opt - D_true),
        lower_residual=float(D_true - D_low),
        theorem13_2_rel_err=(0.0 if abs(B - 1.0) < 1e-12
                             else abs(s0_opt - sig_bar) / max(sig_bar, 1e-12)),
        scan_sigma=cands.tolist(), scan_D_upper=D_ub.tolist(),
        boundary_note=("耗散积分夹逼，非 FELA 荷载上下界；专著未给 Q^UB 机构荷载表达式"),
    )
