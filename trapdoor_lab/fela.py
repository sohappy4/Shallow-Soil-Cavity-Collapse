# -*- coding: utf-8 -*-
"""
fela.py —— 极限分析上下界夹逼（数学规划极限分析 / MP-LA）
=========================================================
把专著 §3.6、§7.7、§20.4.2 的上、下限定理化为您线规划（LP）对偶夹逼，
用 SciPy HiGHS 求解，零新依赖（aarch64 银河麒麟无 cvxpy/SOCP 预编译包）。

理论依据（专著原文，非自创）：
  * 下限定理：满足平衡、边界、不违反屈服准则的静力场，给出的抗力 ≤ 真值。
    下限定理只依赖屈服面，与剪胀角 ψ 无关（§7.7.1）——本模块下限天然 ψ 无关。
  * 上限定理：满足速度边界与关联流动的机构，给出的抗力 ≥ 真值。
    非关联材料按 Davis(1968) 折减强度后回到关联框架，折减上限才严格（§7.3/§7.7）。
  * 二者为对偶数学规划，真值居中，细化时区间单调收窄（§20.4.2、图 3.3）。

实现（对齐专著幂型包线 τ(σ)=A·σ_c·((σ_t+σ)/σ_c)^B，非另立 c/φ）：
  1. 把幂型包线在应力区间上做 K 段线性化。曲线凹(B<1)时：弦线在曲线下方（内接，
     保守，给抗力下界），切线在曲线上方（外接，给抗力上界）；B>1 反之。
  2. 两个 LP 都以「沿滑面抗剪阻力 Ra=∫τ dz」最大化为目标，约束分别为内接/外接
     屈服多面体。由几何包含关系严格有 Ra_LB ≤ Ra_true ≤ Ra_UB；K→∞ 区间收窄。
  3. 主动工况支护压力 q_a = over − Ra/l 是 Ra 的减函数，故夹逼方向反转：
     q_a_LB = over − Ra_UB/l ≤ q_a_true ≤ over − Ra_LB/l = q_a_UB。
  4. Davis 折减单独报告（改材料参数口径，不混入同一问题的数学夹逼区间）。

【诚实边界】
  - 线性化屈服面（多面体）极限分析，不是 SOCP 圆锥规划，不是幂型精确解；
  - 用内接/外接夹持同一幂型包线，K→∞ 才收敛到精确夹逼；
  - 上限的 Davis 折减按专著 §7.3 口径单独报告，不冒充严格荷载上界；
  - N_sat=18.2212 是文献锚点，本模块只做量级对照，不拟合、不推导该常数。
"""
from __future__ import annotations

from typing import Dict, Any, List, Tuple
import math
import numpy as np
from scipy.optimize import linprog

from collapse_kernel import tau_of_z, sigma_of_z, psi0_of_Omega

__all__ = [
    "linearize_strength", "envelope_tau", "davis_reduce", "resistance_bounds",
    "solve_fela", "fela_figure_data", "DEFAULT_FELA_OPTS",
]


DEFAULT_FELA_OPTS: Dict[str, Any] = dict(
    K=12,            # 屈服面线性化段数
    n_points=41,     # 滑面离散点数
    psi_deg=None,    # 剪胀角；给出时报告 Davis 折减上限
)


# ---------------------------------------------------------------------------
# 1. 幂型强度包线的分段线性化
# ---------------------------------------------------------------------------
def linearize_strength(p: Dict[str, float], sigma_a: float, sigma_b: float,
                       K: int = 12, mode: str = "inscribed"):
    """把 τ(σ)=A·σ_c·((σ_t+σ)/σ_c)^B 在 [sigma_a,sigma_b] 上线性化为 K 段。

    返回 (lines, info)，lines 为 (a, b, c) 表示 a·σ + b·τ ≤ c（b 恒为 1）。
      mode="inscribed"     内接：凹函数用弦、凸函数用切点弦中点 → 收缩屈服域；
      mode="circumscribed" 外接：凹函数用切线、凸函数用弦 → 扩张屈服域。
    """
    sc, A, B, st = p["sigma_c"], p["A"], p["B"], p["sigma_t"]
    K = max(2, int(K))
    sig = np.linspace(float(sigma_a), float(sigma_b), K + 1)
    tau = A * sc * ((st + sig) / sc) ** B
    concave = bool(B < 1.0)
    lines: List[Tuple[float, float, float]] = []

    for i in range(K):
        s0, s1 = sig[i], sig[i + 1]
        t0, t1 = tau[i], tau[i + 1]
        slope = (t1 - t0) / (s1 - s0)
        a_chord, c_chord = -slope, t0 - slope * s0
        sm = 0.5 * (s0 + s1)
        tm = A * sc * ((st + sm) / sc) ** B
        der = A * B * ((st + sm) / sc) ** (B - 1.0)
        a_tan, c_tan = -der, tm - der * sm
        if mode == "inscribed":
            chosen = (a_chord, 1.0, c_chord) if concave else (a_tan, 1.0, c_tan)
        else:
            chosen = (a_tan, 1.0, c_tan) if concave else (a_chord, 1.0, c_chord)
        lines.append(chosen)

    info = dict(K=K, sigma_a=float(sigma_a), sigma_b=float(sigma_b),
                concave=concave, mode=mode)
    return lines, info


def envelope_tau(lines, sigma):
    """给定线性化包线族，返回该 σ 处允许的最大 τ（包络）。"""
    sig = np.asarray(sigma, dtype=float)
    env = np.full_like(sig, np.inf)
    for (a, b, c) in lines:
        env = np.minimum(env, (c - a * sig) / b)
    return env


def davis_reduce(p: Dict[str, float], psi_deg: float = None):
    """Davis(1968) 折减：以参考割线摩擦角 φ₀ 与剪胀角 ψ 计算剪阻比。

    对幂型包线等效为 A ← A·ratio（ratio≤1）。参考角用 ψ0_of_Omega(1)，
    即 Ω=1 处的切线角，避免引入专著之外的新参数。
    """
    q = dict(p)
    if psi_deg is None:
        return q, 1.0
    A = float(p["A"])
    phi = float(math.atan(psi0_of_Omega(1.0, A, float(p["B"]))))
    # Davis(1968) 折减只在 psi<=phi 的关联域内有意义；psi>phi 属剪胀超限，
    # 无对应极限定理，这里 clip 到 phi（退化为关联、ratio=1）并标注。
    psi = min(math.radians(float(psi_deg)), phi)
    sinp, cosp = math.sin(psi), math.cos(psi)
    sinf, cosf = math.sin(phi), math.cos(phi)
    # 标准 Davis: tan(phi_red)=sin(phi)*cos(psi)/(1-sin(phi)*sin(psi))
    # 对幂型包线以折减前后割线剪阻比等效折减 A；psi 小折减多，psi→phi 时→1。
    den = max(1.0 - sinf * sinp, 1e-6)
    tan_red = sinf * cosp / den
    ratio = tan_red / max(math.tan(phi), 1e-9)
    ratio = float(np.clip(ratio, 0.2, 1.0))
    # 关联域（psi>=phi）直接返回原参数，避免 A*1.0 的浮点尾差污染对照
    if ratio >= 1.0 - 1e-12:
        return dict(p), 1.0
    q["A"] = A * ratio
    return q, float(ratio)


# ---------------------------------------------------------------------------
# 2. 抗力上下界 LP
# ---------------------------------------------------------------------------
def _resistance_lp(p, sigma, w, lines):
    """在给定线性化屈服域内最大化 Ra=Σ w·τ，返回 Ra 与 τ 分布。"""
    n = len(sigma)
    c_obj = -np.asarray(w, dtype=float)          # 最小化 -Ra → 最大化 Ra
    A_ub, b_ub = [], []
    for i in range(n):
        for (a, b, c) in lines:
            row = np.zeros(n)
            row[i] = b
            A_ub.append(row)
            b_ub.append(c - a * float(sigma[i]))
    bounds = [(0.0, None)] * n
    res = linprog(c_obj, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  bounds=bounds, method="highs")
    if not res.success:
        return None, None, str(res.message)
    return float(np.dot(w, res.x)), res.x, None


def resistance_bounds(p: Dict[str, float], opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """抗力 Ra 的内接/外接 LP 夹逼 + 主动支护压力 q_a 的反序夹逼。"""
    o = dict(DEFAULT_FELA_OPTS)
    if opts:
        o.update({k: v for k, v in opts.items() if v is not None})
    n = max(6, int(o["n_points"]))
    K = max(2, int(o["K"]))
    h = float(p["h"]); gamma = float(p["gamma"]); q2 = float(p.get("q2", 0.0))
    l = float(p["l"])
    z = np.linspace(0.0, h, n)
    sig = sigma_of_z(z, q2, gamma)
    w = np.full(n, h / (n - 1)); w[0] *= 0.5; w[-1] *= 0.5

    sa, sb = max(1e-6, q2), q2 + gamma * h
    lines_in, info_in = linearize_strength(p, sa, sb, K, "inscribed")
    lines_out, info_out = linearize_strength(p, sa, sb, K, "circumscribed")

    Ra_lb, tau_lb, msg_lb = _resistance_lp(p, sig, w, lines_in)
    Ra_ub, tau_ub, msg_ub = _resistance_lp(p, sig, w, lines_out)
    if Ra_lb is None or Ra_ub is None:
        return dict(success=False, message="; ".join(x for x in (msg_lb, msg_ub) if x))

    over = q2 + gamma * h
    # 主动支护压力 q_a = over - Ra/l（Ra 越大，所需支护越小）
    qa_lb = over - Ra_ub / l       # 用最大抗力 → 最小支护 → 支护压力下界
    qa_ub = over - Ra_lb / l       # 用最小抗力 → 最大支护 → 支护压力上界

    # L1 精确 Ra（同参数闭式）作为区间核验参考，不用于构造界
    from collapse_kernel import R_scale
    Ra_exact = R_scale(p["sigma_c"], p["sigma_t"], q2, gamma, h, p["A"], p["B"])
    in_interval = bool(Ra_lb - 1e-6 <= Ra_exact <= Ra_ub + 1e-6)

    # 适用域门控（诚实边界）：专著幂型包线 tau∝sigma^B 在深埋随 sigma 无限增长，
    # Ra∝h^(1+B) 快于自重 gamma*h，q_a 变负且无饱和；真实土体受临界状态封顶，
    # 文献 N_sat 正来自该饱和。本 LP 只线性化同一幂型包线、不含临界状态帽，
    # 故深埋外推无物理意义。按 h/D 门控并显式退化标记（同 L3 的 rate-degraded）。
    D = 2.0 * l
    h_over_D = h / D
    applicable = bool(h_over_D <= 1.4 and qa_lb > 0.0)
    regime = ("shallow-applicable" if applicable
              else "out-of-domain-deep (no critical-state cap; bracket not physical)")

    result = dict(
        success=True,
        applicable=applicable,
        h_over_D=float(h_over_D),
        regime=regime,
        Ra_LB=float(Ra_lb), Ra_UB=float(Ra_ub), Ra_exact=float(Ra_exact),
        Ra_interval_ok=in_interval,
        q_active_LB=float(qa_lb), q_active_UB=float(qa_ub),
        q_active_mid=float(0.5 * (qa_lb + qa_ub)),
        q_active_width=float(qa_ub - qa_lb),
        q_active_rel_width=float((qa_ub - qa_lb) / max(abs(0.5*(qa_lb+qa_ub)), 1e-9)),
        z_m=z.tolist(), sigma_kPa=sig.tolist(),
        tau_LB_kPa=(tau_lb.tolist() if tau_lb is not None else None),
        tau_UB_kPa=(tau_ub.tolist() if tau_ub is not None else None),
        envelope_in_kPa=envelope_tau(lines_in, sig).tolist(),
        envelope_out_kPa=envelope_tau(lines_out, sig).tolist(),
        tau_exact_kPa=tau_of_z(z, p["sigma_c"], p["sigma_t"], q2, gamma,
                               p["A"], p["B"]).tolist(),
        n_points=n, K=K,
        lin_info=dict(inscribed=info_in, circumscribed=info_out),
    )

    # Davis 折减：单独报告折减后的运动上限（非关联工程口径）
    if o.get("psi_deg") is not None:
        p_red, ratio = davis_reduce(p, float(o["psi_deg"]))
        lines_red, _ = linearize_strength(p_red, sa, sb, K, "circumscribed")
        Ra_red, tau_red, msg_red = _resistance_lp(p_red, sig, w, lines_red)
        result["davis"] = dict(
            psi_deg=float(o["psi_deg"]), ratio=ratio,
            A_reduced=float(p_red["A"]),
            Ra_UB_reduced=(float(Ra_red) if Ra_red is not None else None),
            q_active_UB_reduced=(over - Ra_red / l if Ra_red is not None else None),
            note="按专著 §7.3 Davis 折减的关联上限口径；仅工程估计，不混入数学夹逼区间",
        )
    return result


def solve_fela(params: Dict[str, float], opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """FELA 夹逼统一入口。"""
    o = dict(DEFAULT_FELA_OPTS)
    if opts:
        o.update({k: v for k, v in opts.items() if v is not None})
    rb = resistance_bounds(params, o)
    ok = bool(rb.get("success")) and bool(rb.get("Ra_interval_ok")) \
         and bool(rb.get("applicable"))
    return dict(
        solution_type="fela-lp",
        meta=dict(
            method="linearized power-law yield surface + paired LP (inscribed/circumscribed)",
            source="external-mature-LP-limit-analysis (Sloan-type); user authorization 2026-11-24",
            boundary="Linearized polyhedral yield surface, NOT SOCP, NOT exact power-law solution; "
                     "static lower bound is psi-independent; upper bound Davis reduction reported "
                     "separately; bracket from monograph §3.6/§7.7/§20.4.2; "
                     "N_sat=18.2212 is literature anchor only, not derived here; "
                     "out-of-domain deep cases are flagged not physical (no critical-state cap)",
            K=int(rb.get("K", o["K"])),
            n_points=int(rb.get("n_points", o["n_points"])),
            psi_deg=o.get("psi_deg"),
            applicable=rb.get("applicable"),
            regime=rb.get("regime"),
        ),
        resistance=rb,
        bracket=dict(
            ok=ok,
            Ra_LB_kPa=rb.get("Ra_LB"), Ra_exact_kPa=rb.get("Ra_exact"),
            Ra_UB_kPa=rb.get("Ra_UB"),
            q_LB_kPa=rb.get("q_active_LB"), q_UB_kPa=rb.get("q_active_UB"),
            q_mid_kPa=rb.get("q_active_mid"),
        ) if rb.get("success") else dict(ok=False, message=rb.get("message")),
        params=dict(l=float(params["l"]), h=float(params["h"]),
                    gamma=float(params["gamma"]), q2=float(params.get("q2", 0.0)),
                    sigma_c=float(params["sigma_c"]), A=float(params["A"]),
                    B=float(params["B"]), sigma_t=float(params["sigma_t"])),
    )


# ---------------------------------------------------------------------------
# 5. 夹逼带绘图数据（供 Web/报告使用）
# ---------------------------------------------------------------------------
def fela_figure_data(params: Dict[str, float], opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """生成夹逼带图数据：q_LB/q_exact/q_UB 随线性化段数 K 的收窄。

    返回纯数组字典，绘图由 app.py 完成，保持内核无 matplotlib 依赖。
    """
    o = dict(DEFAULT_FELA_OPTS)
    if opts:
        o.update({k: v for k, v in opts.items() if v is not None})
    Ks = [2, 3, 4, 6, 8, 12, 16, 24, 32]
    q_lb, q_ub, q_exact, widths = [], [], [], []
    applicable_flags = []
    for K in Ks:
        r = resistance_bounds(params, dict(o, K=K))
        if not r.get("success"):
            q_lb.append(None); q_ub.append(None)
            q_exact.append(None); widths.append(None)
            applicable_flags.append(False)
            continue
        q_lb.append(r["q_active_LB"])
        q_ub.append(r["q_active_UB"])
        q_exact.append(r["Ra_exact"])
        widths.append(r["q_active_width"])
        applicable_flags.append(bool(r.get("applicable")))

    from collapse_kernel import limit_loads
    q1a = limit_loads(
        float(params["sigma_c"]), float(params["sigma_t"]),
        float(params.get("q2", 0.0)), float(params["gamma"]),
        float(params["l"]), float(params["h"]),
        float(params["A"]), float(params["B"]))[0]

    return dict(
        K=Ks,
        q_LB=q_lb, q_UB=q_ub,
        Ra_exact=q_exact,
        width=widths,
        applicable=applicable_flags,
        q_l1_closed_form=float(q1a),
        note=("主动支护压力 q=over-Ra/l 是抗力 Ra 的减函数："
              "内接多面体给 Ra 下界→q 上界，外接给 Ra 上界→q 下界；"
              "K 增大时夹逼带单调收窄并包住 L1 闭式解"),
    )
