# -*- coding: utf-8 -*-
"""三维轴对称锥形滑体极限平衡求解器（方案 M8 配套，横向引用）。

方法来源与边界：
  横向引用经典三维土拱/活门极限平衡（Terzaghi 活门、Meyerhof 等）的
  轴对称滑体方法；滑体三形态（漏斗 funnel / 瓮形 urn / 井筒 shaft）
  对齐专著的三形态几何命名，但**把二维剖面推广为三维轴对称**属本程序的
  扩展，不是专著原有公式。

  这是**简化轴对称极限平衡**（刚体滑体 + 侧摩阻力积分），不是完整三维
  FEM/FELA，不是 SOCP；与 L1/L2（平面应变）只作**量级对照**，
  不做同一性比较（维度不同、应力路径不同）。

几何约定：
  z 自活门底面向上（0 → h）；活门等效圆半径 r0（由洞半宽 l 换算）；
  r(z) 为滑体半径剖面：
    funnel: r = r0 + z·tanθ                （向上线性扩大，漏斗）
    urn   : r = r0 + 4·a·z(h−z)/h²          （中部最宽，瓮形）
    shaft : r = r0                          （等径，井筒）

极限平衡：
  自重   W = γ∫πr²dz
  侧阻   S = ∫(c + K·γ(h−z)·tanφ)·2πr·sqrt(1+(dr/dz)²) dz
  支护   q = (W − S)/(π r0²)   （q≤0 表示自稳）
  临界埋深 h_cr 由 q(h)=0 二分求得。
"""
from __future__ import annotations

import math
from typing import Any, Callable, Dict, List

import numpy as np

__all__ = [
    "DEFAULT_CONE_PARAMS", "MORPHOLOGIES", "radius_profile",
    "solve_cone", "critical_depth", "selfcheck_cone",
]

DEFAULT_CONE_PARAMS: Dict[str, float] = dict(
    l=2.5,            # 洞半宽 (m)
    h=3.0,            # 覆土 (m)
    gamma=18.0,       # 重度 (kN/m3)
    phi_deg=25.0,     # 内摩擦角
    c=10.0,           # 黏聚力 (kPa)
    theta_deg=20.0,   # 漏斗扩展角（滑面与竖直夹角）
    urn_a=0.8,        # 瓮形中部额外半径系数（×l）
    K=0.6,            # 侧向土压力系数
    r0_factor=1.0,    # 等效圆半径 = r0_factor × l
)

MORPHOLOGIES = ("funnel", "urn", "shaft")


def radius_profile(morph: str, z: np.ndarray, h: float, r0: float,
                   p: Dict[str, float]) -> np.ndarray:
    """滑体半径 r(z) 及其导数（数值中心差分）。"""
    if morph == "funnel":
        return r0 + z * math.tan(math.radians(float(p["theta_deg"])))
    if morph == "urn":
        a = float(p["urn_a"]) * r0
        return r0 + 4.0 * a * z * (h - z) / (h * h)
    if morph == "shaft":
        return np.full_like(z, r0)
    raise ValueError("未知形态：%r（应为 funnel/urn/shaft）" % morph)


def _integrate(morph: str, h: float, r0: float, p: Dict[str, float],
               n: int = 401):
    """沿 z 数值积分自重 W 与侧摩阻 S，返回 (W, S, r_prof, z_prof)。"""
    z = np.linspace(0.0, h, n)
    r = radius_profile(morph, z, h, r0, p)
    dz = z[1] - z[0]
    # dr/dz 中心差分
    drdz = np.gradient(r, dz)
    phi = math.radians(float(p["phi_deg"]))
    gamma = float(p["gamma"])
    c = float(p["c"])
    K = float(p["K"])
    # 侧摩阻：有效正应力 ≈ K·γ·(h−z)（深度线性），加上黏聚力
    tau = c + K * gamma * np.maximum(h - z, 0.0) * math.tan(phi)
    dA_slope = np.sqrt(1.0 + drdz ** 2)      # 斜滑面面积放大
    W = gamma * np.trapz(np.pi * r ** 2, z)
    S = np.trapz(tau * 2.0 * np.pi * r * dA_slope, z)
    return float(W), float(S), r, z


def solve_cone(params: Dict[str, float] = None,
               morph: str = "funnel",
               opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """给定参数与形态，计算三维锥形滑体的支护压力与平衡量。

    返回 dict：meta / result / profile / comparison / boundary。
    """
    p = dict(DEFAULT_CONE_PARAMS)
    if params:
        p.update({k: v for k, v in params.items() if k in p})
    if morph not in MORPHOLOGIES:
        raise ValueError("形态必须是 %s 之一" % (MORPHOLOGIES,))
    o = dict(n=401)
    if opts:
        o.update(opts)

    h, l = float(p["h"]), float(p["l"])
    r0 = float(p["r0_factor"]) * l
    area0 = math.pi * r0 ** 2
    W, S, r_prof, z_prof = _integrate(morph, h, r0, p, int(o["n"]))
    q = (W - S) / area0
    q_simplified = q if q > 0.0 else 0.0     # 自稳时支护需求取 0

    # 与 L1 二维闭式解量级对照（维度不同，只标注量级，不判同一性）
    try:
        from collapse_kernel import limit_loads
        q_l1 = limit_loads(
            float(p.get("sigma_c", 100.0)), float(p.get("sigma_t", 20.0)),
            float(p.get("q2", 0.0)), float(p["gamma"]), l, h,
            float(p.get("A", 0.5)), float(p.get("B", 0.5)))[0]
    except Exception:
        q_l1 = None

    return dict(
        meta=dict(morph=morph, n_integration=int(o["n"]),
                  r0=float(r0), h=float(h), gamma=float(p["gamma"]),
                  phi_deg=float(p["phi_deg"]), c=float(p["c"])),
        result=dict(
            weight_kN=float(W),
            side_resistance_kN=float(S),
            support_pressure_kPa=float(q),
            support_demand_kPa=float(q_simplified),
            self_stable=bool(q <= 0.0),
            mobilized_fraction=float(S / W) if W > 0 else 0.0,
            r_top=float(r_prof[-1]),
            r_max=float(np.max(r_prof)),
        ),
        profile=dict(
            z=[float(v) for v in z_prof],
            r=[float(v) for v in r_prof],
        ),
        comparison=dict(
            l1_2d_pressure_kPa=(float(q_l1) if q_l1 is not None else None),
            note=("三维轴对称极限平衡 vs 二维平面应变闭式解：维度与应力路径不同，"
                  "只作量级对照，不判同一性；q 差异属预期。"),
        ),
        boundary=("简化三维轴对称刚体滑体极限平衡（侧摩阻积分为工程近似），"
                  "横向引用经典土拱方法并按专著三形态命名扩展，非专著原有公式；"
                  "不是三维 FEM/FELA/SOCP；参数 c/φ/K 为独立输入，"
                  "与 L1 的幂型 Hoek–Brown 参数体系不同族，不做同一性对照。"),
    )


def critical_depth(params: Dict[str, float] = None, morph: str = "funnel",
                   h_max: float = 30.0, tol: float = 1e-4) -> Dict[str, Any]:
    """二分求 q(h)=0 的临界埋深 h_cr（自稳埋深）；无根则报告。"""
    p = dict(DEFAULT_CONE_PARAMS)
    if params:
        p.update({k: v for k, v in params.items() if k in p})

    def q_at(h):
        return solve_cone(dict(p, h=h), morph, dict(n=121))["result"]["support_pressure_kPa"]

    q_lo, q_hi = q_at(1e-3), q_at(h_max)
    if q_lo * q_hi > 0:
        return dict(has_root=False, h_cr=None, q_at_low=float(q_lo),
                    q_at_high=float(q_hi), note="搜索区间内 q 不变号，无临界埋深")
    lo, hi = 1e-3, h_max
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if q_at(mid) * q_lo > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    h_cr = 0.5 * (lo + hi)
    return dict(has_root=True, h_cr=float(h_cr),
                q_at_low=float(q_lo), q_at_high=float(q_hi),
                morph=morph, tolerance=float(tol))


def selfcheck_cone() -> Dict[str, Any]:
    """自检：三形态半径剖面、量纲、自稳趋势、临界埋深、边界声明。"""
    checks: Dict[str, Any] = {}
    p = dict(DEFAULT_CONE_PARAMS)
    # 形态剖面单调/有界
    z = np.linspace(0, p["h"], 50)
    r_fun = radius_profile("funnel", z, p["h"], p["l"], p)
    r_urn = radius_profile("urn", z, p["h"], p["l"], p)
    r_shaft = radius_profile("shaft", z, p["h"], p["l"], p)
    checks["funnel_increasing"] = bool(np.all(np.diff(r_fun) > 0))
    checks["urn_peaks_in_middle"] = bool(r_urn[len(r_urn) // 2] > r_urn[0] and
                                         r_urn[len(r_urn) // 2] > r_urn[-1])
    checks["shaft_constant"] = bool(np.allclose(r_shaft, p["l"]))
    checks["all_positive_radius"] = bool(np.all(r_fun > 0) and np.all(r_urn > 0))

    r = solve_cone(p, "funnel")
    checks["weight_positive"] = bool(r["result"]["weight_kN"] > 0)
    checks["resistance_positive"] = bool(r["result"]["side_resistance_kN"] > 0)
    checks["pressure_finite"] = bool(np.isfinite(r["result"]["support_pressure_kPa"]))
    checks["mobilized_in_range"] = bool(0.0 <= r["result"]["mobilized_fraction"] <= 2.0)

    # 埋深增大→自重增长快于侧阻→支护压力上升（对 funnel）
    q1 = solve_cone(dict(p, h=2.0), "funnel")["result"]["support_pressure_kPa"]
    q2 = solve_cone(dict(p, h=8.0), "funnel")["result"]["support_pressure_kPa"]
    checks["pressure_increases_with_depth"] = bool(q2 > q1)

    # 黏聚力/摩擦增大→侧阻增大→支护压力下降
    rq = solve_cone(dict(p, c=40.0), "funnel")["result"]["support_pressure_kPa"]
    checks["cohesion_reduces_pressure"] = bool(rq < r["result"]["support_pressure_kPa"])

    cr = critical_depth(p, "funnel")
    checks["critical_depth_schema"] = bool("has_root" in cr and
                                           (cr["has_root"] == (cr["h_cr"] is not None)))
    checks["boundary_honest"] = ("横向引用" in r["boundary"] and "非专著原有公式" in r["boundary"])
    checks["all_pass"] = all(v for k, v in checks.items() if k != "all_pass")
    return checks


if __name__ == "__main__":
    import json
    print(json.dumps(selfcheck_cone(), indent=2, ensure_ascii=False))
