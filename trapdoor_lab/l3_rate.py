# -*- coding: utf-8 -*-
"""
l3_rate.py —— L3 速率过程层
============================
严格实现专著第 18 章闭式链：
  式(18.1) I = gamma_dot*d / sqrt(p/rho)          惯性数
  式(18.2) mu(I) = mu_s + (mu_d-mu_s)/(1+I0/I)    颗粒流变
  式(18.4) L/h0 ~ sqrt(h0/d)                      铺展尺度律（诊断）
  式(18.5) Q = C_b*rho_b*sqrt(g)*(B_c-k_d*d)^2.5  Beverloo 出流
           v = 4Q/(pi*B_c^2*rho_b)                 烟囱推进速率
  式(18.6) R(T) = (0.75*B_c^2*h/tan(phi_r))^(1/3)  地表漏斗终径
  式(18.7) gamma*D/(xi*tau0) <= 4                 整柱自稳判据（永不自稳当 tau0=0）
  适用域   h/D <= 1.4（偏差<=5%）；h/D=2 偏 14%；h/D=10 偏 174%
  三门槛   h_cr（几何上限）/ 张开支静力全域无解 / 闭合支静力带（不可混用）
超出适用域时：给出速率链的"量级演示"但 solution_type 标记 degraded，
并明确提示应转用静力/极限分析解。
"""
from __future__ import annotations

import math
from typing import Dict, Any, List
import numpy as np

__all__ = ["DEFAULT_RATE_PARAMS", "solve_rate"]

DEFAULT_RATE_PARAMS: Dict[str, float] = dict(
    h=6.0,          # 覆土 (m)
    B_c=1.0,        # 烟囱通道特征宽度 (m)
    d=0.01,         # 平均粒径 (m)
    rho_b=1600.0,   # 散体堆密度 (kg/m^3)
    C_b=0.60,       # Beverloo 出流系数
    k_d=1.5,        # 空隙环修正系数
    phi_r=35.0,     # 休止角 (deg)
    mu_s=None,      # =tan(phi)，None 时由 A 推定
    mu_d=0.70,      # 惯性分支摩擦系数
    I0=0.30,        # 惯性数尺度
    xi=1.0,         # 侧壁摩擦动员系数
    tau0=0.0,       # 表观黏聚截距 (kPa)
    g=9.81,
)


def mu_of_I(I, mu_s, mu_d, I0):
    I = np.asarray(I, dtype=float)
    out = np.full_like(I, mu_s, dtype=float)
    pos = I > 0
    out[pos] = mu_s + (mu_d - mu_s) / (1.0 + I0 / I[pos])
    return out


def solve_rate(params: Dict[str, float], phi_deg: float = None,
               gamma: float = 18.0, nsteps: int = 120) -> Dict[str, Any]:
    """
    计算速率链。params 至少含 h, B_c, d, rho_b, C_b, k_d, phi_r, mu_d, I0, xi, tau0, g。
    phi_deg 为峰摩擦角（用于 mu_s=tan(phi)），缺省 30°。
    """
    p = dict(DEFAULT_RATE_PARAMS)
    p.update({k: v for k, v in params.items() if v is not None})
    for k in ("h", "B_c", "d", "rho_b", "C_b", "k_d", "phi_r", "mu_d", "I0", "xi", "tau0", "g"):
        if k not in p or p[k] is None:
            raise ValueError("缺速率参数 %s" % k)
        p[k] = float(p[k])

    h, Bc, d = p["h"], p["B_c"], p["d"]
    rho, Cb, kd = p["rho_b"], p["C_b"], p["k_d"]
    g = p["g"]
    if min(h, Bc, d, rho, g) <= 0 or Cb <= 0 or kd < 0:
        raise ValueError("速率参数必须为正（k_d>=0）")
    if Bc - kd * d <= 0:
        raise ValueError("B_c - k_d*d <= 0：Beverloo 空隙环修正后无有效孔径")

    phi = math.radians(30.0 if phi_deg is None else phi_deg)
    mu_s = math.tan(phi)
    mu_d, I0 = p["mu_d"], p["I0"]
    if mu_d < mu_s:
        mu_d = mu_s * 1.05

    # --- Beverloo 出流（式 18.5）---
    Q_mass = Cb * rho * math.sqrt(g) * (Bc - kd * d) ** 2.5     # kg/s
    Q_vol = Q_mass / rho                                         # m^3/s
    v = 4.0 * Q_vol / (math.pi * Bc ** 2)                        # m/s
    T = h / v                                                    # 贯通时间 s
    R = (0.75 * Bc ** 2 * h / math.tan(math.radians(p["phi_r"]))) ** (1.0 / 3.0)

    # --- 惯性数与 mu(I)（式 18.1/18.2）---
    # 烟囱流剪切率 gamma_dot ~ v/d，压力 p ~ rho*g*Bc  => I = v/sqrt(g*Bc)
    I_chim = v / math.sqrt(g * Bc)
    mu_chim = float(mu_of_I(I_chim, mu_s, mu_d, I0))

    # --- 铺展尺度律（式 18.4，诊断）---
    L_over_h0 = math.sqrt(h / d) if h > 0 else 0.0

    # --- 适用域：h/D（D 取烟囱同尺度洞径 B_c）---
    ratio = h / Bc
    if ratio <= 1.4:
        domain, dev = "严格适用域", 0.0
    elif ratio <= 2.0:
        domain, dev = "临界（偏差约 14%）", 14.0
    else:
        domain, dev = "超界（偏差可达 174%）", min(174.0, 20.0 * (ratio - 2.0) + 14.0)
    degraded = ratio > 1.4

    # --- 整柱自稳（式 18.7）---
    xi, tau0 = p["xi"], p["tau0"]
    if tau0 <= 0:
        self_stable, stability_note = False, "tau0=0：无黏性土永不自稳（式18.7）"
    else:
        lhs = gamma * (2.0 * Bc) / (xi * tau0)   # D=2B_c
        self_stable = lhs <= 4.0
        stability_note = "gamma*D/(xi*tau0)=%.3f %s 4" % (lhs, "<=" if self_stable else ">")

    # --- 五阶段链与时间线 ---
    stages: List[Dict[str, Any]] = [
        dict(stage=1, name="静力核查", t=0.0, z=0.0,
             note="几何贯通阈值与静力门槛核对（三门槛不可混用）"),
        dict(stage=2, name="烟囱流接管", t=0.0, z=0.0,
             note="越过阈值后由 μ(I) 颗粒流接管，Beverloo 孔口控制"),
    ]
    times = np.linspace(0.0, T, nsteps)
    timeline = []
    for tt in times:
        zz = min(v * tt, h)
        timeline.append(dict(t=float(tt), z=float(zz), v=v, I=I_chim, mu=mu_chim))
    stages.append(dict(stage=3, name="恒速推进", t=T, z=h, note="v=4Q/(pi*rho*B_c^2)，流量与料柱高度无关"))
    stages.append(dict(stage=4, name="贯通", t=T, z=h, note="累计出流体积 pi*B_c^2*h/4"))
    stages.append(dict(stage=5, name="地表漏斗成形", t=T, z=h,
                       note="体积守恒 R=(0.75*B_c^2*h/tan(phi_r))^(1/3)，与 C_b 无关"))

    warnings: List[str] = []
    if degraded:
        warnings.append("h/D=%.2f 超出速率模型严格适用域(h/D<=1.4)：时间线仅作量级演示，"
                        "应转用静力/极限分析解。" % ratio)
    if self_stable:
        warnings.append("整柱自稳判据满足（%s）：全高贯通式坍塌不发生。" % stability_note)
    if I_chim < 1e-2 or I_chim > 1.0:
        warnings.append("惯性数 I=%.3g 不在 mu(I) 典型区间 1e-2~1，流变律外推需谨慎。" % I_chim)

    return dict(
        solution_type="rate-degraded" if degraded else "rate",
        Q_mass=Q_mass, Q_vol=Q_vol, v=v, T=T, R=R,
        I=I_chim, mu_s=mu_s, mu_d=mu_d, mu_chim=mu_chim,
        L_over_h0=L_over_h0, h_over_D=ratio, domain=domain, domain_dev_pct=dev,
        self_stable=self_stable, stability_note=stability_note,
        stages=stages, timeline=timeline, warnings=warnings,
        params=p,
    )
