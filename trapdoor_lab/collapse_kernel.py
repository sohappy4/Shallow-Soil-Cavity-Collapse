# -*- coding: utf-8 -*-
"""
collapse_kernel.py —— 浅埋土洞活门试验 L1 准静态闭式解内核
=================================================================
这是纯计算内核：不画图、不碰界面、不读全局变量。
对外只暴露 solve(params, mode) -> Result。

符号约定（与专著、论文 JRMGE 统一，V2 修订版）：
  l      洞半宽 (m)，D = 2l 为洞全宽
  n      埋深比 = h/(2l) = h/D     —— 注意不是 h/l
  Omega  脆性比 = (sigma_t+q2)/sigma_c
  Lambda 覆强比 = gamma*l/sigma_c
  t      出露坡度参数 = (sigma_t+q2)/(B*S)
  psi    剪胀角 (rad)，psi(z)=arctan[A*B*((sigma_t+q2+gamma*z)/sigma_c)^(B-1)]
  q1a/q1b 主动支护压力 / 被动顶起压力（输出，不是输入）
  B -> 1 走线性极限分支，避免 strength_S 除零

公式来源：
  Eq.6/7 S；Eq.21 t；Eq.28/29 sigma/tau；Eq.30 psi；Eq.32 h_cr；
  Eq.35a/35b q1a/q1b；Eq.36 R_a；Eq.38 psi0(Omega)；Eq.44 Ka；Eq.46 A_psi；
  Eq.42 h_cr=l*cot(phi)；Eq.48 rho=1∓n*tan(phi)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
import math
import numpy as np

__all__ = [
    "ParamsError", "Result", "validate_params", "solve",
    "figure_data", "reproduce_tables", "DEFAULT_PARAMS", "PRESETS",
]

# ---------------------------------------------------------------------------
# 0. 默认参数与预设土性
# ---------------------------------------------------------------------------
DEFAULT_PARAMS: Dict[str, float] = dict(
    sigma_c=100.0,  # 强度参考应力 (kPa)
    A=0.6,          # 幂型强度准则系数 (=tan(phi) 当 B->1)
    B=0.6,          # 幂型强度准则指数 0<B<=1
    sigma_t=5.0,    # 抗拉强度 (kPa)
    gamma=18.0,     # 土体重度 (kN/m^3)
    q2=10.0,        # 地表均布超载 (kPa)
    l=2.5,          # 洞半宽 (m)
    h=3.0,          # 覆盖层厚度 (m)
)

PRESETS: Dict[str, Dict[str, float]] = {
    "密实砂": dict(DEFAULT_PARAMS, B=0.95, A=0.72, sigma_t=1.0, gamma=17.5),
    "松散砂": dict(DEFAULT_PARAMS, B=0.90, A=0.58, sigma_t=0.5, gamma=16.5),
    "软黏土": dict(DEFAULT_PARAMS, B=0.20, A=0.35, sigma_t=12.0, gamma=17.0),
    "结构性土": dict(DEFAULT_PARAMS, B=0.55, A=0.50, sigma_t=18.0, gamma=18.5),
}

LINEAR_EPS = 1e-9          # |B-1| 小于此值即走线性极限分支
ZONE_OMEGA = [0.0196, 0.184, 0.963]   # 论文 Fig.5 / Table 2 分区边界


class ParamsError(ValueError):
    """参数不合法。注意：这是对物理量合法性的校验，不是理论适用域拦截。"""
    pass


@dataclass
class Result:
    """统一结果对象。所有字段均可 JSON 序列化。"""
    params: Dict[str, float]
    mode: str
    solution_type: str                 # "power" | "linear-limit"
    # 派生量
    S: Optional[float]
    t: Optional[float]
    Omega: float
    Lambda: float
    n: float
    eta: float                         # 兼容别名，= n
    psi0_deg: float
    psi_h_deg: float
    h_cr: float
    # 荷载
    q1a: float
    q1b: float
    Ra: float
    over: float
    mirror_error: float
    # 判定
    is_global_collapse: bool
    is_global_applicable: bool         # 被动工况无 h_cr 门控，此值为 False
    judgement_text: str
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        d = self.__dict__.copy()
        return d


# ---------------------------------------------------------------------------
# 1. 输入校验（合法性，不是适用域）
# ---------------------------------------------------------------------------
def validate_params(p: Dict[str, float]) -> None:
    required = ("sigma_c", "A", "B", "sigma_t", "gamma", "q2", "l", "h")
    for k in required:
        if k not in p:
            raise ParamsError("缺少参数: %s" % k)
        try:
            v = float(p[k])
        except (TypeError, ValueError):
            raise ParamsError("参数 %s 不是数字: %r" % (k, p[k]))
        if not math.isfinite(v):
            raise ParamsError("参数 %s 必须是有限数" % k)
    if p["sigma_c"] <= 0:
        raise ParamsError("sigma_c 必须 > 0")
    if p["A"] <= 0:
        raise ParamsError("A 必须 > 0（A = tan(phi) 当 B->1）")
    if not (0 < p["B"] <= 1):
        raise ParamsError("B 必须满足 0 < B <= 1")
    if p["sigma_t"] < 0:
        raise ParamsError("sigma_t 不能为负")
    if p["gamma"] <= 0:
        raise ParamsError("gamma 必须 > 0")
    if p["q2"] < 0:
        raise ParamsError("q2 不能为负")
    if p["l"] <= 0:
        raise ParamsError("l（洞半宽）必须 > 0")
    if p["h"] <= 0:
        raise ParamsError("h（覆土厚度）必须 > 0")


# ---------------------------------------------------------------------------
# 2. 基础公式（Eq.6/7, 21, 28~30, 36, 32, 35a/35b）
# ---------------------------------------------------------------------------
def strength_S(sigma_c: float, A: float, B: float) -> Optional[float]:
    if abs(B - 1.0) < LINEAR_EPS:
        return None                      # B=1 无有限 S，走线性分支
    return sigma_c * A ** (1.0 / (1.0 - B)) * B ** (B / (1.0 - B))


def t_param(sigma_c, sigma_t, q2, A, B, S=None):
    if S is None:
        S = strength_S(sigma_c, A, B)
    if S is None:
        return None
    return (sigma_t + q2) / (B * S)


def sigma_of_z(z, q2, gamma):
    return q2 + gamma * np.asarray(z, dtype=float)


def tau_of_z(z, sigma_c, sigma_t, q2, gamma, A, B):
    return A * sigma_c * ((sigma_t + q2 + gamma * np.asarray(z, dtype=float)) / sigma_c) ** B


def psi_of_z(z, sigma_c, sigma_t, q2, gamma, A, B):
    if abs(B - 1.0) < LINEAR_EPS:
        return np.full_like(np.asarray(z, dtype=float), math.atan(A))
    val = A * B * ((sigma_t + q2 + gamma * np.asarray(z, dtype=float)) / sigma_c) ** (B - 1.0)
    return np.arctan(val)


def R_scale(sigma_c, sigma_t, q2, gamma, h, A, B):
    """Eq.36 深度积分量 R_a。B=1 时用解析极限式。"""
    top = (sigma_t + q2 + gamma * h) ** (1.0 + B)
    bot = (sigma_t + q2) ** (1.0 + B)
    if abs(B - 1.0) < LINEAR_EPS:
        # lim_{B->1}: A*[(st+q2+g h)^2-(st+q2)^2]/(2g)
        return A * ((sigma_t + q2 + gamma * h) ** 2 - (sigma_t + q2) ** 2) / (2.0 * gamma)
    return A * sigma_c ** (1.0 - B) * (top - bot) / ((1.0 + B) * gamma)


def limit_loads(sigma_c, sigma_t, q2, gamma, l, h, A, B):
    Ra = R_scale(sigma_c, sigma_t, q2, gamma, h, A, B)
    over = q2 + gamma * h
    return over - Ra / l, over + Ra / l, Ra, over


def h_critical(sigma_c, sigma_t, q2, gamma, l, A, B):
    """Eq.32；B->1 极限为 l/A（数值已验证，含 sigma_t/q2 亦收敛）。"""
    if abs(B - 1.0) < LINEAR_EPS:
        return l / A
    S = strength_S(sigma_c, A, B)
    t = t_param(sigma_c, sigma_t, q2, A, B, S)
    return (B * S / gamma) * ((gamma * l / S + t ** B) ** (1.0 / B) - t)


def psi0_of_Omega(Omega, A, B):
    if abs(B - 1.0) < LINEAR_EPS:
        return np.full_like(np.asarray(Omega, dtype=float), math.atan(A))
    return np.arctan(A * B * np.asarray(Omega, dtype=float) ** (B - 1.0))


def Ka_local(sigma, tau, psi):
    s, c = np.sin(psi), np.cos(psi)
    return (sigma * c - tau * (1.0 - s)) / (sigma * c + tau * (1.0 + s))


def A_psi_from_psi0(psi0, Omega, B):
    return np.tan(psi0) / (B * np.asarray(Omega, dtype=float) ** (B - 1.0))


# ---------------------------------------------------------------------------
# 3. 滑动面坐标（主动瓮形 / 被动漏斗形）
# ---------------------------------------------------------------------------
def slip_profile(p: Dict[str, float], mode: str = "active", npts: int = 400):
    sc, A, B, st = p["sigma_c"], p["A"], p["B"], p["sigma_t"]
    q2, g, l, h = p["q2"], p["gamma"], p["l"], p["h"]
    if abs(B - 1.0) < LINEAR_EPS:
        # 线性极限：滑动面为直线 f' = -cot(phi) = -1/A（主动）
        phi = math.atan(A)
        dx = h / math.tan(phi) if phi > 0 else h * 1e6
        if mode == "active":
            x2 = max(l - dx, 0.0)
            x = np.linspace(x2, l, npts)
        else:
            x2 = l + dx
            # 与幂型分支保持同一约定：被动支 x 从 x2 递减到 l，
            # F=h-z 从 h 递减到 0。这样关于洞中线 x=l 镜像后恰与主动支重合。
            # 旧写法 np.linspace(l, x2) 方向相反，导致 B=1 工况镜像校验失败。
            x = np.linspace(x2, l, npts)
        z = np.linspace(0.0, h, npts)
        return np.clip(x, 0, None), h - z
    S = strength_S(sc, A, B)
    t = t_param(sc, st, q2, A, B, S)
    z = np.linspace(0.0, h, npts)
    base = t + g * z / (B * S)
    varpi = S / g * (base ** B - t ** B)
    if mode == "active":
        delta = S / g * (base[-1] ** B - t ** B)
        x2 = l - delta
        x = x2 + varpi
    else:
        x2 = l + S / g * (base[-1] ** B - t ** B)
        x = x2 - varpi
    return np.clip(x, 0, None), h - z


# ---------------------------------------------------------------------------
# 4. 统一求解接口
# ---------------------------------------------------------------------------
def solve(params: Dict[str, float], mode: str = "active") -> Result:
    """
    求解一个工况。
    params: 参数字典；mode: 'active'(主动塌落) 或 'passive'(被动顶起)
    返回 Result。任何不合法参数抛 ParamsError。
    """
    validate_params(params)
    if mode not in ("active", "passive"):
        raise ParamsError("mode 只能是 'active' 或 'passive'")

    p = {k: float(v) for k, v in params.items()}
    sc, A, B, st = p["sigma_c"], p["A"], p["B"], p["sigma_t"]
    g, q2, l, h = p["gamma"], p["q2"], p["l"], p["h"]

    linear = abs(B - 1.0) < LINEAR_EPS
    sol_type = "linear-limit" if linear else "power"

    S = None if linear else strength_S(sc, A, B)
    t = None if linear else t_param(sc, st, q2, A, B, S)

    Omega = (st + q2) / sc
    Lambda = g * l / sc
    n = h / (2.0 * l)                    # === 修正点：n=h/(2l)，不是 h/l ===

    psi_surf = float(np.degrees(float(psi_of_z(0.0, sc, st, q2, g, A, B))))
    psi_bot = float(np.degrees(float(psi_of_z(h, sc, st, q2, g, A, B))))

    h_cr = h_critical(sc, st, q2, g, l, A, B)
    q1a, q1b, Ra, over = limit_loads(sc, st, q2, g, l, h, A, B)
    mirror_error = abs((q1b - q1a) - 2.0 * Ra / l)

    warnings: List[str] = []
    if linear:
        warnings.append("B=1：已切换到线性极限分支（h_cr=l/A，q1a/q1b 取 B→1 解析极限），非幂型闭式解。")
    if mirror_error > 1e-8:
        warnings.append("镜面对称校验偏差 %.3e，超出 1e-8。" % mirror_error)
    if not math.isfinite(h_cr) or h_cr <= 0:
        warnings.append("h_cr 计算异常，请检查参数。")

    # === 修正点：h_cr 门控只对主动工况适用（专著：被动出露点随覆土外移，无 h_cr）===
    if mode == "active":
        is_global = h <= h_cr
        applicable = True
        if is_global:
            judgement = ("h=%.3f m <= h_cr=%.3f m：浅埋，滑动面贯通地表，发生整体塌落。" % (h, h_cr))
        else:
            judgement = ("h=%.3f m > h_cr=%.3f m：该覆土厚度下不会发生整体贯通塌落，"
                         "仅给出静力/极限解（局部坍塌数据为二期功能）。" % (h, h_cr))
    else:
        is_global = True                  # 被动不设 h_cr 门控
        applicable = False
        judgement = ("被动顶起工况：专著不设 h_cr 门控，出露点随覆土外移；"
                     "输出为被动极限顶起压力与出露形态。")

    return Result(
        params=p, mode=mode, solution_type=sol_type,
        S=S, t=t, Omega=Omega, Lambda=Lambda, n=n, eta=n,
        psi0_deg=psi_surf, psi_h_deg=psi_bot, h_cr=h_cr,
        q1a=q1a, q1b=q1b, Ra=Ra, over=over, mirror_error=mirror_error,
        is_global_collapse=is_global, is_global_applicable=applicable,
        judgement_text=judgement, warnings=warnings,
    )


# ---------------------------------------------------------------------------
# 5. 六图的数据（纯数组，交给表现层画图）
# ---------------------------------------------------------------------------
def figure_data(p: Dict[str, float], mode: str = "active") -> Dict[str, Any]:
    validate_params(p)
    sc, A, B, st = p["sigma_c"], p["A"], p["B"], p["sigma_t"]
    q2, g, l, h = p["q2"], p["gamma"], p["l"], p["h"]
    z = np.linspace(0.0, h, 200)

    xa, Fa = slip_profile(p, "active")
    xp, Fp = slip_profile(p, "passive")

    sig = sigma_of_z(z, q2, g)
    tau = tau_of_z(z, sc, st, q2, g, A, B)
    psi = np.degrees(psi_of_z(z, sc, st, q2, g, A, B))

    AZ, BZ = 0.6, 0.5
    Om = np.logspace(np.log10(0.005), np.log10(2.0), 300)
    psi0 = np.degrees(psi0_of_Omega(Om, AZ, BZ))

    zz = np.linspace(0.0, 2.0 * h, 300)
    sg = sigma_of_z(zz, q2, g)
    Ka_pow = Ka_local(sg, tau_of_z(zz, sc, st, q2, g, A, B),
                      psi_of_z(zz, sc, st, q2, g, A, B))
    phi, sinf = math.atan(A), math.sin(math.atan(A))
    Ka_linc = (1 - sinf) * (sg - st * sinf) / ((1 + sinf) * (sg + st * sinf))
    Ka_rank = (1 - sinf) / (1 + sinf) * np.ones_like(zz)

    p0 = np.linspace(0.5, 40.0, 200)
    series = []
    for om in [0.02, 0.05, 0.15, 0.3, 0.9]:
        ratio = np.tan(np.radians(p0)) / (B * om ** (B - 1.0))
        series.append({"Omega": om, "ratio": ratio.tolist()})

    def li(a):
        return [float(x) for x in np.asarray(a).tolist()]

    return {
        "mode": mode,
        "fig1_slip_surface": {"active_x": li(xa), "active_F": li(Fa),
                              "passive_x": li(xp), "passive_F": li(Fp),
                              "l": l, "h": h},
        "fig2_stress_field": {"z": li(z), "sigma": li(sig), "tau": li(tau)},
        "fig3_dilatancy": {"z": li(z), "psi": li(psi),
                           "psi_linear": li(np.degrees(np.arctan(A)) * np.ones_like(z))},
        "fig4_zoning": {"Omega": li(Om), "psi0": li(psi0),
                        "zone_omega": ZONE_OMEGA},
        "fig5_Ka_depth": {"z": li(zz), "Ka_power": li(Ka_pow),
                          "Ka_linc": li(Ka_linc), "Ka_rank": li(Ka_rank), "h": h},
        "fig6_nonassociated": {"psi0_deg": li(p0), "series": series},
    }


# ---------------------------------------------------------------------------
# 6. 试验对照复现（Table 3 主动 / Table 4 被动）
# ---------------------------------------------------------------------------
def reproduce_tables() -> Dict[str, Any]:
    def hcr_lin(l, phi):
        return l / math.tan(phi)
    def Ga(eta, phi):
        return 1.0 - eta * math.tan(phi)
    def Gb(eta, phi):
        return 1.0 + eta * math.tan(phi)

    phi_t = math.radians(31.6)
    eta_t, l_t, h_t = 1.0, 0.125, 0.250
    meas_Ga, psi_meas = 0.47, math.radians(24.25)
    a_assoc, a_non = Ga(eta_t, phi_t), Ga(eta_t, psi_meas)
    hcr_a, hcr_n = hcr_lin(l_t, phi_t) * 1000, hcr_lin(l_t, psi_meas) * 1000
    t3 = {
        "psi_deg": 24.25, "h_mm": 250.0,
        "hcr_assoc_mm": hcr_a, "hcr_non_mm": hcr_n,
        "meas_Ga": meas_Ga,
        "assoc_Ga": a_assoc, "assoc_dev_pct": (a_assoc / meas_Ga - 1) * 100,
        "non_Ga": a_non, "non_dev_pct": (a_non / meas_Ga - 1) * 100,
    }
    rows = [(0.5, 57.36, 26.5, 1.48), (1.0, 54.99, 30.0, 1.88),
            (1.5, 53.62, 35.5, 2.35), (2.0, 52.64, 40.5, 2.88),
            (3.0, 51.25, 44.5, 3.86)]
    t4rows, da, dn = [], [], []
    for eta_r, phi_pk, th, meas in rows:
        av, nv = Gb(eta_r, math.radians(phi_pk)), Gb(eta_r, math.radians(th))
        d1, d2 = (av / meas - 1) * 100, (nv / meas - 1) * 100
        da.append(d1); dn.append(d2)
        t4rows.append(dict(n=eta_r, phi_peak_deg=phi_pk, theta_deg=th,
                           meas=meas, assoc=av, assoc_dev_pct=d1,
                           non=nv, non_dev_pct=d2))
    return {
        "table3": t3,
        "table4": {"rows": t4rows,
                   "assoc_mean_dev_pct": float(np.mean(da)),
                   "non_mean_dev_pct": float(np.mean(dn))},
    }
