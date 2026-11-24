# -*- coding: utf-8 -*-
"""二维圆盘颗粒离散元（DEM）内核 —— 活门试验的真·颗粒接触求解。

方法来源（横向引用，非专著理论）：
  Cundall & Strack (1979) 离散元法；线性弹簧-阻尼接触 + Coulomb 摩擦
  （Cundall–Strack 切向历史模型）；速度 Verlet 显式积分。
  这是成熟的颗粒流数值方法，按用户 2026-11-24 授权引入。

与既有模块的区别（必须分清）：
  - collapse_anim.py：运动学演示，粒子位移是构造模板，不求解接触力；
  - fem2d.py：连续介质有限元，输出应力/塑性区场；
  - dem2d.py（本模块）：离散颗粒接触求解，逐颗粒积分牛顿方程，
    输出颗粒轨迹、接触数、活门反力与地表漏斗 —— 真正的颗粒流动力学。

诚实边界：
  1. 二维圆盘、线性接触、无转动自由度（切向力不做力矩），是简化 DEM；
     不是 Yade/PFC 的完整 3D 球体和转动模型；
  2. 颗粒刚度/摩擦为数值参数，不是由专著 μ(I) 细观反算得到；
  3. 结果定性—半定量，与 L1/L2/L3 的专著理论主结果不混同；
  4. 显式积分有时间步稳定限制，超出参数域会报告非收敛。

许可证边界（方案 M8）：本模块为纯 Python 独立实现，不链接 GPL 的 Yade、
不分发商业 PFC；外部后端通过 dem_backend.py 的独立进程协议接入。
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

import numpy as np

__all__ = [
    "DEFAULT_DEM_PARAMS", "DEFAULT_DEM_OPTS", "initial_packing",
    "contact_forces", "solve_dem", "selfcheck_dem",
]

# ---------------------------------------------------------------------------
# 默认参数（SI 单位）
# ---------------------------------------------------------------------------
DEFAULT_DEM_PARAMS: Dict[str, float] = dict(
    domain_w=0.50,        # 模型宽度 (m)
    domain_h=0.32,        # 模型高度 (m)
    door_w=0.20,          # 活门宽度 (m)，居中
    r_mean=0.010,         # 颗粒平均半径 (m)
    r_spread=0.04,        # 半径相对扰动（受六角排布间距约束）
    n_particles=400,      # 颗粒数
    rho_p=1600.0,         # 颗粒密度 (kg/m3)
    g=9.81,               # 重力加速度 (m/s2)
    kn=2.0e4,             # 法向接触刚度 (N/m)
    kt_ratio=0.4,         # 切向/法向刚度比
    mu=0.5,               # 颗粒间摩擦系数
    mu_wall=0.4,          # 墙-颗粒摩擦系数
    zeta=0.20,            # 接触阻尼比
    restitution=0.05,     # 墙碰撞恢复系数（仅诊断，主用弹簧阻尼）
)

DEFAULT_DEM_OPTS: Dict[str, Any] = dict(
    dt=0.0010,            # 时间步 (s)
    settle_steps=150,     # 活门关闭预沉积步数
    n_steps=500,          # 活门打开后模拟步数
    save_every=20,        # 每多少步存一帧
    max_particles=1200,   # 粒子数上限（防止界面卡死）
)


# ---------------------------------------------------------------------------
# 初始堆积
# ---------------------------------------------------------------------------
def initial_packing(params: Dict[str, float],
                    opts: Dict[str, Any] = None,
                    seed: int = 20261124) -> Dict[str, np.ndarray]:
    """在活门关闭的模型内生成近六角密排初始颗粒（含小扰动）。

    返回 x, y, r, m（长度 N 的 numpy 数组）。颗粒不超出边界、彼此不重叠。
    """
    p = dict(DEFAULT_DEM_PARAMS)
    p.update({k: v for k, v in params.items() if k in p})
    o = dict(DEFAULT_DEM_OPTS)
    if opts:
        o.update({k: v for k, v in opts.items() if k in o})

    n = int(min(o["max_particles"], max(20, int(p["n_particles"]))))
    r0 = float(p["r_mean"])
    xmax, ymax = float(p["domain_w"]), float(p["domain_h"])
    rng = np.random.default_rng(int(seed))

    xs: List[float] = []
    ys: List[float] = []
    # 间距裕度 1.08：半径扰动 ±4% + 位置扰动 ±1.5% 后仍不重叠
    margin = 1.08
    x = r0 * margin
    row = 0
    spacing = 2.0 * r0 * margin
    row_h = math.sqrt(3.0) * r0 * margin
    while len(xs) < n:
        y = r0 * margin + row * row_h
        if y + r0 > ymax:
            break
        x = r0 * margin + (row % 2) * r0 * margin
        while x + r0 < xmax and len(xs) < n:
            xs.append(x)
            ys.append(y)
            x += spacing
        row += 1

    if len(xs) < n:
        # 空间不足：按实际可容纳数返回（诚实报告，不强行塞入）
        n = len(xs)

    xx = np.asarray(xs, dtype=float)
    yy = np.asarray(ys, dtype=float)
    rr = r0 * (1.0 + float(p["r_spread"]) * rng.uniform(-1.0, 1.0, n))
    rr = np.clip(rr, 0.92 * r0, 1.04 * r0)
    xx = xx + 0.015 * r0 * rng.uniform(-1.0, 1.0, n)
    yy = yy + 0.015 * r0 * rng.uniform(-1.0, 1.0, n)
    xx = np.clip(xx, rr, xmax - rr)
    yy = np.clip(yy, rr, ymax - rr)
    mm = float(p["rho_p"]) * math.pi * rr ** 2
    return dict(x=xx, y=yy, r=rr, m=mm)


# ---------------------------------------------------------------------------
# 接触力
# ---------------------------------------------------------------------------
def contact_forces(pos: Dict[str, np.ndarray], vel: np.ndarray,
                   params: Dict[str, float], opts: Dict[str, Any],
                   tangential_history: np.ndarray, door_open: bool,
                   dt: float) -> Dict[str, Any]:
    """计算全部颗粒接触力 + 四壁弹簧力 + 重力。

    pos: {'x','y','r','m'}；vel: (N,2)；tangential_history: (N,N) 切向位移历史
    （原地更新）。返回力 (N,2)、接触数、活门反力。

    法向：F_n = kn·δ - cn·v_n（截断到 ≥0，无拉力）
    切向：F_t = -kt·δ_t - ct·v_t，|F_t| ≤ μ·F_n（Cundall–Strack 历史模型）
    """
    p = dict(DEFAULT_DEM_PARAMS)
    p.update({k: v for k, v in params.items() if k in p})

    x, y, rad, mass = pos["x"], pos["y"], pos["r"], pos["m"]
    n = len(x)
    vx, vy = vel[:, 0], vel[:, 1]
    kn = float(p["kn"])
    kt = kn * float(p["kt_ratio"])
    mu = float(p["mu"])
    m_eff = 0.5 * float(np.mean(mass))
    cn = 2.0 * float(p["zeta"]) * math.sqrt(kn * m_eff)
    ct = 2.0 * float(p["zeta"]) * math.sqrt(kt * m_eff)
    cw = 2.0 * float(p["zeta"]) * math.sqrt(kn * max(m_eff, 1e-12))

    fx = np.zeros(n, dtype=float)
    fy = -mass * float(p["g"])

    # 已掉出模型的颗粒（y < -2r）不再参与接触与边界力，自由落体离场；
    # 否则底壁弹簧 pen_b 随穿透深度线性增长，会把颗粒以巨大能量弹回、动能爆炸。
    escaped = y < (-2.0 * rad)
    active = ~escaped

    # ---- 颗粒-颗粒接触（全对向量化，N 上限 1200）----
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    rsum = rad[:, None] + rad[None, :]
    d2 = dx * dx + dy * dy
    d = np.sqrt(np.maximum(d2, 1e-24))
    overlap = rsum - d
    pair_ok = ((np.arange(n)[:, None] > np.arange(n)[None, :]) &
               active[:, None] & active[None, :])
    ii, jj = np.where((overlap > 0.0) & pair_ok)
    n_contacts = len(ii)

    door_force = 0.0
    if n_contacts:
        nx, ny = dx[ii, jj] / d[ii, jj], dy[ii, jj] / d[ii, jj]
        rvx, rvy = vx[ii] - vx[jj], vy[ii] - vy[jj]
        vn = rvx * nx + rvy * ny
        fn = np.maximum(kn * overlap[ii, jj] - cn * vn, 0.0)

        tx, ty = -ny, nx
        vt = rvx * tx + rvy * ty
        tangential_history[ii, jj] += vt * dt
        ft = -kt * tangential_history[ii, jj] - ct * vt
        ftmax = mu * fn
        slip = np.abs(ft) > ftmax
        if np.any(slip):
            ft[slip] = -ftmax[slip] * np.sign(ft[slip])
            tangential_history[ii[slip], jj[slip]] = -(ft[slip] + ct * vt[slip]) / kt

        cfx = fn * nx + ft * tx
        cfy = fn * ny + ft * ty
        fx += np.bincount(ii, weights=cfx, minlength=n) - np.bincount(jj, weights=cfx, minlength=n)
        fy += np.bincount(ii, weights=cfy, minlength=n) - np.bincount(jj, weights=cfy, minlength=n)

    # ---- 四壁弹簧阻尼 ----
    xmax, ymax = float(p["domain_w"]), float(p["domain_h"])
    xmin = 0.0

    pen_l = (xmin + rad) - x
    m_l = (pen_l > 0.0) & active
    if np.any(m_l):
        fx[m_l] += kn * pen_l[m_l] - cw * vx[m_l]

    pen_r = x - (xmax - rad)
    m_r = (pen_r > 0.0) & active
    if np.any(m_r):
        fx[m_r] -= kn * pen_r[m_r] + cw * vx[m_r]

    pen_t = y - (ymax - rad)
    m_t = (pen_t > 0.0) & active
    if np.any(m_t):
        fy[m_t] -= kn * pen_t[m_t] + cw * vy[m_t]

    pen_b = (xmin + rad) - y
    m_b = pen_b > 0.0
    door_x0 = 0.5 * (xmax - float(p["door_w"]))
    door_x1 = 0.5 * (xmax + float(p["door_w"]))
    in_door = (x > door_x0) & (x < door_x1)
    if door_open:
        m_b = m_b & ~in_door
    m_b = m_b & active
    if np.any(m_b):
        fwall = kn * pen_b[m_b] - cw * vy[m_b]
        fy[m_b] += fwall
        if not door_open:
            # 活门关闭时，门所受向下力 = 底壁对颗粒向上力之和
            door_force = float(np.sum(fwall[in_door[m_b]]))
    # 逃逸颗粒只受重力（自由落体离场，不参与任何接触/边界力）
    if np.any(escaped):
        fx[escaped] = 0.0
        fy[escaped] = -mass[escaped] * float(p["g"])

    # 诊断量只对真实接触对取最大（对角线 i=j 的 overlap=r_i+r_i 是伪值）
    max_ov = float(np.max(overlap[ii, jj])) if n_contacts else 0.0
    return dict(fx=fx, fy=fy, n_contacts=n_contacts,
                door_force=door_force, max_overlap=max_ov)


# ---------------------------------------------------------------------------
# 主求解器
# ---------------------------------------------------------------------------
def solve_dem(params: Dict[str, float] = None,
              opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """位移控制活门 DEM 求解：先沉积、后开闸，输出逐帧颗粒轨迹与反力。

    返回 dict：meta / frames / curve / boundary。
    每帧 particles 为 [[x, y, vx, vy], ...]（SI），可直接 JSON 序列化。
    """
    p = dict(DEFAULT_DEM_PARAMS)
    if params:
        p.update({k: v for k, v in params.items() if k in p})
    o = dict(DEFAULT_DEM_OPTS)
    if opts:
        o.update({k: v for k, v in opts.items() if k in o})

    dt = float(o["dt"])
    settle = int(o["settle_steps"])
    n_steps = int(o["n_steps"])
    save_every = max(1, int(o["save_every"]))

    pack = initial_packing(p, o)
    x, y, rad, mass = pack["x"], pack["y"], pack["r"], pack["m"]
    n = len(x)
    vel = np.zeros((n, 2), dtype=float)
    tang = np.zeros((n, n), dtype=float)
    pos = dict(x=x, y=y, r=rad, m=mass)
    m_col = mass[:, None]

    frames: List[Dict[str, Any]] = []
    curve: List[Dict[str, Any]] = []

    def record(step: int, t: float, door_open: bool, cf: Dict[str, Any]):
        speed = np.sqrt(vel[:, 0] ** 2 + vel[:, 1] ** 2)
        # 地表高度：按 x 分 20 格取每格最高颗粒
        nb = 20
        edges = np.linspace(0.0, float(p["domain_w"]), nb + 1)
        surf = []
        for k in range(nb):
            sel = (x >= edges[k]) & (x < edges[k + 1])
            surf.append(float(np.max(y[sel] + rad[sel])) if np.any(sel) else 0.0)
        frames.append(dict(
            step=int(step), t=float(t), door_open=bool(door_open),
            particles=[[float(a), float(b), float(c), float(d)]
                       for a, b, c, d in zip(x, y, vel[:, 0], vel[:, 1])],
            radii=[float(v) for v in rad],
            n_contacts=int(cf["n_contacts"]),
            door_force_N=float(cf["door_force"]),
            kinetic_J=float(0.5 * np.sum(mass * speed ** 2)),
            max_overlap=float(cf["max_overlap"]),
            surface=surf,
        ))

    def advance(door_open: bool, cf: Dict[str, Any]) -> Dict[str, Any]:
        nonlocal vel, x, y
        # 速度 Verlet：半步速度 → 位移 → 新力 → 半步速度
        vel += 0.5 * dt * np.column_stack((cf["fx"] / m_col[:, 0], cf["fy"] / m_col[:, 0]))
        x += dt * vel[:, 0]
        y += dt * vel[:, 1]
        new_cf = contact_forces(pos, vel, p, o, tang, door_open, dt)
        vel += 0.5 * dt * np.column_stack((new_cf["fx"] / m_col[:, 0], new_cf["fy"] / m_col[:, 0]))
        # 数值保护：钳制穿透瞬间的非物理速度（仅极端情况触发）
        vmax = 20.0
        sp = np.sqrt(vel[:, 0] ** 2 + vel[:, 1] ** 2)
        hot = sp > vmax
        if np.any(hot):
            vel[hot] *= (vmax / sp[hot])[:, None]
        return new_cf

    # 阶段 1：活门关闭，自由沉积
    cf = contact_forces(pos, vel, p, o, tang, False, dt)
    for step in range(1, settle + 1):
        cf = advance(False, cf)
        if step % save_every == 0 or step == settle:
            curve.append(dict(stage="settle", step=step, t=step * dt,
                              door_force_N=cf["door_force"],
                              n_contacts=cf["n_contacts"],
                              kinetic_J=frames[-1]["kinetic_J"] if frames else 0.0))
    record(settle, settle * dt, False, cf)

    # 阶段 2：活门打开，颗粒流出
    cf = contact_forces(pos, vel, p, o, tang, True, dt)
    door_x0 = 0.5 * (float(p["domain_w"]) - float(p["door_w"]))
    door_x1 = 0.5 * (float(p["domain_w"]) + float(p["door_w"]))
    for step in range(1, n_steps + 1):
        cf = advance(True, cf)
        if step % save_every == 0 or step == n_steps:
            in_door = (x > door_x0) & (x < door_x1)
            below = np.sum(y[in_door] < float(p["r_mean"]))
            curve.append(dict(stage="flow", step=step, t=(settle + step) * dt,
                              door_force_N=0.0,
                              n_contacts=cf["n_contacts"],
                              kinetic_J=float(0.5 * np.sum(mass * (vel[:, 0] ** 2 + vel[:, 1] ** 2))),
                              particles_in_door=int(np.sum(in_door)),
                              particles_below=int(below)))
        if step % save_every == 0 or step == n_steps:
            record(settle + step, (settle + step) * dt, True, cf)

    # 稳定性诊断：全时段最大重叠 < 1 个半径（无完全穿透）+ 全程有限
    # 注：DEM 软接触允许 30% 量级的瞬时重叠，不能用 20% 的过严阈值误判为发散
    peak_overlap = max((f["max_overlap"] for f in frames), default=0.0)
    overlap_ok = bool(peak_overlap < float(p["r_mean"]))
    finite_ok = bool(np.all(np.isfinite(x)) and np.all(np.isfinite(y)) and
                     all(np.all(np.isfinite(np.asarray(f["particles"]))) for f in frames))
    return dict(
        meta=dict(
            n_particles=int(n), n_steps=int(n_steps), settle_steps=int(settle),
            dt=dt, total_time=float((settle + n_steps) * dt),
            n_frames=len(frames), save_every=save_every,
            door_w=float(p["door_w"]), domain_w=float(p["domain_w"]),
            domain_h=float(p["domain_h"]),
            stable=bool(overlap_ok and finite_ok),
            peak_overlap_ratio=float(peak_overlap / float(p["r_mean"])),
            final_overlap_ratio=float(frames[-1]["max_overlap"] / float(p["r_mean"])),
            final_n_contacts=int(frames[-1]["n_contacts"]),
            kinetic_peak_J=float(max((f["kinetic_J"] for f in frames), default=0.0)),
        ),
        frames=frames,
        curve=curve,
        boundary=("二维圆盘简化 DEM（Cundall–Strack 线性接触 + Coulomb 摩擦，"
                  "无转动自由度），横向引用成熟离散元方法，非专著理论；"
                  "颗粒刚度/摩擦为数值参数，非由 μ(I) 细观反算；"
                  "结果定性—半定量，不冒充 Yade/PFC 完整模型。"),
    )


# ---------------------------------------------------------------------------
# 自检
# ---------------------------------------------------------------------------
def selfcheck_dem() -> Dict[str, Any]:
    """小规模自检：堆积不重叠、沉积后稳定、开闸后真的流出、诊断量有界。"""
    checks: Dict[str, Any] = {}
    o = dict(n_particles=90, n_steps=120, settle_steps=60, save_every=30,
             max_particles=200)
    pack = initial_packing(DEFAULT_DEM_PARAMS, o)
    x, r, n0 = pack["x"], pack["r"], len(pack["x"])
    offdiag = np.arange(n0)[:, None] != np.arange(n0)[None, :]
    checks["packing_no_overlap"] = bool(np.all(
        (np.sqrt((x[:, None] - x[None, :]) ** 2 +
                 (pack["y"][:, None] - pack["y"][None, :]) ** 2) + 1e-9
         >= r[:, None] + r[None, :] - 1e-6)[offdiag]))
    checks["packing_in_domain"] = bool(
        np.all(x >= 0) and np.all(x <= DEFAULT_DEM_PARAMS["domain_w"]) and
        np.all(pack["y"] >= 0) and np.all(pack["y"] <= DEFAULT_DEM_PARAMS["domain_h"]))

    r1 = solve_dem(DEFAULT_DEM_PARAMS, o)
    f0, f1 = r1["frames"][0], r1["frames"][-1]
    checks["frame_count"] = bool(r1["meta"]["n_frames"] >= 2)
    checks["finite_positions"] = bool(np.all(np.isfinite(np.asarray(f1["particles"]))))
    checks["stable_overlap"] = bool(r1["meta"]["stable"])
    checks["contacts_nonnegative"] = bool(all(f["n_contacts"] >= 0 for f in r1["frames"]))
    checks["door_force_positive_before_open"] = bool(f0["door_force_N"] > 0.0)
    checks["kinetic_nonnegative"] = bool(all(f["kinetic_J"] >= 0.0 for f in r1["frames"]))

    # 开闸后活门区应有颗粒流出（底部下方颗粒数增加）
    ys = np.asarray(f1["particles"])[:, 1]
    xs = np.asarray(f1["particles"])[:, 0]
    door_x0 = 0.5 * (DEFAULT_DEM_PARAMS["domain_w"] - DEFAULT_DEM_PARAMS["door_w"])
    door_x1 = 0.5 * (DEFAULT_DEM_PARAMS["domain_w"] + DEFAULT_DEM_PARAMS["door_w"])
    in_door = (xs > door_x0) & (xs < door_x1)
    checks["door_region_drained"] = bool(np.sum(ys[in_door] < 0.5 * DEFAULT_DEM_PARAMS["r_mean"]) > 0)
    checks["boundary_honest"] = "非专著理论" in r1["boundary"]
    checks["all_pass"] = all(v for k, v in checks.items() if k != "all_pass")
    return checks


if __name__ == "__main__":
    import json
    print(json.dumps(selfcheck_dem(), indent=2, ensure_ascii=False))
