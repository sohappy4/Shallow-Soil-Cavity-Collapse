# -*- coding: utf-8 -*-
"""
collapse_anim.py —— 二维坍塌运动学动画（第一层可视化）
======================================================
目标：让"土洞坍塌"真正动起来——土体粒子沿滑动面下落、活门下移、
地表漏斗成形。它把已有三层内核的结果驱动成逐帧几何：

  * 滑体边界 x(z)      ← L1 collapse_kernel.slip_profile（或 L2 打靶）
  * 前锋高度 z_front(t) ← L3 l3_rate.solve_rate 的 timeline（v·t）
  * 漏斗终径 R          ← L3 体积守恒 R=(0.75·B_c²·h/tanφ_r)^(1/3)

【诚实边界】这是运动学演示(kinematic visualization)，不是应力场求解：
  - 粒子位移场是构造出来的运动学模板，满足"滑面内才动、前锋扫过才动、
    体积守恒缩放"，不是 DEM 颗粒接触求解，也不是有限元应力云图；
  - 不输出任何应力/塑性区数值；这些属 FELA/FEM（方案 M7）远期功能；
  - 被动工况按主动镜像显示（顶起隆起），与本书"镜像对称"结论一致。

坐标约定（与 fig1 一致）：x 水平（洞中线 x=l），z 竖直向下（地表 z=0，
洞底 z=h）。返回逐帧归一化几何，前端自行画 SVG。
"""
from __future__ import annotations

from typing import Dict, Any, List
import math
import numpy as np

from collapse_kernel import slip_profile, validate_params
from l3_rate import solve_rate

__all__ = ["build_animation", "DEFAULT_ANIM_PARAMS"]


DEFAULT_ANIM_PARAMS: Dict[str, Any] = dict(
    n_frames=48,        # 帧数
    n_particles_x=26,   # 粒子横向网格数
    n_particles_z=18,   # 粒子纵向网格数
    passive=False,      # True=被动顶起镜像
)


def _shape_weight(x: float, z: float, xl: float, xr: float, h: float) -> float:
    """粒子在滑体内的相对活动权重：滑面附近(边界)为 0，滑体中轴为 1。"""
    if xr - xl <= 1e-12:
        return 0.0
    rel = (x - xl) / (xr - xl)          # 0..1
    w = math.sin(math.pi * max(0.0, min(1.0, rel)))
    return w


def build_animation(params: Dict[str, float],
                    rate: Dict[str, Any] = None,
                    opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """生成逐帧坍塌运动学几何。

    params: L1 参数（含 mode 所需 h、l、B、A 等）
    rate:   可选，solve_rate 的结果；缺省时内部按 params 计算
    opts:   n_frames / n_particles_x / n_particles_z / passive
    返回:
      {
        "meta": {... 边界说明、物理参数、适用性 ...},
        "geometry": {"l":..., "h":..., "x_outcrop":..., "R":...},
        "frames": [
           {"i":0,"t":0,"z_front":0,"door_down":0,"crater_r":0,
            "surface":[[x,z],...], "slip":[[x,z],...],
            "particles":[[x,z,mode_flag],...], "stage":"...", "note":"..."},
           ...
        ]
      }
      particles 第三列 mode_flag: 0=静止 1=运动中（供前端着色）
    """
    validate_params(params)
    o = dict(DEFAULT_ANIM_PARAMS)
    if opts:
        o.update({k: v for k, v in opts.items() if v is not None})
    n_frames = max(8, min(90, int(o["n_frames"])))
    nx = max(8, min(60, int(o["n_particles_x"])))
    nz = max(6, min(40, int(o["n_particles_z"])))
    passive = bool(o["passive"])

    p = {k: float(v) for k, v in params.items()}
    h, l = p["h"], p["l"]
    D = 2.0 * l

    if rate is None:
        rate = solve_rate(p, phi_deg=math.degrees(math.atan(p["A"])),
                          gamma=p["gamma"])
    T = float(rate["T"])
    R_final = float(rate["R"])
    Bc = float(rate["params"]["B_c"])
    timeline = rate.get("timeline", [])

    # 滑面边界（主动支；被动按镜像取用同一几何）
    x_slip, F = slip_profile(p, "active", npts=160)
    # x_slip: 从 x2(地表出露) 到 l(洞底)，F=h-z：F=h→z=0 地表, F=0→z=h 洞底
    z_slip = h - F
    # 转成按 z 升序（地表 0 → 洞底 h）
    order = np.argsort(z_slip)
    zs = z_slip[order]
    xs = x_slip[order]
    x_outcrop = float(xs[0])      # 地表出露半宽
    x_floor = float(xs[-1])       # 洞底半宽（≈l）

    def slip_half_width(z: float) -> float:
        return float(np.interp(z, zs, xs)) if len(zs) > 1 else l

    # 初始土体粒子网格：滑体横向区间 [x_l, x_r]，
    # 左半支 x_l=slip(z)<=l，右半支 x_r=2l-x_l（关于洞中线镜像）。
    particles0: List[tuple] = []
    for iz in range(nz):
        z = h * (iz + 0.5) / nz
        x_l = slip_half_width(z)
        x_r = 2.0 * l - x_l
        if x_r - x_l <= 1e-9:
            continue
        for ix in range(nx):
            x = x_l + (x_r - x_l) * (ix + 0.5) / nx
            particles0.append((float(x), float(z)))

    def front_at(t: float) -> float:
        if timeline:
            for k in range(1, len(timeline)):
                a, b = timeline[k - 1], timeline[k]
                if t <= b["t"]:
                    f = (t - a["t"]) / ((b["t"] - a["t"]) or 1.0)
                    return float(a["z"] + f * (b["z"] - a["z"]))
            return float(timeline[-1]["z"])
        return min(float(rate["v"]) * t, h)

    def stage_at(t: float) -> str:
        name = "静力核查"
        for s in rate.get("stages", []):
            if t >= s["t"] - 1e-12:
                name = s["name"]
        return name

    frames: List[Dict[str, Any]] = []
    for i in range(n_frames):
        frac = i / (n_frames - 1)
        t = frac * T
        zf = front_at(t)

        # 活门下移：与前锋推进同步，限幅 0.35D（视觉可读，非真实位移）
        door_down = min(0.35 * D, (zf / h) * 0.35 * D)
        # 地表漏斗半径：前锋接近地表时成形，∝ frac² 到终径 R
        # 用 frac(=t/T) 而非 zf/h，避免 timeline 插值浮点误差被平方放大，
        # 保证末帧严格等于终径。
        crater_r = R_final * frac ** 2 if h > 0 else 0.0

        pts: List[List[float]] = []
        # zf 为前锋自洞底向上推进的高度（0=洞底, h=地表）。
        # 被扫过区：z > h - zf（靠近洞底一侧先动）；未被扫过：z <= h - zf 静止。
        for (x0, z0) in particles0:
            swept = z0 - (h - zf)       # >0 表示该粒子已被前锋扫过
            # 主动瓮形：滑体横向区间 [x_l, x_r]，左半支 x_l=slip(z)<=l，
            # 右半支 x_r=2l-x_l（关于洞中线镜像），地表宽、洞底收窄至洞宽。
            x_left = slip_half_width(z0)
            x_right = 2.0 * l - x_left
            weight = _shape_weight(x0, z0, min(x_left, x_right),
                                   max(x_left, x_right), h)
            if swept > 0:
                # 被扫过越久（越靠洞底）下落越多，限幅避免穿出洞底过多
                drop = min(0.28 * D, 0.5 * swept) * weight
                # 向洞中线汇聚：越靠洞底收拢越明显
                gather = 0.35 * (z0 / h) * weight
                xn = x0 + (l - x0) * gather
                zn = z0 + drop
                flag = 1
            else:
                # 前锋之上未扫过：静止；处于地表漏斗范围时随漏斗下沉
                if crater_r > 0 and z0 < 0.14 * h and abs(x0 - l) < crater_r:
                    dx = abs(x0 - l)
                    zn = z0 + 0.09 * h * (1.0 - (dx / crater_r) ** 2) * (zf / h)
                    xn = x0
                    flag = 1
                else:
                    xn, zn, flag = x0, z0, 0
            pts.append([round(float(xn), 4), round(float(zn), 4), flag])

        # 地表线（漏斗凹陷）
        surface: List[List[float]] = []
        for k in range(41):
            x = D * 1.6 * k / 40.0
            if crater_r > 0 and abs(x - l) < crater_r:
                z = 0.10 * h * (1.0 - ((x - l) / crater_r) ** 2) * (zf / h)
            else:
                z = 0.0
            surface.append([round(float(x), 4), round(float(z), 4)])

        # 滑面折线（左右两支）
        slip_pts: List[List[float]] = []
        for zz, xx in zip(zs, xs):
            slip_pts.append([round(float(xx), 4), round(float(zz), 4)])
        for zz, xx in zip(zs[::-1], xs[::-1]):
            slip_pts.append([round(float(2.0 * l - xx), 4), round(float(zz), 4)])

        frames.append(dict(
            i=i, t=round(float(t), 6), z_front=round(float(zf), 6),
            door_down=round(float(door_down), 6), crater_r=round(float(crater_r), 8),
            stage=stage_at(t), surface=surface, slip=slip_pts, particles=pts,
        ))

    boundary_note = (
        "运动学演示：粒子位移为满足滑面约束与体积守恒的构造模板，"
        "不是 DEM 颗粒接触求解，也不是有限元应力云图。"
        "本动画不输出应力/塑性区数值；那属 FELA/FEM 远期功能（方案 M7）。"
    )
    return dict(
        meta=dict(
            solution_type="kinematic-animation",
            passive=passive, T=float(T), v=float(rate["v"]), R_final=R_final,
            boundary_note=boundary_note,
            domain=rate.get("domain", ""), domain_dev_pct=rate.get("domain_dev_pct", 0.0),
            warnings=list(rate.get("warnings", [])),
            n_frames=n_frames, n_particles=len(particles0),
        ),
        geometry=dict(l=float(l), h=float(h), D=D, x_outcrop=x_outcrop,
                      x_floor=x_floor, R_final=R_final),
        frames=frames,
    )
