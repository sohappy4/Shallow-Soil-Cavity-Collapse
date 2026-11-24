# -*- coding: utf-8 -*-
"""
fem_anim.py —— 二维 FEM 应力场动态云图 / 演示动图
=================================================
把 solve_fem 的逐位移步场（field_history）做成：

  * frames：逐帧 PNG（base64），供前端滑块播放 —— 方案 M5「应力场动态云图」；
  * gif：GIF 动图（base64），自动循环播放 —— 方案 M6「演示视频」的最小实现。

每帧两面板：偏应力 q 云图（统一色标，避免帧间跳变）+ 塑性区（红=屈服）。
标题显示当前位移步、活门下移量、活门压力，直接对应压力—位移曲线的演化。

【边界】底层是 FEM 数值场（非运动学演示），但仍是横向扩展的定性—半定量结果：
DP 圆化 ≠ 幂型 Hoek–Brown，小变形理想弹塑性，不是严格 FELA 上下界。
"""
from __future__ import annotations

from typing import Dict, Any
import io
import base64

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

__all__ = ["build_fem_animation", "DEFAULT_ANIM_OPTS"]

DEFAULT_ANIM_OPTS: Dict[str, Any] = dict(
    dpi=105,
    fps=4,
    cmap="viridis",
)


def build_fem_animation(result: Dict[str, Any],
                        opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """从 solve_fem 结果生成逐帧 PNG 与 GIF 动图。

    result 必须包含 field_history（solve_fem 默认返回）。
    返回 dict(meta, frames, gif)，其中 frames/gif 均为 base64 字符串。
    """
    o = dict(DEFAULT_ANIM_OPTS)
    if opts:
        o.update({k: v for k, v in opts.items() if v is not None})

    hist = result.get("field_history") or []
    if not hist:
        raise ValueError("结果缺少 field_history，无法生成动态云图")
    curve = result["curve"]
    centers = np.asarray(result["field"]["centers"], dtype=float)
    x, y = centers[:, 0], centers[:, 1]

    # 全场统一色标，保证帧间可比
    qmax = max(float(np.max(f["q"])) for f in hist)
    qmax = max(qmax, 1.0)
    emax = max(float(np.max(f["epsp"])) for f in hist)
    emax = max(emax, 1e-9)

    frames = []
    for f in hist:
        step_no = int(f["step"])
        q = np.asarray(f["q"], dtype=float)
        ep = np.asarray(f["epsp"], dtype=float)
        yld = np.asarray(f["yielded"], dtype=bool)
        c = curve[step_no] if step_no < len(curve) else curve[-1]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 4.2))
        # 全场统一 levels，保证帧间色标可比（不随单帧最大值漂移）
        levels = np.linspace(0.0, qmax, 13)
        try:
            cn = ax1.tricontourf(x, y, q, levels=levels, cmap=o["cmap"])
            fig.colorbar(cn, ax=ax1, label="q_dev (kPa)")
        except Exception:
            sc = ax1.scatter(x, y, c=q, cmap=o["cmap"],
                             vmin=0.0, vmax=qmax, s=24)
            fig.colorbar(sc, ax=ax1, label="q_dev (kPa)")
        ax1.set_title("step %d   door %.3f m   p_door=%.1f kPa"
                      % (step_no, float(f["door_down"]), float(c["q_kPa"])))
        ax1.set_xlabel("x (m)")
        ax1.set_ylabel("y (m)")

        ax2.scatter(x[~yld], y[~yld], c="#9fb0c0", s=16, label="elastic")
        ax2.scatter(x[yld], y[yld], c="#d94b3a", s=20, label="yielded")
        ax2.set_title("plastic zone %.1f%%" % (100.0 * float(f["plastic_frac"])))
        ax2.set_xlabel("x (m)")
        ax2.set_ylabel("y (m)")
        ax2.legend(fontsize=8, loc="upper right")

        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=int(o["dpi"]))
        plt.close(fig)
        frames.append(base64.b64encode(buf.getvalue()).decode("ascii"))

    # GIF 演示动图（Pillow 随 matplotlib 安装）
    from PIL import Image
    imgs = [Image.open(io.BytesIO(base64.b64decode(s))).convert("RGB")
            for s in frames]
    gbuf = io.BytesIO()
    imgs[0].save(gbuf, format="GIF", save_all=True,
                 append_images=imgs[1:],
                 duration=int(1000.0 / max(1, int(o["fps"]))),
                 loop=0, optimize=False)

    return dict(
        meta=dict(
            n_frames=len(frames),
            fps=int(o["fps"]),
            qmax=float(qmax),
            emax=float(emax),
            boundary=result["meta"]["boundary"],
            note="FEM 数值场动态云图；DP 圆化、小变形理想弹塑性、非严格 FELA 上下界",
        ),
        frames=frames,
        gif=base64.b64encode(gbuf.getvalue()).decode("ascii"),
    )
