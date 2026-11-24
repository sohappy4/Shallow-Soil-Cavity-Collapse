# -*- coding: utf-8 -*-
"""plots.py —— 把 collapse_kernel.figure_data() 的数据画成 PNG(base64)。
只依赖 matplotlib（venv 已有），Agg 后端无窗口也能出图。
图中文字用英文，规避麒麟系统缺字体导致方块；中文说明放在网页正文。
"""
from __future__ import annotations

import io
import base64
from typing import Dict, Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ZONE_LABELS = ["I", "II", "III", "IV"]


def _b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def render_all(f: Dict[str, Any], omega: float, psi0_point: float) -> Dict[str, str]:
    """返回 {图名: base64 png}。omega/psi0_point 为当前工况在相图上的落位。"""
    out: Dict[str, str] = {}
    l, h = f["fig1_slip_surface"]["l"], f["fig1_slip_surface"]["h"]

    # fig1 滑动面形态
    d = f["fig1_slip_surface"]
    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    ax.plot(d["active_x"], d["active_F"], "b-", lw=2, label="Active (urn, inward)")
    ax.plot(-__import__("numpy").array(d["active_x"]), d["active_F"], "b-", lw=2)
    ax.plot(d["passive_x"], d["passive_F"], "g-", lw=2, label="Passive (funnel, outward)")
    ax.plot(-__import__("numpy").array(d["passive_x"]), d["passive_F"], "g-", lw=2)
    ax.axhline(h, color="k", lw=1); ax.axhline(0, color="k", lw=1)
    ax.plot([-l, l], [0, 0], "r-", lw=4, label="cavity")
    ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")
    ax.set_title("Fig.1  Slip-surface morphology")
    ax.legend(fontsize=8, loc="upper center"); ax.grid(alpha=0.3)
    fig.tight_layout(); out["fig1"] = _b64(fig)

    # fig2 应力场
    d = f["fig2_stress_field"]
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    ax.plot(d["sigma"], d["z"], "b-", lw=2, label="sigma = q2 + gamma*z (Eq.28)")
    ax.plot(d["tau"], d["z"], "r-", lw=2, label="tau (Eq.29)")
    ax.set_xlabel("stress (kPa)"); ax.set_ylabel("depth below surface (m)")
    ax.invert_yaxis(); ax.set_title("Fig.2  Stress field on slip surface")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); out["fig2"] = _b64(fig)

    # fig3 剪胀角
    d = f["fig3_dilatancy"]
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    ax.plot(d["psi"], d["z"], "b-", lw=2, label="power-type (B<1)")
    ax.plot(d["psi_linear"], d["z"], "k--", lw=1.5, label="linear (B=1), psi const")
    ax.set_xlabel("dilatancy angle psi (deg)"); ax.set_ylabel("depth below surface (m)")
    ax.invert_yaxis(); ax.set_title("Fig.3  Dilatancy angle vs depth (Eq.30)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); out["fig3"] = _b64(fig)

    # fig4 相图 + 当前落位
    d = f["fig4_zoning"]
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.semilogx(d["Omega"], d["psi0"], "b-", lw=2)
    zb = d["zone_omega"]
    for v, lab in zip(zb, ["I|II", "II|III", "III|IV"]):
        ax.axvline(v, color="grey", ls="--", lw=1)
    ax.axvspan(0.005, zb[0], color="0.95")
    ax.axvspan(zb[0], zb[1], color="#ffe08a")
    ax.axvspan(zb[1], zb[2], color="0.95")
    ax.axvspan(zb[2], 2.0, color="0.90")
    ax.text(0.008, 62, "I", fontsize=11, ha="center", weight="bold")
    ax.text(0.075, 62, "II", fontsize=11, ha="center", weight="bold")
    ax.text(0.45, 62, "III", fontsize=11, ha="center", weight="bold")
    ax.text(1.35, 62, "IV", fontsize=11, ha="center", weight="bold")
    ax.plot([omega], [psi0_point], "ro", ms=9, zorder=5, label="current case")
    ax.annotate("current", (omega, psi0_point), textcoords="offset points",
                xytext=(8, 6), fontsize=8, color="r")
    ax.set_xlabel("Omega = (sigma_t+q2)/sigma_c"); ax.set_ylabel("psi0 (deg)")
    ax.set_title("Fig.4  Morphology zoning (Eq.38)"); ax.set_ylim(0, 85)
    ax.set_xlim(0.005, 2.0); ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); out["fig4"] = _b64(fig)

    # fig5 Ka(z)
    d = f["fig5_Ka_depth"]
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    ax.fill_betweenx(d["z"], 0.27, 0.33, color="0.85", label="back-calc band")
    ax.plot(d["Ka_power"], d["z"], "b-", lw=2, label="power-type (Eq.44)")
    ax.plot(d["Ka_linc"], d["z"], "-.", c="g", lw=1.8, label="linear+cohesion")
    ax.plot(d["Ka_rank"], d["z"], "--", c="k", lw=1.5, label="Rankine")
    ax.axhline(d["h"], color="r", lw=1, ls=":", label="ground depth h")
    ax.set_xlabel("Ka(z)"); ax.set_ylabel("depth below surface (m)")
    ax.invert_yaxis(); ax.set_xlim(0, 1.05)
    ax.set_title("Fig.5  Lateral pressure coefficient (Eq.44)")
    ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3)
    fig.tight_layout(); out["fig5"] = _b64(fig)

    # fig6 非关联修正
    d = f["fig6_nonassociated"]
    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    for s in d["series"]:
        ax.plot(d["psi0_deg"], s["ratio"], lw=1.7, label="Omega=%s" % s["Omega"])
    ax.axhline(1.0, color="k", ls="--", lw=1.5, label="associated (A_psi=A)")
    ax.axvspan(2, 15, color="#ffd0a0", alpha=0.6)
    ax.set_xlabel("measured psi0 (deg)"); ax.set_ylabel("A_psi/A (Eq.46)")
    ax.set_title("Fig.6  Non-associated correction"); ax.set_ylim(0, 2)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); out["fig6"] = _b64(fig)

    return out
