# -*- coding: utf-8 -*-
"""
benchmarks.py —— 验证基准集（分期）
=====================================
只录入专著/论文中"有明确数值"的锚点，不虚构逐点试验数据。
一期：Table 3/4（见 collapse_kernel.reproduce_tables）+ 三项闭式基准
二期：1g 出露角、离心机 beta_eff、FELA N_sat —— 均为文献锚点，用于量级核对
"""
from __future__ import annotations

from typing import Dict, Any

__all__ = ["BENCHMARKS", "check_benchmarks", "saturation_anchor"]

BENCHMARKS: Dict[str, Any] = {
    "phase1": {
        "psi0_deg": 37.56,
        "psi_h_deg": 22.67,
        "h_cr_m": 5.389,
        "mirror_error": 0.0,
        "table3_assoc_dev_pct": -18.1,
        "table3_non_dev_pct": 17.0,
        "table4_assoc_mean_pct": 25.4,
        "table4_non_mean_pct": -9.5,
    },
    "phase2_literature": {
        "outcrop_angle_1g_deg": {
            "value": 61.0, "source": "Abdulla & Goodings (1996), 1g 沉陷坑出露角",
            "type": "定性锚点：本书形态谱出露角的定量佐证",
        },
        "beta_eff_centrifuge": {
            "value_min": 4.17, "value_max": 5.26, "scatter_pct": 11.0,
            "source": "Dewoolkar et al. (2007), 密砂 150g 系列, h/D=2~6",
            "type": "单参数 beta 全程不漂移，Janssen 型侧壁摩擦闭式校验",
        },
        "N_sat_FELA": {
            "value": 18.2212, "h_sat_over_D": 4.555, "tolerance_pct": 1.0,
            "source": "深埋圆柱孔静力下界 N=18.22；本书 N_sat=18.2212 四位重合",
            "type": "上下界夹逼至 1%（极限分析常数，无初等闭式）",
        },
        "h_sat_cr_distinction": {
            "value": "不同门槛",
            "source": "专著 §18.2.2 / §20.4.2",
            "type": "h_sat 是力学饱和高度（拱效应层）；h_cr 是几何贯通阈值（上限法）——不可混用",
        },
    },
    "phase3": {
        "FELA_bounds": "待实现：SOCP 上下界夹逼内核",
        "DEM_backend": "待实现：Yade/PFC 独立进程细观标定",
        "three_d_cone": "待实现：三维锥形破坏求解器（专著称进行中）",
    },
}


def check_benchmarks(psi0, psi_h, h_cr, mirror_error, t3, t4) -> Dict[str, Any]:
    """核对一期基准，返回逐项通过标志。"""
    b = BENCHMARKS["phase1"]
    checks = {
        "psi0": abs(psi0 - b["psi0_deg"]) < 0.01,
        "psi_h": abs(psi_h - b["psi_h_deg"]) < 0.01,
        "h_cr": abs(h_cr - b["h_cr_m"]) < 0.01,
        "mirror": mirror_error < 1e-10,
        "t3_assoc": abs(t3["assoc_dev_pct"] - b["table3_assoc_dev_pct"]) < 0.2,
        "t3_non": abs(t3["non_dev_pct"] - b["table3_non_dev_pct"]) < 0.2,
        "t4_assoc": abs(t4["assoc_mean_dev_pct"] - b["table4_assoc_mean_pct"]) < 0.2,
        "t4_non": abs(t4["non_mean_dev_pct"] - b["table4_non_mean_pct"]) < 0.2,
    }
    return {"checks": checks, "all_pass": all(checks.values())}


def saturation_anchor(params: Dict[str, float],
                      h_over_D_list=(1.0, 1.4, 2.0, 4.555, 6.0, 10.0)) -> Dict[str, Any]:
    """文献锚点核验：N_sat=18.2212, h_sat≈4.555 D。

    【定位】这是极限分析常数，专著 §18.2.2 / §20.4.2 明确其无初等闭式，
    由上下界夹逼至 1% 以内；本程序不声称从 L1/L2/L3 推导出它。
    本函数只做三件可核验的事：
      (1) 常数内部自洽：h_sat/D = N_sat/4；
      (2) 区分 h_sat（力学饱和高度）与 h_cr（几何贯通阈值），二者不可混用；
      (3) 给出不同 h/D 下 L1 全高动员闭式外推值与 L3 速率模型适用域判定，
          并明确超出贯通域后闭式外推不具物理意义。
    """
    from collapse_kernel import solve
    b = BENCHMARKS["phase2_literature"]["N_sat_FELA"]
    N_sat = float(b["value"])
    h_sat_D_printed = float(b["h_sat_over_D"])
    h_sat_D_exact = N_sat / 4.0
    D = 2.0 * params["l"]
    h_cr = float(solve(params, "active").h_cr)
    rows = []
    for hd in h_over_D_list:
        pp = dict(params, h=float(hd) * D)
        r = solve(pp, "active")
        rows.append(dict(
            h_over_D=float(hd), h_m=float(pp["h"]),
            q1a_kPa=float(r.q1a),
            globally_through=bool(r.is_global_collapse),
            l3_applicable=float(hd) <= 1.4 + 1e-12,
            note=("L1 闭式有效" if r.is_global_collapse
                  else "超出 h_cr，全高动员闭式外推不具物理意义"),
        ))
    rel = abs(h_sat_D_exact - h_sat_D_printed) / h_sat_D_printed
    return dict(
        N_sat=N_sat, h_sat_over_D_printed=h_sat_D_printed,
        h_sat_over_D_exact=h_sat_D_exact, consistency_rel_err=rel,
        consistency_ok=rel < b["tolerance_pct"] / 100.0,
        h_cr_m=h_cr, h_cr_over_D=h_cr / D,
        thresholds_distinct=abs(h_cr / D - h_sat_D_exact) > 0.1,
        threshold_note="h_sat 为力学饱和高度（拱效应层），h_cr 为几何贯通阈值（上限法）",
        rows=rows,
        source=b["source"],
        boundary_note="N_sat 为文献极限分析常数，本程序不推导、不拟合、不冒充严格上下界",
    )
