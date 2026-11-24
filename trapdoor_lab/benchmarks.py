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

__all__ = [
    "BENCHMARKS", "check_benchmarks", "saturation_anchor",
    "LITERATURE_ANCHORS", "check_literature_anchors",
]

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

# ---------------------------------------------------------------------------
# 二期 1g / 离心机量化锚点集
# ---------------------------------------------------------------------------
# 数据全部逐条摘自专著送审版正文（§10.3–10.5、§11、§18、§20.4），
# 只录入"作者给出了明确数值或区间"的锚点，不虚构逐点试验曲线。
# 每个锚点的 ref 字段指向专著章节，source 指向原始文献。
LITERATURE_ANCHORS: Dict[str, Any] = {
    "one_g": [
        {
            "key": "mcnulty_displacement_sensitivity",
            "label": "McNulty 位移极敏感性",
            "value": 0.0002, "unit": "倍活门直径 D",
            "effect": "活门荷载改变多达 50%",
            "source": "McNulty (1965), An Experimental Study of Arching in Sand [50]",
            "ref": "专著 §10.3.1",
            "kind": "threshold",
        },
        {
            "key": "mcnulty_h_over_D_threshold",
            "label": "McNulty h/D 门槛",
            "value": 2.0, "unit": "h/D",
            "effect": "h/D<2 曲线形状随 h/D 显著变化；h/D>2 几乎不再变化",
            "source": "McNulty (1965) [50]",
            "ref": "专著 §10.3.1，与 Li 等 (2024) 覆盖比 2 分界呼应",
            "kind": "threshold",
        },
        {
            "key": "mcnulty_equal_settlement_axisym",
            "label": "轴对称等沉面高度",
            "value_min": 1.0, "value_max": 1.5, "unit": "倍活门直径 D",
            "effect": "该面以上无拱效应、无竖向应力折减；低于平面应变 1.5–2.5D",
            "source": "McNulty (1965) / Bhandari & Han [46]",
            "ref": "专著 §10.3.1",
            "kind": "range",
        },
        {
            "key": "mcnulty_equal_settlement_plane",
            "label": "平面应变等沉面高度",
            "value_min": 1.5, "value_max": 2.5, "unit": "倍活门宽度",
            "effect": "平面应变高于轴对称",
            "source": "Bhandari & Han [46]（转述太沙基/McNulty）",
            "ref": "专著 §10.3.1",
            "kind": "range",
        },
        {
            "key": "mcnulty_displacement_to_steady_axisym",
            "label": "轴对称位移稳定门槛",
            "value": 3.0, "unit": "% 活门直径",
            "effect": "压力约在 3% 位移处趋于常数（平面应变约 10%）",
            "source": "McNulty (1965) [46]",
            "ref": "专著 §10.3.1",
            "kind": "threshold",
        },
        {
            "key": "mcnulty_displacement_to_steady_plane",
            "label": "平面应变位移稳定门槛",
            "value": 10.0, "unit": "% 活门宽度",
            "effect": "太沙基平面应变约 10% 位移处压力趋稳",
            "source": "Terzaghi (1936) / Bhandari & Han [46]",
            "ref": "专著 §10.3.1",
            "kind": "threshold",
        },
        {
            "key": "abdulla_outcrop_angle",
            "label": "Abdulla 沉陷坑出露角",
            "value": 61.0, "unit": "deg",
            "effect": "主动塌落沉陷坑出露角，本书形态谱的定量佐证",
            "source": "Abdulla & Goodings (1996) [17]",
            "ref": "专著 §2.2.1 / §10.5",
            "kind": "value",
        },
        {
            "key": "han_dilatancy",
            "label": "Han 实测剪胀角",
            "value": 8.0, "unit": "deg",
            "effect": "0.36 m 活门、相对密度 70% 河砂，PIV 测得剪切带倾角即 psi",
            "source": "Han et al. (2017) [55]",
            "ref": "专著 §10.4.1",
            "kind": "value",
        },
        {
            "key": "liang_dilatancy",
            "label": "Liang 实测剪胀角",
            "value": 7.5, "unit": "deg",
            "effect": "0.25 m 活门、相对密度 45% 标准砂，PIV 测得",
            "source": "Liang et al. (2020) [56]",
            "ref": "专著 §10.4.1",
            "kind": "value",
        },
        {
            "key": "li_cover_ratio",
            "label": "Li 覆盖比分界",
            "value": 2.0, "unit": "覆盖比 h/W",
            "effect": "裂缝贯穿/不贯穿分界（h/W=2 对应 h/l≈4）",
            "source": "Li et al. (2024) [47]",
            "ref": "专著 §10.4 / §10.5.1",
            "kind": "threshold",
        },
    ],
    "centrifuge": [
        {
            "key": "dewoolkar_beta_eff_min",
            "label": "Dewoolkar beta_eff 下限",
            "value": 4.17, "unit": "无量纲",
            "effect": "密砂 150g 系列，h/D=2–6 全程不漂移",
            "source": "Dewoolkar, Santichaianant & Ko (2007) [20]",
            "ref": "专著 §18.4 / §20.4.2",
            "kind": "value",
        },
        {
            "key": "dewoolkar_beta_eff_max",
            "label": "Dewoolkar beta_eff 上限",
            "value": 5.26, "unit": "无量纲",
            "effect": "同上；离散约 ±11%",
            "source": "Dewoolkar et al. (2007) [20]",
            "ref": "专著 §18.4 / §20.4.2",
            "kind": "value",
        },
        {
            "key": "dewoolkar_scatter",
            "label": "Dewoolkar 数据离散",
            "value": 11.0, "unit": "%",
            "effect": "beta_eff 全程离散幅度",
            "source": "Dewoolkar et al. (2007) [20]",
            "ref": "专著 §20.4.2",
            "kind": "value",
        },
        {
            "key": "dewoolkar_plateau_h_over_D",
            "label": "Dewoolkar 荷载平台门槛",
            "value": 5.0, "unit": "h/D",
            "effect": "最小门上荷载在 h/D≥5 后进入恒定平台（与上覆压力、模型高度无关）",
            "source": "Dewoolkar et al. (2007) [20]",
            "ref": "专著 §11 / 图 11.4",
            "kind": "threshold",
        },
    ],
    "rate_applicability": [
        {
            "key": "l3_applicability",
            "label": "L3 速率模型严格适用域",
            "value": 1.4, "unit": "h/D",
            "effect": "h/D≲1.4 相对偏差不超过 5%",
            "source": "专著 §18（速率过程层适用域）",
            "ref": "专著 §18 / §20.4.2",
            "kind": "threshold",
        },
        {
            "key": "l3_dev_at_hd2",
            "label": "L3 在 h/D=2 的偏差",
            "value": 14.0, "unit": "%",
            "effect": "h/D=2 时全高动员模型偏危险约 14%",
            "source": "专著速率层适用域讨论",
            "ref": "专著 §20.4.2",
            "kind": "value",
        },
        {
            "key": "l3_dev_at_hd10",
            "label": "L3 在 h/D=10 的偏差",
            "value": 174.0, "unit": "%",
            "effect": "h/D=10 时偏差达 174%，须转静力下限判据或 FELA 曲线",
            "source": "专著速率层适用域讨论",
            "ref": "专著 §20.4.2",
            "kind": "value",
        },
        {
            "key": "self_stable_tol",
            "label": "自稳判据数值吻合度",
            "value": 2.0, "unit": "%",
            "effect": "无黏性土永不自稳；数值极限分析对浅埋域自稳阈值吻合到 ±2% 以内",
            "source": "专著 §20.4.2 引数值极限分析",
            "ref": "专著 §20.4.2",
            "kind": "value",
        },
        {
            "key": "N_sat_consistency",
            "label": "N_sat 四位重合",
            "value": 18.2212, "unit": "无量纲",
            "effect": "本书 N_sat=18.2212 与深埋圆柱孔静力下界 N=18.22 四位重合",
            "source": "Cox et al. (1961) [94] / 本书 §20.4.2",
            "ref": "专著 §18.2.2 / §20.4.2 文献 [94]",
            "kind": "value",
        },
    ],
}


def check_literature_anchors() -> Dict[str, Any]:
    """核验二期文献锚点集的结构完整性与关键数值自洽。

    注意：这些锚点是"文献报告值"，本程序不重新拟合试验曲线，也不声称
    从 L1/L2/L3 推导出它们。本函数做的是可复现的机械核验：
      (1) 分组齐全、每组非空；
      (2) 每个锚点必备字段齐全、数值有限；
      (3) 关键换算自洽：beta_eff 区间宽度与离散度、N_sat/4=h_sat/D、
          L3 适用域门槛与 saturation_anchor 使用值一致。
    """
    import math

    groups = ("one_g", "centrifuge", "rate_applicability")
    errors: list = []
    for g in groups:
        rows = LITERATURE_ANCHORS.get(g)
        if not isinstance(rows, list) or not rows:
            errors.append("分组 '%s' 缺失或为空" % g)
            continue
        for a in rows:
            for field in ("key", "label", "unit", "source", "ref", "kind"):
                if not a.get(field):
                    errors.append("%s.%s 缺字段 %s" % (g, a.get("key", "?"), field))
            vals = [a[k] for k in ("value", "value_min", "value_max") if k in a]
            if not vals:
                errors.append("%s.%s 无数值" % (g, a.get("key", "?")))
            for v in vals:
                if not isinstance(v, (int, float)) or not math.isfinite(v):
                    errors.append("%s.%s 数值非有限" % (g, a.get("key", "?")))

    def by_key(g, k):
        for a in LITERATURE_ANCHORS.get(g, []):
            if a.get("key") == k:
                return a
        return {}

    bmin = by_key("centrifuge", "dewoolkar_beta_eff_min").get("value")
    bmax = by_key("centrifuge", "dewoolkar_beta_eff_max").get("value")
    scatter = by_key("centrifuge", "dewoolkar_scatter").get("value")
    beta_span_ok = False
    if None not in (bmin, bmax, scatter) and bmin:
        # 区间中点相对半宽的百分比，应量级接近所列离散度（±11%）
        mid = (bmin + bmax) / 2.0
        half_pct = (bmax - bmin) / 2.0 / mid * 100.0
        beta_span_ok = abs(half_pct - scatter) <= 2.0

    nsat = by_key("rate_applicability", "N_sat_consistency").get("value")
    hsat_printed = BENCHMARKS["phase2_literature"]["N_sat_FELA"]["h_sat_over_D"]
    # 与 saturation_anchor 同口径：印刷值 4.555 为三位小数，用 1% 相对容差
    nsat_ok = (nsat is not None and
               abs(nsat / 4.0 - hsat_printed) / hsat_printed < 0.01)

    l3_dom = by_key("rate_applicability", "l3_applicability").get("value")
    # saturation_anchor 内部以 1.4 作为 h/D 适用门槛，须与锚点一致
    l3_dom_ok = (l3_dom is not None and abs(l3_dom - 1.4) < 1e-9)

    checks = {
        "groups_complete": not any("分组" in e for e in errors),
        "fields_finite": not any(("缺字段" in e) or ("无数值" in e) or ("非有限" in e) for e in errors),
        "beta_span_consistent": beta_span_ok,
        "Nsat_over_4_eq_hsat": nsat_ok,
        "l3_domain_consistent": l3_dom_ok,
    }
    return {
        "checks": checks,
        "all_pass": all(checks.values()) and not errors,
        "errors": errors,
        "counts": {g: len(LITERATURE_ANCHORS.get(g, [])) for g in groups},
        "beta_half_span_pct": (
            (bmax - bmin) / 2.0 / ((bmin + bmax) / 2.0) * 100.0
            if None not in (bmin, bmax) and bmin else None
        ),
        "note": "文献报告值，程序不拟合、不推导；仅做机械自洽核验",
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
