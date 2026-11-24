# -*- coding: utf-8 -*-
"""五个工程案例模板与复算引擎（方案 M9）。

案例来源：专著《浅埋土洞坍塌》第 19 章（§19.2–19.5），参数与数值结果逐条对齐。
**所有案例参数均为典型取值示例，非特定场地实测值**，仅用于演示判别流程、
公式用法与形态判读；实际工程须按第 10 章场地试验测定后代入，不得直接套用。

复算公式（专著）：
  式(16.32) A_ψ = tanψ₀ / (B·Ω^(B−1))，Ω=(σ_t+q₂)/σ_c
  式(16.33) S_ψ = σ_c·(A_ψ·B)^(B/(1−B))·(A−A_ψ·B)/(1−B)
  式(16.19a) L_c = (S_ψ/γ)·{(p+γh/(B·S_ψ))^B − p^B}，p=(cotψ₀)^(1/(1−B))
  式(19.1)  h_cr = (B·S_ψ/γ)·{(γl/S_ψ+p^B)^(1/B) − p}
  被动顶起（关联流动）：Δ 与 L_c 同式结构，q₁ 按式(19.14) 逐项求和。

界值身份：关联解为严格运动学上界；以实测 ψ 反定 S_ψ 的解是**工程近似**
（非严格界）；静力下限独立校核见 fela.py。深埋 h>h_cr 时 L_c=l 驻值解失效，
按专著 §16.5 改用静力场判定。
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

__all__ = [
    "ENGINEERING_CASES", "CASE_TEMPLATES", "evaluate_case",
    "evaluate_all", "case_summary", "selfcheck_cases",
]

# ---------------------------------------------------------------------------
# 五案例模板（参数逐条对齐专著第19章；expected 为原书公开数值）
# ---------------------------------------------------------------------------
ENGINEERING_CASES: List[Dict[str, Any]] = [
    dict(
        id="karst-cavity",
        name="案例一：覆盖型岩溶土洞塌陷（两土性对照）",
        source="专著 §19.2.1 / 表19.1 / 式(19.2)–(19.6)",
        scenario="岩溶区路基下隐伏土洞，洞跨 4 m、覆土 6 m；同一几何荷载下两种土性对照。",
        mode="active",
        params=dict(sigma_c=500.0, sigma_t=50.0, A=0.6, B=0.5,
                    gamma=19.0, q2=20.0, l=2.0, h=6.0, psi0_deg=12.0),
        expected=dict(Omega=0.140, A_psi=0.159, S_psi=41.39, p=22.13,
                      L_c=1.205, x2=0.795, h_cr=10.33, surface_width=1.59),
        morphology="硬塑黏土瓮形 / 松散粉砂近井筒",
        judgement="硬塑黏土 A：口小 1.59 m、肚大 4 m 的瓮形，存在贯通机构；"
                  "松散粉砂 B：ψ<0 且 σ_t+q₂>0，闭合支不可容许（式19.6 越界 7.67 kPa），"
                  "取 ψ→0 同时性极限（近井筒，地表宽 4 m）为下限、§17.6 渐进模型为上包络。",
        caution="参数为典型示例，非场地实测；粉砂 σ_t=3 kPa 为表观拉伸截距（毛细吸力），"
                "纯无黏粉砂须按 B=1、σ_t=0 重算。",
        variants=[
            dict(label="硬塑黏土 A（瓮形）", mode="active",
                 params=dict(sigma_c=500.0, sigma_t=50.0, A=0.6, B=0.5,
                             gamma=19.0, q2=20.0, l=2.0, h=6.0, psi0_deg=12.0),
                 expected=dict(Omega=0.140, A_psi=0.159, S_psi=41.39, p=22.13,
                               L_c=1.205, x2=0.795, h_cr=10.33, surface_width=1.59)),
            dict(label="松散粉砂 B（近井筒，渐进）", mode="progressive",
                 params=dict(sigma_c=120.0, sigma_t=3.0, A=0.5, B=0.5,
                             gamma=19.0, q2=20.0, l=2.0, h=6.0, psi0_deg=-3.0),
                 expected=dict(Omega=0.192, surface_width=4.00)),
        ],
    ),
    dict(
        id="urban-road",
        name="案例二：城市道路塌陷（管线渗漏空洞）",
        source="专著 §19.2.2 / 式(19.7)",
        scenario="管线渗漏形成空洞，洞跨 4 m、覆土 6 m，路面交通荷载折算 q₂=20 kPa；弱剪胀 ψ=8°。",
        mode="active",
        params=dict(sigma_c=500.0, sigma_t=50.0, A=0.6, B=0.5,
                    gamma=19.0, q2=20.0, l=2.0, h=6.0, psi0_deg=8.0),
        expected=dict(A_psi=0.1052, S_psi=28.79, p=50.63, L_c=0.81,
                      x2=1.19, h_cr=15.55, surface_width=2.37),
        morphology="瓮形（小开口、大深度）",
        judgement="h_cr=15.55 m ≫ h=6 m，贯通机构存在；地表裂缝仅洞宽 59%，"
                  "不得因缝窄判危害小。以空腔充填+注浆加固为主，路面沉降观测+地质雷达复扫。",
        caution="关联流动 ψ=38.7° 会给出 h_cr=2.92 m<h，误判为稳定拱、错过预警——务必实测 ψ。",
    ),
    dict(
        id="goaf",
        name="案例三：采空区覆岩垮落（窄条带采空）",
        source="专著 §19.3.1 / 式(19.8)–(19.9)",
        scenario="窄条带采空区，采空宽 16 m、覆岩 40 m，煤系地层 γ=22，ψ=8°。",
        mode="active",
        params=dict(sigma_c=800.0, sigma_t=60.0, A=0.55, B=0.6,
                    gamma=22.0, q2=0.0, l=8.0, h=40.0, psi0_deg=8.0),
        expected=dict(Omega=0.075, A_psi=0.0831, S_psi=11.14, p=135.05,
                      L_c=4.85, x2=3.15, h_cr=71.55, surface_width=6.31),
        morphology="瓮形（尺度效应显著）",
        judgement="h_cr=71.55 m ≫ h=40 m，贯通机构存在；Λ=0.220 远大于案例一，"
                  "大跨度采空地表沉陷不可避免，与长壁全垮落法工程认知一致。",
        caution="A_ψ、S_ψ 对 B 敏感；大 B 采空区须以现场岩体试验标定。",
    ),
    dict(
        id="shield-tunnel",
        name="案例四：盾构掘进地层损失（中密砂）",
        source="专著 §19.3.2 / 式(19.10)–(19.11)",
        scenario="盾构开挖直径 6 m、覆土 15 m，中密砂 γ=18，地表超载 30 kPa，ψ=10°。",
        mode="active",
        params=dict(sigma_c=300.0, sigma_t=8.0, A=0.55, B=0.5,
                    gamma=18.0, q2=30.0, l=3.0, h=15.0, psi0_deg=10.0),
        expected=dict(Omega=0.127, A_psi=0.126, S_psi=18.35, p=32.16,
                      L_c=2.22, x2=0.78, h_cr=21.43, surface_width=1.56),
        morphology="瓮形（地表影响 0.26 倍洞径）",
        judgement="h_cr=21.43 m > h=15 m，地层损失将反映至地表；地表影响半宽 0.78 m。",
        caution="σ_t=8 kPa 为表观拉伸截距；Peck 沉降槽宽度 i≈K_s·z₀ 为经验式，本书 x₂ 为解析解。",
    ),
    dict(
        id="tunnel-uplift",
        name="案例五：隧道底板顶起（被动工况）",
        source="专著 §19.4 / 表19.2 / 式(19.12)–(19.14)",
        scenario="隧道半宽 3 m、底板至地表 12 m，承压水/高地应力顶起 u<0；关联流动。",
        mode="passive",
        params=dict(sigma_c=600.0, sigma_t=60.0, A=0.6, B=0.5,
                    gamma=20.0, q2=30.0, l=3.0, h=12.0, psi0_deg=None),
        expected=dict(S=108.00, p=1.667, Delta=6.378, x2=9.378, q1=1109.51),
        morphology="外倾漏斗（地表影响 3.13 倍洞跨）",
        judgement="顶起无 h_cr 约束、恒贯通地表；耗散仍为张开支，勿套闭合支。"
                  "q₁≈1.11 MPa，与第16章表16.4 交叉验证一致。",
        caution="u·s<0 锁定外倾 s=+1；与主动塌落几何相反但耗散支相同。",
    ),
]

# 便于按 id 取用
CASE_TEMPLATES: Dict[str, Dict[str, Any]] = {c["id"]: c for c in ENGINEERING_CASES}


# ---------------------------------------------------------------------------
# 复算引擎
# ---------------------------------------------------------------------------
def _active_calc(p: Dict[str, float]) -> Dict[str, Any]:
    """主动工况：非关联修正 + L_c / x₂ / h_cr（专著式16.32/16.33/16.19a/19.1）。"""
    sc, st, A, B = p["sigma_c"], p["sigma_t"], p["A"], p["B"]
    g, q2, l, h = p["gamma"], p["q2"], p["l"], p["h"]
    psi0 = p["psi0_deg"]
    Omega = (st + q2) / sc
    A_psi = math.tan(math.radians(psi0)) / (B * Omega ** (B - 1.0))
    S_psi = sc * (A_psi * B) ** (B / (1.0 - B)) * (A - A_psi * B) / (1.0 - B)
    pparam = (1.0 / math.tan(math.radians(psi0))) ** (1.0 / (1.0 - B))
    Lc = S_psi / g * ((pparam + g * h / (B * S_psi)) ** B - pparam ** B)
    hcr = B * S_psi / g * ((g * l / S_psi + pparam ** B) ** (1.0 / B) - pparam)
    x2 = l - Lc
    # 自洽反验：把 h_cr 代回 L_c 公式应等于 l（专著要求每次必做）
    Lc_at_cr = S_psi / g * ((pparam + g * hcr / (B * S_psi)) ** B - pparam ** B)
    return dict(
        Omega=Omega, A_psi=A_psi, S_psi=S_psi, p=pparam,
        L_c=Lc, x2=x2, h_cr=hcr,
        surface_width=2.0 * x2,
        global_collapse=bool(h <= hcr),
        branch_valid=bool(x2 >= 0.0),
        selfcheck_Lc_at_hcr=Lc_at_cr,
        selfcheck_error=abs(Lc_at_cr - l),
        nonassociated_reduction=float(S_psi / (sc * A * A * B)),
    )


def _passive_calc(p: Dict[str, float]) -> Dict[str, Any]:
    """被动顶起：关联流动，Δ / x₂ / q₁（专著式19.12–19.14）。"""
    sc, st, A, B = p["sigma_c"], p["sigma_t"], p["A"], p["B"]
    g, q2, l, h = p["gamma"], p["q2"], p["l"], p["h"]
    S = sc * A * A * B
    Omega = (st + q2) / sc
    pparam = 4.0 * Omega / (A * A)
    Delta = S / g * ((pparam + g * h / (B * S)) ** 0.5 - math.sqrt(pparam))
    x2 = l + Delta
    q1 = (st * Delta + 0.5 * S * pparam * Delta +
          g * math.sqrt(pparam) * Delta ** 2 +
          g * g * Delta ** 3 / (3.0 * S) +
          g * l * h + q2 * x2) / l
    return dict(
        Omega=Omega, S=S, p=pparam, Delta=Delta, x2=x2,
        surface_width=2.0 * x2, q1=q1,
        global_collapse=True,          # 顶起恒贯通，无 h_cr 约束
        branch_valid=True,
        h_cr=None,
        selfcheck_Lc_at_hcr=None,
        selfcheck_error=None,
    )


def _progressive_calc(p: Dict[str, float]) -> Dict[str, Any]:
    """渐进破坏（ψ<0）：取 ψ→0 同时性极限为下限，报告渐进上包络待场地标定。"""
    # ψ→0 时滑动面近竖直，L_c→0，x₂→l
    l = p["l"]
    Omega = (p["sigma_t"] + p["q2"]) / p["sigma_c"]
    # 式(19.6) 越界诊断（专著：整体同时破坏不可容许的判据）
    st, B, q2 = p["sigma_t"], p["B"], p["q2"]
    sigma_x2 = -(2.0 * st + (1.0 - B) * q2) / (1.0 + B)
    exceed = max(0.0, -sigma_x2 - st)
    return dict(
        Omega=Omega,
        L_c=0.0, x2=l, surface_width=2.0 * l,
        progressive_lower=l,            # 同时性 ψ→0 极限（下限）
        progressive_upper=None,          # §17.6 渐进上包络需 l_eff(h)，本模板不虚构
        closure_sigma_at_surface=sigma_x2,
        closure_exceed_kPa=exceed,
        closure_inadmissible=bool(exceed > 0.0),
        global_collapse=True,
        branch_valid=False,              # 闭合支不可容许
        h_cr=None,
        selfcheck_Lc_at_hcr=None,
        selfcheck_error=None,
    )


def evaluate_case(case: Dict[str, Any]) -> Dict[str, Any]:
    """复算单个案例，返回 {case, params, result, comparison, boundary}。"""
    mode = case.get("mode", "active")
    p = dict(case["params"])
    if mode == "active":
        calc = _active_calc(p)
    elif mode == "passive":
        calc = _passive_calc(p)
    elif mode == "progressive":
        calc = _progressive_calc(p)
    else:
        raise ValueError("未知工况：%r" % mode)

    # 与原书数值逐项对照（只在 expected 中给出的键比较）
    exp = case.get("expected", {})
    comparison = {}
    for key, val in exp.items():
        if key == "morphology":
            continue
        if key in calc and calc[key] is not None:
            comparison[key] = dict(
                computed=float(calc[key]),
                book=float(val),
                abs_error=float(abs(calc[key] - val)),
                rel_error=float(abs(calc[key] - val) / max(abs(val), 1e-12)),
                ok=bool(abs(calc[key] - val) <= max(0.02, 0.01 * abs(val))),
            )
    # 同一案例的多种土性对照（如案例一 A/B），逐一复算
    variants = []
    for v in case.get("variants", []):
        vm = v.get("mode", mode)
        if vm == "active":
            vc = _active_calc(v["params"])
        elif vm == "passive":
            vc = _passive_calc(v["params"])
        elif vm == "progressive":
            vc = _progressive_calc(v["params"])
        else:
            raise ValueError("未知变体工况：%r" % vm)
        vcomp = {}
        for key, val in v.get("expected", {}).items():
            if key in vc and vc[key] is not None:
                vcomp[key] = dict(
                    computed=float(vc[key]), book=float(val),
                    abs_error=float(abs(vc[key] - val)),
                    rel_error=float(abs(vc[key] - val) / max(abs(val), 1e-12)),
                    ok=bool(abs(vc[key] - val) <= max(0.02, 0.01 * abs(val))))
        variants.append(dict(label=v["label"], mode=vm, params=v["params"],
                             result=vc, comparison=vcomp,
                             comparison_ok=all(x["ok"] for x in vcomp.values())))

    boundary = ("案例参数与数值结果来自专著第19章，均为典型取值示例、非场地实测；"
                "以实测 ψ 反定 S_ψ 的解是工程近似而非严格界值；"
                "关联解为运动学上界，静力下限另由 fela.py 校核。")
    return dict(
        case=dict(id=case["id"], name=case["name"], source=case["source"],
                  scenario=case["scenario"], mode=mode,
                  morphology=case.get("morphology", ""),
                  judgement=case.get("judgement", ""),
                  caution=case.get("caution", "")),
        params=p,
        result=calc,
        comparison=comparison,
        comparison_ok=all(v["ok"] for v in comparison.values()),
        variants=variants,
        all_comparison_ok=(all(v["ok"] for v in comparison.values()) and
                           all(v["comparison_ok"] for v in variants)),
        boundary=boundary,
    )


def evaluate_all() -> List[Dict[str, Any]]:
    """复算全部五案例。"""
    return [evaluate_case(c) for c in ENGINEERING_CASES]


def case_summary() -> Dict[str, Any]:
    """汇总：全部复算通过、地表宽/洞跨统计、关键警示。"""
    results = evaluate_all()
    ratios = []
    for r in results:
        p, res = r["params"], r["result"]
        ratios.append(dict(id=r["case"]["id"], mode=r["case"]["mode"],
                           width_over_span=res["surface_width"] / (2.0 * p["l"])))
        for v in r["variants"]:
            ratios.append(dict(id=r["case"]["id"] + "/" + v["label"],
                               mode=v["mode"],
                               width_over_span=(v["result"]["surface_width"] /
                                                (2.0 * v["params"]["l"]))))
    # 按工况分流统计，避免把渐进近井筒变体（宽/跨=1.00）混入“主动瓮形”区间
    active = [x for x in ratios if x["mode"] == "active"]
    progressive = [x for x in ratios if x["mode"] == "progressive"]
    uplift = [x for x in ratios if x["mode"] == "passive"]
    return dict(
        n_cases=len(results),
        all_comparison_ok=all(r["all_comparison_ok"] for r in results),
        selfcheck_ok=all(
            (r["result"].get("selfcheck_error") is None or
             r["result"]["selfcheck_error"] < 1e-6) for r in results),
        width_over_span=ratios,
        active_range=[min(x["width_over_span"] for x in active),
                      max(x["width_over_span"] for x in active)],
        progressive_range=[x["width_over_span"] for x in progressive],
        uplift_range=[x["width_over_span"] for x in uplift],
        key_warning=("地表影响宽度/洞跨：主动瓮形约 0.26–0.59、松散粉砂近井筒 1.00、"
                     "顶起约 3.13；地表坑小不等于地下空腔小，按坑口反推空洞会系统性低估风险。"),
        caveat="案例参数为典型示例，非特定场地实测值；实际工程须按场地试验测定后代入。",
    )


def selfcheck_cases() -> Dict[str, Any]:
    """自检：案例数量/字段完整、五案例数值复算、反验、警示声明。"""
    checks: Dict[str, Any] = {}
    checks["case_count"] = (len(ENGINEERING_CASES) == 5)
    checks["unique_ids"] = (len(CASE_TEMPLATES) == 5)
    checks["required_fields"] = all(
        all(k in c for k in ("id", "name", "source", "scenario", "mode",
                             "params", "expected", "morphology", "judgement"))
        for c in ENGINEERING_CASES)
    results = evaluate_all()
    checks["all_book_comparison_ok"] = all(r["all_comparison_ok"] for r in results)
    # 专著强制反验：L_c(h_cr) = l
    actives = [r["result"] for r in results if r["result"].get("selfcheck_error") is not None]
    checks["selfcheck_Lc_at_hcr"] = all(x["selfcheck_error"] < 1e-6 for x in actives)
    # 关联退化：A_ψ=A 时 S_ψ=S（用案例一参数验证式16.34 退化性）
    p = dict(ENGINEERING_CASES[0]["params"])
    a = _active_calc(p)
    checks["nonassociated_reduces_strength"] = (a["A_psi"] < p["A"] and
                                                a["S_psi"] < p["sigma_c"] * p["A"] ** 2 * p["B"])
    # 顶起恒贯通、无 h_cr
    uplift = results[-1]["result"]
    checks["passive_no_hcr"] = (uplift["h_cr"] is None and uplift["global_collapse"])
    # 渐进闭合支不可容许（案例一第二个土性变体）
    prog = results[0]["variants"][1]["result"]
    checks["progressive_inadmissible"] = (prog["closure_inadmissible"] and
                                          prog["closure_exceed_kPa"] > 7.0)
    checks["boundary_honest"] = all("非场地实测" in r["case"]["caution"] or
                                    "典型" in r["case"]["caution"] or
                                    "非严格界值" in r["boundary"] for r in results)
    checks["all_pass"] = all(v for k, v in checks.items() if k != "all_pass")
    return checks


if __name__ == "__main__":
    import json
    print(json.dumps(selfcheck_cases(), indent=2, ensure_ascii=False))
    print(json.dumps(case_summary(), indent=2, ensure_ascii=False))
