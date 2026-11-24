# -*- coding: utf-8 -*-
"""自检.py —— 命令行一键核对基准数值（不需要浏览器）。"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collapse_kernel import DEFAULT_PARAMS, solve, reproduce_tables
from l2_shooting import solve_shooting
from l3_rate import solve_rate
from gtt import check_table_13_2, mechanism_scan, tangent_params
from benchmarks import check_literature_anchors
from cases import selfcheck_cases, case_summary

r = solve(DEFAULT_PARAMS, "active")
t = reproduce_tables()
print("=" * 64)
print(" 浅埋土洞活门试验数值模拟程序 · 基准自检")
print("=" * 64)
print(" 地表剪胀角 psi0   = %.2f deg    (基准 37.56)" % r.psi0_deg)
print(" 底部剪胀角 psi(h) = %.2f deg    (基准 22.67)" % r.psi_h_deg)
print(" 埋深比 n = h/(2l) = %.4f        (基准 0.6000，不是 h/l)" % r.n)
print(" 临界覆土 h_cr     = %.3f m     (基准 5.389)" % r.h_cr)
print(" 镜面对称误差      = %.2e      (基准 0)" % r.mirror_error)
print(" 整体塌落判定      = %s" % ("是" if r.is_global_collapse else "否"))
print("-" * 64)
print(" Table 3 主动: assoc %+.1f%% / non-assoc %+.1f%%  (基准 -18.1 / +16.9)"
      % (t["table3"]["assoc_dev_pct"], t["table3"]["non_dev_pct"]))
print(" Table 4 被动: assoc %+.1f%% / non-assoc %+.1f%%  (基准 +25.4 / -9.5)"
      % (t["table4"]["assoc_mean_dev_pct"], t["table4"]["non_mean_dev_pct"]))
lin = solve(dict(DEFAULT_PARAMS, B=1.0), "active")
print(" B=1 线性分支: h_cr=%.4f m (理论 %.4f), 类型 %s"
      % (lin.h_cr, DEFAULT_PARAMS["l"] / DEFAULT_PARAMS["A"], lin.solution_type))
l2 = solve_shooting(DEFAULT_PARAMS, "active")
l2p = solve_shooting(DEFAULT_PARAMS, "passive")
l3 = solve_rate({})
print("-" * 64)
print(" L2 打靶: 迭代 %d 次, L_c=%.6f (闭式 %.6f, 误差 %.1e), 三级筛选 %s"
      % (l2["iterations"], l2["L_c"], l2["closed_form_Lc"], l2["lcf_error"],
         "全过" if l2["admissible"] else "未过"))
print(" L2 被动镜像: x2a+x2p=%.6f (2l=%.6f), 可容许 %s"
      % (l2["x2"] + l2p["x2"], 2 * DEFAULT_PARAMS["l"], l2p["admissible"]))
print(" L3 速率链: v=%.3f m/s, T=%.3f s, R=%.3f m, h/D=%.1f, 类型 %s"
      % (l3["v"], l3["T"], l3["R"], l3["h_over_D"], l3["solution_type"]))
gtt_chk = check_table_13_2()
mech = mechanism_scan(DEFAULT_PARAMS, "active")
tan0 = tangent_params(0.0, DEFAULT_PARAMS["sigma_c"], DEFAULT_PARAMS["A"],
                      DEFAULT_PARAMS["sigma_t"], 1.0)
print("-" * 64)
print(" GTT 表13.2: 逐项 %s; 原表印刷舍入差行 sigma_0=%s kPa"
      % ("通过" if gtt_chk["all_pass"] else "未过", gtt_chk["rounding_rows"]))
print(" GTT 定理13.2: sigma0_opt=%.6f, sigma_bar=%.6f kPa, 相对误差 %.2e"
      % (mech["sigma0_opt"], mech["sigma_bar"], mech["theorem13_2_rel_err"]))
print(" GTT 耗散夹逼: D_true=%.3f <= D_UB=%.3f, D_LB=%.3f <= D_true; B=1 c_t=%.3f"
      % (mech["D_true"], mech["D_upper"], mech["D_lower"], tan0["c_t"]))
lit = check_literature_anchors()
print("-" * 64)
print(" 二期文献锚点: 1g %d / 离心机 %d / 速率适用域 %d 项, 自洽核验 %s"
      % (lit["counts"]["one_g"], lit["counts"]["centrifuge"],
         lit["counts"]["rate_applicability"], "通过" if lit["all_pass"] else "未过"))
print("   beta_eff 区间半宽 %.2f%% (标注 ±11%%), 出露角 61°, N_sat=18.2212"
      % lit["beta_half_span_pct"])
csc = selfcheck_cases()
csu = case_summary()
uplift_ratio = [x["width_over_span"] for x in csu["width_over_span"]
                if x["id"] == "tunnel-uplift"]
print("-" * 64)
print(" M9 工程案例: 5 案例复算 %s, 反验 %s, 地表宽/洞跨 主动 %.2f–%.2f / 顶起 %.2f"
      % ("全部吻合" if csc["all_book_comparison_ok"] else "有偏差",
         "通过" if csc["selfcheck_Lc_at_hcr"] else "未过",
         csu["active_range"][0], csu["active_range"][1],
         uplift_ratio[0] if uplift_ratio else float("nan")))
print("=" * 64)
checks = [
    abs(r.psi0_deg - 37.56) < 0.01,
    abs(r.psi_h_deg - 22.67) < 0.01,
    r.mirror_error < 1e-10,
    abs(r.n - 0.6) < 1e-10,
    abs(t["table4"]["assoc_mean_dev_pct"] - 25.4) < 0.1,
    abs(lin.h_cr - DEFAULT_PARAMS["l"] / DEFAULT_PARAMS["A"]) < 1e-9,
    l2["converged"] and l2["admissible"] and l2["lcf_error"] < 1e-8,
    l2p["admissible"] and abs(l2["x2"] + l2p["x2"] - 2 * DEFAULT_PARAMS["l"]) < 1e-9,
    abs(l3["T"] - l3["params"]["h"] / l3["v"]) < 1e-9 and len(l3["stages"]) == 5,
    gtt_chk["all_pass"] and gtt_chk["rounding_rows"] == [50.0],
    mech["refined"] and mech["theorem13_2_rel_err"] < 1e-5
    and mech["upper_residual"] >= 0 and mech["lower_residual"] >= -1e-6,
    tan0["selfcheck_residual"] < 1e-10
    and abs(tan0["c_t"] - DEFAULT_PARAMS["A"] * DEFAULT_PARAMS["sigma_t"]) < 1e-10,
    lit["all_pass"] and lit["counts"]["one_g"] >= 10,
    csc["all_pass"] and csu["all_comparison_ok"] and csu["n_cases"] == 5,
]
print("自检结果: %s (%d/%d 通过)" % ("全部通过" if all(checks) else "存在失败",
                                   sum(checks), len(checks)))
sys.exit(0 if all(checks) else 1)
