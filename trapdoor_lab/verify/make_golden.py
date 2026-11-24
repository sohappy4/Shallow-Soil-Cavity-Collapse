# -*- coding: utf-8 -*-
"""生成 golden_baseline.json —— 从当前内核实跑值写入，避免手抄数字出错。

运行：  ../venv/bin/python verify/make_golden.py
生成：  tests/golden_baseline.json
注意：本脚本只用于「冻结基线」。若内核公式被有意修改，须经人工复核后
      显式重跑本脚本更新 JSON，不允许测试自动改写基准。
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from collapse_kernel import DEFAULT_PARAMS, solve, reproduce_tables
from l2_shooting import solve_shooting
from gtt import mechanism_scan
from benchmarks import saturation_anchor

r = solve(DEFAULT_PARAMS, "active")
t = reproduce_tables()
l2 = solve_shooting(DEFAULT_PARAMS, "active")
mech = mechanism_scan(DEFAULT_PARAMS, "active")
sat = saturation_anchor(DEFAULT_PARAMS)
lin = solve(dict(DEFAULT_PARAMS, B=1.0), "active")

baseline = {
    "_meta": {
        "purpose": "一期/二期冻结基准；由 verify/make_golden.py 实跑生成，勿手工编辑",
        "params": DEFAULT_PARAMS,
        "generated_by": "trapdoor_lab/verify/make_golden.py",
    },
    "baseline": {
        "psi0_deg": {"value": r.psi0_deg, "tol": 0.01},
        "psi_h_deg": {"value": r.psi_h_deg, "tol": 0.01},
        "n": {"value": r.n, "tol": 1e-12},
        "h_cr_m": {"value": r.h_cr, "tol": 0.01},
        "mirror_error": {"value": r.mirror_error, "tol": 1e-10},
        "q1_diff": {"value": r.q1b - r.q1a, "tol": 0.01,
                    "note": "应等于 2*Ra/l"},
        "b1_h_cr_m": {"value": lin.h_cr, "tol": 1e-9,
                      "note": "线性极限 h_cr=l/A，不与 L1 闭式解比"},
        "table3_assoc_dev_pct": {"value": t["table3"]["assoc_dev_pct"], "tol": 0.1},
        "table3_non_dev_pct": {"value": t["table3"]["non_dev_pct"], "tol": 0.1},
        "table4_assoc_mean_dev_pct": {"value": t["table4"]["assoc_mean_dev_pct"], "tol": 0.1},
        "table4_non_mean_dev_pct": {"value": t["table4"]["non_mean_dev_pct"], "tol": 0.1},
        "l2_Lc_m": {"value": l2["L_c"], "tol": 1e-8,
                    "note": "L2 打靶与闭式解相对误差基准 6.9e-12"},
        "gtt_sigma0_opt_kPa": {"value": mech["sigma0_opt"], "tol": 0.01,
                               "note": "定理13.2：应等于滑面正应力加权均值"},
        "gtt_D_true": {"value": mech["D_true"], "tol": 0.01},
        "sat_N_sat": {"value": sat["N_sat"], "tol": 1e-9,
                      "note": "文献常数，非本程序推导"},
        "sat_h_sat_over_D": {"value": sat["h_sat_over_D_exact"], "tol": 1e-9},
    },
}

out = os.path.join(ROOT, "tests", "golden_baseline.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(baseline, f, ensure_ascii=False, indent=2)
    f.write("\n")
print("已写入", out)
