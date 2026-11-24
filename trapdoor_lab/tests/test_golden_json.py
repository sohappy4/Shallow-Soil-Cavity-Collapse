# -*- coding: utf-8 -*-
"""冻结 golden JSON 一致性测试（方案第七章要求）。

基准数值由 verify/make_golden.py 实跑生成并冻结，本测试只读取不改写；
公式若被有意改动，须人工复核后显式重跑生成脚本更新 JSON。
"""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from collapse_kernel import DEFAULT_PARAMS, solve, reproduce_tables
from l2_shooting import solve_shooting
from gtt import mechanism_scan
from benchmarks import saturation_anchor

with open(os.path.join(HERE, "golden_baseline.json"), encoding="utf-8") as f:
    GOLDEN = json.load(f)

B = GOLDEN["baseline"]


class TestFrozenGolden(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = solve(DEFAULT_PARAMS, "active")
        cls.t = reproduce_tables()
        cls.l2 = solve_shooting(DEFAULT_PARAMS, "active")
        cls.mech = mechanism_scan(DEFAULT_PARAMS, "active")
        cls.sat = saturation_anchor(DEFAULT_PARAMS)
        cls.lin = solve(dict(DEFAULT_PARAMS, B=1.0), "active")

    def _actual(self, key):
        return {
            "psi0_deg": self.r.psi0_deg,
            "psi_h_deg": self.r.psi_h_deg,
            "n": self.r.n,
            "h_cr_m": self.r.h_cr,
            "mirror_error": self.r.mirror_error,
            "q1_diff": self.r.q1b - self.r.q1a,
            "b1_h_cr_m": self.lin.h_cr,
            "table3_assoc_dev_pct": self.t["table3"]["assoc_dev_pct"],
            "table3_non_dev_pct": self.t["table3"]["non_dev_pct"],
            "table4_assoc_mean_dev_pct": self.t["table4"]["assoc_mean_dev_pct"],
            "table4_non_mean_dev_pct": self.t["table4"]["non_mean_dev_pct"],
            "l2_Lc_m": self.l2["L_c"],
            "gtt_sigma0_opt_kPa": self.mech["sigma0_opt"],
            "gtt_D_true": self.mech["D_true"],
            "sat_N_sat": self.sat["N_sat"],
            "sat_h_sat_over_D": self.sat["h_sat_over_D_exact"],
        }[key]

    def test_golden_params_match_default(self):
        self.assertEqual(GOLDEN["_meta"]["params"], DEFAULT_PARAMS)

    def test_all_frozen_values(self):
        actual = {k: self._actual(k) for k in B}
        for key, spec in B.items():
            with self.subTest(key=key):
                self.assertLessEqual(
                    abs(actual[key] - spec["value"]), spec["tol"],
                    "%s: 实跑 %.10g，冻结 %.10g，容差 %g"
                    % (key, actual[key], spec["value"], spec["tol"]))

    def test_golden_file_has_expected_fields(self):
        expected = {"psi0_deg", "psi_h_deg", "n", "h_cr_m", "mirror_error",
                    "q1_diff", "b1_h_cr_m", "table3_assoc_dev_pct",
                    "table3_non_dev_pct", "table4_assoc_mean_dev_pct",
                    "table4_non_mean_dev_pct", "l2_Lc_m",
                    "gtt_sigma0_opt_kPa", "gtt_D_true",
                    "sat_N_sat", "sat_h_sat_over_D"}
        self.assertEqual(set(B), expected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
