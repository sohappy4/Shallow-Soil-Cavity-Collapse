# -*- coding: utf-8 -*-
"""L2/L3 + 基准集回归测试。"""
import os
import sys
import math
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collapse_kernel import DEFAULT_PARAMS, solve, reproduce_tables
from collapse_kernel import slip_profile
from l2_shooting import solve_shooting
from l3_rate import solve_rate, mu_of_I
from benchmarks import check_benchmarks, BENCHMARKS


class TestL2(unittest.TestCase):
    def test_active_converges_and_matches_closed_form(self):
        r = solve_shooting(DEFAULT_PARAMS, "active")
        self.assertTrue(r["converged"])
        self.assertLess(r["lcf_error"], 1e-8)
        self.assertLess(abs(r["residual"]), 1e-8)
        self.assertTrue(r["admissible"])

    def test_passive_is_pure_mirror(self):
        a = solve_shooting(DEFAULT_PARAMS, "active")
        p = solve_shooting(DEFAULT_PARAMS, "passive")
        self.assertTrue(p["admissible"])
        self.assertAlmostEqual(p["L_c"], a["L_c"], places=9)
        self.assertAlmostEqual(p["x2"] + a["x2"], 2 * DEFAULT_PARAMS["l"], places=9)

    def test_closed_form_cross_check(self):
        # RK4 打靶与闭式解网格点数不同（401 vs 400），在公共区间插值后逐点比对曲线
        import numpy as np
        r = solve_shooting(DEFAULT_PARAMS, "active")
        xc, fc = slip_profile(DEFAULT_PARAMS, "active")
        xr, fr = np.array(r["x"]), np.array(r["F"])
        lo = max(xr.min(), xc.min())
        hi = min(xr.max(), xc.max())
        xs = np.linspace(lo, hi, 200)
        fr_i = np.interp(xs, xr, fr)
        fc_i = np.interp(xs, xc, fc)
        # 闭式解本身只有 400 点网格，线性插值带来 O(dx^2)≈3e-6 m 的误差（约 3 微米），
        # 不随 RK4 步数下降，属网格离散而非求解分歧；L_c 的解析吻合已达 6.9e-12。
        maxerr = float(np.max(np.abs(fr_i - fc_i)))
        self.assertLess(maxerr, 1e-5)

    def test_derivative_outcrop(self):
        r = solve_shooting(DEFAULT_PARAMS, "active")
        p = r["p"]
        self.assertAlmostEqual(p, 1.929, places=2)

    def test_deep_no_through_going(self):
        r = solve_shooting(dict(DEFAULT_PARAMS, h=8.0), "active")
        self.assertTrue(r["no_through_going"])
        self.assertAlmostEqual(r["h_cr_implied"], 5.389, places=2)

    def test_b1_rejected(self):
        with self.assertRaises(ValueError):
            solve_shooting(dict(DEFAULT_PARAMS, B=1.0), "active")


class TestMuOfI(unittest.TestCase):
    def test_limit(self):
        self.assertAlmostEqual(float(mu_of_I(0.0, 0.6, 0.7, 0.3)), 0.6, places=10)
        self.assertAlmostEqual(float(mu_of_I(1e12, 0.6, 0.7, 0.3)), 0.7, places=6)


class TestL3(unittest.TestCase):
    def test_chain(self):
        r = solve_rate({})
        self.assertGreater(r["v"], 0)
        self.assertAlmostEqual(r["T"], r["params"]["h"] / r["v"], places=9)
        self.assertAlmostEqual(r["R"],
            (0.75 * r["params"]["B_c"] ** 2 * r["params"]["h"]
             / math.tan(math.radians(r["params"]["phi_r"]))) ** (1 / 3), places=9)

    def test_beverloo_v(self):
        r = solve_rate({})
        g = r["params"]["g"]
        eb = r["params"]["B_c"] - r["params"]["k_d"] * r["params"]["d"]
        v = 4 * r["params"]["C_b"] * math.sqrt(g) * eb ** 2.5 / (math.pi * r["params"]["B_c"] ** 2)
        self.assertAlmostEqual(r["v"], v, places=9)

    def test_five_stages(self):
        r = solve_rate({})
        self.assertEqual(len(r["stages"]), 5)
        self.assertEqual([s["stage"] for s in r["stages"]], [1, 2, 3, 4, 5])

    def test_domain_degradation(self):
        shallow = solve_rate(dict(h=1.2, B_c=1.0))
        deep = solve_rate(dict(h=6.0, B_c=1.0))
        self.assertEqual(shallow["solution_type"], "rate")
        self.assertEqual(deep["solution_type"], "rate-degraded")
        self.assertIn("超界", deep["domain"])

    def test_zero_cohesion_never_self_stable(self):
        r = solve_rate(dict(tau0=0.0))
        self.assertFalse(r["self_stable"])
        self.assertIn("永不自稳", r["stability_note"])


class TestBenchmarks(unittest.TestCase):
    def test_phase1_all_pass(self):
        t = reproduce_tables()
        r = solve(DEFAULT_PARAMS, "active")
        res = check_benchmarks(r.psi0_deg, r.psi_h_deg, r.h_cr, r.mirror_error,
                               t["table3"], t["table4"])
        self.assertTrue(res["all_pass"], res["checks"])

    def test_phase2_anchors_present(self):
        p2 = BENCHMARKS["phase2_literature"]
        self.assertAlmostEqual(p2["beta_eff_centrifuge"]["value_min"], 4.17, places=2)
        self.assertAlmostEqual(p2["beta_eff_centrifuge"]["value_max"], 5.26, places=2)
        self.assertAlmostEqual(p2["N_sat_FELA"]["value"], 18.2212, places=4)
        self.assertAlmostEqual(p2["N_sat_FELA"]["h_sat_over_D"], 4.555, places=3)
        self.assertAlmostEqual(p2["outcrop_angle_1g_deg"]["value"], 61.0, places=1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
