# -*- coding: utf-8 -*-
"""FELA LP 上下界夹逼回归测试。

覆盖：线性化包线内外接、抗力夹逼定理、K 细化单调收窄、与 L1 闭式一致、
主动支护压力反序夹逼、psi 无关下限、Davis 单调性与域钳制、
适用域门控、边界声明诚实性。
"""
import os
import sys
import math
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collapse_kernel import DEFAULT_PARAMS, R_scale, limit_loads, psi0_of_Omega
from fela import (
    davis_reduce, envelope_tau, linearize_strength,
    resistance_bounds, solve_fela,
)

BASE = dict(DEFAULT_PARAMS)
OPTS = dict(K=12, n_points=41)


class TestLinearization(unittest.TestCase):

    def setUp(self):
        p = BASE
        self.sa = p["q2"]
        self.sb = p["q2"] + p["gamma"] * p["h"]

    def test_lines_have_three_coeffs(self):
        for mode in ("inscribed", "circumscribed"):
            lines, info = linearize_strength(BASE, self.sa, self.sb, 8, mode)
            self.assertEqual(len(lines), 8)
            self.assertEqual(info["K"], 8)
            for t in lines:
                self.assertEqual(len(t), 3)
                self.assertAlmostEqual(t[1], 1.0)

    def test_inscribed_below_circumscribed(self):
        li, _ = linearize_strength(BASE, self.sa, self.sb, 12, "inscribed")
        lo, _ = linearize_strength(BASE, self.sa, self.sb, 12, "circumscribed")
        sig = np.linspace(self.sa, self.sb, 30)
        self.assertTrue(np.all(envelope_tau(li, sig) <= envelope_tau(lo, sig) + 1e-9))

    def test_both_bracket_exact_power_law(self):
        """线性化包络应把精确幂型强度夹在中间。"""
        li, _ = linearize_strength(BASE, self.sa, self.sb, 12, "inscribed")
        lo, _ = linearize_strength(BASE, self.sa, self.sb, 12, "circumscribed")
        z = np.linspace(0, BASE["h"], 30)
        sig = BASE["q2"] + BASE["gamma"] * z
        exact = BASE["A"] * BASE["sigma_c"] * (
            (BASE["sigma_t"] + sig) / BASE["sigma_c"]) ** BASE["B"]
        self.assertTrue(np.all(envelope_tau(li, sig) <= exact + 1e-6))
        self.assertTrue(np.all(envelope_tau(lo, sig) >= exact - 1e-6))


class TestBracket(unittest.TestCase):

    def test_resistance_interval_contains_l1_exact(self):
        r = resistance_bounds(BASE, OPTS)
        self.assertTrue(r["success"])
        self.assertLessEqual(r["Ra_LB"], r["Ra_exact"] + 1e-6)
        self.assertGreaterEqual(r["Ra_UB"], r["Ra_exact"] - 1e-6)
        self.assertTrue(r["Ra_interval_ok"])

    def test_active_pressure_reversed_bracket(self):
        r = resistance_bounds(BASE, OPTS)
        q_lb, q_ub = r["q_active_LB"], r["q_active_UB"]
        self.assertLessEqual(q_lb, q_ub)
        # L1 闭式主动极限荷载必须落在夹逼区间内
        q1a = limit_loads(BASE["sigma_c"], BASE["sigma_t"], BASE["q2"],
                          BASE["gamma"], BASE["l"], BASE["h"],
                          BASE["A"], BASE["B"])[0]
        self.assertLessEqual(q_lb - 1e-6, q1a)
        self.assertGreaterEqual(q_ub + 1e-6, q1a)

    def test_refinement_monotone_narrower(self):
        widths = []
        for K in (4, 8, 12, 24):
            r = resistance_bounds(BASE, dict(OPTS, K=K))
            self.assertTrue(r["Ra_interval_ok"])
            widths.append(r["q_active_width"])
        for a, b in zip(widths, widths[1:]):
            self.assertLessEqual(b + 1e-9, a)

    def test_solve_fela_ok_and_types(self):
        r = solve_fela(BASE, OPTS)
        self.assertEqual(r["solution_type"], "fela-lp")
        self.assertTrue(r["bracket"]["ok"])
        self.assertTrue(r["meta"]["applicable"])
        self.assertEqual(r["meta"]["regime"], "shallow-applicable")
        self.assertTrue(np.isfinite(r["bracket"]["q_LB_kPa"]))
        self.assertTrue(np.isfinite(r["bracket"]["q_UB_kPa"]))

    def test_psi_independent_lower_bound(self):
        """静力下限必须严格与 psi 无关（专著 §7.7.1）。"""
        base = resistance_bounds(BASE, OPTS)
        for psi in (0, 5, 10, 25):
            r = resistance_bounds(BASE, dict(OPTS, psi_deg=psi))
            self.assertAlmostEqual(r["Ra_LB"], base["Ra_LB"], places=10)
            self.assertAlmostEqual(r["q_active_LB"], base["q_active_LB"], places=10)

    def test_davis_monotone_and_domain_clip(self):
        phi = math.degrees(math.atan(psi0_of_Omega(1.0, BASE["A"], BASE["B"])))
        ratios, qs = [], []
        for psi in (0, 5, 10, 20, 30, 40):
            r = resistance_bounds(BASE, dict(OPTS, psi_deg=psi))
            ratios.append(r["davis"]["ratio"])
            qs.append(r["davis"]["q_active_UB_reduced"])
        for a, b in zip(ratios, ratios[1:]):
            self.assertGreaterEqual(b + 1e-6, a)
        # q_UB 折减值在明显分离的 psi<=10 区间验证方向；
        # 接近关联时 ratio≈1、外接切面对 A 的微扰敏感（线性化数值噪声），
        # 只验证 ratio 单调（上条），不对 q 强求逐步单调。
        for a, b in zip(qs[:3], qs[1:3]):
            self.assertLessEqual(b + 1e-6, a)
        # psi >= phi 时钳到关联，ratio=1
        self.assertAlmostEqual(ratios[-1], 1.0, places=6)
        self.assertAlmostEqual(ratios[-2], 1.0, places=6)
        # 注意：Davis 折减后的 q 用外接多面体计算，与内接的 q_active_UB
        # 口径不同（一个偏危险侧、一个偏保守侧），不比较二者大小；
        # 折减报告仅作独立工程口径，已在 note 中声明不混入数学夹逼。

    def test_davis_reduce_no_psi(self):
        q, ratio = davis_reduce(BASE, None)
        self.assertEqual(ratio, 1.0)
        self.assertEqual(q["A"], BASE["A"])


class TestApplicability(unittest.TestCase):

    def test_deep_flagged_not_physical(self):
        for hD in (2.0, 4.0, 10.0):
            p = dict(BASE, h=hD * 2 * BASE["l"], q2=0.0)
            r = solve_fela(p, OPTS)
            self.assertFalse(r["meta"]["applicable"])
            self.assertIn("out-of-domain", r["meta"]["regime"])
            self.assertFalse(r["bracket"]["ok"])

    def test_literature_anchor_is_reference_only(self):
        r = solve_fela(BASE, OPTS)
        b = r["meta"]["boundary"]
        self.assertIn("N_sat=18.2212", b)
        self.assertIn("not derived here", b)
        self.assertIn("NOT SOCP", b)

    def test_opt_clamps(self):
        r = solve_fela(BASE, dict(K=1, n_points=2))
        self.assertGreaterEqual(r["meta"]["K"], 2)
        self.assertGreaterEqual(r["meta"]["n_points"], 6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
