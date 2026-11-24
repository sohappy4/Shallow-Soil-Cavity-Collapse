# -*- coding: utf-8 -*-
"""一期回归基准：全部数值来自 cavity_collapse_calculator.py 实跑 + 论文 Table 3/4。"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collapse_kernel import DEFAULT_PARAMS, solve, figure_data, reproduce_tables, ParamsError


class TestBaseline(unittest.TestCase):
    def setUp(self):
        self.r = solve(DEFAULT_PARAMS, "active")

    def test_psi_surface(self):
        self.assertAlmostEqual(self.r.psi0_deg, 37.56, places=2)

    def test_psi_bottom(self):
        self.assertAlmostEqual(self.r.psi_h_deg, 22.67, places=2)

    def test_mirror_symmetry(self):
        self.assertLessEqual(self.r.mirror_error, 1e-10)
        self.assertAlmostEqual(self.r.q1b - self.r.q1a, 2 * self.r.Ra / DEFAULT_PARAMS["l"], places=9)

    def test_burial_ratio_definition(self):
        # 修正点：n = h/(2l)，不是 h/l（后者会翻倍）
        self.assertAlmostEqual(self.r.n, DEFAULT_PARAMS["h"] / (2 * DEFAULT_PARAMS["l"]), places=12)
        self.assertAlmostEqual(self.r.n, 0.6, places=12)

    def test_critical_depth(self):
        self.assertAlmostEqual(self.r.h_cr, 5.389, places=2)

    def test_global_collapse_active(self):
        self.assertTrue(self.r.is_global_collapse)  # h=3 < h_cr=5.389
        deep = solve(dict(DEFAULT_PARAMS, h=8.0), "active")
        self.assertFalse(deep.is_global_collapse)
        self.assertIn("不会发生整体贯通塌落", deep.judgement_text)


class TestLinearLimit(unittest.TestCase):
    """初稿 M2 验收错误在此修正：B=1 与线性解比对，不与 L1 闭式解比对。"""

    def test_b1_critical_depth(self):
        p = dict(DEFAULT_PARAMS, B=1.0)
        r = solve(p, "active")
        self.assertEqual(r.solution_type, "linear-limit")
        self.assertAlmostEqual(r.h_cr, p["l"] / p["A"], places=10)
        self.assertAlmostEqual(r.psi0_deg, __import__("math").degrees(__import__("math").atan(p["A"])), places=10)

    def test_b1_no_nan(self):
        p = dict(DEFAULT_PARAMS, B=1.0)
        r = solve(p, "active")
        for k in ("h_cr", "q1a", "q1b", "Ra", "mirror_error"):
            self.assertTrue(__import__("math").isfinite(getattr(r, k)), k)


class TestPassiveMode(unittest.TestCase):
    def test_passive_has_no_hcr_gate(self):
        r = solve(DEFAULT_PARAMS, "passive")
        self.assertFalse(r.is_global_applicable)
        self.assertIn("不设 h_cr 门控", r.judgement_text)


class TestValidation(unittest.TestCase):
    def test_invalid_b(self):
        with self.assertRaises(ParamsError):
            solve(dict(DEFAULT_PARAMS, B=1.5), "active")
        with self.assertRaises(ParamsError):
            solve(dict(DEFAULT_PARAMS, B=0.0), "active")

    def test_invalid_sigma_c(self):
        with self.assertRaises(ParamsError):
            solve(dict(DEFAULT_PARAMS, sigma_c=0.0), "active")

    def test_invalid_mode(self):
        with self.assertRaises(ParamsError):
            solve(DEFAULT_PARAMS, "whatever")


class TestReproduction(unittest.TestCase):
    def test_table3(self):
        t3 = reproduce_tables()["table3"]
        self.assertAlmostEqual(t3["assoc_Ga"], 0.385, places=3)
        self.assertAlmostEqual(t3["assoc_dev_pct"], -18.1, places=1)
        self.assertAlmostEqual(t3["non_Ga"], 0.55, places=2)
        self.assertAlmostEqual(t3["non_dev_pct"], 16.9, places=1)

    def test_table4(self):
        t4 = reproduce_tables()["table4"]
        self.assertAlmostEqual(t4["assoc_mean_dev_pct"], 25.4, places=1)
        self.assertAlmostEqual(t4["non_mean_dev_pct"], -9.5, places=1)


class TestFigureData(unittest.TestCase):
    def test_shapes(self):
        f = figure_data(DEFAULT_PARAMS, "active")
        self.assertEqual(len(f["fig1_slip_surface"]["active_x"]), 400)
        self.assertEqual(len(f["fig2_stress_field"]["z"]), 200)
        self.assertEqual(len(f["fig4_zoning"]["zone_omega"]), 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
