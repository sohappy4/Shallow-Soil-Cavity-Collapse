# -*- coding: utf-8 -*-
"""GTT 广义切线法回归测试：切线层、表 13.2、定理 13.1/13.2、耗散夹逼。"""
import os
import sys
import math
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collapse_kernel import DEFAULT_PARAMS
from gtt import (
    tangent_params, psi_from_sigma, scan_tangent, table_13_2, check_table_13_2,
    chord_lower_params, chord_tau, dissipation_bracket, mechanism_scan,
)


class TestTangent(unittest.TestCase):
    def test_table_13_1_identity(self):
        """定理 13.1：切点倾角 delta_t 恒等于式 (13.12) 的 psi。"""
        for B in (0.4, 0.5, 0.8):
            for s in (0.0, 12.0, 50.0, 130.0):
                r = tangent_params(s, 100.0, 0.6325, 10.0, B)
                self.assertAlmostEqual(r["delta_deg"], psi_from_sigma(s, 100.0, 0.6325, 10.0, B), places=10)
                self.assertLess(r["selfcheck_residual"], 1e-10)
                self.assertLess(r["c_t_formula_residual"], 1e-10)

    def test_b1_mohr_coulomb_limit(self):
        """B=1 应退化为 MC：tan(delta)=A、c_t=A*sigma_t、与 sigma_0 无关。"""
        a = tangent_params(0.0, 100.0, 0.6325, 10.0, 1.0)
        b = tangent_params(80.0, 100.0, 0.6325, 10.0, 1.0)
        self.assertAlmostEqual(a["tan_delta"], 0.6325, places=12)
        self.assertAlmostEqual(a["c_t"], 0.6325 * 10.0, places=12)
        self.assertAlmostEqual(a["c_t"], b["c_t"], places=12)
        self.assertAlmostEqual(a["delta_deg"], b["delta_deg"], places=12)

    def test_scan_residual(self):
        rows = scan_tangent(100.0, 0.6325, 10.0, 0.5, 0.0, 100.0, n=25)
        self.assertEqual(len(rows), 25)
        self.assertTrue(all(r["theorem13_1_residual"] < 1e-9 for r in rows))

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            tangent_params(-20.0, 100.0, 0.6325, 10.0, 0.5)
        with self.assertRaises(ValueError):
            tangent_params(10.0, 100.0, 0.6325, 10.0, 1.2)
        with self.assertRaises(ValueError):
            chord_lower_params(100.0, 0.6325, 10.0, 0.5, [10.0])


class TestTable132(unittest.TestCase):
    def test_table_structure(self):
        d = table_13_2()
        self.assertEqual(len(d["rows"]), 6)
        self.assertIn("selfcheck_residual", d["rows"][0])

    def test_printed_values_with_disclosed_rounding(self):
        r = check_table_13_2()
        self.assertTrue(r["all_pass"])
        self.assertFalse(r["strict_pass"])
        # 原表 B=0.5, sigma_0=50 行的印刷舍入差必须被显式记录，不能静默吞掉
        self.assertEqual(r["rounding_rows"], [50.0])
        row = r["results"][1]
        self.assertAlmostEqual(row["deviation"]["c_t"], 0.0194, places=3)
        self.assertAlmostEqual(row["deviation"]["selfcheck"], 0.0332, places=3)


class TestBracket(unittest.TestCase):
    def test_tangent_upper_chord_lower(self):
        b = dissipation_bracket([0.0, 15.0, 50.0, 120.0], 100.0, 0.6325, 10.0, 0.5, 40.0)
        self.assertGreaterEqual(b["upper_residual"], -1e-10)
        self.assertGreaterEqual(b["lower_residual"], -1e-10)
        self.assertGreater(b["max_gap"], 0.0)

    def test_b1_all_identical(self):
        b = dissipation_bracket([0.0, 15.0, 50.0, 120.0], 100.0, 0.6325, 10.0, 1.0, 40.0)
        self.assertLess(b["max_gap"], 1e-10)

    def test_chord_segments_are_below_envelope(self):
        import numpy as np
        nodes = np.linspace(0.0, 100.0, 11)
        segs = chord_lower_params(100.0, 0.6325, 10.0, 0.5, nodes)
        self.assertEqual(len(segs), 10)
        xs = np.linspace(0.0, 100.0, 501)
        low = chord_tau(xs, 100.0, 0.6325, 10.0, 0.5, nodes)
        true = 0.6325 * 100.0 ** 0.5 * (xs + 10.0) ** 0.5
        self.assertTrue(np.all(low <= true + 1e-10))


class TestTheorem132(unittest.TestCase):
    def test_optimal_tangent_equals_mean_normal_stress(self):
        for mode in ("active", "passive"):
            r = mechanism_scan(DEFAULT_PARAMS, mode)
            self.assertTrue(r["refined"])
            self.assertLess(r["theorem13_2_rel_err"], 1e-5)
            self.assertGreaterEqual(r["upper_residual"], 0.0)
            self.assertGreaterEqual(r["lower_residual"], -1e-6)
            self.assertIn("非 FELA 荷载上下界", r["boundary_note"])

    def test_b1_scan_is_flat(self):
        r = mechanism_scan(dict(DEFAULT_PARAMS, B=1.0), "active")
        self.assertEqual(r["theorem13_2_rel_err"], 0.0)
        self.assertLess(r["upper_residual"], 1e-9)
        self.assertLess(r["lower_residual"], 1e-9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
