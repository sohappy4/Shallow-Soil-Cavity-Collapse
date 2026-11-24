# -*- coding: utf-8 -*-
"""文献饱和高度锚点回归测试：只核验自洽与门槛区分，不声称推导 N_sat。"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collapse_kernel import DEFAULT_PARAMS
from benchmarks import BENCHMARKS, saturation_anchor


class TestSaturationAnchor(unittest.TestCase):
    def test_constants_self_consistent(self):
        r = saturation_anchor(DEFAULT_PARAMS)
        self.assertTrue(r["consistency_ok"])
        self.assertAlmostEqual(r["h_sat_over_D_exact"], 18.2212 / 4.0, places=10)
        self.assertLess(r["consistency_rel_err"], 0.001)

    def test_h_sat_and_h_cr_are_different_thresholds(self):
        r = saturation_anchor(DEFAULT_PARAMS)
        self.assertTrue(r["thresholds_distinct"])
        self.assertLess(r["h_cr_over_D"], 2.0)
        self.assertGreater(r["h_sat_over_D_exact"], 4.0)
        self.assertIn("不可混用", BENCHMARKS["phase2_literature"]["h_sat_cr_distinction"]["type"])

    def test_domain_flags_and_extrapolation_note(self):
        r = saturation_anchor(DEFAULT_PARAMS)
        by = {x["h_over_D"]: x for x in r["rows"]}
        self.assertTrue(by[1.0]["l3_applicable"])
        self.assertTrue(by[1.4]["l3_applicable"])
        self.assertFalse(by[2.0]["l3_applicable"])
        self.assertTrue(by[2.0]["globally_through"] is False)
        self.assertIn("不具物理意义", by[4.555]["note"])
        self.assertIn("不推导", r["boundary_note"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
