# -*- coding: utf-8 -*-
"""主动/被动镜像对照回归测试（防止插值方向错误重现）。"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import Handler
from collapse_kernel import DEFAULT_PARAMS, figure_data


class TestMirror(unittest.TestCase):
    def test_power_law_mirror(self):
        s = figure_data(DEFAULT_PARAMS, "active")["fig1_slip_surface"]
        r = Handler._mirror_check(s)
        self.assertTrue(r["mirror_ok"], r["max_mirror_error_m"])
        self.assertLess(r["max_mirror_error_m"], 1e-5)
        self.assertAlmostEqual(r["x2_sum"], r["two_l"], places=9)

    def test_linear_limit_mirror(self):
        # 基准 h=3、l=2.5、dx=h/A=5>l 时主动出露被裁剪到 x=0，
        # 几何上不再与被动支镜像；取 h=0.6（dx=1<l）核验未裁剪情形。
        p = dict(DEFAULT_PARAMS, B=1.0, h=0.6)
        s = figure_data(p, "active")["fig1_slip_surface"]
        r = Handler._mirror_check(s)
        self.assertTrue(r["mirror_ok"], r["max_mirror_error_m"])
        self.assertAlmostEqual(r["x2_sum"], r["two_l"], places=9)

    def test_linear_limit_clipped_case_is_flagged(self):
        # 深埋线性分支主动出露被裁剪到 x=0，接口应给出说明而不是假装对称
        s = figure_data(dict(DEFAULT_PARAMS, B=1.0), "active")["fig1_slip_surface"]
        r = Handler._mirror_check(s)
        self.assertIn("裁剪", r["note"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
