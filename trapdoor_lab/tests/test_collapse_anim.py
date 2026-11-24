# -*- coding: utf-8 -*-
"""二维坍塌运动学动画回归测试。

核验：帧结构完整、时序单调、粒子运动学约束（前锋未扫过不动）、
数值有限、边界说明如实（非 DEM/应力场）、被动镜像可选。
"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collapse_kernel import DEFAULT_PARAMS
from collapse_anim import build_animation


class TestCollapseAnimation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = build_animation(DEFAULT_PARAMS, opts={"n_frames": 10})

    def test_structure_and_meta(self):
        a = self.a
        self.assertEqual(a["meta"]["solution_type"], "kinematic-animation")
        self.assertEqual(len(a["frames"]), 10)
        self.assertGreater(a["meta"]["n_particles"], 0)
        for key in ("l", "h", "D", "x_outcrop", "x_floor", "R_final"):
            self.assertIn(key, a["geometry"])
        self.assertIn("不是 DEM", a["meta"]["boundary_note"])
        self.assertIn("有限元", a["meta"]["boundary_note"])

    def test_frame_fields(self):
        for f in self.a["frames"]:
            for key in ("i", "t", "z_front", "door_down", "crater_r",
                        "stage", "surface", "slip", "particles"):
                self.assertIn(key, f)
            self.assertGreater(len(f["surface"]), 2)
            self.assertGreater(len(f["slip"]), 2)
            self.assertGreater(len(f["particles"]), 0)

    def test_monotonic_time_front_crater(self):
        fr = self.a["frames"]
        ts = [f["t"] for f in fr]
        zs = [f["z_front"] for f in fr]
        cr = [f["crater_r"] for f in fr]
        dd = [f["door_down"] for f in fr]
        for seq in (ts, zs, cr, dd):
            self.assertTrue(all(seq[i] <= seq[i + 1] for i in range(len(seq) - 1)))
        self.assertAlmostEqual(ts[0], 0.0, places=9)
        self.assertAlmostEqual(ts[-1], self.a["meta"]["T"], places=5)
        self.assertLess(abs(zs[-1] - DEFAULT_PARAMS["h"]), 1e-6)
        self.assertLess(abs(cr[-1] - self.a["meta"]["R_final"]), 1e-7)

    def test_first_frame_no_motion(self):
        f0 = self.a["frames"][0]
        self.assertEqual(f0["z_front"], 0.0)
        self.assertEqual(f0["door_down"], 0.0)
        self.assertEqual(f0["crater_r"], 0.0)
        self.assertTrue(all(p[2] == 0 for p in f0["particles"]))

    def test_particles_move_only_after_front(self):
        a = self.a
        z0 = [p[1] for p in a["frames"][0]["particles"]]
        h = a["geometry"]["h"]
        for f in a["frames"]:
            zf = f["z_front"]
            for i, p in enumerate(f["particles"]):
                swept = z0[i] - (h - zf)
                if swept <= 0 and z0[i] >= 0.14 * h:
                    # 未扫过且不在地表漏斗带：应静止
                    self.assertEqual(p[2], 0)
                    self.assertAlmostEqual(p[1], z0[i], places=4)

    def test_drops_increase_with_depth(self):
        a = self.a
        f = a["frames"][-1]
        z0 = [p[1] for p in a["frames"][0]["particles"]]
        bands = []
        for lo, hi in ((0, 0.75), (1.5, 2.25), (2.25, 3.0)):
            ds = [f["particles"][i][1] - z0[i]
                  for i, z in enumerate(z0) if lo <= z < hi]
            if ds:
                bands.append(float(np.mean(ds)))
        self.assertGreaterEqual(len(bands), 3)
        self.assertTrue(all(bands[i] <= bands[i + 1] + 1e-6
                            for i in range(len(bands) - 1)))

    def test_all_coordinates_finite(self):
        for f in self.a["frames"]:
            vals = []
            for p in f["particles"]:
                vals.extend(p[:2])
            for s in f["surface"]:
                vals.extend(s)
            for s in f["slip"]:
                vals.extend(s)
            self.assertTrue(bool(np.all(np.isfinite(vals))))

    def test_frame_count_clamped(self):
        small = build_animation(DEFAULT_PARAMS, opts={"n_frames": 2})
        big = build_animation(DEFAULT_PARAMS, opts={"n_frames": 999})
        self.assertEqual(len(small["frames"]), 8)   # 下限
        self.assertEqual(len(big["frames"]), 90)    # 上限

    def test_passive_mirror_option(self):
        ap = build_animation(DEFAULT_PARAMS, opts={"n_frames": 6, "passive": True})
        self.assertTrue(ap["meta"]["passive"])
        self.assertEqual(len(ap["frames"]), 8)  # n_frames 下限为 8


if __name__ == "__main__":
    unittest.main(verbosity=2)
