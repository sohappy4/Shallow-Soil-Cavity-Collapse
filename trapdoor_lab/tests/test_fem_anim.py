# -*- coding: utf-8 -*-
"""FEM 应力场动态云图（方案 M5）与演示动图（方案 M6）回归测试。"""
import os
import sys
import json
import base64
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fem2d import solve_fem, DEFAULT_FEM_PARAMS
from fem_anim import build_fem_animation

BASE = dict(DEFAULT_FEM_PARAMS, l=2.5, h=3.0)
SMALL = dict(nx=6, ny=6, n_steps=6)
RESULT = solve_fem(BASE, opts=SMALL)
ANIM = build_fem_animation(RESULT)


class TestFieldHistory(unittest.TestCase):

    def test_history_length_and_keys(self):
        h = RESULT["field_history"]
        self.assertEqual(len(h), SMALL["n_steps"] + 1)
        self.assertEqual(h[0]["step"], 0)
        self.assertEqual(h[-1]["step"], SMALL["n_steps"])
        for f in h:
            for k in ("step", "door_down", "q", "epsp", "yielded",
                      "plastic_frac"):
                self.assertIn(k, f)

    def test_history_json_serializable(self):
        # 逐帧场全部为 Python 标量，不能直接含 numpy 类型
        json.dumps(RESULT["field_history"])

    def test_history_matches_field(self):
        last = RESULT["field_history"][-1]
        field = RESULT["field"]
        np.testing.assert_allclose(last["q"], field["q"], rtol=1e-12)
        np.testing.assert_allclose(last["epsp"], field["epsp"], rtol=1e-12)
        self.assertEqual(last["yielded"], field["yielded"])

    def test_plastic_zone_grows_monotonically(self):
        pf = [f["plastic_frac"] for f in RESULT["field_history"]]
        self.assertAlmostEqual(pf[0], 0.0)
        for a, b in zip(pf, pf[1:]):
            self.assertGreaterEqual(b + 1e-12, a)

    def test_pressure_drops_from_initial(self):
        q = [c["q_kPa"] for c in RESULT["curve"]]
        self.assertGreater(q[0], 2.0 * q[-1])


class TestFemAnimation(unittest.TestCase):

    def test_frames_are_png(self):
        fr = ANIM["frames"]
        self.assertEqual(len(fr), SMALL["n_steps"] + 1)
        for s in fr:
            raw = base64.b64decode(s)
            self.assertEqual(raw[:8], b"\x89PNG\r\n\x1a\n")

    def test_gif_animated(self):
        raw = base64.b64decode(ANIM["gif"])
        self.assertIn(raw[:6], (b"GIF87a", b"GIF89a"))
        # 多帧 GIF 含有多个图像分隔符
        self.assertGreater(raw.count(b"\x2c"), 2)

    def test_meta_bounds_and_honesty(self):
        m = ANIM["meta"]
        self.assertEqual(m["n_frames"], SMALL["n_steps"] + 1)
        self.assertGreater(m["fps"], 0)
        self.assertGreaterEqual(m["qmax"], max(
            max(f["q"]) for f in RESULT["field_history"]))
        self.assertIn("非严格 FELA", m["note"])
        self.assertIn("NOT rigorous FELA", m["boundary"])

    def test_uniform_color_scale(self):
        # qmax 取全场最大，而非仅末帧最大
        qmax_all = max(max(f["q"]) for f in RESULT["field_history"])
        qmax_last = max(RESULT["field_history"][-1]["q"])
        self.assertGreaterEqual(ANIM["meta"]["qmax"], qmax_all)
        if qmax_last < qmax_all:
            self.assertGreater(ANIM["meta"]["qmax"], qmax_last)


if __name__ == "__main__":
    unittest.main(verbosity=2)
