# -*- coding: utf-8 -*-
"""fem2d.py 二维弹塑性 FEM-SRF 求解器回归测试。

覆盖：DP 本构自检、径向返回一致性、拉力截断、平面应变弹性、
小网格求解稳定性、荷载—位移曲线物理趋势、塑性区发展、
强度折减效应、对称性、场量有界性。
"""
import os
import sys
import math
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fem2d import (
    DEFAULT_FEM_PARAMS, build_mesh, dp_params, dp_radial_return,
    dp_yield, selfcheck_dp, solve_fem,
)


SMALL = dict(nx=8, ny=6, n_steps=8)
BASE = dict(DEFAULT_FEM_PARAMS, l=2.5, h=3.0)


class TestDPConstitutive(unittest.TestCase):

    def test_selfcheck_all_pass(self):
        sc = selfcheck_dp()
        self.assertTrue(sc["all_pass"], str(sc["checks"]))
        self.assertLess(abs(sc["f_after"]), 1e-6)

    def test_radial_return_maps_to_surface(self):
        E, nu, c, phi = 20000.0, 0.3, 25.0, 20.0
        m, k, mp = dp_params(c, phi, phi)
        for de in (np.array([0.0, 0.0, 0.02]),
                   np.array([-0.02, -0.01, 0.06]),   # 剪切主导，必屈服
                   np.array([0.02, -0.03, -0.02])):
            sig, dg, yld = dp_radial_return(np.zeros(4), de, E, nu, m, k, mp)
            self.assertTrue(yld)
            self.assertGreater(dg, 0.0)
            self.assertLess(abs(dp_yield(sig, m, k)), 1e-6)

    def test_small_strain_stays_elastic(self):
        E, nu, c, phi = 20000.0, 0.3, 25.0, 20.0
        m, k, mp = dp_params(c, phi, phi)
        sig, dg, yld = dp_radial_return(
            np.zeros(4), np.array([1e-7, 0.0, 0.0]), E, nu, m, k, mp)
        self.assertFalse(yld)
        self.assertAlmostEqual(dg, 0.0)

    def test_tension_cutoff_no_explosion(self):
        """DP 锥无拉断帽：净拉应力下必须靠拉力截断保持有界。"""
        E, nu, c, phi = 20000.0, 0.3, 25.0, 20.0
        m, k, mp = dp_params(c, phi, 0.0)   # 非关联，m_psi=0
        sig, dg, yld = dp_radial_return(
            np.zeros(4), np.array([0.03, 0.0, 0.0]), E, nu, m, k, mp)
        self.assertTrue(np.all(np.isfinite(sig)))
        self.assertLess(np.max(np.abs(sig)), 1e4)
        self.assertTrue(yld)

    def test_dp_params_ucs_positive(self):
        m, k, mp = dp_params(25.0, 20.0, 10.0)
        self.assertGreater(m, 0.0)
        self.assertGreater(k, 0.0)
        self.assertGreater(mp, 0.0)
        self.assertLess(mp, m)


class TestMesh(unittest.TestCase):

    def test_mesh_shape_and_ordering(self):
        nodes, elems = build_mesh(8, 6, -5.0, 15.0, 0.0, 3.0)
        self.assertEqual(nodes.shape, (63, 2))
        self.assertEqual(elems.shape, (48, 4))
        self.assertEqual(len(set(elems.ravel().tolist())), 63)
        # 首单元逆时针
        e = elems[0]
        self.assertTrue(nodes[e[0], 0] < nodes[e[1], 0])
        self.assertTrue(nodes[e[3], 1] > nodes[e[0], 1])  # 单元逆时针：d 在 a 上方


class TestFemSolver(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = solve_fem(BASE, opts=SMALL)

    def test_converged_and_finite(self):
        r = self.result
        self.assertTrue(r["meta"]["converged"])
        self.assertTrue(np.isfinite(r["peak_q_kPa"]))
        self.assertGreater(r["peak_q_kPa"], 0.0)
        self.assertLess(r["peak_q_kPa"], 1e4)

    def test_curve_structure(self):
        r = self.result
        self.assertEqual(len(r["curve"]), SMALL["n_steps"] + 1)
        self.assertEqual(r["curve"][0]["step"], 0)
        self.assertEqual(r["curve"][0]["door_down"], 0.0)
        for c in r["curve"]:
            for key in ("q_kPa", "reaction_kN", "plastic_frac",
                        "iters", "residual"):
                self.assertTrue(math.isfinite(c[key]), "%s 非有限" % key)

    def test_arching_pressure_drop(self):
        """土拱效应：初始静止压力显著高于下沉后残余压力。"""
        q = [c["q_kPa"] for c in self.result["curve"]]
        q0, qres = q[0], float(np.mean(q[3:]))
        self.assertGreater(q0, 0.0)
        self.assertGreater(qres, 0.0)
        self.assertGreater(q0, 2.0 * qres)

    def test_plastic_zone_grows(self):
        pf = [c["plastic_frac"] for c in self.result["curve"]]
        self.assertAlmostEqual(pf[0], 0.0)
        self.assertGreater(pf[-1], pf[1])
        self.assertLessEqual(pf[-1], 1.0)

    def test_field_bounded_and_consistent(self):
        f = self.result["field"]
        n = self.result["meta"]["n_elements"]
        for key in ("sxx", "syy", "txy", "p", "q", "epsp"):
            arr = np.asarray(f[key], dtype=float)
            self.assertEqual(len(arr), n)
            self.assertTrue(np.all(np.isfinite(arr)), key)
            self.assertLess(np.max(np.abs(arr)), 1e5, key)
        self.assertEqual(len(f["centers"]), n)
        self.assertEqual(len(f["yielded"]), n)

    def test_symmetry(self):
        """几何/荷载关于洞中线 x=D/2 对称，塑性区应近似镜像。"""
        r = self.result
        x = np.asarray([c[0] for c in r["field"]["centers"]])
        y = np.asarray([c[1] for c in r["field"]["centers"]])
        ep = np.asarray(r["field"]["epsp"])
        # 按行（同一 y）分组统计左右对称性
        errs = []
        for yy in np.unique(np.round(y, 6)):
            sel = np.abs(y - yy) < 1e-6
            xs, es = x[sel], ep[sel]
            mid = 0.5 * r["params"]["D"]
            left = es[xs < mid]
            right = es[xs > mid]
            if len(left) and len(right):
                errs.append(abs(left.mean() - right.mean()))
        self.assertLess(max(errs), 0.15)

    def test_strength_reduction_increases_damage(self):
        """强度折减 Fs 增大时塑性区占比不应减少。"""
        r1 = solve_fem(BASE, opts=dict(SMALL, sr=1.0))
        r2 = solve_fem(BASE, opts=dict(SMALL, sr=1.4))
        self.assertLessEqual(r1["field"]["plastic_frac"] + 1e-9,
                             r2["field"]["plastic_frac"])
        self.assertLessEqual(r2["meta"]["c_eff"], r1["meta"]["c_eff"])
        self.assertLess(r2["meta"]["phi_eff_deg"], r1["meta"]["phi_eff_deg"])

    def test_meta_boundary_honest(self):
        r = self.result
        self.assertEqual(r["solution_type"], "fem-srf")
        self.assertIn("NOT rigorous FELA", r["meta"]["boundary"])
        self.assertIn("external-mature-FEM", r["meta"]["source"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
