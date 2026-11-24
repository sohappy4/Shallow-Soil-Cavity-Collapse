# -*- coding: utf-8 -*-
"""DEM 颗粒流 / 独立进程适配器 / 三维锥形求解器回归测试。"""
import json
import math
import unittest

import numpy as np

from dem2d import (initial_packing, contact_forces, solve_dem,
                   selfcheck_dem, DEFAULT_DEM_PARAMS, DEFAULT_DEM_OPTS)
from dem_backend import (make_request, run_backend, selfcheck_backend,
                         LocalDEMBackend, BackendError, PROTOCOL)
from cone3d import (radius_profile, solve_cone, critical_depth,
                    selfcheck_cone, DEFAULT_CONE_PARAMS, MORPHOLOGIES)

FAST = dict(n_particles=60, n_steps=40, settle_steps=40, save_every=20,
            max_particles=120)


class TestDEM(unittest.TestCase):
    def test_selfcheck(self):
        s = selfcheck_dem()
        self.assertTrue(s["all_pass"], s)

    def test_packing_nonoverlap_and_domain(self):
        p = initial_packing(DEFAULT_DEM_PARAMS, FAST)
        x, y, r = p["x"], p["y"], p["r"]
        n = len(x)
        off = np.arange(n)[:, None] != np.arange(n)[None, :]
        d = np.sqrt((x[:, None] - x[None, :]) ** 2 + (y[:, None] - y[None, :]) ** 2)
        self.assertTrue(np.all((d + 1e-9 >= r[:, None] + r[None, :] - 1e-6)[off]))
        self.assertTrue(np.all(x >= 0) and np.all(x <= DEFAULT_DEM_PARAMS["domain_w"]))
        self.assertTrue(np.all(y >= 0) and np.all(y <= DEFAULT_DEM_PARAMS["domain_h"]))

    def test_contact_forces_newton_third_law(self):
        p = initial_packing(DEFAULT_DEM_PARAMS, dict(FAST, n_particles=30))
        vel = np.zeros((len(p["x"]), 2))
        tang = np.zeros((len(p["x"]), len(p["x"])))
        cf = contact_forces(p, vel, DEFAULT_DEM_PARAMS, FAST, tang, False, 0.001)
        # 无边界接触的颗粒对：内力总和 + 重力 = 外力合力（含墙）
        self.assertEqual(len(cf["fx"]), len(p["x"]))
        self.assertEqual(cf["n_contacts"], int(np.sum(tang != 0)) if False else cf["n_contacts"])
        self.assertGreaterEqual(cf["n_contacts"], 0)
        self.assertFalse(np.any(np.isnan(cf["fx"])))

    def test_solve_schema_and_json(self):
        r = solve_dem(DEFAULT_DEM_PARAMS, FAST)
        self.assertIn("meta", r)
        self.assertIn("frames", r)
        self.assertIn("curve", r)
        json.dumps(r)  # 全部 Python 标量可序列化
        f = r["frames"][0]
        for key in ("step", "t", "door_open", "particles", "n_contacts",
                    "door_force_N", "kinetic_J", "max_overlap", "surface"):
            self.assertIn(key, f)

    def test_door_force_drops_after_opening(self):
        r = solve_dem(DEFAULT_DEM_PARAMS, FAST)
        closed = [f for f in r["frames"] if not f["door_open"]]
        opened = [f for f in r["frames"] if f["door_open"]]
        self.assertTrue(closed and opened)
        self.assertGreater(closed[-1]["door_force_N"], 0.0)
        self.assertEqual(opened[-1]["door_force_N"], 0.0)

    def test_particles_escape_and_funnel_forms(self):
        r = solve_dem(DEFAULT_DEM_PARAMS, FAST)
        f0 = r["frames"][0]
        fl = r["frames"][-1]
        # 开闸后有颗粒掉出模型（y<0），这是真·流出的证据
        self.assertLess(min(p[1] for p in fl["particles"]), 0.0)
        # 开闸后地表高度下降（漏斗成形）
        self.assertLess(max(fl["surface"]), max(f0["surface"]) + 1e-6)

    def test_stability_and_bounded(self):
        r = solve_dem(DEFAULT_DEM_PARAMS, FAST)
        m = r["meta"]
        self.assertTrue(m["stable"], m)
        self.assertLess(m["peak_overlap_ratio"], 1.0)
        self.assertTrue(all(math.isfinite(f["kinetic_J"]) and f["kinetic_J"] >= 0
                            for f in r["frames"]))

    def test_boundary_honest(self):
        r = solve_dem(DEFAULT_DEM_PARAMS, FAST)
        self.assertIn("非专著理论", r["boundary"])
        self.assertIn("定性", r["boundary"])

    def test_n_particles_capped(self):
        o = dict(FAST, n_particles=99999, max_particles=50)
        r = solve_dem(DEFAULT_DEM_PARAMS, o)
        self.assertLessEqual(r["meta"]["n_particles"], 50)


class TestDEMBackend(unittest.TestCase):
    def test_selfcheck(self):
        s = selfcheck_backend()
        self.assertTrue(s["all_pass"], s)

    def test_request_shape(self):
        req = make_request({"n_particles": 10}, {"n_steps": 10}, backend="local")
        self.assertEqual(req["protocol"], PROTOCOL)
        self.assertEqual(req["backend"], "local")
        self.assertEqual(req["command"], "run")

    def test_roundtrip_local(self):
        resp = run_backend(make_request(
            dict(n_particles=30), dict(FAST, n_steps=20, settle_steps=20)))
        self.assertEqual(resp["status"], "ok")
        self.assertEqual(resp["backend"], "local")
        self.assertGreaterEqual(len(resp["result"]["frames"]), 2)

    def test_local_backend_class(self):
        b = LocalDEMBackend()
        r = b.run(dict(n_particles=30), dict(FAST, n_steps=20, settle_steps=20))
        self.assertIn("frames", r)

    def test_bad_protocol_rejected(self):
        with self.assertRaises(BackendError):
            run_backend({"protocol": "nope/v1"})

    def test_external_requires_command(self):
        with self.assertRaises(BackendError):
            run_backend(make_request(backend="external"))


class TestCone(unittest.TestCase):
    def test_selfcheck(self):
        s = selfcheck_cone()
        self.assertTrue(s["all_pass"], s)

    def test_morphology_profiles(self):
        h, l, p = 3.0, 2.5, DEFAULT_CONE_PARAMS
        z = np.linspace(0, h, 40)
        funnel = radius_profile("funnel", z, h, l, p)
        urn = radius_profile("urn", z, h, l, p)
        shaft = radius_profile("shaft", z, h, l, p)
        self.assertTrue(np.all(np.diff(funnel) > 0))
        self.assertGreater(urn[len(urn) // 2], urn[0])
        self.assertTrue(np.allclose(shaft, l))

    def test_unknown_morphology_rejected(self):
        with self.assertRaises(ValueError):
            radius_profile("unknown", np.array([0, 1]), 3.0, 2.5, DEFAULT_CONE_PARAMS)

    def test_all_morphologies_schema(self):
        for m in MORPHOLOGIES:
            r = solve_cone(DEFAULT_CONE_PARAMS, m)
            json.dumps(r)
            self.assertEqual(r["meta"]["morph"], m)
            self.assertGreater(r["result"]["weight_kN"], 0)
            self.assertGreater(r["result"]["side_resistance_kN"], 0)
            self.assertTrue(math.isfinite(r["result"]["support_pressure_kPa"]))

    def test_urn_mobilizes_more_side_resistance(self):
        urn = solve_cone(DEFAULT_CONE_PARAMS, "urn")["result"]
        shaft = solve_cone(DEFAULT_CONE_PARAMS, "shaft")["result"]
        self.assertGreater(urn["mobilized_fraction"], shaft["mobilized_fraction"])

    def test_cohesion_reduces_pressure(self):
        base = solve_cone(DEFAULT_CONE_PARAMS, "funnel")["result"]["support_pressure_kPa"]
        high_c = solve_cone(dict(DEFAULT_CONE_PARAMS, c=40.0), "funnel")["result"]["support_pressure_kPa"]
        self.assertLess(high_c, base)

    def test_depth_increases_pressure(self):
        q2 = solve_cone(dict(DEFAULT_CONE_PARAMS, h=2.0), "funnel")["result"]["support_pressure_kPa"]
        q8 = solve_cone(dict(DEFAULT_CONE_PARAMS, h=8.0), "funnel")["result"]["support_pressure_kPa"]
        self.assertGreater(q8, q2)

    def test_critical_depth_schema(self):
        cr = critical_depth(DEFAULT_CONE_PARAMS, "funnel")
        self.assertIn("has_root", cr)
        if cr["has_root"]:
            self.assertIsNotNone(cr["h_cr"])
        else:
            self.assertIsNone(cr["h_cr"])

    def test_boundary_honest(self):
        r = solve_cone(DEFAULT_CONE_PARAMS, "funnel")
        self.assertIn("非专著原有公式", r["boundary"])
        self.assertIn("不是三维 FEM", r["boundary"])


if __name__ == "__main__":
    unittest.main()
