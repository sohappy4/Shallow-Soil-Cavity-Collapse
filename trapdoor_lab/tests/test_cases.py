# -*- coding: utf-8 -*-
"""方案 M9：五个工程案例模板与复算引擎回归测试。

覆盖：模板数量/字段、案例一 A/B 变体、原书数值逐项对照、
L_c(h_cr)=l 自洽反验、被动顶起无 h_cr、渐进闭合支不可容许、
以及 /api/cases 与 /api/case 的端到端返回结构。
"""
import json
import math
import threading
import unittest
import urllib.request
from http.server import ThreadingHTTPServer

from app import Handler
from cases import (ENGINEERING_CASES, CASE_TEMPLATES, evaluate_case,
                   evaluate_all, case_summary, selfcheck_cases,
                   _active_calc, _passive_calc, _progressive_calc)

REQUIRED = ("id", "name", "source", "scenario", "mode", "params",
            "expected", "morphology", "judgement")


class TestCaseTemplates(unittest.TestCase):
    def test_count_and_unique_ids(self):
        self.assertEqual(len(ENGINEERING_CASES), 5)
        self.assertEqual(len(CASE_TEMPLATES), 5)
        self.assertEqual(len({c["id"] for c in ENGINEERING_CASES}), 5)

    def test_required_fields(self):
        for c in ENGINEERING_CASES:
            for k in REQUIRED:
                self.assertIn(k, c, "%s 缺字段 %s" % (c["id"], k))
            for k in ("sigma_c", "A", "B", "sigma_t", "gamma", "q2", "l", "h"):
                self.assertIn(k, c["params"], "%s 参数缺 %s" % (c["id"], k))
            self.assertTrue(c["source"].startswith("专著"))
            self.assertTrue(c["expected"])

    def test_ids_exact(self):
        self.assertEqual(
            [c["id"] for c in ENGINEERING_CASES],
            ["karst-cavity", "urban-road", "goaf",
             "shield-tunnel", "tunnel-uplift"])

    def test_karst_two_variants(self):
        c = CASE_TEMPLATES["karst-cavity"]
        self.assertEqual(len(c["variants"]), 2)
        self.assertEqual([v["mode"] for v in c["variants"]],
                         ["active", "progressive"])
        self.assertIn("瓮形", c["variants"][0]["label"])
        self.assertIn("井筒", c["variants"][1]["label"])

    def test_selfcheck_all_pass(self):
        s = selfcheck_cases()
        self.assertTrue(s["all_pass"], s)
        for k, v in s.items():
            if k != "all_pass":
                self.assertTrue(v, "%s 未通过" % k)


class TestCaseRecompute(unittest.TestCase):
    def test_all_book_comparison_ok(self):
        for r in evaluate_all():
            self.assertTrue(r["comparison_ok"], r["case"]["id"])
            self.assertTrue(r["all_comparison_ok"], r["case"]["id"])
        self.assertTrue(case_summary()["all_comparison_ok"])

    def test_karst_active_A_numbers(self):
        r = evaluate_case(CASE_TEMPLATES["karst-cavity"])
        res = r["result"]
        book = r["comparison"]
        for k in ("Omega", "A_psi", "S_psi", "p", "L_c", "x2", "h_cr",
                  "surface_width"):
            self.assertIn(k, book)
            self.assertLess(book[k]["rel_error"], 0.011, k)
        self.assertTrue(res["global_collapse"])   # h=6 < h_cr=10.33
        self.assertTrue(res["branch_valid"])

    def test_karst_progressive_B_closure_inadmissible(self):
        r = evaluate_case(CASE_TEMPLATES["karst-cavity"])
        v = r["variants"][1]
        self.assertTrue(v["comparison_ok"], v["comparison"])
        self.assertTrue(v["result"]["closure_inadmissible"])
        self.assertGreater(v["result"]["closure_exceed_kPa"], 7.0)
        self.assertFalse(v["result"]["branch_valid"])
        self.assertEqual(v["result"]["surface_width"], 4.0)

    def test_selfcheck_Lc_at_hcr(self):
        for r in evaluate_all():
            e = r["result"].get("selfcheck_error")
            if e is not None:
                self.assertLess(e, 1e-6, r["case"]["id"])
        # 独立取一组参数再验一次（不是只依赖模板缓存）
        p = dict(CASE_TEMPLATES["urban-road"]["params"])
        a = _active_calc(p)
        self.assertAlmostEqual(
            _active_calc(dict(p, h=a["h_cr"]))["L_c"], p["l"], places=6)
        self.assertTrue(case_summary()["selfcheck_ok"])

    def test_passive_uplift_no_hcr(self):
        r = evaluate_case(CASE_TEMPLATES["tunnel-uplift"])
        res = r["result"]
        self.assertIsNone(res["h_cr"])
        self.assertTrue(res["global_collapse"])
        self.assertTrue(res["branch_valid"])
        for k in ("S", "p", "Delta", "x2", "q1"):
            self.assertLess(r["comparison"][k]["rel_error"], 0.011, k)
        self.assertGreater(res["q1"], 1000.0)  # ≈1.11 MPa

    def test_nonassociated_reduction(self):
        p = dict(CASE_TEMPLATES["karst-cavity"]["params"])
        a = _active_calc(p)
        self.assertLess(a["A_psi"], p["A"])
        self.assertLess(a["S_psi"], p["sigma_c"] * p["A"] ** 2 * p["B"])
        self.assertLess(a["nonassociated_reduction"], 1.0)

    def test_active_vs_passive_morphology_width(self):
        ratios = {x["id"]: x["width_over_span"]
                  for x in case_summary()["width_over_span"]
                  if "/" not in x["id"]}
        for cid in ("karst-cavity", "urban-road", "goaf", "shield-tunnel"):
            self.assertLess(ratios[cid], 1.0, cid)  # 主动瓮形：口小于洞跨
        self.assertGreater(ratios["tunnel-uplift"], 3.0)  # 顶起外倾

    def test_width_range_by_mode(self):
        s = case_summary()
        self.assertLess(s["active_range"][0], 0.3)
        self.assertLess(s["active_range"][1], 0.6)   # 只含主动瓮形，不含渐进近井筒
        self.assertEqual(s["progressive_range"], [1.0])
        self.assertAlmostEqual(s["uplift_range"][0], 3.1259, places=3)
        self.assertIn("近井筒 1.00", s["key_warning"])
        self.assertTrue(all("mode" in x for x in s["width_over_span"]))

    def test_json_serializable(self):
        for r in evaluate_all():
            json.dumps(r, ensure_ascii=False)
        json.dumps(case_summary(), ensure_ascii=False)
        json.dumps(selfcheck_cases(), ensure_ascii=False)

    def test_pure_calc_helpers_direct(self):
        # 三个分支都能独立调用，不抛异常且字段自洽
        a = _active_calc(dict(CASE_TEMPLATES["urban-road"]["params"]))
        b = _passive_calc(dict(CASE_TEMPLATES["tunnel-uplift"]["params"]))
        c = _progressive_calc(
            dict(CASE_TEMPLATES["karst-cavity"]["variants"][1]["params"]))
        self.assertGreater(a["h_cr"], 0)
        self.assertGreater(b["q1"], 0)
        self.assertTrue(math.isclose(c["x2"], c["progressive_lower"]))


class TestCaseAPI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        cls.port = cls.httpd.server_port
        cls.thread = threading.Thread(target=cls.httpd.serve_forever,
                                      daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()
        cls.thread.join(timeout=5)

    def _get(self, path):
        with urllib.request.urlopen(
                "http://127.0.0.1:%d%s" % (self.port, path), timeout=30) as r:
            return json.loads(r.read().decode("utf-8"))

    def _post(self, path, obj):
        req = urllib.request.Request(
            "http://127.0.0.1:%d%s" % (self.port, path),
            data=json.dumps(obj).encode("utf-8"),
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read().decode("utf-8"))

    def test_cases_index(self):
        d = self._get("/api/cases")
        self.assertEqual(len(d["cases"]), 5)
        for c in d["cases"]:
            for k in ("id", "name", "source", "scenario", "mode",
                      "morphology", "judgement", "caution"):
                self.assertIn(k, c)
        self.assertTrue(d["summary"]["all_comparison_ok"])
        self.assertTrue(d["selfcheck"]["all_pass"])

    def test_case_endpoint(self):
        for cid in CASE_TEMPLATES:
            d = self._post("/api/case", {"id": cid})
            self.assertNotIn("error", d, cid)
            z = d["result"]
            self.assertEqual(z["case"]["id"], cid)
            self.assertTrue(z["all_comparison_ok"], cid)
            self.assertTrue(d["fig"].startswith("iVBOR"), cid)  # PNG base64

    def test_unknown_id(self):
        d = self._post("/api/case", {"id": "not-a-case"})
        self.assertIn("error", d)


if __name__ == "__main__":
    unittest.main(verbosity=2)
