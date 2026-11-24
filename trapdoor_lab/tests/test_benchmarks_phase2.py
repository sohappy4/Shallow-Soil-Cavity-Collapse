# -*- coding: utf-8 -*-
"""二期 1g / 离心机文献锚点回归测试。

只核验：数据齐全、数值有限、关键换算自洽、忠于专著原文数值。
不声称从 L1/L2/L3 推导出文献试验值，不拟合试验曲线，不虚构逐点数据。
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from benchmarks import (
    LITERATURE_ANCHORS, BENCHMARKS, check_literature_anchors,
)


def _by_key(group, key):
    for a in LITERATURE_ANCHORS.get(group, []):
        if a.get("key") == key:
            return a
    return {}


class TestPhase2Anchors(unittest.TestCase):
    def test_all_checks_pass(self):
        r = check_literature_anchors()
        self.assertTrue(r["all_pass"], "核验失败: %s" % r["errors"])
        self.assertTrue(all(r["checks"].values()))
        self.assertEqual(r["errors"], [])

    def test_group_counts(self):
        r = check_literature_anchors()
        self.assertGreaterEqual(r["counts"]["one_g"], 10)
        self.assertGreaterEqual(r["counts"]["centrifuge"], 4)
        self.assertGreaterEqual(r["counts"]["rate_applicability"], 5)

    def test_every_anchor_has_provenance(self):
        for group, rows in LITERATURE_ANCHORS.items():
            for a in rows:
                with self.subTest(group=group, key=a.get("key")):
                    self.assertTrue(a.get("label"))
                    self.assertTrue(a.get("unit"))
                    self.assertTrue(a.get("source"))   # 原始文献
                    self.assertTrue(a.get("ref"))      # 专著章节
                    self.assertIn(a.get("kind"),
                                  ("value", "range", "threshold"))

    def test_mcnulty_values_match_book(self):
        self.assertEqual(_by_key("one_g", "mcnulty_displacement_sensitivity")["value"], 0.0002)
        self.assertEqual(_by_key("one_g", "mcnulty_h_over_D_threshold")["value"], 2.0)
        self.assertEqual(_by_key("one_g", "mcnulty_equal_settlement_axisym")["value_min"], 1.0)
        self.assertEqual(_by_key("one_g", "mcnulty_equal_settlement_axisym")["value_max"], 1.5)
        self.assertEqual(_by_key("one_g", "mcnulty_equal_settlement_plane")["value_min"], 1.5)
        self.assertEqual(_by_key("one_g", "mcnulty_equal_settlement_plane")["value_max"], 2.5)
        self.assertEqual(_by_key("one_g", "mcnulty_displacement_to_steady_axisym")["value"], 3.0)
        self.assertEqual(_by_key("one_g", "mcnulty_displacement_to_steady_plane")["value"], 10.0)

    def test_outcrop_and_dilatancy_values(self):
        self.assertEqual(_by_key("one_g", "abdulla_outcrop_angle")["value"], 61.0)
        self.assertEqual(_by_key("one_g", "han_dilatancy")["value"], 8.0)
        self.assertEqual(_by_key("one_g", "liang_dilatancy")["value"], 7.5)
        self.assertEqual(_by_key("one_g", "li_cover_ratio")["value"], 2.0)

    def test_centrifuge_beta_eff_and_plateau(self):
        self.assertEqual(_by_key("centrifuge", "dewoolkar_beta_eff_min")["value"], 4.17)
        self.assertEqual(_by_key("centrifuge", "dewoolkar_beta_eff_max")["value"], 5.26)
        self.assertEqual(_by_key("centrifuge", "dewoolkar_scatter")["value"], 11.0)
        self.assertEqual(_by_key("centrifuge", "dewoolkar_plateau_h_over_D")["value"], 5.0)

    def test_beta_span_matches_scatter(self):
        r = check_literature_anchors()
        # 4.17–5.26 的半宽/中点为 11.56%，与标注 ±11% 量级一致
        self.assertAlmostEqual(r["beta_half_span_pct"], 11.558854718981971, places=6)
        self.assertLess(abs(r["beta_half_span_pct"] - 11.0), 2.0)

    def test_rate_applicability_values(self):
        self.assertEqual(_by_key("rate_applicability", "l3_applicability")["value"], 1.4)
        self.assertEqual(_by_key("rate_applicability", "l3_dev_at_hd2")["value"], 14.0)
        self.assertEqual(_by_key("rate_applicability", "l3_dev_at_hd10")["value"], 174.0)
        self.assertEqual(_by_key("rate_applicability", "self_stable_tol")["value"], 2.0)
        self.assertEqual(_by_key("rate_applicability", "N_sat_consistency")["value"], 18.2212)

    def test_Nsat_over_4_matches_hsat_printed(self):
        r = check_literature_anchors()
        self.assertTrue(r["checks"]["Nsat_over_4_eq_hsat"])
        printed = BENCHMARKS["phase2_literature"]["N_sat_FELA"]["h_sat_over_D"]
        self.assertLess(abs(18.2212 / 4.0 - printed) / printed, 0.01)

    def test_l3_domain_matches_saturation_anchor(self):
        # saturation_anchor 内部以 1.4 作为 h/D 门槛，锚点必须一致，防止两处漂移
        self.assertEqual(_by_key("rate_applicability", "l3_applicability")["value"], 1.4)
        r = check_literature_anchors()
        self.assertTrue(r["checks"]["l3_domain_consistent"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
