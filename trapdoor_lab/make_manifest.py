# -*- coding: utf-8 -*-
"""生成发布清单 MANIFEST.sha256（可复现校验）。

覆盖：全部 Python 模块、测试、文档、配置、模板与启动脚本；
排除：venv / __pycache__ / build* / dist* / *.pyc / 临时文件。
用法：../venv/bin/python make_manifest.py
输出：项目根目录 MANIFEST.sha256（供发布包校验）
"""
from __future__ import annotations

import hashlib
import os
from typing import List

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXCLUDE_DIRS = {"venv", "__pycache__", "build", "build_new", "build_new2",
                "dist", "dist_new", "dist_new2", ".git", ".workbuddy",
                ".pip_tmp", "__pypackages__"}
EXCLUDE_EXT = {".pyc", ".pyo", ".tmp", ".swp", ".log", ".gz", ".zip", ".rar"}
EXCLUDE_FILES = {"MANIFEST.sha256", "_inline_check.js", "_askpass_tmp.sh"}


def _iter_files():
    for base, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS and
                   not d.startswith(".")]
        for name in sorted(files):
            if name in EXCLUDE_FILES or os.path.splitext(name)[1] in EXCLUDE_EXT:
                continue
            full = os.path.join(base, name)
            rel = os.path.relpath(full, ROOT).replace(os.sep, "/")
            if rel.startswith("trapdoor_lab/verify/_"):
                continue
            yield rel, full


def main() -> None:
    rows: List[str] = []
    for rel, full in _iter_files():
        h = hashlib.sha256()
        with open(full, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        rows.append("%s  %s" % (h.hexdigest(), rel))
    rows.sort(key=lambda s: s.split("  ", 1)[1])
    out = os.path.join(ROOT, "MANIFEST.sha256")
    with open(out, "w", encoding="utf-8") as f:
        f.write("# Shallow-Soil-Cavity-Collapse release manifest\n")
        f.write("# sha256  relative_path\n")
        f.write("\n".join(rows) + "\n")
    print("wrote %s (%d files)" % (out, len(rows)))


if __name__ == "__main__":
    main()
