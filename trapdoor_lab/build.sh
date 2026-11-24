#!/bin/bash
# 可复现打包脚本：PyInstaller --onedir
# 用法：bash build.sh
# 说明：本机 /tmp 是 10MB tmpfs，必须用项目内临时目录，否则 pip/PyInstaller 报 Errno 28。
set -e
cd "$(dirname "$0")" || exit 1

PY="../venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "找不到 ../venv/bin/python，请先创建虚拟环境并安装依赖。"
  exit 1
fi
if ! "$PY" -c "import PyInstaller" 2>/dev/null; then
  echo "缺少 PyInstaller，请先执行："
  echo "  TMPDIR=../.pip_tmp $PY -m pip install --no-cache-dir pyinstaller"
  exit 1
fi
if [ -d dist/trapdoor-lab ]; then
  echo "dist/trapdoor-lab 已存在。批量删除会被环境安全策略拦截，请手动清理后再运行："
  echo "  rm -rf dist/trapdoor-lab"
  exit 1
fi

mkdir -p ../.pip_tmp
TMPDIR="$(cd .. && pwd)/.pip_tmp" "$PY" -m PyInstaller --noconfirm --onedir \
  --name trapdoor-lab --log-level=WARN app.py

# ---- M9：随包分发授权/文档/案例，并生成发布包内清单 ----
DOCS="dist/trapdoor-lab/文档"
mkdir -p "$DOCS"
for f in \
  "../LICENSE" "../THIRD_PARTY_NOTICES.md" "../授权说明.md" "../RELEASE.md" \
  "../README.md" "../MANIFEST.sha256" \
  "运行说明.md" "API说明.md" "参数手册.md" "验收报告.md" \
  "验证基准说明.md" "图号映射表.md" "cases.py" "tests/test_cases.py"; do
  [ -f "$f" ] && cp "$f" "$DOCS/" 2>/dev/null || true
done

# 发布包清单（源码文件 SHA-256 + 打包目录文件数）
if [ -f "$DOCS/MANIFEST.sha256" ]; then
  echo "源码校验清单已随包分发：文档/MANIFEST.sha256"
fi

echo
echo "打包完成：dist/trapdoor-lab/"
echo "启动脚本：dist/trapdoor-lab/启动打包版.sh"
echo "随附文档：dist/trapdoor-lab/文档/（授权、许可证声明、发布说明、五份技术文档）"
echo "注意：改了源码后需重新运行本脚本，打包版才会更新。"
