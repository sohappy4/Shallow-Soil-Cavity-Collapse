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

echo
echo "打包完成：dist/trapdoor-lab/"
echo "启动脚本：dist/trapdoor-lab/启动打包版.sh"
echo "注意：改了源码后需重新运行本脚本，打包版才会更新。"
