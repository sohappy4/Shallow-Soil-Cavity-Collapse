#!/bin/bash
# 浅埋土洞活门试验数值模拟程序 —— 一键启动（银河麒麟 / Linux）
# 用法：在本文件上点右键「运行」，或在终端执行 bash 一键启动.sh
cd "$(dirname "$0")" || exit 1

PY="../venv/bin/python"
if [ ! -x "$PY" ]; then
  echo "找不到项目虚拟环境，改用系统 python3"
  PY="python3"
fi

if ! "$PY" -c "import numpy, matplotlib, scipy" 2>/dev/null; then
  echo "缺少依赖 numpy/matplotlib/scipy。"
  echo "请在终端进入项目目录后执行： pip install -r requirements.txt"
  echo
  read -r -p "按回车关闭窗口..." unused
  exit 1
fi

echo "正在启动，浏览器将自动打开（若 8765 被占用，程序会自动改用 8766 等端口，见下方实际地址）"
echo
"$PY" app.py
code=$?
echo
if [ "$code" -ne 0 ]; then
  echo "程序异常退出（代码 $code）。请把上面的提示信息拍照留存。"
  read -r -p "按回车关闭窗口..." unused
fi
exit "$code"
