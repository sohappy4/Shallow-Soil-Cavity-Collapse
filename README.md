# Shallow-Soil-Cavity-Collapse

浅埋土洞活门试验数值模拟程序。理论依据：专著《浅埋土洞坍塌：形态、临界埋深与统一变分理论》及论文
*"A unified variational limit-analysis theory for active collapse and passive uplift of shallow soil cavities"*（孙海朋）。

程序把「静力层—极限分析层—速率过程层」三层理论产品化：
主动塌落（瓮形）与被动顶起（漏斗形）互为镜像，剪胀角控制形态，临界埋深 `h_cr` 控制是否贯通至地表。

## 快速启动

银河麒麟 V10 SP1 / aarch64 已验证。

**方式一（推荐，源码版）**：双击桌面上的 **活门试验模拟程序** 图标；
或在本目录执行：

```bash
cd trapdoor_lab
bash 一键启动.sh
```

浏览器自动打开 `http://127.0.0.1:8765`。若该端口被占用，程序自动顺延到 8766、8767…
并以终端打印的地址为准。关闭时按 `Ctrl+C`。

**方式二（独立打包版，不需 Python）**：
双击 `trapdoor_lab/dist/trapdoor-lab/启动打包版.sh`。约 157 MB，可整体拷到同为 aarch64 的机器。
脚本会自动探测实际端口并打印浏览器地址。

## 功能

| 层 | 界面按钮 | 内容 |
| --- | --- | --- |
| L1 | 开始计算 | 准静态闭式解：S、t、Ω、Λ、ψ₀、ψ(h)、h_cr、q1a/q1b；六张图；Table 3/4 复现 |
| L2 | L2 打靶求滑动面 | 专著算法 17.1：RK4 + 二分打靶，三级筛选，与闭式解对照 |
| L3 | L3 速率时间线 | 式 (18.1)–(18.7)：Beverloo 出流、v、T、R、μ(I)、五阶段；带时间线播放器 |
| GTT | GTT 切点 / 表13.2 / 定理13.2 | 广义切线法：切线参数、包线夹逼、切点扫描；原表舍入差显式标注 |
| 锚点 | 文献锚点 h_sat / h_cr | N_sat=18.2212、h_sat/D=4.555 自洽核验与门槛区分 |
| 对照 | 主动/被动镜像对照 | 主动滑面关于 x=l 镜像后与被动逐点比对 |
| 动画 | 二维坍塌过程动画 | 滑体土粒逐帧下落、活门下移、地表漏斗成形（运动学演示，非 DEM/应力云图）；样本 `trapdoor_lab/verify/anim_preview.png` |
| 动态 | 参数扫描动态图件 | h / A / B / σ_t 逐帧重算，直接观察 h 跨 h_cr 时结论翻转 |
| 报告 | 导出计算报告 (HTML) | 输入、结果、六图、Table 3/4、可复现性元数据一并导出 |

每次 `/api/solve` 返回参数快照、git commit、代码 sha256、依赖版本、时间戳，随报告导出。

## 验证

```bash
cd trapdoor_lab
../venv/bin/python -m unittest discover -s tests   # 68 项回归
../venv/bin/python 自检.py                          # 13 项基准自检
```

`.github/workflows/test.yml` 在每次 push/PR 上跑全量单测、自检与编译检查。
冻结基准在 `tests/golden_baseline.json`，由 `verify/make_golden.py` 实跑生成，测试只读不改。

## 文档

- `trapdoor_lab/运行说明.md` — 小白操作步骤
- `trapdoor_lab/参数手册.md` — 参数符号、单位、范围、预设
- `trapdoor_lab/API说明.md` — 六个接口与返回字段
- `trapdoor_lab/验证基准说明.md` — 基准分级、容差、理论张力
- `trapdoor_lab/图号映射表.md` — 代码图 / 论文图 / 专著图三级对照
- `trapdoor_lab/验收报告.md` — 实跑结果与边界记录

## 环境

依赖：Python 3.8+，numpy、scipy、matplotlib（见 `requirements.txt`）。
本机 aarch64 无 PySide6/PyQt5 预编译包，因此界面采用 Python 标准库 `http.server` + 浏览器，
零 GUI 框架依赖。重新打包见 `trapdoor_lab/build.sh`。

## 理论边界（必读）

- `n = h/(2l) = h/D`，不是 `h/l`；`q1a`/`q1b` 是输出，不是输入；B=1 走线性极限分支。
- `h_sat`（力学饱和，拱效应层）与 `h_cr`（几何贯通，上限法）是两个门槛，不可混用。
- L3 演示锚点对应 h/D=6，但专著写适用域 h/D≤1.4；程序按公式如实算并标记 `rate-degraded`。
- GTT 只给切线层与固定机构上的**耗散积分夹逼**；专著未给 Q^UB(σ₀) 机构荷载表达式，
  程序不虚构严格荷载上界，也不冒充 FELA 上下界。
- L1 闭式解只有沿滑面的一维场，「参数扫描动态图件」不是二维应力云图；
  二维云图需 FELA/有限元层，未实现。
- 「二维坍塌过程动画」是运动学演示：土粒位移是满足滑面约束与体积守恒的构造模板，
  不是 DEM 颗粒接触求解，也不输出应力/塑性区数值。真·力学场需 FELA/FEM（M7）
  或 Yade/PFC 独立后端（M8）；后者在 aarch64 银河麒麟安装极难。
- 未实现：完整 FELA SOCP、DEM 后端、三维锥形求解器、局部坍塌形态数据。

## 安全

对话中曾明文出现的 GitHub 密码与令牌均已暴露，应到 GitHub 网页撤销并改用 SSH key。
不要把凭据写入仓库、配置文件或对话。

## 引用

如本项目对你的研究有帮助，请引用孙海朋的相关论文与专著。
