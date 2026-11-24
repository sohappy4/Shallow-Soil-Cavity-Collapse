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
| FEM | FEM 二维力学场（应力/塑性区） | 平面应变 Q4 + Drucker–Prager 径向返回 + 位移控制活门；输出偏应力云图、塑性区、压力—位移曲线；样本 `trapdoor_lab/verify/fem_preview.png`；定性—半定量扩展，非 FELA 上下界 |
| FEM | FEM 应力场动态云图 | 逐位移步应力/塑性区演化，全场统一色标，可滑块逐帧查看 + GIF 自动播放（方案 M5/M6）；样本 `trapdoor_lab/verify/fem_demo.gif` |
| FELA | FELA 极限分析上下界夹逼 | 线性化屈服面 + 对偶线性规划：内接多面体给静力下限、外接给运动上限；主动支护压力夹逼并包住 L1 闭式解，Davis 折减独立报告（方案 M7）；样本 `trapdoor_lab/verify/fela_bracket.png` |
| DEM | DEM 二维颗粒流（真接触求解） | Cundall–Strack 线性接触 + Coulomb 摩擦，逐颗粒牛顿积分；活门开闭、颗粒流出、地表漏斗，四阶段速度着色图（方案 M8）；样本 `trapdoor_lab/verify/dem_preview.png` |
| 3D | 三维锥形极限平衡 | 轴对称漏斗/瓮形/井筒三形态滑体，自重与侧阻积分、临界埋深；与 L1 二维解只作量级对照（方案 M8 配套） |
| 动态 | 参数扫描动态图件 | h / A / B / σ_t 逐帧重算，直接观察 h 跨 h_cr 时结论翻转 |
| 报告 | 导出计算报告 (HTML) | 输入、结果、六图、Table 3/4、可复现性元数据一并导出 |

每次 `/api/solve` 返回参数快照、git commit、代码 sha256、依赖版本、时间戳，随报告导出。

## 验证

```bash
cd trapdoor_lab
../venv/bin/python -m unittest discover -s tests   # 128 项回归（含 23 FEM + 13 FELA + 24 DEM/3D）
../venv/bin/python 自检.py                          # 13 项基准自检
```

`.github/workflows/test.yml` 在每次 push/PR 上跑全量单测、自检与编译检查。
冻结基准在 `tests/golden_baseline.json`，由 `verify/make_golden.py` 实跑生成，测试只读不改。

## 文档

- `trapdoor_lab/运行说明.md` — 小白操作步骤
- `trapdoor_lab/参数手册.md` — 参数符号、单位、范围、预设
- `trapdoor_lab/API说明.md` — 全部接口与返回字段（含 `/api/fem`）
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
- L1 闭式解只有沿滑面的一维场；「参数扫描动态图件」不是二维应力云图。
  真·二维力学场由新增 FEM 模块提供（见下条），但两者理论族不同，不做同一性对照。
- 「二维坍塌过程动画」是运动学演示：土粒位移是满足滑面约束与体积守恒的构造模板，
  不是 DEM 颗粒接触求解，也不输出应力/塑性区数值。
- 「FEM 二维力学场」是横向引用的成熟有限元扩展（Q4 + DP 径向返回 + 强度折减框架），
  按用户 2026-11-24 授权引入：应力场/塑性区/压力—位移曲线为数值解，
  但**不是严格 FELA 上下界**，DP 圆化准则与专著幂型 Hoek–Brown 不同族，
  小变形理想弹塑性、网格较粗，仅作定性—半定量演示；L1/L2/L3 仍为专著理论主结果。
- 「FELA 极限分析上下界夹逼」用线性化屈服面的对偶线性规划给出严格数学夹逼
  （内接→静力下限、外接→运动上限，q 夹逼包住 L1 闭式解），但仍**不是 SOCP 圆锥规划**、
  不是幂型精确解；下限与 ψ 无关、Davis 折减单独报告；深埋 h/D>1.4 因幂型包线无
  临界状态封顶，显式标记非物理、不给夹逼结论；N_sat=18.2212 仅为文献锚点，本模块不推导。
- 「DEM 二维颗粒流」是**真·颗粒接触求解**（非运动学演示）：Cundall–Strack 线性接触、
  切向历史与 Coulomb 摩擦、逐颗粒牛顿积分。但为二维圆盘、无转动自由度、刚度/摩擦为
  数值参数（非 μ(I) 细观反算），横向引用成熟离散元方法、非专著理论，定性—半定量；
  Yade/PFC 等外部后端通过 `dem_backend.py` 的独立进程 JSON 协议接入（GPL/商业代码
  不链接、不分发），当前程序包只含自实现的 MIT 参考后端。
- 「三维锥形极限平衡」是轴对称刚体滑体 + 侧摩阻积分，三形态命名对齐专著，但轴对称
  推广为本程序扩展、非专著原有公式；不是三维 FEM/FELA/SOCP，与 L1 二维解只作量级对照。
- 未实现：完整 FELA SOCP 圆锥规划、Yade/PFC 外部 DEM 求解器本身（仅定义接入协议）、
  三维 FEM/FELA、局部坍塌形态数据。

## 安全

对话中曾明文出现的 GitHub 密码与令牌均已暴露，应到 GitHub 网页撤销并改用 SSH key。
不要把凭据写入仓库、配置文件或对话。

## 引用

如本项目对你的研究有帮助，请引用孙海朋的相关论文与专著。
