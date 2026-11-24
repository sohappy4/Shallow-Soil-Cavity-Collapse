# Shallow-Soil-Cavity-Collapse

浅埋土洞坍塌机理的统一变分极限分析理论配套计算程序。
基于论文 *"A unified variational limit-analysis theory for active collapse and
passive uplift of shallow soil cavities"* (Haipeng Sun)。

在一个无需预设滑动面的统一变分框架内，闭式求解浅埋土洞的两种破坏模式
——主动坍塌 (active collapse) 与被动隆起 (passive uplift)——
给出滑动面形态、滑动面应力场、临界覆盖厚度、极限荷载、
剪胀角随深度的演化、侧压力系数演化，以及非关联剪胀修正，
并自动复现论文 Table 3 / Table 4 的试验对比。

## 功能一键概览

运行脚本后：

* 终端打印核心闭式解：强度组合参数 S、出露坡度参数 t、脆性比 Ω、
  覆强比 Λ、地表/底部剪胀角 ψ0 与 ψ(h)、临界覆盖厚度 h_cr、
  主动支护压力 q1a、被动顶起压力 q1b，以及镜面对称校验
  (q1b − q1a = 2·R_a/l)。
* 自动生成 6 张图：

| 文件 | 内容 | 对应论文 |
|------|------|----------|
| fig1_slip_surface.png | 滑动面形态（主动 urn / 被动 funnel） | Fig.1, Eq.25 |
| fig2_stress_field.png | 滑动面应力场 σ、τ 随深度 | Fig.4(a), Eq.28/29 |
| fig3_dilatancy.png    | 剪胀角 ψ 随深度衰减 | Fig.4(b), Eq.30 |
| fig4_zoning.png       | 四形态分区图（Ω–ψ0 平面） | Fig.5 / Table 2, Eq.38 |
| fig5_Ka_depth.png     | 侧压力系数 Ka(z) 随深度演化 | Fig.6, Eq.44 |
| fig6_nonassociated.png| 非关联修正 Aψ/A–ψ0 | Fig.7(a), Eq.46 |

* 终端自动复现 Table 3（主动试验 Chen et al., 2019）与
  Table 4（被动试验 Rui et al., 2025b），打印预测值与偏差。

## 环境配置

```bash
# 创建并激活虚拟环境
python3 -m venv venv
source venv/bin/activate

# 安装依赖
pip install -r requirements.txt
```

> 若创建虚拟环境时报 `ensurepip is not available`，先执行
> `sudo apt install python3-venv python3-pip -y` 再重试。

## 运行

```bash
python3 cavity_collapse_calculator.py
```

脚本所在目录会生成 fig1~fig6 共 6 张 png。

## 换工况 / 改参数

只改脚本顶部 `PARAMS` 这一块即可，其余一行都不用动（单位已注释）：

```python
PARAMS = dict(
    sigma_c = 100.0,   # 强度参考应力 (kPa)
    A       = 0.6,     # 幂型强度准则系数
    B       = 0.6,     # 幂型强度准则指数 (0<B<=1, B=1 退化为 Mohr-Coulomb)
    sigma_t = 5.0,     # 抗拉强度 (kPa)
    gamma   = 18.0,    # 土体重度 (kN/m^3)
    q2      = 10.0,    # 地表均布超载 (kPa)
    l       = 2.5,     # 洞体半宽 (m)
    h       = 3.0,     # 覆盖层厚度 (m)
)
```

默认参数即论文 Fig.4 基准算例，可用于自检：地表剪胀角应约 37.6°、
底部约 22.7°，镜面对称校验误差为 0。

## 适用范围

平面应变、准静态、沿单一连续剪切带破坏的剪胀或等容材料；
主动模式额外要求 h ≤ h_cr。所有极限解为运动学上限，
关联流动下为严格上界，非关联修正后为近似值。

## 引用

如本项目对你的研究有帮助，请引用原论文。
