# 第三方组件与许可证声明

本声明覆盖 `Shallow-Soil-Cavity-Collapse` 仓库及其发布包。主程序代码采用仓库根目录
`LICENSE` 所载 **MIT License**（Copyright © 2026 sohappy4）。本文件只声明**外部依赖**
的许可证与边界，不改动 MIT 主许可。

依赖版本与许可证自本仓库 `venv` 的 `importlib.metadata` 元数据实测读取（2026-11-24），
不是推测。版本锁定见 `requirements.txt`。

## 一、运行时依赖（随程序使用）

| 组件 | 锁定版本 | 许可证 | 用途 | 分发方式 |
|---|---|---|---|---|
| numpy | 1.24.4 | BSD-3-Clause | 数组、线性代数、DEM/FEM/FELA 内核 | 由用户环境/打包工具提供 |
| scipy | 1.10.1 | BSD-3-Clause（含 Enthought 等版权声明） | `linprog`/HiGHS、稀疏矩阵求解 | 同上 |
| matplotlib | 3.7.5 | PSF License（含BSD类组件） | 六图、应力云图、动态图、报告 | 同上 |
| six | 1.17.0 | MIT | matplotlib 间接依赖 | 间接 |
| python-dateutil | 3.0.0 | BSD / Apache 双许可 | matplotlib 间接依赖 | 间接 |
| kiwisolver | 1.4.7 | BSD-3-Clause | matplotlib 布局 | 间接 |
| pyparsing | 3.1.4 | MIT | matplotlib 间接依赖 | 间接 |
| contourpy | 1.1.1 | BSD-3-Clause | matplotlib 等高线（云图） | 间接 |
| cycler | 0.12.1 | BSD 类 | matplotlib 样式循环 | 间接 |
| fonttools | 4.57.0 | MIT | matplotlib 字体 | 间接 |
| pillow | 10.4.0 | HPND | 图像编码（PNG/GIF） | 间接 |

上述均为宽松许可（BSD/MIT/PSF/HPND/Apache），允许随 MIT 主程序分发，
保留各自版权声明即可。完整许可文本以各发行版随附的 `LICENSE`/`METADATA` 为准。

## 二、构建期依赖（不进入运行时）

| 组件 | 版本 | 许可证 | 说明 |
|---|---|---|---|
| PyInstaller | 6.22.3 | GPLv2-or-later **with special exception** | 仅用于生成打包产物；其 bootloader 特殊例外允许打包非 GPL 应用而不传染主程序代码 |

打包产物 `dist/trapdoor-lab/` 含 PyInstaller bootloader（依其特殊例外分发），
**不把 PyInstaller 源码链接进本程序代码**。

## 三、DEM 外部后端（许可证隔离，方案硬约束）

实施方案 §五·许可证边界规定：

> 仓库现有 MIT LICENSE 保留；PySide6 按 LGPL 动态链接使用；DEM 后端
> （Yade GPL / PFC 商业）以独立进程 + 文件/标准输入输出通信，不链接、
> 不随包分发，规避传染。

本程序 **不包含** Yade、PFC 或任何 GPL/商业求解代码。

- `dem2d.py` 是本仓库自实现的 MIT 许可二维 DEM 参考后端（Cundall–Strack 接触）；
- `dem_backend.py` 只定义 `dem/v1` 独立进程 JSON 通信协议与 `ExternalDEMBackend`
  契约基类，由使用者在**自己的环境**提供外部脚本；
- 外部脚本既不与主程序链接，也不随本发布包分发。用户自行接入 Yade/PFC 时，
  须自行遵守相应 GPL/商业许可义务。

## 四、GUI 技术路线说明

实施方案原拟 PySide6（LGPL 动态链接）；因银河麒麟 aarch64 无预编译 wheel，
实际采用 **Python 标准库 `http.server` + Firefox/浏览器** 的本地 Web 界面，
未引入 PySide6/Qt，故无 LGPL 动态链接场景。此变更已在 `README.md` 与
`运行说明.md` 记录。

## 五、数据与理论来源

- 专著《浅埋土洞坍塌》第 19 章五案例的参数与数值结果，属学术引用；
  案例模板已标注"典型取值示例、非场地实测"，不随程序授权转让其著作权。
- 参考文献（Terzaghi、Davis、Cundall–Strack、Hoek–Brown 等）仅在文档与
  代码注释中注明出处，不复制受版权保护的原文图表。

## 六、声明边界

本文件按已安装依赖的元数据如实记录，不声称覆盖使用者后续自行添加的组件；
新增依赖时须同步更新本声明与 `requirements.txt`。商用授权与法务梳理见
`授权说明.md`。
