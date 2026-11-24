# 发布说明（RELEASE）

## 版本

- **版本号**：v1.1.0（M9 工程商用预留版）
- **发布日期**：2026-11-24
- **目标平台**：银河麒麟桌面操作系统 V10 SP1（aarch64）；亦兼容 x86_64 Linux
- **基线提交**：见下方「构建记录」，以实际 `git rev-parse HEAD` 为准

> v1.0 为方案 M1–M6 的基线（tag `v0.1-baseline` 及其后续）；
> v1.1.0 补齐 M7（FELA 上下界）、M8（DEM/三维锥形）与 M9（案例/授权/发布包）。

## 本次发布内容

| 类别 | 内容 |
|---|---|
| 专著理论主线 | L1 闭式解、L2 打靶—积分、L3 速率时间线、GTT 广义切线法、h_cr/h_sat 锚点、主动/被动镜像 |
| 可视化 | 六图、参数扫描动态图、二维坍塌运动学动画、FEM 静态与动态应力云图 |
| 力学场扩展 | Q4-DP 弹塑性 FEM、LP-FELA 上下界夹逼、二维 DEM 颗粒流、三维锥形极限平衡 |
| 工程案例 | 专著第 19 章五案例模板（案例一含两土性变体），参数与数值逐条对齐原书 |
| 工程化 | golden 基准、CI 回归、可复现性元数据、五份技术文档、授权与许可证声明、发布清单 |

## 两套运行方式

### A. 源码运行（推荐，银河麒麟 aarch64 已验证）

```bash
cd trapdoor_lab
bash 一键启动.sh          # 自动检测依赖、启动本地 Web 服务并打开浏览器
```

依赖：`numpy 1.24.4`、`scipy 1.10.1`、`matplotlib 3.7.5`（见 `requirements.txt`）。
服务默认 8765，被占用时自动顺延到 8766–8784，实际地址以终端打印为准。

自检与回归：

```bash
../venv/bin/python 自检.py                       # 13 项基准自检
../venv/bin/python -m unittest discover -s tests # 全量回归
```

### B. 独立打包版（无需 Python 环境）

```bash
cd trapdoor_lab
bash build.sh                       # 需先安装 PyInstaller
cd dist/trapdoor-lab
bash 启动打包版.sh                   # 浏览器打开终端显示的地址
```

注意：打包版是源码的快照，改源码后必须重跑 `build.sh`。旧 `dist/` 存在时
`build.sh` 会提示手动清理（环境安全策略不允许脚本批量删除）。

## 完整性校验

项目根目录 `MANIFEST.sha256` 列出全部交付文件的 SHA-256：

```bash
sha256sum -c MANIFEST.sha256        # 在仓库根目录执行
```

任何文件被改动都会报告 FAILED。清单由 `trapdoor_lab/make_manifest.py` 生成，
新增交付文件后须重新生成并随版本一并提交。

## 已知边界（发布即声明）

1. L1/L2/L3 与 GTT 为专著理论主结果；FEM/FELA/DEM/三维模块为按用户授权横向
   引用的扩展，非专著原有理论，定性—半定量，不冒充严格上下界。
2. 五案例参数为**典型取值示例、非场地实测**，不得直接用于工程设计。
3. DEM 外部后端（Yade GPL / PFC 商业）仅定义独立进程协议，不链接、不分发。
4. 未实现：完整 FELA SOCP、Yade/PFC 求解器本身、三维 FEM/FELA、μ(I) 细观反算、
   局部坍塌形态数据。
5. 程序输出为工程辅助分析结果，不替代勘察、试验与注册工程师签章（MIT "AS IS"）。

## 构建记录（发布时填写）

```
version      : v1.1.0
git commit   : 3fad0a2（M9 功能提交）
python       : 3.8.10 (aarch64 venv) / CI 3.11 (x86_64)
dependencies : numpy 1.24.4 / scipy 1.10.1 / matplotlib 3.7.5
regression   : 146/146（含 18 项工程案例）
self-check   : 14/14
CI           : success，run 45375058015（全量单测/自检/编译/JS 语法四步全过）
```
