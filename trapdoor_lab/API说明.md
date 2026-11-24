# 活门试验数值模拟程序 · API 说明

服务地址：`http://127.0.0.1:8765`（端口被占用时自动顺延，以终端打印的地址为准）。
全部请求/响应均为 UTF-8 JSON。程序只用 Python 标准库 `http.server`，无第三方 Web 框架。

## 一、GET 接口

| 路径 | 说明 | 响应 |
| --- | --- | --- |
| `/` 或 `/index.html` | 内嵌 Web 界面 | HTML |
| `/api/health` | 存活探测 | `{"ok":true}` |
| `/api/tables` | 论文 Table 3 / Table 4 复现数据 | `{"table3":{...},"table4":{...}}` |
| `/api/bench2` | 二期文献基准集（1g / 离心机 / 速率适用域） | `{"anchors":{...},"check":{...}}` |

`/api/bench2` 无需请求体，返回三组共 19 项文献锚点（`anchors`）及其机械自洽核验
（`check`，含 `all_pass`、`checks`、`counts`、`beta_half_span_pct`）。
数据逐条摘自专著，`note` 明确「文献报告值，程序不拟合、不推导」。

## 二、POST 接口

### 1. `/api/solve` —— L1 准静态闭式解

```json
请求：{"params": {"sigma_c":100,"A":0.6,"B":0.6,"sigma_t":5,
                  "gamma":18,"q2":10,"l":2.5,"h":3}, "mode":"active"}
```

`mode` 取 `active`（主动塌落）或 `passive`（被动顶起）。

响应主要字段：

- `result`：求解结果，含 `solution_type`（`power` / `linear-limit`）、
  `psi0_deg`、`psi_h_deg`、`n`、`omega`/`Omega`、`Lambda`、`h_cr`、
  `q1a`、`q1b`、`Ra`、`over`、`is_global_collapse`、`judgement_text`、`warnings` 等。
- `figures`：六张图的 base64 PNG，键为 `fig1`–`fig6`。
- `provenance`：可复现性元数据（参数快照、`git_commit`、`code_sha256`、
  依赖版本、时间戳）。

### 2. `/api/l2` —— L2 打靶—积分

请求同 `/api/solve`。响应 `result` 含 `converged`、`iterations`、`x2`、`L_c`、
`closed_form_Lc`、`lcf_error`、`residual`、`screening`（方向/耗散/边界三级筛选）、
`admissible`；`fig` 为打靶解与闭式解对照图。要求 `0<B<1`，B=1 时返回错误并走 L1 线性分支。

### 3. `/api/l3` —— L3 速率过程层

```json
请求：{"params": {"h":30,"B_c":1,"d":0.01,"rho_b":1600,"C_b":0.6,
                  "k_d":1.5,"phi_r":35,"mu_d":0.70,"I0":0.30,
                  "xi":1.0,"tau0":0.0,"g":9.81},
       "phi_deg": 30.96}
```

响应 `result` 含 `Q_vol`、`v`、`T`、`R`、`I`、`mu_chim`、`h_over_D`、
`domain`、`solution_type`（超界为 `rate-degraded`）、`timeline`、`stages`、
`self_stable`、`warnings`；`fig` 为时间线与五阶段图。

### 4. `/api/gtt` —— 广义切线法

```json
请求：{"params": {同 /api/solve}, "mode":"active", "sigma0": 40.0}
```

`sigma0` 可省略或为 `null`，此时自动取滑面正应力加权均值。
响应 `result` 含：

- `tangent`：`tan_delta`、`delta_deg`(=ψ)、`c_t`、`tau_at_point`、
  `selfcheck_residual`、`c_t_formula_residual`；
- `table_check`：表 13.2 逐项核对，`rounding_rows` 记录原书印刷舍入差；
- `mechanism`：`sigma_bar`、`sigma0_opt`、`theorem13_2_rel_err`、
  `D_true`/`D_upper`/`D_lower`、`boundary_note`、`scan_sigma`/`scan_D_upper`；
- `fig`：包线夹逼图与定理 13.2 扫描图。

**边界**：不输出 Q^UB，专著未给机构荷载表达式，本接口不虚构严格荷载上界。

### 5. `/api/mirror` —— 主动/被动滑面镜像对照

请求 `{"params": {同 /api/solve}}`。响应 `result` 含 `x2_active`、`x2_passive`、
`x2_sum`（应等于 `two_l=2l`）、`max_mirror_error_m`（400 点网格线性插值残差，
约 1e−5 m）、`mirror_ok`、`clipped`、`note`；`fig` 为主动/被动/镜像虚线对照图。

注意：当主动出露点被裁剪到 x=0（线性分支 dx>l）时，几何上不再严格镜像，
`clipped=true`、`mirror_ok=false`，`note` 会说明原因，不假装对称。

### 5b. `/api/anim` —— 二维坍塌过程动画

请求：`{"params": {同 /api/solve}, "n_frames": 48, "passive": false, "phi_deg": null}`。
响应 `result` 含 `meta`（`solution_type="kinematic-animation"`、`T`、`v`、`R_final`、
`boundary_note`）、`geometry`（`l/h/D/x_outcrop/R_final`）、`frames`（逐帧）。

每帧含 `t`、`z_front`（前锋自洞底推进高度）、`door_down`（活门下移）、
`crater_r`（地表漏斗半径）、`stage`、`surface`（地表折线）、`slip`（滑面折线）、
`particles`（`[x,z,flag]`，flag=1 表示该粒子运动中）。前端据此画 SVG 并播放。

**边界**：运动学演示——粒子位移是满足滑面约束与体积守恒的构造模板，
不是 DEM 颗粒接触求解，也不是有限元应力云图；不输出应力/塑性区数值（属 M7 远期）。

### 6. `/api/sat` —— 文献饱和高度锚点

请求 `{"params": {同 /api/solve}}`。
响应含 `N_sat`、`h_sat_over_D_printed`/`h_sat_over_D_exact`、`consistency_ok`、
`h_cr_m`、`h_cr_over_D`、`thresholds_distinct`、`rows`（不同 h/D 的外推与适用域）、
`boundary_note`。该接口**不推导** N_sat，只做常数自洽与门槛区分。

## 三、错误格式

参数非法或内部异常统一返回 `{"error": "错误说明"}`，HTTP 状态码仍为 200，
前端据此在界面上显示红色提示。典型错误：

- `B 必须在 (0,1]`（L2 打靶）；
- `sigma_c 必须 > 0`、`sigma_0 + sigma_t 必须 > 0`；
- `内部错误: ...`（非预期异常，需反馈日志）。

## 四、调用示例（curl）

```bash
curl -s http://127.0.0.1:8765/api/health
curl -s -X POST http://127.0.0.1:8765/api/solve \
  -H 'Content-Type: application/json' \
  -d '{"params":{"sigma_c":100,"A":0.6,"B":0.6,"sigma_t":5,"gamma":18,"q2":10,"l":2.5,"h":3},"mode":"active"}'
```

## 五、可复现性

每次 `/api/solve` 的 `provenance` 随结果返回，并写入「导出计算报告 (HTML)」。
`code_sha256` 覆盖全部计算模块；`git_commit` 读不到时写 `unknown`，不伪造。
