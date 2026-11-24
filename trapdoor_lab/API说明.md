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
| `/api/cases` | 工程案例索引（五案例元信息 + 汇总 + 自检） | `{"cases":[...],"summary":{...},"selfcheck":{...}}` |

`/api/bench2` 无需请求体，返回三组共 19 项文献锚点（`anchors`）及其机械自洽核验
（`check`，含 `all_pass`、`checks`、`counts`、`beta_half_span_pct`）。
数据逐条摘自专著，`note` 明确「文献报告值，程序不拟合、不推导」。

`/api/cases` 无需请求体，供前端下拉框一次性拉取。`cases[]` 每项含
`id/name/source/scenario/mode/morphology/judgement/caution`（不含参数，参数在 `/api/case` 中）；
`summary` 含 `n_cases`、`all_comparison_ok`、`selfcheck_ok`、`width_over_span`、
`active_range`、`key_warning`、`caveat`；`selfcheck` 为 `selfcheck_cases()` 的逐项布尔结果。

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

### 5c. `/api/fem` —— 二维弹塑性 FEM 力学场

请求：
```json
{"params": {"l":2.5,"h":3.0,"gamma":18,"q2":0},
 "fem_params": {"E":20000,"nu":0.3,"c":25,"phi":20,"psi":0},
 "opts": {"nx":12,"ny":9,"n_steps":12,"sr":1.0}}
```

`fem_params` 可省略（用内置默认）；`opts` 可省略（默认 12×9、12 步）。
`sr` 为强度折减系数：`c'=c/sr`，`tanφ'=tanφ/sr`。

响应：
- `result.meta`：`method`、`source`、`boundary`、`converged`、`sr`、`c_eff`、
  `phi_eff_deg`、`grid`、`n_elements`、`n_nodes`；
- `result.curve`：每步 `step`、`door_down`、`reaction_kN`、`q_kPa`、
  `plastic_frac`、`iters`、`residual`；step0 为未扰动静止压力对照点；
- `result.field`：`centers`（单元中心 `[x,y]`）、`sxx`/`syy`/`txy`/`p`/`q`、
  `epsp`、`yielded`、`plastic_frac`；
- `result.peak_q_kPa` / `peak_step`、`params`；
- `fig`：base64 PNG，三面板（偏应力 q 云图 / 塑性区 / 压力—位移曲线）；
- `selfcheck`：DP 本构自检结果。

**边界**：横向引用的成熟有限元（平面应变 Q4 + DP 径向返回 + 拉力截断 + 位移控制活门），
按用户授权引入；不是严格 FELA 上下界，DP 圆化与专著幂型 Hoek–Brown 不同族，
小变形理想弹塑性、网格较粗，定性—半定量演示。截面位移步内以弹性刚度修正 Newton
（默认松弛 0.4、容差 1e-3），`converged` 为工程收敛标志，非严格数学收敛证明。

### 5d. `/api/fem_anim` —— FEM 应力场动态云图

请求与 `/api/fem` 相同；`opts` 默认 `nx=8, ny=6, n_steps=8`，
另可传 `anim_opts`（如 `{"fps":4}`）。

响应：
- `result.frames`：逐帧 PNG 的 base64 列表，长度 `n_steps+1`（步 0…n_steps）；
- `result.gif`：同组帧生成的 GIF 动图（base64），自动循环播放；
- `result.meta`：`n_frames`、`fps`、`qmax`（全场统一色标上限，避免帧间漂移）、
  `boundary`、`note`；
- `curve`、`peak_q_kPa`、`converged`：对应 FEM 求解的压力—位移曲线与工程收敛标志。

每帧两面板：左为偏应力 `q_dev` 云图（全场统一色标），右为塑性区（红=屈服单元），
标题给出步号、活门下移量、门压 `p_door`。前端用滑块逐帧查看，GIF 供自动播放。

**边界**：底层是 FEM 数值场（非运动学演示），但仍是定性—半定量扩展：
DP 圆化、小变形理想弹塑性、非严格 FELA 上下界。

### 5e. `/api/fela` —— FELA 极限分析上下界夹逼（方案 M7）

POST，请求：
```json
{"params": {"l":2.5, "h":3.0, "gamma":18, "q2":0},
 "opts":   {"K":12, "n_points":41, "psi_deg":10}}
```
- `params` 可只传需覆盖的键，其余自动补全自 `DEFAULT_PARAMS`；
- `opts.K`：屈服包线线性化段数（默认 12，范围 2–64）；`n_points`：滑面离散点数；
  `psi_deg`：剪胀角（可选，仅用于 Davis 折减报告，`>φ₀` 自动钳制为关联）。

响应：
- `result.resistance`：`Ra_LB/Ra_exact/Ra_UB`（抗力上下界与 L1 精确值）、
  `Ra_interval_ok`、`q_active_LB/UB/width`、`h_over_D`、`applicable`、`davis`；
- `result.bracket`：`ok`、`q_LB_kPa`、`q_UB_kPa`（主动支护压力夹逼区间）；
- `result.meta`：`regime`（`shallow-applicable` / `out-of-range`）、实际 `K`/`n_points`；
- `band`：`K` 序列及对应 `q_LB`/`q_UB`/`width`，供夹逼带图；
- `fig`：`q` 上下界随 K 收窄 + 区间宽度对数下降双面板 PNG（base64）。

数学口径：内接多面体→静力下限抗力，外接→运动上限抗力；主动支护压力
`q=over−Ra/l` 是 `Ra` 的减函数，故 q 夹逼方向反转，L1 闭式解落在 `[q_LB, q_UB]` 内。

**边界**：线性化多面体极限分析，不是 SOCP 圆锥规划、不是幂型精确解；
下限与 ψ 无关，Davis 折减单独报告（不混入数学夹逼）；深埋 `h/D>1.4`
标记 `out-of-range`，不给出夹逼结论；`N_sat=18.2212` 仅为文献锚点，不由此推导。

### 5f. `/api/dem` —— 二维 DEM 颗粒流（方案 M8）

POST，请求：
```json
{"params": {"l":2.5, "h":3.0, "n_particles":360, "mu":0.5},
 "opts":   {"settle_steps":150, "n_steps":500, "save_every":20, "max_particles":1200}}
```
- `params` 可只传覆盖键，其余补全 `DEFAULT_DEM_PARAMS`（含 `kn`、`mu`、`rho_p` 等）；
- `opts.settle_steps`：活门关闭预沉积步数；`n_steps`：开闸后步数；
  `save_every`：存帧间隔；`max_particles`：粒子上限（防界面卡死，默认 1200）。

响应：
- `result.frames`：逐帧颗粒轨迹。每帧 `{step,t,door_open,particles,radii,n_contacts,
  door_force_N,kinetic_J,max_overlap,surface}`，`particles` 为 `[[x,y,vx,vy],...]`（SI）；
- `result.curve`：沉积/流动两阶段的反力、接触数、动能时间序列；
- `result.meta`：`n_particles`、`stable`、`peak_overlap_ratio`、
  `final_overlap_ratio`、`kinetic_peak_J` 等；
- `selfcheck`：DEM 内核 10 项自检；`fig`：四阶段速度着色 PNG（base64）。

**边界**：二维圆盘 Cundall–Strack 接触，无转动自由度；刚度/摩擦为数值参数，
非 μ(I) 细观反算；真·颗粒接触求解（非运动学模板），但横向引用、非专著理论、定性—半定量。

### 5g. `/api/dem_backend` —— DEM 独立进程后端协议（dem/v1）

POST，请求：`{"params":{...}, "opts":{...}, "backend":"local", "command":[...]}`。
`backend=local` 直接调用自实现 MIT 参考后端；`backend=external` 时必须提供
`command`（外部脚本 argv，脚本从 stdin 读 v1 请求、向 stdout 写 v1 响应）。

响应：`response` 为 v1 响应对象（`protocol/backend/status/result`），`selfcheck` 为 9 项适配器自检。

**许可证边界**：Yade(GPL)/PFC(商业) 仅以独立进程 + 文件/stdio 通信接入，
不链接、不随包分发；本模块只定义契约，不含任何 GPL/商业求解代码。

### 5h. `/api/cone` —— 三维锥形极限平衡（方案 M8 配套）

POST，请求：`{"params":{"l":2.5,"h":3.0,"phi_deg":25,"c":10,"theta_deg":20,"K":0.6}}`。
响应：
- `results`：三形态（funnel/urn/shaft）各自的 `meta`（morph/r0/h…）、
  `result`（`weight_kN`/`side_resistance_kN`/`support_pressure_kPa`/
  `self_stable`/`mobilized_fraction`/`r_top`/`r_max`）、`profile`（z,r）与 `comparison`；
- `critical_depth`：二分求 `q(h)=0` 的临界埋深（无根时 `has_root=false`）；
- `selfcheck`：12 项内核自检；`fig`：三形态剖面 + 支护压力柱状图（base64）。

**边界**：简化三维轴对称刚体滑体极限平衡，非三维 FEM/FELA/SOCP；三形态命名对齐专著，
但轴对称推广为本程序扩展、非专著原有公式；与 L1 二维解只作量级对照，不判同一性。

### 5i. `/api/case` —— 工程案例复算（方案 M9）

POST，请求：`{"id":"karst-cavity"}`。合法 id 为
`karst-cavity`、`urban-road`、`goaf`、`shield-tunnel`、`tunnel-uplift`。
响应：
- `result.case`：案例元信息（id/name/source/scenario/mode/morphology/judgement/caution）；
- `result.params`：模板参数（含 `psi0_deg`；主动/渐进为非关联角，顶起为 null）；
- `result.result`：复算结果。主动含 `Omega/A_psi/S_psi/p/L_c/x2/h_cr/surface_width/
  global_collapse/branch_valid/selfcheck_Lc_at_hcr/selfcheck_error/
  nonassociated_reduction`；顶起含 `S/p/Delta/x2/q1/h_cr=null`；
  渐进含 `closure_sigma_at_surface/closure_exceed_kPa/closure_inadmissible`；
- `result.comparison`：逐键对照，每项 `{computed,book,abs_error,rel_error,ok}`；
- `result.variants`：案例一 A/B 变体，同结构；`result.all_comparison_ok` 为总判定；
- `result.boundary`：适用边界文字；`fig`：复算值 vs 原书值对照条形图（base64 PNG）。

复算公式：式(16.32) `A_ψ`、式(16.33) `S_ψ`、式(16.19a) `L_c`、式(19.1) `h_cr`、
式(19.12)–(19.14) 顶起；每次主动复算强制反验 `L_c(h_cr)=l`。
**边界**：五案例参数均为专著典型示例、非场地实测；以实测 ψ 反定 S_ψ 为工程近似、非严格界值。

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
curl -s http://127.0.0.1:8765/api/cases
curl -s -X POST http://127.0.0.1:8765/api/case \
  -H 'Content-Type: application/json' -d '{"id":"karst-cavity"}'
curl -s -X POST http://127.0.0.1:8765/api/solve \
  -H 'Content-Type: application/json' \
  -d '{"params":{"sigma_c":100,"A":0.6,"B":0.6,"sigma_t":5,"gamma":18,"q2":10,"l":2.5,"h":3},"mode":"active"}'
```

## 五、可复现性

每次 `/api/solve` 的 `provenance` 随结果返回，并写入「导出计算报告 (HTML)」。
`code_sha256` 覆盖全部计算模块；`git_commit` 读不到时写 `unknown`，不伪造。
