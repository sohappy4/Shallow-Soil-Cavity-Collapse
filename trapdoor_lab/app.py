# -*- coding: utf-8 -*-
"""
app.py —— 活门试验数值模拟程序 本地 Web 服务
=================================================
只用 Python 标准库 http.server，不依赖 Flask/PySide6（麒麟 aarch64 装不上）。
浏览器打开 http://127.0.0.1:8765 即可使用。

运行：  python3 app.py
退出：  终端按 Ctrl+C
"""
from __future__ import annotations

import json
import io
import os
import sys
import base64
import hashlib
import platform
import subprocess
import datetime
import webbrowser
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def provenance(params, mode: str = "active") -> dict:
    """可复现性元数据（方案第七章要求）：参数快照 + 代码指纹 + 依赖版本 + 时间。

    git commit 读不到时写 unknown，不伪造；代码指纹取全部计算模块的 sha256。
    """
    modules = ["app.py", "collapse_kernel.py", "plots.py", "l2_shooting.py",
               "l3_rate.py", "gtt.py", "benchmarks.py"]
    digest = hashlib.sha256()
    for name in modules:
        path = os.path.join(HERE, name)
        if os.path.exists(path):
            with open(path, "rb") as f:
                digest.update(name.encode("utf-8"))
                digest.update(f.read())
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=os.path.dirname(HERE),
            stderr=subprocess.DEVNULL).decode("utf-8").strip()
    except Exception:
        commit = "unknown"
    import numpy
    import scipy
    import matplotlib
    return dict(
        timestamp=datetime.datetime.now().isoformat(timespec="seconds"),
        params=dict(params), mode=mode,
        git_commit=commit,
        code_sha256=digest.hexdigest(),
        dependencies=dict(python=platform.python_version(),
                          numpy=numpy.__version__, scipy=scipy.__version__,
                          matplotlib=matplotlib.__version__),
    )


def _fig_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from collapse_kernel import (
    DEFAULT_PARAMS, PRESETS, solve, figure_data, reproduce_tables, ParamsError,
)
from plots import render_all
from l2_shooting import solve_shooting
from l3_rate import solve_rate, DEFAULT_RATE_PARAMS
from gtt import (
    tangent_params, table_13_2, check_table_13_2,
    dissipation_bracket, mechanism_scan, tau_envelope, chord_tau,
)
from benchmarks import (
    saturation_anchor, LITERATURE_ANCHORS, check_literature_anchors,
)

PORT = int(os.environ.get("TRAPDOOR_PORT", "8765"))

HTML = r"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>浅埋土洞活门试验数值模拟程序</title>
<style>
  :root { --bg:#f5f6f8; --card:#fff; --ink:#1c1e22; --muted:#666; --line:#dde1e6;
          --blue:#2f6fb0; --warn:#b36b00; --ok:#1a7f37; }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--bg); color:var(--ink);
         font-family:"Noto Sans CJK SC","Microsoft YaHei",sans-serif; font-size:14px; }
  header { background:#fff; border-bottom:1px solid var(--line); padding:14px 22px; }
  header h1 { margin:0; font-size:19px; }
  header p { margin:4px 0 0; color:var(--muted); font-size:12px; }
  .wrap { display:grid; grid-template-columns:320px 1fr; gap:16px; padding:16px 22px; }
  .card { background:var(--card); border:1px solid var(--line); border-radius:8px; padding:14px 16px; }
  h2 { font-size:15px; margin:0 0 10px; padding-bottom:8px; border-bottom:1px solid var(--line); }
  label { display:block; margin:9px 0 3px; font-size:12px; color:var(--muted); }
  input, select, button { width:100%; padding:7px 9px; border:1px solid var(--line);
                          border-radius:5px; font-size:13px; background:#fff; font-family:inherit; }
  button { cursor:pointer; background:var(--blue); color:#fff; border:0; font-weight:600; margin-top:12px; }
  button.ghost { background:#eef2f7; color:var(--ink); margin-top:6px; }
  .row { display:grid; grid-template-columns:1fr 1fr; gap:8px; }
  .kv { display:grid; grid-template-columns:1fr 1fr; gap:6px 12px; font-size:13px; margin:6px 0; }
  .kv span:nth-child(odd) { color:var(--muted); }
  .banner { padding:10px 12px; border-radius:6px; margin:10px 0; font-size:13px; }
  .banner.ok { background:#e8f5ec; color:var(--ok); }
  .banner.warn { background:#fff4e0; color:var(--warn); }
  .grid { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
  figure { margin:0; background:#fff; border:1px solid var(--line); border-radius:8px; padding:8px; }
  figure img { width:100%; display:block; }
  figcaption { text-align:center; font-size:12px; color:var(--muted); padding:4px 0; }
  table { border-collapse:collapse; width:100%; font-size:12px; margin-top:8px; }
  th,td { border:1px solid var(--line); padding:5px 7px; text-align:right; }
  th { background:#f0f2f5; }
  .note { color:var(--muted); font-size:12px; line-height:1.7; }
  details { margin-top:12px; }
  summary { cursor:pointer; color:var(--blue); }
  #err { color:#b00020; font-size:12px; min-height:16px; }
</style>
</head>
<body>
<header>
  <h1>浅埋土洞活门试验数值模拟程序 · L1 闭式解内核</h1>
  <p>理论依据：专著《浅埋土洞坍塌：形态、临界埋深与统一变分理论》 / JRMGE 配套论文　·　单位：kPa、m、kN/m³</p>
</header>
<div class="wrap">
  <div class="card">
    <h2>工况与参数</h2>
    <label>预设土性</label>
    <select id="preset">
      <option value="">自定义 / 基准算例</option>
      <option>密实砂</option><option>松散砂</option><option>软黏土</option><option>结构性土</option>
    </select>
    <label>工况</label>
    <select id="mode"><option value="active">主动塌落 (u&gt;0)</option>
    <option value="passive">被动顶起 (u&lt;0)</option></select>
    <div class="row">
      <div><label>σ_c 强度 (kPa)</label><input id="sigma_c" type="number" step="any"></div>
      <div><label>A 系数</label><input id="A" type="number" step="any"></div>
      <div><label>B 指数</label><input id="B" type="number" step="any"></div>
      <div><label>σ_t 抗拉 (kPa)</label><input id="sigma_t" type="number" step="any"></div>
      <div><label>γ 重度 (kN/m³)</label><input id="gamma" type="number" step="any"></div>
      <div><label>q₂ 超载 (kPa)</label><input id="q2" type="number" step="any"></div>
      <div><label>l 洞半宽 (m)</label><input id="l" type="number" step="any"></div>
      <div><label>h 覆土 (m)</label><input id="h" type="number" step="any"></div>
    </div>
    <div id="err"></div>
    <button onclick="run()">开始计算</button>
    <button class="ghost" onclick="reset()">恢复基准算例</button>
    <button class="ghost" onclick="exportReport()">导出计算报告 (HTML)</button>
    <details><summary>试验对照（论文 Table 3 / 4）</summary><div id="tables" class="note"></div></details>
    <hr style="border:0;border-top:1px solid var(--line);margin:14px 0">
    <h2>L2 打靶—积分 / L3 速率链</h2>
    <button class="ghost" onclick="runL2()">L2 打靶求滑动面</button>
    <button class="ghost" onclick="runL3()">L3 速率时间线</button>
    <button class="ghost" onclick="toggleRate()">显示/隐藏 L3 参数</button>
    <div id="rateparams" style="display:block;margin-top:8px">
      <label>B_c 烟囱宽度 (m)</label><input id="rate_Bc" type="number" step="any" value="1.0">
      <label>d 粒径 (m)</label><input id="rate_d" type="number" step="any" value="0.01">
      <label>ρ_b 堆密度 (kg/m³)</label><input id="rate_rho" type="number" step="any" value="1600">
      <label>C_b 出流系数</label><input id="rate_Cb" type="number" step="any" value="0.6">
      <label>φ_r 休止角 (°)</label><input id="rate_phir" type="number" step="any" value="35">
    </div>
    <hr style="border:0;border-top:1px solid var(--line);margin:14px 0">
    <h2>GTT 广义切线法</h2>
    <label>切点 σ₀ (kPa，留空则自动取滑面均值)</label><input id="gtt_sigma0" type="number" step="any" placeholder="自动">
    <button class="ghost" onclick="runGTT()">GTT 切点 / 表13.2 / 定理13.2</button>
    <button class="ghost" onclick="runSat()">文献锚点 h_sat / h_cr</button>
    <button class="ghost" onclick="runBench2()">二期文献基准集 (1g/离心机)</button>
    <button class="ghost" onclick="runMirror()">主动/被动镜像对照</button>
    <hr style="border:0;border-top:1px solid var(--line);margin:14px 0">
    <h2>参数扫描动态图件</h2>
    <label>扫描参数</label>
    <select id="scan_param">
      <option value="h">覆土 h (m)</option><option value="A">系数 A</option>
      <option value="B">指数 B</option><option value="sigma_t">抗拉 σ_t (kPa)</option>
    </select>
    <div class="row">
      <div><label>起点</label><input id="scan_min" type="number" step="any" value="1"></div>
      <div><label>终点</label><input id="scan_max" type="number" step="any" value="8"></div>
    </div>
    <label>帧数</label><input id="scan_steps" type="number" value="16" min="2" max="60">
    <button class="ghost" onclick="toggleScan()">播放 / 停止扫描</button>
    <div id="scaninfo" class="note" style="margin-top:6px"></div>
  </div>

  <div class="card">
    <h2>计算结果</h2>
    <div id="banner"></div>
    <div class="kv" id="kv"></div>
    <div class="grid" id="figs"></div>
    <div id="extra" class="note" style="margin-top:14px"></div>
    <div id="player" style="display:none;margin-top:12px;border:1px solid var(--line);border-radius:6px;padding:10px">
      <b>L3 时间线播放</b>
      <div style="display:grid;grid-template-columns:230px 1fr;gap:12px;align-items:center">
        <svg id="l3svg" viewBox="0 0 220 260" width="220" height="260"></svg>
        <div>
          <input id="l3slider" type="range" min="0" max="1000" value="0" oninput="drawL3Frame()" style="padding:0">
          <button class="ghost" onclick="togglePlay()">播放 / 暂停</button>
          <div id="l3frameinfo" class="note" style="margin-top:8px"></div>
        </div>
      </div>
    </div>
    <div class="note" style="margin-top:14px">
      注：埋深比 n = h/(2l) = h/D（洞全宽 D=2l）；q1a/q1b 为程序输出，非输入。
      B=1 自动切换线性极限分支（h_cr = l/A），不与闭式解比较。
      被动工况按专著不设 h_cr 门控。
      L2 为专著算法 17.1 的 RK4 + 二分实现；L3 为式(18.1)–(18.7) 闭式链，
      h/D&gt;1.4 自动标记为超界演示并提示退化到静力解。
      GTT 只给出切线参数、耗散积分夹逼与切点判定；专著未给出机构荷载 Q<sup>UB</sup>(σ₀)
      的表达式，本程序不虚构严格荷载上界，也不把夹逼冒充 FELA 上下界。
    </div>
  </div>
</div>
<script>
const KEYS = ["sigma_c","A","B","sigma_t","gamma","q2","l","h"];
const PRESETS = __PRESETS__;
const DEFAULT = __DEFAULT__;
let LAST = null, LAST_L3 = null, PLAY_TIMER = null;
let SCAN = {timer:null, i:0, busy:false};

function fill(p){ KEYS.forEach(k=>document.getElementById(k).value = p[k]); }
function read(){
  const p={}; for(const k of KEYS){ const v=parseFloat(document.getElementById(k).value);
    if(!isFinite(v)) throw new Error("参数 "+k+" 不是有效数字"); p[k]=v; }
  return p;
}
function reset(){ document.getElementById("preset").value=""; fill(DEFAULT);
  document.getElementById("err").textContent=""; }
window.onload=()=>{ fill(DEFAULT); run(); };

document.getElementById("preset").onchange=e=>{
  const v=e.target.value; if(v && PRESETS[v]) fill(PRESETS[v]); };

async function run(){
  document.getElementById("err").textContent=""; hidePlayer();
  let payload;
  try{ payload={params:read(), mode:document.getElementById("mode").value}; }
  catch(e){ document.getElementById("err").textContent=e.message; return; }
  const r = await fetch("/api/solve",{method:"POST",headers:{"Content-Type":"application/json"},
    body:JSON.stringify(payload)});
  const data = await r.json();
  if(data.error){ document.getElementById("err").textContent=data.error; return; }
  LAST = data;
  const res = data.result;
  const warn = res.warnings.length>0;
  document.getElementById("banner").innerHTML =
    `<div class="banner ${res.is_global_collapse?'ok':'warn'}"><b>判定：</b>${res.judgement_text}</div>` +
    (warn?`<div class="banner warn"><b>提示：</b>${res.warnings.join(" ")}</div>`:"");
  const rows = [
    ["求解器类型", res.solution_type==="linear-limit"?"线性极限 (B=1)":"幂型闭式解"],
    ["组合强度 S", fmt(res.S)], ["出露参数 t", fmt(res.t)],
    ["脆性比 Ω", fmt(res.Omega)], ["覆强比 Λ", fmt(res.Lambda)],
    ["埋深比 n = h/(2l)", fmt(res.n)],
    ["地表剪胀角 ψ₀", fmt(res.psi0_deg,2)+"°"], ["底部剪胀角 ψ(h)", fmt(res.psi_h_deg,2)+"°"],
    ["临界覆土 h_cr", fmt(res.h_cr,3)+" m"],
    ["主动支护 q1a", fmt(res.q1a,3)+" kPa"], ["被动顶起 q1b", fmt(res.q1b,3)+" kPa"],
    ["上覆压力 q₂+γh", fmt(res.over,3)+" kPa"],
    ["镜面对称误差", res.mirror_error.toExponential(2)]
  ];
  document.getElementById("kv").innerHTML = rows.map(x=>`<span>${x[0]}</span><span>${x[1]}</span>`).join("");
  const names={fig1:"Fig.1 滑动面形态",fig2:"Fig.2 滑动面应力场",fig3:"Fig.3 剪胀角演化",
    fig4:"Fig.4 形态分区与工况落位",fig5:"Fig.5 侧压力系数 Ka(z)",fig6:"Fig.6 非关联修正"};
  document.getElementById("figs").innerHTML = Object.keys(names).map(k=>
    `<figure><img src="data:image/png;base64,${data.figures[k]}" alt="${names[k]}">
     <figcaption>${names[k]}</figcaption></figure>`).join("");
}
function fmt(v,d=4){ return v===null||v===undefined?"—":Number(v).toFixed(d); }

// 参数扫描动态图件：复用 /api/solve 与 run() 的渲染，逐帧改一个输入参数。
// 注意 L1 闭式解只有沿滑面的一维场，因此这里是"动态图件"而非二维应力云图；
// 二维云图需要 FELA 层，程序不虚构。
async function scanStep(){
  if(SCAN.busy) return;
  SCAN.busy=true;
  const key=document.getElementById("scan_param").value;
  const lo=parseFloat(document.getElementById("scan_min").value);
  const hi=parseFloat(document.getElementById("scan_max").value);
  const steps=Math.max(2,Math.min(60,parseInt(document.getElementById("scan_steps").value)||16));
  if(!isFinite(lo)||!isFinite(hi)||hi<=lo){
    document.getElementById("err").textContent="扫描起终点无效"; SCAN.busy=false; return; }
  if(key==="B" && hi>1){
    document.getElementById("err").textContent="B 只能在 (0,1] 内扫描，超过 1 无定义"; SCAN.busy=false; return; }
  const frac=SCAN.i/(steps-1);
  const val=lo+(hi-lo)*frac;
  document.getElementById(key).value=val.toFixed(6);
  try{
    await run();
    const res=LAST&&LAST.result;
    const judge=res?(res.is_global_collapse?"整体贯通塌落":"不整体贯通"):"—";
    document.getElementById("scaninfo").innerHTML=
      "第 "+(SCAN.i+1)+"/"+steps+" 帧　"+key+"="+fmt(val,4)+
      "<br>判定："+judge+"　h_cr="+fmt(res?res.h_cr:null,3)+" m";
  }catch(e){
    document.getElementById("err").textContent="扫描出错: "+e.message;
  }finally{
    SCAN.busy=false;
  }
  SCAN.i++;
  if(SCAN.i>=steps) SCAN.i=0;
}

function toggleScan(){
  if(SCAN.timer){
    clearInterval(SCAN.timer); SCAN.timer=null;
    document.getElementById("scaninfo").innerHTML+="<br>（已停止）";
    return;
  }
  SCAN.i=0;
  scanStep();
  SCAN.timer=setInterval(scanStep, 750);
}

async function loadTables(){ const r=await fetch("/api/tables"); const d=await r.json();
  const t3=d.table3, t4=d.table4;
  let s=`<b>Table 3 主动试验 (Chen 2019)</b><br>实测 ρa=${t3.meas_Ga}：关联 ${fmt(t3.assoc_Ga,3)} (${fmt(t3.assoc_dev_pct,1)}%)，非关联 ${fmt(t3.non_Ga,2)} (${fmt(t3.non_dev_pct,1)}%)。<br><br>`;
  s+=`<b>Table 4 被动试验 (Rui 2025b)</b><table><tr><th>n</th><th>φp</th><th>θ</th><th>实测</th><th>关联</th><th>非关联</th></tr>`;
  t4.rows.forEach(x=>{ s+=`<tr><td>${x.n}</td><td>${x.phi_peak_deg}</td><td>${x.theta_deg}</td><td>${x.meas}</td><td>${fmt(x.assoc,2)}</td><td>${fmt(x.non,2)}</td></tr>`; });
  s+=`</table><br>平均偏差：关联 ${fmt(t4.assoc_mean_dev_pct,1)}%，非关联 ${fmt(t4.non_mean_dev_pct,1)}%。`;
  document.getElementById("tables").innerHTML=s;
}
loadTables();
function toggleRate(){ const e=document.getElementById("rateparams");
  e.style.display=e.style.display==="none"?"block":"none"; }

async function runL2(){
  document.getElementById("err").textContent="";
  let payload; try{ payload={params:read(), mode:document.getElementById("mode").value}; }
  catch(e){ document.getElementById("err").textContent=e.message; return; }
  const r=await fetch("/api/l2",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(payload)});
  const d=await r.json();
  if(d.error){ document.getElementById("err").textContent=d.error; return; }
  const s=d.result;
  document.getElementById("banner").innerHTML=`<div class="banner ${s.admissible?'ok':'warn'}"><b>L2 打靶：</b>${s.converged?'已收敛':'未收敛'}，迭代 ${s.iterations} 次，出露点 x₂=${fmt(s.x2,4)} m，水平内缩 L_c=${fmt(s.L_c,4)} m（闭式 ${fmt(s.closed_form_Lc,4)}，相对误差 ${s.lcf_error.toExponential(1)}）。三级筛选：方向 ${s.screening.direction?'✓':'✗'} / 耗散 ${s.screening.dissipation?'✓':'✗'} / 边界 ${s.screening.boundary?'✓':'✗'}，可容许=${s.admissible}</div>`;
  document.getElementById("kv").innerHTML=[["打靶残差",s.residual.toExponential(2)],["临界埋深 h_cr",fmt(s.h_cr_implied,3)+" m"],["Φ(s)",fmt(s.screening.phi,3)]].map(x=>`<span>${x[0]}</span><span>${x[1]}</span>`).join("");
  document.getElementById("figs").innerHTML=`<figure><img src="data:image/png;base64,${d.fig}" alt="L2"><figcaption>L2 打靶滑动面（红线为打靶解，蓝线为闭式解）</figcaption></figure>`;
  document.getElementById("extra").innerHTML="L2 按专著算法 17.1：欧拉方程通解 + RK4 积分 + x₂ 二分打靶；被动支为主动支纯镜像。";
}

async function runL3(){
  document.getElementById("err").textContent="";
  const params=read();
  const q=id=>parseFloat(document.getElementById(id).value);
  const rate={h:params.h, B_c:q("rate_Bc"), d:q("rate_d"), rho_b:q("rate_rho"),
              C_b:q("rate_Cb"), k_d:1.5, phi_r:q("rate_phir"), mu_d:0.70, I0:0.30,
              xi:1.0, tau0:0.0, g:9.81};
  if(Object.values(rate).some(v=>!isFinite(v))){ document.getElementById("err").textContent="L3 参数无效"; return; }
  const r=await fetch("/api/l3",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({params:rate, phi_deg:Math.degrees(Math.atan(params.A))})});
  const d=await r.json();
  if(d.error){ document.getElementById("err").textContent=d.error; return; }
  const s=d.result;
  const isDeg = s.solution_type.indexOf("degraded")>=0;
  document.getElementById("banner").innerHTML=`<div class="banner ${isDeg?'warn':'ok'}"><b>L3 速率链：</b>${s.domain}　h/D=${fmt(s.h_over_D,2)}　类型 ${s.solution_type}</div>`+(s.warnings.length?`<div class="banner warn">${s.warnings.join(" ")}</div>`:"");
  document.getElementById("kv").innerHTML=[["Beverloo Q",fmt(s.Q_vol,3)+" m³/s"],["推进速率 v",fmt(s.v,3)+" m/s"],["贯通时间 T",fmt(s.T,3)+" s"],["漏斗终径 R",fmt(s.R,3)+" m"],["惯性数 I",fmt(s.I,3)],["μ(I)",fmt(s.mu_chim,3)],["整柱自稳",s.self_stable?"是":"否"]].map(x=>`<span>${x[0]}</span><span>${x[1]}</span>`).join("");
  document.getElementById("figs").innerHTML=`<figure><img src="data:image/png;base64,${d.fig}" alt="L3"><figcaption>L3 塌落时间线（深度 z–时间 t）与五阶段</figcaption></figure>`;
  document.getElementById("extra").innerHTML="<b>五阶段：</b>"+s.stages.map(x=>`${x.stage}.${x.name}(t=${fmt(x.t,2)}s)`).join(" → ")+`<br>${s.stability_note}`;
  LAST_L3=s;
  document.getElementById("l3slider").value=0;
  document.getElementById("player").style.display="block";
  drawL3Frame();
}

function hidePlayer(){
  if(PLAY_TIMER){ clearInterval(PLAY_TIMER); PLAY_TIMER=null; }
  const e=document.getElementById("player"); if(e) e.style.display="none";
}

function zAtTime(t){
  const tl=LAST_L3.timeline;
  if(!tl.length) return 0;
  if(t<=tl[0].t) return tl[0].z;
  for(let i=1;i<tl.length;i++){
    if(tl[i].t>=t){
      const a=tl[i-1], b=tl[i];
      const f=(t-a.t)/((b.t-a.t)||1);
      return a.z+f*(b.z-a.z);
    }
  }
  return tl[tl.length-1].z;
}

function stageAtTime(t){
  const st=LAST_L3.stages||[];
  let name="静力核查";
  for(const s of st){ if(t>=s.t) name=s.name; }
  return name;
}

function drawL3Frame(){
  if(!LAST_L3) return;
  const frac=parseFloat(document.getElementById("l3slider").value)/1000;
  const t=frac*LAST_L3.T, z=zAtTime(t);
  const H=LAST_L3.params.h||3;
  const top=20, bottom=235, full=bottom-top;
  const yf=top+(z/H)*full;
  const w=92, x0=64;
  const svg=document.getElementById("l3svg");
  svg.innerHTML=
    '<rect x="'+x0+'" y="'+top+'" width="'+w+'" height="'+full+'" fill="#f3f5f7" stroke="#888"/>'+
    '<rect x="'+x0+'" y="'+top+'" width="'+w+'" height="'+Math.max(0,yf-top)+'" fill="#d9e4f0"/>'+
    '<line x1="'+(x0-8)+'" y1="'+yf+'" x2="'+(x0+w+8)+'" y2="'+yf+'" stroke="#c00" stroke-width="2.5"/>'+
    '<text x="6" y="'+(top+full+14)+'" font-size="10">z=0</text>'+
    '<text x="6" y="'+(top+10)+'" font-size="10">z='+fmt(H,1)+' m</text>'+
    '<text x="'+(x0+w+14)+'" y="'+(yf+4)+'" font-size="11" fill="#c00">'+fmt(z,2)+' m</text>';
  document.getElementById("l3frameinfo").innerHTML=
    "t = "+fmt(t,2)+" s / "+fmt(LAST_L3.T,2)+" s<br>"+
    "前锋深度 z = "+fmt(z,3)+" m<br>"+
    "当前阶段："+stageAtTime(t)+"<br>"+
    "v = "+fmt(LAST_L3.v,3)+" m/s";
}

function togglePlay(){
  if(!LAST_L3) return;
  if(PLAY_TIMER){ clearInterval(PLAY_TIMER); PLAY_TIMER=null; return; }
  PLAY_TIMER=setInterval(()=>{
    const e=document.getElementById("l3slider");
    let v=parseFloat(e.value)+6;
    if(v>1000){ v=0; }
    e.value=v; drawL3Frame();
  },40);
}

async function runMirror(){
  document.getElementById("err").textContent=""; hidePlayer();
  let params; try{ params=read(); }
  catch(e){ document.getElementById("err").textContent=e.message; return; }
  const r=await fetch("/api/mirror",{method:"POST",headers:{"Content-Type":"application/json"},
    body:JSON.stringify({params:params})});
  const d=await r.json();
  if(d.error){ document.getElementById("err").textContent=d.error; return; }
  const s=d.result;
  document.getElementById("banner").innerHTML=
    `<div class="banner ${s.mirror_ok?'ok':'warn'}"><b>镜像对照：</b>主动出露 x₂=${fmt(s.x2_active,4)} m，被动出露 x₂=${fmt(s.x2_passive,4)} m，两者之和 ${fmt(s.x2_sum,4)} m（2l=${fmt(s.two_l,4)} m）；逐点最大误差 ${s.max_mirror_error_m.toExponential(2)} m</div>`;
  document.getElementById("kv").innerHTML=[["洞半宽 l",fmt(s.l,3)+" m"],["覆土 h",fmt(s.h,3)+" m"],["镜像逐点最大误差",s.max_mirror_error_m.toExponential(2)+" m"]].map(x=>`<span>${x[0]}</span><span>${x[1]}</span>`).join("");
  document.getElementById("figs").innerHTML=`<figure><img src="data:image/png;base64,${d.fig}" alt="mirror"><figcaption>主动塌落与被动顶起滑面镜像对照（虚线为主动支关于洞中线 x=l 的镜像）</figcaption></figure>`;
  document.getElementById("extra").innerHTML=s.note;
}

async function runGTT(){
  document.getElementById("err").textContent=""; hidePlayer();
  let params; try{ params=read(); }
  catch(e){ document.getElementById("err").textContent=e.message; return; }
  const s0raw=document.getElementById("gtt_sigma0").value.trim();
  const payload={params:params, mode:document.getElementById("mode").value,
                 sigma0:s0raw===""?null:parseFloat(s0raw)};
  if(payload.sigma0!==null && !isFinite(payload.sigma0)){
    document.getElementById("err").textContent="GTT 切点不是有效数字"; return; }
  const r=await fetch("/api/gtt",{method:"POST",headers:{"Content-Type":"application/json"},
    body:JSON.stringify(payload)});
  const d=await r.json();
  if(d.error){ document.getElementById("err").textContent=d.error; return; }
  const s=d.result, t=s.tangent, m=s.mechanism, chk=s.table_check;
  const round=chk.rounding_rows.length?`（表 13.2 第 ${chk.rounding_rows.join(",")} kPa 行为原表印刷舍入差，已在结果中标出）`:"";
  document.getElementById("banner").innerHTML=
    `<div class="banner ${chk.all_pass?'ok':'warn'}"><b>GTT 切线：</b>tanδ=${fmt(t.tan_delta,4)}，δ=ψ=${fmt(t.delta_deg,2)}°，c_t=${fmt(t.c_t,3)} kPa；表 13.2 ${chk.all_pass?'逐项通过':'存在偏差'} ${round}</div>`+
    `<div class="banner ${m.refined&&m.theorem13_2_rel_err<1e-5?'ok':'warn'}"><b>定理 13.2：</b>耗散上界最小点 σ₀=${fmt(m.sigma0_opt,4)} kPa，滑面正应力加权均值 σ̄=${fmt(m.sigma_bar,4)} kPa，相对误差 ${m.theorem13_2_rel_err.toExponential(2)}</div>`;
  document.getElementById("kv").innerHTML=[
    ["切点 σ₀",fmt(t.sigma_0,3)+" kPa"],["滑面倾角 90°−ψ",fmt(90-t.delta_deg,2)+"°"],
    ["切点剪应力 τ",fmt(t.tau_at_point,3)+" kPa"],["切线构造残差",t.selfcheck_residual.toExponential(2)],
    ["D 真值",fmt(m.D_true,3)],["D 上界",fmt(m.D_upper,3)],["D 下界",fmt(m.D_lower,3)],
    ["上界余量",fmt(m.upper_residual,3)],["下界余量",fmt(m.lower_residual,4)]
  ].map(x=>`<span>${x[0]}</span><span>${x[1]}</span>`).join("");
  document.getElementById("figs").innerHTML=`<figure><img src="data:image/png;base64,${d.fig}" alt="GTT"><figcaption>GTT 幂型包线、切线上界、弦线下界与定理由线扫描</figcaption></figure>`;
  document.getElementById("extra").innerHTML=
    `<b>表 13.2 校核：</b>${chk.results.map(x=>`B=${x.B},σ₀=${x.sigma_0}:${x.strict_ok?'精确':'舍入'}✓`).join("　")}`+
    `<br><b>边界：</b>${m.boundary_note}。定理 13.1 切点倾角等于剪胀角；B=1 退化为 Mohr–Coulomb。`;
}

async function runSat(){
  document.getElementById("err").textContent=""; hidePlayer();
  let params; try{ params=read(); }
  catch(e){ document.getElementById("err").textContent=e.message; return; }
  const r=await fetch("/api/sat",{method:"POST",headers:{"Content-Type":"application/json"},
    body:JSON.stringify({params:params})});
  const d=await r.json();
  if(d.error){ document.getElementById("err").textContent=d.error; return; }
  const s=d.result;
  document.getElementById("banner").innerHTML=
    `<div class="banner ${s.consistency_ok?'ok':'warn'}"><b>文献锚点：</b>N_sat=${fmt(s.N_sat,4)}，h_sat/D=${fmt(s.h_sat_over_D_printed,3)}（=${fmt(s.h_sat_over_D_exact,4)}，相对差 ${s.consistency_rel_err.toExponential(2)}）</div>`+
    `<div class="banner warn"><b>门槛区分：</b>h_cr/D=${fmt(s.h_cr_over_D,3)}（几何贯通）≠ h_sat/D=${fmt(s.h_sat_over_D_exact,3)}（力学饱和）；${s.threshold_note}</div>`;
  document.getElementById("kv").innerHTML=[
    ["h_cr (m)",fmt(s.h_cr_m,3)],["h_cr/D",fmt(s.h_cr_over_D,3)],
    ["h_sat/D (文献)",fmt(s.h_sat_over_D_exact,4)],["常数自洽",s.consistency_ok?"是":"否"]
  ].map(x=>`<span>${x[0]}</span><span>${x[1]}</span>`).join("");
  document.getElementById("figs").innerHTML="";
  document.getElementById("extra").innerHTML=
    `<table><tr><th>h/D</th><th>h (m)</th><th>L1 q1a (kPa)</th><th>整体贯通</th><th>L3 适用</th><th>说明</th></tr>`+
    s.rows.map(x=>`<tr><td>${fmt(x.h_over_D,3)}</td><td>${fmt(x.h_m,3)}</td><td>${fmt(x.q1a_kPa,2)}</td><td>${x.globally_through?'是':'否'}</td><td>${x.l3_applicable?'是':'否'}</td><td>${x.note}</td></tr>`).join("")+
    `</table><br><b>边界：</b>${s.boundary_note}。${s.source}`;
}

async function runBench2(){
  document.getElementById("err").textContent=""; hidePlayer();
  const r=await fetch("/api/bench2");
  const d=await r.json();
  if(d.error){ document.getElementById("err").textContent=d.error; return; }
  const c=d.check, a=d.anchors;
  document.getElementById("banner").innerHTML=
    `<div class="banner ${c.all_pass?'ok':'warn'}"><b>二期文献基准集：</b>自洽核验 ${c.all_pass?'通过':'未过'}；`+
    `1g ${c.counts.one_g} 项 / 离心机 ${c.counts.centrifuge} 项 / 速率适用域 ${c.counts.rate_applicability} 项</div>`+
    `<div class="banner warn"><b>定位：</b>${c.note}。全部数值逐条摘自专著，不做拟合、不推导。</div>`;
  document.getElementById("kv").innerHTML=Object.keys(c.checks)
    .map(k=>`<span>${k}</span><span>${c.checks[k]?'✓':'✗'}</span>`).join("");
  document.getElementById("figs").innerHTML="";
  const grp=(title,rows)=>`<h3>${title}</h3><table><tr><th>锚点</th><th>数值</th><th>单位</th><th>含义 / 效应</th><th>来源</th><th>专著</th></tr>`+
    rows.map(x=>{ const v=(x.value_min!==undefined)?`${x.value_min}–${x.value_max}`:x.value;
      return `<tr><td>${x.label}</td><td>${v}</td><td>${x.unit}</td><td>${x.effect}</td><td>${x.source}</td><td>${x.ref}</td></tr>`;}).join("")+`</table>`;
  document.getElementById("extra").innerHTML=
    grp("1g 活门试验",a.one_g)+"<br>"+grp("离心机试验",a.centrifuge)+"<br>"+grp("速率层适用域 / 极限分析常数",a.rate_applicability);
}

async function exportReport(){
  if(!LAST){ return; }
  let tdata=null;
  try{ const rr=await fetch("/api/tables"); tdata=await rr.json(); }catch(e){}
  const w=window.open("","_blank");
  const res=LAST.result;
  const prov=LAST.provenance?"<h2>可复现性元数据</h2><pre>"+JSON.stringify(LAST.provenance,null,2)+"</pre>":"";
  let tbl="";
  if(tdata && tdata.table3 && tdata.table4){
    const t3=tdata.table3, t4=tdata.table4;
    tbl="<h2>试验对照（论文 Table 3 / Table 4）</h2>"+
      "<p>Table 3 主动：实测 ρa="+t3.meas_Ga+"；关联 "+fmt(t3.assoc_Ga,3)+" ("+fmt(t3.assoc_dev_pct,1)+"%)，非关联 "+fmt(t3.non_Ga,2)+" ("+fmt(t3.non_dev_pct,1)+"%)</p>"+
      "<table><tr><th>n</th><th>φp</th><th>θ</th><th>实测</th><th>关联</th><th>非关联</th></tr>"+
      t4.rows.map(x=>"<tr><td>"+x.n+"</td><td>"+x.phi_peak_deg+"</td><td>"+x.theta_deg+"</td><td>"+x.meas+"</td><td>"+fmt(x.assoc,2)+"</td><td>"+fmt(x.non,2)+"</td></tr>").join("")+
      "</table><p>Table 4 平均偏差：关联 "+fmt(t4.assoc_mean_dev_pct,1)+"%，非关联 "+fmt(t4.non_mean_dev_pct,1)+"%。</p>";
  }
  w.document.write("<meta charset=\"utf-8\"><h1>活门试验数值模拟计算报告</h1>"+
    "<p>生成时间："+new Date().toLocaleString()+"</p>"+
    "<p>工况："+(res.mode==="active"?"主动塌落":"被动顶起")+"　求解器："+res.solution_type+"</p>"+
    "<h2>输入参数</h2><pre>"+JSON.stringify(res.params,null,2)+"</pre>"+
    "<h2>结果</h2><pre>"+JSON.stringify(res,null,2)+"</pre>"+prov+tbl+
    "<h2>图件</h2>"+Object.keys(LAST.figures).map(k=>"<p><img width=\"640\" src=\"data:image/png;base64,"+LAST.figures[k]+"\"></p>").join(""));
  w.document.close();
}
</script>
</body></html>
"""


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt_, *args):
        sys.stderr.write("[web] " + (fmt_ % args) + "\n")

    def _send(self, code: int, body: bytes, ctype: str):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        path = urlparse(self.path).path
        if path in ("/", "/index.html"):
            html = HTML.replace("__PRESETS__", json.dumps(PRESETS, ensure_ascii=False)) \
                       .replace("__DEFAULT__", json.dumps(DEFAULT_PARAMS, ensure_ascii=False))
            self._send(200, html.encode("utf-8"), "text/html; charset=utf-8")
        elif path == "/api/tables":
            self._send(200, json.dumps(reproduce_tables(), ensure_ascii=False).encode("utf-8"),
                       "application/json; charset=utf-8")
        elif path == "/api/health":
            self._send(200, b'{"ok":true}', "application/json")
        elif path == "/api/bench2":
            self._json({"anchors": LITERATURE_ANCHORS,
                        "check": check_literature_anchors()})
        else:
            self._send(404, b"not found", "text/plain")

    def _json(self, obj):
        self._send(200, json.dumps(obj, ensure_ascii=False).encode("utf-8"),
                   "application/json; charset=utf-8")

    def do_POST(self):
        path = urlparse(self.path).path
        if path not in ("/api/solve", "/api/l2", "/api/l3", "/api/gtt",
                        "/api/sat", "/api/mirror"):
            self._send(404, b"not found", "text/plain")
            return
        try:
            n = int(self.headers.get("Content-Length", "0"))
            payload = json.loads(self.rfile.read(n).decode("utf-8"))
            if path == "/api/mirror":
                fdata = figure_data(payload["params"], "active")
                s = fdata["fig1_slip_surface"]
                res = self._mirror_check(s)
                fig = self._plot_mirror(s, payload["params"])
                self._json({"result": res, "fig": fig})
            elif path == "/api/sat":
                self._json({"result": saturation_anchor(payload["params"])})
            elif path == "/api/solve":
                params, mode = payload["params"], payload.get("mode", "active")
                res = solve(params, mode)
                fdata = figure_data(params, mode)
                figures = render_all(fdata, res.Omega, res.psi0_deg)
                self._json({"result": res.to_dict(), "figures": figures,
                            "provenance": provenance(params, mode)})
            elif path == "/api/l2":
                params, mode = payload["params"], payload.get("mode", "active")
                r = solve_shooting(params, mode)
                fig = self._plot_l2(params, mode, r)
                self._json({"result": r, "fig": fig})
            elif path == "/api/l3":
                rate = solve_rate(payload["params"], payload.get("phi_deg"))
                fig = self._plot_l3(rate)
                self._json({"result": rate, "fig": fig})
            else:
                params = payload["params"]
                mode = payload.get("mode", "active")
                s0 = payload.get("sigma0")
                mech = mechanism_scan(params, mode)
                if s0 is None:
                    s0 = mech["sigma_bar"]
                s0 = float(s0)
                tan = tangent_params(s0, params["sigma_c"], params["A"],
                                     params["sigma_t"], params["B"])
                nodes = np.linspace(mech["sigma_min"], mech["sigma_max"], 65)
                br = dissipation_bracket(
                    np.linspace(mech["sigma_min"], mech["sigma_max"], 160),
                    params["sigma_c"], params["A"], params["sigma_t"], params["B"],
                    s0, nodes)
                res = dict(tangent=tan, table=table_13_2(),
                           table_check=check_table_13_2(),
                           bracket=dict(upper_residual=br["upper_residual"],
                                        lower_residual=br["lower_residual"],
                                        max_gap=br["max_gap"]),
                           mechanism=mech)
                fig = self._plot_gtt(params, s0, br, mech)
                self._json({"result": res, "fig": fig})
        except (ParamsError, KeyError, ValueError) as e:
            self._json({"error": str(e)})
        except Exception as e:  # noqa
            self._json({"error": "内部错误: %s" % e})

    @staticmethod
    def _plot_l2(params, mode, r):
        from collapse_kernel import slip_profile
        fig, ax = plt.subplots(figsize=(6.2, 4.5))
        ax.plot(r["x"], r["F"], "r-", lw=2.2, label="L2 shooting (RK4+bisection)")
        xc, Fc = slip_profile(params, mode)
        ax.plot(xc, Fc, "b--", lw=1.6, label="L1 closed form")
        ax.axhline(params["h"], color="k", lw=1); ax.axhline(0, color="k", lw=1)
        ax.plot([-params["l"], params["l"]], [0, 0], "r-", lw=3, alpha=0.5)
        ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")
        ax.set_title("L2 slip surface vs L1 closed form")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        fig.tight_layout()
        return _fig_b64(fig)

    @staticmethod
    def _mirror_check(s):
        """主动滑面关于洞中线 x=l 镜像后，与被动滑面逐点比对。"""
        xa = np.asarray(s["active_x"], dtype=float)
        Fa = np.asarray(s["active_F"], dtype=float)
        xp = np.asarray(s["passive_x"], dtype=float)
        Fp = np.asarray(s["passive_F"], dtype=float)
        xm = 2.0 * s["l"] - xa[::-1]
        Fm = Fa[::-1]
        # 被动支 x 从 x2 递减到 l；np.interp 要求横坐标升序，必须先排序
        order = np.argsort(xp)
        xp_sorted, Fp_sorted = xp[order], Fp[order]
        common_lo = max(xm.min(), xp_sorted.min())
        common_hi = min(xm.max(), xp_sorted.max())
        xs = np.linspace(common_lo, common_hi, 200)
        Fm_i = np.interp(xs, xm, Fm)
        Fp_i = np.interp(xs, xp_sorted, Fp_sorted)
        maxerr = float(np.max(np.abs(Fm_i - Fp_i)))
        clipped = bool(xa.min() <= 1e-12)
        note = ("主动滑面关于 x=l 镜像后与被动滑面重合；"
                "残差来自 400 点网格线性插值，非求解分歧")
        if clipped:
            note += ("；注意：主动出露被裁剪到 x=0（几何边界，dx>l），"
                     "此时与被动支不严格镜像，属内核几何约定")
        return dict(l=s["l"], h=s["h"], clipped=clipped,
                    x2_active=float(xa.min()), x2_passive=float(xp.max()),
                    x2_sum=float(xa.min() + xp.max()), two_l=2.0 * s["l"],
                    max_mirror_error_m=maxerr,
                    # 闭式图仅 400 点网格，线性插值带来 O(dx^2) 误差（约 1e-5 m），
                    # 与 L2 曲线交叉校验同源，不是求解分歧。
                    mirror_ok=(maxerr < 1e-5 and not clipped),
                    note=note)

    @staticmethod
    def _plot_mirror(s, params):
        xa = np.asarray(s["active_x"], dtype=float)
        Fa = np.asarray(s["active_F"], dtype=float)
        xp = np.asarray(s["passive_x"], dtype=float)
        Fp = np.asarray(s["passive_F"], dtype=float)
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        ax.plot(xa, Fa, "b-", lw=2.0, label="active (collapse)")
        ax.plot(xp, Fp, "r-", lw=2.0, label="passive (uplift)")
        ax.plot(2.0 * s["l"] - xa[::-1], Fa[::-1], "k--", lw=1.3,
                label="active mirrored about x=l")
        ax.axvline(s["l"], color="k", ls=":", lw=1.2, label="cavity centre x=l")
        ax.axhline(s["h"], color="k", lw=1)
        ax.axhline(0, color="k", lw=1)
        ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")
        ax.set_title("Active / passive mirror check")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)
        fig.tight_layout()
        return _fig_b64(fig)

    @staticmethod
    def _plot_gtt(params, s0, br, mech):
        sc, A, st, B = params["sigma_c"], params["A"], params["sigma_t"], params["B"]
        sg = np.asarray(br["sigma"], dtype=float)
        true = np.asarray(tau_envelope(sg, sc, A, st, B))
        up = np.asarray(br["tau_upper"])
        low = np.asarray(br["tau_lower"])
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        ax1.plot(sg, true, "k-", lw=2.2, label="power envelope (13.9)")
        ax1.plot(sg, up, "r--", lw=1.6, label="GTT tangent (upper)")
        ax1.plot(sg, low, "b-.", lw=1.6, label="chord (lower)")
        ax1.plot(s0, float(tau_envelope(s0, sc, A, st, B)), "ro", ms=5)
        ax1.set_xlabel("sigma (kPa)"); ax1.set_ylabel("tau (kPa)")
        ax1.set_title("Strength envelope bracket"); ax1.legend(fontsize=7); ax1.grid(alpha=0.3)
        ax2.plot(mech["scan_sigma"], mech["scan_D_upper"], "b-", lw=1.8,
                 label="D_UB(sigma_0)")
        ax2.axvline(mech["sigma_bar"], color="k", ls=":", lw=1.4, label="sigma_bar")
        ax2.plot([mech["sigma0_opt"]], [mech["D_upper"]], "ro", ms=5,
                 label="optimum")
        ax2.axhline(mech["D_true"], color="g", ls="--", lw=1.4, label="D_true")
        ax2.set_xlabel("sigma_0 (kPa)"); ax2.set_ylabel("D (force/length)")
        ax2.set_title("Theorem 13.2: dissipation scan"); ax2.legend(fontsize=7)
        ax2.grid(alpha=0.3)
        fig.tight_layout()
        return _fig_b64(fig)

    @staticmethod
    def _plot_l3(rate):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
        ts = [x["t"] for x in rate["timeline"]]
        zs = [x["z"] for x in rate["timeline"]]
        ax1.plot(ts, zs, "b-", lw=2)
        ax1.set_xlabel("t (s)"); ax1.set_ylabel("chimney depth z (m)")
        ax1.set_title("Collapse timeline"); ax1.grid(alpha=0.3)
        labels = ["static", "takeover", "advance", "through", "funnel"]
        tt = [0, 0, rate["T"] * 0.5, rate["T"], rate["T"]]
        vv = [0, 0, rate["v"], rate["v"], 0]
        ax2.plot(tt, vv, "o-", lw=2)
        ax2.set_xticks(range(5)); ax2.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
        ax2.set_ylabel("v (m/s)"); ax2.set_title("Five stages"); ax2.grid(alpha=0.3)
        fig.tight_layout()
        return _fig_b64(fig)


def _bind_server():
    """端口被上次残留的实例占用时，自动顺延到下一个可用端口。

    以前直接绑定固定 8765，一旦上次窗口没关干净就抛
    OSError: [Errno 98] Address already in use，脚本窗口一闪而过，
    用户看到的就是「启动不了」。现在最多顺延 20 个端口，并在控制台
    明确打印实际地址，浏览器也打开实际地址。
    """
    last_error = None
    for port in range(PORT, PORT + 20):
        try:
            return ThreadingHTTPServer(("127.0.0.1", port), Handler), port
        except OSError as e:
            last_error = e
            print("端口 %d 已被占用，尝试 %d ..." % (port, port + 1), flush=True)
    raise RuntimeError(
        "端口 %d–%d 全部不可用，请关闭已打开的程序窗口后重试。\n原始错误: %s"
        % (PORT, PORT + 19, last_error)
    )


def main():
    try:
        srv, actual_port = _bind_server()
    except RuntimeError as e:
        print("=" * 60, flush=True)
        print("启动失败：%s" % e, flush=True)
        print("=" * 60, flush=True)
        raise
    url = "http://127.0.0.1:%d" % actual_port
    # flush=True 必要：PyInstaller 冻结程序里 stdout 为块缓冲，
    # 不强制刷新时用户在终端看不到实际端口（端口顺延时尤其致命）。
    print("=" * 60, flush=True)
    print(" 浅埋土洞活门试验数值模拟程序 已启动", flush=True)
    print(" 浏览器地址: " + url, flush=True)
    if actual_port != PORT:
        print(" 注意：默认端口 %d 被占用，已自动改用 %d" % (PORT, actual_port), flush=True)
    print(" 关闭程序: 在本窗口按 Ctrl+C", flush=True)
    print("=" * 60, flush=True)
    if os.environ.get("TRAPDOOR_NO_BROWSER") != "1":
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\n已关闭。")
        srv.shutdown()


if __name__ == "__main__":
    main()
