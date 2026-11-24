# -*- coding: utf-8 -*-
"""
浅埋土洞坍塌机理 —— 统一变分极限分析 一键计算与绘图小程序
=================================================================
本程序是论文
  "A unified variational limit-analysis theory for active collapse
   and passive uplift of shallow soil cavities"  (Haipeng Sun)
配套的可复现计算脚本。

功能：
  1. 打印核心闭式解结果（强度组合参数 S、临界覆盖厚度、剪胀角、
     主动/被动极限荷载、镜面对称校验）；
  2. 自动生成 6 张图：
     fig1 滑动面形态 (主动 urn / 被动 funnel)
     fig2 滑动面应力场 (sigma, tau 随深度)
     fig3 剪胀角 psi 随深度衰减
     fig4 四形态分区图 (Fig.5, Omega-psi0 平面)
     fig5 侧压力系数 Ka(z) 随深度演化 (Fig.6, Eq.44)
     fig6 非关联修正 A_psi/A - psi0 (Fig.7a, Eq.46)
  3. 自动复现 Table 3 (主动试验 Chen 2019) 与 Table 4 (被动试验 Rui 2025b)
     并打印预测值与偏差，对照论文。

小白用法：
  1) 只改下面 PARAMS 里的参数 (单位见注释)，其余一行都不用碰；
  2) 终端运行:  python3 cavity_collapse_calculator.py
  3) 当前目录会生成 fig1~fig6 共 6 张 png。

图中文字标签用英文，避免 Linux 缺中文字体出现方块；代码注释用中文方便你看。
"""

import os
import numpy as np

import matplotlib
matplotlib.use("Agg")            # 无窗口环境下也能出图
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.linewidth": 1.0,
                     "figure.dpi": 150, "savefig.dpi": 150})

# 图片输出目录：默认与脚本同目录。可用环境变量 OUTDIR 临时改变。
OUTDIR = os.environ.get("OUTDIR", os.path.dirname(os.path.abspath(__file__)))
os.makedirs(OUTDIR, exist_ok=True)

# ===========================================================================
# 1. 输入参数 (PARAMS) —— 论文 Fig.4 基准算例。只改这一块。
# ===========================================================================
PARAMS = dict(
    sigma_c = 100.0,   # sigma_c  强度参考应力 (kPa)
    A       = 0.6,     # A        幂型强度准则系数 (tan 量级)
    B       = 0.6,     # B        幂型强度准则指数 (0<B<=1, B=1 退化为 Mohr-Coulomb)
    sigma_t = 5.0,     # sigma_t  抗拉强度 (kPa)
    gamma   = 18.0,    # gamma    土体重度 (kN/m^3)
    q2      = 10.0,    # q2       地表均布超载 (kPa)
    l       = 2.5,     # l        洞体半宽 (m)
    h       = 3.0,     # h        覆盖层厚度 (m)
)

# 分区图 (Fig.5) 用论文的代表参数
ZONING = dict(A=0.6, B=0.5)
# 分区竖直虚线 (论文给出的代表值)
ZONE_OMEGA = [0.0196, 0.184, 0.963]


# ===========================================================================
# 2. 核心公式 (每个函数注释对应论文公式编号)
# ===========================================================================
def strength_S(sigma_c, A, B):
    """强度组合参数 S = sigma_c * A^(1/(1-B)) * B^(B/(1-B))   (见 Eq.6/7 下方)"""
    if abs(B - 1.0) < 1e-12:                 # B->1 线性极限单独处理
        return None
    return sigma_c * A**(1.0/(1.0-B)) * B**(B/(1.0-B))


def t_param(sigma_c, sigma_t, q2, A, B):
    """出露坡度参数 t = (sigma_t+q2) / (B*S) = (cot psi0)^(1/(1-B))   (Eq.21)"""
    S = strength_S(sigma_c, A, B)
    return (sigma_t + q2) / (B * S)


def sigma_of_z(z, q2, gamma):
    """滑动面法向应力 sigma(z) = q2 + gamma*z   (Eq.28)"""
    return q2 + gamma*np.asarray(z, dtype=float)


def tau_of_z(z, sigma_c, sigma_t, q2, gamma, A, B):
    """滑动面剪应力 tau(z) = A*sigma_c*((sigma_t+q2+gamma*z)/sigma_c)^B   (Eq.29)"""
    return A*sigma_c*((sigma_t + q2 + gamma*np.asarray(z, dtype=float))/sigma_c)**B


def psi_of_z(z, sigma_c, sigma_t, q2, gamma, A, B):
    """剪胀角 psi(z) = arctan( A*B*((sigma_t+q2+gamma*z)/sigma_c)^(B-1) )   (Eq.30)"""
    val = A*B*((sigma_t + q2 + gamma*np.asarray(z, dtype=float))/sigma_c)**(B-1.0)
    return np.arctan(val)


def R_scale(sigma_c, sigma_t, q2, gamma, h, A, B):
    """剪应力深度积分量 R_a = A*sigma_c^(1-B)*[(sig_t+q2+g h)^(1+B)-(sig_t+q2)^(1+B)]/((1+B)g)  (Eq.36)"""
    top = (sigma_t + q2 + gamma*h)**(1.0+B)
    bot = (sigma_t + q2)**(1.0+B)
    return A*sigma_c**(1.0-B)*(top-bot)/((1.0+B)*gamma)


def limit_loads(sigma_c, sigma_t, q2, gamma, l, h, A, B):
    """主动/被动极限荷载 q1a, q1b   (Eq.35a/35b)"""
    Ra = R_scale(sigma_c, sigma_t, q2, gamma, h, A, B)
    over = q2 + gamma*h
    q1a = over - Ra/l
    q1b = over + Ra/l
    return q1a, q1b, Ra, over


def h_critical(sigma_c, sigma_t, q2, gamma, l, A, B):
    """临界覆盖厚度 h_cr   (Eq.32)"""
    S = strength_S(sigma_c, A, B)
    t = t_param(sigma_c, sigma_t, q2, A, B)
    return (B*S/gamma)*((gamma*l/S + t**B)**(1.0/B) - t)


def psi0_of_Omega(Omega, A, B):
    """地表剪胀角 psi0 = arctan(A*B*Omega^(B-1))，Omega=(sig_t+q2)/sig_c   (Eq.38)"""
    return np.arctan(A*B*np.asarray(Omega, dtype=float)**(B-1.0))


def Ka_local(sigma, tau, psi):
    """侧压力系数 Ka=(sig cos psi - tau(1-sin psi))/(sig cos psi + tau(1+sin psi))   (Eq.44a)"""
    s = np.sin(psi); c = np.cos(psi)
    return (sigma*c - tau*(1.0-s)) / (sigma*c + tau*(1.0+s))


def A_psi_from_psi0(psi0, Omega, B):
    """非关联势斜率参数 A_psi = tan(psi0) / (B*Omega^(B-1))   (Eq.46)"""
    return np.tan(psi0) / (B*np.asarray(Omega, dtype=float)**(B-1.0))


def S_psi(sigma_c, A, B, A_psi):
    """非关联修正参数 S_psi   (Eq.47, 照论文形式: sigma_c*(A_psi*B)^(B/(1-B))*(A-A_psi*B)^(1-B))"""
    return sigma_c*(A_psi*B)**(B/(1.0-B))*(A - A_psi*B)**(1.0-B)


# 线性 Mohr-Coulomb 退化 (用于试验对比 Table 3/4, Eq.42/43/48)
def h_cr_linear(l, phi):
    """h_cr = l*cot(phi)   (Eq.42, 经典自然拱高)"""
    return l/np.tan(phi)

def Ga_linear(eta, phi):
    """主动土压力比 G_a = 1 - eta*tan(phi)   (Eq.48a)"""
    return 1.0 - eta*np.tan(phi)

def Gb_linear(eta, phi):
    """被动土压力比 G_b = 1 + eta*tan(phi)   (Eq.48b)"""
    return 1.0 + eta*np.tan(phi)


# ===========================================================================
# 3. 滑动面坐标求解 (主动 urn / 被动 funnel)，用于 fig1
# ===========================================================================
def slip_profile(p, mode="active", n=400):
    """返回滑动面 (x, F) 列表。F(0)=h(地表), F(l)=0(洞缘)。
       深度形式：z = (B*S/gamma)*[(t^B + gamma*varpi/S)^(1/B) - t]  (Eq.26)
       varpi 为从出露点起的水平投影，主动 x=x2-varpi 向内、被动 x=l+varpi 向外。"""
    sc, A, B, st = p["sigma_c"], p["A"], p["B"], p["sigma_t"]
    q2, g, l, h = p["q2"], p["gamma"], p["l"], p["h"]
    S = strength_S(sc, A, B); t = t_param(sc, st, q2, A, B)

    # 把 z 离散到 [0,h]，反解 varpi，再映射到 x
    z = np.linspace(0.0, h, n)
    base = t + g*z/(B*S)                      # [t^B+g*varpi/S]^(1/B)  (Eq.26 反解)
    varpi = S/g*(base**B - t**B)              # Eq.23/26
    # 出露点水平位置：主动 x2 = l - delta_a(h)，delta_a 取 z=h 处 varpi
    if mode == "active":
        delta_a = S/g*(base[-1]**B - t**B)    # = l - x2  (Eq.31)
        x2 = l - delta_a
        x = x2 + varpi                        # 出露点 x2 -> 洞缘 l
    else:
        x2 = l + S/g*(base[-1]**B - t**B)     # = l + Delta  (Eq.33)
        x = x2 - varpi                        # 出露点 x2 -> 洞缘 l
    F = h - z
    return np.clip(x, 0, None), F


# ===========================================================================
# 4. 主流程
# ===========================================================================
def main():
    p = PARAMS
    sc, A, B, st = p["sigma_c"], p["A"], p["B"], p["sigma_t"]
    q2, g, l, h = p["q2"], p["gamma"], p["l"], p["h"]

    # ---------- 控制台结果 ----------
    S = strength_S(sc, A, B); t = t_param(sc, st, q2, A, B)
    Omega = (st+q2)/sc; Lambda = g*l/sc; eta = h/(2*l)
    q1a, q1b, Ra, over = limit_loads(sc, st, q2, g, l, h, A, B)
    h_cr = h_critical(sc, st, q2, g, l, A, B)
    psi_surf = np.degrees(psi_of_z(0.0, sc, st, q2, g, A, B))
    psi_bot  = np.degrees(psi_of_z(h,   sc, st, q2, g, A, B))

    print("="*66)
    print(" 浅埋土洞坍塌 变分极限分析  —— 核心闭式解结果")
    print("="*66)
    print(f" 强度组合参数   S        = {S:8.3f} kPa")
    print(f" 出露坡度参数   t        = {t:8.3f} (=(cot psi0)^(1/(1-B)))")
    print(f" 脆性比         Omega    = {Omega:8.4f}")
    print(f" 覆强比         Lambda   = {Lambda:8.4f}")
    print(f" 埋深比         eta      = {eta:8.4f}")
    print("-"*66)
    print(f" 地表剪胀角     psi0     = {psi_surf:8.2f} deg   (论文 37.6)")
    print(f" 底部剪胀角     psi(h)   = {psi_bot:8.2f} deg   (论文 22.7)")
    print(f" 临界覆盖厚度   h_cr     = {h_cr:8.3f} m")
    if h <= h_cr:
        print(f"                -> h(={h} m) <= h_cr：浅埋，滑动面贯通地表")
    else:
        print(f"                -> h(={h} m) >  h_cr：深埋，坍塌不贯通地表")
    print("-"*66)
    print(f" 上覆压力       q2+g.h   = {over:8.3f} kPa")
    print(f" 主动支护压力   q1a      = {q1a:8.3f} kPa   (< 上覆, arching)")
    print(f" 被动顶起压力   q1b      = {q1b:8.3f} kPa   (> 上覆)")
    print(f" 镜面对称校验   q1b-q1a  = {q1b-q1a:8.3f}")
    print(f"                 2*Ra/l   = {2*Ra/l:8.3f}   (应相等, 误差 {abs((q1b-q1a)-2*Ra/l):.2e})")
    print("="*66)

    # ---------- fig1 滑动面形态 ----------
    xa, Fa = slip_profile(p, "active")
    xp, Fp = slip_profile(p, "passive")
    fig, ax = plt.subplots(figsize=(7, 5))
    for sx, sF, lab, col in [(xa, Fa, "Active (urn, inward)", "#1f77b4"),
                             (xp, Fp, "Passive (funnel, outward)", "#2ca02c")]:
        # 右半结构 + 镜像成完整剖面
        ax.plot(sx,  sF, color=col, lw=2, label=lab)
        ax.plot(-sx, sF, color=col, lw=2)
    # 地表与洞体
    ax.axhline(h, color="k", lw=1)
    ax.axhline(0, color="k", lw=1)
    ax.plot([-l, l], [0, 0], "r-", lw=4, label="cavity")
    ax.annotate("ground surface z=h", (l*1.05, h), fontsize=9)
    ax.annotate("cavity roof/floor z=0", (l*1.05, 0), fontsize=9)
    ax.set_xlabel("x  (m)"); ax.set_ylabel("z  (m)")
    ax.set_title("Fig.1  Slip-surface morphology")
    ax.legend(loc="upper left"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "fig1_slip_surface.png")); plt.close(fig)

    # ---------- fig2 应力场 ----------
    z = np.linspace(0, h, 200)
    sig = sigma_of_z(z, q2, g)
    tau = tau_of_z(z, sc, st, q2, g, A, B)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot(sig, z, "b-", lw=2, label="sigma = q2 + gamma*z  (Eq.28)")
    ax.plot(tau, z, "r-", lw=2, label="tau  (Eq.29)")
    ax.set_xlabel("stress (kPa)"); ax.set_ylabel("depth below surface (m)")
    ax.invert_yaxis(); ax.set_title("Fig.2  Stress field on the slip surface")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "fig2_stress_field.png")); plt.close(fig)

    # ---------- fig3 剪胀角衰减 ----------
    psi_pow = np.degrees(psi_of_z(z, sc, st, q2, g, A, B))
    psi_lin = np.degrees(np.arctan(A))*np.ones_like(z)   # B=1 时为常数
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.plot(psi_pow, z, "b-", lw=2, label="power-type (B<1), psi decays")
    ax.plot(psi_lin, z, "k--", lw=1.5, label="linear (B=1), psi constant")
    ax.set_xlabel("dilatancy angle psi (deg)"); ax.set_ylabel("depth below surface (m)")
    ax.invert_yaxis(); ax.set_title("Fig.3  Dilatancy angle vs depth (Eq.30)")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "fig3_dilatancy.png")); plt.close(fig)

    # ---------- fig4 四形态分区 (Fig.5) ----------
    AZ, BZ = ZONING["A"], ZONING["B"]
    Om = np.logspace(np.log10(0.005), np.log10(2.0), 400)
    psi0 = np.degrees(psi0_of_Omega(Om, AZ, BZ))
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.semilogx(Om, psi0, "b-", lw=2)
    for om, lab in zip(ZONE_OMEGA, ["I | II", "II | III", "III | IV"]):
        ax.axvline(om, color="grey", ls="--", lw=1)
    # zone labels
    ax.axvspan(0.005, ZONE_OMEGA[0], color="0.95")
    ax.axvspan(ZONE_OMEGA[0], ZONE_OMEGA[1], color="#ffe08a")  # Zone II 关键带
    ax.axvspan(ZONE_OMEGA[1], ZONE_OMEGA[2], color="0.95")
    ax.axvspan(ZONE_OMEGA[2], 2.0, color="0.9")
    ax.text(0.008, 60, "I\nhard rock", fontsize=9, ha="center")
    ax.text(0.05, 40, "II\nstiff clay\n(critical)", fontsize=9, ha="center", color="k")
    ax.text(0.4, 22, "III\nmedium soils", fontsize=9, ha="center")
    ax.text(1.3, 9, "IV\nsoft soil / shaft", fontsize=9, ha="center")
    ax.set_xlabel("Brittleness ratio  Omega = (sigma_t+q2)/sigma_c")
    ax.set_ylabel("Outcrop dilatancy angle psi0 (deg)")
    ax.set_title("Fig.4  Morphology zoning map (Eq.38, A=0.6, B=0.5)")
    ax.set_xlim(0.005, 2.0); ax.set_ylim(0, 85); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "fig4_zoning.png")); plt.close(fig)

    # ---------- fig5 侧压力系数 Ka(z) (Fig.6) ----------
    zz = np.linspace(0, 2*h, 300)              # 超出 h 部分为外推
    sg = sigma_of_z(zz, q2, g)
    # power-type (关联)
    tp_ = tau_of_z(zz, sc, st, q2, g, A, B)
    pp_ = psi_of_z(zz, sc, st, q2, g, A, B)
    Ka_pow = Ka_local(sg, tp_, pp_)
    # linear + cohesion (Eq.45)
    phi = np.arctan(A); sinf = np.sin(phi)
    Ka_linc = (1-sinf)*(sg - st*sinf)/((1+sinf)*(sg + st*sinf))
    # Rankine (B=1, sigma_t=0)
    Ka_rank = (1-sinf)/(1+sinf)*np.ones_like(zz)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.fill_betweenx(zz, 0.27, 0.33, color="0.85",
                     label="back-calculated band (loose sand, 0.30+-0.03)")
    ax.plot(Ka_pow, zz, "b-", lw=2, label="power-type (B=0.6, sig_t=5)")
    ax.plot(Ka_linc, zz, "-.", c="g", lw=1.8, label="linear + cohesion (B=1, sig_t=5)")
    ax.plot(Ka_rank, zz, "--", c="k", lw=1.5, label="Rankine constant (B=1, sig_t=0)")
    ax.axhline(h, color="r", lw=1, ls=":", label="ground depth h=3 m")
    ax.set_xlabel("lateral pressure coefficient  Ka(z)")
    ax.set_ylabel("depth below surface (m)"); ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_title("Fig.5  Lateral pressure coefficient vs depth (Eq.44)")
    ax.legend(fontsize=8, loc="lower right"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "fig5_Ka_depth.png")); plt.close(fig)

    # ---------- fig6 非关联修正 A_psi/A - psi0 (Fig.7a) ----------
    psi0_axis = np.linspace(np.radians(0.5), np.radians(40), 200)
    fig, ax = plt.subplots(figsize=(7, 5))
    for Om in [0.02, 0.05, 0.15, 0.3, 0.9]:
        ratio = A_psi_from_psi0(psi0_axis, Om, B) / A
        ax.plot(np.degrees(psi0_axis), ratio, lw=1.8, label=f"Omega={Om}")
    ax.axhline(1.0, color="k", ls="--", lw=1.5, label="associated-flow limit (A_psi=A)")
    ax.axvspan(2, 15, color="#ffd0a0", alpha=0.6, label="typical measured psi0 (2~15 deg)")
    ax.set_xlabel("measured dilatancy angle psi0 (deg)")
    ax.set_ylabel(r"potential-slope ratio $A_\psi/A$  (Eq.46)")
    ax.set_title("Fig.6  Non-associated correction (Eq.46)")
    ax.set_ylim(0, 2.0); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUTDIR, "fig6_nonassociated.png")); plt.close(fig)

    # ---------- Table 3 复现：主动试验 Chen 2019 ----------
    print()
    print("Table 3  Active trapdoor test (Chen et al., 2019)")
    print("-"*66)
    phi_t = np.radians(31.6); l_t = 0.125; h_t = 0.250; eta_t = 1.0
    meas_Ga = 0.47; psi_meas = np.radians(24.25)
    Ga_assoc = Ga_linear(eta_t, phi_t)
    Ga_non   = Ga_linear(eta_t, psi_meas)
    hcr_assoc = h_cr_linear(l_t, phi_t)*1000   # mm
    hcr_non   = h_cr_linear(l_t, psi_meas)*1000
    print(f"  等效剪胀角 psi   : measured 24.25 | assoc {np.degrees(phi_t):.1f} | non-assoc 24.25")
    print(f"  临界厚度 h_cr    : assoc {hcr_assoc:.0f} mm | non-assoc {hcr_non:.0f} mm  (h=250 mm)")
    print(f"  破坏模式(h=250)  : assoc -> h_cr<h, 深埋 | non-assoc -> h_cr>h, 贯通地表")
    print(f"  主动土压力比 Ga  : measured {meas_Ga} | assoc {Ga_assoc:.3f} "
          f"({(Ga_assoc/meas_Ga-1)*100:+.1f}%) | non-assoc {Ga_non:.2f} "
          f"({(Ga_non/meas_Ga-1)*100:+.1f}%)")
    print("  (论文: assoc 0.385/-18.1%, non-assoc 0.55/+17.0%)")

    # ---------- Table 4 复现：被动试验 Rui 2025b ----------
    print()
    print("Table 4  Passive trapdoor tests (Rui et al., 2025b)")
    print("-"*66)
    #  eta, phi_peak(deg), theta_max(deg), measured_Gb
    rows = [(0.5, 57.36, 26.5, 1.48),
            (1.0, 54.99, 30.0, 1.88),
            (1.5, 53.62, 35.5, 2.35),
            (2.0, 52.64, 40.5, 2.88),
            (3.0, 51.25, 44.5, 3.86)]
    print("  eta   phi_pk  theta   meas | assoc  dev  | non-assoc dev")
    dev_assoc=[]; dev_non=[]
    for eta_r, phi_pk, th, meas in rows:
        assoc = Gb_linear(eta_r, np.radians(phi_pk))
        non   = Gb_linear(eta_r, np.radians(th))
        da = (assoc/meas-1)*100; dn = (non/meas-1)*100
        dev_assoc.append(da); dev_non.append(dn)
        print(f"  {eta_r:3.1f}  {phi_pk:5.1f}  {th:5.1f}  {meas:4.2f} | "
              f"{assoc:5.2f} {da:+5.1f}% |  {non:5.2f} {dn:+5.1f}%")
    print(f"  Mean dev: assoc {np.mean(dev_assoc):+.1f}%  (论文 +25.4)  | "
          f"non-assoc {np.mean(dev_non):+.1f}%  (论文 -9.5)")
    print("="*66)

    print(f"\n图片已保存到: {OUTDIR}")
    for f in ["fig1_slip_surface.png","fig2_stress_field.png","fig3_dilatancy.png",
              "fig4_zoning.png","fig5_Ka_depth.png","fig6_nonassociated.png"]:
        print("  ", f)


if __name__ == "__main__":
    main()
