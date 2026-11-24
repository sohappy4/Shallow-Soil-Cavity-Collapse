# -*- coding: utf-8 -*-
"""
fem2d.py —— 二维平面应变弹塑性有限元（FEM-SRF）
================================================
横向引用成熟有限元方法，扩展专著沿滑面一维场为真正的二维应力/塑性区场。
依据用户 2026-11-24 授权：专著理论不足时可深入扩展与横向引用。

方法（标准岩土有限元，非专著原有内容，已在文档显式标注来源）：
  * 4 节点四边形等参单元（Q4），平面应变，2×2 高斯积分；
  * Drucker–Prager 屈服面 f = q + m·p - k ≤ 0（拉应力为正，压缩 p<0），
    其中 m = 6sinφ/(3-sinφ), k = 6c·cosφ/(3-sinφ)（压缩子午线与 MC 匹配），
    单轴抗压强度严格等于 UCS = 2c·cosφ/(1-sinφ)（见 selfcheck）；
  * 径向返回映射（backward Euler），非关联流动取剪胀角 ψ：
    m_ψ = 6sinψ/(3-sinψ)；Δγ = f_trial/(3G + K·m·m_ψ)；
  * 位移控制活门（底边洞段施加向下位移），修正 Newton 迭代，
    记录活门反力—位移曲线，峰值即极限荷载；
  * 强度折减 Fs: c'=c/Fs, tanφ'=tanφ/Fs，扫描估计临界折减系数。

【诚实边界】
  - 这是 FEM 强度折减（应力场 + 塑性区），不是严格极限分析 FELA 上下界；
  - 本构是 DP 圆化（非幂型 Hoek–Brown），与专著强度准则不同族，量级对照而非同一理论；
  - 小变形、理想弹塑性（无硬化）；网格较粗，仅作定性—半定量演示。
"""
from __future__ import annotations

from typing import Dict, Any, List, Tuple
import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

__all__ = [
    "dp_params", "dp_yield", "dp_radial_return", "selfcheck_dp",
    "build_mesh", "solve_fem", "DEFAULT_FEM_PARAMS",
]


DEFAULT_FEM_PARAMS: Dict[str, float] = dict(
    E=20000.0,      # 变形模量 kPa
    nu=0.30,        # 泊松比
    c=25.0,         # 黏聚力 kPa
    phi=20.0,       # 摩擦角 deg
    psi=0.0,        # 剪胀角 deg（0=非关联无剪胀）
    gamma=18.0,     # 容重 kN/m3
    sigma_v0=0.0,   # 地表超载 kPa
)


# ---------------------------------------------------------------------------
# 1. Drucker–Prager 本构
# ---------------------------------------------------------------------------
def dp_params(c: float, phi_deg: float, psi_deg: float = None):
    """返回 (m, k, m_psi)。拉应力为正，压缩 p<0。"""
    phi = math.radians(float(phi_deg))
    sp_, cp = math.sin(phi), math.cos(phi)
    den = 3.0 - sp_
    m = 6.0 * sp_ / den
    k = 6.0 * float(c) * cp / den
    if psi_deg is None:
        m_psi = m
    else:
        ps = math.radians(float(psi_deg))
        m_psi = 6.0 * math.sin(ps) / (3.0 - math.sin(ps))
    return float(m), float(k), float(m_psi)


def dp_yield(sig: np.ndarray, m: float, k: float) -> float:
    """屈服函数 f = q + m·p - k；f<=0 弹性。sig = [sxx, syy, szz, txy]。"""
    sxx, syy, szz, txy = sig
    p = (sxx + syy + szz) / 3.0
    s1, s2, s3 = sxx - p, syy - p, szz - p
    j2 = 0.5 * (s1 * s1 + s2 * s2 + s3 * s3) + txy * txy
    q = math.sqrt(3.0 * max(j2, 0.0))
    return float(q + m * p - k)   # 转 Python float，避免 numpy 标量污染 JSON


def dp_radial_return(sig: np.ndarray, de: np.ndarray, E: float, nu: float,
                     m: float, k: float, m_psi: float) -> Tuple[np.ndarray, float, bool]:
    """对一个高斯点做弹塑性径向返回。

    sig: 上一收敛应力 [sxx,syy,szz,txy]
    de:  应变增量 [exx,eyy,gxy]（工程剪应变），平面应变 ezz=0
    返回 (新应力, 等效塑性应变增量, 是否屈服)
    """
    G = E / (2.0 * (1.0 + nu))
    lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    exx, eyy, gxy = de
    # 弹性试算（平面应变：szz 由 exx+eyy 决定）
    sxx = sig[0] + (lam + 2 * G) * exx + lam * eyy
    syy = sig[1] + lam * exx + (lam + 2 * G) * eyy
    szz = sig[2] + lam * (exx + eyy)
    txy = sig[3] + G * gxy
    trial = np.array([sxx, syy, szz, txy], dtype=float)

    f = dp_yield(trial, m, k)
    if f <= 0.0:
        return trial, 0.0, False

    p = (sxx + syy + szz) / 3.0
    s1, s2, s3 = sxx - p, syy - p, szz - p
    j2 = 0.5 * (s1 * s1 + s2 * s2 + s3 * s3) + txy * txy
    q = math.sqrt(3.0 * max(j2, 0.0))
    if q < 1e-12:
        return trial, 0.0, False
    K = lam + 2.0 * G / 3.0

    def stress_at(dgam):
        # 用当前试算偏应力的法向（回映射内部固定 n，标准 backward Euler）
        pg = (trial[0] + trial[1] + trial[2]) / 3.0
        s1g, s2g, s3g = trial[0] - pg, trial[1] - pg, trial[2] - pg
        qg = math.sqrt(3.0 * (0.5 * (s1g**2 + s2g**2 + s3g**2) + trial[3] ** 2))
        if qg < 1e-12:
            n = np.zeros(4)
        else:
            # 工程分量表示下，正应力与 txy 共用同一法向系数 1.5/q。
            # 若 txy 误用 3/q，剪应力缩减快一倍，返回路径非径向、
            # f(dgam) 非单调，二分失效并使应力逐帧爆炸。
            n = np.array([1.5 * s1g / qg, 1.5 * s2g / qg,
                          1.5 * s3g / qg, 1.5 * trial[3] / qg])
        cs = K * m_psi * dgam
        return trial - dgam * 2 * G * n - np.array([cs, cs, cs, 0.0])

    # 一致性求根：f(dgam) 随 dgam 单调递减时纯二分稳健。
    dgam0 = f / (3.0 * G + K * m * m_psi)
    lo, hi = 0.0, max(2.0 * dgam0, 1e-4)
    while dp_yield(stress_at(hi), m, k) > 0.0 and hi < 1e3:
        hi *= 2.0

    if dp_yield(stress_at(hi), m, k) > 0.0:
        # DP 锥对静水拉应力无上界：当 m·p > k 时，即使 q→0 也有 f>0，
        # 塑性乘子无解（活门下方土体被拉张时会遇到）。
        # 按岩土有限元惯例施加 Rankine 拉力截断：土体不抗拉，
        # 把正（拉）正应力截断到 0；残余剪应力再二分缩放到屈服面内。
        sig_new = stress_at(min(dgam0, 1.0))
        sig_new = np.array([min(sig_new[0], 0.0),
                            min(sig_new[1], 0.0),
                            min(sig_new[2], 0.0),
                            sig_new[3]], dtype=float)
        f2 = dp_yield(sig_new, m, k)
        if f2 > 1e-8:
            slo, shi = 0.0, 1.0
            for _ in range(60):
                sm = 0.5 * (slo + shi)
                tmp = sig_new.copy()
                tmp[3] *= sm
                if dp_yield(tmp, m, k) > 0.0:
                    shi = sm
                else:
                    slo = sm
                if shi - slo < 1e-12:
                    break
            sig_new[3] *= 0.5 * (slo + shi)
        return sig_new, float(max(dgam0, 0.0)), True

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if dp_yield(stress_at(mid), m, k) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12 * max(1.0, hi):
            break
    dgam = 0.5 * (lo + hi)
    sig_new = stress_at(dgam)
    return sig_new, float(dgam), True


def selfcheck_dp() -> Dict[str, Any]:
    """解析自检：单轴抗压强度 = UCS；屈服后 f≈0；弹性区不屈服。"""
    c, phi, psi = 25.0, 20.0, 10.0
    E, nu = 20000.0, 0.3
    m, k, mp = dp_params(c, phi, psi)
    ucs = 2 * c * math.cos(math.radians(phi)) / (1 - math.sin(math.radians(phi)))
    # 单轴压缩应变：sig=[-ucs,0,0] 需 eps_yy=eps_zz=-nu*ucs/E；这里直接测屈服函数
    sig_uni = np.array([-ucs, 0.0, 0.0, 0.0])  # 近似（平面应力），f 应≈0
    f_uni = dp_yield(sig_uni, m, k)
    # 径向返回：给一个大剪应变试算（剪切不产生侧限围压，必屈服），返回后 f≈0
    sig0 = np.zeros(4)
    de = np.array([0.0, 0.0, 0.02])   # 纯剪 gxy
    new, dgam, yld = dp_radial_return(sig0, de, E, nu, m, k, mp)
    f_after = dp_yield(new, m, k)
    # 小应变弹性
    new_e, dg_e, yld_e = dp_radial_return(sig0, np.array([1e-6, 0, 0]), E, nu, m, k, mp)
    checks = {
        "UCS_formula_positive": ucs > 0,
        "uniaxial_f_near_zero": abs(f_uni) < 0.05 * ucs,
        "return_maps_to_surface": (not yld) or abs(f_after) < 1e-6,
        "return_yields_under_big_strain": yld and dgam > 0,
        "small_strain_elastic": (not yld_e) and abs(dg_e) < 1e-15,
    }
    return dict(checks=checks, all_pass=all(checks.values()),
                ucs=ucs, f_uni=f_uni, f_after=f_after, m=m, k=k, m_psi=mp)


# ---------------------------------------------------------------------------
# 2. 网格与单元
# ---------------------------------------------------------------------------
def build_mesh(nx: int, ny: int, x0: float, x1: float,
               y0: float, y1: float) -> Tuple[np.ndarray, np.ndarray]:
    """结构化四边形网格；节点按行优先，单元逆时针。"""
    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    X, Y = np.meshgrid(xs, ys)
    nodes = np.column_stack([X.ravel(), Y.ravel()])       # (nnode,2)
    elems = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            b = a + 1
            d = (j + 1) * (nx + 1) + i
            cnode = d + 1
            elems.append([a, b, cnode, d])
    return nodes, np.array(elems, dtype=int)


_GAUSS = [(-1 / math.sqrt(3), -1 / math.sqrt(3)),
          (1 / math.sqrt(3), -1 / math.sqrt(3)),
          (1 / math.sqrt(3), 1 / math.sqrt(3)),
          (-1 / math.sqrt(3), 1 / math.sqrt(3))]


def _q4_shape(xi: float, eta: float):
    """返回 N(4), dN/dxi(4), dN/deta(4)。"""
    N = np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                  (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)]) * 0.25
    dNx = np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)]) * 0.25
    dNy = np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]) * 0.25
    return N, dNx, dNy


def _D_plane(E: float, nu: float) -> np.ndarray:
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    D = np.zeros((3, 3))
    D[0, 0] = lam + 2 * G
    D[1, 1] = lam + 2 * G
    D[0, 1] = D[1, 0] = lam
    D[2, 2] = G
    return D


def _b_matrix(dNdx: np.ndarray, dNdy: np.ndarray) -> np.ndarray:
    """平面应变 B (3×8)，应变向量 [exx,eyy,gxy]。"""
    B = np.zeros((3, 8))
    for a in range(4):
        B[0, 2 * a] = dNdx[a]
        B[1, 2 * a + 1] = dNdy[a]
        B[2, 2 * a] = dNdy[a]
        B[2, 2 * a + 1] = dNdx[a]
    return B


# ---------------------------------------------------------------------------
# 3. 有限元求解（修正 Newton + 位移控制活门）
# ---------------------------------------------------------------------------
def solve_fem(params: Dict[str, float], opts: Dict[str, Any] = None) -> Dict[str, Any]:
    """平面应变活门试验有限元求解。

    params: E/nu/c/phi/psi/gamma/sigma_v0，以及 l(半宽)/h(覆土)
    opts:   nx/ny(网格)、n_steps(活门位移步数)、door_frac(终位移/洞宽)、
            sr(强度折减 Fs，默认1)，max_iter/tol
    返回：应力场、塑性区、活门荷载-位移曲线、峰值荷载、收敛信息。
    """
    o = dict(nx=20, ny=16, n_steps=24, door_frac=0.10, sr=1.0,
             max_iter=60, tol=1e-3, relax=0.4)
    if opts:
        o.update({k: v for k, v in opts.items() if v is not None})
    nx, ny = max(6, min(60, int(o["nx"]))), max(6, min(50, int(o["ny"])))
    n_steps = max(4, min(80, int(o["n_steps"])))
    sr = max(0.3, float(o["sr"]))

    E = float(params["E"]) * float(params.get("E_factor", 1.0))
    nu = float(params["nu"])
    c = float(params["c"]) / sr
    phi_deg = math.degrees(math.atan(math.tan(math.radians(float(params["phi"]))) / sr))
    psi_deg = float(params.get("psi", 0.0))
    gamma = float(params.get("gamma", 18.0))
    q2 = float(params.get("sigma_v0", 0.0))
    l = float(params["l"])
    h = float(params["h"])
    D = 2.0 * l
    if nu <= 0 or nu >= 0.49:
        raise ValueError("nu 必须在 (0, 0.49)")
    if phi_deg <= 0 or c <= 0:
        raise ValueError("折减后 c/phi 必须为正")

    m, k, mp = dp_params(c, phi_deg, psi_deg)
    Dm = _D_plane(E, nu)

    # 区域：x∈[-0.8D, 1.8D]，y∈[0, h]（y 向上，地表 y=h，洞底 y=0）
    x0, x1 = -0.8 * D, 1.8 * D
    y0, y1 = 0.0, h
    nodes, elems = build_mesh(nx, ny, x0, x1, y0, y1)
    nnode = nodes.shape[0]
    ndof = 2 * nnode
    nele = elems.shape[0]

    # 单元几何预计算：每个单元 4 高斯点的 B 矩阵与 Jacobian 权重
    Bgp = np.zeros((nele, 4, 3, 8))       # (nele, ngp, 3, 8)
    wgp = np.zeros((nele, 4))
    Ngp = np.zeros((nele, 4, 4))
    for e in range(nele):
        en = nodes[elems[e]]                # 4×2
        for gi, (xi, eta) in enumerate(_GAUSS):
            N, dNx, dNy = _q4_shape(xi, eta)
            J = np.array([[np.dot(dNx, en[:, 0]), np.dot(dNx, en[:, 1])],
                          [np.dot(dNy, en[:, 0]), np.dot(dNy, en[:, 1])]])
            detJ = np.linalg.det(J)
            if detJ <= 0:
                raise ValueError("单元 %d Jacobian 非正" % e)
            Jinv = np.linalg.inv(J)
            dNxy = Jinv @ np.vstack([dNx, dNy])   # 2×4
            Bgp[e, gi] = _b_matrix(dNxy[0], dNxy[1])
            wgp[e, gi] = detJ
            Ngp[e, gi] = N

    # 组装弹性刚度（稀疏）
    rows, cols, vals = [], [], []
    for e in range(nele):
        dofs = (elems[e][:, None] * 2 + np.arange(2)).ravel()
        Ke = np.zeros((8, 8))
        for gi in range(4):
            B = Bgp[e, gi]
            Ke += B.T @ Dm @ B * wgp[e, gi]
        for i in range(8):
            for j in range(8):
                rows.append(dofs[i]); cols.append(dofs[j]); vals.append(Ke[i, j])
    K = sp.coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsc()

    # 体力（自重）与地表超载：等价节点荷载
    Fbody = np.zeros(ndof)
    for e in range(nele):
        dofs = (elems[e][:, None] * 2 + np.arange(2)).ravel()
        fe = np.zeros(8)
        for gi in range(4):
            fe[1::2] += -gamma * Ngp[e, gi] * wgp[e, gi]  # y 向上，重力向下
        Fbody[dofs] += fe
    # 地表超载施加于顶层节点
    topy = np.max(nodes[:, 1])
    for nd in range(nnode):
        if abs(nodes[nd, 1] - topy) < 1e-9:
            Fbody[2 * nd + 1] -= q2 * (x1 - x0) / (nx + 1)

    # 边界条件
    fixed = []
    for nd in range(nnode):
        x, y = nodes[nd]
        if abs(x - x0) < 1e-9 or abs(x - x1) < 1e-9:   # 左右：水平固定
            fixed.append(2 * nd)
        if abs(y - y0) < 1e-9:                          # 底边：全固定
            fixed.append(2 * nd); fixed.append(2 * nd + 1)
    door_nodes = [nd for nd in range(nnode)
                  if abs(nodes[nd, 1] - y0) < 1e-9 and (0.0 - 1e-9) <= nodes[nd, 0] <= D + 1e-9]
    # 活门节点在底边固定集中，改由位移控制：从 fixed 中移除其竖向
    for nd in door_nodes:
        if 2 * nd + 1 in fixed:
            fixed.remove(2 * nd + 1)
    fixed = np.array(sorted(set(fixed)), dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)

    # 地应力初始化：全约束底边（含活门）+弹性解，建立初始重力场
    init_fixed = np.unique(np.concatenate([fixed, np.array([2 * nd + 1 for nd in door_nodes], dtype=int)]))
    init_free = np.setdiff1d(np.arange(ndof), init_fixed)
    u0 = np.zeros(ndof)
    rhs0 = Fbody.copy()
    u0[init_free] = spl.spsolve(K[init_free][:, init_free], rhs0[init_free])
    # 初始应力（弹性）
    sig_gp = np.zeros((nele, 4, 4))     # 每高斯点应力
    for e in range(nele):
        ue = u0[(elems[e][:, None] * 2 + np.arange(2)).ravel()]
        for gi in range(4):
            eps = Bgp[e, gi] @ ue
            sxx = Dm[0, 0] * eps[0] + Dm[0, 1] * eps[1]
            syy = Dm[1, 0] * eps[0] + Dm[1, 1] * eps[1]
            szz = Dm[0, 1] * (eps[0] + eps[1])   # 平面应变：szz = lam·(exx+eyy)
            txy = Dm[2, 2] * eps[2]
            sig_gp[e, gi] = [sxx, syy, szz, txy]

    # 活门位移控制：修正 Newton（初始应力法）
    # ------------------------------------------------------------------
    # 每个位移步保存步初收敛状态 (u_step0, sig_step0, ep_step0)。
    # 迭代中每轮都从同一参考状态出发，按当前候选位移的总增量重算应力，
    # 再以弹性刚度修正不平衡力；塑性应变不跨迭代累加，避免正反馈发散。
    # 活门竖向自由度由位移直接控制，不参与修正。
    door_total = float(o["door_frac"]) * D
    du = door_total / n_steps

    door_dofs = np.array([2 * nd + 1 for nd in door_nodes], dtype=int)
    door_set = set(door_dofs.tolist())
    active_free = np.array([d for d in free if d not in door_set], dtype=int)
    K_ff = K[active_free][:, active_free].tocsc()
    gp_dofs = [(elems[e][:, None] * 2 + np.arange(2)).ravel() for e in range(nele)]
    u = u0.copy()
    epsp_cum = np.zeros((nele, 4))
    curve = []
    converged_all = True
    max_iter, tol = int(o["max_iter"]), float(o["tol"])

    def internal_force(sig):
        F = np.zeros(ndof)
        for e in range(nele):
            fe = np.zeros(8)
            for gi in range(4):
                s = sig[e, gi]
                fe += Bgp[e, gi].T @ np.array([s[0], s[1], s[3]]) * wgp[e, gi]
            F[gp_dofs[e]] += fe
        return F

    def cell_fields(sig, ep):
        """把高斯点场聚合成单元中心场，全部转 Python 标量便于 JSON 序列化。"""
        cs = np.mean(sig, axis=1)
        ce = np.mean(ep, axis=1)
        qq = []
        for s in cs:
            pp = (s[0] + s[1] + s[2]) / 3.0
            j2 = 0.5 * ((s[0]-pp)**2 + (s[1]-pp)**2 + (s[2]-pp)**2) + s[3]**2
            qq.append(math.sqrt(3.0 * max(j2, 0.0)))
        yd = ce > 1e-8
        return dict(q=[float(v) for v in qq],
                    epsp=[float(v) for v in ce],
                    yielded=[bool(v) for v in yd],
                    plastic_frac=float(np.mean(yd)))

    # step0：未扰动初始反力（静止土压力对照点）
    Fint0 = internal_force(sig_gp)
    reaction0 = float(np.sum(Fint0[door_dofs] - Fbody[door_dofs]))
    curve.append(dict(step=0, door_down=0.0, reaction_kN=reaction0,
                      q_kPa=reaction0 / D, plastic_frac=0.0,
                      iters=0, residual=0.0))
    field_history = [dict(step=0, door_down=0.0,
                          **cell_fields(sig_gp, epsp_cum))]

    for step in range(1, n_steps + 1):
        target = -du * step                       # y 向上，向下位移为负
        u_step0 = u.copy()
        sig_step0 = sig_gp.copy()
        ep_step0 = epsp_cum.copy()
        u_it = u_step0.copy()
        u_it[door_dofs] = target
        sig_it = sig_step0.copy()
        unorm = 1.0
        sig_new = sig_step0.copy()
        ep_new = ep_step0.copy()
        for it in range(max_iter):
            Fint = internal_force(sig_it)
            resid = Fint - Fbody
            du_corr = np.zeros(ndof)
            du_corr[active_free] = float(o["relax"]) * spl.spsolve(
                K_ff, -resid[active_free])
            u_it = u_it + du_corr
            u_it[door_dofs] = target
            # 从步初参考状态重算应力（关键：不跨迭代累加塑性应变）
            sig_new = np.empty_like(sig_step0)
            ep_new = ep_step0.copy()
            for e in range(nele):
                ue0, ue1 = u_step0[gp_dofs[e]], u_it[gp_dofs[e]]
                for gi in range(4):
                    deps = Bgp[e, gi] @ (ue1 - ue0)
                    sn, dg, yld = dp_radial_return(
                        sig_step0[e, gi], deps, E, nu, m, k, mp)
                    sig_new[e, gi] = sn
                    ep_new[e, gi] += dg
            unorm = np.linalg.norm(du_corr) / (np.linalg.norm(u_it) + 1e-12)
            sig_it = sig_new
            if unorm < tol:
                break
        u = u_it
        sig_gp, epsp_cum = sig_new, ep_new
        if unorm >= tol * 10:
            converged_all = False
        Fint = internal_force(sig_gp)
        reaction = float(np.sum(Fint[door_dofs] - Fbody[door_dofs]))
        q_door = reaction / D                     # 活门上托压力，kPa（正）
        curve.append(dict(step=step, door_down=float(du * step),
                          reaction_kN=reaction,
                          q_kPa=float(q_door),
                          plastic_frac=float(np.mean(epsp_cum > 1e-8)),
                          iters=int(it + 1), residual=float(unorm)))
        field_history.append(dict(step=step, door_down=float(du * step),
                                  **cell_fields(sig_gp, epsp_cum)))

    # 峰值荷载与最终场
    qs = [p["q_kPa"] for p in curve]
    peak = max(qs) if qs else 0.0
    peak_step = int(np.argmax(qs)) + 1 if qs else 0
    # 单元中心应力（4 高斯点平均）
    cell_sig = np.mean(sig_gp, axis=1)
    cell_ep = np.mean(epsp_cum, axis=1)
    cell_centers = np.mean(nodes[elems], axis=1)
    yielded = cell_ep > 1e-8
    # von Mises 等效应力 q = sqrt(3 J2)，用于云图
    cell_q = []
    for s in cell_sig:
        pp = (s[0] + s[1] + s[2]) / 3.0
        j2 = 0.5 * ((s[0]-pp)**2 + (s[1]-pp)**2 + (s[2]-pp)**2) + s[3]**2
        cell_q.append(math.sqrt(3.0 * max(j2, 0.0)))
    return dict(
        solution_type="fem-srf",
        meta=dict(
            method="FEM plane-strain Q4 + DP radial return + displacement-controlled trapdoor",
            source="external-mature-FEM (not from monograph; per user authorization 2026-11-24)",
            boundary="FEM strength-reduction stress/plastic field; NOT rigorous FELA bounds; "
                     "DP circle (not power-law Hoek-Brown); small-strain perfect plasticity",
            converged=bool(converged_all), sr=sr,
            c_eff=c, phi_eff_deg=phi_deg, n_steps=n_steps,
            grid=(nx, ny), n_elements=nele, n_nodes=nnode,
        ),
        params=dict(E=E, nu=nu, c=c, phi_deg=phi_deg, psi_deg=psi_deg,
                    gamma=gamma, l=l, h=h, D=D, door_total=door_total),
        curve=curve, peak_q_kPa=float(peak), peak_step=peak_step,
        field_history=field_history,
        field=dict(
            centers=cell_centers.tolist(),
            sxx=cell_sig[:, 0].tolist(), syy=cell_sig[:, 1].tolist(),
            txy=cell_sig[:, 3].tolist(),
            p=((cell_sig[:, 0] + cell_sig[:, 1] + cell_sig[:, 2]) / 3.0).tolist(),
            q=cell_q,
            epsp=cell_ep.tolist(), yielded=yielded.tolist(),
            plastic_frac=float(np.mean(yielded)),
        ),
    )
