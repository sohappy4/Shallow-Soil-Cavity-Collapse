# -*- coding: utf-8 -*-
"""
l2_shooting.py —— L2 打靶—积分求解器
=====================================
严格实现专著「算法 17.1　浅埋土洞塌落滑动面数值求解（主动塌落·内倾支）」：
  S = sigma_c A^(1/(1-B)) B^(B/(1-B))
  p = (sigma_t+q2)/(B*S)                       # 式(16.13)
  |f'(x2)| = p^(1-B)
  欧拉方程通解(16.10a): f' = -(gamma*x/S + m)^((1-B)/B)     [主动 s=-1]
  打靶变量 x2：由 f(x2)=h 起积分，令 f(l)=0 为残差，二分收敛
  被动 s=+1 为镜像：x -> 2l-x，方程右端反号
三级筛选：(i) u·s<0  (ii) Φ(s)>=0  (iii) sigma(x2)>=-sigma_t
"""
from __future__ import annotations

import math
from typing import Dict, Any, List
import numpy as np

__all__ = ["solve_shooting"]

TOL = 1e-9
MAX_ITER = 200


def _S(sigma_c, A, B):
    if abs(B - 1.0) < 1e-9:
        return None
    return sigma_c * A ** (1.0 / (1.0 - B)) * B ** (B / (1.0 - B))


def _rk4_integrate(x0: float, x1: float, f0: float, m: float,
                   S: float, B: float, gamma: float, sign: float,
                   nsteps: int = 400):
    """四阶 Runge-Kutta 积分 f' = sign*(gamma*x/S+m)^((1-B)/B)。"""
    hstep = (x1 - x0) / nsteps
    xs = [x0]
    fs = [f0]
    x, f = x0, f0
    expo = (1.0 - B) / B

    def deriv(xx):
        base = gamma * xx / S + m
        if base <= 0:
            return None
        return sign * base ** expo

    for _ in range(nsteps):
        k1 = deriv(x)
        k2 = deriv(x + hstep / 2)
        k3 = deriv(x + hstep / 2)
        k4 = deriv(x + hstep)
        if any(k is None for k in (k1, k2, k3, k4)):
            return None, None
        f += hstep * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        x += hstep
        xs.append(x)
        fs.append(f)
    return np.array(xs), np.array(fs)


def solve_shooting(params: Dict[str, float], mode: str = "active",
                   tol: float = TOL, max_iter: int = MAX_ITER) -> Dict[str, Any]:
    """
    打靶求解滑动面。返回字典（可 JSON 序列化）：
      x, F        滑动面坐标（F 为自洞顶量起的高度，0..h）
      x2, L_c     出露点、水平内缩量
      m, t        欧拉通解参数
      residual    打靶残差 f(l)
      iterations  二分迭代次数
      screening   {'direction':bool,'dissipation':bool,'boundary':bool,'all':bool}
      admissible  三级筛选是否全过
      h_cr_implied 令 L_c=l 反解的临界埋深
      closed_form_Lc 闭式 L_c（交叉校验用）
      lcf_error    与闭式的相对误差
    """
    for k in ("sigma_c", "A", "B", "sigma_t", "gamma", "q2", "l", "h"):
        if k not in params:
            raise ValueError("缺少参数 %s" % k)
    sc, A, B = float(params["sigma_c"]), float(params["A"]), float(params["B"])
    st, g, q2 = float(params["sigma_t"]), float(params["gamma"]), float(params["q2"])
    l, h = float(params["l"]), float(params["h"])
    if not (0 < B < 1):
        raise ValueError("打靶求解要求 0<B<1；B=1 请使用线性极限分支")
    if not (sc > 0 and A > 0 and g > 0 and l > 0 and h > 0 and st >= 0 and q2 >= 0):
        raise ValueError("参数不合法")
    if mode not in ("active", "passive"):
        raise ValueError("mode 只能 active / passive")

    S = _S(sc, A, B)
    p = (st + q2) / (B * S)              # 式(16.13)
    sign = -1.0 if mode == "active" else +1.0
    s = -1 if mode == "active" else +1   # 几何分支
    u = +1 if mode == "active" else -1   # 运动方向

    # 闭式 L_c（用于交叉校验与打靶初值参考）
    L_c_closed = S / g * ((p + g * h / (B * S)) ** B - p ** B)
    h_cr_implied = S * B / g * ((p ** B + g * l / S) ** (1.0 / B) - p)

    # --- 打靶：统一在主动（内倾）支上以 x2a 为变量，残差 r=f(l)。
    #     被动支是主动的纯镜像 x->2l-x（专著算法 17.1 末段），不单独积分。 ---
    def residual(x2a):
        if x2a <= 0 or x2a >= l:
            return None
        m = p ** B - g * x2a / S         # 由 |f'(x2)|=p^(1-B) 定 m
        xs, fs = _rk4_integrate(x2a, l, h, m, S, B, g, -1.0)
        if xs is None:
            return None
        return fs[-1]

    lo, hi = 1e-10 * l, l * (1 - 1e-10)
    r_lo, r_hi = residual(lo), residual(hi)
    # h>h_cr 时无贯通解（残差同号），如实报告
    if r_lo is None or r_hi is None or r_lo * r_hi > 0:
        L_c = L_c_closed             # 仍给闭式形态作参考
        x2 = l - L_c if mode == "active" else l + L_c
        m = p ** B - g * min(x2, 2 * l - x2) / S
        xs, fs = _rk4_integrate(min(x2, 2 * l - x2), l, h, m, S, B, g, sign)
        return dict(
            x=xs.tolist() if xs is not None else [],
            F=fs.tolist() if fs is not None else [],
            x2=x2, L_c=L_c, m=m, t=S * B / g * (p ** B + g * L_c / S) ** (1.0 / B),
            residual=float(fs[-1]) if fs is not None else float("nan"),
            iterations=0, converged=False,
            no_through_going=True,
            screening=dict(direction=(u * s < 0), dissipation=None, boundary=None, all=False),
            admissible=False, h_cr_implied=h_cr_implied,
            closed_form_Lc=L_c_closed, lcf_error=0.0,
        )

    it = 0
    for it in range(1, max_iter + 1):
        mid = 0.5 * (lo + hi)
        rm = residual(mid)
        if rm is None:
            break
        if abs(rm) < tol or (hi - lo) < tol * l:
            lo = hi = mid
            break
        if r_lo * rm <= 0:
            hi, r_hi = mid, rm
        else:
            lo, r_lo = mid, rm

    x2a = 0.5 * (lo + hi)                       # 主动侧出露点（打靶基准）
    L_c = l - x2a
    m = p ** B - g * x2a / S
    xs, fs = _rk4_integrate(x2a, l, h, m, S, B, g, -1.0)   # 主动支积分
    res = float(fs[-1])

    if mode == "passive":
        # 纯镜像输出：x -> 2l-x，F 不变（专著算法 17.1 末段）
        x_out = [2.0 * l - xx for xx in xs[::-1]]
        F_out = fs[::-1].tolist()
        x2 = 2.0 * l - x2a
    else:
        x_out, F_out, x2 = xs.tolist(), fs.tolist(), x2a

    # --- 三级筛选 ---
    # (i) u·s<0
    c1 = (u * s < 0)
    # (ii) Φ(s)=(sgn(u)s*B+1)*Q - sgn(u)s*sigma_t >=0，张开支 Q=S*p
    Q = S * p
    phi = (u * s * B + 1.0) * Q - u * s * st
    c2 = bool(phi >= 0)
    # (iii) 横截条件 sigma(x2)=q2 >= -sigma_t
    c3 = bool(q2 >= -st)
    # 附加几何自检：单调、非负、闭合（积分基准支 fs 从 h 降到 0）
    mono = all(fs[i + 1] <= fs[i] + 1e-7 for i in range(len(fs) - 1))
    closed = abs(res) < 1e-5
    screening = dict(direction=c1, dissipation=c2, boundary=c3,
                     monotonic=bool(mono), closed=bool(closed), phi=float(phi))
    admissible = c1 and c2 and c3 and mono and closed

    lcf_err = abs(L_c - L_c_closed) / L_c_closed if L_c_closed > 0 else 0.0
    return dict(
        x=x_out, F=F_out, x2=x2, L_c=L_c, m=m,
        t=S * B / g * (p ** B + g * L_c / S) ** (1.0 / B),
        residual=res, iterations=it, converged=abs(res) < tol,
        no_through_going=False, screening=screening, admissible=admissible,
        h_cr_implied=h_cr_implied, closed_form_Lc=L_c_closed, lcf_error=float(lcf_err),
        S=S, p=p,
    )
