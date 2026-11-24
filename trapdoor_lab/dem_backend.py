# -*- coding: utf-8 -*-
"""DEM 独立进程后端适配器协议（方案 M8 许可证隔离边界）。

方案对 DEM 后端的硬约束（实施方案 §五·许可证边界）：
  Yade 为 GPL、PFC 为商业软件 —— 必须以**独立进程 + 文件/标准输入输出通信**
  接入，**不链接进主程序、不随包分发**，规避许可证传染。

本模块实现：
  1. 统一的 JSON 请求/响应协议（stdin/stdout 或文件）；
  2. `local` 参考后端 —— 直接调用本仓库自带的纯 Python `dem2d.py`
     （MIT，无传染风险，开箱可用）；
  3. `ExternalDEMBackend` 基类 —— 定义外部进程命令行的调用契约，
     **不包含任何 Yade/PFC 具体代码**，由使用者在自己的环境自行提供后端脚本。

协议版本 v1：
  请求 {"protocol":"dem/v1", "backend":"local|external", "params":{...},
        "opts":{...}, "command":"run"}
  响应 {"protocol":"dem/v1", "backend":..., "status":"ok|error",
        "result":{meta,frames,curve,boundary}, "error":...}

CLI：
  python dem_backend.py --backend local --input request.json --output response.json
  echo '<json>' | python dem_backend.py --backend local     # stdin/stdout 模式
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from typing import Any, Dict, Optional

from dem2d import solve_dem, selfcheck_dem, DEFAULT_DEM_PARAMS, DEFAULT_DEM_OPTS

__all__ = [
    "PROTOCOL", "make_request", "run_backend", "LocalDEMBackend",
    "ExternalDEMBackend", "BackendError", "selfcheck_backend",
]

PROTOCOL = "dem/v1"


class BackendError(RuntimeError):
    """后端协议错误（参数非法、响应格式错误、外部进程失败）。"""


def make_request(params: Dict[str, float] = None,
                 opts: Dict[str, Any] = None,
                 backend: str = "local") -> Dict[str, Any]:
    """构造一个 v1 请求对象。"""
    return dict(protocol=PROTOCOL, backend=backend,
                params=dict(params or {}), opts=dict(opts or {}), command="run")


class LocalDEMBackend:
    """内置参考后端：进程内调用 dem2d（MIT 自实现，无许可证传染）。"""

    name = "local"

    def run(self, params: Dict[str, float] = None,
            opts: Dict[str, Any] = None) -> Dict[str, Any]:
        return solve_dem(params, opts)

    def selfcheck(self) -> Dict[str, Any]:
        return selfcheck_dem()


class ExternalDEMBackend:
    """外部独立进程后端基类（Yade/PFC 等）。

    仅定义**通信契约**，不提供、不链接、不分发任何 GPL/商业求解代码。
    使用者需自行实现 `command` 指定的脚本，使其：
      - 从 stdin 读取 v1 请求 JSON；
      - 向 stdout 写 v1 响应 JSON（frames 格式与 dem2d 一致）；
      - 不依赖本仓库的 Python 环境（完全独立进程）。

    构造函数参数：
      command: 外部脚本命令行（list[str]），如 ["/usr/bin/python3", "/opt/dem/yade_bridge.py"]
      timeout: 子进程超时（秒）
    """

    name = "external"

    def __init__(self, command: list, timeout: float = 600.0):
        if not isinstance(command, (list, tuple)) or not command:
            raise BackendError("外部后端 command 必须是非空 argv 列表")
        self.command = [str(c) for c in command]
        self.timeout = float(timeout)

    def run(self, params: Dict[str, float] = None,
            opts: Dict[str, Any] = None) -> Dict[str, Any]:
        req = make_request(params, opts, backend="external")
        payload = json.dumps(req, ensure_ascii=False).encode("utf-8")
        try:
            proc = subprocess.run(
                self.command, input=payload, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, timeout=self.timeout, check=False)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise BackendError("外部 DEM 后端进程失败：%s" % type(exc).__name__) from exc
        if proc.returncode != 0:
            msg = proc.stderr.decode("utf-8", errors="replace")[-300:]
            raise BackendError("外部 DEM 后端退出码 %d：%s" % (proc.returncode, msg))
        try:
            resp = json.loads(proc.stdout.decode("utf-8"))
        except (ValueError, UnicodeDecodeError) as exc:
            raise BackendError("外部 DEM 后端响应不是合法 JSON") from exc
        return _validate_response(resp, expect_backend="external")

    def selfcheck(self) -> Dict[str, Any]:
        """外部后端自检：只验证协议可用性，不验证其内部物理。"""
        try:
            r = self.run(dict(DEFAULT_DEM_PARAMS, n_particles=30),
                         dict(DEFAULT_DEM_OPTS, n_steps=10, settle_steps=10,
                              max_particles=60))
            return dict(protocol_ok=True, n_frames=len(r.get("frames", [])),
                        boundary=r.get("boundary", ""))
        except BackendError as exc:
            return dict(protocol_ok=False, error=str(exc))


def _validate_response(resp: Dict[str, Any], expect_backend: str = None) -> Dict[str, Any]:
    if not isinstance(resp, dict):
        raise BackendError("响应不是对象")
    if resp.get("protocol") != PROTOCOL:
        raise BackendError("协议版本不匹配：%r" % resp.get("protocol"))
    if resp.get("status") != "ok":
        raise BackendError("后端返回错误：%s" % str(resp.get("error", ""))[:200])
    result = resp.get("result")
    if not isinstance(result, dict) or "frames" not in result or "meta" not in result:
        raise BackendError("响应缺少 result.frames / result.meta")
    if expect_backend and resp.get("backend") not in (None, expect_backend):
        raise BackendError("响应 backend 字段与期望不符")
    return result


def run_backend(request: Dict[str, Any]) -> Dict[str, Any]:
    """按请求分派后端，返回 v1 响应对象。"""
    if not isinstance(request, dict) or request.get("protocol") != PROTOCOL:
        raise BackendError("请求协议必须是 %s" % PROTOCOL)
    name = request.get("backend", "local")
    params, opts = request.get("params", {}), request.get("opts", {})
    if name == "local":
        result = LocalDEMBackend().run(params, opts)
        backend = "local"
    elif name == "external":
        # 外部进程的具体命令行由请求给出（不硬编码任何 GPL/商业路径）
        cmd = request.get("command")
        if not cmd:
            raise BackendError("external 后端必须提供 command（外部脚本 argv）")
        result = ExternalDEMBackend(cmd).run(params, opts)
        backend = "external"
    else:
        raise BackendError("未知后端：%r" % name)
    return dict(protocol=PROTOCOL, backend=backend, status="ok", result=result)


def selfcheck_backend() -> Dict[str, Any]:
    """适配器自检：请求构造、本地后端 v1 往返、错误协议拒绝、外部基类接口存在。"""
    checks: Dict[str, Any] = {}
    req = make_request(dict(n_particles=20),
                       dict(n_steps=20, settle_steps=20, save_every=10, max_particles=50))
    checks["request_protocol"] = (req["protocol"] == PROTOCOL)
    checks["request_shape"] = ("params" in req and "opts" in req and "command" in req)
    resp = run_backend(req)
    checks["response_protocol"] = (resp["protocol"] == PROTOCOL and resp["status"] == "ok")
    checks["response_result"] = ("frames" in resp["result"] and "meta" in resp["result"])
    checks["local_backend_frames"] = (len(resp["result"]["frames"]) >= 2)
    try:
        run_backend(dict(protocol="bad/v9"))
        checks["reject_bad_protocol"] = False
    except BackendError:
        checks["reject_bad_protocol"] = True
    try:
        run_backend(make_request(backend="external"))
        checks["reject_external_without_command"] = False
    except BackendError:
        checks["reject_external_without_command"] = True
    checks["external_class_contract"] = hasattr(ExternalDEMBackend, "run") and hasattr(ExternalDEMBackend, "selfcheck")
    checks["license_boundary_declared"] = ("不链接" in __doc__ and "不随包分发" in __doc__)
    checks["all_pass"] = all(v for k, v in checks.items() if k != "all_pass")
    return checks


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="DEM v1 独立进程后端桥")
    ap.add_argument("--backend", default="local", choices=["local", "external"])
    ap.add_argument("--input", help="请求 JSON 文件（缺省从 stdin 读）")
    ap.add_argument("--output", help="响应 JSON 文件（缺省写到 stdout）")
    ap.add_argument("--selfcheck", action="store_true")
    args = ap.parse_args(argv)

    if args.selfcheck:
        print(json.dumps(selfcheck_backend(), indent=2, ensure_ascii=False))
        return 0
    try:
        raw = open(args.input, "r", encoding="utf-8").read() if args.input else sys.stdin.read()
        req = json.loads(raw)
        if args.backend != "local":
            req["backend"] = args.backend
        resp = run_backend(req)
    except Exception as exc:
        resp = dict(protocol=PROTOCOL, backend=args.backend, status="error",
                    error="%s: %s" % (type(exc).__name__, str(exc)[:200]))
    text = json.dumps(resp, ensure_ascii=False)
    if args.output:
        open(args.output, "w", encoding="utf-8").write(text)
    else:
        print(text)
    return 0 if resp.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(_main())
