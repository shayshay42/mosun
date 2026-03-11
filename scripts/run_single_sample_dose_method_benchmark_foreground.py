#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import threading
from pathlib import Path


OBJECTIVES = ("simple", "clinical", "tracking")


def stream_pipe(prefix: str, pipe, log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8", newline="") as log:
        for line in iter(pipe.readline, ""):
            msg = f"[{prefix}] {line}"
            sys.stdout.write(msg)
            sys.stdout.flush()
            log.write(line)
            log.flush()
    pipe.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument(
        "--out-root",
        default=None,
        help="Benchmark output directory. Defaults to generated/figures/optimization/single_sample_dose_method_benchmark_20260311 under repo root.",
    )
    parser.add_argument("--threads-per-job", type=int, default=4)
    parser.add_argument("--solver", default="rodas4p")
    parser.add_argument("--saveat-dt", default="1.0")
    parser.add_argument("--rs-evals", default="180")
    parser.add_argument("--mc-steps", default="260")
    parser.add_argument("--nm-iters", default="60")
    parser.add_argument("--lbfgs-iters", default="25")
    parser.add_argument("--print-every", default="5")
    parser.add_argument("--objectives", nargs="+", default=list(OBJECTIVES))
    args = parser.parse_args()

    root = Path(args.root).resolve()
    out_root = (
        Path(args.out_root).resolve()
        if args.out_root
        else root / "generated" / "figures" / "optimization" / "single_sample_dose_method_benchmark_20260311"
    )
    out_root.mkdir(parents=True, exist_ok=True)

    procs: list[tuple[str, subprocess.Popen[str]]] = []
    threads: list[threading.Thread] = []

    for obj in args.objectives:
        env = os.environ.copy()
        env["JULIA_NUM_THREADS"] = str(args.threads_per_job)
        env["SINGLE_BENCH_OBJECTIVE"] = obj
        env["SINGLE_BENCH_SOLVER"] = str(args.solver)
        env["SINGLE_BENCH_SAVEAT_DT"] = str(args.saveat_dt)
        env["SINGLE_BENCH_RS_EVALS"] = str(args.rs_evals)
        env["SINGLE_BENCH_MC_STEPS"] = str(args.mc_steps)
        env["SINGLE_BENCH_NM_ITERS"] = str(args.nm_iters)
        env["SINGLE_BENCH_LBFGS_ITERS"] = str(args.lbfgs_iters)
        env["SINGLE_BENCH_PRINT_EVERY"] = str(args.print_every)
        env["SINGLE_BENCH_OUT_ROOT"] = str(out_root)
        cmd = ["julia", "--project=./julia", "julia/benchmark_single_sample_dose_methods_mosun.jl"]
        proc = subprocess.Popen(
            cmd,
            cwd=root,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        procs.append((obj, proc))
        log_path = out_root / f"{obj}.stdout.log"
        t = threading.Thread(target=stream_pipe, args=(obj, proc.stdout, log_path), daemon=True)
        t.start()
        threads.append(t)

    returncodes: dict[str, int] = {}
    for obj, proc in procs:
        returncodes[obj] = proc.wait()

    for t in threads:
        t.join()

    failed = {k: v for k, v in returncodes.items() if v != 0}
    if failed:
        for obj, code in failed.items():
            print(f"[launcher] objective={obj} exited with code={code}", file=sys.stderr)
        return 1

    plot_cmd = [sys.executable, str(root / "scripts" / "plot_single_sample_dose_method_benchmark.py")]
    plot_env = os.environ.copy()
    plot_env["SINGLE_BENCH_RESULTS_ROOT"] = str(out_root)
    plot_proc = subprocess.run(plot_cmd, cwd=root, env=plot_env, text=True)
    return plot_proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
