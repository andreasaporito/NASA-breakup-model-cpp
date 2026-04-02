#!/usr/bin/env python3
"""
benchmark.py — Run breakup model simulations across a sweep of
minimalCharacteristicLength values and record wall-clock time.

Usage
-----
# Benchmark the current branch (no git checkout, no rebuild):
python benchmark.py \
    --repo   /home/andrea/LSMS_project/NASA-breakup-model-cpp \
    --data   /home/andrea/LSMS_project/NASA-breakup-model-cpp/example_config/CTwrtMinCT/data.yaml \
    --output results_branch_A.csv \
    --label  "branch_A"

# Benchmark a specific branch (will checkout + cmake + make):
python benchmark.py \
    --repo    /home/andrea/LSMS_project/NASA-breakup-model-cpp \
    --data    /home/andrea/LSMS_project/NASA-breakup-model-cpp/example_config/CTwrtMinCT/data.yaml \
    --output  results_branch_B.csv \
    --label   "branch_B" \
    --branch  my-other-branch \
    --rebuild          # pass this flag to trigger cmake + make
"""
from __future__ import annotations
import argparse
import csv
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional

import yaml  # pip install pyyaml


# ---------------------------------------------------------------------------
# Sweep values for minimalCharacteristicLength
# ---------------------------------------------------------------------------
MIN_CL_VALUES = [
    0.0010,
    0.0025,
    0.0050,
    0.0075,
    0.0100,
    0.0125,
    0.0150,
    0.0175,
    0.0200,
    0.0300,
    0.0400,
    0.0500,
    0.0750,
    0.1000,
    0.2000
]

# How many times to repeat each simulation (results are averaged).
REPEATS = 3

# Regex patterns to extract info from spdlog output
_RE_TIME = re.compile(r"The simulation took\s+(\d+)\s+ms")
_RE_FRAGS = re.compile(r"The simulation produced\s+(\d+)\s+fragments")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def git_checkout(repo: Path, branch: str) -> None:
    print(f"\n[git] Checking out branch '{branch}' …")
    subprocess.run(["git", "checkout", branch], cwd=repo, check=True)


def cmake_build(repo: Path, jobs: int = 4) -> None:
    build_dir = repo / "build"
    build_dir.mkdir(exist_ok=True)
    print(f"\n[cmake] Configuring in {build_dir} …")
    subprocess.run(["cmake", ".."], cwd=build_dir, check=True)
    print(f"[make]  Building with -j{jobs} …")
    subprocess.run(["make", f"-j{jobs}"], cwd=build_dir, check=True)


def write_config(tmpdir: Path, min_cl: float, data_yaml: Path) -> Path:
    """Write a minimal config.yaml for one simulation run."""
    cfg = {
        "simulation": {
            "minimalCharacteristicLength": float(min_cl),
            "simulationType": "COLLISION",
            "inputSource": [str(data_yaml.resolve())],
            "enforceMassConservation": True,
        },
        # Route all output to throwaway files so I/O doesn't pollute timing.
        "inputOutput": {
            "target": [str(tmpdir / "input.csv")],
            "kepler": False,
        },
        "resultOutput": {
            "target": [str(tmpdir / "output.csv")],
            "kepler": False,
        },
    }
    path = tmpdir / f"config_{min_cl:.6f}.yaml"
    with open(path, "w") as fh:
        yaml.dump(cfg, fh, default_flow_style=False)
    return path


def run_once(executable: Path, config: Path) -> dict:
    """
    Run the executable once and return a dict with
    'ms' (int | None) and 'fragments' (int | None).
    """
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(executable), str(config)],
        capture_output=True,
        text=True,
    )
    wall_s = time.perf_counter() - t0

    output = proc.stdout + proc.stderr

    ms_match = _RE_TIME.search(output)
    fr_match = _RE_FRAGS.search(output)

    if proc.returncode != 0 or ms_match is None:
        print(f"  [!] Run failed or produced no timing. stderr:\n{proc.stderr[:400]}")
        return {"ms": None, "fragments": None, "wall_s": wall_s}

    return {
        "ms": int(ms_match.group(1)),
        "fragments": int(fr_match.group(1)) if fr_match else None,
        "wall_s": wall_s,
    }


# ---------------------------------------------------------------------------
# Main benchmark loop
# ---------------------------------------------------------------------------

def benchmark(
    repo: Path,
    data_yaml: Path,
    output_csv: Path,
    label: str,
    # branch can be a string or None. If None, we benchmark the current branch without git checkout.
    # TypeError: unsupported operand type(s) for |: 'type' and 'NoneType'
    branch: Optional[str] = None,
    rebuild: bool = False,
    repeats: int = REPEATS,
    min_cl_values: list[float] = MIN_CL_VALUES,
) -> None:

    if branch:
        git_checkout(repo, branch)

    if rebuild:
        cmake_build(repo)

    executable = repo / "build" / "breakupModel"
    if not executable.exists():
        sys.exit(
            f"[error] Executable not found at {executable}.\n"
            "        Build the project first, or pass --rebuild."
        )

    rows: list[dict] = []

    with tempfile.TemporaryDirectory(prefix="breakup_bench_") as _tmp:
        tmpdir = Path(_tmp)

        for min_cl in min_cl_values:
            config_path = write_config(tmpdir, min_cl, data_yaml)
            print(f"\n  minCL = {min_cl:.4f}  ({repeats} run(s))")

            times_ms: list[int] = []
            frags_val: int | None = None

            for rep in range(repeats):
                result = run_once(executable, config_path)
                if result["ms"] is not None:
                    times_ms.append(result["ms"])
                    frags_val = result["fragments"]
                    print(f"    run {rep + 1}: {result['ms']} ms  "
                          f"({result['fragments']} fragments)")
                else:
                    print(f"    run {rep + 1}: FAILED")

            if not times_ms:
                print(f"  [!] All runs failed for minCL={min_cl}. Skipping.")
                continue

            rows.append({
                "label":     label,
                "min_cl":    min_cl,
                "mean_ms":   sum(times_ms) / len(times_ms),
                "min_ms":    min(times_ms),
                "max_ms":    max(times_ms),
                "std_ms":    _std(times_ms),
                "fragments": frags_val,
                "n_runs":    len(times_ms),
            })

    # Write CSV
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["label", "min_cl", "mean_ms", "min_ms", "max_ms",
                  "std_ms", "fragments", "n_runs"]
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n✓ Results saved to {output_csv}")
    _print_table(rows)


def _std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return (sum((v - mean) ** 2 for v in values) / (len(values) - 1)) ** 0.5


def _print_table(rows: list[dict]) -> None:
    print(f"\n{'minCL':>10}  {'mean_ms':>10}  {'min_ms':>8}  "
          f"{'max_ms':>8}  {'fragments':>10}")
    print("-" * 56)
    for r in rows:
        print(f"{r['min_cl']:>10.4f}  {r['mean_ms']:>10.1f}  "
              f"{r['min_ms']:>8}  {r['max_ms']:>8}  "
              f"{str(r['fragments']):>10}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Benchmark breakupModel across minimalCharacteristicLength values."
    )
    p.add_argument("--repo",    required=True,  type=Path,
                   help="Root of the git repository.")
    p.add_argument("--data",    required=True,  type=Path,
                   help="Path to data.yaml (satellite inputs).")
    p.add_argument("--output",  required=True,  type=Path,
                   help="Output CSV file path, e.g. results_main.csv")
    p.add_argument("--label",   default="default",
                   help="Label for this run (appears in CSV and plots).")
    p.add_argument("--branch",  default=None,
                   help="Git branch to checkout before benchmarking.")
    p.add_argument("--rebuild", action="store_true",
                   help="Run cmake + make before benchmarking.")
    p.add_argument("--repeats", type=int, default=REPEATS,
                   help=f"Number of repetitions per minCL value (default: {REPEATS}).")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    benchmark(
        repo=args.repo,
        data_yaml=args.data,
        output_csv=args.output,
        label=args.label,
        branch=args.branch,
        rebuild=args.rebuild,
        repeats=args.repeats,
    )
