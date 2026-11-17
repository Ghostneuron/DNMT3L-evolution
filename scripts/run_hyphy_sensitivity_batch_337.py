from __future__ import annotations

import argparse
import subprocess
import sys
import time
from collections import deque
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RUNNER = PROJECT / "scripts" / "run_hyphy_selection_337.py"
LOG_DIR = PROJECT / "logs" / "selection_reanalysis_337" / "batch_launchers"
ALIGNERS = ["MACSE", "MAFFT", "PRANK", "MUSCLE5", "COBALT"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--parallel", type=int, default=2)
    parser.add_argument("--threads-per-job", type=int, default=4)
    parser.add_argument(
        "--models",
        nargs="+",
        choices=["FUBAR", "MEME_STANDARD"],
        default=["FUBAR", "MEME_STANDARD"],
    )
    args = parser.parse_args()
    if args.parallel < 1 or args.parallel > 3:
        raise ValueError("--parallel must be between 1 and 3")
    if args.parallel * args.threads_per_job > 8:
        raise ValueError("The batch safety cap is eight assigned threads")

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    queue = deque((analysis, model) for model in args.models for analysis in ALIGNERS)
    active: dict[subprocess.Popen, tuple[str, str, object]] = {}
    failures: list[str] = []

    while queue or active:
        while queue and len(active) < args.parallel:
            analysis, model = queue.popleft()
            label = f"{analysis}_{model}"
            launcher_log = (LOG_DIR / f"{label}.launcher.log").open("w")
            command = [
                sys.executable,
                str(RUNNER),
                "--analysis",
                analysis,
                "--model",
                model,
                "--threads",
                str(args.threads_per_job),
            ]
            process = subprocess.Popen(
                command,
                cwd=PROJECT,
                stdout=launcher_log,
                stderr=subprocess.STDOUT,
                text=True,
            )
            active[process] = (analysis, model, launcher_log)
            print(f"STARTED {label} launcher_PID={process.pid}", flush=True)

        time.sleep(5)
        for process, (analysis, model, launcher_log) in list(active.items()):
            exit_code = process.poll()
            if exit_code is None:
                continue
            launcher_log.close()
            label = f"{analysis}_{model}"
            print(f"FINISHED {label} exit={exit_code}", flush=True)
            if exit_code != 0:
                failures.append(label)
            del active[process]

    if failures:
        print(f"FAILED: {', '.join(failures)}", file=sys.stderr)
        sys.exit(1)
    print("All requested sensitivity jobs completed", flush=True)


if __name__ == "__main__":
    main()
