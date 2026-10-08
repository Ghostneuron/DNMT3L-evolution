from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
SUMMARIES = PROJECT / "analysis_v2" / "summaries"
CTERMINAL_START = 178
CTERMINAL_END = 380
REGIONS = {
    "switching_helix": (341, 355),
    "switching_helix_adjacent_loop": (356, 360),
}


def bh_adjust(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[order] = np.minimum(adjusted, 1.0)
    return result


def main() -> None:
    sites = pd.read_csv(SUMMARIES / "amino_acid_constraint_by_site.csv")
    stable = sites.loc[
        (sites["Aligners_Retaining_Site"] == 3)
        & sites["Alignment_Stable"].astype(str).str.lower().isin({"true", "1"})
        & sites["Human_Position"].between(CTERMINAL_START, CTERMINAL_END)
    ].copy()

    rows: list[dict[str, object]] = []
    for scope, frame in stable.groupby("Scope"):
        entropy = frame.set_index("Human_Position")["Median_AA_Entropy"]
        for region, (start, end) in REGIONS.items():
            length = end - start + 1
            focal_positions = list(range(start, end + 1))
            focal_complete = all(position in entropy.index for position in focal_positions)
            focal_median = float(entropy.loc[focal_positions].median()) if focal_complete else np.nan

            windows: list[float] = []
            for window_start in range(CTERMINAL_START, CTERMINAL_END - length + 2):
                positions = list(range(window_start, window_start + length))
                if all(position in entropy.index for position in positions):
                    windows.append(float(entropy.loc[positions].median()))

            at_least_as_variable = (
                int(sum(value >= focal_median for value in windows))
                if focal_complete and windows
                else np.nan
            )
            empirical_p = (
                (1 + at_least_as_variable) / (1 + len(windows))
                if focal_complete and windows
                else np.nan
            )
            rows.append(
                {
                    "Scope": scope,
                    "Structural_Region": region,
                    "Region_Start": start,
                    "Region_End": end,
                    "Region_Length": length,
                    "Focal_Window_Complete": focal_complete,
                    "Focal_Median_Entropy": focal_median,
                    "Eligible_Comparable_Windows": len(windows),
                    "Windows_At_Least_As_Variable": at_least_as_variable,
                    "Empirical_Upper_Tail_p": empirical_p,
                }
            )

    result = pd.DataFrame(rows).sort_values(["Scope", "Structural_Region"])
    result["BH_q_across_evaluable_tests"] = bh_adjust(result["Empirical_Upper_Tail_p"])
    result.to_csv(SUMMARIES / "switching_helix_window_sensitivity.csv", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
