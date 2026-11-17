from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Seq import Seq


PROJECT = Path(__file__).resolve().parents[1]
ROOT = PROJECT / "analysis" / "selection_reanalysis_337"
RUNS = ROOT / "runs"
CROSSWALK = PROJECT / "results" / "selection_337_human_coordinate_crosswalk.tsv"
FUNCTIONAL = PROJECT / "metadata" / "functional_site_manifest.tsv"
ALL_SITES_OUT = PROJECT / "results" / "selection_337_all_site_summary.tsv"
ANALYSIS_OUT = PROJECT / "results" / "selection_337_analysis_summary.tsv"
FUNCTIONAL_OUT = PROJECT / "results" / "selection_337_functional_site_summary.tsv"
DECISION_OUT = PROJECT / "results" / "selection_337_preliminary_decisions.tsv"
FINAL_DECISION_OUT = PROJECT / "results" / "selection_337_final_decisions.tsv"
MULTIHIT_OUT = PROJECT / "results" / "selection_337_functional_multihit_summary.tsv"
TOPOLOGY_OUT = PROJECT / "results" / "selection_337_topology_sensitivity_summary.tsv"
FUBAR_DECISION_OUT = PROJECT / "results" / "selection_337_fubar_candidate_decisions.tsv"
STATE_OUT = PROJECT / "results" / "selection_337_fubar_candidate_state_counts.tsv"
ANALYSES = ["CONSENSUS", "MACSE", "MAFFT", "PRANK", "MUSCLE5", "COBALT"]
RAW_ALIGNERS = ["MACSE", "MAFFT", "PRANK", "MUSCLE5", "COBALT"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> dict | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    with path.open() as handle:
        return json.load(handle)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def bh(values: list[float]) -> list[float]:
    count = len(values)
    order = sorted(range(count), key=values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[original_index] * count / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def main() -> None:
    crosswalk_rows = read_tsv(CROSSWALK)
    crosswalk: dict[str, dict[int, dict[str, str]]] = defaultdict(dict)
    for row in crosswalk_rows:
        crosswalk[row["analysis_id"]][int(row["filtered_codon_site"])] = row

    all_rows: list[dict[str, object]] = []
    analysis_rows: list[dict[str, object]] = []
    by_analysis_position: dict[tuple[str, int], dict[str, object]] = {}

    for analysis in ANALYSES:
        fubar_path = RUNS / f"{analysis}_FUBAR" / f"{analysis}_FUBAR.json"
        meme_path = RUNS / f"{analysis}_MEME_STANDARD" / f"{analysis}_MEME_STANDARD.json"
        fubar = load_json(fubar_path)
        meme = load_json(meme_path)
        fubar_rows = fubar["MLE"]["content"]["0"] if fubar else []
        meme_rows = meme["MLE"]["content"]["0"] if meme else []
        if fubar_rows and len(fubar_rows) != len(crosswalk[analysis]):
            raise ValueError(f"{analysis}: FUBAR/crosswalk length mismatch")
        if meme_rows and len(meme_rows) != len(crosswalk[analysis]):
            raise ValueError(f"{analysis}: MEME/crosswalk length mismatch")
        meme_q = bh([float(row[6]) for row in meme_rows]) if meme_rows else []

        for site, mapping in sorted(crosswalk[analysis].items()):
            row: dict[str, object] = {
                "analysis_id": analysis,
                "filtered_codon_site": site,
                "human_position": int(mapping["human_position"]),
                "source_alignment_codon_column": int(mapping["source_alignment_codon_column"]),
                "pre_filter_occupancy": float(mapping["pre_filter_occupancy"]),
                "post_filter_occupancy": float(mapping["post_filter_occupancy"]),
                "FUBAR_alpha": "",
                "FUBAR_beta": "",
                "FUBAR_posterior_purifying": "",
                "FUBAR_posterior_diversifying": "",
                "FUBAR_bayes_factor_diversifying": "",
                "MEME_alpha": "",
                "MEME_beta_positive": "",
                "MEME_positive_weight": "",
                "MEME_LRT": "",
                "MEME_p": "",
                "MEME_BH_q": "",
                "MEME_approx_branches_EBF_ge_100": "",
            }
            if fubar_rows:
                frow = fubar_rows[site - 1]
                row.update(
                    {
                        "FUBAR_alpha": frow[0],
                        "FUBAR_beta": frow[1],
                        "FUBAR_posterior_purifying": frow[3],
                        "FUBAR_posterior_diversifying": frow[4],
                        "FUBAR_bayes_factor_diversifying": frow[5],
                    }
                )
            if meme_rows:
                mrow = meme_rows[site - 1]
                row.update(
                    {
                        "MEME_alpha": mrow[0],
                        "MEME_beta_positive": mrow[3],
                        "MEME_positive_weight": mrow[4],
                        "MEME_LRT": mrow[5],
                        "MEME_p": mrow[6],
                        "MEME_BH_q": meme_q[site - 1],
                        "MEME_approx_branches_EBF_ge_100": mrow[7],
                    }
                )
            all_rows.append(row)
            by_analysis_position[(analysis, int(mapping["human_position"]))] = row

        analysis_rows.append(
            {
                "analysis_id": analysis,
                "retained_sites": len(crosswalk[analysis]),
                "FUBAR_complete": bool(fubar_rows),
                "FUBAR_diversifying_PP_ge_0.90": sum(float(row[4]) >= 0.90 for row in fubar_rows),
                "FUBAR_purifying_PP_ge_0.90": sum(float(row[3]) >= 0.90 for row in fubar_rows),
                "MEME_complete": bool(meme_rows),
                "MEME_nominal_p_lt_0.05": sum(float(row[6]) < 0.05 for row in meme_rows),
                "MEME_BH_q_lt_0.05": sum(value < 0.05 for value in meme_q),
            }
        )

    all_fields = [
        "analysis_id", "filtered_codon_site", "human_position", "source_alignment_codon_column",
        "pre_filter_occupancy", "post_filter_occupancy", "FUBAR_alpha", "FUBAR_beta",
        "FUBAR_posterior_purifying", "FUBAR_posterior_diversifying",
        "FUBAR_bayes_factor_diversifying", "MEME_alpha", "MEME_beta_positive",
        "MEME_positive_weight", "MEME_LRT", "MEME_p", "MEME_BH_q",
        "MEME_approx_branches_EBF_ge_100",
    ]
    write_tsv(ALL_SITES_OUT, all_rows, all_fields)
    write_tsv(ANALYSIS_OUT, analysis_rows, list(analysis_rows[0]))

    feature_rows = [row for row in read_tsv(FUNCTIONAL) if row["feature_type"] == "site"]
    features_by_position: dict[int, list[str]] = defaultdict(list)
    for row in feature_rows:
        features_by_position[int(row["human_np_position"])].append(row["feature_id"])
    functional_rows: list[dict[str, object]] = []
    for position in sorted(features_by_position):
        for analysis in ANALYSES:
            result = by_analysis_position.get((analysis, position), {})
            functional_rows.append(
                {
                    "human_position": position,
                    "feature_ids": ";".join(features_by_position[position]),
                    "analysis_id": analysis,
                    "retained": bool(result),
                    "filtered_codon_site": result.get("filtered_codon_site", ""),
                    "post_filter_occupancy": result.get("post_filter_occupancy", ""),
                    "FUBAR_posterior_purifying": result.get("FUBAR_posterior_purifying", ""),
                    "FUBAR_posterior_diversifying": result.get("FUBAR_posterior_diversifying", ""),
                    "MEME_p": result.get("MEME_p", ""),
                    "MEME_BH_q": result.get("MEME_BH_q", ""),
                }
            )
    write_tsv(FUNCTIONAL_OUT, functional_rows, list(functional_rows[0]))

    multihit_path = (
        RUNS
        / "CONSENSUS_MEME_MULTIHIT_FUNCTIONAL"
        / "CONSENSUS_MEME_MULTIHIT_FUNCTIONAL.json"
    )
    multihit = load_json(multihit_path)
    multihit_rows: list[dict[str, object]] = []
    if multihit:
        content = multihit["MLE"]["content"]["0"]
        for position in sorted(features_by_position):
            primary = by_analysis_position.get(("CONSENSUS", position), {})
            if not primary:
                continue
            site = int(primary["filtered_codon_site"])
            row = content[site - 1]
            multihit_rows.append(
                {
                    "human_position": position,
                    "feature_ids": ";".join(features_by_position[position]),
                    "filtered_codon_site": site,
                    "standard_MEME_p": primary.get("MEME_p", ""),
                    "standard_MEME_BH_q": primary.get("MEME_BH_q", ""),
                    "DoubleTriple_MEME_LRT": row[5],
                    "DoubleTriple_MEME_p": row[6],
                    "DoubleTriple_approx_branches_EBF_ge_100": row[7],
                    "relative_rate_two_nucleotide_hits": row[13],
                    "relative_rate_three_nucleotide_hits": row[14],
                }
            )
    if multihit_rows:
        write_tsv(MULTIHIT_OUT, multihit_rows, list(multihit_rows[0]))

    taxonomy_fubar = load_json(
        RUNS / "CONSENSUS_FUBAR_NCBI_TAXONOMY" / "CONSENSUS_FUBAR_NCBI_TAXONOMY.json"
    )
    taxonomy_meme = load_json(
        RUNS
        / "CONSENSUS_MEME_STANDARD_FUNCTIONAL_NCBI_TAXONOMY"
        / "CONSENSUS_MEME_STANDARD_FUNCTIONAL_NCBI_TAXONOMY.json"
    )
    topology_rows: list[dict[str, object]] = []
    if taxonomy_fubar:
        fubar_content = taxonomy_fubar["MLE"]["content"]["0"]
        meme_content = taxonomy_meme["MLE"]["content"]["0"] if taxonomy_meme else []
        for site, mapping in sorted(crosswalk["CONSENSUS"].items()):
            position = int(mapping["human_position"])
            primary = by_analysis_position[("CONSENSUS", position)]
            topology_rows.append(
                {
                    "human_position": position,
                    "filtered_codon_site": site,
                    "prespecified_functional_site": position in features_by_position,
                    "gene_tree_FUBAR_posterior_diversifying": primary["FUBAR_posterior_diversifying"],
                    "ncbi_taxonomy_FUBAR_posterior_diversifying": fubar_content[site - 1][4],
                    "gene_tree_MEME_p": primary["MEME_p"],
                    "gene_tree_MEME_BH_q": primary["MEME_BH_q"],
                    "ncbi_taxonomy_focal_MEME_p": (
                        meme_content[site - 1][6]
                        if meme_content and position in features_by_position
                        else ""
                    ),
                }
            )
        write_tsv(TOPOLOGY_OUT, topology_rows, list(topology_rows[0]))

        candidate_rows: list[dict[str, object]] = []
        for row in topology_rows:
            position = int(row["human_position"])
            primary_pp = float(row["gene_tree_FUBAR_posterior_diversifying"])
            taxonomy_pp = float(row["ncbi_taxonomy_FUBAR_posterior_diversifying"])
            raw_values = [
                float(by_analysis_position[(analysis, position)]["FUBAR_posterior_diversifying"])
                for analysis in RAW_ALIGNERS
                if by_analysis_position.get((analysis, position), {}).get("FUBAR_posterior_diversifying", "") != ""
            ]
            raw_direction = sum(value > 0.5 for value in raw_values)
            raw_threshold = sum(value >= 0.90 for value in raw_values)
            if max([primary_pp, taxonomy_pp, *raw_values]) < 0.90:
                continue
            if primary_pp >= 0.90 and raw_direction >= 4 and taxonomy_pp >= 0.90:
                decision = "PERVASIVE_SUPPORT_ROBUST_TO_ALIGNMENT_AND_TOPOLOGY"
            elif primary_pp >= 0.90 and raw_direction >= 4 and taxonomy_pp < 0.90:
                decision = "PRIMARY_PERVASIVE_SUPPORT_BUT_TOPOLOGY_SENSITIVE"
            elif primary_pp < 0.90:
                decision = "SENSITIVITY_ONLY_NOT_PRIMARY_PERVASIVE_SUPPORT"
            else:
                decision = "ALIGNMENT_SENSITIVE_PERVASIVE_SUPPORT"
            candidate_rows.append(
                {
                    "human_position": position,
                    "prespecified_functional_site": position in features_by_position,
                    "gene_tree_consensus_posterior": primary_pp,
                    "ncbi_taxonomy_posterior": taxonomy_pp,
                    "raw_aligners_same_direction_count": raw_direction,
                    "raw_aligners_posterior_ge_0.90_count": raw_threshold,
                    "decision": decision,
                }
            )
        write_tsv(FUBAR_DECISION_OUT, candidate_rows, list(candidate_rows[0]))

        consensus_fasta = read_fasta(
            ROOT / "inputs" / "Placentalia_CONSENSUS_human_coordinate_occupancy80_deduplicated.fasta"
        )
        state_rows: list[dict[str, object]] = []
        for candidate in candidate_rows:
            position = int(candidate["human_position"])
            site = int(by_analysis_position[("CONSENSUS", position)]["filtered_codon_site"])
            codon_counts: Counter[str] = Counter()
            amino_acid_counts: Counter[str] = Counter()
            for sequence in consensus_fasta.values():
                codon = sequence[(site - 1) * 3 : site * 3]
                codon_counts[codon] += 1
                amino_acid = str(Seq(codon).translate()) if set(codon) <= set("ACGT") else "X"
                amino_acid_counts[amino_acid] += 1
            human_codon = consensus_fasta["Homo_sapiens"][(site - 1) * 3 : site * 3]
            state_rows.append(
                {
                    "human_position": position,
                    "filtered_codon_site": site,
                    "human_codon": human_codon,
                    "human_amino_acid": str(Seq(human_codon).translate()),
                    "resolved_sequence_count": sum(
                        count for state, count in amino_acid_counts.items() if state != "X"
                    ),
                    "unresolved_sequence_count": amino_acid_counts.get("X", 0),
                    "amino_acid_counts": ";".join(
                        f"{state}:{count}" for state, count in sorted(amino_acid_counts.items())
                    ),
                    "codon_counts": ";".join(
                        f"{state}:{count}" for state, count in sorted(codon_counts.items())
                    ),
                    "decision": candidate["decision"],
                }
            )
        write_tsv(STATE_OUT, state_rows, list(state_rows[0]))

    decision_positions = sorted({int(row["human_position"]) for row in all_rows})
    decision_rows: list[dict[str, object]] = []
    all_sensitivities_complete = all(
        by_analysis_position.get((analysis, position), {}).get("MEME_p", "") != ""
        for analysis in RAW_ALIGNERS
        for position in decision_positions
        if (analysis, position) in by_analysis_position
    ) and all((RUNS / f"{analysis}_MEME_STANDARD" / f"{analysis}_MEME_STANDARD.json").exists() for analysis in RAW_ALIGNERS)
    for position in decision_positions:
        primary = by_analysis_position.get(("CONSENSUS", position), {})
        primary_q = primary.get("MEME_BH_q", "")
        nominal_count = sum(
            float(by_analysis_position[(analysis, position)]["MEME_p"]) < 0.05
            for analysis in RAW_ALIGNERS
            if by_analysis_position.get((analysis, position), {}).get("MEME_p", "") != ""
        )
        corrected_aligners = [
            analysis
            for analysis in RAW_ALIGNERS
            if by_analysis_position.get((analysis, position), {}).get("MEME_BH_q", "") != ""
            and float(by_analysis_position[(analysis, position)]["MEME_BH_q"]) < 0.05
        ]
        fubar_direction_count = sum(
            float(by_analysis_position[(analysis, position)]["FUBAR_posterior_diversifying"]) >= 0.5
            for analysis in RAW_ALIGNERS
            if by_analysis_position.get((analysis, position), {}).get("FUBAR_posterior_diversifying", "") != ""
        )
        if primary_q == "":
            label = "PENDING_PRIMARY_MEME"
        elif float(primary_q) >= 0.05:
            label = "NO_CORRECTED_EPISODIC_SUPPORT"
        elif not all_sensitivities_complete:
            label = "PRIMARY_Q_LT_0.05; ALIGNER/MULTIHIT/TOPOLOGY_SENSITIVITIES_PENDING"
        elif nominal_count >= 4:
            label = "PRIMARY_AND_ALIGNER_CRITERIA_MET; MULTIHIT/TOPOLOGY_PENDING"
        else:
            label = "PRIMARY_ONLY_OR_MODEL_SENSITIVE"
        decision_rows.append(
            {
                "human_position": position,
                "prespecified_functional_site": position in features_by_position,
                "primary_MEME_p": primary.get("MEME_p", ""),
                "primary_MEME_BH_q": primary_q,
                "raw_aligners_nominal_p_lt_0.05": nominal_count,
                "raw_aligners_BH_q_lt_0.05": len(corrected_aligners),
                "raw_aligners_BH_q_lt_0.05_names": ";".join(corrected_aligners),
                "primary_FUBAR_diversifying_posterior": primary.get("FUBAR_posterior_diversifying", ""),
                "raw_aligners_same_FUBAR_direction": fubar_direction_count,
                "preliminary_decision": label,
            }
        )
    write_tsv(DECISION_OUT, decision_rows, list(decision_rows[0]))
    write_tsv(FINAL_DECISION_OUT, decision_rows, list(decision_rows[0]))

    print(f"Wrote {len(all_rows)} analysis-site rows")
    for row in analysis_rows:
        print(
            f"{row['analysis_id']}: FUBAR={row['FUBAR_complete']}, MEME={row['MEME_complete']}, "
            f"MEME q<0.05={row['MEME_BH_q_lt_0.05']}"
        )


if __name__ == "__main__":
    main()
