# HyPhy focal site filter audit

Date: 6 October 2026

## Scope

This note records a pre-output audit of the focal-site filter used by `scripts/run_conditional_sensitivities.py` for the conditional MEME and FEL Double+Triple analyses. No conditional sensitivity result had been generated when the issue was identified.

The MEME focal set consists of strict alignment-robust episodic candidates at `q <= 0.05`; the FEL set consists of strict alignment-robust pervasive candidates at the same threshold. Historical h358 is included in every scope for both models as a diagnostic audit case. Conditional results can qualify but cannot create a primary candidate.

## Finding

In the installed HyPhy 2.5.94 code, `selection.io.sitelist_matches_pattern` tests the user-supplied `limit-to-sites` string with a comma-prefixed integer expression. A direct test of that installed helper showed that:

- a range such as `1-2,5-6` did not select the intended sites;
- the first item in an ordinary string such as `1,2,5,6` was not recognized; and
- an explicit leading comma, as in `,1,2,5,6`, ensured that every intended entry was selected; and
- because the installed regular expression does not require a boundary after the index, a requested index can also trigger decimal prefixes (for example, requesting `354` also evaluates `3` and `35`).

The earlier unexecuted script prepared merged ranges and was therefore unsuitable for this installed HyPhy version.

## Correction

The script now supplies a sorted, explicit, leading-comma list of one-based filtered-site indices. It does not use ranges. Exact requested rows, rather than every row that HyPhy happened to evaluate, are extracted for the focal summary and multiplicity correction. After each run, it requires all of the following before a focal result can be written:

1. the JSON-recorded site filter exactly matches the requested string;
2. the JSON contains one row per filtered alignment site;
3. the JSON row length matches the reported header length; and
4. every requested row is classified either as evaluated or as a zero-information row with `p = 1`; zero-information targets remain in the correction family as conservative non-supporting tests; and
5. every additional evaluated row is a decimal prefix of a requested index, with the complete evaluated set, prefix-expanded set, and unevaluable target set recorded in the summary table.

Header names, rather than fixed column offsets, identify each p value and the estimated double- and triple-nucleotide relative rates.

## Consequence

The missing-first-item problem was found before the reported conditional models were run. The decimal-prefix behavior was identified during release inspection of the completed JSON rows. It causes extra site models to be computed but does not change which prespecified rows are extracted, corrected, or interpreted. A requested constant or otherwise zero-information site is retained with `p = 1` instead of being dropped. The release validator independently reconstructs the evaluated rows from each JSON file, verifies every requested row or conservative zero-information default, and rejects any extra row that is not a decimal prefix of a requested index.
