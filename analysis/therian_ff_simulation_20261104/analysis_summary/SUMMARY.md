# Phase-1 molecular-dynamics summary

## Validation

- Six canonical replicas passed the completion and analysis gates (3 Y/E/E and 3 F/T/D).
- Each replica contributed 500 production frames from a 5.0-ns trajectory; the aggregate sampled production time was 30.0 ns.
- Each Y/E/E trajectory contained 142,702 atoms and each F/T/D trajectory contained 142,689 atoms.
- The preserved duplicate-corrupted directory was not admitted.

## Main findings

- The B261-A728 aromatic contact was present in 100% of frames in every replica. Replica-mean distances were 0.343 +/- 0.004 nm for Y/E/E and 0.337 +/- 0.009 nm for F/T/D (mean +/- SD across three replica means). The paired direction was mixed (2/3 lower in F/T/D).
- B291-A725 and B293-A725 distances varied among replicas, without a common direction across all three pairs. Sub-0.5-nm contacts at these sites were infrequent and replica-specific.
- B293-A722 was shorter in F/T/D in all three matched replica indices, but the mean difference was small (-0.050 nm) relative to between-replica variation, and neither endpoint formed a persistent sub-0.5-nm contact.
- Protein-backbone RMSD averaged 0.213 +/- 0.004 nm for Y/E/E and 0.201 +/- 0.014 nm for F/T/D.
- Mean DNMT3A active-loop CA RMSF overlapped (0.514 +/- 0.095 nm versus 0.449 +/- 0.054 nm), with a mixed paired direction.

## Interpretation boundary

The six short trajectories support stable maintenance of the common B261-A728 aromatic contact and do not show a replica-consistent global-stability or active-loop-mobility shift between the reconstructed states. The B293-A722 distance tendency is exploratory. These simulations do not estimate binding free energy, catalytic activity, imprinting, organismal fitness, or adaptation. Frames are repeated observations within trajectories; the biological comparison uses three independent replica summaries per endpoint.
