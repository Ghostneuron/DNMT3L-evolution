# Run status

- Analysis ID: `DNMT3L_FF_THERIAN_MD_20261104_V1`
- Status: `PHASE1_COMPLETE_AND_ANALYZED`
- Started: 2026-11-04 America/Los_Angeles
- Working directory: this directory
- Source structure: `inputs/2QRV.pdb`
- Source SHA-256: `cebae059fc4045e72d97deeb1137bba0ac5e644ec74033b675daa208c7ba8150`
- Primary variants: `amniote_YEE`, `therian_FTD`
- Production replicas: three 5-ns replicas per endpoint completed in
  `production_hmr/`; all six reached step 1,275,000
- Completed outputs: matched 2-fs and 4-fs smoke tests plus six validated
  scientific production trajectories (30 ns aggregate production time)
- Environment: `${HOME}/dnmt3l_ff_openmm_20261109`
- OpenMM: `8.6.1`
- Available platforms: Reference, CPU, OpenCL
- Validated platform: CPU; the advertised OpenCL platform reports no compatible
  device on this Apple M4 system, including outside the sandbox.
- Preserved failed attempt: `.conda/` is an incomplete environment from a conda
  transaction stopped after external-volume linking remained unsuitable for
  execution; it must not be used.
- Current action: production is complete. All six trajectories have been
  analyzed with the prespecified contact, RMSD, and active-loop RMSF readouts.
- Authoritative live status: `production_hmr/BATCH_STATUS.json`
- No simulation or analysis process is expected to remain active.

## Recovery after interrupted session (2026-11-06)

- Four checkpoint files at step 415,000 were preserved in their replica
  directories and backed up under `~/dnmt3l_checkpoint_backup/`.
- The batch resumed at 2026-11-06 14:47 America/Los_Angeles using the original
  `${TMPDIR}/dnmt3l_ff_openmm_20261104/bin/python` environment. The resumed
  driver PID was 49198 and the four stage-1 replica PIDs were 49199-49202.
- The controller writes stdout/stderr to `~/dnmt3l_batch_resume.log`; each
  replica continues to append to its existing `driver.log`, `state.tsv`, and
  `trajectory.dcd`.
- Before the first resumed report, every replica had 166 state-data rows and
  156 DCD production frames, consistent with 10 equilibration reports plus 156
  production reports through step 415,000.
- At 2026-11-06 15:01 America/Los_Angeles, all four replicas had reached step
  417,500. This verifies successful checkpoint loading and append-safe
  continuation beyond the interruption point without duplicate steps.
- The batch controller still owns all four running replicas and retains the
  two stage-2 jobs (seeds 1103 and 2103) as pending; it will launch that matched
  pair automatically after stage 1 completes.

### Runtime persistence repair (2026-11-09)

- macOS temporary-file cleanup removed Python standard-library and package
  files from `${TMPDIR}/dnmt3l_ff_openmm_20261104` while the already-loaded
  stage-1 processes continued to run. Without repair, new stage-2 processes
  would have failed during Python initialization.
- A persistent replacement environment was created at
  `${HOME}/dnmt3l_ff_openmm_20261109` from the local conda package cache,
  using Python 3.11.18 and OpenMM 8.6.1 with the prespecified dependencies.
- All four backed-up step-415,000 checkpoints loaded successfully in the
  replacement environment with matching potential energies.
- The original executable was preserved as
  `${TMPDIR}/dnmt3l_ff_openmm_20261104/bin/python3.11.cleaned_20261109`, and
  the original controller path now symlinks to the persistent Python 3.11.18
  executable. This preserves the controller's frozen stage-2 command without
  altering simulation parameters, seeds, inputs, or outputs.
- Post-repair checks confirmed that the original path imports OpenMM 8.6.1,
  exposes the CPU platform, executes `run_replica.py --help`, and leaves all
  four stage-1 processes active under controller PID 49198.

### Stage transition (2026-11-10)

- All four stage-1 replicas completed successfully at step 1,275,000 with
  `COMPLETE.json` records and exit code 0.
- The first stage-2 launch attempt exited before simulation initialization
  because the detached controller passed an invalid inherited stdin file
  descriptor to its children. No stage-2 checkpoint, trajectory, or state data
  was created by that failed attempt.
- `run_phase1_batch.py` now passes `stdin=subprocess.DEVNULL` explicitly to
  every replica process. This is an operational recovery fix only; commands,
  seeds, thread counts, timesteps, and scientific parameters are unchanged.
- The controller was relaunched as macOS launchd job
  `org.openai.codex.dnmt3l.phase2.20261110`, using the persistent environment
  at `${HOME}/dnmt3l_ff_openmm_20261109`. It skipped all four completed
  stage-1 jobs and launched only seeds 1103 and 2103.
- At 2026-11-10 08:42 America/Los_Angeles, both stage-2 replicas had reached
  step 2,500 and written their first state/checkpoint records, verifying that
  the repaired launch is functioning independently of Codex.

### Stage-2 seed-1103 checkpoint recovery (2026-11-10)

- The amniote seed-1103 worker exited with signal 11 after its last durable
  report at step 37,500. The log contained no OpenMM exception or NaN; the
  therian seed-2103 worker continued normally.
- The step-37,500 checkpoint was preserved as
  `checkpoint.step37500.pre_resume.chk` and loaded successfully in the
  persistent OpenMM 8.6.1 environment. Its 15 state rows and five production
  DCD frames were consistent with the checkpoint boundary.
- Seed 1103 was relaunched alone from that checkpoint under launchd job
  `org.openai.codex.dnmt3l.seed1103.resume.20261110`, using the unchanged seed,
  five CPU threads, 4-fs timestep, inputs, output directory, and reporting
  intervals. Seed 2103 was not restarted or duplicated.

### Stage-2 duplicate-writer repair (2026-11-11)

- `launchctl submit` inferred keep-alive behavior for both the controller and
  the standalone seed-1103 job. After abnormal exits, it relaunched the
  controller and temporarily created two concurrent seed-1103 writers. The
  keep-alive controller and its duplicate child were stopped; seed 2103 was
  not running at the time of repair.
- The complete affected seed-1103 directory is preserved as
  `production_hmr/amniote_YEE_seed1103_duplicate_corrupted_20261111/`.
- The canonical seed-1103 directory was rebuilt from its readable step-285,000
  checkpoint by retaining the first state row and corresponding DCD frame for
  each reported step. This removed exactly 52 duplicate state rows and 52
  duplicate DCD frames. `RECOVERY_PROVENANCE.json` records the transformation;
  `recover_duplicate_dcd.py` implements the guarded reconstruction.
- After repair, seed 1103 had 114 unique state rows and 104 production frames,
  exactly matching step 285,000. Seed 2103 remained clean with 85 unique state
  rows and 75 production frames at step 212,500. Both checkpoints loaded in
  OpenMM 8.6.1 with finite energies.
- Both replicas now run as explicit one-shot LaunchAgents with `KeepAlive`
  false: `org.openai.codex.dnmt3l.seed1103.oneshot.20261111` and
  `org.openai.codex.dnmt3l.seed2103.oneshot.20261111`. No batch controller is
  active, so it cannot create duplicate replicas after a failure.
- At 2026-11-11 08:52 America/Los_Angeles, seed 1103 had advanced cleanly to
  step 287,500 and seed 2103 to step 215,000. State rows were unique and DCD
  frame counts matched both checkpoint boundaries.

### Recurrent five-thread OpenMM CPU abort (2026-11-11)

- The one-shot seed-1103 worker stopped at its durable step-312,500 checkpoint
  with `SIGABRT` while seed 2103 continued normally.
- macOS diagnostic report
  `python3.11-2026-11-11-104521.ips` identifies an invalid-pointer free inside
  `OpenMM::CpuNeighborList::threadComputeNeighborList` on an OpenMM CPU worker
  thread. The main thread was executing `LangevinMiddleIntegrator::step`.
- This is a native multithreaded OpenMM CPU neighbor-list failure, not a NaN,
  Python exception, invalid checkpoint, or output-format error. The repeated
  five-thread aborts make another unchanged five-thread retry unsafe without a
  resource-setting decision.
- At 2026-11-11 12:35 America/Los_Angeles, seed 2103 remained active at step
  277,500. Seed 1103 was stopped with an intact checkpoint at step 312,500.

### Approved two-thread seed-1103 recovery (2026-11-11)

- The user authorized reducing seed 1103 from five to two CPU threads to avoid
  the recurrent OpenMM CPU neighbor-list memory abort. This changes only host
  computational parallelism; the checkpoint, force field, integrator,
  timestep, temperature, pressure, seed, inputs, and output cadence are
  unchanged.
- The step-312,500 checkpoint was preserved as
  `checkpoint.step312500.pre_threads2_resume.chk` and loaded successfully in a
  two-thread OpenMM CPU context. Before restart, the canonical output contained
  125 unique state rows and 115 DCD frames, exactly matching the checkpoint.
- Seed 1103 now runs as non-respawning LaunchAgent
  `org.openai.codex.dnmt3l.seed1103.threads2.oneshot.20261111` (PID 95271 at
  launch). Seed 2103 remains in its existing five-thread one-shot job.
- At 2026-11-11 12:56 America/Los_Angeles, seed 1103 reached step 315,000 with
  126 unique state rows and 116 DCD frames. Seed 2103 was active at step
  282,500 with 113 unique state rows and 103 DCD frames.

### Seed-2103 five-thread abort (2026-11-11/12)

- Seed 2103 stopped with `SIGABRT` after its durable step-462,500 checkpoint.
  macOS report `python3.11-2026-11-11-232838.ips` shows the same invalid-pointer
  free in `OpenMM::CpuNeighborList::threadComputeNeighborList` seen for seed
  1103. No NaN or Python exception was reported.
- At 2026-11-12 07:51 America/Los_Angeles, the approved two-thread seed-1103
  recovery remained healthy at step 525,000. Seed 2103 was stopped at step
  462,500 pending authorization for the matching two-thread recovery.

### Approved two-thread seed-2103 recovery (2026-11-12)

- The user authorized reducing seed 2103 from five to two CPU threads after
  the matching native OpenMM CPU neighbor-list abort. This changes only host
  computational parallelism; the checkpoint, force field, integrator,
  timestep, temperature, pressure, seed, inputs, and output cadence are
  unchanged.
- Before restart, seed 2103 contained 185 unique state rows through step
  462,500 and 175 production DCD frames, exactly matching the checkpoint
  boundary. The checkpoint was preserved as
  `checkpoint.step462500.pre_threads2_resume.chk` and loaded successfully in
  an OpenMM CPU context with `Threads=2` and a finite potential energy.
- Seed 2103 now runs as non-respawning LaunchAgent
  `org.openai.codex.dnmt3l.seed2103.threads2.oneshot.20261112` (PID 43718 at
  launch). The inactive five-thread job was unloaded before launch, and the
  process inventory contained exactly one writer for each stage-2 replica and
  no batch controller.
- At 2026-11-12 08:13 America/Los_Angeles, seed 2103 reached step 465,000 with
  186 unique state rows and 176 DCD frames. Seed 1103 was active at step
  527,500 with 211 unique state rows and 201 DCD frames. Both frame counts
  exactly matched their report boundaries; aggregate progress was 6,092,500
  of 7,650,000 steps (79.641%).

### Continued two-thread monitoring (2026-11-12)

- At 13:01 America/Los_Angeles, the process inventory still contained exactly
  the two intended stage-2 workers (seed 1103 PID 95271 and seed 2103 PID
  43718), both with `--threads 2`; no batch controller or duplicate writer was
  active. The external volume had 347 GiB free.
- Seed 1103 had an unusually long report gap after step 550,000 while still
  consuming two CPU cores. A read-only five-second process sample showed both
  OpenMM worker threads actively performing Langevin integration, PME/direct
  nonbonded calculations, and neighbor-list work rather than a deadlock. The
  sample is preserved as
  `production_hmr/diagnostics/dnmt3l_seed1103_sample_20261112_1301.txt`
  (SHA-256
  `72b50b98a6565f9c0867560bdb8a93ef088ba3477c3a19ec4ff840a2914610cf`).
- No process intervention was made. At 13:02, seed 1103 advanced normally to
  step 552,500. Its new record remained physically plausible (302.264 K,
  1,419.276 nm^3, 1.01917 g/mL), with 221 unique state rows and 211 DCD frames,
  exactly matching the report boundary.
- At the same validation point, seed 2103 remained healthy at step 487,500
  with 195 unique state rows and 185 DCD frames. Aggregate progress was
  6,140,000 of 7,650,000 steps (80.261%). Both non-respawning LaunchAgents
  remained active on their original PIDs.

### External-drive disconnect recovery (2026-11-12)

- At approximately 16:40 America/Los_Angeles, `${EXTERNAL_VOLUME}` disappeared
  from the mount table while both OpenMM workers were still consuming CPU.
  Because their open output handles were no longer trustworthy, both one-shot
  LaunchAgents were unloaded before any checkpoint-based relaunch. Process and
  launchd checks then confirmed that no simulation worker or batch controller
  remained, eliminating duplicate-writer risk.
- After the drive remounted, the authoritative filesystem outputs were intact
  through step 585,000 for seed 1103 and step 520,000 for seed 2103. Both state
  tables ended with complete newlines and contained unique, monotonic rows;
  their 224 and 198 production DCD frames exactly matched those boundaries.
  Temperatures, volumes, and densities were physically plausible.
- The recovered checkpoints were preserved as
  `checkpoint.step585000.pre_disconnect_resume.chk` and
  `checkpoint.step520000.pre_disconnect_resume.chk`. Each loaded at its exact
  saved step with a finite energy in OpenMM 8.6.1 using the CPU platform and
  `Threads=2`.
- Re-registering the old LaunchAgents produced `EX_CONFIG` before Python
  started because their stale post-remount configuration included an
  unnecessary external-volume `WorkingDirectory`. No output was changed by
  that failed registration. Fresh non-respawning jobs omit only that launchd
  setting while retaining the identical executable, inputs, seeds, output
  directories, integrator, timestep, reporting intervals, and two-thread CPU
  resource setting:
  `org.openai.codex.dnmt3l.seed1103.threads2.disconnectresume.20261112` and
  `org.openai.codex.dnmt3l.seed2103.threads2.disconnectresume.20261112`.
- At 18:24 America/Los_Angeles, seed 1103 (PID 2835 at launch) had advanced to
  step 587,500 and seed 2103 (PID 2836 at launch) to step 522,500. Their 235
  and 209 state rows were unique and their 225 and 199 DCD frames exactly
  matched the new boundaries. Aggregate progress was 6,210,000 of 7,650,000
  steps (81.176%).

### Recurrent external-drive I/O recovery (2026-11-14)

- Both non-respawning workers exited with code 1 after `DCDReporter` received
  `OSError: [Errno 5] Input/output error` while seeking in the trajectory file.
  This was an external-storage interruption, not an OpenMM integration error,
  NaN, or checkpoint incompatibility.
- The authoritative outputs remained clean through step 1,085,000 for seed
  1103 and step 1,020,000 for seed 2103. State rows were complete, unique, and
  monotonic; DCD frame counts exactly matched both boundaries. The checkpoints
  were preserved as `checkpoint.step1085000.pre_io_resume.chk` and
  `checkpoint.step1020000.pre_io_resume.chk`.
- Both checkpoints loaded at their exact steps with finite energies in the
  unchanged OpenMM 8.6.1 CPU environment using `Threads=2`. With no simulation
  or controller process present, the existing non-respawning LaunchAgents were
  restarted as run 2 (PIDs 33007 and 33008 at launch).
- At 19:00 America/Los_Angeles, seed 1103 had advanced to step 1,087,500 and
  seed 2103 to step 1,022,500. Their 435 and 409 state rows remained unique,
  and their 425 and 399 DCD frames exactly matched the new boundaries.
  Aggregate progress was 7,210,000 of 7,650,000 steps (94.248%).

### Seed-1103 completion (2026-11-15)

- Seed 1103 completed normally at step 1,275,000 with launchd exit code 0,
  `COMPLETE.json`, `final_state.xml`, 510 unique state rows, and exactly 500
  production DCD frames. Its completion record confirms the unchanged seed,
  4-fs timestep, 5.0-ns production length, CPU platform, and two-thread
  resource setting.
- At 16:59 America/Los_Angeles, seed 2103 remained the sole active simulation
  worker (PID 33008) at step 1,207,500. It had 483 unique state rows and 473
  production DCD frames, exactly matching the checkpoint boundary, with
  67,500 steps remaining. Aggregate progress was 7,582,500 of 7,650,000 steps
  (99.118%).

### Phase-1 completion and analysis (2026-11-15)

- Seed 2103 completed normally at step 1,275,000. The final six-replica audit
  confirmed 510 unique monotonic state rows, 500 production DCD frames,
  `final_state.xml`, `checkpoint.chk`, and `COMPLETE.json` for every canonical
  replica. Aggregate production time is 30 ns.
- `scripts/analyze_replica.py` read all six canonical trajectories and produced
  500-row per-frame metric tables plus active-loop RMSF profiles. The preserved
  `amniote_YEE_seed1103_duplicate_corrupted_20261111` directory was excluded.
- `scripts/summarize_phase1_md.py` validated replica dimensions and wrote the
  canonical aggregate under `analysis_summary/`. The inference unit is the
  replica summary (three per state), not an individual trajectory frame.
- The DNMT3L-261/DNMT3A-F728 contact was present in every sampled frame in both
  states. Replica-mean distance, protein-backbone RMSD, and active-loop RMSF
  overlapped between states. The DNMT3L-293/DNMT3A1-G726 distance was shorter
  in the F/T/D state in all three matched replica indices, but the mean shift
  was small, did not cross the 0.5-nm contact threshold, and remains exploratory.
- No simulation or trajectory-analysis process was left running after the
  summary completed.

## Deprecated preparation

- `prepared/amniote_YEE` and `prepared/therian_FTD` passed mutation and energy
  validation but used independent 1.0-nm padding boxes (11.3968 versus 11.3419
  nm; 142,570 versus 140,297 atoms). They are retained for provenance and must
  not be used for the matched comparison.
- Primary preparation uses a fixed cubic 11.4-nm box for both variants.

## Completion gates

- Environment imports and platform inventory pass.
- Both prepared structures contain the requested residues at 261/291/293.
- Both minimized systems have finite energy and identical atom counts except
  for the expected mutation-dependent atom difference.
- Both smoke runs terminate normally with parseable trajectory, state log, and
  checkpoint.
- Production trajectories meet the declared step count and replicate count.
- Analysis reads every expected trajectory and reports replicate-level values.

## Validated smoke gates

- The original 2-fs matched smoke runs completed for both endpoints.
- The 3-amu protein-H HMR derivation preserved potential energies to numerical
  precision and did not repartition water hydrogen masses.
- Concurrent 4-fs smoke runs completed normally for both endpoints. Final
  temperatures were 299.715 and 299.006 K; final densities were 1.0027 and
  0.9997 g/mL for amniote and therian states, respectively.
- At five CPU threads per trajectory, paired 4-fs runs achieved approximately
  2.3 ns/day each (about 4.6 ns/day aggregate).

## Rejected branch

- GBn2/NoCutoff implicit solvent is rejected for this comparison. The therian
  endpoint produced a NaN coordinate during minimization. The amniote job was
  terminated after the matched endpoint failed; no implicit trajectory was
  launched and no implicit result may be interpreted.

