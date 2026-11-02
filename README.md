# MDDIPIC MATLAB paper code

MATLAB implementations of MDDIPIC and DRAP for multi-antenna OTFS sensing and channel estimation, with DDAFR, MDDAFR, OMP and NOMP baselines. A small reproducible demonstration covers sensing, communication, model-order stopping, computational options and LMMSE equalization.

## Requirements

MATLAB R2025b is the reference release. Required toolboxes are Communications Toolbox (`qammod`, `qamdemod`, optional `awgn`) and Phased Array System Toolbox (`delayseq`). Parallel Computing Toolbox is optional. MATLAB Java support is used for source/checkpoint signatures.

## Quick start

Clone or download this repository and make its root the current MATLAB folder. In a fresh MATLAB session, run:

```matlab
setup;
report = run_capsule_demo;
```

The demo uses a 32×16 full-DD QPSK pilot, eight receive antennas, three sensing targets, five TDL-C communication taps, and pilot SNRs of 0 and 20 dB. All six estimators share the same observations. The run also compares stopping criteria, direct/Schur/cache implementations and DRAP refinement cycles, then equalizes one 16-QAM frame.

This preset uses reduced dimensions and trial counts. It demonstrates the algorithms rather than regenerating the manuscript's Monte Carlo figures. The equalization example uses TDL-C; the manuscript BER study uses its separately specified TDL-B configuration.

## Options

```matlab
report = run_capsule_demo('Trials', 2, 'SNRdB', [0 10 20]);
report = run_capsule_demo('UseParallel', true);
report = run_capsule_demo('ShowFigures', false, 'OutputDir', 'results/batch_check');
checks = run_release_checks;
```

Serial execution is the default. Resume reuses checkpoints only when the source, configuration and MATLAB environment match. After editing source or configuration, use a new output directory. `run_release_checks` exercises numerical contracts and aggregation without fitting estimators.

## Results

Generated outputs go to `results/capsule_demo/`:

| Files | Contents |
|---|---|
| `run_config.json`, `environment.json`, `console.txt` | Configuration, environment and log |
| `case_*.mat` | Per-observation checkpoints |
| `estimator_metrics.csv`, `summary.csv` | Per-fit errors/status, linear mean NMSE and median time |
| `stopping_histories.csv`, `stopping_decisions.csv` | Acquisition statistics and model-order decisions |
| `computational_ablation.csv` | Score/cache and refinement-cycle comparisons |
| `equalization.csv` | BER and solver convergence |
| `replacement_checks.csv`, `numerical_checks.csv`, `summary_checks.csv` | Numerical and reporting checks |
| `demo_summary.pdf`, `demo_summary.png` | Demo plots |
| `review_report.mat`, `completion.json` | Combined report and completion status |

The demo reports failures explicitly. Parallel timings and single-trial results should be interpreted within this small workload. See [LIBRARY_STRUCTURE.md](LIBRARY_STRUCTURE.md) for interfaces, equations, dimensions and conventions. Numerical contract checks passed on MATLAB R2025b/Windows; use `run_release_checks` to repeat them.

The delay reference backend builds a dense matrix from installed `delayseq` responses. It uses quadratic storage per path and is intended for small examples; its timing does not characterize an FFT implementation.

## License and citation

Copyright © 2026 Koorosh Ahmadi. Original software is available under GPL-3.0-only with the additional MATLAB permission in [LICENSES/MATLAB_EXCEPTION.txt](LICENSES/MATLAB_EXCEPTION.txt). MATLAB and toolbox licenses are obtained separately.

Please cite the associated manuscript and the software version used: [CITATION.md](CITATION.md), [CITATION.cff](CITATION.cff). Algorithm and dependency attribution appears in [REFERENCES.bib](REFERENCES.bib) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
