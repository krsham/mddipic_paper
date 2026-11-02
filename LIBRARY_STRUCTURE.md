# Library structure and numerical conventions

The library has one physical path model and several estimators. `run_capsule_demo.m` orchestrates a small experiment; it does not contain an alternative estimator. The public classes retain the estimator interfaces used by the manuscript.

## Folder map

| Location | Responsibility |
|---|---|
| `run_release_checks.m` | Focused replacement and algebraic checks, no estimator experiment |
| `setup.m` | Add explicit package folders; check dependencies and path resolution |
| `run_capsule_demo.m` | Serial/parallel job plan, resumable checkpoints, output tables and completion status |
| `demo/` | Small configuration, observations, estimator factory, stopping/ablation/equalization exercises and plots |
| `src/common/` | Physical Doppler/delay/gain/array model and end-to-end frame helpers |
| `src/proposed/` | MDDIPIC and DRAP numerical implementations |
| `src/baselines/` | DDAFR, MDDAFR, fixed-grid OMP and 3-D NOMP |
| `src/ber/` | Finite-delay forward/adjoint operators and iterative SIMO LMMSE equalization |
| `src/scenario/` | Sensing geometry, TDL profile tables and communication truth generation |
| `src/evaluation/` | One-to-one association and reconstruction/coordinate scoring |
| `tests/` | Small algebraic contracts invoked by the main runner |
| `results/` | Generated local outputs; empty in the supplied ZIP |

No other checkout, base-workspace variable, personal directory or exploratory branch is needed.

## Array layout and units

| Quantity | Shape | Meaning |
|---|---|---|
| `pilot` / `XDD` | M × N | Complex DD symbols; delay index varies down rows |
| `tx` | L × 1 | Serialized TD frame; L=MN in the no-padding configuration |
| `rxTD` | L × Na | One receive antenna per column |
| `YDD` | M × N × Na | Demodulated received grids |
| `yDD` | MN·Na × 1 | `YDD(:)`, contiguous MN-sample antenna blocks |
| `tau` | P × 1 | Seconds |
| `nu` | P × 1 | Hz |
| `theta` | P × 1 | Radians, broadside-referenced receive angle |
| `gain` | P × 1 | Complex effective amplitude, including the single-transmit beam response |
| `varianceTD` | scalar | E[abs(w)^2] per complex TD sample and receive antenna |
| `varianceDD` | scalar | M²N · varianceTD for the included no-padding transform |

Noise standard deviation per real/imaginary component is sqrt(varianceTD/2). The demo draws one unit-variance complex noise array for each truth scene and scales it with SNR. It does not normalize each noise draw to a fixed realized norm.

## Physical model

`delayDopplerGainAngle` forms one path as

```
t = (0:L-1)' / fs;
Yp = gain(p) * delayseq(tx .* exp(1j*2*pi*nu(p)*t), tau(p), fs) ...
     * steeringVector(theta(p),Na).';
```

The order is Doppler, finite delay, gain, receive steering. The final transpose is `.'`, not `'`: the array response is not conjugated when generating the receive waveform. The steering convention is

```
a(theta) = exp(1j*pi*sin(theta) * ((0:Na-1)' - Na/2));
```

Finite-frame delay/truncation is used; this is not a circular full-MIMO channel. `otfsMimoChannel` modulates the pilot, applies the path model, adds explicitly supplied noise, and demodulates each antenna. `otfsChannel` is its single-channel counterpart used while preparing the sensing OMP dictionary. `applyMimoDelayDopplerAngleDelayseq_v1` supplies the corresponding TD forward operation for data equalization.

`helperOTFSmod` uses `vec(ifft(XDD.').'/M)` without padding. `helperOTFSdemod` performs the inverse map. Therefore TD covariance for independent unit-power DD symbols is I/(M²N). The scale enters both estimator noise units and equalizer regularization. The compatibility names delegate to the new project transforms `otfs_modulate` and `otfs_demodulate`; see THIRD_PARTY_NOTICES.md.

## Estimator interfaces

Use `capsule_make_estimator` as a complete construction example. Construct → set non-tunable properties → `setup` → `reset` before a fresh observation → `step` → `release` is the lifecycle. A changed pilot, grid or rho generally requires release/reconstruction and setup. The demo creates an estimator per observation so its noise-variance setting is explicit.

| Class | Input to `step` | Noise argument | Role |
|---|---|---|---|
| `MDDIPICEstimator` | Stacked DD vector | DD variance | Public joint 3-D search, concentrated active-support score, joint gains and path revisits |
| `AggressiveMDDIPICEstimator` | Stacked DD vector | DD variance | Implementation superclass; `MDDIPICEstimator` fixes joint refinement |
| `DRAPMDDIPICEstimator` | Stacked DD vector | DD variance | Residual DD/spatial projections, alternating fine refinement and joint gains |
| `DDAFRBenchmarkEstimator` | Stacked DD vector | DD variance | Coarse/fine joint dictionary search, scalar LS gain and sequential cancellation |
| `MDDAFRBenchmarkEstimator` | Stacked DD vector | DD variance | MUSIC angular candidates followed by local DD search and cancellation |
| `Classic3DOMPSensingEstimator` | Stacked DD vector | DD variance | Fixed-grid sensing OMP, joint LS gain refits |
| `SrivastavaOMPBenchmarkEstimator` | TD array/vector | TD variance | Fixed-grid communication OMP, joint LS refits |
| `Li3DNOMPCommonHostEstimator` | TD array/vector | TD variance | Continuous angle/Doppler Newton updates, fractional-delay refinement and cyclic revisits |

All expose `[tau,nu,theta,gain,count] = step(estimator, samples, knownOrder, variance)`. Some classes additionally return an information struct. `knownOrder=0` has estimator-specific semantics; the demonstration uses supplied order for six-method accuracy comparisons and evaluates MDDIPIC stopping separately. MDDAFR known-order MUSIC needs P<Na.

The proposed classes do **not** infer their `rho` property from the third `step` argument. The factory explicitly sets `rho=varianceDD` for the demo. Native baselines retain their LS/scalar-gain rules. Physical inputs are shared, but refinement budgets and objective functions differ.

### MDDIPIC score and computational switches

For the active dictionary R, the concentrated score is

```
y' * R * ((R'*R + rho*I) \ (R'*y)).
```

Each atom factors as `kron(a(theta),h(tau,nu))`. The Schur realization updates the concentrated score without a fresh full support solve for every candidate. `UseSchurScoreUpdate=false` selects the direct augmented solve; `true` selects Schur. `EnableDDResponseCache=true` stores temporal responses only for the current observation. `capsule_ablation_demo` compares reconstructed responses for direct, Schur and cached Schur at the same grids and rho, and requires relative differences below 1e-8.

DRAP replaces the full support-conditioned candidate search with normalized residual projections. `NumAlternatingRefinements` is the number of complete DD/angle cycles. Its estimates are not expected to equal MDDIPIC, or to remain identical when the cycle count changes.

### Stopping

The diagnostic histories are `NewPathCoarsePeakScoreHistory` and `ResidualChangeHistory`. Only new-path acquisition peaks enter the drop statistic; IPI peaks are excluded. The rule is

```
drop_dB(p) = 10*log10(q(p-1)/q(p));   p >= 2.
```

A crossing rejects candidate p and retains p-1 paths. The native peak-drop branch tests after local refinement and before appending the candidate. The demonstration checks it against a decision reconstructed from the full-order acquisition history. It also reconstructs the paper's candidate-reject residual reference. The native residual branch accepts a provisional candidate before checking its residual-change statistic, so its returned count is not labeled equivalent to that reconstructed reference.

Thresholds 6 dB and residual scale 0.1 are fixed demonstration settings from the paper setup, not newly calibrated small-grid thresholds. The demo covers nonempty scenes and is not a zero-target detection test.

## Channel truth

`generate_regime_traffic_v2` places point targets at a common longitudinal coordinate near 90 m, with the requested lateral spacing, a 5 m receiver height and 0.8 m target height. Round-trip delay/Doppler follow geometry; effective gains have independent random phases and unit magnitude.

`tdl_abc_profile_v1` contains the 3GPP normalized TDL-A/B/C tables. `generate_comm_channel_v1` supplies the repository's Rayleigh/Jakes/receive-angle construction; `generate_tdl_abc_channel_v1` applies strongest-five profile selection. The demo selects C. The 634 ns scale refers to the full profile before truncation; the retained five-tap RMS spread is different. AoA clustering is a project extension, not an angle model specified by the scalar TDL table.

`make_isac_scenario_config_v1` contains the original scenario defaults. `capsule_demo_config` explicitly overrides waveform dimensions, receive antennas and demo budgets. Inspect the saved `run_config.json` for the effective settings.

## Reconstruction, association and aggregation

`score_isac_raw_estimate_v2` synthesizes a clean response from the estimated paths and computes `norm(Yclean-Yhat)^2/norm(Yclean)^2`. It requires a valid estimate with the supplied path count. `associate_isac_estimates_v2` uses a one-to-one dynamic-programming assignment over normalized delay, Doppler and angle errors. Communication rows are weighted by realized truth tap power; sensing weights are equal. The assignment implementation is intended for P<=20.

The demo averages per-trial NMSE in the linear domain before displaying dB. It retains failed fits in the detailed table and withholds a group aggregate if any requested fit failed. Use larger trial counts for statistical comparisons.

## Equalization and adjoint

`makeDelayseqReferenceKernel_v1` builds D by applying the installed `delayseq` to separate complex unit vectors. Its Hermitian adjoint is D'. Storage is 16*L^2 bytes per path; construction requires L toolbox calls per path and application costs O(L^2). This dense reference implementation is intended for small frames; timings should identify the backend. Applying `delayseq` with negative delay is not generally that adjoint.

`makeMimoDelayseqOperator_v1` combines these path operators with Doppler, gains and steering. It supplies `op.forward(x)` and `op.adjoint(Y)`. Construction checks installed-delayseq forward agreement and the adjoint inner-product identity. `mimoOtfsLMMSEEqualizerDelayseq_v1` solves

```
(A'*A + lambda*I) * tHat = A'*y;
lambda = varianceTD * M^2 * N / symbolVarianceDD;
```

PCG is used first, with an augmented-LSQR fallback when needed. The recovered TD waveform is demodulated to DD symbols. Solver diagnostics are kept; failed solves do not receive a BER value in the demo. `check_capsule_numerics` compares this iterative solution with a small independently basis-probed dense system. Each basis column is propagated separately; delaying `eye(L)(:)` and reshaping is not used.

## File-by-file entry points

| File | Use |
|---|---|
| `capsule_demo_config` | Reduced workload and search budgets |
| `capsule_observation` | Truth, QPSK pilot, paired AWGN and TD/DD variance scaling |
| `capsule_make_estimator` | Public class configuration for all six methods |
| `capsule_run_case` | Known-order fitting, timing, response scoring and failure records |
| `capsule_stopping_demo` | Full-order histories, fixed thresholds and native peak-drop comparison |
| `capsule_ablation_demo` | Exact score/cache checks and DRAP cycle comparison |
| `capsule_equalization_demo` | One paired 16-QAM frame with all estimated channels and perfect CSI |
| `capsule_plot_summary` | Save linear-averaged NMSE and median-time plots |
| `capsule_source_signature` | Prevent checkpoint reuse after source edits |
| `check_capsule_numerics` | Waveform, physical operator, adjoint and dense LMMSE checks |

`apply_farrow_delay_1d` is retained solely for the optional legacy delay-method argument in the common helpers. Its active coefficients are three sinc taps, despite the historical file name. It is not used by the paper-aligned demo's `delayseq` path. 


## Release checks and checkpoint boundary

`run_release_checks` calls the direct-DFT transform/padding tests, installed-delay forward/adjoint tests, small physical-channel/dense-LMMSE checks and synthetic aggregation regression. It writes four compact evidence files without estimator fitting. `run_capsule_demo` invokes the same contract checks before its experiment.

`capsule_summarize_cases` preallocates a typed record schema. `check_capsule_summary` covers linear NMSE averaging, median runtime, partly failed groups and empty groups. The five aggregation contracts passed on MATLAB R2025b.

Checkpoint reuse requires matching source, configuration and environment signatures.

Original software ownership: Koorosh Ahmadi. Manuscript authorship and external method attribution are separate; see COPYRIGHT.txt.
