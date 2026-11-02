function report = capsule_stopping_demo(config, result)
    % CAPSULE_STOPPING_DEMO Inspect full-order histories and candidate rejection.
    %   Both reported decisions are reconstructed from a common acquisition
    %   sequence. This matches the paper figure's offline decision convention.
    %   The small demo does not recalibrate the thresholds or estimate a ROC.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    obs = result.observation;
    estimator = capsule_make_estimator('MDDIPIC', config, obs.pilot, ...
                                       obs.varianceDD, obs.domain);
    guard = onCleanup(@() release(estimator)); %#ok<NASGU>
    setup(estimator);
    step(estimator, obs.noisyDD(:), config.demo.maximumOrder, obs.varianceDD);
    peaks = estimator.NewPathCoarsePeakScoreHistory(:);
    changes = estimator.ResidualChangeHistory(:);
    drops = [NaN; 10 * log10(max(peaks(1:end - 1), realmin) ./ max(peaks(2:end), realmin))];
    residualThreshold = config.demo.residualScale * config.array.Nr * ...
        config.waveform.M * config.waveform.N * obs.varianceDD;
    peakCrossing = find((1:numel(peaks))' > 1 & drops >= config.demo.peakDrop_dB, 1);
    residualCrossing = find((1:numel(changes))' > 1 & changes < residualThreshold, 1);
    peakOrder = config.demo.maximumOrder;
    residualOrder = config.demo.maximumOrder;
    if ~isempty(peakCrossing)
        peakOrder = peakCrossing - 1;
    end
    if ~isempty(residualCrossing)
        residualOrder = residualCrossing - 1;
    end
    histories = table((1:numel(peaks))', peaks, drops, changes, ...
                      'VariableNames', {'candidate', 'acquisitionPeak', 'drop_dB', 'residualChangeEnergy'});
    decisions = table(["peakdrop"; "residual"], [peakOrder; residualOrder], ...
                      repmat(obs.truth.P, 2, 1), ...
                      'VariableNames', {'rule', 'estimatedOrder', 'trueOrder'});

    % The native residual branch accepts its provisional candidate before
    % testing. Show that convention separately rather than equating counts.
    native = capsule_make_estimator('MDDIPIC', config, obs.pilot, obs.varianceDD, obs.domain);
    nativeGuard = onCleanup(@() release(native)); %#ok<NASGU>
    native.AutoStopRule = 'peakdrop';
    setup(native);
    [~, ~, ~, ~, nativePeakOrder] = step(native, obs.noisyDD(:), 0, obs.varianceDD);
    assert(nativePeakOrder == peakOrder, 'Native and reconstructed peak-drop decisions differ.');
    report = struct('histories', histories, 'decisions', decisions, ...
                    'nativePeakOrder', nativePeakOrder, 'residualThreshold', residualThreshold);
end
