function output = capsule_ablation_demo(config, result)
    % CAPSULE_ABLATION_DEMO Compare exact scoring/cache options and DRAP cycles.
    %   Score and cache variants must reconstruct the same response. DRAP cycle
    %   counts are separate searches; equality with MDDIPIC is not expected.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    obs = result.observation;
    labels = ["direct", "schur", "schur_cache", "drap_1", "drap_2"];
    responses = cell(numel(labels), 1);
    elapsed = zeros(numel(labels), 1);
    setupTime = zeros(numel(labels), 1);
    for k = 1:numel(labels)
        if k <= 3
            estimator = capsule_make_estimator('MDDIPIC', config, obs.pilot, ...
                                               obs.varianceDD, obs.domain, 'Schur', k > 1, 'Cache', k == 3);
        else
            estimator = capsule_make_estimator('DRAP', config, obs.pilot, ...
                                               obs.varianceDD, obs.domain, 'Cycles', k - 3);
        end
        guard = onCleanup(@() release(estimator));
        timer = tic;
        setup(estimator);
        setupTime(k) = toc(timer);
        timer = tic;
        [tau, nu, theta, gain] = step(estimator, obs.noisyDD(:), obs.truth.P, obs.varianceDD);
        elapsed(k) = toc(timer);
        response = otfsMimoChannel(obs.pilot, config.waveform.M, config.waveform.fs_Hz, ...
                                   tau, nu, gain, theta, config.array.Nr, 'none', 'NONE', 0);
        responses{k} = response(:);
        clear guard;
    end
    nmse = cellfun(@(x) norm(x - obs.cleanDD(:))^2 / norm(obs.cleanDD(:))^2, responses);
    relativeDifference = cellfun(@(x) norm(x - responses{1}) / max(norm(responses{1}), realmin), responses);
    assert(all(relativeDifference(1:3) < 1e-8), ...
           'Direct/Schur/cache response identity exceeded 1e-8.');
    output = table(labels', setupTime, elapsed, nmse, relativeDifference, ...
                   'VariableNames', {'variant', 'setup_s', 'estimate_s', 'responseNMSE', 'differenceFromDirect'});
end
