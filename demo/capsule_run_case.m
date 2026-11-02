function result = capsule_run_case(config, domain, trial, snrDB)
    % CAPSULE_RUN_CASE Fit all six estimators to one identical noisy observation.
    %   Records setup and estimation time separately, with explicit fit status.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    obs = capsule_observation(config, domain, trial, snrDB);
    names = config.demo.methodNames;
    rows = cell(numel(names), 1);
    estimates = cell(numel(names), 1);
    for k = 1:numel(names)
        row = struct('domain', string(domain), 'trial', trial, 'snr_dB', snrDB, ...
                     'method', names(k), 'trueOrder', obs.truth.P, 'estimatedOrder', NaN, ...
                     'responseNMSE', NaN, 'responseNMSE_dB', NaN, ...
                     'delayRMSE_bins', NaN, 'dopplerRMSE_bins', NaN, 'angleRMSE_deg', NaN, ...
                     'setup_s', NaN, 'estimate_s', NaN, 'success', false, 'message', "");
        estimator = [];
        try
            estimator = capsule_make_estimator(names(k), config, obs.pilot, obs.varianceDD, domain);
            guard = onCleanup(@() release(estimator));
            timer = tic;
            setup(estimator);
            row.setup_s = toc(timer);
            reset(estimator);
            if names(k) == "NOMP" || (names(k) == "OMP" && string(domain) == "communication")
                samples = obs.noisyTD(:);
                variance = obs.varianceTD;
            else
                samples = obs.noisyDD(:);
                variance = obs.varianceDD;
            end
            timer = tic;
            [tau, nu, theta, gain, count] = step(estimator, samples, obs.truth.P, variance);
            row.estimate_s = toc(timer);
            estimate = struct('tau_sec', tau(:), 'nu_Hz', nu(:), ...
                              'theta_rad', theta(:), 'gain', gain(:), ...
                              'P_est', count, 'success', true);
            sys = struct('M', config.waveform.M, 'fs', config.waveform.fs_Hz, ...
                         'Na', config.array.Nr, 'padtype', 'NONE', 'padlen', 0);
            score = score_isac_raw_estimate_v2(obs.truth, estimate, config, domain, ...
                                               obs.cleanDD, obs.pilot, sys);
            row.estimatedOrder = count;
            row.success = score.valid;
            row.message = string(score.failureReason);
            if score.valid
                row.responseNMSE = score.responseNMSE;
                row.responseNMSE_dB = score.responseNMSE_dB;
                % Match the paper: equal target weights for sensing; realized
                % tap-power weights for communication coordinate errors.
                weights = score.association.truthWeights;
                row.delayRMSE_bins = sqrt(sum(weights .* score.delayError_bins.^2));
                row.dopplerRMSE_bins = sqrt(sum(weights .* score.dopplerError_bins.^2));
                row.angleRMSE_deg = sqrt(sum(weights .* score.angleError_deg.^2));
            end
            estimates{k} = estimate;
            clear guard;
        catch exception
            row.message = string(exception.identifier) + ": " + string(exception.message);
            if exist('guard', 'var')
                clear guard;
            elseif ~isempty(estimator)
                release(estimator);
            end
        end
        rows{k} = row;
    end
    result = struct('rows', vertcat(rows{:}), 'estimates', {estimates}, ...
                    'observation', obs, 'schema', 'MDDIPIC_REVIEW_CASE_V1');
end
