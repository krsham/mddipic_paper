function score = score_isac_raw_estimate_v2( ...
                                            truth, est, cfg, domain, Yclean, x_dd, sys, varargin)
    % SCORE_ISAC_RAW_ESTIMATE_V2 Metrics with optional full-response override.
    %
    % ResponseEstimate permits structured estimators to be scored from their
    % complete retained response rather than an approximate path-proxy collapse.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    p = inputParser;
    p.addParameter('ResponseEstimate', [], @(x)isnumeric(x) || isempty(x));
    p.addParameter('ParameterMetricsValid', true, @(x)islogical(x) && isscalar(x));
    p.addParameter('ParameterMetricNote', '', @(x)ischar(x) || isstring(x));
    p.parse(varargin{:});
    opt = p.Results;

    score = empty_score();
    score.schemaVersion = 'ISAC_SCORE_V2';
    score.parameterMetricsValid = opt.ParameterMetricsValid;
    score.parameterMetricNote = char(opt.ParameterMetricNote);

    P = truth.P;
    if ~est.success || est.P_est ~= P
        score.failureReason = 'estimator failure or path-count mismatch';
        return;
    end

    if opt.ParameterMetricsValid
        assoc = associate_isac_estimates_v2(truth, est.tau_sec, est.nu_Hz, ...
                                          est.theta_rad, cfg, domain);
        score.association = assoc;
        if ~assoc.valid
            score.failureReason = 'association invalid';
            return;
        end

        q = assoc.truthToEstimate;
        score.truthToEstimate = q(:);
        tau = est.tau_sec(q);
        nu = est.nu_Hz(q);
        th = est.theta_rad(q);
        g = est.gain(q);
        dt = tau(:) - truth.delays_sec(:);
        dn = nu(:) - truth.dopplers_Hz(:);
        da = wrap_pi(th(:) - truth.angles_rad(:));
        dg = g(:) - truth.gains(:);

        score.delayError_s = dt;
        score.delayError_bins = dt / cfg.waveform.delayBin_s;
        score.dopplerError_Hz = dn;
        score.dopplerError_bins = dn / cfg.waveform.dopplerBin_Hz;
        score.angleError_rad = da;
        score.angleError_deg = rad2deg(da);
        score.gainError = dg;
        score.gainRelativeError = abs(dg) ./ max(abs(truth.gains(:)), realmin);

        if strcmpi(domain, 'communication')
            % LoS-specific metrics are only meaningful when the truth model
            % explicitly contains a designated LoS path. TDL-B is NLOS.
            if isfield(truth, 'hasLoS') && truth.hasLoS
                score.losDelayError_bins = score.delayError_bins(1);
                score.losDopplerError_bins = score.dopplerError_bins(1);
                score.losAngleError_deg = score.angleError_deg(1);
            end
            if isfield(truth, 'pathPowers') && numel(truth.pathPowers) == P
                w = truth.pathPowers(:);
            else
                w = abs(truth.gains(:)).^2;
            end
            w = w / max(sum(w), realmin);
            score.powerWeightedDelaySquaredError_bins2 = ...
                sum(w .* abs(score.delayError_bins).^2);
            score.powerWeightedDopplerSquaredError_bins2 = ...
                sum(w .* abs(score.dopplerError_bins).^2);
            score.powerWeightedAngleSquaredError_deg2 = ...
                sum(w .* abs(score.angleError_deg).^2);
        end

        if strcmpi(domain, 'sensing')
            c = cfg.constants.c_mps;
            fc = cfg.waveform.fc_Hz;
            score.rangeError_m = c * dt / 2;
            score.radialVelocityError_mps = c * dn / (2 * fc);
            score.radialVelocityError_kmh = 3.6 * score.radialVelocityError_mps;
        end
    end

    if isempty(opt.ResponseEstimate)
        [Yhat, ~, ~, ~] = otfsMimoChannel(x_dd, sys.M, sys.fs, ...
                                     est.tau_sec(:), est.nu_Hz(:), est.gain(:), est.theta_rad(:), ...
                                     sys.Na, 'none', sys.padtype, sys.padlen);
        score.responseScoringRule = 'reconstructed_physical_path_response';
    else
        Yhat = opt.ResponseEstimate(:);
        score.responseScoringRule = 'complete_retained_estimator_response';
    end

    if numel(Yhat) ~= numel(Yclean) || any(~isfinite(Yhat(:)))
        score.failureReason = 'invalid response estimate';
        return;
    end
    nmse = norm(Yclean(:) - Yhat(:))^2 / max(norm(Yclean(:))^2, realmin);
    score.responseNMSE = nmse;
    score.responseNMSE_dB = 10 * log10(max(nmse, realmin));
    score.valid = isfinite(nmse);
    score.failureReason = '';
end

function score = empty_score()
    score = struct('schemaVersion', '', 'valid', false, ...
                 'parameterMetricsValid', false, 'parameterMetricNote', '', ...
                 'association', struct(), 'truthToEstimate', [], ...
                 'responseNMSE', Inf, 'responseNMSE_dB', Inf, ...
                 'responseScoringRule', '', ...
                 'delayError_s', [], 'delayError_bins', [], ...
                 'dopplerError_Hz', [], 'dopplerError_bins', [], ...
                 'angleError_rad', [], 'angleError_deg', [], ...
                 'gainError', [], 'gainRelativeError', [], ...
                 'losDelayError_bins', [], 'losDopplerError_bins', [], ...
                 'losAngleError_deg', [], ...
                 'powerWeightedDelaySquaredError_bins2', [], ...
                 'powerWeightedDopplerSquaredError_bins2', [], ...
                 'powerWeightedAngleSquaredError_deg2', [], ...
                 'rangeError_m', [], 'radialVelocityError_mps', [], ...
                 'radialVelocityError_kmh', [], 'failureReason', '');
end

function x = wrap_pi(x)
    x = mod(x + pi, 2 * pi) - pi;
end
