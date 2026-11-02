function report = capsule_equalization_demo(config, result)
    % CAPSULE_EQUALIZATION_DEMO Equalize one full-load 16-QAM frame per channel estimate.
    %   The same data bits and noise are used for every method and perfect CSI.
    %   Uses the five-tap TDL-C demo configuration and the finite-frame equalizer.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    obs = result.observation;
    old = rng;
    cleanup = onCleanup(@() rng(old)); %#ok<NASGU>
    rng(obs.seed + 500, 'twister');
    M = config.waveform.M;
    N = config.waveform.N;
    bits = randi([0 1], 4 * M * N, 1);
    symbols = qammod(bits, 16, 'InputType', 'bit', 'UnitAveragePower', true);
    tx = helperOTFSmod(reshape(symbols, M, N), 0, 'NONE');
    truth = obs.truth;
    clean = applyMimoDelayDopplerAngleDelayseq_v1(tx, truth.delays_sec, ...
                                                  truth.dopplers_Hz, truth.gains, truth.angles_rad, config.waveform.fs_Hz, config.array.Nr);
    dataSNR_dB = 10;
    variance = mean(abs(clean(:)).^2) / 10^(dataSNR_dB / 10);
    noisy = clean + sqrt(variance / 2) * (randn(size(clean)) + 1j * randn(size(clean)));
    names = ["PerfectCSI", config.demo.methodNames];
    ber = nan(numel(names), 1);
    flags = nan(numel(names), 1);
    messages = strings(numel(names), 1);
    for k = 1:numel(names)
        if k == 1
            tau = truth.delays_sec;
            nu = truth.dopplers_Hz;
            theta = truth.angles_rad;
            gain = truth.gains;
        else
            estimate = result.estimates{k - 1};
            if isempty(estimate)
                messages(k) = "Pilot estimator failed; equalization skipped.";
                continue
            end
            tau = estimate.tau_sec;
            nu = estimate.nu_Hz;
            theta = estimate.theta_rad;
            gain = estimate.gain;
        end
        try
            operator = makeMimoDelayseqOperator_v1(M * N, config.waveform.fs_Hz, ...
                                                   tau, nu, gain, theta, config.array.Nr, 'ValidationCount', 2);
            [estimated, info] = mimoOtfsLMMSEEqualizerDelayseq_v1(noisy, M, N, ...
                                                                  operator, variance, 1, 'Tol', 1e-7, 'MaxIter', 250);
            flags(k) = info.flag;
            if info.flag ~= 0
                messages(k) = "Iterative solver did not converge; BER withheld.";
                continue
            end
            detected = qamdemod(estimated(:), 16, 'OutputType', 'bit', 'UnitAveragePower', true);
            ber(k) = mean(bits ~= detected);
        catch exception
            messages(k) = string(exception.identifier) + ": " + string(exception.message);
        end
    end
    report = table(names', ber, flags, messages, repmat(dataSNR_dB, numel(names), 1), ...
                   'VariableNames', {'method', 'BER', 'solverFlag', 'message', 'dataSNR_dB'});
end
