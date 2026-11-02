function observation = capsule_observation(config, domain, trial, snrDB)
    % CAPSULE_OBSERVATION Generate paired pilot observations with explicit noise units.
    %   The same domain/trial has the same truth, pilot and unit-noise samples
    %   at every SNR. Only the noise standard deviation changes with snrDB.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    seed = 2601100 + 10000 * double(string(domain) == "communication") + trial;
    previous = rng;
    cleanup = onCleanup(@() rng(previous)); %#ok<NASGU>
    rng(seed, 'twister');
    M = config.waveform.M;
    N = config.waveform.N;
    pilot = qammod(randi([0 3], M, N), 4, 'UnitAveragePower', true);
    if string(domain) == "sensing"
        truth = generate_regime_traffic_v2(config, config.demo.sensingPaths, ...
                                           config.demo.spacing_m, seed + 100);
    else
        truth = generate_tdl_abc_channel_v1(config, 'C', seed + 100);
    end
    [cleanDD, tx, cleanTD] = otfsMimoChannel(pilot, M, config.waveform.fs_Hz, ...
                                             truth.delays_sec, truth.dopplers_Hz, truth.gains, truth.angles_rad, ...
                                             config.array.Nr, 'none', 'NONE', 0);

    % The transform is scaled, not unitary: sigmaDD^2 = M^2*N*sigmaTD^2.
    varianceTD = mean(abs(cleanTD(:)).^2) / 10^(snrDB / 10);
    unitNoise = (randn(size(cleanTD)) + 1j * randn(size(cleanTD))) / sqrt(2);
    noisyTD = cleanTD + sqrt(varianceTD) * unitNoise;
    noisyDD = complex(zeros(size(cleanDD)));
    for antenna = 1:config.array.Nr
        noisyDD(:, :, antenna) = helperOTFSdemod(noisyTD(:, antenna), M, 0, 0, 'NONE');
    end
    truthParameters = [truth.delays_sec(:), truth.dopplers_Hz(:), truth.angles_rad(:)];
    assert(all(isfinite(truthParameters(:))), 'Nonfinite channel parameters.');
    assert(max(truth.delays_sec) <= max(config.demo.tauGrid), ...
           'Demo coarse delay support is too short for this realization.');
    observation = struct('truth', truth, 'pilot', pilot, 'tx', tx, ...
                         'cleanDD', cleanDD, 'cleanTD', cleanTD, 'noisyDD', noisyDD, 'noisyTD', noisyTD, ...
                         'varianceTD', varianceTD, 'varianceDD', M^2 * N * varianceTD, ...
                         'domain', string(domain), 'trial', trial, 'snrDB', snrDB, 'seed', seed);
end
