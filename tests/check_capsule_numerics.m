function checks = check_capsule_numerics()
    % CHECK_CAPSULE_NUMERICS Verify waveform normalization and finite-channel adjoints.
    %   Tests use small deterministic arrays. They check algebraic contracts,
    %   not estimator rankings or the statistical claims in the manuscript.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    previous = rng;
    cleanup = onCleanup(@() rng(previous)); %#ok<NASGU>
    rng(2601102, 'twister');
    M = 8;
    N = 4;
    fs = 10e6;
    antennas = 4;
    grid = (randn(M, N) + 1j * randn(M, N)) / sqrt(2);
    signal = helperOTFSmod(grid, 0, 'NONE');
    recovered = helperOTFSdemod(signal, M, 0, 0, 'NONE');
    roundTrip = norm(grid(:) - recovered(:)) / norm(grid(:));
    energyScaling = abs(norm(grid(:))^2 / norm(signal)^2 - M^2 * N) / (M^2 * N);
    assert(roundTrip < 1e-12 && energyScaling < 1e-12, ...
           'OTFS round-trip or transform-energy scaling failed.');

    tau = [0.37; 2.18] / fs;
    nu = [-740; 1300];
    gain = [1; 0.4 * exp(0.7j)];
    theta = deg2rad([-15; 20]);
    direct = delayDopplerGainAngle(signal, tau, nu, gain, theta, fs, antennas);
    independent = applyMimoDelayDopplerAngleDelayseq_v1(signal, tau, nu, gain, theta, fs, antennas);
    channelParity = norm(direct(:) - independent(:)) / norm(direct(:));
    assert(channelParity < 1e-12, 'The two physical forward implementations disagree.');

    operator = makeMimoDelayseqOperator_v1(M * N, fs, tau, nu, gain, theta, antennas);
    % Probe each TD basis vector independently. Flattening an identity before
    % delayseq would connect columns at their boundaries and is not equivalent.
    dense = complex(zeros(M * N * antennas, M * N));
    for k = 1:M * N
        basis = zeros(M * N, 1);
        basis(k) = 1;
        response = operator.forward(basis);
        dense(:, k) = response(:);
    end
    forwardError = norm(dense * signal - direct(:)) / norm(direct(:));
    assert(forwardError < 2e-11, 'Dense and reference forward operators disagree.');

    varianceTD = 1e-5;
    [estimated, info] = mimoOtfsLMMSEEqualizerDelayseq_v1(direct, M, N, operator, ...
                                                          varianceTD, 1, 'Tol', 1e-11, 'MaxIter', 300);
    lambda = varianceTD * M^2 * N;
    exactTD = (dense' * dense + lambda * eye(M * N)) \ (dense' * direct(:));
    exactDD = helperOTFSdemod(exactTD, M, 0, 0, 'NONE');
    equalizerError = norm(estimated(:) - exactDD(:)) / norm(exactDD(:));
    assert(info.flag == 0 && equalizerError < 1e-8, ...
           'Iterative LMMSE does not agree with the small dense solution.');
    checks = table(["OTFS round trip"; "OTFS energy scaling"; "physical forward parity"; ...
                    "operator forward parity"; "operator adjoint"; "dense LMMSE parity"], ...
                   [roundTrip; energyScaling; channelParity; forwardError; ...
                    operator.info.adjointRelativeError; equalizerError], ...
                   'VariableNames', {'check', 'relativeError'});
end
