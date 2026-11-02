function truth = generate_regime_traffic_v2(cfg, P, spacing, seed)
    % GENERATE_REGIME_TRAFFIC_V2 Highway point-target geometry with variable count/spacing.
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    validateattributes(P, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(spacing, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    previous = rng;
    cleanup = onCleanup(@()rng(previous)); %#ok<NASGU>
    rng(seed, 'twister');
    receiver = [0 0 5];
    jitter = 2 * (rand - 0.5);
    positions = [((0:P - 1)' - (P - 1) / 2) * spacing, repmat(90 + jitter, P, 1), repmat(0.8, P, 1)];
    velocities = repmat([0 60 / 3.6 0], P, 1);
    offset = positions - receiver;
    ranges = vecnorm(offset, 2, 2);
    radial = sum(offset .* velocities, 2) ./ ranges;
    gains = exp(1j * 2 * pi * rand(P, 1));
    truth = struct('P', P, 'delays_sec', 2 * ranges / cfg.constants.c_mps, ...
                   'dopplers_Hz', 2 * cfg.waveform.fc_Hz * radial / cfg.constants.c_mps, ...
                   'angles_rad', asin(offset(:, 1) ./ ranges), 'gains', gains, 'pathPowers', abs(gains).^2, ...
                   'targetMask', true(P, 1), 'objectGroups', (1:P)', 'objectCount', P, ...
                   'positions_m', positions, 'velocities_mps', velocities, 'receiver_m', receiver, ...
                   'scene', "traffic", 'layout', "spacing_" + string(spacing));
end
