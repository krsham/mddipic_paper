function Y = applyMimoDelayDopplerAngleDelayseq_v1(x, tau, nu, g, theta, fs, Na)
    % APPLYMIMODELAYDOPPLERANGLEDELAYSEQ_V1 Exact project physical atom using delayseq.
    %
    % Y = applyMimoDelayDopplerAngleDelayseq_v1(x,tau,nu,g,theta,fs,Na)
    %
    % Implements, path by path,
    %   x(t) -> Doppler -> delayseq delay -> complex gain -> ULA steering
    %
    % using the same physical ordering as delayDopplerGainAngle.
    %
    % Inputs
    %   x      [L x 1] complex TD waveform
    %   tau    [P x 1] seconds
    %   nu     [P x 1] Hz
    %   g      [P x 1] complex path gains
    %   theta  [P x 1] radians
    %   fs     scalar sampling rate [Hz]
    %   Na     scalar receive antennas
    %
    % Output
    %   Y      [L x Na] receive-array TD samples

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    x = x(:);
    tau = tau(:);
    nu = nu(:);
    g = g(:);
    theta = theta(:);
    P = numel(tau);
    assert(numel(nu) == P && numel(g) == P && numel(theta) == P, ...
           'Path parameter vectors must have equal length.');
    L = numel(x);
    t = (0:L - 1).' / fs;
    Y = complex(zeros(L, Na));

    for p = 1:P
        xdop = x .* exp(1j * 2 * pi * nu(p) * t);
        xdel = delayseq(xdop, tau(p), fs);
        a = steeringVector(theta(p), Na);
        Y = Y + (g(p) * xdel) * a.';
    end
end
