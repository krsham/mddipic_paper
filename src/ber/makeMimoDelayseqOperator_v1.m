function op = makeMimoDelayseqOperator_v1(L, fs, tau, nu, g, theta, Na, varargin)
    % MAKEMIMODELAYSEQOPERATOR_V1 SIMO delayseq operator with a basis-probed adjoint.
    %
    % Path ordering:
    %   Doppler -> exact MATLAB delayseq finite operator -> gain -> steering.
    %
    % Each delay matrix is built by calling the installed delayseq on basis
    % vectors. Storage is O(P*L^2); application is O(P*L^2 + P*L*Na).
    % The reference backend supports the capsule demonstration, not BER timing.
    %
    % Name-value:
    %   'DelayTolerance'  default 2e-12
    %   'ValidationCount' default 4
    %   'Seed'            default 20260917

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    p = inputParser;
    p.addParameter('DelayTolerance', 2e-12, @(x)isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('ValidationCount', 4, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('Seed', 20260917, @(x)isnumeric(x) && isscalar(x));
    p.parse(varargin{:});
    o = p.Results;

    tau = tau(:);
    nu = nu(:);
    g = g(:);
    theta = theta(:);
    P = numel(tau);
    assert(numel(nu) == P && numel(g) == P && numel(theta) == P, ...
           'Path parameter vectors must have equal length.');
    assert(L >= 1 && Na >= 1 && fs > 0);

    t = (0:L - 1).' / fs;

    DK = cell(P, 1);
    delayForwardErr = zeros(P, 1);
    delayAdjErr = zeros(P, 1);
    phaseDop = cell(P, 1);
    avec = cell(P, 1);

    for ip = 1:P
        DK{ip} = makeDelayseqReferenceKernel_v1( ...
                                             L, tau(ip), fs, ...
                                             'Tolerance', o.DelayTolerance, ...
                                             'ValidationCount', o.ValidationCount, ...
                                             'Seed', o.Seed + 97 * ip);

        delayForwardErr(ip) = DK{ip}.info.maxForwardRelativeError;
        delayAdjErr(ip) = DK{ip}.info.adjointRelativeError;
        phaseDop{ip} = exp(1j * 2 * pi * nu(ip) * t);
        avec{ip} = steeringVector(theta(ip), Na);
    end

    % Whole-operator adjoint audit.
    old = rng;
    cleanup = onCleanup(@()rng(old)); %#ok<NASGU>
    rng(o.Seed + 7919, 'twister');

    xt = (randn(L, 1) + 1j * randn(L, 1)) / sqrt(2);
    Z = (randn(L, Na) + 1j * randn(L, Na)) / sqrt(2);

    FY = mimo_forward_exact(xt, L, Na, P, g, phaseDop, avec, DK);
    AZ = mimo_adjoint_exact(Z, L, Na, P, g, phaseDop, avec, DK);

    lhs = sum(conj(FY) .* Z, 'all');
    rhs = xt' * AZ;
    adjErr = abs(lhs - rhs) / max([abs(lhs), abs(rhs), realmin]);

    if ~(isfinite(adjErr) && adjErr <= 2e-10)
        error('MimoDelayseqOperator:AdjointFailure', ...
              'Whole SIMO operator adjoint identity failed: %.3e.', adjErr);
    end

    op = struct();
    op.forward = @(x)mimo_forward_exact(x, L, Na, P, g, phaseDop, avec, DK);
    op.adjoint = @(Y)mimo_adjoint_exact(Y, L, Na, P, g, phaseDop, avec, DK);

    op.L = L;
    op.Na = Na;
    op.P = P;
    op.fs = fs;
    op.tau = tau;
    op.nu = nu;
    op.g = g;
    op.theta = theta;

    op.info = struct( ...
                   'schemaVersion', "MIMO_DELAYSEQ_REFERENCE_OPERATOR_V1", ...
                   'delayEngine', "Installed delayseq basis responses", ...
                   'ordering', "Doppler -> delayseq -> gain -> steering", ...
                   'adjointMode', "basis-matrix-conjugate-transpose", ...
                   'maxDelayForwardRelativeError', max(delayForwardErr), ...
                   'maxDelayAdjointRelativeError', max(delayAdjErr), ...
                   'adjointRelativeError', adjErr, ...
                   'pathDelayForwardErrors', delayForwardErr, ...
                   'pathDelayAdjointErrors', delayAdjErr);
end

function Y = mimo_forward_exact(x, L, Na, P, g, phaseDop, avec, DK)
    x = x(:);
    assert(numel(x) == L, 'Forward input must contain L samples.');

    Y = complex(zeros(L, Na));
    for jp = 1:P
        q = g(jp) * DK{jp}.forward(phaseDop{jp} .* x);
        Y = Y + q * avec{jp}.';
    end
end

function x = mimo_adjoint_exact(Y, L, Na, P, g, phaseDop, avec, DK)
    assert(isequal(size(Y), [L Na]), 'Adjoint input must be L x Na.');

    x = complex(zeros(L, 1));
    for jp = 1:P
        % q*a.' -> adjoint Y*conj(a)
        z = conj(g(jp)) * (Y * conj(avec{jp}));
        z = DK{jp}.adjoint(z);
        x = x + conj(phaseDop{jp}) .* z;
    end
end
