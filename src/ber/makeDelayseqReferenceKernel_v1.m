function kernel = makeDelayseqReferenceKernel_v1(L, tau, fs, varargin)
    % MAKEDELAYSEQREFERENCEKERNEL_V1 Finite complex-baseband delay and its adjoint.
    %   The installed DELAYSEQ is evaluated on each unit vector separately.
    %   Its responses form D(:,k); forward and adjoint operations are D*x and
    %   D'*y. Requires Phased Array System Toolbox.
    %
    %   Each probe is purely imaginary so DELAYSEQ takes its complex-input path.
    %   Dividing its response by 1j recovers the complex-linear basis response.
    %   Random complex-vector parity is checked before returning the operator.
    %   A negative delay is not used as an approximation to the adjoint.
    %
    %   L samples; tau seconds; fs Hz. Storage is 16*L^2 bytes per path (double
    %   complex), construction takes L DELAYSEQ calls, and application is O(L^2).
    %   This reference backend targets the small capsule, not paper BER timing.
    %   Tolerance (2e-12), ValidationCount (4), Seed (20260917) are name-values.
    %   Validation restores the caller's random-number state.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    validateattributes(L, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
    validateattributes(tau, {'numeric'}, {'scalar', 'real', 'finite'});
    validateattributes(fs, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    parser = inputParser;
    parser.addParameter('Tolerance', 2e-12, @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x > 0);
    parser.addParameter('ValidationCount', 4, @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    parser.addParameter('Seed', 20260917, @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 0 && x == floor(x));
    parser.parse(varargin{:});
    options = parser.Results;

    delayMatrix = complex(zeros(L, L));
    for column = 1:L
        probe = complex(zeros(L, 1));
        probe(column) = 1j;
        response = delayseq(probe, tau, fs);
        assert(isequal(size(response), [L 1]), 'MDDIPIC:DelayShape', ...
               'Installed delayseq returned an unexpected response shape.');
        delayMatrix(:, column) = response / 1j;
    end

    previous = rng;
    cleanup = onCleanup(@() rng(previous)); %#ok<NASGU>
    rng(options.Seed, 'twister');
    forwardErrors = zeros(options.ValidationCount, 1);
    for trial = 1:options.ValidationCount
        x = (randn(L, 1) + 1j * randn(L, 1)) / sqrt(2);
        reference = delayseq(x, tau, fs);
        forwardErrors(trial) = norm(delayMatrix * x - reference) / max(norm(reference), realmin);
    end
    y = (randn(L, 1) + 1j * randn(L, 1)) / sqrt(2);
    lhs = (delayMatrix * x)' * y;
    rhs = x' * (delayMatrix' * y);
    adjointError = abs(lhs - rhs) / max([abs(lhs), abs(rhs), realmin]);
    maxForwardError = max(forwardErrors);
    assert(isfinite(maxForwardError) && maxForwardError <= options.Tolerance && ...
           isfinite(adjointError) && adjointError <= 2e-12, ...
           'MDDIPIC:DelayValidation', 'Installed delayseq parity %.3e; adjoint %.3e.', ...
           maxForwardError, adjointError);

    kernel.forward = @(x) apply_delay_matrix(delayMatrix, x, L);
    kernel.adjoint = @(y) apply_delay_matrix(delayMatrix', y, L);
    kernel.info = struct('schemaVersion', "DELAYSEQ_REFERENCE_KERNEL_V1", ...
                         'backend', "Installed delayseq basis responses", ...
                         'delaySamples', tau * fs, 'matrixBytes', 16 * L^2, ...
                         'validationRelativeErrors', forwardErrors, ...
                         'maxForwardRelativeError', maxForwardError, ...
                         'adjointRelativeError', adjointError, 'pass', true);
end

function output = apply_delay_matrix(matrix, input, L)
    assert(isvector(input) && numel(input) == L, 'MDDIPIC:DelayInput', ...
           'Delay input must be a vector containing L samples.');
    output = matrix * input(:);
end
