function y = apply_farrow_delay_1d(x, delay_sec, fs)
    % APPLY_FARROW_DELAY_1D Legacy three-tap fractional-delay approximation.
    %
    % Y = apply_farrow_delay_1d(X,delay_sec,fs) returns a vector of length numel(X).
    % The integer delay prepends zeros. The fractional part uses three sinc taps
    % and compensates their one-sample causal-filter shift. Despite the compatibility
    % name, this helper uses sinc taps. The demo uses delayseq.
    % This fallback is retained because it is exposed by delayDopplerGainAngle;
    % it is not interchangeable with the finite FFT delay operator.
    %
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    x = x(:); % Ensure column vector
    N = length(x);
    D_total = delay_sec * fs;
    D_int = floor(D_total);
    d = D_total - D_int;

    %% 1. Manual Integer Shifting (Causal, prepend zeros)
    if D_int >= N
        y = zeros(N, 1);
        return;
    elseif D_int > 0
        x_shifted = [zeros(D_int, 1); x(1:N - D_int)];
    else
        x_shifted = x;
    end

    %% 2. Three-tap fractional-delay filter

    B = sinc((-1:1) - d);
    %% 3. GROUP DELAY COMPENSATION
    % The mathematical filter uses n-1, n, n+1 (non-causal).
    % MATLAB's filter() assumes n, n-1, n-2 (causal), adding a 1-sample delay.
    % To fix this, we pad the input with 1 zero at the end, filter, and
    % shift the output left by 1 sample.
    x_pad = [x_shifted; 0];          % Pad input to evaluate the non-causal tap at n=N
    y_pad = filter(B, 1, x_pad);     % Apply causal filter (length N+1)
    y = y_pad(2:end);                % Shift left by 1 to cancel the 1-sample bulk delay (length N)
end
