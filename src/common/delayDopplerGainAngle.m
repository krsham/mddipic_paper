function [y_sum, y_cell] = delayDopplerGainAngle(signal, delays_sec, dopplers_Hz, gains, angles_rad, fs, Na, delay_method)
    % DELAYDOPPLERGAINANGLE Apply Doppler, delay, gain, and ULA steering.
    %
    % For path p, the temporal signal follows
    %   y_p(t) = alpha_p x(t-tau_p) exp(j*2*pi*nu_p*(t-tau_p)),
    % and is spatially expanded using steeringVector(theta_p, Na).
    %
    % Inputs:
    %   signal        - N x 1 time-domain signal vector
    %   delays_sec    - P x 1 path delays [s]
    %   dopplers_Hz   - P x 1 path Doppler shifts [Hz]
    %   gains         - P x 1 complex path gains
    %   angles_rad    - P x 1 AoA/AoD angles [rad]
    %   fs            - Sampling frequency [Hz]
    %   Na            - Number of array elements
    %   delay_method  - Optional delay implementation:
    %                     'delayseq' (DEFAULT)
    %                     'farrow'   -> apply_farrow_delay_1d
    %
    % Outputs:
    %   y_sum         - N x Na summed received array signal
    %   y_cell        - Optional P x 1 cell array; y_cell{p} is N x Na
    %
    % Examples:
    %   % Default: MATLAB delayseq
    %   Y = delayDopplerGainAngle(x, tau, nu, alpha, theta, fs, Na);
    %
    %   % Explicit delayseq
    %   Y = delayDopplerGainAngle(x, tau, nu, alpha, theta, fs, Na, 'delayseq');
    %
    %   % Existing fractional-delay helper
    %   Y = delayDopplerGainAngle(x, tau, nu, alpha, theta, fs, Na, 'farrow');
    %
    % Notes:
    %   - delayseq is the default delay backend.
    %   - The spatial response uses a nonconjugating transpose:
    %       y_p = sig_gain * steeringVector(theta_p,Na).'
    %     i.e., non-conjugating transpose on the steering vector.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    if nargin < 8 || isempty(delay_method)
        delay_method = 'delayseq';
    end

    delay_method = normalize_delay_method(delay_method);

    % Ensure column vectors.
    signal      = signal(:);
    delays_sec  = delays_sec(:);
    dopplers_Hz = dopplers_Hz(:);
    gains       = gains(:);
    angles_rad  = angles_rad(:);

    P = numel(delays_sec);
    if numel(dopplers_Hz) ~= P || numel(gains) ~= P || numel(angles_rad) ~= P
        error('delayDopplerGainAngle:PathCountMismatch', ...
              ['delays_sec, dopplers_Hz, gains, and angles_rad must contain ' ...
               'the same number of paths.']);
    end
    if ~isscalar(fs) || ~isfinite(fs) || fs <= 0
        error('delayDopplerGainAngle:InvalidFs', ...
              'fs must be a positive finite scalar.');
    end
    if ~isscalar(Na) || ~isfinite(Na) || Na < 1 || Na ~= floor(Na)
        error('delayDopplerGainAngle:InvalidNa', ...
              'Na must be a positive integer scalar.');
    end

    N = numel(signal);
    t_vec = (0:N - 1).' / fs;

    y_cell = cell(P, 1);
    y_sum  = zeros(N, Na);

    for p = 1:P
        %% 1) Doppler first
        phase_shift = exp(1j * 2 * pi * dopplers_Hz(p) * t_vec);
        sig_dop = phase_shift .* signal;

        %% 2) Delay second
        switch delay_method
            case 'delayseq'
                % delayseq(x, delay_seconds, fs)
                sig_del = delayseq(sig_dop, delays_sec(p), fs);
                sig_del = sig_del(:);  % Preserve project column convention.

            case 'farrow'
                sig_del = apply_farrow_delay_1d(sig_dop, delays_sec(p), fs);

            otherwise
                % normalize_delay_method prevents this branch.
                error('delayDopplerGainAngle:InternalDelayMethod', ...
                      'Unexpected normalized delay method ''%s''.', delay_method);
        end

        %% 3) Complex path gain
        sig_gain = gains(p) * sig_del;

        %% 4) Spatial steering
        a_theta = steeringVector(angles_rad(p), Na);  % Na x 1

        % Nonconjugating spatial expansion:
        % (N x 1) * (1 x Na) -> N x Na using non-conjugating transpose.
        y_p = sig_gain * a_theta.';

        y_cell{p} = y_p;
        y_sum = y_sum + y_p;
    end

    if nargout < 2
        y_cell = [];
    end
end

function method = normalize_delay_method(method)
    % NORMALIZE_DELAY_METHOD Normalize supported delay-method aliases.

    if ~(ischar(method) || (isstring(method) && isscalar(method)))
        error('delayDopplerGainAngle:InvalidDelayMethod', ...
              'delay_method must be ''delayseq'' or ''farrow''.');
    end

    method = lower(strtrim(char(method)));

    switch method
        case {'delayseq', 'delay-seq'}
            method = 'delayseq';

        case {'farrow', 'apply_farrow', 'apply_farrow_delay_1d'}
            method = 'farrow';

        otherwise
            error('delayDopplerGainAngle:InvalidDelayMethod', ...
                  ['Unsupported delay_method ''%s''. Use ''delayseq'' ' ...
                   '(default) or ''farrow''.'], method);
    end
end
