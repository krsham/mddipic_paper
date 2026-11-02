function [y_sum, y_cell] = delayDopplerGain(signal, delays_sec, dopplers_Hz, gains, fs, delay_method)
    % DELAYDOPPLERGAIN Apply Doppler, delay, and gain to a time-domain signal.
    %
    % Physical model for path p:
    %   y_p(t) = alpha_p x(t-tau_p) exp(j*2*pi*nu_p*(t-tau_p))
    %
    % The implementation applies Doppler first and delay second:
    %   x_dop(t)      = x(t) exp(j*2*pi*nu_p*t)
    %   x_dd(t)       = x_dop(t-tau_p)
    %                 = x(t-tau_p) exp(j*2*pi*nu_p*(t-tau_p))
    %   y_p(t)        = alpha_p x_dd(t)
    %
    % Inputs:
    %   signal        - N x 1 time-domain signal vector
    %   delays_sec    - P x 1 path delays [s]
    %   dopplers_Hz   - P x 1 path Doppler shifts [Hz]
    %   gains         - P x 1 complex path gains
    %   fs            - Sampling frequency [Hz]
    %   delay_method  - Optional delay implementation:
    %                     'delayseq' (DEFAULT)
    %                     'farrow'   -> apply_farrow_delay_1d
    %
    % Outputs:
    %   y_sum         - N x 1 sum of all received path signals
    %   y_cell        - Optional P x 1 cell array; y_cell{p} is N x 1
    %
    % Examples:
    %   % Default: MATLAB delayseq
    %   y = delayDopplerGain(x, tau, nu, alpha, fs);
    %
    %   % Explicit delayseq
    %   y = delayDopplerGain(x, tau, nu, alpha, fs, 'delayseq');
    %
    %   % Existing fractional-delay helper
    %   y = delayDopplerGain(x, tau, nu, alpha, fs, 'farrow');
    %
    % Notes:
    %   - delayseq is the default delay backend.
    %   - The optional three-tap approximation is apply_farrow_delay_1d.m.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    if nargin < 6 || isempty(delay_method)
        delay_method = 'delayseq';
    end

    delay_method = normalize_delay_method(delay_method);

    % Ensure column vectors.
    signal      = signal(:);
    delays_sec  = delays_sec(:);
    dopplers_Hz = dopplers_Hz(:);
    gains       = gains(:);

    P = numel(delays_sec);
    if numel(dopplers_Hz) ~= P || numel(gains) ~= P
        error('delayDopplerGain:PathCountMismatch', ...
              ['delays_sec, dopplers_Hz, and gains must contain the same ' ...
               'number of paths.']);
    end
    if ~isscalar(fs) || ~isfinite(fs) || fs <= 0
        error('delayDopplerGain:InvalidFs', ...
              'fs must be a positive finite scalar.');
    end

    N = numel(signal);
    t_vec = (0:N - 1).' / fs;

    y_cell = cell(P, 1);
    y_sum  = zeros(N, 1);

    for p = 1:P
        %% 1) Doppler first
        phase_shift = exp(1j * 2 * pi * dopplers_Hz(p) * t_vec);
        x_dop = phase_shift .* signal;

        %% 2) Delay second
        switch delay_method
            case 'delayseq'
                % delayseq(x, delay_seconds, fs)
                x_dd = delayseq(x_dop, delays_sec(p), fs);
                x_dd = x_dd(:);  % Preserve project column-vector convention.

            case 'farrow'
                x_dd = apply_farrow_delay_1d(x_dop, delays_sec(p), fs);

            otherwise
                % normalize_delay_method prevents this branch.
                error('delayDopplerGain:InternalDelayMethod', ...
                      'Unexpected normalized delay method ''%s''.', delay_method);
        end

        %% 3) Complex path gain
        y_p = gains(p) * x_dd;

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
        error('delayDopplerGain:InvalidDelayMethod', ...
              'delay_method must be ''delayseq'' or ''farrow''.');
    end

    method = lower(strtrim(char(method)));

    switch method
        case {'delayseq', 'delay-seq'}
            method = 'delayseq';

        case {'farrow', 'apply_farrow', 'apply_farrow_delay_1d'}
            method = 'farrow';

        otherwise
            error('delayDopplerGain:InvalidDelayMethod', ...
                  ['Unsupported delay_method ''%s''. Use ''delayseq'' ' ...
                   '(default) or ''farrow''.'], method);
    end
end
