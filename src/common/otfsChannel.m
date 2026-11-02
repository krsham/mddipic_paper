function [y_dd, tx_sig, rx_sig_clean, rx_sig_noisy] = otfsChannel(x_dd, M, fs, delays, dopplers, gains, noise_opt, padtype, padlen)
    % OTFSCHANNEL Simulates the end-to-end OTFS transceiver over a delay-Doppler channel.
    %
    % Inputs:
    %   x_dd        - M x N matrix of Delay-Doppler domain symbols.
    %   M           - Number of subcarriers (Delay bins).
    %   delta_f     - Subcarrier spacing in Hz.
    %   fs          - Sampling frequency in Hz (typically M * delta_f).
    %   delays      - Px1 vector of path delays in seconds.
    %   dopplers    - Px1 vector of path Doppler shifts in Hz.
    %   gains       - Px1 vector of complex path gains.
    %   noise_opt   - Noise configuration. Can be:
    %                   - Scalar (SNR in dB): Adds AWGN using MATLAB's awgn() with 'measured' power.
    %                   - Vector (Nx1): Adds the provided external noise vector directly.
    %                   - String 'none': No noise added.
    %   padtype     - OTFS padding type: 'CP', 'ZP', 'RCP', 'RZP', or 'NONE'.
    %   padlen      - Length of the CP/ZP in samples.
    %
    % Outputs:
    %   y_dd        - M x N matrix of received Delay-Doppler domain symbols.
    %   tx_sig      - Transmitted time-domain signal (column vector).
    %   rx_sig_clean- Received time-domain signal before noise (column vector).
    %   rx_sig_noisy- Received time-domain signal after noise (column vector).

    %% 1. OTFS Modulation
    % Modulate the M x N DD domain grid into the time domain
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    [tx_sig, ~] = helperOTFSmod(x_dd, padlen, padtype);

    %% 2. Time-Domain Channel Application
    % Apply Doppler, finite-frame delay and path gain.
    rx_sig_clean = delayDopplerGain(tx_sig, delays, dopplers, gains, fs);

    %% 3. Noise Addition
    if ischar(noise_opt) && strcmpi(noise_opt, 'none')
        % No noise
        rx_sig_noisy = rx_sig_clean;
    elseif isnumeric(noise_opt)
        % noise_opt is SNR in dB. Use awgn with 'measured' to handle complex baseband signals correctly.
        rx_sig_noisy = awgn(rx_sig_clean, noise_opt, 'measured');
    elseif isvector(noise_opt)
        % noise_opt is an external noise vector
        if length(noise_opt) ~= length(rx_sig_clean)
            error('External noise vector length (%d) must match signal length (%d).', ...
                  length(noise_opt), length(rx_sig_clean));
        end
        rx_sig_noisy = rx_sig_clean + noise_opt(:);

    else
        error('Invalid noise_opt. Provide SNR in dB, an external noise vector, or ''none''.');
    end

    %% 4. OTFS Demodulation
    % Determine the correct sample offset for the demodulator based on padtype
    switch upper(padtype)
        case {'CP', 'RCP'}
            offset = padlen; % Skip the prepended Cyclic Prefix
        case {'ZP', 'RZP', 'NONE'}
            offset = 0;      % Zero padding is at the end, or no padding exists
        otherwise
            error('Invalid padtype. Must be CP, ZP, RCP, RZP, or NONE.');
    end

    % Demodulate the noisy time-domain signal back to the M x N DD domain
    [y_dd, ~] = helperOTFSdemod(rx_sig_noisy, M, padlen, offset, padtype);
end
