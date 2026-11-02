function [y_dd_mimo, tx_sig, rx_sig_clean, rx_sig_noisy] = otfsMimoChannel(x_dd, M, fs, delays_sec, dopplers_Hz, gains, angles_rad, Na, noise_opt, padtype, padlen)
    % OTFSMIMOCHANNEL Modulate, propagate and demodulate a receive-array OTFS frame.
    %
    % [YDD,tx,cleanTD,noisyTD] = otfsMimoChannel(XDD,M,fs,tau,nu,g,theta,Na,noise,padtype,padlen)
    %   XDD: M-by-N symbols. tau [s], nu [Hz], g and theta [rad] have P elements.
    %   noise: 'none', scalar per-antenna SNR [dB], or an L-by-Na noise matrix.
    %   YDD: M-by-N-by-Na. cleanTD/noisyTD: L-by-Na. tx: L-by-1.
    %
    % Each path applies Doppler, fractional delay, gain, then nonconjugated
    % receive steering. The demo uses padtype='NONE', padlen=0.
    % An explicit noise matrix keeps every estimator on identical observations.
    % The scalar-SNR mode uses MathWorks awgn; see THIRD_PARTY_NOTICES.md.
    %
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    [tx_sig, ~] = helperOTFSmod(x_dd, padlen, padtype);
    N_tx = length(tx_sig);

    %% 2. Time-Domain MIMO Channel Application
    rx_sig_clean = delayDopplerGainAngle(tx_sig, delays_sec, dopplers_Hz, ...
                                         gains, angles_rad, fs, Na);

    %% 3. Noise addition
    rx_sig_noisy = zeros(N_tx, Na);

    % Handle text before matrices because MATLAB character arrays are matrices.
    if ischar(noise_opt) || isstring(noise_opt)
        if strcmpi(noise_opt, 'none')
            rx_sig_noisy = rx_sig_clean;
        else
            error('Invalid string for noise_opt. Use ''none''.');
        end

    elseif isnumeric(noise_opt) && isscalar(noise_opt)
        % noise_opt is a scalar SNR in dB. Apply AWGN independently per antenna.
        for a = 1:Na
            rx_sig_noisy(:, a) = awgn(rx_sig_clean(:, a), noise_opt, 'measured');
        end

    elseif isnumeric(noise_opt) && ~isscalar(noise_opt)
        % noise_opt is an external noise matrix
        if ~isequal(size(noise_opt), [N_tx, Na])
            error('External noise matrix size [%d x %d] must match signal size [%d x %d].', ...
                  size(noise_opt, 1), size(noise_opt, 2), N_tx, Na);
        end
        rx_sig_noisy = rx_sig_clean + noise_opt;

    else
        error('Invalid noise_opt. Provide SNR in dB, an N_tx x Na noise matrix, or ''none''.');
    end

    %% 4. OTFS Demodulation (Time Domain -> DD Domain)
    switch upper(padtype)
        case {'CP', 'RCP'}
            offset = padlen;
        case {'ZP', 'RZP', 'NONE'}
            offset = 0;
        otherwise
            error('Invalid padtype. Must be CP, ZP, RCP, RZP, or NONE.');
    end

    N_sym = size(x_dd, 2);
    y_dd_mimo = zeros(M, N_sym, Na);

    for a = 1:Na
        [y_dd_mimo(:, :, a), ~] = helperOTFSdemod(rx_sig_noisy(:, a), M, padlen, offset, padtype);
    end
end
