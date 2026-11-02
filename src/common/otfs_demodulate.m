function [grid, timeFrequency] = otfs_demodulate(samples, M, paddingLength, offset, paddingMode)
    % OTFS_DEMODULATE Rectangular-pulse Zak transform after removing padding.
    %   SAMPLES is one serialized frame, M is the number of delay rows.
    %   OFFSET counts samples skipped at the start of each padded subsymbol
    %   (CP/ZP), or at the start of the full frame (RCP/RZP). Use the prefix
    %   length for CP/RCP and zero for ZP/RZP. NONE requires offset zero.
    %   GRID = M*fft(timeGrid,[],2); TIMEFREQUENCY = fft(timeGrid,[],1).
    %   The explicit FFT dimensions also support a single delay/Doppler bin.
    %   See THIRD_PARTY_NOTICES.md for references and mathematical provenance.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    if nargin < 5
        paddingMode = 'CP';
    end
    validateattributes(samples, {'double', 'single'}, {'vector', 'nonempty', 'finite'}, mfilename, 'samples');
    validateattributes(M, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'M');
    validateattributes(paddingLength, {'numeric'}, ...
                       {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, mfilename, 'paddingLength');
    validateattributes(offset, {'numeric'}, ...
                       {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, mfilename, 'offset');
    paddingMode = validatestring(paddingMode, {'CP', 'ZP', 'RCP', 'RZP', 'NONE'});
    samples = samples(:);
    switch paddingMode
        case {'CP', 'ZP'}
            N = numel(samples) / (M + paddingLength);
            assert(N >= 1 && N == floor(N) && offset <= paddingLength, ...
                   'MDDIPIC:FrameShape', 'Invalid padded frame length or offset.');
            padded = reshape(samples, M + paddingLength, N);
            timeGrid = padded(offset + (1:M), :);
        case {'RCP', 'RZP'}
            N = (numel(samples) - paddingLength) / M;
            assert(N >= 1 && N == floor(N) && offset <= paddingLength, ...
                   'MDDIPIC:FrameShape', 'Invalid reduced-padding frame length or offset.');
            timeGrid = reshape(samples(offset + (1:M * N)), M, N);
        case 'NONE'
            N = numel(samples) / M;
            assert(N >= 1 && N == floor(N) && offset == 0, ...
                   'MDDIPIC:FrameShape', 'Unpadded frame length must be M*N and offset must be zero.');
            timeGrid = reshape(samples, M, N);
    end
    grid = M * fft(timeGrid, [], 2);
    if nargout > 1
        timeFrequency = fft(timeGrid, [], 1);
    end
end
