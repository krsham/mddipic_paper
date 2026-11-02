function [samples, timeFrequency] = otfs_modulate(grid, paddingLength, paddingMode)
    % OTFS_MODULATE Rectangular-pulse inverse Zak transform with frame padding.
    %   GRID is M-by-N: delay rows, Doppler columns. SAMPLES is a column vector.
    %   The normalization is T = ifft(GRID,[],2)/M; therefore, without padding,
    %   norm(GRID,'fro')^2 = M^2*N*norm(SAMPLES)^2.
    %   TIMEFREQUENCY is fft(T,[],1), before padding.
    %
    %   Modes: CP/ZP add a prefix/zeros to each M-sample subsymbol; RCP/RZP
    %   add one prefix/zero block to the full serialized frame; NONE adds none.
    %   Default mode is CP. MATLAB column order serializes successive columns.
    %   These equations follow the rectangular-pulse inverse Zak transform;
    %   see THIRD_PARTY_NOTICES.md for mathematical references and provenance.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    if nargin < 3
        paddingMode = 'CP';
    end
    validateattributes(grid, {'double', 'single'}, {'2d', 'nonempty', 'finite'}, mfilename, 'grid');
    validateattributes(paddingLength, {'numeric'}, ...
                       {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, mfilename, 'paddingLength');
    paddingMode = validatestring(paddingMode, {'CP', 'ZP', 'RCP', 'RZP', 'NONE'});
    [M, ~] = size(grid);
    timeGrid = ifft(grid, [], 2) / M;
    if nargout > 1
        timeFrequency = fft(timeGrid, [], 1);
    end
    switch paddingMode
        case 'CP'
            assert(paddingLength <= M, 'MDDIPIC:PaddingLength', 'CP length must not exceed M.');
            padded = [timeGrid(M - paddingLength + 1:M, :); timeGrid];
            samples = padded(:);
        case 'ZP'
            padded = [timeGrid; zeros(paddingLength, size(grid, 2), 'like', timeGrid)];
            samples = padded(:);
        case 'RCP'
            frame = timeGrid(:);
            assert(paddingLength <= numel(frame), 'MDDIPIC:PaddingLength', ...
                   'RCP length must not exceed the unpadded frame length.');
            samples = [frame(end - paddingLength + 1:end); frame];
        case 'RZP'
            samples = [timeGrid(:); zeros(paddingLength, 1, 'like', timeGrid)];
        case 'NONE'
            samples = timeGrid(:);
    end
end
