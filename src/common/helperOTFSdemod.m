function [grid, timeFrequency] = helperOTFSdemod(samples, M, paddingLength, offset, varargin)
    % HELPEROTFSDEMOD Compatibility name used by the paper estimators.
    %   See OTFS_DEMODULATE for the original transform implementation and contract.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    [grid, timeFrequency] = otfs_demodulate(samples, M, paddingLength, offset, varargin{:});
end
