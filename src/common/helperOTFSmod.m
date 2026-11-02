function [samples, timeFrequency] = helperOTFSmod(grid, paddingLength, varargin)
    % HELPEROTFSMOD Compatibility name used by the paper estimators.
    %   Delegates to the project implementation; no MathWorks example is bundled.
    %   See OTFS_MODULATE for normalization, dimensions and padding conventions.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    [samples, timeFrequency] = otfs_modulate(grid, paddingLength, varargin{:});
end
