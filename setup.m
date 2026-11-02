function environment = setup()
    % SETUP Add this package to the path and check the required MATLAB features.
    %   ENVIRONMENT = SETUP() adds only the source and demo folders below this
    %   file. It does not clear your workspace, change folders, download code,
    %   start a parallel pool, or run a simulation.
    %
    %   MATLAB R2025b is the reference release. Communications Toolbox supplies
    %   qammod/qamdemod; Phased Array System Toolbox supplies delayseq.
    %   Other releases are accepted, but the finite-delay adjoint check must pass.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    root = fileparts(mfilename('fullpath'));
    folders = {'src/common', 'src/proposed', 'src/baselines', 'src/ber', ...
               'src/scenario', 'src/evaluation', 'demo', 'tests'};
    addpath(root);
    for k = 1:numel(folders)
        addpath(fullfile(root, folders{k}));
    end

    required = {'delayseq', 'qammod', 'qamdemod', 'pcg', 'lsqr', ...
                'helperOTFSmod', 'helperOTFSdemod'};
    missing = required(cellfun(@(name) isempty(which(name)), required));
    if ~isempty(missing)
        error('MDDIPIC:MissingDependency', ...
              'Missing functions: %s. See README.md.', strjoin(missing, ', '));
    end
    if isempty(meta.class.fromName('matlab.System'))
        error('MDDIPIC:MissingSystemObjects', 'MATLAB System objects are required.');
    end

    % A second checkout on the path can silently select a different estimator.
    ownFunctions = {'helperOTFSmod', 'helperOTFSdemod', 'MDDIPICEstimator', ...
                    'DRAPMDDIPICEstimator', 'Li3DNOMPCommonHostEstimator'};
    for k = 1:numel(ownFunctions)
        resolved = which(ownFunctions{k});
        if ~startsWith(resolved, [root filesep], 'IgnoreCase', ispc)
            error('MDDIPIC:PathConflict', ...
                  '%s resolves outside this package: %s', ownFunctions{k}, resolved);
        end
    end

    environment = struct('root', root, 'matlabVersion', version, ...
                         'release', version('-release'), 'platform', computer, ...
                         'toolboxes', ver, 'validation', 'runtime checks pending');
end
