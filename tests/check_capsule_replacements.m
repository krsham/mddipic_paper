function checks = check_capsule_replacements()
    % CHECK_CAPSULE_REPLACEMENTS Check the replacement transforms and delay backend.
    %   Direct DFT sums provide an independent transform reference. Padding modes,
    %   one-row/one-column grids, integer/fractional signed delays, zero response,
    %   and adjoints are covered without running an estimator experiment.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    previous = rng;
    cleanup = onCleanup(@() rng(previous)); %#ok<NASGU>
    rng(2601103, 'twister');
    names = strings(0, 1);
    errors = zeros(0, 1);
    modes = {'NONE', 'CP', 'ZP', 'RCP', 'RZP'};
    shapes = [8 4; 1 4; 4 1];
    for shape = 1:size(shapes, 1)
        M = shapes(shape, 1);
        N = shapes(shape, 2);
        grid = randn(M, N) + 1j * randn(M, N);
        % Positive-exponent sum, divided by N, is MATLAB's inverse DFT.
        inverseDFT = exp(2j * pi * (0:N - 1).' * (0:N - 1) / N) / N;
        directTime = (grid * inverseDFT) / M;
        forwardDFT = exp(-2j * pi * (0:M - 1).' * (0:M - 1) / M);
        directTF = forwardDFT * directTime;
        for mode = modes
            for paddingLength = unique([0 min(2, M)])
                offset = 0;
                if any(strcmp(mode{1}, {'CP', 'RCP'}))
                    offset = paddingLength;
                end
                % Inspect the transmitted padding, not only what demodulation keeps.
                switch mode{1}
                    case 'CP'
                        packed = [directTime(M - paddingLength + 1:M, :); directTime];
                        expectedSignal = packed(:);
                    case 'ZP'
                        packed = [directTime; zeros(paddingLength, N)];
                        expectedSignal = packed(:);
                    case 'RCP'
                        frame = directTime(:);
                        expectedSignal = [frame(end - paddingLength + 1:end); frame];
                    case 'RZP'
                        expectedSignal = [directTime(:); zeros(paddingLength, 1)];
                    case 'NONE'
                        expectedSignal = directTime(:);
                end
                [signal, tfMod] = helperOTFSmod(grid, paddingLength, mode{1});
                [recovered, tfDemod] = helperOTFSdemod(signal, M, paddingLength, offset, mode{1});
                assert(isequal(size(signal), size(expectedSignal)), ...
                       'MDDIPIC:PaddingShape', 'Transmitted padding has the wrong size.');
                err = max([norm(signal - expectedSignal) / norm(expectedSignal), ...
                           norm(recovered - grid, 'fro') / norm(grid, 'fro'), ...
                           norm(tfMod - directTF, 'fro') / norm(directTF, 'fro'), ...
                           norm(tfDemod - directTF, 'fro') / norm(directTF, 'fro')]);
                if strcmp(mode{1}, 'NONE')
                    err = max(err, norm(signal - directTime(:)) / norm(directTime(:)));
                end
                assert(err < 2e-12, 'MDDIPIC:TransformCheck', ...
                       'Transform/padding check failed for %dx%d %s pad %d: %.3e.', ...
                       M, N, mode{1}, paddingLength, err);
                names(end + 1, 1) = sprintf('OTFS %dx%d %s pad %d', M, N, mode{1}, paddingLength); %#ok<AGROW>
                errors(end + 1, 1) = err; %#ok<AGROW>
            end
        end
    end

    % Prime frame length and near-boundary delays exercise finite-frame behavior.
    L = 17;
    fs = 10e6;
    delays = [0 2 -2 0.37 -0.37 2.18 -2.18 L + 2 -L - 2];
    for delaySamples = delays
        kernel = makeDelayseqReferenceKernel_v1(L, delaySamples / fs, fs);
        err = max(kernel.info.maxForwardRelativeError, kernel.info.adjointRelativeError);
        names(end + 1, 1) = sprintf('delayseq %.3g samples', delaySamples); %#ok<AGROW>
        errors(end + 1, 1) = err; %#ok<AGROW>
    end
    checks = table(names, errors, 'VariableNames', {'check', 'relativeError'});
end
