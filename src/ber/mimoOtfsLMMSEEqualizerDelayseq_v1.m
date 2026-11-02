function [Xhat, info] = mimoOtfsLMMSEEqualizerDelayseq_v1( ...
                                                         rxTD, M, N, op, sigma2TD, Px, varargin)
    % MIMOOTFSLMMSEEQUALIZERDELAYSEQ_V1 Iterative SIMO OTFS LMMSE equalizer.
    %
    % Solves full-load LMMSE in time-domain coordinates using the finite-frame
    % delayseq physical operator:
    %
    %   min_t ||y-A t||^2 + lambda ||t||^2,
    %   lambda = sigma2TD * (M^2*N) / Px.
    %
    % Since helperOTFSmod with PADTYPE='NONE' obeys
    %   E[t t^H] = Px/(M^2*N) I,
    % this is equivalent to DD-domain LMMSE for unit-power DD symbols.  The
    % recovered TD waveform is finally mapped back to DD by helperOTFSdemod.
    %
    % Inputs
    %   rxTD      [MN x Na]
    %   M,N       OTFS dimensions
    %   op        makeMimoDelayseqOperator_v1 output
    %   sigma2TD  per-complex-sample AWGN variance at each antenna
    %   Px        DD symbol variance (1 for UnitAveragePower QAM)
    %
    % Name-value
    %   Tol       1e-7
    %   MaxIter   250
    %   X0TD      [] optional warm start
    %   Solver    'pcg' (default) or 'lsqr'
    %
    % Output
    %   Xhat      [M x N] estimated DD symbols
    %   info      solver diagnostics

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    p = inputParser;
    p.addParameter('Tol', 1e-7, @(x)isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('MaxIter', 250, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('X0TD', [], @(x)isnumeric(x));
    p.addParameter('Solver', 'pcg', @(x)ischar(x) || isstring(x));
    p.parse(varargin{:});
    o = p.Results;

    L = M * N;
    assert(op.L == L, 'Operator length must equal M*N.');
    assert(size(rxTD, 1) == L && size(rxTD, 2) == op.Na, 'rxTD size mismatch.');
    assert(sigma2TD >= 0 && Px > 0);

    scale = M^2 * N;
    lambda = sigma2TD * scale / Px;
    rhs = op.adjoint(rxTD);

    if isempty(o.X0TD)
        x0 = complex(zeros(L, 1));
    else
        x0 = o.X0TD(:);
        assert(numel(x0) == L);
    end

    normal = @(x) op.adjoint(op.forward(x)) + lambda * x;
    solver = lower(string(o.Solver));

    switch solver
        case "pcg"
            [tHat, flag, relres, iter, resvec] = pcg(normal, rhs, o.Tol, o.MaxIter, [], [], x0);
            usedFallback = false;
            primaryFlag = flag;

            % If PCG does not finish, solve the equivalent augmented LS problem.
            if flag ~= 0
                b = [rxTD(:); complex(zeros(L, 1))];
                [t2, flag2, relres2, iter2, resvec2] = lsqr( ...
                                                      @(x, tr) augop(x, tr, op, L, lambda), b, o.Tol, max(2 * o.MaxIter, 400));
                if flag2 == 0 || relres2 < relres
                    tHat = t2;
                    flag = flag2;
                    relres = relres2;
                    iter = iter2;
                    resvec = resvec2;
                    solver = "lsqr-fallback";
                    usedFallback = true;
                end
            end

        case "lsqr"
            b = [rxTD(:); complex(zeros(L, 1))];
            [tHat, flag, relres, iter, resvec] = lsqr( ...
                                                @(x, tr) augop(x, tr, op, L, lambda), b, o.Tol, o.MaxIter);
            primaryFlag = flag;
            usedFallback = false;

        otherwise
            error('Unsupported Solver "%s".', solver);
    end

    Xhat = helperOTFSdemod(tHat, M, 0, 0, 'NONE');

    nr = norm(normal(tHat) - rhs) / max(norm(rhs), realmin);
    info = struct();
    info.schemaVersion = "MIMO_OTFS_LMMSE_DELAYSEQ_V1";
    info.solver = solver;
    info.flag = flag;
    info.primaryFlag = primaryFlag;
    info.usedFallback = usedFallback;
    info.relres = relres;
    info.iter = iter;
    info.resvec = resvec;
    info.finalNormalResidual = nr;
    info.lambda = lambda;
    info.sigma2TD = sigma2TD;
    info.sigma2DD = scale * sigma2TD;
    info.xTD = tHat;
    info.operatorInfo = op.info;
end

function y = augop(x, transp_flag, op, L, lambda)
    if strcmp(transp_flag, 'notransp')
        A = op.forward(x);
        y = [A(:); sqrt(lambda) * x(:)];
    else
        nObs = L * op.Na;
        Y = reshape(x(1:nObs), L, op.Na);
        r = x(nObs + 1:end);
        y = op.adjoint(Y) + sqrt(lambda) * r;
    end
end
