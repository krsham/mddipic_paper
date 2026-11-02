classdef DRAPMDDIPICEstimator < matlab.System
    % DRAPMDDIPICESTIMATOR Dimension-reduced alternating-projection OTFS estimator.
    %
    % [tau,nu,theta,g,P] = step(est,yDD,PExact,sigma2DD)
    %   yDD is (M*N*Na)-by-1 with contiguous antenna blocks.
    %   PExact>0 supplies the order; PExact=0 enables the selected stopping rule.
    %   Output coordinates are seconds, Hz and radians; g contains complex gains.
    %
    % For each coarse DD atom h, h'*Y reduces the MN-by-Na observation to an
    % Na-element spatial vector. An angle search is then applied to that vector.
    % Fine refinement alternates DD estimation after beamforming with spatial
    % estimation after temporal projection. NumAlternatingRefinements counts
    % complete DD/angle cycles; every cycle restarts the same grid schedule.
    %
    % Candidate scores use the residual and full atom norm. Accepted gains use
    % the joint fit controlled by rho; set rho=sigma2DD before setup for LMMSE.
    % This differs from MDDIPIC's support-conditioned candidate score. All active
    % paths are revisited after a multipath acquisition. DD-response caching is
    % local to one observation and changes computation, not the search.
    %
    % Dependencies: helperOTFSmod, helperOTFSdemod, delayseq, steeringVector.
    % See LIBRARY_STRUCTURE.md for the projection equations and interface.
    %
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    properties
        M = 16
        N = 16
        Na = 8
        fs = 150e6
        padtype = 'NONE'
        padlen = 0
        x_dd = []

        tau_grid_coarse = []
        nu_grid_coarse = []
        theta_grid_coarse = []
        theta_grid_fine = []

        m_tau = 10
        n_nu = 10
        max_iter_fine = 5
        % Compatibility fields; a center tie is not DD convergence.
        % All max_iter_fine levels are always evaluated.
        eps_tau = 1e-12
        eps_nu = 1e-12

        P_max = 10
        rho = 0

        % Automatic stopping. 'residual' uses residual change;
        % 'peakdrop' rejects the first candidate whose coarse DRAP score
        % drops by PeakDropThreshold_dB relative to the preceding path.
        stop_scale = 0.01
        AutoStopRule = 'residual'
        PeakDropThreshold_dB = 6

        % Number L of complete fine DD -> fine AoA alternating cycles.
        % L=1 performs one DD/angle refinement cycle.
        NumAlternatingRefinements = 1

        % Implementation-only ablation; disabling this does not alter the
        % candidate grid, scores, stopping rule, or generated responses.
        EnableDDResponseCache = true
    end

    properties (SetAccess = private)
        P_detected = 0
        LastResidualNorm = NaN
        LastStopMetric = NaN
        CandidateEvaluations = 0

        % Structural instrumentation used by unit/acceptance tests.
        PathSearchCalls = 0
        CoarseSearchCalls = 0
        FineDDSearchCalls = 0
        FineAngleSearchCalls = 0
        AlternatingCyclesExecuted = 0

        LastCoarseDDProjectionSize = [0 0]
        LastCoarseAngleScoreSize = [0 0]
        LastFineBeamformedSize = [0 0]
        LastFineSpatialProjectionSize = [0 0]
        LastCoarseDDIndex = 0
        LastCoarseThetaIndex = 0
        NewPathCoarsePeakScoreHistory = []
        DDResponseCacheHits = 0
        DDResponseCacheMisses = 0
    end

    properties (Access = private)
        MN
        d_tau
        d_nu

        txSig
        timeVec
        demodOffset

        PsiCoarse
        denDDCoarse
        AThetaCoarse
        denThetaCoarse
        AThetaFine
        denThetaFine

        Gtau
        Gnu
        Gdd
        GthetaCoarse
        GthetaFine
        DDResponseCache = []
    end

    methods

        function obj = DRAPMDDIPICEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.validateConfiguration();

            obj.MN = obj.M * obj.N;
            obj.d_tau = 1 / obj.fs;
            obj.d_nu = obj.fs / (obj.M * obj.N);

            [obj.txSig, ~] = helperOTFSmod( ...
                                          obj.x_dd, obj.padlen, upper(char(obj.padtype)));
            obj.txSig = obj.txSig(:);
            obj.timeVec = (0:numel(obj.txSig) - 1).' / obj.fs;
            obj.demodOffset = obj.getDemodOffset();

            obj.Gtau = numel(obj.tau_grid_coarse);
            obj.Gnu = numel(obj.nu_grid_coarse);
            obj.Gdd = obj.Gtau * obj.Gnu;
            obj.GthetaCoarse = numel(obj.theta_grid_coarse);
            obj.GthetaFine = numel(obj.theta_grid_fine);

            obj.PsiCoarse = complex(zeros(obj.MN, obj.Gdd));
            obj.denDDCoarse = zeros(obj.Gdd, 1);

            idx = 1;
            for j = 1:obj.Gnu
                for i = 1:obj.Gtau
                    h = obj.generateDDResponse( ...
                                               obj.tau_grid_coarse(i), obj.nu_grid_coarse(j));
                    obj.PsiCoarse(:, idx) = h;
                    obj.denDDCoarse(idx) = real(h' * h);
                    idx = idx + 1;
                end
            end

            obj.AThetaCoarse = complex(zeros(obj.Na, obj.GthetaCoarse));
            obj.denThetaCoarse = zeros(1, obj.GthetaCoarse);

            for k = 1:obj.GthetaCoarse
                a = steeringVector(obj.theta_grid_coarse(k), obj.Na);
                obj.AThetaCoarse(:, k) = a;
                obj.denThetaCoarse(k) = real(a' * a);
            end

            obj.AThetaFine = complex(zeros(obj.Na, obj.GthetaFine));
            obj.denThetaFine = zeros(1, obj.GthetaFine);

            for k = 1:obj.GthetaFine
                a = steeringVector(obj.theta_grid_fine(k), obj.Na);
                obj.AThetaFine(:, k) = a;
                obj.denThetaFine(k) = real(a' * a);
            end
        end

        function [tauEst, nuEst, thetaEst, gEst, Pest] = ...
                stepImpl(obj, yVec, PExact, noiseVarianceDD)

            yVec = yVec(:);

            if numel(yVec) ~= obj.MN * obj.Na
                error('DRAPMDDIPICEstimator:ObservationSize', ...
                      'yVec must contain M*N*Na = %d samples.', obj.MN * obj.Na);
            end

            if PExact < 0 || PExact ~= floor(PExact)
                error('DRAPMDDIPICEstimator:PExact', ...
                      'PExact must be a nonnegative integer.');
            end

            if noiseVarianceDD < 0
                error('DRAPMDDIPICEstimator:NoiseVariance', ...
                      'noiseVarianceDD must be nonnegative.');
            end

            if PExact > 0
                Plimit = PExact;
                autoStop = false;
            else
                Plimit = obj.P_max;
                autoStop = true;
            end

            Y = reshape(yVec, obj.MN, obj.Na);
            % The IPI sweep revisits many identical fractional DD points.
            % Keep generated responses only for this observation, so the
            % cache cannot leak across trials or change the search grid.
            if obj.EnableDDResponseCache
                obj.DDResponseCache = containers.Map( ...
                                                     'KeyType', 'char', 'ValueType', 'any');
            else
                obj.DDResponseCache = [];
            end
            obj.DDResponseCacheHits = 0;
            obj.DDResponseCacheMisses = 0;

            tauEst = zeros(Plimit, 1);
            nuEst = zeros(Plimit, 1);
            thetaEst = zeros(Plimit, 1);
            gEst = complex(zeros(Plimit, 1));

            Hactive = complex(zeros(obj.MN, Plimit));
            Aactive = complex(zeros(obj.Na, Plimit));

            obj.resetDiagnostics(yVec);
            obj.NewPathCoarsePeakScoreHistory = nan(Plimit, 1);
            residualPrev = yVec;

            for p = 1:Plimit
                if p == 1
                    Yres = Y;
                else
                    Yres = Y - obj.reconstructMatrix( ...
                                                     Hactive(:, 1:p - 1), Aactive(:, 1:p - 1), gEst(1:p - 1));
                end

                [tauP, nuP, thetaP, hP, aP, coarsePeak] = obj.estimateOnePath(Yres);
                obj.NewPathCoarsePeakScoreHistory(p) = coarsePeak;

                if autoStop && strcmpi(obj.AutoStopRule, 'peakdrop') && p > 1
                    qPrev = obj.NewPathCoarsePeakScoreHistory(p - 1);
                    qCurr = obj.NewPathCoarsePeakScoreHistory(p);
                    drop_dB = 10 * log10(max(qPrev, realmin) / max(qCurr, realmin));
                    if isfinite(drop_dB) && drop_dB >= obj.PeakDropThreshold_dB
                        obj.P_detected = p - 1;
                        break;
                    end
                end

                tauEst(p) = tauP;
                nuEst(p) = nuP;
                thetaEst(p) = thetaP;
                Hactive(:, p) = hP;
                Aactive(:, p) = aP;

                gActive = obj.estimateGains( ...
                                            Y, Hactive(:, 1:p), Aactive(:, 1:p));
                gEst(1:p) = gActive;

                residualCurrent = yVec - reshape(obj.reconstructMatrix( ...
                                                                       Hactive(:, 1:p), Aactive(:, 1:p), gActive), [], 1);

                obj.P_detected = p;
                obj.LastResidualNorm = norm(residualCurrent);

                if p > 1
                    obj.LastStopMetric = norm( ...
                                              residualCurrent - residualPrev)^2;
                else
                    obj.LastStopMetric = NaN;
                end

                if autoStop && strcmpi(obj.AutoStopRule, 'residual') && p > 1
                    epsilon = ...
                        obj.stop_scale * obj.Na * obj.MN * noiseVarianceDD;

                    if obj.LastStopMetric < epsilon
                        break;
                    end
                end

                if p > 1
                    for z = 1:p
                        keep = (1:p) ~= z;

                        YresZ = Y - obj.reconstructMatrix( ...
                                                          Hactive(:, keep), Aactive(:, keep), gEst(keep));

                        [tauZ, nuZ, thetaZ, hZ, aZ] = ...
                            obj.estimateOnePath(YresZ);

                        tauEst(z) = tauZ;
                        nuEst(z) = nuZ;
                        thetaEst(z) = thetaZ;
                        Hactive(:, z) = hZ;
                        Aactive(:, z) = aZ;

                        gActive = obj.estimateGains( ...
                                                    Y, Hactive(:, 1:p), Aactive(:, 1:p));
                        gEst(1:p) = gActive;
                    end

                    residualCurrent = yVec - reshape(obj.reconstructMatrix( ...
                                                                           Hactive(:, 1:p), Aactive(:, 1:p), gEst(1:p)), [], 1);

                    obj.LastResidualNorm = norm(residualCurrent);
                end

                residualPrev = residualCurrent;

                if ~autoStop && p == PExact
                    break;
                end
            end

            Pest = obj.P_detected;

            tauEst = tauEst(1:Pest);
            nuEst = nuEst(1:Pest);
            thetaEst = thetaEst(1:Pest);
            gEst = gEst(1:Pest);
        end

        function resetImpl(obj)
            obj.P_detected = 0;
            obj.LastResidualNorm = NaN;
            obj.LastStopMetric = NaN;
            obj.CandidateEvaluations = 0;

            obj.PathSearchCalls = 0;
            obj.CoarseSearchCalls = 0;
            obj.FineDDSearchCalls = 0;
            obj.FineAngleSearchCalls = 0;
            obj.AlternatingCyclesExecuted = 0;

            obj.LastCoarseDDProjectionSize = [0 0];
            obj.LastCoarseAngleScoreSize = [0 0];
            obj.LastFineBeamformedSize = [0 0];
            obj.LastFineSpatialProjectionSize = [0 0];
            obj.LastCoarseDDIndex = 0;
            obj.LastCoarseThetaIndex = 0;
            obj.NewPathCoarsePeakScoreHistory = [];
        end

    end

    methods (Access = private)

        function resetDiagnostics(obj, yVec)
            obj.P_detected = 0;
            obj.LastResidualNorm = norm(yVec);
            obj.LastStopMetric = NaN;
            obj.CandidateEvaluations = 0;

            obj.PathSearchCalls = 0;
            obj.CoarseSearchCalls = 0;
            obj.FineDDSearchCalls = 0;
            obj.FineAngleSearchCalls = 0;
            obj.AlternatingCyclesExecuted = 0;

            obj.LastCoarseDDProjectionSize = [0 0];
            obj.LastCoarseAngleScoreSize = [0 0];
            obj.LastFineBeamformedSize = [0 0];
            obj.LastFineSpatialProjectionSize = [0 0];
            obj.LastCoarseDDIndex = 0;
            obj.LastCoarseThetaIndex = 0;
            obj.NewPathCoarsePeakScoreHistory = [];
        end

        function [tauFinal, nuFinal, thetaFinal, hFinal, aFinal, coarsePeakScore] = ...
                estimateOnePath(obj, Yres)

            obj.PathSearchCalls = obj.PathSearchCalls + 1;

            [tauCurrent, nuCurrent, thetaCurrent, aCurrent, coarsePeakScore] = ...
                obj.coarseConditionalSearch(Yres);

            hCurrent = [];

            % Every alternating cycle restarts the SAME shrinking schedule.
            %
            % This is intentional for the L=1/2/3 ablation:
            % repeating the cycle must not implicitly increase the allowed
            % delay/Doppler resolution.
            for ell = 1:obj.NumAlternatingRefinements
                denA = max(real(aCurrent' * aCurrent), realmin);

                yBeamformed = ...
                    Yres * conj(aCurrent) / denA;

                obj.LastFineBeamformedSize = ...
                    size(yBeamformed);

                [tauCurrent, nuCurrent, hCurrent] = ...
                    obj.fineDDSingle( ...
                                     yBeamformed, tauCurrent, nuCurrent);

                [thetaCurrent, aCurrent] = ...
                    obj.fineAngleFromDDProjection( ...
                                                  Yres, hCurrent);

                obj.AlternatingCyclesExecuted = ...
                    obj.AlternatingCyclesExecuted + 1;
            end

            tauFinal = tauCurrent;
            nuFinal = nuCurrent;
            thetaFinal = thetaCurrent;
            hFinal = hCurrent;
            aFinal = aCurrent;
        end

        function [bestTau, bestNu, bestTheta, bestA, bestScore] = ...
                coarseConditionalSearch(obj, Yres)

            obj.CoarseSearchCalls = ...
                obj.CoarseSearchCalls + 1;

            % First reduce MN x Na to Gdd x Na:
            %
            % q_i^T = h_i^H Y_res.
            Qdd = ...
                obj.PsiCoarse' * Yres;

            obj.LastCoarseDDProjectionSize = ...
                size(Qdd);

            % Only after DD projection perform the AoA-grid search.
            spatialResponse = ...
                Qdd * conj(obj.AThetaCoarse);

            denom = ...
                obj.denDDCoarse * obj.denThetaCoarse;

            scores = ...
                abs(spatialResponse).^2 ./ ...
                max(denom, realmin);

            obj.LastCoarseAngleScoreSize = ...
                size(scores);

            % Conditional AoA ML for EACH DD candidate.
            [conditionalMax, thetaIdxPerDD] = ...
                max(scores, [], 2);

            % Then compare DD-candidate maxima.
            [bestScore, ddIdx] = ...
                max(conditionalMax);

            thetaIdx = ...
                thetaIdxPerDD(ddIdx);

            obj.CandidateEvaluations = ...
                obj.CandidateEvaluations + numel(scores);

            obj.LastCoarseDDIndex = ...
                ddIdx;

            obj.LastCoarseThetaIndex = ...
                thetaIdx;

            [tauIdx, nuIdx] = ...
                ind2sub( ...
                        [obj.Gtau, obj.Gnu], ...
                        ddIdx);

            bestTau = ...
                obj.tau_grid_coarse(tauIdx);

            bestNu = ...
                obj.nu_grid_coarse(nuIdx);

            bestTheta = ...
                obj.theta_grid_coarse(thetaIdx);

            bestA = ...
                obj.AThetaCoarse(:, thetaIdx);
        end

        function [tauFine, nuFine, hFine] = ...
                fineDDSingle(obj, yBeamformed, tauInit, nuInit)

            obj.FineDDSearchCalls = ...
                obj.FineDDSearchCalls + 1;

            tauFine = tauInit;
            nuFine = nuInit;

            if obj.max_iter_fine == 0
                hFine = ...
                    obj.getDDResponse(tauFine, nuFine);

                return;
            end

            hFine = [];

            for s = 1:obj.max_iter_fine
                % SAME resolution schedule in every alternating cycle.
                wTau = ...
                    obj.d_tau / (obj.m_tau^(s - 1));

                wNu = ...
                    obj.d_nu / (obj.n_nu^(s - 1));

                Gamma = ...
                    -floor(obj.m_tau / 2):floor(obj.m_tau / 2);

                Lambda = ...
                    -floor(obj.n_nu / 2):floor(obj.n_nu / 2);

                tauGrid = ...
                    tauFine + Gamma * wTau;

                tauGrid = ...
                    tauGrid(tauGrid >= 0);

                if isempty(tauGrid)
                    tauGrid = 0;
                end

                nuGrid = ...
                    nuFine + Lambda * wNu;

                [Psi, denPsi] = ...
                    obj.generateDDGrid( ...
                                       tauGrid, nuGrid);

                score = ...
                    abs(Psi' * yBeamformed).^2 ./ ...
                    max(denPsi, realmin);

                obj.CandidateEvaluations = ...
                    obj.CandidateEvaluations + ...
                    numel(score);

                [~, bestIdx] = ...
                    max(score);

                [tauIdx, nuIdx] = ...
                    ind2sub( ...
                            [numel(tauGrid), numel(nuGrid)], ...
                            bestIdx);

                tauNew = ...
                    tauGrid(tauIdx);

                nuNew = ...
                    nuGrid(nuIdx);

                hNew = ...
                    Psi(:, bestIdx);

                tauFine = tauNew;
                nuFine = nuNew;
                hFine = hNew;
                % An unchanged coarse-grid maximum is not convergence:
                % the next shrinking level can reveal a better candidate.
            end
        end

        function [thetaFine, aFine] = ...
                fineAngleFromDDProjection(obj, Yres, hFine)

            obj.FineAngleSearchCalls = ...
                obj.FineAngleSearchCalls + 1;

            denH = ...
                max(real(hFine' * hFine), realmin);

            % MN-dimensional DD result -> Na-dimensional spatial statistic.
            qSpatial = ...
                (hFine' * Yres) / denH;

            obj.LastFineSpatialProjectionSize = ...
                size(qSpatial);

            % Fine AoA search therefore operates only in Na dimensions.
            e = ...
                qSpatial * conj(obj.AThetaFine);

            score = ...
                abs(e).^2 ./ ...
                max(obj.denThetaFine, realmin);

            obj.CandidateEvaluations = ...
                obj.CandidateEvaluations + ...
                numel(score);

            [~, idx] = ...
                max(score);

            thetaFine = ...
                obj.theta_grid_fine(idx);

            aFine = ...
                obj.AThetaFine(:, idx);
        end

        function g = estimateGains(obj, Y, H, A)
            p = size(H, 2);

            if p == 0
                g = complex(zeros(0, 1));
                return;
            end

            if obj.rho == 0
                B = obj.buildDictionary(H, A);
                y = Y(:);
                g = lsqminnorm(B, y);
            else
                Gram = ...
                    (H' * H) .* (A' * A) + ...
                    obj.rho * eye(p);

                HY = ...
                    H' * Y;

                rhs = ...
                    sum(HY .* A', 2);

                Gram = ...
                    (Gram + Gram') / 2;

                if rcond(Gram) > 1e-12
                    g = Gram \ rhs;
                else
                    g = pinv(Gram) * rhs;
                end
            end
        end

        function Yhat = reconstructMatrix(~, H, A, g)
            if isempty(g)
                Yhat = ...
                    complex(zeros(size(H, 1), size(A, 1)));

                return;
            end

            Yhat = ...
                H * (g(:) .* A.');
        end

        function B = buildDictionary(~, H, A)
            p = size(H, 2);

            B = ...
                complex( ...
                        zeros( ...
                          size(H, 1) * size(A, 1), ...
                          p));

            for k = 1:p
                B(:, k) = ...
                    kron(A(:, k), H(:, k));
            end
        end

        function [Psi, denPsi] = ...
                generateDDGrid(obj, tauGrid, nuGrid)

            GtauLocal = ...
                numel(tauGrid);

            GnuLocal = ...
                numel(nuGrid);

            G = ...
                GtauLocal * GnuLocal;

            Psi = ...
                complex(zeros(obj.MN, G));

            denPsi = ...
                zeros(G, 1);

            idx = 1;

            for j = 1:GnuLocal
                for i = 1:GtauLocal
                    [isCoarse, coarseIdx] = ...
                        obj.lookupCoarse( ...
                                         tauGrid(i), ...
                                         nuGrid(j));

                    if isCoarse
                        h = ...
                            obj.PsiCoarse(:, coarseIdx);

                        den = ...
                            obj.denDDCoarse(coarseIdx);
                    else
                        h = ...
                            obj.generateDDResponse( ...
                                                   tauGrid(i), ...
                                                   nuGrid(j));

                        den = ...
                            real(h' * h);
                    end

                    Psi(:, idx) = h;
                    denPsi(idx) = den;

                    idx = idx + 1;
                end
            end
        end

        function h = getDDResponse(obj, tau, nu)
            [isCoarse, coarseIdx] = ...
                obj.lookupCoarse(tau, nu);

            if isCoarse
                h = ...
                    obj.PsiCoarse(:, coarseIdx);
            else
                h = ...
                    obj.generateDDResponse(tau, nu);
            end
        end

        function [tf, idx] = ...
                lookupCoarse(obj, tau, nu)

            tolTau = ...
                max(1e-14, 1e-10 * obj.d_tau);

            tolNu = ...
                max(1e-10, 1e-10 * obj.d_nu);

            [dt, i] = ...
                min(abs(obj.tau_grid_coarse - tau));

            [dn, j] = ...
                min(abs(obj.nu_grid_coarse - nu));

            tf = ...
                dt <= tolTau && ...
                dn <= tolNu;

            if tf
                idx = ...
                    sub2ind( ...
                            [obj.Gtau, obj.Gnu], ...
                            i, j);
            else
                idx = 0;
            end
        end

        function h = generateDDResponse(obj, tau, nu)
            useCache = isa(obj.DDResponseCache, 'containers.Map');
            if useCache
                key = sprintf('%.17g|%.17g', tau, nu);
                if isKey(obj.DDResponseCache, key)
                    h = obj.DDResponseCache(key);
                    obj.DDResponseCacheHits = obj.DDResponseCacheHits + 1;
                    return;
                end
            end
            phase = ...
                exp(1j * 2 * pi * nu * obj.timeVec);

            xDop = ...
                obj.txSig .* phase;

            xDel = ...
                delayseq( ...
                         xDop, tau, obj.fs);

            [yDD, ~] = ...
                helperOTFSdemod( ...
                                xDel, ...
                                obj.M, ...
                                obj.padlen, ...
                                obj.demodOffset, ...
                                upper(char(obj.padtype)));

            h = ...
                yDD(:);
            if useCache
                obj.DDResponseCache(key) = h;
                obj.DDResponseCacheMisses = obj.DDResponseCacheMisses + 1;
            end
        end

        function offset = getDemodOffset(obj)
            switch upper(char(obj.padtype))
                case {'CP', 'RCP'}
                    offset = obj.padlen;

                case {'ZP', 'RZP', 'NONE'}
                    offset = 0;

                otherwise
                    error( ...
                          'DRAPMDDIPICEstimator:PadType', ...
                          'Unsupported padtype %s.', ...
                          obj.padtype);
            end
        end

        function validateConfiguration(obj)
            if isempty(obj.x_dd) || ...
                    size(obj.x_dd, 1) ~= obj.M || ...
                    size(obj.x_dd, 2) ~= obj.N

                error( ...
                      'DRAPMDDIPICEstimator:PilotSize', ...
                      'x_dd must be an M-by-N matrix.');
            end

            if isempty(obj.tau_grid_coarse) || ...
                    isempty(obj.nu_grid_coarse)

                error( ...
                      'DRAPMDDIPICEstimator:DDGrid', ...
                      ['tau_grid_coarse and ' ...
                       'nu_grid_coarse must be nonempty.']);
            end

            if isempty(obj.theta_grid_coarse) || ...
                    isempty(obj.theta_grid_fine)

                error( ...
                      'DRAPMDDIPICEstimator:AngleGrid', ...
                      ['Both coarse and fine angle ' ...
                       'grids must be nonempty.']);
            end

            if obj.m_tau < 1 || ...
                    obj.n_nu < 1 || ...
                    obj.max_iter_fine < 0 || ...
                    obj.P_max < 1

                error( ...
                      'DRAPMDDIPICEstimator:Configuration', ...
                      ['Invalid refinement or ' ...
                       'path-count configuration.']);
            end

            if obj.rho < 0
                error( ...
                      'DRAPMDDIPICEstimator:Regularization', ...
                      'rho must be nonnegative.');
            end

            if obj.NumAlternatingRefinements < 1 || ...
                    obj.NumAlternatingRefinements ~= ...
                    floor(obj.NumAlternatingRefinements)

                error( ...
                      'DRAPMDDIPICEstimator:AlternatingRefinements', ...
                      ['NumAlternatingRefinements ' ...
                       'must be a positive integer.']);
            end

            validStop = ismember(lower(string(obj.AutoStopRule)), ...
                                 ["residual", "peakdrop"]);
            if ~isscalar(validStop) || ~validStop
                error( ...
                      'DRAPMDDIPICEstimator:AutoStopRule', ...
                      'AutoStopRule must be ''residual'' or ''peakdrop''.');
            end
            if ~isscalar(obj.PeakDropThreshold_dB) || ...
                    ~isfinite(obj.PeakDropThreshold_dB) || ...
                    obj.PeakDropThreshold_dB <= 0
                error( ...
                      'DRAPMDDIPICEstimator:PeakDropThreshold', ...
                      'PeakDropThreshold_dB must be a positive finite scalar.');
            end
        end

    end
end
