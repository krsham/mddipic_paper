classdef AggressiveMDDIPICEstimator < matlab.System
    % AGGRESSIVEMDDIPICESTIMATOR Joint support-conditioned OTFS path estimator.
    %
    % Use MDDIPICEstimator for the paper's joint three-dimensional refinement.
    % This implementation class retains its repository name for compatibility.
    %
    % [tau,nu,theta,g,P] = step(est,yDD,PExact,sigma2DD)
    %   yDD      : (M*N*Na)-by-1, antenna blocks stacked in column order.
    %   PExact   : positive path count; zero selects automatic stopping.
    %   sigma2DD : variance per complex DD sample, used by residual stopping.
    %   tau,nu   : estimated delays [s] and Dopplers [Hz].
    %   theta,g  : receive angles [rad] and complex path gains.
    %
    % Set rho=sigma2DD before setup for the paper's identity-covariance LMMSE
    % policy. The third step input does not set rho automatically. Coarse and
    % fine candidates use the same support-conditioned concentrated objective;
    % accepted path gains are refitted jointly after each coordinate update.
    % The IPI sweep includes the newest path and uses the full observation.
    %
    % UseSchurScoreUpdate selects an exact Schur or direct augmented solve.
    % EnableDDResponseCache reuses temporal responses within one observation.
    % Neither option changes the candidate grids or the objective. The physical
    % atom is kron(a(theta),h(tau,nu)); its inner products factor into spatial
    % and temporal products. See LIBRARY_STRUCTURE.md for the score derivation.
    %
    % Dependencies: helperOTFSmod, helperOTFSdemod, delayseq, steeringVector.
    % See REFERENCES.bib: muppaneniChannelRadarParameter2023 and the manuscript.
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
        % Angular refinement factor/cardinality parameter. In joint-3D
        % mode, Gamma_theta = -floor(m_theta/2):floor(m_theta/2), while
        % the angular step shrinks by m_theta^(s-1), mirroring DDIPIC.
        m_theta = 4
        max_iter_fine = 5

        % Joint 3-D shrinking-grid refinement; false selects split refinement.
        Joint3DRefinement = true

        % Cache fractional DD responses within one observation.
        EnableDDResponseCache = false

        % Evaluate the same concentrated score with Schur updates (true)
        % or a direct augmented normal-equation solve (false).
        UseSchurScoreUpdate = true

        % Compatibility fields; a center tie is not DD convergence.
        % All max_iter_fine levels are always evaluated.
        eps_tau = 1e-12
        eps_nu = 1e-12

        P_max = 10
        rho = 0

        % Automatic stopping. Default uses residual change; 'peakdrop' uses the calibrated relative
        % drop between consecutive new-path coarse Schur peaks.
        stop_scale = 0.01
        AutoStopRule = 'residual'
        PeakDropThreshold_dB = 6

        % Numerical tolerance used only to reject nearly dependent
        % Schur-complement candidates.
        schur_tolerance = 1e-12
    end

    properties (SetAccess = private)
        P_detected = 0
        LastResidualNorm = NaN
        LastStopMetric = NaN
        CandidateEvaluations = 0
        Fine3DSearchCalls = 0
        FineDDSearchCalls = 0
        FineAngleSearchCalls = 0
        DDResponseCacheHits = 0
        DDResponseCacheMisses = 0
        SchurScoreBatches = 0
        DirectCandidateSolves = 0

        % Model-order diagnostics. These histories are observational only:
        % they do not alter candidate ranking or known-P estimator outputs.
        ResidualChangeHistory = []
        ResidualNormHistory = []
        NewPathCoarsePeakScoreHistory = []
        NewPathCoarseBackgroundHistory = []
        NewPathCoarsePBRHistory = []
    end

    properties (Access = private)
        MN
        d_tau
        d_nu
        d_theta

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

        % Capture only the coarse search that acquires each new path.
        % IPI/revisit searches are excluded from the stopping statistic.
        CaptureNewPathCoarseStatistic = false
        CurrentNewPathIndex = 0
    end

    methods

        function obj = AggressiveMDDIPICEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.validateConfiguration();

            obj.MN = obj.M * obj.N;
            obj.d_tau = 1 / obj.fs;
            obj.d_nu = obj.fs / (obj.M * obj.N);
            obj.d_theta = obj.inferAngularSpacing();

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
                error('AggressiveMDDIPICEstimator:ObservationSize', ...
                      'yVec must contain M*N*Na = %d samples.', obj.MN * obj.Na);
            end
            if PExact < 0 || PExact ~= floor(PExact)
                error('AggressiveMDDIPICEstimator:PExact', ...
                      'PExact must be a nonnegative integer.');
            end
            if noiseVarianceDD < 0
                error('AggressiveMDDIPICEstimator:NoiseVariance', ...
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

            % Optional implementation-only cache. The lifetime is one
            % observation so no state can leak across trials.
            if obj.EnableDDResponseCache
                obj.DDResponseCache = containers.Map( ...
                                                     'KeyType', 'char', 'ValueType', 'any');
            else
                obj.DDResponseCache = [];
            end
            obj.DDResponseCacheHits = 0;
            obj.DDResponseCacheMisses = 0;
            obj.SchurScoreBatches = 0;
            obj.DirectCandidateSolves = 0;

            obj.ResidualChangeHistory = nan(Plimit, 1);
            obj.ResidualNormHistory = nan(Plimit, 1);
            obj.NewPathCoarsePeakScoreHistory = nan(Plimit, 1);
            obj.NewPathCoarseBackgroundHistory = nan(Plimit, 1);
            obj.NewPathCoarsePBRHistory = nan(Plimit, 1);
            obj.CaptureNewPathCoarseStatistic = false;
            obj.CurrentNewPathIndex = 0;

            tauEst = zeros(Plimit, 1);
            nuEst = zeros(Plimit, 1);
            thetaEst = zeros(Plimit, 1);
            gEst = complex(zeros(Plimit, 1));

            Hactive = complex(zeros(obj.MN, Plimit));
            Aactive = complex(zeros(obj.Na, Plimit));

            obj.P_detected = 0;
            obj.LastResidualNorm = norm(yVec);
            obj.LastStopMetric = NaN;
            obj.CandidateEvaluations = 0;
            obj.Fine3DSearchCalls = 0;
            obj.FineDDSearchCalls = 0;
            obj.FineAngleSearchCalls = 0;

            residualPrev = yVec;

            for p = 1:Plimit
                Hprev = Hactive(:, 1:p - 1);
                Aprev = Aactive(:, 1:p - 1);

                obj.CaptureNewPathCoarseStatistic = true;
                obj.CurrentNewPathIndex = p;
                [tauP, nuP, thetaP, hP, aP] = ...
                    obj.estimateOnePath(Y, Hprev, Aprev);
                obj.CaptureNewPathCoarseStatistic = false;

                % Optional unknown-order peak-drop stopping. Candidate p is
                % the first rejected component when its coarse support-aware
                % acquisition peak drops by the calibrated threshold relative
                % to candidate p-1. Known-P execution is unaffected.
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

                residualCurrent = yVec - obj.buildDictionary( ...
                                                             Hactive(:, 1:p), Aactive(:, 1:p)) * gActive;

                obj.P_detected = p;
                obj.LastResidualNorm = norm(residualCurrent);

                if p > 1
                    obj.LastStopMetric = norm( ...
                                              residualCurrent - residualPrev)^2;
                    obj.ResidualChangeHistory(p) = obj.LastStopMetric;
                else
                    obj.LastStopMetric = NaN;
                    obj.ResidualChangeHistory(p) = NaN;
                end

                if autoStop && strcmpi(obj.AutoStopRule, 'residual') && p > 1
                    epsilon = ...
                        obj.stop_scale * obj.Na * obj.MN * noiseVarianceDD;
                    if obj.LastStopMetric < epsilon
                        break;
                    end
                end

                % Sequential IPI / coordinate refinement. The newest path
                % is included. Search is performed on the ORIGINAL full
                % observation Y while the other active atoms remain inside
                % the concentrated multipath score. No leave-one-out
                % residual is formed for parameter search.
                %
                % residualCurrent above is retained only for diagnostics
                % and automatic model-order stopping when PExact==0.
                % Gains are refreshed after every path.
                if p > 1
                    for z = 1:p
                        keep = (1:p) ~= z;
                        Hothers = Hactive(:, keep);
                        Aothers = Aactive(:, keep);

                        obj.CaptureNewPathCoarseStatistic = false;
                        [tauZ, nuZ, thetaZ, hZ, aZ] = ...
                            obj.estimateOnePath(Y, Hothers, Aothers);

                        tauEst(z) = tauZ;
                        nuEst(z) = nuZ;
                        thetaEst(z) = thetaZ;
                        Hactive(:, z) = hZ;
                        Aactive(:, z) = aZ;

                        gActive = obj.estimateGains( ...
                                                    Y, Hactive(:, 1:p), Aactive(:, 1:p));
                        gEst(1:p) = gActive;
                    end

                    residualCurrent = yVec - obj.buildDictionary( ...
                                                                 Hactive(:, 1:p), Aactive(:, 1:p)) * gEst(1:p);
                    obj.LastResidualNorm = norm(residualCurrent);
                end

                obj.ResidualNormHistory(p) = norm(residualCurrent);
                residualPrev = residualCurrent;

                % Manual mode stops only after the final IPI sweep.
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
            obj.Fine3DSearchCalls = 0;
            obj.FineDDSearchCalls = 0;
            obj.FineAngleSearchCalls = 0;
            obj.DDResponseCacheHits = 0;
            obj.DDResponseCacheMisses = 0;
            obj.SchurScoreBatches = 0;
            obj.DirectCandidateSolves = 0;
            obj.ResidualChangeHistory = [];
            obj.ResidualNormHistory = [];
            obj.NewPathCoarsePeakScoreHistory = [];
            obj.NewPathCoarseBackgroundHistory = [];
            obj.NewPathCoarsePBRHistory = [];
            obj.CaptureNewPathCoarseStatistic = false;
            obj.CurrentNewPathIndex = 0;
        end

    end

    methods (Access = private)

        function [tauFinal, nuFinal, thetaFinal, hFinal, aFinal] = ...
                estimateOnePath(obj, Y, Hprev, Aprev)

            ctx = obj.makeContext(Y, Hprev, Aprev);

            [tauCoarse, nuCoarse, thetaCoarse] = ...
                obj.coarseSearch(Y, Hprev, Aprev, ctx);

            tauFine = tauCoarse;
            nuFine = nuCoarse;
            thetaFine = thetaCoarse;
            hFine = obj.generateDDResponse(tauFine, nuFine);
            aFine = steeringVector(thetaFine, obj.Na);

            if obj.Joint3DRefinement
                for s = 1:obj.max_iter_fine
                    wTau = obj.d_tau / (obj.m_tau^(s - 1));
                    wNu = obj.d_nu / (obj.n_nu^(s - 1));
                    wTheta = obj.d_theta / (obj.m_theta^(s - 1));

                    Gamma = -floor(obj.m_tau / 2):floor(obj.m_tau / 2);
                    Lambda = -floor(obj.n_nu / 2):floor(obj.n_nu / 2);
                    Omega = -floor(obj.m_theta / 2):floor(obj.m_theta / 2);

                    tauGrid = tauFine + Gamma * wTau;
                    tauGrid = tauGrid(tauGrid >= 0);
                    if isempty(tauGrid)
                        tauGrid = 0;
                    end
                    nuGrid = nuFine + Lambda * wNu;

                    thetaGrid = thetaFine + Omega * wTheta;
                    thetaMin = min(obj.theta_grid_fine(:));
                    thetaMax = max(obj.theta_grid_fine(:));
                    thetaGrid = thetaGrid( ...
                                          thetaGrid >= thetaMin & thetaGrid <= thetaMax);
                    thetaGrid = unique(thetaGrid, 'stable');
                    if isempty(thetaGrid)
                        thetaGrid = min(max(thetaFine, thetaMin), thetaMax);
                    end

                    [tauFine, nuFine, thetaFine, hFine, aFine] = ...
                        obj.fine3DSearch( ...
                                         Y, Hprev, Aprev, ctx, tauGrid, nuGrid, thetaGrid);
                end

                tauFinal = tauFine;
                nuFinal = nuFine;
                thetaFinal = thetaFine;
                hFinal = hFine;
                aFinal = aFine;
                return;
            end

            % Split refinement when Joint3DRefinement=false.
            for s = 1:obj.max_iter_fine
                wTau = obj.d_tau / (obj.m_tau^(s - 1));
                wNu = obj.d_nu / (obj.n_nu^(s - 1));

                Gamma = -floor(obj.m_tau / 2):floor(obj.m_tau / 2);
                Lambda = -floor(obj.n_nu / 2):floor(obj.n_nu / 2);

                tauGrid = tauFine + Gamma * wTau;
                tauGrid = tauGrid(tauGrid >= 0);
                if isempty(tauGrid)
                    tauGrid = 0;
                end
                nuGrid = nuFine + Lambda * wNu;

                [tauFine, nuFine, hFine] = obj.fineDDSearch( ...
                                                          Y, Hprev, Aprev, ctx, thetaCoarse, tauGrid, nuGrid);
            end

            [thetaFinal, aFinal] = obj.fineAngleSearch( ...
                                                      Y, Hprev, Aprev, ctx, tauFine, nuFine, hFine);

            tauFinal = tauFine;
            nuFinal = nuFine;
            hFinal = hFine;
        end

        function [bestTau, bestNu, bestTheta] = ...
                coarseSearch(obj, Y, Hprev, Aprev, ctx)

            E = obj.PsiCoarse' * Y * conj(obj.AThetaCoarse);
            D = obj.denDDCoarse * obj.denThetaCoarse + obj.rho;

            if ctx.numActive == 0
                scores = abs(E).^2 ./ max(D, realmin);
            else
                Hcorr = Hprev' * obj.PsiCoarse;
                Acorr = Aprev' * obj.AThetaCoarse;

                Ccross = complex(zeros( ...
                                       ctx.numActive, obj.Gdd * obj.GthetaCoarse));

                for k = 1:obj.GthetaCoarse
                    cols = (k - 1) * obj.Gdd + (1:obj.Gdd);
                    Ccross(:, cols) = Hcorr .* Acorr(:, k);
                end

                scoresVec = obj.multipathScores( ...
                                                ctx, E(:).', D(:).', Ccross);
                scores = reshape( ...
                                 scoresVec, obj.Gdd, obj.GthetaCoarse);
            end

            obj.CandidateEvaluations = ...
                obj.CandidateEvaluations + numel(scores);

            if obj.CaptureNewPathCoarseStatistic
                validScores = real(scores(isfinite(scores) & scores >= 0));
                if isempty(validScores)
                    peakScore = NaN;
                    bgScore = NaN;
                    pbr = NaN;
                else
                    peakScore = max(validScores);
                    medScore = median(validScores);
                    bgScore = medScore / log(2);
                    pbr = peakScore / max(bgScore, realmin);
                end
                idxNew = obj.CurrentNewPathIndex;
                if idxNew >= 1 && idxNew <= numel(obj.NewPathCoarsePBRHistory)
                    obj.NewPathCoarsePeakScoreHistory(idxNew) = peakScore;
                    obj.NewPathCoarseBackgroundHistory(idxNew) = bgScore;
                    obj.NewPathCoarsePBRHistory(idxNew) = pbr;
                end
            end

            [~, linearIdx] = max(scores(:));
            [ddIdx, thetaIdx] = ind2sub(size(scores), linearIdx);
            [tauIdx, nuIdx] = ind2sub( ...
                                     [obj.Gtau, obj.Gnu], ddIdx);

            bestTau = obj.tau_grid_coarse(tauIdx);
            bestNu = obj.nu_grid_coarse(nuIdx);
            bestTheta = obj.theta_grid_coarse(thetaIdx);
        end

        function [bestTau, bestNu, bestTheta, bestH, bestA] = ...
                fine3DSearch(obj, Y, Hprev, Aprev, ctx, tauGrid, nuGrid, thetaGrid)
            [Psi, denPsi] = obj.generateDDGrid(tauGrid, nuGrid);

            Gtheta = numel(thetaGrid);
            A = complex(zeros(obj.Na, Gtheta));
            denTheta = zeros(1, Gtheta);
            for k = 1:Gtheta
                a = steeringVector(thetaGrid(k), obj.Na);
                A(:, k) = a;
                denTheta(k) = real(a' * a);
            end

            E = Psi' * Y * conj(A);
            D = denPsi * denTheta + obj.rho;

            if ctx.numActive == 0
                scores = abs(E).^2 ./ max(D, realmin);
            else
                Hcorr = Hprev' * Psi;
                Acorr = Aprev' * A;

                GddLocal = size(Psi, 2);
                Ccross = complex(zeros( ...
                                       ctx.numActive, GddLocal * Gtheta));
                for k = 1:Gtheta
                    cols = (k - 1) * GddLocal + (1:GddLocal);
                    Ccross(:, cols) = Hcorr .* Acorr(:, k);
                end

                scoresVec = obj.multipathScores( ...
                                                ctx, E(:).', D(:).', Ccross);
                scores = reshape(scoresVec, GddLocal, Gtheta);
            end

            obj.CandidateEvaluations = ...
                obj.CandidateEvaluations + numel(scores);
            obj.Fine3DSearchCalls = obj.Fine3DSearchCalls + 1;

            [~, linearIdx] = max(scores(:));
            [ddIdx, thetaIdx] = ind2sub(size(scores), linearIdx);
            [tauIdx, nuIdx] = ind2sub( ...
                                     [numel(tauGrid), numel(nuGrid)], ddIdx);

            bestTau = tauGrid(tauIdx);
            bestNu = nuGrid(nuIdx);
            bestTheta = thetaGrid(thetaIdx);
            bestH = Psi(:, ddIdx);
            bestA = A(:, thetaIdx);
        end

        function [bestTau, bestNu, bestH] = fineDDSearch( ...
                                                       obj, Y, Hprev, Aprev, ctx, thetaFixed, tauGrid, nuGrid)

            [Psi, denPsi] = obj.generateDDGrid(tauGrid, nuGrid);
            a = steeringVector(thetaFixed, obj.Na);
            denA = real(a' * a);

            yProjected = Y * conj(a);
            e = Psi' * yProjected;
            d = denPsi * denA + obj.rho;

            if ctx.numActive == 0
                scores = abs(e).^2 ./ max(d, realmin);
            else
                Hcorr = Hprev' * Psi;
                Acorr = Aprev' * a;
                Ccross = Hcorr .* Acorr;
                scores = obj.multipathScores( ...
                                             ctx, e(:).', d(:).', Ccross).';
            end

            obj.CandidateEvaluations = ...
                obj.CandidateEvaluations + numel(scores);
            obj.FineDDSearchCalls = obj.FineDDSearchCalls + 1;

            [~, bestIdx] = max(scores);
            [tauIdx, nuIdx] = ind2sub( ...
                                     [numel(tauGrid), numel(nuGrid)], bestIdx);

            bestTau = tauGrid(tauIdx);
            bestNu = nuGrid(nuIdx);
            bestH = Psi(:, bestIdx);
        end

        function [bestTheta, bestA] = fineAngleSearch( ...
                                                     obj, Y, Hprev, Aprev, ctx, tauFixed, nuFixed, hFixed)

            if nargin < 8 || isempty(hFixed)
                hFixed = obj.generateDDResponse(tauFixed, nuFixed);
            end

            rowProjection = hFixed' * Y;
            e = rowProjection * conj(obj.AThetaFine);
            denH = real(hFixed' * hFixed);
            d = denH * obj.denThetaFine + obj.rho;

            if ctx.numActive == 0
                scores = abs(e).^2 ./ max(d, realmin);
            else
                Hcorr = Hprev' * hFixed;
                Acorr = Aprev' * obj.AThetaFine;
                Ccross = Hcorr .* Acorr;
                scores = obj.multipathScores(ctx, e, d, Ccross);
            end

            obj.CandidateEvaluations = ...
                obj.CandidateEvaluations + numel(scores);
            obj.FineAngleSearchCalls = obj.FineAngleSearchCalls + 1;

            [~, idx] = max(scores);
            bestTheta = obj.theta_grid_fine(idx);
            bestA = obj.AThetaFine(:, idx);
        end

        function ctx = makeContext(obj, Y, Hprev, Aprev)
            q = size(Hprev, 2);

            ctx.numActive = q;
            ctx.C = [];
            ctx.q = [];
            ctx.u = [];

            if q == 0
                return;
            end

            GramH = Hprev' * Hprev;
            GramA = Aprev' * Aprev;
            C = GramH .* GramA + obj.rho * eye(q);
            C = (C + C') / 2;

            HY = Hprev' * Y;
            qvec = sum(HY .* Aprev', 2);

            ctx.C = C;
            ctx.q = qvec;
            ctx.u = obj.safeSolve(C, qvec);
        end

        function scores = multipathScores(obj, ctx, e, d, Ccross)
            % MULTIPATHSCORES Exact alternative implementations of the
            % same concentrated multipath score.
            if obj.UseSchurScoreUpdate
                obj.SchurScoreBatches = obj.SchurScoreBatches + 1;
                scores = obj.schurScores(ctx, e, d, Ccross);
            else
                scores = obj.directAugmentedScores(ctx, e, d, Ccross);
            end
        end

        function scores = directAugmentedScores(obj, ctx, e, d, Ccross)
            % DIRECTAUGMENTEDSCORES Reference per-candidate augmented solve.
            %
            % For each candidate, solve
            %   [C c; c^H d] [x] = [q; e]
            % and evaluate the concentrated objective increment relative to
            % the fixed active support. The Kronecker/factored inner
            % products used to construct C,c,d are shared with Schur mode;
            % only the block-system evaluation differs.
            %
            % A candidate is subjected to the same near-dependence
            % admissibility guard as Schur mode so numerical pathology does
            % not masquerade as an algorithm difference.

            K = numel(e);
            scores = -inf(size(e));

            % Numerical-admissibility guard only; candidate scores below are
            % still computed by independent augmented solves.
            V = obj.safeSolve(ctx.C, Ccross);
            schur = real(d - sum(conj(Ccross) .* V, 1));
            scale = max(1, max(abs(d)));
            tol = obj.schur_tolerance * scale;
            valid = schur > tol;

            base = real(ctx.q' * ctx.u);

            qActive = size(ctx.C, 1);
            G = complex(zeros(qActive + 1, qActive + 1));
            G(1:qActive, 1:qActive) = ctx.C;
            qplus = complex(zeros(qActive + 1, 1));
            qplus(1:qActive) = ctx.q;

            for k = 1:K
                if ~valid(k)
                    continue;
                end

                c = Ccross(:, k);
                G(1:qActive, end) = c;
                G(end, 1:qActive) = c';
                G(end, end) = d(k);
                qplus(end) = e(k);

                x = obj.safeSolve((G + G') / 2, qplus);
                scores(k) = real(qplus' * x) - base;
                obj.DirectCandidateSolves = ...
                    obj.DirectCandidateSolves + 1;
            end
        end

        function scores = schurScores(obj, ctx, e, d, Ccross)
            V = obj.safeSolve(ctx.C, Ccross);

            cCinvC = sum(conj(Ccross) .* V, 1);
            cCinvy = sum(conj(Ccross) .* ctx.u, 1);

            schur = real(d - cCinvC);

            scale = max(1, max(abs(d)));
            tol = obj.schur_tolerance * scale;

            numer = abs(e - cCinvy).^2;
            scores = -inf(size(numer));

            valid = schur > tol;
            scores(valid) = numer(valid) ./ schur(valid);
        end

        function X = safeSolve(~, A, B)
            if isempty(A)
                X = zeros(0, size(B, 2), 'like', B);
                return;
            end

            rc = rcond(A);

            if isfinite(rc) && rc > 1e-12
                X = A \ B;
            else
                X = pinv(A) * B;
            end
        end

        function g = estimateGains(obj, Y, H, A)
            p = size(H, 2);

            if p == 0
                g = complex(zeros(0, 1));
                return;
            end

            B = obj.buildDictionary(H, A);
            y = Y(:);

            if obj.rho == 0
                % Numerically robust ordinary LS.
                g = lsqminnorm(B, y);
            else
                Gram = (H' * H) .* (A' * A) + obj.rho * eye(p);
                HY = H' * Y;
                rhs = sum(HY .* A', 2);
                g = obj.safeSolve((Gram + Gram') / 2, rhs);
            end
        end

        function B = buildDictionary(~, H, A)
            p = size(H, 2);
            B = complex(zeros(size(H, 1) * size(A, 1), p));

            for k = 1:p
                B(:, k) = kron(A(:, k), H(:, k));
            end
        end

        function [Psi, denPsi] = generateDDGrid(obj, tauGrid, nuGrid)
            GtauLocal = numel(tauGrid);
            GnuLocal = numel(nuGrid);
            G = GtauLocal * GnuLocal;

            Psi = complex(zeros(obj.MN, G));
            denPsi = zeros(G, 1);

            idx = 1;

            for j = 1:GnuLocal
                for i = 1:GtauLocal

                    [isCoarse, coarseIdx] = obj.lookupCoarse( ...
                                                            tauGrid(i), nuGrid(j));

                    if isCoarse
                        h = obj.PsiCoarse(:, coarseIdx);
                        den = obj.denDDCoarse(coarseIdx);
                    else
                        h = obj.generateDDResponse( ...
                                                   tauGrid(i), nuGrid(j));
                        den = real(h' * h);
                    end

                    Psi(:, idx) = h;
                    denPsi(idx) = den;
                    idx = idx + 1;
                end
            end
        end

        function [tf, idx] = lookupCoarse(obj, tau, nu)
            tolTau = max(1e-14, 1e-10 * obj.d_tau);
            tolNu = max(1e-10, 1e-10 * obj.d_nu);

            [dt, i] = min(abs(obj.tau_grid_coarse - tau));
            [dn, j] = min(abs(obj.nu_grid_coarse - nu));

            tf = dt <= tolTau && dn <= tolNu;

            if tf
                idx = sub2ind([obj.Gtau, obj.Gnu], i, j);
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
                    obj.DDResponseCacheHits = ...
                        obj.DDResponseCacheHits + 1;
                    return;
                end
            end

            phase = exp(1j * 2 * pi * nu * obj.timeVec);
            xDop = obj.txSig .* phase;

            xDel = delayseq( ...
                            xDop, tau, obj.fs);

            [yDD, ~] = helperOTFSdemod( ...
                                      xDel, obj.M, obj.padlen, obj.demodOffset, ...
                                      upper(char(obj.padtype)));

            h = yDD(:);

            if useCache
                obj.DDResponseCache(key) = h;
                obj.DDResponseCacheMisses = ...
                    obj.DDResponseCacheMisses + 1;
            end
        end

        function offset = getDemodOffset(obj)
            switch upper(obj.padtype)
                case {'CP', 'RCP'}
                    offset = obj.padlen;

                case {'ZP', 'RZP', 'NONE'}
                    offset = 0;

                otherwise
                    error('AggressiveMDDIPICEstimator:PadType', ...
                          'Unsupported padtype %s.', obj.padtype);
            end
        end

        function dTheta = inferAngularSpacing(obj)
            coarse = sort(unique(obj.theta_grid_coarse(:)));
            if numel(coarse) > 1
                delta = diff(coarse);
                dTheta = median(delta(delta > 0));
            else
                fine = sort(unique(obj.theta_grid_fine(:)));
                if numel(fine) > 1
                    delta = diff(fine);
                    dTheta = median(delta(delta > 0));
                else
                    error('AggressiveMDDIPICEstimator:AngleGrid', ...
                          'At least two distinct AoA grid points are required.');
                end
            end
            if isempty(dTheta) || ~isfinite(dTheta) || dTheta <= 0
                error('AggressiveMDDIPICEstimator:AngleGrid', ...
                      'Unable to infer a positive angular grid spacing.');
            end
        end

        function validateConfiguration(obj)
            if isempty(obj.x_dd) || ...
                    size(obj.x_dd, 1) ~= obj.M || ...
                    size(obj.x_dd, 2) ~= obj.N

                error('AggressiveMDDIPICEstimator:PilotSize', ...
                      'x_dd must be an M-by-N matrix.');
            end

            if isempty(obj.tau_grid_coarse) || ...
                    isempty(obj.nu_grid_coarse)

                error('AggressiveMDDIPICEstimator:DDGrid', ...
                      ['tau_grid_coarse and nu_grid_coarse ' ...
                       'must be nonempty.']);
            end

            if isempty(obj.theta_grid_coarse) || ...
                    isempty(obj.theta_grid_fine)

                error('AggressiveMDDIPICEstimator:AngleGrid', ...
                      ['Both coarse and fine angle grids ' ...
                       'must be nonempty.']);
            end

            if obj.m_tau < 1 || obj.n_nu < 1 || obj.m_theta < 1 || ...
                    obj.m_theta ~= floor(obj.m_theta) || ...
                    obj.max_iter_fine < 0 || obj.P_max < 1

                error('AggressiveMDDIPICEstimator:Configuration', ...
                      'Invalid refinement or path-count configuration.');
            end
            if ~islogical(obj.Joint3DRefinement) || ...
                    ~isscalar(obj.Joint3DRefinement)
                error('AggressiveMDDIPICEstimator:Joint3DRefinement', ...
                      'Joint3DRefinement must be a logical scalar.');
            end

            if ~islogical(obj.EnableDDResponseCache) || ...
                    ~isscalar(obj.EnableDDResponseCache)
                error('AggressiveMDDIPICEstimator:DDResponseCache', ...
                      'EnableDDResponseCache must be a logical scalar.');
            end

            if ~islogical(obj.UseSchurScoreUpdate) || ...
                    ~isscalar(obj.UseSchurScoreUpdate)
                error('AggressiveMDDIPICEstimator:SchurScoreUpdate', ...
                      'UseSchurScoreUpdate must be a logical scalar.');
            end

            validStop = ismember(lower(string(obj.AutoStopRule)), ["residual", "peakdrop"]);
            if ~isscalar(validStop) || ~validStop
                error('AggressiveMDDIPICEstimator:AutoStopRule', ...
                      'AutoStopRule must be ''residual'' or ''peakdrop''.');
            end
            if ~isscalar(obj.PeakDropThreshold_dB) || ...
                    ~isfinite(obj.PeakDropThreshold_dB) || obj.PeakDropThreshold_dB <= 0
                error('AggressiveMDDIPICEstimator:PeakDropThreshold', ...
                      'PeakDropThreshold_dB must be a positive finite scalar.');
            end

            if obj.rho < 0
                error('AggressiveMDDIPICEstimator:Regularization', ...
                      'rho must be nonnegative.');
            end
        end

    end
end
