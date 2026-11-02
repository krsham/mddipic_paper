classdef MDDAFRBenchmarkEstimator < matlab.System
    % MDDAFRBENCHMARKESTIMATOR MUSIC-assisted DD search with sequential cancellation.
    %
    % [tau,nu,theta,g,P] = step(est,yDD,PExact,sigma2DD)
    %   yDD: stacked M*N*Na DD samples. Coordinates: seconds, Hz and radians.
    %   Known-order operation requires 0<PExact<Na so a noise subspace remains.
    %
    % Implements the structure of Zacharia and Vani Devi, Algorithm 2. MUSIC
    % first selects receive-angle candidates from Y.'*conj(Y), where Y is MN-by-Na.
    % The DD search is crossed with those fixed angles and followed by a single
    % local fractional DD grid. A scalar LS gain and sequential cancellation are
    % used by the paper configuration. Earlier coordinates are not revisited.
    %
    % The physical atom uses this package's finite-frame delayseq convention;
    % the single-transmit response is absorbed into gain. Normalized correlation
    % keeps atom scaling consistent with the gain estimate. This is an adaptation
    % of the source circular full-MIMO model to the shared sampled observation.
    % Project implementation of Zacharia and Vani Devi;
    % see REFERENCES.bib: zachariaDDAFR2026.
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

        % Used only if music_theta_grid is empty.
        theta_grid_fine = []

        % High-resolution physical-angle grid for Eq. (25).
        music_theta_grid = []

        P_max = 10
        eta = 0

        FineDelayPoints = 0
        FineDopplerPoints = 0

        % Canonical reconstruction assumes unit-normalized source columns.
        NormalizeProjection = true
        JointGainUpdate = false

        % 0 -> use Na-PExact smallest-eigenvalue eigenvectors.
        NumNoiseEigenvectors = 0
    end

    properties (SetAccess = private)
        P_detected = 0
        CandidateEvaluations = 0
        MusicSpectrumEvaluations = 0
        CoarseSearchCalls = 0
        FineSearchCalls = 0
        LastResidualNorm = NaN

        LastMusicAngles = []
        LastMusicSpectrum = []
        LastNoiseEigenvalues = []
        LastCoarseIndex = [0 0]
        LastFineIndex = [0 0]
        LastFineGridSize = [0 0]
        LastSelectedMusicIndex = 0
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

        MusicA
        MusicDenA
        MusicGrid

        Gtau
        Gnu
        Gdd
        Gmusic
    end

    methods

        function obj = MDDAFRBenchmarkEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

        function [Pmusic, thetaGrid, thetaPeaks, eigenvalues] = ...
                musicSpectrum(obj, yVec, PExact)
            % Public Eq.-(25) diagnostic.
            yVec = yVec(:);
            obj.validateObservation(yVec);
            obj.validateKnownP(PExact);

            [Pmusic, thetaPeaks, eigenvalues] = ...
                obj.computeMusic(yVec, PExact);
            thetaGrid = obj.MusicGrid;
        end

        function [tauFine, nuFine] = fineGridAround(obj, tau0, nu0)
            [tauFine, nuFine] = obj.makeFineDDGrids(tau0, nu0);
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

            obj.PsiCoarse = complex(zeros(obj.MN, obj.Gdd));
            obj.denDDCoarse = zeros(obj.Gdd, 1);

            idx = 1;
            for j = 1:obj.Gnu
                for i = 1:obj.Gtau
                    h = obj.generateDDResponse( ...
                                               obj.tau_grid_coarse(i), ...
                                               obj.nu_grid_coarse(j));
                    obj.PsiCoarse(:, idx) = h;
                    obj.denDDCoarse(idx) = max(real(h' * h), realmin);
                    idx = idx + 1;
                end
            end

            if isempty(obj.music_theta_grid)
                obj.MusicGrid = obj.theta_grid_fine(:).';
            else
                obj.MusicGrid = obj.music_theta_grid(:).';
            end

            obj.Gmusic = numel(obj.MusicGrid);
            obj.MusicA = complex(zeros(obj.Na, obj.Gmusic));
            obj.MusicDenA = zeros(1, obj.Gmusic);

            for k = 1:obj.Gmusic
                a = steeringVector(obj.MusicGrid(k), obj.Na);
                obj.MusicA(:, k) = a;
                obj.MusicDenA(k) = max(real(a' * a), realmin);
            end
        end

        function [tauEst, nuEst, thetaEst, gEst, Pest] = ...
                stepImpl(obj, yVec, PExact, noiseVarianceDD) %#ok<INUSD>

            yVec = yVec(:);
            obj.validateObservation(yVec);
            obj.validateKnownP(PExact);
            obj.resetDiagnostics();

            [Pmusic, thetaMusic, eigenvalues] = obj.computeMusic(yVec, PExact);
            obj.LastMusicSpectrum = Pmusic;
            obj.LastMusicAngles = thetaMusic(:);
            obj.LastNoiseEigenvalues = eigenvalues(:);
            obj.MusicSpectrumEvaluations = obj.Gmusic;

            Ared = complex(zeros(obj.Na, PExact));
            denAred = zeros(1, PExact);
            for k = 1:PExact
                a = steeringVector(thetaMusic(k), obj.Na);
                Ared(:, k) = a;
                denAred(k) = max(real(a' * a), realmin);
            end

            tauEst = zeros(PExact, 1);
            nuEst = zeros(PExact, 1);
            thetaEst = zeros(PExact, 1);
            gEst = complex(zeros(PExact, 1));
            B = complex(zeros(obj.MN * obj.Na, PExact));

            residual = yVec;
            Pest = 0;

            for p = 1:PExact
                Yres = reshape(residual, obj.MN, obj.Na);

                Q = obj.PsiCoarse' * Yres;
                E = Q * conj(Ared);

                if obj.NormalizeProjection
                    score = abs(E) ./ sqrt(max(obj.denDDCoarse * denAred, realmin));
                else
                    score = abs(E);
                end

                obj.CandidateEvaluations = obj.CandidateEvaluations + numel(score);
                obj.CoarseSearchCalls = obj.CoarseSearchCalls + 1;

                [~, idxLinear] = max(score(:));
                [idd, iMusic] = ind2sub(size(score), idxLinear);
                [itau, inu] = ind2sub([obj.Gtau, obj.Gnu], idd);

                tau0 = obj.tau_grid_coarse(itau);
                nu0 = obj.nu_grid_coarse(inu);
                theta0 = thetaMusic(iMusic);

                obj.LastCoarseIndex = [itau inu];
                obj.LastSelectedMusicIndex = iMusic;

                [tauP, nuP, hP, fineIndex, fineSize] = ...
                    obj.fineDDSearch(Yres, tau0, nu0, theta0);

                obj.FineSearchCalls = obj.FineSearchCalls + 1;
                obj.LastFineIndex = fineIndex;
                obj.LastFineGridSize = fineSize;

                aP = Ared(:, iMusic);
                atom = kron(aP, hP);

                tauEst(p) = tauP;
                nuEst(p) = nuP;
                thetaEst(p) = theta0;
                B(:, p) = atom;

                if obj.JointGainUpdate
                    gEst(1:p) = lsqminnorm(B(:, 1:p), yVec);
                    residual = yVec - B(:, 1:p) * gEst(1:p);
                else
                    den = max(real(atom' * atom), realmin);
                    gEst(p) = (atom' * residual) / den;
                    residual = residual - atom * gEst(p);
                end

                Pest = p;
                obj.P_detected = p;
                obj.LastResidualNorm = norm(residual, 2);

                % In known-P mode eta is diagnostic only; do not stop early.
            end
        end

        function resetImpl(obj)
            obj.resetDiagnostics();
        end

    end

    methods (Access = private)

        function [spectrum, thetaPeaks, eigenvalues] = ...
                computeMusic(obj, yVec, PExact)

            Ydd = reshape(yVec, obj.MN, obj.Na);

            % Project convention: Ydd = sum h_p g_p a_p^T.
            % The cited formulation arranges observations as Nr x MN.
            Ysource = Ydd.';
            R = Ysource * Ysource';
            R = (R + R') / 2;

            [V, D] = eig(R, 'vector');
            [eigenvalues, ord] = sort(real(D), 'ascend');
            V = V(:, ord);

            if obj.NumNoiseEigenvectors == 0
                nNoise = obj.Na - PExact;
            else
                nNoise = obj.NumNoiseEigenvectors;
            end

            if nNoise < 1 || nNoise >= obj.Na
                error('MDDAFRBenchmarkEstimator:NoiseSubspace', ...
                      ['Noise-subspace dimension must satisfy 1 <= e < Na. ', ...
                       'Current e=%d, Na=%d.'], nNoise, obj.Na);
            end

            Vn = V(:, 1:nNoise);
            denom = sum(abs(Vn' * obj.MusicA).^2, 1);
            spectrum = 1 ./ max(real(denom), realmin);

            peakIdx = obj.selectMusicPeaks(spectrum, PExact);
            thetaPeaks = obj.MusicGrid(peakIdx).';
        end

        function peakIdx = selectMusicPeaks(obj, spectrum, P)
            spectrum = spectrum(:).';
            n = numel(spectrum);
            isPeak = false(1, n);

            if n == 1
                isPeak(1) = true;
            else
                isPeak(1) = spectrum(1) >= spectrum(2);
                isPeak(n) = spectrum(n) >= spectrum(n - 1);
                for k = 2:n - 1
                    isPeak(k) = spectrum(k) >= spectrum(k - 1) && ...
                        spectrum(k) >= spectrum(k + 1) && ...
                        (spectrum(k) > spectrum(k - 1) || ...
                         spectrum(k) > spectrum(k + 1));
                end
            end

            localIdx = find(isPeak);
            [~, ord] = sort(spectrum(localIdx), 'descend');
            localIdx = localIdx(ord);

            peakIdx = zeros(0, 1);
            for k = 1:numel(localIdx)
                peakIdx(end + 1, 1) = localIdx(k); %#ok<AGROW>
                if numel(peakIdx) == P
                    return;
                end
            end

            % Source does not specify a fallback. This deterministic fallback
            % is used only when the sampled pseudospectrum contains fewer
            % than P local maxima.
            [~, globalOrd] = sort(spectrum, 'descend');
            for k = 1:numel(globalOrd)
                idx = globalOrd(k);
                if ~any(peakIdx == idx)
                    peakIdx(end + 1, 1) = idx; %#ok<AGROW>
                    if numel(peakIdx) == P
                        return;
                    end
                end
            end

            if numel(peakIdx) < P
                error('MDDAFRBenchmarkEstimator:MusicPeaks', ...
                      'Unable to identify %d MUSIC peak locations.', P);
            end
        end

        function [tauBest, nuBest, hBest, idxBest, gridSize] = ...
                fineDDSearch(obj, Yres, tau0, nu0, theta)

            [tauFine, nuFine] = obj.makeFineDDGrids(tau0, nu0);
            nTau = numel(tauFine);
            nNu = numel(nuFine);
            nDD = nTau * nNu;

            H = complex(zeros(obj.MN, nDD));
            denH = zeros(nDD, 1);

            idx = 1;
            for j = 1:nNu
                for i = 1:nTau
                    h = obj.generateDDResponse(tauFine(i), nuFine(j));
                    H(:, idx) = h;
                    denH(idx) = max(real(h' * h), realmin);
                    idx = idx + 1;
                end
            end

            a = steeringVector(theta, obj.Na);
            denA = max(real(a' * a), realmin);
            E = H' * Yres * conj(a);

            if obj.NormalizeProjection
                score = abs(E) ./ sqrt(max(denH * denA, realmin));
            else
                score = abs(E);
            end

            obj.CandidateEvaluations = obj.CandidateEvaluations + numel(score);

            [~, idd] = max(score);
            [itau, inu] = ind2sub([nTau, nNu], idd);

            tauBest = tauFine(itau);
            nuBest = nuFine(inu);
            hBest = H(:, idd);
            idxBest = [itau inu];
            gridSize = [nTau nNu];
        end

        function [tauFine, nuFine] = makeFineDDGrids(obj, tau0, nu0)
            nTau = obj.FineDelayPoints;
            if nTau == 0
                nTau = floor(obj.M / 2) + 1;
            end

            nNu = obj.FineDopplerPoints;
            if nNu == 0
                nNu = floor(obj.N / 2) + 1;
            end

            tauFine = tau0 + linspace(-0.5, 0.5, nTau) * obj.d_tau;
            tauFine = tauFine(tauFine >= 0);
            if isempty(tauFine)
                tauFine = 0;
            end

            nuFine = nu0 + linspace(-0.5, 0.5, nNu) * obj.d_nu;

            tauFine = tauFine(:).';
            nuFine = nuFine(:).';
        end

        function h = generateDDResponse(obj, tau, nu)
            phase = exp(1j * 2 * pi * nu * obj.timeVec);
            xDop = obj.txSig .* phase;
            xDel = delayseq(xDop, tau, obj.fs);

            [yDD, ~] = helperOTFSdemod( ...
                                      xDel, obj.M, obj.padlen, obj.demodOffset, ...
                                      upper(char(obj.padtype)));
            h = yDD(:);
        end

        function offset = getDemodOffset(obj)
            switch upper(char(obj.padtype))
                case {'CP', 'RCP'}
                    offset = obj.padlen;
                case {'ZP', 'RZP', 'NONE'}
                    offset = 0;
                otherwise
                    error('MDDAFRBenchmarkEstimator:PadType', ...
                          'Unsupported padtype %s.', obj.padtype);
            end
        end

        function validateObservation(obj, yVec)
            if isempty(obj.MN)
                error('MDDAFRBenchmarkEstimator:NotSetup', ...
                      'Call setup(est) or step(est,...) before diagnostics.');
            end
            expected = obj.MN * obj.Na;
            if numel(yVec) ~= expected
                error('MDDAFRBenchmarkEstimator:ObservationSize', ...
                      'yVec must contain M*N*Na = %d samples.', expected);
            end
        end

        function validateKnownP(obj, PExact)
            if ~isscalar(PExact) || PExact < 1 || PExact ~= floor(PExact)
                error('MDDAFRBenchmarkEstimator:KnownPRequired', ...
                      ['MDDAFR benchmark requires known PExact>=1 because ', ...
                       'the source does not specify an unknown-order MUSIC ', ...
                       'model-selection rule.']);
            end
            if PExact > obj.P_max
                error('MDDAFRBenchmarkEstimator:PExact', ...
                      'PExact=%d exceeds P_max=%d.', PExact, obj.P_max);
            end
            if PExact >= obj.Na && obj.NumNoiseEigenvectors == 0
                error('MDDAFRBenchmarkEstimator:MusicRank', ...
                      ['Default MUSIC requires PExact < Na so that the ', ...
                       'noise subspace is nonempty.']);
            end
        end

        function validateConfiguration(obj)
            if obj.M < 1 || obj.N < 1 || obj.Na < 2 || obj.fs <= 0
                error('MDDAFRBenchmarkEstimator:Configuration', ...
                      'M,N,fs must be positive and Na must be at least 2.');
            end

            if isempty(obj.x_dd) || ~isequal(size(obj.x_dd), [obj.M obj.N])
                error('MDDAFRBenchmarkEstimator:PilotSize', ...
                      'x_dd must be an M-by-N matrix.');
            end

            if isempty(obj.tau_grid_coarse) || isempty(obj.nu_grid_coarse)
                error('MDDAFRBenchmarkEstimator:DDGrid', ...
                      'Coarse delay and Doppler grids must be nonempty.');
            end

            if isempty(obj.music_theta_grid) && isempty(obj.theta_grid_fine)
                error('MDDAFRBenchmarkEstimator:MusicGrid', ...
                      'music_theta_grid or theta_grid_fine must be nonempty.');
            end

            if any(obj.tau_grid_coarse < 0)
                error('MDDAFRBenchmarkEstimator:DelayGrid', ...
                      'Delay candidates must be nonnegative.');
            end

            if obj.P_max < 1 || obj.P_max ~= floor(obj.P_max)
                error('MDDAFRBenchmarkEstimator:Pmax', ...
                      'P_max must be a positive integer.');
            end

            values = [obj.FineDelayPoints obj.FineDopplerPoints];
            if any(values < 0) || any(values ~= floor(values)) || any(values == 1)
                error('MDDAFRBenchmarkEstimator:FineGrid', ...
                      ['Fine point counts must be 0 (source default) or ', ...
                       'integers >=2.']);
            end

            if obj.NumNoiseEigenvectors < 0 || ...
                    obj.NumNoiseEigenvectors ~= floor(obj.NumNoiseEigenvectors)
                error('MDDAFRBenchmarkEstimator:NoiseEigenvectors', ...
                      'NumNoiseEigenvectors must be a nonnegative integer.');
            end

            if obj.NumNoiseEigenvectors >= obj.Na
                error('MDDAFRBenchmarkEstimator:NoiseEigenvectors', ...
                      'NumNoiseEigenvectors must be smaller than Na.');
            end

            if ~islogical(obj.NormalizeProjection) || ~isscalar(obj.NormalizeProjection) || ...
                    ~islogical(obj.JointGainUpdate) || ~isscalar(obj.JointGainUpdate)
                error('MDDAFRBenchmarkEstimator:LogicalProperty', ...
                      'NormalizeProjection and JointGainUpdate must be scalar logicals.');
            end
        end

        function resetDiagnostics(obj)
            obj.P_detected = 0;
            obj.CandidateEvaluations = 0;
            obj.MusicSpectrumEvaluations = 0;
            obj.CoarseSearchCalls = 0;
            obj.FineSearchCalls = 0;
            obj.LastResidualNorm = NaN;
            obj.LastMusicAngles = [];
            obj.LastMusicSpectrum = [];
            obj.LastNoiseEigenvalues = [];
            obj.LastCoarseIndex = [0 0];
            obj.LastFineIndex = [0 0];
            obj.LastFineGridSize = [0 0];
            obj.LastSelectedMusicIndex = 0;
        end

    end
end
