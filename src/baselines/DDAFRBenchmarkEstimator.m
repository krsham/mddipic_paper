classdef DDAFRBenchmarkEstimator < matlab.System
    % DDAFRBENCHMARKESTIMATOR Coarse/fine DD-angle search with sequential cancellation.
    %
    % [tau,nu,theta,g,P] = step(est,yDD,PExact,sigma2DD)
    %   yDD: stacked M*N*Na DD samples. Coordinates: seconds, Hz and radians.
    %   PExact>0 supplies the order. Automatic stopping requires an explicit eta.
    %
    % Implements the search/cancellation structure of Zacharia and Vani Devi,
    % Algorithm 1. The single-transmit receive-array atom is formed with this
    % package's finite-frame Doppler-before-delay operator and centered ULA.
    % The transmit response is absorbed into the complex gain. This numerical
    % atom differs from the source's circular full-MIMO dictionary.
    %
    % Default normalized projection uses |r'*e|^2/(r'*r), followed by the scalar
    % LS gain (r'*e)/(r'*r). One local fractional grid is searched per acquisition.
    % Previously accepted coordinates are held fixed. The paper configuration
    % sets JointGainUpdate=false and NormalizeProjection=true explicitly.
    %
    % FineDelayPoints=0 and FineDopplerPoints=0 request M/2+1 and N/2+1 points
    % across a half-bin window on either side of the coarse estimate. The angular
    % fine search is restricted to its coarse cell.
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
        theta_grid_coarse = []
        theta_grid_fine = []

        P_max = 10

        % Source stopping threshold eta. Used only when PExact == 0.
        eta = 0

        % 0 selects the source defaults M/2+1 and N/2+1.
        FineDelayPoints = 0
        FineDopplerPoints = 0

        % Canonical reconstruction assumes unit-normalized source columns,
        % as implied by Algorithm-1 line 9's projection-as-gain update.
        NormalizeProjection = true

        % false = printed Algorithm-1 sequential single-path update.
        % true  = MODIFIED active-dictionary joint LS after every addition.
        JointGainUpdate = false
    end

    properties (SetAccess = private)
        P_detected = 0
        CandidateEvaluations = 0
        CoarseSearchCalls = 0
        FineSearchCalls = 0
        LastResidualNorm = NaN
        LastCoarseIndex = [0 0 0]
        LastFineIndex = [0 0 0]
        LastFineGridSize = [0 0 0]
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

        Gtau
        Gnu
        Gdd
        Gtheta
    end

    methods

        function obj = DDAFRBenchmarkEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

        function [scores, cube] = coarseScoreMap(obj, yVec)
            % Public diagnostic for Algorithm-1 line 3.
            yVec = yVec(:);
            obj.validateObservation(yVec);

            Y = reshape(yVec, obj.MN, obj.Na);
            scores = obj.computeCoarseScores(Y);
            cube = reshape(scores, obj.Gtau, obj.Gnu, obj.Gtheta);
        end

        function [tauFine, nuFine, thetaFine] = fineGridAround( ...
                                                             obj, tau0, nu0, theta0)
            % Public diagnostic exposing the source Fig.-2 local grid.
            thetaIdx = obj.nearestCoarseThetaIndex(theta0);
            [tauFine, nuFine, thetaFine] = obj.makeFineGrids( ...
                                                           tau0, nu0, theta0, thetaIdx);
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
            obj.Gtheta = numel(obj.theta_grid_coarse);

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

            obj.AThetaCoarse = complex(zeros(obj.Na, obj.Gtheta));
            obj.denThetaCoarse = zeros(1, obj.Gtheta);

            for k = 1:obj.Gtheta
                a = steeringVector(obj.theta_grid_coarse(k), obj.Na);
                obj.AThetaCoarse(:, k) = a;
                obj.denThetaCoarse(k) = max(real(a' * a), realmin);
            end
        end

        function [tauEst, nuEst, thetaEst, gEst, Pest] = ...
                stepImpl(obj, yVec, PExact, noiseVarianceDD) %#ok<INUSD>

            yVec = yVec(:);
            obj.validateObservation(yVec);
            obj.validatePathRequest(PExact);
            obj.resetDiagnostics();

            if PExact > 0
                nMax = PExact;
            else
                if obj.eta <= 0
                    error('DDAFRBenchmarkEstimator:EtaRequired', ...
                          ['PExact=0 requires a positive source stopping ', ...
                           'threshold eta. The paper does not prescribe ', ...
                           'a universal eta calibration.']);
                end
                nMax = obj.P_max;
            end

            tauEst = zeros(nMax, 1);
            nuEst = zeros(nMax, 1);
            thetaEst = zeros(nMax, 1);
            gEst = complex(zeros(nMax, 1));
            B = complex(zeros(obj.MN * obj.Na, nMax));

            residual = yVec;
            Pest = 0;

            for p = 1:nMax
                Yres = reshape(residual, obj.MN, obj.Na);

                coarseScores = obj.computeCoarseScores(Yres);
                obj.CandidateEvaluations = obj.CandidateEvaluations + ...
                    numel(coarseScores);
                obj.CoarseSearchCalls = obj.CoarseSearchCalls + 1;

                [~, idxLinear] = max(coarseScores(:));
                [idd, itheta] = ind2sub(size(coarseScores), idxLinear);
                [itau, inu] = ind2sub([obj.Gtau, obj.Gnu], idd);

                tau0 = obj.tau_grid_coarse(itau);
                nu0 = obj.nu_grid_coarse(inu);
                theta0 = obj.theta_grid_coarse(itheta);
                obj.LastCoarseIndex = [itau inu itheta];

                [tauP, nuP, thetaP, hP, aP, fineIndex, fineSize] = ...
                    obj.fineSearch(Yres, tau0, nu0, theta0, itheta);

                obj.FineSearchCalls = obj.FineSearchCalls + 1;
                obj.LastFineIndex = fineIndex;
                obj.LastFineGridSize = fineSize;

                tauEst(p) = tauP;
                nuEst(p) = nuP;
                thetaEst(p) = thetaP;

                atom = kron(aP, hP);
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

                if PExact == 0 && obj.LastResidualNorm <= obj.eta
                    break;
                end
            end

            tauEst = tauEst(1:Pest);
            nuEst = nuEst(1:Pest);
            thetaEst = thetaEst(1:Pest);
            gEst = gEst(1:Pest);
        end

        function resetImpl(obj)
            obj.resetDiagnostics();
        end

    end

    methods (Access = private)

        function scores = computeCoarseScores(obj, Y)
            % Separable exact evaluation of |Psi^H r| on the common model.
            Q = obj.PsiCoarse' * Y;
            E = Q * conj(obj.AThetaCoarse);

            if obj.NormalizeProjection
                den = obj.denDDCoarse * obj.denThetaCoarse;
                scores = abs(E) ./ sqrt(max(den, realmin));
            else
                scores = abs(E);
            end
        end

        function [tauBest, nuBest, thetaBest, hBest, aBest, idxBest, gridSize] = ...
                fineSearch(obj, Yres, tau0, nu0, theta0, thetaCoarseIdx)

            [tauFine, nuFine, thetaFine] = obj.makeFineGrids( ...
                                                           tau0, nu0, theta0, thetaCoarseIdx);

            nTau = numel(tauFine);
            nNu = numel(nuFine);
            nTheta = numel(thetaFine);
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

            A = complex(zeros(obj.Na, nTheta));
            denA = zeros(1, nTheta);
            for k = 1:nTheta
                a = steeringVector(thetaFine(k), obj.Na);
                A(:, k) = a;
                denA(k) = max(real(a' * a), realmin);
            end

            E = (H' * Yres) * conj(A);

            if obj.NormalizeProjection
                score = abs(E) ./ sqrt(max(denH * denA, realmin));
            else
                score = abs(E);
            end

            obj.CandidateEvaluations = obj.CandidateEvaluations + numel(score);

            [~, idxLinear] = max(score(:));
            [idd, itheta] = ind2sub(size(score), idxLinear);
            [itau, inu] = ind2sub([nTau, nNu], idd);

            tauBest = tauFine(itau);
            nuBest = nuFine(inu);
            thetaBest = thetaFine(itheta);
            hBest = H(:, idd);
            aBest = A(:, itheta);
            idxBest = [itau inu itheta];
            gridSize = [nTau nNu nTheta];
        end

        function [tauFine, nuFine, thetaFine] = ...
                makeFineGrids(obj, tau0, nu0, theta0, thetaCoarseIdx)

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

            if obj.Gtheta == 1
                thetaFine = obj.theta_grid_fine(:).';
            else
                if thetaCoarseIdx <= 1
                    cellWidth = abs(obj.theta_grid_coarse(2) - ...
                                    obj.theta_grid_coarse(1));
                elseif thetaCoarseIdx >= obj.Gtheta
                    cellWidth = abs(obj.theta_grid_coarse(end) - ...
                                    obj.theta_grid_coarse(end - 1));
                else
                    left = abs(obj.theta_grid_coarse(thetaCoarseIdx) - ...
                               obj.theta_grid_coarse(thetaCoarseIdx - 1));
                    right = abs(obj.theta_grid_coarse(thetaCoarseIdx + 1) - ...
                                obj.theta_grid_coarse(thetaCoarseIdx));
                    cellWidth = min(left, right);
                end

                halfWidth = 0.5 * cellWidth;
                mask = abs(obj.theta_grid_fine - theta0) <= halfWidth + 1e-12;
                thetaFine = obj.theta_grid_fine(mask);
            end

            if isempty(thetaFine)
                thetaFine = theta0;
            end

            tauFine = tauFine(:).';
            nuFine = nuFine(:).';
            thetaFine = thetaFine(:).';
        end

        function idx = nearestCoarseThetaIndex(obj, theta)
            [~, idx] = min(abs(obj.theta_grid_coarse - theta));
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
                    error('DDAFRBenchmarkEstimator:PadType', ...
                          'Unsupported padtype %s.', obj.padtype);
            end
        end

        function validateObservation(obj, yVec)
            if isempty(obj.MN)
                error('DDAFRBenchmarkEstimator:NotSetup', ...
                      'Call setup(est) or step(est,...) before diagnostics.');
            end
            expected = obj.MN * obj.Na;
            if numel(yVec) ~= expected
                error('DDAFRBenchmarkEstimator:ObservationSize', ...
                      'yVec must contain M*N*Na = %d samples.', expected);
            end
        end

        function validatePathRequest(obj, PExact)
            if ~isscalar(PExact) || PExact < 0 || PExact ~= floor(PExact)
                error('DDAFRBenchmarkEstimator:PExact', ...
                      'PExact must be a nonnegative integer.');
            end
            if PExact > obj.P_max
                error('DDAFRBenchmarkEstimator:PExact', ...
                      'PExact=%d exceeds P_max=%d.', PExact, obj.P_max);
            end
        end

        function validateConfiguration(obj)
            if obj.M < 1 || obj.N < 1 || obj.Na < 1 || obj.fs <= 0
                error('DDAFRBenchmarkEstimator:Configuration', ...
                      'M, N, Na, and fs must be positive.');
            end

            if isempty(obj.x_dd) || ~isequal(size(obj.x_dd), [obj.M obj.N])
                error('DDAFRBenchmarkEstimator:PilotSize', ...
                      'x_dd must be an M-by-N matrix.');
            end

            if isempty(obj.tau_grid_coarse) || isempty(obj.nu_grid_coarse) || ...
                    isempty(obj.theta_grid_coarse) || isempty(obj.theta_grid_fine)
                error('DDAFRBenchmarkEstimator:Grid', ...
                      'All coarse grids and theta_grid_fine must be nonempty.');
            end

            if any(obj.tau_grid_coarse < 0)
                error('DDAFRBenchmarkEstimator:DelayGrid', ...
                      'Delay candidates must be nonnegative.');
            end

            if obj.P_max < 1 || obj.P_max ~= floor(obj.P_max)
                error('DDAFRBenchmarkEstimator:Pmax', ...
                      'P_max must be a positive integer.');
            end

            values = [obj.FineDelayPoints obj.FineDopplerPoints];
            if any(values < 0) || any(values ~= floor(values))
                error('DDAFRBenchmarkEstimator:FineGrid', ...
                      'FineDelayPoints/FineDopplerPoints must be nonnegative integers.');
            end
            if any(values == 1)
                error('DDAFRBenchmarkEstimator:FineGrid', ...
                      'Fine grid point counts must be 0 (source default) or >=2.');
            end

            if ~islogical(obj.NormalizeProjection) || ~isscalar(obj.NormalizeProjection) || ...
                    ~islogical(obj.JointGainUpdate) || ~isscalar(obj.JointGainUpdate)
                error('DDAFRBenchmarkEstimator:LogicalProperty', ...
                      'NormalizeProjection and JointGainUpdate must be scalar logicals.');
            end
        end

        function resetDiagnostics(obj)
            obj.P_detected = 0;
            obj.CandidateEvaluations = 0;
            obj.CoarseSearchCalls = 0;
            obj.FineSearchCalls = 0;
            obj.LastResidualNorm = NaN;
            obj.LastCoarseIndex = [0 0 0];
            obj.LastFineIndex = [0 0 0];
            obj.LastFineGridSize = [0 0 0];
        end

    end
end
