classdef SrivastavaOMPBenchmarkEstimator < matlab.System
    % SRIVASTAVAOMPBENCHMARKESTIMATOR Fixed-grid OMP on array time-domain samples.
    %
    % [tau,nu,theta,g,P] = step(est,yTD,PExact,sigma2TD)
    %   yTD: L-by-Na samples or their column-major vector, L=M*N without padding.
    %   sigma2TD: noise variance in time-domain units.
    %   Outputs: seconds, Hz, radians, complex gains and selected path count.
    %
    % Temporal and receive-array dictionaries represent the same finite-frame
    % physical atom as the proposed estimators. Selected coordinates stay on the
    % grid; active gains are refitted jointly by LS. The single-transmit host
    % absorbs the transmit response into gain. The communication demo uses
    % this class; sensing uses Classic3DOMPSensingEstimator on DD samples.
    % Single-transmit adaptation of the sparse OTFS baseline;
    % see REFERENCES.bib: srivastava4DSparse2022.
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

        tau_grid = []
        nu_grid = []
        theta_grid = []

        P_max = 10
        NormalizeSelection = true
        AtomNormFloor = 1e-14
    end

    properties (SetAccess = private)
        P_detected = 0
        LastResidualNorm = NaN
        LastResidualEnergy = NaN
        LastNoiseEnergyThreshold = NaN
        CandidateEvaluations = 0
        SelectedLinearIndices = []
        SelectedDDIndices = []
        SelectedThetaIndices = []
    end

    properties (Access = private)
        txSig
        timeVec
        L

        PsiTD
        denTD
        ATheta
        denTheta

        tauOfDD
        nuOfDD
        Gtau
        Gnu
        Gtheta
        Gdd
    end

    methods

        function obj = SrivastavaOMPBenchmarkEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.validateConfiguration();

            [obj.txSig, ~] = helperOTFSmod( ...
                                          obj.x_dd, obj.padlen, upper(char(obj.padtype)));
            obj.txSig = obj.txSig(:);
            obj.L = numel(obj.txSig);
            obj.timeVec = (0:obj.L - 1).' / obj.fs;

            obj.Gtau = numel(obj.tau_grid);
            obj.Gnu = numel(obj.nu_grid);
            obj.Gtheta = numel(obj.theta_grid);
            obj.Gdd = obj.Gtau * obj.Gnu;

            % Temporal dictionary. Loop order is Doppler outer, delay inner.
            obj.PsiTD = complex(zeros(obj.L, obj.Gdd));
            obj.denTD = zeros(obj.Gdd, 1);
            obj.tauOfDD = zeros(obj.Gdd, 1);
            obj.nuOfDD = zeros(obj.Gdd, 1);

            idx = 1;
            for j = 1:obj.Gnu
                for i = 1:obj.Gtau
                    h = obj.generateTDResponse( ...
                                               obj.tau_grid(i), obj.nu_grid(j));
                    obj.PsiTD(:, idx) = h;
                    obj.denTD(idx) = real(h' * h);
                    obj.tauOfDD(idx) = obj.tau_grid(i);
                    obj.nuOfDD(idx) = obj.nu_grid(j);
                    idx = idx + 1;
                end
            end

            obj.ATheta = complex(zeros(obj.Na, obj.Gtheta));
            obj.denTheta = zeros(1, obj.Gtheta);
            for k = 1:obj.Gtheta
                a = steeringVector(obj.theta_grid(k), obj.Na);
                obj.ATheta(:, k) = a;
                obj.denTheta(k) = real(a' * a);
            end

            if any(obj.denTD <= obj.AtomNormFloor)
                error('SrivastavaOMPBenchmarkEstimator:ZeroTemporalAtom', ...
                      ['At least one delay-Doppler candidate has essentially ', ...
                       'zero temporal energy. Restrict the search grid.']);
            end
            if any(obj.denTheta <= obj.AtomNormFloor)
                error('SrivastavaOMPBenchmarkEstimator:ZeroSpatialAtom', ...
                      'At least one spatial candidate has zero energy.');
            end
        end

        function [tauEst, nuEst, thetaEst, gEst, Pest] = ...
                stepImpl(obj, yTD, PExact, noiseVarianceTD)

            if PExact < 0 || PExact ~= floor(PExact)
                error('SrivastavaOMPBenchmarkEstimator:PExact', ...
                      'PExact must be a nonnegative integer.');
            end
            if noiseVarianceTD < 0 || ~isscalar(noiseVarianceTD)
                error('SrivastavaOMPBenchmarkEstimator:NoiseVariance', ...
                      'noiseVarianceTD must be a nonnegative scalar.');
            end

            yVec = obj.canonicalizeObservation(yTD);
            Y = reshape(yVec, obj.L, obj.Na);

            if PExact > 0
                Plimit = PExact;
                autoStop = false;
            else
                Plimit = obj.P_max;
                autoStop = true;
            end

            tauEst = zeros(Plimit, 1);
            nuEst = zeros(Plimit, 1);
            thetaEst = zeros(Plimit, 1);
            gEst = complex(zeros(Plimit, 1));

            PhiActive = complex(zeros(obj.L * obj.Na, Plimit));
            selectedMask = false(obj.Gdd, obj.Gtheta);

            residual = yVec;
            obj.resetDiagnostics();
            obj.LastNoiseEnergyThreshold = ...
                noiseVarianceTD * obj.L * obj.Na;
            obj.LastResidualNorm = norm(residual);
            obj.LastResidualEnergy = real(residual' * residual);

            % Paper Algorithm 1 checks the residual/noise-power condition
            % before entering the next OMP iteration.
            if autoStop && obj.LastResidualEnergy < ...
                    obj.LastNoiseEnergyThreshold
                Pest = 0;
                tauEst = zeros(0, 1);
                nuEst = zeros(0, 1);
                thetaEst = zeros(0, 1);
                gEst = complex(zeros(0, 1));
                return;
            end

            for p = 1:Plimit
                R = reshape(residual, obj.L, obj.Na);

                % Exact factorization of all full-atom correlations:
                % vec(h*a^T)^H vec(R) = h^H R conj(a).
                Q = obj.PsiTD' * R;                 % Gdd x Na
                C = Q * conj(obj.ATheta);           % Gdd x Gtheta

                if obj.NormalizeSelection
                    den = sqrt(obj.denTD * obj.denTheta);
                    scores = abs(C) ./ max(den, obj.AtomNormFloor);
                else
                    scores = abs(C);
                end

                scores(selectedMask) = -Inf;
                [bestScore, lin] = max(scores(:));
                if ~isfinite(bestScore)
                    break;
                end

                [idd, itheta] = ind2sub([obj.Gdd, obj.Gtheta], lin);
                selectedMask(idd, itheta) = true;

                h = obj.PsiTD(:, idd);
                a = obj.ATheta(:, itheta);
                atom = kron(a, h);  % vec(h*a^T)
                PhiActive(:, p) = atom;

                % Algorithm 1, Step 5: LS on the accumulated support.
                % Backslash computes the same least-squares solution as the
                % Moore-Penrose expression for a full-column-rank support,
                % while being numerically preferable to explicit inversion.
                gActive = PhiActive(:, 1:p) \ yVec;
                residual = yVec - PhiActive(:, 1:p) * gActive;

                tauEst(p) = obj.tauOfDD(idd);
                nuEst(p) = obj.nuOfDD(idd);
                thetaEst(p) = obj.theta_grid(itheta);
                gEst(1:p) = gActive;

                obj.P_detected = p;
                obj.LastResidualNorm = norm(residual);
                obj.LastResidualEnergy = real(residual' * residual);
                obj.CandidateEvaluations = obj.CandidateEvaluations + ...
                    obj.Gdd * obj.Gtheta;
                obj.SelectedLinearIndices(p, 1) = lin;
                obj.SelectedDDIndices(p, 1) = idd;
                obj.SelectedThetaIndices(p, 1) = itheta;

                if autoStop && obj.LastResidualEnergy < ...
                        obj.LastNoiseEnergyThreshold
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
            obj.resetDiagnostics();
        end

    end

    methods (Access = private)

        function h = generateTDResponse(obj, tau, nu)
            % Exact unit-gain temporal atom under the current project model.
            phase = exp(1j * 2 * pi * nu * obj.timeVec);
            xDop = phase .* obj.txSig;
            h = delayseq(xDop, tau, obj.fs);
            h = h(:);
        end

        function yVec = canonicalizeObservation(obj, yTD)
            if ismatrix(yTD) && ~isvector(yTD) && ...
                    isequal(size(yTD), [obj.L, obj.Na])
                yVec = yTD(:);
                return;
            end

            yVec = yTD(:);
            expected = obj.L * obj.Na;
            if numel(yVec) ~= expected
                error('SrivastavaOMPBenchmarkEstimator:ObservationSize', ...
                      ['Observation must be L-by-Na = %d-by-%d or contain ', ...
                       '%d stacked samples.'], obj.L, obj.Na, expected);
            end
        end

        function validateConfiguration(obj)
            if isempty(obj.x_dd) || ~isequal(size(obj.x_dd), [obj.M, obj.N])
                error('SrivastavaOMPBenchmarkEstimator:PilotSize', ...
                      'x_dd must be a nonempty M-by-N known pilot grid.');
            end
            if isempty(obj.tau_grid) || isempty(obj.nu_grid) || ...
                    isempty(obj.theta_grid)
                error('SrivastavaOMPBenchmarkEstimator:EmptyGrid', ...
                      'tau_grid, nu_grid, and theta_grid must be nonempty.');
            end
            if any(obj.tau_grid(:) < 0)
                error('SrivastavaOMPBenchmarkEstimator:NegativeDelay', ...
                      'tau_grid must contain nonnegative delays.');
            end
            if obj.Na < 1 || obj.Na ~= floor(obj.Na)
                error('SrivastavaOMPBenchmarkEstimator:Na', ...
                      'Na must be a positive integer.');
            end
            if obj.M < 1 || obj.N < 1 || ...
                    obj.M ~= floor(obj.M) || obj.N ~= floor(obj.N)
                error('SrivastavaOMPBenchmarkEstimator:FrameSize', ...
                      'M and N must be positive integers.');
            end
            if obj.fs <= 0
                error('SrivastavaOMPBenchmarkEstimator:Fs', ...
                      'fs must be positive.');
            end
            if obj.P_max < 1 || obj.P_max ~= floor(obj.P_max)
                error('SrivastavaOMPBenchmarkEstimator:Pmax', ...
                      'P_max must be a positive integer.');
            end
        end

        function resetDiagnostics(obj)
            obj.P_detected = 0;
            obj.LastResidualNorm = NaN;
            obj.LastResidualEnergy = NaN;
            obj.LastNoiseEnergyThreshold = NaN;
            obj.CandidateEvaluations = 0;
            obj.SelectedLinearIndices = zeros(0, 1);
            obj.SelectedDDIndices = zeros(0, 1);
            obj.SelectedThetaIndices = zeros(0, 1);
        end

    end
end
