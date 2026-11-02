classdef Classic3DOMPSensingEstimator < matlab.System
    % CLASSIC3DOMPSENSINGESTIMATOR Fixed-grid three-dimensional OMP for sensing.
    %
    % [tau,nu,theta,g,P,info] = step(est,yDD,PExact,sigma2DD)
    %   yDD: M*N*Na DD samples in antenna-block order.
    %   Outputs: delay [s], Doppler [Hz], angle [rad], complex gain and order.
    %
    % Precomputed temporal and spatial dictionaries avoid storing every full
    % Kronecker atom. Greedy selection uses normalized residual correlation;
    % all accepted gains are jointly refitted by LS. Coordinates remain on the
    % specified fixed grids. The demo factory uses fourfold grid oversampling.
    % See REFERENCES.bib: srivastava4DSparse2022 for the sparse OTFS baseline
    % context. This receive-array implementation is adapted to the common model.
    %
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    properties
        M (1, 1) double = 32
        N (1, 1) double = 16
        Na (1, 1) double = 6
        fs (1, 1) double = 4.8e6

        padtype (1, :) char = 'NONE'
        padlen (1, 1) double = 0
        x_dd double = []

        tau_grid double = []      % seconds
        nu_grid double = []       % Hz
        theta_grid double = []    % radians

        P_max (1, 1) double = 4

        % Normalized matched-filter selection:
        % |d^H r|^2 / ||d||^2. This is invariant to atom scaling.
        NormalizeSelection (1, 1) logical = true

        % Noise-energy stopping for blind model-order operation:
        % stop if ||r||^2 <= NoiseStopScale * L * sigma2DD.
        NoiseStopScale (1, 1) double = 1.0

        % If true, retain the selected atom matrix in info for diagnostics.
        ReturnSelectedDictionary (1, 1) logical = false
    end

    properties (SetAccess = private)
        CandidateEvaluations (1, 1) double = 0
        IterationsExecuted (1, 1) double = 0
        FinalResidualNorm (1, 1) double = NaN
    end

    properties (Access = private)
        MN (1, 1) double
        Lobs (1, 1) double
        Gtau (1, 1) double
        Gnu (1, 1) double
        Gtheta (1, 1) double
        Gdd (1, 1) double
        Gtotal (1, 1) double

        % Hdd(:,j) = h(tau_i,nu_k), j = sub2ind([Gtau,Gnu],i,k)
        Hdd double
        HddEnergy double

        % A(:,q) = a(theta_q)
        A double
        AEnergy double

        SelectedFullIndex double
        SelectedDDIndex double
        SelectedThetaIndex double
        ResidualEnergyHistory double
        MaxScoreHistory double
    end

    methods

        function obj = Classic3DOMPSensingEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.validateConfiguration();

            obj.MN = obj.M * obj.N;
            obj.Lobs = obj.MN * obj.Na;
            obj.Gtau = numel(obj.tau_grid);
            obj.Gnu = numel(obj.nu_grid);
            obj.Gtheta = numel(obj.theta_grid);
            obj.Gdd = obj.Gtau * obj.Gnu;
            obj.Gtotal = obj.Gdd * obj.Gtheta;

            % Precompute all DD atoms.  This is the reusable dictionary
            % component for a fixed pilot and fixed (tau,nu) grid.
            obj.Hdd = complex(zeros(obj.MN, obj.Gdd));
            obj.HddEnergy = zeros(obj.Gdd, 1);

            for inu = 1:obj.Gnu
                for itau = 1:obj.Gtau
                    j = sub2ind([obj.Gtau, obj.Gnu], itau, inu);
                    [yDD, ~, ~, ~] = otfsChannel( ...
                                              obj.x_dd, obj.M, obj.fs, ...
                                              obj.tau_grid(itau), obj.nu_grid(inu), 1, ...
                                              'none', obj.padtype, obj.padlen);
                    h = yDD(:);
                    obj.Hdd(:, j) = h;
                    obj.HddEnergy(j) = real(h' * h);
                end
            end

            obj.A = complex(zeros(obj.Na, obj.Gtheta));
            obj.AEnergy = zeros(obj.Gtheta, 1);
            for ia = 1:obj.Gtheta
                a = steeringVector(obj.theta_grid(ia), obj.Na);
                obj.A(:, ia) = a;
                obj.AEnergy(ia) = real(a' * a);
            end

            if any(obj.HddEnergy <= eps)
                error('Classic3DOMP:ZeroDDAtom', ...
                      'At least one delay-Doppler atom has near-zero energy.');
            end
            if any(obj.AEnergy <= eps)
                error('Classic3DOMP:ZeroSpatialAtom', ...
                      'At least one spatial atom has near-zero energy.');
            end

            obj.resetDiagnostics();
        end

        function [tau_est, nu_est, theta_est, g_est, P_est, info] = ...
                stepImpl(obj, y_vec, P_exact, sigma2DD)

            y = y_vec(:);
            if numel(y) ~= obj.Lobs
                error('Classic3DOMP:ObservationSize', ...
                      'Expected M*N*Na=%d observation samples, got %d.', ...
                      obj.Lobs, numel(y));
            end
            if ~all(isfinite(y))
                error('Classic3DOMP:NonFiniteObservation', ...
                      'Observation contains NaN or Inf.');
            end

            if nargin < 3 || isempty(P_exact)
                P_exact = 0;
            end
            if nargin < 4 || isempty(sigma2DD)
                sigma2DD = NaN;
            end

            if P_exact < 0 || abs(P_exact - round(P_exact)) > 0
                error('Classic3DOMP:InvalidPExact', ...
                      'P_exact must be a nonnegative integer.');
            end

            if P_exact > 0
                maxIter = min(round(P_exact), round(obj.P_max));
            else
                maxIter = round(obj.P_max);
            end

            obj.resetDiagnostics();

            r = y;
            supportFull = zeros(maxIter, 1);
            supportDD = zeros(maxIter, 1);
            supportTheta = zeros(maxIter, 1);
            B = complex(zeros(obj.Lobs, maxIter));
            g = complex(zeros(maxIter, 1));

            residualEnergy = nan(maxIter + 1, 1);
            residualEnergy(1) = real(r' * r);
            maxScore = nan(maxIter, 1);

            noiseThreshold = NaN;
            if P_exact == 0 && isfinite(sigma2DD) && sigma2DD >= 0
                noiseThreshold = obj.NoiseStopScale * obj.Lobs * sigma2DD;
                if residualEnergy(1) <= noiseThreshold
                    tau_est = zeros(0, 1);
                    nu_est = zeros(0, 1);
                    theta_est = zeros(0, 1);
                    g_est = complex(zeros(0, 1));
                    P_est = 0;
                    info = obj.makeInfo(y, r, B(:, []), g_est, ...
                                        supportFull([]), supportDD([]), supportTheta([]), ...
                                        residualEnergy(1), maxScore([]), noiseThreshold);
                    return;
                end
            end

            P_est = 0;

            for t = 1:maxIter
                % Reshape according to project convention:
                % Y is (MN)-by-Na, columns are antenna observations.
                Yres = reshape(r, obj.MN, obj.Na);

                % EXACT factorized implementation of full 3-D correlations:
                % C(j,q) = h_j^H Yres conj(a_q)
                C = obj.Hdd' * Yres * conj(obj.A);

                if obj.NormalizeSelection
                    Den = obj.HddEnergy * obj.AEnergy.';
                    score = abs(C).^2 ./ max(Den, realmin('double'));
                else
                    score = abs(C).^2;
                end

                % OMP cannot reselect a previously active grid atom.
                if t > 1
                    score(supportFull(1:t - 1)) = -Inf;
                end

                [bestScore, idxFull] = max(score(:));
                if ~isfinite(bestScore)
                    break;
                end

                [jDD, ia] = ind2sub([obj.Gdd, obj.Gtheta], idxFull);
                supportFull(t) = idxFull;
                supportDD(t) = jDD;
                supportTheta(t) = ia;
                maxScore(t) = bestScore;

                % Add selected physical atom.
                B(:, t) = kron(obj.A(:, ia), obj.Hdd(:, jDD));

                % Orthogonal projection / joint LS coefficient update.
                Bt = B(:, 1:t);
                G = Bt' * Bt;
                if rcond(G) < 1e-12 && exist('lsqminnorm', 'file') == 2
                    g(1:t) = lsqminnorm(Bt, y);
                else
                    g(1:t) = Bt \ y;
                end

                r = y - Bt * g(1:t);
                residualEnergy(t + 1) = real(r' * r);

                obj.CandidateEvaluations = obj.CandidateEvaluations + obj.Gtotal;
                obj.IterationsExecuted = t;
                P_est = t;

                if P_exact == 0 && isfinite(noiseThreshold) && ...
                        residualEnergy(t + 1) <= noiseThreshold
                    break;
                end
            end

            supportFull = supportFull(1:P_est);
            supportDD = supportDD(1:P_est);
            supportTheta = supportTheta(1:P_est);
            B = B(:, 1:P_est);
            g_est = g(1:P_est);
            residualEnergy = residualEnergy(1:P_est + 1);
            maxScore = maxScore(1:P_est);

            tau_est = zeros(P_est, 1);
            nu_est = zeros(P_est, 1);
            theta_est = zeros(P_est, 1);
            for t = 1:P_est
                [itau, inu] = ind2sub([obj.Gtau, obj.Gnu], supportDD(t));
                tau_est(t) = obj.tau_grid(itau);
                nu_est(t) = obj.nu_grid(inu);
                theta_est(t) = obj.theta_grid(supportTheta(t));
            end

            obj.SelectedFullIndex = supportFull;
            obj.SelectedDDIndex = supportDD;
            obj.SelectedThetaIndex = supportTheta;
            obj.ResidualEnergyHistory = residualEnergy;
            obj.MaxScoreHistory = maxScore;
            obj.FinalResidualNorm = norm(r);

            info = obj.makeInfo(y, r, B, g_est, supportFull, supportDD, ...
                                supportTheta, residualEnergy, maxScore, noiseThreshold);
        end

        function resetImpl(obj)
            obj.resetDiagnostics();
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
                names = {'MN', 'Lobs', 'Gtau', 'Gnu', 'Gtheta', 'Gdd', 'Gtotal', ...
                         'Hdd', 'HddEnergy', 'A', 'AEnergy', ...
                         'SelectedFullIndex', 'SelectedDDIndex', 'SelectedThetaIndex', ...
                         'ResidualEnergyHistory', 'MaxScoreHistory'};
                for k = 1:numel(names)
                    s.(names{k}) = obj.(names{k});
                end
                s.CandidateEvaluations = obj.CandidateEvaluations;
                s.IterationsExecuted = obj.IterationsExecuted;
                s.FinalResidualNorm = obj.FinalResidualNorm;
            end
        end

        function loadObjectImpl(obj, s, wasLocked)
            if wasLocked
                names = {'MN', 'Lobs', 'Gtau', 'Gnu', 'Gtheta', 'Gdd', 'Gtotal', ...
                         'Hdd', 'HddEnergy', 'A', 'AEnergy', ...
                         'SelectedFullIndex', 'SelectedDDIndex', 'SelectedThetaIndex', ...
                         'ResidualEnergyHistory', 'MaxScoreHistory'};
                for k = 1:numel(names)
                    obj.(names{k}) = s.(names{k});
                end
                obj.CandidateEvaluations = s.CandidateEvaluations;
                obj.IterationsExecuted = s.IterationsExecuted;
                obj.FinalResidualNorm = s.FinalResidualNorm;
            end
            loadObjectImpl@matlab.System(obj, s, wasLocked);
        end

    end

    methods (Access = private)

        function validateConfiguration(obj)
            if isempty(obj.x_dd) || ~isequal(size(obj.x_dd), [obj.M obj.N])
                error('Classic3DOMP:PilotSize', ...
                      'x_dd must be M-by-N = %d-by-%d.', obj.M, obj.N);
            end
            if isempty(obj.tau_grid) || isempty(obj.nu_grid) || isempty(obj.theta_grid)
                error('Classic3DOMP:EmptyGrid', ...
                      'tau_grid, nu_grid, and theta_grid must be nonempty.');
            end
            if any(~isfinite(obj.tau_grid)) || any(obj.tau_grid < 0)
                error('Classic3DOMP:TauGrid', ...
                      'tau_grid must contain finite nonnegative delays [s].');
            end
            if any(~isfinite(obj.nu_grid))
                error('Classic3DOMP:NuGrid', ...
                      'nu_grid must contain finite Dopplers [Hz].');
            end
            if any(~isfinite(obj.theta_grid))
                error('Classic3DOMP:ThetaGrid', ...
                      'theta_grid must contain finite angles [rad].');
            end
            if obj.P_max < 1 || abs(obj.P_max - round(obj.P_max)) > 0
                error('Classic3DOMP:Pmax', 'P_max must be a positive integer.');
            end
        end

        function resetDiagnostics(obj)
            obj.CandidateEvaluations = 0;
            obj.IterationsExecuted = 0;
            obj.FinalResidualNorm = NaN;
            obj.SelectedFullIndex = zeros(0, 1);
            obj.SelectedDDIndex = zeros(0, 1);
            obj.SelectedThetaIndex = zeros(0, 1);
            obj.ResidualEnergyHistory = zeros(0, 1);
            obj.MaxScoreHistory = zeros(0, 1);
        end

        function info = makeInfo(obj, y, r, B, g, sFull, sDD, sTheta, ...
                                 residualEnergy, maxScore, noiseThreshold)
            info = struct;
            info.algorithm = 'conventional_3d_omp_common_model';
            info.gridSize = [obj.Gtau obj.Gnu obj.Gtheta];
            info.dictionaryAtomCount = obj.Gtotal;
            info.observationLength = obj.Lobs;
            info.candidateEvaluations = obj.CandidateEvaluations;
            info.iterations = obj.IterationsExecuted;
            info.supportFullIndex = sFull(:);
            info.supportDDIndex = sDD(:);
            info.supportThetaIndex = sTheta(:);
            info.residualEnergy = residualEnergy(:);
            info.maxScore = maxScore(:);
            info.finalResidualNorm = norm(r);
            info.noiseEnergyThreshold = noiseThreshold;
            info.normalizeSelection = obj.NormalizeSelection;
            info.yHat = B * g;
            info.responseResidualNMSE = real(r' * r) / max(real(y' * y), realmin);

            if size(B, 2) >= 2
                d = sqrt(max(real(diag(B' * B)), realmin));
                Gn = (B' * B) ./ (d * d.');
                Gn(1:size(Gn, 1) + 1:end) = 0;
                info.maxSelectedAtomCoherence = max(abs(Gn), [], 'all');
            else
                info.maxSelectedAtomCoherence = 0;
            end

            if obj.ReturnSelectedDictionary
                info.selectedDictionary = B;
            end
        end

    end
end
