classdef Li3DNOMPCommonHostEstimator < matlab.System
    % LI3DNOMPCOMMONHOSTESTIMATOR Three-dimensional NOMP on time-domain samples.
    %
    % [tau,nu,theta,g,P,info] = step(est,yTD,PExact,sigma2TD)
    %   yTD: L-by-Na samples, or their column-major vector, L=M*N without padding.
    %   sigma2TD: variance per complex time-domain sample, not DD variance.
    %   Outputs: seconds, Hz, radians, complex gains, order and diagnostics.
    %
    % The paper uses AtomConvention='project': Doppler multiplication precedes
    % finite-frame delayseq, and the receive response uses steeringVector.
    % Coarse acquisition is followed by a local angle/Doppler search, Rs Newton
    % updates, bounded fractional-delay refinement and at most Rc cyclic sweeps.
    % Backtracking protects the single-path score; worsening cyclic steps roll
    % back. Gains are fitted jointly by ordinary LS.
    %
    % The fractional-delay extension and centered array convention adapt Li's
    % 3-D construction to the same physical waveform as MDDIPIC. The optional
    % 'li-source' convention retains integer delay and an uncentered ULA; it is
    % not used by run_capsule_demo.
    % Project implementation of the NOMP construction and Li 3-D adaptation;
    % see REFERENCES.bib: mamandipoorNOMP2016 and liPDMA2021.
    %
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    properties
        M = 64
        N = 16
        Na = 16
        fs = 10e6
        padtype = 'NONE'
        padlen = 0
        x_dd = []
        tau_grid = []
        nu_coarse_grid = []
        theta_coarse_grid = []
        RhoNu = 4
        RhoTheta = 4
        Rs = 4
        Rc = 5
        P_max = 10
        Pfa = 1e-3
        AtomConvention = 'project'
        EnableFractionalDelayRefinement = true
        DelaySearchHalfWidthBins = 0.5
        DelayTolerance = 1e-12
        DelayMaxFunEvals = 60
        HessianRcondFloor = 1e-10
        AtomNormFloor = 1e-14
        NewtonMaxBacktracks = 8
    end

    properties (SetAccess = private)
        P_detected = 0
        LastResidualEnergy = NaN
        LastStoppingThreshold = NaN
        CandidateEvaluations = 0
        NewtonUpdates = 0
        CoarseSearchEvaluations = 0
        PreciseSearchEvaluations = 0
        DelayRefinementEvaluations = 0
        DelayRefinementAccepted = 0
        ResidualEnergyByPath = []
    end

    properties (Access = private)
        txSig
        timeVec
        spatialIndex
        L
        PsiTDCoarse
        denTDCoarse
        AThetaCoarse
        denThetaCoarse
        thetaMin
        thetaMax
        nuMin
        nuMax
        tauMin
        tauMax
        dThetaCoarse
        dNuCoarse
        dTauCoarse
    end

    methods

        function obj = Li3DNOMPCommonHostEstimator(varargin)
            setProperties(obj, nargin, varargin{:});
        end

        function p = debugAtom(obj, tau, nu, theta)
            p = obj.atom(tau, nu, theta);
        end

        function [p, pTheta, pNu, pThetaTheta, pThetaNu, pNuNu] = ...
                debugAtomDerivatives(obj, tau, nu, theta)
            [p, pTheta, pNu, pThetaTheta, pThetaNu, pNuNu] = ...
                obj.atomWithDerivatives(tau, nu, theta);
        end

        function [J, grad, H] = debugObjectiveDerivatives( ...
                                                        obj, y, tau, nu, theta, g)
            [p, pT, pN, pTT, pTN, pNN] = ...
                obj.atomWithDerivatives(tau, nu, theta);
            [J, grad, H] = obj.objectiveDerivatives(y(:), g, p, pT, pN, pTT, pTN, pNN);
        end

        function [tau, nu, theta, score] = debugCoarseSearch(obj, residual)
            [tau, nu, theta, score] = obj.coarseSearch(residual(:));
        end

        function [nu, theta, g, score] = ...
                debugPreciseSearch(obj, residual, tau, nu0, theta0)
            [nu, theta, g, score] = ...
                obj.preciseSearch(residual(:), tau, nu0, theta0);
        end

        function [theta, nu, g, history] = ...
                debugNewtonRefine(obj, y, tau, theta, nu, g, nIter)
            [theta, nu, g, history] = ...
                obj.newtonRefineSingle(y(:), tau, theta, nu, g, nIter);
        end

        function [tau, g, scoreBefore, scoreAfter] = ...
                debugDelayRefine(obj, y, tau, nu, theta)
            [tau, g, scoreBefore, scoreAfter] = ...
                obj.refineDelay(y(:), tau, nu, theta);
        end

        function g = debugJointGain(obj, y, tau, nu, theta)
            g = obj.jointGainUpdate(y(:), tau, nu, theta);
        end

        function Phi = debugActiveDictionary(obj, tau, nu, theta)
            Phi = obj.activeDictionary(tau, nu, theta);
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
            if strcmpi(obj.AtomConvention, 'li-source')
                obj.spatialIndex = (0:obj.Na - 1).';
            else
                obj.spatialIndex = (0:obj.Na - 1).' - obj.Na / 2;
            end

            obj.thetaMin = min(obj.theta_coarse_grid);
            obj.thetaMax = max(obj.theta_coarse_grid);
            obj.nuMin = min(obj.nu_coarse_grid);
            obj.nuMax = max(obj.nu_coarse_grid);
            obj.tauMin = min(obj.tau_grid);
            obj.tauMax = max(obj.tau_grid);
            obj.dThetaCoarse = obj.gridSpacing(obj.theta_coarse_grid);
            obj.dNuCoarse = obj.gridSpacing(obj.nu_coarse_grid);
            obj.dTauCoarse = obj.gridSpacing(obj.tau_grid);

            Gtau = numel(obj.tau_grid);
            Gnu = numel(obj.nu_coarse_grid);
            obj.PsiTDCoarse = complex(zeros(obj.L, Gtau * Gnu));
            obj.denTDCoarse = zeros(Gtau * Gnu, 1);
            idx = 1;
            for inu = 1:Gnu
                for itau = 1:Gtau
                    h = obj.temporalAtom( ...
                                         obj.tau_grid(itau), obj.nu_coarse_grid(inu));
                    obj.PsiTDCoarse(:, idx) = h;
                    obj.denTDCoarse(idx) = real(h' * h);
                    idx = idx + 1;
                end
            end

            obj.AThetaCoarse = complex(zeros(obj.Na, numel(obj.theta_coarse_grid)));
            obj.denThetaCoarse = zeros(1, numel(obj.theta_coarse_grid));
            for ith = 1:numel(obj.theta_coarse_grid)
                a = obj.spatialAtom(obj.theta_coarse_grid(ith));
                obj.AThetaCoarse(:, ith) = a;
                obj.denThetaCoarse(ith) = real(a' * a);
            end
        end

        function [tauEst, nuEst, thetaEst, gEst, Pest, info] = ...
                stepImpl(obj, yTD, PExact, noiseVarianceTD)
            obj.validateStepInputs(PExact, noiseVarianceTD);
            y = obj.canonicalizeObservation(yTD);
            if PExact > 0
                Ptarget = PExact;
                autoStop = false;
            else
                Ptarget = obj.P_max;
                autoStop = true;
            end
            tauEst = zeros(Ptarget, 1);
            nuEst = zeros(Ptarget, 1);
            thetaEst = zeros(Ptarget, 1);
            gEst = complex(zeros(Ptarget, 1));
            residual = y;
            obj.resetDiagnostics();
            if noiseVarianceTD > 0
                Kobs = obj.Na * obj.L;
                obj.LastStoppingThreshold = -noiseVarianceTD * ...
                    log(1 - (1 - obj.Pfa)^(1 / Kobs));
            else
                obj.LastStoppingThreshold = 0;
            end

            newtonHistory = cell(Ptarget, 1);
            cycleResidualHistory = cell(Ptarget, 1);
            residualHistory = zeros(Ptarget + 1, 1);
            residualHistory(1) = real(y' * y);
            for ip = 1:Ptarget
                [tau0, nu0, theta0, bestScore] = obj.coarseSearch(residual);
                if autoStop && bestScore < obj.LastStoppingThreshold
                    break;
                end
                [nu1, theta1, g1] = obj.preciseSearch( ...
                                                    residual, tau0, nu0, theta0);
                [theta1, nu1, g1, h1] = obj.newtonRefineSingle( ...
                                                            residual, tau0, theta1, nu1, g1, obj.Rs);

                tau1 = tau0;
                if obj.EnableFractionalDelayRefinement
                    [tau1, g1] = obj.refineDelay( ...
                                                residual, tau0, nu1, theta1);
                    % MODIFIED alternating cleanup: one theta/nu Newton
                    % update after delay movement, then one final delay fit.
                    [theta1, nu1, g1, h2] = obj.newtonRefineSingle( ...
                                                                residual, tau1, theta1, nu1, g1, 1);
                    [tau1, g1] = obj.refineDelay( ...
                                                residual, tau1, nu1, theta1);
                    h1 = [h1; h2]; %#ok<AGROW>
                end

                tauEst(ip) = tau1;
                nuEst(ip) = nu1;
                thetaEst(ip) = theta1;
                gEst(ip) = g1;
                gEst(1:ip) = obj.jointGainUpdate( ...
                                                 y, tauEst(1:ip), nuEst(1:ip), thetaEst(1:ip));

                for ic = 1:obj.Rc
                    paramsChanged = false;
                    residualBeforeCycle = norm(y - obj.activeDictionary( ...
                                                                      tauEst(1:ip), nuEst(1:ip), thetaEst(1:ip)) * gEst(1:ip));
                    oldTau = tauEst(1:ip);
                    oldNu = nuEst(1:ip);
                    oldTheta = thetaEst(1:ip);
                    oldGain = gEst(1:ip);
                    for jp = 1:ip
                        yMinus = y;
                        for kp = 1:ip
                            if kp ~= jp
                                yMinus = yMinus - gEst(kp) * obj.atom( ...
                                                                  tauEst(kp), nuEst(kp), thetaEst(kp));
                            end
                        end
                        p0 = obj.atom(tauEst(jp), nuEst(jp), thetaEst(jp));
                        gj = (p0' * yMinus) / max(real(p0' * p0), obj.AtomNormFloor);
                        [thNew, nuNew, gj] = obj.newtonRefineSingle( ...
                                                                  yMinus, tauEst(jp), thetaEst(jp), nuEst(jp), gj, 1);
                        tauNew = tauEst(jp);
                        if obj.EnableFractionalDelayRefinement
                            [tauNew, gj] = obj.refineDelay( ...
                                                          yMinus, tauNew, nuNew, thNew);
                        end
                        paramsChanged = paramsChanged || ...
                            abs(tauNew - tauEst(jp)) > obj.DelayTolerance || ...
                            abs(nuNew - nuEst(jp)) > 1e-9 || ...
                            abs(thNew - thetaEst(jp)) > 1e-13;
                        tauEst(jp) = tauNew;
                        nuEst(jp) = nuNew;
                        thetaEst(jp) = thNew;
                    end
                    gEst(1:ip) = obj.jointGainUpdate( ...
                                                     y, tauEst(1:ip), nuEst(1:ip), thetaEst(1:ip));
                    residualAfterCycle = norm(y - obj.activeDictionary( ...
                                                                     tauEst(1:ip), nuEst(1:ip), thetaEst(1:ip)) * gEst(1:ip));
                    if residualAfterCycle > residualBeforeCycle * (1 + 1e-12)
                        tauEst(1:ip) = oldTau;
                        nuEst(1:ip) = oldNu;
                        thetaEst(1:ip) = oldTheta;
                        gEst(1:ip) = oldGain;
                        cycleResidualHistory{ip}(end + 1, :) = ...
                            [residualBeforeCycle residualBeforeCycle]; %#ok<AGROW>
                        break;
                    end
                    cycleResidualHistory{ip}(end + 1, :) = ...
                        [residualBeforeCycle residualAfterCycle]; %#ok<AGROW>
                    if ~paramsChanged,
                        break;
                    end
                end

                Phi = obj.activeDictionary( ...
                                           tauEst(1:ip), nuEst(1:ip), thetaEst(1:ip));
                residual = y - Phi * gEst(1:ip);
                obj.P_detected = ip;
                obj.LastResidualEnergy = real(residual' * residual);
                obj.ResidualEnergyByPath(ip, 1) = obj.LastResidualEnergy;
                residualHistory(ip + 1) = obj.LastResidualEnergy;
                newtonHistory{ip} = h1;
            end
            Pest = obj.P_detected;
            tauEst = tauEst(1:Pest);
            nuEst = nuEst(1:Pest);
            thetaEst = thetaEst(1:Pest);
            gEst = gEst(1:Pest);
            Phi = obj.activeDictionary(tauEst, nuEst, thetaEst);
            yHat = Phi * gEst;
            info = struct('yHat', yHat, 'residual', y - yHat, ...
                          'residualEnergyByPath', obj.ResidualEnergyByPath, ...
                          'residualEnergyHistory', residualHistory(1:Pest + 1), ...
                          'newtonHistory', {newtonHistory(1:Pest)}, ...
                          'cycleResidualHistory', {cycleResidualHistory(1:Pest)}, ...
                          'newtonUpdates', obj.NewtonUpdates, ...
                          'delayRefinementEvaluations', obj.DelayRefinementEvaluations, ...
                          'delayRefinementAccepted', obj.DelayRefinementAccepted, ...
                          'candidateEvaluations', obj.CandidateEvaluations, ...
                          'coarseSearchEvaluations', obj.CoarseSearchEvaluations, ...
                          'preciseSearchEvaluations', obj.PreciseSearchEvaluations, ...
                          'atomConvention', obj.AtomConvention, ...
                          'fractionalDelayRefinement', obj.EnableFractionalDelayRefinement);
        end

        function resetImpl(obj)
            obj.resetDiagnostics();
        end

    end

    methods (Access = private)

        function [tauBest, nuBest, thetaBest, bestScore] = ...
                coarseSearch(obj, residual)
            R = reshape(residual, obj.L, obj.Na);
            C = (obj.PsiTDCoarse' * R) * conj(obj.AThetaCoarse);
            den = obj.denTDCoarse * obj.denThetaCoarse;
            score = abs(C).^2 ./ max(den, obj.AtomNormFloor);
            [bestScore, lin] = max(score(:));
            [idd, ith] = ind2sub(size(score), lin);
            Gtau = numel(obj.tau_grid);
            inu = floor((idd - 1) / Gtau) + 1;
            itau = mod(idd - 1, Gtau) + 1;
            tauBest = obj.tau_grid(itau);
            nuBest = obj.nu_coarse_grid(inu);
            thetaBest = obj.theta_coarse_grid(ith);
            nEval = numel(score);
            obj.CandidateEvaluations = obj.CandidateEvaluations + nEval;
            obj.CoarseSearchEvaluations = obj.CoarseSearchEvaluations + nEval;
        end

        function [nuBest, thetaBest, gBest, bestScore] = ...
                preciseSearch(obj, residual, tau, nu0, theta0)
            nuCand = obj.localFineGrid(nu0, obj.dNuCoarse, obj.RhoNu, obj.nuMin, obj.nuMax);
            thCand = obj.localFineGrid(theta0, obj.dThetaCoarse, obj.RhoTheta, obj.thetaMin, obj.thetaMax);
            bestScore = -Inf;
            nuBest = nu0;
            thetaBest = theta0;
            gBest = 0;
            for inu = 1:numel(nuCand)
                h = obj.temporalAtom(tau, nuCand(inu));
                for ith = 1:numel(thCand)
                    p = kron(obj.spatialAtom(thCand(ith)), h);
                    den = real(p' * p);
                    z = p' * residual;
                    score = abs(z)^2 / max(den, obj.AtomNormFloor);
                    if score > bestScore
                        bestScore = score;
                        nuBest = nuCand(inu);
                        thetaBest = thCand(ith);
                        gBest = z / max(den, obj.AtomNormFloor);
                    end
                end
            end
            nEval = numel(nuCand) * numel(thCand);
            obj.CandidateEvaluations = obj.CandidateEvaluations + nEval;
            obj.PreciseSearchEvaluations = obj.PreciseSearchEvaluations + nEval;
        end

        function [theta, nu, g, history] = newtonRefineSingle( ...
                                                           obj, y, tau, theta, nu, g, nIter)
            history = zeros(nIter + 1, 1);
            history(1) = obj.concentratedScore(y, obj.atom(tau, nu, theta));
            used = 0;
            for ir = 1:nIter
                [p, pT, pN, pTT, pTN, pNN] = obj.atomWithDerivatives(tau, nu, theta);
                [~, grad, H] = obj.objectiveDerivatives(y, g, p, pT, pN, pTT, pTN, pNN);
                if any(~isfinite(H), 'all') || any(~isfinite(grad)),
                    break;
                end
                if rcond(H) >= obj.HessianRcondFloor
                    delta = H \ grad;
                else
                    delta = pinv(H) * grad;
                end
                if any(~isfinite(delta)),
                    break;
                end
                oldScore = obj.concentratedScore(y, p);
                accepted = false;
                scale = 1;
                for ib = 0:obj.NewtonMaxBacktracks
                    thNew = min(max(theta - scale * delta(1), obj.thetaMin), obj.thetaMax);
                    nuNew = min(max(nu - scale * delta(2), obj.nuMin), obj.nuMax);
                    pNew = obj.atom(tau, nuNew, thNew);
                    newScore = obj.concentratedScore(y, pNew);
                    if newScore + 1e-12 * max(1, oldScore) >= oldScore
                        den = real(pNew' * pNew);
                        g = (pNew' * y) / max(den, obj.AtomNormFloor);
                        theta = thNew;
                        nu = nuNew;
                        accepted = true;
                        break;
                    end
                    scale = scale / 2;
                end
                if ~accepted,
                    break;
                end
                used = used + 1;
                obj.NewtonUpdates = obj.NewtonUpdates + 1;
                history(used + 1) = newScore;
            end
            history = history(1:used + 1);
        end

        function [tauBest, gBest, score0, scoreBest] = ...
                refineDelay(obj, y, tau0, nu, theta)
            p0 = obj.atom(tau0, nu, theta);
            score0 = obj.concentratedScore(y, p0);
            tauBest = tau0;
            scoreBest = score0;
            den0 = real(p0' * p0);
            gBest = (p0' * y) / max(den0, obj.AtomNormFloor);
            if ~obj.EnableFractionalDelayRefinement || ...
                    ~isfinite(obj.dTauCoarse) || obj.dTauCoarse <= 0
                return;
            end
            lo = max(obj.tauMin, tau0 - obj.DelaySearchHalfWidthBins * obj.dTauCoarse);
            hi = min(obj.tauMax, tau0 + obj.DelaySearchHalfWidthBins * obj.dTauCoarse);
            if hi - lo <= obj.DelayTolerance,
                return;
            end
            fun = @(tau)-obj.concentratedScore(y, obj.atom(tau, nu, theta));
            options = optimset('Display', 'off', 'TolX', obj.DelayTolerance, ...
                             'MaxFunEvals', obj.DelayMaxFunEvals);
            [candidate, fval, ~, output] = fminbnd(fun, lo, hi, options);
            obj.DelayRefinementEvaluations = obj.DelayRefinementEvaluations + output.funcCount;
            candidateScore = -fval;
            if candidateScore + 1e-12 * max(1, score0) >= score0
                tauBest = candidate;
                scoreBest = candidateScore;
                p = obj.atom(tauBest, nu, theta);
                gBest = (p' * y) / max(real(p' * p), obj.AtomNormFloor);
                if abs(tauBest - tau0) > obj.DelayTolerance
                    obj.DelayRefinementAccepted = obj.DelayRefinementAccepted + 1;
                end
            end
        end

        function [J, grad, H] = objectiveDerivatives( ...
                                                   ~, y, g, p, pT, pN, pTT, pTN, pNN)
            e = y - g * p;
            J = 2 * real(g * (y' * p)) - abs(g)^2 * real(p' * p);
            grad = [2 * real(g * (e' * pT)); 2 * real(g * (e' * pN))];
            H11 = 2 * real(g * (e' * pTT) - abs(g)^2 * (pT' * pT));
            H12 = 2 * real(g * (e' * pTN) - abs(g)^2 * (pN' * pT));
            H21 = 2 * real(g * (e' * pTN) - abs(g)^2 * (pT' * pN));
            H22 = 2 * real(g * (e' * pNN) - abs(g)^2 * (pN' * pN));
            H = [H11 H12; H21 H22];
        end

        function s = concentratedScore(obj, y, p)
            s = abs(p' * y)^2 / max(real(p' * p), obj.AtomNormFloor);
        end

        function g = jointGainUpdate(obj, y, tau, nu, theta)
            Phi = obj.activeDictionary(tau, nu, theta);
            g = lsqminnorm(Phi, y);
        end

        function Phi = activeDictionary(obj, tau, nu, theta)
            Phi = complex(zeros(obj.L * obj.Na, numel(tau)));
            for k = 1:numel(tau),
                Phi(:, k) = obj.atom(tau(k), nu(k), theta(k));
            end
        end

        function p = atom(obj, tau, nu, theta)
            p = kron(obj.spatialAtom(theta), obj.temporalAtom(tau, nu));
        end

        function [p, pT, pN, pTT, pTN, pNN] = atomWithDerivatives(obj, tau, nu, theta)
            [h, hN, hNN] = obj.temporalAtomDerivatives(tau, nu);
            n = obj.spatialIndex;
            a = exp(1j * pi * sin(theta) * n);
            aT = (1j * pi * n * cos(theta)) .* a;
            aTT = (-1j * pi * n * sin(theta) - (pi * n * cos(theta)).^2) .* a;
            p = kron(a, h);
            pT = kron(aT, h);
            pN = kron(a, hN);
            pTT = kron(aTT, h);
            pTN = kron(aT, hN);
            pNN = kron(a, hNN);
        end

        function a = spatialAtom(obj, theta)
            a = exp(1j * pi * sin(theta) * obj.spatialIndex);
        end

        function h = temporalAtom(obj, tau, nu)
            if strcmpi(obj.AtomConvention, 'li-source')
                delayed = delayseq(obj.txSig, tau, obj.fs);
                delayed = delayed(:);
                h = exp(1j * 2 * pi * nu * obj.timeVec) .* delayed;
            else
                h = delayseq(exp(1j * 2 * pi * nu * obj.timeVec) .* obj.txSig, tau, obj.fs);
                h = h(:);
            end
        end

        function [h, hN, hNN] = temporalAtomDerivatives(obj, tau, nu)
            phase = exp(1j * 2 * pi * nu * obj.timeVec);
            if strcmpi(obj.AtomConvention, 'li-source')
                delayed = delayseq(obj.txSig, tau, obj.fs);
                delayed = delayed(:);
                h = phase .* delayed;
                hN = (1j * 2 * pi * obj.timeVec) .* h;
                hNN = -(2 * pi * obj.timeVec).^2 .* h;
            else
                x0 = phase .* obj.txSig;
                x1 = (1j * 2 * pi * obj.timeVec) .* x0;
                x2 = -(2 * pi * obj.timeVec).^2 .* x0;
                h = delayseq(x0, tau, obj.fs);
                h = h(:);
                hN = delayseq(x1, tau, obj.fs);
                hN = hN(:);
                hNN = delayseq(x2, tau, obj.fs);
                hNN = hNN(:);
            end
        end

        function g = localFineGrid(~, x0, d, rho, xmin, xmax)
            if ~isfinite(d) || d <= 0 || rho < 1,
                g = min(max(x0, xmin), xmax);
                return;
            end
            step = d / rho;
            lo = max(xmin, x0 - d / 2);
            hi = min(xmax, x0 + d / 2);
            g = lo:step:hi;
            if isempty(g) || abs(g(end) - hi) > 1e-12 * max(1, abs(hi)),
                g = [g hi];
            end
            if all(abs(g - x0) > 1e-12 * max(1, abs(x0))),
                g = sort([g x0]);
            end
        end

        function d = gridSpacing(~, g)
            g = sort(unique(g(:)));
            if numel(g) < 2,
                d = Inf;
            else,
                d = median(diff(g));
            end
        end

        function y = canonicalizeObservation(obj, yTD)
            y = yTD(:);
            if numel(y) ~= obj.L * obj.Na
                error('Li3DNOMPCommonHostEstimator:ObservationSize', ...
                      'Expected %d stacked TD samples.', obj.L * obj.Na);
            end
        end

        function validateConfiguration(obj)
            if isempty(obj.x_dd) || ~isequal(size(obj.x_dd), [obj.M, obj.N])
                error('Li3DNOMPCommonHostEstimator:PilotSize', 'x_dd must be M-by-N.');
            end
            if isempty(obj.tau_grid) || isempty(obj.nu_coarse_grid) || isempty(obj.theta_coarse_grid)
                error('Li3DNOMPCommonHostEstimator:Grid', 'All grids must be nonempty.');
            end
            if ~any(strcmpi(obj.AtomConvention, {'project', 'li-source'}))
                error('Li3DNOMPCommonHostEstimator:AtomConvention', ...
                      'AtomConvention must be project or li-source.');
            end
            if obj.Rs < 0 || obj.Rc < 0 || obj.Rs ~= floor(obj.Rs) || obj.Rc ~= floor(obj.Rc)
                error('Li3DNOMPCommonHostEstimator:Iterations', 'Rs/Rc must be integers.');
            end
            if obj.DelaySearchHalfWidthBins <= 0 || obj.DelayTolerance <= 0
                error('Li3DNOMPCommonHostEstimator:DelayRefinement', 'Invalid delay settings.');
            end
        end

        function validateStepInputs(~, PExact, noiseVarianceTD)
            if ~isscalar(PExact) || PExact < 0 || PExact ~= floor(PExact)
                error('Li3DNOMPCommonHostEstimator:PExact', 'PExact must be a nonnegative integer.');
            end
            if ~isscalar(noiseVarianceTD) || noiseVarianceTD < 0
                error('Li3DNOMPCommonHostEstimator:NoiseVariance', 'Invalid noise variance.');
            end
        end

        function resetDiagnostics(obj)
            obj.P_detected = 0;
            obj.LastResidualEnergy = NaN;
            obj.LastStoppingThreshold = NaN;
            obj.CandidateEvaluations = 0;
            obj.NewtonUpdates = 0;
            obj.CoarseSearchEvaluations = 0;
            obj.PreciseSearchEvaluations = 0;
            obj.DelayRefinementEvaluations = 0;
            obj.DelayRefinementAccepted = 0;
            obj.ResidualEnergyByPath = [];
        end

    end
end
