function assoc = associate_isac_estimates_v2(truth, tauHat, nuHat, thetaHat, cfg, domain)
    % ASSOCIATE_ISAC_ESTIMATES_V2 Truth-power-aware one-to-one assignment.
    %
    % Communication paths use truth-power-weighted rows in the normalized
    % geometric assignment cost. This protects dominant channel components
    % from being exchanged with much weaker paths during one-to-one matching.
    % The rule is valid for both the former Rician truth and the TDL-B Rayleigh
    % truth; it does not assume that a designated LoS component exists.
    % Sensing retains the equal-target V1 geometric policy.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    tauT = truth.delays_sec(:);
    nuT = truth.dopplers_Hz(:);
    thT = truth.angles_rad(:);
    tauH = tauHat(:);
    nuH = nuHat(:);
    thH = thetaHat(:);

    Pt = numel(tauT);
    Ph = numel(tauH);
    assoc = struct('schemaVersion', 'ISAC_ASSOCIATION_V2_TRUTH_POWER_AWARE', ...
                 'domain', char(domain), 'truthCount', Pt, 'estimateCount', Ph, ...
                 'valid',( Pt == Ph && Ph > 0), 'truthToEstimate', [], ...
                 'totalCost', Inf, 'costMatrix', [], 'truthWeights', []);
    if ~assoc.valid
        return;
    end

    dt = (tauT - tauH.') / cfg.waveform.delayBin_s;
    dn = (nuT - nuH.') / cfg.waveform.dopplerBin_Hz;
    da = rad2deg(wrap_pi(thT - thH.')) / 5;
    C = dt.^2 + dn.^2 + da.^2;

    if strcmpi(domain, 'communication')
        if isfield(truth, 'pathPowers') && numel(truth.pathPowers) == Pt
            w = truth.pathPowers(:);
        else
            w = abs(truth.gains(:)).^2;
        end
        w = w / max(sum(w), realmin);
        C = C .* w;
    else
        w = ones(Pt, 1) / Pt;
    end

    [assign, totalCost] = dp_assignment(C);
    assoc.truthToEstimate = assign(:);
    assoc.totalCost = totalCost;
    assoc.costMatrix = C;
    assoc.truthWeights = w;
end

function [assign, bestCost] = dp_assignment(C)
    n = size(C, 1);
    if size(C, 2) ~= n
        error('Assignment cost matrix must be square.');
    end
    if n > 20
        error('DP assignment intended for small sparse path counts (n<=20).');
    end
    nStates = 2^n;
    dp = inf(n + 1, nStates);
    parentMask = zeros(n + 1, nStates, 'uint32');
    parentCol = zeros(n + 1, nStates, 'uint16');
    dp(1, 1) = 0;
    for i = 1:n
        for mask = 0:nStates - 1
            base = dp(i, mask + 1);
            if ~isfinite(base) || sum(bitget(mask, 1:n)) ~= i - 1
                continue;
            end
            for j = 1:n
                if ~bitget(mask, j)
                    newMask = bitset(mask, j, 1);
                    val = base + C(i, j);
                    if val < dp(i + 1, newMask + 1)
                        dp(i + 1, newMask + 1) = val;
                        parentMask(i + 1, newMask + 1) = uint32(mask);
                        parentCol(i + 1, newMask + 1) = uint16(j);
                    end
                end
            end
        end
    end
    fullMask = nStates - 1;
    bestCost = dp(n + 1, fullMask + 1);
    assign = zeros(n, 1);
    mask = fullMask;
    for i = n:-1:1
        j = double(parentCol(i + 1, mask + 1));
        if j == 0,
            error('Assignment backtracking failed.');
        end
        assign(i) = j;
        mask = double(parentMask(i + 1, mask + 1));
    end
end

function x = wrap_pi(x)
    x = mod(x + pi, 2 * pi) - pi;
end
