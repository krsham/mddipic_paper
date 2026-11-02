function truth = generate_tdl_abc_channel_v1(cfg, model, seed)
    % GENERATE_TDL_ABC_CHANNEL_V1 Five strongest mean-power taps, common random draws.
    % Reuse canonical B geometry/Jakes/AoA/Rayleigh draws; change PDP only.
    % The 634-ns parameter scales the FULL profile; truncated RMS is reported.
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    assert(~cfg.comm.normalizeTotalPower, 'TDLABC:Normalization', ...
           'Canonical instantaneous Rayleigh amplitude variation must be retained.');
    profile = tdl_abc_profile_v1(model);
    truth = generate_comm_channel_v1(cfg, 'NumPaths', 5, 'Seed', seed);
    [~, order] = sort(profile.power_dB, 'descend');
    selected = order(1:5);
    [~, order] = sort(profile.delayNorm(selected), 'ascend');
    selected = selected(order);
    power = 10.^(profile.power_dB(selected) / 10);
    power = power / sum(power);
    if profile.model ~= "B"
        standardGain = truth.gains ./ sqrt(truth.averagePathPowers);
        truth.gains = sqrt(power) .* standardGain;
        truth.excessDelays_sec = profile.delayNorm(selected) * cfg.comm.rmsDelaySpread_s;
        truth.delays_sec = truth.referenceRange_m / cfg.constants.c_mps + truth.excessDelays_sec;
        truth.pathPowers = abs(truth.gains).^2;
        truth.averagePathPowers = power;
        truth.channelProfile = char("3GPP_TDL_" + profile.model + "_P5_PROJECT_ADAPTED");
        truth.selectedTDLTapIndices = selected;
        truth.tdlNormalizedDelays = profile.delayNorm(selected);
        truth.tdlAveragePowers_dB = profile.power_dB(selected);
        truth.tdlFullTapCount = numel(profile.delayNorm);
        truth.meta.status = 'ADAPTED: strongest-five truncation; project Jakes/AoA extension';
        truth.meta.sourcePaper = '';
        truth.meta.sourceStandard = profile.source;
        truth.meta.delayBins = truth.delays_sec / cfg.waveform.delayBin_s;
        truth.meta.excessDelayBins = truth.excessDelays_sec / cfg.waveform.delayBin_s;
    end
    assert(all(truth.delays_sec < cfg.waveform.maxUnambiguousDelay_s), ...
           'TDLABC:Delay', 'Retained delay exceeds waveform unambiguous support.');
    meanDelay = sum(power .* truth.excessDelays_sec);
    truth.retainedPdpRms_s = sqrt(sum(power .* (truth.excessDelays_sec - meanDelay).^2));
    truth.tdlModel = profile.model;
end
