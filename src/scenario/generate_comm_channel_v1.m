function ch = generate_comm_channel_v1(cfg, varargin)
    % GENERATE_COMM_CHANNEL_V1 3GPP TDL-B Rayleigh vehicular SIMO-OTFS truth.
    %
    % Rayleigh TDL-B construction adapted from:
    %
    %   A. Mohebbi, W.-P. Zhu, and M. O. Ahmad,
    %   "Two-Dimensional DPS-BEM Based Channel Estimation for MIMO OTFS Systems,"
    %   IEEE Trans. Veh. Technol., vol. 75, no. 3, pp. 4352-4368, Mar. 2026,
    %   doi: 10.1109/TVT.2025.361045.
    %
    % Mohebbi et al. use the 3GPP TR 38.901 V17.0.0 TDL-B NLOS profile:
    %   - 23 Rayleigh taps;
    %   - RMS delay spread = 634 ns, giving max normalized-TDL delay
    %     4.7834*634 ns = 3.0327 us;
    %   - Jakes Doppler nu_p = nu_max*cos(phi_p);
    %   - vehicle speed = 120 km/h.
    %
    % Source standard:
    %   3GPP TR 38.901 V17.0.0 (Release 17), Sec. 7.7.2, Table 7.7.2-2.
    %
    % Project adaptation
    % ------------------
    % This project retains its existing waveform, carrier, SIMO array, absolute
    % one-way propagation delay, and known-model-order P sweeps. Therefore:
    %   1) the full 23-tap TDL-B power-delay profile is stored below;
    %   2) NumPaths=P selects the P strongest TDL-B components, then sorts them
    %      by delay;
    %   3) TDL-B itself does not prescribe a physical AoA per tap, so this SIMO
    %      adaptation retains the project's continuous clustered-AoA extension.
    %
    % Gains are zero-mean circular complex Gaussian with the selected TDL-B
    % average tap powers; magnitudes are therefore Rayleigh distributed.
    %
    % Usage
    %   ch = generate_comm_channel_v1(cfg)
    %   ch = generate_comm_channel_v1(cfg,'NumPaths',8,'Seed',4)
    %
    % Name-value inputs
    %   NumPaths         : integer 1..23
    %   Seed             : RNG seed or []
    %   Speed_kmh        : default cfg.comm.speed_kmh (120 km/h)
    %   LinkRange_m      : absolute one-way reference range [min max] or scalar
    %   LateralOffset_m  : reference lateral offset [min max] or scalar
    %   RmsDelaySpread_s : TDL delay scale, default 634 ns
    %   NlosAoASpread_deg: project SIMO angular spread around reference AoA
    %
    % Backward-compatibility inputs RiceK_dB/RmsExcessDelay_s/MaxExcessDelay_s
    % are accepted but intentionally ignored on this Rayleigh TDL-B branch.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    validate_cfg(cfg);

    p = inputParser;
    p.addParameter('NumPaths', cfg.comm.numPathsNominal, ...
                   @(x)isnumeric(x) && isscalar(x) && x >= 1 && x <= 23 && mod(x, 1) == 0);
    p.addParameter('Seed', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('Speed_kmh', cfg.comm.speed_kmh, ...
                   @(x)isnumeric(x) && isscalar(x) && x >= 0);
    p.addParameter('LinkRange_m', cfg.comm.linkRange_m, @isnumeric);
    p.addParameter('LateralOffset_m', cfg.comm.lateralOffset_m, @isnumeric);
    p.addParameter('RmsDelaySpread_s', cfg.comm.rmsDelaySpread_s, ...
                   @(x)isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('NlosAoASpread_deg', cfg.comm.nlosAoASpread_deg, ...
                   @(x)isnumeric(x) && isscalar(x) && x >= 0);

    % Accepted only so older experiment wrappers do not fail on parsing.
    p.addParameter('RiceK_dB', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('RmsExcessDelay_s', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));
    p.addParameter('MaxExcessDelay_s', [], @(x)isempty(x) || (isnumeric(x) && isscalar(x)));
    p.parse(varargin{:});
    opt = p.Results;

    if ~isempty(opt.Seed)
        rng(opt.Seed, 'twister');
    end

    P = opt.NumPaths;
    c = cfg.constants.c_mps;
    fc = cfg.waveform.fc_Hz;

    % 3GPP TR 38.901 V17.0.0, Table 7.7.2-2: TDL-B.
    tdlDelayNorm = [ ...
                    0.0000; 0.1072; 0.2155; 0.2095; 0.2870; 0.2986; ...
                    0.3752; 0.5055; 0.3681; 0.3697; 0.5700; 0.5283; ...
                    1.1021; 1.2756; 1.5474; 1.7842; 2.0169; 2.8294; ...
                    3.0219; 3.6187; 4.1067; 4.2790; 4.7834];
    tdlPower_dB = [ ...
                   0.0; -2.2; -4.0; -3.2; -9.8; -1.2; -3.4; -5.2; ...
                   -7.6; -3.0; -8.9; -9.0; -4.8; -5.7; -7.5; -1.9; ...
                   -7.6; -12.2; -9.8; -11.4; -14.9; -9.2; -11.3];

    % Select the P strongest mean-power taps.
    [~, strengthOrder] = sort(tdlPower_dB, 'descend');
    selected = strengthOrder(1:P);
    [~, delayOrder] = sort(tdlDelayNorm(selected), 'ascend');
    selected = selected(delayOrder);

    delayNorm = tdlDelayNorm(selected);
    avgPower = 10.^(tdlPower_dB(selected) / 10);
    avgPower = avgPower / sum(avgPower);

    % Existing project absolute-delay convention: TDL delays are excess
    % delays on top of a geometry-derived one-way reference propagation time.
    R0 = draw_uniform_or_fixed(opt.LinkRange_m, 'LinkRange_m');
    x0 = draw_uniform_or_fixed(opt.LateralOffset_m, 'LateralOffset_m');
    if abs(x0) >= R0
        error('LateralOffset_m must satisfy |offset| < LinkRange_m.');
    end
    y0 = sqrt(max(R0^2 - x0^2, 0));
    theta0 = atan2(x0, y0);
    tau0 = R0 / c;

    excess = delayNorm * opt.RmsDelaySpread_s;
    delays = tau0 + excess;
    if any(delays >= cfg.waveform.maxUnambiguousDelay_s)
        error(['Generated TDL-B communication delay exceeds 1/DeltaF = %.3g s. ', ...
               'Reduce reference range/delay spread or change waveform.'], ...
              cfg.waveform.maxUnambiguousDelay_s);
    end

    % Mohebbi / TDL-B vehicular Doppler construction.
    vmax = opt.Speed_kmh / 3.6;
    nuMax = fc * vmax / c;
    jakesAngles = 2 * pi * rand(P, 1);
    dopplers = nuMax * cos(jakesAngles);

    % TDL-B specifies NLOS Rayleigh fading but no per-tap physical AoA.
    % Retain the project's SIMO angular extension without introducing LoS.
    angles = theta0 + deg2rad(opt.NlosAoASpread_deg) * randn(P, 1);
    maxAbs = deg2rad(cfg.comm.maxAbsAoA_deg);
    angles = min(max(angles, -maxAbs), maxAbs);

    % CN(0,sigma_p^2) complex coefficients -> Rayleigh magnitudes.
    gains = sqrt(avgPower / 2) .* (randn(P, 1) + 1j * randn(P, 1));
    if cfg.comm.normalizeTotalPower
        gains = gains / max(norm(gains), realmin);
    end

    ch = struct();
    ch.P = P;
    ch.delays_sec = delays(:);
    ch.excessDelays_sec = excess(:);
    ch.dopplers_Hz = dopplers(:);
    ch.gains = gains(:);
    ch.pathPowers = abs(gains(:)).^2;
    ch.averagePathPowers = avgPower(:);
    ch.angles_rad = angles(:);
    ch.angles_deg = rad2deg(angles(:));

    ch.hasLoS = false;
    ch.channelProfile = '3GPP_TDL_B_RAYLEIGH_MOHEBBI_ADAPTED';
    ch.selectedTDLTapIndices = selected(:);
    ch.tdlNormalizedDelays = delayNorm(:);
    ch.tdlAveragePowers_dB = tdlPower_dB(selected);
    ch.tdlFullTapCount = 23;
    ch.rmsDelaySpread_s = opt.RmsDelaySpread_s;

    % Compatibility/reference geometry metadata. No path is designated LoS.
    ch.losRange_m = R0;
    ch.losPosition_m = [x0, y0];
    ch.losAngle_rad = theta0;
    ch.referenceRange_m = R0;
    ch.referenceAngle_rad = theta0;

    ch.speed_kmh = opt.Speed_kmh;
    ch.nuMax_Hz = nuMax;
    ch.riceK_dB = -Inf;

    ch.meta = struct();
    ch.meta.status = ['ADAPTED 3GPP TDL-B Rayleigh vehicular channel; ', ...
                      'Mohebbi 2026 / 3GPP TR 38.901 V17.0.0; strongest-P truncation ', ...
                      'and clustered AoA are project adaptations'];
    ch.meta.sourcePaper = 'Mohebbi et al., IEEE TVT, 2026, doi:10.1109/TVT.2025.361045';
    ch.meta.sourceStandard = '3GPP TR 38.901 V17.0.0, Sec. 7.7.2, Table 7.7.2-2';
    ch.meta.delayBins = ch.delays_sec / cfg.waveform.delayBin_s;
    ch.meta.excessDelayBins = ch.excessDelays_sec / cfg.waveform.delayBin_s;
    ch.meta.dopplerBins = ch.dopplers_Hz / cfg.waveform.dopplerBin_Hz;
    ch.meta.jakesAngles_rad = jakesAngles;
end

function validate_cfg(cfg)
    req = {'constants', 'waveform', 'comm'};
    for i = 1:numel(req)
        if ~isfield(cfg, req{i})
            error('cfg.%s is required.', req{i});
        end
    end
end

function v = draw_uniform_or_fixed(x, name)
    if isscalar(x)
        v = x;
    elseif numel(x) == 2 && x(2) >= x(1)
        v = x(1) + (x(2) - x(1)) * rand;
    else
        error('%s must be a scalar or [min max].', name);
    end
end
