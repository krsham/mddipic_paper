function cfg = make_isac_scenario_config_v1()
    % MAKE_ISAC_SCENARIO_CONFIG_V1 Canonical scenario-generation configuration.
    %
    % Waveform, array and channel-generation defaults. The small demo
    % overrides selected values in capsule_demo_config.
    %
    % Model sources and conventions
    % --------------------------
    % EXPLICIT (Wang et al., 2025 numerical setup):
    %   fc = 24.25 GHz, M = 64, N = 16, Delta_f = 156.25 kHz,
    %   nominal Nt = Nr = 16, half-wavelength ULA.
    %
    % EXPLICIT / 3GPP-informed (TR 38.901 Rel-19 automotive sensing, referring
    % to TR 37.885 Highway Option A):
    %   Highway vehicle type: passenger vehicle (Type 2), baseline speed
    %   140 km/h, lane-based motion, no clustered dropping. The standard also
    %   permits a 10 m fixed minimum sensing-target separation option for
    %   ISAC-Automotive (TR 38.901 V19.4.0, Table 7.9.1-2).
    %
    % EXPLICIT / ADAPTED (Mohebbi et al., IEEE TVT 2026 + 3GPP TR 38.901):
    %   communication uses the TDL-B NLOS Rayleigh power-delay profile;
    %   the full source profile has 23 taps, RMS delay spread = 634 ns,
    %   Jakes-type per-path Doppler, and vehicle speed = 120 km/h.
    %
    % ASSUMED / project choices:
    %   - sensing: three straight highway lanes, 4 m lane spacing, nearest
    %     lane center 10 m from the RSU, targets at y=20--120 m, and
    %     equal-magnitude random-phase effective sensing gains;
    %   - communication: retain the project 24.25-GHz waveform and absolute
    %     one-way reference range; for known-order P sweeps select the strongest
    %     P TDL-B taps; retain a 10-deg clustered receive-AoA extension because
    %     the scalar TDL-B table does not assign a physical AoA to each tap.
    %
    % Coordinate convention
    % ---------------------
    % RSU at [0,0] m. ULA axis is x. Broadside is +y. Highway runs along y.
    % The project steering vector uses exp(j*pi*sin(theta)*n), therefore
    % theta = atan2(x,y) is consistent with theta=0 at broadside.
    %
    % Output
    %   cfg : nested struct consumed by generate_sensing_targets_v1 and
    %         generate_comm_channel_v1.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    c = 299792458;

    %% Waveform and array configuration
    cfg.constants.c_mps = c;

    cfg.waveform.fc_Hz       = 24.25e9;
    cfg.waveform.M           = 64;
    cfg.waveform.N           = 16;
    cfg.waveform.deltaF_Hz   = 156.25e3;
    cfg.waveform.fs_Hz       = cfg.waveform.M * cfg.waveform.deltaF_Hz;
    cfg.waveform.bandwidth_Hz = cfg.waveform.M * cfg.waveform.deltaF_Hz;
    cfg.waveform.delayBin_s  = 1 / cfg.waveform.bandwidth_Hz;
    cfg.waveform.dopplerBin_Hz = cfg.waveform.deltaF_Hz / cfg.waveform.N;
    cfg.waveform.maxUnambiguousDelay_s = 1 / cfg.waveform.deltaF_Hz;

    cfg.array.Nt = 16;
    cfg.array.Nr = 16;
    cfg.array.spacing_lambda = 0.5;

    %% Sensing truth generator: 3GPP-informed Highway Option-A style
    cfg.sensing.numTargetsNominal = 3;
    cfg.sensing.profile = 'highway_option_a_informed';

    cfg.sensing.rsuPosition_m = [0, 0];
    cfg.sensing.numLanes = 3;
    cfg.sensing.laneSpacing_m = 4;
    cfg.sensing.nearestLaneCenter_m = 10;
    cfg.sensing.laneCentersX_m = cfg.sensing.nearestLaneCenter_m + ...
        (0:cfg.sensing.numLanes - 1) * cfg.sensing.laneSpacing_m;

    cfg.sensing.longitudinalRange_m = [20, 120];
    cfg.sensing.minTargetSeparation_m = 10;
    cfg.sensing.baselineSpeed_kmh = 140;
    cfg.sensing.motionDirectionSign = -1; % vehicles approach RSU along -y
    cfg.sensing.vehicleType = 2;          % passenger vehicle, metadata only
    cfg.sensing.maxPlacementAttempts = 5000;

    % Main estimator-comparison gain model: controlled, not radar-equation.
    cfg.sensing.gainMode = 'equal_random_phase';
    cfg.sensing.gainMagnitude = 1;

    %% Communication truth generator: adapted 3GPP TDL-B Rayleigh channel
    % Source: Mohebbi et al., IEEE TVT 2026, doi:10.1109/TVT.2025.361045,
    % using 3GPP TR 38.901 V17.0.0 TDL-B for vehicular communication.
    %
    % The standard TDL-B profile contains 23 Rayleigh taps. The existing
    % project known-order sweeps are retained by selecting the strongest P
    % TDL-B components in generate_comm_channel_v1.
    cfg.comm.profile = '3gpp_tdl_b_rayleigh_mohebbi_adapted';
    cfg.comm.numPathsNominal = 5;
    cfg.comm.tdlFullTapCount = 23;
    cfg.comm.riceK_dB = -Inf; % no deterministic LoS component in TDL-B
    cfg.comm.speed_kmh = 120;
    cfg.comm.linkRange_m = [20, 120];
    cfg.comm.lateralOffset_m = [10, 18];

    % Mohebbi simulation: desired TDL-B RMS delay spread = 634 ns.
    % With max normalized delay 4.7834 this gives 3.0327 us excess delay.
    cfg.comm.rmsDelaySpread_s = 634e-9;
    cfg.comm.rmsExcessDelay_s = cfg.comm.rmsDelaySpread_s; % compatibility
    cfg.comm.maxExcessDelay_s = 4.7834 * cfg.comm.rmsDelaySpread_s;

    % TDL-B itself does not define a physical AoA per tap. For this SIMO
    % path-parameter study, retain the existing continuous clustered-AoA
    % extension so all estimators receive the same physical-angle truth.
    cfg.comm.nlosAoASpread_deg = 10;
    cfg.comm.maxAbsAoA_deg = 80;

    % Keep instantaneous Rayleigh amplitude variation. The MC runner sets
    % noise from each realized receive power, so per-scene renormalization
    % is unnecessary and would distort the Rayleigh amplitude law.
    cfg.comm.normalizeTotalPower = false;

    %% Metadata / traceability
    cfg.provenance.waveform = 'EXPLICIT: Wang et al. 2025 numerical configuration';
    cfg.provenance.sensingMobility = ['3GPP-INFORMED: TR 38.901 Rel-19 automotive sensing; ', ...
                                      'TR 37.885 Highway Option A baseline speed/lane motion'];
    cfg.provenance.sensingGeometry = ['ASSUMED: straight 3-lane road, 4 m lane spacing, ', ...
                                      '10 m nearest lane center, 20-120 m longitudinal region'];
    cfg.provenance.comm = ['ADAPTED: Mohebbi et al. 2026 vehicular setup using ', ...
                           '3GPP TR 38.901 V17.0.0 TDL-B Rayleigh PDP, 634 ns RMS delay spread, ', ...
                           '120 km/h, and Jakes per-path Doppler; project waveform retained'];
    cfg.provenance.commSpreads = ['SOURCE: TDL-B normalized delay/power profile and ', ...
                                  '634 ns delay scaling from Mohebbi/TR 38.901; ASSUMED extension: ', ...
                                  '10 deg clustered per-path AoA for SIMO parameter estimation'];
end
