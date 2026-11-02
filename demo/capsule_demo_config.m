function config = capsule_demo_config()
    % CAPSULE_DEMO_CONFIG Small waveform and bounded searches for the small demo.
    %   The physical bandwidth and carrier follow the paper. M, antenna count,
    %   path count in sensing, and refinement budgets are reduced for inspection.
    %   This preset demonstrates execution; its curves are not manuscript results.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    config = make_isac_scenario_config_v1();
    config.waveform.M = 32;
    config.waveform.N = 16;
    config.waveform.fs_Hz = 10e6;
    config.waveform.bandwidth_Hz = 10e6;
    config.waveform.deltaF_Hz = config.waveform.fs_Hz / config.waveform.M;
    config.waveform.delayBin_s = 1 / config.waveform.fs_Hz;
    config.waveform.dopplerBin_Hz = config.waveform.fs_Hz / ...
        (config.waveform.M * config.waveform.N);
    config.waveform.maxUnambiguousDelay_s = config.waveform.M / config.waveform.fs_Hz;
    config.array.Nr = 8;
    config.array.Nt = 1;

    config.demo.sensingPaths = 3;
    config.demo.spacing_m = 6;
    config.demo.tauGrid = (0:20) * config.waveform.delayBin_s;
    config.demo.nuGrid = (-1:1) * config.waveform.dopplerBin_Hz;
    config.demo.angleGrid = deg2rad(-60:10:60);
    config.demo.fineAngleGrid = deg2rad(-60:1:60);
    config.demo.refinementLevels = 3;
    config.demo.shrinkFactor = 4;
    config.demo.maximumOrder = 7;
    config.demo.peakDrop_dB = 6;
    config.demo.residualScale = 0.1;
    config.demo.methodNames = ["MDDIPIC", "DRAP", "DDAFR", "MDDAFR", "OMP", "NOMP"];
end
