function checks = check_capsule_summary()
%CHECK_CAPSULE_SUMMARY Check schema, linear averaging and failed/empty groups.
%   Uses fixed synthetic table rows; no channel simulation is performed.
%
%   SPDX-License-Identifier: GPL-3.0-only
%   Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
%   Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
%   Attribution and paper citation: CITATION.md.

    domain = ["sensing"; "sensing"; "sensing"; "sensing"; "communication"];
    method = ["A"; "A"; "B"; "B"; "A"];
    snr_dB = [0; 0; 0; 0; 0];
    success = [true; true; true; false; true];
    responseNMSE = [1; 9; 2; NaN; 0.2];
    estimate_s = [0.1; 0.3; 0.2; NaN; 0.4];
    metrics = table(domain, method, snr_dB, success, responseNMSE, estimate_s);
    result = capsule_summarize_cases(metrics, ["A", "B"], ...
                                    ["sensing", "communication"], [0, 20]);
    assert(height(result) == 8, 'Summary must retain every requested group.');
    good = result(result.domain == "sensing" & result.method == "A" & result.snr_dB == 0, :);
    assert(good.requestedTrials == 2 && good.successfulTrials == 2);
    assert(abs(good.meanResponseNMSE - 5) < 1e-12, 'NMSE must be averaged linearly.');
    assert(abs(good.medianEstimate_s - 0.2) < 1e-12);
    failed = result(result.domain == "sensing" & result.method == "B" & result.snr_dB == 0, :);
    assert(failed.requestedTrials == 2 && failed.successfulTrials == 1);
    assert(isnan(failed.meanResponseNMSE) && isnan(failed.medianEstimate_s));
    empty = result(result.snr_dB == 20, :);
    assert(all(empty.requestedTrials == 0) && all(isnan(empty.meanResponseNMSE)));
    checks = table(["summary schema"; "linear NMSE mean"; "median time"; ...
                    "failed group retained"; "empty group retained"], true(5, 1), ...
                   'VariableNames', {'check', 'passed'});
end
