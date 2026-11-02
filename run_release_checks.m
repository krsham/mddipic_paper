function report = run_release_checks()
    % RUN_RELEASE_CHECKS Check replacement numerical contracts without estimator fits.
    %   Writes small CSV/JSON evidence to results/release_checks. For the complete
    %   sensing/communication demonstration use RUN_CAPSULE_DEMO instead.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    environment = setup();
    output = fullfile(environment.root, 'results', 'release_checks');
    if ~isfolder(output)
        mkdir(output);
    end
    completion = fullfile(output, 'completion.json');
    if isfile(completion)
        delete(completion);
    end
    report.replacements = check_capsule_replacements();
    report.numerics = check_capsule_numerics();
    report.summary = check_capsule_summary();
    writetable(report.replacements, fullfile(output, 'replacement_checks.csv'));
    writetable(report.numerics, fullfile(output, 'numerical_checks.csv'));
    writetable(report.summary, fullfile(output, 'summary_checks.csv'));
    % No personal filesystem path is written into this compact handoff record.
    report.completion = struct('complete', true, ...
                               'sourceSignature', capsule_source_signature(environment.root), ...
                               'release', environment.release, 'platform', environment.platform, ...
                               'matlabVersion', environment.matlabVersion, 'publicationEvidence', false);
    file = fopen(completion, 'w');
    assert(file >= 0, 'MDDIPIC:WriteEvidence', 'Cannot write completion.json.');
    guard = onCleanup(@() fclose(file)); %#ok<NASGU>
    fprintf(file, '%s\n', jsonencode(report.completion, 'PrettyPrint', true));
    fprintf('Release checks passed. Evidence: %s\n', output);
end
