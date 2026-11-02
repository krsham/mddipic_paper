function report = run_capsule_demo(varargin)
    % RUN_CAPSULE_DEMO Check the paper library on a small, reproducible experiment.
    %   RUN_CAPSULE_DEMO runs numerical checks, six-method sensing/communication
    %   fits, a stopping comparison, a computational ablation and one BER frame.
    %   Results are written beneath results/capsule_demo. Nothing is uploaded.
    %
    %   Name-value options:
    %     Trials       Number of independent scenes per domain/SNR (default 1).
    %     SNRdB        Pilot SNR vector (default [0 20]).
    %     UseParallel  Parallelize independent observations (default false).
    %     Resume       Reuse matching per-case MAT checkpoints (default true).
    %     ShowFigures  Keep the summary figure open (default true).
    %     OutputDir    Absolute or relative results directory.
    %     ValidateOnly Run MATLAB's static checkcode, then return (default false).
    %
    %   Examples:
    %     report = run_capsule_demo;
    %     report = run_capsule_demo('Trials',2,'UseParallel',true);
    %     report = run_capsule_demo('ShowFigures',false,'OutputDir','/results');
    %
    %   This preset is for checking the package. It does not regenerate the
    %   manuscript's Monte Carlo figures or recalibrate stopping thresholds.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    environment = setup();
    parser = inputParser;
    parser.addParameter('Trials', 1, @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
    parser.addParameter('SNRdB', [0 20], @(x) isnumeric(x) && isreal(x) && isvector(x) && ~isempty(x) && all(isfinite(x)));
    parser.addParameter('UseParallel', false, @(x) islogical(x) && isscalar(x));
    parser.addParameter('Resume', true, @(x) islogical(x) && isscalar(x));
    parser.addParameter('ShowFigures', true, @(x) islogical(x) && isscalar(x));
    parser.addParameter('OutputDir', fullfile(environment.root, 'results', 'capsule_demo'), ...
                        @(x) ischar(x) || (isstring(x) && isscalar(x)));
    parser.addParameter('ValidateOnly', false, @(x) islogical(x) && isscalar(x));
    parser.parse(varargin{:});
    options = parser.Results;
    assert(numel(unique(options.SNRdB)) == numel(options.SNRdB), ...
           'SNRdB values must be unique to avoid duplicate checkpoint jobs.');
    if options.ValidateOnly
        sources = dir(fullfile(environment.root, '**', '*.m'));
        diagnostics = cell(numel(sources), 1);
        for k = 1:numel(sources)
            diagnostics{k} = checkcode(fullfile(sources(k).folder, sources(k).name), '-id');
        end
        report = struct('environment', environment, 'checkcode', {diagnostics});
        fprintf('MATLAB checkcode completed for %d files. Inspect report.checkcode.\n', numel(sources));
        return
    end

    output = char(options.OutputDir);
    if ~isfolder(output)
        mkdir(output);
    end
    % Resolve a relative results path before workers start.
    [ok, attributes] = fileattrib(output);
    assert(ok, 'Could not resolve the output directory.');
    output = attributes.Name;
    diary(fullfile(output, 'console.txt'));
    diaryGuard = onCleanup(@() diary('off')); %#ok<NASGU>
    oldRandom = rng;
    randomGuard = onCleanup(@() rng(oldRandom)); %#ok<NASGU>
    config = capsule_demo_config();
    sourceSignature = capsule_source_signature(environment.root);
    plan = struct('schema', 'MDDIPIC_REVIEW_PLAN_V1', 'sourceSignature', sourceSignature, ...
                  'config', config, 'trials', options.Trials, 'snrDB', options.SNRdB(:)', ...
                  'parallel', options.UseParallel, 'matlabRelease', environment.release, ...
                  'platform', environment.platform, 'matlabVersion', environment.matlabVersion);
    encoded = jsonencode(plan, 'PrettyPrint', true);
    planFile = fullfile(output, 'run_config.json');
    processingSourceSignature = sourceSignature;
    if isfile(planFile) && ~isequaln(jsondecode(fileread(planFile)), jsondecode(encoded))
        error('MDDIPIC:CheckpointMismatch', ...
              'Source/configuration changed. Choose a new OutputDir for this package.');
    end
    write_text(planFile, encoded);
    reporting = struct('caseSourceSignature', sourceSignature, ...
                       'processingSourceSignature', processingSourceSignature);
    write_text(fullfile(output, 'reporting_version.json'), jsonencode(reporting, 'PrettyPrint', true));
    write_text(fullfile(output, 'environment.json'), jsonencode(environment, 'PrettyPrint', true));
    completionFile = fullfile(output, 'completion.json');
    if isfile(completionFile)
        delete(completionFile);
    end

    fprintf('MDDIPIC package check: %d scenes, pilot SNRs %s dB.\n', ...
            options.Trials, mat2str(options.SNRdB));
    fprintf('Waveform %dx%d, %d receive antennas. Results: %s\n', ...
            config.waveform.M, config.waveform.N, config.array.Nr, output);
    summaryChecks = check_capsule_summary();
    writetable(summaryChecks, fullfile(output, 'summary_checks.csv'));
    replacementChecks = check_capsule_replacements();
    writetable(replacementChecks, fullfile(output, 'replacement_checks.csv'));
    numericalChecks = check_capsule_numerics();
    writetable(numericalChecks, fullfile(output, 'numerical_checks.csv'));
    domains = ["sensing", "communication"];
    jobs = struct('domain', {}, 'trial', {}, 'snrDB', {});
    for domain = domains
        for snrDB = options.SNRdB(:)'
            for trial = 1:options.Trials
                jobs(end + 1) = struct('domain', domain, 'trial', trial, 'snrDB', snrDB); %#ok<AGROW>
            end
        end
    end
    cases = cell(numel(jobs), 1);
    if options.UseParallel
        if ~license('test', 'Distrib_Computing_Toolbox')
            error('MDDIPIC:ParallelToolbox', 'UseParallel requires Parallel Computing Toolbox.');
        end
        parfor k = 1:numel(jobs)
            cases{k} = run_job(jobs(k), config, output, options.Resume, sourceSignature);
        end
    else
        for k = 1:numel(jobs)
            fprintf('  %s, trial %d, pilot SNR %g dB\n', ...
                    jobs(k).domain, jobs(k).trial, jobs(k).snrDB);
            cases{k} = run_job(jobs(k), config, output, options.Resume, sourceSignature);
        end
    end
    rows = cellfun(@(item) item.rows, cases, 'UniformOutput', false);
    metrics = struct2table(vertcat(rows{:}));
    writetable(metrics, fullfile(output, 'estimator_metrics.csv'));
    summary = capsule_summarize_cases(metrics, config.demo.methodNames, domains, options.SNRdB);
    writetable(summary, fullfile(output, 'summary.csv'));

    % Use the highest-SNR scene for the focused diagnostics, independently
    % of the order in which SNR values were supplied.
    sensingIndex = find([jobs.domain] == "sensing" & [jobs.trial] == 1 & ...
                        [jobs.snrDB] == max(options.SNRdB), 1);
    commIndex = find([jobs.domain] == "communication" & [jobs.trial] == 1 & ...
                     [jobs.snrDB] == max(options.SNRdB), 1);
    stopping = capsule_stopping_demo(config, cases{sensingIndex});
    ablation = capsule_ablation_demo(config, cases{sensingIndex});
    ber = capsule_equalization_demo(config, cases{commIndex});
    writetable(stopping.histories, fullfile(output, 'stopping_histories.csv'));
    writetable(stopping.decisions, fullfile(output, 'stopping_decisions.csv'));
    writetable(ablation, fullfile(output, 'computational_ablation.csv'));
    writetable(ber, fullfile(output, 'equalization.csv'));
    capsule_plot_summary(summary, output, options.ShowFigures);
    success = all(metrics.success) && all(isfinite(ber.BER));
    completion = struct('complete', success, 'sourceSignature', sourceSignature, ...
                        'processingSourceSignature', processingSourceSignature, ...
                        'caseCount', numel(jobs), 'failedFits', sum(~metrics.success), ...
                        'failedEqualizers', sum(~isfinite(ber.BER)), 'publicationEvidence', false);
    write_text(completionFile, jsonencode(completion, 'PrettyPrint', true));
    report = struct('outputDir', output, 'metrics', metrics, 'summary', summary, ...
                    'stopping', stopping, 'ablation', ablation, 'ber', ber, ...
                    'summaryChecks', summaryChecks, 'numericalChecks', numericalChecks, 'completion', completion);
    save(fullfile(output, 'review_report.mat'), 'report', '-v7');
    disp(summary);
    if ~success
        error('MDDIPIC:IncompleteDemo', ...
              'Some fits/solvers failed. See estimator_metrics.csv and equalization.csv; no successful completion recorded.');
    end
    fprintf('Package check completed. This is demonstration output, not publication evidence.\n');
end

function result = run_job(job, config, output, resume, signature)
    file = fullfile(output, sprintf('case_%s_%03d_%+.17g.mat', job.domain, job.trial, job.snrDB));
    if resume && isfile(file)
        saved = load(file, 'result', 'sourceSignature');
        assert(strcmp(saved.sourceSignature, signature), 'Checkpoint source mismatch.');
        result = saved.result;
        return
    end
    result = capsule_run_case(config, job.domain, job.trial, job.snrDB);
    sourceSignature = signature; %#ok<NASGU>
    temporary = [tempname(output) '.mat'];
    save(temporary, 'result', 'sourceSignature', '-v7');
    movefile(temporary, file, 'f');
end

function write_text(path, contents)
    file = fopen(path, 'w', 'n', 'UTF-8');
    if file < 0
        error('MDDIPIC:WriteFailed', 'Cannot open %s for writing.', path);
    end
    guard = onCleanup(@() fclose(file)); %#ok<NASGU>
    fprintf(file, '%s\n', contents);
end
