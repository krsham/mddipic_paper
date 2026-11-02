function summary = capsule_summarize_cases(metrics, methods, domains, snrs)
%CAPSULE_SUMMARIZE_CASES Aggregate paired observations without dropping failures.
%   NMSE is averaged in linear units. Timing is the median. Empty or partly
%   failed groups retain NaN aggregates and explicit requested/success counts.
%
%   SPDX-License-Identifier: GPL-3.0-only
%   Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
%   Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
%   Attribution and paper citation: CITATION.md.

    methods = reshape(string(methods), 1, []);
    domains = reshape(string(domains), 1, []);
    snrs = sort(unique(snrs(:)'));
    template = struct('domain', "", 'method', "", 'snr_dB', NaN, ...
                      'requestedTrials', 0, 'successfulTrials', 0, ...
                      'meanResponseNMSE', NaN, 'medianEstimate_s', NaN);
    % Preallocate one record per domain, method and SNR group.
    records = repmat(template, numel(domains) * numel(methods) * numel(snrs), 1);
    index = 0;
    for domain = domains
        for method = methods
            for snrDB = snrs
                group = metrics(metrics.domain == domain & metrics.method == method & ...
                                metrics.snr_dB == snrDB, :);
                good = group.success;
                item = template;
                item.domain = domain;
                item.method = method;
                item.snr_dB = snrDB;
                item.requestedTrials = height(group);
                item.successfulTrials = sum(good);
                if ~isempty(good) && all(good)
                    item.meanResponseNMSE = mean(group.responseNMSE);
                    item.medianEstimate_s = median(group.estimate_s);
                end
                index = index + 1;
                records(index) = item;
            end
        end
    end
    summary = struct2table(records);
end
