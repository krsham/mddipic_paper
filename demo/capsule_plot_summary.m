function capsule_plot_summary(summary, outputFolder, showFigures)
    % CAPSULE_PLOT_SUMMARY Save a compact reconstruction and estimation-time overview.
    %   NMSE is averaged in the linear domain before conversion to dB. Timings
    %   are medians. No smoothing or extrapolation is applied to the demo points.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    visibility = 'off';
    if showFigures
        visibility = 'on';
    end
    figureHandle = figure('Color', 'white', 'Visible', visibility);
    guard = onCleanup(@() close_if_hidden(figureHandle, showFigures)); %#ok<NASGU>
    layout = tiledlayout(figureHandle, 2, 2, 'TileSpacing', 'compact');
    domains = ["sensing", "communication"];
    colors = lines(6);
    for d = 1:numel(domains)
        for panel = 1:2
            axesHandle = nexttile(layout, 2 * (d - 1) + panel);
            hold(axesHandle, 'on');
            methods = unique(summary.method, 'stable');
            for m = 1:numel(methods)
                rows = summary(summary.domain == domains(d) & summary.method == methods(m), :);
                if panel == 1
                    values = 10 * log10(max(rows.meanResponseNMSE, realmin));
                    ylabel(axesHandle, 'Response NMSE (dB)');
                else
                    values = 1000 * rows.medianEstimate_s;
                    ylabel(axesHandle, 'Estimation time (ms)');
                end
                plot(axesHandle, rows.snr_dB, values, '-o', 'Color', colors(m, :), ...
                     'DisplayName', char(methods(m)));
            end
            xlabel(axesHandle, 'Pilot SNR (dB)');
            title(axesHandle, domains(d));
            grid(axesHandle, 'on');
            if d == 1 && panel == 1
                legend(axesHandle, 'Location', 'best');
            end
        end
    end
    title(layout, 'Small-scale package check - not publication curves');
    exportgraphics(figureHandle, fullfile(outputFolder, 'demo_summary.png'), 'Resolution', 180);
    exportgraphics(figureHandle, fullfile(outputFolder, 'demo_summary.pdf'), 'ContentType', 'vector');
end

function close_if_hidden(handle, keepOpen)
    if ~keepOpen && isgraphics(handle)
        close(handle);
    end
end
