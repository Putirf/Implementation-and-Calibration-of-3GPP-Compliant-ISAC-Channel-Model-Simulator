function plot_ISAC_spatial_consistency_result(result_input, calibration_root)
%PLOT_ISAC_SPATIAL_CONSISTENCY_RESULT Plot ISAC target spatial consistency calibration.
%
%   plot_ISAC_spatial_consistency_result(result_folder)
%   plot_ISAC_spatial_consistency_result(result_struct)
%   plot_ISAC_spatial_consistency_result(..., calibration_root)

if nargin < 1 || isempty(result_input)
    result_input = fullfile(pwd, 'results', 'ISAC_spatial_consistency');
end
if nargin < 2 || isempty(calibration_root)
    project_folder = fileparts(mfilename('fullpath'));
    calibration_root = fullfile(project_folder, 'ISAC calibration', ...
        'Spatial consistency-UrbanGridandInO');
end

results = localLoadResults(result_input);
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
save_root = fullfile(pwd, 'results', 'ISAC_spatial_consistency_figures', timestamp);
if ~exist(save_root, 'dir')
    mkdir(save_root);
end
for result_idx = 1:numel(results)
    result = results{result_idx};
    benchmark = localReadBenchmark(result, calibration_root);
    localPlotResult(result, benchmark, save_root);
end
fprintf('Saved figures to %s\n', save_root);
end

function results = localLoadResults(result_input)
if isstruct(result_input)
    results = num2cell(result_input);
    return;
end

if ischar(result_input) || isstring(result_input)
    result_path = char(result_input);
    if isfolder(result_path)
        files = dir(fullfile(result_path, 'ISAC_SC_*.mat'));
        results = {};
        for file_idx = 1:numel(files)
            loaded = load(fullfile(files(file_idx).folder, files(file_idx).name), 'result');
            results{end + 1} = loaded.result; %#ok<AGROW>
        end
        if isempty(results)
            error('No ISAC_SC_*.mat result files found in %s.', result_path);
        end
        return;
    elseif isfile(result_path)
        loaded = load(result_path, 'result');
        results = {loaded.result};
        return;
    end
end

error('result_input must be a result struct, result .mat file, or result folder.');
end

function benchmark = localReadBenchmark(result, calibration_root)
xlsx_file = localBenchmarkFile(result, calibration_root);
if isempty(xlsx_file) || ~isfile(xlsx_file)
    warning('Spatial consistency benchmark file not found for %s.', result.case_name);
    benchmark = struct();
    return;
end

sheet = result.sheet;
if ~localHasSheet(xlsx_file, sheet)
    warning('Sheet "%s" not found in "%s".', sheet, xlsx_file);
    benchmark = struct();
    return;
end
fprintf('Loaded benchmark: %s (sheet=%s)\n', xlsx_file, sheet);

benchmark = struct();
benchmark.delay = localReadBenchmarkMetric(xlsx_file, sheet, 'xcorr - delay');
benchmark.aoa = localReadBenchmarkMetric(xlsx_file, sheet, 'xcorr - AoA');
benchmark.los = localReadBenchmarkMetric(xlsx_file, sheet, 'xcorr - LOS');
benchmark.file = xlsx_file;
benchmark.sheet = sheet;
end

function xlsx_file = localBenchmarkFile(result, calibration_root)
case_name = lower(string(result.case_name));
project_folder = fileparts(mfilename('fullpath'));
if contains(case_name, 'urbangrid') || contains(lower(string(result.scenario)), 'urbangrid')
    candidates = {
        fullfile(project_folder, 'benchmark', 'UrbanGrid.xlsx')
        fullfile(calibration_root, ...
            'SpatialConsistencyCalibration_Urban_grid_PedestrianUT_OneScatteringPoint_v004_CATT_CICTCI_xiaomi.xlsx')
        };
elseif contains(case_name, 'inh') || contains(lower(string(result.scenario)), 'inh')
    candidates = {
        fullfile(project_folder, 'benchmark', 'InH.xlsx')
        fullfile(calibration_root, ...
            'SpatialConsistencyCalibration_InH_TerrestrialUT_v002_xiaomi_IDCC.xlsx')
        };
else
    xlsx_file = '';
    return;
end

xlsx_file = '';
for candidate_idx = 1:numel(candidates)
    if isfile(candidates{candidate_idx})
        xlsx_file = candidates{candidate_idx};
        return;
    end
end
end

function metric = localReadBenchmarkMetric(xlsx_file, sheet, metric_name)
row_source = 25;
row_data_start = 29;
row_data_end = 159;

source = readcell(xlsx_file, 'Sheet', sheet, 'Range', sprintf('A%d:EH%d', row_source, row_source));
source = string(source(1, :));
source(ismissing(source)) = "";

block = localMetricBlock(metric_name);
if isempty(block)
    metric = struct('distance', [], 'company', string.empty(1, 0), 'x', []);
    return;
end

distance = readmatrix(xlsx_file, 'Sheet', sheet, ...
    'Range', sprintf('A%d:A%d', row_data_start, row_data_end));
data = readmatrix(xlsx_file, 'Sheet', sheet, ...
    'Range', sprintf('%s%d:%s%d', localColLetter(block(1)), row_data_start, ...
    localColLetter(block(2)), row_data_end));

company = source(block(1):block(2));
valid = company ~= "" & lower(company) ~= "0" & lower(company) ~= "<missing>";
if size(data, 2) > numel(valid)
    data = data(:, 1:numel(valid));
end

metric = struct();
metric.distance = distance(:);
metric.company = company(valid);
metric.x = data(:, valid);
end

function block = localMetricBlock(metric_name)
switch lower(string(metric_name))
    case "xcorr - delay"
        block = [2, 31];
    case "xcorr - aoa"
        block = [33, 62];
    case "xcorr - los"
        block = [64, 93];
    otherwise
        block = [];
end
end

function localPlotResult(result, benchmark, save_root)
metric_keys = localMetricKeysForResult(result);
metric_labels = localMetricLabels(metric_keys);

figure_name = sprintf('%s %s %.0fGHz Spatial Consistency', ...
    result.case_name, result.procedure, result.fc/1e9);
main_fig = figure('Name', figure_name);
tiledlayout(1, numel(metric_keys), 'TileSpacing', 'compact', 'Padding', 'compact');

for metric_idx = 1:numel(metric_keys)
    nexttile;
    hold on;
    grid on;

    key = metric_keys{metric_idx};
    ours = result.xcorr(metric_idx, :);
    plot(result.distance, ours * 100, 'r-', 'LineWidth', 2);
    legend_entries = "ours equal-drop avg";

    has_benchmark = isfield(benchmark, key) && ~isempty(benchmark.(key).x);
    if has_benchmark
        blk = benchmark.(key);
        colors = lines(max(1, numel(blk.company)));
        for col_idx = 1:numel(blk.company)
            plot(blk.distance, blk.x(:, col_idx), '-', ...
                'Color', colors(col_idx, :), 'LineWidth', 0.75);
        end
        mean_curve = mean(blk.x, 2, 'omitnan');
        plot(blk.distance, mean_curve, 'k-', 'LineWidth', 1.5);
        legend([legend_entries, blk.company, "benchmark mean"], 'Location', 'northeast');
    else
        legend(legend_entries, 'Location', 'northeast');
    end

    xlabel('Distance (m)');
    ylabel('Cross-correlation (%)');
    title(metric_labels(metric_idx));
    ylim([-10, 105]);
    xlim([min(result.distance), max(result.distance)]);
end
localSaveFigure(main_fig, save_root, figure_name);
localPlotPairCount(result, save_root);
localPlotLspRawDebug(result, save_root);
localPlotSspRawDebug(result, save_root);
end

function localPlotPairCount(result, save_root)
if ~isfield(result, 'pair_count_per_bin') || isempty(result.pair_count_per_bin)
    return;
end

figure_name = sprintf('%s %s %.0fGHz Pair Count', ...
    result.case_name, result.procedure, result.fc/1e9);
fig = figure('Name', figure_name);
bar(result.distance, result.pair_count_per_bin, 1);
grid on;
xlabel('Distance (m)');
ylabel('Pair count');
title('pair count per distance group');
xlim([min(result.distance), max(result.distance)]);
localSaveFigure(fig, save_root, figure_name);
end

function localPlotLspRawDebug(result, save_root)
if ~isfield(result, 'lsp_raw_xcorr_pooled') || isempty(result.lsp_raw_xcorr_pooled)
    return;
end

metric_names = string(result.lsp_raw_metric_name);
if contains(lower(string(result.case_name) + " " + string(result.scenario)), "urbangrid")
    keep_metric = ~contains(metric_names, "NLOS");
    metric_names = metric_names(keep_metric);
    lsp_equal = result.lsp_raw_xcorr_equal_drop_average(keep_metric, :);
    lsp_pooled = result.lsp_raw_xcorr_pooled(keep_metric, :);
    lsp_theory = result.lsp_raw_theory(keep_metric, :);
else
    lsp_equal = result.lsp_raw_xcorr_equal_drop_average;
    lsp_pooled = result.lsp_raw_xcorr_pooled;
    lsp_theory = result.lsp_raw_theory;
end

num_metrics = numel(metric_names);
if num_metrics == 0
    return;
end
num_cols = min(4, num_metrics);
num_rows = ceil(num_metrics / num_cols);

figure_name = sprintf('%s %s %.0fGHz LSP Raw Spatial Correlation', ...
    result.case_name, result.procedure, result.fc/1e9);
fig = figure('Name', figure_name);
tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

for metric_idx = 1:num_metrics
    nexttile;
    hold on;
    grid on;

    if isfield(result, 'lsp_raw_xcorr_equal_drop_average')
        plot(result.distance, lsp_equal(metric_idx, :) * 100, ...
            'r-', 'LineWidth', 1.2);
    end
    plot(result.distance, lsp_pooled(metric_idx, :) * 100, ...
        'r--', 'LineWidth', 1.2);
    if isfield(result, 'lsp_raw_theory')
        plot(result.distance, lsp_theory(metric_idx, :) * 100, ...
            'k-', 'LineWidth', 1.2);
    end

    title(metric_names(metric_idx), 'Interpreter', 'none');
    xlabel('Distance (m)');
    ylabel('Correlation (%)');
    ylim([-10, 105]);
    xlim([min(result.distance), max(result.distance)]);
    legend(["equal-drop avg", "pooled", "theory"], 'Location', 'northeast');
end
localSaveFigure(fig, save_root, figure_name);
end

function localPlotSspRawDebug(result, save_root)
if ~isfield(result, 'ssp_raw_xcorr_pooled') || isempty(result.ssp_raw_xcorr_pooled)
    return;
end

metric_names = string(result.ssp_raw_metric_name);
if contains(lower(string(result.case_name) + " " + string(result.scenario)), "urbangrid")
    keep_metric = ~contains(metric_names, "NLOS");
    metric_names = metric_names(keep_metric);
    ssp_equal = result.ssp_raw_xcorr_equal_drop_average(keep_metric, :);
    ssp_pooled = result.ssp_raw_xcorr_pooled(keep_metric, :);
    ssp_theory = result.ssp_raw_theory(keep_metric, :);
else
    ssp_equal = result.ssp_raw_xcorr_equal_drop_average;
    ssp_pooled = result.ssp_raw_xcorr_pooled;
    ssp_theory = result.ssp_raw_theory;
end

num_metrics = numel(metric_names);
if num_metrics == 0
    return;
end
num_cols = min(3, num_metrics);
num_rows = ceil(num_metrics / num_cols);

figure_name = sprintf('%s %s %.0fGHz SSP Raw Spatial Correlation', ...
    result.case_name, result.procedure, result.fc/1e9);
fig = figure('Name', figure_name);
tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

for metric_idx = 1:num_metrics
    nexttile;
    hold on;
    grid on;

    plot(result.distance, ssp_equal(metric_idx, :) * 100, ...
        'r-', 'LineWidth', 1.2);
    plot(result.distance, ssp_pooled(metric_idx, :) * 100, ...
        'r--', 'LineWidth', 1.2);
    plot(result.distance, ssp_theory(metric_idx, :) * 100, ...
        'k-', 'LineWidth', 1.2);

    title(metric_names(metric_idx), 'Interpreter', 'none');
    xlabel('Distance (m)');
    ylabel('Correlation (%)');
    ylim([-10, 105]);
    xlim([min(result.distance), max(result.distance)]);
    legend(["equal-drop avg", "pooled", "theory"], 'Location', 'northeast');
end
localSaveFigure(fig, save_root, figure_name);
end

function localSaveFigure(fig, save_root, figure_name)
if nargin < 2 || isempty(save_root) || ~ishandle(fig)
    return;
end

file_stem = regexprep(figure_name, '[^\w.-]', '_');
savefig(fig, fullfile(save_root, [file_stem, '.fig']));
exportgraphics(fig, fullfile(save_root, [file_stem, '.png']), 'Resolution', 200);
end

function metric_keys = localMetricKeysForResult(result)
scenario = lower(string(result.case_name) + " " + string(result.scenario));
if contains(scenario, "urbangrid")
    metric_keys = {'delay', 'aoa'};
else
    metric_keys = {'delay', 'aoa', 'los'};
end
end

function metric_labels = localMetricLabels(metric_keys)
metric_labels = strings(1, numel(metric_keys));
for metric_idx = 1:numel(metric_keys)
    switch metric_keys{metric_idx}
        case 'delay'
            metric_labels(metric_idx) = "xcorr - delay";
        case 'aoa'
            metric_labels(metric_idx) = "xcorr - AoA";
        case 'los'
            metric_labels(metric_idx) = "xcorr - LOS";
    end
end
end

function letters = localColLetter(col)
letters = '';
while col > 0
    r = mod(col - 1, 26);
    letters = [char('A' + r), letters]; %#ok<AGROW>
    col = floor((col - 1) / 26);
end
end

function tf = localHasSheet(xlsx_file, sheet)
try
    tf = any(strcmp(sheetnames(xlsx_file), sheet));
catch
    tf = false;
end
end
