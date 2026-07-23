function Plot_Comm_channel_calibration_fullScale_config2(selected_scenario)
if nargin < 1 || strlength(string(selected_scenario)) == 0
    selected_scenario = 'UMi';
end

script_dir = fileparts(mfilename('fullpath'));
result_root = fullfile(script_dir, 'results', 'comm_full_scale_config2');
benchmark_file = fullfile(script_dir, 'Phase2Config2Calibration_v28_CMCC.xlsx');

set(0, 'defaultAxesFontName', 'Times');
set(0, 'defaultTextFontName', 'Times');
set(0, 'defaultAxesFontSize', 13);
set(0, 'defaultTextFontSize', 13);
set(groot, 'defaultFigureWindowStyle', 'docked');

% Map each metric to its saved field name and benchmark block index.
metric_configs = {
    localMetricConfig('CouplingLoss', 'Coupling loss (dB)', 'CouplingLoss', 1)
    localMetricConfig('GeometrySIR', 'Geometry SIR (dB)', 'Geometry_SIR', 2)
    % The CMCC workbook block order follows the original Flexible plot usage:
    % figure 3 uses block 3 for the largest singular value and
    % figure 4 uses block 4 for the smallest singular value.
    localMetricConfig('LargestSingularValue', '10log10(largest singular value)', 'Largest_Singular_Value', 3)
    localMetricConfig('SmallestSingularValue', '10log10(smallest singular value)', 'Smallest_Singular_Value', 4)
    localMetricConfig('RatioSingularValue', '10log10(ratio singular value)', 'Ratio_Singular_Value', 5)
};

scenario_config = localPlotScenarioConfig(selected_scenario);
frequency_tags = {'6GHz', '30GHz', '60GHz', '70GHz'};
color_map = [
    0 0.447058823529412 0.741176470588235
    0.850980392156863 0.325490196078431 0.098039215686274
    0.4 0.6 0.12
    0.494117647058824 0.184313725490196 0.556862745098039
];
line_styles = {'-', '--', ':', '-.'};
marker_map = {'o', 's', 'x', '+'};

% Load one MAT file per carrier frequency before starting the plot loop.
results = cell(numel(frequency_tags), 1);
for freq_idx = 1:numel(frequency_tags)
    fname = fullfile(result_root, sprintf('%s_%s.mat', scenario_config.ResultLabel, frequency_tags{freq_idx}));
    if ~isfile(fname)
        error('The result file does not exist: %s', fname);
    end
    loaded = load(fname, 'result');
    results{freq_idx} = loaded.result;
end

% Draw each figure in three layers: gray benchmark runs, benchmark mean, then sim markers.
for metric_idx = 1:numel(metric_configs)
    figure_tag = sprintf('%s_%s', scenario_config.ResultLabel, metric_configs{metric_idx}.FileTag);
    fig = figure( ...
        'InvertHardcopy', 'off', ...
        'Color', [1 1 1], ...
        'Renderer', 'painters', ...
        'WindowStyle', 'docked', ...
        'Name', figure_tag, ...
        'NumberTitle', 'off');
    hold on;
    grid on;

    benchmark_handles = gobjects(numel(frequency_tags), 1);
    sim_handles = gobjects(numel(frequency_tags), 1);
    legend_handles = gobjects(2 * numel(frequency_tags), 1);
    legend_labels = cell(2 * numel(frequency_tags), 1);
    benchmark_blocks = cell(numel(frequency_tags), 1);

    for freq_idx = 1:numel(frequency_tags)
        sheet_name = sprintf('%s-%s', scenario_config.BenchmarkPrefix, frequency_tags{freq_idx});
        benchmark_block = localReadBenchmarkBlock(benchmark_file, sheet_name, metric_configs{metric_idx}.MetricIndex);
        benchmark_blocks{freq_idx} = benchmark_block;

        if ~isempty(benchmark_block)
            gray_handles = plot(benchmark_block.Values, benchmark_block.CDF, ...
                'Color', [0.85 0.85 0.85], ...
                'LineWidth', 1.5);
            set(gray_handles, 'HandleVisibility', 'off');
        end
    end

    for freq_idx = 1:numel(frequency_tags)
        benchmark_block = benchmark_blocks{freq_idx};
        if ~isempty(benchmark_block)
            benchmark_handles(freq_idx) = plot( ...
                benchmark_block.Mean, benchmark_block.CDF, ...
                'LineWidth', 2.0, ...
                'LineStyle', line_styles{freq_idx}, ...
                'Color', color_map(freq_idx, :));
        else
            benchmark_handles(freq_idx) = plot(nan, nan, ...
                'LineWidth', 2.0, ...
                'LineStyle', line_styles{freq_idx}, ...
                'Color', color_map(freq_idx, :));
        end
    end

    for freq_idx = 1:numel(frequency_tags)
        sim_data = results{freq_idx}.(metric_configs{metric_idx}.FieldName);
        sim_cdf = results{freq_idx}.Pr;
        marker_indices = unique(round(linspace(5, numel(sim_cdf), min(10, numel(sim_cdf)))));
        sim_handles(freq_idx) = plot(sim_data, sim_cdf, ...
            'LineStyle', 'none', ...
            'LineWidth', 1.8, ...
            'Color', color_map(freq_idx, :), ...
            'Marker', marker_map{freq_idx}, ...
            'MarkerEdgeColor', color_map(freq_idx, :), ...
            'MarkerFaceColor', [1 1 1], ...
            'MarkerIndices', marker_indices, ...
            'MarkerSize', 7);
    end

    for freq_idx = 1:numel(frequency_tags)
        legend_handles(freq_idx) = benchmark_handles(freq_idx);
        legend_labels{freq_idx} = sprintf('%s, benchmark', localFormatFrequencyLabel(frequency_tags{freq_idx}));
        legend_handles(numel(frequency_tags) + freq_idx) = sim_handles(freq_idx);
        legend_labels{numel(frequency_tags) + freq_idx} = sprintf('%s, sim.', localFormatFrequencyLabel(frequency_tags{freq_idx}));
    end

    xlabel(metric_configs{metric_idx}.AxisLabel);
    ylabel('CDF (%)');
    ylim([0, 100]);
    yticks(0:10:100);
    set(gca, 'GridLineStyle', '--', 'GridColor', [0.75, 0.75, 0.75], 'GridAlpha', 1);
    title(metric_configs{metric_idx}.AxisLabel, 'Interpreter', 'none');
    legend(legend_handles, legend_labels, 'Location', 'southeast', 'FontSize', 12);
    saveas(fig, fullfile(result_root, [figure_tag, '.png']));
end
end

function config = localMetricConfig(field_name, axis_label, file_tag, metric_index)
config = struct();
config.FieldName = field_name;
config.AxisLabel = axis_label;
config.FileTag = file_tag;
config.MetricIndex = metric_index;
end

function scenario_config = localPlotScenarioConfig(selected_scenario)
selected_key = lower(string(selected_scenario));
switch selected_key
    case "uma"
        scenario_config = struct('ResultLabel', 'UMa', 'BenchmarkPrefix', 'UMa');
    case "umi"
        scenario_config = struct('ResultLabel', 'UMi', 'BenchmarkPrefix', 'UMi');
    case "inh"
        scenario_config = struct('ResultLabel', 'InH', 'BenchmarkPrefix', 'InH');
    otherwise
        error('Unsupported selected_scenario: %s', string(selected_scenario));
end
end

function label = localFormatFrequencyLabel(frequency_tag)
label = strrep(frequency_tag, 'GHz', ' GHz');
end

function benchmark_block = localReadBenchmarkBlock(benchmark_file, sheet_name, metric_index)
% Extract the benchmark runs and their mean from the matching workbook block.
benchmark_block = [];

if ~isfile(benchmark_file)
    warning('Benchmark workbook not found: %s', benchmark_file);
    return;
end

try
    workbook_sheets = sheetnames(benchmark_file);
catch
    warning('Unable to read workbook sheets from %s.', benchmark_file);
    return;
end

if ~any(strcmp(workbook_sheets, sheet_name))
    warning('Sheet %s was not found in %s.', sheet_name, benchmark_file);
    return;
end

raw_matrix = readmatrix(benchmark_file, 'Sheet', sheet_name, 'Range', 'B29:EY129');
if isempty(raw_matrix)
    return;
end

block_offset = 31 * (metric_index - 1);
source_cols = (1:19) + block_offset;
mean_col = 30 + block_offset;

if size(raw_matrix, 2) < mean_col
    return;
end

value_matrix = raw_matrix(:, source_cols);
mean_values = raw_matrix(:, mean_col);
valid_rows = ~(isnan(mean_values) | all(isnan(value_matrix), 2));

benchmark_block = struct();
benchmark_block.CDF = (0:size(raw_matrix, 1)-1).';
benchmark_block.CDF = benchmark_block.CDF(valid_rows);
benchmark_block.Values = value_matrix(valid_rows, :);
benchmark_block.Mean = mean_values(valid_rows);
end
