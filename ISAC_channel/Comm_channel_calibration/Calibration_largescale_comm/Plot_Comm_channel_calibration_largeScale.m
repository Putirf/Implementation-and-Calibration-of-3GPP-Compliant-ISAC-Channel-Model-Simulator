function Plot_Comm_channel_calibration_largeScale(selected_scenario)
if nargin < 1 || strlength(string(selected_scenario)) == 0
    selected_scenario = 'UMi';
end

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(fileparts(script_dir));
addpath(project_dir);
result_root = fullfile(script_dir, 'results', 'comm_large_scale');
benchmark_file = fullfile(script_dir, 'Phase1Calibration_v42_CMCC.xlsx');

set(0, 'defaultAxesFontName', 'Times');
set(0, 'defaultTextFontName', 'Times');
set(0, 'defaultAxesFontSize', 13);
set(0, 'defaultTextFontSize', 13);
set(groot, 'defaultFigureWindowStyle', 'docked');

% Map each metric to its saved field name and benchmark workbook columns.
metric_configs = {
    localMetricConfig('CouplingLoss', 'Coupling loss (dB)', 'CouplingLoss', ...
        'Coupling Loss', 'A', 'AE')
    localMetricConfig('GeometrySINR', 'Geometry SINR (dB)', 'Geometry SINR', ...
        'Geometry SINR - with white noise added', 'AG', 'BJ')
    localMetricConfig('GeometrySIR', 'Geometry SIR (dB)', 'Geometry SIR', ...
        'Geometry SIR - without white noise added', 'BL', 'CO')
};

scenario_config = localPlotScenarioConfig(selected_scenario);
frequency_tags = {'6GHz', '30GHz', '70GHz'};
color_map = [
    0 0.447058823529412 0.741176470588235
    0.850980392156863 0.325490196078431 0.098039215686274
    0.4 0.6 0.12
];
line_styles = {'-', '--', '-.'};
marker_map = {'o', 's', '^'};

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
        benchmark_block = localReadBenchmarkBlock(benchmark_file, sheet_name, metric_configs{metric_idx});
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
            benchmark_mean = mean(benchmark_block.Values, 2, 'omitnan');
            valid_benchmark = ~(isnan(benchmark_mean) | isnan(benchmark_block.CDF));
            benchmark_handles(freq_idx) = plot( ...
                benchmark_mean(valid_benchmark), benchmark_block.CDF(valid_benchmark), ...
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

function config = localMetricConfig(field_name, axis_label, file_tag, benchmark_header, block_start_col, block_end_col)
config = struct();
config.FieldName = field_name;
config.AxisLabel = axis_label;
config.FileTag = file_tag;
config.BenchmarkHeader = benchmark_header;
config.BlockStartCol = block_start_col;
config.BlockEndCol = block_end_col;
end

function scenario_config = localPlotScenarioConfig(selected_scenario)
selected_key = lower(string(selected_scenario));
switch selected_key
    case "uma"
        scenario_config = struct('ResultLabel', 'UMa', 'BenchmarkPrefix', 'UMa');
    case "umi"
        scenario_config = struct('ResultLabel', 'UMi', 'BenchmarkPrefix', 'UMi');
    case "inh"
        scenario_config = struct('ResultLabel', 'InH', 'BenchmarkPrefix', 'Indoor');
    otherwise
        error('Unsupported selected_scenario: %s', string(selected_scenario));
end
end

function label = localFormatFrequencyLabel(frequency_tag)
label = strrep(frequency_tag, 'GHz', ' GHz');
end

function benchmark_block = localReadBenchmarkBlock(benchmark_file, sheet_name, metric_config)
% Extract the benchmark runs for the requested metric block and keep only valid rows.
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

header_cell = readcell(benchmark_file, 'Sheet', sheet_name, ...
    'Range', sprintf('%s28:%s28', metric_config.BlockStartCol, metric_config.BlockStartCol));
header_text = "";
if ~isempty(header_cell)
    cell_value = header_cell{1};
    if ischar(cell_value) || isstring(cell_value)
        header_text = string(strtrim(char(cell_value)));
    end
end

if header_text == "" || ~strcmpi(header_text, metric_config.BenchmarkHeader)
    warning('Metric header %s was not found at %s28 in sheet %s.', ...
        metric_config.BenchmarkHeader, metric_config.BlockStartCol, sheet_name);
    return;
end

block_start_col = localColNumber(metric_config.BlockStartCol);
source_row = readcell(benchmark_file, 'Sheet', sheet_name, ...
    'Range', sprintf('%s25:%s25', metric_config.BlockStartCol, metric_config.BlockEndCol));

company_cols = [];
for offset = 1:numel(source_row)
    cell_value = source_row{offset};
    if ischar(cell_value) || isstring(cell_value)
        cell_text = strtrim(char(cell_value));
        if ~isempty(cell_text) && ~strcmpi(cell_text, 'Source') && ~strcmpi(cell_text, 'Mean')
            absolute_col = block_start_col + offset - 1;
            company_cols(end + 1) = absolute_col; %#ok<AGROW>
        end
    end
end

if isempty(company_cols)
    warning('No benchmark source columns were found for %s in sheet %s.', ...
        metric_config.BenchmarkHeader, sheet_name);
    return;
end

cdf_values = readmatrix(benchmark_file, 'Sheet', sheet_name, 'Range', 'A29:A129');
value_matrix = readmatrix(benchmark_file, 'Sheet', sheet_name, ...
    'Range', sprintf('%s29:%s129', localColLetter(company_cols(1)), localColLetter(company_cols(end))));
col_selector = company_cols - company_cols(1) + 1;
value_matrix = value_matrix(:, col_selector);

valid_rows = ~all(isnan(value_matrix), 2) & ~isnan(cdf_values(:));

benchmark_block = struct();
benchmark_block.CDF = cdf_values(valid_rows);
benchmark_block.Values = value_matrix(valid_rows, :);
end

function letters = localColLetter(col)
letters = '';
while col > 0
    remainder = mod(col - 1, 26);
    letters = [char('A' + remainder), letters]; %#ok<AGROW>
    col = floor((col - 1) / 26);
end
end

function col = localColNumber(col_letters)
col = 0;
col_letters = char(col_letters);
for idx = 1:numel(col_letters)
    col = col * 26 + double(col_letters(idx)) - double('A') + 1;
end
end
