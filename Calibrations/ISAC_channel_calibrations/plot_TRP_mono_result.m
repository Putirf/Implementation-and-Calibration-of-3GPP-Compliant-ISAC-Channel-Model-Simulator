function plot_TRP_mono_result(TRP_mono_target_result, calibrationRoot)
%PLOT_TRP_MONO_RESULT Plot multi-frequency TRP calibration CDF figures.
%
%   plot_TRP_mono_result(TRP_mono_target_result)
%   plot_TRP_mono_result(TRP_mono_target_result, calibrationRoot)
%
% Input may be a result array (for example, 6 GHz and 30 GHz).
% Each result is plotted separately, with frequency in titles and filenames.
%   .fc
%   .option
%   .scenario
%   .sen_type
%   .link_type (optional, e.g. 'TRPmo(t)', 'TRPmo(b)', 'TRP-TRP', 'TRP-UE')
%   .LS.CouplingLoss
%   .Full.CouplingLoss
%   .Full.DS
%   .Full.SA

paths = get_ISAC_calibration_paths();
addpath(paths.isacRoot);

if nargin < 2 || isempty(calibrationRoot)
    calibrationRoot = paths.referenceRoot;
end

figureOutputRoot = paths.resultsRoot;
if ~isfolder(figureOutputRoot)
    mkdir(figureOutputRoot);
end

fprintf("===== start plot TRP mono target result =====\n");

results = TRP_mono_target_result;
if isempty(results)
    warning('No calibration result was supplied.');
    return;
end

for ir = 1:numel(results)
    result = results(ir);
    if ~isfield(result, 'fc')
        error('TRP_mono_target_result(%d).fc is required.', ir);
    end
end

[~, frequencyOrder] = sort([results.fc]);
results = results(frequencyOrder);
entries = repmat(struct( ...
    'result', struct(), 'caseInfo', struct(), 'linkType', "", ...
    'sheet', "", 'benchLS', [], 'benchFull', []), 1, numel(results));

for ir = 1:numel(results)
    result = results(ir);

    caseInfo = localResolveCase(result);
    [xlsxLS, xlsxFull] = localGetCalibrationFiles(caseInfo, calibrationRoot);
    linkType = localGetLinkType(result);
    sheet = localGetSheetName(result.fc, linkType);

    fprintf('Plot result %d: scenario=%s, sen_type=%s, link_type=%s, fc=%.0f GHz, sheet=%s\n', ...
        ir, string(result.scenario), string(result.sen_type), linkType, result.fc/1e9, sheet);

    entries(ir).result = result;
    entries(ir).caseInfo = caseInfo;
    entries(ir).linkType = linkType;
    entries(ir).sheet = sheet;
    entries(ir).benchLS = localReadBenchmark(xlsxLS, sheet);
    entries(ir).benchFull = localReadBenchmark(xlsxFull, sheet);
end

for plotIndex = 1:numel(entries)
    entry = entries(plotIndex);

if localHasMetric(entry.result, 'LS', 'CouplingLoss')
    localPlotMetric(entry, 'LS', 'CouplingLoss', figureOutputRoot);
end

if localHasMetric(entry.result, 'Full', 'CouplingLoss')
    localPlotMetric(entry, 'Full', 'CouplingLoss', figureOutputRoot);
end
if localHasMetric(entry.result, 'Full', 'DS')
    localPlotMetric(entry, 'Full', 'DS', figureOutputRoot);
end
if localHasMetric(entry.result, 'Full', 'SA')
    localPlotSAMetrics(entry, figureOutputRoot);
end
end

end

function tf = localHasMetric(result, scaleName, metricName)
tf = isfield(result, scaleName) && isfield(result.(scaleName), metricName);
end

function caseInfo = localResolveCase(result)
scenario = lower(string(result.scenario));
senType = lower(string(result.sen_type));

caseInfo = struct();
caseInfo.option = 2;
if isfield(result, 'option') && ~isempty(result.option)
    caseInfo.option = result.option;
end

if contains(scenario, "uma") && contains(senType, "uav")
    caseInfo.caseCombination = 1;
    caseInfo.commCaseNum = 1;
    caseInfo.senCaseNum = 1;
    caseInfo.folder = "UAV-UMa-AV";
elseif contains(scenario, "uma") && contains(senType, "human")
    caseInfo.caseCombination = 2;
    caseInfo.commCaseNum = 1;
    caseInfo.senCaseNum = 2;
    caseInfo.folder = "Human Outdoor-UMa";
elseif contains(scenario, "umi") && contains(senType, "human")
    caseInfo.caseCombination = 3;
    caseInfo.commCaseNum = 2;
    caseInfo.senCaseNum = 2;
    caseInfo.folder = "Human Outdoor-UMi";
elseif contains(scenario, "inh") && contains(senType, "human")
    caseInfo.caseCombination = 4;
    caseInfo.commCaseNum = 3;
    caseInfo.senCaseNum = 2;
    caseInfo.folder = "Human Indoor-InH";
elseif contains(scenario, "inf") && contains(scenario, "sh") && contains(senType, "human")
    caseInfo.caseCombination = 5;
    caseInfo.commCaseNum = 4;
    caseInfo.senCaseNum = 2;
    caseInfo.folder = "Human Indoor-InF-SH";
elseif contains(scenario, "inf") && contains(scenario, "sh") && contains(senType, "agv")
    caseInfo.caseCombination = 6;
    caseInfo.commCaseNum = 4;
    caseInfo.senCaseNum = 4;
    caseInfo.folder = "AGV-InF-SH";
elseif contains(scenario, "urbangrid") && contains(senType, "vehicle")
    caseInfo.caseCombination = 7;
    caseInfo.commCaseNum = 5;
    caseInfo.senCaseNum = 3;
    caseInfo.folder = "Auto-Urban Grid";
else
    error('Unsupported calibration case: scenario="%s", sen_type="%s".', ...
        string(result.scenario), string(result.sen_type));
end

caseInfo.scenario = string(result.scenario);
caseInfo.senType = string(result.sen_type);
end

function [xlsxLS, xlsxFull] = localGetCalibrationFiles(caseInfo, calibrationRoot)
folder = fullfile(calibrationRoot, char(caseInfo.folder));

switch caseInfo.caseCombination
    case 1
        xlsxLS = fullfile(folder, 'LargeScaleCalibration_UMa-AV_TerrestrialUT_v038_Ericsson_Ericsson2.xlsx');
        if caseInfo.option == 1
            xlsxFull = fullfile(folder, 'FullCalibration_UMa-AV_ConcatenationOption1_TerrestrialUT_v008_mod.xlsx');
        else
            xlsxFull = fullfile(folder, 'FullCalibration_UMa-AV_ConcatenationOption2_TerrestrialUT_v037_Ericsson_mod.xlsx');
        end
    case 2
        xlsxLS = fullfile(folder, 'HumanOutdoor_LargeScaleCalibration_UMa_TerrestrialUT_v015_MTK_BUPT.xlsx');
        xlsxFull = fullfile(folder, 'HumanOutdoor_FullCalibration_UMa_ConcatenationOption2_TerrestrialUT_v010_Apple_mod.xlsx');
    case 3
        xlsxLS = fullfile(folder, 'HumanOutdoorLargeScaleCalibration_UMi_TerrestrialUT_v006_ITRI_BUPT.xlsx');
        if caseInfo.option == 1
            xlsxFull = fullfile(folder, 'HumanOutdoorFullCalibration_UMi_ConcatenationOption1_TerrestrialUT_v003_mod.xlsx');
        else
            xlsxFull = fullfile(folder, 'HumanOutdoorFullCalibration_UMi_ConcatenationOption2_TerrestrialUT_v005_Apple_mod.xlsx');
        end
    case 4
        xlsxLS = fullfile(folder, 'LargeScaleCalibration_InH_TerrestrialUT_v015_Spreadtrum_OPPO.xlsx');
        if caseInfo.option == 1
            xlsxFull = fullfile(folder, 'FullCalibration_InH_ConcatenationOption1_TerrestrialUT_v003_mod.xlsx');
        else
            xlsxFull = fullfile(folder, 'FullCalibration_InH_ConcatenationOption2_TerrestrialUT_v013_Apple_mod.xlsx');
        end
    case 5
        xlsxLS = fullfile(folder, 'Human_LargeScaleCalibration_InF-SH_TerrestrialUT_v009_CATT_ITRI.XLSX');
        xlsxFull = fullfile(folder, 'FullCalibration_InF-SH_ConcatenationOption2_TerrestrialUT_v005_CATT_mod.xlsx');
    case 6
        xlsxLS = fullfile(folder, 'AGVLargeScaleCalibration_InF-SH_TerrestrialUT_OneScatteringPoint_v006_BUPT_IDCC.xlsx');
        xlsxFull = fullfile(folder, 'AGVFullCalibration_InF-SH_ConcatenationOption2_TerrestrialUT_OneScatteringPoint_v006_IDCC_mod.xlsx');
    case 7
        xlsxLS = fullfile(folder, 'AutoLargeScaleCalibration_UrbanGrid_PedestrianUT_OneScatteringPoint_v014_BUPT_mod.xlsx');
        xlsxFull = fullfile(folder, 'AutoFullCalibration_UrbanGrid_ConcatenationOption2_PedestrianUT_OneScatteringPoint_v013_mod.xlsx');
    otherwise
        error('Unsupported case combination: %d.', caseInfo.caseCombination);
end

if ~isfile(xlsxLS)
    error('Large-scale calibration file not found: %s', xlsxLS);
end
if ~isfile(xlsxFull)
    error('Full calibration file not found: %s', xlsxFull);
end
end

function linkType = localGetLinkType(result)
if isfield(result, 'link_type') && ~isempty(result.link_type)
    linkType = string(result.link_type);
else
    linkType = "TRPmo(t)";
end
end

function sheet = localGetSheetName(fc, linkType)
if fc == 30e9
    tag = '30GHz';
else
    tag = '6GHz';
end

sheet = sprintf('%s-%s', linkType, tag);
end

function bench = localReadBenchmark(xlsxFile, sheet)
persistent CACHE
if isempty(CACHE)
    CACHE = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

key = [char(xlsxFile), '||', char(sheet)];
if isKey(CACHE, key)
    bench = CACHE(key);
    return;
end

if ~localHasSheet(xlsxFile, sheet)
    warning('Sheet "%s" not found in "%s". Skip.', sheet, xlsxFile);
    bench = [];
    return;
end

try
    rowLabel = 28;
    rowCompany = 25;
    rowDataStart = 29;
    rowDataEnd = 129;

    headerRow = readcell(xlsxFile, 'Sheet', sheet, ...
        'Range', sprintf('A%d:GF%d', rowLabel, rowLabel));
    headerRow = headerRow(1, :);

    startCol.CouplingLoss = localFindCol(headerRow, 'Coupling loss') + 1;
    startCol.DS = localFindCol(headerRow, 'Delay Spread');
    startCol.ASD = localFindCol(headerRow, 'ASD');
    startCol.ZSD = localFindCol(headerRow, 'ZSD');
    startCol.ASA = localFindCol(headerRow, 'ASA');
    startCol.ZSA = localFindCol(headerRow, 'ZSA');

    bench = struct();
    bench.CouplingLoss = localReadBlock(xlsxFile, sheet, startCol.CouplingLoss, rowCompany, rowDataStart, rowDataEnd);
    bench.DS = localReadBlock(xlsxFile, sheet, startCol.DS, rowCompany, rowDataStart, rowDataEnd);
    bench.ASD = localReadBlock(xlsxFile, sheet, startCol.ASD, rowCompany, rowDataStart, rowDataEnd);
    bench.ZSD = localReadBlock(xlsxFile, sheet, startCol.ZSD, rowCompany, rowDataStart, rowDataEnd);
    bench.ASA = localReadBlock(xlsxFile, sheet, startCol.ASA, rowCompany, rowDataStart, rowDataEnd);
    bench.ZSA = localReadBlock(xlsxFile, sheet, startCol.ZSA, rowCompany, rowDataStart, rowDataEnd);

    CACHE(key) = bench;
catch ME
    warning('Read benchmark failed: file="%s", sheet="%s". (%s)', xlsxFile, sheet, ME.message);
    bench = [];
end
end

function col = localFindCol(headerRow, key)
col = [];
for c = 1:numel(headerRow)
    v = headerRow{c};
    if ischar(v) || isstring(v)
        if contains(string(v), key, 'IgnoreCase', true)
            col = c;
            return;
        end
    end
end
end

function blk = localReadBlock(xlsxFile, sheet, startCol, rowCompany, rowDataStart, rowDataEnd)
if isempty(startCol)
    blk = struct('company', string.empty(1, 0), 'x', [], 'f', []);
    return;
end

cP = 1;
c1 = startCol;
c2 = startCol + 28;

rngCompany = sprintf('%s%d:%s%d', localColLetter(c1), rowCompany, localColLetter(c2), rowCompany);
compRow = readcell(xlsxFile, 'Sheet', sheet, 'Range', rngCompany);
compRow = compRow(1, :);

compStr = strings(1, numel(compRow));
for i = 1:numel(compRow)
    v = compRow{i};
    if ismissing(v)
        compStr(i) = "";
    elseif isstring(v)
        compStr(i) = v;
    elseif ischar(v)
        compStr(i) = string(v);
    else
        compStr(i) = "";
    end
end

compStr = strtrim(compStr);
valid = compStr ~= "" & lower(compStr) ~= "<missing>";

if ~any(valid)
    blk = struct('company', string.empty(1, 0), 'x', [], 'f', []);
    return;
end

company = compStr(valid);

rngP = sprintf('%s%d:%s%d', localColLetter(cP), rowDataStart, localColLetter(cP), rowDataEnd);
p = readmatrix(xlsxFile, 'Sheet', sheet, 'Range', rngP);
f = p(:) / 100;

validIdx = find(valid);
cStart = c1 + validIdx(1) - 1;
cEnd = c1 + validIdx(end) - 1;

rngXall = sprintf('%s%d:%s%d', localColLetter(cStart), rowDataStart, localColLetter(cEnd), rowDataEnd);
xAll = readmatrix(xlsxFile, 'Sheet', sheet, 'Range', rngXall);
x = xAll(:, validIdx - validIdx(1) + 1);

blk = struct();
blk.company = company;
blk.x = x;
blk.f = f;
end

function letters = localColLetter(col)
letters = '';
while col > 0
    r = mod(col - 1, 26);
    letters = [char('A' + r), letters]; %#ok<AGROW>
    col = floor((col - 1) / 26);
end
end

function localPlotMetric(entry, scaleName, metricName, figureOutputRoot)
titleKey = sprintf('%s %s', scaleName, localMetricLabel(metricName));
figName = sprintf('%s %s %s %s %gGHz', entry.caseInfo.scenario, ...
    entry.caseInfo.senType, entry.linkType, titleKey, entry.result.fc/1e9);
fig = figure('Name', figName, 'Color', 'w');
hold on;
grid on;
localPlotFrequencyCurves(entry, scaleName, metricName, []);
xlabel(localMetricAxisLabel(metricName));
ylabel('CDF');
title(figName, 'Interpreter', 'none');
ylim([0 1]);
set(gca, 'FontSize', 12);
localSaveFigure(fig, figName, figureOutputRoot);
end

function localPlotSAMetrics(entry, figureOutputRoot)
saNames = ["ASD", "ZSD", "ASA", "ZSA"];
figName = sprintf('%s %s %s Full SA %gGHz', entry.caseInfo.scenario, ...
    entry.caseInfo.senType, entry.linkType, entry.result.fc/1e9);
fig = figure('Name', figName, 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for isa = 1:numel(saNames)
    metricName = char(saNames(isa));
    nexttile;
    hold on;
    grid on;
    localPlotFrequencyCurves(entry, 'Full', metricName, isa);
    xlabel(localMetricAxisLabel(metricName));
    ylabel('CDF');
    title(metricName);
    ylim([0 1]);
    set(gca, 'FontSize', 12);
end
sgtitle(figName, 'Interpreter', 'none');
localSaveFigure(fig, figName, figureOutputRoot);
end

function localPlotFrequencyCurves(entry, scaleName, metricName, saIndex)
legendHandles = gobjects(0);
legendLabels = strings(0);
bench = localBenchmarkForScale(entry, scaleName);
if ~isempty(bench) && isfield(bench, metricName)
    blk = bench.(metricName);
    if ~isempty(blk.x)
        colors = hsv(numel(blk.company) + 1);
        colors = 0.8 * colors(2:end, :);
        for i = 1:numel(blk.company)
            [benchmarkX, benchmarkF] = localSanitizeBenchmarkCurve( ...
                blk.x(:, i), blk.f, metricName, entry.linkType);
            h = plot(benchmarkX, benchmarkF, 'Color', colors(i, :), 'LineWidth', 1);
            legendHandles(end + 1) = h; %#ok<AGROW>
            legendLabels(end + 1) = blk.company(i); %#ok<AGROW>
        end
    end
end
oursData = localSimulationData(entry.result, scaleName, metricName, saIndex);
oursData = oursData(:);
oursData = oursData(isfinite(oursData));
if ~isempty(oursData)
    xOurs = sort(oursData);
    fOurs = (1:numel(xOurs)).' / numel(xOurs);
    hOurs = plot(xOurs, fOurs, '-', 'Color', [1 0 0], 'LineWidth', 2);
    legendHandles(end + 1) = hOurs;
    legendLabels(end + 1) = "Simulation";
end

if ~isempty(legendHandles)
    legendCount = numel(legendHandles);
    isSAMetric = ~isempty(saIndex);
    if isSAMetric && legendCount > 12
        % A 2-by-2 SA figure has much less vertical room per axes.  Use two
        % compact columns so the legend stays below the subplot title.
        legendFontSize = 6;
        legendTokenSize = [10, 6];
        legendColumns = 2;
    elseif legendCount > 16
        legendFontSize = 8;
        legendTokenSize = [14, 8];
        legendColumns = 1;
    elseif legendCount > 8
        legendFontSize = 9;
        legendTokenSize = [16, 9];
        legendColumns = 1;
    else
        legendFontSize = 10;
        legendTokenSize = [20, 10];
        legendColumns = 1;
    end

    lgd = legend(legendHandles, legendLabels, ...
        'Location', 'southeast', ...
        'Interpreter', 'none', ...
        'FontSize', legendFontSize, ...
        'NumColumns', legendColumns);
    lgd.ItemTokenSize = legendTokenSize;
end
end

function [x, f] = localSanitizeBenchmarkCurve(x, f, metricName, linkType)
% Remove invalid spreadsheet cells without changing simulation samples.
x = x(:);
f = f(:);
valid = isfinite(x) & isfinite(f);
x = x(valid);
f = f(valid);

if strcmp(metricName, 'DS') && strcmp(string(linkType), "UEmo(t)") && ...
        numel(x) >= 2
    % Some supplied calibration workbooks contain several corrupt DS cells
    % at the upper tail (for example 3552.64 ns followed by 2.25709e11 ns).
    % Discard the first point after a four-decade jump and the remaining tail
    % so multiple consecutive corrupt values cannot evade the check.
    previous_x = x(1:end-1);
    next_x = x(2:end);
    positive_pair = previous_x > 0 & next_x > 0;
    jump_ratio = zeros(numel(x) - 1, 1);
    jump_ratio(positive_pair) = next_x(positive_pair) ./ previous_x(positive_pair);
    corrupt_start = find(jump_ratio > 1e4, 1, 'first') + 1;
    if ~isempty(corrupt_start)
        x(corrupt_start:end) = [];
        f(corrupt_start:end) = [];
    end
end
end

function bench = localBenchmarkForScale(entry, scaleName)
if strcmp(scaleName, 'LS')
    bench = entry.benchLS;
else
    bench = entry.benchFull;
end
end

function data = localSimulationData(result, scaleName, metricName, saIndex)
data = [];
if strcmp(metricName, 'ASD') || strcmp(metricName, 'ZSD') || ...
        strcmp(metricName, 'ASA') || strcmp(metricName, 'ZSA')
    if localHasMetric(result, scaleName, 'SA')
        data = squeeze(result.(scaleName).SA(saIndex, :, :));
    end
elseif localHasMetric(result, scaleName, metricName)
    data = result.(scaleName).(metricName);
end
end

function label = localMetricLabel(metricName)
switch metricName
    case 'CouplingLoss'
        label = 'Coupling loss';
    otherwise
        label = metricName;
end
end

function label = localMetricAxisLabel(metricName)
switch metricName
    case 'CouplingLoss'
        label = 'Coupling loss (dB)';
    case 'DS'
        label = 'DS (ns)';
    otherwise
        label = sprintf('%s (degree)', metricName);
end
end

function localSaveFigure(fig, figName, figureOutputRoot)
safeName = regexprep(char(figName), '[<>:"/\\|?*]', '_');
outputFile = fullfile(figureOutputRoot, [safeName, '.png']);
figOutputFile = fullfile(figureOutputRoot, [safeName, '.fig']);
drawnow;
exportgraphics(fig, outputFile, 'Resolution', 300);
savefig(fig, figOutputFile);
fprintf('Saved figures: %s, %s\n', outputFile, figOutputFile);
end

function tf = localHasSheet(xlsxFile, sheet)
try
    sn = sheetnames(xlsxFile);
    tf = any(strcmp(sn, sheet));
catch
    tf = false;
end
end
