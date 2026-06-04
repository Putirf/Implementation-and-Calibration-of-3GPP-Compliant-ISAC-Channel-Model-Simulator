function plot_TRP_mono_result(TRP_mono_target_result, calibrationRoot)
%PLOT_TRP_MONO_RESULT Plot TRP monostatic target calibration CDF figures.
%
%   plot_TRP_mono_result(TRP_mono_target_result)
%   plot_TRP_mono_result(TRP_mono_target_result, calibrationRoot)
%
% Input is expected to contain:
%   .fc
%   .option
%   .scenario
%   .sen_type
%   .link_type (optional, e.g. 'TRPmo(t)', 'TRPmo(b)', 'TRP-TRP', 'TRP-UE')
%   .LS.CouplingLoss
%   .Full.CouplingLoss
%   .Full.DS
%   .Full.SA

if nargin < 2 || isempty(calibrationRoot)
    thisDir = fileparts(mfilename('fullpath'));
    calibrationRoot = fullfile(fileparts(thisDir), 'ISAC calibration');
end

fprintf("===== start plot TRP mono target result =====\n");

for ir = 1:numel(TRP_mono_target_result)
    result = TRP_mono_target_result(ir);
    if ~isfield(result, 'fc')
        error('TRP_mono_target_result(%d).fc is required.', ir);
    end

    caseInfo = localResolveCase(result);
    [xlsxLS, xlsxFull] = localGetCalibrationFiles(caseInfo, calibrationRoot);
    linkType = localGetLinkType(result);
    sheet = localGetSheetName(result.fc, linkType);

    fprintf('Plot result %d: scenario=%s, sen_type=%s, link_type=%s, fc=%.0f GHz, sheet=%s\n', ...
        ir, string(result.scenario), string(result.sen_type), linkType, result.fc/1e9, sheet);

    if isfield(result, 'LS') && isfield(result.LS, 'CouplingLoss')
        benchLS = localReadBenchmark(xlsxLS, sheet);
        if ~isempty(benchLS)
            localPlotMetric(benchLS, result.LS.CouplingLoss, 'CouplingLoss', ...
                result.fc, sheet, 'LS', linkType, caseInfo);
        end
    end

    if isfield(result, 'Full')
        benchFull = localReadBenchmark(xlsxFull, sheet);
        if ~isempty(benchFull)
            if isfield(result.Full, 'CouplingLoss')
                localPlotMetric(benchFull, result.Full.CouplingLoss, 'CouplingLoss', ...
                    result.fc, sheet, 'Full', linkType, caseInfo);
            end

            if isfield(result.Full, 'DS')
                localPlotMetric(benchFull, result.Full.DS, 'DS', ...
                    result.fc, sheet, 'Full', linkType, caseInfo);
            end

            if isfield(result.Full, 'SA')
                saName = ["ASD", "ZSD", "ASA", "ZSA"];
                saData = cell(1, 4);
                for isa = 1:4
                    data = squeeze(result.Full.SA(isa, :, :));
                    saData{isa} = data;
                end
                localPlotSAMetrics(benchFull, saData, saName, ...
                    result.fc, sheet, 'Full', linkType, caseInfo);
            end
        end
    end
end

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

function localPlotMetric(bench, oursData, metricName, fc, sheet, caliType, linkType, caseInfo)
titleKey = sprintf('%s %s', caliType, metricName);
titleKey = strrep(titleKey, '_', ' ');

figName = sprintf('%s %s %s %.0fGHz', caseInfo.scenario, caseInfo.senType, titleKey, fc/1e9);
figure('Name', figName);
hold on;
grid on;

blk = bench.(metricName);
if ~isempty(blk.x)
    colors = lines(numel(blk.company));
    for i = 1:numel(blk.company)
        plot(blk.x(:, i), blk.f, 'Color', colors(i, :), 'LineWidth', 0.5);
    end
end

oursData = oursData(:);
oursData = oursData(~isnan(oursData));
if ~isempty(oursData)
    xOurs = sort(oursData);
    n = numel(xOurs);
    fOurs = (1:n).' / n;
    plot(xOurs, fOurs, 'Color', [0.85 0 0], 'LineWidth', 2);
end

if ~isempty(blk.x)
    legend([blk.company, "ours"], 'Location', 'southeast');
else
    legend("ours", 'Location', 'southeast');
end

xlabel(strrep(metricName, '_', ' '));
ylabel('CDF');
title(sprintf('%s %s %s', titleKey, sheet, linkType));
set(gca, 'FontSize', 12);
end

function localPlotSAMetrics(bench, saData, saName, fc, sheet, caliType, linkType, caseInfo)
figName = sprintf('%s %s %s SA %.0fGHz', caseInfo.scenario, caseInfo.senType, caliType, fc/1e9);
figure('Name', figName);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for isa = 1:4
    metricName = char(saName(isa));
    nexttile;
    hold on;
    grid on;

    blk = bench.(metricName);
    if ~isempty(blk.x)
        colors = lines(numel(blk.company));
        for i = 1:numel(blk.company)
            plot(blk.x(:, i), blk.f, 'Color', colors(i, :), 'LineWidth', 0.5);
        end
    end

    oursData = saData{isa};
    oursData = oursData(:);
    oursData = oursData(~isnan(oursData));
    if ~isempty(oursData)
        xOurs = sort(oursData);
        n = numel(xOurs);
        fOurs = (1:n).' / n;
        plot(xOurs, fOurs, 'Color', [0.85 0 0], 'LineWidth', 2);
    end

    if ~isempty(blk.x)
        legend([blk.company, "ours"], 'Location', 'southeast');
    else
        legend("ours", 'Location', 'southeast');
    end

    xlabel(metricName);
    ylabel('CDF');
    title(sprintf('%s %s %s', metricName, sheet, linkType));
    set(gca, 'FontSize', 12);
end
end

function tf = localHasSheet(xlsxFile, sheet)
try
    sn = sheetnames(xlsxFile);
    tf = any(strcmp(sn, sheet));
catch
    tf = false;
end
end
