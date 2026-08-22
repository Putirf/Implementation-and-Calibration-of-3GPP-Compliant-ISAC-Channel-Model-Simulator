clc;
clearvars;

%% ISAC spatial consistency calibration controls
selected_case_ids = [1 2];
frequency_presets = {'ISAC_FR1', 'ISAC_FR2'};
num_ue_per_drop = 120;
num_drops = 500;
full_channel = true;
plot_layout = true;
result_root = fullfile(pwd, 'results', 'ISAC_spatial_consistency');

rng(1);

calibration_cases = localCalibrationCases();
if ~exist(result_root, 'dir')
    mkdir(result_root);
end

if isempty(selected_case_ids)
    case_list = 1:numel(calibration_cases);
else
    case_list = selected_case_ids;
end

results_to_plot = [];
for case_idx = case_list
    for preset_idx = 1:numel(frequency_presets)
        result = localRunCase( ...
            calibration_cases(case_idx), frequency_presets{preset_idx}, ...
            num_ue_per_drop, num_drops, full_channel, plot_layout);

        file_name = sprintf('ISAC_SC_%s_%s.mat', ...
            result.case_name, result.frequency_preset);
        file_name = regexprep(file_name, '[^\w.-]', '_');
        save(fullfile(result_root, file_name), 'result');
        if isempty(results_to_plot)
            results_to_plot = result;
        else
            results_to_plot(end + 1) = result; %#ok<SAGROW>
        end
        fprintf('Saved %s\n', fullfile(result_root, file_name));
    end
end

plot_ISAC_spatial_consistency_result(results_to_plot);

function calibration_cases = localCalibrationCases()
calibration_cases = struct('name', {}, 'scenarioFcn', {}, 'sensingFcn', {});

calibration_cases(end + 1) = localCase( ...
    'UrbanGrid', ...
    @() comm_scenario.UrbanGrid(500, 2), ...
    @() sensing_types.Vehicle());

calibration_cases(end + 1) = localCase( ...
    'InH', ...
    @() comm_scenario.InH, ...
    @() sensing_types.Human());
end

function case_config = localCase(name, scenario_fcn, sensing_fcn)
case_config = struct();
case_config.name = name;
case_config.scenarioFcn = scenario_fcn;
case_config.sensingFcn = sensing_fcn;
end

function result = localRunCase(case_config, frequency_preset, num_ue, num_drops, full_channel, plot_layout)
fprintf('Run ISAC spatial consistency: case=%s, preset=%s\n', ...
    case_config.name, string(frequency_preset));

scenario = case_config.scenarioFcn();
scenario.applyIsacFrequencyPreset(frequency_preset);
scenario.spatial_consistency_enable = true;
scenario.spatial_consistency_procedure = 'B';

if strcmp(scenario.name, 'UrbanGrid')
    scenario = localConfigureUrbanGridCalibration(scenario);
end

sensing_type = case_config.sensingFcn();
[BS_list, ~] = network_layout.Drop_BaseStation(scenario, plot_layout);

max_group = localMaxDistanceGroup(scenario);
distance_axis = 0:max_group;
metric_keys = localMetricKeys(scenario);
xcorr_per_drop = nan(numel(metric_keys), numel(distance_axis), num_drops);
pair_count_per_drop = zeros(num_drops, numel(distance_axis));
pooled_stats = localEmptyPooledStats(numel(metric_keys), numel(distance_axis));
[lsp_metric_keys, lsp_metric_labels, lsp_theory] = localLspDebugMetricConfig(scenario, distance_axis);
lsp_xcorr_per_drop = nan(numel(lsp_metric_keys), numel(distance_axis), num_drops);
lsp_pooled_stats = localEmptyPooledStats(numel(lsp_metric_keys), numel(distance_axis));
[ssp_metric_keys, ssp_metric_labels, ssp_theory] = localSspDebugMetricConfig(scenario, distance_axis);
ssp_xcorr_per_drop = nan(numel(ssp_metric_keys), numel(distance_axis), num_drops);
ssp_pooled_stats = localEmptyPooledStats(numel(ssp_metric_keys), numel(distance_axis));

for drop_idx = 1:num_drops
    fprintf('  Drop %d/%d\n', drop_idx, num_drops);
    drop_plot_layout = plot_layout && drop_idx == 1;
    drop_data = localSingleDrop(BS_list, scenario, sensing_type, num_ue, max_group, full_channel, drop_plot_layout);
    [xcorr_per_drop(:, :, drop_idx), pair_count_per_drop(drop_idx, :), drop_pooled_stats] = localPairDistanceCorrelation( ...
        drop_data, distance_axis, metric_keys, strcmp(scenario.name, 'UrbanGrid'));
    pooled_stats = localMergePooledStats(pooled_stats, drop_pooled_stats);

    lsp_drop_data = drop_data;
    lsp_drop_data.metrics = drop_data.lsp_raw_metrics;
    [lsp_xcorr_per_drop(:, :, drop_idx), ~, lsp_drop_pooled_stats] = localPairDistanceCorrelation( ...
        lsp_drop_data, distance_axis, lsp_metric_keys, false);
    lsp_pooled_stats = localMergePooledStats(lsp_pooled_stats, lsp_drop_pooled_stats);

    ssp_drop_data = drop_data;
    ssp_drop_data.metrics = drop_data.ssp_raw_metrics;
    [ssp_xcorr_per_drop(:, :, drop_idx), ~, ssp_drop_pooled_stats] = localPairDistanceCorrelation( ...
        ssp_drop_data, distance_axis, ssp_metric_keys, false);
    ssp_pooled_stats = localMergePooledStats(ssp_pooled_stats, ssp_drop_pooled_stats);
end

xcorr_metrics = mean(xcorr_per_drop, 3, 'omitnan');
xcorr_pooled = localPooledXcorr(pooled_stats);
pair_count_per_bin = sum(pair_count_per_drop, 1);
localPrintXcorrComparison(metric_keys, distance_axis, xcorr_metrics, xcorr_pooled, pair_count_per_bin);
lsp_xcorr_equal_avg = mean(lsp_xcorr_per_drop, 3, 'omitnan');
lsp_xcorr_pooled = localPooledXcorr(lsp_pooled_stats);
localPrintLspDebugSummary(lsp_metric_labels, lsp_xcorr_pooled, lsp_theory);
ssp_xcorr_equal_avg = mean(ssp_xcorr_per_drop, 3, 'omitnan');
ssp_xcorr_pooled = localPooledXcorr(ssp_pooled_stats);
localPrintLspDebugSummary(ssp_metric_labels, ssp_xcorr_pooled, ssp_theory);

result = struct();
result.case_name = case_config.name;
result.scenario = [scenario.name, scenario.subname];
result.frequency_preset = char(frequency_preset);
result.fc = scenario.frequency;
result.sheet = localSheetName(scenario.frequency);
result.link_type = 'TRP-UE';
result.procedure = scenario.spatial_consistency_procedure;
result.distance = distance_axis;
result.metric_name = localMetricNames(metric_keys);
result.xcorr = xcorr_metrics;
result.xcorr_equal_drop_average = xcorr_metrics;
result.xcorr_per_drop = xcorr_per_drop;
result.pair_count_per_bin = pair_count_per_bin;
result.pair_count_per_drop = pair_count_per_drop;
result.lsp_raw_metric_name = lsp_metric_labels;
result.lsp_raw_xcorr_equal_drop_average = lsp_xcorr_equal_avg;
result.lsp_raw_xcorr_pooled = lsp_xcorr_pooled;
result.lsp_raw_xcorr_per_drop = lsp_xcorr_per_drop;
result.lsp_raw_theory = lsp_theory;
result.lsp_raw_pair_count = lsp_pooled_stats.count;
result.ssp_raw_metric_name = ssp_metric_labels;
result.ssp_raw_xcorr_equal_drop_average = ssp_xcorr_equal_avg;
result.ssp_raw_xcorr_pooled = ssp_xcorr_pooled;
result.ssp_raw_xcorr_per_drop = ssp_xcorr_per_drop;
result.ssp_raw_theory = ssp_theory;
result.ssp_raw_pair_count = ssp_pooled_stats.count;
result.num_ue_per_drop = num_ue;
result.num_drops = num_drops;
result.max_distance_group_m = max_group;
los_metric_idx = find(strcmp(metric_keys, 'los'), 1);
result.los_state_xcorr_equal_drop_average = [];
result.los_state_xcorr_pooled = [];
result.los_state_pair_count = [];
if ~isempty(los_metric_idx)
    result.los_state_xcorr_equal_drop_average = xcorr_metrics(los_metric_idx, :);
    result.los_state_xcorr_pooled = xcorr_pooled(los_metric_idx, :);
    result.los_state_pair_count = pair_count_per_bin;
end
result.note = ['Multi-drop target-level spatial consistency calibration with equal-weight averaging across drops. ', ...
    'ST and UT positions follow TR 38.901 Table 7.9.6.3-1 distribution. ', ...
    'Distance groups are evaluated up to twice the maximum target-level correlation distance.'];
end

function metric_keys = localMetricKeys(scenario)
if strcmp(scenario.name, 'InH')
    metric_keys = {'delay', 'aoa', 'los'};
else
    metric_keys = {'delay', 'aoa'};
end
end

function metric_names = localMetricNames(metric_keys)
metric_names = cell(size(metric_keys));
for idx = 1:numel(metric_keys)
    switch metric_keys{idx}
        case 'delay'
            metric_names{idx} = 'delay';
        case 'aoa'
            metric_names{idx} = 'AoA';
        case 'los'
            metric_names{idx} = 'LOS';
        otherwise
            error('Unsupported metric key: %s.', metric_keys{idx});
    end
end
end

function scenario = localConfigureUrbanGridCalibration(scenario)
scenario.ISD = 500;
scenario.BS_layer_num = 2;
scenario.grid_layer_num = 3;
scenario.grid_dx = 250;
scenario.grid_dy = 433;
scenario.x_range = [-scenario.grid_dx * (2 * scenario.grid_layer_num - 1) / 2, ...
    scenario.grid_dx * (2 * scenario.grid_layer_num - 1) / 2];
scenario.y_range = [-scenario.grid_dy * (2 * scenario.grid_layer_num - 1) / 2, ...
    scenario.grid_dy * (2 * scenario.grid_layer_num - 1) / 2];
end

function max_group = localMaxDistanceGroup(scenario)
if strcmp(scenario.name, 'UrbanGrid')
    max_group = 80;    % LOS target-UT correlation distance 40 m, NOTE 2 => 2*40 m.
elseif strcmp(scenario.name, 'InH')
    max_group = 20;    % Indoor office target-UT correlation distance 10 m.
else
    max_group = 100;
end
end

function drop_data = localSingleDrop(BS_list, scenario, sensing_type, num_ue, max_group, full_channel, plot_layout)
if plot_layout
    figure(1);
    hold on;
end
[ST_list, st_position, sc_info] = network_layout.drop_ST_sc_calibration(BS_list, scenario, sensing_type, plot_layout, max_group);
[UE_list, ue_positions] = network_layout.drop_UE_ISAC_sc_calibration(BS_list, scenario, ST_list, sc_info, plot_layout, num_ue, max_group);
tx = localNearestEquipment(BS_list, st_position);
if scenario.spatial_consistency_enable
    localApplyTargetLosStateSpatialConsistency(tx, UE_list, ST_list, scenario);
    localApplyTargetLspSpatialConsistency(tx, UE_list, ST_list, scenario);
    localApplyTargetProcedureBSpatialConsistency(UE_list, ST_list, scenario);
end
lsp_raw_metrics = localCollectTargetUeLspRawMetrics(UE_list, ST_list, scenario);
ssp_raw_metrics = localCollectTargetUeProcedureBRawMetrics(UE_list, ST_list, scenario);
metrics = localEmptyMetrics(num_ue);
st_id = ST_list(1).ID;

for ue_idx = 1:num_ue
    rx = UE_list(ue_idx);
    link = channel.Target_channel(tx, rx, ST_list, scenario, full_channel, 0);
    metrics = localStoreMetrics(metrics, ue_idx, link, rx, st_id, scenario);
end

drop_data = struct();
drop_data.st_position = st_position;
drop_data.ue_positions = ue_positions;
drop_data.metrics = metrics;
drop_data.lsp_raw_metrics = lsp_raw_metrics;
drop_data.ssp_raw_metrics = ssp_raw_metrics;
end

function equipment = localNearestEquipment(equipment_list, position)
pos = reshape([equipment_list.Position], 3, []).';
[~, idx] = min(sum((pos(:, 1:2) - position(1:2)).^2, 2));
equipment = equipment_list(idx);
end

function localApplyTargetLosStateSpatialConsistency(tx, UE_list, ST_list, scenario)
% Generate spatially correlated LOS/NLOS state random variables in ST.rand_LoS.
if ~ismember(scenario.name, {'UrbanGrid', 'InH'})
    error('LOS/NLOS state spatial consistency calibration is implemented only for UrbanGrid and InH.');
end

d_cor = localLosStateCorrelationDistance(scenario);
for st_idx = 1:numel(ST_list)
    st = ST_list(st_idx);

    if strcmp(scenario.name, 'UrbanGrid')
        st = localAssignRandLos(st, tx.ID, 0);
        for ue_idx = 1:numel(UE_list)
            st = localAssignRandLos(st, UE_list(ue_idx).ID, 0);
        end
        continue;
    end

    % STX-SPST is a separate link set under TR 38.901 7.9.5. In this
    % calibration it has one link, so the correlated process reduces to one
    % uniform sample.
    st = localAssignRandLos(st, tx.ID, rand);

    if isempty(UE_list)
        continue;
    end

    % Table 7.9.6.3-1 requires the LOS/NLOS-status correlation metric for
    % indoor office. Note 3 excludes the NLOSv blockage state; it does not
    % make the ordinary binary LOS/NLOS state constant. Preserve spatially
    % correlated LOS/NLOS variation so that this metric is well-defined.
    ue_positions = reshape([UE_list.Position], 3, []).';
    corr_uniform = localCorrelatedUniformRandom(ue_positions(:, 1:2), d_cor);
    for ue_idx = 1:numel(UE_list)
        st = localAssignRandLos(st, UE_list(ue_idx).ID, corr_uniform(ue_idx));
    end
end
end

function d_cor = localLosStateCorrelationDistance(scenario)
if strcmp(scenario.name, 'UrbanGrid')
    d_cor = 50;   % Table 7.6.3.1-2 UMa LOS/NLOS state correlation distance.
elseif strcmp(scenario.name, 'InH')
    d_cor = 10;   % Table 7.6.3.1-2 Indoor LOS/NLOS state correlation distance.
else
    error('Unsupported scenario for LOS/NLOS state spatial consistency: %s.', scenario.name);
end
end

function uniform_values = localCorrelatedUniformRandom(positions_2d, d_cor)
num_pos = size(positions_2d, 1);
if num_pos == 0
    uniform_values = [];
    return;
end

corr_matrix = localSpatialCorrelationMatrix(positions_2d, d_cor);
sqrt_corr = chol(corr_matrix, 'lower');
gaussian_values = sqrt_corr * randn(num_pos, 1);
uniform_values = 0.5 * (1 + erf(gaussian_values / sqrt(2)));
uniform_values = min(max(uniform_values, eps), 1 - eps);
end

function st = localAssignRandLos(st, object_id, value)
if isempty(st.rand_LoS)
    st.rand_LoS = nan(1, object_id);
elseif numel(st.rand_LoS) < object_id
    st.rand_LoS(end+1:object_id) = nan;
end
st.rand_LoS(object_id) = value;
end

function localApplyTargetLspSpatialConsistency(tx, UE_list, ST_list, scenario)
% Generate spatially correlated LSP Gaussian variables before channel construction.
if ~ismember(scenario.name, {'UrbanGrid', 'InH'})
    error('Spatial consistency LSP calibration is implemented only for UrbanGrid and InH.');
end

for st_idx = 1:numel(ST_list)
    st = ST_list(st_idx);
    tx_state = localTargetLinkState(tx, st, scenario);
    tx_lsp = localCorrelatedLspRaw(st.Position, tx_state, scenario);
    localAssignLspRaw(tx, tx_state, st.ID, tx_lsp);

    rx_states = strings(numel(UE_list), 1);
    rx_positions = nan(numel(UE_list), 3);
    for ue_idx = 1:numel(UE_list)
        ue = UE_list(ue_idx);
        rx_states(ue_idx) = localTargetLinkState(ue, st, scenario);
        rx_positions(ue_idx, :) = ue.Position;
    end

    for state = ["LOS", "NLOS", "O2I"]
        group_idx = find(rx_states == state);
        if isempty(group_idx)
            continue;
        end

        group_lsp = localCorrelatedLspRaw(rx_positions(group_idx, :), state, scenario);
        for local_idx = 1:numel(group_idx)
            ue = UE_list(group_idx(local_idx));
            localAssignLspRaw(ue, state, st.ID, group_lsp(:, local_idx));
        end
    end
end
end

function localApplyTargetProcedureBSpatialConsistency(UE_list, ST_list, scenario)
% Generate spatially correlated Procedure B uniforms for target-UE links.
if ~ismember(scenario.name, {'UrbanGrid', 'InH'})
    error('Spatial consistency Procedure B random fields are implemented only for UrbanGrid and InH.');
end

if isempty(UE_list) || isempty(ST_list)
    return;
end

for st_idx = 1:numel(ST_list)
    st = ST_list(st_idx);
    rx_states = strings(numel(UE_list), 1);
    rx_positions = nan(numel(UE_list), 3);
    for ue_idx = 1:numel(UE_list)
        ue = UE_list(ue_idx);
        rx_states(ue_idx) = localTargetLinkState(ue, st, scenario);
        rx_positions(ue_idx, :) = ue.Position;
    end

    for state = ["LOS", "NLOS"]
        group_idx = find(rx_states == state);
        if isempty(group_idx)
            continue;
        end

        procB_raw = localCorrelatedProcedureBRaw(rx_positions(group_idx, :), state, scenario);
        for local_idx = 1:numel(group_idx)
            ue = UE_list(group_idx(local_idx));
            localAssignProcedureBRaw(ue, st.ID, procB_raw(local_idx));
        end
    end
end
end

function procB_raw = localCorrelatedProcedureBRaw(positions, state, scenario)
max_clusters = 20;
uniform_field_names = {'tau', 'AOA', 'AOD', 'ZOA', 'ZOD'};
procB_raw = repmat(localEmptyProcedureBRaw(max_clusters), size(positions, 1), 1);
if isempty(positions)
    return;
end

for field_idx = 1:numel(uniform_field_names)
    field_name = uniform_field_names{field_idx};
    d_cor = localProcedureBCorrelationDistance(scenario, state, field_name);
    corr_uniform = localCorrelatedUniformMatrix(positions(:, 1:2), d_cor, max_clusters);
    for pos_idx = 1:size(positions, 1)
        procB_raw(pos_idx).(field_name) = corr_uniform(pos_idx, :);
    end
end

d_cor_shadow = localProcedureBCorrelationDistance(scenario, state, 'shadow');
corr_shadow = localCorrelatedGaussianMatrix(positions(:, 1:2), d_cor_shadow, max_clusters);
for pos_idx = 1:size(positions, 1)
    procB_raw(pos_idx).shadow = corr_shadow(pos_idx, :);
end
end

function raw = localEmptyProcedureBRaw(num_clusters)
raw = struct();
raw.tau = nan(1, num_clusters);
raw.AOA = nan(1, num_clusters);
raw.AOD = nan(1, num_clusters);
raw.ZOA = nan(1, num_clusters);
raw.ZOD = nan(1, num_clusters);
raw.shadow = nan(1, num_clusters);
end

function uniform_values = localCorrelatedUniformMatrix(positions_2d, d_cor, num_fields)
num_pos = size(positions_2d, 1);
if num_pos == 0
    uniform_values = [];
    return;
end

corr_matrix = localSpatialCorrelationMatrix(positions_2d, d_cor);
sqrt_corr = chol(corr_matrix, 'lower');
gaussian_values = sqrt_corr * randn(num_pos, num_fields);
uniform_values = 0.5 * (1 + erf(gaussian_values / sqrt(2)));
uniform_values = min(max(uniform_values, eps), 1 - eps);
end

function gaussian_values = localCorrelatedGaussianMatrix(positions_2d, d_cor, num_fields)
num_pos = size(positions_2d, 1);
if num_pos == 0
    gaussian_values = [];
    return;
end

corr_matrix = localSpatialCorrelationMatrix(positions_2d, d_cor);
sqrt_corr = chol(corr_matrix, 'lower');
gaussian_values = sqrt_corr * randn(num_pos, num_fields);
end

function d_cor = localProcedureBCorrelationDistance(scenario, state, field_name)
if strcmp(field_name, 'shadow')
    d_cor = localLspCorrelationDistance(scenario, state, 1);
    return;
end
if any(strcmp(field_name, {'AOA', 'ZOA'}))
    d_cor = 50;
    return;
end

[mu_lgDS, sigma_lgDS] = localProcedureBDsStats(scenario, state);
d_cor = 3e8 * 2 * 10^(mu_lgDS + sigma_lgDS);
end

function [mu_lgDS, sigma_lgDS] = localProcedureBDsStats(scenario, state)
fc_eval = max(scenario.frequency/1e9, 6);
switch scenario.name
    case 'UrbanGrid'
        if state == "LOS"
            mu_lgDS = -7.067 - 0.0794*log10(1 + fc_eval);
            sigma_lgDS = 0.57 + 0.026*log10(fc_eval);
        else
            mu_lgDS = -6.47 - 0.134*log10(fc_eval);
            sigma_lgDS = 0.39;
        end
    case 'InH'
        if state == "LOS"
            mu_lgDS = -0.01*log10(1 + fc_eval) - 7.692;
            sigma_lgDS = 0.18;
        else
            mu_lgDS = -0.28*log10(1 + fc_eval) - 7.173;
            sigma_lgDS = 0.1*log10(fc_eval + 1) + 0.055;
        end
    otherwise
        error('Unsupported scenario for Procedure B DS statistics: %s.', scenario.name);
end
end

function localAssignProcedureBRaw(equipment, st_id, raw)
if isempty(equipment.SC_procB_raw)
    proc_table = repmat(localEmptyProcedureBRaw(numel(raw.tau)), 1, st_id);
else
    proc_table = equipment.SC_procB_raw;
end
if numel(proc_table) < st_id
    proc_table(end+1:st_id) = repmat(localEmptyProcedureBRaw(numel(raw.tau)), 1, st_id - numel(proc_table));
end
proc_table(st_id) = raw;
equipment.SC_procB_raw = proc_table;
end

function state = localTargetLinkState(equipment, st, scenario)
if isprop(st, 'Indoor') && st.Indoor
    state = "O2I";
    return;
end

if strcmp(scenario.name, 'UrbanGrid')
    state = "LOS";
    return;
end

if strcmp(scenario.name, 'InH')
    d2d_in = norm(st.Position(1:2) - equipment.Position(1:2));
    switch scenario.subname
        case {'B', 'open_office'}
            pr_los = min(max(exp(-(d2d_in - 5)/70.8), 0.54*exp(-(d2d_in - 49)/211.7)), 1);
        case 'mixed_office'
            pr_los = min(max(exp(-(d2d_in - 1.2)/4.7), 0.32*exp(-(d2d_in - 6.5)/32.6)), 1);
        otherwise
            pr_los = min(max(exp(-(d2d_in - 5)/70.8), 0.54*exp(-(d2d_in - 49)/211.7)), 1);
    end
elseif strcmp(scenario.name, 'UrbanGrid')
    d2d = norm(st.Position(1:2) - equipment.Position(1:2));
    if scenario.frequency < 30e9
        pr_los = min(1, 1.05 * exp(-0.0114 * d2d));
    elseif d2d <= 18
        pr_los = 1;
    else
        pr_los = 18/d2d + exp(-d2d/63)*(1 - 18/d2d);
    end
else
    error('Unsupported scenario for target LSP spatial consistency: %s.', scenario.name);
end

rand_los = localRandLos(st, equipment.ID);
if rand_los < pr_los
    state = "LOS";
else
    state = "NLOS";
end
end

function value = localRandLos(st, id)
if isempty(st.rand_LoS) || id < 1 || id > numel(st.rand_LoS)
    value = rand;
else
    value = st.rand_LoS(id);
end
end

function lsp_raw = localCorrelatedLspRaw(positions, state, scenario)
if isempty(positions)
    lsp_raw = [];
    return;
end

vector_len = localLspVectorLength(state);
lsp_raw = nan(vector_len, size(positions, 1));
for param_idx = 1:vector_len
    d_cor = localLspCorrelationDistance(scenario, state, param_idx);
    corr_matrix = localSpatialCorrelationMatrix(positions(:, 1:2), d_cor);
    sqrt_corr = chol(corr_matrix, 'lower');
    lsp_raw(param_idx, :) = (sqrt_corr * randn(size(positions, 1), 1)).';
end
end

function n = localLspVectorLength(state)
if state == "LOS"
    n = 7;
else
    n = 6;
end
end

function corr_matrix = localSpatialCorrelationMatrix(positions_2d, d_cor)
num_pos = size(positions_2d, 1);
corr_matrix = eye(num_pos);
for row = 1:num_pos
    for col = row+1:num_pos
        distance = norm(positions_2d(row, :) - positions_2d(col, :));
        corr_value = exp(-distance / d_cor);
        corr_matrix(row, col) = corr_value;
        corr_matrix(col, row) = corr_value;
    end
end
corr_matrix = corr_matrix + 1e-10*eye(num_pos);
end

function d_cor = localLspCorrelationDistance(scenario, state, param_idx)
% LSP order: LOS [SF K DS ASD ASA ZSD ZSA], NLOS/O2I [SF DS ASD ASA ZSD ZSA].
switch scenario.name
    case 'UrbanGrid'
        d_cor = localUrbanGridCorrelationDistance(state, param_idx);
    case 'InH'
        d_cor = localInHCorrelationDistance(state, param_idx);
    otherwise
        error('Unsupported scenario for LSP spatial consistency: %s.', scenario.name);
end
end

function d_cor = localUrbanGridCorrelationDistance(state, param_idx)
if state == "LOS"
    distances = [37, 12, 30, 18, 15, 15, 15];
elseif state == "NLOS"
    distances = [50, 40, 50, 50, 50, 50];
else
    distances = [7, 10, 11, 17, 25, 25];
end
d_cor = distances(param_idx);
end

function d_cor = localInHCorrelationDistance(state, param_idx)
if state == "LOS"
    distances = [10, 4, 8, 7, 5, 4, 4];
elseif state == "NLOS"
    distances = [6, 5, 3, 3, 4, 4];
else
    error('InH O2I LSP spatial consistency is not defined for this calibration.');
end
d_cor = distances(param_idx);
end

function [metric_keys, metric_labels, theory] = localLspDebugMetricConfig(scenario, distance_axis)
states = ["LOS", "NLOS"];
metric_keys = {};
metric_labels = strings(1, 0);
theory = [];

for state_idx = 1:numel(states)
    state = states(state_idx);
    param_names = localLspParamNames(state);
    for param_idx = 1:numel(param_names)
        metric_keys{end + 1} = localLspDebugFieldName(state, param_names(param_idx)); %#ok<AGROW>
        metric_labels(end + 1) = state + "_" + param_names(param_idx); %#ok<AGROW>
        d_cor = localLspCorrelationDistance(scenario, state, param_idx);
        theory(end + 1, :) = exp(-distance_axis / d_cor); %#ok<AGROW>
    end
end
end

function [metric_keys, metric_labels, theory] = localSspDebugMetricConfig(scenario, distance_axis)
states = ["LOS", "NLOS"];
field_names = ["tau", "AOA", "AOD", "ZOA", "ZOD", "shadow"];
metric_keys = {};
metric_labels = strings(1, 0);
theory = [];

for state_idx = 1:numel(states)
    state = states(state_idx);
    for field_idx = 1:numel(field_names)
        field_name = field_names(field_idx);
        metric_keys{end + 1} = localSspDebugFieldName(state, field_name); %#ok<AGROW>
        metric_labels(end + 1) = state + "_" + field_name; %#ok<AGROW>
        d_cor = localProcedureBCorrelationDistance(scenario, state, char(field_name));
        theory(end + 1, :) = exp(-distance_axis / d_cor); %#ok<AGROW>
    end
end
end

function field_name = localSspDebugFieldName(state, proc_field)
field_name = char("ssp_" + state + "_" + proc_field);
end

function param_names = localLspParamNames(state)
if state == "LOS"
    param_names = ["SF", "K", "DS", "ASD", "ASA", "ZSD", "ZSA"];
else
    param_names = ["SF", "DS", "ASD", "ASA", "ZSD", "ZSA"];
end
end

function field_name = localLspDebugFieldName(state, param_name)
field_name = char("lsp_" + state + "_" + param_name);
end

function metrics = localCollectTargetUeLspRawMetrics(UE_list, ST_list, scenario)
[metric_keys, ~, ~] = localLspDebugMetricConfig(scenario, 0);
num_ue = numel(UE_list);
metrics = struct();
for metric_idx = 1:numel(metric_keys)
    metrics.(metric_keys{metric_idx}) = nan(num_ue, 1);
end

if isempty(ST_list)
    return;
end

st_id = ST_list(1).ID;
for ue_idx = 1:num_ue
    ue = UE_list(ue_idx);
    for state = ["LOS", "NLOS"]
        source_field = char("LSP_raw_" + state);
        if isempty(ue.(source_field)) || size(ue.(source_field), 2) < st_id
            continue;
        end
        raw_values = ue.(source_field)(:, st_id);
        param_names = localLspParamNames(state);
        for param_idx = 1:min(numel(param_names), numel(raw_values))
            field_name = localLspDebugFieldName(state, param_names(param_idx));
            metrics.(field_name)(ue_idx) = raw_values(param_idx);
        end
    end
end
end

function metrics = localCollectTargetUeProcedureBRawMetrics(UE_list, ST_list, scenario)
[metric_keys, ~, ~] = localSspDebugMetricConfig(scenario, 0);
num_ue = numel(UE_list);
metrics = struct();
for metric_idx = 1:numel(metric_keys)
    metrics.(metric_keys{metric_idx}) = nan(num_ue, 1);
end

if isempty(ST_list)
    return;
end

st_id = ST_list(1).ID;
for ue_idx = 1:num_ue
    ue = UE_list(ue_idx);
    if isempty(ue.SC_procB_raw) || numel(ue.SC_procB_raw) < st_id
        continue;
    end

    state = localTargetLinkState(ue, ST_list(1), scenario);
    proc_raw = ue.SC_procB_raw(st_id);
    for proc_field = ["tau", "AOA", "AOD", "ZOA", "ZOD", "shadow"]
        if ~isfield(proc_raw, char(proc_field)) || isempty(proc_raw.(char(proc_field)))
            continue;
        end
        field_name = localSspDebugFieldName(state, proc_field);
        values = proc_raw.(char(proc_field));
        if numel(values) >= 3
            metrics.(field_name)(ue_idx) = values(3);
        end
    end
end
end

function localAssignLspRaw(equipment, state, object_id, lsp_raw)
field_name = char("LSP_raw_" + state);
if isempty(equipment.(field_name))
    lsp_table = nan(size(lsp_raw, 1), object_id);
else
    lsp_table = equipment.(field_name);
end
if size(lsp_table, 1) < size(lsp_raw, 1)
    lsp_table(end+1:size(lsp_raw, 1), :) = nan;
end
if size(lsp_table, 2) < object_id
    lsp_table(:, end+1:object_id) = nan;
end
lsp_table(1:size(lsp_raw, 1), object_id) = lsp_raw(:);
equipment.(field_name) = lsp_table;
end

function metrics = localEmptyMetrics(n)
metrics = struct();
metrics.delay = nan(n, 1);
metrics.aoa = nan(n, 1);
metrics.los = nan(n, 1);
end

function metrics = localStoreMetrics(metrics, idx, link, ue, st_id, scenario)
metrics.delay(idx) = localThirdClusterDelayMetric(link);
metrics.aoa(idx) = localRawThirdClusterAoaMetric(ue, st_id, scenario);
% Table 7.9.6.3-1 LOS/NLOS-status metric is evaluated for the ST-UT leg.
metrics.los(idx) = double(numel(link.bLOS) >= 2 && logical(link.bLOS(2)));
end

function value = localThirdClusterDelayMetric(link)
txrx = 2;          % SPST-SRX, i.e. target-UE side.
cluster_idx = 3;   % Table 7.9.6.3-1 metric: third cluster.
if isempty(link.tau_n_keep) || size(link.tau_n_keep, 1) < txrx || size(link.tau_n_keep, 2) < cluster_idx
    value = nan;
    return;
end
value = link.tau_n_keep(txrx, cluster_idx);
end

function value = localRawThirdClusterAoaMetric(ue, st_id, scenario)
cluster_idx = 3;   % Table 7.9.6.3-1 metric: third cluster.
value = nan;
if isempty(ue.SC_procB_raw) || numel(ue.SC_procB_raw) < st_id
    return;
end

proc_raw = ue.SC_procB_raw(st_id);
if ~isfield(proc_raw, 'AOA') || numel(proc_raw.AOA) < cluster_idx
    return;
end

normalized_aoa = 2*proc_raw.AOA(cluster_idx) - 1;
if strcmp(scenario.name, 'InH')
    % Diagnostic: keep the original Procedure-B raw AOA sample, but treat
    % its normalized [-1,1] range as a circular [-pi,pi] coordinate.
    value = exp(1j*pi*normalized_aoa);
else
    value = normalized_aoa;
end
end

function [xcorr_metrics, pair_count_per_bin, pooled_stats] = localPairDistanceCorrelation(drop_data, distance_axis, metric_keys, los_only)
metrics = drop_data.metrics;
ue_pos = drop_data.ue_positions(:, 1:2);
num_ue = size(ue_pos, 1);
pair_count = num_ue * (num_ue - 1) / 2;
pair_distance = nan(pair_count, 1);
pair_i = nan(pair_count, 1);
pair_j = nan(pair_count, 1);
pair_idx = 0;

for i = 1:num_ue-1
    for j = i+1:num_ue
        pair_idx = pair_idx + 1;
        pair_i(pair_idx) = i;
        pair_j(pair_idx) = j;
        pair_distance(pair_idx) = norm(ue_pos(i, :) - ue_pos(j, :));
    end
end

if los_only
    los_pair = metrics.los(pair_i) == 1 & metrics.los(pair_j) == 1;
else
    los_pair = true(size(pair_distance));
end

xcorr_metrics = nan(numel(metric_keys), numel(distance_axis));
pair_count_per_bin = zeros(1, numel(distance_axis));
pooled_stats = localEmptyPooledStats(numel(metric_keys), numel(distance_axis));
for bin_idx = 1:numel(distance_axis)
    bin_start = distance_axis(bin_idx);
    if bin_idx == numel(distance_axis)
        in_bin = pair_distance >= bin_start & pair_distance <= bin_start + 1;
    else
        in_bin = pair_distance >= bin_start & pair_distance < bin_start + 1;
    end
    in_bin = in_bin & los_pair;
    pair_count_per_bin(bin_idx) = nnz(in_bin);

    for metric_idx = 1:numel(metric_keys)
        metric_values = metrics.(metric_keys{metric_idx});
        x = metric_values(pair_i(in_bin));
        y = metric_values(pair_j(in_bin));
        xcorr_metrics(metric_idx, bin_idx) = localPairXcorrFromValues(x, y);
        pooled_stats = localAccumulatePooledStats(pooled_stats, metric_idx, bin_idx, x, y);
    end
end
end

function value = localPairXcorrFromValues(x, y)
if numel(x) < 2
    value = nan;
    return;
end
value = localXcorr(x, y);
end

function value = localXcorr(x, y)
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if numel(x) < 2 || std(x) == 0 || std(y) == 0
    value = nan;
    return;
end
if isreal(x) && isreal(y)
    value = (mean(x.*y)-mean(x)*mean(y)) / sqrt(mean(x.^2)-mean(x)^2) / sqrt(mean(y.^2)-mean(y)^2);
else
    value = abs((mean(x.*conj(y))-mean(x)*conj(mean(y))) / sqrt(mean(x.*conj(x))-mean(x)*conj(mean(x))) / sqrt(mean(y.*conj(y))-mean(y)*conj(mean(y))));
end
end

function stats = localEmptyPooledStats(num_metrics, num_bins)
stats = struct();
stats.count = zeros(num_metrics, num_bins);
stats.sum_x = zeros(num_metrics, num_bins);
stats.sum_y = zeros(num_metrics, num_bins);
stats.sum_x2 = zeros(num_metrics, num_bins);
stats.sum_y2 = zeros(num_metrics, num_bins);
stats.sum_xy = zeros(num_metrics, num_bins);
end

function stats = localAccumulatePooledStats(stats, metric_idx, bin_idx, x, y)
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if isempty(x)
    return;
end
stats.count(metric_idx, bin_idx) = stats.count(metric_idx, bin_idx) + numel(x);
stats.sum_x(metric_idx, bin_idx) = stats.sum_x(metric_idx, bin_idx) + sum(x);
stats.sum_y(metric_idx, bin_idx) = stats.sum_y(metric_idx, bin_idx) + sum(y);
if isreal(x) && isreal(y)
    stats.sum_x2(metric_idx, bin_idx) = stats.sum_x2(metric_idx, bin_idx) + sum(x.^2);
    stats.sum_y2(metric_idx, bin_idx) = stats.sum_y2(metric_idx, bin_idx) + sum(y.^2);
    stats.sum_xy(metric_idx, bin_idx) = stats.sum_xy(metric_idx, bin_idx) + sum(x.*y);
else
    stats.sum_x2(metric_idx, bin_idx) = stats.sum_x2(metric_idx, bin_idx) + sum(x.*conj(x));
    stats.sum_y2(metric_idx, bin_idx) = stats.sum_y2(metric_idx, bin_idx) + sum(y.*conj(y));
    stats.sum_xy(metric_idx, bin_idx) = stats.sum_xy(metric_idx, bin_idx) + sum(x.*conj(y));
end
end

function stats = localMergePooledStats(stats, new_stats)
fields = {'count', 'sum_x', 'sum_y', 'sum_x2', 'sum_y2', 'sum_xy'};
for field_idx = 1:numel(fields)
    field_name = fields{field_idx};
    stats.(field_name) = stats.(field_name) + new_stats.(field_name);
end
end

function xcorr_pooled = localPooledXcorr(stats)
xcorr_pooled = nan(size(stats.count));
for metric_idx = 1:size(stats.count, 1)
    for bin_idx = 1:size(stats.count, 2)
        n = stats.count(metric_idx, bin_idx);
        if n < 2
            continue;
        end
        mean_x = stats.sum_x(metric_idx, bin_idx) / n;
        mean_y = stats.sum_y(metric_idx, bin_idx) / n;
        is_complex_metric = ~isreal(mean_x) || ~isreal(mean_y) || ...
            ~isreal(stats.sum_xy(metric_idx, bin_idx));
        if is_complex_metric
            var_x = real(stats.sum_x2(metric_idx, bin_idx) / n - mean_x*conj(mean_x));
            var_y = real(stats.sum_y2(metric_idx, bin_idx) / n - mean_y*conj(mean_y));
        else
            var_x = stats.sum_x2(metric_idx, bin_idx) / n - mean_x^2;
            var_y = stats.sum_y2(metric_idx, bin_idx) / n - mean_y^2;
        end
        if var_x <= 0 || var_y <= 0
            continue;
        end
        if is_complex_metric
            cov_xy = stats.sum_xy(metric_idx, bin_idx) / n - mean_x * conj(mean_y);
            xcorr_pooled(metric_idx, bin_idx) = abs(cov_xy) / sqrt(var_x) / sqrt(var_y);
        else
            cov_xy = stats.sum_xy(metric_idx, bin_idx) / n - mean_x * mean_y;
            xcorr_pooled(metric_idx, bin_idx) = cov_xy / sqrt(var_x) / sqrt(var_y);
        end
    end
end
end

function localPrintXcorrComparison(metric_keys, distance_axis, xcorr_equal_avg, xcorr_pooled, pair_count_per_bin)
fprintf('  Correlation aggregation comparison:\n');
for metric_idx = 1:numel(metric_keys)
    valid = isfinite(xcorr_equal_avg(metric_idx, :)) & isfinite(xcorr_pooled(metric_idx, :));
    if any(valid)
        mae_pct = mean(abs(xcorr_equal_avg(metric_idx, valid) - xcorr_pooled(metric_idx, valid)), 'omitnan') * 100;
    else
        mae_pct = nan;
    end
    fprintf('    %s: mean |equal-drop - pooled| = %.2f %%\n', metric_keys{metric_idx}, mae_pct);
end

print_bins = 1:min(10, numel(distance_axis));
fprintf('    first bins distance/count: ');
fprintf('%.0fm/%d ', [distance_axis(print_bins); pair_count_per_bin(print_bins)]);
fprintf('\n');
end

function localPrintLspDebugSummary(metric_labels, lsp_xcorr_pooled, lsp_theory)
fprintf('  LSP raw spatial correlation check against exp(-d/d_cor):\n');
for metric_idx = 1:numel(metric_labels)
    valid = isfinite(lsp_xcorr_pooled(metric_idx, :)) & isfinite(lsp_theory(metric_idx, :));
    if any(valid)
        mae_pct = mean(abs(lsp_xcorr_pooled(metric_idx, valid) - lsp_theory(metric_idx, valid)), 'omitnan') * 100;
        near_idx = find(valid, 1, 'first');
        fprintf('    %s: MAE=%.2f %%, first-bin pooled/theory=%.1f/%.1f %%\n', ...
            metric_labels(metric_idx), mae_pct, ...
            lsp_xcorr_pooled(metric_idx, near_idx) * 100, lsp_theory(metric_idx, near_idx) * 100);
    else
        fprintf('    %s: no valid pairs\n', metric_labels(metric_idx));
    end
end
end

function sheet = localSheetName(fc)
if fc == 30e9
    sheet = 'TRP-UE-30GHz';
else
    sheet = 'TRP-UE-6GHz';
end
end
