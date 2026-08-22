clc;
clearvars -except use_isac_frequency_preset isac_frequency_preset custom_frequency_config custom_bs_position custom_ue_position custom_st_position
close all;

%% Basic ISAC channel scaffold
enable_full_channel = true;
enable_rp = true;
current_time = 0;
plot_controller = true;
if ~exist('use_isac_frequency_preset', 'var')
    use_isac_frequency_preset = true;
end
if ~exist('isac_frequency_preset', 'var')
    isac_frequency_preset = 'ISAC_FR1';
end
if ~exist('custom_frequency_config', 'var')
    custom_frequency_config = struct();
end
if ~exist('custom_bs_position', 'var')
    custom_bs_position = [];
end
if ~exist('custom_ue_position', 'var')
    custom_ue_position = [];
end
if ~exist('custom_st_position', 'var')
    custom_st_position = [];
end

%% Choose communication scenario: UMa, UMi, RMa, InF(SL/DL/SH/DH/HH), InH, UrbanGrid, HighWay
scenario = comm_scenario.UrbanGrid(500, 2);
if use_isac_frequency_preset
    scenario.applyIsacFrequencyPreset(isac_frequency_preset, custom_frequency_config);
end

%% Choose sensing type: UAV(scenario), Human, Vehicle, AGV
sensing_type = sensing_types.Human;

%% Drop network nodes
[BS_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller, custom_bs_position);
[UE_list, UE_pos_list, ~] = network_layout.Drop_UE_ISAC(BS_list, scenario, plot_controller, custom_ue_position);
st_reference_ue_pos = UE_pos_list(1:min(scenario.UE_per_sec, size(UE_pos_list, 1)), :);
[ST_list, ~] = network_layout.Drop_ST( ...
    BS_list, st_reference_ue_pos, scenario, sensing_type, plot_controller, custom_st_position);

if scenario.spatial_consistency_enable && strcmpi(scenario.spatial_consistency_procedure, 'A')
    localApplyProcedureAXnSpatialConsistency(BS_list, UE_list, ST_list, scenario);
end

%% Build initial channel links
link_list_UE = [];
link_list_ST = [];
link_list_RP = [];

for BS_num = 1:numel(BS_list)
    STX = BS_list(BS_num);
    SRX = STX;

    for UE_num = 1:numel(UE_list)
        link_UE = channel.Comm_channel(STX, UE_list(UE_num), scenario, enable_full_channel, current_time);
        link_list_UE = [link_list_UE; link_UE]; %#ok<AGROW>
    end

    for ST_num = 1:numel(ST_list)
        ST = ST_list(ST_num);
        link_ST = channel.Target_channel(STX, SRX, ST, scenario, enable_full_channel, current_time);
        link_list_ST = [link_list_ST; link_ST]; %#ok<AGROW>
    end

    if enable_rp
        attached_equipment = STX;
        [RP_list, ~] = network_layout.Drop_RP(attached_equipment, scenario);
        attached_equipment.RP = RP_list;

        for rp_num = 1:scenario.RP_per_equipment
            link_RP = channel.Comm_channel(attached_equipment, RP_list(rp_num), scenario, enable_full_channel, current_time);
            link_list_RP = [link_list_RP; link_RP]; %#ok<AGROW>
        end
    end
end

%% Save generated link lists
output_dir = fullfile(pwd, 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_file = fullfile(output_dir, ['ISAC_channel_link_lists_' timestamp '.mat']);
save(output_file, 'link_list_UE', 'link_list_ST', 'link_list_RP');

function localApplyProcedureAXnSpatialConsistency(BS_list, UE_list, ST_list, scenario)
max_clusters = 20;
d_cor = localProcedureAXnCorrelationDistance(scenario);

for tx_idx = 1:numel(BS_list)
    tx = BS_list(tx_idx);
    rx_objects = num2cell(UE_list(:));
    localAssignCommProcedureAXn(tx, rx_objects, scenario, d_cor, max_clusters);
    localAssignTargetTxProcedureAXn(tx, ST_list, scenario, d_cor, max_clusters);
end

target_rx_objects = [num2cell(BS_list(:)); num2cell(UE_list(:))];
localAssignTargetRxProcedureAXn(target_rx_objects, ST_list, scenario, d_cor, max_clusters);
end

function localAssignCommProcedureAXn(tx, rx_objects, scenario, d_cor, max_clusters)
if isempty(rx_objects)
    return;
end
if ~iscell(rx_objects)
    rx_objects = num2cell(rx_objects(:));
end

states = strings(numel(rx_objects), 1);
positions = nan(numel(rx_objects), 3);
for obj_idx = 1:numel(rx_objects)
    rx = rx_objects{obj_idx};
    states(obj_idx) = localCommLinkState(tx, rx, scenario);
    positions(obj_idx, :) = rx.Position;
end

for state = unique(states).'
    group_idx = find(states == state);
    x_n = localCorrelatedProcedureAXn(positions(group_idx, 1:2), d_cor, max_clusters);
    for local_idx = 1:numel(group_idx)
        rx = rx_objects{group_idx(local_idx)};
        localAssignXnTable(tx, 'SC_procA_comm_Xn', rx.ID, x_n(local_idx, :));
    end
end
end

function localAssignTargetTxProcedureAXn(tx, ST_list, scenario, d_cor, max_clusters)
if isempty(ST_list)
    return;
end

states = strings(numel(ST_list), 1);
positions = nan(numel(ST_list), 3);
for st_idx = 1:numel(ST_list)
    st = ST_list(st_idx);
    states(st_idx) = localTargetLinkState(tx, st, scenario);
    positions(st_idx, :) = st.Position;
end

for state = unique(states).'
    group_idx = find(states == state);
    x_n = localCorrelatedProcedureAXn(positions(group_idx, 1:2), d_cor, max_clusters);
    for local_idx = 1:numel(group_idx)
        st = ST_list(group_idx(local_idx));
        localAssignXnTable(tx, 'SC_procA_target_tx_Xn', st.ID, x_n(local_idx, :));
    end
end
end

function localAssignTargetRxProcedureAXn(rx_objects, ST_list, scenario, d_cor, max_clusters)
if isempty(rx_objects) || isempty(ST_list)
    return;
end
if ~iscell(rx_objects)
    rx_objects = num2cell(rx_objects(:));
end

for st_idx = 1:numel(ST_list)
    st = ST_list(st_idx);
    states = strings(numel(rx_objects), 1);
    positions = nan(numel(rx_objects), 3);
    for obj_idx = 1:numel(rx_objects)
        rx = rx_objects{obj_idx};
        states(obj_idx) = localTargetLinkState(rx, st, scenario);
        positions(obj_idx, :) = rx.Position;
    end

    for state = unique(states).'
        group_idx = find(states == state);
        x_n = localCorrelatedProcedureAXn(positions(group_idx, 1:2), d_cor, max_clusters);
        for local_idx = 1:numel(group_idx)
            rx = rx_objects{group_idx(local_idx)};
            localAssignXnTable(rx, 'SC_procA_target_rx_Xn', st.ID, x_n(local_idx, :));
        end
    end
end
end

function x_n = localCorrelatedProcedureAXn(positions_2d, d_cor, max_clusters)
num_pos = size(positions_2d, 1);
if num_pos == 0
    x_n = [];
    return;
end

corr_matrix = localBinarySignGaussianCorrelationMatrix(positions_2d, d_cor);
sqrt_corr = chol(corr_matrix, 'lower');
gaussian_values = sqrt_corr * randn(num_pos, max_clusters);
x_n = ones(num_pos, max_clusters);
x_n(gaussian_values < 0) = -1;
end

function localAssignXnTable(owner, field_name, object_id, x_n)
if isempty(owner.(field_name))
    x_n_table = nan(object_id, numel(x_n));
else
    x_n_table = owner.(field_name);
end
if size(x_n_table, 1) < object_id
    x_n_table(end+1:object_id, :) = nan;
end
if size(x_n_table, 2) < numel(x_n)
    x_n_table(:, end+1:numel(x_n)) = nan;
end
x_n_table(object_id, 1:numel(x_n)) = x_n;
owner.(field_name) = x_n_table;
end

function state = localCommLinkState(tx, rx, scenario)
if isprop(rx, 'Indoor') && rx.Indoor
    state = "O2I";
    return;
end

if localRandLos(rx, tx.ID) < localLosProbability(tx, rx, scenario)
    state = "LOS";
else
    state = "NLOS";
end
end

function state = localTargetLinkState(equipment, st, scenario)
if localRandLos(st, equipment.ID) < localLosProbability(equipment, st, scenario)
    state = "LOS";
else
    state = "NLOS";
end
end

function value = localRandLos(object, id)
if ~isprop(object, 'rand_LoS') || isempty(object.rand_LoS) || id < 1 || id > numel(object.rand_LoS)
    value = rand;
else
    value = object.rand_LoS(id);
end
end

function pr_los = localLosProbability(tx, rx, scenario)
d2d = norm(rx.Position(1:2) - tx.Position(1:2));
if isprop(rx, 'height') && ~isempty(rx.height)
    rx_height = rx.height;
else
    rx_height = rx.Position(3);
end

switch scenario.name
    case {'UMa', 'UrbanGrid'}
        if strcmp(scenario.name, 'UrbanGrid') && scenario.frequency < 30e9
            pr_los = min(1, 1.05 * exp(-0.0114 * d2d));
        elseif d2d <= 18
            pr_los = 1;
        else
            c_tmp = max((rx_height - 13)/10, 0)^1.5;
            pr_los = (18/d2d + exp(-d2d/63)*(1 - 18/d2d)) * ...
                (1 + c_tmp*(5/4)*((d2d/100)^3)*exp(-d2d/150));
        end
    case 'UMi'
        pr_los = min(18/d2d, 1)*(1 - exp(-d2d/36)) + exp(-d2d/36);
    case 'RMa'
        pr_los = min(exp(-(d2d - 10)/1000), 1);
    case 'InH'
        switch scenario.subname
            case {'B', 'open_office'}
                pr_los = min(max(exp(-(d2d - 5)/70.8), 0.54*exp(-(d2d - 49)/211.7)), 1);
            case 'mixed_office'
                pr_los = min(max(exp(-(d2d - 1.2)/4.7), 0.32*exp(-(d2d - 6.5)/32.6)), 1);
            otherwise
                pr_los = min(max(exp(-(d2d - 5)/70.8), 0.54*exp(-(d2d - 49)/211.7)), 1);
        end
    otherwise
        pr_los = 1;
end
end

function d_cor = localProcedureAXnCorrelationDistance(scenario)
switch scenario.name
    case 'RMa'
        d_cor = 60;
    case {'UMa', 'UrbanGrid'}
        d_cor = 50;
    case 'UMi'
        d_cor = 15;
    case {'InH', 'InF'}
        d_cor = 10;
    otherwise
        d_cor = 50;
end
end

function corr_matrix = localBinarySignGaussianCorrelationMatrix(positions_2d, d_cor)
num_pos = size(positions_2d, 1);
corr_matrix = eye(num_pos);
for row = 1:num_pos
    for col = row+1:num_pos
        distance = norm(positions_2d(row, :) - positions_2d(col, :));
        target_corr = exp(-distance / d_cor);
        gaussian_corr = sin((pi/2) * target_corr);
        corr_matrix(row, col) = gaussian_corr;
        corr_matrix(col, row) = gaussian_corr;
    end
end
corr_matrix = corr_matrix + 1e-10*eye(num_pos);
end
