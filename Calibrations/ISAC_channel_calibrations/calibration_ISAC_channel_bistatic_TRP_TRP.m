% TRP-TRP bistatic snapshot calibration against the 3GPP company curves.
% Best-N selection uses unique physical TRP pairs by default so reciprocal
clc;
% Start from the case configuration in this script.  Keeping old control
% variables from the base workspace can silently rerun a previous case.
clearvars;

% close all;

calibration_dir = fileparts(mfilename('fullpath'));
addpath(calibration_dir);
paths = get_ISAC_calibration_paths();
addpath(paths.isacRoot);

%% Calibration controls
% Select exactly one calibration case from localCalibrationCases() below.
% 1: UAV-UMa-AV
% 2: Human Outdoor-UMa
% 3: Human Outdoor-UMi
% 4: Human Indoor-InH
% 5: Human Indoor-InF-SH
% 6: AGV-InF-SH
% 7: Automotive-Urban Grid
case_id = 5;
if ~exist('full_cali', 'var')
    full_cali = true;
end
if ~exist('do_background', 'var')
    do_background = false;
    
    
end
if ~exist('plot_controller', 'var')
    plot_controller = true;
end
calibration_root = paths.referenceRoot;
if ~exist('current_time', 'var')
    current_time = 0;
end

if ~exist('use_isac_frequency_preset', 'var')
    use_isac_frequency_preset = true;
end
if ~exist('isac_frequency_preset', 'var')
    isac_frequency_preset = 'ISAC_FR2';
end
if ~exist('custom_frequency_config', 'var')
    custom_frequency_config = struct();
end

if ~exist('run_all_calibration_cases', 'var')
    run_all_calibration_cases = false;
end
if ~exist('calibration_simulation_times', 'var')
    calibration_simulation_times = [];
end
if ~exist('calibration_frequency_presets', 'var')
    calibration_frequency_presets = {'ISAC_FR1', 'ISAC_FR2'};
end
if ~exist('trp_pair_selection_mode', 'var')
    trp_pair_selection_mode = 'unique_direction_balanced';
end
calibration_cases = localCalibrationCases();
case_id_list = localCaseIdList(case_id, run_all_calibration_cases, numel(calibration_cases));
frequency_preset_list = localFrequencyPresetList(isac_frequency_preset, ...
    run_all_calibration_cases, calibration_frequency_presets);

for case_pos = 1:numel(case_id_list)
    active_case_id = case_id_list(case_pos);
    for preset_idx = 1:numel(frequency_preset_list)
        active_preset = frequency_preset_list{preset_idx};
        localRunCalibrationCase( ...
            calibration_cases(active_case_id), active_case_id, active_preset, use_isac_frequency_preset, ...
            custom_frequency_config, full_cali, do_background, current_time, plot_controller, ...
            calibration_root, trp_pair_selection_mode, calibration_simulation_times);
    end
end

function case_id_list = localCaseIdList(case_id, run_all_calibration_cases, num_cases)
if run_all_calibration_cases
    case_id_list = 1:num_cases;
else
    case_id_list = case_id(:).';
end

if ~isnumeric(case_id_list) || any(case_id_list ~= fix(case_id_list)) ...
        || any(case_id_list < 1) || any(case_id_list > num_cases)
    error('case_id must contain integer case IDs in the range 1 to %d.', num_cases);
end
end

function frequency_preset_list = localFrequencyPresetList( ...
    isac_frequency_preset, run_all_calibration_cases, calibration_frequency_presets)
if run_all_calibration_cases
    frequency_preset_list = calibration_frequency_presets;
else
    frequency_preset_list = localCellstrList(isac_frequency_preset, 'isac_frequency_preset');
end
end

function value_list = localCellstrList(value, value_name)
if ischar(value)
    value_list = cellstr(value);
elseif isstring(value)
    value_list = cellstr(value(:).');
elseif iscell(value)
    value_list = value(:).';
else
    error('%s must be a char, string, or cell array of character vectors.', value_name);
end

for value_idx = 1:numel(value_list)
    if isstring(value_list{value_idx})
        value_list{value_idx} = char(value_list{value_idx});
    end
    if ~ischar(value_list{value_idx}) || isempty(value_list{value_idx})
        error('%s must contain only nonempty character vectors.', value_name);
    end
end
end

function calibration_cases = localCalibrationCases()
calibration_cases = struct( ...
    'name', {}, ...
    'scenarioFcn', {}, ...
    'sensingFcn', {});

calibration_cases(end + 1) = localCase( ...
    'UAV-UMa-AV', ...
    @() comm_scenario.UMa, ...
    @(scenario) sensing_types.UAV(scenario));

calibration_cases(end + 1) = localCase( ...
    'Human Outdoor-UMa', ...
    @() comm_scenario.UMa, ...
    @(scenario) sensing_types.Human());

calibration_cases(end + 1) = localCase( ...
    'Human Outdoor-UMi', ...
    @() comm_scenario.UMi, ...
    @(scenario) sensing_types.Human());

calibration_cases(end + 1) = localCase( ...
    'Human Indoor-InH', ...
    @() comm_scenario.InH, ...
    @(scenario) sensing_types.Human());

calibration_cases(end + 1) = localCase( ...
    'Human Indoor-InF-SH', ...
    @() comm_scenario.InF('SH'), ...
    @(scenario) sensing_types.Human());

calibration_cases(end + 1) = localCase( ...
    'AGV-InF-SH', ...
    @() comm_scenario.InF('SH'), ...
    @(scenario) sensing_types.AGV());

calibration_cases(end + 1) = localCase( ...
    'Auto-Urban Grid', ...
    @() comm_scenario.UrbanGrid, ...
    @(scenario) sensing_types.Vehicle());
end

function case_config = localCase(case_name, scenario_fcn, sensing_fcn)
case_config = struct();
case_config.name = case_name;
case_config.scenarioFcn = scenario_fcn;
case_config.sensingFcn = sensing_fcn;
end

function [scenario, sensing_type] = localCreateCalibrationObjects(case_config)
scenario = case_config.scenarioFcn();
sensing_type = case_config.sensingFcn(scenario);
end

function localRunCalibrationCase( ...
    case_config, case_id, isac_frequency_preset, use_isac_frequency_preset, ...
    custom_frequency_config, full_cali, do_background, current_time, plot_controller, calibration_root, ...
    trp_pair_selection_mode, calibration_simulation_times)

[scenario, sensing_type] = localCreateCalibrationObjects(case_config);
if use_isac_frequency_preset
    scenario.applyIsacFrequencyPreset(isac_frequency_preset, custom_frequency_config);
end
if ~isempty(calibration_simulation_times)
    validateattributes(calibration_simulation_times, {'numeric'}, ...
        {'scalar', 'integer', 'positive'}, mfilename, 'calibration_simulation_times');
    scenario.simulation_times = calibration_simulation_times;
end
scenario.spatial_consistency_enable = false;
if localIsAgvInFCalibration(scenario, sensing_type)
    % Select distinct physical TRP pairs and balance their STX/SRX roles.
    trp_pair_selection_mode = 'unique_direction_balanced';
elseif localIsAutomotiveUrbanGridCalibration(scenario, sensing_type)
    % TR 38.901 Table 7.9.6.1-3: Urban Grid Automotive calibration uses
    % one UT type per calibration and the general Best-N = 4 rule.  The
    % bundled company reference curves are the pedestrian-UT calibration.
    scenario.user_case = "Pedestrian";
    scenario.best_N = 4;
    trp_pair_selection_mode = 'unique_direction_balanced';
    localValidateAutomotiveUrbanGridScenario(scenario);
end

fprintf('Run TRP-TRP bistatic calibration case %d: %s, preset=%s, pair selection=%s\n', ...
    case_id, case_config.name, string(isac_frequency_preset), string(trp_pair_selection_mode));

sim_data = localRunCalibrationLayout(scenario, sensing_type, full_cali, do_background, current_time, ...
    plot_controller, trp_pair_selection_mode);
[TRP_TRP_bistatic_target_result, selected_link_id] = ...
    localBuildTargetCalibrationResult(sim_data, full_cali);
plot_TRP_mono_result(TRP_TRP_bistatic_target_result, calibration_root);
if do_background
    TRP_TRP_bistatic_background_result = ...
        localBuildBackgroundCalibrationResult(sim_data, full_cali, selected_link_id);
    plot_TRP_mono_result(TRP_TRP_bistatic_background_result, calibration_root);
end
end

function sim_data = localRunCalibrationLayout(scenario, sensing_type, full_cali, do_background, current_time, ...
    plot_controller, trp_pair_selection_mode)
[BS_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller);
if localIsHumanIndoorCalibration(scenario, sensing_type)
    % Tables 7.9.6.1-1/2: Human calibration inherits the single
    % dual-polarized isotropic TRP antenna assumption.
    BS_list = tools.configureCalibrationAntennas( ...
        BS_list, [45, -45], 'TRP');
elseif localIsAgvInFCalibration(scenario, sensing_type)
    % Table 7.9.6.1-4 inherits unspecified antenna parameters from
    % Table 7.9.6.1-1: use one dual-polarized isotropic TRP antenna.
    BS_list = tools.configureCalibrationAntennas( ...
        BS_list, [45, -45], 'AGV TRP');
    BS_list = localConfigureAgvCalibrationTrps(BS_list);
end
trp_pair_idx = localTrpPairs(numel(BS_list), trp_pair_selection_mode);

link_list_ST = [];
link_list_background = [];
UE_list = [];
ST_list = [];
target_links_per_sim = zeros(scenario.simulation_times, 1);

simulation_times = scenario.simulation_times;
original_ST_per_cell = scenario.ST_per_cell;
scenario.ST_per_cell = 1;
st_count_guard = onCleanup(@() setScenarioStPerCell(scenario, original_ST_per_cell));

total_target_links = size(trp_pair_idx, 1) * simulation_times;
progress_bar = waitbar(0, 'creating TRP-TRP bistatic links...');
progress_guard = onCleanup(@() close(progress_bar));
created_target_links = 0;

for sim_idx = 1:simulation_times
    plot_this_drop = plot_controller && sim_idx == 1;
    if localIsHumanIndoorCalibration(scenario, sensing_type)
        % Tables 7.9.6.1-2 and 7.8-1/2: use 20 indoor UTs at
        % 1 m height, uniformly distributed with zero minimum distance.
        ue_pos = localDropHumanIndoorCalibrationUts(scenario, 20);
        [UE_list_drop, ~, ~] = network_layout.Drop_UE_ISAC( ...
            BS_list, scenario, plot_this_drop, ue_pos);
        % For TRP-TRP sensing the target minimum distance is specified
        % only relative to STX/SRX (0 m for indoor Human).  Do not impose
        % a UE-target exclusion region on the convex-hull target drop.
        [ST_list_drop, ~] = network_layout.Drop_ST( ...
            BS_list, [], scenario, sensing_type, plot_this_drop);
    elseif localIsAgvInFCalibration(scenario, sensing_type)
        % Table 7.9.6.1-4: 30 indoor UTs and one AGV.  Let Drop_ST select
        % the target-drop geometry, consistently with the other indoor
        % calibration cases.
        ue_pos = localDropAgvCalibrationUts(BS_list, scenario, 30);
        [UE_list_drop, UE_pos_list, ~] = network_layout.Drop_UE_ISAC( ...
            BS_list, scenario, plot_this_drop, ue_pos);
        [ST_list_drop, ~] = network_layout.Drop_ST( ...
            BS_list, UE_pos_list, scenario, sensing_type, plot_this_drop);
    else
        [UE_list_drop, UE_pos_list, ~] = network_layout.Drop_UE_ISAC( ...
            BS_list, scenario, plot_this_drop);
        if localIsAutomotiveUrbanGridCalibration(scenario, sensing_type)
            % The TRP-TRP minimum-distance rule applies only to STX/SRX and
            % the target.  Pedestrian UEs must not bias the target drop.
            st_reference_ue_pos = [];
        else
            st_reference_ue_pos = UE_pos_list( ...
                1:min(scenario.UE_per_sec, size(UE_pos_list, 1)), :);
        end
        [ST_list_drop, ~] = network_layout.Drop_ST( ...
            BS_list, st_reference_ue_pos, scenario, sensing_type, plot_this_drop);
    end
    if localIsHumanIndoorCalibration(scenario, sensing_type) ...
            || localIsAgvInFCalibration(scenario, sensing_type)
        UE_list_drop = tools.configureCalibrationAntennas( ...
            UE_list_drop, [90, 0], 'UT');
    end
    ST_list_drop(1).ID = sim_idx;

    UE_list = [UE_list; UE_list_drop]; %#ok<AGROW>
    ST_list = [ST_list; ST_list_drop]; %#ok<AGROW>
    ST = ST_list_drop(1);
    target_leg_step24_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    target_leg_small_scale_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

    for pair_num = 1:size(trp_pair_idx, 1)
        tx_idx = trp_pair_idx(pair_num, 1);
        rx_idx = trp_pair_idx(pair_num, 2);
        STX = BS_list(tx_idx);
        SRX = BS_list(rx_idx);

        link_ST = channel.Target_channel(STX, SRX, ST, scenario, full_cali, current_time, ...
            target_leg_step24_cache, target_leg_small_scale_cache);
        link_list_ST = [link_list_ST; link_ST]; %#ok<AGROW>
        if do_background
            link_background = channel.TRP_TRP_background_channel( ...
                STX, SRX, scenario, full_cali, current_time);
            link_list_background = [link_list_background; link_background]; %#ok<AGROW>
        end
        target_links_per_sim(sim_idx) = target_links_per_sim(sim_idx) + 1;

        created_target_links = created_target_links + 1;
        progress_ratio = created_target_links / total_target_links;
        waitbar(progress_ratio, progress_bar, localProgressText(scenario, progress_ratio));
    end
end

sim_data = struct();
sim_data.scenario = scenario;
sim_data.sensing_type = sensing_type;
sim_data.BS_list = BS_list;
sim_data.UE_list = UE_list;
sim_data.ST_list = ST_list;
sim_data.link_list_ST = link_list_ST;
sim_data.link_list_background = link_list_background;
sim_data.target_links_per_sim = target_links_per_sim;
sim_data.trp_pair_idx = trp_pair_idx;
sim_data.trp_pair_selection_mode = trp_pair_selection_mode;
end

function pair_idx = localTrpPairs(num_trp, trp_pair_selection_mode)
if num_trp < 2
    error('TRP-TRP bistatic calibration requires at least two TRPs.');
end

switch char(trp_pair_selection_mode)
    case 'unique_unordered'
        % One representative direction per physical TRP pair.  For LS
        % calibration the reciprocal links have the same PL/SF/RCS score.
        pair_idx = nchoosek(1:num_trp, 2);
    case {'unique_direction_balanced', 'ordered', 'ordered_all'}
        % The balanced mode needs both directions during link generation.
        % Best-N selection later treats reciprocal directions as one
        % physical pair and chooses only one direction from each pair.
        pair_idx = zeros(num_trp * (num_trp - 1), 2);
        row = 0;
        for tx_idx = 1:num_trp
            for rx_idx = 1:num_trp
                if tx_idx == rx_idx
                    continue;
                end
                row = row + 1;
                pair_idx(row, :) = [tx_idx, rx_idx];
            end
        end
    otherwise
        error('Unsupported TRP pair selection mode: %s.', ...
            string(trp_pair_selection_mode));
end
end

function tf = localIsAgvInFCalibration(scenario, sensing_type)
tf = strcmp(scenario.name, 'InF') && strcmp(sensing_type.sensing_type, 'AGV');
end

function tf = localIsHumanIndoorCalibration(scenario, sensing_type)
tf = ismember(scenario.name, {'InH', 'InF'}) ...
    && strcmp(sensing_type.sensing_type, 'Human');
end

function ue_pos = localDropHumanIndoorCalibrationUts(scenario, number_of_uts)
% Tables 7.9.6.1-2 and 7.8-1/2: indoor Human UTs are uniformly
% distributed in the room with zero 2D minimum distance and 1 m height.
ue_pos = zeros(number_of_uts, 3);
ue_pos(:, 1) = scenario.x_range(1) + diff(scenario.x_range) * rand(number_of_uts, 1);
ue_pos(:, 2) = scenario.y_range(1) + diff(scenario.y_range) * rand(number_of_uts, 1);
ue_pos(:, 3) = 1;
end

function tf = localIsAutomotiveUrbanGridCalibration(scenario, sensing_type)
tf = strcmp(scenario.name, 'UrbanGrid') ...
    && strcmp(sensing_type.sensing_type, 'Vehicle');
end

function localValidateAutomotiveUrbanGridScenario(scenario)
% TR 38.901 Tables 7.9.6.1-3 and 7.9.6.2-3.
if scenario.BS_height ~= 25
    error('Urban Grid Automotive calibration requires a 25 m TRP height.');
end
if scenario.frequency == 6e9
    expected_isd = 500;
elseif scenario.frequency == 30e9
    expected_isd = 250;
else
    error('Urban Grid Automotive calibration supports only 6 GHz and 30 GHz.');
end
if scenario.ISD ~= expected_isd
    error('Urban Grid Automotive calibration at %.0f GHz requires ISD=%d m.', ...
        scenario.frequency/1e9, expected_isd);
end
if scenario.best_N ~= 4
    error('Urban Grid Automotive TRP-TRP bistatic calibration requires Best-N=4.');
end
end

function BS_list = localConfigureAgvCalibrationTrps(BS_list)
% Table 7.9.6.1-4: local X-axis and antenna boresight point downward.
for bs_idx = 1:numel(BS_list)
    for sector_idx = 1:numel(BS_list(bs_idx).sector)
        BS_list(bs_idx).sector(sector_idx).antenna.beta = 90;
    end
end
end

function ue_pos = localDropAgvCalibrationUts(BS_list, scenario, number_of_uts)
% Table 7.8-7: uniform indoor UT drop with 1 m minimum 2D TRP distance.
trp_pos = reshape([BS_list.Position], 3, []).';
ue_pos = zeros(number_of_uts, 3);
for ue_idx = 1:number_of_uts
    while true
        candidate = [ ...
            scenario.x_range(1) + diff(scenario.x_range) * rand, ...
            scenario.y_range(1) + diff(scenario.y_range) * rand];
        d_2d = vecnorm(trp_pos(:, 1:2) - candidate, 2, 2);
        if all(d_2d >= 1)
            ue_pos(ue_idx, :) = [candidate, 1.5];
            break;
        end
    end
end
end

function setScenarioStPerCell(scenario, st_per_cell)
scenario.ST_per_cell = st_per_cell;
end

function text = localProgressText(scenario, progress_ratio)
percent = floor(progress_ratio * 1000) / 10;
text = sprintf('(%s, %.0fGHz) %.1f %%', scenario.name, scenario.frequency / 1e9, percent);
end

function [target_result, selected_link_id] = localBuildTargetCalibrationResult(sim_data, full_cali)
scenario = sim_data.scenario;
sensing_type = sim_data.sensing_type;
link_list_ST = sim_data.link_list_ST;
trp_pair_idx = sim_data.trp_pair_idx;
trp_pair_selection_mode = sim_data.trp_pair_selection_mode;
best_N = scenario.best_N;
simulation_times = scenario.simulation_times;
target_links_per_sim = sim_data.target_links_per_sim;
if any(target_links_per_sim ~= target_links_per_sim(1))
    error('TRP-TRP bistatic calibration generated a different number of target links across simulation drops.');
end
links_per_sim = target_links_per_sim(1);
if links_per_sim < best_N
    error('TRP-TRP bistatic calibration generated %d links per drop, but best_N is %d.', links_per_sim, best_N);
end

[coupling_loss_ls, coupling_loss_full] = arrayfun( ...
    @(link) link.calc_coupling_loss(), link_list_ST, 'UniformOutput', false);

coupling_loss_ls = cell2mat(coupling_loss_ls);
coupling_loss_ls = reshape(coupling_loss_ls, links_per_sim, simulation_times).';
selected_link_id = localSelectTrpPairLinks(coupling_loss_ls, trp_pair_idx, best_N, ...
    trp_pair_selection_mode);

selected_link_count = size(selected_link_id, 2);
selected_ls_loss = zeros(simulation_times, selected_link_count);
for ST_num = 1:simulation_times
    selected_ls_loss(ST_num, :) = coupling_loss_ls(ST_num, selected_link_id(ST_num, :));
end

target_result = localResultHeader(scenario, sensing_type, 'TRP-TRP');
target_result.trp_pair_selection_mode = trp_pair_selection_mode;
target_result.best_N = selected_link_count;
if strcmp(char(trp_pair_selection_mode), 'ordered_all')
    target_result.trp_pair_selection_basis = 'all_directed_STX_not_equal_SRX_pairs';
elseif strcmp(char(trp_pair_selection_mode), 'unique_direction_balanced')
    target_result.trp_pair_selection_basis = ...
        'best_N_distinct_physical_pairs_with_balanced_STX_SRX_direction';
else
    target_result.trp_pair_selection_basis = 'largest_STX_ST_SRX_target_channel_power';
end
target_result.LS.CouplingLoss = sort(selected_ls_loss);

if full_cali
    coupling_loss_full = cell2mat(coupling_loss_full);
    coupling_loss_full = reshape(coupling_loss_full, links_per_sim, simulation_times).';

    delay_spread = cell2mat(arrayfun( ...
        @(link) link.calc_delay_spread(), link_list_ST, 'UniformOutput', false));
    delay_spread = reshape(delay_spread, links_per_sim, simulation_times).';

    sa_aod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_path_AOD(1, :)), link_list_ST, 'UniformOutput', false)).';
    sa_zod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_path_ZOD(1, :)), link_list_ST, 'UniformOutput', false)).';
    sa_aoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_path_AOA(2, :)), link_list_ST, 'UniformOutput', false)).';
    sa_zoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_path_ZOA(2, :)), link_list_ST, 'UniformOutput', false)).';
    angular_spread = reshape([sa_aod; sa_zod; sa_aoa; sa_zoa], 4, links_per_sim, simulation_times);
    angular_spread = permute(angular_spread, [1, 3, 2]);

    selected_full_loss = zeros(simulation_times, selected_link_count);
    selected_delay_spread = zeros(simulation_times, selected_link_count);
    selected_angular_spread = zeros(4, simulation_times, selected_link_count);
    for ST_num = 1:simulation_times
        selected_full_loss(ST_num, :) = coupling_loss_full(ST_num, selected_link_id(ST_num, :));
        selected_delay_spread(ST_num, :) = delay_spread(ST_num, selected_link_id(ST_num, :));
        selected_angular_spread(:, ST_num, :) = angular_spread(:, ST_num, selected_link_id(ST_num, :));
    end

    target_result.Full.CouplingLoss = sort(selected_full_loss);
    target_result.Full.DS = selected_delay_spread * 1e9;
    target_result.Full.SA = selected_angular_spread;
end
end

function background_result = localBuildBackgroundCalibrationResult( ...
    sim_data, full_cali, selected_link_id)
scenario = sim_data.scenario;
sensing_type = sim_data.sensing_type;
link_list_background = sim_data.link_list_background;
simulation_times = scenario.simulation_times;
links_per_sim = sim_data.target_links_per_sim(1);
best_N = size(selected_link_id, 2);

if numel(link_list_background) ~= simulation_times * links_per_sim
    error('TRP-TRP background link count does not match the target-link layout.');
end

[coupling_loss_ls, coupling_loss_full] = arrayfun( ...
    @(link) link.calc_coupling_loss(), link_list_background, 'UniformOutput', false);
coupling_loss_ls = reshape(cell2mat(coupling_loss_ls), ...
    links_per_sim, simulation_times).';

selected_ls_loss = zeros(simulation_times, best_N);
for sim_idx = 1:simulation_times
    selected_ls_loss(sim_idx, :) = ...
        coupling_loss_ls(sim_idx, selected_link_id(sim_idx, :));
end

background_result = localResultHeader(scenario, sensing_type, 'TRP-TRP(b)');
background_result.trp_pair_selection_mode = sim_data.trp_pair_selection_mode;
background_result.best_N = best_N;
background_result.trp_pair_selection_basis = ...
    'same_pairs_selected_by_STX_ST_SRX_target_channel_power';
background_result.LS.CouplingLoss = sort(selected_ls_loss);

if full_cali
    coupling_loss_full = reshape(cell2mat(coupling_loss_full), ...
        links_per_sim, simulation_times).';
    delay_spread = reshape(cell2mat(arrayfun( ...
        @(link) link.calc_delay_spread(), link_list_background, ...
        'UniformOutput', false)), links_per_sim, simulation_times).';

    sa_aod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_n_m_AOD, link.phi_LOS_AOD), ...
        link_list_background, 'UniformOutput', false)).';
    sa_zod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_n_m_ZOD, link.theta_LOS_ZOD), ...
        link_list_background, 'UniformOutput', false)).';
    sa_aoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_n_m_AOA, link.phi_LOS_AOA), ...
        link_list_background, 'UniformOutput', false)).';
    sa_zoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_n_m_ZOA, link.theta_LOS_ZOA), ...
        link_list_background, 'UniformOutput', false)).';
    angular_spread = reshape([sa_aod; sa_zod; sa_aoa; sa_zoa], ...
        4, links_per_sim, simulation_times);
    angular_spread = permute(angular_spread, [1, 3, 2]);

    selected_full_loss = zeros(simulation_times, best_N);
    selected_delay_spread = zeros(simulation_times, best_N);
    selected_angular_spread = zeros(4, simulation_times, best_N);
    for sim_idx = 1:simulation_times
        pair_ids = selected_link_id(sim_idx, :);
        selected_full_loss(sim_idx, :) = coupling_loss_full(sim_idx, pair_ids);
        selected_delay_spread(sim_idx, :) = delay_spread(sim_idx, pair_ids);
        selected_angular_spread(:, sim_idx, :) = angular_spread(:, sim_idx, pair_ids);
    end

    background_result.Full.CouplingLoss = sort(selected_full_loss);
    background_result.Full.DS = selected_delay_spread * 1e9;
    background_result.Full.SA = selected_angular_spread;
end
end

function selected_link_id = localSelectTrpPairLinks(pair_selection_score, trp_pair_idx, best_N, ...
    trp_pair_selection_mode)
if ~ismember(char(trp_pair_selection_mode), ...
        {'unique_unordered', 'unique_direction_balanced', 'ordered', 'ordered_all'})
    error('Unsupported TRP pair selection mode: %s.', ...
        string(trp_pair_selection_mode));
end

if size(pair_selection_score, 2) ~= size(trp_pair_idx, 1)
    error('TRP pair score count does not match the generated pair count.');
end

if strcmp(char(trp_pair_selection_mode), 'ordered_all')
    selected_link_id = repmat(1:size(pair_selection_score, 2), ...
        size(pair_selection_score, 1), 1);
elseif strcmp(char(trp_pair_selection_mode), 'unique_direction_balanced')
    selected_link_id = localSelectBalancedPhysicalPairs( ...
        pair_selection_score, trp_pair_idx, best_N);
else
    [~, link_id_sort] = sort(pair_selection_score, 2, 'descend');
    selected_link_id = link_id_sort(:, 1:best_N);
end
end

function selected_link_id = localSelectBalancedPhysicalPairs( ...
    pair_selection_score, trp_pair_idx, best_N)
% Rank reciprocal links as one physical pair, then choose one direction.
canonical_pairs = sort(trp_pair_idx, 2);
[physical_pairs, ~, physical_pair_id] = unique(canonical_pairs, 'rows');
num_physical_pairs = size(physical_pairs, 1);
if num_physical_pairs < best_N
    error('Only %d distinct physical TRP pairs are available, but best_N is %d.', ...
        num_physical_pairs, best_N);
end

forward_link_id = zeros(num_physical_pairs, 1);
reverse_link_id = zeros(num_physical_pairs, 1);
for pair_num = 1:num_physical_pairs
    member_ids = find(physical_pair_id == pair_num);
    if numel(member_ids) ~= 2
        error('Balanced selection requires both directions of every physical TRP pair.');
    end
    pair = physical_pairs(pair_num, :);
    is_forward = trp_pair_idx(member_ids, 1) == pair(1) ...
        & trp_pair_idx(member_ids, 2) == pair(2);
    forward_link_id(pair_num) = member_ids(is_forward);
    reverse_link_id(pair_num) = member_ids(~is_forward);
end

% LS reciprocal links should have equal scores. Averaging makes the
% physical-pair ranking independent of whichever direction is later used.
physical_pair_score = (pair_selection_score(:, forward_link_id) ...
    + pair_selection_score(:, reverse_link_id)) / 2;
[~, physical_pair_rank] = sort(physical_pair_score, 2, 'descend');
selected_physical_pair = physical_pair_rank(:, 1:best_N);

num_drops = size(pair_selection_score, 1);
selected_link_id = zeros(num_drops, best_N);
for drop_idx = 1:num_drops
    for rank_idx = 1:best_N
        pair_num = selected_physical_pair(drop_idx, rank_idx);
        % Alternate direction by rank within each drop and reverse the
        % pattern on the next drop. For Best-N=4 this is exactly 2/2.
        if mod(drop_idx + rank_idx, 2) == 0
            selected_link_id(drop_idx, rank_idx) = forward_link_id(pair_num);
        else
            selected_link_id(drop_idx, rank_idx) = reverse_link_id(pair_num);
        end
    end
end
end

function result = localResultHeader(scenario, sensing_type, link_type)
result = struct();
result.link_type = link_type;
result.fc = scenario.frequency;
result.option = scenario.option;
result.scenario = [scenario.name, scenario.subname];
result.sen_type = [sensing_type.sensing_type, sensing_type.sensing_sub_type];
end
