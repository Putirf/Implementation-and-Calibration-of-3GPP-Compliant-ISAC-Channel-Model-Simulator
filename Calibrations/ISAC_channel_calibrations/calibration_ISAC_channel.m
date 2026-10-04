clc;
clearvars;
clear classes;
% close all;

calibration_dir = fileparts(mfilename('fullpath'));
cd(calibration_dir);
addpath(calibration_dir);
paths = get_ISAC_calibration_paths();
addpath(paths.isacRoot, '-begin');

%% Calibration controls
% 1: UAV-UMa-AV (only for TRP_mono) 
% 2: Human Outdoor-UMa
% 3: Human Outdoor-UMi
% 4: Human Indoor-InH
% 5: Human Indoor-InF-SH
% 6: AGV-InF-SH
% 7: Automotive-Urban Grid

case_id = 7;
full_cali = true;
do_background = false;  % false: target only; true: also run RP background links
current_time = 0;
plot_controller = true;
calibration_root = paths.referenceRoot;

sensing_mode = 'UE_mono';    % sensing_mode options: 'TRP_mono', 'UE_mono'
sensing_mode = validatestring(sensing_mode, {'TRP_mono','UE_mono'});

use_isac_frequency_preset = true;
isac_frequency_preset = 'ISAC_FR1'; % 'ISAC_FR1' or 'ISAC_FR2'
custom_frequency_config = struct();

run_all_calibration_cases = false;
run_all_frequency_presets = false; 
calibration_frequency_presets = {'ISAC_FR1', 'ISAC_FR2'};


calibration_cases = localCalibrationCases();
if run_all_frequency_presets
    frequency_preset_list = calibration_frequency_presets; %#ok<UNRCH>
else
    frequency_preset_list = {isac_frequency_preset};
end

if run_all_calibration_cases
    for case_idx = 1:numel(calibration_cases) %#ok<UNRCH>
        target_results = cell(1, numel(frequency_preset_list));
        rp_results = cell(1, numel(frequency_preset_list));
        for preset_idx = 1:numel(frequency_preset_list)
            active_preset = frequency_preset_list{preset_idx};
            [target_results{preset_idx}, rp_results{preset_idx}] = localRunCalibrationCase( ...
                calibration_cases(case_idx), case_idx, active_preset, use_isac_frequency_preset, ...
                custom_frequency_config, full_cali, do_background, current_time, ...
                plot_controller, sensing_mode);
        end
        plot_TRP_mono_result([target_results{:}], calibration_root);
        if do_background
            plot_TRP_mono_result([rp_results{:}], calibration_root);
        end
    end
else
    target_results = cell(1, numel(frequency_preset_list));
    rp_results = cell(1, numel(frequency_preset_list));
    for preset_idx = 1:numel(frequency_preset_list)
        active_preset = frequency_preset_list{preset_idx};
        [target_results{preset_idx}, rp_results{preset_idx}] = localRunCalibrationCase( ...
            calibration_cases(case_id), case_id, active_preset, use_isac_frequency_preset, ...
            custom_frequency_config, full_cali, do_background, current_time, ...
            plot_controller, sensing_mode);
    end
    plot_TRP_mono_result([target_results{:}], calibration_root);
    if do_background
        plot_TRP_mono_result([rp_results{:}], calibration_root); %#ok<UNRCH>
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

function [target_result, rp_result] = localRunCalibrationCase( ...
    case_config, case_id, isac_frequency_preset, use_isac_frequency_preset, ...
    custom_frequency_config, full_cali, do_background, current_time, ...
    plot_controller, sensing_mode)

fprintf('Run calibration case %d: %s, mode=%s, preset=%s\n', ...
    case_id, case_config.name, sensing_mode, string(isac_frequency_preset));

[scenario, sensing_type] = localCreateCalibrationObjects(case_config);
if strcmp(sensing_mode, 'UE_mono')
    localValidateUeMonostaticCase(scenario, sensing_type);
end
if strcmp(sensing_mode, 'UE_mono') && ...
        strcmp(sensing_type.sensing_type, 'Human') && ...
        ismember(scenario.name, {'UMa', 'UMi'})
    % Use ten sensing UEs per cell so Best-4 UE-monostatic selection has
    % the UMa/UMi calibration candidate population specified by 3GPP.
    % The mixed indoor/outdoor drop and O2I handling remain enabled.
    scenario.UE_per_sec = 10;
end
if strcmp(sensing_type.sensing_type, 'Human') && ...
        ismember(scenario.name, {'InH', 'InF'})
    % TR 38.901 Table 7.9.6.1-2: indoor Human calibration has
    % 20 UTs in total at 1 m for both TRP- and UE-monostatic modes,
    % not one UT per BS/sector.
    scenario.UE_height = 1;
end
if use_isac_frequency_preset
    scenario.applyIsacFrequencyPreset(isac_frequency_preset, custom_frequency_config);
end
if strcmp(sensing_mode, 'UE_mono') && ...
        strcmp(sensing_type.sensing_type, 'Vehicle')
    % Tables 7.9.6.1-3/7.9.6.2-3: the bundled UE-monostatic
    % Automotive references use one pedestrian-UT calibration at a time.
    scenario.user_case = "Pedestrian";
    scenario.best_N = 4;
    localValidateAutomotiveUrbanGridScenario(scenario);
end
if full_cali
    % TR 38.901 Table 7.9.6.2: full calibration uses -40 dB path
    % dropping after concatenation, independently of coupling Option 1/2.
    scenario.target_path_drop_threshold_db = -40;
end
sim_data = localRunCalibrationLayout(scenario, sensing_type, full_cali, ...
    do_background, current_time, plot_controller, sensing_mode);

target_result = localBuildTargetCalibrationResult(sim_data, full_cali);

if do_background
    rp_result = localBuildRpCalibrationResult(sim_data, full_cali);
else
    rp_result = struct([]);
end
end

function sim_data = localRunCalibrationLayout(scenario, sensing_type, full_cali, ...
        do_background, current_time, plot_controller, sensing_mode)
[BS_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller);
if localIsHumanIndoorCalibration(scenario, sensing_type)
    % Tables 7.9.6.1-1/2: parameters not repeated for Human calibration
    % inherit the single dual-polarized isotropic antenna assumption.
    BS_list = tools.configureCalibrationAntennas( ...
        BS_list, [45, -45], 'TRP');
elseif localIsAgvInFCalibration(scenario, sensing_type)
    % Table 7.9.6.1-4 inherits the single dual-polarized isotropic antenna
    % assumption from Table 7.9.6.1-1.
    BS_list = tools.configureCalibrationAntennas( ...
        BS_list, [45, -45], 'AGV TRP');
    BS_list = localConfigureAgvCalibrationTrps(BS_list);
end
[ue_layouts, UE_list, ~, ST_list, ST_pos_list] = ...
    localDropCalibrationLayouts( ...
    BS_list, scenario, sensing_type, sensing_mode, plot_controller);
if strcmp(sensing_mode, 'UE_mono')
    ue_count_per_drop = cellfun(@numel, ue_layouts);
    if any(ue_count_per_drop ~= ue_count_per_drop(1))
        error('ISACCalibration:InconsistentUeCountPerDrop', ...
            'Every UE_mono drop must contain the same number of sensing UEs.');
    end
    target_links_per_drop = ue_count_per_drop(1);
    sensing_node_list = UE_list;
else
    target_links_per_drop = [];
    sensing_node_list = BS_list;
end
if plot_controller && ~isempty(ST_pos_list)
    localPlotFirstCalibrationTarget(ST_pos_list(1, :));
end

link_list_UE = [];
link_list_ST = [];
link_list_RP = [];

if strcmp(sensing_mode, 'UE_mono')
    total_target_links = target_links_per_drop * numel(ST_list);
else
    total_target_links = numel(sensing_node_list) * numel(ST_list);
end
progress_bar = waitbar(0, 'creating all links...');
progress_guard = onCleanup(@() close(progress_bar));

if do_background
    total_background_links = numel(BS_list) * numel(UE_list);
    completed_background_links = 0;
    for BS_num = 1:numel(BS_list)
        communication_tx = BS_list(BS_num);
        for UE_num = 1:numel(UE_list)
            link_UE = channel.Comm_channel(communication_tx, UE_list(UE_num), ...
                scenario, full_cali, current_time);
            link_list_UE = [link_list_UE; link_UE]; %#ok<AGROW>

            completed_background_links = completed_background_links + 1;
            background_progress = completed_background_links / ...
                total_background_links;
            waitbar(background_progress, progress_bar, sprintf( ...
                '(%s, %.0fGHz) BS-UE background %.1f %%', ...
                scenario.name, scenario.frequency / 1e9, ...
                100 * background_progress));
            drawnow limitrate;
        end
    end
end

if strcmp(sensing_mode, 'UE_mono')
    completed_target_links = 0;
    for ST_num = 1:numel(ST_list)
        ST = ST_list(ST_num);
        current_ue_layout = ue_layouts{ST_num};
        for node_num = 1:numel(current_ue_layout)
            STX = current_ue_layout(node_num);
            link_ST = channel.Target_channel( ...
                STX, STX, ST, scenario, full_cali, current_time);
            link_list_ST = [link_list_ST; link_ST]; %#ok<AGROW>

            completed_target_links = completed_target_links + 1;
            progress_ratio = completed_target_links / total_target_links;
            waitbar(progress_ratio, progress_bar, ...
                localProgressText(scenario, progress_ratio));
            drawnow limitrate;

            if do_background
                attached_equipment = STX;
                [RP_list, ~] = network_layout.Drop_RP(attached_equipment, scenario);
                attached_equipment.RP = RP_list;
                for rp_num = 1:scenario.RP_per_equipment
                    link_RP = channel.Comm_channel(attached_equipment, ...
                        RP_list(rp_num), scenario, full_cali, current_time);
                    link_list_RP = [link_list_RP; link_RP]; %#ok<AGROW>
                end
            end
        end
    end
else
    completed_target_links = 0;
    for node_num = 1:numel(sensing_node_list)
        STX = sensing_node_list(node_num);
        for ST_num = 1:numel(ST_list)
            ST = ST_list(ST_num);
            link_ST = channel.Target_channel( ...
                STX, STX, ST, scenario, full_cali, current_time);
            link_list_ST = [link_list_ST; link_ST]; %#ok<AGROW>

            completed_target_links = completed_target_links + 1;
            progress_ratio = completed_target_links / total_target_links;
            waitbar(progress_ratio, progress_bar, ...
                localProgressText(scenario, progress_ratio));
            drawnow limitrate;

            if do_background
                attached_equipment = STX;
                [RP_list, ~] = network_layout.Drop_RP(attached_equipment, scenario);
                attached_equipment.RP = RP_list;
                for rp_num = 1:scenario.RP_per_equipment
                    link_RP = channel.Comm_channel(attached_equipment, ...
                        RP_list(rp_num), scenario, full_cali, current_time);
                    link_list_RP = [link_list_RP; link_RP]; %#ok<AGROW>
                end
            end
        end
    end
end

sim_data = struct();
sim_data.scenario = scenario;
sim_data.sensing_type = sensing_type;
sim_data.BS_list = BS_list;
sim_data.UE_list = UE_list;
sim_data.ST_list = ST_list;
sim_data.link_list_UE = link_list_UE;
sim_data.link_list_ST = link_list_ST;
sim_data.link_list_RP = link_list_RP;
sim_data.sensing_mode = sensing_mode;
sim_data.sensing_node_list = sensing_node_list;
sim_data.target_links_per_drop = target_links_per_drop;
end

function [ue_layouts, UE_list, UE_pos_list, ST_list, ST_pos_list] = ...
        localDropCalibrationLayouts( ...
        BS_list, scenario, sensing_type, sensing_mode, plot_controller)
drop_count = scenario.ST_per_cell;
ue_layouts = cell(drop_count, 1);
UE_list = [];
UE_pos_list = [];
ST_list = [];
ST_pos_list = [];

original_st_per_cell = scenario.ST_per_cell;
restore_st_count = onCleanup(@() localRestoreTargetDropCount( ...
    scenario, original_st_per_cell));
scenario.ST_per_cell = 1;

for drop_idx = 1:drop_count
    plot_this_layout = plot_controller && drop_idx == 1;
    [drop_ue_list, drop_ue_positions] = localDropCalibrationUes( ...
        BS_list, scenario, sensing_type, sensing_mode, plot_this_layout);
    if localIsHumanCalibration(sensing_type) ...
            || localIsAgvInFCalibration(scenario, sensing_type)
        % Tables 7.9.6.1-1/2: Human calibration UTs use one
        % dual-polarized isotropic antenna for both indoor and outdoor
        % scenarios.  The UT polarization pair is 0/90 degrees.
        drop_ue_list = tools.configureCalibrationAntennas( ...
            drop_ue_list, [90, 0], 'UT');
    end
    st_exclusion_ue_positions = localTargetExclusionUePositions( ...
        drop_ue_positions, scenario, sensing_mode);
    if strcmp(sensing_mode, 'UE_mono')
        % UE-monostatic calibration has no BS-ST exclusion, but preserves
        % the scenario's UE-ST minimum distance for the current UE layout.
        [drop_st, drop_st_position] = network_layout.Drop_ST( ...
            BS_list, st_exclusion_ue_positions, scenario, ...
            sensing_type, false, [], 0);
    else
        [drop_st, drop_st_position] = network_layout.Drop_ST( ...
            BS_list, st_exclusion_ue_positions, scenario, ...
            sensing_type, false);
    end
    drop_st.ID = drop_idx;

    ue_layouts{drop_idx} = drop_ue_list;
    UE_list = [UE_list; drop_ue_list]; %#ok<AGROW>
    UE_pos_list = [UE_pos_list; drop_ue_positions]; %#ok<AGROW>
    ST_list = [ST_list; drop_st]; %#ok<AGROW>
    ST_pos_list = [ST_pos_list; drop_st_position]; %#ok<AGROW>
end
clear restore_st_count;
end

function tf = localIsHumanIndoorCalibration(scenario, sensing_type)
tf = ismember(scenario.name, {'InH', 'InF'}) ...
    && strcmp(sensing_type.sensing_type, 'Human');
end

function tf = localIsHumanCalibration(sensing_type)
tf = strcmp(sensing_type.sensing_type, 'Human');
end

function tf = localIsAgvInFCalibration(scenario, sensing_type)
tf = strcmp(scenario.name, 'InF') ...
    && strcmp(sensing_type.sensing_type, 'AGV');
end

function BS_list = localConfigureAgvCalibrationTrps(BS_list)
% Table 7.9.6.1-4: local X-axis and antenna boresight point downward.
for bs_idx = 1:numel(BS_list)
    for sector_idx = 1:numel(BS_list(bs_idx).sector)
        BS_list(bs_idx).sector(sector_idx).antenna.beta = 90;
    end
end
end

function [UE_list, UE_pos_list] = localDropCalibrationUes( ...
        BS_list, scenario, sensing_type, sensing_mode, plot_controller)
if strcmp(sensing_type.sensing_type, 'Human') && ...
        ismember(scenario.name, {'InH', 'InF'})
    % Table 7.9.6.1-2 specifies 20 UTs in total for indoor Human
    % calibration in TRP- and UE-monostatic modes, uniformly distributed
    % over the indoor scenario area.
    indoor_ue_count = 20;
    ue_xy_range = [scenario.x_range; scenario.y_range];
    indoor_ue_positions = zeros(indoor_ue_count, 2);
    for ue_idx = 1:indoor_ue_count
        indoor_ue_positions(ue_idx, :) = ...
            network_layout.drop_InH_ISAC(ue_xy_range);
    end
    [UE_list, UE_pos_list, ~] = network_layout.Drop_UE_ISAC( ...
        BS_list, scenario, plot_controller, indoor_ue_positions, 'indoor');
elseif strcmp(sensing_mode, 'UE_mono') && ...
        strcmp(sensing_type.sensing_type, 'AGV') && ...
        strcmp(scenario.name, 'InF')
    % Tables 7.9.6.1-4 and 7.8-7: 30 indoor UTs at 1.5 m,
    % uniformly dropped with a 1 m minimum 2D distance to every TRP.
    agv_ue_positions = localDropAgvCalibrationUts(BS_list, scenario, 30);
    [UE_list, UE_pos_list, ~] = network_layout.Drop_UE_ISAC( ...
        BS_list, scenario, plot_controller, agv_ue_positions, 'indoor');
else
    [UE_list, UE_pos_list, ~] = network_layout.Drop_UE_ISAC( ...
        BS_list, scenario, plot_controller, [], ...
        localUeEnvironment(scenario, sensing_type, sensing_mode));
end
end

function localValidateUeMonostaticCase(scenario, sensing_type)
target_type = char(sensing_type.sensing_type);
is_supported = strcmp(target_type, 'Human') ...
    || (strcmp(target_type, 'AGV') && strcmp(scenario.name, 'InF') ...
        && strcmp(scenario.subname, 'SH')) ...
    || (strcmp(target_type, 'Vehicle') && strcmp(scenario.name, 'UrbanGrid'));
if ~is_supported
    error('ISACCalibration:UnsupportedUeMonostaticTarget', ...
        'Unsupported UE_mono calibration: target=%s, scenario=%s-%s.', ...
        target_type, scenario.name, scenario.subname);
end
end

function ue_positions = localDropAgvCalibrationUts(BS_list, scenario, ue_count)
trp_positions = reshape([BS_list.Position], 3, []).';
ue_positions = zeros(ue_count, 3);
for ue_idx = 1:ue_count
    while true
        candidate_xy = [ ...
            scenario.x_range(1) + diff(scenario.x_range) * rand, ...
            scenario.y_range(1) + diff(scenario.y_range) * rand];
        if all(vecnorm(trp_positions(:, 1:2) - candidate_xy, 2, 2) >= 1)
            ue_positions(ue_idx, :) = [candidate_xy, 1.5];
            break;
        end
    end
end
end

function localValidateAutomotiveUrbanGridScenario(scenario)
if scenario.BS_height ~= 25
    error('Automotive Urban Grid UE_mono requires 25 m TRP height.');
end
if scenario.frequency == 6e9
    expected_isd = 500;
elseif scenario.frequency == 30e9
    expected_isd = 250;
else
    error('Automotive Urban Grid UE_mono supports only 6 GHz and 30 GHz.');
end
if scenario.ISD ~= expected_isd
    error('Automotive Urban Grid at %.0f GHz requires ISD=%d m.', ...
        scenario.frequency/1e9, expected_isd);
end
end

function localRestoreTargetDropCount(scenario, original_st_per_cell)
scenario.ST_per_cell = original_st_per_cell;
end

function text = localProgressText(scenario, progress_ratio)
percent = floor(progress_ratio * 1000) / 10;
text = sprintf('(%s, %.0fGHz) %.1f %%', scenario.name, scenario.frequency / 1e9, percent);
end

function localPlotFirstCalibrationTarget(st_position)
% Match the bistatic calibration view: show only the first Monte Carlo ST.
figure(1);
plot3(st_position(1), st_position(2), st_position(3), '*', ...
    'Color', 'g', 'MarkerSize', 8, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
axis equal;
view(0, 90);
end

function ue_environment = localUeEnvironment(scenario, sensing_type, sensing_mode)
if strcmp(sensing_mode, 'UE_mono') ...
        && ismember(scenario.name, {'UMa', 'UMi'}) ...
        && strcmp(sensing_type.sensing_type, 'Human')
    ue_environment = 'outdoor';
else
    ue_environment = 'mixed';
end
end

function ue_positions = localTargetExclusionUePositions(UE_pos_list, scenario, sensing_mode)
if strcmp(sensing_mode, 'UE_mono')
    % In UE-monostatic calibration every UE is a sensing endpoint. Enforce
    % the scenario's UE-ST minimum distance against every sensing UE.
    ue_positions = UE_pos_list;
else
    % Preserve the existing target-drop behaviour of all other calibration
    % modes, including TRP_mono.
    ue_positions = UE_pos_list(1:scenario.UE_per_sec, :);
end
end

function target_result = localBuildTargetCalibrationResult(sim_data, full_cali)
scenario = sim_data.scenario;
sensing_type = sim_data.sensing_type;
ST_list = sim_data.ST_list;
link_list_ST = sim_data.link_list_ST;
best_N = scenario.best_N;

[coupling_loss_ls, coupling_loss_full] = arrayfun( ...
    @(link) link.calc_coupling_loss(), link_list_ST, 'UniformOutput', false);

coupling_loss_ls = cell2mat(coupling_loss_ls);
coupling_loss_ls = localReshapeTargetLinks(coupling_loss_ls, sim_data);
[~, link_id_sort] = sort(coupling_loss_ls, 2, 'descend');
selected_link_id = link_id_sort(:, 1:best_N);

selected_ls_loss = zeros(numel(ST_list), best_N);
for ST_num = 1:numel(ST_list)
    selected_ls_loss(ST_num, :) = coupling_loss_ls(ST_num, selected_link_id(ST_num, :));
end

if strcmp(sim_data.sensing_mode, 'UE_mono')
    link_type = 'UEmo(t)';
else
    link_type = 'TRPmo(t)';
end
target_result = localResultHeader(scenario, sensing_type, link_type);
target_result.LS.CouplingLoss = sort(selected_ls_loss);

if full_cali
    coupling_loss_full = cell2mat(coupling_loss_full);
    coupling_loss_full = localReshapeTargetLinks(coupling_loss_full, sim_data);
    % Tables 7.9.6.1-1 and 7.9.6.2-1/2: Best-N STX/SRX selection is based
    % on the target-channel large-scale power scaling factor (PL, SF and
    % sigma_M).  Full coupling loss, DS and angular spreads must therefore
    % use the same links selected above, without re-ranking on instantaneous
    % small-scale RCS, polarization or fast fading.
    selected_full_link_id = selected_link_id;

    delay_spread = cell2mat(arrayfun( ...
        @(link) link.calc_delay_spread(), link_list_ST, 'UniformOutput', false));
    delay_spread = localReshapeTargetLinks(delay_spread, sim_data);

    sa_aod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_path_AOD(1, :)), link_list_ST, 'UniformOutput', false)).';
    sa_zod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_path_ZOD(1, :)), link_list_ST, 'UniformOutput', false)).';
    sa_aoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_path_AOA(2, :)), link_list_ST, 'UniformOutput', false)).';
    sa_zoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_path_ZOA(2, :)), link_list_ST, 'UniformOutput', false)).';
    sa_aod = localReshapeTargetLinks(sa_aod, sim_data);
    sa_zod = localReshapeTargetLinks(sa_zod, sim_data);
    sa_aoa = localReshapeTargetLinks(sa_aoa, sim_data);
    sa_zoa = localReshapeTargetLinks(sa_zoa, sim_data);
    angular_spread = cat(1, ...
        reshape(sa_aod, 1, size(sa_aod, 1), size(sa_aod, 2)), ...
        reshape(sa_zod, 1, size(sa_zod, 1), size(sa_zod, 2)), ...
        reshape(sa_aoa, 1, size(sa_aoa, 1), size(sa_aoa, 2)), ...
        reshape(sa_zoa, 1, size(sa_zoa, 1), size(sa_zoa, 2)));

    selected_full_loss = zeros(numel(ST_list), best_N);
    selected_delay_spread = zeros(numel(ST_list), best_N);
    selected_angular_spread = zeros(4, numel(ST_list), best_N);
    for ST_num = 1:numel(ST_list)
        full_ids = selected_full_link_id(ST_num, :);
        selected_full_loss(ST_num, :) = coupling_loss_full(ST_num, full_ids);
        selected_delay_spread(ST_num, :) = delay_spread(ST_num, full_ids);
        selected_angular_spread(:, ST_num, :) = angular_spread(:, ST_num, full_ids);
    end

    target_result.Full.CouplingLoss = sort(selected_full_loss);
    target_result.Full.DS = selected_delay_spread * 1e9;
    target_result.Full.SA = selected_angular_spread;
end
end

function values_by_drop = localReshapeTargetLinks(values, sim_data)
values = values(:);
target_count = numel(sim_data.ST_list);
if strcmp(sim_data.sensing_mode, 'UE_mono')
    links_per_drop = sim_data.target_links_per_drop;
    expected_count = target_count * links_per_drop;
    if numel(values) ~= expected_count
        error('ISACCalibration:UnexpectedTargetLinkCount', ...
            'Expected %d target links but received %d.', ...
            expected_count, numel(values));
    end
    values_by_drop = reshape(values, links_per_drop, target_count).';
else
    values_by_drop = reshape(values, target_count, []);
end
end

function rp_result = localBuildRpCalibrationResult(sim_data, full_cali)
scenario = sim_data.scenario;
sensing_type = sim_data.sensing_type;
sensing_node_list = sim_data.sensing_node_list;
link_list_RP = sim_data.link_list_RP;

[coupling_loss_ls, coupling_loss_full] = arrayfun( ...
    @(link) link.calc_coupling_loss(), link_list_RP, 'UniformOutput', false);

coupling_loss_ls = cell2mat(coupling_loss_ls);
coupling_loss_ls = reshape(coupling_loss_ls, numel(sensing_node_list), []);

if strcmp(sim_data.sensing_mode, 'UE_mono')
    link_type = 'UEmo(b)';
else
    link_type = 'TRPmo(b)';
end
rp_result = localResultHeader(scenario, sensing_type, link_type);
rp_result.LS.CouplingLoss = sort(coupling_loss_ls);

if full_cali
    coupling_loss_full = cell2mat(coupling_loss_full);
    coupling_loss_full = reshape(coupling_loss_full, scenario.RP_per_equipment, []).';
    combined_full_loss = pow2db(sum(db2pow(coupling_loss_full), 2));

    delay_spread = cell2mat(arrayfun( ...
        @(link) link.calc_delay_spread(), link_list_RP, 'UniformOutput', false));
    delay_spread = reshape(delay_spread, numel(sensing_node_list), []);

    sa_aod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_n_m_AOD, link.phi_LOS_AOD), link_list_RP, 'UniformOutput', false)).';
    sa_zod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_n_m_ZOD, link.theta_LOS_ZOD), link_list_RP, 'UniformOutput', false)).';
    sa_aoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_n_m_AOA, link.phi_LOS_AOA), link_list_RP, 'UniformOutput', false)).';
    sa_zoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_n_m_ZOA, link.theta_LOS_ZOA), link_list_RP, 'UniformOutput', false)).';
    angular_spread = reshape([sa_aod; sa_zod; sa_aoa; sa_zoa], ...
        4, numel(sensing_node_list), []);

    rp_result.Full.CouplingLoss = sort(combined_full_loss);
    rp_result.Full.DS = delay_spread * 1e9;
    rp_result.Full.SA = angular_spread;
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
