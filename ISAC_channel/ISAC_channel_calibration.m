clc;
clearvars;
% close all;

%% Calibration controls
case_id = 7;
full_cali = true;
do_rp = true;
current_time = 0;
plot_controller = true;
calibration_root = [];

use_isac_frequency_preset = true;
isac_frequency_preset = 'ISAC_FR1';
custom_frequency_config = struct();

run_all_calibration_cases = true;
calibration_frequency_presets = {'ISAC_FR1', 'ISAC_FR2'};


calibration_cases = localCalibrationCases();

if run_all_calibration_cases
    for case_idx = 1:numel(calibration_cases) %#ok<UNRCH>
        for preset_idx = 1:numel(calibration_frequency_presets)
            active_preset = calibration_frequency_presets{preset_idx};
            localRunCalibrationCase( ...
                calibration_cases(case_idx), case_idx, active_preset, use_isac_frequency_preset, ...
                custom_frequency_config, full_cali, do_rp, current_time, plot_controller, calibration_root);
        end
    end
else
    localRunCalibrationCase( ...
        calibration_cases(case_id), case_id, isac_frequency_preset, use_isac_frequency_preset, ...
        custom_frequency_config, full_cali, do_rp, current_time, plot_controller, calibration_root);
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
    custom_frequency_config, full_cali, do_rp, current_time, plot_controller, calibration_root)

fprintf('Run calibration case %d: %s, preset=%s\n', ...
    case_id, case_config.name, string(isac_frequency_preset));

[scenario, sensing_type] = localCreateCalibrationObjects(case_config);
if use_isac_frequency_preset
    scenario.applyIsacFrequencyPreset(isac_frequency_preset, custom_frequency_config);
end
sim_data = localRunCalibrationLayout(scenario, sensing_type, full_cali, do_rp, current_time, plot_controller);

TRP_mono_target_result = localBuildTargetCalibrationResult(sim_data, full_cali);
plot_TRP_mono_result(TRP_mono_target_result, calibration_root);

if do_rp
    TRP_mono_RP_result = localBuildRpCalibrationResult(sim_data, full_cali);
    plot_TRP_mono_result(TRP_mono_RP_result, calibration_root);
end
end

function sim_data = localRunCalibrationLayout(scenario, sensing_type, full_cali, do_rp, current_time, plot_controller)
[BS_list, ~] = network_layout.Drop_BaseStation(scenario);
[UE_list, UE_pos_list, ~] = network_layout.Drop_UE_ISAC(BS_list, scenario, plot_controller);
[ST_list, ~] = network_layout.Drop_ST( ...
    BS_list, UE_pos_list(1:scenario.UE_per_sec, :), scenario, sensing_type, plot_controller);

link_list_UE = [];
link_list_ST = [];
link_list_RP = [];

total_target_links = numel(BS_list) * numel(ST_list);
progress_bar = waitbar(0, 'creating all links...');
progress_guard = onCleanup(@() close(progress_bar));

for BS_num = 1:numel(BS_list)
    STX = BS_list(BS_num);
    SRX = STX;

    for UE_num = 1:numel(UE_list)
        link_UE = channel.Comm_channel(STX, UE_list(UE_num), scenario, full_cali, current_time);
        link_list_UE = [link_list_UE; link_UE]; %#ok<AGROW>
    end

    for ST_num = 1:numel(ST_list)
        ST = ST_list(ST_num);
        link_ST = channel.Target_channel(STX, SRX, ST, scenario, full_cali, current_time);
        link_list_ST = [link_list_ST; link_ST]; %#ok<AGROW>

        progress_ratio = ((BS_num - 1) * numel(ST_list) + ST_num) / total_target_links;
        waitbar(progress_ratio, progress_bar, localProgressText(scenario, progress_ratio));

        if do_rp
            attached_equipment = STX;
            [RP_list, ~] = network_layout.Drop_RP(attached_equipment, scenario);
            attached_equipment.RP = RP_list;

            for rp_num = 1:scenario.RP_per_equipment
                link_RP = channel.Comm_channel(attached_equipment, RP_list(rp_num), scenario, full_cali, current_time);
                link_list_RP = [link_list_RP; link_RP]; %#ok<AGROW>
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
end

function text = localProgressText(scenario, progress_ratio)
percent = floor(progress_ratio * 1000) / 10;
text = sprintf('(%s, %.0fGHz) %.1f %%', scenario.name, scenario.frequency / 1e9, percent);
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
coupling_loss_ls = reshape(coupling_loss_ls, numel(ST_list), []);
[~, link_id_sort] = sort(coupling_loss_ls, 2, 'descend');
selected_link_id = link_id_sort(:, 1:best_N);

selected_ls_loss = zeros(numel(ST_list), best_N);
for ST_num = 1:numel(ST_list)
    selected_ls_loss(ST_num, :) = coupling_loss_ls(ST_num, selected_link_id(ST_num, :));
end

target_result = localResultHeader(scenario, sensing_type, 'TRPmo(t)');
target_result.LS.CouplingLoss = sort(selected_ls_loss);

if full_cali
    coupling_loss_full = cell2mat(coupling_loss_full);
    coupling_loss_full = reshape(coupling_loss_full, numel(ST_list), []);

    delay_spread = cell2mat(arrayfun( ...
        @(link) link.calc_delay_spread(), link_list_ST, 'UniformOutput', false));
    delay_spread = reshape(delay_spread, numel(ST_list), []);

    sa_aod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_path_AOD(1, :)), link_list_ST, 'UniformOutput', false)).';
    sa_zod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_path_ZOD(1, :)), link_list_ST, 'UniformOutput', false)).';
    sa_aoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_path_AOA(2, :)), link_list_ST, 'UniformOutput', false)).';
    sa_zoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_path_ZOA(2, :)), link_list_ST, 'UniformOutput', false)).';
    angular_spread = reshape([sa_aod; sa_zod; sa_aoa; sa_zoa], 4, numel(ST_list), []);

    selected_full_loss = zeros(numel(ST_list), best_N);
    selected_delay_spread = zeros(numel(ST_list), best_N);
    selected_angular_spread = zeros(4, numel(ST_list), best_N);
    for ST_num = 1:numel(ST_list)
        selected_full_loss(ST_num, :) = coupling_loss_full(ST_num, selected_link_id(ST_num, :));
        selected_delay_spread(ST_num, :) = delay_spread(ST_num, selected_link_id(ST_num, :));
        selected_angular_spread(:, ST_num, :) = angular_spread(:, ST_num, selected_link_id(ST_num, :));
    end

    target_result.Full.CouplingLoss = sort(selected_full_loss);
    target_result.Full.DS = selected_delay_spread * 1e9;
    target_result.Full.SA = selected_angular_spread;
end
end

function rp_result = localBuildRpCalibrationResult(sim_data, full_cali)
scenario = sim_data.scenario;
sensing_type = sim_data.sensing_type;
BS_list = sim_data.BS_list;
link_list_RP = sim_data.link_list_RP;

[coupling_loss_ls, coupling_loss_full] = arrayfun( ...
    @(link) link.calc_coupling_loss(), link_list_RP, 'UniformOutput', false);

coupling_loss_ls = cell2mat(coupling_loss_ls);
coupling_loss_ls = reshape(coupling_loss_ls, numel(BS_list), []);

rp_result = localResultHeader(scenario, sensing_type, 'TRPmo(b)');
rp_result.LS.CouplingLoss = sort(coupling_loss_ls);

if full_cali
    coupling_loss_full = cell2mat(coupling_loss_full);
    coupling_loss_full = reshape(coupling_loss_full, scenario.RP_per_equipment, []).';
    combined_full_loss = pow2db(sum(db2pow(coupling_loss_full), 2));

    delay_spread = cell2mat(arrayfun( ...
        @(link) link.calc_delay_spread(), link_list_RP, 'UniformOutput', false));
    delay_spread = reshape(delay_spread, numel(BS_list), []);

    sa_aod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_n_m_AOD, link.phi_LOS_AOD), link_list_RP, 'UniformOutput', false)).';
    sa_zod = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_n_m_ZOD, link.theta_LOS_ZOD), link_list_RP, 'UniformOutput', false)).';
    sa_aoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.phi_n_m_AOA, link.phi_LOS_AOA), link_list_RP, 'UniformOutput', false)).';
    sa_zoa = cell2mat(arrayfun( ...
        @(link) link.calc_angular_spreads(link.theta_n_m_ZOA, link.theta_LOS_ZOA), link_list_RP, 'UniformOutput', false)).';
    angular_spread = reshape([sa_aod; sa_zod; sa_aoa; sa_zoa], 4, numel(BS_list), []);

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
