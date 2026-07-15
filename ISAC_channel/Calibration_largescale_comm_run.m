function Calibration_largescale_comm_run(selected_scenario)
if nargin < 1 || strlength(string(selected_scenario)) == 0
    selected_scenario = 'UMi';
end

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(fileparts(script_dir));
addpath(project_dir);
original_dir = pwd;
cleanup_dir = onCleanup(@() cd(original_dir)); %#ok<NASGU>
cd(script_dir);

plot_controller = false;
current_time = 0;
enable_full_channel = false;
result_root = fullfile(script_dir, 'results', 'comm_large_scale');
if ~exist(result_root, 'dir')
    mkdir(result_root);
end

% Build the selectable scenario list so the same runner can be reused for UMa/UMi/InH.
scenario_configs = {
    localScenarioConfig('UMa', @() comm_scenario.UMa(), 'UMa')
    localScenarioConfig('UMi', @() comm_scenario.UMi(), 'UMi')
    localScenarioConfig('InH', @() localCreateInHScenario(), 'InH')
};
scenario_config = localSelectScenarioConfig(scenario_configs, selected_scenario);

% Large-scale calibration uses the communication-only 6/30/70 GHz comparison set.
frequency_configs = {
    localFrequencyConfig(6e9, 20e6, '6GHz', 'ISAC_FR1')
    localFrequencyConfig(30e9, 100e6, '30GHz', 'ISAC_FR2')
    localFrequencyConfig(70e9, 100e6, '70GHz', 'CUSTOM')
};

% Run one calibration case per carrier frequency and store the sorted metrics.
for freq_idx = 1:numel(frequency_configs)
    fprintf('Running %s %s...\n', ...
        scenario_config.Name, frequency_configs{freq_idx}.Tag);
    result = localRunPureCommCase( ...
        scenario_config, frequency_configs{freq_idx}, ...
        current_time, enable_full_channel, plot_controller);

    save(fullfile(result_root, ...
        sprintf('%s_%s.mat', result.ScenarioLabel, result.FrequencyTag)), 'result');
end
end

function config = localScenarioConfig(name, factory, label)
config = struct();
config.Name = name;
config.Factory = factory;
config.Label = label;
end

function scenario_config = localSelectScenarioConfig(all_configs, selected_scenario)
selected_key = lower(string(selected_scenario));
for idx = 1:numel(all_configs)
    config = all_configs{idx};
    if lower(string(config.Name)) == selected_key || lower(string(config.Label)) == selected_key
        scenario_config = config;
        return;
    end
end

error('Unsupported selected_scenario: %s', string(selected_scenario));
end

function scenario = localCreateInHScenario()
scenario = comm_scenario.InH();
scenario.subname = 'open_office';
end

function config = localFrequencyConfig(frequency_hz, bandwidth_hz, tag, source_mode)
config = struct();
config.FrequencyHz = frequency_hz;
config.BandwidthHz = bandwidth_hz;
config.Tag = tag;
config.SourceMode = source_mode;
end

function result = localRunPureCommCase( ...
    scenario_config, frequency_config, current_time, enable_full_channel, plot_controller)

% Instantiate the scenario and overwrite only the calibration-specific assumptions.
scenario = scenario_config.Factory();
scenario = localApplyLargeScaleCommSettings(scenario, frequency_config);

% Rebuild BS/UE objects with the calibration antenna settings before collecting metrics.
[BS_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller);
[BS_list, ~] = localApplyBsCalibrationConfiguration(BS_list, scenario);
[UE_list, ~] = localDropCalibrationUEs(BS_list, scenario, plot_controller);
close all;

noise_mw = 10^((10 * log10(db2pow(-174) * scenario.BW) + 9) / 10);

ue_count = numel(UE_list);
coupling_loss = nan(ue_count, 1);
geometry_sir = nan(ue_count, 1);
geometry_sinr = nan(ue_count, 1);

progress_bar = waitbar(0, sprintf('Running %s %s...', ...
    scenario.name, frequency_config.Tag));
progress_guard = onCleanup(@() close(progress_bar)); %#ok<NASGU>

% For each UE, compare all candidate sectors and keep the serving-cell large-scale metrics.
for ue_idx = 1:ue_count
    current_ue = UE_list(ue_idx);
    sector_metric = [];
    sector_rsrp = [];

    for bs_idx = 1:numel(BS_list)
        current_bs = BS_list(bs_idx);
        link = channel.Comm_channel(current_bs, current_ue, scenario, enable_full_channel, current_time);
        base_coupling_loss = localCalibrationCouplingLossDb(link);

        for sector_idx = 1:numel(current_bs.sector)
            sector = current_bs.sector(sector_idx);
            port_gain = localSectorVirtualizedPortGain(sector, link);
            sector_metric(end + 1, 1) = base_coupling_loss + port_gain; %#ok<AGROW>
            sector_rsrp(end + 1, 1) = sector_metric(end) + current_bs.Power; %#ok<AGROW>
        end
    end

    [coupling_loss(ue_idx), serving_sector_idx] = max(sector_metric);
    serving_rsrp = sector_rsrp(serving_sector_idx);
    interference_rsrp = sector_rsrp;
    interference_rsrp(serving_sector_idx) = -inf;

    interference_linear = sum(10.^(interference_rsrp / 10));
    geometry_sir(ue_idx) = serving_rsrp - 10 * log10(interference_linear);
    geometry_sinr(ue_idx) = serving_rsrp - 10 * log10(interference_linear + noise_mw);

    waitbar(ue_idx / ue_count, progress_bar, sprintf('(%s, %s) %.1f %%', ...
        scenario.name, frequency_config.Tag, ue_idx / ue_count * 100));
end

result = struct();
result.ScenarioName = scenario.name;
result.ScenarioLabel = scenario_config.Label;
result.FrequencyTag = frequency_config.Tag;
result.FrequencyHz = frequency_config.FrequencyHz;
result.BandwidthHz = frequency_config.BandwidthHz;
result.FrequencySource = frequency_config.SourceMode;
result.CouplingLoss = sort(coupling_loss);
result.GeometrySIR = sort(geometry_sir);
result.GeometrySINR = sort(geometry_sinr);
result.Pr = (1:ue_count).' / ue_count * 100;
end

function coupling_loss_db = localCalibrationCouplingLossDb(link)
coupling_loss_db = link.calc_coupling_loss();
end

function pathloss_db = localCalibrationLosPathlossDb(link)
tx_height = link.TX.height;
rx_height = link.RX.height;
fc_ghz = link.fc / 1e9;

switch link.scenario.name
    case 'UMa'
        rx_breakpoint_height = rx_height;
        if rx_height > 22.5
            rx_breakpoint_height = 1.5;
        end

        c_term = 0;
        if rx_breakpoint_height > 13 && rx_breakpoint_height <= 22.5
            if link.d_2D <= 18
                g_term = 0;
            else
                g_term = 5 / 4 * ((link.d_2D / 100) ^ 3) * exp(-link.d_2D / 150);
            end
            c_term = ((rx_breakpoint_height - 13) / 10) ^ 1.5 * g_term;
        end

        if rand < 1 / (1 + c_term) || link.isRP
            environment_height = 1;
        else
            environment_height = 3 * randi([4, max(4, link.RX.n_fl)], 1);
        end

        effective_tx_height = tx_height - environment_height;
        effective_rx_height = rx_breakpoint_height - environment_height;
        breakpoint_distance = 4 * effective_tx_height * effective_rx_height * link.fc / 3e8;

        if link.d_2D >= 10 && link.d_2D <= breakpoint_distance
            los_pathloss_db = 22 * log10(link.d_3D) + 28 + 20 * log10(fc_ghz);
        elseif link.d_2D > breakpoint_distance && link.d_2D < 5000
            los_pathloss_db = 40 * log10(link.d_3D) + 28 + 20 * log10(fc_ghz) ...
                - 9 * log10(breakpoint_distance ^ 2 + (tx_height - rx_height) ^ 2);
        else
            los_pathloss_db = 22 * log10(link.d_3D) + 28 + 20 * log10(fc_ghz);
        end

        pathloss_db = los_pathloss_db;
        if link.O2I
            pathloss_db = pathloss_db + localCalibrationO2ILossDb(link);
        end

    case 'UMi'
        effective_tx_height = tx_height - 1;
        effective_rx_height = rx_height - 1;
        breakpoint_distance = 4 * effective_tx_height * effective_rx_height * link.fc / 3e8;

        if link.d_2D <= breakpoint_distance
            los_pathloss_db = 32.4 + 21 * log10(link.d_3D) + 20 * log10(fc_ghz);
        elseif link.d_2D > breakpoint_distance && link.d_2D < 5000
            los_pathloss_db = 32.4 + 40 * log10(link.d_3D) + 20 * log10(fc_ghz) ...
                - 9.5 * log10(breakpoint_distance ^ 2 + (tx_height - rx_height) ^ 2);
        else
            los_pathloss_db = 32.4 + 21 * log10(link.d_3D) + 20 * log10(fc_ghz);
        end

        pathloss_db = los_pathloss_db;
        if link.O2I
            pathloss_db = pathloss_db + localCalibrationO2ILossDb(link);
        end

    case 'InH'
        pathloss_db = 32.4 + 17.3 * log10(link.d_3D) + 20 * log10(fc_ghz);

    otherwise
        error('Unsupported scenario for calibration LOS pathloss: %s', link.scenario.name);
end
end

function loss_db = localCalibrationO2ILossDb(link)
[penetration_loss_db, ~] = localO2IPenetrationLossDb(link.scenario.name, link.fc / 1e9, link.RX.O2IPL);
loss_db = penetration_loss_db + 0.5 * link.d_2D_in;
end

function [penetration_loss_db, penetration_sigma_db] = localO2IPenetrationLossDb(scenario_name, fc_ghz, penetration_level)
if fc_ghz < 6
    penetration_loss_db = 20;
    penetration_sigma_db = 0;
    return;
end

if strcmp(scenario_name, 'RMa')
    penetration_level = 'low';
elseif strcmp(scenario_name, 'InF')
    penetration_level = 'high';
end

material_loss_db = [2 + 0.2 * fc_ghz; 23 + 0.3 * fc_ghz; 5 + 4 * fc_ghz; 4.85 + 0.12 * fc_ghz];
if strcmp(penetration_level, 'low')
    penetration_loss_db = 5 - 10 * log10( ...
        0.3 * 10 ^ (-material_loss_db(1) / 10) + 0.7 * 10 ^ (-material_loss_db(3) / 10));
    penetration_sigma_db = 4.4;
elseif strcmp(penetration_level, 'high')
    penetration_loss_db = 5 - 10 * log10( ...
        0.7 * 10 ^ (-material_loss_db(2) / 10) + 0.3 * 10 ^ (-material_loss_db(3) / 10));
    penetration_sigma_db = 6.5;
else
    error('Unsupported penetration level: %s', string(penetration_level));
end
end

function scenario = localApplyLargeScaleCommSettings(scenario, frequency_config)
% Apply the 7.8-1 large-scale assumptions that differ by scenario/frequency.
scenario.BS_sec_num = 3;
scenario.UE_sec_num = 1;
scenario.UE_per_sec = 50;
scenario.frequency = frequency_config.FrequencyHz;
scenario.BW = frequency_config.BandwidthHz;
scenario.BS_noise_figure = 5;
scenario.UE_noise_figure = 9;

switch scenario.name
    case 'UMa'
        if frequency_config.FrequencyHz == 6e9
            scenario.BS_Tx_power = 49;
        else
            scenario.BS_Tx_power = 35;
        end
    case 'UMi'
        if frequency_config.FrequencyHz == 6e9
            scenario.BS_Tx_power = 44;
        else
            scenario.BS_Tx_power = 35;
        end
    case 'InH'
        scenario.BS_Tx_power = 24;
        scenario.BS_UE_min_d = 0;
end
end

function [bs_list, sector_angles] = localApplyBsCalibrationConfiguration(bs_list, scenario)
% Replace the default BS antenna with the large-scale calibration array and sector angles.
bs_cfg = localBuildBsAntennaConfig(scenario.name);
sector_angles = localSectorAngles();

for bs_idx = 1:numel(bs_list)
    current_bs = bs_list(bs_idx);
    current_bs.antenna_params = bs_cfg;
    current_bs.Power = scenario.BS_Tx_power;
    current_bs.BW = scenario.BW;
    current_bs.frequency = scenario.frequency;
    current_bs.sector_num = scenario.BS_sec_num;

    for sector_idx = 1:numel(current_bs.sector)
        sector = current_bs.sector(sector_idx);
        sector.frequency = scenario.frequency;
        sector.BW = scenario.BW;
        sector.boresight = sector_angles(sector_idx);

        ang = struct('alpha', sector.boresight, 'beta', 0, 'gamma', 0);
        sector.antenna = antennas.antenna_array(bs_cfg, ang);
        sector.antenna.attachedDevice = sector;
        sector.antenna.attachedType = 'BS';
        localApplyElementModels(sector.antenna, 'dipole', bs_cfg.pol_model);
    end
end
end

function [ue_list, ue_pos_list] = localDropCalibrationUEs(bs_list, scenario, plot_controller)
% Drop UEs with the calibration indoor/outdoor rules and attach isotropic UT arrays.
ut_cfg = localBuildUtAntennaConfig();
ue_list = [];
ue_pos_list = [];
ue_id = 0;
lowhigh = {'low', 'high'};
total_bs_sector_num = scenario.total_BS_sector_num;

if ismember(scenario.name, {'UMa', 'UMi'})
    for ue_idx_per_sector = 1:scenario.UE_per_sec
        for bs_idx = 1:numel(bs_list)
            current_bs = bs_list(bs_idx);
            for sector_idx = 1:numel(current_bs.sector)
                current_sector = current_bs.sector(sector_idx);
                current_ue = elements.Equipment(scenario, 'UE');
                ue_id = ue_id + 1;
                current_ue.ID = ue_id;
                current_ue.antenna_params = ut_cfg;
                current_ue.sector_num = 1;
                current_ue.fcin = localRandomizedCenterFrequency(scenario);

                if rand < 0.8
                    current_ue.Indoor = true;
                    current_ue.d_2D_in = scenario.UE_max_d_2D_indoor * rand(1);
                    current_ue.Position = [ ...
                        network_layout.drop_in_hexagonUE( ...
                            current_bs.Position(1:2), scenario.R, ...
                            scenario.BS_UE_min_d + current_ue.d_2D_in, ...
                            current_sector.boresight, current_bs.sector_num), ...
                        0];
                    n_fl_max = randi([4, 8], 1);
                    current_ue.n_fl = randi([1, n_fl_max], 1);
                    current_ue.O2IPL = lowhigh{round(rand) + 1};
                else
                    current_ue.Indoor = false;
                    current_ue.d_2D_in = 0;
                    current_ue.Position = [ ...
                        network_layout.drop_in_hexagonUE( ...
                            current_bs.Position(1:2), scenario.R, ...
                            scenario.BS_UE_min_d, current_sector.boresight, ...
                            current_bs.sector_num), ...
                        0];
                    current_ue.n_fl = 1;
                end

                current_ue.height = 3 * (current_ue.n_fl - 1) + 1.5;
                current_ue.Position(3) = current_ue.height;
                current_ue.inital_Position = current_ue.Position;
                current_ue.rand_LoS = rand(1, numel(bs_list));
                current_ue.O2Isigma = randn;
                current_ue.carPL = [9, 5];
                current_ue = localAttachUtAntenna(current_ue, ut_cfg);

                ue_pos_list = [ue_pos_list; current_ue.Position]; %#ok<AGROW>
                ue_list = [ue_list; current_ue]; %#ok<AGROW>
            end
        end
    end
elseif strcmp(scenario.name, 'InH')
    for bs_idx = 1:numel(bs_list)
        current_bs = bs_list(bs_idx);
        for sector_idx = 1:numel(current_bs.sector)
            for ue_idx_per_sector = 1:scenario.UE_per_sec
                current_ue = elements.Equipment(scenario, 'UE');
                ue_id = ue_id + 1;
                current_ue.ID = ue_id;
                current_ue.antenna_params = ut_cfg;
                current_ue.sector_num = 1;
                current_ue.fcin = localRandomizedCenterFrequency(scenario);
                current_ue.d_2D_in = 0;
                current_ue.n_fl = 1;
                current_ue.height = 1.0;
                current_ue.Position = [network_layout.drop_InH_ISAC([scenario.x_range; scenario.y_range]), current_ue.height];
                current_ue.inital_Position = current_ue.Position;
                current_ue.rand_LoS = rand(1, numel(bs_list));
                current_ue.O2Isigma = 0;
                current_ue.carPL = [0, 0];
                current_ue = localAttachUtAntenna(current_ue, ut_cfg);

                ue_pos_list = [ue_pos_list; current_ue.Position]; %#ok<AGROW>
                ue_list = [ue_list; current_ue]; %#ok<AGROW>
            end
        end
    end
else
    error('Unsupported scenario for calibration UE dropping: %s', scenario.name);
end

if plot_controller && ~isempty(ue_pos_list)
    figure(1);
    plot3(ue_pos_list(:, 1), ue_pos_list(:, 2), ue_pos_list(:, 3), '.', ...
        'Color', [0 0.447058823529412 0.741176470588235], 'HandleVisibility', 'off');
    hold on;
    axis equal;
    view(0, 90);
end
end

function ue = localAttachUtAntenna(ue, ut_cfg)
ue.sector = elements.Sector(ue);
ang = struct('alpha', rand * 360 - 180, 'beta', ut_cfg.beta, 'gamma', 0);
ue.sector.antenna = antennas.antenna_array(ut_cfg, ang);
ue.sector.antenna.attachedDevice = ue;
ue.sector.antenna.attachedType = 'UE';
localApplyElementModels(ue.sector.antenna, 'isotropic', 'model-2');
end

function out = localRandomizedCenterFrequency(scenario)
fc = scenario.frequency;
bw = scenario.BW;
out = (bw / 20) * (floor(20 * rand) - (19 / 2)) + fc;
end

function cfg = localBuildBsAntennaConfig(scenario_name)
cfg = struct();
cfg.array = struct('Mg', 1, 'Ng', 1, 'dg_H', 2.5, 'dg_V', 2.5);
cfg.panel = struct( ...
    'M', 10, ...
    'N', 1, ...
    'Kv', 10, ...
    'Kh', 1, ...
    'd_H', 0.5, ...
    'd_V', 0.5, ...
    'P', 1, ...
    'X_pol', 0, ...
    'ele_downtilt', localBsDowntilt(scenario_name), ...
    'ele_panning', 0);
cfg.pol_model = 'model-2';
end

function cfg = localBuildUtAntennaConfig()
cfg = struct();
cfg.array = struct('Mg', 1, 'Ng', 1, 'dg_H', 2.5, 'dg_V', 2.5);
cfg.panel = struct( ...
    'M', 1, ...
    'N', 1, ...
    'Kv', 1, ...
    'Kh', 1, ...
    'd_H', 0.5, ...
    'd_V', 0.5, ...
    'P', 1, ...
    'X_pol', 0, ...
    'ele_downtilt', 90, ...
    'ele_panning', 0);
cfg.beta = 0;
end

function out = localBsDowntilt(scenario_name)
if strcmp(scenario_name, 'InH')
    out = 110;
else
    out = 102;
end
end

function out = localSectorAngles()
out = [30, 150, 270];
end

function localApplyElementModels(array_obj, antenna_model, pol_model)
if iscell(array_obj.panel)
    panel_cells = array_obj.panel;
else
    panel_cells = {array_obj.panel};
end

for panel_idx = 1:numel(panel_cells)
    panel_obj = panel_cells{panel_idx};
    for element_idx = 1:numel(panel_obj.element_list)
        panel_obj.element_list(element_idx).antenna_model = antenna_model;
        panel_obj.element_list(element_idx).pol_model = pol_model;
    end
end
end

function gain_db = localSectorVirtualizedPortGain(sector, link)
% Read the virtualized CRS-port gain seen by this UE for serving/interference comparison.
gain_db = 0;
if isempty(sector) || isempty(link)
    return;
end

shim_link = localBuildShimLink(sector, link);

if iscell(sector.antenna.panel)
    panel_obj = sector.antenna.panel{1};
else
    panel_obj = sector.antenna.panel(1);
end

gain_values = panel_obj.panel_gain(shim_link);
if ~isempty(gain_values)
    gain_db = gain_values(1);
end
end

function shim_link = localBuildShimLink(sector, link)
shim_link = struct();
shim_link.sector = struct( ...
    'frequency', sector.frequency, ...
    'h_BS', sector.attached_equipment.height);
shim_link.BS_pos_wrap = link.TX_pos_wrap(1:2);
shim_link.theta_LOS_ZOD = link.theta_LOS_ZOD;
shim_link.phi_LOS_AOD = link.phi_LOS_AOD;
shim_link.theta_LOS_ZOA = link.theta_LOS_ZOA;
shim_link.phi_LOS_AOA = link.phi_LOS_AOA;
shim_link.UE = struct( ...
    'pos', link.RX.Position(1:2), ...
    'h_UT', link.RX.height, ...
    'antenna', link.RX.sector.antenna);
end
