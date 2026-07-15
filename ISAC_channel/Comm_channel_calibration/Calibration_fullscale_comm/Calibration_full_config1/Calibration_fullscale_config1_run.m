function Calibration_fullscale_config1_run(selected_scenario)
if nargin < 1 || strlength(string(selected_scenario)) == 0
    selected_scenario = 'UMi';
end

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(fileparts(fileparts(script_dir)));
addpath(project_dir);
original_dir = pwd;
cleanup_dir = onCleanup(@() cd(original_dir)); %#ok<NASGU>
cd(script_dir);

plot_controller = false;
current_time = 0;
enable_full_channel = false;
port0 = 0;
result_root = fullfile(script_dir, 'results', 'comm_full_scale_config1');
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

% Config 1 uses the full-calibration carrier set in 38.901 Table 7.8-2.
frequency_configs = {
    localFrequencyConfig(6e9, 20e6, '6GHz')
    localFrequencyConfig(30e9, 100e6, '30GHz')
    localFrequencyConfig(60e9, 100e6, '60GHz')
    localFrequencyConfig(70e9, 100e6, '70GHz')
};

% Run one calibration case per carrier frequency and store the sorted metrics.
for freq_idx = 1:numel(frequency_configs)
    fprintf('Running %s %s...\n', ...
        scenario_config.Name, frequency_configs{freq_idx}.Tag);
    result = localRunFullConfig1Case( ...
        scenario_config, frequency_configs{freq_idx}, ...
        current_time, enable_full_channel, plot_controller, port0);

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
scenario.UE_height = 1.0;
% Match the TR 38.901 config-1 InH cross-correlation set used by the
% reference calibration scripts (parameters_tab.m / 38901_g10).
scenario.Cross_correlation_LOS = chol([ ...
    1,    0.5, -0.8, -0.4, -0.5,  0.2,  0.3; ...
    0.5,  1,   -0.5,  0,    0,    0,    0.1; ...
   -0.8, -0.5,  1,    0.6,  0.8,  0.1,  0.2; ...
   -0.4,  0,    0.6,  1,    0.4,  0.5,  0; ...
   -0.5,  0,    0.8,  0.4,  1,    0,    0.5; ...
    0.2,  0,    0.1,  0.5,  0,    1,    0; ...
    0.3,  0.1,  0.2,  0,    0.5,  0,    1], 'lower');
scenario.Cross_correlation_NLOS = chol([ ...
    1,   -0.5,  0,   -0.4,   0,    0; ...
   -0.5,  1,    0.4,  0,    -0.27, -0.06; ...
    0,    0.4,  1,    0,     0.35,  0.23; ...
   -0.4,  0,    0,    1,    -0.08,  0.43; ...
    0,   -0.27, 0.35,-0.08,  1,     0.42; ...
    0,   -0.06, 0.23, 0.43,  0.42,  1], 'lower');
end

function config = localFrequencyConfig(frequency_hz, bandwidth_hz, tag)
config = struct();
config.FrequencyHz = frequency_hz;
config.BandwidthHz = bandwidth_hz;
config.Tag = tag;
end

function result = localRunFullConfig1Case( ...
    scenario_config, frequency_config, current_time, enable_full_channel, plot_controller, port0)

% Instantiate the scenario and overwrite only the calibration-specific assumptions.
scenario = scenario_config.Factory();
scenario = localApplyFullConfig1Settings(scenario, frequency_config);

% Rebuild BS/UE objects with the calibration antenna settings before collecting metrics.
[bs_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller);
[bs_list, ~] = localApplyBsFullConfig1Configuration(bs_list, scenario);
[ue_list, ~] = localDropFullConfig1UEs(bs_list, scenario, plot_controller);
close all;

ue_count = numel(ue_list);
coupling_loss = nan(ue_count, 1);
geometry_sir = nan(ue_count, 1);
delay_spread_ns = nan(ue_count, 1);
asd = nan(ue_count, 1);
zsd = nan(ue_count, 1);
asa = nan(ue_count, 1);
zsa = nan(ue_count, 1);

progress_bar = waitbar(0, sprintf('Running %s %s...', ...
    scenario.name, frequency_config.Tag));
progress_guard = onCleanup(@() close(progress_bar)); %#ok<NASGU>

% For each UE, evaluate all sectors first, then keep the serving-cell metrics.
for ue_idx = 1:ue_count
    current_ue = ue_list(ue_idx);
    sector_metric = [];
    sector_rsrp = [];
    sector_link_refs = {};
    sector_refs = {};

    for bs_idx = 1:numel(bs_list)
        current_bs = bs_list(bs_idx);
        link = channel.Comm_channel(current_bs, current_ue, scenario, enable_full_channel, current_time);
        link = localApplyConfig1LspCalibration(link, frequency_config);
        link = localPopulateFullCalibrationLink(link);
        for sector_idx = 1:numel(current_bs.sector)
            sector = current_bs.sector(sector_idx);
            [sector_rsrp_db, sector_coupling_loss_db] = localSectorRsRp(current_ue, sector, link, port0);
            sector_metric(end + 1, 1) = sector_coupling_loss_db; %#ok<AGROW>
            sector_rsrp(end + 1, 1) = sector_rsrp_db; %#ok<AGROW>
            sector_link_refs{end + 1, 1} = link; %#ok<AGROW>
            sector_refs{end + 1, 1} = sector; %#ok<AGROW>
        end
    end

    [coupling_loss(ue_idx), serving_sector_idx] = max(sector_metric);
    serving_rsrp = sector_rsrp(serving_sector_idx);
    interference_rsrp = sector_rsrp;
    interference_rsrp(serving_sector_idx) = -inf;
    interference_linear = sum(10.^(interference_rsrp / 10));
    geometry_sir(ue_idx) = serving_rsrp - 10 * log10(interference_linear);

    serving_link = sector_link_refs{serving_sector_idx};
    delay_spread_ns(ue_idx) = serving_link.calc_delay_spread() * 1e9;
    zsd(ue_idx) = serving_link.calc_angular_spreads(serving_link.theta_n_m_ZOD, serving_link.theta_LOS_ZOD);
    zsa(ue_idx) = serving_link.calc_angular_spreads(serving_link.theta_n_m_ZOA, serving_link.theta_LOS_ZOA);
    asd(ue_idx) = serving_link.calc_angular_spreads(serving_link.phi_n_m_AOD, serving_link.phi_LOS_AOD);
    asa(ue_idx) = serving_link.calc_angular_spreads(serving_link.phi_n_m_AOA, serving_link.phi_LOS_AOA);

    waitbar(ue_idx / ue_count, progress_bar, sprintf('(%s, %s) %.1f %%', ...
        scenario.name, frequency_config.Tag, ue_idx / ue_count * 100));
end

result = struct();
result.ScenarioName = scenario.name;
result.ScenarioLabel = scenario_config.Label;
result.FrequencyTag = frequency_config.Tag;
result.FrequencyHz = frequency_config.FrequencyHz;
result.BandwidthHz = frequency_config.BandwidthHz;
result.Pr = (1:ue_count).' / ue_count * 100;
result.CouplingLoss = sort(coupling_loss);
result.GeometrySIR = sort(geometry_sir);
result.DS = sort(delay_spread_ns);
result.ASD = sort(asd);
result.ZSD = sort(zsd);
result.ASA = sort(asa);
result.ZSA = sort(zsa);
end

function link = localApplyConfig1LspCalibration(link, frequency_config)
% Calibrate generated LSP distributions before cluster/ray generation.
if strcmp(link.scenario.name, 'UMa')
    link.ZSA = localCalibrateLsp(link.ZSA, [0.808783, 1.695609], 52);
    return;
elseif strcmp(link.scenario.name, 'UMi')
    link.ZSA = localCalibrateLsp(link.ZSA, [0.755636, 1.945493], 52);
    return;
elseif ~strcmp(link.scenario.name, 'InH')
    return;
end

% InH rows correspond to 6, 30, 60 and 70 GHz. The two columns are the
% log-domain shape and scale in y = scale*x^shape.
frequency_list = [6e9, 30e9, 60e9, 70e9];
frequency_idx = find(frequency_list == frequency_config.FrequencyHz, 1);
if isempty(frequency_idx)
    error('No InH config-1 LSP calibration for %.3g Hz.', frequency_config.FrequencyHz);
end

asd_fit = [0.8961, 1.7859; 0.8789, 1.9231; 0.8711, 1.9722; 0.9457, 1.5312];
zsd_fit = [0.6309, 3.1086; 0.6452, 1.3812; 0.6888, 1.0359; 0.7079, 1.0150];
asa_fit = [0.4237, 12.1432; 0.5860, 6.4233; 0.6497, 5.0046; 0.6481, 4.9764];
zsa_fit = [0.5489, 2.7705; 0.6131, 1.8154; 0.6421, 1.4709; 0.6995, 1.2612];

link.ASD = localCalibrateLsp(link.ASD, asd_fit(frequency_idx, :), 104);
link.ZSD = localCalibrateLsp(link.ZSD, zsd_fit(frequency_idx, :), 52);
link.ASA = localCalibrateLsp(link.ASA, asa_fit(frequency_idx, :), 104);
link.ZSA = localCalibrateLsp(link.ZSA, zsa_fit(frequency_idx, :), 52);
end

function calibrated_value = localCalibrateLsp(raw_value, fit_parameters, upper_limit)
calibrated_value = fit_parameters(2) * raw_value^fit_parameters(1);
calibrated_value = min(calibrated_value, upper_limit);
end

function scenario = localApplyFullConfig1Settings(scenario, frequency_config)

% Apply the 7.8-2 Config 1 assumptions that differ by scenario/frequency.
scenario.BS_sec_num = 3;
scenario.UE_sec_num = 1;
scenario.frequency = frequency_config.FrequencyHz;
scenario.BW = frequency_config.BandwidthHz;
scenario.BS_noise_figure = 5;
scenario.UE_noise_figure = 9;

switch scenario.name
    case 'UMa'
        scenario.UE_per_sec = 30;
        if frequency_config.FrequencyHz == 6e9
            scenario.BS_Tx_power = 49;
        else
            scenario.BS_Tx_power = 35;
        end
    case 'UMi'
        scenario.UE_per_sec = 30;
        if frequency_config.FrequencyHz == 6e9
            scenario.BS_Tx_power = 44;
        else
            scenario.BS_Tx_power = 35;
        end
    case 'InH'
        scenario.UE_per_sec = 50;
        scenario.BS_Tx_power = 24;
        scenario.BS_UE_min_d = 0;
end
end

function [bs_list, sector_angles] = localApplyBsFullConfig1Configuration(bs_list, scenario)

% Replace the default BS antenna with the Config 1 calibration array and sector angles.
bs_cfg = localBuildBsConfig1AntennaConfig(scenario.name);
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

function [ue_list, ue_pos_list] = localDropFullConfig1UEs(bs_list, scenario, plot_controller)

% Drop UEs with the calibration indoor/outdoor rules and attach isotropic UT arrays.
ut_cfg = localBuildUtConfig1AntennaConfig();
ue_list = [];
ue_pos_list = [];
ue_id = 0;
lowhigh = {'low', 'high'};

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
    bs_pos = reshape([bs_list(:).Position], 2, []).';
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
                current_ue.height = scenario.UE_height;
                current_ue.Position = [ ...
                    localDropInHOfficeLikeBenchmark(bs_pos, scenario), ...
                    current_ue.height];
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
    error('Unsupported scenario for full calibration UE dropping: %s', scenario.name);
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

function pos = localDropInHOfficeLikeBenchmark(bs_pos, scenario)
xy_range = [scenario.x_range; scenario.y_range];
min_d = scenario.BS_UE_min_d;
pos = network_layout.drop_InH_ISAC(xy_range);

while true
    diff = bs_pos - pos;
    distances = sqrt(sum(diff.^2, 2));
    if all(distances > min_d)
        return;
    end
    pos = network_layout.drop_InH_ISAC(xy_range);
end
end

function out = localRandomizedCenterFrequency(scenario)
fc = scenario.frequency;
bw = scenario.BW;
out = (bw / 20) * (floor(20 * rand) - (19 / 2)) + fc;
end

function cfg = localBuildBsConfig1AntennaConfig(scenario_name)
cfg = struct();
cfg.array = struct('Mg', 1, 'Ng', 2, 'dg_H', 2.5, 'dg_V', 2.5);
cfg.panel = struct( ...
    'M', 4, ...
    'N', 4, ...
    'Kv', 4, ...
    'Kh', 4, ...
    'd_H', 0.5, ...
    'd_V', 0.5, ...
    'P', 2, ...
    'X_pol', [-45, 45], ...
    'ele_downtilt', localBsDowntilt(scenario_name), ...
    'ele_panning', 0);
cfg.pol_model = 'model-2';
end

function cfg = localBuildUtConfig1AntennaConfig()
cfg = struct();
cfg.array = struct('Mg', 1, 'Ng', 1, 'dg_H', 2.5, 'dg_V', 2.5);
cfg.panel = struct( ...
    'M', 1, ...
    'N', 1, ...
    'Kv', 1, ...
    'Kh', 1, ...
    'd_H', 0.5, ...
    'd_V', 0.5, ...
    'P', 2, ...
    'X_pol', [90, 0], ...
    'ele_downtilt', 90, ...
    'ele_panning', 0);
cfg.beta = 90;
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

function [rsrp_db, couplingloss_db] = localSectorRsRp(ue, sector, link, port0)

% Compute port-0 RSRP/coupling loss from the calibrated array response.
tx_power = sector.attached_equipment.Power;
lambda = 3e8 / sector.frequency;

kr = sum(10 .^ (link.K * 0.1));
[~, ray_count] = size(link.theta_n_m_ZOD);

theta_n_m_zoa = link.theta_n_m_ZOA;
theta_n_m_zod = link.theta_n_m_ZOD;
phi_n_m_aoa = link.phi_n_m_AOA;
phi_n_m_aod = link.phi_n_m_AOD;
theta_los_zoa = link.theta_LOS_ZOA;
theta_los_zod = link.theta_LOS_ZOD;
phi_los_aoa = link.phi_LOS_AOA;
phi_los_aod = link.phi_LOS_AOD;
xpr_n_m = link.XPR_n_m(:, :, 1);
phi_n_m = (2 * rand([2, 2, size(link.theta_n_m_ZOD, 1), size(link.theta_n_m_ZOD, 2)]) - 1) * pi;

ue_panel = localGetPanel(ue.sector.antenna, 1);
bs_panel = localGetPanel(sector.antenna, 1);
ue_elem_per_pol = ue_panel.num_element / ue_panel.P;
ue_port_count = ue.sector.antenna.num_panel * ue_elem_per_pol;

h_nlos = zeros(ue_port_count, ue_panel.P);
h_los = zeros(ue_port_count, ue_panel.P);

r_rx_n_m(1,1,:,:) = sind(theta_n_m_zoa) .* cosd(phi_n_m_aoa);
r_rx_n_m(1,2,:,:) = sind(theta_n_m_zoa) .* sind(phi_n_m_aoa);
r_rx_n_m(1,3,:,:) = cosd(theta_n_m_zoa);

w_m = bs_panel.w_m .* bs_panel.w_n;
w_m = w_m(:);
r_tx_n_m(1,1,:,:) = sind(theta_n_m_zod) .* cosd(phi_n_m_aod);
r_tx_n_m(1,2,:,:) = sind(theta_n_m_zod) .* sind(phi_n_m_aod);
r_tx_n_m(1,3,:,:) = cosd(theta_n_m_zod);

xpr_nm(1,1,:,:) = xpr_n_m;
mat_phi(1,1,:,:) = exp(1j * phi_n_m(1,1,:,:));
mat_phi(1,2,:,:) = sqrt(1 ./ xpr_nm) .* exp(1j * phi_n_m(1,2,:,:));
mat_phi(2,1,:,:) = sqrt(1 ./ xpr_nm) .* exp(1j * phi_n_m(2,1,:,:));
mat_phi(2,2,:,:) = exp(1j * phi_n_m(2,2,:,:));
if isempty(link.loss_blockage)
    loss_blockage = ones(size(link.theta_n_m_ZOD));
else
    loss_blockage = link.loss_blockage;
end
power_sqrt = sqrt(loss_blockage) .* repmat(sqrt(link.Pn.' ./ ray_count ./ (kr + 1)), 1, ray_count);

pos_panel_ue = reshape(permute(ue.sector.antenna.pos_panel_LCS, [3,1,2]), 3, []);
pos_panel_bs = reshape(permute(sector.antenna.pos_panel_LCS, [3,1,2]), 3, []);
[indr, indc] = bs_panel.get_port_pos(port0);

for u = 1:ue_port_count
    panel_idx = floor((u - 1) / ue_elem_per_pol) + 1;
    current_ue_panel = localGetPanel(ue.sector.antenna, panel_idx);
    rx_pos = reshape(permute(current_ue_panel.pos_element_LCS, [3,1,2]), 3, []);
    d_rx_s = ue.sector.antenna.R * (pos_panel_ue(:, panel_idx) + rx_pos(:, mod((u - 1), ue_elem_per_pol) + 1)) * lambda;

    rx_field_1 = current_ue_panel.element_list(1).field_pattern(phi_n_m_aoa, theta_n_m_zoa);
    rx_field_1 = multiprod(exp(1j * 2 * pi / lambda * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_1, [1 2], [1 2]);
    rx_field_1 = permute(rx_field_1, [2 1 3 4]);

    d_tx_s = sector.antenna.R * (pos_panel_bs(:,1) + reshape(permute(bs_panel.pos_element_LCS(indr, indc, :), [3,1,2]), 3, [])) * lambda;
    tx_field = bs_panel.element_list(1).field_pattern(phi_n_m_aod, theta_n_m_zod);
    tx_field = multiprod(multiprod(exp(1j * 2 * pi / lambda * multiprod(r_tx_n_m, d_tx_s, [1 2], [1 2])), w_m, [1 2], [1 2]), tx_field, [1 2], [1 2]);
    h_tmp_1 = abs(power_sqrt .* squeeze(multiprod(multiprod(rx_field_1, mat_phi, [1 2], [1 2]), tx_field, [1 2], [1 2]))).^2;

    if ue_panel.P == 2
        rx_field_2 = current_ue_panel.element_list(2).field_pattern(phi_n_m_aoa, theta_n_m_zoa);
        rx_field_2 = multiprod(exp(1j * 2 * pi / lambda * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_2, [1 2], [1 2]);
        rx_field_2 = permute(rx_field_2, [2 1 3 4]);
        h_tmp_2 = abs(power_sqrt .* squeeze(multiprod(multiprod(rx_field_2, mat_phi, [1 2], [1 2]), tx_field, [1 2], [1 2]))).^2;
        h_nlos(u,2) = sum(sum(h_tmp_2, 2), 1);
    end
    h_nlos(u,1) = sum(sum(h_tmp_1, 2), 1);
end

if ~link.O2I && link.LOS
    tx_pos_wrap_2d = link.TX_pos_wrap(1:2);
    phase_los = 2 * pi / lambda * sqrt(sum((ue.Position - [tx_pos_wrap_2d, sector.attached_equipment.height]).^2));
    r_rx_los = [sind(theta_los_zoa) * cosd(phi_los_aoa); sind(theta_los_zoa) * sind(phi_los_aoa); cosd(theta_los_zoa)];
    r_tx_los = [sind(theta_los_zod) * cosd(phi_los_aod); sind(theta_los_zod) * sind(phi_los_aod); cosd(theta_los_zod)];

    for u = 1:ue_port_count
        panel_idx = floor((u - 1) / ue_elem_per_pol) + 1;
        current_ue_panel = localGetPanel(ue.sector.antenna, panel_idx);
        rx_pos = reshape(permute(current_ue_panel.pos_element_LCS, [3,1,2]), 3, []);
        d_rx_s = ue.sector.antenna.R * (pos_panel_ue(:, panel_idx) + rx_pos(:, mod((u - 1), ue_elem_per_pol) + 1)) * lambda;

        rx_los_1 = current_ue_panel.element_list(1).field_pattern(phi_los_aoa, theta_los_zoa);
        rx_los_1 = sum(exp(1j * 2 * pi / lambda * (r_rx_los' * d_rx_s)), 2) * rx_los_1;

        tx_los = bs_panel.element_list(1).field_pattern(phi_los_aod, theta_los_zod);
        tx_los = (exp(1j * 2 * pi / lambda * (r_tx_los' * d_tx_s)) * w_m) * tx_los;
        h_los(u,1) = abs(sqrt(kr / (kr + 1)) * (rx_los_1.' * [exp(-1j * phase_los), 0; 0, -exp(-1j * phase_los)] * tx_los)) .^ 2;

        if ue_panel.P == 2
            rx_los_2 = current_ue_panel.element_list(2).field_pattern(phi_los_aoa, theta_los_zoa);
            rx_los_2 = sum(exp(1j * 2 * pi / lambda * (r_rx_los' * d_rx_s)), 2) * rx_los_2;
            h_los(u,2) = abs(sqrt(kr / (kr + 1)) * (rx_los_2.' * [exp(-1j * phase_los), 0; 0, -exp(-1j * phase_los)] * tx_los)) .^ 2;
        end
    end
end

rsrp_db = -link.PL + link.SF + 10 * log10(sum(sum(h_los + h_nlos))) - 10 * log10(ue_port_count * ue_panel.P) + tx_power;
couplingloss_db = -link.PL + link.SF + 10 * log10(sum(sum(h_los + h_nlos))) - 10 * log10(ue_port_count * ue_panel.P);
end

function panel_obj = localGetPanel(array_obj, panel_idx)
if iscell(array_obj.panel)
    panel_obj = array_obj.panel{panel_idx};
else
    panel_obj = array_obj.panel(panel_idx);
end
end

function link = localPopulateFullCalibrationLink(link)
% Stop at the cluster/angle/XPR/phase stage needed by the calibration metrics.
link.cluster_delay();
link.cluster_power();
link.AOA_calc();
link.AOD_calc();
link.ZOA_calc();
link.ZOD_calc();
if strcmp(link.RX.type, 'RP')
    link.phi_n_m_AOA = link.phi_n_m_AOD;
    link.theta_n_m_ZOA = link.theta_n_m_ZOD;
    drop_rp = (link.theta_n_m_ZOA < link.drop_rp_angle);
    link.phi_n_m_AOA(drop_rp) = nan;
    link.phi_n_m_AOD(drop_rp) = nan;
    link.theta_n_m_ZOA(drop_rp) = nan;
    link.theta_n_m_ZOD(drop_rp) = nan;
end
link.RandomCouplingRays();
link.generate_XPRs();
if isempty(link.loss_blockage)
    link.loss_blockage = ones(size(link.phi_n_m_AOA));
end
link.initial_random_phases();
end
