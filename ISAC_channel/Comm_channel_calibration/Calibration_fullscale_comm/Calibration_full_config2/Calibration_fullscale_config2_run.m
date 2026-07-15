function Calibration_fullscale_config2_run(selected_scenario, selected_frequency_tags)
if nargin < 1 || strlength(string(selected_scenario)) == 0
    selected_scenario = 'UMi';
end
if nargin < 2
    selected_frequency_tags = [];
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
result_root = fullfile(script_dir, 'results', 'comm_full_scale_config2');
if ~exist(result_root, 'dir')
    mkdir(result_root);
end

benchmark_file = fullfile(script_dir, 'Phase2Config2Calibration_v28_CMCC.xlsx');
if ~isfile(benchmark_file)
    error('Config 2 benchmark workbook was not found: %s', benchmark_file);
end

% Build the selectable scenario list so the same runner can be reused for UMa/UMi/InH.
scenario_configs = {
    localScenarioConfig('UMa', @() comm_scenario.UMa(), 'UMa')
    localScenarioConfig('UMi', @() comm_scenario.UMi(), 'UMi')
    localScenarioConfig('InH', @() localCreateInHScenario(), 'InH')
};
scenario_config = localSelectScenarioConfig(scenario_configs, selected_scenario);

% Config 2 uses the full-calibration carrier set, with optional frequency filtering.
frequency_configs = {
    localFrequencyConfig(6e9, 20e6, '6GHz')
    localFrequencyConfig(30e9, 100e6, '30GHz')
    localFrequencyConfig(60e9, 100e6, '60GHz')
    localFrequencyConfig(70e9, 100e6, '70GHz')
};
frequency_configs = localSelectFrequencyConfigs(frequency_configs, selected_frequency_tags);

% Run one calibration case per carrier frequency and store the sorted metrics.
for freq_idx = 1:numel(frequency_configs)
    fprintf('Running Config2 %s %s...\n', ...
        scenario_config.Name, frequency_configs{freq_idx}.Tag);
    result = localRunFullConfig2Case( ...
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

function selected_configs = localSelectFrequencyConfigs(all_configs, selected_frequency_tags)
if isempty(selected_frequency_tags)
    selected_configs = all_configs;
    return;
end

if ischar(selected_frequency_tags) || isstring(selected_frequency_tags)
    selected_keys = cellstr(string(selected_frequency_tags));
elseif iscell(selected_frequency_tags)
    selected_keys = cellfun(@(x) char(string(x)), selected_frequency_tags, 'UniformOutput', false);
else
    error('selected_frequency_tags must be a string, char vector, or cell array of strings.');
end

selected_configs = {};
for idx = 1:numel(all_configs)
    config = all_configs{idx};
    if any(strcmpi(config.Tag, selected_keys))
        selected_configs{end + 1, 1} = config; %#ok<AGROW>
    end
end

if isempty(selected_configs)
    error('No matching frequency tags were selected.');
end
end

function scenario = localCreateInHScenario()
% InH uses an explicit benchmark-like correlation setup for Config 2.
scenario = comm_scenario.InH();
scenario.subname = 'open_office';
scenario.UE_height = 1.0;
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

function result = localRunFullConfig2Case( ...
    scenario_config, frequency_config, current_time, enable_full_channel, plot_controller, port0)

% Instantiate the scenario and overwrite only the calibration-specific assumptions.
scenario = scenario_config.Factory();
scenario = localApplyFullConfig2Settings(scenario, frequency_config);

% Rebuild BS/UE objects with the Config 2 antenna settings before collecting metrics.
[bs_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller);
[bs_list, ~] = localApplyBsFullConfig2Configuration(bs_list, scenario);
[ue_list, ~] = localDropFullConfig2UEs(bs_list, scenario, plot_controller);
close all;

ue_count = numel(ue_list);
coupling_loss = nan(ue_count, 1);
geometry_sir = nan(ue_count, 1);
largest_sv = nan(ue_count, 1);
smallest_sv = nan(ue_count, 1);
ratio_sv = nan(ue_count, 1);

progress_bar = waitbar(0, sprintf('Running Config2 %s %s...', ...
    scenario.name, frequency_config.Tag));
progress_guard = onCleanup(@() close(progress_bar)); %#ok<NASGU>

% For each UE, evaluate all sectors first, then extract serving-cell singular-value metrics.
for ue_idx = 1:ue_count
    current_ue = ue_list(ue_idx);
    sector_metric = [];
    sector_rsrp = [];
    sector_link_refs = {};
    sector_refs = {};
    sector_phase_refs = {};

    for bs_idx = 1:numel(bs_list)
        current_bs = bs_list(bs_idx);
        link = channel.Comm_channel(current_bs, current_ue, scenario, enable_full_channel, current_time);
        link = localPopulateFullCalibrationLink(link);
        for sector_idx = 1:numel(current_bs.sector)
            sector = current_bs.sector(sector_idx);
            [sector_rsrp_db, sector_coupling_loss_db, sector_phi_n_m] = localSectorRsRp(current_ue, sector, link, port0);
            sector_metric(end + 1, 1) = sector_coupling_loss_db; %#ok<AGROW>
            sector_rsrp(end + 1, 1) = sector_rsrp_db; %#ok<AGROW>
            sector_link_refs{end + 1, 1} = link; %#ok<AGROW>
            sector_refs{end + 1, 1} = sector; %#ok<AGROW>
            sector_phase_refs{end + 1, 1} = sector_phi_n_m; %#ok<AGROW>
        end
    end

    [serving_rsrp, serving_sector_idx] = max(sector_rsrp);
    coupling_loss(ue_idx) = sector_metric(serving_sector_idx);
    interference_rsrp = sector_rsrp;
    interference_rsrp(serving_sector_idx) = -inf;
    interference_linear = sum(10.^(interference_rsrp / 10));
    geometry_sir(ue_idx) = serving_rsrp - 10 * log10(interference_linear);

    serving_link = sector_link_refs{serving_sector_idx};
    serving_sector = sector_refs{serving_sector_idx};
    serving_link.PHI_n_m = sector_phase_refs{serving_sector_idx};
    [largest_sv_db, smallest_sv_db, ratio_sv_db] = ...
        localComputeSingularValues(current_ue, serving_sector, serving_link, current_time);
    largest_sv(ue_idx) = largest_sv_db;
    smallest_sv(ue_idx) = smallest_sv_db;
    ratio_sv(ue_idx) = ratio_sv_db;

    waitbar(ue_idx / ue_count, progress_bar, sprintf('(Config2 %s, %s) %.1f %%', ...
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
result.LargestSingularValue = sort(largest_sv);
result.SmallestSingularValue = sort(smallest_sv);
result.RatioSingularValue = sort(ratio_sv);
end

function scenario = localApplyFullConfig2Settings(scenario, frequency_config)
% Apply the 7.8-2 Config 2 assumptions that differ by scenario/frequency.
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

function [bs_list, sector_angles] = localApplyBsFullConfig2Configuration(bs_list, scenario)
% Replace the default BS antenna with the Config 2 calibration array and sector angles.
bs_cfg = localBuildBsConfig2AntennaConfig(scenario.name);
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

function [ue_list, ue_pos_list] = localDropFullConfig2UEs(bs_list, scenario, plot_controller)
% Drop UEs with the calibration indoor/outdoor rules and attach isotropic UT arrays.
ut_cfg = localBuildUtConfig2AntennaConfig();
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
                current_ue.height = 1.0;
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

function cfg = localBuildBsConfig2AntennaConfig(scenario_name)
cfg = struct();
cfg.array = struct('Mg', 1, 'Ng', 1, 'dg_H', 2.5, 'dg_V', 2.5);
cfg.panel = struct( ...
    'M', 2, ...
    'N', 2, ...
    'Kv', 1, ...
    'Kh', 1, ...
    'd_H', 0.5, ...
    'd_V', 0.5, ...
    'P', 1, ...
    'X_pol', 0, ...
    'ele_downtilt', localBsDowntilt(scenario_name), ...
    'ele_panning', 0);
cfg.pol_model = 'model-2';
end

function cfg = localBuildUtConfig2AntennaConfig()
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

function [rsrp_db, couplingloss_db, phi_n_m] = localSectorRsRp(ue, sector, link, port0)
% Compute port-0 RSRP/coupling loss and keep the random phases for later channel rebuild.
tx_power = sector.attached_equipment.Power;
lambda = 3e8 / sector.frequency;
if ~isempty(ue.fcin)
    lambda_mt = 3e8 / ue.fcin;
else
    lambda_mt = [];
end

kr = sum(10 .^ (link.K * 0.1));
[cluster_count, ray_count] = size(link.theta_n_m_ZOD);

theta_n_m_zoa = link.theta_n_m_ZOA;
theta_n_m_zod = link.theta_n_m_ZOD;
phi_n_m_aoa = link.phi_n_m_AOA;
phi_n_m_aod = link.phi_n_m_AOD;
theta_los_zoa = link.theta_LOS_ZOA;
theta_los_zod = link.theta_LOS_ZOD;
phi_los_aoa = link.phi_LOS_AOA;
phi_los_aod = link.phi_LOS_AOD;
xpr_n_m = link.XPR_n_m(:, :, 1);
phi_n_m = (2 * rand([2, 2, cluster_count, ray_count]) - 1) * pi;

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
power_sqrt = repmat(sqrt(link.Pn.' ./ ray_count ./ (kr + 1)), 1, ray_count);

pos_panel_ue = reshape(permute(ue.sector.antenna.pos_panel_LCS, [3,1,2]), 3, []);
pos_panel_bs = reshape(permute(sector.antenna.pos_panel_LCS, [3,1,2]), 3, []);
[indr, indc] = bs_panel.get_port_pos(port0);

for u = 1:ue_port_count
    panel_idx = floor((u - 1) / ue_elem_per_pol) + 1;
    current_ue_panel = localGetPanel(ue.sector.antenna, panel_idx);
    rx_pos = reshape(permute(current_ue_panel.pos_element_LCS, [3,1,2]), 3, []);
    d_rx_s = ue.sector.antenna.R * (pos_panel_ue(:, panel_idx) + rx_pos(:, mod((u - 1), ue_elem_per_pol) + 1)) * lambda;

    rx_field_1 = current_ue_panel.element_list(1).field_pattern(phi_n_m_aoa, theta_n_m_zoa);
    if isempty(lambda_mt)
        rx_field_1 = multiprod(exp(1j * 2 * pi / lambda * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_1, [1 2], [1 2]);
    else
        rx_field_1 = multiprod(exp(1j * 2 * pi / lambda_mt * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_1, [1 2], [1 2]);
    end
    rx_field_1 = permute(rx_field_1, [2 1 3 4]);

    d_tx_s = sector.antenna.R * (pos_panel_bs(:,1) + reshape(permute(bs_panel.pos_element_LCS(indr, indc, :), [3,1,2]), 3, [])) * lambda;
    tx_field = bs_panel.element_list(1).field_pattern(phi_n_m_aod, theta_n_m_zod);
    if isempty(lambda_mt)
        tx_field = multiprod(multiprod(exp(1j * 2 * pi / lambda * multiprod(r_tx_n_m, d_tx_s, [1 2], [1 2])), w_m, [1 2], [1 2]), tx_field, [1 2], [1 2]);
    else
        tx_field = multiprod(multiprod(exp(1j * 2 * pi / lambda_mt * multiprod(r_tx_n_m, d_tx_s, [1 2], [1 2])), w_m, [1 2], [1 2]), tx_field, [1 2], [1 2]);
    end
    h_tmp_1 = abs(power_sqrt .* squeeze(multiprod(multiprod(rx_field_1, mat_phi, [1 2], [1 2]), tx_field, [1 2], [1 2]))).^2;

    if ue_panel.P == 2
        rx_field_2 = current_ue_panel.element_list(2).field_pattern(phi_n_m_aoa, theta_n_m_zoa);
        if isempty(lambda_mt)
            rx_field_2 = multiprod(exp(1j * 2 * pi / lambda * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_2, [1 2], [1 2]);
        else
            rx_field_2 = multiprod(exp(1j * 2 * pi / lambda_mt * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_2, [1 2], [1 2]);
        end
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
        if isempty(lambda_mt)
            rx_los_1 = sum(exp(1j * 2 * pi / lambda * (r_rx_los' * d_rx_s)), 2) * rx_los_1;
        else
            rx_los_1 = sum(exp(1j * 2 * pi / lambda_mt * (r_rx_los' * d_rx_s)), 2) * rx_los_1;
        end

        tx_los = bs_panel.element_list(1).field_pattern(phi_los_aod, theta_los_zod);
        if isempty(lambda_mt)
            tx_los = (exp(1j * 2 * pi / lambda * (r_tx_los' * d_tx_s)) * w_m) * tx_los;
        else
            tx_los = (exp(1j * 2 * pi / lambda_mt * (r_tx_los' * d_tx_s)) * w_m) * tx_los;
        end
        h_los(u,1) = abs(sqrt(kr / (kr + 1)) * (rx_los_1.' * [exp(-1j * phase_los), 0; 0, -exp(-1j * phase_los)] * tx_los)) .^ 2;

        if ue_panel.P == 2
            rx_los_2 = current_ue_panel.element_list(2).field_pattern(phi_los_aoa, theta_los_zoa);
            if isempty(lambda_mt)
                rx_los_2 = sum(exp(1j * 2 * pi / lambda * (r_rx_los' * d_rx_s)), 2) * rx_los_2;
            else
                rx_los_2 = sum(exp(1j * 2 * pi / lambda_mt * (r_rx_los' * d_rx_s)), 2) * rx_los_2;
            end
            h_los(u,2) = abs(sqrt(kr / (kr + 1)) * (rx_los_2.' * [exp(-1j * phase_los), 0; 0, -exp(-1j * phase_los)] * tx_los)) .^ 2;
        end
    end
end

rsrp_db = -link.PL + link.SF + 10 * log10(sum(sum(h_los + h_nlos))) - 10 * log10(ue_port_count * ue_panel.P) + tx_power;
couplingloss_db = -link.PL + link.SF + 10 * log10(sum(sum(h_los + h_nlos))) - 10 * log10(ue_port_count * ue_panel.P);
end

function [largest_sv_db, smallest_sv_db, ratio_sv_db] = localComputeSingularValues(ue, sector, link, current_time)
% Build the serving-link covariance and reduce it to the Config 2 singular-value metrics.
[tau, channel_coef] = localBuildServingChannel(ue, sector, link, current_time);
fc = sector.frequency;
frequency_bins = fc;
expf = permute(repmat(exp(-1j * 2 * pi * tau * frequency_bins), 1, 1, size(channel_coef, 1), size(channel_coef, 2)), [3, 4, 1, 2]);
channel = repmat(channel_coef, 1, 1, 1, numel(frequency_bins));
channel_freq = permute(sum(channel .* expf, 3), [1, 2, 4, 3]);
mean_covariance = sum(multiprod(channel_freq, permute(conj(channel_freq), [2, 1, 3]), [1, 2], [1, 2]), 3);
singular_linear = svd(mean_covariance) / numel(frequency_bins);
singular_linear = sort(real(singular_linear(:)), 'ascend');

if numel(singular_linear) < 2
    singular_linear(2, 1) = singular_linear(1);
end

smallest_sv_db = 10 * log10(max(singular_linear(1), realmin));
largest_sv_db = 10 * log10(max(singular_linear(end), realmin));
ratio_sv_db = 10 * log10(max(singular_linear(end), realmin) / max(singular_linear(1), realmin));
end

function [tau, channel_coef] = localBuildServingChannel(ue, sector, link, current_time)
% Rebuild the serving-link channel coefficients at element level for the singular-value calculation.
fc = sector.frequency;
lambda = 3e8 / fc;
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
phi_n_m = link.PHI_n_m;
tau_n = link.tau_n_LOS;
cluster_idx = 1:length(tau_n);
strong_cluster_id = link.strong_cluster_id;
map_delay = link.map_delay;
delay_shifts = unique(map_delay);
delay_map = zeros(length(delay_shifts), length(map_delay));
for delay_idx = 1:length(delay_shifts)
    delay_map(delay_idx, :) = map_delay == delay_shifts(delay_idx);
end

tau_main = tau_n(cluster_idx ~= strong_cluster_id(1) & cluster_idx ~= strong_cluster_id(2));
tau_sub = tau_n(strong_cluster_id) + delay_shifts.';
tau = [tau_main(:); tau_sub(:)];
[tau, sort_idx] = sort(tau);

ue_panel = localGetPanel(ue.sector.antenna, 1);
bs_panel = localGetPanel(sector.antenna, 1);
ue_pol_count = ue_panel.P;
ue_elem_per_pol = ue_panel.num_element / ue_pol_count;
ue_port_count = ue.sector.antenna.num_panel * ue_elem_per_pol;

bs_pol_count = bs_panel.P;
bs_elem_per_pol = bs_panel.num_element / bs_pol_count;
bs_port_count = sector.antenna.num_panel * bs_elem_per_pol;

channel_matrix = complex(zeros(ue_port_count * ue_pol_count, bs_port_count * bs_pol_count, numel(tau)));
los_matrix = complex(zeros(ue_port_count * ue_pol_count, bs_port_count * bs_pol_count, 1));

r_rx_n_m(1,1,:,:) = sind(theta_n_m_zoa) .* cosd(phi_n_m_aoa);
r_rx_n_m(1,2,:,:) = sind(theta_n_m_zoa) .* sind(phi_n_m_aoa);
r_rx_n_m(1,3,:,:) = cosd(theta_n_m_zoa);

rx_velocity = ue.velocity * [sind(ue.theta_v) * cosd(ue.phi_v); sind(ue.theta_v) * sind(ue.phi_v); cosd(ue.theta_v)];
phase_vnm = squeeze(exp(1j * 2 * pi / lambda * multiprod(r_rx_n_m, rx_velocity, [1, 2], [1, 2]) * current_time));

r_tx_n_m(1,1,:,:) = sind(theta_n_m_zod) .* cosd(phi_n_m_aod);
r_tx_n_m(1,2,:,:) = sind(theta_n_m_zod) .* sind(phi_n_m_aod);
r_tx_n_m(1,3,:,:) = cosd(theta_n_m_zod);

xpr_nm(1,1,:,:) = xpr_n_m;
mat_phi(1,1,:,:) = exp(1j * phi_n_m(1,1,:,:));
mat_phi(1,2,:,:) = sqrt(1 ./ xpr_nm) .* exp(1j * phi_n_m(1,2,:,:));
mat_phi(2,1,:,:) = sqrt(1 ./ xpr_nm) .* exp(1j * phi_n_m(2,1,:,:));
mat_phi(2,2,:,:) = exp(1j * phi_n_m(2,2,:,:));
power_sqrt = repmat(sqrt(link.Pn.' ./ ray_count ./ (kr + 1)), 1, ray_count);

pos_panel_ue = reshape(permute(ue.sector.antenna.pos_panel_LCS, [3,1,2]), 3, []);
pos_panel_bs = reshape(permute(sector.antenna.pos_panel_LCS, [3,1,2]), 3, []);

for u = 0:ue_port_count - 1
    ue_panel_idx = floor(u / ue_elem_per_pol) + 1;
    current_ue_panel = localGetPanel(ue.sector.antenna, ue_panel_idx);
    rx_pos = reshape(permute(current_ue_panel.pos_element_LCS, [3,1,2]), 3, []);
    d_rx_s = ue.sector.antenna.R * (pos_panel_ue(:, ue_panel_idx) + rx_pos(:, u + 1)) * lambda;

    rx_field_1 = current_ue_panel.element_list(1).field_pattern(phi_n_m_aoa, theta_n_m_zoa);
    rx_field_1 = multiprod(exp(1j * 2 * pi / lambda * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_1, [1 2], [1 2]);
    rx_field_1 = permute(rx_field_1, [2 1 3 4]);
    if ue_pol_count == 2
        rx_field_2 = current_ue_panel.element_list(2).field_pattern(phi_n_m_aoa, theta_n_m_zoa);
        rx_field_2 = multiprod(exp(1j * 2 * pi / lambda * multiprod(r_rx_n_m, d_rx_s, [1 2], [1 2])), rx_field_2, [1 2], [1 2]);
        rx_field_2 = permute(rx_field_2, [2 1 3 4]);
    end

    for p = 0:bs_port_count - 1
        bs_panel_idx = floor(p / bs_elem_per_pol) + 1;
        current_bs_panel = localGetPanel(sector.antenna, bs_panel_idx);
        port_idx = mod(p, bs_elem_per_pol);
        [indr, indc] = current_bs_panel.get_element_pos(port_idx);
        d_tx_s = sector.antenna.R * (pos_panel_bs(:, bs_panel_idx) + reshape(permute(current_bs_panel.pos_element_LCS(indr, indc, :), [3,1,2]), 3, [])) * lambda;

        tx_field_1 = current_bs_panel.element_list(1).field_pattern(phi_n_m_aod, theta_n_m_zod);
        tx_field_1 = multiprod( ...
            exp(1j * 2 * pi / lambda * multiprod(r_tx_n_m, d_tx_s, [1 2], [1 2])), ...
            tx_field_1, [1 2], [1 2]);
        h_tmp_1 = power_sqrt .* squeeze(multiprod(multiprod(rx_field_1, mat_phi, [1 2], [1 2]), tx_field_1, [1 2], [1 2])) .* phase_vnm;
        h_tmp_11 = sum(h_tmp_1((cluster_idx ~= strong_cluster_id(1) & cluster_idx ~= strong_cluster_id(2)), :), 2);
        h_tmp_12 = h_tmp_1(strong_cluster_id, :);
        h_tmp_12 = (h_tmp_12 * delay_map.').';
        h_tmp_1 = [h_tmp_11; h_tmp_12(:)];
        channel_matrix((u * ue_pol_count + 1), p * bs_pol_count + 1, :) = h_tmp_1(sort_idx);

        if ue_pol_count == 2
            h_tmp_2 = power_sqrt .* squeeze(multiprod(multiprod(rx_field_2, mat_phi, [1 2], [1 2]), tx_field_1, [1 2], [1 2])) .* phase_vnm;
            h_tmp_21 = sum(h_tmp_2((cluster_idx ~= strong_cluster_id(1) & cluster_idx ~= strong_cluster_id(2)), :), 2);
            h_tmp_22 = h_tmp_2(strong_cluster_id, :);
            h_tmp_22 = (h_tmp_22 * delay_map.').';
            h_tmp_2 = [h_tmp_21; h_tmp_22(:)];
            channel_matrix((u + 1) * ue_pol_count, p * bs_pol_count + 1, :) = h_tmp_2(sort_idx);
        end

        if bs_pol_count == 2
            tx_field_2 = current_bs_panel.element_list(2).field_pattern(phi_n_m_aod, theta_n_m_zod);
            tx_field_2 = multiprod( ...
                exp(1j * 2 * pi / lambda * multiprod(r_tx_n_m, d_tx_s, [1 2], [1 2])), ...
                tx_field_2, [1 2], [1 2]);
            h_tmp_3 = power_sqrt .* squeeze(multiprod(multiprod(rx_field_1, mat_phi, [1 2], [1 2]), tx_field_2, [1 2], [1 2])) .* phase_vnm;
            h_tmp_31 = sum(h_tmp_3((cluster_idx ~= strong_cluster_id(1) & cluster_idx ~= strong_cluster_id(2)), :), 2);
            h_tmp_32 = h_tmp_3(strong_cluster_id, :);
            h_tmp_32 = (h_tmp_32 * delay_map.').';
            h_tmp_3 = [h_tmp_31; h_tmp_32(:)];
            channel_matrix((u * ue_pol_count + 1), (p + 1) * bs_pol_count, :) = h_tmp_3(sort_idx);

            if ue_pol_count == 2
                h_tmp_4 = power_sqrt .* squeeze(multiprod(multiprod(rx_field_2, mat_phi, [1 2], [1 2]), tx_field_2, [1 2], [1 2])) .* phase_vnm;
                h_tmp_41 = sum(h_tmp_4((cluster_idx ~= strong_cluster_id(1) & cluster_idx ~= strong_cluster_id(2)), :), 2);
                h_tmp_42 = h_tmp_4(strong_cluster_id, :);
                h_tmp_42 = (h_tmp_42 * delay_map.').';
                h_tmp_4 = [h_tmp_41; h_tmp_42(:)];
                channel_matrix((u + 1) * ue_pol_count, (p + 1) * bs_pol_count, :) = h_tmp_4(sort_idx);
            end
        end
    end
end

if ~link.O2I && link.LOS
    tx_pos_wrap_2d = link.TX_pos_wrap(1:2);
    d_los = sqrt(sum((ue.Position - [tx_pos_wrap_2d, sector.attached_equipment.height]).^2));
    phi_los = 2 * pi / lambda * d_los;
    r_rx_los = [sind(theta_los_zoa) * cosd(phi_los_aoa); sind(theta_los_zoa) * sind(phi_los_aoa); cosd(theta_los_zoa)];
    r_tx_los = [sind(theta_los_zod) * cosd(phi_los_aod); sind(theta_los_zod) * sind(phi_los_aod); cosd(theta_los_zod)];
    phase_los = exp(1j * 2 * pi / lambda * (r_rx_los.' * rx_velocity) * current_time);

    for u = 0:ue_port_count - 1
        ue_panel_idx = floor(u / ue_elem_per_pol) + 1;
        current_ue_panel = localGetPanel(ue.sector.antenna, ue_panel_idx);
        rx_pos = reshape(permute(current_ue_panel.pos_element_LCS, [3,1,2]), 3, []);
        d_rx_s = ue.sector.antenna.R * (pos_panel_ue(:, ue_panel_idx) + rx_pos(:, u + 1)) * lambda;
        rx_los_1 = current_ue_panel.element_list(1).field_pattern(phi_los_aoa, theta_los_zoa);
        rx_los_1 = (exp(1j * 2 * pi / lambda * (r_rx_los' * d_rx_s))) * rx_los_1;
        if ue_pol_count == 2
            rx_los_2 = current_ue_panel.element_list(2).field_pattern(phi_los_aoa, theta_los_zoa);
            rx_los_2 = (exp(1j * 2 * pi / lambda * (r_rx_los' * d_rx_s))) * rx_los_2;
        end

        for p = 0:bs_port_count - 1
            bs_panel_idx = floor(p / bs_elem_per_pol) + 1;
            current_bs_panel = localGetPanel(sector.antenna, bs_panel_idx);
            port_idx = mod(p, bs_elem_per_pol);
            [indr, indc] = current_bs_panel.get_element_pos(port_idx);
            d_tx_s = sector.antenna.R * (pos_panel_bs(:, bs_panel_idx) + reshape(permute(current_bs_panel.pos_element_LCS(indr, indc, :), [3,1,2]), 3, [])) * lambda;
            tx_los_1 = current_bs_panel.element_list(1).field_pattern(phi_los_aod, theta_los_zod);
            tx_los_1 = (exp(1j * 2 * pi / lambda * (r_tx_los' * d_tx_s))) * tx_los_1;
            los_matrix((u * ue_pol_count + 1), p * bs_pol_count + 1, 1) = ...
                (rx_los_1.' * [exp(-1j * phi_los), 0; 0, -exp(-1j * phi_los)] * tx_los_1) .* phase_los;

            if ue_pol_count == 2
                los_matrix((u + 1) * ue_pol_count, p * bs_pol_count + 1, 1) = ...
                    (rx_los_2.' * [exp(-1j * phi_los), 0; 0, -exp(-1j * phi_los)] * tx_los_1) .* phase_los;
            end

            if bs_pol_count == 2
                tx_los_2 = current_bs_panel.element_list(2).field_pattern(phi_los_aod, theta_los_zod);
                tx_los_2 = (exp(1j * 2 * pi / lambda * (r_tx_los' * d_tx_s))) * tx_los_2;
                los_matrix((u * ue_pol_count + 1), (p + 1) * bs_pol_count, 1) = ...
                    (rx_los_1.' * [exp(-1j * phi_los), 0; 0, -exp(-1j * phi_los)] * tx_los_2) .* phase_los;

                if ue_pol_count == 2
                    los_matrix((u + 1) * ue_pol_count, (p + 1) * bs_pol_count, 1) = ...
                        (rx_los_2.' * [exp(-1j * phi_los), 0; 0, -exp(-1j * phi_los)] * tx_los_2) .* phase_los;
                end
            end
        end
    end

    channel_matrix(:, :, 1) = channel_matrix(:, :, 1) + sqrt(kr / (kr + 1)) * los_matrix;
end

channel_coef = channel_matrix;
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
link.initial_random_phases();
end
