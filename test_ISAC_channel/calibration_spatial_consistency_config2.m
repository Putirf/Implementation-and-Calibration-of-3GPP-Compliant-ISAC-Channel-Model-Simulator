function results = calibration_spatial_consistency_config2(selected_procedure, options)
%CALIBRATION_SPATIAL_CONSISTENCY_CONFIG2 Config2 mobility calibration.
% Based on TR 38.901 V19.1.0 Table 7.8-5:
% UMi, 30 GHz, outdoor moving UTs, 30 km/h, metrics 1/2/7/8/9.
%
% Examples:
%   calibration_spatial_consistency_config2('A')
%   calibration_spatial_consistency_config2('B')
%   calibration_spatial_consistency_config2({'A','B'})
%
% The default run is a small smoke test. Increase num_ue and duration_s
% through OPTIONS for a production calibration run.

if nargin < 1 || isempty(selected_procedure)
    selected_procedure = 'A';
end
if nargin < 2 || isempty(options)
    options = struct();
end
procedures = localNormalizeProcedures(selected_procedure);
options = localDefaultOptions(options);

results = repmat(localEmptyResult(), 1, numel(procedures));
for procedure_idx = 1:numel(procedures)
    procedure = procedures{procedure_idx};
    rng(options.rng_seed, 'twister');
    results(procedure_idx) = localRunProcedure(procedure, options, procedure_idx == 1);
    if options.save_results
        if ~exist(options.output_dir, 'dir')
            mkdir(options.output_dir);
        end
        result = results(procedure_idx); %#ok<NASGU>
        save(fullfile(options.output_dir, ['Config2_Proc' procedure '.mat']), 'result');
    end
end
if options.plot_results
    localPlotResults(results, options);
end
if options.plot_time_series
    localPlotTimeSeries(results, options);
end
end

function result = localRunProcedure(procedure, options, plot_scene)
scenario = comm_scenario.UMi;
scenario.layer_num = 1;
scenario.BS_sec_num = 3;
scenario.UE_sec_num = 1;
scenario.UE_per_sec = options.num_ue;
scenario.ST_per_cell = 0;
scenario.RP_per_equipment = 0;
scenario.applyIsacFrequencyPreset('ISAC_FR2');
% 3GPP TR 38.901 calibration assumptions for 30 GHz UMi.
scenario.BS_Tx_power = 35; % dBm
scenario.BW = 100e6;       % Hz
scenario.spatial_consistency_enable = true;
scenario.spatial_consistency_procedure = procedure;

time_nodes = 0:options.sample_dt:options.duration_s;
if time_nodes(end) < options.duration_s
    time_nodes(end+1) = options.duration_s; %#ok<AGROW>
end

[BS_list, BS_sector_list] = localDropSingleBs(scenario, [0 0 scenario.BS_height]);
BS = BS_list(1);

[scenario.total_BS_sector_num] = numel(BS_sector_list);
[UE_list, ~] = localDropOutdoorUmiUes( ...
    BS, BS_sector_list, scenario, options.duration_s);

if plot_scene && options.plot_scene
    localPlotScenario(BS_list, UE_list, scenario, procedure, options);
end

if strcmp(procedure,'A')
    % Procedure-A X_n uses the UMi-specific 15 m decorrelation distance.
    BS.SC_procA_comm_Xn = localProcedureAXnTable(UE_list, 15, 20);
    [BS.SC_procA_comm_raw_LOS, BS.SC_procA_comm_raw_NLOS, ...
        BS.SC_procA_comm_raw_O2I, initial_lsp] = ...
        localProcedureAInitialFields(UE_list,20);
    BS.LSP_raw_LOS = initial_lsp.LOS;
    BS.LSP_raw_NLOS = initial_lsp.NLOS;
    BS.LSP_raw_O2I = initial_lsp.O2I;
else
    [BS.SC_procB_raw, lsp_path] = localProcedureBPaths( ...
        UE_list, time_nodes);
    procB_path = BS.SC_procB_raw;
end

num_ue = numel(UE_list);
num_samples = numel(time_nodes);
total_steps = num_samples*num_ue;
progress_step = max(1,floor(total_steps/100));
completed_steps = 0;
progress_bar = localCreateProgressBar(options,scenario);
progress_cleanup = onCleanup(@()localCloseProgressBar(progress_bar)); %#ok<NASGU>
power3 = nan(num_samples, num_ue);
delay3 = nan(num_samples, num_ue);
aoa3 = nan(num_samples, num_ue);
zoa3 = nan(num_samples, num_ue);
absolute_delay3 = nan(num_samples, num_ue);
coupling_loss = nan(1, num_ue);
sir_db = nan(1, num_ue);
sinr_db = nan(1, num_ue);
serving_sector_idx = nan(1, num_ue);
cluster3_active_idx = nan(1, num_ue);
cluster3_original_id = nan(1, num_ue);
los_flag = false(1, num_ue);
previous_los_state = false(1,num_ue);
procedure_a_regenerated_any = false(1,num_ue);
procedure_a_regeneration_count = zeros(1,num_ue);
los_nlos_transition_any = false(1,num_ue);
los_nlos_transition_count = zeros(1,num_ue);
initial_k_db = nan(1,num_ue);
initial_zsa_deg = nan(1,num_ue);
initial_theta_los_zoa_deg = nan(1,num_ue);
initial_cluster3_power_for_angle = nan(1,num_ue);

for time_idx = 1:num_samples
    t_k = time_nodes(time_idx);
    if strcmp(procedure, 'B')
        localAssignProcedureBTimeSlice(BS, procB_path, lsp_path, time_idx);
    end
    for ue_idx = 1:num_ue
        UE = UE_list(ue_idx);
        UE.Position = localMoveFromInitial(UE, t_k);
        UE.height = UE.Position(3);
        if time_idx == 1
            [link, serving_sector_idx(ue_idx), coupling_loss(ue_idx), ...
                sir_db(ue_idx), sinr_db(ue_idx)] = localSelectServingSector( ...
                BS, BS_sector_list, UE, scenario, t_k);
            [cluster3_active_idx(ue_idx),cluster3_original_id(ue_idx)] = ...
                localThirdActiveClusterIndex(link);
            los_flag(ue_idx) = link.LOS && ~link.O2I;
            previous_los_state(ue_idx) = link.LOS;
        else
            BS.sector = BS_sector_list(serving_sector_idx(ue_idx));
            link = channel.Comm_channel(BS, UE, scenario, true, t_k);
            regenerated_now = link.procedureA_regenerated;
            los_transition_now = link.LOS ~= previous_los_state(ue_idx);
            procedure_a_regenerated_any(ue_idx) = ...
                procedure_a_regenerated_any(ue_idx) || regenerated_now;
            procedure_a_regeneration_count(ue_idx) = ...
                procedure_a_regeneration_count(ue_idx) + double(regenerated_now);
            los_nlos_transition_any(ue_idx) = ...
                los_nlos_transition_any(ue_idx) || los_transition_now;
            los_nlos_transition_count(ue_idx) = ...
                los_nlos_transition_count(ue_idx) + double(los_transition_now);
            previous_los_state(ue_idx) = link.LOS;
            if link.procedureA_regenerated
                [cluster3_active_idx(ue_idx),cluster3_original_id(ue_idx)] = ...
                    localThirdActiveClusterIndex(link);
            else
                cluster3_active_idx(ue_idx) = localTrackedClusterIndex( ...
                    link,cluster3_original_id(ue_idx));
            end
        end
        [power3(time_idx, ue_idx), delay3(time_idx, ue_idx), ...
            aoa3(time_idx, ue_idx),zoa3(time_idx,ue_idx), ...
            absolute_delay3(time_idx,ue_idx)] = localThirdClusterValues( ...
            link, cluster3_active_idx(ue_idx));
        if time_idx == 1
            active_idx = cluster3_active_idx(ue_idx);
            if ~isempty(link.K)
                initial_k_db(ue_idx) = link.K(1);
            end
            initial_zsa_deg(ue_idx) = link.ZSA;
            initial_theta_los_zoa_deg(ue_idx) = link.theta_LOS_ZOA;
            if isfinite(active_idx) && numel(link.Pn_LOS) >= active_idx
                initial_cluster3_power_for_angle(ue_idx) = ...
                    link.Pn_LOS(round(active_idx));
            end
        end
        completed_steps = completed_steps + 1;
        if mod(completed_steps,progress_step) == 0 || ...
                completed_steps == total_steps
            localUpdateProgressBar(progress_bar,completed_steps,total_steps,scenario);
        end
    end
end

% Metrics 7--9: compute the sample standard deviation inside every 100 ms
% collecting window, then average the window standard deviations per UE.
% Unwrap AoA along time before forming the windows.
power_db = 10*log10(power3);
delay_ns = 1e9*delay3;
unwrapped_aoa_deg = rad2deg(unwrap(deg2rad(aoa3),[],1));
varying_power = localMeanWindowStd(power_db,time_nodes,options.window_dt);
varying_delay_ns = localMeanWindowStd(delay_ns,time_nodes,options.window_dt);
varying_aoa_deg = localMeanWindowStd( ...
    unwrapped_aoa_deg,time_nodes,options.window_dt);

result = localEmptyResult();
result.procedure = procedure;
result.scenario = scenario.name;
result.frequency = scenario.frequency;
result.speed_kmh = 30;
result.time_nodes = time_nodes;
result.ue_initial_position = reshape([UE_list.inital_Position],3,[]).';
result.ue_direction_deg = [UE_list.phi_v];
result.ue_min_2d_distance_m = arrayfun(@(ue)localSegmentMinimumDistance( ...
    ue.inital_Position(1:2),ue.velocity*[cosd(ue.phi_v),sind(ue.phi_v)], ...
    options.duration_s,BS.Position(1:2)),UE_list);
result.coupling_loss_db = coupling_loss;
result.sir_db = sir_db;
result.sinr_db = sinr_db;
result.serving_sector_idx = serving_sector_idx;
result.cluster3_active_idx = cluster3_active_idx;
result.cluster3_original_id = cluster3_original_id;
result.los_flag = los_flag;
result.initial_k_db = initial_k_db;
result.initial_zsa_deg = initial_zsa_deg;
result.initial_theta_los_zoa_deg = initial_theta_los_zoa_deg;
result.initial_cluster3_power_for_angle = ...
    initial_cluster3_power_for_angle;
result.cluster3_power_db = 10*log10(power3);
result.cluster3_power_linear = power3;
result.cluster3_delay_ns = 1e9*delay3;
result.cluster3_aoa_deg = aoa3;
result.cluster3_zoa_deg = zoa3;
result.cluster3_absolute_delay_ns = 1e9*absolute_delay3;
result.avg_varying_power = varying_power;
result.avg_varying_power_db = varying_power; % legacy field name
result.avg_varying_delay_ns = varying_delay_ns;
result.avg_varying_aoa_deg = varying_aoa_deg;
result.metric9_diagnostic = localMetric9Diagnostic(UE_list, ...
    varying_aoa_deg,procedure_a_regenerated_any, ...
    procedure_a_regeneration_count,los_nlos_transition_any, ...
    los_nlos_transition_count);
result.options = options;
result.implementation = struct( ...
    'main_file',which(mfilename), ...
    'channel_file',which('channel.Comm_channel'), ...
    'varying_rate_definition','mean of per-window sample standard deviations', ...
    'varying_metric_normalization','N_minus_1', ...
    'varying_metric_window_s',options.window_dt, ...
    'angle_delay_mode','previous_unnormalized_tau_tilde', ...
    'angle_delay_expression','state.tau_absolute');
result.validation = localBenchmarkValidation(result,options);
fprintf('  Angle delay     %s (%s)\n', ...
    result.implementation.angle_delay_mode, ...
    result.implementation.angle_delay_expression);
localPrintValidation(result.validation,procedure);
localPrintMetric9Diagnostic(result.metric9_diagnostic);
localWarnMissingBenchmark(result.validation,procedure,options.benchmark_xlsx);
end

function [raw_los,raw_nlos,raw_o2i,lsp] = ...
    localProcedureAInitialFields(UE_list,n_cluster)
% Clause 7.6.3.1 spatial fields at t_0.  The cluster-specific variables
% are generated before delay sorting. UMi correlation distances are 12 m
% for LOS and 15 m for NLOS/O2I.
positions = reshape([UE_list.Position],3,[]).';
raw_los = localProcedureAInitialRaw(positions,12,n_cluster);
raw_nlos = localProcedureAInitialRaw(positions,15,n_cluster);
raw_o2i = localProcedureAInitialRaw(positions,15,n_cluster);

max_id = max([UE_list.ID]);
lsp = struct('LOS',nan(7,max_id),'NLOS',nan(6,max_id), ...
    'O2I',nan(6,max_id));
% Clause 7.5 Step 4 requires LSPs of different BS-UT links to be
% uncorrelated. Co-sited sectors still share these values because they
% use the same BS-site/UT link state in localSelectServingSector().
ue_ids = [UE_list.ID];
lsp.LOS(:,ue_ids) = randn(7,numel(ue_ids));
lsp.NLOS(:,ue_ids) = randn(6,numel(ue_ids));
lsp.O2I(:,ue_ids) = randn(6,numel(ue_ids));
end

function table = localProcedureAInitialRaw(positions,d_cor,n_cluster)
gaussian_fields = {'shadow','AOA_jitter','AOD_jitter', ...
    'ZOA_jitter','ZOD_jitter'};
uniform_fields = {'tau','AOA','AOD','ZOA','ZOD'};
table = struct();
for idx = 1:numel(gaussian_fields)
    table.(gaussian_fields{idx}) = localCorrelatedGaussian( ...
        positions(:,1:2),d_cor,n_cluster);
end
for idx = 1:numel(uniform_fields)
    gaussian = localCorrelatedGaussian( ...
        positions(:,1:2),d_cor,n_cluster);
    values = 0.5*(1+erf(gaussian/sqrt(2)));
    table.(uniform_fields{idx}) = min(max(values,eps),1-eps);
end
end

function [serving_link, serving_idx, coupling_loss, sir_db, sinr_db] = ...
    localSelectServingSector(BS, BS_sector_list, UE, scenario, t_k)
num_sector = numel(BS_sector_list);
rsrp_dbm = -inf(1,num_sector);
coupling_db = nan(1,num_sector);
sector_phase = cell(1,num_sector);

% One BS-site/UE link supplies the common Steps 1--10 channel state for
% every co-sited sector.  Only the sector antenna response is regenerated,
% matching Flexible's Link38901 + per-sector RSRP_calc workflow.
BS.sector = BS_sector_list(1);
shared_link = channel.Comm_channel(BS, UE, scenario, true, t_k);
for sector_idx = 1:num_sector
    BS.sector = BS_sector_list(sector_idx);
    % TR 38.901 requires co-sited sectors to share the outcome of Steps
    % 1--9. Step-10 initial phases are sector/link specific.
    shared_link.initial_random_phases();
    sector_phase{sector_idx} = shared_link.PHI_n_m;
    shared_link.generate_channel();
    coupling_db(sector_idx) = localScalarFullCouplingLoss( ...
        shared_link);
    if isfinite(coupling_db(sector_idx))
        rsrp_dbm(sector_idx) = BS.Power + coupling_db(sector_idx);
    end
end

[~, serving_idx] = max(rsrp_dbm);
BS.sector = BS_sector_list(serving_idx);
shared_link.PHI_n_m = sector_phase{serving_idx};
if strcmpi(scenario.spatial_consistency_procedure,'A')
    % The constructor saved the pre-selection Step-10 phase. Replace it
    % immediately with the phase of the selected serving sector so that
    % the first mobility update restores the same realization.
    shared_link.syncProcedureAPhaseState();
end
shared_link.generate_channel();
serving_link = shared_link;
coupling_loss = coupling_db(serving_idx);

interference_idx = true(1,num_sector);
interference_idx(serving_idx) = false;
interference_mw = sum(db2pow(rsrp_dbm(interference_idx)),'omitnan');
sir_db = rsrp_dbm(serving_idx) - 10*log10(interference_mw);
noise_dbm = scenario.N0_dBm + 10*log10(scenario.BW) + ...
    scenario.UE_noise_figure;
sinr_db = rsrp_dbm(serving_idx) - ...
    10*log10(interference_mw + db2pow(noise_dbm));
end

function coupling_db = localScalarFullCouplingLoss(link)
pl_db = mean(link.PL(:),'omitnan');
sf_db = mean(link.SF(:),'omitnan');
if ~(isfinite(pl_db) && isfinite(sf_db))
    coupling_db = nan;
    return;
end

% Flexible's coupling loss uses the port-0 channel gain and normalizes by
% the number of receive elements/polarizations.
if ~isempty(link.H_port0)
    channel_power = sum(link.H_port0(:),'omitnan');
    if isfinite(channel_power) && channel_power > 0
        channel_power_db = pow2db(channel_power);
    else
        channel_power_db = 0;
    end
else
    % Fallback for links created before H_port0 was available.
    channel_power = sum(link.H_fastfading(:),'omitnan');
    if ~(isfinite(channel_power) && channel_power > 0)
        channel_power_db = 0;
    else
        channel_power_db = pow2db(channel_power);
    end
end

rx_elements = max(1, mean(link.U(:),'omitnan'));
rx_polarizations = max(1, link.RX.antenna_params.panel.P);
coupling_db = -pl_db + sf_db + channel_power_db - ...
    10*log10(rx_elements * rx_polarizations);
coupling_db = double(coupling_db(1));
end

function options = localDefaultOptions(options)
defaults = struct('num_ue',100,'duration_s',1.0,'sample_dt',0.01, ...
    'window_dt',0.1,'rng_seed',7,'plot_results',true,'save_results',true, ...
    'plot_scene',true,'plot_time_series',false,'show_progress',true, ...
    'debug_ue_idx',1, ...
    'output_dir',fullfile(pwd,'results','Config2_ProcA_ProcB'), ...
    'benchmark_xlsx',fullfile(fileparts(mfilename('fullpath')), ...
        'Phase3SpatialConsistency_v15_Ericsson.xlsx'));
names = fieldnames(defaults);
for idx = 1:numel(names)
    name = names{idx};
    if ~isfield(options,name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
end

if options.sample_dt <= 0 || options.window_dt < options.sample_dt
    error('Require 0 < sample_dt <= window_dt.');
end
end

function progress_bar = localCreateProgressBar(options,scenario)
progress_bar = [];
if ~options.show_progress
    return;
end
try
    progress_bar = waitbar(0,localProgressMessage(scenario,0));
catch
    progress_bar = [];
end
end

function localUpdateProgressBar(progress_bar,completed_steps,total_steps,scenario)
if isempty(progress_bar) || ~isgraphics(progress_bar)
    return;
end
fraction = completed_steps/max(total_steps,1);
waitbar(fraction,progress_bar,localProgressMessage(scenario,100*fraction));
end

function message = localProgressMessage(scenario,percent)
message = sprintf('(Config 2, 3D-UMi, %.0f GHz)  %.1f %%', ...
    scenario.frequency/1e9,percent);
end

function localCloseProgressBar(progress_bar)
if ~isempty(progress_bar) && isgraphics(progress_bar)
    delete(progress_bar);
end
end

function localPlotScenario(BS_list, UE_list, scenario, procedure, options)
ue_xy = reshape([UE_list.Position],3,[]).';
bs_xy = reshape([BS_list.Position],3,[]).';

if ismember(scenario.name,{'UMi','UMa','RMa'})
    % Use the scenario's own cell/sector geometry.
    scenario.plot_BS_pos(bs_xy);
    fig = gcf;
    fig.Name = ['Config2 scene - Procedure ' procedure];
    hold on; grid on; axis equal; view(0,90);
else
    figure('Name',['Config2 scene - Procedure ' procedure]);
    hold on; grid on; axis equal;
end

scatter(ue_xy(:,1),ue_xy(:,2),18,[0.85 0.15 0.15],'filled', ...
    'DisplayName','UE');
scatter(bs_xy(:,1),bs_xy(:,2),110,'kp','filled', ...
    'DisplayName','BS site');
for bs_idx = 1:size(bs_xy,1)
    text(bs_xy(bs_idx,1),bs_xy(bs_idx,2), ...
        sprintf('  BS%d',bs_idx),'FontWeight','bold', ...
        'VerticalAlignment','bottom');
end
xlabel('x position (m)');
ylabel('y position (m)');
title(sprintf('Code-generated Config2 scene: %s, %d BS site, %d sectors, %d UE', ...
    scenario.name,size(bs_xy,1),numel(BS_list(1).sector),size(ue_xy,1)));
legend('Location','best');

if options.save_results
    if ~exist(options.output_dir,'dir')
        mkdir(options.output_dir);
    end
    saveas(gcf,fullfile(options.output_dir, ...
        ['Config2_scene_Proc' procedure '.png']));
end
end

function procedures = localNormalizeProcedures(value)
if ischar(value) || isstring(value)
    value = {char(value)};
end
procedures = cellfun(@(x)upper(char(x)),value,'UniformOutput',false);
if any(~ismember(procedures,{'A','B'}))
    error('Procedure must be A, B, or {''A'',''B''}.');
end
end

function [BS_list, BS_sector_list] = localDropSingleBs(scenario, position)
BS = elements.Equipment(scenario, 'BS');
BS.ID = 1;
BS.inital_Position = position;
BS.Position = position;
BS.cluster_wrapped = [0 0];
BS.antenna_params = localBsAntennaParams();
scenario.total_BS_sector_num = 1;

sector = elements.Sector(BS);
sector.ID = [1 1 1];
sector.equipment_type = 'BS';
sector.boresight = 0;
ang.alpha = 0; ang.beta = 0; ang.gamma = 0;
sector.antenna = antennas.antenna_array(BS.antenna_params, ang);
sector.antenna.attachedDevice = sector;
sector.antenna.attachedType = 'BS';
BS.sector = [];
BS_sector_list = [];
for sector_idx = 1:BS.sector_num
    sector = elements.Sector(BS);
    sector.ID = [1 sector_idx sector_idx];
    sector.equipment_type = 'BS';
    sector.boresight = -90 + (360/BS.sector_num)*(sector_idx-1);
    ang.alpha = sector.boresight;
    sector.antenna = antennas.antenna_array(BS.antenna_params, ang);
    sector.antenna.attachedDevice = sector;
    sector.antenna.attachedType = 'BS';
    BS.sector = [BS.sector; sector];
    BS_sector_list = [BS_sector_list; sector]; %#ok<AGROW>
end
BS_list = BS;
end

function [UE_list, ue_pos_list] = localDropOutdoorUmiUes( ...
    BS, BS_sector_list, scenario, duration_s)
% Config2 requires 100% outdoor UTs.  Generate them as outdoor from the
% beginning instead of using Drop_UE_ISAC's normal 80% indoor mixture and
% changing only the Indoor flag afterwards.
UE_list = [];
ue_pos_list = zeros(numel(BS_sector_list)*scenario.UE_per_sec,3);
ue_id = 0;
for sector_idx = 1:numel(BS_sector_list)
    sector = BS_sector_list(sector_idx);
    for ue_in_sector = 1:scenario.UE_per_sec
        ue_id = ue_id + 1;
        UE = elements.Equipment(scenario,'UE');
        UE.ID = ue_id;
        UE.fcin = (scenario.BW/20)*(floor(20*rand)-(19/2)) + ...
            scenario.frequency;
        UE.Indoor = false;
        UE.d_2D_in = 0;
        UE.n_fl = 1;
        UE.O2IPL = 'low';
        UE.O2Isigma = randn;
        UE.carPL = [9,5];
        % The UMi 10 m BS-UT minimum distance is a deployment constraint,
        % so it must hold over the complete Config2 trajectory, not only at
        % t=0. Draw position and direction jointly and reject a trajectory
        % that enters the exclusion region during the one-second run.
        trajectory_is_valid = false;
        while ~trajectory_is_valid
            xy = network_layout.drop_in_hexagonUE( ...
                BS.Position(1:2), scenario.R, scenario.BS_UE_min_d, ...
                sector.boresight, BS.sector_num);
            phi_v = 360*rand();
            velocity_xy = (30/3.6)*[cosd(phi_v),sind(phi_v)];
            min_distance = localSegmentMinimumDistance( ...
                xy,velocity_xy,duration_s,BS.Position(1:2));
            trajectory_is_valid = min_distance >= scenario.BS_UE_min_d;
        end
        UE.height = scenario.UE_height;
        UE.Position = [xy UE.height];
        UE.inital_Position = UE.Position;
        UE.velocity = 30/3.6;
        UE.theta_v = 90;
        UE.phi_v = phi_v;
        UE.rand_LoS = rand(1,scenario.total_BS_sector_num);

        ang.alpha = rand*360-180;
        ang.beta = UE.antenna_params.beta;
        ang.gamma = 0;
        UE.sector.antenna = antennas.antenna_array( ...
            UE.antenna_params,ang);
        UE.sector.antenna.attachedDevice = UE;
        UE.sector.antenna.attachedType = 'UE';

        UE_list = [UE_list; UE]; %#ok<AGROW>
        ue_pos_list(ue_id,:) = UE.Position;
    end
end
end

function min_distance = localSegmentMinimumDistance(start_xy,velocity_xy, ...
    duration_s,center_xy)
relative_start = start_xy-center_xy;
speed_squared = dot(velocity_xy,velocity_xy);
if duration_s <= 0 || speed_squared <= 0
    closest_time = 0;
else
    closest_time = -dot(relative_start,velocity_xy)/speed_squared;
    closest_time = min(max(closest_time,0),duration_s);
end
closest_point = relative_start+closest_time*velocity_xy;
min_distance = norm(closest_point);
end

function params = localBsAntennaParams()
params = struct();
params.array.Mg = 1; params.array.Ng = 2;
params.array.dg_H = 2.5; params.array.dg_V = 2.5;
params.panel.M = 4; params.panel.N = 4;
params.panel.Kv = 4; params.panel.Kh = 4;
params.panel.d_H = 0.5; params.panel.d_V = 0.5;
params.panel.P = 2; params.panel.X_pol = [-45 45];
params.panel.ele_downtilt = 102; params.panel.ele_panning = 0;
params.panel.antenna_model = 'directional';
params.pol_model = 'model-2';
end

function pos = localMoveFromInitial(node,t)
v = node.velocity*[cosd(node.phi_v),sind(node.phi_v),cosd(node.theta_v)];
pos = node.inital_Position + v*t;
end

function [active_idx,original_id] = localThirdActiveClusterIndex(link)
active_idx = nan;
original_id = nan;
% Metrics 7--9 use cluster n = 3 after Step-5 delay sorting and Step-6
% weak-cluster removal, not pre-sort random-draw identity 3.
active_ids = localActiveOriginalIds(link);
if numel(link.Pn) >= 3 && numel(active_ids) >= 3
    active_idx = 3;
    original_id = active_ids(3);
end
end

function active_idx = localTrackedClusterIndex(link,original_id)
active_idx = nan;
if ~(isscalar(original_id) && isfinite(original_id))
    return;
end
active_ids = localActiveOriginalIds(link);
match = find(active_ids == original_id,1,'first');
if ~isempty(match) && match <= numel(link.Pn)
    active_idx = match;
end
end

function active_ids = localActiveOriginalIds(link)
active_ids = link.tau_order(:).';
keep = logical(link.keep(:).');
if numel(keep) == numel(active_ids)
    active_ids = active_ids(keep);
end
if numel(active_ids) > numel(link.Pn)
    active_ids = active_ids(1:numel(link.Pn));
end
end

function mean_window_std = localMeanWindowStd(samples,time_nodes,window_dt)
tol = max(1e-12,eps(max(abs(time_nodes)))*16);
window_starts = time_nodes(1):window_dt:(time_nodes(end)-window_dt+tol);
if isempty(window_starts)
    error('Metric window_dt exceeds the simulated duration.');
end

window_std = nan(numel(window_starts),size(samples,2));
for window_idx = 1:numel(window_starts)
    window_start = window_starts(window_idx);
    in_window = time_nodes >= window_start-tol & ...
        time_nodes <= window_start+window_dt+tol;
    if nnz(in_window) < 2
        error('Each metric window must contain at least two time samples.');
    end
    window_std(window_idx,:) = std(samples(in_window,:),0,1);
end

mean_window_std = mean(window_std,1,'omitnan');
% A UE contributes only when the same tracked cluster is present for the
% complete trajectory and every collecting window has a valid deviation.
incomplete = any(~isfinite(samples),1) | any(~isfinite(window_std),1);
mean_window_std(incomplete) = nan;
end

function diagnostic = localMetric9Diagnostic(UE_list,aoa_metric,regenerated_any, ...
        regeneration_count,transition_any,transition_count)
ue_id = [UE_list.ID];
valid = isfinite(aoa_metric);
valid_idx = find(valid);
[~,order] = sort(aoa_metric(valid_idx),'descend');
top_idx = valid_idx(order(1:min(10,numel(order))));
diagnostic = struct( ...
    'ue_id',ue_id, ...
    'procedureA_regenerated_any',logical(regenerated_any), ...
    'procedureA_regeneration_count',regeneration_count, ...
    'los_nlos_transition_any',logical(transition_any), ...
    'los_nlos_transition_count',transition_count, ...
    'top_ue_id',ue_id(top_idx), ...
    'top_avg_varying_aoa_deg',aoa_metric(top_idx), ...
    'top_procedureA_regenerated',logical(regenerated_any(top_idx)), ...
    'top_los_nlos_transition',logical(transition_any(top_idx)));
if any(valid & transition_any)
    diagnostic.transition_aoa_quantiles = prctile( ...
        aoa_metric(valid & transition_any),[50 90 99]);
else
    diagnostic.transition_aoa_quantiles = nan(1,3);
end
if any(valid & ~transition_any)
    diagnostic.no_transition_aoa_quantiles = prctile( ...
        aoa_metric(valid & ~transition_any),[50 90 99]);
else
    diagnostic.no_transition_aoa_quantiles = nan(1,3);
end
end

function localPrintMetric9Diagnostic(diagnostic)
fprintf('\nMetric 9 tail diagnostic\n');
fprintf('  UE ID       avg AoA std [deg]   regenerated   LOS<->NLOS\n');
for idx = 1:numel(diagnostic.top_ue_id)
    fprintf('  %-11d %-19.6g %-13d %d\n', ...
        diagnostic.top_ue_id(idx), ...
        diagnostic.top_avg_varying_aoa_deg(idx), ...
        diagnostic.top_procedureA_regenerated(idx), ...
        diagnostic.top_los_nlos_transition(idx));
end
fprintf('  transition    P50/P90/P99 = [%.6g %.6g %.6g] deg\n', ...
    diagnostic.transition_aoa_quantiles);
fprintf('  no transition P50/P90/P99 = [%.6g %.6g %.6g] deg\n', ...
    diagnostic.no_transition_aoa_quantiles);
end

function [p,d,a,z,absolute_delay] = localThirdClusterValues(link, active_idx)
p = nan; d = nan; a = nan; z = nan; absolute_delay = nan;
if ~(isscalar(active_idx) && isfinite(active_idx) && active_idx >= 1)
    return;
end
active_idx = round(active_idx);
% Calibration metric 7 refers to the cluster power after Step 6
% normalization, not the unnormalized pre-normalization weight.
if numel(link.Pn) >= active_idx
    p = link.Pn(active_idx);
end
% Metric 8 uses the channel-model cluster delay (the LOS-scaled excess
% delay for a LOS link), consistent with Step 5 and Step 11.
if link.LOS && ~link.O2I && numel(link.tau_n_LOS) >= active_idx
    d = link.tau_n_LOS(active_idx);
elseif numel(link.tau_n) >= active_idx
    d = link.tau_n(active_idx);
end
% Metric 9 uses the nominal cluster-level AOA directly. Do not derive it
% from a circular mean of the cluster's 20 ray AOAs.
if numel(link.phi_n_AOA_cluster) >= active_idx && ...
        isfinite(link.phi_n_AOA_cluster(active_idx))
    a = link.phi_n_AOA_cluster(active_idx);
end
if size(link.theta_n_m_ZOA,1) >= active_idx
    ray_zoa = link.theta_n_m_ZOA(active_idx,:);
    ray_zoa = ray_zoa(isfinite(ray_zoa));
    if ~isempty(ray_zoa)
        z = mean(ray_zoa);
    end
end
if numel(link.tau_absolute) >= active_idx
    absolute_delay = link.tau_absolute(active_idx);
end
end

function table = localProcedureAXnTable(UE_list,d_cor,n_cluster)
max_id = max([UE_list.ID]);
positions = reshape([UE_list.Position],3,[]).';
g = localCorrelatedGaussian(positions(:,1:2),d_cor,n_cluster);
table = nan(max_id,n_cluster);
table([UE_list.ID],:) = 2*(g >= 0)-1;
end

function [proc_table,lsp_path] = localProcedureBPaths(UE_list,time_nodes)
max_id = max([UE_list.ID]);
n_time = numel(time_nodes);
n_cluster = 20;
empty = localEmptyProcedureBRaw(n_time,n_cluster);
proc_table = repmat(empty,1,max_id);
lsp_path = struct('LOS',nan(n_time,7,max_id), ...
    'NLOS',nan(n_time,6,max_id),'O2I',nan(n_time,6,max_id));
for idx = 1:numel(UE_list)
    ue = UE_list(idx);
    path = ue.inital_Position(1:2) + time_nodes(:)*ue.velocity* ...
        [cosd(ue.phi_v),sind(ue.phi_v)];
    proc_table(ue.ID).tau = localUniformPath(path,15,n_cluster);
    proc_table(ue.ID).AOA = localUniformPath(path,15,n_cluster);
    proc_table(ue.ID).AOD = localUniformPath(path,15,n_cluster);
    proc_table(ue.ID).ZOA = localUniformPath(path,15,n_cluster);
    proc_table(ue.ID).ZOD = localUniformPath(path,15,n_cluster);
    proc_table(ue.ID).shadow = localGaussianPath(path,15,n_cluster);
    lsp_path.LOS(:,:,ue.ID) = localGaussianPath(path,50,7);
    lsp_path.NLOS(:,:,ue.ID) = localGaussianPath(path,50,6);
    lsp_path.O2I(:,:,ue.ID) = localGaussianPath(path,10,6);
end
end

function localAssignProcedureBTimeSlice(BS,procB_path,lsp_path,time_idx)
for id = 1:numel(procB_path)
    raw = procB_path(id);
    raw.tau = raw.tau(time_idx,:); raw.AOA = raw.AOA(time_idx,:);
    raw.AOD = raw.AOD(time_idx,:); raw.ZOA = raw.ZOA(time_idx,:);
    raw.ZOD = raw.ZOD(time_idx,:); raw.shadow = raw.shadow(time_idx,:);
    BS.SC_procB_raw(id) = raw;
end
BS.LSP_raw_LOS = squeeze(lsp_path.LOS(time_idx,:,:));
BS.LSP_raw_NLOS = squeeze(lsp_path.NLOS(time_idx,:,:));
BS.LSP_raw_O2I = squeeze(lsp_path.O2I(time_idx,:,:));
end

function raw = localEmptyProcedureBRaw(n_time,n_cluster)
raw = struct('tau',nan(n_time,n_cluster),'AOA',nan(n_time,n_cluster), ...
    'AOD',nan(n_time,n_cluster),'ZOA',nan(n_time,n_cluster), ...
    'ZOD',nan(n_time,n_cluster),'shadow',nan(n_time,n_cluster));
end

function values = localUniformPath(path,d_cor,n_value)
values = 0.5*(1+erf(localGaussianPath(path,d_cor,n_value)/sqrt(2)));
values = min(max(values,eps),1-eps);
end

function values = localGaussianPath(path,d_cor,n_value)
dx = path(:,1)-path(:,1).'; dy = path(:,2)-path(:,2).';
R = exp(-sqrt(dx.^2+dy.^2)/d_cor) + 1e-10*eye(size(path,1));
values = chol(R,'lower')*randn(size(path,1),n_value);
end

function values = localCorrelatedGaussian(positions,d_cor,n_value)
dx = positions(:,1)-positions(:,1).'; dy = positions(:,2)-positions(:,2).';
R = exp(-sqrt(dx.^2+dy.^2)/d_cor) + 1e-10*eye(size(positions,1));
values = chol(R,'lower')*randn(size(positions,1),n_value);
end

function result = localEmptyResult()
result = struct('procedure','','scenario','','frequency',nan,'speed_kmh',nan, ...
    'time_nodes',[],'ue_initial_position',[],'ue_direction_deg',[], ...
    'ue_min_2d_distance_m',[],'coupling_loss_db',[],'sir_db',[],'sinr_db',[], ...
    'serving_sector_idx',[],'cluster3_active_idx',[], ...
    'cluster3_original_id',[],'los_flag',[], ...
    'initial_k_db',[],'initial_zsa_deg',[], ...
    'initial_theta_los_zoa_deg',[], ...
    'initial_cluster3_power_for_angle',[], ...
    'cluster3_power_db',[],'cluster3_power_linear',[], ...
    'cluster3_delay_ns',[],'cluster3_aoa_deg',[], ...
    'cluster3_zoa_deg',[],'cluster3_absolute_delay_ns',[], ...
    'avg_varying_power',[],'avg_varying_power_db',[], ...
    'avg_varying_delay_ns',[],'avg_varying_aoa_deg',[], ...
    'metric9_diagnostic',struct(), ...
    'options',struct(),'implementation',struct(),'validation',struct([]));
end

function localPlotResults(results, options)
procedure_names = {results.procedure};
if numel(procedure_names) == 1
    figure_title = ['Config2 Procedure ' procedure_names{1} ' calibration'];
    output_suffix = ['Proc' procedure_names{1}];
else
    figure_title = ['Config2 Procedures ' strjoin(procedure_names, ' and ') ...
        ' calibration'];
    output_suffix = ['Proc' strjoin(procedure_names, '_Proc')];
end
fig = figure('Name',figure_title);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
fields = {'coupling_loss_db','sinr_db','avg_varying_power', ...
    'avg_varying_delay_ns','avg_varying_aoa_deg'};
labels = {'coupling loss','Wideband SINR (dB)', ...
    'avg varying rate -- power (dB)', ...
    'avg varying rate -- delay (ns)', ...
    'avg varying rate -- AoA (degree)'};
for metric_idx = 1:numel(fields)
    nexttile; hold on; grid on;
    for result_idx = 1:numel(results)
        benchmark = localLoadBenchmark(options.benchmark_xlsx, ...
            results(result_idx).procedure, metric_idx);
        if benchmark.available
            for curve_idx = 1:size(benchmark.curves,2)
                plot(benchmark.curves(:,curve_idx), benchmark.cdf, ...
                    'Color',[0.75 0.75 0.75], 'LineWidth',0.7, ...
                    'HandleVisibility','off');
            end
            plot(benchmark.mean_curve, benchmark.cdf, 'k-', ...
                'LineWidth',2.0, 'DisplayName','Benchmark mean');
        end
    end
    for result_idx = 1:numel(results)
        x = results(result_idx).(fields{metric_idx});
        x = x(isfinite(x));
        if ~isempty(x)
            % Plot the actual empirical CDF. With many UEs, forcing the
            % simulation onto only 101 integer-percentile points makes the
            % segment from P99 to one rare maximum look like a 1% tail.
            x = sort(x(:));
            cdf_percent = (1:numel(x)).'/numel(x)*100;
            if strcmp(results(result_idx).procedure,'A')
                procedure_color = [0.85 0.00 0.00];
            else
                procedure_color = [0.00 0.35 0.85];
            end
            plot(x,cdf_percent,'LineWidth',1.2, ...
                'Color',procedure_color, ...
                'DisplayName','Config2 simulation');
        end
    end
    xlabel(labels{metric_idx}); ylabel('CDF (%)'); legend('Location','best');
    if metric_idx == 3
        power_limits = xlim(gca);
        power_axis_max = max(1, ceil(power_limits(2)));
        xlim([0 power_axis_max]);
        xticks(0:power_axis_max);
    end
end
if options.save_results
    if ~exist(options.output_dir,'dir')
        mkdir(options.output_dir);
    end
    exportgraphics(fig,fullfile(options.output_dir, ...
        ['Config2_spatial_consistency_calibration_' output_suffix '.png']), ...
        'Resolution',180);
end
end

function validation = localBenchmarkValidation(result,options)
fields = {'coupling_loss_db','sinr_db','avg_varying_power', ...
    'avg_varying_delay_ns','avg_varying_aoa_deg'};
names = {'coupling loss','Wideband SINR','varying power', ...
    'varying delay','varying AoA'};
validation = repmat(struct('metric','','sim_quantiles',nan(1,3), ...
    'benchmark_quantiles',nan(1,3),'median_error',nan, ...
    'sim_p99',nan,'sim_max',nan,'benchmark_p99',nan,'benchmark_max',nan, ...
    'benchmark_available',false),1,numel(fields));
for metric_idx = 1:numel(fields)
    benchmark = localLoadBenchmark(options.benchmark_xlsx, ...
        result.procedure,metric_idx);
    sim_values = result.(fields{metric_idx});
    sim_values = sim_values(isfinite(sim_values));
    validation(metric_idx).metric = names{metric_idx};
    validation(metric_idx).benchmark_available = benchmark.available;
    if ~isempty(sim_values)
        validation(metric_idx).sim_quantiles = ...
            prctile(sim_values,[10 50 90]);
        validation(metric_idx).sim_p99 = prctile(sim_values,99);
        validation(metric_idx).sim_max = max(sim_values);
    end
    if benchmark.available
        validation(metric_idx).benchmark_quantiles = interp1( ...
            benchmark.cdf,benchmark.mean_curve,[10 50 90], ...
            'linear','extrap');
        validation(metric_idx).median_error = ...
            validation(metric_idx).sim_quantiles(2) - ...
            validation(metric_idx).benchmark_quantiles(2);
        validation(metric_idx).benchmark_p99 = interp1( ...
            benchmark.cdf,benchmark.mean_curve,99,'linear','extrap');
        validation(metric_idx).benchmark_max = benchmark.mean_curve(end);
    end
end
aoa_idx = find(strcmp({validation.metric},'varying AoA'),1);
if ~isempty(aoa_idx) && validation(aoa_idx).benchmark_available
    fprintf(['  AoA tail        sim P99/max [%9.4g %9.4g]  ' ...
        'benchmark P99/end [%9.4g %9.4g]\n'], ...
        validation(aoa_idx).sim_p99,validation(aoa_idx).sim_max, ...
        validation(aoa_idx).benchmark_p99,validation(aoa_idx).benchmark_max);
end
end

function localWarnMissingBenchmark(validation,procedure,xlsx_path)
if all([validation.benchmark_available])
    return;
end
missing = {validation(~[validation.benchmark_available]).metric};
warning('Config2:MissingBenchmarkData', ...
    ['No numeric Procedure %s benchmark data for: %s. ' ...
     'Checked sheet Config2-Proc%s in %s. Simulation curves are still plotted.'], ...
    procedure,strjoin(missing,', '),procedure,xlsx_path);
end

function localPrintValidation(validation,procedure)
fprintf('\nConfig2 Procedure %s benchmark check (P10 / P50 / P90)\n', ...
    procedure);
for idx = 1:numel(validation)
    fprintf('  %-16s sim [%9.4g %9.4g %9.4g]  benchmark [%9.4g %9.4g %9.4g]\n', ...
        validation(idx).metric,validation(idx).sim_quantiles, ...
        validation(idx).benchmark_quantiles);
end
end

function localPlotTimeSeries(results, options)
% Show the third-cluster state for one UE over the simulated trajectory.
ue_idx = max(1, round(options.debug_ue_idx));
figure('Name',sprintf('Config2 third-cluster diagnostic - UE %d',ue_idx));
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for result_idx = 1:numel(results)
    result = results(result_idx);
    if isempty(result.time_nodes) || ue_idx > size(result.cluster3_power_linear,2)
        continue;
    end
    t_ms = 1000*result.time_nodes(:);
    name = ['Procedure ' result.procedure];

    nexttile(1); hold on; grid on;
    plot(t_ms,result.cluster3_power_linear(:,ue_idx), ...
        'DisplayName',name,'LineWidth',1.1);
    ylabel('P_n (cluster 3)');

    nexttile(2); hold on; grid on;
    plot(t_ms,result.cluster3_delay_ns(:,ue_idx), ...
        'DisplayName',name,'LineWidth',1.1);
    ylabel('delay (ns)');

    nexttile(3); hold on; grid on;
    plot(t_ms,result.cluster3_aoa_deg(:,ue_idx), ...
        'DisplayName',name,'LineWidth',1.1);
    ylabel('AOA (degree)');
    xlabel('time (ms)');
end

nexttile(1); legend('Location','best');
nexttile(2); legend('Location','best');
nexttile(3); legend('Location','best');
end

function benchmark = localLoadBenchmark(xlsx_path, procedure, metric_idx)
benchmark = struct('available',false,'cdf',[],'curves',[],'mean_curve',[]);
if isempty(xlsx_path) || ~isfile(xlsx_path)
    return;
end

% Procedure A and B must use their own contributed calibration sheet.
sheet_name = ['Config2-Proc' upper(procedure)];
% Percentile axis is A29:A129. The workbook's metric blocks are:
% coupling B:AD, SIR AG:BI, power BL:CN, delay CQ:DS, and AOA DV:EX.
blocks = {'B29:AD129','AG29:BI129','BL29:CN129', ...
    'CQ29:DS129','DV29:EX129'};
try
    cdf = readmatrix(xlsx_path,'Sheet',sheet_name,'Range','A29:A129');
    curves = readmatrix(xlsx_path,'Sheet',sheet_name, ...
        'Range',blocks{metric_idx});
catch
    return;
end

valid_rows = isfinite(cdf);
valid_cols = sum(isfinite(curves),1) >= 2;
if ~any(valid_rows) || ~any(valid_cols)
    return;
end
cdf = cdf(valid_rows);
curves = curves(valid_rows,valid_cols);
% Convert the workbook's ps and millidegree values to the Table 7.8-5
% output units ns and degree. Coupling, SIR and power are already in dB.
if metric_idx == 4 || metric_idx == 5
    curves = curves/1e3;
end
benchmark.available = true;
benchmark.cdf = cdf;
benchmark.curves = curves;
benchmark.mean_curve = mean(curves,2,'omitnan');
end
