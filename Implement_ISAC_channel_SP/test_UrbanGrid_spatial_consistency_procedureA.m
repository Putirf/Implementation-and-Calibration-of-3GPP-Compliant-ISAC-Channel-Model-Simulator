clc;
clearvars;
close all;

rng(7);

%% UrbanGrid Procedure A smoke test
simulation_time = 3;       % seconds
requested_delta_t = 1;     % seconds
enable_full_channel = true;
plot_controller = false;

scenario = comm_scenario.UrbanGrid(500, 2);
scenario.applyIsacFrequencyPreset('ISAC_FR1');
scenario.spatial_consistency_enable = true;
scenario.spatial_consistency_procedure = 'A';
scenario.UE_per_sec = 1;
scenario.ST_per_cell = 1;
scenario.RP_per_equipment = 0;

sensing_type = sensing_types.Vehicle();

custom_bs_position = [0, 0, scenario.BS_height];
custom_ue_position = [65, 4.5, scenario.UE_height];
custom_st_position = [90, 0, sensing_type.height];

[BS_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller, custom_bs_position);
[UE_list, ~, ~] = network_layout.Drop_UE_ISAC(BS_list, scenario, plot_controller, custom_ue_position);
[ST_list, ~] = network_layout.Drop_ST(BS_list, custom_ue_position, scenario, sensing_type, plot_controller, custom_st_position);

UE = UE_list(1);
ST = ST_list(1);
BS = BS_list(1);

UE.Indoor = false;
UE.d_2D_in = 0;
UE.O2IPL = 'low';
UE.velocity = 3/3.6;
UE.theta_v = 90;
UE.phi_v = 0;
ST.velocity = 3/3.6;
ST.theta_v = 90;
ST.phi_v = 90;

UE.rand_LoS(:) = 0;
ST.rand_LoS(:) = 0;
BS.LSP_raw_LOS = randn(7, max(UE.ID, ST.ID));

max_velocity = max([UE.velocity, ST.velocity]);
delta_t = min(requested_delta_t, 0.99 / max_velocity);
time_nodes = 0:delta_t:simulation_time;
if time_nodes(end) < simulation_time
    time_nodes(end + 1) = simulation_time; %#ok<SAGROW>
end

comm_delay3 = nan(size(time_nodes));
comm_aoa3 = nan(size(time_nodes));
comm_h = nan(size(time_nodes));
target_delay3 = nan(size(time_nodes));
target_aoa3 = nan(size(time_nodes));
target_h = nan(size(time_nodes));
ue_track = nan(numel(time_nodes), 3);
st_track = nan(numel(time_nodes), 3);

for time_idx = 1:numel(time_nodes)
    t_k = time_nodes(time_idx);
    UE.Position = localMoveFromInitial(UE, t_k);
    UE.height = UE.Position(3);
    ST.Position = localMoveFromInitial(ST, t_k);
    ST.height = ST.Position(3);

    link_comm = channel.Comm_channel(BS, UE, scenario, enable_full_channel, t_k);
    link_target = channel.Target_channel(BS, BS, ST, scenario, enable_full_channel, t_k);

    comm_delay3(time_idx) = localVectorValue(link_comm.tau_n, 3);
    comm_aoa3(time_idx) = localVectorValue(link_comm.phi_prime_AOA, 3);
    comm_h(time_idx) = localFirstMagnitude(link_comm.H_full);

    target_delay3(time_idx) = localMatrixValue(link_target.tau_n_keep, 2, 3);
    target_aoa3(time_idx) = localMatrixValue(link_target.phi_prime_AOA, 2, 3);
    target_h(time_idx) = localFirstMagnitude(link_target.H_full);

    ue_track(time_idx, :) = UE.Position;
    st_track(time_idx, :) = ST.Position;
end

figure('Name', 'UrbanGrid Procedure A trajectories');
plot(ue_track(:,1), ue_track(:,2), 'bo-', 'LineWidth', 1.5);
hold on;
grid on;
plot(st_track(:,1), st_track(:,2), 'rs-', 'LineWidth', 1.5);
plot(BS.Position(1), BS.Position(2), 'k^', 'MarkerSize', 8, 'LineWidth', 1.5);
axis equal;
xlabel('x (m)');
ylabel('y (m)');
legend('UE', 'ST', 'BS', 'Location', 'best');
title('UrbanGrid Procedure A moving nodes');

figure('Name', 'UrbanGrid Procedure A channel update');
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(time_nodes, comm_delay3, 'bo-', 'LineWidth', 1.2);
grid on;
xlabel('t_k (s)');
ylabel('delay (s)');
title('Comm delay cluster 3');

nexttile;
plot(time_nodes, comm_aoa3, 'bo-', 'LineWidth', 1.2);
grid on;
xlabel('t_k (s)');
ylabel('AoA offset (deg)');
title('Comm raw AoA cluster 3');

nexttile;
plot(time_nodes, comm_h, 'bo-', 'LineWidth', 1.2);
grid on;
xlabel('t_k (s)');
ylabel('|H|');
title('Comm first channel coeff');

nexttile;
plot(time_nodes, target_delay3, 'rs-', 'LineWidth', 1.2);
grid on;
xlabel('t_k (s)');
ylabel('delay (s)');
title('Target delay cluster 3');

nexttile;
plot(time_nodes, target_aoa3, 'rs-', 'LineWidth', 1.2);
grid on;
xlabel('t_k (s)');
ylabel('AoA offset (deg)');
title('Target raw AoA cluster 3');

nexttile;
plot(time_nodes, target_h, 'rs-', 'LineWidth', 1.2);
grid on;
xlabel('t_k (s)');
ylabel('|H|');
title('Target first channel coeff');

fprintf('Procedure A test completed. time nodes: ');
fprintf('%.3f ', time_nodes);
fprintf('\n');

function pos = localMoveFromInitial(node, t_k)
velocity_vec = node.velocity * [cosd(node.phi_v), sind(node.phi_v), cosd(node.theta_v)];
pos = node.inital_Position + velocity_vec * t_k;
end

function value = localVectorValue(values, idx)
value = nan;
if ~isempty(values) && numel(values) >= idx
    value = values(idx);
end
end

function value = localMatrixValue(values, row, col)
value = nan;
if ~isempty(values) && size(values, 1) >= row && size(values, 2) >= col
    value = values(row, col);
end
end

function value = localFirstMagnitude(values)
value = nan;
if ~isempty(values)
    value = abs(values(1));
end
end
