clc;
clearvars;
close all;

rng(11);

%% UrbanGrid Procedure A example: 1 BS, 2 UE, 2 targets, t = 5 s
simulation_time = 5;
requested_delta_t = 1;
enable_full_channel = true;
plot_controller = false;

scenario = comm_scenario.UrbanGrid(500, 2);
scenario.applyIsacFrequencyPreset('ISAC_FR1');
scenario.spatial_consistency_enable = true;
scenario.spatial_consistency_procedure = 'A';
scenario.UE_per_sec = 2;
scenario.ST_per_cell = 2;
scenario.RP_per_equipment = 0;

sensing_type = sensing_types.Vehicle();

custom_bs_position = [0, 0, scenario.BS_height];
custom_ue_position = [
    55,  4.5, scenario.UE_height;
    78,  4.5, scenario.UE_height];
custom_st_position = [
    95,  0.0, sensing_type.height;
    128, 0.0, sensing_type.height];

[BS_list, ~] = network_layout.Drop_BaseStation(scenario, plot_controller, custom_bs_position);
[UE_list, ~, ~] = network_layout.Drop_UE_ISAC(BS_list, scenario, plot_controller, custom_ue_position);
[ST_list, ~] = network_layout.Drop_ST(BS_list, custom_ue_position, scenario, sensing_type, plot_controller, custom_st_position);

BS = BS_list(1);
num_ue = numel(UE_list);
num_st = numel(ST_list);
for st_idx = 1:num_st
    ST_list(st_idx).ID = num_ue + st_idx;
end
max_peer_id = max([UE_list.ID, ST_list.ID]);

for ue_idx = 1:num_ue
    UE_list(ue_idx).Indoor = false;
    UE_list(ue_idx).d_2D_in = 0;
    UE_list(ue_idx).O2IPL = 'low';
    UE_list(ue_idx).velocity = 3/3.6;
    UE_list(ue_idx).theta_v = 90;
    UE_list(ue_idx).phi_v = 0;
    UE_list(ue_idx).rand_LoS(:) = 0;
end

for st_idx = 1:num_st
    ST_list(st_idx).velocity = 5/3.6;
    ST_list(st_idx).theta_v = 90;
    ST_list(st_idx).phi_v = 90 + 10*(st_idx - 1);
    ST_list(st_idx).rand_LoS(:) = 0;
end

BS.LSP_raw_LOS = randn(7, max_peer_id);
BS.SC_procA_comm_Xn = localCorrelatedXn(reshape([UE_list.Position], 3, []).', max_peer_id, [UE_list.ID]);
BS.SC_procA_target_tx_Xn = localCorrelatedXn(reshape([ST_list.Position], 3, []).', max_peer_id, [ST_list.ID]);
BS.SC_procA_target_rx_Xn = localCorrelatedXn(BS.Position, max_peer_id, [ST_list.ID]);

max_velocity = max([[UE_list.velocity], [ST_list.velocity]]);
delta_t = min(requested_delta_t, 0.99 / max_velocity);
time_nodes = 0:delta_t:simulation_time;
if time_nodes(end) < simulation_time
    time_nodes(end + 1) = simulation_time; %#ok<SAGROW>
end

ue_track = nan(numel(time_nodes), 3, num_ue);
st_track = nan(numel(time_nodes), 3, num_st);
link_list_UE = cell(numel(time_nodes), num_ue);
link_list_ST = cell(numel(time_nodes), num_st);

for time_idx = 1:numel(time_nodes)
    t_k = time_nodes(time_idx);

    for ue_idx = 1:num_ue
        UE_list(ue_idx).Position = localMoveFromInitial(UE_list(ue_idx), t_k);
        UE_list(ue_idx).height = UE_list(ue_idx).Position(3);
        ue_track(time_idx, :, ue_idx) = UE_list(ue_idx).Position;

        link_comm = channel.Comm_channel(BS, UE_list(ue_idx), scenario, enable_full_channel, t_k);
        link_list_UE{time_idx, ue_idx} = link_comm;
    end

    for st_idx = 1:num_st
        ST_list(st_idx).Position = localMoveFromInitial(ST_list(st_idx), t_k);
        ST_list(st_idx).height = ST_list(st_idx).Position(3);
        st_track(time_idx, :, st_idx) = ST_list(st_idx).Position;

        link_target = channel.Target_channel(BS, BS, ST_list(st_idx), scenario, enable_full_channel, t_k);
        link_list_ST{time_idx, st_idx} = link_target;
    end
end

figure('Name', 'UrbanGrid Procedure A 2UE 2ST trajectories');
hold on;
grid on;
plot(BS.Position(1), BS.Position(2), 'k^', 'MarkerSize', 8, 'LineWidth', 1.5);
for ue_idx = 1:num_ue
    plot(ue_track(:,1,ue_idx), ue_track(:,2,ue_idx), 'o-', 'LineWidth', 1.2);
end
for st_idx = 1:num_st
    plot(st_track(:,1,st_idx), st_track(:,2,st_idx), 's-', 'LineWidth', 1.2);
end
axis equal;
xlabel('x (m)');
ylabel('y (m)');
legend('BS', 'UE 1', 'UE 2', 'Target 1', 'Target 2', 'Location', 'best');
title('UrbanGrid Procedure A trajectories, t = 0 to 5 s');

fprintf('Procedure A 2UE/2ST example completed. Final t = %.3f s.\n', time_nodes(end));
fprintf('Generated %d comm channels and %d target channels.\n', numel(link_list_UE), numel(link_list_ST));
fprintf('Time nodes: ');
fprintf('%.3f ', time_nodes);
fprintf('\n');

function table = localCorrelatedXn(positions, max_peer_id, peer_ids)
max_clusters = 20;
d_cor = 50;
if size(positions, 1) == 1 && numel(peer_ids) > 1
    positions = repmat(positions, numel(peer_ids), 1);
end
corr_matrix = eye(size(positions, 1));
for row = 1:size(positions, 1)
    for col = row+1:size(positions, 1)
        target_corr = exp(-norm(positions(row,1:2) - positions(col,1:2)) / d_cor);
        gaussian_corr = sin((pi/2) * target_corr);
        corr_matrix(row, col) = gaussian_corr;
        corr_matrix(col, row) = gaussian_corr;
    end
end
sqrt_corr = chol(corr_matrix + 1e-10*eye(size(corr_matrix)), 'lower');
x_n = ones(size(positions, 1), max_clusters);
x_n(sqrt_corr * randn(size(positions, 1), max_clusters) < 0) = -1;
table = nan(max_peer_id, max_clusters);
for idx = 1:numel(peer_ids)
    table(peer_ids(idx), :) = x_n(idx, :);
end
end

function pos = localMoveFromInitial(node, t_k)
velocity_vec = node.velocity * [cosd(node.phi_v), sind(node.phi_v), cosd(node.theta_v)];
pos = node.inital_Position + velocity_vec * t_k;
end
