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
scenario = comm_scenario.UMa;
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
