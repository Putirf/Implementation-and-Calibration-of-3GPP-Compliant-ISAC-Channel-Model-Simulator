function [ST_list, ST_pos_list, sc_info] = drop_ST_sc_calibration(BS_list, scenario, sensing_type, varargin)
%DROP_ST_SC_CALIBRATION Drop one target for ISAC spatial consistency calibration.

[plot_controller, max_group] = localParseInputs(varargin{:});

if strcmp(scenario.name, 'UrbanGrid')
    [st_position, sc_info] = localUrbanGridTargetPosition(scenario, sensing_type.height, max_group);
else
    [st_position, sc_info] = localIndoorTargetPosition(scenario, sensing_type.height);
end

ST_list = elements.Target();
ST_list.ID = 1;
ST_list.type = sensing_type;
ST_list.SP = 1;
ST_list.height = st_position(3);
ST_list.Position = st_position;
ST_list.inital_Position = st_position;
ST_list.rand_LoS = rand(1, numel(BS_list) + 200);
ST_list.is_single_STSP = true;

ST_pos_list = ST_list.Position;

if plot_controller
    figure(1);
    plot3(ST_pos_list(:, 1), ST_pos_list(:, 2), ST_pos_list(:, 3), ...
        'p', 'Color', [0 0.55 0], 'MarkerFaceColor', [0 0.75 0], ...
        'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', 'SC target');
    hold on;
    axis equal;
    view(0, 90);
end
end

function [plot_controller, max_group] = localParseInputs(varargin)
plot_controller = false;
max_group = [];

if nargin >= 1 && ~isempty(varargin{1})
    plot_controller = logical(varargin{1});
end
if nargin >= 2
    max_group = varargin{2};
end
end

function [st_position, sc_info] = localUrbanGridTargetPosition(scenario, height, max_group)
if isempty(max_group)
    max_group = 80;
end

lane_width = scenario.Lanewidth;
sidewalk_width = scenario.Sidewalkwidth;
block_half_width = (scenario.grid_dx - (4*lane_width + 2*sidewalk_width)) / 2;
side_sign = 2 * randi([0, 1]) - 1;

outer_lane_x = side_sign * (block_half_width + sidewalk_width + 1.5*lane_width);
sidewalk_center_x = side_sign * (block_half_width + sidewalk_width/2);

st_y = -scenario.grid_dy/2 + rand * scenario.grid_dy;

st_position = [outer_lane_x, st_y, height];
sc_info = struct();
sc_info.side_sign = side_sign;
sc_info.block_half_width = block_half_width;
sc_info.sidewalk_width = sidewalk_width;
sc_info.sidewalk_center_x = sidewalk_center_x;
sc_info.window_min_y = -scenario.grid_dy/2;
sc_info.window_max_y = scenario.grid_dy/2;
sc_info.max_group = max_group;
sc_info.ue_distribution_length = scenario.grid_dy;
end

function [st_position, sc_info] = localIndoorTargetPosition(scenario, height)
x = scenario.x_range(1) + rand * diff(scenario.x_range);
y = scenario.y_range(1) + rand * diff(scenario.y_range);
st_position = [x, y, height];

sc_info = struct();
sc_info.st_position = st_position;
end
