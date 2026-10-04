function [UE_list, ue_pos_list] = drop_UE_ISAC_sc_calibration(BS_list, scenario, ST_list, sc_info, varargin)
%DROP_UE_ISAC_SC_CALIBRATION Drop UTs for ISAC spatial consistency calibration.

[plot_controller, num_ue] = localParseInputs(varargin{:});

if strcmp(scenario.name, 'UrbanGrid')
    ue_pos_list = localUrbanGridUePositions(scenario, ST_list, sc_info, num_ue);
elseif strcmp(scenario.name, 'InH')
    ue_pos_list = localIndoorUePositions(BS_list, scenario, num_ue);
else
    ue_pos_list = localUniformScenarioPositions(scenario, num_ue);
end

UE_list = [];
for ue_idx = 1:size(ue_pos_list, 1)
    ue = elements.Equipment(scenario, 'UE');
    ue.ID = numel(BS_list) + ue_idx;
    ue.fcin = scenario.frequency;
    ue.d_2D_in = localIndoorDistance(scenario, ue_pos_list(ue_idx, :));
    ue.Position = ue_pos_list(ue_idx, :);
    ue.inital_Position = ue.Position;
    ue.height = ue.Position(3);
    ue.rand_LoS = rand(1, numel(BS_list) + ue_idx + 19);
    ue.n_fl = 1;
    ue.O2Isigma = randn;
    ue.carPL = [0, 0];

    ang.alpha = rand * 360 - 180;
    ang.beta = ue.antenna_params.beta;
    ang.gamma = 0;
    ue.sector.antenna = antennas.antenna_array(ue.antenna_params, ang);
    ue.sector.antenna.attachedDevice = ue;
    ue.sector.antenna.attachedType = 'UE';

    UE_list = [UE_list; ue]; %#ok<AGROW>
end

if plot_controller
    figure(1);
    plot3(ue_pos_list(:, 1), ue_pos_list(:, 2), ue_pos_list(:, 3), ...
        '.', 'Color', [0 0.35 0.9], 'MarkerSize', 12, 'DisplayName', 'SC UE');
    hold on;
    axis equal;
    view(0, 90);
end
end

function [plot_controller, num_ue] = localParseInputs(varargin)
plot_controller = false;
num_ue = 120;

if nargin >= 1 && ~isempty(varargin{1})
    plot_controller = logical(varargin{1});
end
if nargin >= 2 && ~isempty(varargin{2})
    num_ue = varargin{2};
end
end

function ue_pos_list = localUrbanGridUePositions(scenario, ST_list, sc_info, num_ue)
required_fields = {'side_sign', 'block_half_width', 'window_min_y', 'window_max_y'};
for field_idx = 1:numel(required_fields)
    if ~isfield(sc_info, required_fields{field_idx})
        error('UrbanGrid sc_info must contain %s.', required_fields{field_idx});
    end
end

% Table 7.9.6.3-1: UTs are in the same-side sidewalk with a uniform
% lateral distance in [0,1] m and longitudinally span the 433 m center grid.
ue_lateral_offset = rand(num_ue, 1);
ue_y = sc_info.window_min_y + rand(num_ue, 1) * (sc_info.window_max_y - sc_info.window_min_y);
ue_x = sc_info.side_sign * (sc_info.block_half_width + ue_lateral_offset);
ue_pos_list = [ue_x, ...
    ue_y, scenario.UE_height * ones(num_ue, 1)];

if isempty(ST_list) || isempty(ST_list.Position)
    return;
end
end

function ue_pos_list = localIndoorUePositions(BS_list, scenario, num_ue)
% Table 7.9.6.3-1: drop 10 UTs uniformly per cell (120 UTs total).
% For other num_ue values, distribute UTs as evenly as possible so that
% reduced-size diagnostic runs still cover all indoor cells.
bs_positions = reshape([BS_list.Position], 3, []).';
x_centers = unique(bs_positions(:, 1), 'sorted');
y_centers = unique(bs_positions(:, 2), 'sorted');
num_x_cells = numel(x_centers);
num_y_cells = numel(y_centers);
num_cells = num_x_cells * num_y_cells;

if num_cells ~= numel(BS_list)
    error('InH calibration requires a rectangular TRP grid with one TRP per cell.');
end

x_edges = [scenario.x_range(1); ...
    (x_centers(1:end-1) + x_centers(2:end)) / 2; scenario.x_range(2)];
y_edges = [scenario.y_range(1); ...
    (y_centers(1:end-1) + y_centers(2:end)) / 2; scenario.y_range(2)];

cell_counts = floor(num_ue / num_cells) * ones(num_cells, 1);
num_remainder = mod(num_ue, num_cells);
if num_remainder > 0
    remainder_cells = randperm(num_cells, num_remainder);
    cell_counts(remainder_cells) = cell_counts(remainder_cells) + 1;
end

ue_pos_list = nan(num_ue, 3);
ue_idx = 0;
cell_idx = 0;
for x_idx = 1:num_x_cells
    for y_idx = 1:num_y_cells
        cell_idx = cell_idx + 1;
        num_cell_ue = cell_counts(cell_idx);
        if num_cell_ue == 0
            continue;
        end

        row_idx = ue_idx + (1:num_cell_ue);
        ue_pos_list(row_idx, 1) = x_edges(x_idx) + ...
            rand(num_cell_ue, 1) * (x_edges(x_idx + 1) - x_edges(x_idx));
        ue_pos_list(row_idx, 2) = y_edges(y_idx) + ...
            rand(num_cell_ue, 1) * (y_edges(y_idx + 1) - y_edges(y_idx));
        ue_pos_list(row_idx, 3) = scenario.UE_height;
        ue_idx = ue_idx + num_cell_ue;
    end
end
end

function ue_pos_list = localUniformScenarioPositions(scenario, num_ue)
ue_x = scenario.x_range(1) + rand(num_ue, 1) * diff(scenario.x_range);
ue_y = scenario.y_range(1) + rand(num_ue, 1) * diff(scenario.y_range);
ue_pos_list = [ue_x, ue_y, scenario.UE_height * ones(num_ue, 1)];
end

function d2Din = localIndoorDistance(scenario, position)
if strcmp(scenario.name, 'InH')
    d2Din = sqrt(sum(position(1:2).^2));
else
    d2Din = 0;
end
end
