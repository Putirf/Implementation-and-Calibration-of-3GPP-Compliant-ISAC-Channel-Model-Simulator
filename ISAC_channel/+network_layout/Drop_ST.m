function  [ST_list, ST_pos_list] = Drop_ST(BS_list, UE_pos_list, scenario, type, varargin)
ST_list = []; ST_pos_list = [];
id      = 0;
[plot_controller, St_pos, bs_st_min_d_override] = localParseInputs(varargin{:});

BS_ST_min_d = scenario.BS_ST_min_d;
if ~isempty(bs_st_min_d_override)
    BS_ST_min_d = bs_st_min_d_override;
end
UE_ST_min_d = scenario.UE_ST_min_d;

if isempty(St_pos)
    for st = 1:scenario.ST_per_cell
        BS = BS_list(1);
        St      = elements.Target();
        id      = id + 1;
        St.ID   = id;
        St.type = type;
        St.SP = 1;
        if ismember(scenario.name, {'RMa', 'UMa', 'UMi'})
            R = scenario.R;
            ST_pos_2D = network_layout.drop_in_hexagonST(BS.Position(1:2), UE_pos_list, R, BS_ST_min_d, UE_ST_min_d);
        elseif ismember(scenario.name, {'InH', 'InF'})
            bs_pos_list = zeros(numel(BS_list), 3);
            for m = 1:numel(BS_list)
                bs_pos_list(m, :) = BS_list(m).Position;
            end 

            % convex hull or whole area
            % ST_pos_2D = network_layout.drop_in_Indoor_convex(bs_pos_list(:, 1:2), UE_pos_list, [scenario.x_range; scenario.y_range], BS_ST_min_d, UE_ST_min_d);
              ST_pos_2D = network_layout.drop_in_Indoor_uniform(bs_pos_list(:, 1:2), UE_pos_list, [scenario.x_range; scenario.y_range], BS_ST_min_d, UE_ST_min_d);

        elseif strcmp(scenario.name, 'UrbanGrid')
            % Table 7.9.6.1-3: drop the Automotive target uniformly in the
            % centre road grid.  STX/SRX minimum-distance eligibility is a
            % property of each candidate TRP pair and is applied by the
            % TRP-TRP calibration driver after this geometry is generated.
            [ST_pos_2D, side, ~] = network_layout.drop_Vehicle_ST(scenario);
            St.phi_v = (side - 1) * 90;
            St.velocity = 60 / 3.6; % Table 7.9.6.1-3: 60 km/h in all lanes.
            St.theta_v = 90;
            St.SP = 1;
        end
        if strcmp(scenario.name, 'UrbanGrid') ...
                && strcmp(type.sensing_type, 'Vehicle') && St.SP == 1
            % Table 7.9.6.1-3, Note 2: one Vehicle SPST is at (0,0,H/2)
            % in the target LCS.  Vehicle type 2 has H=1.6 m.
            St.height = type.height / 2;
        else
            St.height = type.height;
        end
        ST_pos_3D = [ST_pos_2D, St.height];
        St.inital_Position = ST_pos_3D;
        St.Position = ST_pos_3D;
        St.rand_LoS = rand(1, numel(BS_list) * 19);

        ST_pos_list = [ST_pos_list; St.Position]; %#ok<AGROW>
        ST_list     = [ST_list; St]; %#ok<AGROW>
        St.is_single_STSP = true;
    end
else
    default_st_height = type.height;
    if strcmp(scenario.name, 'UrbanGrid') && strcmp(type.sensing_type, 'Vehicle')
        default_st_height = type.height / 2;
    end
    St_pos = localPosition3D(St_pos, default_st_height, 'custom_st_position');
    for st = 1:size(St_pos, 1)
        St      = elements.Target();
        id      = id + 1;
        St.ID   = id;
        St.type = type;
        St.SP = 1;
        if strcmp(scenario.name, 'UrbanGrid')
            % TEMP UrbanGrid calibration: use the single scattering point option from Table 7.9.6.1-3 note 2.
            % Original temporary state: St.SP = 5;
            St.SP = 1;
        end
        St.height = St_pos(st, 3);
        St.Position = St_pos(st, :);
        St.inital_Position = St.Position;
        St.rand_LoS = rand(1, numel(BS_list) * 19);
        St.is_single_STSP = true;

        ST_pos_list = [ST_pos_list; St.Position]; %#ok<AGROW>
        ST_list     = [ST_list; St]; %#ok<AGROW>
    end
end
if plot_controller
    figure(1);
    for i = 1:numel(BS_list)
        indx = i:numel(BS_list):numel(ST_list);
        plot3(ST_pos_list(indx, 1), ST_pos_list(indx, 2), ST_pos_list(indx, 3), ...
            '*', 'color', 'g', 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        hold on;
    end

end
axis equal; view(0, 90);
end

function [plot_controller, custom_st_position, bs_st_min_d_override] = localParseInputs(varargin)
plot_controller = false;
custom_st_position = [];
bs_st_min_d_override = [];

if nargin >= 1 && ~isempty(varargin{1})
    if islogical(varargin{1}) || (isnumeric(varargin{1}) && isscalar(varargin{1}))
        plot_controller = logical(varargin{1});
    else
        custom_st_position = varargin{1};
    end
end

if nargin >= 2
    custom_st_position = varargin{2};
end
if nargin >= 3 && ~isempty(varargin{3})
    bs_st_min_d_override = varargin{3};
    if ~isnumeric(bs_st_min_d_override) || ~isscalar(bs_st_min_d_override) || ...
            ~isreal(bs_st_min_d_override) || ~isfinite(bs_st_min_d_override) || ...
            bs_st_min_d_override < 0
        error('BS-ST minimum-distance override must be a nonnegative finite scalar.');
    end
end
end

function pos_3d = localPosition3D(pos, default_height, label)
if ~isnumeric(pos) || isempty(pos) || size(pos, 2) < 2 || size(pos, 2) > 3
    error('%s must be an N-by-2 or N-by-3 numeric matrix.', label);
end

if size(pos, 2) == 2
    pos_3d = [pos, default_height * ones(size(pos, 1), 1)];
else
    pos_3d = pos;
end
end
