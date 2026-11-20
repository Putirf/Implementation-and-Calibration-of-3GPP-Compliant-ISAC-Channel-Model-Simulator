function [UE_list, positions] = Drop_UrbanGridSensingUE(BS_list, scenario, plot_controller, target, center_only)
if nargin < 4, target = []; end
if nargin < 5, center_only = false; end
% UE-monostatic Automotive UTs, TR 38.901 Table 7.9.6.1-3.
ut_type = validatestring(char(scenario.user_case), {'RSU', 'Vehicle'});
if strcmp(ut_type, 'RSU')
    positions = network_layout.drop_RSU_UT(scenario);
    positions(:, 3) = 5;
    headings = zeros(size(positions, 1), 1);
    speed = 0;
else
    % TR 37.885 6.1.2, Urban Option A: type-2 vehicles, 60 km/h.
    speed = 60 / 3.6;
    centers = scenario.Centerposition;
    x_roads = unique([centers(:, 1)-scenario.grid_dx/2; ...
        centers(:, 1)+scenario.grid_dx/2]);
    y_roads = unique([centers(:, 2)-scenario.grid_dy/2; ...
        centers(:, 2)+scenario.grid_dy/2]);
    lane_offsets = [-1.5, -0.5, 0.5, 1.5] * scenario.Lanewidth;
    positions = zeros(0, 3);
    headings = zeros(0, 1);
    for orientation = 1:2
        if orientation == 1
            roads = y_roads;
            limits = [min(x_roads), max(x_roads)];
        else
            roads = x_roads;
            limits = [min(y_roads), max(y_roads)];
        end
        for road = roads(:).'
            for lane = 1:4
                anchor = [];
                if ~isempty(target)
                    target_horizontal = mod(round(target.phi_v/90),2)==0;
                    perpendicular_axis = 1+(orientation==1);
                    if target_horizontal == (orientation==1) ...
                            && abs(target.Position(perpendicular_axis) ...
                            -road-lane_offsets(lane)) < 1e-7
                        anchor = target.Position(3-perpendicular_axis);
                    end
                end
                along = localLanePositions(limits, speed, anchor);
                across = repmat(road+lane_offsets(lane), numel(along), 1);
                if orientation == 1
                    xy = [along, across];
                    heading = 180 * (lane > 2);
                else
                    xy = [across, along];
                    heading = 90 + 180 * (lane <= 2);
                end
                positions = [positions; xy, 1.6*ones(numel(along), 1)]; %#ok<AGROW>
                headings = [headings; repmat(heading, numel(along), 1)]; %#ok<AGROW>
            end
        end
    end
end
if center_only && strcmp(ut_type, 'Vehicle')
    eligible = abs(positions(:,1)) <= scenario.grid_dx/2 ...
        & abs(positions(:,2)) <= scenario.grid_dy/2;
    positions = positions(eligible,:);
    headings = headings(eligible);
end
[UE_list, positions, ~] = network_layout.Drop_UE_ISAC( ...
    BS_list, scenario, plot_controller, positions, 'outdoor');
for k = 1:numel(UE_list)
    UE_list(k).velocity = speed;
    UE_list(k).phi_v = headings(k);
    UE_list(k).theta_v = 90;
    UE_list(k).sector.antenna.attachedType = ['UE-', ut_type];
end
UE_list = tools.configureCalibrationAntennas(UE_list, [90, 0], 'Urban UT');
end

function along = localLanePositions(limits, speed, anchor)
% Center spacing = 5 m body + max(2 m, exponential(mean = speed*2 s)).
% Start with the equilibrium residual, so the road boundary does not
% introduce an artificial concentration of vehicles at its entrance.
mu = speed * 2;
if ~isempty(anchor)
    along = zeros(0,1);
    for direction = [-1,1]
        x = anchor;
        while true
            x = x+direction*(5+max(2,-mu*log(max(realmin,rand))));
            if x < limits(1) || x > limits(2), break; end
            along(end+1,1) = x;
        end
    end
    along = sort(along);
    return;
end
q = exp(-2/mu);
mean_spacing = 7 + mu*q;
if rand < 7*(1-q)/mean_spacing
    first_interval = 7;
elseif rand < 7/(7+mu)
    first_interval = 7 - mu*log(max(realmin, rand));
else
    first_interval = 7 - mu*log(max(realmin, rand)) ...
        - mu*log(max(realmin, rand));
end
x = limits(1) + rand*first_interval;
along = zeros(0, 1);
while x <= limits(2)
    along(end+1, 1) = x; %#ok<AGROW>
    x = x + 5 + max(2, -mu*log(max(realmin, rand)));
end
end
