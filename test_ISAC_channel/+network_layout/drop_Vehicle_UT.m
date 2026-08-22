function [veh_pos, side, lane_id, s_new, d_min] = drop_Vehicle_UT(scenario, vehicle_list, speed, center_xy, lane_id,isST)
% drop_vehicle_UT (gap-based placement, leave-last-slot rule)
%
% Place ONE vehicle on a rectangular ring road around centre region, with 2 lanes.
% Vehicle is placed uniformly on feasible gaps of the chosen lane, enforcing headway:
%   circular distance to nearest neighbors on same lane >= d_min
% with:
%   d_min = max(2, Exp(mean = speed*2))
%
% Strategy (no rejection sampling over whole perimeter):
%   1) Compute gaps between existing vehicles (including wrap-around)
%   2) If any gap >= 2*d_min, feasible interval exists: [s_i+d_min, s_{i+1}-d_min]
%   3) Choose among feasible intervals with probability proportional to interval length
%      and sample uniformly inside it.
%   4) If lane has no feasible gap, try the other lane.
%   5) If both lanes have no feasible gap, return empty (keep last slot empty).
%
% Inputs:
%   scenario     : has fields ISD, Lanewidth, Sidewalkwidth
%   vehicle_list : struct('s',[],'lane_id',[],'pos',[])  (fields can be empty)
%   speed        : scalar m/s
%   center_xy    : [cx cy] (default [0 0])
%   lane_id      : optional preferred lane (1/2). If empty -> random lane.
%
% Outputs:
%   veh_pos : 1x3 [x y z] or [] if cannot place (leave last slot)
%   side    : 1..4 (1=up,2=right,3=down,4=left) or []
%   lane_id : 1 or 2 (final lane used) or []
%   s_new   : arc-length coordinate in [0,P) or []
%   d_min   : sampled minimum headway used (always returned unless error in input)

    if nargin < 4 || isempty(center_xy), center_xy = [0 0]; end
    if nargin < 5, lane_id = []; end
    if isempty(speed) || ~isscalar(speed) || speed < 0
        error('speed must be a non-negative scalar (m/s).');
    end

    % ---- choose lane preference ----
    if isempty(lane_id)
        lane_id = randi(2);
    end
    if ~(lane_id==1 || lane_id==2)
        error('lane_id must be 1 or 2.');
    end

    % ---- geometry (centre-region size) ----
    % TEMP UrbanGrid calibration: vehicle UE road width is independent from TRP ISD.
    % Original state: ISD = scenario.ISD;
    ISD   = scenario.grid_dx;
    laneW = scenario.Lanewidth;

    % dy = ISD * 433/250;
    dy = scenario.grid_dy;

    roadTerm = 4*laneW;
    Lx = ISD - roadTerm;
    Ly = dy  - roadTerm;
    if Lx <= 0 || Ly <= 0
        error('Invalid center-region size: Lx=%.3f, Ly=%.3f. Check ISD/laneW/swW.', Lx, Ly);
    end

    cx = center_xy(1); cy = center_xy(2);

    % ---- headway threshold ----
    meanD = speed * 2;           % mean distance (m)
    d_min = max(2, localExprnd(meanD));

    % ---- parse vehicle_list ----
    [s_exist, lane_exist] = localParseVehicleList(vehicle_list);

    % ---- try preferred lane; if fail, try other lane; if fail, return empty ----
    try_order = [lane_id, 3-lane_id];

    % ---------- keepLast policy (only centre grid, and not ST) ----------
    keepLast = (~isST) && isequal(center_xy(:).', [0 0]);

    if keepLast
        s1 = s_exist(lane_exist == 1);
        s2 = s_exist(lane_exist == 2);

        P1 = localPerimeterForLane(1, cx, cy, Lx, Ly, laneW);  % helper below
        P2 = localPerimeterForLane(2, cx, cy, Lx, Ly, laneW);

        nF1 = localCountFeasibleGaps(s1, P1, d_min);
        nF2 = localCountFeasibleGaps(s2, P2, d_min);

        if (nF1 + nF2) <= 1
            veh_pos = [];
            side = [];
            lane_id = [];
            s_new = [];
            warning('Center grid: keep last slot empty. Stop at n-1 vehicles.');
            return;
        end
    end

    % ---------- try place in preferred lane, then the other lane ----------
    for lane_try = try_order(:).'
        [xL, xR, yB, yT, P] = localLaneGeom(lane_try, cx, cy, Lx, Ly, laneW);

        s_lane = s_exist(lane_exist == lane_try);

        [ok, s_new_try] = localTryPlaceInLane(s_lane, P, d_min);
        if ok
            lane_id = lane_try;
            s_new   = s_new_try;
            [x, y, side] = localPerimPointAndSide(xL, xR, yB, yT, s_new);
            veh_pos = [x, y, 1.5];
            return;
        end
    end

    % ---- both lanes fail: keep last slot empty ----
    veh_pos = [];
    side = [];
    lane_id = [];
    s_new = [];
    warning('drop_vehicle_UT:NoSlot', ...
        'Both lanes have no available gap >= 2*d_min (d_min=%.3f). Keep last slot empty (max n-1 vehicles).', d_min);
end

% ================= helpers =================

function r = localExprnd(mu)
% exponential random variable with mean mu
    if mu <= 0
        r = 0;
        return;
    end
    u = rand;
    r = -mu * log(max(u, realmin));
end

function [s, lane_id] = localParseVehicleList(vehicle_list)
% Support your struct('s',[],'lane_id',[],'pos',[])
    if isempty(vehicle_list)
        s = zeros(0,1);
        lane_id = zeros(0,1);
        return;
    end
    if ~isstruct(vehicle_list) || ~isfield(vehicle_list,'s') || ~isfield(vehicle_list,'lane_id')
        error('vehicle_list must be a struct with fields: s and lane_id.');
    end
    s = vehicle_list.s(:);
    lane_id = vehicle_list.lane_id(:);
    if numel(s) ~= numel(lane_id)
        error('vehicle_list.s and vehicle_list.lane_id must have same length.');
    end
    if any(~(lane_id==1 | lane_id==2))
        error('vehicle_list.lane_id must contain only 1 or 2.');
    end
end

function [ok, s_new] = localTryPlaceInLane(s_lane, P, d_min)
% Gap-based placement on a circular lane of length P.
% Existing vehicles at s_lane (any order), in [0,P).
% A gap between consecutive vehicles must be >= 2*d_min to fit a new vehicle
% such that distances to both neighbors >= d_min.
%
% Choose among feasible gaps proportional to free length (gap - 2*d_min),
% then sample uniformly inside [s_i+d_min, s_{i+1}-d_min] (with wrap).

    s_new = [];

    if isempty(s_lane)
        ok = true;
        s_new = rand * P;
        return;
    end

    s = sort(mod(s_lane(:), P));
    K = numel(s);

    % gaps including wrap-around
    gap = zeros(K,1);
    gap(1:K-1) = s(2:K) - s(1:K-1);
    gap(K)     = (s(1) + P) - s(K);

    feasible = gap >= 2*d_min;
    if ~any(feasible)
        ok = false;
        return;
    end

    freeLen = zeros(K,1);
    freeLen(feasible) = gap(feasible) - 2*d_min;

    totalFree = sum(freeLen);
    if totalFree <= 0
        ok = false;
        return;
    end

    r = rand * totalFree;
    cum = cumsum(freeLen);
    idx = find(cum >= r, 1, 'first');

    s_start = s(idx);
    if idx < K
        s_end = s(idx+1);
    else
        s_end = s(1) + P;
    end

    a = s_start + d_min;
    b = s_end   - d_min;

    % numerical safety
    if b < a
        ok = false;
        return;
    end

    s_pick = a + (b-a) * rand;
    s_new = mod(s_pick, P);
    ok = true;
end

function [x, y, side] = localPerimPointAndSide(xL, xR, yB, yT, s)
% Clockwise, start at bottom-left (xL,yB):
% bottom -> right -> top -> left
% side: 1=up(top), 2=right, 3=down(bottom), 4=left

    Lx = xR - xL;
    Ly = yT - yB;

    if s <= Lx
        x = xL + s;  y = yB;  side = 3;
    elseif s <= Lx + Ly
        x = xR;      y = yB + (s - Lx);  side = 2;
    elseif s <= 2*Lx + Ly
        x = xR - (s - (Lx + Ly));  y = yT;  side = 1;
    else
        x = xL;      y = yT - (s - (2*Lx + Ly));  side = 4;
    end
end

function nFeasible = localCountFeasibleGaps(s_lane, P, d_min)
    if isempty(s_lane)
        nFeasible = 1;
        return;
    end
    s = sort(mod(s_lane(:), P));
    K = numel(s);
    gap = zeros(K,1);
    gap(1:K-1) = s(2:K) - s(1:K-1);
    gap(K)     = (s(1) + P) - s(K);
    nFeasible = nnz(gap >= 2*d_min);
end

function P = localPerimeterForLane(lane_try, cx, cy, Lx, Ly, laneW)
    [~, ~, ~, ~, P] = localLaneGeom(lane_try, cx, cy, Lx, Ly, laneW);
end

function [xL, xR, yB, yT, P] = localLaneGeom(lane_try, cx, cy, Lx, Ly, laneW)
    offset = (lane_try - 0.5) * laneW;

    xL = cx - Lx/2 - offset;
    xR = cx + Lx/2 + offset;
    yB = cy - Ly/2 - offset;
    yT = cy + Ly/2 + offset;

    P  = 2 * ((xR - xL) + (yT - yB));
end
