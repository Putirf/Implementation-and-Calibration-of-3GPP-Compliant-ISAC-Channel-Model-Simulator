function UE_pos = drop_Pedestrian_UT(scenario, N, center_xy)
% drop_Pedestrian_UT
% Return all pedestrian UE positions (16*N UEs).
%
% Inputs:
%   scenario   : with fields ISD, Lanewidth, Sidewalkwidth
%   N          : number of road grids (default = 1)
%   center_xy  : [cx cy] (default = [0 0])
%
% Output:
%   UE_pos     : (16*N) x 3 array [x y z]

    if nargin < 2 || isempty(N), N = 1; end
    if nargin < 3 || isempty(center_xy), center_xy = [0 0]; end

    % ---- geometry ----
    % TEMP UrbanGrid calibration: road-grid width is independent from TRP ISD.
    % Original state: ISD = scenario.ISD;
    ISD   = scenario.grid_dx;
    laneW = scenario.Lanewidth;
    swW   = scenario.Sidewalkwidth;

    % dy = ISD * 433/250;
    dy = scenario.grid_dy;
    % TR 38.901 Table 7.9.6.1-3 uses {(250m-17m)+(433m-17m)}*2*N.
    % With TR 37.885 urban road dimensions, 4*laneW + swW = 4*3.5 + 3 = 17 m.
    roadTerm = 4*laneW + swW;

    Lx = ISD - roadTerm;
    Ly = dy  - roadTerm;
    if Lx <= 0 || Ly <= 0
        error('Invalid geometry: Lx=%.3f, Ly=%.3f.', Lx, Ly);
    end

    % total sidewalk length across N grids
    A = (Lx + Ly) * 2 * N;

    % spacing
    totalUE = 16 * N;
    X = A / totalUE;

    perim1 = 2*(Lx + Ly);   % one grid perimeter

    % ---- generate all positions ----
    UE_pos = zeros(totalUE,3);

    cx = center_xy(1);
    cy = center_xy(2);

    xL = cx - Lx/2;
    xR = cx + Lx/2;
    yB = cy - Ly/2;
    yT = cy + Ly/2;
    randstart = rand * A;
    for k = 1:totalUE
        s  = (k-1) * X + randstart;
        s1 = mod(s, perim1);     % wrap to one grid

        [x, y] = localPerimPoint(xL, xR, yB, yT, Lx, Ly, s1);

        UE_pos(k,:) = [x, y, 1.5];
    end
end


function [x, y] = localPerimPoint(xL, xR, yB, yT, Lx, Ly, s)
% Map s in [0, 2*(Lx+Ly)) to rectangle perimeter (clockwise)
% start at bottom-left corner.

    if s <= Lx
        x = xL + s;  y = yB;                      % bottom
    elseif s <= Lx + Ly
        x = xR;      y = yB + (s - Lx);           % right
    elseif s <= 2*Lx + Ly
        x = xR - (s - (Lx + Ly));  y = yT;        % top
    else
        x = xL;      y = yT - (s - (2*Lx + Ly));  % left
    end
end
