function [veh_pos, side, lane_id] = drop_Vehicle_ST(scenario)

center_xy = [0 0];
lane_id = randi(1);
% lane_id = 1;
% ---- geometry (centre-region size) ----
ISD   = scenario.ISD;
laneW = scenario.Lanewidth;
swW   = scenario.Sidewalkwidth;

dy = ISD * 433/250;

roadTerm = 4*laneW + 2*swW;
Lx = ISD - roadTerm;
Ly = dy  - roadTerm;

cx = center_xy(1); cy = center_xy(2);

offset = (lane_id - 0.5) * laneW + swW;
xL = cx - Lx/2 - offset;
xR = cx + Lx/2 + offset;
yB = cy - Ly/2 - offset;
yT = cy + Ly/2 + offset;
P  = 2 * ((xR - xL) + (yT - yB));

s = rand * P;

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

veh_pos = [x, y];

end
