function RSU_pos = drop_RSU_UT(scenario)
% drop_RSU_UT
% RSU at all vertices of every grid.
%
% Inputs:
%   scenario: must contain
%       - scenario.ISD
%       - scenario.Centerposition  (Gx2 or Gx3), each row is a grid center [x y (z)]
%
% Output:
%   RSU_pos: (4*G) x 3, each grid contributes 4 vertices (z=0)


    C = scenario.Centerposition;
    if size(C,2) < 2
        error('scenario.Centerposition must be Gx2 or Gx3.');
    end
    C = C(:,1:2);   % only x,y

    dx = scenario.ISD;           % short side
    dy = scenario.ISD * 433/250; % long side (left-right sides are the long edges)

    hx = dx/2;
    hy = dy/2;

    G = size(C,1);
    RSU_pos = zeros(4*G, 3);

    % Vertex offsets (counter-clockwise): bottom-left, bottom-right, top-right, top-left
    offs = [-hx -hy;
             hx -hy;
             hx  hy;
            -hx  hy];

    for g = 1:G
        base = (g-1)*4;
        verts_xy = C(g,:) + offs;           % 4x2
        RSU_pos(base+1:base+4, 1:2) = verts_xy;
        RSU_pos(base+1:base+4, 3)   = 0;    % RSU height if you want, change here
    end

    % optional: remove duplicates (shared vertices between adjacent grids)
    RSU_pos = unique(round(RSU_pos, 10), 'rows');  % keep if you want unique RSUs
end