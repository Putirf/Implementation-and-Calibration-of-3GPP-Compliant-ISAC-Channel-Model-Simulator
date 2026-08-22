function pos = drop_in_Indoor_uniform(trp_pos_list, ue_pos_list, xy_range, min_d_trp, min_d_ue)
% drop_in_Indoor_rectReject
% Uniformly drop 1 target in a rectangular area (xy_range),
% then reject and redraw if outside the convex hull of TRP deployment.
%
% INPUTS:
%   trp_pos_list : [Ntrp x 2] TRP (x,y)
%   ue_pos_list  : [Nue  x 2] UE  (x,y), can be []
%   xy_range     : [xmin xmax; ymin ymax]
%   min_d_trp    : min 2D distance to any TRP (default 0)
%   min_d_ue     : min 2D distance to any UE  (default 0)
%
% OUTPUT:
%   pos          : [1 x 2] target position (x,y)

    if nargin < 4 || isempty(min_d_trp), min_d_trp = 0; end
    if nargin < 5 || isempty(min_d_ue),  min_d_ue  = 0; end
    if isempty(ue_pos_list), ue_pos_list = zeros(0,2); end

    x = trp_pos_list(:,1); y = trp_pos_list(:,2);
    k = convhull(x, y);
    hx = x(k); hy = y(k);

    xmin = xy_range(1,1); xmax = xy_range(1,2);
    ymin = xy_range(2,1); ymax = xy_range(2,2);

    if ~(xmax > xmin && ymax > ymin)
        error('xy_range must be [xmin xmax; ymin ymax] with xmax>xmin and ymax>ymin.');
    end

    while true
        Px = xmin + (xmax - xmin) * rand;
        Py = ymin + (ymax - ymin) * rand;
        P  = [Px, Py];

        if ~inpolygon(Px, Py, hx, hy)
            continue;
        end

        if min_d_trp > 0
            d_trp = sqrt(sum((trp_pos_list(:,1:2) - P).^2, 2));
            if any(d_trp < min_d_trp)
                continue;
            end
        end

        if min_d_ue > 0 && ~isempty(ue_pos_list)
            d_ue = sqrt(sum((ue_pos_list(:,1:2) - P).^2, 2));
            if any(d_ue < min_d_ue)
                continue;
            end
        end

        pos = P;
        return;
    end
end
