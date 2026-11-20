function pos = drop_in_hexagonST(center, ue_pos_list, R, min_d, min_d_ue, ...
        height, endpoints_3d, min_distance_3d)
    if nargin < 8
        height = 0;
        endpoints_3d = zeros(0, 3);
        min_distance_3d = 0;
    end
    distance_to_center = -inf;
    distance_to_ue = -inf(size(ue_pos_list, 1), 1);
    distance_3d = -inf(size(endpoints_3d, 1), 1);
    center = center(1:2);
    while distance_to_center < min_d || any(distance_to_ue < min_d_ue) ...
            || any(distance_3d < min_distance_3d)
        hex_basis_a = R*exp(1j*pi/3)*rand;
        hex_basis_b = R*exp(-1j*pi/3)*rand;
        hex_sector_turn = randsrc(1,1,[0, 1, 2]);
        local_hex_point_complex = (hex_basis_a+hex_basis_b)*exp(1j*2*pi/3*hex_sector_turn);
        local_hex_point = ([real(local_hex_point_complex); imag(local_hex_point_complex)])';
        pos = local_hex_point + center;
        distance_to_center = sqrt(sum((pos-center).^2));
        distance_to_ue = sqrt(sum((pos - ue_pos_list(:, 1:2)).^2, 2));
        distance_3d = vecnorm([pos, height] - endpoints_3d, 2, 2);
    end
end

