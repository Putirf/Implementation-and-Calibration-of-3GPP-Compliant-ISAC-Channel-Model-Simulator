function pos = drop_in_hexagonUE(center, R, min_d, boresight,sector_num)
distance_to_center = 0;
center = center(1:2);
if sector_num==1
    sector_rotation = [1,0;0,1];
    sector_center = center;
else
    R = R/sqrt(3);
    sector_center = center + R*[cosd(boresight),sind(boresight)];
    sector_rotation = [cosd(boresight),-sind(boresight);sind(boresight),cosd(boresight)];
    
end

while distance_to_center <= min_d
    hex_basis_a = R*exp(1j*pi/3)*rand;
    hex_basis_b = R*exp(-1j*pi/3)*rand;
    hex_sector_turn = randsrc(1,1,[0, 1, 2]);
    local_hex_point_complex = (hex_basis_a+hex_basis_b)*exp(1j*2*pi/3*hex_sector_turn);
    local_hex_point = (sector_rotation*[real(local_hex_point_complex); imag(local_hex_point_complex)])';
    pos = local_hex_point + sector_center;
    distance_to_center = sqrt(sum((pos-center).^2));
end
end

