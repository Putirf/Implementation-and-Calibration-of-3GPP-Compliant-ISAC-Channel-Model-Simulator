function pos = drop_in_hexagonRP(center, uk,rp_angle,a_d,b_d,c_d)
rp_angle = rp_angle + 120 * (uk-1);
yt = (gamrnd(a_d, 1/b_d,1,1) + c_d) *[cosd(rp_angle),sind(rp_angle)];
pos = yt + center(1:2);
end

