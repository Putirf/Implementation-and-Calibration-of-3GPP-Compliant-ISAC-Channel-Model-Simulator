function pos = drop_InH_ISAC(xy_range)
        pos = [(xy_range(1,2)-xy_range(1,1))*rand+xy_range(1,1),(xy_range(2,2)-xy_range(2,1))*rand+xy_range(2,1)];
end