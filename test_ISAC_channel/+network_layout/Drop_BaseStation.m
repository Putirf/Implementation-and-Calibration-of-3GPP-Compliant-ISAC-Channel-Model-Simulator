function [BS_list,BS_sector_list] = Drop_BaseStation(scenario, varargin)

    [plot_controller, custom_bs_position] = localParseInputs(varargin{:});
    if isempty(custom_bs_position)
        BS_pos_list = scenario.BSposition;
    else
        BS_pos_list = localPosition3D(custom_bs_position, scenario.BS_height, 'custom_bs_position');
    end

    if plot_controller
        scenario.plot_BS_pos(BS_pos_list);
    end


    BS_list = [];
    BS_sector_list = [];
    id     = 0;
    ang.beta  = 0; 
    ang.gamma = 0;
    for pos_idx = 1:size(BS_pos_list,1)
        BS_pos                  = BS_pos_list(pos_idx,:);
        BS                      = elements.Equipment(scenario,'BS'); 
        BS.ID                   = pos_idx;
        BS.inital_Position      = BS_pos;
        BS.Position             = BS_pos;
    
        for n = 1:BS.sector_num
            id                          = id + 1;
            sector                      = elements.Sector(BS);
            sector.ID                   = [pos_idx, n, id];
            sector.equipment_type       = 'BS';
            
            if BS.sector_num == 1            
                sector.boresight        = 0;
            else
                sector.boresight        =(-90+(360/BS.sector_num)*(n-1));
            end
            ang.alpha                   = sector.boresight; 
            sector.antenna              = antennas.antenna_array(BS.antenna_params, ang);
            sector.antenna.attachedDevice = sector;
            sector.antenna.attachedType   = 'BS';
            if ismember(scenario.name,{'InF','InH'})
                sector.boresight          = 0;
                sector.antenna.attachedType   = 'BSInF';
            end
            BS_sector_list              = [BS_sector_list; sector]; %#ok<AGROW>
            BS.sector                   = [BS.sector; sector];
        end
        BS_list = [BS_list;BS]; %#ok<AGROW>
    end
    scenario.total_BS_sector_num = id;
end

function [plot_controller, custom_bs_position] = localParseInputs(varargin)
plot_controller = true;
custom_bs_position = [];

if nargin >= 1 && ~isempty(varargin{1})
    if islogical(varargin{1}) || (isnumeric(varargin{1}) && isscalar(varargin{1}))
        plot_controller = logical(varargin{1});
    else
        custom_bs_position = varargin{1};
    end
end

if nargin >= 2
    custom_bs_position = varargin{2};
end
end

function pos_3d = localPosition3D(pos, default_height, label)
if ~isnumeric(pos) || isempty(pos) || size(pos, 2) < 2 || size(pos, 2) > 3
    error('%s must be an N-by-2 or N-by-3 numeric matrix.', label);
end

if size(pos, 2) == 2
    pos_3d = [pos, default_height * ones(size(pos, 1), 1)];
else
    pos_3d = pos;
end
end
