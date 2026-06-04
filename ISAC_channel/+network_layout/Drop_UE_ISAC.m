function [UE_list, ue_pos_list,vehicle_list_center] = Drop_UE_ISAC(BS_list, scenario,varargin)
UE_list = []; ue_pos_list = [];vehicle_list_center = [];
lowhigh = {'low','high'};
id      = 0;
indoor = false;
[plot_controller, Ue_pos] = localParseInputs(varargin{:});

fc      = scenario.frequency;
BW      = scenario.BW;
min_d   = scenario.BS_UE_min_d;
total_BS_sector_num = scenario.total_BS_sector_num;



if isempty(Ue_pos)
    if strcmp(scenario.name,'UrbanGrid')
        % Pedstrain UT
        Pedstrain_UE_pos = network_layout.drop_Pedstrain_UT(scenario,1,[0,0]);
        for uk = 1:16
            Ue      = elements.Equipment(scenario,'UE');
            id      = id + 1;
            Ue.ID   = id;
            Ue.fcin = (BW/20)*(floor(20*rand)-(19/2)) + fc;
            Ue.d_2D_in = zeros(1);
            Ue.height     = 1.5;
            Ue_pos_2D      = Pedstrain_UE_pos(uk,:);
            Ue_pos_3D = [Ue_pos_2D, Ue.height];
            Ue.inital_Position  = Ue_pos_3D;
            Ue.Position = Ue_pos_3D;
            Ue.rand_LoS = rand(1,total_BS_sector_num);

            Ue.n_fl     = 1;
            Ue.O2Isigma = 0;
            Ue.carPL = [0,0]; % for metallized car windows, 米 = 20 can be used.
            ang.alpha   = rand*360-180;
            ang.beta = Ue.antenna_params.beta;
            ang.gamma = 0;
            Ue.sector.antenna  = antennas.antenna_array(Ue.antenna_params, ang);
            Ue.sector.antenna.attachedDevice = Ue;
            Ue.sector.antenna.attachedType   = 'UE-Pedstrain';
            ue_pos_list = [ue_pos_list; Ue.Position]; %#ok<AGROW>
            UE_list     = [UE_list;Ue]; %#ok<AGROW>
        end

        % RSU UT
        RSU_pos = network_layout.drop_RSU_UT(scenario);
        for uk = 1:(sqrt(numel(BS_list)/2)+1)^2
            Ue      = elements.Equipment(scenario,'UE');
            id      = id + 1;
            Ue.ID   = id;
            Ue.fcin = (BW/20)*(floor(20*rand)-(19/2)) + fc;
            Ue.d_2D_in = zeros(1);
            Ue_pos_2D     = RSU_pos(uk,:);
            Ue.height     = 5;
            Ue_pos_3D = [Ue_pos_2D, Ue.height];
            Ue.inital_Position  = Ue_pos_3D;
            Ue.Position = Ue_pos_3D;
            Ue.rand_LoS = rand(1,total_BS_sector_num);

            Ue.n_fl     = 1;
            Ue.O2Isigma = 0;
            Ue.carPL = [0,0]; % for metallized car windows, 米 = 20 can be used.
            
            ang.alpha   = rand*360-180;
            ang.beta = Ue.antenna_params.beta;
            ang.gamma = 0;
            Ue.sector.antenna  = antennas.antenna_array(Ue.antenna_params, ang);
            Ue.sector.antenna.attachedDevice = Ue;
            Ue.sector.antenna.attachedType   = 'UE-RSU';
            ue_pos_list = [ue_pos_list; Ue.Position]; %#ok<AGROW>
            UE_list     = [UE_list;Ue]; %#ok<AGROW>
        end


        % vehicle UT
        for m = 1:numel(scenario.Centerposition)/2
            vehicle_list = struct('s', [], 'lane_id', [],'pos',[],'grid_id',[]);
            for uk = 1:scenario.vehicle_per_sec
                Ue      = elements.Equipment(scenario,'UE');
                id      = id + 1;
                Ue.ID   = id;
                Ue.fcin = (BW/20)*(floor(20*rand)-(19/2)) + fc;
                Ue.d_2D_in = zeros(1);
                Ue.velocity = 60*3600/1000; % 60km/hr
                [pos,side,lane_id, s]      = network_layout.drop_Vehicle_UT(scenario, vehicle_list, Ue.velocity, scenario.Centerposition(m,:),[],false);
                if isempty(pos)   
                    break;   
                end
                vehicle_list.lane_id(uk,1) = lane_id;
                vehicle_list.s(uk,1) = s;
                vehicle_list.pos(uk,:) = pos;
                vehicle_list.grid_id(uk,1) = m;
                Ue.height     = 1.6;
                Ue_pos_3D = [pos, Ue.height];
                Ue.phi_v = (side-1)*90;
                Ue.inital_Position  = Ue_pos_3D;
                Ue.Position = Ue_pos_3D;
                Ue.rand_LoS = rand(1,total_BS_sector_num);
                Ue.n_fl     = 1;
                Ue.O2Isigma = 0;
                Ue.carPL = [0,0]; % for metallized car windows, 米 = 20 can be used.

                ang.alpha   = rand*360-180;
                ang.beta = Ue.antenna_params.beta;
                ang.gamma = 0;
                Ue.sector.antenna  = antennas.antenna_array(Ue.antenna_params, ang);
                Ue.sector.antenna.attachedDevice = Ue;
                Ue.sector.antenna.attachedType   = 'UE-Vehicle';
                ue_pos_list = [ue_pos_list; Ue.Position]; %#ok<AGROW>
                UE_list     = [UE_list;Ue]; %#ok<AGROW>
            end
            if isequal(scenario.Centerposition(m,:), [0 0])
                vehicle_list_center = vehicle_list;
            end
        end

    else

        for m = 1:numel(BS_list)
            BS = BS_list(m);
            for s = 1 : BS.sector_num
                sector = BS.sector(s);
                for uk = 1:scenario.UE_per_sec
                    Ue      = elements.Equipment(scenario,'UE');
                    id      = id + 1;
                    Ue.ID   = id;
                    if ismember(scenario.name,{'RMa','UMa','UMi'})
                        R = scenario.R;
                        Ue.fcin = (BW/20)*(floor(20*rand)-(19/2)) + fc;
                        if rand < 0.8  || indoor
                            % indoor
                            Ue.Indoor  = true;
                            Ue.d_2D_in = scenario.UE_max_d_2D_indoor*rand(1);
                            Ue_pos      = network_layout.drop_in_hexagonUE(BS.Position(1:2), R, min_d+Ue.d_2D_in, sector.boresight,BS.sector_num);
                            N_fl        = randi([4, 8],1);
                            Ue.n_fl     = randi([1, N_fl],1);
                            Ue.O2IPL    = lowhigh{round(rand)+1};
                        else
                            % outdoor
                            Ue.d_2D_in = zeros(1);
                            Ue_pos      = network_layout.drop_in_hexagonUE(BS.Position(1:2), R, min_d, sector.boresight,BS.sector_num);
                            Ue.n_fl     = 1;
                        end

                        Ue.O2Isigma = randn;
                        Ue.carPL = [9,5];               % for metallized car windows, m = 20 can be used.

                    elseif ismember(scenario.name,{'InH','InF'})
                        Ue.fcin     = fc;
                        Ue_pos      = network_layout.drop_InH_ISAC([scenario.x_range;scenario.y_range]);
                    end
                    Ue.rand_LoS = rand(1,total_BS_sector_num);
                    Ue.height     = 1.5;
                    Ue.Position = [Ue_pos, Ue.height];
                    Ue.inital_Position  = [Ue_pos, Ue.height];

                    ang.alpha   = rand*360-180;
                    ang.beta = Ue.antenna_params.beta;
                    ang.gamma = 0;

                    Ue.sector.antenna  = antennas.antenna_array(Ue.antenna_params, ang);
                    Ue.sector.antenna.attachedDevice = Ue;
                    Ue.sector.antenna.attachedType   = 'UE';
                    ue_pos_list = [ue_pos_list; Ue.Position]; %#ok<AGROW>
                    UE_list     = [UE_list;Ue]; %#ok<AGROW>
                end
            end
        end
    end
else
    Ue_pos = localPosition3D(Ue_pos, scenario.UE_height, 'custom_ue_position');
    for uk = 1:size(Ue_pos,1)
        Ue      = elements.Equipment(scenario,'UE');
        id      = id + 1;
        Ue.ID   = id;
        Ue.fcin = (BW/20)*(floor(20*rand)-(19/2)) + fc;
        if rand < 0.8
            Ue.Indoor  = true;
            Ue.d_2D_in = 25*rand(1);
        else
            Ue.d_2D_in = zeros(1);
            Ue.n_fl     = 1;
        end
        Ue.Position = Ue_pos(uk,:);
        Ue.O2Isigma = randn;
        Ue.carPL = [9,5]; % for metallized car windows, 米 = 20 can be used.
        Ue.rand_LoS = rand(1,total_BS_sector_num);
        Ue.height     = Ue_pos(uk,3);
        Ue.inital_Position  = Ue.Position;
        ang.alpha   = rand*360-180;
        ang.beta = Ue.antenna_params.beta;
        ang.gamma = 0;
        Ue.sector.antenna  = antennas.antenna_array(Ue.antenna_params, ang);
        Ue.sector.antenna.attachedDevice = Ue;
        Ue.sector.antenna.attachedType   = 'UE';
        ue_pos_list = [ue_pos_list; Ue.Position]; %#ok<AGROW>
        UE_list     = [UE_list;Ue]; %#ok<AGROW>
    end

end
if plot_controller
    figure(1);
    for i = 1:total_BS_sector_num
        indx = i:total_BS_sector_num:numel(UE_list);
        plot3(ue_pos_list(indx,1),ue_pos_list(indx,2),ue_pos_list(indx,3),'.','color','b', 'HandleVisibility','off');hold on;
    end
end
axis equal;view(0,90);
end

function [plot_controller, custom_ue_position] = localParseInputs(varargin)
plot_controller = false;
custom_ue_position = [];

if nargin >= 1 && ~isempty(varargin{1})
    if islogical(varargin{1}) || (isnumeric(varargin{1}) && isscalar(varargin{1}))
        plot_controller = logical(varargin{1});
    else
        custom_ue_position = varargin{1};
    end
end

if nargin >= 2
    custom_ue_position = varargin{2};
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

