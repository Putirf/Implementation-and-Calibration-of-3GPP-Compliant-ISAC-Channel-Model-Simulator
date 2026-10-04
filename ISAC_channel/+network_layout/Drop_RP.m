function [RP_list, RP_pos_list] = Drop_RP(attached_equipment,scenario, varargin)
RP_list = []; RP_pos_list = [];
id      = 0;

rp_angle = rand*360-180;
[a_d, b_d, c_d, a_h, b_h, c_h, rp_mode] = ...
    localRpParameters(attached_equipment, scenario);

for rp = 1:scenario.RP_per_equipment
    RP      = elements.Equipment(scenario,'RP');
    id      = id + 1;
    RP.ID   = id;
    RP.d_2D_in = 0;
    RP_pos_2D      = network_layout.drop_in_hexagonRP(attached_equipment.Position,rp,rp_angle,a_d,b_d,c_d);
    RP.height     = gamrnd(a_h, 1/b_h,1,1) + c_h;
    RP_pos_3D   = [RP_pos_2D, RP.height];
    RP.inital_Position  = RP_pos_3D;
    RP.Position = RP_pos_3D;
    RP.rand_LoS = 1;
    RP.RP = struct('mode', rp_mode, 'c_d', c_d, 'c_h', c_h);

    RP.antenna_params  = attached_equipment.antenna_params;
    RP.velocity = attached_equipment.velocity;
    RP.theta_v = attached_equipment.theta_v;
    RP.phi_v = attached_equipment.phi_v;

    RP.antenna_params.attachedDevice = RP;
    RP.antenna_params.attachedType   = 'RP';
    RP_pos_list = [RP_pos_list; RP_pos_3D]; %#ok<AGROW>
    RP_list     = [RP_list;RP]; %#ok<AGROW>
end
% if ~isempty(varargin) && varargin{1}
%     figure(1);
%     for i = 1:numel(BSsector_list)
%         indx = i:numel(BSsector_list):numel(RP_list);
%         plot3(rp_pos_list(indx,1),rp_pos_list(indx,2),rp_pos_list(indx,3),'+','color',BSsector_list(i).lineColor,'MarkerSize', 5,'LineWidth', 1.5, 'HandleVisibility','off');hold on;
%     end
% end
% axis equal;view(0,90);
end



function [a_d, b_d, c_d, a_h, b_h, c_h, rp_mode] = ...
        localRpParameters(attached_equipment, scenario)
if strcmp(attached_equipment.type, 'UE')
    % TR 38.901 Table 7.9.4.2-2 Part 1: UT-monostatic background.
    rp_mode = 'UE_mono';
    switch scenario.name
        case 'UMi'
            values = [10.0220, 1.2522, 11.0040, 3.0487, 1.9128, 0.1785];
        case {'UMa', 'UrbanGrid'}
            values = [2.9072, 0.1031, 3.8471, 1.6640, 1.6215, -1.4205];
        case 'RMa'
            values = [10.2421, 0.0526, 3.3131, 0.3175, 1.4150, 1.5906];
        case 'InH'
            values = [4.3733, 0.4457, 4.6302, 0.2974, 0.4103, 2.9711];
        case 'InF'
            values = [0.231418, 0.128133, 2.004903, ...
                0.462968, 0.281526, -16.921515];
        otherwise
            error('DropRP:UnsupportedUeMonostaticScenario', ...
                'UE_mono RP parameters are not available for scenario "%s".', ...
                scenario.name);
    end
else
    % Existing TRP-monostatic configuration remains unchanged.
    rp_mode = 'TRP_mono';
    values = [scenario.a_d, scenario.b_d, scenario.c_d, ...
        scenario.a_h, scenario.b_h, scenario.c_h];
end

a_d = values(1);
b_d = values(2);
c_d = values(3);
a_h = values(4);
b_h = values(5);
c_h = values(6);
end


