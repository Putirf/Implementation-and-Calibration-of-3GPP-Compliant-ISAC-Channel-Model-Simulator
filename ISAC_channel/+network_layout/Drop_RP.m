function [RP_list, RP_pos_list] = Drop_RP(attached_equipment,scenario, varargin)
RP_list = []; RP_pos_list = [];
id      = 0;

rp_angle = rand*360-180;
a_d = scenario.a_d;
b_d = scenario.b_d;
c_d = scenario.c_d;
a_h = scenario.a_h;
b_h = scenario.b_h;
c_h = scenario.c_h;

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



