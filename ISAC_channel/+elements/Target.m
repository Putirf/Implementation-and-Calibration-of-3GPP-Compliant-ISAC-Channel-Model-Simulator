classdef Target < handle
    properties
        ID
        type
        sub_type
        rand_LoS
        inital_Position
        Position
        height        
        velocity               % 3 km/h
        theta_v
        phi_v        
        RCS
        
        % link_list
        % best_link_id
        % couplingloss_TRP_mono_LS
        % couplingloss_TRP_TRP_LS
        % couplingloss_TRP_UE_LS
        % couplingloss_UE_UE_LS
        % couplingloss_UE_mono_LS
        % couplingloss_TRP_mono_Full
        % couplingloss_TRP_TRP_Full
        % couplingloss_TRP_UE_Full
        % couplingloss_UE_UE_Full
        % couplingloss_UE_mono_Full
        % DS_TRP_mono
        % SA4_TRP_mono

        is_single_STSP              % boolen
        SP
    end
    properties(Dependent)
        pos3D
        STSP_pos
        STSP_velocity
    end

    methods
        function obj = Target(varargin)
            if isempty(varargin)
                obj.velocity       = 3/3.6;   % 3 km/h
                obj.theta_v = 90;
                obj.phi_v   = rand*360-180;
                obj.height    = 1.5;
            else
                params = varargin{1};
                if isstruct(params)
                    names = fieldnames(params);
                    for name_idx = 1:numel(names)
                        obj.(names{name_idx}) = params.(names{name_idx});
                    end
                end
            end
        end

        function value = get.pos3D(obj)
            value = [obj.Position';obj.height];
        end
        
        function value = get.STSP_pos(obj)
            value = obj.Position;
        end

        function value = get.STSP_velocity(obj)
            value = obj.velocity * [cosd(obj.phi_v), sind(obj.phi_v), cosd(obj.theta_v)];
        end
    
    end
end
