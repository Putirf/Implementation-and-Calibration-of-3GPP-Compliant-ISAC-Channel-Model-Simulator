classdef Sector < handle
    %BS_SECTOR 
    
    properties
        ID           % the ID of sector [equipment id, local ID, global ID]
        equipment_type 
        frequency    % the RF in Hz
        BW           % the bandwidth
        boresight    % the direction of the sector
        antenna      % the configurational antenna
        attached_equipment  % the attached BaseStation of current sector
        PHI_n_m
        link_params
    end
    % properties(Dependent)
    %     num_UE       % the number of attached UE
    % end
    % 
    methods
        function obj = Sector(varargin)
            if length(varargin) >= 1
                obj.attached_equipment = varargin{1};
            end
            if length(varargin) == 2
                obj.frequency = varargin{2};
            else 
                obj.frequency = obj.attached_equipment.frequency;
                obj.BW = obj.attached_equipment.BW;
            end
        end
        
        % function plot_links(obj)
        %     pos = zeros(6,obj.num_UE);
        %     pos(1:2:5,:) = [obj.attached_UE(:).pos3D];
        %     pos(2:2:6,:) = [reshape([obj.attached_UE(:).attachBS_pos],2,[]); obj.attached_equipment.h_BS*ones(1,obj.num_UE)];
        %     figure(1); hold on;
        %     plot3(pos(1:2,:),pos(3:4,:),pos(5:6,:),'color',obj.lineColor); hold on;
        % end
        % 
        % %% get and set function
        % function out = get.num_UE(obj)
        %     out = length(obj.attached_UE);
        % end
    end
end

