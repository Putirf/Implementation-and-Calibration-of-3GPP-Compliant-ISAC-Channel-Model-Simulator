classdef InH < comm_scenario.Comm_Scenario
    properties
        BS_num
    end

    properties(Dependent)
        BSposition
    end
    
    methods
        function obj = InH(varargin)
            obj.name    = 'InH';
            if numel(varargin)>0
                obj.subname    = varargin;
            else
                obj.subname    = 'open_office';    % open_office , mixed_office
            end
            obj.x_range = [-60, 60];
            obj.y_range = [-25, 25];
            obj.BS_UE_min_d    = 0;
            obj.BS_ST_min_d    = 0;
            obj.BS_height   = 3;
            obj.BS_Tx_power = 24;    % dB
            obj.BS_num    = [6,2];
            obj.ISD         = 20;


            obj.a_d = 4.236;
            obj.b_d = 0.19255;
            obj.c_d = 4.99;
            obj.a_h = 1.3293;
            obj.b_h = 0.1442;
            obj.c_h = -13.19 ;


            obj.Cross_correlation_LOS = chol([    1,    0.5, -0.8, -0.4, -0.5,   0.2,   -0.1;          % SF/K/DS/ASD/ASA/ZSD/ZSA
                                  0.5,     1, -0.5,    0,    0,     0,    0.1;
                                 -0.8,  -0.5,    1,  0.6,  0.8,   0.1,    0.2;
                                 -0.4,     0,  0.6,    1,  0.4,   0.2,    0.2;
                                 -0.5,     0,  0.8,  0.4,    1,   0.1,    0.3;
                                  0.2,     0,  0.1,  0.2,  0.1,     1,    0.2;
                                 -0.1,   0.1,  0.2,  0.2,  0.3,   0.2,      1],'lower'); 
            obj.Cross_correlation_NLOS = chol([    1,  -0.5,     0,   -0.4,     0,    -0.1;     % SF/DS/ASD/ASA/ZSD/ZSA
                                  -0.5,     1,   0.4,     0,  -0.1,    -0.1;
                                     0,   0.4,     1,     0,   0.3,     0.2;
                                  -0.4,     0,     0,     1,   0.1,       0;
                                     0,  -0.1,   0.3,   0.1,     1,     0.4;
                                  -0.1,  -0.1,   0.2,     0,   0.4,      1],'lower'); 
            % if ~isempty(varargin)
            %     obj.BS_num   = varargin{1};
            %     if numel(varargin)>1
            %         obj.InH_case = varargin{2};
            %     end
            % else
            %     obj.BS_num = 2;
            % end
            % if obj.BS_num == 2
            %     obj.BSposition = [-30 0;30 0];
            % elseif obj.BS_num == 3
            %     obj.BSposition = [0 0; -40 0;40 0];
            % elseif obj.BS_num == 6 && strcmp(obj.InH_case,'A')
            %     obj.BSposition = [-10 0; 10 0; -30 0; 30 0; -50 0;50 0];
            % elseif obj.BS_num == 6 && strcmp(obj.InH_case,'B')
            %     obj.BSposition = [0 10; -40 10;40 10; 0 -10; -40 -10;40 -10];
            % elseif obj.BS_num == 12 
            %     obj.BSposition = [-10 10; 10 10; -30 10; 30 10; -50 10;50 10; -10 -10; 10 -10; -30 -10; 30 -10; -50 -10;50 -10];
            % end
        end

        function BS_pos_list = get.BSposition(obj)
            BS_pos_list = [];
            Dx_side = (sum(abs(obj.x_range)) - (obj.BS_num(1) - 1)*obj.ISD)/2;
            Dy_side = (sum(abs(obj.y_range)) - (obj.BS_num(2) - 1)*obj.ISD)/2;
            for x = 1:obj.BS_num(1)
                position_x = obj.x_range(1) + Dx_side + obj.ISD * (x-1);
                for y = 1:obj.BS_num(2)
                    position_y = obj.y_range(1) + Dy_side + obj.ISD * (y-1);
                    position = [position_x ,position_y , obj.BS_height];                    
                    BS_pos_list = [BS_pos_list;position]; %#ok<AGROW>
                end
            end
        end

        function plot_BS_pos(obj,BS_pos_list)
            figure(1); hold off;
            plot3(BS_pos_list(:,1),BS_pos_list(:,2),BS_pos_list(:,3),'ro','markersize',5,'linewidth',2, 'HandleVisibility','off');hold on;

            around = [obj.x_range(1), obj.y_range(1);obj.x_range(2), obj.y_range(1);obj.x_range(2), obj.y_range(2);obj.x_range(1), obj.y_range(2);obj.x_range(1), obj.y_range(1)];
            plot(around(:,1),around(:,2),'k','linewidth',2); grid on;
        end
        
    end
end

