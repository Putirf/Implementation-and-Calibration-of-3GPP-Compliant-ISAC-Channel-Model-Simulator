classdef UMi < comm_scenario.Comm_Scenario
    properties
        layer_num
    end
    properties(Dependent)
        R   % radius
        BSposition
    end
    
    methods
        function obj = UMi
            obj.name       = 'UMi';
            obj.ISD        = 200;            % 500 m (option: 200m) for 3D-UMa; 200 m for 3D-UMi.
            obj.BS_height  = 10;             % 25 m for 3D-UMa; 10 m for 3D-UMi.
            obj.BS_Tx_power   = 56;          % dBm
            obj.BS_UE_min_d      = 10;          % the min 2D distance between BS and UE, 35 m for 3D-UMa; 10 m for 3D-UMi.
            obj.BS_ST_min_d      = 10;          % the min 2D distance between BS and UE, 35 m for 3D-UMa; 10 m for 3D-UMi.
            obj.layer_num  = 3;              % the number of layer of cells.
            obj.x_range    = [-ceil(obj.R*(1.5*obj.layer_num-0.5)), ceil(obj.R*(1.5*obj.layer_num-0.5))];
            obj.y_range    = [-ceil(obj.ISD*(obj.layer_num-0.5)), ceil(obj.ISD*(obj.layer_num-0.5))];
            obj.UE_max_d_2D_indoor = 25;  

            obj.Cross_correlation_LOS = chol([ 1  , 0.5, -0.4, -0.5, -0.4,   0,  0;
                                 0.5,   1, -0.7, -0.2, -0.3,   0,  0; % SF/K/DS/ASD/ASA/ZSD/ZSA
                                -0.4,-0.7,    1,  0.5,  0.8,   0,0.2;
                                -0.5,-0.2,  0.5,    1,  0.4, 0.5,0.3;
                                -0.4,-0.3,  0.8,  0.4,    1,   0,  0;
                                   0,   0,    0,  0.5,    0,   1,  0;
                                   0,   0,  0.2,  0.3,    0,   0,  1],'lower');
            obj.Cross_correlation_NLOS = chol([   1,  -0.7,  0,  -0.4,   0,   0;
                                 -0.7,     1,  0,   0.4,-0.5,   0;
                                    0,     0,  1,     0, 0.5, 0.5;
                                 -0.4,   0.4,  0,     1,   0, 0.2;
                                    0,  -0.5,0.5,     0,   1,   0;
                                    0,     0,0.5,   0.2,   0,   1],'lower');% SF/DS/ASD/ASA/ZSD/ZSA
            obj.Cross_correlation_O2I = chol([    1.0000   -0.5000    0.0000    0.5300    0.0000    0.4000;
                                   -0.5000    1.0000    0.4000    0.0000   0.0000   -0.5300;
                                    0.0000    0.4000    1.0000    0.0000   -0.2000    0.0000;
                                    0.5300    0.0000    0.0000    1.0000    0.0000    0.4200;
                                    0.0000    0.0000   -0.2000    0.0000    1.0000    0.0000;
                                    0.4000   -0.5300    0.0000    0.4200    0.0000    1.0000],'lower');

            % RP
            obj.a_d = 6.1996;
            obj.b_d = 0.1558;
            obj.c_d = 15.2697;
            obj.a_h = 12.0487;
            obj.b_h = 2.3261;
            obj.c_h = 0.0157;
        end
        
        function value = get.R(obj)
            value = obj.ISD/sqrt(3);
        end

        % function BS_pos_list = get.BSposition(obj)
        %     if ~isempty(obj.BSposition_override)
        %         BS_pos_list = obj.BSposition_override;
        %         return;
        %     end
        %     BS_pos_list = [0 0];
        %     for layer = 2:obj.layer_num
        %         vertex = obj.ISD*(layer-1).*exp(1j*[30 90 150 210 270 330 30]./180*pi);
        %         for n = 1:6
        %             d = (vertex(n+1)-vertex(n))/(layer-1);
        %             for m = 1:layer-1
        %                 pos = vertex(n)+d*m;
        %                 BS_pos_list = [BS_pos_list;[real(pos), imag(pos)]];
        %             end
        %         end
        %     end
        %     BS_pos_list(abs(BS_pos_list)<10^-7) = 0;
        % end
        function BS_pos_list = get.BSposition(obj)
            BS_pos_list = [0 0 obj.BS_height];
            for layer = 2:obj.layer_num
                vertex = obj.ISD*(layer-1).*exp(1j*[30 90 150 210 270 330 30]./180*pi);
                for n = 1:6
                    d = (vertex(n+1)-vertex(n))/(layer-1);
                    for m = 1:layer-1
                        pos = vertex(n)+d*m;
                        BS_pos_list = [BS_pos_list;[real(pos), imag(pos),obj.BS_height]]; %#ok<AGROW>
                    end
                end
            end
            BS_pos_list(abs(BS_pos_list)<10^-7) = 0;
        end

        function plot_BS_pos(obj,BS_pos_list)
            figure(1); hold off;
            % The calibration scene draws the BS as a black star.  Do not
            % draw an additional red circle at the same coordinates.
            hold on;
            around_umi = obj.R*[1,0;0.5,sqrt(3)/2;-0.5,sqrt(3)/2;-1,0;-0.5,-sqrt(3)/2;0.5,-sqrt(3)/2;1,0];   % umi sector
            % 3 sector
            around_3_sector = obj.R/sqrt(3)*[sqrt(3)/2,-0.5;0,0;0,1;sqrt(3)/2,1.5;sqrt(3),1;sqrt(3),0;sqrt(3)/2,-0.5;sqrt(3)/2,-1.5;0,-2;-sqrt(3)/2,-1.5;-sqrt(3)/2,-0.5;-sqrt(3),0;-sqrt(3),1;-sqrt(3)/2,1.5;0,1;0,0;-sqrt(3)/2,-0.5]; 
            if obj.BS_sec_num == 3
                around = around_3_sector;
            else
                around = around_umi;
            end
            for i = 1:size(BS_pos_list,1)
                bspos = BS_pos_list(i,1:2) + around;
                plot(bspos(:,1),bspos(:,2),'k', 'HandleVisibility','off');hold on;
            end
        end

        
        % function setCustomBSposition(obj, BS_pos)
        %     obj.BSposition_override = BS_pos;
        % end
        % 
        % function clearCustomBSposition(obj)
        %     obj.BSposition_override = [];
        % end
    end
end
