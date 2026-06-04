classdef UMa < comm_scenario.Comm_Scenario
    properties
        layer_num
    end
    properties(Dependent)
        R   % radius
        BSposition
        
    end
    
    methods
        function obj = UMa
            obj.name       = 'UMa';
            obj.ISD        = 500;            % [200-500] m for 3D-UMa; 200 m for 3D-UMi.
            obj.BS_height  = 25;             % 25 m for 3D-UMa; 10 m for 3D-UMi.
            obj.BS_Tx_power   = 56;          % dBm
            obj.BS_UE_min_d      = 35;          % the min 2D distance between BS and UE, 35 m for 3D-UMa
            obj.BS_ST_min_d      = 35;          % the min 2D distance between BS and UE, 35 m for 3D-UMa
            obj.layer_num  = 3;              % the number of layer of cells.
            obj.x_range    = [-ceil(obj.R*(1.5*obj.layer_num-0.5)), ceil(obj.R*(1.5*obj.layer_num-0.5))];
            obj.y_range    = [-ceil(obj.ISD*(obj.layer_num-0.5)), ceil(obj.ISD*(obj.layer_num-0.5))];
            obj.UE_max_d_2D_indoor = 25;  

            obj.Cross_correlation_LOS = chol([    1.0000    0.0000   -0.4000   -0.5000   -0.5000    0.0000   -0.8000;
                                    0.0000    1.0000   -0.4000    0.0000   -0.2000    0.0000    0.0000;
                                   -0.4000   -0.4000    1.0000    0.4000    0.8000   -0.2000    0.0000;
                                   -0.5000    0.0000    0.4000    1.0000    0.0000    0.5000    0.0000;
                                   -0.5000   -0.2000    0.8000    0.0000    1.0000   -0.3000    0.4000;
                                    0.0000    0.0000   -0.2000    0.5000   -0.3000    1.0000    0.0000;
                                   -0.8000    0.0000    0.0000    0.0000    0.4000    0.0000    1.0000],'lower');
            obj.Cross_correlation_NLOS = chol([   1.0000   -0.4000   -0.6000    0.0000    0.0000   -0.4000;
                                   -0.4000    1.0000    0.4000    0.6000   -0.5000    0.0000;
                                   -0.6000    0.4000    1.0000    0.4000    0.5000   -0.1000;
                                    0.0000    0.6000    0.4000    1.0000    0.0000    0.0000;
                                    0.0000   -0.5000    0.5000    0.0000    1.0000    0.0000;
                                   -0.4000    0.0000   -0.1000    0.0000    0.0000    1.0000],'lower');
            obj.Cross_correlation_O2I = chol([    1.0000   -0.5000    0.0000    0.5300    0.0000    0.4000;
                                   -0.5000    1.0000    0.4000    0.0000   0.0000   -0.5300;
                                    0.0000    0.4000    1.0000    0.0000   0.0000    0.4200;
                                    0.5300    0.0000    0.0000    1.0000    0.0000    0.0000;
                                    0.0000    0.0000    0.0000    0.0000    1.0000    0.0000;
                                    0.4000   -0.5300    0.4200    0.0000    0.0000    1.0000],'lower');

            % RP
            obj.a_d = 10.3370;
            obj.b_d = 0.1317;
            obj.c_d = 68.7778;
            obj.a_h = 16.2253;
            obj.b_h = 1.9218;
            obj.c_h = 2.6142;  
        end
        
        function value = get.R(obj)
            value = obj.ISD/sqrt(3);
        end

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
            plot3(BS_pos_list(:,1),BS_pos_list(:,2),BS_pos_list(:,3),'ro','markersize',5,'linewidth',2, 'HandleVisibility','off');hold on;
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
    end
end

