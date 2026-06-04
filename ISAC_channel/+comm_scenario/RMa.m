classdef RMa < comm_scenario.Comm_Scenario
    properties
        W                   % building width
        h                   % building height
        layer_num
    end
    properties(Dependent)
        R   % radius
        BSposition
    end

    methods
        function obj = RMa
            obj.name       = 'RMa';
            obj.W          = 20;
            obj.h          = 5;
            obj.ISD        = 1732;           % 1732 or 5000 m
            obj.BS_height  = 35;             % 25 m for 3D-UMa; 10 m for 3D-UMi.
            obj.BS_Tx_power   = 46;          % dBm
            obj.BS_min_d      = 35;          % the min 2D distance between BS and UE, 35 m for 3D-UMa; 10 m for 3D-UMi.
            obj.layer_num  = 3;              % the number of layer of cells.
            obj.x_range    = [-ceil(obj.R*(1.5*obj.layer_num-0.5)), ceil(obj.R*(1.5*obj.layer_num-0.5))];
            obj.y_range    = [-ceil(obj.ISD*(obj.layer_num-0.5)), ceil(obj.ISD*(obj.layer_num-0.5))];

            obj.Cross_correlation_LOS = chol([    1,     0, -0.5,    0,    0,      0,   -0.8;          % SF/K/DS/ASD/ASA/ZSD/ZSA
                0,     1,    0,    0,    0,     0,      0;
                -0.5,     0,    1,    0,    0,     0,      0;
                0,     0,    0,    1,    0,   0.5,      0;
                0,     0,    0,    0,    1,     0,      0;
                0,     0,    0,  0.5,    0,     1,      0;
                -0.8,     0,    0,    0,    0,     0,      1],'lower');
            obj.Cross_correlation_NLOS = chol([    1,  -0.5,   0.6,    0,       0,   -0.4;     % SF/DS/ASD/ASA/ZSD/ZSA
                -0.5,     1,  -0.4,    0,   -0.5,      0;
                0.6,  -0.4,     1,    0,    0.5,   -0.1;
                0,     0,     0,    1,      0,      0;
                0,  -0.5,   0.5,    0,      1,      0;
                -0.4,     0,  -0.1,    0,      0,      1],'lower');
            obj.Cross_correlation_O2I = chol([1, 0,    0,    0,     0,     0;     % SF/DS/ASD/ASA/ZSD/ZSA
                0, 1,    0,    0,     0,     0;
                0, 0,    1, -0.7,  0.66,  0.47;
                0, 0, -0.7,    1, -0.55, -0.22;
                0, 0, 0.66,-0.55,     1,     0;
                0, 0, 0.47,-0.22,     0,    1],'lower');


        end

        function value = get.R(obj)
            value = obj.ISD/sqrt(3);
        end

        function BS_pos_list = get.BSposition(obj)
            BS_pos_list = [0 0];
            for layer = 2:obj.layer_num
                vertex = obj.ISD*(layer-1).*exp(1j*[30 90 150 210 270 330 30]./180*pi);
                for n = 1:6
                    d = (vertex(n+1)-vertex(n))/(layer-1);
                    for m = 1:layer-1
                        pos = vertex(n)+d*m;
                        BS_pos_list = [BS_pos_list;[real(pos), imag(pos)]]; %#ok<AGROW>
                    end
                end
            end
            BS_pos_list(abs(BS_pos_list)<10^-7) = 0;
        end

        function plot_BS_pos(obj)
            figure(1); hold off;
            plot3(obj.BS_pos_list(:,1),obj.BS_pos_list(:,2),obj.BS_pos_list(:,3),'ro','markersize',5,'linewidth',2, 'HandleVisibility','off');hold on;
            around_umi = scenario.R*[1,0;0.5,sqrt(3)/2;-0.5,sqrt(3)/2;-1,0;-0.5,-sqrt(3)/2;0.5,-sqrt(3)/2;1,0];   % umi sector
            % 3 sector
            around_3_sector = scenario.R/sqrt(3)*[sqrt(3)/2,-0.5;0,0;0,1;sqrt(3)/2,1.5;sqrt(3),1;sqrt(3),0;sqrt(3)/2,-0.5;sqrt(3)/2,-1.5;0,-2;-sqrt(3)/2,-1.5;-sqrt(3)/2,-0.5;-sqrt(3),0;-sqrt(3),1;-sqrt(3)/2,1.5;0,1;0,0;-sqrt(3)/2,-0.5];
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
