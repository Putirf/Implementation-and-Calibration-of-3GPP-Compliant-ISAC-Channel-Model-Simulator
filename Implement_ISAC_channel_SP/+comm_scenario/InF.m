classdef InF < comm_scenario.Comm_Scenario
    % INH_3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        % SL DL SH DH HH
        r               % Clutter density
        hc              % Effective clutter height,  < Ceiling height, 0-10 m
        d_cluster       % Typical clutter size

        L               % Room length 
        W               % Room width 
        H               % Room height 
        V               % Room size  Rectangular: 20-160000 m^2
        S               % Room surface area

        BS_num

        d_subsce
        p_subsce
    end

    properties(Dependent)
        BSposition
    end
    
    methods
        function obj = InF(varargin)
            obj.name    = 'InF';
            if numel(varargin)>0
                obj.subname = varargin{1};
            else
                obj.subname = 'SH';
            end

            if strcmp(obj.subname,'SL') || strcmp(obj.subname,'DH')
                obj.L       = 120;
                obj.W       = 60;
                obj.ISD     = 20;
            else
                obj.L       = 300;
                obj.W       = 150;
                obj.ISD     = 50;
            end
            obj.x_range = [-obj.L/2, obj.L/2];
            obj.y_range = [-obj.W/2, obj.W/2];

            if strcmp(obj.subname,'SL') || strcmp(obj.subname,'SH')
                obj.r  = 0.2;               % Low clutter density (<40%)
                obj.hc = 2;
                obj.d_cluster = 10;
                obj.H = 10;
            else
                obj.r  = 0.6;               % High clutter density (≥40%) 
                obj.hc = 6;
                obj.d_cluster = 2;
                obj.H = 15;
            end

            if strcmp(obj.subname,'SL') || strcmp(obj.subname,'DL')
                obj.BS_height  = 1.5;  % 8
            else
                obj.BS_height  = 8;  % 8
            end
            if strcmp(obj.subname,'DH')
                obj.d_subsce = 1;
                obj.p_subsce = 0.6;
            else
                obj.d_subsce = 0;
                obj.p_subsce = 1;
            end


            obj.V           = obj.L*obj.W*obj.H;
            obj.S           = 2*(obj.L*obj.W+obj.L*obj.H + obj.W*obj.H);
            obj.BS_UE_min_d    = 0;
            obj.BS_ST_min_d    = 0;
            obj.BS_Tx_power = 24;    % dB
            obj.BS_num      = [6,3];


            % RP
            obj.a_d = 0.039836;
            obj.b_d = 0.179783;
            obj.c_d = 1.13002 ;
            obj.a_h = 0.283447;
            obj.b_h = 0.435965;
            obj.c_h = -17.043530;


            obj.Cross_correlation_LOS = chol([             1,    0,    0,    0,    0,    0,    0;          % SF/K/DS/ASD/ASA/ZSD/ZSA   
                                            0,    1, -0.7, -0.5,    0,    0,    0;   
                                            0, -0.7,    1,    0,    0,    0,    0;  
                                            0, -0.5,    0,    1,    0,    0,    0;   
                                            0,    0,    0,    0,    1,    0,    0;   
                                            0,    0,    0,    0,    0,    1,    0;   
                                            0,    0,    0,    0,    0,    0,    1],'lower'); 
            obj.Cross_correlation_NLOS = chol([    1, 0, 0, 0, 0, 0;   
                                0, 1, 0, 0, 0, 0;   
                                0, 0, 1, 0, 0, 0;   
                                0, 0, 0, 1, 0, 0;   
                                0, 0, 0, 0, 1, 0;   
                                0, 0, 0, 0, 0, 1],'lower'); 
        end

        function BS_pos_list = get.BSposition(obj)
            BS_pos_list = [];
            Dx_side = (sum(abs(obj.x_range)) - (obj.BS_num(1) - 1)*obj.ISD)/2;
            Dy_side = (sum(abs(obj.y_range)) - (obj.BS_num(2) - 1)*obj.ISD)/2;
            for x = 1: obj.BS_num(1)
                position_x = obj.x_range(1) + Dx_side + obj.ISD * (x-1);
                for y = 1: obj.BS_num(2)
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

