classdef UrbanGrid < comm_scenario.Comm_Scenario
    properties
        layer_num
        Lanewidth
        Sidewalkwidth
        vehicle_per_sec
        user_case
        BS_layer_num
        grid_layer_num
        grid_dx
        grid_dy
    end
    properties(Dependent)
        BSposition
        Centerposition
    end

    methods
        function obj = UrbanGrid(varargin)
            if ~isempty(varargin)
                obj.ISD   = varargin{1};
                if numel(varargin)>1
                    obj.layer_num = varargin{2};
                else
                    obj.layer_num = 2;
                end
            else
                obj.ISD   = 500;
                obj.layer_num = 2;
            end

            obj.name       = 'UrbanGrid';
            % Keep the TRP layout independent of the 250 m x 433 m road
            % grid.  FR1 uses the TR 36.885/37.885 500 m macro layout;
            % FR2 uses the 18-BS layout specified by TR 38.901.
            obj.BS_layer_num = obj.layer_num;
            obj.grid_layer_num = obj.layer_num;
            obj.grid_dx = 250;
            obj.grid_dy = 433;
            obj.BS_height  = 25;
            % obj.x_range    = [-obj.ISD*(2*obj.layer_num-1)/2,obj.ISD*(2*obj.layer_num-1)/2];
            % obj.y_range    = [-(obj.ISD*433/250)*(2*obj.layer_num-1)/2,(obj.ISD*433/250)*(2*obj.layer_num-1)/2];
            obj.x_range    = [-obj.grid_dx*(2*obj.grid_layer_num-1)/2,obj.grid_dx*(2*obj.grid_layer_num-1)/2];
            obj.y_range    = [-obj.grid_dy*(2*obj.grid_layer_num-1)/2,obj.grid_dy*(2*obj.grid_layer_num-1)/2];
            obj.Lanewidth = 3.5;
            obj.Sidewalkwidth = 3;
            % Automotive Urban Grid does not use an additional fixed
            % equipment-to-target rejection distance.  Separation follows
            % the TR 37.885 road, lane, pedestrian and TRP geometry.
            obj.BS_ST_min_d = 0;
            obj.UE_ST_min_d = 0;
            obj.BS_Tx_power = 56;       % dB
            obj.vehicle_per_sec = 3;
            obj.user_case = "Pedestrian";       % Pedestrian, RSU, Vehicle

            % RP
            obj.a_d = 10.3370;
            obj.b_d = 0.1317;
            obj.c_d = 68.7778;
            obj.a_h = 16.2253;
            obj.b_h = 1.9218;
            obj.c_h = 2.6142;


            obj.Cross_correlation_LOS = chol([     1,   0.5, -0.4, -0.5, -0.4,     0,      0;          % SF/K/DS/ASD/ASA/ZSD/ZSA
                                        0.5,     1, -0.7, -0.2, -0.3,     0,      0;
                                       -0.4,  -0.7,    1,  0.5,  0.8,     0,    0.2;
                                       -0.5,  -0.2,  0.5,    1,  0.4,   0.5,    0.3;
                                       -0.4,  -0.3,  0.8,  0.4,    1,     0,      0;
                                          0,     0,    0,  0.5,    0,     1,      0;
                                          0,     0,  0.2,  0.3,    0,     0,      1],'lower'); 
            obj.Cross_correlation_NLOS = chol([     1,  -0.7,     0,   -0.4,    0,       0;          % SF/DS/ASD/ASA/ZSD/ZSA
                                        -0.7,     1,     0,    0.4, -0.5,       0;
                                           0,     0,     1,     0,   0.5,     0.5;
                                        -0.4,   0.4,     0,     1,     0,     0.2;
                                           0,  -0.5,   0.5,     0,     1,       0;
                                           0,     0,   0.5,   0.2,     0,      1],'lower'); 
        end

        function BS_pos_list = get.BSposition(obj)
            if obj.ISD == 500
                % TR 36.885 Figure A.1.3-1 / TR 37.885 baseline Urban
                % macro layout: one centre site and its six first-tier
                % neighbours.  Calibration does not model interference,
                % so no geographical-distance wrap-around replicas are
                % added here.
                dy = sqrt(3) * obj.ISD / 2;
                BS_pos_list = [ ...
                     0,             0,  obj.BS_height; ...
                    -obj.ISD,       0,  obj.BS_height; ...
                     obj.ISD,       0,  obj.BS_height; ...
                    -obj.ISD/2,    dy,  obj.BS_height; ...
                     obj.ISD/2,    dy,  obj.BS_height; ...
                    -obj.ISD/2,   -dy,  obj.BS_height; ...
                     obj.ISD/2,   -dy,  obj.BS_height];
                return;
            end

            % TR 38.901 Table 7.9.6.1-3 special Urban Grid layout for
            % ISD=250 m: two BSs per road-grid anchor, 18 BSs in total.
            L = obj.BS_layer_num;
            n = 2*L - 1;

            dx = obj.ISD;
            % dy = obj.ISD*433/250;
            % TEMP UrbanGrid calibration: BS anchor row spacing follows ISD relative to road-grid width.
            % Original temporary state: dy = 2 * obj.grid_dy;
            % Original state before split: dy = obj.ISD*433/250;
            dy = (obj.ISD / obj.grid_dx) * obj.grid_dy;

            x_list = ((-(n-1)/2):((n-1)/2)) * dx;
            y_list = ((-(n-1)/2):((n-1)/2)) * dy;

            [X, Y] = meshgrid(x_list, y_list);   % X: dx grid, Y: dy grid
            BS_pos = [X(:), Y(:)];

            % TRP locations follow the configured ISD. The centre road grid
            % used for vehicle dropping remains 250 m by 433 m.
            offset1 = [10, dy/2 - 10];
            offset2 = [-dx/2+10, -10];

            BS_pos_list = zeros(n^2*2,3);
            BS_pos_list(1:2:end,1:2)  = BS_pos + offset1;
            BS_pos_list(2:2:end,1:2)  = BS_pos + offset2;
            BS_pos_list(:,3) = obj.BS_height;
        end

        function Center_pos = get.Centerposition(obj)
            % TEMP UrbanGrid calibration: road grid centers can expand independently from BS count.
            % Original state: L = obj.layer_num;
            L = obj.grid_layer_num;
            n = 2*L - 1;

            % Original state: dx = obj.ISD;
            dx = obj.grid_dx;
            % dy = obj.ISD*433/250;
            dy = obj.grid_dy;

            x_list = ((-(n-1)/2):((n-1)/2)) * dx;
            y_list = ((-(n-1)/2):((n-1)/2)) * dy;
            [X, Y] = meshgrid(x_list, y_list);   % X: dx grid, Y: dy grid

            Center_pos = [X(:), Y(:)];
        end


        function plot_BS_pos(obj, BS_pos_list,varargin)
            % plotUrbanGrid(scn)
            % scn: scenarios.UrbanGrid_3D instance
            %
            % Optional name-value:
            %   'lanesPerDir' (default 2)  
            %   'usePatch'    (default true)

            figure(1); hold off;
            plot3(BS_pos_list(:,1),BS_pos_list(:,2),BS_pos_list(:,3),'ro','markersize',5,'linewidth',2, 'HandleVisibility','off');hold on;
            
            p = inputParser;
            addParameter(p,'lanesPerDir',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
            addParameter(p,'usePatch',true,@(x)islogical(x)&&isscalar(x));
            parse(p,varargin{:});
            lanesPerDir = p.Results.lanesPerDir;

            % TEMP UrbanGrid calibration: plot the road-grid layer, not the BS layer.
            % Original state: L = obj.layer_num;
            L  = obj.grid_layer_num;
            n  = 2*L - 1;

            % Original state: dx = obj.ISD;
            dx = obj.grid_dx;
            % dy = obj.ISD*433/250;
            dy = obj.grid_dy;

            laneW = obj.Lanewidth;
            swW   = obj.Sidewalkwidth;

            streetW = 2*lanesPerDir*laneW + 2*swW;

            blkW = dx - streetW;
            blkH = dy - streetW;
            if blkW <= 0 || blkH <= 0
                error('streetW too large: dx-streetW=%.2f, dy-streetW=%.2f', blkW, blkH);
            end

            xMin = -dx*(n-1)/2 - dx/2;
            xMax =  dx*(n-1)/2 + dx/2;
            yMin = -dy*(n-1)/2 - dy/2;
            yMax =  dy*(n-1)/2 + dy/2;

            figure(1); hold on; axis equal;
            xlim([xMin xMax]); ylim([yMin yMax]);
            grid on;

            xCenters = ((-(n-1)/2):((n-1)/2)) * dx;
            yCenters = ((-(n-1)/2):((n-1)/2)) * dy;

            for ix = 1:numel(xCenters)
                for iy = 1:numel(yCenters)
                    cx = xCenters(ix);
                    cy = yCenters(iy);

                    x0 = cx - blkW/2;  x1 = cx + blkW/2;
                    y0 = cy - blkH/2;  y1 = cy + blkH/2;

                    obj.plotRect([x0 y0 blkW blkH],'k',1);

                    bx0 = x0 - swW;  bx1 = x1 + swW;
                    by0 = y0 - swW;  by1 = y1 + swW;
                    obj.plotRect([bx0 by0 bx1-bx0 by1-by0],'k',1);
                end
            end

            for kx = 0:n
                for ky = 0:n
                    x = xMin + kx*dx;
                    y = yMin + ky*dy;
                    plot([x x],[yMin+(2+2*max(0,ky-1))*laneW+(dy-2*laneW)*max(0,ky-1) yMin+(2+4*max(0,ky-1))*laneW+(dy-4*laneW)*ky],'k','LineWidth',1); % vertical road centerline
                    plot([x+laneW x+laneW],[yMin+(2+2*max(0,ky-1))*laneW+(dy-2*laneW)*max(0,ky-1) yMin+(2+4*max(0,ky-1))*laneW+(dy-4*laneW)*ky],'k--','LineWidth',0.8); % vertical road centerline
                    plot([x-laneW x-laneW],[yMin+(2+2*max(0,ky-1))*laneW+(dy-2*laneW)*max(0,ky-1) yMin+(2+4*max(0,ky-1))*laneW+(dy-4*laneW)*ky],'k--','LineWidth',0.8); % vertical road centerline
                    plot([xMin+(2+2*max(0,kx-1))*laneW+(dx-2*laneW)*max(0,kx-1) xMin+(2+4*max(0,kx-1))*laneW+(dx-4*laneW)*kx],[y y],'k','LineWidth',1); % horizontal road centerline
                    plot([xMin+(2+2*max(0,kx-1))*laneW+(dx-2*laneW)*max(0,kx-1) xMin+(2+4*max(0,kx-1))*laneW+(dx-4*laneW)*kx],[y+laneW y+laneW],'k--','LineWidth',0.8); % horizontal road centerline
                    plot([xMin+(2+2*max(0,kx-1))*laneW+(dx-2*laneW)*max(0,kx-1) xMin+(2+4*max(0,kx-1))*laneW+(dx-4*laneW)*kx],[y-laneW y-laneW],'k--','LineWidth',0.8); % horizontal road centerline
                end
            end

            xlabel('x (m)'); ylabel('y (m)');

        end

        function plotRect(~,rect, color, lw)
            % rect = [x y w h]
            x = rect(1); y = rect(2); w = rect(3); h = rect(4);
            around = [x y; x+w y; x+w y+h; x y+h; x y];
            plot(around(:,1), around(:,2), color, 'LineWidth', lw);
        end

    end
end
