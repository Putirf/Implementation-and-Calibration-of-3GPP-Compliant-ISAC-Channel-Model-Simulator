classdef antenna_array < handle
    properties
        Mg
        Ng
        dg_H
        dg_V

        alpha
        beta
        gamma

        attachedType
        attachedDevice

        panel
    end

    properties(Dependent)
        num_panel
        pos_panel_LCS
        R
    end

    methods
        function obj = antenna_array(varargin)
            if nargin == 0
                obj.useDefaultLayout();
                return;
            end

            cfg = varargin{1};
            ang = varargin{2};
            obj.applyArrayConfig(cfg.array, ang);
            obj.panel = obj.makePanelGrid(cfg.panel);
        end

        function gain = array_gain(obj, link)
            gain = cell(obj.Mg, obj.Ng);
            for m = 1:obj.Mg
                for n = 1:obj.Ng
                    panelObj = obj.getPanel(m, n);
                    gain{m, n} = panelObj.panel_gain(link);
                end
            end
        end

        function set_num_panel(obj, mn)
            if ~(isequal(size(mn), [1, 2]) || isequal(size(mn), [2, 1]))
                error('Please input the parameters with the size of [2 x 1] or [1 x 2]');
            end
            obj.Mg = mn(1);
            obj.Ng = mn(2);
        end

        function out = get.num_panel(obj)
            out = obj.Mg * obj.Ng;
        end

        function out = get.pos_panel_LCS(obj)
            out = obj.gridPositions(obj.Mg, obj.Ng, obj.dg_H, obj.dg_V);
        end

        function out = get.R(obj)
            out = obj.eulerRotation(obj.alpha, obj.beta, obj.gamma).';
        end
    end

    methods(Access = private)
        function useDefaultLayout(obj)
            obj.Mg = 1;
            obj.Ng = 1;
            obj.dg_H = 2.5;
            obj.dg_V = 2.5;
            obj.alpha = 0;
            obj.beta = 0;
            obj.gamma = 0;

            obj.panel = antennas.antenna_panel();
            obj.panel.attachedArray = obj;
            obj.panel.ID = [1, 1];
        end

        function applyArrayConfig(obj, arrayCfg, angleCfg)
            obj.Mg = arrayCfg.Mg;
            obj.Ng = arrayCfg.Ng;
            obj.dg_H = arrayCfg.dg_H;
            obj.dg_V = arrayCfg.dg_V;

            obj.alpha = angleCfg.alpha;
            obj.beta = angleCfg.beta;
            obj.gamma = angleCfg.gamma;
        end

        function grid = makePanelGrid(obj, panelCfg)
            grid = cell(obj.Mg, obj.Ng);
            for row = 1:obj.Mg
                for col = 1:obj.Ng
                    grid{row, col} = obj.newPanel(panelCfg, row, col);
                end
            end
        end

        function panelObj = newPanel(obj, panelCfg, row, col)
            panelObj = antennas.antenna_panel(panelCfg);
            panelObj.attachedArray = obj;
            panelObj.ID = [row, col];
        end

        function panelObj = getPanel(obj, row, col)
            if iscell(obj.panel)
                panelObj = obj.panel{row, col};
            else
                panelObj = obj.panel(row, col);
            end
        end
    end

    methods(Static, Access = private)
        function pos = gridPositions(numRow, numCol, spacingH, spacingV)
            [rowIdx, colIdx] = meshgrid(0:numRow-1, 0:numCol-1);
            pos = zeros(numRow, numCol, 3);
            pos(:, :, 2) = (colIdx' - (numCol - 1) / 2) * spacingH;
            pos(:, :, 3) = (rowIdx' - (numRow - 1) / 2) * spacingV;
        end

        function rot = eulerRotation(alpha, beta, gamma)
            ca = cosd(alpha); sa = sind(alpha);
            cb = cosd(beta);  sb = sind(beta);
            cg = cosd(gamma); sg = sind(gamma);

            rot = [ca*cb,                  sa*cb,                  -sb;
                   ca*sb*sg - sa*cg,       sa*sb*sg + ca*cg,       cb*sg;
                   ca*sb*cg + sa*sg,       sa*sb*cg - ca*sg,       cb*cg];
        end
    end
end
