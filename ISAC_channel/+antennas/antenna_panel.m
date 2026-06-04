classdef antenna_panel < handle
    properties
        ID
        M
        N
        Kv
        Kh
        d_H
        d_V
        P
        X_pol
        ele_downtilt
        ele_panning

        attachedArray
        element_list
    end

    properties(Dependent)
        alpha
        beta
        gamma
        attachedType
        num_element
        pos_LCS
        pos_element_LCS
        w_m
        w_n
        num_port
    end

    methods
        function obj = antenna_panel(varargin)
            if nargin == 0
                obj.loadDefaultPanel();
            else
                obj.loadPanelConfig(varargin{1});
            end
            obj.element_list = obj.createElements(varargin);
        end

        function gain = panel_gain(obj, link)
            lambda = 3e8 / link.sector.frequency;
            txPos = [link.BS_pos_wrap, link.sector.h_BS]';
            rxPos = [link.UE.pos, link.UE.h_UT]';
            txRay = obj.losUnitVector(link.theta_LOS_ZOD, link.phi_LOS_AOD);

            phaseLos = 2 * pi / lambda * sqrt(sum((rxPos - txPos).^2));
            rxPanel = obj.firstPanel(link.UE.antenna.panel);
            rxField = rxPanel.element_list(1).field_pattern(link.phi_LOS_AOA, link.theta_LOS_ZOA);

            gainLos = zeros(1, obj.N);
            losMatrix = [exp(-1j*phaseLos), 0; 0, -exp(-1j*phaseLos)];
            for col = 1:obj.N
                txField = obj.weightedTxField(link, lambda, txRay, col);
                gainLos(col) = sum(abs(rxField.' * losMatrix * txField).^2);
            end

            gain = 10 * log10(gainLos);
        end

        function set_ele_downtilt(obj, new_ele_downtilt)
            obj.ele_downtilt = new_ele_downtilt;
        end

        function [indr, indc] = get_port_pos(obj, port_indx)
            obj.assertIndexInRange(port_indx, obj.num_port, 'port');
            [baseRow, baseCol] = obj.gridIndex(mod(port_indx, obj.M*obj.N/(obj.Kv*obj.Kh)), obj.M/obj.Kv);
            indr = (baseRow*obj.Kv + 1):(baseRow + 1)*obj.Kv;
            indc = (baseCol*obj.Kh + 1):(baseCol + 1)*obj.Kh;
        end

        function [indr, indc] = get_element_pos(obj, element_indx)
            obj.assertIndexInRange(element_indx, obj.num_element, 'element');
            [baseRow, baseCol] = obj.gridIndex(mod(element_indx, obj.M*obj.N), obj.M);
            indr = baseRow + 1;
            indc = baseCol + 1;
        end

        function out = get.alpha(obj)
            out = obj.attachedArray.alpha;
        end

        function out = get.beta(obj)
            out = obj.attachedArray.beta;
        end

        function out = get.gamma(obj)
            out = obj.attachedArray.gamma;
        end

        function out = get.attachedType(obj)
            out = obj.attachedArray.attachedType;
        end

        function out = get.num_element(obj)
            out = obj.M * obj.N * obj.P;
        end

        function out = get.w_m(obj)
            weights = obj.steeringWeights(obj.Kv, obj.d_V, cosd(obj.ele_downtilt), true);
            out = obj.expandWeights(weights, 1, obj.Kh);
        end

        function out = get.w_n(obj)
            weights = obj.steeringWeights(obj.Kh, obj.d_H, sind(obj.ele_panning), false);
            out = obj.expandWeights(weights, obj.Kv, 1);
        end

        function out = get.num_port(obj)
            out = obj.M * obj.N / (obj.Kv * obj.Kh) * obj.P;
        end

        function out = get.pos_element_LCS(obj)
            out = obj.elementPositions(obj.M, obj.N, obj.d_H, obj.d_V);
        end

        function out = get.pos_LCS(obj)
            out = squeeze(obj.attachedArray.pos_panel_LCS(obj.ID(1), obj.ID(2), :));
        end
    end

    methods(Access = private)
        function loadDefaultPanel(obj)
            obj.M = 2;
            obj.N = 2;
            obj.Kv = 1;
            obj.Kh = 1;
            obj.d_H = 0.5;
            obj.d_V = 0.5;
            obj.P = 1;
            obj.X_pol = 0;
            obj.ele_downtilt = 102;
            obj.ele_panning = 0;
        end

        function loadPanelConfig(obj, cfg)
            fields = {'M', 'N', 'Kv', 'Kh', 'd_H', 'd_V', 'P', 'X_pol', 'ele_downtilt', 'ele_panning'};
            for idx = 1:numel(fields)
                obj.(fields{idx}) = cfg.(fields{idx});
            end
        end

        function elements = createElements(obj, args)
            cfg = [];
            if ~isempty(args)
                cfg = args{1};
            end

            if obj.P == 1
                obj.X_pol = 0;
            elseif obj.P ~= 2
                error('Only one or two polarizations are supported.');
            end

            elements(1, obj.P) = antennas.antenna_element();
            for idx = 1:obj.P
                elements(idx) = antennas.antenna_element(obj.X_pol(idx));
                elements(idx).attachedPanel = obj;
                if ~isempty(cfg) && isfield(cfg, 'pol_model')
                    elements(idx).pol_model = cfg.pol_model;
                end
            end
        end

        function field = weightedTxField(obj, link, lambda, txRay, col)
            rowList = 1:obj.Kv;
            elemPos = permute(obj.pos_element_LCS(rowList, col, :), [3, 1, 2]);
            txOffset = obj.attachedArray.R * (obj.pos_LCS + elemPos) * lambda;
            baseField = obj.element_list(1).field_pattern(link.phi_LOS_AOD, link.theta_LOS_ZOD);
            field = (exp(1j * 2*pi/lambda * (txRay' * txOffset)) * obj.w_m) * baseField;
        end

        function panelObj = firstPanel(~, panelGrid)
            if iscell(panelGrid)
                panelObj = panelGrid{1};
            else
                panelObj = panelGrid(1);
            end
        end

        function assertIndexInRange(~, idx, maxCount, label)
            if idx < 0 || idx >= maxCount
                error('the %s %d does not exist.', label, idx);
            end
        end
    end

    methods(Static, Access = private)
        function vec = losUnitVector(theta, phi)
            vec = [sind(theta)*cosd(phi); sind(theta)*sind(phi); cosd(theta)];
        end

        function [row, col] = gridIndex(linearIdx, rowsPerColumn)
            row = mod(linearIdx, rowsPerColumn);
            col = floor(linearIdx / rowsPerColumn);
        end

        function out = steeringWeights(count, spacing, trigValue, asColumn)
            idx = 0:count-1;
            if asColumn
                idx = idx.';
            end
            out = exp(-1j * 2*pi * idx * spacing * trigValue) / sqrt(count);
        end

        function out = expandWeights(weights, repeatRows, repeatCols)
            out = weights;
            if repeatRows ~= 1 || repeatCols ~= 1
                out = repmat(weights, repeatRows, repeatCols);
            end
        end

        function pos = elementPositions(numRow, numCol, spacingH, spacingV)
            [rowIdx, colIdx] = meshgrid(0:numRow-1, 0:numCol-1);
            pos = zeros(numRow, numCol, 3);
            pos(:, :, 2) = (colIdx' - (numCol - 1) / 2) * spacingH;
            pos(:, :, 3) = (rowIdx' - (numRow - 1) / 2) * spacingV;
        end
    end
end
