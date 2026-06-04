classdef HighWay < comm_scenario.Comm_Scenario
    properties
        Lanewidth
        Roadlength
    end

    properties(Dependent)
        BSposition
    end

    methods
        function obj = HighWay(varargin)
            if ~isempty(varargin)
                obj.ISD   = varargin{1};
                if numel(varargin)>1
                    obj.BS_height  = varargin{2};
                    if numel(varargin)>2
                        obj.Roadlength = varargin{3};
                    else
                        obj.Roadlength = 2000;
                    end
                else
                    obj.BS_height  = 35;
                    obj.Roadlength = 2000;
                end
            else
                obj.ISD   = 1732;
                obj.BS_height  = 35;
                obj.Roadlength = 2000;
            end
            
            obj.name            = 'HighWay';
            obj.Lanewidth       = 4;
            obj.x_range         = [-obj.Roadlength/2,obj.Roadlength/2];
            obj.y_range         = [-3*obj.Lanewidth,3*obj.Lanewidth];            
        end

        function BS_pos_list = get.BSposition(obj)
            x_pos = obj.x_range(1):obj.ISD:obj.x_range(2);
            if x_pos(end) < obj.x_range(2)
                x_pos = [x_pos, obj.x_range(2)];
            end
            BS_pos_list = [x_pos(:), zeros(numel(x_pos), 1), obj.BS_height * ones(numel(x_pos), 1)];
        end 


        
    end
end
