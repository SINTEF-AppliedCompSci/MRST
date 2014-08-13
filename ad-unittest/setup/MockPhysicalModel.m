classdef MockPhysicalModel < PhysicalModel
    
    methods
        function model = MockPhysicalModel()
            model = model@PhysicalModel([]);
        end
        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case 'x1'
                    fn = 'x';
                    index = 1;
                case 'x2'
                    fn = 'x';
                    index = 2;
                case 'x'
                    fn = 'x';
                    index = ':';
                case 'y'
                    fn = 'y';
                    index = 1;
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@PhysicalModel(model, name);
            end
        end
    end
end