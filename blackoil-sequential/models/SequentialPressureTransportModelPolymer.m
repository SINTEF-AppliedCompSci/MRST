classdef SequentialPressureTransportModelPolymer < SequentialPressureTransportModel
    properties
        polymer
    end
    
    methods
        function model = SequentialPressureTransportModelPolymer(...
                pressureModel, transportModel, varargin)
            
            % Call parent
            model = model@SequentialPressureTransportModel(...
                pressureModel, transportModel, varargin{:});
            
            if isprop(model.transportModel, 'polymer')
               model.polymer = model.transportModel.polymer;
            end
        end
    end
end
