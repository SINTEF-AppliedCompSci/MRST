classdef SequentialPressureTransportModelDG < SequentialPressureTransportModel
    
    methods
        % ----------------------------------------------------------------%
        function model = SequentialPressureTransportModelDG(pressureModel, transportModel, varargin)
            model = model@SequentialPressureTransportModel(pressureModel, transportModel, varargin{:});
        end
        
        % ----------------------------------------------------------------%
        function state = validateState(model, state)
            state = model.pressureModel.validateState(state);
            state = assignDofFromState(model.transportModel.disc, state);
        end
    end
    
end