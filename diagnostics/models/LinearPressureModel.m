classdef LinearPressureModel < PressureModel
    properties
        useAverageMobility = true;
    end
    
    methods
        function model = LinearPressureModel(parent, varargin)
            model = model@PressureModel(parent);
            model.useIncTol = false;
            model.pressureTol = 1e-4;
            model = merge_options(model, varargin{:});
        end
        
        function model = validateModel(model, forces, varargin)
            model = validateModel@PressureModel(model, forces, varargin{:});
            if model.useAverageMobility
                fd = model.parentModel.FlowDiscretization;
                fd = fd.setStateFunction('FaceComponentMobility', FaceComponentMobilityAverage(model.parentModel, []));
                fd = fd.setStateFunction('FaceMobility', FaceMobilityAverage(model.parentModel, []));
                model.parentModel.FlowDiscretization = fd;
            end
        end
        
        function [state, report] = updateAfterConvergence(model, varargin)
            [state, report] = updateAfterConvergence@PressureModel(model, varargin{:});
            mob = model.parentModel.getProp(state, 'mobility');
            state.mob = horzcat(mob{:});
            if isfield(state, 'statePressure')
                state = rmfield(state, 'statePressure');
            end
        end
        
        function [eqs, names, types, state, src] = getModelEquations(pmodel, state0, state, dt, drivingForces)
            model = pmodel.parentModel;
            ncomp = model.getNumberOfComponents();
            [eqs, names, types, state] = model.getModelEquations(state0, state, dt, drivingForces);
            ceqs = eqs(1:ncomp);
            
            w = model.getProp(state, 'PressureReductionFactors');
            % Assemble equations and add in sources
            pressure_equation = 0;
            for i = 1:numel(ceqs)
                pressure_equation = pressure_equation + w{i}.*ceqs{i};
            end
            subs = (ncomp+1):numel(eqs);
            eqs = [{pressure_equation}, eqs(subs)];
            names = ['pressure', names(subs)];
            types = ['cell', types(subs)];
            src = model.FacilityModel.getComponentSources(state);
        end
        
    end
end