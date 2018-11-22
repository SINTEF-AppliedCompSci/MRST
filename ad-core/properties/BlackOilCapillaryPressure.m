classdef BlackOilCapillaryPressure < GridProperty
    properties
    end
    
    methods
        function pc = evaluateOnGrid(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            pc = cell(1, nph);
            
            f = model.fluid;
            if model.water && model.oil && isfield(f, 'pcOW')
                sW = model.getProp(state, 'sw');
                % Note sign! Water is always first
                pc{phInd == 1} = -prop.evaluateFunctionOnGrid(f.pcOW, sW);
            end
            
            if model.gas && model.oil && isfield(f, 'pcOG')
                sG = model.getProp(state, 'sg');
                pc{phInd == 3} = prop.evaluateFunctionOnGrid(f.pcOG, sG);
            end
        end
    end
end