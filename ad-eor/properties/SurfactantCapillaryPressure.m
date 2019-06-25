classdef SurfactantCapillaryPressure < BlackOilCapillaryPressure
% Implementation only for 2 phases oil-water
    properties
    end
    
    methods
        function prop = SurfactantCapillaryPressure(varargin)
            prop@BlackOilCapillaryPressure(varargin{:});
            prop = prop.dependsOn({'sW', 'surfactant'}, 'state');
        end
        
        function pc = evaluateOnDomain(prop, model, state)
            
            fluid = model.fluid;
            if ~isfield(fluid, 'pcOW')
                pc = evaluateOnDomain@BlackOilCapillaryPressure(prop, model, ...
                                                                state);
            else
                [act, phInd] = model.getActivePhases();
                nph = sum(act);
                pc = cell(1, nph);
                c = model.getProps(state, 'surfactant');
                sW = model.getProps(state, 'sW');
                pcow = prop.evaluateFunctionOnDomainWithArguments(fluid.pcOW, sW);
                pcow = pcow.*fluid.ift(c)/fluid.ift(0);
                % Note sign! Water is always first
                pc{phInd == 1} = -pcow;
            end
        end
    end
end