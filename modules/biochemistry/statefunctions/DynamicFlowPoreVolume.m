classdef DynamicFlowPoreVolume < PoreVolume
    % Effective pore-volume after pressure-multiplier
    properties
    end
    
    methods
        function gp = DynamicFlowPoreVolume(model, varargin)
            gp@PoreVolume(model, varargin{:});
            if model.bacteriamodel   
                gp = gp.dependsOn('nbact', 'state');
            end
            gp = gp.dependsOn('pressure', 'state');
            assert(isfield(model.fluid, 'pvMultR'), 'pvMultR missing from fluid.');
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for rock-compressibility
            pv = evaluateOnDomain@PoreVolume(prop, model, state);
            p = model.getProp(state, 'pressure');
            if model.bacteriamodel
                nbact = model.getProp(state, 'nbact');
                pvMult = prop.evaluateFluid(model, 'pvMultR',p, nbact);
                pv = pv.*pvMult;
            end
        end
    end
end