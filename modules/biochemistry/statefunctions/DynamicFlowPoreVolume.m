classdef DynamicFlowPoreVolume < PoreVolume
    % Effective pore-volume after pressure-multiplier
    properties
    end
    
    methods
        function gp = DynamicFlowPoreVolume(model, varargin)
            gp@PoreVolume(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'nbact'}, 'state');
            assert(isfield(model.fluid, 'pvMultR'), 'pvMultR missing from fluid.');
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for rock-compressibility
            pv = evaluateOnDomain@PoreVolume(prop, model, state);
            p = model.getProp(state, 'pressure');

            nbact = model.getProp(state, 'nbact');

            pvMult = prop.evaluateFluid(model, 'pvMultR',p, nbact);
            pv = pv.*pvMult;
        end
    end
end