classdef MultipliedPoreVolumeDG < MultipliedPoreVolume
    % Effective pore-volume after pressure-multiplier
    properties
    end
    
    methods
        function gp = MultipliedPoreVolumeDG(model, varargin)
            gp@MultipliedPoreVolume(model, varargin{:});
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for possible
            % rock-compressibility
            f = model.fluid;
            pv = model.operators.pv;
            pv = pv(state.cells);
            if isfield(f, 'pvMultR')
                p = model.getProp(state, 'pressure');
                pvMult = prop.evaluateFunctionOnDomainWithArguments(f.pvMultR, p);
                pvMult = pvMult(state.cells);
                pv = pv.*pvMult;
            end
        end
    end
end