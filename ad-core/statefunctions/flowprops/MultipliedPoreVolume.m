classdef MultipliedPoreVolume < StateFunction
    % Effective pore-volume after pressure-multiplier
    properties
    end
    
    methods
        function gp = MultipliedPoreVolume(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn({'pressure'}, 'state');
            end
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for possible
            % rock-compressibility
            f = model.fluid;
            pv = model.operators.pv;
            if isfield(f, 'pvMultR')
                p = model.getProp(state, 'pressure');
                pvMult = prop.evaluateFunctionOnDomainWithArguments(f.pvMultR, p);
                pv = pv.*pvMult;
            end
        end
    end
end