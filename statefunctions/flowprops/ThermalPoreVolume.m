classdef ThermalPoreVolume < StateFunction
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function gp = ThermalPoreVolume(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn('Temperature');
                gp = gp.dependsOn({'pressure'}, 'state');
            end
            gp.label = '\phi';
        end
        
        %-----------------------------------------------------------------%
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for possible
            % rock-compressibility-dilatation as a function of p and T
            f = model.fluid;
            pv = model.operators.pv;
            if isfield(f, 'pvMultR')
                [p, T] = model.getProps(state, 'pressure', 'Temperature');
                pvMult = prop.evaluateFluid(model, 'pvMultR', p, T);
                pv     = pv.*pvMult;
            end
        end
    end
end