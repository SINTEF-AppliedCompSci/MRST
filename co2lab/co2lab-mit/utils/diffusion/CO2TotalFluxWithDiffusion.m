classdef CO2TotalFluxWithDiffusion < ComponentTotalFlux
    properties
        D
        componentDiffusion
        faceAverage = false;
    end
    
    methods
        function cf = CO2TotalFluxWithDiffusion(model)
            cf@ComponentTotalFlux(model);
            ncomp = model.getNumberOfComponents();
            cf.D = getFaceDiffusivity(model.G, model.rock);
            cf.componentDiffusion = ones(1, ncomp);
        end

        function v = evaluateOnDomain(prop, model, state)
            op = model.operators;
            T = prop.D(op.internalConn);
            faceAvg = prop.faceAverage;
            
            % Get mass fraction of CO2
            cellMass = model.getProps(state, 'ComponentPhaseMass');
            waterInBrine = cellMass{1, 1};
            Co2InBrine   = cellMass{2, 1};
            massBrine    = waterInBrine + Co2InBrine;
            if isa(waterInBrine, 'GenericAD')
                X_co2_val = Co2InBrine.val ./ massBrine.val;
                X_co2_val = min(X_co2_val, 0.999);
                X_co2 = double2GenericAD(X_co2_val, waterInBrine);
            else
                X_co2 = Co2InBrine ./ massBrine;
                X_co2 = min(X_co2, 0.999);
            end
            
            % Add diffuse flux
            v = evaluateOnDomain@ComponentTotalFlux(prop, model, state);
            C = prop.componentDiffusion;
            cellDens = model.getProps(state, 'Density');
            s = state.s;
            c = 2;
            grad_xi = op.Grad(X_co2);
            if faceAvg
                faceDens = op.faceAvg(cellDens{1});
                faceS = op.faceAvg(s{1});
            else
                flag = value(grad_xi) < 0;
                faceDens = op.faceUpstr(flag, cellDens{1});
                faceS = op.faceUpstr(flag, s{1});
            end
            Tc = C(c).*T;
            diff_flux = -faceDens.*faceS.*Tc.*grad_xi;  % kg/s/m^2
            v{2} = v{2} + diff_flux;
         
        end
    end
end