classdef CompositionalViscosityLV < GridProperty
    properties
    end
    
    methods
        function gp = CompositionalViscosityLV(model, varargin)
            gp@GridProperty(model, varargin{:});
            gp = gp.dependsOn({'PhasePressures', 'PhaseCompressibilityFactors'});
        end
        
        function mu = evaluateOnDomain(prop, model, state)
            [act, phInd] = model.getActivePhases();
            nph = sum(act);
            mu = cell(1, nph);
            
            [p, T, x, y] = model.getProps(state, 'pressure', 'T', 'x', 'y');
            Z = prop.getEvaluatedDependencies(state, 'PhaseCompressibilityFactors');
            if model.water
                p_phase = prop.getEvaluatedDependencies(state, 'PhasePressures');
                f = model.fluid;
                wix = phInd == 1;
                pw = p_phase{wix};
                mu{wix} = prop.evaluateFunctionOnGrid(f.muW, pw);
            end
            pm = model.EOSModel.PropertyModel;
            oix = phInd == 2;
            mu{oix} = pm.computeViscosity(p, x, Z{oix}, T, true);
            gix = phInd == 3;
            mu{gix} = pm.computeViscosity(p, y, Z{oix}, T, false);
        end
    end
end