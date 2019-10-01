classdef CompositionalViscosityLV < StateFunction
    properties
        useCompactEvaluation = true;
    end
    
    methods
        function gp = CompositionalViscosityLV(model, varargin)
            gp@StateFunction(model, varargin{:});
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
                mu{wix} = prop.evaluateFunctionOnDomainWithArguments(f.muW, pw);
            end
            pm = model.EOSModel.PropertyModel;
            oix = phInd == 2;
            mu{oix} = pm.computeViscosity(p, x, Z{oix}, T, true);
            gix = phInd == 3;
            
            twoPhase = model.getTwoPhaseFlag(state);
            if prop.useCompactEvaluation && ~all(twoPhase)
                muV = mu{oix};
                if any(twoPhase)
                    if iscell(y)
                        y = cellfun(@(x) x(twoPhase), y, 'uniformoutput', false);
                    else
                        y = y(twoPhase, :);
                    end
                    muV(twoPhase) = pm.computeViscosity(p(twoPhase), y, Z{gix}(twoPhase), T(twoPhase), false);
                end
            else
                muV = pm.computeViscosity(p, y, Z{gix}, T, false);
            end
            mu{gix} = muV;
        end
    end
end