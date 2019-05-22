classdef GravityPotentialDifference < AutoDiffFunction
    properties
        saturationWeighting = false;
    end
    
    methods
        function gp = GravityPotentialDifference(varargin)
            gp@AutoDiffFunction(varargin{:});
            gp = gp.dependsOn('Density', 'FlowPropertyFunctions');
            if gp.saturationWeighting
                gp = gp.dependsOn('s', 'state');
            end
        end
        function gRhoDz = evaluateOnDomain(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            
            gRhoDz = cell(1, nph);
            gdz = model.getGravityGradient();
            if norm(model.gravity) > 0
                nm = model.getPhaseNames();
                rho = model.getProp(state, 'Density');
                avg = model.operators.faceAvg;
                for i = 1:nph
                    if prop.saturationWeighting
                        s = model.getProp(state, ['s', nm(i)]);
                        rhof = avg(s.*rho{i})./max(avg(s), 1e-8);
                    else
                        rhof = avg(rho{i});
                    end
                    gRhoDz{i} = - rhof.*gdz;
                end
            else
                [gRhoDz{:}] = deal(gdz);
            end
        end
    end
end