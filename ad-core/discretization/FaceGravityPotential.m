classdef FaceGravityPotential < GridProperty
    properties

    end
    
    methods
        function gRhoDz = evaluateOnGrid(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            
            gRhoDz = cell(1, nph);
            if norm(model.gravity) > 0
                gdz = model.getGravityGradient();
                rho = model.getProp(state, 'Density');
                for i = 1:nph
                    rhof = model.operators.faceAvg(rho{i});
                    gRhoDz{i} = - rhof.*gdz;
                end
            end
        end

    end
end