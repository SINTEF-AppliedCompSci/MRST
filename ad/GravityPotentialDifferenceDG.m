classdef GravityPotentialDifferenceDG < GravityPotentialDifference
    
    methods
        function gp = GravityPotentialDifferenceDG(varargin)
            gp = gp@GravityPotentialDifference(varargin{:});
        end
        function gRhoDz = evaluateOnDomain(prop, model, state)
            act = model.getActivePhases();
            nph = sum(act);
            
            gRhoDz = cell(1, nph);
            
            nf = model.G.faces.num;
            gdz = zeros(nf,1);
            gdz(model.operators.internalConn) = model.getGravityGradient();
%             [W, x, cellNo, faceNo] = model.disc.getCubature((1:model.G.cells.num)', 'surface');
            [W, x, cellNo, faceNo] = model.disc.getCubature(find(model.operators.internalConn), 'face');
            gdz = gdz(faceNo);
            
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